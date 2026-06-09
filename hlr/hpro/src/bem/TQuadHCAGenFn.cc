//
// Project     : HLIBpro
// File        : TQuadHCAGenFn.cc
// Description : class providing HCA functionality for BEM bilinear forms
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <set>

#include "unordered_map.hh"

#include "hpro/bem/TGaussQuad.hh"

#include "hpro/bem/TQuadHCAGenFn.hh"

namespace Hpro
{

namespace
{

//
// local types
//

using  tri_set_t = std::set< idx_t >;
using  val_map_t = std::unordered_map< idx_t, idx_t >;

}// namespace anonymous

//!
//! constructor
//! - \a order defines (maximal) quadrature order for
//!   evaluating the integrals 
//!
template < typename  ansatzsp_t,
           typename  testsp_t,
           typename  value_t >
TQuadHCAGenFn< ansatzsp_t, testsp_t, value_t >::TQuadHCAGenFn ( const ansatzsp_t *    ansatzsp,
                                                                const testsp_t *      testsp,
                                                                const uint            quad_order,
                                                                const TPermutation *  row_perm_i2e,
                                                                const TPermutation *  col_perm_i2e,
                                                                stat_t *              stat )
        : TPermHCAGeneratorFn< value_t >( row_perm_i2e, col_perm_i2e )
        , _ansatz_sp( ansatzsp )
        , _test_sp( testsp )
        , _quad_order( quad_order )
        , _stat( stat )
{
    if (( _ansatz_sp == nullptr ) || ( _test_sp == nullptr ))
        HERROR( ERR_ARG, "(TQuadHCAGenFn)", "function space is nullptr" );

    //
    // build quadrature points for triangles based on
    // 1D Gauss rules with tensor product ansatz and
    // Duffy transformation
    //

    TGaussQuad   gauss;
    
    _quad_rules_cache.resize( quad_order+1 );
    
    for ( uint  n = 1 ; n <= quad_order; ++n )
    {
        std::vector< double >  pts_1d, wghts_1d;

        gauss.build( n, pts_1d, wghts_1d );

        // actual and padded number of quadrature points
        const size_t  npts        = n*n;
        const size_t  padded_npts = CFG::Mach::simd_padded_size< real_t >( npts );

        _quad_rules_cache[n].npts = npts;
        _quad_rules_cache[n].x.resize( padded_npts );
        _quad_rules_cache[n].y.resize( padded_npts );
        _quad_rules_cache[n].w.resize( padded_npts );

        size_t  quad_idx = 0;

        // define "real" quadrature points
        for ( uint j = 0; j < n; ++j )
        {
            for ( uint i = 0; i < n; ++i, ++quad_idx )
            {
                const T2Point  p( pts_1d[i] * ( 1.0 - pts_1d[j] ),
                                  pts_1d[i] * pts_1d[j] );

                _quad_rules_cache[n].x[ quad_idx ] = p.x();
                _quad_rules_cache[n].y[ quad_idx ] = p.y();
                _quad_rules_cache[n].w[ quad_idx ] = wghts_1d[i] * wghts_1d[j] * pts_1d[i];
            }// for
        }// for

        // fill rest with zero (only weights to avoid division by zero in coordinate computations)
        for ( ; quad_idx < padded_npts ; ++quad_idx )
        {
            _quad_rules_cache[n].x[ quad_idx ] = _quad_rules_cache[n].x[ 0 ];
            _quad_rules_cache[n].y[ quad_idx ] = _quad_rules_cache[n].y[ 0 ];
            _quad_rules_cache[n].w[ quad_idx ] = 0.0;
        }// for
    }// for
}
    
//!
//! Evaluate \f$ \int \phi_i(x) D_x \gamma(x, y_{l}) dx\f$ for i ∈ \a is
//! and points \f$ y_{l} \f$ defined by \a pts.
//! Store results in \a matrix at index (i,l).
//!
template < typename  ansatzsp_t,
           typename  testsp_t,
           typename  value_t >
void
TQuadHCAGenFn< ansatzsp_t, testsp_t, value_t >::integrate_dx_perm  ( const std::vector< idx_t > &    idxs,
                                                                     const std::vector< T3Point > &  pts,
                                                                     BLAS::Matrix< value_t > &       matrix ) const
{
    const auto  quad_rule = this->get_quad_rule();
    const auto  npts      = quad_rule->npts;
    const auto  fnspace   = this->ansatz_space();
    const auto  grid      = fnspace->grid();
    
    //
    // go over all points and indices and compute integral
    //
    
    std::vector< value_t >  quad_vals( CFG::Mach::simd_padded_size< value_t >( npts ) );
    
    for ( size_t  pts_i = 0; pts_i < pts.size(); pts_i++ )
    {
        const T3Point  y( pts[pts_i] );
        size_t         idx_i = 0;

        for ( const auto  idx : idxs )
        {
            //
            // integrate over all triangles in support
            //

            value_t  value   = value_t(0);

            for ( auto  tri_idx : fnspace->support( idx ) )
            {
                //
                // do quadrature
                //

                const TGrid::triangle_t  tri = grid->triangle( tri_idx );
                value_t                  val = value_t(0);

                this->eval_dx( tri_idx, y, *quad_rule, quad_vals );

                if ( _stat != nullptr )
                    _stat->n_eval_dx++;
                
                for ( size_t  k = 0; k < npts; ++k )
                    val += ( real_t(quad_rule->w[k]) *
                             fnspace->eval_basis_unit( idx, quad_rule->x[k], quad_rule->y[k], tri.vtx ) *
                             quad_vals[k] );
                    
                value += real_t(grid->tri_size( tri_idx )) * val;
            }// for
        
            matrix( idx_t( idx_i ), idx_t(pts_i) ) = value;

            ++idx_i;
        }// for
    }// for
}

//!
//! Evaluate \f$ \int \phi_j(y) D_y \gamma(x_{l}, y) dy\f$ for j ∈ \a is 
//! and points \f$ x_{l} \f$ defined by \a pts. 
//! Store results in \a matrix at index (j,l).
//!
template < typename  ansatzsp_t,
           typename  testsp_t,
           typename  value_t >
void
TQuadHCAGenFn< ansatzsp_t, testsp_t, value_t >::integrate_dy_perm  ( const std::vector< idx_t > &    idxs,
                                                                     const std::vector< T3Point > &  pts,
                                                                     BLAS::Matrix< value_t > &       matrix ) const
{
    const auto  quad_rule = this->get_quad_rule();
    const auto  npts      = quad_rule->npts;
    const auto  fnspace   = this->test_space();
    const auto  grid      = fnspace->grid();
    
    //
    // go over all points and indices and compute integral
    //

    std::vector< value_t >  quad_vals( CFG::Mach::simd_padded_size< value_t >( npts ) );
    
    for ( size_t  pts_i = 0; pts_i < pts.size(); pts_i++ )
    {
        const T3Point  x( pts[pts_i] );
        size_t         idx_i = 0;

        for ( const auto  idx : idxs )
        {
            //
            // integrate over all triangles in support
            //

            value_t  value   = value_t(0);

            for ( auto  tri_idx : fnspace->support( idx ) )
            {
                //
                // do quadrature on triangle
                //

                const auto  tri = grid->triangle( tri_idx );
                value_t     val = value_t(0);
                
                this->eval_dy( x, tri_idx, *quad_rule, quad_vals );

                if ( _stat != nullptr )
                    _stat->n_eval_dy++;
                
                for ( size_t  k = 0; k < npts; ++k )
                    val += ( real_t(quad_rule->w[k]) *
                             fnspace->eval_basis_unit( idx, quad_rule->x[k], quad_rule->y[k], tri.vtx ) *
                             quad_vals[k] );
                    
                value += real_t(grid->tri_size( tri_idx )) * val;
            }// for
        
            matrix( idx_t(idx_i), idx_t(pts_i) ) = value;

            ++idx_i;
        }// for
    }// for
}

//!
//! constructor
//!
template < typename  ansatzsp_t,
           typename  testsp_t,
           typename  value_t >
TInvarBasisQuadHCAGenFn< ansatzsp_t, testsp_t, value_t >::TInvarBasisQuadHCAGenFn ( const ansatzsp_t *    ansatzsp,
                                                                                    const testsp_t *      testsp,
                                                                                    const uint            quad_order,
                                                                                    const TPermutation *  row_perm_i2e,
                                                                                    const TPermutation *  col_perm_i2e,
                                                                                    stat_t *              stat )
: TQuadHCAGenFn< ansatzsp_t, testsp_t, value_t >( ansatzsp, testsp,
                                                  quad_order,
                                                  row_perm_i2e, col_perm_i2e,
                                                  stat )
{
    //
    // precompute ansatz values
    //
    
    const size_t  n_ansatz_tri_bases = this->ansatz_space()->n_unit_bases();
    
    _ansatz_val.resize( n_ansatz_tri_bases );

    for ( uint  nbasis = 0; nbasis < n_ansatz_tri_bases; ++nbasis )
    {
        _ansatz_val[ nbasis ].resize( this->_quad_order+1 );
        
        for ( uint  order = 1; order <= this->_quad_order; ++order )
        {
            const auto    rule        = this->get_quad_rule( order );
            const size_t  npts        = rule->x.size();
            const size_t  padded_npts = CFG::Mach::simd_padded_size< value_t >( npts );
            size_t        i           = 0;
            
            _ansatz_val[ nbasis ][ order ].resize( padded_npts );
            
            for ( ; i < npts; ++i )
            {
                _ansatz_val[ nbasis ][ order ][ i ] = this->ansatz_space()->eval_basis_unit( nbasis,
                                                                                             rule->x[i],
                                                                                             rule->y[i] );
            }// for

            // zero pad
            for ( ; i < padded_npts; ++i )
                _ansatz_val[ nbasis ][ order ][ i ] = ansatz_value_t(0);
        }// for
    }// for
    
    //
    // precompute test values
    //
    
    const size_t  n_test_tri_bases = this->test_space()->n_unit_bases();
    
    _test_val.resize( n_test_tri_bases );

    for ( uint  nbasis = 0; nbasis < n_test_tri_bases; ++nbasis )
    {
        _test_val[ nbasis ].resize( this->_quad_order+1 );
        
        for ( uint  order = 1; order <= this->_quad_order; ++order )
        {
            const auto    rule        = this->get_quad_rule( order );
            const size_t  npts        = rule->x.size();
            const size_t  padded_npts = CFG::Mach::simd_padded_size< value_t >( npts );
            size_t        i           = 0;
            
            _test_val[ nbasis ][ order ].resize( padded_npts );
            
            for ( ; i < npts; ++i )
            {
                _test_val[ nbasis ][ order ][ i ] = this->test_space()->eval_basis_unit( nbasis,
                                                                                         rule->x[i],
                                                                                         rule->y[i] );
            }// for

            // zero pad
            for ( ; i < padded_npts; ++i )
                _test_val[ nbasis ][ order ][ i ] = test_value_t(0);
        }// for
    }// for
    
}
    
//!
//! Evaluate \f$ \int \phi_i(x) D_x \gamma(x, y_{l}) dx\f$ for i ∈ \a is
//! and points \f$ y_{l} \f$ defined by \a pts.
//! Store results in \a matrix at index (i,l).
//!
template < typename  ansatzsp_t,
           typename  testsp_t,
           typename  value_t >
void
TInvarBasisQuadHCAGenFn< ansatzsp_t, testsp_t, value_t >::integrate_dx_perm  ( const std::vector< idx_t > &    idxs,
                                                                               const std::vector< T3Point > &  pts,
                                                                               BLAS::Matrix< value_t > &       matrix ) const
{
    //
    // store position of all indices for indirect access of <matrix>
    //
    
    const size_t  nidx = idxs.size();
    val_map_t     is_map;

    for ( size_t  i = 0; i < nidx; ++i )
        is_map[ idxs[i] ] = idx_t( i );
    
    //
    // collect all triangles in support of basis functions
    //

    const ansatzsp_t *  fnspace = this->ansatz_space();
    const TGrid *       grid    = fnspace->grid();
    tri_set_t           triangles;

    for ( const auto  idx : idxs )
    {
        for ( auto  tri_idx : fnspace->support( idx ) )
            triangles.insert( tri_idx );
    }// for

    //
    // loop over each triangle in function space and integrate function
    //

    const auto              quad_rule = this->get_quad_rule();
    const size_t            npts      = quad_rule->npts;
    std::vector< value_t >  quad_vals( CFG::Mach::simd_padded_size< value_t >( npts ) );
    
    for ( size_t  pts_i = 0; pts_i < pts.size(); pts_i++ )
    {
        const T3Point  y( pts[pts_i] );
            
        for ( const auto  tri_idx : triangles )
        {
            const TGrid::triangle_t  tri      = grid->triangle( tri_idx );
            const auto               tri_size = real_t( grid->tri_size( tri_idx ) );
            
            // evaluate dx-integral
            this->eval_dx( tri_idx, y, *quad_rule, quad_vals );

            if ( this->_stat != nullptr )
                this->_stat->n_eval_dx++;
                
            //
            // compute final result for each basis function by combining
            // quadrature values with precomputed basis function values
            //

            for ( auto  idx : fnspace->triangle_indices( tri_idx ) )
            {
                // exclude indices not in <is>
                if ( is_map.find( idx ) == is_map.end() )
                    continue;

                const idx_t                            idx_i      = is_map[ idx ];
                const std::vector< ansatz_value_t > *  basis_vals = ansatz_val( idx, tri, this->_quad_order );

                //
                // combine quadrature values and basis function values
                //
                
                value_t  val = value_t(0);

                for ( size_t  k = 0; k < npts; ++k )
                    val += ( real_t(quad_rule->w[k]) * (*basis_vals)[k] * quad_vals[k] );
                    
                matrix( idx_t(idx_i), idx_t(pts_i) ) += tri_size * val;
            }// for
        }// for
    }// for
}

//!
//! Evaluate \f$ \int \phi_j(y) D_y \gamma(x_{l}, y) dy\f$ for j ∈ \a is 
//! and points \f$ x_{l} \f$ defined by \a pts. 
//! Store results in \a matrix at index (j,l).
//!
template < typename  ansatzsp_t,
           typename  testsp_t,
           typename  value_t >
void
TInvarBasisQuadHCAGenFn< ansatzsp_t, testsp_t, value_t >::integrate_dy_perm  ( const std::vector< idx_t > &    idxs,
                                                                               const std::vector< T3Point > &  pts,
                                                                               BLAS::Matrix< value_t > &       matrix ) const
{
    //
    // store position of all indices for indirect access of <matrix>
    //
    
    const size_t  nidx = idxs.size();
    val_map_t     is_map;

    for ( size_t  i = 0; i < nidx; ++i )
        is_map[ idxs[i] ] = idx_t( i );
    
    //
    // collect all triangles in support of basis functions
    //

    const testsp_t *  fnspace = this->test_space();
    const TGrid *     grid    = fnspace->grid();
    tri_set_t         triangles;

    for ( const auto  idx : idxs )
    {
        for ( auto  tri_idx : fnspace->support( idx ) )
            triangles.insert( tri_idx );
    }// for

    //
    // loop over each triangle in function space and integrate function
    //

    const auto              quad_rule = this->get_quad_rule();
    const size_t            npts      = quad_rule->npts;
    std::vector< value_t >  quad_vals( CFG::Mach::simd_padded_size< value_t >( npts ) );
    
    for ( size_t  pts_i = 0; pts_i < pts.size(); pts_i++ )
    {
        const T3Point  x( pts[pts_i] );
            
        for (const auto  tri_idx : triangles )
        {
            const TGrid::triangle_t  tri      = grid->triangle( tri_idx );
            const auto               tri_size = real_t( grid->tri_size( tri_idx ) );
            
            // evaluate dx-integral
            this->eval_dy( x, tri_idx, *quad_rule, quad_vals );

            if ( this->_stat != nullptr )
                this->_stat->n_eval_dy++;
                
            //
            // compute final result for each basis function by combining
            // quadrature values with precomputed basis function values
            //

            for ( auto  idx : fnspace->triangle_indices( tri_idx ) )
            {
                // exclude indices not in <is>
                if ( is_map.find( idx ) == is_map.end() )
                    continue;

                const idx_t                          idx_i      = is_map[ idx ];
                const std::vector< test_value_t > *  basis_vals = test_val( idx, tri, this->_quad_order );

                //
                // combine quadrature values and basis function values
                //
                
                value_t  val = value_t(0);

                for ( size_t  k = 0; k < npts; ++k )
                    val += ( real_t(quad_rule->w[k]) * (*basis_vals)[k] * quad_vals[k] );
                    
                matrix( idx_t(idx_i), idx_t(pts_i) ) += tri_size * val;
            }// for
        }// for
    }// for
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

#define INST_ALL( type1, type2 )                                              \
    template class TQuadHCAGenFn< TConstFnSpace< type1 >,  TConstFnSpace< type1 >,  type2 >; \
    template class TQuadHCAGenFn< TConstFnSpace< type1 >,  TLinearFnSpace< type1 >, type2 >; \
    template class TQuadHCAGenFn< TLinearFnSpace< type1 >, TConstFnSpace< type1 >,  type2 >; \
    template class TQuadHCAGenFn< TLinearFnSpace< type1 >, TLinearFnSpace< type1 >, type2 >; \
                                                                        \
    template class TInvarBasisQuadHCAGenFn< TConstFnSpace< type1 >,  TConstFnSpace< type1 >,  type2 >; \
    template class TInvarBasisQuadHCAGenFn< TConstFnSpace< type1 >,  TLinearFnSpace< type1 >, type2 >; \
    template class TInvarBasisQuadHCAGenFn< TLinearFnSpace< type1 >, TConstFnSpace< type1 >,  type2 >; \
    template class TInvarBasisQuadHCAGenFn< TLinearFnSpace< type1 >, TLinearFnSpace< type1 >, type2 >;

INST_ALL( float,  float )
INST_ALL( double, double )
INST_ALL( float,  std::complex< float > )
INST_ALL( double, std::complex< double > )

}// namespace Hpro
