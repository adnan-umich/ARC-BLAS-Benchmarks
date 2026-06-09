//
// Project     : HLIBpro
// File        : TBEMBF.cc
// Description : classes for bilinearforms in BEM-applications
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <set>

#include "unordered_map.hh"

#include "hpro/base/error.hh"
#include "hpro/base/config.hh"

#include "hpro/bem/TQuadBEMBF.hh"

namespace Hpro
{

namespace B = BLAS;

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// TInvarBasisQuadBEMBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template < typename T_ansatzsp, typename T_testsp, typename T_value >
TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, T_value >::TInvarBasisQuadBEMBF ( const T_ansatzsp *  aansatzsp,
                                                                              const T_testsp *    atestsp,
                                                                              const uint          aorder,
                                                                              const bool          dist_ada,
                                                                              const real_t        ada_err )
        : TQuadBEMBF< T_ansatzsp, T_testsp, T_value >( aansatzsp, atestsp, aorder, dist_ada, ada_err )
{
    //
    // compute ansatz and test function values at quadrature points for
    // all local basis functions in unit triangle
    //

    compute_basis_func();
}

//////////////////////////////////////
//
// evaluate basis functions
//

//
// compute ansatz and test basis functions for all quadrature points
//
template < typename T_ansatzsp, typename T_testsp, typename T_value >
void
TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, T_value >::compute_basis_func ()
{
    const size_t  n_ansatz_tri_bases = this->ansatz_space()->n_unit_bases();
    
    _ansatz_val.resize( n_ansatz_tri_bases );

    for ( uint  nbasis = 0; nbasis < n_ansatz_tri_bases; ++nbasis )
    {
        _ansatz_val[ nbasis ].resize( 4 );
        
        for ( uint  ncommon = 0; ncommon < 4; ++ncommon )
        {
            _ansatz_val[ nbasis ][ ncommon ].resize( this->_quad_order+1 );
            
            for ( uint  order = 1; order <= this->_quad_order; ++order )
            {
                const auto    rule        = this->quad_rule( ncommon, order );
                const size_t  npts        = rule->npts;
                const size_t  padded_npts = CFG::Mach::simd_padded_size< ansatz_value_t >( npts );
                size_t        i           = 0;
                
                _ansatz_val[ nbasis ][ ncommon ][ order ].resize( padded_npts );

                for ( ; i < npts; ++i )
                {
                    const ansatz_value_t  val = this->ansatz_space()->eval_basis_unit( nbasis,
                                                                                       rule->pts1[i].x(),
                                                                                       rule->pts1[i].y() );
                    
                    _ansatz_val[ nbasis ][ ncommon ][ order ][ i ] = val;
                }// for

                // zero padding
                for ( ; i < padded_npts; ++i )
                    _ansatz_val[ nbasis ][ ncommon ][ order ][ i ] = ansatz_value_t(0);
            }// for
        }// for
    }// for
    
    const size_t  n_test_tri_bases = this->test_space()->n_unit_bases();

    _test_val.resize( n_test_tri_bases );
    
    for ( uint  nbasis = 0; nbasis < n_test_tri_bases; ++nbasis )
    {
        _test_val[ nbasis ].resize( 4 );
        
        for ( uint  ncommon = 0; ncommon < 4; ++ncommon )
        {
            _test_val[ nbasis ][ ncommon ].resize( this->_quad_order+1 );
            
            for ( uint  order = 1; order <= this->_quad_order; ++order )
            {
                const auto    rule        = this->quad_rule( ncommon, order );
                const size_t  npts        = rule->npts;
                const size_t  padded_npts = CFG::Mach::simd_padded_size< test_value_t >( npts );
                size_t        i           = 0;
                
                _test_val[ nbasis ][ ncommon ][ order ].resize( padded_npts );
                
                for ( ; i < npts; ++i )
                {
                    const test_value_t  val = this->test_space()->eval_basis_unit( nbasis,
                                                                                   rule->pts2[i].x(),
                                                                                   rule->pts2[i].y() );
                    
                    _test_val[ nbasis ][ ncommon ][ order ][ i ] = val;
                }// for

                // zero padding
                for ( ; i < padded_npts; ++i )
                    _test_val[ nbasis ][ ncommon ][ order ][ i ] = test_value_t(0);
            }// for
        }// for
    }// for
}
    
//////////////////////////////////////
//
// evaluate bilinearform
//

//
// evaluate subblock defined by \a row_ind × \a col_ind; the indices
// in \a row_ind and \a col_ind can be arbitrary, e.g. must not be
// contiguous
//
template < typename ansatzsp_t, typename testsp_t, typename value_t >
void
TInvarBasisQuadBEMBF< ansatzsp_t, testsp_t, value_t >::eval  ( const std::vector< idx_t > &  row_ind,
                                                               const std::vector< idx_t > &  col_ind,
                                                               BLAS::Matrix< value_t > &     values ) const
{
    eval( row_ind, col_ind, values, this->_ada_error );
}

template < typename ansatzsp_t, typename testsp_t, typename value_t >
void
TInvarBasisQuadBEMBF< ansatzsp_t, testsp_t, value_t >::eval  ( const std::vector< idx_t > &  row_ind,
                                                               const std::vector< idx_t > &  col_ind,
                                                               BLAS::Matrix< value_t > &     values,
                                                               const double                  eps ) const
{
    using  real_t = real_type_t< value_t >;
    
    //
    // local types for storing triangle sets and for
    // mapping indices to position in \a values
    //
    
    using  tri_set_t = std::set< idx_t >;
    using  val_map_t = std::unordered_map< idx_t, idx_t >;

    //
    // store position of all indices in <rowind> and <colind>
    // for indirect access of <values>
    //
    
    const size_t  nrows = row_ind.size();
    const size_t  ncols = col_ind.size();
    val_map_t     rowmap, colmap;

    for ( size_t  i = 0; i < nrows; ++i )
        rowmap[ row_ind[i] ] = idx_t( i );

    for ( size_t  i = 0; i < ncols; ++i )
        colmap[ col_ind[i] ] = idx_t( i );
    
    //
    // collect all triangles in support of ansatz functions
    //

    const ansatzsp_t *  ansatz_sp   = this->ansatz_space();
    const TGrid *       ansatz_grid = ansatz_sp->grid();
    tri_set_t           ansatz_triangles;

    for ( size_t  i = 0; i < nrows; ++i )
    {
        for ( auto  tri_idx : ansatz_sp->support( row_ind[i] ) )
            ansatz_triangles.insert( tri_idx );
    }// for

    //
    // collect all triangles in support of test functions
    //

    const testsp_t *   test_sp   = this->test_space();
    const TGrid *      test_grid = test_sp->grid();
    tri_set_t          test_triangles;
    
    for ( size_t  i = 0; i < ncols; ++i )
    {
        for ( auto  tri_idx : test_sp->support( col_ind[i] ) )
            test_triangles.insert( tri_idx );
    }// for

    //
    // loop over each pair of triangles in ansatz and test space
    // and integrate kernel function
    //

    std::vector< value_t >  kernel_values;

    B::fill( value_t(0), values );

    for ( const auto  tri0idx : ansatz_triangles )
    {
        TGrid::triangle_t  tri0 = ansatz_grid->triangle( tri0idx );
        const real_t       J0   = real_t( ansatz_grid->tri_size( tri0idx ) );
        
        for ( const auto  tri1idx : test_triangles )
        {
            TGrid::triangle_t  tri1 = test_grid->triangle( tri1idx );
            const real_t       J1   = real_t( test_sp->grid()->tri_size( tri1idx ) );
        
            //
            // decide about relative position of triangles, choose quadrature points
            // and determine points of triangles s.t. order is as expected
            //

            uint  ncommon = this->reorder_common( tri0.vtx, tri1.vtx );
            uint  torder  = this->_quad_order;

            if ( ncommon == 0 )
                torder = this->adjust_order( tri0.vtx, tri1.vtx, torder );
            
            auto  rule = this->quad_rule( ncommon, torder );
    
            if ( rule == nullptr )
                HERROR( ERR_NULL, "(TInvarBasisQuadBEMBF) eval", "quadrature rule is nullptr" );

            auto  weights = rule->template weights< real_t >();
            
            //
            // compute kernel at quadrature points
            //

            const size_t  npts        = rule->npts;
            const size_t  padded_npts = CFG::Mach::simd_padded_size< value_t >( npts );

            if ( kernel_values.size() < padded_npts )
                kernel_values.resize( padded_npts );
                
            this->eval_kernel( tri0idx, tri1idx, tri0, tri1, rule, kernel_values, eps );
            
            //
            // compute final result for each basis function by combining
            // kernel values with precomputed basis function values
            //

            for ( auto  idx0 : ansatz_sp->triangle_indices( tri0idx ) )
            {
                // exclude indices not in <rowind>
                if ( rowmap.find( idx0 ) == rowmap.end() )
                    continue;

                const idx_t  pos0     = rowmap[ idx0 ];
                const auto   phi0_val = ansatz_val( idx0, tri0, ncommon, torder );
                
                for ( auto  idx1 : test_sp->triangle_indices( tri1idx ) )
                {
                    // exclude indices not in <colind>
                    if ( colmap.find( idx1 ) == colmap.end() )
                        continue;
                
                    const idx_t  pos1     = colmap[ idx1 ];
                    const auto   phi1_val = test_val( idx1, tri1, ncommon, torder );
                    value_t      value    = value_t(0);

                    for ( size_t  iv = 0; iv < npts; ++iv )
                        value += value_t( weights[iv] * (*phi0_val)[ iv ] * (*phi1_val)[ iv ] ) * kernel_values[ iv ];

                    values( pos0, pos1 ) += value * J0 * J1;
                }// for
            }// for
        }// for
    }// for
}

//
// specialization for constant basis functions: it is more efficient to
// handle this case per index instead of per triangle as in the general
// implementation above
// (only for real case; in complex case it is slower?)
//
template <>
void
TInvarBasisQuadBEMBF< TConstFnSpace< float >,
                      TConstFnSpace< float >,
                      float >::eval  ( const std::vector< idx_t > &  row_ind,
                                       const std::vector< idx_t > &  col_ind,
                                       BLAS::Matrix< float > &       values,
                                       const double                  eps ) const
{
    const size_t        nrows     = row_ind.size();
    const size_t        ncols     = col_ind.size();
    const ansatzsp_t *  ansatz_sp = this->ansatz_space();
    const testsp_t *    test_sp   = this->test_space();
        
    //
    // compute integral over support triangle pair for each (i,j)
    //
                
    if ( this->_quad_dist_adaptive && ( eps > 0 ) )
    {
        for ( size_t  col = 0; col < ncols; col++ )
        {
            const idx_t        j       = col_ind[ col ];
            const idx_t        tri1idx = *(test_sp->support(j).begin());
            const auto         J1      = value_t( test_sp->grid()->tri_size( tri1idx ) );
            TGrid::triangle_t  tri1    = test_sp->grid()->triangle( tri1idx );
        
            for ( size_t  row = 0; row < nrows; row++ )
            {
                const idx_t        i       = row_ind[ row ];
                const idx_t        tri0idx = *(ansatz_sp->support(i).begin());
                const auto         J0      = value_t( ansatz_sp->grid()->tri_size( tri0idx ) );
                TGrid::triangle_t  tri0    = ansatz_sp->grid()->triangle( tri0idx );

                //
                // increase quadrature order until error criterion is met
                //
            
                auto    value   = value_t(0);
                bool    init    = true;
                auto    ncommon = reorder_common( tri0.vtx, tri1.vtx );

                for ( uint  torder = 1; torder <= this->_quad_order; ++torder )
                {
                    auto  rule = quad_rule( ncommon, torder );

                    if ( rule == nullptr )
                        HERROR( ERR_NULL, "(TInvarBasisQuadBEMBF) eval", "quadrature rule is nullptr" );

                    const size_t  npts        = rule->npts;
                    auto          weights     = rule->weights< real_t >();
                    const size_t  padded_npts = CFG::Mach::simd_padded_size< value_t >( npts );
                    auto          kvalues     = std::vector< value_t >( padded_npts );

                    eval_kernel( tri0idx, tri1idx, tri0, tri1, rule, kvalues, eps );

                    auto  val = value_t( 0 );
                
                    for ( size_t  iv = 0; iv < npts; ++iv )
                        val += value_t( weights[iv] ) * kvalues[ iv ];
                
                    const auto  diff = std::abs( value - val );

                    value = val;
                
                    if ( init )
                        init = false;
                    else
                    {
                        if ( diff / value <= eps )
                            break;
                    }// else
                }// for
            
                values( idx_t( row ), idx_t( col ) ) = J0 * J1 * value;
            }// for
        }// for
    }// if
    else
    {
        std::vector< value_t >  kernel_values;
    
        for ( size_t  col = 0; col < ncols; col++ )
        {
            const idx_t        j       = col_ind[ col ];
            const idx_t        tri1idx = *(test_sp->support(j).begin());
            const auto         J1      = value_t( test_sp->grid()->tri_size( tri1idx ) );
            TGrid::triangle_t  tri1    = test_sp->grid()->triangle( tri1idx );
        
            for ( size_t  row = 0; row < nrows; row++ )
            {
                const idx_t        i       = row_ind[ row ];
                const idx_t        tri0idx = *(ansatz_sp->support(i).begin());
                const auto         J0      = value_t( ansatz_sp->grid()->tri_size( tri0idx ) );
                TGrid::triangle_t  tri0    = ansatz_sp->grid()->triangle( tri0idx );
        
                //
                // decide about relative position of triangles, choose quadrature points
                // and determine points of triangles s.t. order is as expected
                //

                uint  ncommon = reorder_common( tri0.vtx, tri1.vtx );
                uint  torder  = this->_quad_order;

                if ( this->_quad_dist_adaptive && ( ncommon == 0 ))
                    torder = adjust_order( tri0.vtx, tri1.vtx, torder );
                
                auto  rule = quad_rule( ncommon, torder );
                
                if ( rule == nullptr )
                    HERROR( ERR_NULL, "(TInvarBasisQuadBEMBF) eval", "quadrature rule is nullptr" );

                //
                // compute kernel at quadrature points
                //

                value_t       value       = value_t(0);
                const size_t  npts        = rule->npts;
                auto          weights     = rule->weights< real_t >();
                const size_t  padded_npts = CFG::Mach::simd_padded_size< value_t >( npts );

                if ( kernel_values.size() < padded_npts )
                    kernel_values.resize( padded_npts );
            
                eval_kernel( tri0idx, tri1idx, tri0, tri1, rule, kernel_values, eps );
            
                for ( size_t  iv = 0; iv < npts; ++iv )
                    value += value_t( weights[iv] ) * kernel_values[ iv ];

                values( idx_t( row ), idx_t( col ) ) = J0 * J1 * value;
            }// for
        }// for
    }// if
}

template <>
void
TInvarBasisQuadBEMBF< TConstFnSpace< float >,
                      TConstFnSpace< float >,
                      std::complex< float > >::eval  ( const std::vector< idx_t > &             row_ind,
                                                       const std::vector< idx_t > &             col_ind,
                                                       BLAS::Matrix< std::complex< float > > &  values,
                                                       const double                             eps ) const
{
    const size_t        nrows     = row_ind.size();
    const size_t        ncols     = col_ind.size();
    const ansatzsp_t *  ansatz_sp = this->ansatz_space();
    const testsp_t *    test_sp   = this->test_space();
        
    //
    // compute integral over support triangle pair for each (i,j)
    //
                
    if ( this->_quad_dist_adaptive && ( eps > 0 ) )
    {
        for ( size_t  col = 0; col < ncols; col++ )
        {
            const idx_t        j       = col_ind[ col ];
            const idx_t        tri1idx = *(test_sp->support(j).begin());
            const auto         J1      = value_t( test_sp->grid()->tri_size( tri1idx ) );
            TGrid::triangle_t  tri1    = test_sp->grid()->triangle( tri1idx );
        
            for ( size_t  row = 0; row < nrows; row++ )
            {
                const idx_t        i       = row_ind[ row ];
                const idx_t        tri0idx = *(ansatz_sp->support(i).begin());
                const auto         J0      = value_t( ansatz_sp->grid()->tri_size( tri0idx ) );
                TGrid::triangle_t  tri0    = ansatz_sp->grid()->triangle( tri0idx );

                //
                // increase quadrature order until error criterion is met
                //
            
                auto    value   = value_t(0);
                bool    init    = true;
                auto    ncommon = reorder_common( tri0.vtx, tri1.vtx );

                for ( uint  torder = 1; torder <= this->_quad_order; ++torder )
                {
                    auto  rule = quad_rule( ncommon, torder );

                    if ( rule == nullptr )
                        HERROR( ERR_NULL, "(TInvarBasisQuadBEMBF) eval", "quadrature rule is nullptr" );

                    const size_t  npts        = rule->npts;
                    auto          weights     = rule->weights< real_t >();
                    const size_t  padded_npts = CFG::Mach::simd_padded_size< value_t >( npts );
                    auto          kvalues     = std::vector< value_t >( padded_npts );

                    eval_kernel( tri0idx, tri1idx, tri0, tri1, rule, kvalues, eps );

                    auto  val = value_t( 0 );
                
                    for ( size_t  iv = 0; iv < npts; ++iv )
                        val += value_t( weights[iv] ) * kvalues[ iv ];
                
                    const auto  diff = std::abs( value - val );

                    value = val;
                
                    if ( init )
                        init = false;
                    else
                    {
                        if ( diff / std::abs( value ) <= eps )
                            break;
                    }// else
                }// for
            
                values( idx_t( row ), idx_t( col ) ) = J0 * J1 * value;
            }// for
        }// for
    }// if
    else
    {
        std::vector< value_t >  kernel_values;
    
        for ( size_t  col = 0; col < ncols; col++ )
        {
            const idx_t        j       = col_ind[ col ];
            const idx_t        tri1idx = *(test_sp->support(j).begin());
            const auto         J1      = value_t( test_sp->grid()->tri_size( tri1idx ) );
            TGrid::triangle_t  tri1    = test_sp->grid()->triangle( tri1idx );
        
            for ( size_t  row = 0; row < nrows; row++ )
            {
                const idx_t        i       = row_ind[ row ];
                const idx_t        tri0idx = *(ansatz_sp->support(i).begin());
                const auto         J0      = value_t( ansatz_sp->grid()->tri_size( tri0idx ) );
                TGrid::triangle_t  tri0    = ansatz_sp->grid()->triangle( tri0idx );
        
                //
                // decide about relative position of triangles, choose quadrature points
                // and determine points of triangles s.t. order is as expected
                //

                uint  ncommon = reorder_common( tri0.vtx, tri1.vtx );
                uint  torder  = this->_quad_order;

                if ( this->_quad_dist_adaptive && ( ncommon == 0 ))
                    torder = adjust_order( tri0.vtx, tri1.vtx, torder );
                
                auto  rule = quad_rule( ncommon, torder );
                
                if ( rule == nullptr )
                    HERROR( ERR_NULL, "(TInvarBasisQuadBEMBF) eval", "quadrature rule is nullptr" );

                //
                // compute kernel at quadrature points
                //

                value_t       value       = value_t(0);
                const size_t  npts        = rule->npts;
                auto          weights     = rule->weights< real_t >();
                const size_t  padded_npts = CFG::Mach::simd_padded_size< value_t >( npts );

                if ( kernel_values.size() < padded_npts )
                    kernel_values.resize( padded_npts );
            
                eval_kernel( tri0idx, tri1idx, tri0, tri1, rule, kernel_values, eps );
            
                for ( size_t  iv = 0; iv < npts; ++iv )
                    value += value_t( weights[iv] ) * kernel_values[ iv ];

                values( idx_t( row ), idx_t( col ) ) = J0 * J1 * value;
            }// for
        }// for
    }// if
}

template <>
void
TInvarBasisQuadBEMBF< TConstFnSpace< double >,
                      TConstFnSpace< double >,
                      double >::eval  ( const std::vector< idx_t > &  row_ind,
                                        const std::vector< idx_t > &  col_ind,
                                        BLAS::Matrix< double > &      values,
                                        const double                  eps ) const
{
    const size_t        nrows     = row_ind.size();
    const size_t        ncols     = col_ind.size();
    const ansatzsp_t *  ansatz_sp = this->ansatz_space();
    const testsp_t *    test_sp   = this->test_space();
        
    //
    // compute integral over support triangle pair for each (i,j)
    //
                
    if ( this->_quad_dist_adaptive && ( eps > 0 ) )
    {
        for ( size_t  col = 0; col < ncols; col++ )
        {
            const idx_t        j       = col_ind[ col ];
            const idx_t        tri1idx = *(test_sp->support(j).begin());
            const auto         J1      = value_t( test_sp->grid()->tri_size( tri1idx ) );
            TGrid::triangle_t  tri1    = test_sp->grid()->triangle( tri1idx );
        
            for ( size_t  row = 0; row < nrows; row++ )
            {
                const idx_t        i       = row_ind[ row ];
                const idx_t        tri0idx = *(ansatz_sp->support(i).begin());
                const auto         J0      = value_t( ansatz_sp->grid()->tri_size( tri0idx ) );
                TGrid::triangle_t  tri0    = ansatz_sp->grid()->triangle( tri0idx );

                //
                // increase quadrature order until error criterion is met
                //
            
                auto    value   = value_t(0);
                bool    init    = true;
                auto    ncommon = reorder_common( tri0.vtx, tri1.vtx );

                for ( uint  torder = 1; torder <= this->_quad_order; ++torder )
                {
                    auto  rule = quad_rule( ncommon, torder );

                    if ( rule == nullptr )
                        HERROR( ERR_NULL, "(TInvarBasisQuadBEMBF) eval", "quadrature rule is nullptr" );

                    const size_t  npts        = rule->npts;
                    auto          weights     = rule->weights< real_t >();
                    const size_t  padded_npts = CFG::Mach::simd_padded_size< value_t >( npts );
                    auto          kvalues     = std::vector< value_t >( padded_npts );

                    eval_kernel( tri0idx, tri1idx, tri0, tri1, rule, kvalues, eps );

                    auto  val = value_t( 0 );
                
                    for ( size_t  iv = 0; iv < npts; ++iv )
                        val += value_t( weights[iv] ) * kvalues[ iv ];
                
                    const auto  diff = std::abs( value - val );

                    value = val;
                
                    if ( init )
                        init = false;
                    else
                    {
                        if ( diff / value <= eps )
                            break;
                    }// else
                }// for
            
                values( idx_t( row ), idx_t( col ) ) = J0 * J1 * value;
            }// for
        }// for
    }// if
    else
    {
        std::vector< value_t >  kernel_values;
    
        for ( size_t  col = 0; col < ncols; col++ )
        {
            const idx_t        j       = col_ind[ col ];
            const idx_t        tri1idx = *(test_sp->support(j).begin());
            const auto         J1      = value_t( test_sp->grid()->tri_size( tri1idx ) );
            TGrid::triangle_t  tri1    = test_sp->grid()->triangle( tri1idx );
        
            for ( size_t  row = 0; row < nrows; row++ )
            {
                const idx_t        i       = row_ind[ row ];
                const idx_t        tri0idx = *(ansatz_sp->support(i).begin());
                const auto         J0      = value_t( ansatz_sp->grid()->tri_size( tri0idx ) );
                TGrid::triangle_t  tri0    = ansatz_sp->grid()->triangle( tri0idx );

                //
                // decide about relative position of triangles, choose quadrature points
                // and determine points of triangles s.t. order is as expected
                //

                uint  ncommon = reorder_common( tri0.vtx, tri1.vtx );
                uint  torder  = this->_quad_order;

                if ( this->_quad_dist_adaptive && ( ncommon == 0 ))
                    torder = adjust_order( tri0.vtx, tri1.vtx, torder );
                
                auto  rule = quad_rule( ncommon, torder );
                
                if ( rule == nullptr )
                    HERROR( ERR_NULL, "(TInvarBasisQuadBEMBF) eval", "quadrature rule is nullptr" );

                //
                // compute kernel at quadrature points
                //

                value_t       value       = value_t(0);
                const size_t  npts        = rule->npts;
                auto          weights     = rule->weights< real_t >();
                const size_t  padded_npts = CFG::Mach::simd_padded_size< value_t >( npts );

                if ( kernel_values.size() < padded_npts )
                    kernel_values.resize( padded_npts );
            
                eval_kernel( tri0idx, tri1idx, tri0, tri1, rule, kernel_values, eps );
            
                for ( size_t  iv = 0; iv < npts; ++iv )
                    value += value_t( weights[iv] ) * kernel_values[ iv ];

                values( idx_t( row ), idx_t( col ) ) = J0 * J1 * value;
            }// for
        }// for
    }// else
}

template <>
void
TInvarBasisQuadBEMBF< TConstFnSpace< double >,
                      TConstFnSpace< double >,
                      std::complex< double > >::eval  ( const std::vector< idx_t > &              row_ind,
                                                        const std::vector< idx_t > &              col_ind,
                                                        BLAS::Matrix< std::complex< double > > &  values,
                                                        const double                              eps ) const
{
    const size_t        nrows     = row_ind.size();
    const size_t        ncols     = col_ind.size();
    const ansatzsp_t *  ansatz_sp = this->ansatz_space();
    const testsp_t *    test_sp   = this->test_space();
        
    //
    // compute integral over support triangle pair for each (i,j)
    //
                
    if ( this->_quad_dist_adaptive && ( eps > 0 ) )
    {
        for ( size_t  col = 0; col < ncols; col++ )
        {
            const idx_t        j       = col_ind[ col ];
            const idx_t        tri1idx = *(test_sp->support(j).begin());
            const auto         J1      = value_t( test_sp->grid()->tri_size( tri1idx ) );
            TGrid::triangle_t  tri1    = test_sp->grid()->triangle( tri1idx );
        
            for ( size_t  row = 0; row < nrows; row++ )
            {
                const idx_t        i       = row_ind[ row ];
                const idx_t        tri0idx = *(ansatz_sp->support(i).begin());
                const auto         J0      = value_t( ansatz_sp->grid()->tri_size( tri0idx ) );
                TGrid::triangle_t  tri0    = ansatz_sp->grid()->triangle( tri0idx );

                //
                // increase quadrature order until error criterion is met
                //
            
                auto    value   = value_t(0);
                bool    init    = true;
                auto    ncommon = reorder_common( tri0.vtx, tri1.vtx );

                for ( uint  torder = 1; torder <= this->_quad_order; ++torder )
                {
                    auto  rule = quad_rule( ncommon, torder );

                    if ( rule == nullptr )
                        HERROR( ERR_NULL, "(TInvarBasisQuadBEMBF) eval", "quadrature rule is nullptr" );

                    const size_t  npts        = rule->npts;
                    auto          weights     = rule->weights< real_t >();
                    const size_t  padded_npts = CFG::Mach::simd_padded_size< value_t >( npts );
                    auto          kvalues     = std::vector< value_t >( padded_npts );

                    eval_kernel( tri0idx, tri1idx, tri0, tri1, rule, kvalues, eps );

                    auto  val = value_t( 0 );
                
                    for ( size_t  iv = 0; iv < npts; ++iv )
                        val += value_t( weights[iv] ) * kvalues[ iv ];
                
                    const auto  diff = std::abs( value - val );

                    value = val;
                
                    if ( init )
                        init = false;
                    else
                    {
                        if ( diff / std::abs( value ) <= eps )
                            break;
                    }// else
                }// for
            
                values( idx_t( row ), idx_t( col ) ) = J0 * J1 * value;
            }// for
        }// for
    }// if
    else
    {
        std::vector< value_t >  kernel_values;
    
        for ( size_t  col = 0; col < ncols; col++ )
        {
            const idx_t        j       = col_ind[ col ];
            const idx_t        tri1idx = *(test_sp->support(j).begin());
            const auto         J1      = value_t( test_sp->grid()->tri_size( tri1idx ) );
            TGrid::triangle_t  tri1    = test_sp->grid()->triangle( tri1idx );
        
            for ( size_t  row = 0; row < nrows; row++ )
            {
                const idx_t        i       = row_ind[ row ];
                const idx_t        tri0idx = *(ansatz_sp->support(i).begin());
                const auto         J0      = value_t( ansatz_sp->grid()->tri_size( tri0idx ) );
                TGrid::triangle_t  tri0    = ansatz_sp->grid()->triangle( tri0idx );

                //
                // decide about relative position of triangles, choose quadrature points
                // and determine points of triangles s.t. order is as expected
                //

                uint  ncommon = reorder_common( tri0.vtx, tri1.vtx );
                uint  torder  = this->_quad_order;

                if ( this->_quad_dist_adaptive && ( ncommon == 0 ))
                    torder = adjust_order( tri0.vtx, tri1.vtx, torder );
                
                auto  rule = quad_rule( ncommon, torder );
                
                if ( rule == nullptr )
                    HERROR( ERR_NULL, "(TInvarBasisQuadBEMBF) eval", "quadrature rule is nullptr" );

                //
                // compute kernel at quadrature points
                //

                value_t       value       = value_t(0);
                const size_t  npts        = rule->npts;
                auto          weights     = rule->weights< real_t >();
                const size_t  padded_npts = CFG::Mach::simd_padded_size< value_t >( npts );

                if ( kernel_values.size() < padded_npts )
                    kernel_values.resize( padded_npts );
            
                eval_kernel( tri0idx, tri1idx, tri0, tri1, rule, kernel_values, eps );
            
                for ( size_t  iv = 0; iv < npts; ++iv )
                    value += value_t( weights[iv] ) * kernel_values[ iv ];

                values( idx_t( row ), idx_t( col ) ) = J0 * J1 * value;
            }// for
        }// for
    }// else
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

#define INST_ALL( type1, type2 )                                        \
    template class TInvarBasisQuadBEMBF< TConstFnSpace< type1 >,  TConstFnSpace< type1 >,  type2 >; \
    template class TInvarBasisQuadBEMBF< TConstFnSpace< type1 >,  TLinearFnSpace< type1 >, type2 >; \
    template class TInvarBasisQuadBEMBF< TLinearFnSpace< type1 >, TConstFnSpace< type1 >,  type2 >; \
    template class TInvarBasisQuadBEMBF< TLinearFnSpace< type1 >, TLinearFnSpace< type1 >, type2 >;

INST_ALL( float,  float )
INST_ALL( double, double )
INST_ALL( float,  std::complex< float > )
INST_ALL( double, std::complex< double > )

}// namespace Hpro
