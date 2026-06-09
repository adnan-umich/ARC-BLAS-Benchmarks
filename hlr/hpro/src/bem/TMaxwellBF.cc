//
// Project     : HLIBpro
// File        : TMaxwellBF.cc
// Description : bilinear forms for Maxwell operator
// Author      : Ronald Kriemann, Jonas Ballani
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include <set>

#include "unordered_map.hh"

#include "hpro/base/packed.hh"

#include "hpro/bem/TGaussQuad.hh"
#include "hpro/bem/TMaxwellBF.hh"
#include "hpro/bem/TConstEdgeFnSpace.hh"

namespace Hpro
{

namespace B = BLAS;

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// external function definitions
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//
// Helmholtz SLP kernels to be used in TMaxwellEFIEBF
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  value_t >
void
helmholtz_slp_flt ( const TGrid::triangle_t &    tri0,
                    const TGrid::triangle_t &    tri1,
                    const tripair_quad_rule_t *  rule,
                    const value_t                ikappa,
                    const T_ansatzsp *           ansatz_sp,
                    const T_testsp *             test_sp,
                    std::vector< value_t > &     values );

#define  HELMHOLTZ_SLP( suffix )                                        \
    template < typename  value_t,                                       \
               typename  T_ansatzsp,                                    \
               typename  T_testsp,                                      \
               typename  T_packed >                                     \
    void                                                                \
    helmholtz_slp_##suffix ( const TGrid::triangle_t &    tri0,         \
                             const TGrid::triangle_t &    tri1,         \
                             const tripair_quad_rule_t *  rule,         \
                             const value_t                ikappa,       \
                             const T_ansatzsp *           ansatz_sp,    \
                             const T_testsp *             test_sp,      \
                             std::vector< value_t > &     values )

HELMHOLTZ_SLP( simd );
HELMHOLTZ_SLP( re_simd );
HELMHOLTZ_SLP( im_simd );

//
// Helmholtz DLP kernels to be used in TMaxwellMFIEBF
//
#define  HELMHOLTZ_DLP_SIMD( suffix )                                   \
    template < typename  T_ansatzsp,                                    \
               typename  T_testsp,                                      \
               typename  T_packed >                                     \
    void                                                                \
    helmholtz_dlp_wo_normal_##suffix ( const TGrid::triangle_t &                tri0, \
                                       const TGrid::triangle_t &                tri1, \
                                       const tripair_quad_rule_t *              rule, \
                                       const std::complex< double >             ikappa, \
                                       const T_ansatzsp *                       ansatz_sp, \
                                       const T_testsp *                         test_sp, \
                                       std::vector< std::complex< double > > &  values )

HELMHOLTZ_DLP_SIMD( simd );
HELMHOLTZ_DLP_SIMD( re_simd );
HELMHOLTZ_DLP_SIMD( im_simd );

namespace
{

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// local constants
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

const double ONE_OVER_4PI = double(1) / (double(4) * Math::pi< double >());

}// namespace anonymous

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// TMaxwellEFIEBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template < typename  ansatzsp_t,
           typename  testsp_t,
           typename  value_t >
auto
efie_kernel_fn ( const value_t  ikappa )
{
    using  real_t = real_type_t< value_t >;

    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            using  packed_t = packed< real_t, ISA_AVX512F >;
            
            HINFO( "(TMaxwellEFIEBF) using AVX512F kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic )
        {
            using  packed_t = packed< real_t, ISA_MIC >;
            
            HINFO( "(TMaxwellEFIEBF) using MIC kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            using  packed_t = packed< real_t, ISA_AVX2 >;
            
            HINFO( "(TMaxwellEFIEBF) using AVX2 kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            using  packed_t = packed< real_t, ISA_AVX >;
            
            HINFO( "(TMaxwellEFIEBF) using AVX kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            using  packed_t = packed< real_t, ISA_SSE3 >;
            
            HINFO( "(TMaxwellEFIEBF) using SSE3 kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            using  packed_t = packed< real_t, ISA_VSX >;
            
            HINFO( "(TMaxwellEFIEBF) using VSX kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            using  packed_t = packed< real_t, ISA_NEON >;
            
            HINFO( "(TMaxwellEFIEBF) using NEON kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else
        {
            HINFO( "(TMaxwellEFIEBF) using standard kernel" );
            
            return helmholtz_slp_flt< ansatzsp_t, testsp_t, value_t >;
        }// else
    }// if
    else
    {
        HINFO( "(TMaxwellEFIEBF) using standard kernel" );
            
        return helmholtz_slp_flt< ansatzsp_t, testsp_t, value_t >;
    }// else
}

//
// ctor
// 
template < typename  T_ansatzsp,
           typename  T_testsp >
TMaxwellEFIEBF< T_ansatzsp, T_testsp >::TMaxwellEFIEBF ( const value_t       kappa,
                                                         const real_t        eta,
                                                         const ansatzsp_t *  aansatzsp,
                                                         const testsp_t *    atestsp,
                                                         const uint          quad_order )
        : TQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, quad_order )
        , _ikappa( value_t( 0, 1 ) * kappa )
        , _eta( eta )
        , _ikappa_times_eta( _ikappa * _eta )
        , _eta_over_ikappa( _eta / _ikappa )
        , _kernel_fn( efie_kernel_fn< T_ansatzsp, T_testsp, value_t >( _ikappa ) )
{}

template < typename  T_ansatzsp,
           typename  T_testsp >
TMaxwellEFIEBF< T_ansatzsp, T_testsp >::TMaxwellEFIEBF ( const value_t       kappa,
                                                         const real_t        eta,
                                                         const ansatzsp_t *  aansatzsp,
                                                         const testsp_t *    atestsp,
                                                         const real_t        quad_error )
        : TQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, 10, true, quad_error )
        , _ikappa( value_t( 0, 1 ) * kappa )
        , _eta( eta )
        , _ikappa_times_eta( _ikappa * _eta )
        , _eta_over_ikappa( _eta / _ikappa )
        , _kernel_fn( efie_kernel_fn< T_ansatzsp, T_testsp, value_t >( _ikappa ) )
{}

template < typename  T_ansatzsp,
           typename  T_testsp >
matform_t
TMaxwellEFIEBF< T_ansatzsp, T_testsp >::format () const
{
    if ( ( reinterpret_cast< const void * >( this->ansatz_space() ) == reinterpret_cast< const void * >( this->test_space() ) ) ||
         (( this->ansatz_space()->type() == this->test_space()->type() ) &&
          ( this->ansatz_space()->grid() == this->test_space()->grid() )) )
        return symmetric;
    else
        return unsymmetric;
}

//
// evaluate subblock defined by \a row_ind × \a col_ind; the indices
// in \a row_ind and \a col_ind can be arbitrary, e.g. must not be
// contiguous
//
template < typename  T_ansatzsp,
           typename  T_testsp >
void
TMaxwellEFIEBF< T_ansatzsp, T_testsp >::eval  ( const std::vector< idx_t > &  row_ind,
                                                const std::vector< idx_t > &  col_ind,
                                                BLAS::Matrix< value_t > &     values ) const
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
    // for indirect acces of <values>
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
        for ( auto  triangle : ansatz_sp->support( row_ind[i] ) )
            ansatz_triangles.insert( triangle );
    }// for

    //
    // collect all triangles in support of test functions
    //

    const testsp_t *   test_sp   = this->test_space();
    const TGrid *      test_grid = test_sp->grid();
    tri_set_t          test_triangles;

    for ( size_t  i = 0; i < ncols; ++i )
    {
        for ( auto  triangle : test_sp->support( col_ind[i] ) )
            test_triangles.insert( triangle );
    }// for

    //
    // loop over each pair of triangles in ansatz and test space
    // and integrate kernel function
    //

    std::vector< value_t >  kernel_values;
    std::vector< T3Point >  phi0_values;
    std::vector< T3Point >  phi1_values;

    B::fill( value_t(0), values );

    for ( auto  tri0idx : ansatz_triangles )
    {
        TGrid::triangle_t  tri0    = ansatz_grid->triangle( tri0idx );
        const real_t       J0      = real_t( ansatz_grid->tri_size( tri0idx ) );

        for ( auto  tri1idx : test_triangles )
        {
            TGrid::triangle_t  tri1    = test_grid->triangle( tri1idx );
            const real_t       J1      = real_t( test_sp->grid()->tri_size( tri1idx ) );

            //
            // decide about relative position of triangles, choose quadrature points
            // and determine points of triangles s.t. order is as expected
            //

            uint  ncommon = this->reorder_common( tri0.vtx, tri1.vtx );
            uint  torder  = this->_quad_order;

            if ( ncommon == 0 )
                torder = this->adjust_order( tri0.vtx, tri1.vtx, torder );

            auto  rule    = this->quad_rule( ncommon, torder );
            auto  weights = rule -> template weights< real_t >();

            if ( rule == nullptr )
                HERROR( ERR_NULL, "(TMaxwellEFIEBF) eval", "quadrature rule is nullptr" );

            //
            // compute kernel at quadrature points
            //

            const size_t  npts        = rule->npts;
            const size_t  padded_npts = CFG::Mach::simd_padded_size< value_t >( npts );

            if ( kernel_values.size() < padded_npts )
            {
                kernel_values.resize( padded_npts );
                phi0_values.resize( padded_npts );
                phi1_values.resize( padded_npts );
            }// if

            eval_kernel( tri0idx, tri1idx, tri0, tri1, rule, kernel_values );

            //
            // compute final result for each basis function by combining
            // kernel values with precomputed basis function values
            //

            for ( auto  idx0 : ansatz_sp->triangle_indices( tri0idx ) )
            {
                // exclude indices not in <rowind>
                if ( rowmap.find( idx0 ) == rowmap.end() )
                    continue;

                const idx_t   pos0         = rowmap[ idx0 ];
                const real_t  scale_phi0   = ansatz_sp->scaling_factor_basis( idx0, tri0idx );
                const bool    phi0_in_supp = ansatz_sp->eval_basis( idx0, tri0, rule->pts1, phi0_values );

                for ( auto  idx1 : test_sp->triangle_indices( tri1idx ) )
                {
                    // exclude indices not in <colind>
                    if ( colmap.find( idx1 ) == colmap.end() )
                        continue;

                    const idx_t   pos1         = colmap[ idx1 ];
                    const real_t  scale_phi1   = test_sp->scaling_factor_basis( idx1, tri1idx );
                    const bool    phi1_in_supp = test_sp->eval_basis( idx1, tri1, rule->pts2, phi1_values );
                    
                    value_t       value        = value_t(0);
                    
                    if ( phi0_in_supp && phi1_in_supp )
                    {
                        const real_t   scale_phi = scale_phi0 * scale_phi1;
                        const real_t   div_phi   = 4 * scale_phi;
                        const value_t  eik_div   = _eta_over_ikappa * div_phi;
                        
                        for ( size_t  iv = 0; iv < npts; ++iv )
                        {
                            const real_t  dot_phi = scale_phi * real_t( dot( phi0_values[iv], phi1_values[iv] ) );
                            
                            value += real_t( weights[iv] ) * kernel_values[ iv ] * ( -_ikappa_times_eta * dot_phi - eik_div );
                        }// for
                    }// if

                    values( pos0, pos1 ) += value * J0 * J1;
                }// for
            }// for
        }// for
    }// for
}

//
// eval kernel function at quadrature points
//
template < typename  T_ansatzsp,
           typename  T_testsp >
void
TMaxwellEFIEBF< T_ansatzsp, T_testsp >::eval_kernel ( const                        idx_t,
                                                      const                        idx_t,
                                                      const TGrid::triangle_t &    tri0,
                                                      const TGrid::triangle_t &    tri1,
                                                      const tripair_quad_rule_t *  rule,
                                                      std::vector< value_t > &     values,
                                                      const real_t                 eps ) const
{
    _kernel_fn( tri0, tri1, rule, _ikappa, this->ansatz_space(), this->test_space(), values );
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// TMaxwellEFIEMassBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//
// ctor
//
template < typename T_ansatzsp, typename T_testsp >
TMaxwellEFIEMassBF< T_ansatzsp, T_testsp >::TMaxwellEFIEMassBF ( const ansatzsp_t *  aansatzsp,
                                                                 const testsp_t *    atestsp,
                                                                 const uint          aorder )
        : TBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp )
{
    //
    // build quad-points and weights
    //

    TTriGaussQuad  gaussquad;

    gaussquad.build( std::max( 1u, aorder ), _quad_pts, _quad_wghts );
}
    
template < typename  T_ansatzsp,
           typename  T_testsp >
matform_t
TMaxwellEFIEMassBF< T_ansatzsp, T_testsp >::format () const
{
    return unsymmetric;
}

//
// evaluate bilinearform
//
template < typename T_ansatzsp, typename T_testsp >
void
TMaxwellEFIEMassBF< T_ansatzsp, T_testsp >::eval ( const std::vector< idx_t > &  row_ind,
                                                   const std::vector< idx_t > &  col_ind,
                                                   BLAS::Matrix< value_t > &     values ) const
{
    const size_t  nrows = row_ind.size();
    const size_t  ncols = col_ind.size();

    for ( size_t  col = 0; col < ncols; col++ )
    {
        const idx_t  j = col_ind[col];
        
        for ( size_t  row = 0; row < nrows; row++ )
        {
            const idx_t  i = row_ind[row];
        
            //
            // integrate over test space
            //

            auto                    value     = value_t(0);
            const auto              ansatz_sp = this->ansatz_space();
            const auto              test_sp   = this->test_space();
            auto                    support_i = ansatz_sp->support(i);
            auto                    support_j = test_sp->support(j);
            std::vector< T3Point >  phi0_values;
            std::vector< T3Point >  phi1_values;

            if ( ansatz_sp->grid() != test_sp->grid() )
                HERROR( ERR_NOT_IMPL, "(TMaxwellEFIEMassBF) eval", "unequal grids not supported" );

            for ( auto  tri1 : support_j )
            {
                const auto  J1 = value_t( test_sp->grid()->tri_size( tri1 ) );

                //
                // test if support of ansatz function is disjoint
                // and immediately go one if so
                //

                bool  intersects = false;

                for ( auto  tri0 : support_i )
                {
                    if ( tri0 == tri1 ) // only valid for equal grids!
                    {
                        intersects = true;
                        break;
                    }// if
                }// for

                if ( ! intersects )
                    continue;

                //
                // evaluate basis functions at quadrature points
                //

                TGrid::triangle_t  t1;
                const size_t       npts   = _quad_pts.size();

                if ( phi0_values.size() < npts )
                    phi0_values.resize( npts );
                if ( phi1_values.size() < npts )
                    phi1_values.resize( npts );

                t1 = test_sp->grid()->triangle( tri1 );

                const bool  phi0_in_supp = ansatz_sp->eval_basis( i, t1, _quad_pts, phi0_values );
                const bool  phi1_in_supp = test_sp->eval_basis(   j, t1, _quad_pts, phi1_values );
                auto        tvalue       = value_t(0);

                if ( phi0_in_supp && phi1_in_supp )
                {
                    const value_t  scale_phi0 = ansatz_sp->scaling_factor_basis( i, tri1 );
                    const value_t  scale_phi1 = test_sp->scaling_factor_basis(   j, tri1 );
                    const value_t  scale_phi  = scale_phi0 * scale_phi1;

                    for ( size_t  k = 0; k < npts; ++k )
                    {
                        const value_t  dot_phi = scale_phi * value_t( dot( phi0_values[k], phi1_values[k] ) );
                        
                        tvalue += value_t( _quad_wghts[k] ) * dot_phi;
                    }// for
                }// if

                value += J1 * tvalue;
            }// for

            values( idx_t( row ), idx_t( col ) ) = value;
        }// for
    }// for
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// TMaxwellMFIEBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//
// evaluation of Helmholtz DLP kernel without normal directory
//
template < typename  T_ansatzsp,
           typename  T_testsp >
void
helmholtz_dlp_wo_normal_flt ( const TGrid::triangle_t &                tri0,
                              const TGrid::triangle_t &                tri1,
                              const tripair_quad_rule_t *              rule,
                              const std::complex< double >             ikappa,
                              const T_ansatzsp *                       ansatz_sp,
                              const T_testsp *                         test_sp,
                              std::vector< std::complex< double > > &  values )
{
    using real_t = double;
    
    const size_t   n_pts = rule->npts;
    const auto     x1    = rule -> template x1< real_t >();
    const auto     y1    = rule -> template y1< real_t >();
    const auto     x2    = rule -> template x2< real_t >();
    const auto     y2    = rule -> template y2< real_t >();

    const T3Point  v00( ansatz_sp->grid()->vertex( tri0.vtx[0] ) );
    const T3Point  v01( ansatz_sp->grid()->vertex( tri0.vtx[1] ) );
    const T3Point  v02( ansatz_sp->grid()->vertex( tri0.vtx[2] ) );
    const T3Point  v10( test_sp->grid()->vertex( tri1.vtx[0] ) );
    const T3Point  v11( test_sp->grid()->vertex( tri1.vtx[1] ) );
    const T3Point  v12( test_sp->grid()->vertex( tri1.vtx[2] ) );

    #pragma omp simd
    for ( size_t  i = 0; i < n_pts; i++ )
    {
        const double  a1 = x1[i];
        const double  a2 = y1[i];
        const double  a0 = 1.0 - a1 - a2;
        
        const double  b1 = x2[i];
        const double  b2 = y2[i];
        const double  b0 = 1.0 - b1 - b2;

        // u = y - x
        const T3Point  x      = a0 * v00 + a1 * v01 + a2 * v02;
        const T3Point  y      = b0 * v10 + b1 * v11 + b2 * v12;
        const T3Point  ymx    = y - x;
        const auto     sqdist = real_t( dot( ymx, ymx ) );       // |y-x|²

        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n,y-x>
        // ────────────────────────────────────
        //                |x-y|³
        //

        const auto  dist   = Math::sqrt( sqdist );     // |y-x|
        const auto  cudist = sqdist * dist;            // |y-x|³
        const auto  ikymx  = ikappa * dist;            // i · κ · |x-y|
        
        values[ i ] = ONE_OVER_4PI * Math::exp( ikymx ) * ( ikymx - real_t(1) ) / cudist;
    }// for
}

template < typename  ansatzsp_t,
           typename  testsp_t,
           typename  value_t >
auto
mfie_kernel_fn ( const value_t  ikappa )
{
    using  real_t = real_type_t< value_t >;

    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            using  packed_t = packed< real_t, ISA_AVX512F >;

            HINFO( "(TMaxwellMFIEBF) using AVX512F kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_re_simd< ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_im_simd< ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_wo_normal_simd< ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic )
        {
            using  packed_t = packed< real_t, ISA_MIC >;

            HINFO( "(TMaxwellMFIEBF) using MIC kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_re_simd< ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_im_simd< ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_wo_normal_simd< ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            using  packed_t = packed< real_t, ISA_AVX2 >;

            HINFO( "(TMaxwellMFIEBF) using AVX2 kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_re_simd< ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_im_simd< ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_wo_normal_simd< ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            using  packed_t = packed< real_t, ISA_AVX >;

            HINFO( "(TMaxwellMFIEBF) using AVX kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_re_simd< ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_im_simd< ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_wo_normal_simd< ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            using  packed_t = packed< real_t, ISA_SSE3 >;

            HINFO( "(TMaxwellMFIEBF) using SSE3 kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_re_simd< ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_im_simd< ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_wo_normal_simd< ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            using  packed_t = packed< real_t, ISA_VSX >;

            HINFO( "(TMaxwellMFIEBF) using VSX kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_re_simd< ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_im_simd< ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_wo_normal_simd< ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            using  packed_t = packed< real_t, ISA_NEON >;

            HINFO( "(TMaxwellMFIEBF) using NEON kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_re_simd< ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_wo_normal_im_simd< ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_wo_normal_simd< ansatzsp_t, testsp_t, packed_t >;
        }// if
        else
        {
            HINFO( "(TMaxwellMFIEBF) using standard kernel" );
            
            return helmholtz_dlp_wo_normal_flt< ansatzsp_t, testsp_t >;
        }// else
    }// if
    else
    {
        HINFO( "(TMaxwellMFIEBF) using standard kernel" );
            
        return helmholtz_dlp_wo_normal_flt< ansatzsp_t, testsp_t >;
    }// else
}

//
// ctor
//
template < typename  T_ansatzsp,
           typename  T_testsp >
TMaxwellMFIEBF< T_ansatzsp, T_testsp >::TMaxwellMFIEBF ( const value_t       kappa,
                                                         const ansatzsp_t *  aansatzsp,
                                                         const testsp_t *    atestsp,
                                                         const uint          quad_order )
        : TQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, quad_order )
        , _ikappa( value_t( 0, 1 ) * kappa )
        , _kernel_fn( mfie_kernel_fn< T_ansatzsp, T_testsp, value_t >( _ikappa ) )
{}

template < typename  T_ansatzsp,
           typename  T_testsp >
TMaxwellMFIEBF< T_ansatzsp, T_testsp >::TMaxwellMFIEBF ( const value_t       kappa,
                                                         const ansatzsp_t *  aansatzsp,
                                                         const testsp_t *    atestsp,
                                                         const real_t        quad_error )
        : TQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, 10, true, quad_error )
        , _ikappa( value_t( 0, 1 ) * kappa )
        , _kernel_fn( mfie_kernel_fn< T_ansatzsp, T_testsp, value_t >( _ikappa ) )
{}

template < typename  T_ansatzsp,
           typename  T_testsp >
matform_t
TMaxwellMFIEBF< T_ansatzsp, T_testsp >::format () const
{
    return unsymmetric;
}

//
// evaluate subblock defined by \a row_ind × \a col_ind; the indices
// in \a row_ind and \a col_ind can be arbitrary, e.g. must not be
// contiguous
//
template < typename  T_ansatzsp,
           typename  T_testsp >
void
TMaxwellMFIEBF< T_ansatzsp, T_testsp >::eval  ( const std::vector< idx_t > &  row_ind,
                                                const std::vector< idx_t > &  col_ind,
                                                BLAS::Matrix< value_t > &     values ) const
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
    // for indirect acces of <values>
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
        for ( auto  triangle : ansatz_sp->support( row_ind[i] ) )
            ansatz_triangles.insert( triangle );
    }// for

    //
    // collect all triangles in support of test functions
    //

    const testsp_t *   test_sp   = this->test_space();
    const TGrid *      test_grid = test_sp->grid();
    tri_set_t          test_triangles;

    for ( size_t  i = 0; i < ncols; ++i )
    {
        for ( auto  triangle : test_sp->support( col_ind[i] ) )
            test_triangles.insert( triangle );
    }// for

    //
    // loop over each pair of triangles in ansatz and test space
    // and integrate kernel function
    //

    std::vector< value_t >  kernel_values;
    std::vector< T3Point >  phi0_values;
    std::vector< T3Point >  phi1_values;
    std::vector< T3Point >  quad_points0;
    std::vector< T3Point >  quad_points1;

    B::fill( value_t(0), values );

    for ( auto  tri0idx : ansatz_triangles )
    {
        TGrid::triangle_t  tri0    = ansatz_grid->triangle( tri0idx );
        const real_t       J0      = real_t( ansatz_grid->tri_size( tri0idx ) );

        for ( auto  tri1idx : test_triangles )
        {
            // integral vanishes when integrating over identical triangles
            if ( tri1idx == tri0idx )
                continue;

            TGrid::triangle_t  tri1 = test_grid->triangle( tri1idx );
            const real_t       J1   = real_t( test_sp->grid()->tri_size( tri1idx ) );
            const T3Point      n1   = test_grid->tri_normal( tri1idx );

            //
            // decide about relative position of triangles, choose quadrature points
            // and determine points of triangles s.t. order is as expected
            //

            uint  ncommon = this->reorder_common( tri0.vtx, tri1.vtx );
            uint  torder  = this->_quad_order;

            if ( ncommon == 0 )
                torder = this->adjust_order( tri0.vtx, tri1.vtx, torder );

            auto  rule    = this->quad_rule( ncommon, torder );
            auto  weights = rule -> template weights< real_t >();

            if ( rule == nullptr )
                HERROR( ERR_NULL, "(TMaxwellMFIEBF) eval", "quadrature rule is nullptr" );

            //
            // compute kernel at quadrature points
            //

            const size_t  npts        = rule->npts;
            const size_t  padded_npts = CFG::Mach::simd_padded_size< value_t >( npts );

            if ( kernel_values.size() < padded_npts )
            {
                kernel_values.resize( padded_npts );
                phi0_values.resize( padded_npts );
                phi1_values.resize( padded_npts );
                quad_points0.resize( padded_npts );
                quad_points1.resize( padded_npts );
            }// if

            eval_kernel( tri0idx, tri1idx, tri0, tri1, rule, kernel_values );

            ansatz_sp->ref_points2triangle( tri0, rule->pts1, quad_points0 );
            test_sp->ref_points2triangle( tri1, rule->pts2, quad_points1 );
            
            //
            // compute final result for each basis function by combining
            // kernel values with precomputed basis function values
            //

            for ( auto  idx0 : ansatz_sp->triangle_indices( tri0idx ) )
            {
                // exclude indices not in <rowind>
                if ( rowmap.find( idx0 ) == rowmap.end() )
                    continue;

                const idx_t   pos0         = rowmap[ idx0 ];
                const bool    phi0_in_supp = ansatz_sp->eval_basis( idx0, tri0, quad_points0, phi0_values );
                const real_t  scale_phi0   = ansatz_sp->scaling_factor_basis( idx0, tri0idx );

                for ( auto  idx1 : test_sp->triangle_indices( tri1idx ) )
                {
                    // exclude indices not in <colind>
                    if ( colmap.find( idx1 ) == colmap.end() )
                        continue;

                    const idx_t   pos1         = colmap[ idx1 ];
                    const bool    phi1_in_supp = test_sp->eval_basis( idx1, tri1, quad_points1, phi1_values );
                    const real_t  scale_phi1   = test_sp->scaling_factor_basis( idx1, tri1idx );
                    
                    value_t       value        = value_t(0);

                    if ( phi0_in_supp && phi1_in_supp )
                    {
                        for ( size_t  iv = 0; iv < npts; ++iv )
                        {
                            const T3Point  gradg_x_phi0 = cross( quad_points1[iv] - quad_points0[iv], phi0_values[iv] );
                            const T3Point  n1_x_phi1    = cross( n1, phi1_values[iv] );
                            
                            value += real_t( weights[iv] ) * kernel_values[ iv ] * real_t( dot( gradg_x_phi0, n1_x_phi1 ) );
                        }// for
                    }// if
                    
                    values( pos0, pos1 ) += value * scale_phi0 * scale_phi1 * J0 * J1;
                }// for
            }// for
        }// for
    }// for
}

//
// eval kernel function at quadrature points
//
template < typename  T_ansatzsp,
           typename  T_testsp >
void
TMaxwellMFIEBF< T_ansatzsp, T_testsp >::eval_kernel ( const                        idx_t,
                                                      const                        idx_t,
                                                      const TGrid::triangle_t &    tri0,
                                                      const TGrid::triangle_t &    tri1,
                                                      const tripair_quad_rule_t *  rule,
                                                      std::vector< value_t > &     values,
                                                      const real_t                 eps ) const
{
    _kernel_fn( tri0, tri1, rule, _ikappa, this->ansatz_space(), this->test_space(), values );
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// TMaxwellMFIEMassBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//
// ctor
//
template < typename T_ansatzsp, typename T_testsp >
TMaxwellMFIEMassBF< T_ansatzsp, T_testsp >::TMaxwellMFIEMassBF ( const ansatzsp_t *  aansatzsp,
                                                                 const testsp_t *    atestsp,
                                                                 const uint          aorder )
        : TBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp )
{
    //
    // build quad-points and weights
    //

    TTriGaussQuad  gaussquad;

    gaussquad.build( std::max( 1u, aorder ), _quad_pts, _quad_wghts );
}
    
template < typename  T_ansatzsp,
           typename  T_testsp >
matform_t
TMaxwellMFIEMassBF< T_ansatzsp, T_testsp >::format () const
{
    return unsymmetric;
}

//
// evaluate bilinearform
//

template < typename T_ansatzsp, typename T_testsp >
void
TMaxwellMFIEMassBF< T_ansatzsp, T_testsp >::eval ( const std::vector< idx_t > &  row_ind,
                                                   const std::vector< idx_t > &  col_ind,
                                                   BLAS::Matrix< value_t > &     values ) const
{
    const size_t  nrows = row_ind.size();
    const size_t  ncols = col_ind.size();

    for ( size_t  col = 0; col < ncols; col++ )
    {
        const idx_t  j = col_ind[col];
        
        for ( size_t  row = 0; row < nrows; row++ )
        {
            const idx_t  i = row_ind[row];
        
            //
            // integrate over test space
            //

            auto                    value     = value_t(0);
            auto                    ansatz_sp = this->ansatz_space();
            auto                    test_sp   = this->test_space();
            auto                    support_i = ansatz_sp->support(i);
            auto                    support_j = test_sp->support(j);
            std::vector< T3Point >  phi0_values;
            std::vector< T3Point >  phi1_values;

            if ( ansatz_sp->grid() != test_sp->grid() )
                HERROR( ERR_NOT_IMPL, "(TMaxwellMFIEMassBF) eval", "unequal grids not supported" );

            for ( auto  tri1 : support_j )
            {
                const auto  J1 = value_t( test_sp->grid()->tri_size( tri1 ) );

                //
                // test if support of ansatz function is disjoint
                // and immediately go one if so
                //

                bool  intersects = false;

                for ( auto  tri0 : support_i )
                {
                    if ( tri0 == tri1 ) // only valid for equal grids!
                    {
                        intersects = true;
                        break;
                    }// if
                }// for

                if ( ! intersects )
                    continue;

                //
                // evaluate basis functions at quadrature points
                //

                TGrid::triangle_t  t1;
                const size_t       npts   = _quad_pts.size();
                auto               tvalue = value_t(0);

                if ( phi0_values.size() < npts )
                {
                    phi0_values.resize( npts );
                    phi1_values.resize( npts );
                }// if

                t1 = test_sp->grid()->triangle( tri1 );

                ansatz_sp->eval_basis( i, t1, _quad_pts, phi0_values );
                test_sp->eval_basis(   j, t1, _quad_pts, phi1_values );

                const value_t  scale_phi0 = ansatz_sp->scaling_factor_basis( i, tri1 );
                const value_t  scale_phi1 = test_sp->scaling_factor_basis(   j, tri1 );
                const value_t  scale_phi  = scale_phi0 * scale_phi1;
                const T3Point  n1         = test_sp->grid()->tri_normal( tri1 );

                for ( size_t  k = 0; k < npts; ++k )
                {
                    const T3Point  n1_x_phi0 = cross( n1, phi0_values[k] );
                    const T3Point  n1_x_phi1 = cross( n1, phi1_values[k] );

                    tvalue += value_t( _quad_wghts[k] * dot( n1_x_phi0, n1_x_phi1 ) );
                }// for
                tvalue *= scale_phi;

                value += J1 * tvalue;
            }// for

            values( idx_t( row ), idx_t( col ) ) = value;
        }// for
    }// for
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instantiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template class TMaxwellEFIEBF< TConstEdgeFnSpace, TConstEdgeFnSpace >;
template class TMaxwellEFIEMassBF< TConstEdgeFnSpace, TConstEdgeFnSpace >;

template class TMaxwellMFIEBF< TConstEdgeFnSpace, TConstEdgeFnSpace >;
template class TMaxwellMFIEMassBF< TConstEdgeFnSpace, TConstEdgeFnSpace >;

}// namespace
