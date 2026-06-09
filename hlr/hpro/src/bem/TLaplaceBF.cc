//
// Project     : HLIBpro
// File        : TLaplaceBF.cc
// Description : bilinearforms for Laplace operator
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/base/packed.hh"
#include "hpro/blas/lapack.hh"

#include "hpro/bem/TLaplaceBF.hh"

/////////////////////////////////////////////////////////////

namespace Hpro
{

namespace B = BLAS;

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// external functions
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename value_t,
           typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
laplace_slp_simd  ( const TGrid::triangle_t &    tri0,
                    const TGrid::triangle_t &    tri1,
                    const tripair_quad_rule_t *  rule,
                    std::vector< value_t > &     values,
                    const T_ansatzsp *           ansatz_sp,
                    const T_testsp *             test_sp );

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
laplace_dlp_simd  ( const idx_t                                tri_id,
                    const bool                                 adjoint,
                    const TGrid::triangle_t &                  tri0,
                    const TGrid::triangle_t &                  tri1,
                    const tripair_quad_rule_t *                rule,
                    std::vector< value_type_t< T_packed > > &  values,
                    const T_ansatzsp *                         ansatz_sp,
                    const T_testsp *                           test_sp );

template < typename T_packed >
void
laplace_slp_eval_dx_simd  ( const tri_quad_rule_t< value_type_t< T_packed > > & quad_rule,
                            const T3Point                                       x[3],
                            const T3Point &                                     y,
                            std::vector< value_type_t< T_packed > > &           values );

template < typename T_packed >
void
laplace_dlp_eval_dy_simd  ( const tri_quad_rule_t< value_type_t< T_packed > > & quad_rule,
                            const T3Point &                                     x,
                            const T3Point                                       y[3],
                            const T3Point &                                     normal,
                            std::vector< value_type_t< T_packed > > &           values );

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
// Laplace SLP
//
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// machine dependent implementations
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename T_ansatzsp,
           typename T_testsp,
           typename value_t >
void
laplace_slp_flt ( const TGrid::triangle_t &    tri0,
                  const TGrid::triangle_t &    tri1,
                  const tripair_quad_rule_t *  rule,
                  std::vector< value_t > &     values,
                  const T_ansatzsp *           ansatz_sp,
                  const T_testsp *             test_sp )
{
    using real_t  = real_type_t< value_t >;
    
    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const auto     n_pts = rule -> npts;
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
    T3Point        x, y, u;

    #pragma omp simd
    for ( size_t  i = 0; i < n_pts; ++i )
    {
        const value_t  a1 = x1[i];
        const value_t  a2 = y1[i];
        const value_t  a0 = value_t(1) - a1 - a2;
        
        const value_t  b1 = x2[i];
        const value_t  b2 = y2[i];
        const value_t  b0 = value_t(1) - b1 - b2;

        // u = x - y
        x = a0 * v00 + a1 * v01 + a2 * v02;
        y = b0 * v10 + b1 * v11 + b2 * v12;
        u = x - y;

        //
        //   1
        // ─────
        // |x-y|
        //
        
        values[ i ] = ( ONE_OVER_4PI * Math::rsqrt( value_t( dot( u, u ) ) ) );
    }// for
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// TLaplaceSLPBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template < typename  ansatzsp_t,
           typename  testsp_t,
           typename  value_t >
auto
laplace_slp_kernel_fn ()
{
    using  real_t = real_type_t< value_t >;
    
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f()  && CFG::BEM::use_simd_avx512f  )
        {
            HINFO( "(TLaplaceSLPBF) using AVX512F kernel" );
            
            return laplace_slp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_AVX512F > >;
        }// if
        else if ( CFG::Mach::has_mic()  && CFG::BEM::use_simd_mic  )
        {
            HINFO( "(TLaplaceSLPBF) using MIC kernel" );
            
            return laplace_slp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_MIC > >;
        }// if
        else if ( CFG::Mach::has_avx2()  && CFG::BEM::use_simd_avx2  )
        {
            HINFO( "(TLaplaceSLPBF) using AVX2 kernel" );
            
            return laplace_slp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_AVX2 > >;
        }// if
        else if ( CFG::Mach::has_avx()  && CFG::BEM::use_simd_avx  )
        {
            HINFO( "(TLaplaceSLPBF) using AVX kernel" );
            
            return laplace_slp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_AVX > >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            HINFO( "(TLaplaceSLPBF) using SSE3 kernel" );
            
            return laplace_slp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_SSE3 > >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            HINFO( "(TLaplaceSLPBF) using VSX kernel" );
            
            return laplace_slp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_VSX > >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            HINFO( "(TLaplaceSLPBF) using NEON kernel" );
            
            return laplace_slp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_NEON > >;
        }// if
        else
        {
            HINFO( "(TLaplaceSLPBF) using standard kernel" );
            
            return laplace_slp_flt< ansatzsp_t, testsp_t, value_t >;
        }// else
    }// if
    else
    {
        HINFO( "(TLaplaceSLPBF) using standard kernel" );
            
        return laplace_slp_flt< ansatzsp_t, testsp_t, value_t >;
    }// else
}

template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
TLaplaceSLPBF< T_ansatzsp, T_testsp, T_value >::TLaplaceSLPBF ( const ansatzsp_t *  aansatzsp,
                                                                const testsp_t *    atestsp,
                                                                const uint          quad_order )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, quad_order, false )
        , _kernel_fn_dbl( laplace_slp_kernel_fn< T_ansatzsp, T_testsp, double >() )
        , _kernel_fn_sgl( laplace_slp_kernel_fn< T_ansatzsp, T_testsp, float >() )
{}

template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
TLaplaceSLPBF< T_ansatzsp, T_testsp, T_value >::TLaplaceSLPBF ( const ansatzsp_t *  aansatzsp,
                                                                const testsp_t *    atestsp,
                                                                const real_t        quad_error )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, 10, true, quad_error )
        , _kernel_fn_dbl( laplace_slp_kernel_fn< T_ansatzsp, T_testsp, double >() )
        , _kernel_fn_sgl( laplace_slp_kernel_fn< T_ansatzsp, T_testsp, float >() )
{}

//
// eval kernel function at quadrature points
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
void
TLaplaceSLPBF< T_ansatzsp, T_testsp, T_value >::eval_kernel ( const                        idx_t,
                                                              const                        idx_t,
                                                              const TGrid::triangle_t &    tri0,
                                                              const TGrid::triangle_t &    tri1,
                                                              const tripair_quad_rule_t *  rule,
                                                              std::vector< value_t > &     values,
                                                              const real_t                 eps ) const
{
    if constexpr ( std::is_same_v< value_t, double > )
    {
        if (( eps > 1e-7 ) && ( rule->npts >= 50 ))
        {
            //
            // evaluate in FP32 and convert to FP64
            //

            auto  fp32_values = std::vector< float >( values.size() );

            _kernel_fn_sgl( tri0, tri1, rule, fp32_values, this->ansatz_space(), this->test_space() );

            std::copy( fp32_values.cbegin(), fp32_values.cend(), values.begin() );
        }// if
        else
        {
            _kernel_fn_dbl( tri0, tri1, rule, values, this->ansatz_space(), this->test_space() );
        }// else
    }// if

    if constexpr ( std::is_same_v< value_t, float > )
        _kernel_fn_sgl( tri0, tri1, rule, values, this->ansatz_space(), this->test_space() );
    
    ADD_FLOPS( rule->npts * ( 15 + 15 + 3 + 8 ) );
}

//
// return format of bilinear form, e.g. symmetric
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
matform_t
TLaplaceSLPBF< T_ansatzsp, T_testsp, T_value >::format () const
{
    if ( ( reinterpret_cast< const void * >( this->ansatz_space() ) == reinterpret_cast< const void * >( this->test_space() ) ) ||
         (( this->ansatz_space()->type() == this->test_space()->type() ) &&
          ( this->ansatz_space()->grid() == this->test_space()->grid() )) )
        return symmetric;
    else
        return unsymmetric;
}



/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
// Laplace DLP
//
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// machine dependent implementations
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename T_ansatzsp,
           typename T_testsp,
           typename value_t >
void
laplace_dlp_flt ( const idx_t                  tri_id,
                  const bool                   adjoint,
                  const TGrid::triangle_t &    tri0,
                  const TGrid::triangle_t &    tri1,
                  const tripair_quad_rule_t *  rule,
                  std::vector< value_t > &     values,
                  const T_ansatzsp *           ansatz_sp,
                  const T_testsp *             test_sp )
{
    using real_t  = real_type_t< value_t >;
    
    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const size_t   n_pts = rule -> npts;
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
    T3Point        x, y, u, n;
    const bool     has_vtx_normal = ( adjoint ? ansatz_sp->grid()->has_vtx_normal() : test_sp->grid()->has_vtx_normal() );

    if ( ! has_vtx_normal )
        n = ( adjoint ? ansatz_sp->grid()->tri_normal( tri_id ) : test_sp->grid()->tri_normal( tri_id ) );
            
    #pragma omp simd
    for ( size_t  i = 0; i < n_pts; i++ )
    {
        const value_t  a1 = x1[i];
        const value_t  a2 = y1[i];
        const value_t  a0 = 1.0 - a1 - a2;
        
        const value_t  b1 = x2[i];
        const value_t  b2 = y2[i];
        const value_t  b0 = 1.0 - b1 - b2;

        // u = x - y
        x = a0 * v00 + a1 * v01 + a2 * v02;
        y = b0 * v10 + b1 * v11 + b2 * v12;

        u = ( adjoint ? x - y : y - x );

        //
        // ⟨n,x-y⟩
        // ───────
        //  |x-y|³
        //

        if ( has_vtx_normal )
            n = ( adjoint
                  ? ansatz_sp->grid()->tri_normal( tri_id, tri0, a1, a2 )
                  : test_sp->grid()->tri_normal( tri_id, tri1, b1, b2 ) );
        
        const auto  inxy = Math::rsqrt( value_t( dot( u, u ) ) );
        
        values[ i ] = ONE_OVER_4PI * value_t( dot( n, u ) ) * inxy * inxy * inxy;
    }// for
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// TLaplaceDLPBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template < typename  ansatzsp_t,
           typename  testsp_t,
           typename  value_t >
auto
laplace_dlp_kernel_fn ()
{
    using  real_t = real_type_t< value_t >;
    
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f()  && CFG::BEM::use_simd_avx512f  )
        {
            HINFO( "(TLaplaceDLPBF) using AVX512F kernel" );
            
            return laplace_dlp_simd< ansatzsp_t, testsp_t, packed< real_t, ISA_AVX512F > >;
        }// if
        else if ( CFG::Mach::has_mic()  && CFG::BEM::use_simd_mic  )
        {
            HINFO( "(TLaplaceDLPBF) using MIC kernel" );
            
            return laplace_dlp_simd< ansatzsp_t, testsp_t, packed< real_t, ISA_MIC > >;
        }// if
        else if ( CFG::Mach::has_avx2()  && CFG::BEM::use_simd_avx2  )
        {
            HINFO( "(TLaplaceDLPBF) using AVX2 kernel" );
            
            return laplace_dlp_simd< ansatzsp_t, testsp_t, packed< real_t, ISA_AVX2 > >;
        }// if
        else if ( CFG::Mach::has_avx()  && CFG::BEM::use_simd_avx  )
        {
            HINFO( "(TLaplaceDLPBF) using AVX kernel" );
            
            return laplace_dlp_simd< ansatzsp_t, testsp_t, packed< real_t, ISA_AVX > >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            HINFO( "(TLaplaceDLPBF) using SSE3 kernel" );
            
            return laplace_dlp_simd< ansatzsp_t, testsp_t, packed< real_t, ISA_SSE3 > >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            HINFO( "(TLaplaceDLPBF) using VSX kernel" );
            
            return laplace_dlp_simd< ansatzsp_t, testsp_t, packed< real_t, ISA_VSX > >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            HINFO( "(TLaplaceDLPBF) using NEON kernel" );
            
            return laplace_dlp_simd< ansatzsp_t, testsp_t, packed< real_t, ISA_NEON > >;
        }// if
        else
        {
            HINFO( "(TLaplaceDLPBF) using FPU kernel" );
            
            return laplace_dlp_flt< ansatzsp_t, testsp_t, value_t >;
        }// else
    }// if
    else
    {
        HINFO( "(TLaplaceDLPBF) using FPU kernel" );
            
        return laplace_dlp_flt< ansatzsp_t, testsp_t, value_t >;
    }// else
}

//
// ctor
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
TLaplaceDLPBF< T_ansatzsp, T_testsp, T_value >::TLaplaceDLPBF ( const ansatzsp_t *  aansatzsp,
                                                                const testsp_t *    atestsp,
                                                                const bool          adjoint,
                                                                const uint          quad_order )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, quad_order )
        , _adjoint( adjoint )
        , _kernel_fn( laplace_dlp_kernel_fn< T_ansatzsp, T_testsp, T_value >() )
{}

template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
TLaplaceDLPBF< T_ansatzsp, T_testsp, T_value >::TLaplaceDLPBF ( const ansatzsp_t *  aansatzsp,
                                                                const testsp_t *    atestsp,
                                                                const bool          adjoint,
                                                                const real_t        quad_error )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, 10, true, quad_error )
        , _adjoint( adjoint )
        , _kernel_fn( laplace_dlp_kernel_fn< T_ansatzsp, T_testsp, T_value >() )
{}

//
// eval kernel function at quadrature points
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
void
TLaplaceDLPBF< T_ansatzsp, T_testsp, T_value >::eval_kernel ( const idx_t                  tri0idx,
                                                              const idx_t                  tri1idx,
                                                              const TGrid::triangle_t &    tri0,
                                                              const TGrid::triangle_t &    tri1,
                                                              const tripair_quad_rule_t *  rule,
                                                              std::vector< value_t > &     values,
                                                              const real_t                 eps ) const
{
    if ( _adjoint )
        _kernel_fn( tri0idx, _adjoint, tri0, tri1, rule, values, this->ansatz_space(), this->test_space() );
    else
        _kernel_fn( tri1idx, _adjoint, tri0, tri1, rule, values, this->ansatz_space(), this->test_space() );
    
    ADD_FLOPS( rule->npts * ( 15 + 15 + 3 + 7 + 9 ) );
}


//
// return format of bilinear form, e.g. symmetric
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
matform_t
TLaplaceDLPBF< T_ansatzsp, T_testsp, T_value >::format () const
{
    return unsymmetric;
}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// auxiliary functions for SLP and DLP generators
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// evaluate dx-integral of SLP
//
template < typename value_t >
void
laplace_slp_eval_dx_flt ( const tri_quad_rule_t< value_t > & quad_rule,
                          const T3Point                      vx[3],
                          const T3Point &                    vy,
                          std::vector< value_t > &           values )
{
    using real_t  = real_type_t< value_t >;
    
    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const size_t  npts  = quad_rule.npts;
    
    #pragma omp simd
    for ( size_t  k = 0; k < npts; ++k )
    {
        const value_t  a2 = quad_rule.y[k];
        const value_t  a1 = quad_rule.x[k];
        const value_t  a0 = 1.0 - a1 - a2;

        //
        // evaluate SLP kernel
        //

        // d = y - x
        const T3Point  x = a0 * vx[0]  +  a1 * vx[1]  +  a2 * vx[2];
        const T3Point  d = vy - x;

        values[k] = ONE_OVER_4PI * Math::rsqrt( value_t( dot( d, d ) ) );
    }// for
}

//
// evaluate dy-integral of DLP
//
template < typename value_t >
void
laplace_dlp_eval_dy_flt ( const tri_quad_rule_t< value_t > & quad_rule,
                          const T3Point &                    vx,
                          const T3Point                      vy[3],
                          const T3Point &                    normal,
                          std::vector< value_t > &           values )
{
    using real_t  = real_type_t< value_t >;
    
    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());

    const size_t  npts  = quad_rule.npts;
    
    //
    // fall back or handle remaining
    //

    #pragma omp simd
    for ( size_t  k = 0; k < npts; k++ )
    {
        const value_t  a2 = quad_rule.y[k];
        const value_t  a1 = quad_rule.x[k];
        const value_t  a0 = 1.0 - a1 - a2;

        //
        // evaluate DLP kernel
        //
                
        // d = y - x
        const T3Point  y = a0 * vy[0]  +  a1 * vy[1]  +  a2 * vy[2];
        const T3Point  d = y - vx;

        // r = 1 / |d| = 1 / sqrt( <d,d> )
        const auto     r = Math::rsqrt( value_t( d.dot( d ) ) );
                    
        // store < n, d > / |d|^3
        values[k] = ONE_OVER_4PI * value_t( normal.dot( d ) ) * (r*r*r);
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// TLaplaceSLPGenFn
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// constructor
//
template < typename T_ansatzsp,
           typename T_testsp,
           typename T_value >
TLaplaceSLPGenFn< T_ansatzsp, T_testsp, T_value >::TLaplaceSLPGenFn ( const ansatzsp_t *    ansatzsp,
                                                                      const testsp_t *      testsp,
                                                                      const TPermutation *  row_perm_i2e,
                                                                      const TPermutation *  col_perm_i2e,
                                                                      const uint            quad_order )
        : TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, value_t >( ansatzsp,
                                                                    testsp,
                                                                    quad_order,
                                                                    row_perm_i2e,
                                                                    col_perm_i2e )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            HINFO( "(TLaplaceSLPGenFn) using AVX512F kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_AVX512F > >;
        }// if
        else if ( CFG::Mach::has_mic()  && CFG::BEM::use_simd_mic  )
        {
            HINFO( "(TLaplaceSLPGenFn) using MIC kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_MIC > >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            HINFO( "(TLaplaceSLPGenFn) using AVX2 kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_AVX2 > >;
        }// if
        else if ( CFG::Mach::has_avx()  && CFG::BEM::use_simd_avx  )
        {
            HINFO( "(TLaplaceSLPGenFn) using AVX kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_AVX > >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            HINFO( "(TLaplaceSLPGenFn) using SSE3 kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_SSE3 > >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            HINFO( "(TLaplaceSLPGenFn) using VSX kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_VSX > >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            HINFO( "(TLaplaceSLPGenFn) using NEON kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_NEON > >;
        }// if
        else
        {
            HINFO( "(TLaplaceSLPGenFn) using FPU kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_flt< value_t >;
        }// else
    }// if
    else
    {
        HINFO( "(TLaplaceSLPGenFn) using FPU kernel" );
            
        _eval_dxy_impl = laplace_slp_eval_dx_flt< value_t >;
    }// else
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// TLaplaceDLPGenFn
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// constructor
//
template < typename T_ansatzsp,
           typename T_testsp,
           typename T_value >
TLaplaceDLPGenFn< T_ansatzsp, T_testsp, T_value >::TLaplaceDLPGenFn ( const ansatzsp_t *    ansatzsp,
                                                                      const testsp_t *      testsp,
                                                                      const TPermutation *  row_perm_i2e,
                                                                      const TPermutation *  col_perm_i2e,
                                                                      const uint            quad_order )
        : TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, value_t >( ansatzsp,
                                                                    testsp,
                                                                    quad_order,
                                                                    row_perm_i2e,
                                                                    col_perm_i2e )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            HINFO( "(TLaplaceDLPGenFn) using AVX512F kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_AVX512F > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real_t, ISA_AVX512F > >;
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic )
        {
            HINFO( "(TLaplaceDLPGenFn) using MIC kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_MIC > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real_t, ISA_MIC > >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            HINFO( "(TLaplaceDLPGenFn) using AVX2 kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_AVX2 > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real_t, ISA_AVX2 > >;
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            HINFO( "(TLaplaceDLPGenFn) using AVX kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_AVX > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real_t, ISA_AVX > >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            HINFO( "(TLaplaceDLPGenFn) using SSE3 kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_SSE3 > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real_t, ISA_SSE3 > >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            HINFO( "(TLaplaceDLPGenFn) using VSX kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_VSX > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real_t, ISA_VSX > >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            HINFO( "(TLaplaceDLPGenFn) using NEON kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real_t, ISA_NEON > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real_t, ISA_NEON > >;
        }// if
        else
        {
            HINFO( "(TLaplaceDLPGenFn) using standard kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_flt< value_t >;
            _eval_dy_impl = laplace_dlp_eval_dy_flt< value_t >;
        }// else
    }// if
    else
    {
        HINFO( "(TLaplaceDLPGenFn) using standard kernel" );
            
        _eval_dx_impl = laplace_slp_eval_dx_flt< value_t >;
        _eval_dy_impl = laplace_dlp_eval_dy_flt< value_t >;
    }// if
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

#define INST_ALL( type )                                                \
    template class TLaplaceSLPBF< TConstFnSpace< type >,  TConstFnSpace< type >, type >; \
    template class TLaplaceSLPBF< TConstFnSpace< type >,  TLinearFnSpace< type >, type >; \
    template class TLaplaceSLPBF< TLinearFnSpace< type >, TConstFnSpace< type >, type >; \
    template class TLaplaceSLPBF< TLinearFnSpace< type >, TLinearFnSpace< type >, type >; \
                                                                        \
    template class TLaplaceDLPBF< TConstFnSpace< type >,  TConstFnSpace< type >, type >; \
    template class TLaplaceDLPBF< TConstFnSpace< type >,  TLinearFnSpace< type >, type >; \
    template class TLaplaceDLPBF< TLinearFnSpace< type >, TConstFnSpace< type >, type >; \
    template class TLaplaceDLPBF< TLinearFnSpace< type >, TLinearFnSpace< type >, type >; \
                                                                        \
    template class TLaplaceSLPGenFn< TConstFnSpace< type >,  TConstFnSpace< type >, type >; \
    template class TLaplaceSLPGenFn< TConstFnSpace< type >,  TLinearFnSpace< type >, type >; \
    template class TLaplaceSLPGenFn< TLinearFnSpace< type >, TConstFnSpace< type >, type >; \
    template class TLaplaceSLPGenFn< TLinearFnSpace< type >, TLinearFnSpace< type >, type >; \
                                                                        \
    template class TLaplaceDLPGenFn< TConstFnSpace< type >,  TConstFnSpace< type >, type >; \
    template class TLaplaceDLPGenFn< TConstFnSpace< type >,  TLinearFnSpace< type >, type >; \
    template class TLaplaceDLPGenFn< TLinearFnSpace< type >, TConstFnSpace< type >, type >; \
    template class TLaplaceDLPGenFn< TLinearFnSpace< type >, TLinearFnSpace< type >, type >;

INST_ALL( float )
INST_ALL( double )

}// namespace Hpro
