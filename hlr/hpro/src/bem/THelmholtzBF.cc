//
// Project     : HLIBpro
// File        : THelmholtzBF.cc
// Description : bilinear forms for Helmholtz operator
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>
#include <iostream> // DEBUG

#include "hpro/base/packed.hh"

#include "hpro/bem/TConstEdgeFnSpace.hh"

#include "hpro/bem/THelmholtzBF.hh"

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

#define  HELMHOLTZ_DLP( suffix )                                        \
    template < typename  value_t,                                       \
               typename  T_ansatzsp,                                    \
               typename  T_testsp,                                      \
               typename  T_packed >                                     \
    void                                                                \
    helmholtz_dlp_##suffix ( const idx_t                  tri1_id,      \
                             const bool                   adjoint,      \
                             const TGrid::triangle_t &    tri0,         \
                             const TGrid::triangle_t &    tri1,         \
                             const tripair_quad_rule_t *  rule,         \
                             const value_t                ikappa,       \
                             const T_ansatzsp *           ansatz_sp,    \
                             const T_testsp *             test_sp,      \
                             std::vector< value_t > &     values )

HELMHOLTZ_DLP( simd );
HELMHOLTZ_DLP( re_simd );
HELMHOLTZ_DLP( im_simd );

#define  HELMHOLTZ_SLP_HCA( suffix )                                    \
    template < typename value_t, typename  T_packed >                   \
    void                                                                \
    helmholtz_slp_eval_dx_##suffix ( const value_t &           ikappa,  \
                                     const tri_quad_rule_t< real_type_t< value_t > > &  quad_rule, \
                                     const T3Point             x[3],    \
                                     const T3Point &           y,       \
                                     std::vector< value_t > &  values )

HELMHOLTZ_SLP_HCA( simd );
HELMHOLTZ_SLP_HCA( re_simd );
HELMHOLTZ_SLP_HCA( im_simd );

#define  HELMHOLTZ_DLP_HCA( suffix )                                    \
    template < typename value_t, typename  T_packed >                   \
    void                                                                \
    helmholtz_dlp_eval_dy_##suffix ( const value_t &           ikappa, \
                                     const tri_quad_rule_t< real_type_t< value_t > > &  quad_rule, \
                                     const T3Point &           x, \
                                     const T3Point             y[3], \
                                     const T3Point &           normal, \
                                     std::vector< value_t > &  values )

HELMHOLTZ_DLP_HCA( simd );
HELMHOLTZ_DLP_HCA( re_simd );
HELMHOLTZ_DLP_HCA( im_simd );

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// THelmholtzSLPBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//
// standard evaluation of Helmholtz SLP kernel
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
                    std::vector< value_t > &     values )
{
    using real_t  = real_type_t< value_t >;

    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
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

        // transform coordinates
        const T3Point  x    = a0 * v00 + a1 * v01 + a2 * v02;
        const T3Point  y    = b0 * v10 + b1 * v11 + b2 * v12;
        const auto     xmy  = x - y;
        const auto     dist = real_t( dot( xmy, xmy ) );

        // e^( i * kappa * |x-y| ) / |x-y|
        values[ i ] = ( ONE_OVER_4PI * Math::exp( ikappa * Math::sqrt( dist ) ) * Math::rsqrt( dist ) );
    }// for
}

template < typename  ansatzsp_t,
           typename  testsp_t,
           typename  value_t >
auto
helmholtz_slp_kernel_fn ( const value_t  ikappa )
{
    using  real_t = real_type_t< value_t >;

    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            using  packed_t = packed< real_t, ISA_AVX512F >;
            
            HINFO( "(THelmholtzSLPBF) using AVX512F kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic )
        {
            using  packed_t = packed< real_t, ISA_MIC >;
            
            HINFO( "(THelmholtzSLPBF) using MIC kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            using  packed_t = packed< real_t, ISA_AVX2 >;
            
            HINFO( "(THelmholtzSLPBF) using AVX2 kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            using  packed_t = packed< real_t, ISA_AVX >;
            
            HINFO( "(THelmholtzSLPBF) using AVX kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            using  packed_t = packed< real_t, ISA_SSE3 >;
            
            HINFO( "(THelmholtzSLPBF) using SSE3 kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            using  packed_t = packed< real_t, ISA_VSX >;
            
            HINFO( "(THelmholtzSLPBF) using VSX kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            using  packed_t = packed< real_t, ISA_NEON >;
            
            HINFO( "(THelmholtzSLPBF) using NEON kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_slp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_slp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_slp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else
        {
            HINFO( "(THelmholtzSLPBF) using standard kernel" );
            
            return helmholtz_slp_flt< ansatzsp_t, testsp_t, value_t >;
        }// else
    }// if
    else
    {
        HINFO( "(THelmholtzSLPBF) using standard kernel" );
            
        return helmholtz_slp_flt< ansatzsp_t, testsp_t, value_t >;
    }// else
}

//
// ctor
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
THelmholtzSLPBF< T_ansatzsp, T_testsp, T_value >::THelmholtzSLPBF ( const value_t       kappa,
                                                                    const ansatzsp_t *  aansatzsp,
                                                                    const testsp_t *    atestsp,
                                                                    const uint          quad_order )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, quad_order )
        , _ikappa( value_t( 0, 1 ) * kappa )
        , _kernel_fn( helmholtz_slp_kernel_fn< T_ansatzsp, T_testsp, T_value >( _ikappa ) )
{}

template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
THelmholtzSLPBF< T_ansatzsp, T_testsp, T_value >::THelmholtzSLPBF ( const value_t       kappa,
                                                                    const ansatzsp_t *  aansatzsp,
                                                                    const testsp_t *    atestsp,
                                                                    const real_t        quad_error )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, 10, true, quad_error )
        , _ikappa( value_t( 0, 1 ) * kappa )
        , _kernel_fn( helmholtz_slp_kernel_fn< T_ansatzsp, T_testsp, T_value >( _ikappa ) )
{}

//
// eval kernel function at quadrature points
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
void
THelmholtzSLPBF< T_ansatzsp, T_testsp, T_value >::eval_kernel ( const idx_t                  /* tri0idx */,
                                                                const idx_t                  /* tri1idx */,
                                                                const TGrid::triangle_t &    tri0,
                                                                const TGrid::triangle_t &    tri1,
                                                                const tripair_quad_rule_t *  rule,
                                                                std::vector< value_t > &     values,
                                                                const real_t                 eps ) const
{
    _kernel_fn( tri0, tri1, rule, _ikappa, this->ansatz_space(), this->test_space(), values );
}

//
// return format of bilinear form, e.g. symmetric
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
matform_t
THelmholtzSLPBF< T_ansatzsp, T_testsp, T_value >::format () const
{
    if ( ( reinterpret_cast< const void * >( this->ansatz_space() ) == reinterpret_cast< const void * >( this->test_space() ) ) ||
         (( this->ansatz_space()->type() == this->test_space()->type() ) &&
          ( this->ansatz_space()->grid() == this->test_space()->grid() )) )
        return symmetric;
    else
        return unsymmetric;
}




/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// THelmholtzDLPBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//
// standard evaluation of Helmholtz DLP kernel
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  value_t >
void
helmholtz_dlp_flt ( const idx_t                  tri_id,
                    const bool                   adjoint,
                    const TGrid::triangle_t &    tri0,
                    const TGrid::triangle_t &    tri1,
                    const tripair_quad_rule_t *  rule,
                    const value_t                ikappa,
                    const T_ansatzsp *           ansatz_sp,
                    const T_testsp *             test_sp,
                    std::vector< value_t > &     values )
{
    using real_t  = real_type_t< value_t >;

    constexpr real_t  ONE          = real_t(1);
    constexpr real_t  ONE_OVER_4PI = ONE / (real_t(4) * Math::pi< real_t >());

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
        const real_t  a1 = x1[i];
        const real_t  a2 = y1[i];
        const real_t  a0 = 1.0 - a1 - a2;
        
        const real_t  b1 = x2[i];
        const real_t  b2 = y2[i];
        const real_t  b0 = 1.0 - b1 - b2;

        // u = y - x
        const auto    x      = a0 * v00 + a1 * v01 + a2 * v02;
        const auto    y      = b0 * v10 + b1 * v11 + b2 * v12;
        const auto    ymx    = ( adjoint ? x - y : y - x );
        const auto    sqdist = real_t( dot( ymx, ymx ) );       // |y-x|²

        const auto    n      = ( adjoint
                                 ? ansatz_sp->grid()->tri_normal( tri_id, tri0, a1, a2 )
                                 : test_sp->grid()->tri_normal( tri_id, tri1, b1, b2 ) );
        
        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n,y-x>
        // ────────────────────────────────────
        //                |x-y|³
        //

        const real_t  dist      = Math::sqrt( sqdist );       // |y-x|
        const auto    cudist    = real_t( sqdist * dist );    // |y-x|³
        const auto    n_dot_ymx = real_t( dot( n, ymx ) );    // <n(y),y-x>
        const auto    ikymx     = ikappa * real_t( dist );    // i · κ · |x-y|

        if ( adjoint )
            values[ i ] = ONE_OVER_4PI * Math::exp( ikymx ) * ( ONE - ikymx ) * n_dot_ymx / cudist;
        else
            values[ i ] = ONE_OVER_4PI * Math::exp( ikymx ) * ( ikymx - ONE ) * n_dot_ymx / cudist;
    }// for
}

template < typename  ansatzsp_t,
           typename  testsp_t,
           typename  value_t >
auto
helmholtz_dlp_kernel_fn ( const value_t  ikappa )
{
    using  real_t = real_type_t< value_t >;

    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            using  packed_t = packed< real_t, ISA_AVX512F >;
            
            HINFO( "(THelmholtzDLPBF) using AVX512F kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic  )
        {
            using  packed_t = packed< real_t, ISA_MIC >;
            
            HINFO( "(THelmholtzDLPBF) using MIC kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            using  packed_t = packed< real_t, ISA_AVX2 >;
            
            HINFO( "(THelmholtzDLPBF) using AVX2 kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            using  packed_t = packed< real_t, ISA_AVX >;
            
            HINFO( "(THelmholtzDLPBF) using AVX kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            using  packed_t = packed< real_t, ISA_SSE3 >;
            
            HINFO( "(THelmholtzDLPBF) using SSE3 kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            using  packed_t = packed< real_t, ISA_VSX >;
            
            HINFO( "(THelmholtzDLPBF) using VSX kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            using  packed_t = packed< real_t, ISA_NEON >;
            
            HINFO( "(THelmholtzDLPBF) using NEON kernel" );
            
            if      ( std::real( ikappa ) == real_t(0) ) return helmholtz_dlp_re_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else if ( std::imag( ikappa ) == real_t(0) ) return helmholtz_dlp_im_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
            else                                         return helmholtz_dlp_simd< value_t, ansatzsp_t, testsp_t, packed_t >;
        }// if
        else
        {
            HINFO( "(THelmholtzDLPBF) using standard kernel" );
            
            return helmholtz_dlp_flt< ansatzsp_t, testsp_t, value_t >;
        }// else
    }// if
    else
    {
        HINFO( "(THelmholtzDLPBF) using standard kernel" );
            
        return helmholtz_dlp_flt< ansatzsp_t, testsp_t, value_t >;
    }// else
}

//
// ctor
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
THelmholtzDLPBF< T_ansatzsp, T_testsp, T_value >::THelmholtzDLPBF ( const value_t       kappa,
                                                                    const ansatzsp_t *  aansatzsp,
                                                                    const testsp_t *    atestsp,
                                                                    const bool          adjoint,
                                                                    const uint          quad_order )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, quad_order )
        , _adjoint( adjoint )
        , _ikappa( value_t( 0, 1 ) * kappa )
        , _kernel_fn( helmholtz_dlp_kernel_fn< T_ansatzsp, T_testsp, T_value >( _ikappa ) )
{}

template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
THelmholtzDLPBF< T_ansatzsp, T_testsp, T_value >::THelmholtzDLPBF ( const value_t       kappa,
                                                                    const ansatzsp_t *  aansatzsp,
                                                                    const testsp_t *    atestsp,
                                                                    const bool          adjoint,
                                                                    const real_t        quad_error )
: TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, 10, true, quad_error )
        , _adjoint( adjoint )
        , _ikappa( value_t( 0, 1 ) * kappa )
        , _kernel_fn( helmholtz_dlp_kernel_fn< T_ansatzsp, T_testsp, T_value >( _ikappa ) )
{}

//
// eval kernel function at quadrature points
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
void
THelmholtzDLPBF< T_ansatzsp, T_testsp, T_value >::eval_kernel ( const idx_t                  tri0idx,
                                                                const idx_t                  tri1idx,
                                                                const TGrid::triangle_t &    tri0,
                                                                const TGrid::triangle_t &    tri1,
                                                                const tripair_quad_rule_t *  rule,
                                                                std::vector< value_t > &     values,
                                                                const real_t                 eps ) const
{
    if ( _adjoint )
        _kernel_fn( tri0idx, _adjoint, tri0, tri1, rule, _ikappa, this->ansatz_space(), this->test_space(), values );
    else
        _kernel_fn( tri1idx, _adjoint, tri0, tri1, rule, _ikappa, this->ansatz_space(), this->test_space(), values );
}

//
// return format of bilinear form, e.g. symmetric
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
matform_t
THelmholtzDLPBF< T_ansatzsp, T_testsp, T_value >::format () const
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
// computes x-integral of SLP
//
template < typename value_t >
void
helmholtz_slp_eval_dx_flt ( const value_t &                                   ikappa,
                            const tri_quad_rule_t< real_type_t< value_t > > & quad_rule,
                            const T3Point                                     vx[3],
                            const T3Point &                                   vy,
                            std::vector< value_t > &                          values )
{
    using real_t  = real_type_t< value_t >;

    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());

    const size_t  npts = quad_rule.npts;
    
    #pragma omp simd
    for ( size_t  k = 0; k < npts; ++k )
    {
        const double  a2 = quad_rule.y[k];
        const double  a1 = quad_rule.x[k];
        const double  a0 = 1.0 - a1 - a2;

        //
        // evaluate SLP kernel
        //

        // d = y - x
        const T3Point  xmy( vy - ( a0 * vx[0] + a1 * vx[1] + a2 * vx[2] ) );
        const auto     ddotd = real_t( Math::sqrt( dot( xmy, xmy ) ) );

        // e^( i * kappa * | x - y | ) / | x - y |
        values[k] = ONE_OVER_4PI * Math::exp( ikappa * ddotd ) / ddotd;
    }// for
}

//
// computes y-integral of DLP
//
template < typename value_t >
void
helmholtz_dlp_eval_dy_flt ( const value_t &                                   ikappa,
                            const tri_quad_rule_t< real_type_t< value_t > > & quad_rule,
                            const T3Point &                                   vx,
                            const T3Point                                     vy[3],
                            const T3Point &                                   normal,
                            std::vector< value_t > &                          values )
{
    using real_t  = real_type_t< value_t >;

    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());

    const size_t  npts = quad_rule.npts;
    
    #pragma omp simd
    for ( size_t  k = 0; k < npts; ++k )
    {
        const double  a2 = quad_rule.y[k];
        const double  a1 = quad_rule.x[k];
        const double  a0 = 1.0 - a1 - a2;

        //
        // evaluate DLP kernel
        //
                
        // d = y - x
        const T3Point  ymx( ( a0 * vy[0] + a1 * vy[1] + a2 * vy[2] ) - vx );
        const auto     dist2  = real_t( dot( ymx, ymx ) );
        const auto     dist   = Math::sqrt( dist2 );
        const auto     dist3  = dist * dist2;
        const auto     ikdist = ikappa * dist;
                    
        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n,y-x>
        // ────────────────────────────────────
        //                |x-y|³
        //

        values[k] = ( ONE_OVER_4PI * ( Math::exp( ikdist ) * ( ikdist - real_t(1) ) * real_t(dot( normal, ymx )) )
                      / dist3 );
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// THelmholtzSLPGenFn
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// constructor
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
THelmholtzSLPGenFn< T_ansatzsp, T_testsp, T_value >::THelmholtzSLPGenFn ( const value_t         kappa,
                                                                          const ansatzsp_t *    ansatzsp,
                                                                          const testsp_t *      testsp,
                                                                          const TPermutation *  row_perm_i2e,
                                                                          const TPermutation *  col_perm_i2e,
                                                                          const uint            quad_order )
        : TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, value_t >( ansatzsp, testsp, quad_order, row_perm_i2e, col_perm_i2e )
        , _ikappa( value_t( 0, 1 ) * kappa )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            using  packed_t = packed< real_t, ISA_AVX512F >;
            
            HINFO( "(THelmholtzSLPGenFn) using AVX512F kernel" );
            
            if      ( std::real( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
            else if ( std::imag( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
            else                                          _eval_dxy_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic )
        {
            using  packed_t = packed< real_t, ISA_MIC >;
            
            HINFO( "(THelmholtzSLPGenFn) using MIC kernel" );
            
            if      ( std::real( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
            else if ( std::imag( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
            else                                          _eval_dxy_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            using  packed_t = packed< real_t, ISA_AVX2 >;
            
            HINFO( "(THelmholtzSLPGenFn) using AVX2 kernel" );
            
            if      ( std::real( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
            else if ( std::imag( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
            else                                          _eval_dxy_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            using  packed_t = packed< real_t, ISA_AVX >;
            
            HINFO( "(THelmholtzSLPGenFn) using AVX kernel" );
            
            if      ( std::real( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
            else if ( std::imag( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
            else                                          _eval_dxy_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            using  packed_t = packed< real_t, ISA_SSE3 >;
            
            HINFO( "(THelmholtzSLPGenFn) using SSE3 kernel" );
            
            if      ( std::real( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
            else if ( std::imag( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
            else                                          _eval_dxy_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            using  packed_t = packed< real_t, ISA_VSX >;
            
            HINFO( "(THelmholtzSLPGenFn) using VSX kernel" );
            
            if      ( std::real( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
            else if ( std::imag( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
            else                                          _eval_dxy_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            using  packed_t = packed< real_t, ISA_NEON >;
            
            HINFO( "(THelmholtzSLPGenFn) using NEON kernel" );
            
            if      ( std::real( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
            else if ( std::imag( _ikappa ) == real_t(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
            else                                          _eval_dxy_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
        }// if
        else
        {
            HINFO( "(THelmholtzSLPGenFn) using standard kernel" );
            
            _eval_dxy_impl = helmholtz_slp_eval_dx_flt< value_t >;
        }// else
    }// if
    else
    {
        HINFO( "(THelmholtzSLPGenFn) using standard kernel" );
            
        _eval_dxy_impl = helmholtz_slp_eval_dx_flt< value_t >;
    }// else
}
        
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// THelmholtzDLPGenFn
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// constructor
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
THelmholtzDLPGenFn< T_ansatzsp, T_testsp, T_value >::THelmholtzDLPGenFn ( const value_t         kappa,
                                                                          const ansatzsp_t *    ansatzsp,
                                                                          const testsp_t *      testsp,
                                                                          const TPermutation *  row_perm_i2e,
                                                                          const TPermutation *  col_perm_i2e,
                                                                          const uint            quad_order )
        : TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, value_t >( ansatzsp, testsp, quad_order, row_perm_i2e, col_perm_i2e )
        , _ikappa( value_t( 0, 1 ) * kappa )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            using  packed_t = packed< real_t, ISA_AVX512F >;
            
            HINFO( "(THelmholtzDLPGenFn) using AVX512F kernel" );
            
            if ( std::real( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< value_t, packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< value_t, packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< value_t, packed_t >;
            }// else
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic )
        {
            using  packed_t = packed< real_t, ISA_MIC >;
            
            HINFO( "(THelmholtzDLPGenFn) using MIC kernel" );
            
            if ( std::real( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< value_t, packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< value_t, packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< value_t, packed_t >;
            }// else
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            using  packed_t = packed< real_t, ISA_AVX2 >;
            
            HINFO( "(THelmholtzDLPGenFn) using AVX2 kernel" );
            
            if ( std::real( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< value_t, packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< value_t, packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< value_t, packed_t >;
            }// else
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            using  packed_t = packed< real_t, ISA_AVX >;
            
            HINFO( "(THelmholtzDLPGenFn) using AVX kernel" );
            
            if ( std::real( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< value_t, packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< value_t, packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< value_t, packed_t >;
            }// else
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            using  packed_t = packed< real_t, ISA_SSE3 >;
            
            HINFO( "(THelmholtzDLPGenFn) using SSE3 kernel" );
            
            if ( std::real( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< value_t, packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< value_t, packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< value_t, packed_t >;
            }// else
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            using  packed_t = packed< real_t, ISA_VSX >;
            
            HINFO( "(THelmholtzDLPGenFn) using VSX kernel" );
            
            if ( std::real( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< value_t, packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< value_t, packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< value_t, packed_t >;
            }// else
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            using  packed_t = packed< real_t, ISA_NEON >;
            
            HINFO( "(THelmholtzDLPGenFn) using NEON kernel" );
            
            if ( std::real( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< value_t, packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real_t(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< value_t, packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< value_t, packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< value_t, packed_t >;
            }// else
        }// if
        else
        {
            HINFO( "(THelmholtzDLPGenFn) using standard kernel" );
            
            _eval_dx_impl = helmholtz_slp_eval_dx_flt< value_t >;
            _eval_dy_impl = helmholtz_dlp_eval_dy_flt< value_t >;
        }// else
    }// if
    else
    {
        HINFO( "(THelmholtzDLPGenFn) using standard kernel" );
            
        _eval_dx_impl = helmholtz_slp_eval_dx_flt< value_t >;
        _eval_dy_impl = helmholtz_dlp_eval_dy_flt< value_t >;
    }// if
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

#define INST_ALL( type1, type2 )                                        \
    template class THelmholtzSLPBF< TConstFnSpace< type1 >,  TConstFnSpace< type1 >, type2 >; \
    template class THelmholtzSLPBF< TConstFnSpace< type1 >,  TLinearFnSpace< type1 >, type2 >; \
    template class THelmholtzSLPBF< TLinearFnSpace< type1 >, TConstFnSpace< type1 >, type2 >; \
    template class THelmholtzSLPBF< TLinearFnSpace< type1 >, TLinearFnSpace< type1 >, type2 >; \
                                                                        \
    template class THelmholtzDLPBF< TConstFnSpace< type1 >,  TConstFnSpace< type1 >, type2 >; \
    template class THelmholtzDLPBF< TConstFnSpace< type1 >,  TLinearFnSpace< type1 >, type2 >; \
    template class THelmholtzDLPBF< TLinearFnSpace< type1 >, TConstFnSpace< type1 >, type2 >; \
    template class THelmholtzDLPBF< TLinearFnSpace< type1 >, TLinearFnSpace< type1 >, type2 >; \
                                                                        \
    template class THelmholtzSLPGenFn< TConstFnSpace< type1 >,  TConstFnSpace< type1 >, type2 >; \
    template class THelmholtzSLPGenFn< TConstFnSpace< type1 >,  TLinearFnSpace< type1 >, type2 >; \
    template class THelmholtzSLPGenFn< TLinearFnSpace< type1 >, TConstFnSpace< type1 >, type2 >; \
    template class THelmholtzSLPGenFn< TLinearFnSpace< type1 >, TLinearFnSpace< type1 >, type2 >; \
                                                                        \
    template class THelmholtzDLPGenFn< TConstFnSpace< type1 >,  TConstFnSpace< type1 >, type2 >; \
    template class THelmholtzDLPGenFn< TConstFnSpace< type1 >,  TLinearFnSpace< type1 >, type2 >; \
    template class THelmholtzDLPGenFn< TLinearFnSpace< type1 >, TConstFnSpace< type1 >, type2 >; \
    template class THelmholtzDLPGenFn< TLinearFnSpace< type1 >, TLinearFnSpace< type1 >, type2 >;

INST_ALL( float, std::complex< float > )
INST_ALL( double, std::complex< double > )

//
// instantiate Helmholtz for TConstEdgeFnSpace for MaxwellEFIE
//
template
void
helmholtz_slp_flt< TConstEdgeFnSpace, TConstEdgeFnSpace, std::complex< double > > (
    const TGrid::triangle_t &                 tri0,
    const TGrid::triangle_t &                 tri1,
    const tripair_quad_rule_t *               rule,
    const std::complex< double >              ikappa,
    const TConstEdgeFnSpace *                 ansatz_sp,
    const TConstEdgeFnSpace *                 test_sp,
    std::vector< std::complex< double > > &   values );

}// namespace Hpro
