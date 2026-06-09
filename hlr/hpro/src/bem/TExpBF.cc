//
// Project     : HLIBpro
// File        : TExpBF.cc
// Description : bilinearforms for exponential kernel
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/base/packed.hh"

#include "hpro/bem/TExpBF.hh"

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
           typename packed_t >
void
exp_simd  ( const TGrid::triangle_t &    tri0,
            const TGrid::triangle_t &    tri1,
            const tripair_quad_rule_t *  rule,
            std::vector< value_t > &     values,
            const T_ansatzsp *           ansatz_sp,
            const T_testsp *             test_sp );

template < typename value_t,
           typename packed_t >
void
exp_eval_dx_simd  ( const tri_quad_rule_t< real_type_t< value_t > > & quad_rule,
                    const T3Point                                     x[3],
                    const T3Point &                                   y,
                    std::vector< value_t > &                          values );

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
exp_flt ( const TGrid::triangle_t &    tri0,
          const TGrid::triangle_t &    tri1,
          const tripair_quad_rule_t *  rule,
          std::vector< value_t > &     values,
          const T_ansatzsp *           ansatz_sp,
          const T_testsp *             test_sp )
{
    using real_t  = real_type_t< value_t >;

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
    T3Point        x, y, u;

    #pragma omp simd
    for ( size_t  i = 0; i < n_pts; ++i )
    {
        const auto  a1 = x1[i];
        const auto  a2 = y1[i];
        const auto  a0 = real_t(1) - a1 - a2;
        
        const auto  b1 = x2[i];
        const auto  b2 = y2[i];
        const auto  b0 = real_t(1) - b1 - b2;

        // u = x - y
        x = a0 * v00 + a1 * v01 + a2 * v02;
        y = b0 * v10 + b1 * v11 + b2 * v12;
        u = x - y;

        // exp( - |x-y| )
        values[ i ] = Math::exp( - norm2( u ) );

        // x_0 * exp( - |x-y| )
        // values[ i ] = x[0] * Math::exp( - norm2( u ) );
    }// for
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// TExpBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template < typename  ansatzsp_t,
           typename  testsp_t,
           typename  value_t >
auto
exp_kernel_fn ()
{
    using  real_t = real_type_t< value_t >;

    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f()  && CFG::BEM::use_simd_avx512f  )
        {
            HINFO( "(TExpBF) using AVX512F kernel" );
            
            return exp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_AVX512F > >;
        }// if
        else if ( CFG::Mach::has_mic()  && CFG::BEM::use_simd_mic  )
        {
            HINFO( "(TExpBF) using MIC kernel" );
            
            return exp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_MIC > >;
        }// if
        else if ( CFG::Mach::has_avx2()  && CFG::BEM::use_simd_avx2  )
        {
            HINFO( "(TExpBF) using AVX2 kernel" );
            
            return exp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_AVX2 > >;
        }// if
        else if ( CFG::Mach::has_avx()  && CFG::BEM::use_simd_avx  )
        {
            HINFO( "(TExpBF) using AVX kernel" );
            
            return exp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_AVX > >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            HINFO( "(TExpBF) using SSE3 kernel" );
            
            return exp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_SSE3 > >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            HINFO( "(TExpBF) using VSX kernel" );
            
            return exp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_VSX > >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            HINFO( "(TExpBF) using NEON kernel" );
            
            return exp_simd< value_t, ansatzsp_t, testsp_t, packed< real_t, ISA_NEON > >;
        }// if
        else
        {
            HINFO( "(TExpBF) using FPU kernel" );
            
            return exp_flt< ansatzsp_t, testsp_t, value_t >;
        }// else
    }// if
    else
    {
        HINFO( "(TExpBF) using FPU kernel" );
            
        return exp_flt< ansatzsp_t, testsp_t, value_t >;
    }// else
}

template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
TExpBF< T_ansatzsp, T_testsp, T_value >::TExpBF ( const ansatzsp_t *  aansatzsp,
                                                  const testsp_t *    atestsp,
                                                  const uint          quad_order )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, quad_order )
        , _kernel_fn( exp_kernel_fn< T_ansatzsp, T_testsp, T_value >() )
{}

template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
TExpBF< T_ansatzsp, T_testsp, T_value >::TExpBF ( const ansatzsp_t *  aansatzsp,
                                                  const testsp_t *    atestsp,
                                                  const real_t        quad_error )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp, 10, true, quad_error )
        , _kernel_fn( exp_kernel_fn< T_ansatzsp, T_testsp, T_value >() )
{}

//
// eval kernel function at quadrature points
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
void
TExpBF< T_ansatzsp, T_testsp, T_value >::eval_kernel ( const                        idx_t,
                                                       const                        idx_t,
                                                       const TGrid::triangle_t &    tri0,
                                                       const TGrid::triangle_t &    tri1,
                                                       const tripair_quad_rule_t *  rule,
                                                       std::vector< value_t > &     values,
                                                       const real_t                 eps ) const
{
    _kernel_fn( tri0, tri1, rule, values, this->ansatz_space(), this->test_space() );
}

//
// return format of bilinear form, e.g. symmetric
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
matform_t
TExpBF< T_ansatzsp, T_testsp, T_value >::format () const
{
    if ( ( cptrcast( this->ansatz_space(), TFnSpace< value_t > ) == cptrcast( this->test_space(), TFnSpace< value_t > ) ) ||
         (( this->ansatz_space()->type() == this->test_space()->type() ) &&
          ( this->ansatz_space()->grid() == this->test_space()->grid() )) )
        return symmetric;
    else
        return unsymmetric;
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// TExpGenFn
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// computes x-integral of SLP
//
template < typename value_t >
void
exp_eval_dx_flt ( const tri_quad_rule_t< real_type_t< value_t > > & quad_rule,
                  const T3Point                                     vx[3],
                  const T3Point &                                   vy,
                  std::vector< value_t > &                          values )
{
    using real_t = real_type_t< value_t >;

    const size_t  npts = quad_rule.npts;
    
    #pragma omp simd
    for ( size_t  k = 0; k < npts; ++k )
    {
        const auto  a2 = quad_rule.y[k];
        const auto  a1 = quad_rule.x[k];
        const auto  a0 = real_t(1) - a1 - a2;

        //
        // evaluate SLP kernel
        //

        // d = y - x
        const T3Point  xmy( vy - ( a0 * vx[0] + a1 * vx[1] + a2 * vx[2] ) );
        const auto     ddotd = value_t( Math::sqrt( dot( xmy, xmy ) ) );

        // e^( - | x - y | )
        values[k] = Math::exp( - ddotd );
    }// for
}

//
// constructor
//
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
TExpGenFn< T_ansatzsp, T_testsp, T_value >::TExpGenFn ( const ansatzsp_t *    ansatzsp,
                                                        const testsp_t *      testsp,
                                                        const TPermutation *  row_perm_i2e,
                                                        const TPermutation *  col_perm_i2e,
                                                        const uint            quad_order )
        : TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, value_t >( ansatzsp, testsp, quad_order,
                                                                    row_perm_i2e, col_perm_i2e )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            HINFO( "(TExpGenFn) using AVX512F kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< value_t, packed< real_t, ISA_AVX512F > >;
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic )
        {
            HINFO( "(TExpGenFn) using MIC kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< value_t, packed< real_t, ISA_MIC > >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            HINFO( "(TExpGenFn) using AVX2 kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< value_t, packed< real_t, ISA_AVX2 > >;
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            HINFO( "(TExpGenFn) using AVX kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< value_t, packed< real_t, ISA_AVX > >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            HINFO( "(TExpGenFn) using SSE3 kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< value_t, packed< real_t, ISA_SSE3 > >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            HINFO( "(TExpGenFn) using VSX kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< value_t, packed< real_t, ISA_VSX > >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            HINFO( "(TExpGenFn) using NEON kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< value_t, packed< real_t, ISA_NEON > >;
        }// if
        else
        {
            HINFO( "(TExpGenFn) using standard kernel" );
            
            _eval_dxy_impl = exp_eval_dx_flt;
        }// else
    }// if
    else
    {
        HINFO( "(TExpGenFn) using standard kernel" );
            
        _eval_dxy_impl = exp_eval_dx_flt;
    }// else
}
        
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

#define INST_ALL( type )                                                \
    template class TExpBF< TConstFnSpace< type >,  TConstFnSpace< type >, type >; \
    template class TExpBF< TConstFnSpace< type >,  TLinearFnSpace< type >, type >; \
    template class TExpBF< TLinearFnSpace< type >, TConstFnSpace< type >, type >; \
    template class TExpBF< TLinearFnSpace< type >, TLinearFnSpace< type >, type >; \
                                                                        \
    template class TExpGenFn< TConstFnSpace< type >,  TConstFnSpace< type >, type >; \
    template class TExpGenFn< TConstFnSpace< type >,  TLinearFnSpace< type >, type >; \
    template class TExpGenFn< TLinearFnSpace< type >, TConstFnSpace< type >, type >; \
    template class TExpGenFn< TLinearFnSpace< type >, TLinearFnSpace< type >, type >;

INST_ALL( float )
INST_ALL( double )

}// namespace Hpro
