//
// Project     : HLIBpro
// File        : maxwell_simd.inc
// Description : Maxwell kernels using SIMD functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#if !defined(SIMD_ISA)
#  error "no SIMD ISA defined"
#endif

#include <vector>

#include "hpro/base/packed.hh"

#include "hpro/bem/TConstEdgeFnSpace.hh"
#include "hpro/bem/THelmholtzBF.hh"

namespace Hpro
{

namespace
{

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// local constants
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

const double  ONE_OVER_4PI = double(1) / (double(4) * Math::pi< double >());

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// local functions
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// form local quadrature vectors x and y and compute  y-x (Attention: y-x !!!)
//
template < typename T_packed >
inline
void
comp_diff_simd ( const double *  x1,
                 const double *  y1,
                 const double *  x2,
                 const double *  y2,
                 const T_packed  tri0coo[3][3],
                 const T_packed  tri1coo[3][3],
                 T_packed        diff[3] )
{
    using vpacked = T_packed;
    using value_t = typename vpacked::value_t;
    
    const vpacked  vONE( value_t(1) );

    // load quadrature points for first triangle
    const vpacked  a1 = load< vpacked >( x1 );
    const vpacked  a2 = load< vpacked >( y1 );
    const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );

    // load quadrature points for second triangle
    const vpacked  b1 = load< vpacked >( x2 );
    const vpacked  b2 = load< vpacked >( y2 );
    const vpacked  b0 = sub( sub( vONE, b1 ), b2 );

    //
    // compute y-x
    //

    vpacked  tmp1;
    
    diff[0] = sub( muladd( b2, tri1coo[0][2], muladd( b1, tri1coo[0][1], mul( b0, tri1coo[0][0] ) ) ),
                   muladd( a2, tri0coo[0][2], muladd( a1, tri0coo[0][1], mul( a0, tri0coo[0][0] ) ) ) );
                                                                                                     
    diff[1] = sub( muladd( b2, tri1coo[1][2], muladd( b1, tri1coo[1][1], mul( b0, tri1coo[1][0] ) ) ),
                   muladd( a2, tri0coo[1][2], muladd( a1, tri0coo[1][1], mul( a0, tri0coo[1][0] ) ) ) );
                                                                                                     
    diff[2] = sub( muladd( b2, tri1coo[2][2], muladd( b1, tri1coo[2][1], mul( b0, tri1coo[2][0] ) ) ),
                   muladd( a2, tri0coo[2][2], muladd( a1, tri0coo[2][1], mul( a0, tri0coo[2][0] ) ) ) );
}

//
// form local quadrature vectors x and y and compute  |x-y|²
//
template < typename T_packed >
inline
T_packed
comp_sqdist_simd ( const double *  x1,
                   const double *  y1,
                   const double *  x2,
                   const double *  y2,
                   const T_packed  tri0coo[3][3],
                   const T_packed  tri1coo[3][3] )
{
    using vpacked = T_packed;
    using value_t = typename vpacked::value_t;

    const vpacked  vONE( value_t(1) );

    // load quadrature points for first triangle
    const vpacked  a1 = load< vpacked >( x1 );
    const vpacked  a2 = load< vpacked >( y1 );
    const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );

    // load quadrature points for second triangle
    const vpacked  b1 = load< vpacked >( x2 );
    const vpacked  b2 = load< vpacked >( y2 );
    const vpacked  b0 = sub( sub( vONE, b1 ), b2 );

    //
    // compute |y-x|²
    //

    vpacked  dot, tmp1;
    
    // u0-v0
    tmp1 = sub( muladd( b2, tri1coo[0][2], muladd( b1, tri1coo[0][1], mul( b0, tri1coo[0][0] ) ) ),
                muladd( a2, tri0coo[0][2], muladd( a1, tri0coo[0][1], mul( a0, tri0coo[0][0] ) ) ) );
    dot  = mul( tmp1, tmp1 );                                                              
                                                                                                     
    // u1-v1                                                                               
    tmp1 = sub( muladd( b2, tri1coo[1][2], muladd( b1, tri1coo[1][1], mul( b0, tri1coo[1][0] ) ) ),
                muladd( a2, tri0coo[1][2], muladd( a1, tri0coo[1][1], mul( a0, tri0coo[1][0] ) ) ) );
    dot  = muladd( tmp1, tmp1, dot )   ;                                                   
                                                                                                     
    // u2-v2                                                                               
    tmp1 = sub( muladd( b2, tri1coo[2][2], muladd( b1, tri1coo[2][1], mul( b0, tri1coo[2][0] ) ) ),
                muladd( a2, tri0coo[2][2], muladd( a1, tri0coo[2][1], mul( a0, tri0coo[2][0] ) ) ) );
    dot  = muladd( tmp1, tmp1, dot );

    return dot;
}

}// namespace anonymous

///////////////////////////////////////////////////////////////
//
// some standard definitions used in all functions
//
///////////////////////////////////////////////////////////////

#define DEFINE_LOCAL_TYPE                               \
    using vpacked = T_packed;                           \
    using value_t = typename vpacked::value_t;          \
    const size_t   VECTOR_SIZE = vpacked::vector_size

#define LOAD_TRIANGLE_COORD                                             \
    const vpacked  t0[3][3] = { \
        { ansatz_sp->grid()->vertex( tri0.vtx[0] )[0],                  \
          ansatz_sp->grid()->vertex( tri0.vtx[1] )[0],                  \
          ansatz_sp->grid()->vertex( tri0.vtx[2] )[0] },                \
        { ansatz_sp->grid()->vertex( tri0.vtx[0] )[1],                  \
          ansatz_sp->grid()->vertex( tri0.vtx[1] )[1],                  \
          ansatz_sp->grid()->vertex( tri0.vtx[2] )[1] },                \
        { ansatz_sp->grid()->vertex( tri0.vtx[0] )[2],                  \
          ansatz_sp->grid()->vertex( tri0.vtx[1] )[2],                  \
          ansatz_sp->grid()->vertex( tri0.vtx[2] )[2] } };              \
                                                                        \
    const vpacked  t1[3][3] = { \
        { test_sp->grid()->vertex( tri1.vtx[0] )[0],                    \
          test_sp->grid()->vertex( tri1.vtx[1] )[0],                    \
          test_sp->grid()->vertex( tri1.vtx[2] )[0] },                  \
        { test_sp->grid()->vertex( tri1.vtx[0] )[1],                    \
          test_sp->grid()->vertex( tri1.vtx[1] )[1],                    \
          test_sp->grid()->vertex( tri1.vtx[2] )[1] },                  \
        { test_sp->grid()->vertex( tri1.vtx[0] )[2],                    \
          test_sp->grid()->vertex( tri1.vtx[1] )[2],                    \
          test_sp->grid()->vertex( tri1.vtx[2] )[2] } }

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// Helmholtz DLP kernel without normal direction
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
helmholtz_dlp_wo_normal_simd ( const TGrid::triangle_t &                tri0,
                               const TGrid::triangle_t &                tri1,
                               const tripair_quad_rule_t *              rule,
                               const std::complex< double >             ikappa,
                               const T_ansatzsp *                       ansatz_sp,
                               const T_testsp *                         test_sp,
                               std::vector< std::complex< double > > &  values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;

    const vpacked  vONE( value_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    //
    // eval VECTOR_SIZE  add the same time
    //

    const vpacked  vikappa_re( ikappa.real() );
    const vpacked  vikappa_im( ikappa.imag() );
    value_t        res_re[ VECTOR_SIZE ];
    value_t        res_im[ VECTOR_SIZE ];
    
    const size_t   n_pts = rule->npts;
    const auto     x1    = rule -> template x1< double >();
    const auto     y1    = rule -> template y1< double >();
    const auto     x2    = rule -> template x2< double >();
    const auto     y2    = rule -> template y2< double >();

    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute y-x
        //
        
        vpacked  diff[3];
        
        comp_diff_simd( x1 + i, y1 + i,
                        x2 + i, y2 + i,
                        t0, t1, diff );

        // |y-x|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );

        // |y-x|
        const vpacked  vdist  = sqrt( vdist2 );

        // |y-x|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1)
        // ────────────────────────────
        //           |x-y|³
        //

        //
        // compute exp( i·κ·|x-y| )
        //         = polar( exp( re ), im )
        //         = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        vpacked        vikd_re = mul( vikappa_re, vdist );
        
        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist );
        
        // exp( re(·) )
        // (also: multiply with 1/(4·π) / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( vONE_OVER_4PI, vdist3 ), exp( vikd_re ) ); 

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        vpacked  vc, vs;

        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()   (here we would multiply with 1/(4·π) / |x-y|³)
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        // i·κ·|x-y| - 1
        vikd_re = sub( vikd_re, vONE );
        
        // compute exp() · (i·κ·|x-y| - 1) and store final values
        store( mulsub( vc, vikd_re, mul( vs, vikd_im ) ), res_re );
        store( muladd( vc, vikd_im, mul( vs, vikd_re ) ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ i+j ] = std::complex< double >( res_re[j], res_im[j] );
    }// for
}

///////////////////////////////////////////////////////////////
//
// Helmholtz DLP kernel; double wavenumber
//

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
helmholtz_dlp_wo_normal_re_simd ( const TGrid::triangle_t &                tri0,
                                  const TGrid::triangle_t &                tri1,
                                  const tripair_quad_rule_t *              rule,
                                  const std::complex< double >             ikappa,
                                  const T_ansatzsp *                       ansatz_sp,
                                  const T_testsp *                         test_sp,
                                  std::vector< std::complex< double > > &  values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;

    const vpacked  vMINUS_ONE( value_t(-1) );  // used for: re(i·κ·|x-y|)-1 since re() = 0
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    //
    // eval VECTOR_SIZE  add the same time
    //

    const vpacked  vikappa_im( ikappa.imag() );
    value_t        res_re[ VECTOR_SIZE ];
    value_t        res_im[ VECTOR_SIZE ];
    
    const size_t   n_pts = rule->npts;
    const auto     x1    = rule -> template x1< double >();
    const auto     y1    = rule -> template y1< double >();
    const auto     x2    = rule -> template x2< double >();
    const auto     y2    = rule -> template y2< double >();

    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute y-x
        //
        
        vpacked  diff[3];
        
        comp_diff_simd( x1 + i, y1 + i,
                        x2 + i, y2 + i,
                        t0, t1, diff );

        // |y-x|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        
        // |y-x|
        const vpacked  vdist  = sqrt( vdist2 );

        // |y-x|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1)
        // ────────────────────────────
        //           |x-y|³
        //

        //
        // compute exp( i·κ·|x-y| )
        //         = polar( exp( re ), im )
        //         = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist );
        
        // exp( re(·) ) = 1 since re() = 0
        // (also: multiply with 1/(4·π) / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = div( vONE_OVER_4PI, vdist3 ); 

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        vpacked  vc, vs;

        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()   (here we would multiply with 1/(4·π) / |x-y|³)
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        // compute exp() · (i·κ·|x-y| - 1) and store final values
        store( mulsub( vc, vMINUS_ONE, mul( vs, vikd_im    ) ), res_re );
        store( muladd( vc, vikd_im,    mul( vs, vMINUS_ONE ) ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ i+j ] = std::complex< double >( res_re[j], res_im[j] );
    }// for
}

///////////////////////////////////////////////////////////////
//
// Helmholtz DLP kernel; imaginary wavenumber
//

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
helmholtz_dlp_wo_normal_im_simd ( const TGrid::triangle_t &                tri0,
                                  const TGrid::triangle_t &                tri1,
                                  const tripair_quad_rule_t *              rule,
                                  const std::complex< double >             ikappa,
                                  const T_ansatzsp *                       ansatz_sp,
                                  const T_testsp *                         test_sp,
                                  std::vector< std::complex< double > > &  values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;

    const vpacked  vONE( value_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    //
    // eval VECTOR_SIZE  add the same time
    //

    const vpacked  vikappa_re( ikappa.real() );
    value_t        res_re[ VECTOR_SIZE ];
    
    const size_t   n_pts = rule->npts;
    const auto     x1    = rule -> template x1< double >();
    const auto     y1    = rule -> template y1< double >();
    const auto     x2    = rule -> template x2< double >();
    const auto     y2    = rule -> template y2< double >();

    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute y-x
        //
        
        vpacked  diff[3];
        
        comp_diff_simd( x1 + i, y1 + i,
                        x2 + i, y2 + i,
                        t0, t1, diff );

        // |y-x|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        
        // |y-x|
        const vpacked  vdist  = sqrt( vdist2 );

        // |y-x|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1)
        // ────────────────────────────
        //           |x-y|³
        //

        //
        // compute exp( i·κ·|x-y| )
        //         = polar( exp( re ), im )
        //         = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        vpacked        vikd_re = mul( vikappa_re, vdist );
        
        // exp( re(·) )
        // (also: multiply with 1/(4·π) / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( vONE_OVER_4PI, vdist3 ), exp( vikd_re ) ); 

        // exp(re)·cos(), exp(re)·sin() = ( exp(re), 0 ) since im = 0

        // i·κ·|x-y| - 1
        vikd_re = sub( vikd_re, vONE );
        
        // compute exp() · (i·κ·|x-y| - 1) ( only double numbers ! )
        vikd_re = mul( vexp_ikd_re, vikd_re );
        
        // store final values
        store( vikd_re, res_re );
        
        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ i+j ] = std::complex< double >( res_re[j], double(0) );
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// explicit template instantiations
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#define  INST_HELMHOLTZ_DLP_WON_SIMD( T_ansatzsp, T_testsp, suffix )    \
    template void                                                       \
    helmholtz_dlp_wo_normal_##suffix< T_ansatzsp, T_testsp, packed< double, SIMD_ISA > > ( \
        const TGrid::triangle_t &                tri0,                  \
        const TGrid::triangle_t &                tri1,                  \
        const tripair_quad_rule_t *              rule,                  \
        const std::complex< double >             ikappa,                \
        const T_ansatzsp *                       ansatz_sp,             \
        const T_testsp *                         test_sp,               \
        std::vector< std::complex< double > > &  values )

INST_HELMHOLTZ_DLP_WON_SIMD( TConstEdgeFnSpace, TConstEdgeFnSpace, simd );
INST_HELMHOLTZ_DLP_WON_SIMD( TConstEdgeFnSpace, TConstEdgeFnSpace, re_simd );
INST_HELMHOLTZ_DLP_WON_SIMD( TConstEdgeFnSpace, TConstEdgeFnSpace, im_simd );

}// namespace Hpro

// Local Variables:
// mode: c++
// End:
