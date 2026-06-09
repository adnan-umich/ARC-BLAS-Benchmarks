//
// Project     : HLIBpro
// File        : helmholtz_simd.inc
// Description : Helmholtz kernels using SIMD functions
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
// local functions
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// form local quadrature vectors x and y and compute  y-x (Attention: y-x !!!)
//
template < typename packed_t >
inline
void
comp_diff_simd ( const value_type_t< packed_t > * x1,
                 const value_type_t< packed_t > * y1,
                 const value_type_t< packed_t > * x2,
                 const value_type_t< packed_t > * y2,
                 const packed_t                   tri0coo[3][3],
                 const packed_t                   tri1coo[3][3],
                 packed_t                         diff[3],
                 const bool                       adjoint = false )
{
    using vpacked = packed_t;
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

    if ( adjoint )
    {
        //
        // compute x-y
        //

        vpacked  tmp1;
    
        diff[0] = sub( muladd( a2, tri0coo[0][2], muladd( a1, tri0coo[0][1], mul( a0, tri0coo[0][0] ) ) ),
                       muladd( b2, tri1coo[0][2], muladd( b1, tri1coo[0][1], mul( b0, tri1coo[0][0] ) ) ) );
        
        diff[1] = sub( muladd( a2, tri0coo[1][2], muladd( a1, tri0coo[1][1], mul( a0, tri0coo[1][0] ) ) ),
                       muladd( b2, tri1coo[1][2], muladd( b1, tri1coo[1][1], mul( b0, tri1coo[1][0] ) ) ) );
        
        diff[2] = sub( muladd( a2, tri0coo[2][2], muladd( a1, tri0coo[2][1], mul( a0, tri0coo[2][0] ) ) ),
                       muladd( b2, tri1coo[2][2], muladd( b1, tri1coo[2][1], mul( b0, tri1coo[2][0] ) ) ) );
    }// if
    else
    {
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
    }// else
}

//
// form local quadrature vectors x and y and compute  |x-y|²
//
template < typename packed_t >
inline
packed_t
comp_sqdist_simd ( const value_type_t< packed_t > *  x1,
                   const value_type_t< packed_t > *  y1,
                   const value_type_t< packed_t > *  x2,
                   const value_type_t< packed_t > *  y2,
                   const packed_t                    tri0coo[3][3],
                   const packed_t                    tri1coo[3][3] )
{
    using vpacked = packed_t;
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

//
// get normal direction for all quadrature points in vector data
//
template < typename  T_ansatzsp,
           typename  T_packed >
inline
void
get_normal ( const T_packed             va1,
             const T_packed             va2,
             const idx_t                tri_id,
             const TGrid::triangle_t &  tri,
             const T_ansatzsp *         func_sp,
             T_packed                   vn[3] )
{
    using vpacked = T_packed;
    using value_t = typename vpacked::value_t;

    const size_t   VECTOR_SIZE = vpacked::vector_size;
    
    //
    // convert vector data to array
    //
    
    value_t  a1[ VECTOR_SIZE ], a2[ VECTOR_SIZE ];
    value_t  x[ VECTOR_SIZE ],  y[ VECTOR_SIZE ],  z[ VECTOR_SIZE ];

    store( va1, a1 );
    store( va2, a2 );
    
    for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
    {
        const auto  n = func_sp->grid()->tri_normal( tri_id, tri, a1[j], a2[j] );

        x[j] = n.x();
        y[j] = n.y();
        z[j] = n.z();
    }// for

    vn[0] = load< vpacked >( x );
    vn[1] = load< vpacked >( y );
    vn[2] = load< vpacked >( z );
}

}// namespace anonymous

///////////////////////////////////////////////////////////////
//
// some standard definitions used in all functions
//
///////////////////////////////////////////////////////////////

#define DEFINE_LOCAL_TYPE                               \
    using vpacked = packed_t;                           \
    using real_t  = real_type_t< value_t >;             \
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
// Helmholtz SLP kernel
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename value_t,
           typename T_ansatzsp,
           typename T_testsp,
           typename packed_t >
void
helmholtz_slp_simd ( const TGrid::triangle_t &    tri0,
                     const TGrid::triangle_t &    tri1,
                     const tripair_quad_rule_t *  rule,
                     const value_t                ikappa,
                     const T_ansatzsp *           ansatz_sp,
                     const T_testsp *             test_sp,
                     std::vector< value_t > &     values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;
    
    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    const vpacked  vikappa_re( ikappa.real() );
    const vpacked  vikappa_im( ikappa.imag() );

    const size_t   n_pts = rule->npts;
    const auto     x1    = rule -> template x1< real_t >();
    const auto     y1    = rule -> template y1< real_t >();
    const auto     x2    = rule -> template x2< real_t >();
    const auto     y2    = rule -> template y2< real_t >();

    //
    // eval VECTOR_SIZE add the same time
    //

    real_t  res_re[ VECTOR_SIZE ];
    real_t  res_im[ VECTOR_SIZE ];

    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute and store |x-y|
        //
        
        const vpacked  vdist = sqrt( comp_sqdist_simd< vpacked >( x1 + i, y1 + i,
                                                                  x2 + i, y2 + i,
                                                                  t0, t1 ) );

        //
        //   ( i·κ·|x-y| )
        // e
        // ───────────────
        //      |x-y|
        //

        //
        // compute exp( i·κ·|x-y| )
        //         = polar( exp( re ), im )
        //         = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //
        
        // exp( re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) ) )
        // (also: multiply with 1/(4·π) / |x-y| as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( vONE_OVER_4PI, vdist ), exp( mul( vikappa_re, vdist ) ) ); 

        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist ); 

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        vpacked        vc, vs;

        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        // combine double and imaginary parts
        store( vc, res_re );
        store( vs, res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; j++ )
            values[ i+j ] = value_t( res_re[j], res_im[j] );
    }// for
}

///////////////////////////////////////////////////////////////
//
// Helmholtz SLP kernel; double wave number
//

template < typename value_t,
           typename T_ansatzsp,
           typename T_testsp,
           typename packed_t >
void
helmholtz_slp_re_simd ( const TGrid::triangle_t &    tri0,
                        const TGrid::triangle_t &    tri1,
                        const tripair_quad_rule_t *  rule,
                        const value_t                ikappa,
                        const T_ansatzsp *           ansatz_sp,
                        const T_testsp *             test_sp,
                        std::vector< value_t > &     values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;
    
    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    const vpacked  vikappa_im( ikappa.imag() );

    const size_t   n_pts = rule->npts;
    const auto     x1    = rule -> template x1< real_t >();
    const auto     y1    = rule -> template y1< real_t >();
    const auto     x2    = rule -> template x2< real_t >();
    const auto     y2    = rule -> template y2< real_t >();

    //
    // eval VECTOR_SIZE add the same time
    //

    real_t  res_re[ VECTOR_SIZE ];
    real_t  res_im[ VECTOR_SIZE ];

    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute and store |x-y|
        //
        
        const vpacked  vdist = sqrt( comp_sqdist_simd< vpacked >( x1 + i, y1 + i,
                                                                  x2 + i, y2 + i,
                                                                  t0, t1 ) );

        //
        //   ( i·κ·|x-y| )
        // e
        // ───────────────
        //      |x-y|
        //

        //
        // compute exp( i·κ·|x-y| )
        //         = polar( exp( re ), im )
        //         = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //
        
        // exp( re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) ) ) = 1 since re = 0
        // (also: multiply with 1/(4·π) / |x-y| as this would be done later anyway)
        const vpacked  vexp_ikd_re = div( vONE_OVER_4PI, vdist ); 

        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im     = mul( vikappa_im, vdist ); 
        vpacked        vc, vs;

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        // combine double and imaginary parts
        store( vc, res_re );
        store( vs, res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; j++ )
            values[ i+j ] = value_t( res_re[j], res_im[j] );
    }// for
}

///////////////////////////////////////////////////////////////
//
// Helmholtz SLP kernel; imaginary wave number
//

template < typename value_t,
           typename T_ansatzsp,
           typename T_testsp,
           typename packed_t >
void
helmholtz_slp_im_simd ( const TGrid::triangle_t &    tri0,
                        const TGrid::triangle_t &    tri1,
                        const tripair_quad_rule_t *  rule,
                        const value_t                ikappa,
                        const T_ansatzsp *           ansatz_sp,
                        const T_testsp *             test_sp,
                        std::vector< value_t > &     values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;
    
    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const vpacked   vONE_OVER_4PI( ONE_OVER_4PI );
    const vpacked   vikappa_re( ikappa.real() );

    const size_t   n_pts = rule->npts;
    const auto     x1    = rule -> template x1< real_t >();
    const auto     y1    = rule -> template y1< real_t >();
    const auto     x2    = rule -> template x2< real_t >();
    const auto     y2    = rule -> template y2< real_t >();
    
    //
    // eval VECTOR_SIZE add the same time
    //

    real_t  res_re[ VECTOR_SIZE ];
    
    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute and store |x-y|
        //
        
        const vpacked  vdist = sqrt( comp_sqdist_simd< vpacked >( x1 + i, y1 + i,
                                                                  x2 + i, y2 + i,
                                                                  t0, t1 ) );

        //
        //   ( i·κ·|x-y| )
        // e
        // ───────────────
        //      |x-y|
        //

        //
        // compute exp( i·κ·|x-y| )
        //         = polar( exp( re ), im )
        //         = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //
        
        // exp( re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) ) )
        // (also: multiply with 1/(4·π) / |x-y| as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( vONE_OVER_4PI, vdist ), exp( mul( vikappa_re, vdist ) ) ); 

        // exp(re)·cos(im), exp(re)·sin(im) = ( exp(re), 0 ) since im = 0
        // combine double and imaginary parts
        store( vexp_ikd_re, res_re );

        for ( size_t  j = 0; j < VECTOR_SIZE; j++ )
            values[ i+j ] = value_t( res_re[j], real_t(0) );
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// Helmholtz DLP kernel
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename value_t,
           typename T_ansatzsp,
           typename T_testsp,
           typename packed_t >
void
helmholtz_dlp_simd ( const idx_t                  tri_id,
                     const bool                   adjoint,
                     const TGrid::triangle_t &    tri0,
                     const TGrid::triangle_t &    tri1,
                     const tripair_quad_rule_t *  rule,
                     const value_t                ikappa,
                     const T_ansatzsp *           ansatz_sp,
                     const T_testsp *             test_sp,
                     std::vector< value_t > &     values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;

    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const vpacked  vONE( real_t(1) );
    const vpacked  vMINUSONE( real_t(-1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    //
    // eval VECTOR_SIZE  add the same time
    //

    const vpacked  vikappa_re( ikappa.real() );
    const vpacked  vikappa_im( ikappa.imag() );
    real_t         res_re[ VECTOR_SIZE ];
    real_t         res_im[ VECTOR_SIZE ];
    
    const size_t   n_pts = rule->npts;
    const auto     x1    = rule -> template x1< real_t >();
    const auto     y1    = rule -> template y1< real_t >();
    const auto     x2    = rule -> template x2< real_t >();
    const auto     y2    = rule -> template y2< real_t >();

    // set normal direction
    const bool     has_vtx_normal = ( adjoint ? ansatz_sp->grid()->has_vtx_normal() : test_sp->grid()->has_vtx_normal() );
    const auto     n_vec          = ( adjoint ? ansatz_sp->grid()->tri_normal( tri_id ) : test_sp->grid()->tri_normal( tri_id ) );
    vpacked        n[3]           = { n_vec[0], n_vec[1], n_vec[2] };

    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute y-x
        //
        
        vpacked  diff[3];
        
        // load quadrature points for first triangle
        const vpacked  a1 = load< vpacked >( x1 + i );
        const vpacked  a2 = load< vpacked >( y1 + i );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );

        // load quadrature points for second triangle
        const vpacked  b1 = load< vpacked >( x2 + i );
        const vpacked  b2 = load< vpacked >( y2 + i );
        const vpacked  b0 = sub( sub( vONE, b1 ), b2 );

        if ( adjoint )
        {
            //
            // compute x-y
            //

            vpacked  tmp1;
    
            diff[0] = sub( muladd( a2, t0[0][2], muladd( a1, t0[0][1], mul( a0, t0[0][0] ) ) ),
                           muladd( b2, t1[0][2], muladd( b1, t1[0][1], mul( b0, t1[0][0] ) ) ) );
        
            diff[1] = sub( muladd( a2, t0[1][2], muladd( a1, t0[1][1], mul( a0, t0[1][0] ) ) ),
                           muladd( b2, t1[1][2], muladd( b1, t1[1][1], mul( b0, t1[1][0] ) ) ) );
        
            diff[2] = sub( muladd( a2, t0[2][2], muladd( a1, t0[2][1], mul( a0, t0[2][0] ) ) ),
                           muladd( b2, t1[2][2], muladd( b1, t1[2][1], mul( b0, t1[2][0] ) ) ) );
        }// if
        else
        {
            //
            // compute y-x
            //

            vpacked  tmp1;
    
            diff[0] = sub( muladd( b2, t1[0][2], muladd( b1, t1[0][1], mul( b0, t1[0][0] ) ) ),
                           muladd( a2, t0[0][2], muladd( a1, t0[0][1], mul( a0, t0[0][0] ) ) ) );
        
            diff[1] = sub( muladd( b2, t1[1][2], muladd( b1, t1[1][1], mul( b0, t1[1][0] ) ) ),
                           muladd( a2, t0[1][2], muladd( a1, t0[1][1], mul( a0, t0[1][0] ) ) ) );
        
            diff[2] = sub( muladd( b2, t1[2][2], muladd( b1, t1[2][1], mul( b0, t1[2][0] ) ) ),
                           muladd( a2, t0[2][2], muladd( a1, t0[2][1], mul( a0, t0[2][0] ) ) ) );
        }// else
        // comp_diff_simd( & rule->x1[i], & rule->y1[i],
        //                 & rule->x2[i], & rule->y2[i],
        //                 t0, t1, diff, adjoint );

        // |y-x|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |y-x|
        const vpacked  vdist  = sqrt( vdist2 );

        // |y-x|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        // get normal direction
        if ( has_vtx_normal )
        {
            if ( adjoint )
                get_normal( a1, a2, tri_id, tri0, ansatz_sp, n );
            else
                get_normal( b1, b2, tri_id, tri1, test_sp, n );
        }// if
    
        // < n, y-x >
        const vpacked  vndu = ( adjoint
                                ? muladd( diff[0], n[0], muladd( diff[1], n[1], mul( diff[2], n[2] ) ) )
                                : muladd( n[0], diff[0], muladd( n[1], diff[1], mul( n[2], diff[2] ) ) ) );

        
        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n(y),y-x>
        // ───────────────────────────────────────
        //                |x-y|³
        //
        // or
        //
        //   (i·κ·|x-y|)
        // e            (1 - i·κ·|x-y|) <n(x),x-y>
        // ───────────────────────────────────────
        //                |x-y|³
        //
        
        //
        // compute exp( i·κ·|x-y| )
        //         = polar( exp( re ), im )
        //         = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        vpacked        vikd_re = mul( vikappa_re, vdist );
        
        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        vpacked        vikd_im = mul( vikappa_im, vdist );
        
        // exp( re(·) )
        // (also: multiply with 1/(4·π) <n,y-x> / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( mul( vONE_OVER_4PI, vndu ), vdist3 ), exp( vikd_re ) ); 

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        vpacked  vc, vs;

        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()   (here we would multiply with 1/(4·π) <n,y-x> / |x-y|³)
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        if ( adjoint )
        {
            // 1 - i·κ·|x-y|
            vikd_re = sub( vONE, vikd_re );
            vikd_im = mul( vMINUSONE, vikd_im );
        }// if
        else
        {
            // i·κ·|x-y| - 1
            vikd_re = sub( vikd_re, vONE );
        }// else
        
        // compute exp() · (i·κ·|x-y| - 1) and store final values
        store( mulsub( vc, vikd_re, mul( vs, vikd_im ) ), res_re );
        store( muladd( vc, vikd_im, mul( vs, vikd_re ) ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ i+j ] = value_t( res_re[j], res_im[j] );
    }// for
}

///////////////////////////////////////////////////////////////
//
// Helmholtz DLP kernel; double wavenumber
//

template < typename value_t,
           typename T_ansatzsp,
           typename T_testsp,
           typename packed_t >
void
helmholtz_dlp_re_simd ( const idx_t                  tri_id,
                        const bool                   adjoint,
                        const TGrid::triangle_t &    tri0,
                        const TGrid::triangle_t &    tri1,
                        const tripair_quad_rule_t *  rule,
                        const value_t                ikappa,
                        const T_ansatzsp *           ansatz_sp,
                        const T_testsp *             test_sp,
                        std::vector< value_t > &     values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;

    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());

    const vpacked  vONE( real_t(1) );         // used for: re(i·κ·|x-y|)-1 since re() = 0
    const vpacked  vMINUS_ONE( real_t(-1) );  // used for: re(i·κ·|x-y|)-1 since re() = 0
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    //
    // eval VECTOR_SIZE  add the same time
    //

    const vpacked  vikappa_im( ikappa.imag() );
    real_t         res_re[ VECTOR_SIZE ];
    real_t         res_im[ VECTOR_SIZE ];
    
    const size_t   n_pts = rule->npts;
    const auto     x1    = rule -> template x1< real_t >();
    const auto     y1    = rule -> template y1< real_t >();
    const auto     x2    = rule -> template x2< real_t >();
    const auto     y2    = rule -> template y2< real_t >();

    // set normal direction
    const bool     has_vtx_normal = ( adjoint ? ansatz_sp->grid()->has_vtx_normal() : test_sp->grid()->has_vtx_normal() );
    const auto     n_vec          = ( adjoint ? ansatz_sp->grid()->tri_normal( tri_id ) : test_sp->grid()->tri_normal( tri_id ) );
    vpacked        n[3]           = { n_vec[0], n_vec[1], n_vec[2] };

    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute y-x
        //
        
        vpacked  diff[3];
        
        // load quadrature points for first triangle
        const vpacked  a1 = load< vpacked >( x1 + i );
        const vpacked  a2 = load< vpacked >( y1 + i );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );

        // load quadrature points for second triangle
        const vpacked  b1 = load< vpacked >( x2 + i );
        const vpacked  b2 = load< vpacked >( y2 + i );
        const vpacked  b0 = sub( sub( vONE, b1 ), b2 );

        if ( adjoint )
        {
            //
            // compute x-y
            //

            vpacked  tmp1;
    
            diff[0] = sub( muladd( a2, t0[0][2], muladd( a1, t0[0][1], mul( a0, t0[0][0] ) ) ),
                           muladd( b2, t1[0][2], muladd( b1, t1[0][1], mul( b0, t1[0][0] ) ) ) );
        
            diff[1] = sub( muladd( a2, t0[1][2], muladd( a1, t0[1][1], mul( a0, t0[1][0] ) ) ),
                           muladd( b2, t1[1][2], muladd( b1, t1[1][1], mul( b0, t1[1][0] ) ) ) );
        
            diff[2] = sub( muladd( a2, t0[2][2], muladd( a1, t0[2][1], mul( a0, t0[2][0] ) ) ),
                           muladd( b2, t1[2][2], muladd( b1, t1[2][1], mul( b0, t1[2][0] ) ) ) );
        }// if
        else
        {
            //
            // compute y-x
            //

            vpacked  tmp1;
    
            diff[0] = sub( muladd( b2, t1[0][2], muladd( b1, t1[0][1], mul( b0, t1[0][0] ) ) ),
                           muladd( a2, t0[0][2], muladd( a1, t0[0][1], mul( a0, t0[0][0] ) ) ) );
        
            diff[1] = sub( muladd( b2, t1[1][2], muladd( b1, t1[1][1], mul( b0, t1[1][0] ) ) ),
                           muladd( a2, t0[1][2], muladd( a1, t0[1][1], mul( a0, t0[1][0] ) ) ) );
        
            diff[2] = sub( muladd( b2, t1[2][2], muladd( b1, t1[2][1], mul( b0, t1[2][0] ) ) ),
                           muladd( a2, t0[2][2], muladd( a1, t0[2][1], mul( a0, t0[2][0] ) ) ) );
        }// else
        // comp_diff_simd( & rule->x1[i], & rule->y1[i],
        //                 & rule->x2[i], & rule->y2[i],
        //                 t0, t1, diff, adjoint );

        // |y-x|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |y-x|
        const vpacked  vdist  = sqrt( vdist2 );

        // |y-x|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        // get normal direction
        if ( has_vtx_normal )
        {
            if ( adjoint )
                get_normal( a1, a2, tri_id, tri0, ansatz_sp, n );
            else
                get_normal( b1, b2, tri_id, tri1, test_sp, n );
        }// if

        // < n, y-x >
        const vpacked  vndu = ( adjoint
                                ? muladd( diff[0], n[0], muladd( diff[1], n[1], mul( diff[2], n[2] ) ) )
                                : muladd( n[0], diff[0], muladd( n[1], diff[1], mul( n[2], diff[2] ) ) ) );
        
        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n,y-x>
        // ────────────────────────────────────
        //                |x-y|³
        //

        //
        // compute exp( i·κ·|x-y| )
        //         = polar( exp( re ), im )
        //         = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist );
        
        // exp( re(·) ) = 1 since re() = 0
        // (also: multiply with 1/(4·π) <n,y-x> / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = div( mul( vONE_OVER_4PI, vndu ), vdist3 ); 

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        vpacked  vc, vs;

        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()   (here we would multiply with 1/(4·π) <n,y-x> / |x-y|³)
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        if ( adjoint )
        {
            // compute exp() · (1 - i·κ·|x-y|) and store final values
            store( muladd(    vs, vikd_im, vc ), res_re );
            store( negmuladd( vc, vikd_im, vs ), res_im );
        }// if
        else
        {
            // compute exp() · (i·κ·|x-y| - 1) and store final values
            store( mulsub( vc, vMINUS_ONE, mul( vs, vikd_im    ) ), res_re );
            store( muladd( vc, vikd_im,    mul( vs, vMINUS_ONE ) ), res_im );
        }// else

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ i+j ] = value_t( res_re[j], res_im[j] );
    }// for
}

///////////////////////////////////////////////////////////////
//
// Helmholtz DLP kernel; imaginary wavenumber
//

template < typename value_t,
           typename T_ansatzsp,
           typename T_testsp,
           typename packed_t >
void
helmholtz_dlp_im_simd ( const idx_t                  tri_id,
                        const bool                   adjoint,
                        const TGrid::triangle_t &    tri0,
                        const TGrid::triangle_t &    tri1,
                        const tripair_quad_rule_t *  rule,
                        const value_t                ikappa,
                        const T_ansatzsp *           ansatz_sp,
                        const T_testsp *             test_sp,
                        std::vector< value_t > &     values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;

    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const vpacked  vONE( real_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    //
    // eval VECTOR_SIZE  add the same time
    //

    const vpacked  vikappa_re( ikappa.real() );
    real_t         res_re[ VECTOR_SIZE ];
    
    const size_t   n_pts = rule->npts;
    const auto     x1    = rule -> template x1< real_t >();
    const auto     y1    = rule -> template y1< real_t >();
    const auto     x2    = rule -> template x2< real_t >();
    const auto     y2    = rule -> template y2< real_t >();

    // set normal direction
    const bool     has_vtx_normal = ( adjoint ? ansatz_sp->grid()->has_vtx_normal() : test_sp->grid()->has_vtx_normal() );
    const auto     n_vec          = ( adjoint ? ansatz_sp->grid()->tri_normal( tri_id ) : test_sp->grid()->tri_normal( tri_id ) );
    vpacked        n[3]           = { n_vec[0], n_vec[1], n_vec[2] };

    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute y-x
        //
        
        vpacked  diff[3];
        
        // load quadrature points for first triangle
        const vpacked  a1 = load< vpacked >( x1 + i );
        const vpacked  a2 = load< vpacked >( y1 + i );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );

        // load quadrature points for second triangle
        const vpacked  b1 = load< vpacked >( x2 + i );
        const vpacked  b2 = load< vpacked >( y2 + i );
        const vpacked  b0 = sub( sub( vONE, b1 ), b2 );

        if ( adjoint )
        {
            //
            // compute x-y
            //

            vpacked  tmp1;
    
            diff[0] = sub( muladd( a2, t0[0][2], muladd( a1, t0[0][1], mul( a0, t0[0][0] ) ) ),
                           muladd( b2, t1[0][2], muladd( b1, t1[0][1], mul( b0, t1[0][0] ) ) ) );
        
            diff[1] = sub( muladd( a2, t0[1][2], muladd( a1, t0[1][1], mul( a0, t0[1][0] ) ) ),
                           muladd( b2, t1[1][2], muladd( b1, t1[1][1], mul( b0, t1[1][0] ) ) ) );
        
            diff[2] = sub( muladd( a2, t0[2][2], muladd( a1, t0[2][1], mul( a0, t0[2][0] ) ) ),
                           muladd( b2, t1[2][2], muladd( b1, t1[2][1], mul( b0, t1[2][0] ) ) ) );
        }// if
        else
        {
            //
            // compute y-x
            //

            vpacked  tmp1;
    
            diff[0] = sub( muladd( b2, t1[0][2], muladd( b1, t1[0][1], mul( b0, t1[0][0] ) ) ),
                           muladd( a2, t0[0][2], muladd( a1, t0[0][1], mul( a0, t0[0][0] ) ) ) );
        
            diff[1] = sub( muladd( b2, t1[1][2], muladd( b1, t1[1][1], mul( b0, t1[1][0] ) ) ),
                           muladd( a2, t0[1][2], muladd( a1, t0[1][1], mul( a0, t0[1][0] ) ) ) );
        
            diff[2] = sub( muladd( b2, t1[2][2], muladd( b1, t1[2][1], mul( b0, t1[2][0] ) ) ),
                           muladd( a2, t0[2][2], muladd( a1, t0[2][1], mul( a0, t0[2][0] ) ) ) );
        }// else
        // comp_diff_simd( & rule->x1[i], & rule->y1[i],
        //                 & rule->x2[i], & rule->y2[i],
        //                 t0, t1, diff, adjoint );

        // |y-x|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |y-x|
        const vpacked  vdist  = sqrt( vdist2 );

        // |y-x|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        // get normal direction
        if ( has_vtx_normal )
        {
            if ( adjoint )
                get_normal( a1, a2, tri_id, tri0, ansatz_sp, n );
            else
                get_normal( b1, b2, tri_id, tri1, test_sp, n );
        }// if

        // < n, y-x >
        const vpacked  vndu = ( adjoint
                                ? muladd( diff[0], n[0], muladd( diff[1], n[1], mul( diff[2], n[2] ) ) )
                                : muladd( n[0], diff[0], muladd( n[1], diff[1], mul( n[2], diff[2] ) ) ) );

        
        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n,y-x>
        // ────────────────────────────────────
        //                |x-y|³
        //

        //
        // compute exp( i·κ·|x-y| )
        //         = polar( exp( re ), im )
        //         = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        vpacked        vikd_re = mul( vikappa_re, vdist );
        
        // exp( re(·) )
        // (also: multiply with 1/(4·π) <n,y-x> / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( mul( vONE_OVER_4PI, vndu ), vdist3 ), exp( vikd_re ) ); 

        // exp(re)·cos(), exp(re)·sin() = ( exp(re), 0 ) since im = 0

        if ( adjoint )
        {
            // 1 - i·κ·|x-y|
            vikd_re = sub( vONE, vikd_re );
        }// if
        else
        {
            // i·κ·|x-y| - 1
            vikd_re = sub( vikd_re, vONE );
        }// else
        
        // compute exp() · (i·κ·|x-y| - 1) ( only double numbers ! )
        vikd_re = mul( vexp_ikd_re, vikd_re );
        
        // store final values
        store( vikd_re, res_re );
        
        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ i+j ] = value_t( res_re[j], real_t(0) );
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// Helmholtz HCA functions
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename value_t,
           typename packed_t >
void
helmholtz_slp_eval_dx_simd ( const value_t &           ikappa,
                             const tri_quad_rule_t< real_type_t< value_t > > &  quad_rule,
                             const T3Point             vx[3],
                             const T3Point &           vy,
                             std::vector< value_t > &  values )
{
    DEFINE_LOCAL_TYPE;

    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const vpacked  vONE( real_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    const vpacked  vikappa_re( ikappa.real() );
    const vpacked  vikappa_im( ikappa.imag() );
    const size_t   npts  = quad_rule.npts;
    const vpacked  y[3]  = { vy.x(),    vy.y(),    vy.z() };
    const vpacked  x0[3] = { vx[0].x(), vx[0].y(), vx[0].z() };
    const vpacked  x1[3] = { vx[1].x(), vx[1].y(), vx[1].z() };
    const vpacked  x2[3] = { vx[2].x(), vx[2].y(), vx[2].z() };
    real_t         res_re[ VECTOR_SIZE ];
    real_t         res_im[ VECTOR_SIZE ];

    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        // load quadrature points and weight for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // compute d = y - x; <d,d>
        //

        vpacked  diff[3];

        diff[0] = sub( y[0], muladd( a0, x0[0], muladd( a1, x1[0], mul( a2, x2[0] ) ) ) );
        diff[1] = sub( y[1], muladd( a0, x0[1], muladd( a1, x1[1], mul( a2, x2[1] ) ) ) );
        diff[2] = sub( y[2], muladd( a0, x0[2], muladd( a1, x1[2], mul( a2, x2[2] ) ) ) );

        // |x-y|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |x-y|
        const vpacked  vdist  = sqrt( vdist2 );

        //
        //   ( i·κ·|x-y| )
        // e
        // ───────────────
        //      |x-y|
        //
        // compute exp( i·κ·|x-y| ) = polar( exp( re ), im )
        //                          = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // exp( re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) ) )
        // (also: multiply with (1/4·π) / |x-y| as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( vONE_OVER_4PI, vdist ), exp( mul( vikappa_re, vdist ) ) );

        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist ); 
        vpacked        vc, vs;

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()
        // combine double and imaginary parts
        store( mul( vexp_ikd_re, vc ), res_re );
        store( mul( vexp_ikd_re, vs ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ k+j ] = value_t( res_re[j], res_im[j] );
    }// for
}

template < typename value_t,
           typename packed_t >
void
helmholtz_slp_eval_dx_re_simd ( const value_t &           ikappa,
                                const tri_quad_rule_t< real_type_t< value_t > > &            quad_rule,
                                const T3Point             vx[3],
                                const T3Point &           vy,
                                std::vector< value_t > &  values )
{
    DEFINE_LOCAL_TYPE;

    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const vpacked  vONE( real_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    const vpacked  vikappa_im( ikappa.imag() );  // κ double -> i·κ imaginary
    const size_t   npts  = quad_rule.npts;
    const vpacked  y[3]  = { vy.x(),    vy.y(),    vy.z() };
    const vpacked  x0[3] = { vx[0].x(), vx[0].y(), vx[0].z() };
    const vpacked  x1[3] = { vx[1].x(), vx[1].y(), vx[1].z() };
    const vpacked  x2[3] = { vx[2].x(), vx[2].y(), vx[2].z() };
    real_t         res_re[ VECTOR_SIZE ];
    real_t         res_im[ VECTOR_SIZE ];
    
    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        // load quadrature points and weight for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // compute d = y - x; <d,d>
        //

        vpacked  diff[3];

        diff[0] = sub( y[0], muladd( a0, x0[0], muladd( a1, x1[0], mul( a2, x2[0] ) ) ) );
        diff[1] = sub( y[1], muladd( a0, x0[1], muladd( a1, x1[1], mul( a2, x2[1] ) ) ) );
        diff[2] = sub( y[2], muladd( a0, x0[2], muladd( a1, x1[2], mul( a2, x2[2] ) ) ) );

        // |x-y|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |x-y|
        const vpacked  vdist  = sqrt( vdist2 );

        //
        //   ( i·κ·|x-y| )
        // e
        // ───────────────
        //      |x-y|
        //
        // compute exp( i·κ·|x-y| ) = polar( exp( re ), im )
        //                          = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // exp( re( ( i·κ·d[0], i·κ·d[1] ) ) ) = 1 since re() = 0
        // (also: multiply with 1/(4·π) / |x-y| as this would be done later anyway)
        const vpacked  vexp_ikd_re = div( vONE_OVER_4PI, vdist ); 

        // im( ( i·κ·d[0], i·κ·d[1] ) )
        const vpacked  vikd_im     = mul( vikappa_im, vdist ); 

        // cos( i·κ·d[0..1] ), sin( i·κ·d[0..1] )
        vpacked        vc, vs;
        
        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()
        // combine double and imaginary parts
        store( mul( vexp_ikd_re, vc ), res_re );
        store( mul( vexp_ikd_re, vs ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ k+j ] = value_t( res_re[j], res_im[j] );
    }// for
}

template < typename value_t,
           typename packed_t >
void
helmholtz_slp_eval_dx_im_simd ( const value_t &           ikappa,
                                const tri_quad_rule_t< real_type_t< value_t > > &  quad_rule,
                                const T3Point             vx[3],
                                const T3Point &           vy,
                                std::vector< value_t > &  values )
{
    DEFINE_LOCAL_TYPE;

    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const vpacked  vONE( real_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    const vpacked  vikappa_re( ikappa.real() ); // κ imaginary -> i·κ double
    const size_t   npts  = quad_rule.npts;
    const vpacked  y[3]  = { vy.x(),    vy.y(),    vy.z() };
    const vpacked  x0[3] = { vx[0].x(), vx[0].y(), vx[0].z() };
    const vpacked  x1[3] = { vx[1].x(), vx[1].y(), vx[1].z() };
    const vpacked  x2[3] = { vx[2].x(), vx[2].y(), vx[2].z() };
    real_t         res_re[ VECTOR_SIZE ];
    
    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        // load quadrature points and weight for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // compute d = y - x; <d,d>
        //

        vpacked  diff[3];

        diff[0] = sub( y[0], muladd( a0, x0[0], muladd( a1, x1[0], mul( a2, x2[0] ) ) ) );
        diff[1] = sub( y[1], muladd( a0, x0[1], muladd( a1, x1[1], mul( a2, x2[1] ) ) ) );
        diff[2] = sub( y[2], muladd( a0, x0[2], muladd( a1, x1[2], mul( a2, x2[2] ) ) ) );

        // |x-y|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |x-y|
        const vpacked  vdist  = sqrt( vdist2 );

        //
        //   ( i·κ·|x-y| )
        // e
        // ───────────────
        //      |x-y|
        //
        // compute exp( i·κ·|x-y| ) = polar( exp( re ), im )
        //                          = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // exp( re( ( i·κ·d[0], i·κ·d[1] ) ) )
        // (also: multiply with 1/(4·π) / |x-y| as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( vONE_OVER_4PI, vdist ), exp( mul( vikappa_re, vdist ) ) ); 

        // exp(re)·cos(im), exp(re)·sin(im) = ( exp(re), 0 ) since im = 0
        // (here we would multiply with 1/(4·π) / |x-y|)

        // combine double and imaginary parts
        store( vexp_ikd_re, res_re );
        
        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ k+j ] = value_t( res_re[j], real_t(0) );
    }// for
}

template < typename value_t,
           typename packed_t >
void
helmholtz_dlp_eval_dy_simd ( const value_t &           ikappa,
                             const tri_quad_rule_t< real_type_t< value_t > > &            quad_rule,
                             const T3Point &           vx,
                             const T3Point             vy[3],
                             const T3Point &           normal,
                             std::vector< value_t > &  values )
{
    DEFINE_LOCAL_TYPE;

    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const vpacked  vONE( real_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    const vpacked  vikappa_re( ikappa.real() );
    const vpacked  vikappa_im( ikappa.imag() );
    const size_t   npts  = quad_rule.npts;
    const vpacked  x[3]  = { vx.x(),     vx.y(),     vx.z() };
    const vpacked  y0[3] = { vy[0].x(),  vy[0].y(),  vy[0].z() };
    const vpacked  y1[3] = { vy[1].x(),  vy[1].y(),  vy[1].z() };
    const vpacked  y2[3] = { vy[2].x(),  vy[2].y(),  vy[2].z() };
    const vpacked  n[3]  = { normal.x(), normal.y(), normal.z() };
    real_t         res_re[ VECTOR_SIZE ];
    real_t         res_im[ VECTOR_SIZE ];

    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        // load quadrature points and weight for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // compute d = y - x; <d,d>; <n,d>
        //

        vpacked  diff[3];

        diff[0] = sub( muladd( a0, y0[0], muladd( a1, y1[0], mul( a2, y2[0] ) ) ), x[0] );
        diff[1] = sub( muladd( a0, y0[1], muladd( a1, y1[1], mul( a2, y2[1] ) ) ), x[1] );
        diff[2] = sub( muladd( a0, y0[2], muladd( a1, y1[2], mul( a2, y2[2] ) ) ), x[2] );

        // |x-y|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |x-y|
        const vpacked  vdist  = sqrt( vdist2 );

        // |x-y|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        // < n, y-x >
        const vpacked  vndu   = muladd( n[0], diff[0], muladd( n[1], diff[1], mul( n[2], diff[2] ) ) );

        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n,y-x>
        // ────────────────────────────────────
        //                |x-y|³
        //
        // compute exp( i·κ·|x-y| ) = polar( exp( re ), im )
        //                          = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        vpacked        vikd_re = mul( vikappa_re, vdist );
        
        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist );
        
        // exp( re(·) )
        // (also: multiply with (1/4·π) <n,y-x> / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( mul( vONE_OVER_4PI, vndu ), vdist3 ), exp( vikd_re ) ); 

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        vpacked  vc, vs;

        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()   (here we would multiply with <n,y-x> / |x-y|³)
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        // i·κ·|x-y| - 1
        vikd_re = sub( vikd_re, vONE );
        
        // compute exp() · (i·κ·|x-y| - 1) and store final values
        store( mulsub( vc, vikd_re, mul( vs, vikd_im ) ), res_re );
        store( muladd( vc, vikd_im, mul( vs, vikd_re ) ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ k+j ] = value_t( res_re[j], res_im[j] );
    }// for
}

template < typename value_t,
           typename packed_t >
void
helmholtz_dlp_eval_dy_re_simd ( const value_t &           ikappa,
                                const tri_quad_rule_t< real_type_t< value_t > > & quad_rule,
                                const T3Point &           vx,
                                const T3Point             vy[3],
                                const T3Point &           normal,
                                std::vector< value_t > &  values )
{
    DEFINE_LOCAL_TYPE;

    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const vpacked  vONE( real_t(1) );
    const vpacked  vMINUS_ONE( real_t(-1) );  // used for: re(i·κ·|x-y|)-1 since re() = 0
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    const vpacked  vikappa_im( ikappa.imag() );  // κ double -> i·κ imaginary
    const size_t   npts  = quad_rule.npts;
    const vpacked  x[3]  = { vx.x(),     vx.y(),     vx.z() };
    const vpacked  y0[3] = { vy[0].x(),  vy[0].y(),  vy[0].z() };
    const vpacked  y1[3] = { vy[1].x(),  vy[1].y(),  vy[1].z() };
    const vpacked  y2[3] = { vy[2].x(),  vy[2].y(),  vy[2].z() };
    const vpacked  n[3]  = { normal.x(), normal.y(), normal.z() };
    real_t         res_re[ VECTOR_SIZE ];
    real_t         res_im[ VECTOR_SIZE ];

    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        // load quadrature points and weight for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // compute d = y - x; <d,d>; <n,d>
        //

        vpacked  diff[3];

        diff[0] = sub( muladd( a0, y0[0], muladd( a1, y1[0], mul( a2, y2[0] ) ) ), x[0] );
        diff[1] = sub( muladd( a0, y0[1], muladd( a1, y1[1], mul( a2, y2[1] ) ) ), x[1] );
        diff[2] = sub( muladd( a0, y0[2], muladd( a1, y1[2], mul( a2, y2[2] ) ) ), x[2] );

        // |x-y|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |x-y|
        const vpacked  vdist  = sqrt( vdist2 );

        // |x-y|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        // < n, y-x >
        const vpacked  vndu   = muladd( n[0], diff[0], muladd( n[1], diff[1], mul( n[2], diff[2] ) ) );

        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n,y-x>
        // ────────────────────────────────────
        //                |x-y|³
        //
        // compute exp( i·κ·|x-y| ) = polar( exp( re ), im )
        //                          = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // im( ( i·κ·d[0], i·κ·d[1] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist );
        
        // exp( re(·) ) = 1 since re() = 0
        // (also: multiply with 1/(4·π) <n,y-x> / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = div( mul( vONE_OVER_4PI, vndu ), vdist3 ); 

        // cos( i·κ·d[0..1] ), sin( i·κ·d[0..1] )
        vpacked  vc, vs;

        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()   (here we would multiply with 1/(4·π) <n,y-x> / |x-y|³)
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        // compute exp() · (i·κ·|x-y| - 1) and store final values
        store( mulsub( vc, vMINUS_ONE, mul( vs, vikd_im    ) ), res_re );
        store( muladd( vc, vikd_im,    mul( vs, vMINUS_ONE ) ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ k+j ] = value_t( res_re[j], res_im[j] );
    }// for
}

template < typename value_t,
           typename packed_t >
void
helmholtz_dlp_eval_dy_im_simd ( const value_t &           ikappa,
                                const tri_quad_rule_t< real_type_t< value_t > > & quad_rule,
                                const T3Point &           vx,
                                const T3Point             vy[3],
                                const T3Point &           normal,
                                std::vector< value_t > &  values )
{
    DEFINE_LOCAL_TYPE;

    constexpr real_t  ONE_OVER_4PI = real_t(1) / (real_t(4) * Math::pi< real_t >());
    
    const vpacked  vONE( real_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    const vpacked  vikappa_re( ikappa.real() ); // κ imaginary -> i·κ double
    const size_t   npts  = quad_rule.npts;
    const vpacked  x[3]  = { vx.x(),     vx.y(),     vx.z() };
    const vpacked  y0[3] = { vy[0].x(),  vy[0].y(),  vy[0].z() };
    const vpacked  y1[3] = { vy[1].x(),  vy[1].y(),  vy[1].z() };
    const vpacked  y2[3] = { vy[2].x(),  vy[2].y(),  vy[2].z() };
    const vpacked  n[3]  = { normal.x(), normal.y(), normal.z() };
    real_t         res_re[ VECTOR_SIZE ];

    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        // load quadrature points and weight for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // compute d = y - x; <d,d>; <n,d>
        //

        vpacked  diff[3];

        diff[0] = sub( muladd( a0, y0[0], muladd( a1, y1[0], mul( a2, y2[0] ) ) ), x[0] );
        diff[1] = sub( muladd( a0, y0[1], muladd( a1, y1[1], mul( a2, y2[1] ) ) ), x[1] );
        diff[2] = sub( muladd( a0, y0[2], muladd( a1, y1[2], mul( a2, y2[2] ) ) ), x[2] );

        // |x-y|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |x-y|
        const vpacked  vdist  = sqrt( vdist2 );

        // |x-y|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        // < n, y-x >
        const vpacked  vndu   = muladd( n[0], diff[0], muladd( n[1], diff[1], mul( n[2], diff[2] ) ) );

        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n,y-x>
        // ────────────────────────────────────
        //                |x-y|³
        //
        // compute exp( i·κ·|x-y| ) = polar( exp( re ), im )
        //                          = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // re( ( i·κ·d[0], i·κ·d[1] ) )
        vpacked        vikd_re = mul( vikappa_re, vdist );
        
        // exp( re(·) )
        // (also: multiply with <n,y-x> / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( mul( vONE_OVER_4PI, vndu ), vdist3 ),
                                          exp( vikd_re ) );

        // exp(re)·cos(im), exp(re)·sin(im) = ( exp(re), 0 ) since im = 0 
        // (here we would multiply with 1/(4·π) <n,y-x> / |x-y|³)

        // i·κ·|x-y| - 1
        vikd_re = sub( vikd_re, vONE );
        
        // compute exp() · (i·κ·|x-y| - 1)  ( only double numbers ! )
        vikd_re = mul( vexp_ikd_re, vikd_re );

        // compute exp() · (i·κ·|x-y| - 1) and store final values
        store( vikd_re, res_re );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ k+j ] = value_t( res_re[j], real_t(0) );
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// explicit template instantiations
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#define  INST_HELMHOLTZ_SLP_SIMD( value_t, ansatzsp_t, testsp_t, suffix ) \
    template void                                                       \
    helmholtz_slp_##suffix< value_t, ansatzsp_t, testsp_t, packed< real_type_t< value_t >, SIMD_ISA > > ( \
        const TGrid::triangle_t &    tri0,                       \
        const TGrid::triangle_t &    tri1,                       \
        const tripair_quad_rule_t *  rule,                       \
        const value_t                ikappa,                     \
        const ansatzsp_t *           ansatz_sp,                  \
        const testsp_t *             test_sp,                    \
        std::vector< value_t > &     values )

#define  INST_SLP( type1, type2 )                                             \
    INST_HELMHOLTZ_SLP_SIMD( type1, TConstFnSpace< type2 >,     TConstFnSpace< type2 >,     simd ); \
    INST_HELMHOLTZ_SLP_SIMD( type1, TLinearFnSpace< type2 >,    TConstFnSpace< type2 >,     simd ); \
    INST_HELMHOLTZ_SLP_SIMD( type1, TConstFnSpace< type2 >,     TLinearFnSpace< type2 >,    simd ); \
    INST_HELMHOLTZ_SLP_SIMD( type1, TLinearFnSpace< type2 >,    TLinearFnSpace< type2 >,    simd ); \
                                                                        \
    INST_HELMHOLTZ_SLP_SIMD( type1, TConstFnSpace< type2 >,     TConstFnSpace< type2 >,     re_simd ); \
    INST_HELMHOLTZ_SLP_SIMD( type1, TLinearFnSpace< type2 >,    TConstFnSpace< type2 >,     re_simd ); \
    INST_HELMHOLTZ_SLP_SIMD( type1, TConstFnSpace< type2 >,     TLinearFnSpace< type2 >,    re_simd ); \
    INST_HELMHOLTZ_SLP_SIMD( type1, TLinearFnSpace< type2 >,    TLinearFnSpace< type2 >,    re_simd ); \
                                                                        \
    INST_HELMHOLTZ_SLP_SIMD( type1, TConstFnSpace< type2 >,     TConstFnSpace< type2 >,     im_simd ); \
    INST_HELMHOLTZ_SLP_SIMD( type1, TLinearFnSpace< type2 >,    TConstFnSpace< type2 >,     im_simd ); \
    INST_HELMHOLTZ_SLP_SIMD( type1, TConstFnSpace< type2 >,     TLinearFnSpace< type2 >,    im_simd ); \
    INST_HELMHOLTZ_SLP_SIMD( type1, TLinearFnSpace< type2 >,    TLinearFnSpace< type2 >,    im_simd );

INST_SLP( std::complex< float >, float )
INST_SLP( std::complex< double >, double )

INST_HELMHOLTZ_SLP_SIMD( std::complex< double >, TConstEdgeFnSpace, TConstEdgeFnSpace, simd );
INST_HELMHOLTZ_SLP_SIMD( std::complex< double >, TConstEdgeFnSpace, TConstEdgeFnSpace, re_simd );
INST_HELMHOLTZ_SLP_SIMD( std::complex< double >, TConstEdgeFnSpace, TConstEdgeFnSpace, im_simd );

#define  INST_HELMHOLTZ_DLP_SIMD( value_t, ansatzsp_t, testsp_t, suffix ) \
    template void                                                       \
    helmholtz_dlp_##suffix< value_t, ansatzsp_t, testsp_t, packed< real_type_t< value_t >, SIMD_ISA > > ( \
        const idx_t                  tri_id,                            \
        const bool                   adjoint,                           \
        const TGrid::triangle_t &    tri0,                              \
        const TGrid::triangle_t &    tri1,                              \
        const tripair_quad_rule_t *  rule,                              \
        const value_t                ikappa,                            \
        const ansatzsp_t *           ansatz_sp,                         \
        const testsp_t *             test_sp,                           \
        std::vector< value_t > &     values )

#define  INST_DLP( type1, type2 )                                             \
    INST_HELMHOLTZ_DLP_SIMD( type1, TConstFnSpace< type2 >,  TConstFnSpace< type2 >,  simd ); \
    INST_HELMHOLTZ_DLP_SIMD( type1, TLinearFnSpace< type2 >, TConstFnSpace< type2 >,  simd ); \
    INST_HELMHOLTZ_DLP_SIMD( type1, TConstFnSpace< type2 >,  TLinearFnSpace< type2 >, simd ); \
    INST_HELMHOLTZ_DLP_SIMD( type1, TLinearFnSpace< type2 >, TLinearFnSpace< type2 >, simd ); \
                                                                        \
    INST_HELMHOLTZ_DLP_SIMD( type1, TConstFnSpace< type2 >,  TConstFnSpace< type2 >,  re_simd ); \
    INST_HELMHOLTZ_DLP_SIMD( type1, TLinearFnSpace< type2 >, TConstFnSpace< type2 >,  re_simd ); \
    INST_HELMHOLTZ_DLP_SIMD( type1, TConstFnSpace< type2 >,  TLinearFnSpace< type2 >, re_simd ); \
    INST_HELMHOLTZ_DLP_SIMD( type1, TLinearFnSpace< type2 >, TLinearFnSpace< type2 >, re_simd ); \
                                                                        \
    INST_HELMHOLTZ_DLP_SIMD( type1, TConstFnSpace< type2 >,  TConstFnSpace< type2 >,  im_simd ); \
    INST_HELMHOLTZ_DLP_SIMD( type1, TLinearFnSpace< type2 >, TConstFnSpace< type2 >,  im_simd ); \
    INST_HELMHOLTZ_DLP_SIMD( type1, TConstFnSpace< type2 >,  TLinearFnSpace< type2 >, im_simd ); \
    INST_HELMHOLTZ_DLP_SIMD( type1, TLinearFnSpace< type2 >, TLinearFnSpace< type2 >, im_simd );

INST_DLP( std::complex< float >, float )
INST_DLP( std::complex< double >, double )

#define  INST_HELMHOLTZ_SLP_HCA_SIMD( type, suffix )                    \
    template void                                                       \
    helmholtz_slp_eval_dx_##suffix< std::complex< type >, packed< type, SIMD_ISA > > ( \
        const std::complex< type > &           ikappa,                  \
        const tri_quad_rule_t< real_type_t< std::complex< type > > > &  quad_rule, \
        const T3Point                          vx[3],                   \
        const T3Point &                        vy,                      \
        std::vector< std::complex< type > > &  values );

INST_HELMHOLTZ_SLP_HCA_SIMD( float, simd )
INST_HELMHOLTZ_SLP_HCA_SIMD( float, re_simd )
INST_HELMHOLTZ_SLP_HCA_SIMD( float, im_simd )

INST_HELMHOLTZ_SLP_HCA_SIMD( double, simd )
INST_HELMHOLTZ_SLP_HCA_SIMD( double, re_simd )
INST_HELMHOLTZ_SLP_HCA_SIMD( double, im_simd )

#define  INST_HELMHOLTZ_DLP_HCA_SIMD( type, suffix )                    \
    template void                                                       \
    helmholtz_dlp_eval_dy_##suffix< std::complex< type >, packed< type, SIMD_ISA > > ( \
        const std::complex< type > &           ikappa,                  \
        const tri_quad_rule_t< real_type_t< std::complex< type > > > &  quad_rule, \
        const T3Point &                        vx,                      \
        const T3Point                          vy[3],                   \
        const T3Point &                        normal,                  \
        std::vector< std::complex< type > > &  values )

INST_HELMHOLTZ_DLP_HCA_SIMD( float, simd );
INST_HELMHOLTZ_DLP_HCA_SIMD( float, re_simd );
INST_HELMHOLTZ_DLP_HCA_SIMD( float, im_simd );

INST_HELMHOLTZ_DLP_HCA_SIMD( double, simd );
INST_HELMHOLTZ_DLP_HCA_SIMD( double, re_simd );
INST_HELMHOLTZ_DLP_HCA_SIMD( double, im_simd );

}// namespace Hpro

// Local Variables:
// mode: c++
// End:
