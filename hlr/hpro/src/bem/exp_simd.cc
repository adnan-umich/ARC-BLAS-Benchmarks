//
// Project     : HLIBpro
// File        : expsimd.cc
// Description : exponential kernel using SIMD functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#if !defined(SIMD_ISA)
#  error "no SIMD ISA defined"
#endif

#include <vector>

#include "hpro/base/packed.hh"

#include "hpro/bem/TExpBF.hh"

namespace Hpro
{

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// Laplace SLP kernel
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename value_t,
           typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
exp_simd ( const TGrid::triangle_t &    tri0,
           const TGrid::triangle_t &    tri1,
           const tripair_quad_rule_t *  rule,
           std::vector< value_t > &     values,
           const T_ansatzsp *           ansatz_sp,
           const T_testsp *             test_sp )
{
    using vpacked = T_packed;
    using real_t  = real_type_t< value_t >;

    const size_t   VECTOR_SIZE = vpacked::vector_size;
    
    const vpacked  vZERO( value_t(0) );
    const vpacked  vONE(  value_t(1) );

    // load triangle coordinates, all coefficients in t0[i][j]
    // contain coord j (x, y, or z) of vertex i
    const vpacked  t0[3][3] = { { ansatz_sp->grid()->vertex( tri0.vtx[0] )[0],
                                  ansatz_sp->grid()->vertex( tri0.vtx[0] )[1],
                                  ansatz_sp->grid()->vertex( tri0.vtx[0] )[2] },
                                { ansatz_sp->grid()->vertex( tri0.vtx[1] )[0],
                                  ansatz_sp->grid()->vertex( tri0.vtx[1] )[1],
                                  ansatz_sp->grid()->vertex( tri0.vtx[1] )[2] },
                                { ansatz_sp->grid()->vertex( tri0.vtx[2] )[0],
                                  ansatz_sp->grid()->vertex( tri0.vtx[2] )[1],
                                  ansatz_sp->grid()->vertex( tri0.vtx[2] )[2] } };
    
    const vpacked  t1[3][3] = { { test_sp->grid()->vertex( tri1.vtx[0] )[0],
                                  test_sp->grid()->vertex( tri1.vtx[0] )[1],
                                  test_sp->grid()->vertex( tri1.vtx[0] )[2] },
                                { test_sp->grid()->vertex( tri1.vtx[1] )[0],
                                  test_sp->grid()->vertex( tri1.vtx[1] )[1],
                                  test_sp->grid()->vertex( tri1.vtx[1] )[2] },
                                { test_sp->grid()->vertex( tri1.vtx[2] )[0],
                                  test_sp->grid()->vertex( tri1.vtx[2] )[1],
                                  test_sp->grid()->vertex( tri1.vtx[2] )[2] } };

    //
    // eval VECTOR_SIZE add the same time
    //

    const size_t  n_pts = rule->npts;
    const auto    x1    = rule -> template x1< real_t >();
    const auto    y1    = rule -> template y1< real_t >();
    const auto    x2    = rule -> template x2< real_t >();
    const auto    y2    = rule -> template y2< real_t >();
    value_t       res[ VECTOR_SIZE ];

    #pragma ivdep
    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute 1/|x-y|
        //
        
        // load quadrature points for first triangle
        const vpacked  a1 = load< vpacked >( x1 + i );  // ( x1[0], x1[1], ... x1[N] )
        const vpacked  a2 = load< vpacked >( y1 + i );  // ( y1[0], y1[1], ... y1[N] )
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );

        // load quadrature points for second triangle
        const vpacked  b1 = load< vpacked >( x2 + i );  // ( x2[0], x2[1], ... x2[N] )
        const vpacked  b2 = load< vpacked >( y2 + i );  // ( y2[0], y2[1], ... y2[N] )
        const vpacked  b0 = sub( sub( vONE, b1 ), b2 );

        //
        // compute |x-y|²
        //

        vpacked  dot( value_t(0) );
        vpacked  tmp1;
        vpacked  x;

        // u0-v0
        tmp1 =            muladd( a0, t0[0][0], muladd( a1, t0[1][0], mul( a2, t0[2][0] ) ) );
        x    = tmp1;
        tmp1 = sub( tmp1, muladd( b0, t1[0][0], muladd( b1, t1[1][0], mul( b2, t1[2][0] ) ) ) );
        dot  = mul( tmp1, tmp1 );

        // u1-v1
        tmp1 =            muladd( a0, t0[0][1], muladd( a1, t0[1][1], mul( a2, t0[2][1] ) ) );
        tmp1 = sub( tmp1, muladd( b0, t1[0][1], muladd( b1, t1[1][1], mul( b2, t1[2][1] ) ) ) );
        dot  = muladd( tmp1, tmp1, dot );

        // u2-v2
        tmp1 =            muladd( a0, t0[0][2], muladd( a1, t0[1][2], mul( a2, t0[2][2] ) ) );
        tmp1 = sub( tmp1, muladd( b0, t1[0][2], muladd( b1, t1[1][2], mul( b2, t1[2][2] ) ) ) );
        dot  = muladd( tmp1, tmp1, dot );

        // compute exp( - |x-y|) and store in <values>
        store( exp( sub( vZERO, sqrt( dot ) ) ), res );

        // compute x_0 * exp( - |x-y|) and store in <values>
        // store( mul( x, exp( sub( vZERO, sqrt( dot ) ) ) ), res );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[i+j] = value_t(res[j]);
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// HCA functions
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename value_t,
           typename T_packed >
void
exp_eval_dx_simd ( const tri_quad_rule_t< real_type_t< value_t > > & quad_rule,
                   const T3Point                                     vx[3],
                   const T3Point &                                   vy,
                   std::vector< value_t > &                          values )
{
    using vpacked = T_packed;

    const size_t   VECTOR_SIZE = vpacked::vector_size;
    const vpacked  vONE(  value_t(1) );
    const vpacked  vmONE( value_t(-1) );
    
    const size_t   npts  = quad_rule.npts;
    const vpacked  y[3]  = { vy.x(),    vy.y(),    vy.z() };
    const vpacked  x0[3] = { vx[0].x(), vx[0].y(), vx[0].z() };
    const vpacked  x1[3] = { vx[1].x(), vx[1].y(), vx[1].z() };
    const vpacked  x2[3] = { vx[2].x(), vx[2].y(), vx[2].z() };
    value_t        res[ VECTOR_SIZE ];

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
        //   ( -|x-y| )
        // e
        //

        const vpacked  vexp = exp( mul( vmONE, vdist ) );

        store( vexp, res );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ k+j ] = value_t( res[j] );
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// explicit template instantiations
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#define  INST_EXP_SIMD( type, T_ansatzsp, T_testsp )                    \
    template void                                                       \
    exp_simd< type, T_ansatzsp, T_testsp, packed< type, SIMD_ISA > > (   \
        const TGrid::triangle_t &    tri0, \
        const TGrid::triangle_t &    tri1, \
        const tripair_quad_rule_t *  rule, \
        std::vector< type > &        values, \
        const T_ansatzsp *           ansatz_sp, \
        const T_testsp *             test_sp );

#define INST_EXP_ALL( type )                                        \
    INST_EXP_SIMD( type, TConstFnSpace< type >,  TConstFnSpace< type > ) \
    INST_EXP_SIMD( type, TLinearFnSpace< type >, TConstFnSpace< type > )  \
    INST_EXP_SIMD( type, TConstFnSpace< type >,  TLinearFnSpace< type > ) \
    INST_EXP_SIMD( type, TLinearFnSpace< type >, TLinearFnSpace< type > )

INST_EXP_ALL( float )
INST_EXP_ALL( double )

template
void
exp_eval_dx_simd< float, packed< float, SIMD_ISA > > ( const tri_quad_rule_t< float > &  quad_rule,
                                                       const T3Point                     vx[3],
                                                       const T3Point &                   vy,
                                                       std::vector< float > &            values );
template
void
exp_eval_dx_simd< double, packed< double, SIMD_ISA > > ( const tri_quad_rule_t< double > &  quad_rule,
                                                         const T3Point                      vx[3],
                                                         const T3Point &                    vy,
                                                         std::vector< double > &            values );

}// namespace Hpro

// Local Variables:
// mode: c++
// End:
