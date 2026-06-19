#ifndef __HPRO_PACKED_MIC_HH
#define __HPRO_PACKED_MIC_HH
//
// Project     : HLIBpro
// File        : packed_mic.hh
// Description : datatype for packed (vector) operations using MIC
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

/////////////////////////////////////////////////////////////
//
// MIC version of packed type
//

#if defined(__MIC__)

#include <immintrin.h>

namespace Hpro
{

template <>
struct simd_traits< double, ISA_MIC >
{
    // SIMD base type
    using  value_t  = double;

    // SIMD vector type
    using  packed_t = __m512d;

    // SIMD instruction set
    enum { isa = ISA_MIC };

    // SIMD vector size
    enum { vector_size = 8 };

    //
    // SIMD functions
    //

    static packed_t  zero ()                   { return _mm512_setzero_pd(); }
    static packed_t  fill ( const value_t  f ) { return _mm512_set1_pd( f ); }
           
    static packed_t  load  ( const value_t *  f ) { return _mm512_load_pd( f ); }
    static void      store ( const packed_t   a,
                             value_t *        f ) { _mm512_store_pd( f, a ); }
           
    static packed_t  add   ( const packed_t  x,
                             const packed_t  y ) { return _mm512_add_pd( x, y ); }
    static packed_t  sub   ( const packed_t  x,
                             const packed_t  y ) { return _mm512_sub_pd( x, y ); }
    static packed_t  mul   ( const packed_t  x,
                             const packed_t  y ) { return _mm512_mul_pd( x, y ); }
    static packed_t  div   ( const packed_t  x,
                             const packed_t  y ) { return _mm512_div_pd( x, y ); }
           
    static packed_t  muladd     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return _mm512_fmadd_pd( x, y, z ); }
    static packed_t  negmuladd  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return _mm512_fnmadd_pd( x, y, z ); }
    static packed_t  mulsub     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return _mm512_fmsub_pd( x, y, z ); }
    static packed_t  negmulsub  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return _mm512_fnmsub_pd( x, y, z ); }
           
    static packed_t  sqrt   ( const packed_t  x ) { return _mm512_sqrt_pd( x ); }
    static packed_t  rsqrt  ( const packed_t  x ) { return _mm512_invsqrt_pd( x ); }
    static packed_t  exp    ( const packed_t  x ) { return _mm512_exp_pd( x ); }
           
    static void      sincos ( const packed_t   f,
                              packed_t &       s,
                              packed_t &       c )
    {
        s = _mm512_sin_pd( f );
        c = _mm512_cos_pd( f );
    }
};

template <>
struct packed< double, ISA_MIC >
{
    // SIMD base type
    using  value_t  = typename simd_traits< double, ISA_MIC >::value_t;

    // SIMD vector type
    using  packed_t = typename simd_traits< double, ISA_MIC >::packed_t;

    // SIMD instruction set
    enum { isa = simd_traits< double, ISA_MIC >::isa };

    // SIMD vector size
    enum { vector_size = simd_traits< double, ISA_MIC >::vector_size };

    // SIMD data
    packed_t   x;

    // ctors
    packed ()                     : x( simd_traits< double, ISA_MIC >::zero() )    {}
    packed ( const packed_t   y ) : x( y )                                         {}
    packed ( const value_t    f ) : x( simd_traits< double, ISA_MIC >::fill( f ) ) {}
    packed ( const value_t    a,
             const value_t    b,
             const value_t    c,
             const value_t    d,
             const value_t    e,
             const value_t    f,
             const value_t    g,
             const value_t    h ) : x( _mm512_setr_pd( a, b, c, d, e, f, g, h ) ) {}
    packed ( const value_t *  v ) : x( simd_traits< double, ISA_MIC >::load( v ) ) {}
};

}// namespace Hpro

#endif // __MIC__

#endif // __HPRO_PACKED_MIC_HH
