#ifndef __HPRO_PACKED_SSE3_HH
#define __HPRO_PACKED_SSE3_HH
//
// Project     : HLIBpro
// File        : packed_sse3.hh
// Description : datatype for packed (vector) operations using SSE3
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

/////////////////////////////////////////////////////////////
//
// SSE3 version of packed type
//

#if defined(__SSE3__)

#include <pmmintrin.h>

#if HPRO_USE_AMDLIBM == 1
extern "C" __m128d  amd_vrf4_expf    ( __m128 );
extern "C" void     amd_vrf4_sincosf ( __m128, __m128 *, __m128 * );
extern "C" __m128d  amd_vrd2_exp     ( __m128d );
extern "C" void     amd_vrd2_sincos  ( __m128d, __m128d *, __m128d * );
#endif

#if HPRO_USE_ACML == 1
extern "C" __m128d  __vrs4_expf      ( __m128d );
extern "C" void     __vrs4_sincosf   ( __m128d, __m128d *, __m128d * );
extern "C" __m128d  __vrd2_exp       ( __m128d );
extern "C" void     __vrd2_sincos    ( __m128d, __m128d *, __m128d * );
#endif

#if HPRO_USE_SVML == 1
extern "C"  __m128   _mm_exp_ps      ( __m128 );
extern "C"  __m128   _mm_sincos_ps   ( __m128 *, __m128 );
extern "C"  __m128d  _mm_exp_pd      ( __m128d );
extern "C"  __m128d  _mm_sincos_pd   ( __m128d *, __m128d );
#endif

#if HPRO_USE_LIBMVEC == 1
extern "C" __m128   _ZGVbN4v_expf ( __m128   x );
extern "C" __m128   _ZGVbN4v_sinf ( __m128   x );
extern "C" __m128   _ZGVbN4v_cosf ( __m128   x );
extern "C" __m128d  _ZGVbN2v_exp  ( __m128d  x );
extern "C" __m128d  _ZGVbN2v_sin  ( __m128d  x );
extern "C" __m128d  _ZGVbN2v_cos  ( __m128d  x );
#endif

namespace Hpro
{

//
// SSE3 functions
//
template <>
struct simd_traits< float, ISA_SSE3 >
{
    // SIMD base type
    using  value_t  = float;

    // SIMD vector type
    using  packed_t = __m128;

    // SIMD instruction set
    enum { isa = ISA_SSE3 };

    // SIMD vector size
    enum { vector_size = 4 };

    //
    // SIMD functions
    //

    static packed_t  zero ()                   { return _mm_set1_ps( 0.0f ); }
    static packed_t  fill ( const value_t  f ) { return _mm_set1_ps( f ); }
    
    static packed_t  load  ( const value_t *  f ) { return _mm_loadu_ps( f ); }
    static void      store ( const packed_t   a,
                             value_t *        f ) { _mm_storeu_ps( f, a ); }
    
    static packed_t  add   ( const packed_t  x,
                             const packed_t  y ) { return _mm_add_ps( x, y ); }
    static packed_t  sub   ( const packed_t  x,
                             const packed_t  y ) { return _mm_sub_ps( x, y ); }
    static packed_t  mul   ( const packed_t  x,
                             const packed_t  y ) { return _mm_mul_ps( x, y ); }
    static packed_t  div   ( const packed_t  x,
                             const packed_t  y ) { return _mm_div_ps( x, y ); }

    static packed_t  muladd     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return add( mul( x, y ), z ); }
    static packed_t  negmuladd  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return sub( z, mul( x, y ) ); }
    static packed_t  mulsub     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return sub( mul( x, y ), z ); }
    static packed_t  negmulsub  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return sub( sub( zero(), mul( x, y ) ), z ); }

    static packed_t  sqrt   ( const packed_t  x ) { return _mm_sqrt_ps( x ); }
    
    static packed_t  rsqrt  ( const packed_t  x )
    {
        static const packed_t  vthree = fill( 3.0f );
        static const packed_t  vhalf  = fill( 0.5f );
        packed_t               res;
    
        // r = 1/√(x) in single precision
        res = _mm_rsqrt_ps( x );

        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res,  mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );
    
        // // r = r · 1/2 (3.0 - r·r·x)
        // res = mul( res,  mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );

        return res;
    }
    
    static packed_t  exp    ( const packed_t  x )
    {
        #if HPRO_USE_SVML == 1
    
        return _mm_exp_ps( x );

        #elif HPRO_USE_LIBMVEC == 1
    
        return _ZGVbN4v_expf( x );
    
        #elif HPRO_USE_AMDLIBM == 1
    
        return amd_vrf4_exp( x );

        #elif HPRO_USE_ACML == 1
    
        return __vrf4_exp( x );

        #else
    
        //
        // fall back to standard floating point functions
        //
    
        float  sx[4];

        store( x, sx );

        sx[0] = std::exp( sx[0] );
        sx[1] = std::exp( sx[1] );
        sx[2] = std::exp( sx[2] );
        sx[3] = std::exp( sx[3] );

        return load( sx );
    
        #endif
    }

    static void  sincos ( const packed_t   a,
                          packed_t &       s,
                          packed_t &       c )
    {
        #if HPRO_USE_SVML == 1

        s = _mm_sincos_ps( & c, a );
    
        #elif HPRO_USE_LIBMVEC == 1

        s = _ZGVbN4v_sinf( a );
        c = _ZGVbN4v_cosf( a );
    
        #elif HPRO_USE_AMDLIBM == 1
    
        amd_vrf4_sincos( a, & s, & c );
    
        #elif HPRO_USE_ACML == 1
    
        __vrf4_sincos( a, & s, & c );
    
        #else

        //
        // fall back to standard floating point functions
        //
    
        value_t  sa[4], ss[4], sc[4];

        store( a, sa );

        #  if HPRO_HAS_SINCOS == 1    

        Math::sincos( sa[0], ss[0], sc[0] );
        Math::sincos( sa[1], ss[1], sc[1] );
        Math::sincos( sa[2], ss[2], sc[2] );
        Math::sincos( sa[3], ss[3], sc[3] );

        #  else

        ss[0] = std::sin( sa[0] );
        ss[1] = std::sin( sa[1] );
        ss[2] = std::sin( sa[2] );
        ss[3] = std::sin( sa[3] );
        sc[0] = std::cos( sa[0] );
        sc[1] = std::cos( sa[1] );
        sc[2] = std::cos( sa[2] );
        sc[3] = std::cos( sa[3] );

        #  endif
    
        s = load( ss );
        c = load( sc );
    
        #endif
    }
};

template <>
struct packed< float, ISA_SSE3 >
{
    // SIMD base type
    using  value_t  = typename simd_traits< float, ISA_SSE3 >::value_t;

    // SIMD vector type
    using  packed_t = typename simd_traits< float, ISA_SSE3 >::packed_t;

    // SIMD vector size
    enum { isa = simd_traits< float, ISA_SSE3 >::isa };

    // SIMD vector size
    enum { vector_size = simd_traits< float, ISA_SSE3 >::vector_size };

    // SIMD data
    packed_t   x;

    // ctors
    packed ()                   : x( simd_traits< float, ISA_SSE3 >::zero() )    {}
    packed ( packed_t       y ) : x( y )                                         {}
    packed ( const value_t  f ) : x( simd_traits< float, ISA_SSE3 >::fill( f ) ) {}
    packed ( const value_t  a,
             const value_t  b,
             const value_t  c,
             const value_t  d ) : x( _mm_setr_ps( a, b, c, d ) ) {}
};

//
// yields { a[0], b[0] }
//
inline
packed< float, ISA_SSE3 >
unpacklo ( const packed< float, ISA_SSE3 >  a,
           const packed< float, ISA_SSE3 >  b )
{
    return _mm_unpacklo_ps( a.x, b.x );
}

//
// yields { a[1], b[1] }
//
inline
packed< float, ISA_SSE3 >
unpackhi ( const packed< float, ISA_SSE3 >  a,
           const packed< float, ISA_SSE3 >  b )
{
    return _mm_unpackhi_ps( a.x, b.x );
}

//
// a·b with complex a,b
//
inline
packed< float, ISA_SSE3 >
zmul ( const packed< float, ISA_SSE3 >  a,
       const packed< float, ISA_SSE3 >  b )
{
    return _mm_addsub_ps( _mm_mul_ps( _mm_unpacklo_ps( a.x, a.x ), b.x ),
                          _mm_mul_ps( _mm_unpackhi_ps( a.x, a.x ), _mm_shuffle_ps( b.x, b.x, 1 ) ) );
}

//
// a/b with complex a,b
//
inline
packed< float, ISA_SSE3 >
zdiv ( const packed< float, ISA_SSE3 >  a,
       const packed< float, ISA_SSE3 >  b )
{
    const packed< float, ISA_SSE3 >  t1( _mm_addsub_ps( _mm_mul_ps( _mm_shuffle_ps( a.x, a.x, 1 ),
                                                                    _mm_unpacklo_ps( b.x, b.x ) ),
                                                         _mm_mul_ps( a.x, _mm_unpackhi_ps( b.x, b.x ) ) ) );
    const packed< float, ISA_SSE3 >  t2( _mm_mul_ps( b.x, b.x ) );

    return _mm_div_ps( _mm_shuffle_ps( t1.x, t1.x, 1 ), _mm_hadd_ps( t2.x, t2.x ) );
}



//
// SSE3 functions
//
template <>
struct simd_traits< double, ISA_SSE3 >
{
    // SIMD base type
    using  value_t  = double;

    // SIMD vector type
    using  packed_t = __m128d;

    // SIMD instruction set
    enum { isa = ISA_SSE3 };

    // SIMD vector size
    enum { vector_size = 2 };

    //
    // SIMD functions
    //

    static packed_t  zero ()                   { return _mm_set1_pd( 0.0 ); }
    static packed_t  fill ( const value_t  f ) { return _mm_set1_pd( f ); }
    
    static packed_t  load  ( const value_t *  f ) { return _mm_loadu_pd( f ); }
    static void      store ( const packed_t   a,
                             value_t *        f ) { _mm_storeu_pd( f, a ); }
    
    static packed_t  add   ( const packed_t  x,
                             const packed_t  y ) { return _mm_add_pd( x, y ); }
    static packed_t  sub   ( const packed_t  x,
                             const packed_t  y ) { return _mm_sub_pd( x, y ); }
    static packed_t  mul   ( const packed_t  x,
                             const packed_t  y ) { return _mm_mul_pd( x, y ); }
    static packed_t  div   ( const packed_t  x,
                             const packed_t  y ) { return _mm_div_pd( x, y ); }

    static packed_t  muladd     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return add( mul( x, y ), z ); }
    static packed_t  negmuladd  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return sub( z, mul( x, y ) ); }
    static packed_t  mulsub     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return sub( mul( x, y ), z ); }
    static packed_t  negmulsub  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return sub( sub( zero(), mul( x, y ) ), z ); }

    static packed_t  sqrt   ( const packed_t  x ) { return _mm_sqrt_pd( x ); }
    
    static packed_t  rsqrt  ( const packed_t  x )
    {
        static const packed_t  vthree = fill( 3.0 );
        static const packed_t  vhalf  = fill( 0.5 );
        packed_t               res;
    
        // r = 1/√(x) in single precision
        res = _mm_cvtps_pd( _mm_rsqrt_ps( _mm_cvtpd_ps( x ) ) );

        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res,  mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );
    
        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res,  mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );

        return res;
    }
    
    static packed_t  exp    ( const packed_t  x )
    {
        #if HPRO_USE_SVML == 1
    
        return _mm_exp_pd( x );

        #elif HPRO_USE_LIBMVEC == 1
    
        return _ZGVbN2v_exp( x );
    
        #elif HPRO_USE_AMDLIBM == 1
    
        return amd_vrd2_exp( x );

        #elif HPRO_USE_ACML == 1
    
        return __vrd2_exp( x );

        #else
    
        //
        // fall back to standard floating point functions
        //
    
        double  sx[2];

        store( x, sx );

        sx[0] = std::exp( sx[0] );
        sx[1] = std::exp( sx[1] );

        return load( sx );
    
        #endif
    }

    static void  sincos ( const packed_t   a,
                          packed_t &       s,
                          packed_t &       c )
    {
        #if HPRO_USE_SVML == 1

        s = _mm_sincos_pd( & c, a );
    
        #elif HPRO_USE_LIBMVEC == 1

        s = _ZGVbN2v_sin( a );
        c = _ZGVbN2v_cos( a );
    
        #elif HPRO_USE_AMDLIBM == 1
    
        amd_vrd2_sincos( a, & s, & c );
    
        #elif HPRO_USE_ACML == 1
    
        __vrd2_sincos( a, & s, & c );
    
        #else

        //
        // fall back to standard floating point functions
        //
    
        value_t  sa[2], ss[2], sc[2];

        store( a, sa );

        #  if HPRO_HAS_SINCOS == 1    

        Math::sincos( sa[0], ss[0], sc[0] );
        Math::sincos( sa[1], ss[1], sc[1] );

        #  else

        ss[0] = std::sin( sa[0] );
        ss[1] = std::sin( sa[1] );
        sc[0] = std::cos( sa[0] );
        sc[1] = std::cos( sa[1] );

        #  endif
    
        s = load( ss );
        c = load( sc );
    
        #endif
    }
};

template <>
struct packed< double, ISA_SSE3 >
{
    // SIMD base type
    using  value_t  = typename simd_traits< double, ISA_SSE3 >::value_t;

    // SIMD vector type
    using  packed_t = typename simd_traits< double, ISA_SSE3 >::packed_t;

    // SIMD vector size
    enum { isa = simd_traits< double, ISA_SSE3 >::isa };

    // SIMD vector size
    enum { vector_size = simd_traits< double, ISA_SSE3 >::vector_size };

    // SIMD data
    packed_t   x;

    // ctors
    packed ()                     : x( simd_traits< double, ISA_SSE3 >::zero() )   {}
    packed ( packed_t         y ) : x( y )          {}
    packed ( const value_t    f ) : x( simd_traits< double, ISA_SSE3 >::fill( f ) ) {}
    packed ( const value_t    a,
             const value_t    b ) : x( _mm_setr_pd( a, b ) ) {}
    packed ( const value_t *  v ) : x( simd_traits< double, ISA_SSE3 >::load( v ) ) {}
};

//
// yields { a[0], b[0] }
//
inline
packed< double, ISA_SSE3 >
unpacklo ( const packed< double, ISA_SSE3 >  a,
           const packed< double, ISA_SSE3 >  b )
{
    return _mm_unpacklo_pd( a.x, b.x );
}

//
// yields { a[1], b[1] }
//
inline
packed< double, ISA_SSE3 >
unpackhi ( const packed< double, ISA_SSE3 >  a,
           const packed< double, ISA_SSE3 >  b )
{
    return _mm_unpackhi_pd( a.x, b.x );
}

//
// a·b with complex a,b
//
inline
packed< double, ISA_SSE3 >
zmul ( const packed< double, ISA_SSE3 >  a,
       const packed< double, ISA_SSE3 >  b )
{
    return _mm_addsub_pd( _mm_mul_pd( _mm_unpacklo_pd( a.x, a.x ), b.x ),
                          _mm_mul_pd( _mm_unpackhi_pd( a.x, a.x ), _mm_shuffle_pd( b.x, b.x, 1 ) ) );
}

//
// a/b with complex a,b
//
inline
packed< double, ISA_SSE3 >
zdiv ( const packed< double, ISA_SSE3 >  a,
       const packed< double, ISA_SSE3 >  b )
{
    const packed< double, ISA_SSE3 >  t1( _mm_addsub_pd( _mm_mul_pd( _mm_shuffle_pd( a.x, a.x, 1 ),
                                                                     _mm_unpacklo_pd( b.x, b.x ) ),
                                                         _mm_mul_pd( a.x, _mm_unpackhi_pd( b.x, b.x ) ) ) );
    const packed< double, ISA_SSE3 >  t2( _mm_mul_pd( b.x, b.x ) );

    return _mm_div_pd( _mm_shuffle_pd( t1.x, t1.x, 1 ), _mm_hadd_pd( t2.x, t2.x ) );
}

}// namespace Hpro

#endif // __SSE3__

#endif // __HPRO_PACKED_SSE3_HH

