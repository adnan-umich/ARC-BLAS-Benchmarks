#ifndef __HPRO_PACKED_AVX_HH
#define __HPRO_PACKED_AVX_HH
//
// Project     : HLIBpro
// File        : packed_sse3.hh
// Description : datatype for packed (vector) operations using AVX
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

/////////////////////////////////////////////////////////////
//
// AVX version of packed type
//

#if defined(__AVX__)

#include <immintrin.h>

#if HPRO_USE_AMDLIBM == 1
extern "C" void  amd_vrda_exp     ( int  n, double *  f, double *  res );
extern "C" void  amd_vrda_sincos  ( int  n, double *  f, double *  res_s, double *  res_c );
#endif

#if HPRO_USE_ACML == 1
extern "C" void  vrda_exp         ( int  n, double *  f, double *  res );
extern "C" void  vrda_sincos      ( int  n, double *  f, double *  res_s, double *  res_c );
#endif

#if HPRO_USE_SVML == 1
extern "C"  __m256d  _mm256_exp_pd     ( __m256d );
extern "C"  __m256d  _mm256_sincos_pd  ( __m256d *, __m256d );
#endif

#if HPRO_USE_LIBMVEC == 1
extern "C" __m256   _ZGVcN8v_expf ( __m256   x );
extern "C" __m256   _ZGVcN8v_sinf ( __m256   x );
extern "C" __m256   _ZGVcN8v_cosf ( __m256   x );
extern "C" __m256d  _ZGVcN4v_exp  ( __m256d  x );
extern "C" __m256d  _ZGVcN4v_sin  ( __m256d  x );
extern "C" __m256d  _ZGVcN4v_cos  ( __m256d  x );
#endif

namespace Hpro
{

template <>
struct simd_traits< float, ISA_AVX >
{
    // SIMD base type
    using  value_t  = float;

    // SIMD vector type
    using  packed_t = __m256;

    // SIMD instruction set
    enum { isa = ISA_AVX };

    // SIMD vector size
    enum { vector_size = 8 };

    //
    // SIMD functions
    //

    static packed_t  zero ()                   { return _mm256_set1_ps( 0.0f ); }
    static packed_t  fill ( const value_t  f ) { return _mm256_set1_ps( f ); }
    
    static packed_t  load  ( const value_t *  f ) { return _mm256_loadu_ps( f ); }
    static void      store ( const packed_t   a,
                             value_t *        f ) { _mm256_storeu_ps( f, a ); }
    
    static packed_t  add   ( const packed_t  x,
                             const packed_t  y ) { return _mm256_add_ps( x, y ); }
    static packed_t  sub   ( const packed_t  x,
                             const packed_t  y ) { return _mm256_sub_ps( x, y ); }
    static packed_t  mul   ( const packed_t  x,
                             const packed_t  y ) { return _mm256_mul_ps( x, y ); }
    static packed_t  div   ( const packed_t  x,
                             const packed_t  y ) { return _mm256_div_ps( x, y ); }

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

    static packed_t  sqrt   ( const packed_t  x ) { return _mm256_sqrt_ps( x ); }
    
    static packed_t  rsqrt  ( const packed_t  x )
    {
        static const packed_t  vthree = fill( 3.0f );
        static const packed_t  vhalf  = fill( 0.5f );
        packed_t               res;
    
        // r = 1/√(x) in single precision
        res = _mm256_rsqrt_ps( x );

        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res, mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );
    
        // // r = r · 1/2 (3.0 - r·r·x)
        // res = mul( res, mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );

        return res;
    }
    
    static packed_t  exp    ( const packed_t  x )
    {
        #if HPRO_USE_SVML == 1
    
        return _mm256_exp_ps( x );
    
        #elif HPRO_USE_LIBMVEC == 1
    
        return _ZGVcN8v_expf( x );
    
        #elif HPRO_USE_AMDLIBM == 1

        packed_t  res;
    
        amd_vrda_exp( 8,
                      reinterpret_cast< float * >( const_cast< packed_t * >( & x ) ),
                      reinterpret_cast< float * >( & res ) );

        return res;
    
        #elif HPRO_USE_ACML == 1

        packed_t  res;
    
        vrda_exp( 8,
                  reinterpret_cast< float * >( const_cast< packed_t * >( & x ) ),
                  reinterpret_cast< float * >( & res ) );

        return res;
    
        #else

        //
        // fall back to standard floating point functions
        //
    
        value_t  sx[8];

        store( x, sx );

        sx[0] = std::exp( sx[0] );
        sx[1] = std::exp( sx[1] );
        sx[2] = std::exp( sx[2] );
        sx[3] = std::exp( sx[3] );
        sx[4] = std::exp( sx[4] );
        sx[5] = std::exp( sx[5] );
        sx[6] = std::exp( sx[6] );
        sx[7] = std::exp( sx[7] );

        return load( sx );
    
        #endif
    }

    static void      sincos ( const packed_t   a,
                              packed_t &       s,
                              packed_t &       c )
    {
        #if HPRO_USE_SVML == 1

        s = _mm256_sincos_ps( & c, a );
    
        #elif HPRO_USE_LIBMVEC == 1

        s = _ZGVcN8v_sinf( a );
        c = _ZGVcN8v_cosf( a );
    
        #elif HPRO_USE_AMDLIBM == 1

        amd_vrda_sincos( 8,
            reinterpret_cast< float * >( const_cast< packed_t * >( & a ) ),
            reinterpret_cast< float * >( & s ),
            reinterpret_cast< float * >( & c ) );
    
        #elif HPRO_USE_ACML == 1

        vrda_sincos( 8,
                     reinterpret_cast< float * >( const_cast< packed_t * >( & a ) ),
                     reinterpret_cast< float * >( & s ),
                     reinterpret_cast< float * >( & c ) );
    
        #else

        //
        // fall back to standard floating point functions
        //
    
        value_t  sa[8], ss[8], sc[8];

        store( a, sa );

        #  if HPRO_HAS_SINCOS == 1    

        Math::sincos( sa[0], ss[0], sc[0] );
        Math::sincos( sa[1], ss[1], sc[1] );
        Math::sincos( sa[2], ss[2], sc[2] );
        Math::sincos( sa[3], ss[3], sc[3] );
        Math::sincos( sa[4], ss[4], sc[4] );
        Math::sincos( sa[5], ss[5], sc[5] );
        Math::sincos( sa[6], ss[6], sc[6] );
        Math::sincos( sa[7], ss[7], sc[7] );

        #  else

        ss[0] = std::sin( sa[0] );
        ss[1] = std::sin( sa[1] );
        ss[2] = std::sin( sa[2] );
        ss[3] = std::sin( sa[3] );
        ss[4] = std::sin( sa[4] );
        ss[5] = std::sin( sa[5] );
        ss[6] = std::sin( sa[6] );
        ss[7] = std::sin( sa[7] );

        sc[0] = std::cos( sa[0] );
        sc[1] = std::cos( sa[1] );
        sc[2] = std::cos( sa[2] );
        sc[3] = std::cos( sa[3] );
        sc[4] = std::cos( sa[4] );
        sc[5] = std::cos( sa[5] );
        sc[6] = std::cos( sa[6] );
        sc[7] = std::cos( sa[7] );

        #  endif
    
        s = load( ss );
        c = load( sc );

        #endif
    }
};

template <>
struct packed< float, ISA_AVX >
{
    // SIMD base type
    using  value_t  = typename simd_traits< float, ISA_AVX >::value_t;

    // SIMD vector type
    using  packed_t = typename simd_traits< float, ISA_AVX >::packed_t;

    // SIMD vector size
    enum { isa = simd_traits< float, ISA_AVX >::isa };

    // SIMD vector size
    enum { vector_size = simd_traits< float, ISA_AVX >::vector_size };

    // SIMD data
    packed_t   x;

    // ctors
    packed ()                    : x( simd_traits< float, ISA_AVX >::zero() )    {}
    packed ( const packed_t  y ) : x( y )                                        {}
    packed ( const value_t   f ) : x( simd_traits< float, ISA_AVX >::fill( f ) ) {}
    packed ( const value_t   a,
             const value_t   b,
             const value_t   c,
             const value_t   d,
             const value_t   e,
             const value_t   f,
             const value_t   g,
             const value_t   h ) : x( _mm256_setr_ps( a, b, c, d, e, f, g, h ) ) {}
};



template <>
struct simd_traits< double, ISA_AVX >
{
    // SIMD base type
    using  value_t  = double;

    // SIMD vector type
    using  packed_t = __m256d;

    // SIMD instruction set
    enum { isa = ISA_AVX };

    // SIMD vector size
    enum { vector_size = 4 };

    //
    // SIMD functions
    //

    static packed_t  zero ()                   { return _mm256_set1_pd( 0.0 ); }
    static packed_t  fill ( const value_t  f ) { return _mm256_set1_pd( f ); }
    
    static packed_t  load  ( const value_t *  f ) { return _mm256_loadu_pd( f ); }
    static void      store ( const packed_t   a,
                             value_t *        f ) { _mm256_storeu_pd( f, a ); }
    
    static packed_t  add   ( const packed_t  x,
                             const packed_t  y ) { return _mm256_add_pd( x, y ); }
    static packed_t  sub   ( const packed_t  x,
                             const packed_t  y ) { return _mm256_sub_pd( x, y ); }
    static packed_t  mul   ( const packed_t  x,
                             const packed_t  y ) { return _mm256_mul_pd( x, y ); }
    static packed_t  div   ( const packed_t  x,
                             const packed_t  y ) { return _mm256_div_pd( x, y ); }

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

    static packed_t  sqrt   ( const packed_t  x ) { return _mm256_sqrt_pd( x ); }
    
    static packed_t  rsqrt  ( const packed_t  x )
    {
        static const packed_t  vthree = fill( 3.0 );
        static const packed_t  vhalf  = fill( 0.5 );
        packed_t               res;
    
        // r = 1/√(x) in single precision
        res = _mm256_cvtps_pd( _mm_rsqrt_ps( _mm256_cvtpd_ps( x ) ) );

        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res, mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );
    
        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res, mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );

        return res;
    }
    
    static packed_t  exp    ( const packed_t  x )
    {
        #if HPRO_USE_SVML == 1
    
        return _mm256_exp_pd( x );
    
        #elif HPRO_USE_LIBMVEC == 1
    
        return _ZGVcN4v_exp( x );
    
        #elif HPRO_USE_AMDLIBM == 1

        packed_t  res;
    
        amd_vrda_exp( 4,
                      reinterpret_cast< double * >( const_cast< packed_t * >( & x ) ),
                      reinterpret_cast< double * >( & res ) );

        return res;
    
        #elif HPRO_USE_ACML == 1

        packed_t  res;
    
        vrda_exp( 4,
                  reinterpret_cast< double * >( const_cast< packed_t * >( & x ) ),
                  reinterpret_cast< double * >( & res ) );

        return res;
    
        #else

        //
        // fall back to standard floating point functions
        //
    
        value_t  sx[4];

        store( x, sx );

        sx[0] = std::exp( sx[0] );
        sx[1] = std::exp( sx[1] );
        sx[2] = std::exp( sx[2] );
        sx[3] = std::exp( sx[3] );

        return load( sx );
    
        #endif
    }

    static void      sincos ( const packed_t   a,
                              packed_t &       s,
                              packed_t &       c )
    {
        #if HPRO_USE_SVML == 1

        s = _mm256_sincos_pd( & c, a );
    
        #elif HPRO_USE_LIBMVEC == 1

        s = _ZGVcN4v_sin( a );
        c = _ZGVcN4v_cos( a );
    
        #elif HPRO_USE_AMDLIBM == 1

        amd_vrda_sincos( 4,
            reinterpret_cast< double * >( const_cast< packed_t * >( & a ) ),
            reinterpret_cast< double * >( & s ),
            reinterpret_cast< double * >( & c ) );
    
        #elif HPRO_USE_ACML == 1

        vrda_sincos( 4,
                     reinterpret_cast< double * >( const_cast< packed_t * >( & a ) ),
                     reinterpret_cast< double * >( & s ),
                     reinterpret_cast< double * >( & c ) );
    
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
struct packed< double, ISA_AVX >
{
    // SIMD base type
    using  value_t  = typename simd_traits< double, ISA_AVX >::value_t;

    // SIMD vector type
    using  packed_t = typename simd_traits< double, ISA_AVX >::packed_t;

    // SIMD vector size
    enum { isa = simd_traits< double, ISA_AVX >::isa };

    // SIMD vector size
    enum { vector_size = simd_traits< double, ISA_AVX >::vector_size };

    // SIMD data
    packed_t   x;

    // ctors
    packed ()                    : x( simd_traits< double, ISA_AVX >::zero() )    {}
    packed ( const packed_t  y ) : x( y )                                         {}
    packed ( const value_t   f ) : x( simd_traits< double, ISA_AVX >::fill( f ) ) {}
    packed ( const value_t   a,
             const value_t   b,
             const value_t   c,
             const value_t   d ) : x( _mm256_setr_pd( a, b, c, d ) ) {}
    packed ( const value_t *  v ) : x( simd_traits< double, ISA_AVX >::load( v ) ) {}
};

}// namespace Hpro

#endif // __AVX__

#endif // __HPRO_PACKED_AVX_HH
