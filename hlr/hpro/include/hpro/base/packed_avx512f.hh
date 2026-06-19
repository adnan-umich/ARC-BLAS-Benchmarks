#ifndef __HPRO_PACKED_AVX512F_HH
#define __HPRO_PACKED_AVX512F_HH
//
// Project     : HLIBpro
// File        : packed_avx512f.hh
// Description : datatype for packed (vector) operations using AVX512F
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

/////////////////////////////////////////////////////////////
//
// AVX512 version of packed type
//

#if defined(__AVX512F__)

#include <immintrin.h>

#if HPRO_USE_AMDLIBM == 1
extern "C" void  amd_vrsa_expf    ( int  n, float *   f, float *   res );
extern "C" void  amd_vrsa_sincosf ( int  n, float *   f, float *   res_s, float *   res_c );
extern "C" void  amd_vrda_exp     ( int  n, double *  f, double *  res );
extern "C" void  amd_vrda_sincos  ( int  n, double *  f, double *  res_s, double *  res_c );
#endif

#if HPRO_USE_ACML == 1
extern "C" void  vrsa_expf        ( int  n, float *   f, float *   res );
extern "C" void  vrsa_sincosf     ( int  n, float *   f, float *   res_s, float *   res_c );
extern "C" void  vrda_exp         ( int  n, double *  f, double *  res );
extern "C" void  vrda_sincos      ( int  n, double *  f, double *  res_s, double *  res_c );
#endif

#if HPRO_USE_SVML == 1
extern "C"  __m512   _mm512_exp_ps     ( __m512 );
extern "C"  __m512   _mm512_sincos_ps  ( __m512 *, __m512 );
extern "C"  __m512d  _mm512_exp_pd     ( __m512d );
extern "C"  __m512d  _mm512_sincos_pd  ( __m512d *, __m512d );
#endif

#if HPRO_USE_LIBMVEC == 1
extern "C" __m512   _ZGVeN16v_expf  ( __m512   x );
extern "C" __m512   _ZGVeN16v_sinf  ( __m512   x );
extern "C" __m512   _ZGVeN16v_cosf  ( __m512   x );
extern "C" __m512d  _ZGVeN8v_exp    ( __m512d  x );
extern "C" __m512d  _ZGVeN8v_sin    ( __m512d  x );
extern "C" __m512d  _ZGVeN8v_cos    ( __m512d  x );
#endif

namespace Hpro
{

template <>
struct simd_traits< float, ISA_AVX512F >
{
    // SIMD base type
    using  value_t  = float;

    // SIMD vector type
    using  packed_t = __m512;

    // SIMD instruction set
    enum { isa = ISA_AVX512F };

    // SIMD vector size
    enum { vector_size = 16 };

    //
    // SIMD functions
    //

    static packed_t  zero ()                   { return _mm512_setzero_ps(); }
    static packed_t  fill ( const value_t  f ) { return _mm512_set1_ps( f ); }
           
    static packed_t  load  ( const value_t *  f ) { return _mm512_loadu_ps( f ); }
    static void      store ( const packed_t   a,
                             value_t *        f ) { _mm512_storeu_ps( f, a ); }
           
    static packed_t  add   ( const packed_t  x,
                             const packed_t  y ) { return _mm512_add_ps( x, y ); }
    static packed_t  sub   ( const packed_t  x,
                             const packed_t  y ) { return _mm512_sub_ps( x, y ); }
    static packed_t  mul   ( const packed_t  x,
                             const packed_t  y ) { return _mm512_mul_ps( x, y ); }
    static packed_t  div   ( const packed_t  x,
                             const packed_t  y ) { return _mm512_div_ps( x, y ); }
           
    static packed_t  muladd     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return _mm512_fmadd_ps( x, y, z ); }
    static packed_t  negmuladd  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return _mm512_fnmadd_ps( x, y, z ); }
    static packed_t  mulsub     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return _mm512_fmsub_ps( x, y, z ); }
    static packed_t  negmulsub  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return _mm512_fnmsub_ps( x, y, z ); }
           
    static packed_t  sqrt   ( const packed_t  x ) { return _mm512_sqrt_ps( x ); }

    static packed_t  rsqrt  ( const packed_t  x )
    {
        #if HPRO_USE_SVML == 1

        return _mm512_invsqrt_ps( x );

        #else
        
        static const packed_t  vthree = fill( 3.0f );
        static const packed_t  vhalf  = fill( 0.5f );
        packed_t               res;
    
        // r = 1/√(x) in single precision
        res = _mm512_rsqrt14_ps( x );

        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res, mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );
    
        // // r = r · 1/2 (3.0 - r·r·x)
        // res = mul( res, mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );

        return res;

        #endif
    }
    
    static packed_t  exp    ( const packed_t  x )
    {
        #if HPRO_USE_SVML == 1
    
        return _mm512_exp_ps( x );
    
        #elif HPRO_USE_LIBMVEC == 1
    
        return _ZGVeN16v_expf( x );
    
        #elif HPRO_USE_AMDLIBM == 1

        packed_t  res;
    
        amd_vrsa_expf( vector_size,
                       reinterpret_cast< value_t * >( const_cast< packed_t * >( & x ) ),
                       reinterpret_cast< value_t * >( & res ) );

        return res;
    
        #elif HPRO_USE_ACML == 1

        packed_t  res;
    
        vrsa_expf( vector_size,
                   reinterpret_cast< value_t * >( const_cast< packed_t * >( & x ) ),
                   reinterpret_cast< value_t * >( & res ) );

        return res;
    
        #else

        //
        // fall back to standard floating point functions
        //
    
        value_t  sx[ vector_size ];

        store( x, sx );

        for ( int i = 0; i < vector_size; ++i )
            sx[i] = std::exp( sx[i] );

        return load( sx );
    
        #endif
    }
           
    static void      sincos ( const packed_t   a,
                              packed_t &       s,
                              packed_t &       c )
    {
        #if HPRO_USE_SVML == 1

        s = _mm512_sincos_ps( & c, a );
    
        #elif HPRO_USE_LIBMVEC == 1

        s = _ZGVeN16v_sinf( a );
        c = _ZGVeN16v_cosf( a );
    
        #elif HPRO_USE_AMDLIBM == 1

        amd_vrda_sincos( vector_size,
                         reinterpret_cast< value_t * >( const_cast< packed_t * >( & a ) ),
                         reinterpret_cast< value_t * >( & s ),
                         reinterpret_cast< value_t * >( & c ) );
    
        #elif HPRO_USE_ACML == 1

        vrda_sincos( vector_size,
                     reinterpret_cast< value_t * >( const_cast< packed_t * >( & a ) ),
                     reinterpret_cast< value_t * >( & s ),
                     reinterpret_cast< value_t * >( & c ) );
    
        #else

        //
        // fall back to standard floating point functions
        //
    
        value_t  sa[ vector_size ];
        value_t  ss[ vector_size ];
        value_t  sc[ vector_size ];

        store( a, sa );

        #  if HPRO_HAS_SINCOS == 1    

        for ( int i = 0; i < vector_size; ++i )
            Math::sincos( sa[i], ss[i], sc[i] );

        #  else

        for ( int i = 0; i < vector_size; ++i )
            ss[i] = std::sin( sa[i] );

        for ( int i = 0; i < vector_size; ++i )
            sc[i] = std::cos( sa[i] );

        #  endif
    
        s = load( ss );
        c = load( sc );

        #endif
    }
};

template <>
struct packed< float, ISA_AVX512F >
{
    // SIMD base type
    using  value_t  = typename simd_traits< float, ISA_AVX512F >::value_t;

    // SIMD vector type
    using  packed_t = typename simd_traits< float, ISA_AVX512F >::packed_t;

    // SIMD instruction set
    enum { isa = simd_traits< float, ISA_AVX512F >::isa };

    // SIMD vector size
    enum { vector_size = simd_traits< float, ISA_AVX512F >::vector_size };

    // SIMD data
    packed_t   x;

    // ctors
    packed ()                    : x( simd_traits< float, ISA_AVX512F >::zero() )    {}
    packed ( const packed_t  y ) : x( y )                                            {}
    packed ( const value_t   f ) : x( simd_traits< float, ISA_AVX512F >::fill( f ) ) {}
    // packed ( const value_t   a,
    //          const value_t   b,
    //          const value_t   c,
    //          const value_t   d,
    //          const value_t   e,
    //          const value_t   f,
    //          const value_t   g,
    //          const value_t   h ) : x( _mm512_setr_ps( a, b, c, d, e, f, g, h ) ) {}
};



template <>
struct simd_traits< double, ISA_AVX512F >
{
    // SIMD base type
    using  value_t  = double;

    // SIMD vector type
    using  packed_t = __m512d;

    // SIMD instruction set
    enum { isa = ISA_AVX512F };

    // SIMD vector size
    enum { vector_size = 8 };

    //
    // SIMD functions
    //

    static packed_t  zero ()                   { return _mm512_setzero_pd(); }
    static packed_t  fill ( const value_t  f ) { return _mm512_set1_pd( f ); }
           
    static packed_t  load  ( const value_t *  f ) { return _mm512_loadu_pd( f ); }
    static void      store ( const packed_t   a,
                             value_t *        f ) { _mm512_storeu_pd( f, a ); }
           
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

    static packed_t  rsqrt  ( const packed_t  x )
    {
        #if HPRO_USE_SVML == 1

        return _mm512_invsqrt_pd( x );

        #else
        
        static const packed_t  vthree = fill( 3.0 );
        static const packed_t  vhalf  = fill( 0.5 );
        packed_t               res;
    
        // r = 1/√(x) in single precision
        res = _mm512_rsqrt14_pd( x );

        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res, mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );
    
        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res, mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );

        return res;

        #endif
    }
    
    static packed_t  exp    ( const packed_t  x )
    {
        #if HPRO_USE_SVML == 1
    
        return _mm512_exp_pd( x );
    
        #elif HPRO_USE_LIBMVEC == 1
    
        return _ZGVeN8v_exp( x );
    
        #elif HPRO_USE_AMDLIBM == 1

        packed_t  res;
    
        amd_vrda_exp( vector_size,
                      reinterpret_cast< value_t * >( const_cast< packed_t * >( & x ) ),
                      reinterpret_cast< value_t * >( & res ) );

        return res;
    
        #elif HPRO_USE_ACML == 1

        packed_t  res;
    
        vrda_exp( vector_size,
                  reinterpret_cast< value_t * >( const_cast< packed_t * >( & x ) ),
                  reinterpret_cast< value_t * >( & res ) );

        return res;
    
        #else

        //
        // fall back to standard floating point functions
        //
    
        value_t  sx[ vector_size ];

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
           
    static void  sincos ( const packed_t   a,
                          packed_t &       s,
                          packed_t &       c )
    {
        #if HPRO_USE_SVML == 1

        s = _mm512_sincos_pd( & c, a );
    
        #elif HPRO_USE_LIBMVEC == 1

        s = _ZGVeN8v_sin( a );
        c = _ZGVeN8v_cos( a );
    
        #elif HPRO_USE_AMDLIBM == 1

        amd_vrda_sincos( vector_size,
                         reinterpret_cast< value_t * >( const_cast< packed_t * >( & a ) ),
                         reinterpret_cast< value_t * >( & s ),
                         reinterpret_cast< value_t * >( & c ) );
    
        #elif HPRO_USE_ACML == 1

        vrda_sincos( vector_size,
                     reinterpret_cast< value_t * >( const_cast< packed_t * >( & a ) ),
                     reinterpret_cast< value_t * >( & s ),
                     reinterpret_cast< value_t * >( & c ) );
    
        #else

        //
        // fall back to standard floating point functions
        //
    
        value_t  sa[ vector_size ];
        value_t  ss[ vector_size ];
        value_t  sc[ vector_size ];

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
struct packed< double, ISA_AVX512F >
{
    // SIMD base type
    using  value_t  = typename simd_traits< double, ISA_AVX512F >::value_t;

    // SIMD vector type
    using  packed_t = typename simd_traits< double, ISA_AVX512F >::packed_t;

    // SIMD instruction set
    enum { isa = simd_traits< double, ISA_AVX512F >::isa };

    // SIMD vector size
    enum { vector_size = simd_traits< double, ISA_AVX512F >::vector_size };

    // SIMD data
    packed_t   x;

    // ctors
    packed ()                    : x( simd_traits< double, ISA_AVX512F >::zero() )    {}
    packed ( const packed_t  y ) : x( y )                                            {}
    packed ( const value_t   f ) : x( simd_traits< double, ISA_AVX512F >::fill( f ) ) {}
    packed ( const value_t   a,
             const value_t   b,
             const value_t   c,
             const value_t   d,
             const value_t   e,
             const value_t   f,
             const value_t   g,
             const value_t   h ) : x( _mm512_setr_pd( a, b, c, d, e, f, g, h ) ) {}
    packed ( const value_t *  v ) : x( simd_traits< double, ISA_AVX512F >::load( v ) ) {}
};

}// namespace Hpro

#endif // __AVX512F__

#endif // __HPRO_PACKED_AVX512F_HH
