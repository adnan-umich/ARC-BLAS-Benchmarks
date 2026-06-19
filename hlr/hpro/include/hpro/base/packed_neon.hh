#ifndef __HPRO_PACKED_NEON_HH
#define __HPRO_PACKED_NEON_HH
//
// Project     : HLIBpro
// File        : packed_neon.hh
// Description : datatype for packed (vector) operations using NEON
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

/////////////////////////////////////////////////////////////
//
// NEON version of packed type
//

#if defined(__ARM_NEON__)

#include <arm_neon.h>

namespace Hpro
{

//
// NEON functions
//
template <>
struct simd_traits< float, ISA_NEON >
{
    // SIMD base type
    using  value_t  = float;

    // SIMD vector type
    using  packed_t = float32x4_t;

    // SIMD instruction set
    enum { isa = ISA_NEON };

    // SIMD vector size
    enum { vector_size = 4 };

    //
    // SIMD functions
    //

    static packed_t  zero ()                   { return vmovq_n_f32( value_t(0) ); }
    static packed_t  fill ( const value_t  f ) { return vmovq_n_f32( f ); }
    
    static packed_t  load  ( const value_t *  f ) { return vld1q_f32( f ); }
    static void      store ( const packed_t   a,
                             value_t *        f ) { vst1q_f32( f, a ); }
    
    static packed_t  add   ( const packed_t  x,
                             const packed_t  y ) { return vaddq_f32( x, y ); }
    static packed_t  sub   ( const packed_t  x,
                             const packed_t  y ) { return vsubq_f32( x, y ); }
    static packed_t  mul   ( const packed_t  x,
                             const packed_t  y ) { return vmulq_f32( x, y ); }
    static packed_t  div   ( const packed_t  x,
                             const packed_t  y ) { return vdivq_f32( x, y ); }

    static packed_t  muladd     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return vfmaq_f32( x, y, z ); }
    static packed_t  negmuladd  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return sub( z, mul( x, y ) ); }
    static packed_t  mulsub     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return vfmsq_f32( x, y, z ); }
    static packed_t  negmulsub  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return sub( sub( zero(), mul( x, y ) ), z ); }

    static packed_t  sqrt   ( const packed_t  x ) { return vsqrtq_f32( x ); }
    
    static packed_t  rsqrt  ( const packed_t  x )
    {
        static const packed_t  vthree = fill( value_t(3) );
        static const packed_t  vhalf  = fill( value_t(0.5) );
        packed_t               res;
    
        // r = 1/√(x) in single precision
        res = vrsqrteq_f32( x );

        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res,  mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );
    
        // // r = r · 1/2 (3.0 - r·r·x)
        // res = mul( res,  mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );

        return res;
    }
    
    static packed_t  exp    ( const packed_t  x )
    {
        //
        // fall back to standard floating point functions
        //
    
        value_t  sx[ vector_size ];

        store( x, sx );

        sx[0] = std::exp( sx[0] );
        sx[1] = std::exp( sx[1] );
        sx[2] = std::exp( sx[2] );
        sx[3] = std::exp( sx[3] );

        return load( sx );
    }

    static void  sincos ( const packed_t   a,
                          packed_t &       s,
                          packed_t &       c )
    {
        value_t  sa[ vector_size ];
        value_t  ss[ vector_size ];
        value_t  sc[ vector_size ];

        store( a, sa );

        #if HPRO_HAS_SINCOS == 1    

        Math::sincos( sa[0], ss[0], sc[0] );
        Math::sincos( sa[1], ss[1], sc[1] );
        Math::sincos( sa[2], ss[2], sc[2] );
        Math::sincos( sa[3], ss[3], sc[3] );

        #else

        ss[0] = std::sin( sa[0] );
        ss[1] = std::sin( sa[1] );
        ss[2] = std::sin( sa[2] );
        ss[3] = std::sin( sa[3] );
        sc[0] = std::cos( sa[0] );
        sc[1] = std::cos( sa[1] );
        sc[2] = std::cos( sa[2] );
        sc[3] = std::cos( sa[3] );

        #endif
    
        s = load( ss );
        c = load( sc );
    }
};

template <>
struct packed< float, ISA_NEON >
{
    // SIMD base type
    using  value_t  = typename simd_traits< float, ISA_NEON >::value_t;

    // SIMD vector type
    using  packed_t = typename simd_traits< float, ISA_NEON >::packed_t;

    // SIMD vector size
    enum { isa = simd_traits< float, ISA_NEON >::isa };

    // SIMD vector size
    enum { vector_size = simd_traits< float, ISA_NEON >::vector_size };

    // SIMD data
    packed_t   x;

    // ctors
    packed ()                     : x( simd_traits< float, ISA_NEON >::zero() )    {}
    packed ( packed_t         y ) : x( y )                                         {}
    packed ( const value_t    f ) : x( simd_traits< float, ISA_NEON >::fill( f ) ) {}
    packed ( const value_t *  v ) : x( simd_traits< float, ISA_NEON >::load( v ) ) {}
};

}// namespace Hpro

#endif // __ARM_NEON__

#endif // __HPRO_PACKED_NEON_HH

