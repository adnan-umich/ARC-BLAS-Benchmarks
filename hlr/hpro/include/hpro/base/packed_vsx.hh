#ifndef __HPRO_PACKED_VSX_HH
#define __HPRO_PACKED_VSX_HH
//
// Project     : HLIBpro
// File        : packed_vsx.hh
// Description : datatype for packed (vector) operations using VSX
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

/////////////////////////////////////////////////////////////
//
// VSX version of packed type
//

#if defined(__VSX__)

#include <altivec.h>

// VSX exp/sincos from libmass
#define HPRO_USE_MASS  1

#if HPRO_USE_MASS == 1
extern "C" vector double expd2    ( vector double vx );
extern "C" void          sincosd2 ( vector double vx, vector double * vs, vector double * vc );
#endif

//
// VSX functions
//
template <>
struct simd_traits< double, ISA_VSX >
{
    // SIMD base type
    using  value_t  = double;

    // SIMD vector type
    using  packed_t = vector double;

    // SIMD instruction set
    enum { isa = ISA_VSX };

    // SIMD vector size
    enum { vector_size = 2 };

    //
    // SIMD functions
    //

    static packed_t  zero ()                   { return packed_t{ 0, 0 }; }
    static packed_t  fill ( const value_t  f ) { return packed_t{ f, f }; }
    
    static packed_t  load  ( const value_t *  f ) { return packed_t{ f[0], f[1] }; }
    static void      store ( const packed_t   a,
                             value_t *        f ) { f[0] = a[0]; f[1] = a[1]; }
    
    static packed_t  add   ( const packed_t  x,
                             const packed_t  y ) { return vec_add( x, y ); }
    static packed_t  sub   ( const packed_t  x,
                             const packed_t  y ) { return vec_sub( x, y ); }
    static packed_t  mul   ( const packed_t  x,
                             const packed_t  y ) { return vec_mul( x, y ); }
    static packed_t  div   ( const packed_t  x,
                             const packed_t  y ) { return vec_div( x, y ); }

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

    static packed_t  sqrt   ( const packed_t  x ) { return vec_sqrt( x ); }
    
    static packed_t  rsqrt  ( const packed_t  x )
    {
        static const packed_t  vthree = fill( 3.0 );
        static const packed_t  vhalf  = fill( 0.5 );
        packed_t               res;
    
        // r = 1/√(x) in single precision
        res = vec_rsqrte( x );

        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res,  mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );
    
        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res,  mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );

        return res;
    }
    
    static packed_t  exp    ( const packed_t  x )
    {
        #if HPRO_USE_MASS == 1
    
        return expd2( x );

        #else
    
        //
        // fall back to standard floating point functions
        //
    
        value_t  sx[2];

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
        #if HPRO_USE_MASS == 1

        sincosd2( a, & s, & c );
    
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
struct packed< double, ISA_VSX >
{
    // SIMD base type
    using  value_t  = typename simd_traits< double, ISA_VSX >::value_t;

    // SIMD vector type
    using  packed_t = typename simd_traits< double, ISA_VSX >::packed_t;

    // SIMD vector size
    enum { isa = simd_traits< double, ISA_VSX >::isa };

    // SIMD vector size
    enum { vector_size = simd_traits< double, ISA_VSX >::vector_size };

    // SIMD data
    packed_t   x;

    // ctors
    packed ()                     : x( simd_traits< double, ISA_VSX >::zero() )   {}
    packed ( packed_t         y ) : x( y )          {}
    packed ( const value_t    f ) : x( simd_traits< double, ISA_VSX >::fill( f ) ) {}
    packed ( const value_t    a,
             const value_t    b ) : x{ a, b } {}
    packed ( const value_t *  v ) : x( simd_traits< double, ISA_VSX >::load( v ) ) {}
};

//
// yields { a[0], b[0] }
//
inline
packed< double, ISA_VSX >
unpacklo ( const packed< double, ISA_VSX >  a,
           const packed< double, ISA_VSX >  b )
{
    vector unsigned char perm_lo{  0,  1,  2,  3,  4,  5,  6,  7,
                                  16, 17, 18, 19, 20, 21, 22, 23 };
    return vec_perm( a.x, b.x, perm_lo );
}

//
// yields { a[1], b[1] }
//
inline
packed< double, ISA_VSX >
unpackhi ( const packed< double, ISA_VSX >  a,
           const packed< double, ISA_VSX >  b )
{
    vector unsigned char perm_hi{  8,  9, 10, 11, 12, 13, 14, 15,
                                  24, 25, 26, 27, 28, 29, 30, 31 };
    return vec_perm( a.x, b.x, perm_hi );
}

//
// yields { a[1], b[0] }
//
inline
packed< double, ISA_VSX >
unpackhilo ( const packed< double, ISA_VSX >  a,
             const packed< double, ISA_VSX >  b )
{
    vector unsigned char perm_hilo{  8,  9, 10, 11, 12, 13, 14, 15,
                                    16, 17, 18, 19, 20, 21, 22, 23 };
    return vec_perm( a.x, b.x, perm_hilo );
}

//
// a·b with complex a,b
//
inline
packed< double, ISA_VSX >
zmul ( const packed< double, ISA_VSX >  a,
       const packed< double, ISA_VSX >  b )
{
    const packed< double, ISA_VSX >  t1( mul( a, b ) );
    const packed< double, ISA_VSX >  t2( mul( a, unpackhilo( b, b ) ) );
     
    return unpacklo( sub( t1, unpackhilo( t1, t1 ) ),
                     add( t2, unpackhilo( t2, t2 ) ) );
}

//
// a/b with complex a,b
//
inline
packed< double, ISA_VSX >
zdiv ( const packed< double, ISA_VSX >  a,
       const packed< double, ISA_VSX >  b )
{
    const packed< double, ISA_VSX >  t1( mul( a, b ) );
    const packed< double, ISA_VSX >  t2( mul( unpackhilo( a, a ), b ) );
    const packed< double, ISA_VSX >  t3( mul( b, b ) );
        
    return div( unpacklo( add( t1, unpackhilo( t1, t1 ) ),
                          sub( t2, unpackhilo( t2, t2 ) ) ),
                add( t3, unpackhilo( t3, t3 ) ) );
}

#endif // __VSX__

#endif // __HPRO_PACKED_VSX_HH

