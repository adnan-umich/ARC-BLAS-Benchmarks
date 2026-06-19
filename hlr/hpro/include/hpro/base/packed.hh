#ifndef __HPRO_PACKED_HH
#define __HPRO_PACKED_HH
//
// Project     : HLIBpro
// File        : packed.hh
// Description : datatype for packed (vector) operations
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <cmath>

#include "hpro/config.h"
#include "hpro/base/basetypes.hh"
#include "hpro/base/System.hh"

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//
// datatype for vector operations
//
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

namespace Hpro
{

//
// defines instruction set 
//
enum isa_t
{
    ISA_FLT  = 0,  // standard, floating point instructions
    ISA_SSE2,      // SSE2 instructions
    ISA_SSE3,      // SSE3 instructions
    ISA_AVX,       // AVX instructions
    ISA_AVX2,      // AVX2 (FMA) instructions
    ISA_MIC,       // MIC instructions
    ISA_AVX512F,   // AVX512F instructions
    ISA_VSX,       // VSX instructions
    ISA_NEON       // NEON instructions
};

inline
const char *
isa_to_string ( isa_t  isa )
{
    switch ( isa )
    {
        case ISA_FLT     : return "float";
        case ISA_SSE2    : return "SSE2";
        case ISA_SSE3    : return "SSE3";
        case ISA_AVX     : return "AVX";
        case ISA_AVX2    : return "AVX2";
        case ISA_MIC     : return "MIC";
        case ISA_AVX512F : return "AVX512F";
        case ISA_VSX     : return "VSX";
        case ISA_NEON    : return "NEON";
        default          : return "unknown";
    }// switch
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//
// auxiliary traits type for actual type and function set
//
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//
// For each ISA, define simd_traits with corresponding functions
//
template < typename T,
           int      ISA >
struct simd_traits
{
    // SIMD base type
    using  value_t  = T;

    // SIMD vector type
    using  packed_t = T;

    // SIMD ISA
    enum { isa = ISA };

    // SIMD vector size
    enum { vector_size = 1 };

    //
    // SIMD functions
    //

    // return zero value
    static packed_t  zero ()                   { return value_t(0); }

    // return variable filled with \a f
    static packed_t  fill ( const value_t  f ) { return f; }

    // return variable with content f[0],f[1], ...
    static packed_t  load  ( const value_t *  f ) { return f[0]; }

    // store variable at array f[0], f[1], ...
    static void      store ( const packed_t  a,
                             value_t *       f )  { f[0] = a; }

    // return x+y
    static packed_t  add   ( const packed_t  x,
                             const packed_t  y ) { return x+y; }
    // return x-y
    static packed_t  sub   ( const packed_t  x,
                             const packed_t  y ) { return x-y; }
    // return x·y
    static packed_t  mul   ( const packed_t  x,
                             const packed_t  y ) { return x*y; }
    // return x/y
    static packed_t  div   ( const packed_t  x,
                             const packed_t  y ) { return x/y; }

    // return x·y + z
    static packed_t  muladd     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return x*y + z; }
    // return -x·y + z
    static packed_t  negmuladd  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return -x*y + z; }
    // return x·y - z
    static packed_t  mulsub     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return x*y - z; }
    // return -x·y - z
    static packed_t  negmulsub  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return -x*y - z; }

    // return sqrt(x)
    static packed_t  sqrt   ( const packed_t  x ) { return std::sqrt( x ); }

    // return 1/sqrt(x)
    static packed_t  rsqrt  ( const packed_t  x ) { return value_t(1) / std::sqrt( x ); }

    // return exp(x)
    static packed_t  exp    ( const packed_t  x ) { return std::exp( x ); }

    // return compute s = sin(f) and c = cos(f)
    static void      sincos ( const packed_t   f,
                              packed_t &       s,
                              packed_t &       c )
    {
        s = std::sin( f );
        c = std::cos( f );
    }
};

//
// packed data type for SIMD operations
//
template < typename T,
           int      ISA >
struct packed
{
    // SIMD base type
    using  value_t  = typename simd_traits<T,ISA>::value_t;

    // SIMD vector type
    using  packed_t = typename simd_traits<T,ISA>::packed_t;

    // SIMD vector size
    enum { vector_size = simd_traits<T,ISA>::vector_size };

    // SIMD instruction set
    enum { isa = simd_traits<T,ISA >::isa };

    // SIMD data
    packed_t  x;
    
    // ctors
    packed ()              : x( value_t(0) ) {}
    packed ( packed_t  y ) : x( y )          {}
};

//
// SIMD functions
// - implemented using functions from simd_traits< T_packed >
//

// yields { f[0], ..  f[?] }
template < typename T_packed >
T_packed
load  ( const typename T_packed::value_t *  f )
{
    return simd_traits< typename T_packed::value_t, T_packed::isa >::load( f );
}

// store { a[0], .. a[?] } into f[0..?]
template < typename T_packed >
void
store ( const T_packed                a,
        typename T_packed::value_t *  f )
{
    simd_traits< typename T_packed::value_t, T_packed::isa >::store( a.x, f );
}

// return sqrt(a)
template < typename T_packed >
T_packed
sqrt ( const T_packed  a )
{
    return simd_traits< typename T_packed::value_t, T_packed::isa >::sqrt( a.x );
}

// return 1/sqrt(a)
template < typename T_packed >
T_packed
rsqrt ( const T_packed  a )
{
    return simd_traits< typename T_packed::value_t, T_packed::isa >::rsqrt( a.x );
}

// return exp(a)
template < typename T_packed >
T_packed
exp ( const T_packed  a )
{
    return simd_traits< typename T_packed::value_t, T_packed::isa >::exp( a.x );
}


// return a+b
template < typename T_packed >
T_packed
add ( const T_packed  a,
      const T_packed  b )
{
    return simd_traits< typename T_packed::value_t, T_packed::isa >::add( a.x, b.x );
}

// return a-b
template < typename T_packed >
T_packed
sub ( const T_packed  a,
      const T_packed  b )
{
    return simd_traits< typename T_packed::value_t, T_packed::isa >::sub( a.x, b.x );
}

// return a*b
template < typename T_packed >
T_packed
mul ( const T_packed  a,
      const T_packed  b )
{
    return simd_traits< typename T_packed::value_t, T_packed::isa >::mul( a.x, b.x );
}

// return a/b
template < typename T_packed >
T_packed
div ( const T_packed  a,
      const T_packed  b )
{
    return simd_traits< typename T_packed::value_t, T_packed::isa >::div( a.x, b.x );
}

// return a*b + c
template < typename T_packed >
T_packed
muladd ( const T_packed  a,
         const T_packed  b,
         const T_packed  c )
{
    return simd_traits< typename T_packed::value_t, T_packed::isa >::muladd( a.x, b.x, c.x );
}

// return -a*b + c  ( = c - a*b )
template < typename T_packed >
T_packed
negmuladd ( const T_packed  a,
            const T_packed  b,
            const T_packed  c )
{
    return simd_traits< typename T_packed::value_t, T_packed::isa >::negmuladd( a.x, b.x, c.x );
}

// return a*b - c
template < typename T_packed >
T_packed
mulsub ( const T_packed  a,
         const T_packed  b,
         const T_packed  c )
{
    return simd_traits< typename T_packed::value_t, T_packed::isa >::mulsub( a.x, b.x, c.x );
}

// return - a*b - c ( = - (a*b + c))
template < typename T_packed >
T_packed
negmulsub ( const T_packed  a,
            const T_packed  b,
            const T_packed  c )
{
    return simd_traits< typename T_packed::value_t, T_packed::isa >::negmulsub( a.x, b.x, c.x );
}

// compute s = sin(f) and c = cos(f)
template < typename T_packed >
void
sincos ( const T_packed  f,
         T_packed &      s,
         T_packed &      c )
{
    simd_traits< typename T_packed::value_t, T_packed::isa >::sincos( f.x, s.x, c.x );
}

}// namespace Hpro

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//
// include special version for different ISAs
//
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

#include "hpro/base/packed_sse2.hh"
#include "hpro/base/packed_sse3.hh"
#include "hpro/base/packed_avx.hh"
#include "hpro/base/packed_avx2.hh"
#include "hpro/base/packed_mic.hh"
#include "hpro/base/packed_avx512f.hh"
#include "hpro/base/packed_vsx.hh"
#include "hpro/base/packed_neon.hh"

#endif  // HPRO_PACKED_HH
