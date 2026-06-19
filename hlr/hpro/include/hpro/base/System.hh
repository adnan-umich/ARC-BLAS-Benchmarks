#ifndef __HPRO_SYSTEM_HH
#define __HPRO_SYSTEM_HH
//
// Project     : HLIBpro
// File        : System.hh
// Description : module containing basic system routines
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <cmath>
#include <string>
#include <limits>

#include "hpro/base/types.hh"
#include "hpro/base/String.hh"

namespace Hpro
{
    
///////////////////////////////////////////////////
//
// mathematical functions
//
///////////////////////////////////////////////////

#if HPRO_USE_AMDLIBM == 1

extern "C" void  amd_sincos  ( double, double *, double * );

#endif

#if HPRO_USE_ACML == 1

extern "C" void  fastsincos  ( double, double *, double * );

#endif

namespace Math
{
//
// constants
//

template < typename T >  constexpr T  pi () { return T(3.141592653589793238462643383279502884L); }

//
// minimum, maximum, intervals
//

// return minimum/maximum of \a f1 and \a f2
template < typename T >  constexpr T  min     ( const T  f1, const T  f2 ) noexcept { return std::min( f1, f2 ); }
template < typename T >  constexpr T  max     ( const T  f1, const T  f2 ) noexcept { return std::max( f1, f2 ); }

// limit \a f to interval [f_min,f_max]
template < typename T >  constexpr T  limit   ( const T  f_min,
                                                const T  f,
                                                const T  f_max ) noexcept { return std::min( f_max, std::max( f_min, f ) ); }

// return true if lb < f < ub
template < typename T >  constexpr bool  inside ( const T  lb,
                                                  const T  f,
                                                  const T  ub ) noexcept 
{
    return (( lb < f ) && ( f < ub ));
}

// return true if lb ≤ f ≤ ub
template < typename T >  constexpr bool  inside_eq ( const T  lb,
                                                     const T  f,
                                                     const T  ub ) noexcept 
{
    return (( lb <= f ) && ( f <= ub ));
}

//
// modulus, sign and roots
//

// return conjugate value
template < typename T >  constexpr T                  conj   ( const T                  f ) { return f; }
template < typename T >  constexpr std::complex< T >  conj   ( const std::complex< T >  f ) { return std::conj( f ); }

// return absolute value of argument
template < typename T >  constexpr T                  abs    ( const T                  f ) { return std::abs( f ); }
template < typename T >  constexpr T                  abs    ( const std::complex< T >  f ) { return std::abs( f ); }

// return signum of argument
template < typename T >  constexpr T                  sign   ( const T                  f ) { return (f == T(0) ? T(0) : (f > T(0) ? T(1) : T(-1))); }
template < typename T >  constexpr std::complex< T >  sign   ( const std::complex< T >  f ) { return f / abs( f ); }

// return square of argument
template < typename T >  constexpr T                  square ( const T                  f ) { return f*f; }

// return square root of argument
template < typename T >  constexpr T                  sqrt   ( const T                  f ) { return std::sqrt(f); }

// return reciprocal square root of argument
template < typename T >  constexpr T                  rsqrt  ( const T                  f ) { return T(1) / sqrt(f); }

//
// power and logarithms
//

// compute x raised to the power of y
template < typename T1,
           typename T2 > constexpr auto               pow    ( const T1                 x,
                                                               const T2                 y ) { return std::pow( x, y ); }

// compute e raised to the power of x
template < typename T >  constexpr T                  exp    ( const T                  f ) { return std::exp( f ); }

// compute base 2 logarithm of integers
uint  log2  ( const uint  n );

// compute natural logarithm
template < typename T >  constexpr T                  log    ( const T                  f ) { return std::log( f ); }

// compute base 10 logarithm
template < typename T >  constexpr T                  log10  ( const T                  f ) { return std::log10( f ); }

//
// rounding
//

// round argument down to the nearest integer
template < typename T >  constexpr T                  floor  ( const T                  f ) { return std::floor( f ); }
    
// round argument up to the nearest integer
template < typename T >  constexpr T                  ceil   ( const T                  f ) { return std::ceil( f ); }
    
//
// trigonometry
//

// return sine and cosing of argument
template < typename T >  constexpr T                  sin   ( const T                  f ) { return std::sin( f ); }
template < typename T >  constexpr T                  cos   ( const T                  f ) { return std::cos( f ); }

// return arc sine and cosing of argument
template < typename T >  constexpr T                  asin  ( const T                  f ) { return std::asin( f ); }
template < typename T >  constexpr T                  acos  ( const T                  f ) { return std::acos( f ); }

// simultaneously compute sine and cosine
inline void
sincos  ( const float  f,
          float &      s,
          float &      c )
{
    #if HPRO_USE_SVML == 1
        
    ::sincosf( f, & s, & c );
        
    #elif HPRO_USE_AMDLIBM == 1
        
    amd_sincos( f, & s, & c );
        
    #elif HPRO_USE_ACML == 1
        
    fastsincos( f, & s, & c );
        
    #elif HPRO_HAS_SINCOS == 1
        
    ::sincosf( f, & s, & c );
        
    #else
    s = sin( f );
    c = cos( f );
    #endif
}

// simultaneously compute sine and cosine
inline void
sincos  ( const double f,
          double &     s,
          double &     c )
{
    #if HPRO_USE_SVML == 1
        
    ::sincos( f, & s, & c );
        
    #elif HPRO_USE_AMDLIBM == 1
        
    amd_sincos( f, & s, & c );
        
    #elif HPRO_USE_ACML == 1
        
    fastsincos( f, & s, & c );
        
    #elif HPRO_HAS_SINCOS == 1
        
    ::sincos( f, & s, & c );
        
    #else
    s = sin( f );
    c = cos( f );
    #endif
}

//
// check of values
//

// return true if given value contains Inf
template <typename T> bool is_inf ( const T val );

// return true if given value contains NaN
template <typename T> bool is_nan ( const T val );

//
// misc.
//

// compute givens rotation (cs,sn) for given <a>, <b>
void givens ( const float  a, const float  b, float &  cs, float &  sn );
void givens ( const double a, const double b, double & cs, double & sn );

void givens ( const std::complex<float> &   a,
              const std::complex<float> &   b,
              std::complex<float> &         cs,
              std::complex<float> &         sn );
void givens ( const std::complex<double> &  a,
              const std::complex<double> &  b,
              std::complex<double> &        cs,
              std::complex<double> &        sn );

}// namespace Math

///////////////////////////////////////////////////
//
// information about types
//
///////////////////////////////////////////////////

namespace Limits
{
    
// minimal representable value
template <class T>  constexpr T  min      () { return std::numeric_limits< T >::min(); }
    
// maximal representable value
template <class T>  constexpr T  max      () { return std::numeric_limits< T >::max(); }
    
// return difference between 1 and smallest value bigger than 1
template <class T>  constexpr T  epsilon  () { return std::numeric_limits< T >::epsilon(); }

// return NaN value
template <class T>            T  nan      ();

}// namespace Limits

///////////////////////////////////////////////////
//
// Time related types and functions
//
///////////////////////////////////////////////////

namespace Time
{
//
// helper function to print timings
//
void  __autotimer_print ( const std::string &  format,
                          double               sec );

//
// stores time durations
//
struct TDuration
{
    double  value;
    
    // default ctor
    TDuration ( const double  val = 0.0 ) : value(val) {}
    
    // convert to various units
    double  hour     () const { return value / 3600.0; }
    double  minute   () const { return value / 60.0; }
    double  seconds  () const { return value; }
    double  millisec () const { return value * 1000.0; }
    double  microsec () const { return value * 1000000.0; }
    
    // convert to double (seconds)
    operator double () const { return value; }
    
    // convert to string
    std::string  to_string     () const;

    // convert to string (extended version)
    std::string  to_string_ext () const;
    
    //! stream output
    friend std::ostream & operator << ( std::ostream & os, const TDuration & t )
    {
        return os << t.to_string();
    }
};
    
inline
TDuration
operator + ( const TDuration  d1,
             const TDuration  d2 )
{
    return d1.value + d2.value;
}
    
//
// stores time points
//
template < int TIME_TYPE >
struct TBaseTimePoint
{
    double  value;
        
    // default ctor
    TBaseTimePoint ( const double  val = 0.0 ) : value(val) {}
};
    
template < int TIME_TYPE >
inline TDuration
operator - ( const TBaseTimePoint< TIME_TYPE >  t1,
             const TBaseTimePoint< TIME_TYPE >  t2 )
{
    return t1.value - t2.value;
}
    
enum TTimeType
{
    PROCESS_TIME,
    THREAD_TIME,
    WALL_TIME
};
    
//! return CPU time since start of program
double  cpu_time         ();

//! return CPU time of current thread since thread start
double  cpu_time_thread  ();

//! return current wall time in seconds
double  wall_time        ();

//! measure and print time spent in current block
#define DEF_AUTO_TIMER                                      \
    class TAutoTimer                                        \
    {                                                       \
    private:                                                \
        std::string  _format;                               \
        TTimePoint   _tic;                                  \
                                                            \
    public:                                                 \
        TAutoTimer ( const std::string &  aformat )         \
            : _format( aformat )                            \
            , _tic( now() )                                 \
        {}                                                  \
                                                            \
        ~TAutoTimer ()                                      \
        {                                                   \
            auto  toc = since( _tic );                      \
                                                            \
            __autotimer_print( _format, double( toc ) );    \
        }                                                   \
    }


namespace Process
{
using TTimePoint = TBaseTimePoint< PROCESS_TIME >;

//! return process time since start of program
TTimePoint  now ();

//! return process time since \a t
TDuration   since ( const TTimePoint  t );

DEF_AUTO_TIMER;

}// namespace CPU

namespace Thread
{
using TTimePoint = TBaseTimePoint< THREAD_TIME >;

//! return thread time since start of program
TTimePoint  now ();

//! return thread time since \a t
TDuration   since ( const TTimePoint  t );

DEF_AUTO_TIMER;

}// namespace Thread

namespace Wall
{
using TTimePoint = TBaseTimePoint< WALL_TIME >;

//! return current wall time
TTimePoint  now ();

//! return wall time since \a t
TDuration   since ( const TTimePoint  t );

DEF_AUTO_TIMER;

}// namespace Wall

}// namespace Time

///////////////////////////////////////////////////
//
// cryptographic functions
//
///////////////////////////////////////////////////

namespace Crypt
{

//! compute CRC32 checksum of byte stream stored in \a data with \a size entries
uint crc32 ( const size_t size, const void * data );

//! template version if \see crc32 for automatic type handling
template <typename T>
uint crc32 ( const size_t nelem, const T * data )
{
    return crc32( sizeof(T) * nelem, reinterpret_cast< const void * >( data ) );
}

//! compute Adler32 checksum of byte stream stored in \a data with \a size entries
ulong adler32 ( const size_t size, const void * data );

}// namespace Crypt

///////////////////////////////////////////////////
//
// Support for RTTI (independent of C++)
//
///////////////////////////////////////////////////

namespace RTTI
{

// convert a typename to a unique id
uint         type_to_id        ( const std::string &  type );

// return typename to a given id
std::string  id_to_type        ( const uint           id );

// register type in RTTI system (returns id)
uint         register_type     ( const char *         type );

// print all registered types with ids
void         print_registered  ();

}// namespace RTTI

///////////////////////////////////////////////////
//
// memory management
//
///////////////////////////////////////////////////

namespace Mem
{

//
// return current memory usage in bytes
//
size_t       usage     ();

//
// return maximal memory usage of application so far in bytes
//
size_t       max_usage ();

//
// return string containing pretty printed current memory usage
// with base <base> (base = 2 or 10, e.g. MiB vs MB)
//
std::string  to_string ( const uint   base = 10 );

//
// convert given number of bytes to human readable format
// with base <base> (base = 2 or 10, e.g. MiB vs MB)
//
std::string  to_string ( const size_t bytes,
                         const uint   base = 10 );

//
// wrapper around malloc for statistics (only for debugging)
//
void *       alloc     ( const size_t  size );
    
//
// wrapper around free for statistics (only for debugging)
//
void         free      ( void *        ptr );
    
}// namespace Mem

///////////////////////////////////////////////////
//
// machine properties
//
///////////////////////////////////////////////////

namespace Mach
{

//
// return string with list of CPUs associated to process
//
std::string  cpuset   ();

//
// return hostname
//
std::string  hostname ();

}// namespace Mach

}// namespace Hpro

#endif  // __HPRO_SYSTEM_HH
