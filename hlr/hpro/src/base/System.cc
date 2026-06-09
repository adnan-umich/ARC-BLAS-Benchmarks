//
// Project     : HLIBpro
// File        : System.cc
// Description : module containing basic system routines
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#ifdef __INTEL_COMPILER
#  include <mathimf.h>  // not compatible with math.h, hence include before
#endif

#include "hpro/base/types.hh"

#include <ios>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>
#include <unordered_map>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <string.h>

#include "hpro/base/types.hh"
#include "hpro/base/error.hh"
#include "hpro/base/System.hh"
#include "hpro/parallel/TMutex.hh"

///////////////////////////////////////////////////
//
// mathematical functions
//
///////////////////////////////////////////////////

#include <cmath>

#ifdef SUNOS
#  include <sunmath.h>
#endif

#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
#  include <float.h>
#endif

namespace Hpro
{

namespace Math
{

//
// compute base 2 logarithm of integers
//
uint
log2 ( const uint val ) 
{
    uint  n = val;
    uint  l = 0;

    while ( n > 1 )
    {
        l++;
        n >>= 1;
    }// while
    
    return l;
}

//
// check of values
//

//
// return true if given value contains Inf
//

#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
template <> bool is_inf ( const float   val ) { return ! _finite(val) && ! _isnan(val); }
template <> bool is_inf ( const double  val ) { return ! _finite(val) && ! _isnan(val); }
#elif defined(SUNOS) || defined(__INTEL_COMPILER)
template <> bool is_inf ( const float   val ) { return std::isinf( val ); }
template <> bool is_inf ( const double  val ) { return std::isinf( val ); }
#else
template <> bool is_inf ( const float   val ) { return std::isinf( val ); }
template <> bool is_inf ( const double  val ) { return std::isinf( val ); }
#endif
template <> bool is_inf ( const std::complex<float>  val ) { return ( is_inf<float>(  std::real( val ) ) ||
                                                                      is_inf<float>(  std::imag( val ) ) ); }
template <> bool is_inf ( const std::complex<double> val ) { return ( is_inf<double>( std::real( val ) ) ||
                                                                      is_inf<double>( std::imag( val ) ) ); }

//
// return true if given value contains NaN
//

#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
template <> bool is_nan ( const float   val ) { return _isnan( val ) != 0; }
template <> bool is_nan ( const double  val ) { return _isnan( val ) != 0; }
#elif defined(SUNOS) || defined(__INTEL_COMPILER)
template <> bool is_nan ( const float   val ) { return std::isnan( val ); }
template <> bool is_nan ( const double  val ) { return std::isnan( val ); }
#else
template <> bool is_nan ( const float   val ) { return std::isnan( val ); }
template <> bool is_nan ( const double  val ) { return std::isnan( val ); }
#endif
template <> bool is_nan ( const std::complex<float>  val ) { return ( is_nan<float>(  std::real( val ) ) ||
                                                                      is_nan<float>(  std::imag( val ) ) ); }
template <> bool is_nan ( const std::complex<double> val ) { return ( is_nan<double>( std::real( val ) ) ||
                                                                      is_nan<double>( std::imag( val ) ) ); }

//
// compute givens rotation (cs,sn;-sn,cs) for given <a>, <b>
//

extern "C"
{
// LAPACK routines for Givens rotations
int slartg_ ( float  * f, float  * g, float  * cs, float  * sn, float  * r );
int dlartg_ ( double * f, double * g, double * cs, double * sn, double * r );
int clartg_ ( std::complex<float> *  f,  std::complex<float> *  g, float * cs,
              std::complex<float> *  sn, std::complex<float> *  r );
int zlartg_ ( std::complex<double> * f,  std::complex<double> * g, double * cs,
              std::complex<double> * sn, std::complex<double> * r );
}
    
void givens ( const float a, const float b, float & cs, float & sn )
{
    float  r = 0.0;  // will not be used;
            
    slartg_( const_cast< float * >( & a ), const_cast< float * >( & b ), & cs, & sn, & r );
}

void givens ( const double a, const double b, double & cs, double & sn )
{
    double  r = 0.0;  // will not be used;
            
    dlartg_( const_cast< double * >( & a ), const_cast< double * >( & b ), & cs, & sn, & r );
}

void givens ( const std::complex<float> &  a,
              const std::complex<float> &  b,
              std::complex<float> &        cs,
              std::complex<float> &        sn )
{
    std::complex<float>  r  = 0.0;  // will not be used;
    float           tc = 0.0;
         
    clartg_( const_cast< std::complex<float> * >( & a ), const_cast< std::complex<float> * >( & b ),
             & tc, & sn, & r );

    cs = tc;
}

void givens ( const std::complex<double> &  a,
              const std::complex<double> &  b,
              std::complex<double> &        cs,
              std::complex<double> &        sn )
{
    std::complex<double>  r  = 0.0;  // will not be used;
    double                tc = 0.0;
         
    zlartg_( const_cast< std::complex<double> * >( & a ), const_cast< std::complex<double> * >( & b ),
             & tc, & sn, & r );

    cs = tc;
}

}// namespace Math

}// namespace Hpro

///////////////////////////////////////////////////
//
// information about types
//
///////////////////////////////////////////////////

namespace Hpro
{

namespace Limits
{

//
// return NaN
//

template <> float   nan<float>   () { return std::numeric_limits<float>::quiet_NaN(); }
template <> double  nan<double>  () { return std::numeric_limits<double>::quiet_NaN(); }

template <> std::complex< float >  nan< std::complex< float > >  () { return std::complex< float >(  nan<float>(),  nan<float>() ); }
template <> std::complex< double > nan< std::complex< double > > () { return std::complex< double >( nan<double>(), nan<double>() ); }

}// namespace Limits

}// namespace Hpro
    
///////////////////////////////////////////////////
//
// Time related functions
//
///////////////////////////////////////////////////

#if HPRO_HAS_GETRUSAGE == 1
#include <sys/time.h>
#include <sys/resource.h>
#endif

#if HPRO_HAS_LWPINFO == 1
#include <sys/lwp.h>
#endif

#if HPRO_HAS_GETTIMEOFDAY == 1
#include <sys/time.h>
#include <time.h>
#endif

namespace Hpro
{

namespace Time
{

//
// helper function to print timings
//
void
__autotimer_print ( const std::string &  format,
                    double               sec )
{
    std::cout << boost::format( format.c_str() ) % sec << std::endl;
}

//
// basetype storing time values
//
std::string
TDuration::to_string () const
{
    return str( boost::format( "%.2fs" ) % std::max( 0.0, seconds() ) );
}

//
// basetype storing time values
//
std::string
TDuration::to_string_ext () const
{
    const auto  sec = std::max( 0.0, seconds() );

    if ( sec < 1e-6 )
        return str( boost::format( "%.2fns" ) % ( sec * 1e9 ) );
    else if ( sec < 1e-3 )
        return str( boost::format( "%.2fµs" ) % ( sec * 1e6 ) );
    else if ( sec < 0.1 )
        return str( boost::format( "%.2fms" ) % ( sec * 1e3 ) );
    else if ( sec < 60.0 )
        return str( boost::format( "%.2fs" ) % sec );
    else if ( sec < 3600.0 )
    {
        const auto  tmin = int( std::floor( sec / 60.0 ) );
        const auto  tsec = sec - ( tmin * 60.0 );
        
            
        return str( boost::format( "%d:%.2f" ) % tmin % tsec );
    }// if
    else
    {
        const auto  thour = int( std::floor( sec / 3600.0 ) );
        const auto  rhour = sec - (thour * 3600.0);
        const auto  tmin  = int( std::floor( rhour / 60.0 ) );
        const auto  tsec  = rhour - ( tmin * 60.0 );
        
        return str( boost::format( "%d:%d:%.2f" ) % thour % tmin % tsec );
    }// else
}

//
// return CPU time since start of program
//
double
cpu_time ()
{
#if HPRO_HAS_LWPINFO == 1
            
    struct lwpinfo  lwpinfo_data;
    double          sec, nsec;

    _lwp_info( & lwpinfo_data );
            
    sec   = lwpinfo_data.lwp_utime.tv_sec  + lwpinfo_data.lwp_stime.tv_sec;
    nsec  = lwpinfo_data.lwp_utime.tv_nsec + lwpinfo_data.lwp_stime.tv_nsec;
            
    return sec + 1e-9 * nsec;
            
#elif HPRO_HAS_CLOCKGETTIME == 1
            
    struct timespec ts;

    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, & ts );

    return double(ts.tv_sec) + 1e-9 * double(ts.tv_nsec);
            
#elif HPRO_HAS_GETRUSAGE == 1
            
    struct rusage  rusage_data;
    double         sec, usec;
            
    getrusage( RUSAGE_SELF, & rusage_data );
            
    sec  = double( rusage_data.ru_utime.tv_sec  + rusage_data.ru_stime.tv_sec );
    usec = double( rusage_data.ru_utime.tv_usec + rusage_data.ru_stime.tv_usec );
            
    return sec + 1e-6 * usec;
            
#elif HPRO_HAS_GETPROCESSTIMES == 1
            
    static const HANDLE  pid = GetCurrentProcess();
    FILETIME             createtm, exittm, kerneltm, usertm;
            
    if ( GetProcessTimes( pid, & createtm, & exittm, & kerneltm, & usertm ) )
    {
        // convert kernel and user time from 100 ns to seconds
        // (both are stored in 2 32-bit integers)
        double  sec = 0.0;
                
        sec += double(kerneltm.dwLowDateTime) + double(kerneltm.dwHighDateTime) * 4294967296.0;
        sec += double(usertm.dwLowDateTime)   + double(usertm.dwHighDateTime)   * 4294967296.0;
        sec *= 1e-7;  // convert from 100ns to 1sec
                
        return sec;
    }// if
    else
        return 0.0;
            
#elif HPRO_HAS_GETTICKCOUNT == 1

    //
    // fall back to GetTickCount
    //
            
    const ulong ticks = GetTickCount();
            
    return double(ticks) / 1000.0;

#else

    return 0.0;

#endif
}

//
// return CPU time of current thread since thread start
//
double
cpu_time_thread ()
{
#if HPRO_HAS_LWPINFO == 1
            
    ... // TODO
            
#elif HPRO_HAS_CLOCKGETTIME == 1
            
    struct timespec ts;

    clock_gettime( CLOCK_THREAD_CPUTIME_ID, & ts );

    return double(ts.tv_sec) + 1e-9 * double(ts.tv_nsec);
            
#elif HPRO_HAS_GETRUSAGE == 1
            
    struct rusage  rusage_data;
    double         sec, usec;
            
    getrusage( RUSAGE_THREAD, & rusage_data );
            
    sec  = double( rusage_data.ru_utime.tv_sec  + rusage_data.ru_stime.tv_sec );
    usec = double( rusage_data.ru_utime.tv_usec + rusage_data.ru_stime.tv_usec );
            
    return sec + 1e-6 * usec;
            
#elif HPRO_HAS_GETPROCESSTIMES == 1

    static const HANDLE  pid = GetCurrentThread();
    FILETIME             createtm, exittm, kerneltm, usertm;
            
    if ( GetThreadTimes( pid, & createtm, & exittm, & kerneltm, & usertm ) )
    {
        // convert kernel and user time from 100 ns to seconds
        // (both are stored in 2 32-bit integers)
        double  sec = 0.0;
                
        sec += double(kerneltm.dwLowDateTime) + double(kerneltm.dwHighDateTime) * 4294967296.0;
        sec += double(usertm.dwLowDateTime)   + double(usertm.dwHighDateTime)   * 4294967296.0;
        sec *= 1e-7;  // convert from 100ns to 1sec
                
        return sec;
    }// if
    else
        return 0.0;
            
#elif HPRO_HAS_GETTICKCOUNT == 1

    //
    // fall back to GetTickCount
    //
            
    const ulong ticks = GetTickCount();
            
    return double(ticks) / 1000.0;

#else

    return 0.0;

#endif
}

//
// return current wall time in seconds
//
double
wall_time ()
{
#if HPRO_HAS_CLOCKGETTIME == 1
            
    struct timespec ts;

    clock_gettime( CLOCK_REALTIME, & ts );

    return double(ts.tv_sec) + 1e-9 * double(ts.tv_nsec);
            
#elif HPRO_HAS_GETTIMEOFDAY == 1
            
    struct timeval  timeval_data;
            
    gettimeofday( & timeval_data, nullptr );

    return double(timeval_data.tv_sec) + 1e-6 * double(timeval_data.tv_usec);

#elif HPRO_HAS_GETTICKCOUNT == 1
            
    ulong ticks = GetTickCount();

    return double(ticks) / 1000.0;
            
#else
    return 0.0;
#endif
}

namespace Process
{

// return process time since start of program
TTimePoint  now () { return cpu_time(); }

// return process time since \a t
TDuration   since ( const TTimePoint  t ) { return now() - t; }

}// namespace CPU

namespace Thread
{

// return thread time since start of program
TTimePoint  now () { return cpu_time_thread(); }

// return thread time since \a t
TDuration   since ( const TTimePoint  t ) { return now() - t; }

}// namespace Thread

namespace Wall
{

// return wall time since start of program
TTimePoint  now () { return wall_time(); }

// return wall time since \a t
TDuration   since ( const TTimePoint  t ) { return now() - t; }

}// namespace Wall

}// namespace Time

}// namespace Hpro
    
///////////////////////////////////////////////////
//
// cryptographic functions
//
///////////////////////////////////////////////////

namespace Hpro
{

namespace Crypt
{

//
// table for creating CRC32 hash value
//

static uint crc32_table[256] =  {
    0x00000000, 0x77073096, 0xEE0E612C, 0x990951BA,
    0x076DC419, 0x706AF48F, 0xE963A535, 0x9E6495A3,
    0x0EDB8832, 0x79DCB8A4, 0xE0D5E91E, 0x97D2D988,
    0x09B64C2B, 0x7EB17CBD, 0xE7B82D07, 0x90BF1D91,
    0x1DB71064, 0x6AB020F2, 0xF3B97148, 0x84BE41DE,
    0x1ADAD47D, 0x6DDDE4EB, 0xF4D4B551, 0x83D385C7,
    0x136C9856, 0x646BA8C0, 0xFD62F97A, 0x8A65C9EC,
    0x14015C4F, 0x63066CD9, 0xFA0F3D63, 0x8D080DF5,
    0x3B6E20C8, 0x4C69105E, 0xD56041E4, 0xA2677172,
    0x3C03E4D1, 0x4B04D447, 0xD20D85FD, 0xA50AB56B,
    0x35B5A8FA, 0x42B2986C, 0xDBBBC9D6, 0xACBCF940,
    0x32D86CE3, 0x45DF5C75, 0xDCD60DCF, 0xABD13D59,
    0x26D930AC, 0x51DE003A, 0xC8D75180, 0xBFD06116,
    0x21B4F4B5, 0x56B3C423, 0xCFBA9599, 0xB8BDA50F,
    0x2802B89E, 0x5F058808, 0xC60CD9B2, 0xB10BE924,
    0x2F6F7C87, 0x58684C11, 0xC1611DAB, 0xB6662D3D,

    0x76DC4190, 0x01DB7106, 0x98D220BC, 0xEFD5102A,
    0x71B18589, 0x06B6B51F, 0x9FBFE4A5, 0xE8B8D433,
    0x7807C9A2, 0x0F00F934, 0x9609A88E, 0xE10E9818,
    0x7F6A0DBB, 0x086D3D2D, 0x91646C97, 0xE6635C01,
    0x6B6B51F4, 0x1C6C6162, 0x856530D8, 0xF262004E,
    0x6C0695ED, 0x1B01A57B, 0x8208F4C1, 0xF50FC457,
    0x65B0D9C6, 0x12B7E950, 0x8BBEB8EA, 0xFCB9887C,
    0x62DD1DDF, 0x15DA2D49, 0x8CD37CF3, 0xFBD44C65,
    0x4DB26158, 0x3AB551CE, 0xA3BC0074, 0xD4BB30E2,
    0x4ADFA541, 0x3DD895D7, 0xA4D1C46D, 0xD3D6F4FB,
    0x4369E96A, 0x346ED9FC, 0xAD678846, 0xDA60B8D0,
    0x44042D73, 0x33031DE5, 0xAA0A4C5F, 0xDD0D7CC9,
    0x5005713C, 0x270241AA, 0xBE0B1010, 0xC90C2086,
    0x5768B525, 0x206F85B3, 0xB966D409, 0xCE61E49F,
    0x5EDEF90E, 0x29D9C998, 0xB0D09822, 0xC7D7A8B4,
    0x59B33D17, 0x2EB40D81, 0xB7BD5C3B, 0xC0BA6CAD,

    0xEDB88320, 0x9ABFB3B6, 0x03B6E20C, 0x74B1D29A,
    0xEAD54739, 0x9DD277AF, 0x04DB2615, 0x73DC1683,
    0xE3630B12, 0x94643B84, 0x0D6D6A3E, 0x7A6A5AA8,
    0xE40ECF0B, 0x9309FF9D, 0x0A00AE27, 0x7D079EB1,
    0xF00F9344, 0x8708A3D2, 0x1E01F268, 0x6906C2FE,
    0xF762575D, 0x806567CB, 0x196C3671, 0x6E6B06E7,
    0xFED41B76, 0x89D32BE0, 0x10DA7A5A, 0x67DD4ACC,
    0xF9B9DF6F, 0x8EBEEFF9, 0x17B7BE43, 0x60B08ED5,
    0xD6D6A3E8, 0xA1D1937E, 0x38D8C2C4, 0x4FDFF252,
    0xD1BB67F1, 0xA6BC5767, 0x3FB506DD, 0x48B2364B,
    0xD80D2BDA, 0xAF0A1B4C, 0x36034AF6, 0x41047A60,
    0xDF60EFC3, 0xA867DF55, 0x316E8EEF, 0x4669BE79,
    0xCB61B38C, 0xBC66831A, 0x256FD2A0, 0x5268E236,
    0xCC0C7795, 0xBB0B4703, 0x220216B9, 0x5505262F,
    0xC5BA3BBE, 0xB2BD0B28, 0x2BB45A92, 0x5CB36A04,
    0xC2D7FFA7, 0xB5D0CF31, 0x2CD99E8B, 0x5BDEAE1D,

    0x9B64C2B0, 0xEC63F226, 0x756AA39C, 0x026D930A,
    0x9C0906A9, 0xEB0E363F, 0x72076785, 0x05005713,
    0x95BF4A82, 0xE2B87A14, 0x7BB12BAE, 0x0CB61B38,
    0x92D28E9B, 0xE5D5BE0D, 0x7CDCEFB7, 0x0BDBDF21,
    0x86D3D2D4, 0xF1D4E242, 0x68DDB3F8, 0x1FDA836E,
    0x81BE16CD, 0xF6B9265B, 0x6FB077E1, 0x18B74777,
    0x88085AE6, 0xFF0F6A70, 0x66063BCA, 0x11010B5C,
    0x8F659EFF, 0xF862AE69, 0x616BFFD3, 0x166CCF45,
    0xA00AE278, 0xD70DD2EE, 0x4E048354, 0x3903B3C2,
    0xA7672661, 0xD06016F7, 0x4969474D, 0x3E6E77DB,
    0xAED16A4A, 0xD9D65ADC, 0x40DF0B66, 0x37D83BF0,
    0xA9BCAE53, 0xDEBB9EC5, 0x47B2CF7F, 0x30B5FFE9,
    0xBDBDF21C, 0xCABAC28A, 0x53B39330, 0x24B4A3A6,
    0xBAD03605, 0xCDD70693, 0x54DE5729, 0x23D967BF,
    0xB3667A2E, 0xC4614AB8, 0x5D681B02, 0x2A6F2B94,
    0xB40BBE37, 0xC30C8EA1, 0x5A05DF1B, 0x2D02EF8D,
};

//
// compute CRC32 checksum of given data
//
uint
crc32 ( const size_t bsize, const void * data )
{
    uint          crc = 0xffffffff;
    const uchar * str = static_cast< const uchar * >( data );

    for ( size_t i = 0; i < bsize; i++ )
        crc = ((crc) >> 8) ^ crc32_table[(str[i]) ^ ((crc) & 0x000000FF)];

    crc =~ crc;

    return crc;
}

//
// compute Adler32 checksum of given data
// (see RFC 1950)
//
ulong
adler32 ( const size_t size, const void * data )
{
    const ulong  BASE = 65521; // largest prime smaller than 65536
    const char * buf  = static_cast< const char * >( data );
    ulong        s1   = 1L & 0xffff;
    ulong        s2   = (1L >> 16) & 0xffff;

    for ( size_t i = 0; i < size; i++ )
    {
        s1 = (s1 + buf[i]) % BASE;
        s2 = (s2 + s1)     % BASE;
    }// for

    return (s2 << 16) + s1;
}

}// namespace Crypt

}// namespace Hpro
    
///////////////////////////////////////////////////
//
// Support for RTTI (independent of C++)
//
///////////////////////////////////////////////////

namespace Hpro
{

namespace RTTI
{

namespace
{
        
//
// database for id -> name conversion
//
struct id_name_t
{
    uint         id;
    std::string  name; 
};

//
// get instance of Type-DB
//
std::list< id_name_t > &
get_type_db ()
{
    static  std::list< id_name_t >  _type_db;

    return _type_db;
}

}// namespace anonymous
    
//
// convert a typename to a unique id
//
uint
type_to_id ( const std::string &  type )
{
    return Crypt::crc32<char>( type.size(), type.c_str() );
    // return Crypt::adler32( type.size(), type.c_str() );
}

//
// return typename to a given id
//
std::string
id_to_type ( const uint id )
{
    for ( std::list< id_name_t >::const_iterator  iter = get_type_db().begin();
          iter != get_type_db().end();
          ++iter )
    {
        if ( (*iter).id == id )
            return (*iter).name;
    }// for
            
    return "UNKOWN_ID";
}

//
// register type in RTTI system (returns id)
//
uint
register_type ( const char *  type )
{
    //
    // check in db if same ID with different name is present
    //

    const uint     id      = type_to_id( type );
    bool           present = false;
    static TMutex  mutex;

    mutex.lock();
    
    for ( auto  entry : get_type_db() )
    {
        if ( entry.id == id )
        {
            present = true;

            if ( std::string( type ) != entry.name )
                std::cerr << "type \"" << type << "\" and type \"" << entry.name << "\" have same id" << std::endl;
            
            break;
        }// if
    }// for

    //
    // insert into db if not present
    //

    if ( ! present )
    {
        id_name_t  entry;

        entry.id   = id;
        entry.name = type;
        
        get_type_db().push_back( entry );
    }// if

    mutex.unlock();
                
    return id;
}

//
// print all registered types in type db
//
void
print_registered ()
{
    //
    // determine maximal typename length and id length for pretty printing
    // also sort entries by id
    //

    uint                     name_width = 0;
    uint                     id_width   = 0;
    std::list< id_name_t >   sorted_db;
    
    for ( auto  entry : get_type_db() )
    {
        const auto  nlen = uint( entry.name.length() );
        const auto  nid  = uint( std::ceil( std::log10( entry.id ) ) );
        
        name_width = std::max( name_width, nlen );
        id_width   = std::max( id_width,   nid  );

        sorted_db.push_back( entry );
    }// for

    sorted_db.sort( [] ( const id_name_t &  e1, const id_name_t &  e2 ) { return e1.id < e2.id; } );

    //
    // print types/ids as a table
    //
    
    std::cout << std::setw( name_width ) << std::left << "Type Name"
              << " │ "
              << std::setw( id_width )   << std::left << "Type Id"
              << std::endl;

    for ( uint  i = 0; i < name_width; ++i ) std::cout << "─";
    std::cout << "─┼─";
    for ( uint  i = 0; i < id_width;   ++i ) std::cout << "─";
    std::cout << std::endl;
    
    for ( auto  entry : sorted_db )
    {
        std::cout << std::setw( name_width ) << std::left  << entry.name
                  << " │ "
                  << std::setw( id_width )   << std::right << entry.id
                  << std::endl;
    }// for
}

}// namespace RTTI

}// namespace Hpro
    
/////////////////////////////////////////////////////////////////
//
// Memory management
//
/////////////////////////////////////////////////////////////////

#include <hpro/base/alloc.hh>

#ifdef LINUX
#include <stdio.h>
#include <unistd.h>
#endif

#ifdef SUNOS
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#ifdef TRU64
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/procfs.h>
#endif

#ifdef AIX
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/procfs.h>
#endif

#ifdef HPUX
#include <unistd.h>
#include <sys/param.h>
#include <sys/pstat.h>
#endif

#ifdef DARWIN
#include <unistd.h>
#include <mach/task.h>
#include <mach/task_info.h>
#include <mach/mach_init.h>
#endif

#if 0
#include <sys/resource.h>
#endif

namespace Hpro
{

namespace Mem
{
//
// return byte usage
//
    
static size_t                                max_mem_usage = 0;
static std::unordered_map< void *, size_t >  mem_sizes;
size_t                                       mem_used = 0;
TMutex                                       mem_mutex;

//
// wrapper around malloc for statistics
//
void *
alloc ( const size_t  size )
{
    auto  alloc = allocator< uint8_t >();
    auto  p     = alloc.allocate( size );

    {
        TScopedLock  lock( mem_mutex );
        
        mem_sizes[ p ] = size;
        mem_used      += size;
    }

    return p;
}

//
// wrapper around free for statistics
//
void
free ( void *        ptr,
       const size_t  size )
{
    {
        TScopedLock  lock( mem_mutex );
        
        auto  s = mem_sizes[ ptr ];

        mem_used        -= s;
        mem_sizes[ ptr ] = 0;
    }

    auto  alloc = allocator< uint8_t >();
    
    alloc.deallocate( reinterpret_cast< uint8_t * >( ptr ), size );
}

//
// return byte usage
//
size_t
usage ()
{
    // uncomment for debugging malloc/free
    // return mem_used;
    
    size_t mem_size = 0;

#if defined(LINUX)

    //
    // try to read from "/proc/<pid>/smaps"
    //

    bool  success = false;
            
    {
        std::ostringstream  filename;
                
        filename << "/proc/" << getpid() << "/smaps";

        if ( boost::filesystem::exists( filename.str() ) )
        {
            std::ifstream  file( filename.str().c_str() );
                
            if ( ! file )
                HERROR( ERR_FOPEN, "(TMemory) used", filename.str().c_str() );

            std::string  buf;
            size_t       size = 0;
                    
            buf.reserve( 256 );
                    
            while ( file )
            {
                getline( file, buf );
                        
                if ( buf.substr( 0, 7 ) == "Private" )
                {
                    std::istringstream  iss( buf.substr( 16 ) );
                    size_t              lsize = 0;
                            
                    iss >> lsize;
                    size += lsize;
                }// if
            }// while
                    
            if ( size > 0 )
            {
                mem_size = size * 1024;
                success  = true;
            }// if
        }
    }
                
    //
    // then read size as first number in "/proc/<pid>/statm"
    //

    if ( ! success )
    {
        std::ostringstream  filename;
        std::ifstream       file;
                
        filename << "/proc/" << getpid() << "/statm";
                
        file.open( filename.str().c_str() );
                
        if ( ! file )
            HERROR( ERR_FOPEN, "(TMemory) used", filename.str().c_str() );
                
        size_t  size, rss, share, trs, drs, lrs, dt;
                
        file >> size >> rss >> share >> trs >> drs >> lrs >> dt;
                
        mem_size = size * size_t( getpagesize() );
    }
        
#elif defined(SUNOS)
#  if defined (__sparc_v9__) || defined (__sparcv9)
    
    //
    // read size of file "/proc/<pid>/as"
    //
    struct stat64  stat_buf;
    char           buf[256];
        
    snprintf( buf, sizeof(buf), "/proc/%d/as", uint(getpid()) );

    if ( stat64( buf, & stat_buf ) < 0 )
        HERROR( ERR_FOPEN, "(TMemory) used", std::string( "in stat64 to " ) + buf + " : " + syserror() );
        
    mem_size = stat_buf.st_size;

#  else

    //
    // read size of file "/proc/<pid>/as"
    //
    struct stat  stat_buf;
    char         buf[256];
        
    snprintf( buf, sizeof(buf), "/proc/%d/as", uint(getpid()) );
        
    if (stat( buf, & stat_buf ) < 0)
        HERROR( ERR_FOPEN, "(TMemory) used", std::string( "in stat to " ) + buf + " : " + syserror() );
        
    mem_size = stat_buf.st_size;

#  endif
#elif defined(TRU64)
    
    //
    // read size in "/proc/<pid>"
    //
    
    char              buf[256];
    int               file;
    struct prpsinfo * prpsinfo;
    struct prpsinfo   arrayprpsinfo[32];
    
    prpsinfo = arrayprpsinfo;
    snprintf( buf, sizeof(buf), "/proc/%d", (unsigned int) getpid() );
    
    if ((file = open( buf, O_RDONLY, 0555 )) < 0)
        HERROR( ERR_FOPEN, "(TMemory) used", buf );
    
    if ( ioctl( file, PIOCPSINFO, prpsinfo ) == -1 )
    {
        close( file );
        HERROR( ERR_PERM, "(TMemory) used", "in ioctl : " + syserror() );
    }// if
    
    close( file );
    mem_size = prpsinfo->pr_size * getpagesize();
    
#elif defined(AIX)

    //
    // read size from "/proc/<pid>/psinfo"
    // (returns only memory allocated via brk, not mmap)
    //
      
    char           buf[256];
    int            file;
    struct psinfo  psi;

    snprintf( buf, sizeof(buf), "/proc/%d/psinfo", (unsigned int) getpid() );

    if ((file = open( buf, O_RDONLY, 0555 )) < 0)
        HERROR( ERR_FOPEN, "(TMemory) used", "open" );
    
    if ( read( file, & psi, sizeof(psi) ) < 0 )
    {
        close( file );
        HERROR( ERR_FREAD, "(TMemory) used", "read" );
    }// if
    
    close( file );

    mem_size = psi.pr_size;

#elif defined(HPUX)

    //
    // use pstat functions to gather info
    //

    struct pst_status pst;

    if ( pstat_getproc( & pst, sizeof(pst), 0, getpid() ) == -1 )
        HERROR( ERR_PERM, "(TMemory) used", "in pstat_getproc : " + syserror() );

    mem_size  = pst.pst_dsize + pst.pst_tsize + pst.pst_ssize + pst.pst_mmsize;
    mem_size *= getpagesize();

#elif defined(DARWIN)

    struct task_basic_info   tinfo;
    task_t                   task   = mach_task_self();
    mach_msg_type_number_t   icount = TASK_BASIC_INFO_COUNT;

    if ( task_info( task, TASK_BASIC_INFO, (task_info_t) & tinfo, & icount ) != KERN_SUCCESS )
        HERROR( ERR_PERM, "(TMemory) used", "in task_info : " + syserror() );

    mem_size = tinfo.resident_size;
    // mem_size = tinfo.virtual_size;

#endif

    static TMutex  mutex;
            
    mutex.lock();
    max_mem_usage = std::max( max_mem_usage, mem_size );
    mutex.unlock();
            
    return mem_size;
}

//
// return maximal memory usage of application so far in bytes
//
size_t
max_usage ()
{
    return max_mem_usage;
}

//
// return string containing pretty printed memory usage
//
std::string
to_string ( const uint  base )
{
    return to_string( usage(), base );
}

//
// convert given number of bytes to human readable format
//
std::string
to_string ( const size_t byteval,
            const uint   base )
{
    std::ostringstream  str;
    size_t              bytes     = byteval;
    const char *        units10[] = { "B", "kB", "MB", "GB", "TB", "PB", "EB" };
    const char *        units2[]  = { "B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB" };
    const char **       units     = ( base == 10 ? units10 : units2 );
    const size_t        mult      = ( base == 10 ? 1000 : 1024 );
    
    if ( bytes < mult )
        str << bytes << ' ' << units[0];
    else
    {
        size_t  ofs = mult;
        uint    idx = 1;
            
        while (( bytes > ofs * mult ) && ( idx < 6 ))
        {
            ofs *= mult;
            ++idx;
        }// while
        
        const size_t  uval = bytes / (ofs);
        const size_t  lval = size_t( ((double(bytes) / double(ofs)) - double(uval)) * 100.0 );

        str << uval;
                
        if ( lval < 10 ) str << ".0";
        else             str << '.' ;

        str << lval << ' ' << units[idx];
    }// else

    return str.str();
}
        
}// namespace Memory

}// namespace Hpro

///////////////////////////////////////////////////
//
// machine properties
//
///////////////////////////////////////////////////

#ifdef LINUX
#  ifndef _GNU_SOURCE
#    define _GNU_SOURCE
#  endif
#  include <sched.h>
#endif

#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
#  include <Winsock2.h>
#endif

namespace Hpro
{

namespace Mach
{

//
// return string with list of CPUs associated to process
//
std::string
cpuset ()
{
    #if defined(LINUX)
    
    cpu_set_t  cset;

    CPU_ZERO( & cset );
    sched_getaffinity( 0, sizeof(cset), & cset );

    uint                ncores  = 0;
    std::ostringstream  out;
    int                 first   = -1;
    int                 last    = -1;
    bool                comma   = false;
    auto                prn_set = [&] ()
    {
        if ( comma ) out << ',';
                
        if      ( first == last    ) out << first;
        else if ( last  == first+1 ) out << first << ',' << last;
        else                         out << first << '-' << last;
    };

    for ( int  i = 0; i < 1024; ++i )
    {
        if ( CPU_ISSET( i, & cset ) )
        {
            ++ncores;
            
            // first initialization
            if ( first == -1 )
            {
                first = i;
                last  = i;
            }// if

            // new interval 
            if ( last < i-1 )
            {
                prn_set();
                
                first = i;
                last  = i;

                comma = true;
            }// if
            
            last = i;
        }// if
    }// for

    // finish expr
    prn_set();

    // add number of cores
    out << " (#" << ncores << ')';
                
    return out.str();

    #else

    return "";

    #endif
}

//
// return hostname
//
std::string
hostname ()
{
    char  buf[256];

    gethostname( buf, sizeof(buf) );

    return buf;
}

}// namespace Mach

}// namespace Hpro
