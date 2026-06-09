//
// Project     : HLIBpro
// File        : error.cc
// Description : error handling in HLIBpro
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/config.h"

#include <iostream>
#include <cstdlib>
#include <cstring>

#include <stdarg.h>
#include <errno.h>
#include <time.h>

#if HPRO_HAS_BACKTRACE == 1
#include <execinfo.h>
#include <dlfcn.h>
#endif

#if HPRO_HAS_CXXDEMANGLE == 1
#include <cxxabi.h>
#endif

#include <fstream>
#include <vector>

#include "hpro/base/System.hh"
#include "hpro/parallel/TMutex.hh"

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TBlockMatrix.hh"
#include "hpro/matrix/structure.hh"

#include "hpro/base/error.hh"

//
// replace "localtime_r" by "localtime_s" on Windows systems
//
#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
#define localtime_r( time, timeval )  localtime_s( timeval, time )
#endif

namespace Hpro
{

namespace DBG
{

// see definition below
void exception ( const std::string &  msg );

}// namespace DBG

namespace
{

#if HPRO_HAS_BACKTRACE == 1

std::string
demangle ( char *      symbol,
           const bool  only_funcname = false )
{
    #if HPRO_HAS_CXXDEMANGLE == 1
    
    char    no_output    = '\0';
    char *  begin_name   = nullptr;
    char *  begin_offset = nullptr;
    char *  end_offset   = nullptr;
    char *  begin_addr   = nullptr;
    char *  end_addr     = nullptr;

    // find parentheses and +address offset surrounding the mangled name:
    // ./module(function+0x15c) [0x8048a6d]
    for ( char * p = symbol; *p; ++p )
    {
        if ( *p == '(' )
        {
            begin_name = p;
        }// if
        else if ( *p == '+' )
        {
            begin_offset = p;
        }// if
        else if ( *p == ')' && ( begin_offset != nullptr ))
        {
            end_offset = p;
        }// if
        else if ( *p == '[' )
        {
            begin_addr = p;
        }// if
        else if ( *p == ']' )
        {
            end_addr = p;
        }// if
    }// for

    if (( begin_name   != nullptr ) &&
        ( begin_offset != nullptr ) &&
        ( end_offset   != nullptr ) &&
        ( begin_name   <  begin_offset ))
    {
        size_t               funcnamesize = 256;
        std::vector< char >  funcname( funcnamesize );
        
        *begin_name++   = '\0';
        *begin_offset++ = '\0';
        *end_offset     = '\0';

        // mangled name is now in [begin_name, begin_offset) and caller
        // offset in [begin_offset, end_offset). now apply
        // __cxa_demangle():

        int     status;
        char *  ret = abi::__cxa_demangle( begin_name, funcname.data(), & funcnamesize, & status );

        if ( status == 0 )
        {
            // funcname = ret; // use possibly realloc()-ed string

            if ( only_funcname )
                return ret;
            
            if (( begin_addr != nullptr ) && ( end_addr != nullptr ))
            {
                end_addr = & no_output;
                return to_string( "%s : %s %s", symbol, ret, begin_addr );
            }// if
            else
                return to_string( "%s : %s+%s", symbol, ret, begin_offset );
        }// if
        else
        {
            // demangling failed. Output function name as a C function with
            // no arguments.
            return to_string( "%s : %s()+%s", symbol, begin_name, begin_offset );
        }// else
    }// if
    else
    {
        // if function name requested but none available: return empty string
        if ( only_funcname )
            return "";
        
        // couldn't parse the line : print the whole line.
        return to_string( "%s", symbol );
    }// else

    #else
    
    return std::string( symbol );
    
    #endif
}

#endif

//
// inquire environment variable and return true if 
// available or <false> otherwise
//
bool
get_environ ( const std::string &  envname,
              std::string &        content )
{
    char  * val = nullptr;
    
#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)

    size_t  val_len  = 0;

    if ( ! _dupenv_s( & val, & val_len, envname.c_str() ) && ( val != nullptr ))
    {
        content = val;
        free( val );
    }// if
#else
    if (( val = ::getenv( envname.c_str() ) ) != nullptr )
        content = val;
#endif 
    else
    {
        content = "";
        return false;
    }// else

    return true;
}

}// namespace

//////////////////////////////////////////////////////////////////
//
// error exception storing filename, linenumber, function name,
// error code and error message
//
//////////////////////////////////////////////////////////////////

//
// constructor and destructor
//
Error::Error ()
{
    _errno  = NO_ERROR;
    _lineno = 0;

    _error_string = to_string();
}
               
Error::Error ( const std::string & filename,
               const uint          lineno,
               const std::string & fnname,
               const uint          err_no,
               const std::string & errmsg )
{
    _file   = filename;
    _lineno = lineno;
    _fnname = fnname;
    _errno  = err_no;
    _errmsg = errmsg;

    _error_string = to_string();

    // if ( verbose(2) && ( _errno != NO_ERROR ))
    //     print();
    
    // error breakpoint for debugging
    if ( _errno != NO_ERROR )
        DBG::exception( _error_string );
}

Error::Error ( const std::string & fnname,
               const uint          err_no,
               const std::string & errmsg )
{
    _file   = "";
    _lineno = 0;
    _fnname = fnname;
    _errno  = err_no;
    _errmsg = errmsg;

    _error_string = to_string();
    
    // if ( verbose(2) && ( _errno != NO_ERROR ))
    //     print();
    
    // error breakpoint for debugging
    if ( _errno != NO_ERROR )
        DBG::exception( _error_string );
}

Error::Error ( const Error & e )
        : std::exception()
{
    *this = e;
}

Error::~Error () throw ()
{
}
    
//
// reset error
//
void
Error::reset ()
{
    _file   = "";
    _lineno = 0;
    _fnname = "";
    _errno  = NO_ERROR;
    _errmsg = "";
}
    
//
// copy operator
//
Error &
Error::operator = ( const Error & e )
{
    _file   = e._file;
    _lineno = e._lineno;
    _fnname = e._fnname;
    _errno  = e._errno;
    _errmsg = e._errmsg;

    _error_string = e._error_string;
    
    return *this;
}

//
// print error message
//
void
Error::print () const
{
    std::cerr << to_string() << std::endl;

    if ( verbose(5) )
    {
        std::cout << "backtrace:" << std::endl;
        
        DBG::backtrace();
    }// if
}

//
// convert error to string
//
std::string
Error::to_string () const
{
#if 0
    
    std::string  msg;

    msg  = "ERROR " + Hpro::to_string( "%d", _errno );

    if ( _fnname != "" )
        msg += " in function \"" + _fnname + "\"";

    msg += '\n';
    msg += "    -> " + strerror( _errno );
    
    if ( _errmsg != "" )
        msg += " (" + _errmsg + ")";
    
    msg += '\n';
    msg += "       (in file \"" + _file + "\" at line " + Hpro::to_string( "%d", _lineno ) + ')';
    
#else
    
    std::string  msg;

    if ( _fnname != "" )
    {
        msg = " in \"" + _fnname + "\" at \"" + _file + Hpro::to_string( ":%d", _lineno ) + "\"\n";
        msg += "    Error: " + strerror( _errno );
    
        if ( _errmsg != "" )
            msg += " (" + _errmsg + ")";
    }// if
    else
    {
        msg = " at \"" + _file + Hpro::to_string( ":%d", _lineno ) + "\"\n";
        msg += "    Error: " + strerror( _errno );
    
        if ( _errmsg != "" )
            msg += " (" + _errmsg + ")";
    }// else
    
    msg += '\n';
    
#endif
    
    return msg;
}

//
// standard description function for exceptions
//
const char *
Error::what() const throw()
{
    return _error_string.c_str();
    // return "HLIBpro error, see error code";
}

//////////////////////////////////////////////////////////////////
//
// different level of information
//
//////////////////////////////////////////////////////////////////

namespace LOG
{

// indentation level
int  indent_block::__INDENT_LEVEL = 0;

namespace
{

//
// logfile for all LOG output
// (default: std::cout)
//
std::ostream * logfile = nullptr;

//
// general logging function
//
void
prnlog ( const std::string &  prefix,
         const std::string &  amsg)
{
    //
    // guard logging by mutex
    //
    
    static TMutex  mutex;

    TScopedLock    lock( mutex );

    //
    // get current date
    //
    
    time_t  t            = time( nullptr );
    char    timestr[256] = "\0";

#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64) || HPRO_HAS_LOCALTIME_R == 1
    struct tm   timeval;

    localtime_r( & t, & timeval );
    strftime( timestr, sizeof(timestr), "%Y-%m-%d %H:%M:%S", & timeval );
#else
    struct tm * timeval;
    
    timeval = localtime( & t );
    strftime( timestr, sizeof(timestr), "%F %T", timeval );
#endif

    //
    // compose final message
    //

    if ( prefix == "" )
    {
        if ( logfile != nullptr )
            (*logfile) << '[' << timestr << "] " << amsg << std::endl;
        else
            std::cout << '[' << timestr << "] " << amsg << std::endl;
    }// if
    else
    {
        if ( logfile != nullptr )
            (*logfile) << '[' << timestr << "] " << prefix << " : " << amsg << std::endl;
        else
            std::cout << '[' << timestr << "] " << prefix << " : " << amsg << std::endl;
    }// else
}

}// namespace anonymous

//
// intialise logging
//
void
init ()
{
    //
    // read environment variable for nthreads
    //

    std::string  env_logfile;

    if ( get_environ( "HPRO_LOGFILE", env_logfile ) )
        logfile = new std::ofstream( env_logfile.c_str() );
}

//
// finish logging
//
void
done ()
{
    if ( logfile != nullptr )
        delete logfile;
}

//
// write message to log device
//
void
print ( const std::string &  msg )
{
    LOG::prnlog( "", msg );
}

void
printf ( const char * fmt, ... )
{
    if ( fmt == nullptr )
        LOG::print( "" );

    va_list ap;
    char    buffer[1024];

    va_start( ap, fmt );

#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
    vsnprintf_s( buffer, sizeof(buffer), _TRUNCATE, fmt, ap );
#else
    vsnprintf( buffer, sizeof(buffer), fmt, ap );
#endif
    
    va_end(ap);
    
    LOG::print( buffer );
}

void
error ( const char *        filename,
        const uint          lineno,
        const std::string & msg )
{
    LOG::prnlog( "[error]", to_string( "%s(%d) ", filename, lineno ) + msg );
}
    
void
warning ( const char *        filename,
          const uint          lineno,
          const std::string & msg )
{
    LOG::prnlog( "[warn]", to_string( "%s(%d) ", filename, lineno ) + msg );
}
    
void
note ( const char *        filename,
       const uint          lineno,
       const std::string & msg )
{
    LOG::prnlog( "[notice]", to_string( "%s(%d) ", filename, lineno ) + msg );
}

void
info ( const char *        filename,
       const uint          lineno,
       const std::string & msg )
{
    LOG::prnlog( "[info]", to_string( "%s(%d) ", filename, lineno ) + msg );
}

void
debug ( const char *        filename,
        const uint          lineno,
        const std::string & msg )
{
    LOG::prnlog( "[debug]", to_string( "%s(%d) ", filename, lineno ) + msg );
}

}// namespace LOG

///////////////////////////////////////////////////////
//
// C-style error reporting
//
///////////////////////////////////////////////////////

void
report_error ( const char *       file,
               const unsigned int lineno,
               const char *       fnname,
               const unsigned int code,
               const char *       msg )
{
    throw Error( file, lineno, fnname, code, msg );
}
    
//
// return string for a given error code
//
std::string
strerror ( const uint ec )
{
    switch ( ec )
    {
        case NO_ERROR            : return "no error";
        
        case ERR_INIT            : return "not initialised";
        case ERR_LICENSE         : return "invalid license";
        case ERR_NOT_IMPL        : return "functionality not implemented";
        case ERR_CONSISTENCY     : return "detected inconsistency";
        case ERR_COMM            : return "error during communication";
        case ERR_PERM            : return "permission denied";

        case ERR_REAL            : return "data is real valued";
        case ERR_NREAL           : return "data is not real valued";
        case ERR_COMPLEX         : return "data is complex valued";
        case ERR_NCOMPLEX        : return "data is not complex valued";
        case ERR_REAL_CMPLX      : return "invalid mixing of real and complex data"; 
        case ERR_DIV_ZERO        : return "division by zero";
        case ERR_NEG_SQRT        : return "square root of negative number";
        case ERR_INF             : return "<INFINITY> occured as a number-value";
        case ERR_NAN             : return "<NOT-A-NUMBER> occured as a number-value";
        case ERR_NCONVERGED      : return "iteration did not converge";

        case ERR_ARG             : return "invalid argument";
        case ERR_MEM             : return "insufficient memory available";
        case ERR_NULL            : return "unexpected NULL pointer";
        case ERR_SIZE            : return "incorrect size of data";
        case ERR_INDEXSET        : return "incorrect index set";
        case ERR_DIM             : return "invalid or incompatible dimension";
        case ERR_ARR_BOUND       : return "out-of-bound error in array";
        case ERR_DIAG_ENTRY      : return "entry is not on diagonal";

        case ERR_COORD_INVALID   : return "invalid coordinates";
        
        case ERR_CT_INVALID      : return "invalid cluster tree";
        case ERR_CT_TYPE         : return "invalid type of cluster tree";
        case ERR_CT_STRUCT       : return "invalid structure of cluster tree";
        case ERR_CT_INCOMP       : return "given cluster trees are incompatible";
        case ERR_CT_SPARSE       : return "missing sparse matrix for given cluster tree";
        case ERR_CT_DEPTH        : return "depth of cluster tree too large";
        
        case ERR_BCT_INVALID     : return "invalid block cluster tree";
        case ERR_BCT_STRUCT      : return "invalid block cluster tree structure";
        
        case ERR_VEC_INVALID     : return "invalid vector";
        case ERR_VEC_TYPE        : return "wrong vector type";
        case ERR_VEC_STRUCT      : return "vector with invalid structure";
        case ERR_VEC_SIZE        : return "invalid size of vector";
        case ERR_VEC_INCOMP      : return "vector with incompatible dimension";
        case ERR_VEC_NSCALAR     : return "vector is not a scalar vector";
        
        case ERR_MAT_TYPE        : return "invalid matrix type";
        case ERR_MAT_STRUCT      : return "invalid (block-) matrix structure";
        case ERR_MAT_SIZE        : return "invalid matrix size";
        case ERR_MAT_SINGULAR    : return "singular matrix detected";
        case ERR_MAT_NSPARSE     : return "matrix not a sparse matrix";
        case ERR_MAT_NDENSE      : return "matrix not a dense matrix";
        case ERR_MAT_NHMAT       : return "matrix not an H-matrix";
        case ERR_MAT_INCOMP_TYPE : return "matrices with incompatible type";
        case ERR_MAT_INCOMP_CT   : return "matrices with incompatible cluster tree";
        case ERR_MAT_INVALID     : return "invalid matrix";
        case ERR_MAT_NSYM        : return "matrix not symmetric";
        case ERR_MAT_NHERM       : return "matrix not hermitian";
        case ERR_MAT_NPOSDEF     : return "matrix not positive definite";
        
        case ERR_FMT_UNKNOWN     : return "detected unknown file format";
        case ERR_FMT_HFORMAT     : return "error while parsing HLIBpro format";
        case ERR_FMT_SAMG        : return "error while parsing SAMG format";
        case ERR_FMT_MATLAB      : return "error while parsing Matlab format";
        case ERR_FMT_PLTMG       : return "error while parsing PLTMG format";
        case ERR_FMT_HB          : return "error while parsing Harwell-Boeing format";
        case ERR_FMT_CST         : return "error while parsing CST format";
        
        case ERR_GRID_FORMAT     : return "invalid format of grid file";
        case ERR_GRID_DATA       : return "invalid data in grid file";

        case ERR_FOPEN           : return "could not open file";
        case ERR_FCLOSE          : return "could not close file";
        case ERR_FWRITE          : return "could not write to file";
        case ERR_FREAD           : return "could not read from file";
        case ERR_FSEEK           : return "could not seek in file";
        case ERR_FNEXISTS        : return "file does not exist";
        
        case ERR_BS_SIZE         : return "size of bytestream too small";
        case ERR_BS_WRITE        : return "error while writing to bytestream";
        case ERR_BS_READ         : return "error while reading from bytestream";
        case ERR_BS_TYPE         : return "wrong type in bytestream";
        case ERR_BS_DATA         : return "wrong or unexpected data in bytestream";
        
        case ERR_NOZLIB          : return "no zlib support compiled in";
        case ERR_ZLIB_UNZIP      : return "error during zlib uncompression";
        case ERR_NOMETIS         : return "no METIS support compiled in";
        case ERR_NOSCOTCH        : return "no Scotch support compiled in";
        case ERR_SCOTCH          : return "error in call to Scotch function";
        case ERR_NOCHACO         : return "no Chaco support compiled in";
        case ERR_NOLIBGRAPH      : return "no libGraph support compiled in";
        case ERR_NOFFTW3         : return "no FFTW3 support compiled in";
        case ERR_NOCAIRO         : return "no Cairo support compiled in";
        case ERR_NOHDF5          : return "no HDF5 support compiled in";
        case ERR_NOMONGOOSE      : return "no Mongoose support compiled in";
        case ERR_NONETCDF        : return "no NetCDF support compiled in";
        case ERR_NETCDF          : return "error during NetCDF I/O";
        case ERR_NOCGAL          : return "no CGAL support compiled in";
        case ERR_CUDA            : return "error within CUDA function";
        case ERR_HIP             : return "error within HIP function";

        case ERR_MPI             : return "error in call to MPI function";

        case ERR_SOLVER_INVALID  : return "invalid solver";
        case ERR_LRAPX_INVALID   : return "invalid low-rank approximation type";
        case ERR_GRID_INVALID    : return "invalid grid";
        case ERR_FNSPACE_INVALID : return "invalid function space";

        case ERR_THR_CREATE      : return "error while creating thread";
        case ERR_THR_JOIN        : return "error while joining with thread";
        case ERR_THR_DETACH      : return "error while detaching thread";
        case ERR_THR_CANCEL      : return "error during thread cancelation";
        
    }// switch

    return "unknown error";
}

//
// return corresponding error description for current
// system error defined in <errno>
//
std::string
syserror ( const int errcode )
{
    char  buffer[1024];

#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
    strerror_s( buffer, sizeof(buffer), errcode );
#else
#  if HPRO_HAS_STRERROR_R == 1
    strerror_r( errcode, buffer, sizeof(buffer) );
#  else
    return ::strerror( errcode );
#  endif
#endif
    
    return buffer;
}

std::string
syserror ()
{
    return syserror( errno );
}


namespace DBG
{

uint  indent_ofs = 0;
    
//
// debug hook for exceptions
//
void
exception ( const std::string &  /* msg */ )
{
    // if ( verbose( LOG_ERROR ) )
    //     std::cerr << msg << std::endl;
}

//
// write message to log device
//
void
print ( const std::string &  msg )
{
    if ( indent_ofs > 0 )
        for ( uint  i = 0; i < indent_ofs; ++i )
            std::cout << "  ";
    
    std::cout << msg << std::endl;
}

void
printf ( const char * fmt, ... )
{
    if ( fmt == nullptr )
        LOG::print( "" );

    va_list ap;
    char    buffer[1024];

    va_start( ap, fmt );

#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
    vsnprintf_s( buffer, sizeof(buffer), _TRUNCATE, fmt, ap );
#else
    vsnprintf( buffer, sizeof(buffer), fmt, ap );
#endif
    
    va_end(ap);
    
    DBG::print( buffer );
}

void
indent ( const uint  ofs )
{
    if ( ofs > 0 )
        indent_ofs += ofs;
    else
    {
        if ( indent_ofs >= ofs )
            indent_ofs -= ofs;
        else
            indent_ofs =  0;
    }// else
}

//
// print function backtrace up to n entries (n=0 : all available entries)
//
void
backtrace  ( const int  max_entries,
             const int  skip )
{
    #if HPRO_HAS_BACKTRACE == 1
    
    void *   callstack[1024];
    int      nsym    = ::backtrace( callstack, sizeof(callstack) / sizeof(callstack[0]) );
    char **  symbols = backtrace_symbols( callstack, nsym );

    if ( symbols != nullptr )
    {
        int  nentries = 0;
        
        for ( int  i = skip; i < nsym; ++i )
        {
            const auto  res = demangle( symbols[i], true );

            if ( res.size() > 0 )
            {
                std::cout << res << std::endl;
                nentries++;
            }// if

            if ( nentries > max_entries )
                break;
        }// for
        
        free( symbols );
    }// if
    
    #endif
}

//
// breakpoint for debugging
//
void
breakpoint ()
{
}

}// namespace DBG

}// namespace Hpro

