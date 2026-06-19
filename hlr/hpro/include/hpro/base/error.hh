#ifndef __HPRO_ERROR_HH
#define __HPRO_ERROR_HH
//
// Project     : HLIBpro
// File        : error.hh
// Description : error handling in HLIBpro
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <exception>
#include <string>

#include "hpro/base/types.hh"
#include "hpro/base/config.hh"

namespace Hpro
{

// forward decl. for REFERENCE below
template < typename value_t >
class TMatrix;

//////////////////////////////////////////////////////////////////
//
// error, warning, etc. reporting
//
//////////////////////////////////////////////////////////////////

enum loglevel_t
{
    LOG_ERROR    = 0,  // normal error, must abort operation
    LOG_WARNING  = 2,  // warning
    LOG_NOTICE   = 4,  // notice
    LOG_INFO     = 6,  // normal information
    LOG_DEBUG    = 8   // debug information
};

#define HERROR(   code, fnname, msg ) throw Hpro::Error( __FILE__, __LINE__, (fnname), (code), (msg) )
// #define HERROR(   code, fnname, msg ) Hpro::LOG::error( __FILE__, __LINE__, (msg) )
#define HWARNING( msg )               if ( Hpro::CFG::verbose( Hpro::LOG_WARNING ) ) Hpro::LOG::warning( __FILE__, __LINE__, (msg) )
#define HNOTICE(  msg )               if ( Hpro::CFG::verbose( Hpro::LOG_NOTICE  ) ) Hpro::LOG::note(    __FILE__, __LINE__, (msg) )
#define HINFO(    msg )               if ( Hpro::CFG::verbose( Hpro::LOG_INFO    ) ) Hpro::LOG::info(    __FILE__, __LINE__, (msg) )
#define HDEBUG(   msg )               if ( Hpro::CFG::verbose( Hpro::LOG_DEBUG   ) ) Hpro::LOG::debug(   __FILE__, __LINE__, (msg) )

#if HPRO_DEBUG == 0
#define HASSERT( cond, code, fnname, msg )
#else
#define HASSERT( cond, code, fnname, msg ) { if ( ! ( cond ) ) HERROR( code, fnname, msg ); }
#endif

//////////////////////////////////////////////////////////////////
//
// error exception storing filename, linenumber, function name,
// error code and error message
//
//////////////////////////////////////////////////////////////////

class Error : public std::exception
{
private:
    
    // position where error occured
    std::string  _file, _fnname;
    uint         _lineno;

    // error code and message
    uint         _errno;
    std::string  _errmsg;

    // full error string
    std::string  _error_string;
    
public:
    ///////////////////////////////////////////
    //
    // constructor
    //

    Error ();
    
    Error ( const std::string &  filename,
            const uint           lineno,
            const std::string &  fnname,
            const uint           err_no,
            const std::string &  errmsg );

    Error ( const std::string &  fnname,
            const uint           err_no,
            const std::string &  errmsg );

    Error ( const Error & e );
    
    virtual ~Error () throw ();
    
    ///////////////////////////////////////////
    //
    // give access to data
    //

    const std::string & file        () const { return _file;   }
    uint                line_number () const { return _lineno; }
    const std::string & function    () const { return _fnname; }
    uint                error_code  () const { return _errno;  }
    const std::string & error_msg   () const { return _errmsg; }

    // reset error
    virtual void reset ();
    
    // copy operator
    Error & operator = ( const Error & e );
    
    // print error message
    virtual void print () const;

    // convert error to string
    virtual std::string to_string () const;

    // standard description function for exceptions
    virtual const char * what() const throw();
};

//////////////////////////////////////////////////////////////////
//
// different levels of non-critical information
//
//////////////////////////////////////////////////////////////////

namespace LOG
{

// initialise and finish logging
void init ();
void done ();

// write message to log device
void print   ( const std::string &  msg );
void printf  ( const char *         fmt, ... );

// different logging levels
void error   ( const char *  filename, const uint  lineno, const std::string &  msg );
void warning ( const char *  filename, const uint  lineno, const std::string &  msg );
void note    ( const char *  filename, const uint  lineno, const std::string &  msg );
void info    ( const char *  filename, const uint  lineno, const std::string &  msg );
void debug   ( const char *  filename, const uint  lineno, const std::string &  msg );

//
// indentation for logging
//

// indentation block (change indentation for subblock)
struct indent_block
{
    // global counter
    static int  __INDENT_LEVEL;

    indent_block ()
    {
        __INDENT_LEVEL = std::max( 0, __INDENT_LEVEL + 2 );
    }

    ~indent_block ()
    {
        __INDENT_LEVEL = std::max( 0, __INDENT_LEVEL - 2 );
    }
};

// stream manipulator
struct indent
{
    template < class T_char, class T_traits >
    friend
    std::basic_ostream< T_char, T_traits > &
    operator <<  ( std::basic_ostream< T_char, T_traits > &  os,
                   const indent &                               )
    {
        for ( int  i = 0; i < indent_block::__INDENT_LEVEL; ++i )
            os.put( os.widen( ' ' ) );
        
        os.flush();

        return os;
    }
};

}// namespace

//////////////////////////////////////////////////////////////////
//
// C-style error handling
//
//////////////////////////////////////////////////////////////////

// return string for a given error code
std::string  strerror ( const uint ec );

// return corresponding error description for <ec> or 
// current system error defined in <errno>
std::string  syserror ( const uint ec );
std::string  syserror ();

//////////////////////////////////////////////////////////////////
//
// (long) list of error codes
//
//////////////////////////////////////////////////////////////////

// undefine NO_ERROR from <Windows.h>
#ifdef NO_ERROR
#  undef NO_ERROR
#endif

enum { NO_ERROR            = 0,         // no error occured

       ERR_INIT            = 100,       // not initialised
       ERR_LICENSE         = 101,       // invalid license
       ERR_NOT_IMPL        = 102,       // functionality not implemented
       ERR_CONSISTENCY     = 103,       // general consistency error
       ERR_COMM            = 104,       // communication error
       ERR_PERM            = 105,       // permission denied
       ERR_EXCEPTION       = 106,       // general exception error

       ERR_REAL            = 200,       // data is real valued
       ERR_NREAL           = 201,       // data is not real valued
       ERR_COMPLEX         = 202,       // data is complex valued
       ERR_NCOMPLEX        = 203,       // data is not complex valued
       ERR_REAL_CMPLX      = 204,       // invalid mixing of real and complex data
       ERR_DIV_ZERO        = 205,       // division by zero
       ERR_NEG_SQRT        = 206,       // sqrt of negative number
       ERR_INF             = 207,       // infinity occured
       ERR_NAN             = 208,       // not-a-number occured
       ERR_NCONVERGED      = 209,       // iteration did not converge
       
       ERR_ARG             = 300,       // error with argument
       ERR_MEM             = 301,       // insufficient memory available
       ERR_NULL            = 302,       // null pointer encountered
       ERR_SIZE            = 303,       // size of data incorrect
       ERR_INDEXSET        = 304,       // wrong index set
       ERR_DIM             = 305,       // invalid or incompatible dimension
       ERR_ARR_BOUND       = 306,       // out-of-bound error in array
       ERR_DIAG_ENTRY      = 307,       // entry is not on diagonal

       ERR_COORD_INVALID   = 400,       // invalid coordinates
       
       ERR_CT_INVALID      = 500,       // invalid cluster tree
       ERR_CT_TYPE         = 501,       // wrong type of cluster tree
       ERR_CT_STRUCT       = 502,       // invalid structure of cluster tree
       ERR_CT_INCOMP       = 503,       // given cluster trees are incompatible
       ERR_CT_SPARSE       = 504,       // missing sparse matrix for given cluster tree
       ERR_CT_DEPTH        = 505,       // depth of cluster tree too large
       
       ERR_BCT_INVALID     = 600,       // invalid block cluster tree
       ERR_BCT_STRUCT      = 601,       // invalid block cluster tree structure
       
       ERR_VEC_INVALID     = 700,       // invalid vector
       ERR_VEC_TYPE        = 701,       // wrong vector type
       ERR_VEC_STRUCT      = 702,       // invalid vector structure
       ERR_VEC_SIZE        = 703,       // invalid size of vector
       ERR_VEC_INCOMP      = 704,       // vector with incompatible dimension
       ERR_VEC_NSCALAR     = 705,       // vector is not a scalar vector
       
       ERR_MAT_TYPE        = 800,       // invalid matrix type
       ERR_MAT_STRUCT      = 801,       // invalid structure of matrix
       ERR_MAT_SIZE        = 802,       // invalid size of matrix
       ERR_MAT_SINGULAR    = 803,       // singular matrix detected
       ERR_MAT_NSPARSE     = 804,       // matrix not a sparse matrix
       ERR_MAT_NDENSE      = 805,       // matrix not a dense matrix
       ERR_MAT_NHMAT       = 806,       // matrix not an H-matrix
       ERR_MAT_INCOMP_TYPE = 807,       // matrices with incompatible type
       ERR_MAT_INCOMP_CT   = 808,       // matrices with incompatible cluster tree
       ERR_MAT_INVALID     = 809,       // invalid matrix
       ERR_MAT_NSYM        = 810,       // matrix not symmetric
       ERR_MAT_NHERM       = 811,       // matrix not hermitian
       ERR_MAT_NPOSDEF     = 812,       // matrix not positive definite

       ERR_FMT_UNKNOWN     = 900,       // unknown file format detected
       ERR_FMT_HFORMAT     = 901,       // error while parsing HLIBpro format
       ERR_FMT_SAMG        = 902,       // error while parsing SAMG format
       ERR_FMT_MATLAB      = 903,       // error while parsing Matlab format
       ERR_FMT_PLTMG       = 904,       // error while parsing PLTMG format
       ERR_FMT_HB          = 905,       // error while parsing Harwell Boeing format
       ERR_FMT_MTX         = 906,       // error while parsing Matrix Market format
       ERR_FMT_PLY         = 907,       // error while parsing Ply format
       ERR_FMT_CST         = 908,       // error while parsing CST format
      
       ERR_GRID_FORMAT     = 1000,      // invalid format of grid file
       ERR_GRID_DATA       = 1001,      // invalid data in grid file

       ERR_FOPEN           = 1100,      // could not open file
       ERR_FCLOSE          = 1101,      // could not close file
       ERR_FWRITE          = 1102,      // could not write to file
       ERR_FREAD           = 1103,      // could not read from file
       ERR_FSEEK           = 1104,      // could not seek in file
       ERR_FNEXISTS        = 1105,      // file does not exists

       ERR_BS_SIZE         = 1200,      // size of bytestream too small
       ERR_BS_WRITE        = 1201,      // error while writing to bytestream
       ERR_BS_READ         = 1202,      // error while reading from bytestream
       ERR_BS_TYPE         = 1203,      // type error in bytestream
       ERR_BS_DATA         = 1204,      // general data error in bytestream
       
       ERR_NOZLIB          = 1300,      // no zlib support available
       ERR_ZLIB_UNZIP      = 1301,      // error during zlib uncompression
       ERR_NOMETIS         = 1302,      // no METIS support available
       ERR_NOSCOTCH        = 1303,      // no Scotch support available
       ERR_SCOTCH          = 1304,      // error in call to Scotch function
       ERR_NOCHACO         = 1305,      // no Chaco support available
       ERR_NOLIBGRAPH      = 1306,      // no libGraph support available
       ERR_NOFFTW3         = 1307,      // no FFTW3 support available
       ERR_NOCAIRO         = 1308,      // no Cairo support available
       ERR_NOHDF5          = 1309,      // no HDF5 support available
       ERR_NOMONGOOSE      = 1310,      // no Mongoose support available
       ERR_NONETCDF        = 1311,      // no NetCDF support available
       ERR_NETCDF          = 1312,      // error during NetCDF I/O
       ERR_NOCGAL          = 1313,      // no CGAL support available
       ERR_CUDA            = 1314,      // error within CUDA function
       ERR_HIP             = 1315,      // error within HIP function
       
       ERR_MPI             = 1400,      // error in call to MPI function
       
       ERR_SOLVER_INVALID  = 1500,      // invalid solver
       ERR_LRAPX_INVALID   = 1501,      // invalid low-rank approximation type
       ERR_GRID_INVALID    = 1502,      // invalid grid
       ERR_FNSPACE_INVALID = 1503,      // invalid function space

       ERR_THR_CREATE      = 1600,      // error while creating thread
       ERR_THR_JOIN        = 1601,      // error while joining with thread
       ERR_THR_DETACH      = 1602,      // error while detaching thread
       ERR_THR_CANCEL      = 1603,      // error during thread cancelation

       ERR_PAR_PROCSET     = 1700,      // invalid or wrong processor set
       
       //
       // just for convenience
       //

       ERR_LAST_ERROR
};

//////////////////////////////////////////////////////////////////
//
// debugging
//
//////////////////////////////////////////////////////////////////

namespace DBG
{

// debug messages
void print      ( const std::string &  msg );
void printf     ( const char *         fmt, ... );

// add/remove indent in output
void indent     ( const uint  ofs = 1 );

// print function backtrace up to n entries (n=0 : all available entries)
void backtrace  ( const int  n    = 0,
                  const int  skip = 1  );

//
// breakpoint for debugging
//
void
breakpoint ();

}// namespace DBG

}// namespace Hpro

#endif  // __HPRO_ERROR_HH
