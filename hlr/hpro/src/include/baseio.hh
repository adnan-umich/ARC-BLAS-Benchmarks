#ifndef __HPRO_BASEIO_HH
#define __HPRO_BASEIO_HH
//
// Project     : HLIBpro
// File        : baseio.hh
// Description : basic IO related functions and classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <istream>
#include <ostream>
#include <string>
#include <cstdint>

#include "hpro/base/types.hh"
#include "hpro/base/traits.hh"

namespace Hpro
{

//////////////////////////////////////////////////////////////
//
// system independent types
//
//////////////////////////////////////////////////////////////

using  int8_t   = std::int8_t;
using  uint8_t  = std::uint8_t;
                
using  int16_t  = std::int16_t;
using  uint16_t = std::uint16_t;
                
using  int32_t  = std::int32_t;
using  uint32_t = std::uint32_t;
                
using  int64_t  = std::int64_t;
using  uint64_t = std::uint64_t;

///////////////////////////////////////////////////
// 
// byte swapping (little <-> big endian)
//

template <class T>
void
swap_bytes  ( T &  data );

///////////////////////////////////////////////////
// 
// open file with suitable stream, e.g. for
// compression
//

std::unique_ptr< std::istream >
open_read  ( const std::string &  filename );

std::unique_ptr< std::ostream >
open_write ( const std::string &  filename );

///////////////////////////////////////////////////
// 
// return extension without .gz etc.
//

std::string
extension ( const std::string &  filename );

///////////////////////////////////////////////////
// 
// file format detection
//

enum fmt_t {
        FMT_UNKNOWN,
        FMT_HPRO,
        FMT_HPRO_GRID,
        FMT_OCTAVE,
        FMT_SAMG,
        FMT_MATLAB,
        FMT_HDF5,
        FMT_HB,
        FMT_MTX,
        FMT_PLTMG,
        FMT_PLY,
        FMT_SURFMESH,
        FMT_GMSH
};

//
// guess file format of file \a filename
// - if \a filename exists, content of it
//   may be read
//
fmt_t
guess_format ( const std::string &  filename );


///////////////////////////////////////////////////
// 
// lexical cast with default value
//

long
str_to_int ( const std::string &  str,
             const long           def_val = 0 );

long
str_to_int ( const char *         begin,
             const char *         end,
             const long           def_val = 0 );

double
str_to_dbl ( const std::string &  str,
             const double         def_val = 0 );

double
str_to_dbl ( const char *         begin,
             const char *         end,
             const double         def_val = 0 );

///////////////////////////////////////////////////
//
// append file suffix if not already present
//

std::string
add_extension ( const std::string &  filename,
                const std::string &  suffix );

///////////////////////////////////////////////////
//
// compose any of float, double, std::complex< float >
// or std::complex< double > without known what to expect
//

template < typename value_t >
struct get_value
{
    static value_t compose ( real_type_t< value_t >  re,
                             real_type_t< value_t >  /* im */ )
    {
        return re;
    }
};

template <>
struct get_value< std::complex< float > >
{
    static std::complex< float > compose ( float  re, float  im )
    {
        return std::complex< float >( re, im );
    }
};

template <>
struct get_value< std::complex< double > >
{
    static std::complex< double > compose ( double  re, double  im )
    {
        return std::complex< double >( re, im );
    }
};

}// namespace Hpro

#endif  // __HPRO_BASEIO_HH
