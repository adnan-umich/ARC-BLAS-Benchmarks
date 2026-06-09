//
// Project     : HLIBpro
// File        : String.cc
// Description : module containing a class for a string
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <stdarg.h>

#include <boost/algorithm/string.hpp>

#include "hpro/base/String.hh"

namespace Hpro
{

//
// convert various arguments to string
//
std::string
to_string ( const char * style, ... )
{
    if ( style == nullptr )
        return "";

    va_list ap;
    char    buffer[1024];

    va_start( ap, style );

#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
    vsnprintf_s( buffer, sizeof(buffer), _TRUNCATE, style, ap );
#else
    vsnprintf( buffer, sizeof(buffer), style, ap );
#endif
    
    va_end(ap);
    
    return buffer;
}

/////////////////////////////////////////////////////////////////
//
// functions for std::string
//
/////////////////////////////////////////////////////////////////

#define BASETYPE_TOSTRING_IMPL( type ) template <>      \
    std::string                                         \
    to_string< type > ( const type & n )                \
    {                                                   \
        return std::to_string( n );                     \
    }

BASETYPE_TOSTRING_IMPL( bool )
BASETYPE_TOSTRING_IMPL( char )
BASETYPE_TOSTRING_IMPL( unsigned char )
BASETYPE_TOSTRING_IMPL( short )
BASETYPE_TOSTRING_IMPL( unsigned short )
BASETYPE_TOSTRING_IMPL( int )
BASETYPE_TOSTRING_IMPL( unsigned int )
BASETYPE_TOSTRING_IMPL( long )
BASETYPE_TOSTRING_IMPL( unsigned long )
BASETYPE_TOSTRING_IMPL( long long )
BASETYPE_TOSTRING_IMPL( unsigned long long )
BASETYPE_TOSTRING_IMPL( float )
BASETYPE_TOSTRING_IMPL( double )

template <>
std::string
to_string< std::complex< float > > ( const std::complex< float > &  f )
{
    std::ostringstream  out;
    const char *        imag_suffix = ( CFG::IO::use_matlab_syntax ? "*i" : "i" );
    const float         im          = f.imag();
    
    if ( im != float(0) )
    {
        if ( im < float(0) )
            out << f.real() << " - " << std::abs( f.imag() ) << imag_suffix;
        else
            out << f.real() << " + " << f.imag() << imag_suffix;
    }// if
    else
        out << f.real();

    return out.str();
}

template <>
std::string
to_string< std::complex< double > > ( const std::complex< double > &  f )
{
    std::ostringstream  out;
    const char *        imag_suffix = ( CFG::IO::use_matlab_syntax ? "*i" : "i" );
    const double        im          = f.imag();
    
    if ( im != double(0) )
    {
        if ( im < double(0) )
            out << f.real() << " - " << std::abs( f.imag() ) << imag_suffix;
        else
            out << f.real() << " + " << f.imag() << imag_suffix;
    }// if
    else
        out << f.real();

    return out.str();
}

//
// split \a str into sub strings seperated by characters in \a delim
// and put results into \a parts
//
void
split  ( const std::string &           str,
         const std::string &           delim,
         std::vector< std::string > &  parts )
{
    boost::split( parts, str, boost::is_any_of( delim ), boost::token_compress_on );
}

std::vector< std::string >
split  ( const std::string &           str,
         const std::string &           delim )
{
    std::vector< std::string >   parts;

    split( str, delim, parts );

    return  parts;
}

//
// convert to lowercase characters
//
void
to_lower ( std::string &  str )
{
    boost::to_lower( str );
}

//
// convert to lowercase characters
//
void
to_upper ( std::string &  str )
{
    boost::to_upper( str );
}

}// namespace Hpro
