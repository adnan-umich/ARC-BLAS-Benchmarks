#ifndef __HPRO_STRING_HH
#define __HPRO_STRING_HH
//
// Project     : HLIBpro
// File        : String.hh
// Description : module containing a class for a string
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <string>
#include <vector>

#include "hpro/base/config.hh"

namespace Hpro
{
    
/////////////////////////////////////////////////////////////////
//
// convert various arguments to string
//
/////////////////////////////////////////////////////////////////

std::string
to_string ( const char * style, ... );

template <class T>
std::string
to_string ( const T & obj )
{
    return obj.to_string();
}

#define BASETYPE_TOSTRING( type ) template <>           \
    std::string                                         \
    to_string< type > ( const type & n );

BASETYPE_TOSTRING( bool )
BASETYPE_TOSTRING( char )
BASETYPE_TOSTRING( unsigned char )
BASETYPE_TOSTRING( short )
BASETYPE_TOSTRING( unsigned short )
BASETYPE_TOSTRING( int )
BASETYPE_TOSTRING( unsigned int )
BASETYPE_TOSTRING( long )
BASETYPE_TOSTRING( unsigned long )
BASETYPE_TOSTRING( long long )
BASETYPE_TOSTRING( unsigned long long )
BASETYPE_TOSTRING( float )
BASETYPE_TOSTRING( double )
BASETYPE_TOSTRING( size_t )

template <> std::string  to_string< std::complex< float  > > ( const std::complex< float > &   f );
template <> std::string  to_string< std::complex< double > > ( const std::complex< double > &  f );

/////////////////////////////////////////////////////////////////
//
// functions for std::string
//
/////////////////////////////////////////////////////////////////

//! split \a str into sub strings seperated by characters in \a delim
//! and put results into \a parts
void
split     ( const std::string &           str,
            const std::string &           delim,
            std::vector< std::string > &  parts );

//! split \a str into sub strings seperated by characters in \a delim
std::vector< std::string >
split     ( const std::string &           str,
            const std::string &           delim );

}// namespace Hpro

#endif  // __HPRO_STRING_HH
