//
// Project     : HLIBpro
// File        : TByteStream.hh
// Description : class for a byte stream
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <iostream>
#include <cstring>
#include <memory>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "hpro/base/error.hh"
#include "hpro/base/System.hh"
#include "hpro/base/TByteStream.hh"

#include "baseio.hh"

namespace Hpro
{

using std::unique_ptr;

////////////////////////////////////////////////////
//
// constructor and destructor
//

TByteStream::TByteStream ( const TByteStream &  str )
{
    if ( str._extern )
    {
        _data   = str._data;
        _size   = str._size;
        _pos    = str._pos;
        _extern = true;
    }// if
    else
    {
        set_size( str._size );

        memcpy( _data, str._data, _size );
        
        _pos    = str._pos;
        _extern = false;
    }// else
}

////////////////////////////////////////////////////
//
// set stream to hold external data
//

TByteStream &
TByteStream::set_size ( const size_t  n )
{
    if ( ! _extern && ( _data != nullptr ))
        delete _data;

    _data   = ( n > 0 ? new uchar[n] : nullptr );
    _size   = n;
    _extern = false;
    _pos    = 0;

    return *this;
}

//
// set internal stream to data and size WITHOUT copying
//
TByteStream &
TByteStream::set_stream ( void *        stream_data,
                          const size_t  stream_size )
{
    set_size( 0 );

    _size   = stream_size;
    _pos    = 0;
    _data   = static_cast< uchar * >( stream_data );
    _extern = true;

    return *this;
}
    
//
// set internal stream to data and size WITH copying
//
TByteStream &
TByteStream::copy_stream ( void *        stream_data,
                           const size_t  stream_size )
{
    set_size( stream_size );

    memcpy( _data, stream_data, stream_size );

    _pos    = 0;
    _extern = false;

    return *this;
}
    
////////////////////////////////////////////////////
//
// access-methods
//

//
// put n bytes from buffer into stream and increase position
// return 0 on success, else < 0
//
void
TByteStream::put ( const void *  buf,
                   const size_t  n )
{
    if ( _pos + n > size() )
        HERROR( ERR_BS_SIZE, "(TByteStream) put", "" );

    memcpy( _data + _pos, buf, n );
    _pos += n;
}

//
// copy n bytes from stream into buffer (increase position)
//
void
TByteStream::get ( void *        buf,
                   const size_t  n )
{
    if ( _pos + n > size() )
        HERROR( ERR_BS_SIZE, "(TByteStream) get", "" );

    memcpy( buf, _data + _pos, n );
    _pos += n;
}

//
// write stream into file
//
void
TByteStream::save ( const std::string &  filename ) const
{
    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();
    uint64_t                    t   = _size;
    
    out.write( reinterpret_cast< const char * >( & t ), sizeof(t) );
    out.write( reinterpret_cast< const char * >( _data ), _size );
}

//
// load stream from file
//
void
TByteStream::load ( const std::string & filename )
{
    if ( ! boost::filesystem::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TByteStream) read", filename );

    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();
    uint64_t                    t;
    
    in.read( reinterpret_cast< char * >( & t ), sizeof(t) );

    set_size( t );
    
    in.read( reinterpret_cast< char * >( _data ), _size );
}

////////////////////////////////////////////////////
//
// debugging
//

//
// return checksum of content
//
uint
TByteStream::checksum () const
{
    return Crypt::crc32( _size, _data );
}

//
// print information/content of stream
//
void
TByteStream::print ( const bool  show_content ) const
{
    if ( show_content )
        std::cout << "TByteStream ( size = " << _size << " )" << std::endl;
    
    if ( show_content )
    {
        for ( size_t  i = 0; i < _size; )
        {
            std::cout << boost::format( "%8d: " ) % i;
        
            for ( size_t  j = 0; j < 24 && i < _size; j++, i++ )
                std::cout << boost::format( "%02x " ) % _data[i];

            std::cout << std::endl;
        }// for
    }// if

    const uint crc = checksum();
    
    if ( show_content )
        std::cout << boost::format( "  crc = 0x%08x" ) % crc;
    else
        std::cout << "TByteStream " << boost::format( "( size = %ld , crc = 0x%08x )" ) % _size % crc;
}

}// namespace
