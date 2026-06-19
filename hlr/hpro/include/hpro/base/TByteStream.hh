#ifndef __HPRO_TBYTESTREAM_HH
#define __HPRO_TBYTESTREAM_HH
//
// Project     : HLIBpro
// File        : TByteStream.hh
// Description : class for a byte stream
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <stddef.h>

#include "hpro/base/error.hh"
#include "hpro/base/types.hh"

namespace Hpro
{

//!
//! \class TByteStream
//! \brief implements a stream of bytes for storage purposes
//!
class TByteStream
{
protected:
    //! array holding the actual stream
    uchar  * _data;
    size_t   _size;

    //! current position in stream
    size_t   _pos;

    //! indicates, that stream memory is handled outside of class
    bool     _extern;
    
public:

    ////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct bytestream capable of holding \a s bytes
    TByteStream ()
            : _data(NULL)
            , _size(0)
            , _pos(0)
            , _extern(false)
    {}

    //! construct bytestream capable of holding \a s bytes
    TByteStream ( const size_t  s )
            : _data(NULL)
            , _size(0)
            , _pos(0)
            , _extern(false)
    { set_size(s); }

    //! construct bytestream by given array (external storage)
    TByteStream ( const size_t  asize,
                  void *        adata )
            : _data(NULL)
            , _size(0)
            , _pos(0)
            , _extern(false)
    { set_stream( adata, asize ); }

    //! copy constructor
    TByteStream ( const TByteStream & str );

    //! destructor
    ~TByteStream () { set_size(0); }

    ////////////////////////////////////////////////////
    //
    // access internal data
    //

    //! return current position in byte stream
    size_t  position () const { return _pos; }

    //! return capacity of bytestream
    size_t  size     () const { return _size; }

    //! set new bytestream capacity without keeping old data
    TByteStream & set_size ( const size_t n );

    //! set new position in bytestream
    TByteStream & set_pos  ( const size_t p ) { _pos = p; return *this; }

    //! set internal stream to \a data and \a size <b>without</b> copying
    TByteStream & set_stream  ( void * data, const size_t size );

    //! set internal stream to \a data and \a size <b>with</b> copying
    TByteStream & copy_stream ( void * data, const size_t size );

    //! set internal pointer to start of stream
    TByteStream & to_start () { _pos = 0;      return *this; }

    //! set internal pointer to end of stream
    TByteStream & to_end   () { _pos = size(); return *this; }

    //! directly access internal stream data
    uchar *       data ()       { return _data; }
    const uchar * data () const { return _data; }
    
    ////////////////////////////////////////////////////
    //
    // stream-functionality
    //

    //! put \a n bytes to current position in stream
    void put  ( const void * buf, const size_t n );

    //! directly put given data into stream
    template <typename T>
    void put  ( const T &  buf )                     { put( & buf, sizeof(T) ); }

    //! directly put \a nelem elements of data into stream
    template <typename T>
    void dput ( const T * buf, const size_t  nelem ) { put( buf, sizeof(T) * nelem ); }
    
    //! get \a n bytes from current position in stream
    void get  ( void * buf, const size_t n );

    //! directly get data from stream
    template <typename T>
    void get  ( T &  buf )                     { get( & buf, sizeof(T) ); }

    //! directly get \a nelem elements of data from stream
    template <typename T>
    void dget ( T * buf, const size_t  nelem ) { get( buf, sizeof(T) * nelem ); }
    
    //! write stream into file \a filename
    void save ( const std::string & filename ) const;

    //! load stream from file \a filename
    void load ( const std::string & filename );

    ////////////////////////////////////////////////////
    //
    // debugging
    //

    //! return checksum of content
    uint checksum () const;
    
    //! print information/content of stream
    void print ( const bool show_content = false ) const;
};

}// namespace

#endif // __HPRO_TBYTESTREAM_HH
