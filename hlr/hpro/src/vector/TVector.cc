//
// Project     : HLIBpro
// File        : TVector.cc
// Description : baseclass for all vector-classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <iostream>

#include "hpro/base/error.hh"
#include "hpro/parallel/NET.hh"
#include "hpro/io/TVectorIO.hh"

#include "hpro/vector/TVector.hh"

namespace Hpro
{

using std::unique_ptr;
using std::make_unique;

////////////////////////////////////////////////
//
// access entries
//

template < typename value_t >
value_t
TVector< value_t >::entry  ( const idx_t ) const
{
    HERROR( ERR_NOT_IMPL, "(TVector) entry", "" );
}

template < typename value_t >
void
TVector< value_t >::set_entry  ( const idx_t, const value_t )
{
    HERROR( ERR_NOT_IMPL, "(TVector) set_entry", "" );
}

////////////////////////////////////////////////
//
// memory consumption
//

//
// return size in bytes used by this distributed object,
// i.e. of all distributed sub matrices
//
template < typename value_t >
size_t
TVector< value_t >::global_byte_size () const
{
    return byte_size();
}

////////////////////////////////////////////////
//
// restriction
//

//
// create vector restricted to real/imaginary part of coefficients
//
template < typename value_t >
auto
TVector< value_t >::restrict_re () const -> std::unique_ptr< TVector< real_t > >
{
    HERROR( ERR_NOT_IMPL, "(TVector) restrict_re", "" );
}

template < typename value_t >
auto
TVector< value_t >::restrict_im () const -> std::unique_ptr< TVector< real_t > >
{
    HERROR( ERR_NOT_IMPL, "(TVector) restrict_im", "" );
}
    
template < typename value_t >
auto
TVector< value_t >::restrict_abs () const -> std::unique_ptr< TVector< real_t > >
{
    HERROR( ERR_NOT_IMPL, "(TVector) restrict_abs", "" );
}
    
////////////////////////////////////////////////
//
// serialisation
//

template < typename value_t >
void
TVector< value_t >::read  ( TByteStream & s )
{
    typeid_t  t;

    TStreamable::read( s );
    
    s.get( t );

    if ( t != type() )
        HERROR( ERR_BS_TYPE, "(TVector) read", "found " + RTTI::id_to_type( t ) + ","
               + ", expected " + typestr() );
    
    s.get( _ofs );
}

template < typename value_t >
void
TVector< value_t >::write ( TByteStream & s ) const
{
    typeid_t  t = type();

    TStreamable::write( s );
    
    s.put( t );
    s.put( _ofs );
}

//
// returns size of object in bytestream
//
template < typename value_t >
size_t
TVector< value_t >::bs_size () const
{
    return TStreamable::bs_size() + sizeof(typeid_t) + sizeof(_ofs);
}

//
// parallel methods
//

//
// sum up nparts parallel copies
// (if bs != nullptr it will be used)
//
template < typename value_t >
void
TVector< value_t >::sum ( const TProcSet &  ps,
                          const uint        pid,
                          const uint        parts,
                          TByteStream *     bs )
{
    unique_ptr< TVector< value_t > >  v;
    unique_ptr< TByteStream >         tbs;
    uint                              nparts = parts;
    
    if ( bs == nullptr )
    {
        tbs = make_unique< TByteStream >();
        bs  = tbs.get();
    }// if

    while ( nparts > 1 )
    {
        std::vector< TProcSet >  psets;
        int                      idx = -1;

        // create processor sets for 
        ps.split( nparts, psets );

        for ( uint i = 0; i < nparts; i++ )
            if ( pid == psets[i].master() )
            {
                idx = i;
                break;
            }// if

        nparts /= 2;
                
        // if not involved in communication, continue
        if ( idx == -1 )
            continue;

        // send data at all "odd" processors
        if ( idx % 2 == 1 )
        {
            bs->set_size( bs_size() );
            
            write( * bs );
            
            NET::dsend( psets[idx-1].master(), bs->data(), bs->size() );
        }// if
        
        // and receive it at all "even" procs.
        if ( idx % 2 == 0 )
        {
            bs->set_size( NET::dprobe( psets[idx+1].master() ) );
            
            NET::drecv( psets[idx+1].master(), bs->data(), bs->size() );
            
            if ( v.get() == nullptr )
                v = unique_ptr< TVector< value_t > >( create() );
            
            v->read( * bs );

            axpy( value_t(1), v.get() );
        }// if
    }// while
}

//
// sum up all vectors in parallel machine defined by <procs>
//
template < typename value_t >
void
TVector< value_t >::sum ( const TProcSet & procs )
{
    unique_ptr< TVector< value_t > >  v;
    uint                              nparts = procs.size();
    TByteStream                       bs;
    const uint                        pid    = NET::pid();

    LOG::print( "sum()" );
    
    while ( nparts > 1 )
    {
        std::vector< TProcSet >  psets;
        int                      idx = -1;

        // create processor sets for 
        procs.split( nparts, psets );

        for ( uint i = 0; i < nparts; i++ )
            if ( pid == psets[i].master() )
            {
                idx = i;
                break;
            }// if

        nparts /= 2;
                
        // if not involved in communication, continue
        if ( idx == -1 )
            continue;

        if ( idx % 2 == 1 )
        {
            //
            // send data at all "odd" processors
            //
            
            bs.set_size( bs_size() );
            
            write( bs );
            
            NET::dsend( psets[idx-1].master(), bs.data(), bs.size() );
        }// if
        else // if ( idx % 2 == 0 )
        {
            //
            // and receive it at all "even" procs.
            //
            
            bs.set_size( NET::dprobe( psets[idx+1].master() ) );
            
            NET::drecv( psets[idx+1].master(), bs.data(), bs.size() );

            if ( v.get() == nullptr )
                v = unique_ptr< TVector< value_t > >( create() );
            
            v->read( bs );

            axpy( value_t(1), v.get() );
        }// if
    }// while

    //
    // finally, scatter result back from master to all other
    // nodes in local processor group
    //

    LOG::print( "scatter" );
    
    // scatter( procs );
}

//
// stream output
//
template < typename value_t >
void
TVector< value_t >::print ( const uint offset ) const
{
    for ( uint i = 0; i < offset; i++ )
        std::cout << ' ';

    std::cout << typestr() << std::endl;
}

template class TVector< float >;
template class TVector< double >;
template class TVector< std::complex< float > >;
template class TVector< std::complex< double > >;

//////////////////////////////////////////////////////////
//
// debug helpers
//

namespace DBG
{

//
// write vector to file
//
template < typename value_t >
void
write ( const TVector< value_t > *  v,
        const std::string &         filename,
        const std::string &         vecname )
{
    TMatlabVectorIO  vio;

    vio.write( v, filename, vecname );
}

#define INST_WRITE( type ) \
    template void write< type > ( const TVector< type > *, \
                                  const std::string &,     \
                                  const std::string & );

INST_WRITE( float )
INST_WRITE( double )
INST_WRITE( std::complex< float > )
INST_WRITE( std::complex< double > )

}// namespace DBG

}// namespace Hpro
