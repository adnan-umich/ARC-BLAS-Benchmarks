//
// Project     : HLIBpro
// File        : TCoordinate.cc
// Description : class for encapsulating coordinate data
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/error.hh"

#include "hpro/cluster/TCoordinate.hh"

namespace Hpro
{

////////////////////////////////////////////////////////
//
// constructor and destructor
//

TCoordinate::TCoordinate ( const std::vector< double * > &  acoord,
                           const uint                       adim,
                           const coord_data_t               coord_data )
        : _dim( adim )
        , _coord_data( coord_data )
{
    const size_t  n = acoord.size();

    if ( n == 0 )
        HERROR( ERR_ARG, "(TCoordinate) ctor", "number of coordinates is 0" );

    if ( adim == 0 )
        HERROR( ERR_ARG, "(TCoordinate) ctor", "dimension of coordinates is 0" );

    _coord.resize( n );

    if ( coord_data == copy_coord_data )
    {
        for ( size_t  i = 0; i < n; i++ )
        {
            _coord[i] = new double[ _dim ];

            for ( uint  j = 0; j < _dim; ++j )
                _coord[i][j] = acoord[i][j];
        }// for
    }// if
    else
    {
        for ( size_t  i = 0; i < n; i++ )
            _coord[i] = acoord[i];
    }// else
}
                  
TCoordinate::TCoordinate ( const std::vector< double * > &  acoord,
                           const uint                       adim,
                           const std::vector< double * > &  abbmin,
                           const std::vector< double * > &  abbmax,
                           const coord_data_t               coord_data )
        : _dim( adim )
        , _coord_data( coord_data )
{
    if (( acoord.size() != abbmin.size() ) ||
        ( acoord.size() != abbmax.size() ))
        HERROR( ERR_ARG, "(TCoordinate)",
                "bounding box arrays and coord. array differ in size" );
    
    const size_t  n = acoord.size();

    if ( n == 0 )
        HERROR( ERR_ARG, "(TCoordinate) ctor", "number of coordinates is 0" );

    if ( adim == 0 )
        HERROR( ERR_ARG, "(TCoordinate) ctor", "dimension of coordinates is 0" );

    _coord.resize( n );
    _bbmin.resize( n );
    _bbmax.resize( n );

    if ( coord_data == copy_coord_data )
    {
        for ( size_t  i = 0; i < n; i++ )
        {
            _coord[i] = new double[ _dim ];
            _bbmin[i] = new double[ _dim ];
            _bbmax[i] = new double[ _dim ];

            for ( uint  j = 0; j < _dim; ++j )
            {
                _coord[i][j] = acoord[i][j];
                _bbmin[i][j] = abbmin[i][j];
                _bbmax[i][j] = abbmax[i][j];
            }// for
        }// for
    }// if
    else
    {
        for ( size_t  i = 0; i < n; i++ )
        {
            _coord[i] = acoord[i];
            _bbmin[i] = abbmin[i];
            _bbmax[i] = abbmax[i];
        }// for
    }// else
}

//
// special versions for TPoint (always copy data)
//
TCoordinate::TCoordinate ( const std::vector< TPoint > &  acoord )
        : _dim( acoord.size() > 0 ? acoord[0].dim() : 0 )
        , _coord_data( copy_coord_data )
{
    const size_t  n = acoord.size();

    if ( n == 0 )
        HERROR( ERR_ARG, "(TCoordinate) ctor", "number of coordinates is 0" );

    _coord.resize( n );

    for ( size_t  i = 0; i < n; i++ )
    {
        _coord[i] = new double[ _dim ];

        if ( acoord[i].dim() != _dim )
            HERROR( ERR_CONSISTENCY, "(TCoordinate) ctor", "coordinates have differing dimension" );
        
        for ( uint  j = 0; j < _dim; j++ )
            _coord[i][j] = acoord[i][j];
    }// for
}

TCoordinate::TCoordinate ( const std::vector< TPoint > &   acoord,
                           const std::vector< TBBox > &    abbox )
        : _dim( acoord.size() > 0 ? acoord[0].dim() : 0 )
        , _coord_data( copy_coord_data )
{
    if ( acoord.size() != abbox.size() )
        HERROR( ERR_ARG, "(TCoordinate)",
                "bounding box array and coord. array differ in size" );
    
    const size_t  n = acoord.size();

    if ( n == 0 )
        HERROR( ERR_ARG, "(TCoordinate) ctor", "number of coordinates is 0" );

    if ( abbox[0].dim() != 2 )
        HERROR( ERR_ARG, "(TCoordinate) ctor", "bbox has not spatial dimension 2" );

    _coord.resize( n );
    _bbmin.resize( n );
    _bbmax.resize( n );

    for ( size_t  i = 0; i < n; i++ )
    {
        _coord[i] = new double[ _dim ];
        _bbmin[i] = new double[ _dim ];
        _bbmax[i] = new double[ _dim ];
        
        if ( acoord[i].dim() != _dim )
            HERROR( ERR_CONSISTENCY, "(TCoordinate) ctor", "coordinates have differing dimension" );
        
        if ( abbox[i].dim() != _dim )
            HERROR( ERR_CONSISTENCY, "(TCoordinate) ctor", "coordinates have differing dimension" );
        
        for ( uint  j = 0; j < _dim; j++ )
        {
            _coord[i][j] = acoord[i][j];
            _bbmin[i][j] = abbox[i].min()[j];
            _bbmax[i][j] = abbox[i].max()[j];
        }// for
    }// for
}

//
// special versions for T2Point (always copy data)
//
TCoordinate::TCoordinate ( const std::vector< T2Point > &   acoord )
        : _dim( 2 )
        , _coord_data( copy_coord_data )
{
    const size_t  n = acoord.size();

    if ( n == 0 )
        HERROR( ERR_ARG, "(TCoordinate) ctor", "number of coordinates is 0" );

    _coord.resize( n );

    for ( size_t  i = 0; i < n; i++ )
    {
        _coord[i] = new double[ _dim ];
        
        _coord[i][0] = acoord[i][0];
        _coord[i][1] = acoord[i][1];
    }// for
}

TCoordinate::TCoordinate ( const std::vector< T2Point > &   acoord,
                           const std::vector< TBBox > &     abbox )
        : _dim( 2 )
        , _coord_data( copy_coord_data )
{
    if ( acoord.size() != abbox.size() )
        HERROR( ERR_ARG, "(TCoordinate)",
                "bounding box array and coord. array differ in size" );
    
    const size_t  n = acoord.size();

    if ( n == 0 )
        HERROR( ERR_ARG, "(TCoordinate) ctor", "number of coordinates is 0" );

    if ( abbox[0].dim() != 2 )
        HERROR( ERR_ARG, "(TCoordinate) ctor", "bbox has not spatial dimension 2" );

    _coord.resize( n );
    _bbmin.resize( n );
    _bbmax.resize( n );

    for ( size_t  i = 0; i < n; i++ )
    {
        _coord[i] = new double[ _dim ];
        _bbmin[i] = new double[ _dim ];
        _bbmax[i] = new double[ _dim ];
        
        _coord[i][0] = acoord[i][0];
        _coord[i][1] = acoord[i][1];

        _bbmin[i][0] = abbox[i].min()[0];
        _bbmin[i][1] = abbox[i].min()[1];

        _bbmax[i][0] = abbox[i].max()[0];
        _bbmax[i][1] = abbox[i].max()[1];
    }// for
}

//
// special versions for T3Point (always copy data)
//
TCoordinate::TCoordinate ( const std::vector< T3Point > &   acoord )
        : _dim( 3 )
        , _coord_data( copy_coord_data )
{
    const size_t  n = acoord.size();

    if ( n == 0 )
        HERROR( ERR_ARG, "(TCoordinate) ctor", "number of coordinates is 0" );

    _coord.resize( n );

    for ( size_t  i = 0; i < n; i++ )
    {
        _coord[i] = new double[ _dim ];
        
        _coord[i][0] = acoord[i][0];
        _coord[i][1] = acoord[i][1];
        _coord[i][2] = acoord[i][2];
    }// for
}
    
TCoordinate::TCoordinate ( const std::vector< T3Point > &   acoord,
                           const std::vector< TBBox > &     abbox )
        : _dim( 3 )
        , _coord_data( copy_coord_data )
{
    if ( acoord.size() != abbox.size() )
        HERROR( ERR_ARG, "(TCoordinate)",
                "bounding box array and coord. array differ in size" );
    
    const size_t  n = acoord.size();

    if ( n == 0 )
        HERROR( ERR_ARG, "(TCoordinate) ctor", "number of coordinates is 0" );

    if ( abbox[0].dim() != 3 )
        HERROR( ERR_ARG, "(TCoordinate) ctor", "bbox has not spatial dimension 3" );

    _coord.resize( n );
    _bbmin.resize( n );
    _bbmax.resize( n );

    for ( size_t  i = 0; i < n; i++ )
    {
        _coord[i] = new double[ _dim ];
        _bbmin[i] = new double[ _dim ];
        _bbmax[i] = new double[ _dim ];
        
        _coord[i][0] = acoord[i][0];
        _coord[i][1] = acoord[i][1];
        _coord[i][2] = acoord[i][2];

        _bbmin[i][0] = abbox[i].min()[0];
        _bbmin[i][1] = abbox[i].min()[1];
        _bbmin[i][2] = abbox[i].min()[2];

        _bbmax[i][0] = abbox[i].max()[0];
        _bbmax[i][1] = abbox[i].max()[1];
        _bbmax[i][2] = abbox[i].max()[2];
    }// for
}

//
// dtor
//
TCoordinate::~TCoordinate ()
{
    if ( _coord_data == copy_coord_data )
    {
        const size_t  n = ncoord();

        if ( has_bbox() )
        {
            for ( size_t  i = 0; i < n; i++ )
            {
                delete[] _bbmin[i];
                delete[] _bbmax[i];
            }// for
        }// if
        
        for ( size_t  i = 0; i < n; i++ )
            delete[] _coord[i];
    }// if
}

////////////////////////////////////////////////////////
//
// misc.
//

//
// return boundning box of coordinate set
//
TBBox
TCoordinate::bounding_box () const
{
    TBBox  bbox;

    if ( ncoord() == 0 )
        return bbox;
    
    if ( has_bbox() )
    {
        //
        // bounding box of coordinate bounding boxes
        //
        
        bbox.min() = TPoint( dim(), bbmin( 0 ) );
        bbox.max() = TPoint( dim(), bbmax( 0 ) );

        for ( size_t  i = 1; i < ncoord(); ++i )
        {
            const TPoint  bmin( dim(), bbmin( uint(i) ) );
            const TPoint  bmax( dim(), bbmax( uint(i) ) );

            bbox.join( TBBox( bmin, bmax ) );
        }// for
    }// if
    else
    {
        //
        // bounding box of coordinates itself
        //

        bbox.min() = TPoint( dim(), coord( 0 ) );
        bbox.max() = bbox.min();

        for ( idx_t  i = 1; i < idx_t(ncoord()); ++i )
        {
            const TPoint  v( dim(), coord( i ) );

            bbox.join( TBBox( v, v ) );
        }// for
    }// else

    return  bbox;
}

//
// return memory consumption
//
size_t
TCoordinate::byte_size () const
{
    return ( sizeof(_coord) + sizeof(double*) * _coord.size() +
             _coord.size() * sizeof(double) * dim() +
             sizeof(_bbmin) + sizeof(double*) * _bbmin.size() +
             _bbmin.size() * sizeof(double) * dim() +
             sizeof(_bbmax) + sizeof(double*) * _bbmax.size() +
             _bbmax.size() * sizeof(double) * dim() +
             sizeof(_dim) +
             _period.byte_size() +
             sizeof( _coord_data ) );
}

}// namespace
