//
// Project     : HLIBpro
// File        : TFnSpace.cc
// Description : implements function spaces over grids
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/bem/TFnSpace.hh"

namespace Hpro
{

using std::unique_ptr;
using std::make_unique;

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// TFnSpace
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
//
// constructor and destructor
//

template < typename value_t >
TFnSpace< value_t >::TFnSpace ( const TGrid * agrid )
        : _grid(agrid)
{
    if ( agrid == nullptr )
        HERROR( ERR_ARG, "(TFnSpace)", "given grid is nullptr" );
}

template < typename value_t >
TFnSpace< value_t >::~TFnSpace ()
{
}

////////////////////////////////////////////////////////
//
// misc.
//

//
// return coordinate info for grid
//
template < typename value_t >
unique_ptr< TCoordinate >
TFnSpace< value_t >::build_coord ( const bool  with_bbox ) const
{
    const size_t            nindices = _indices.size();
    std::vector< T3Point >  coord( nindices );

    if ( with_bbox )
    {
        std::vector< TBBox >  bbox( nindices );

        for ( size_t  i = 0; i < nindices; ++i )
        {
            //
            // start with position of index
            //
                
            coord[i] = _indices[i];

            TBBox  bbox_i( coord[i], coord[i] );
                
            //
            // and enhance by support, e.g. coordinates of
            // triangles in support
            //
                
            for ( auto  tri : support( idx_t(i) ) )
            {
                // loop over the three triangle vertices
                for ( uint k = 0; k < 3; k++ )
                {
                    const T3Point  pos( _grid->vertex( _grid->triangle( tri ).vtx[k] ) );
                    const TBBox    pos_bbox( pos, pos );
                        
                    bbox_i.join( pos_bbox );
                }// for
            }// for

            bbox[i] = bbox_i;
        }// for

        return make_unique< TCoordinate >( coord, bbox );
    }// if
    else 
    {
        for ( size_t  i = 0; i < nindices; i++ )
        {
            coord[i] = _indices[i];
        }// for

        return make_unique< TCoordinate >( coord );
    }// if
}

//
// return size of grid in bytes
//
template < typename value_t >
size_t
TFnSpace< value_t >::byte_size () const
{
    size_t  size = sizeof(TGrid*);
    
    size += sizeof( _supp_list_ptr ) + _supp_list_ptr.size() * sizeof( idx_t );
    size += sizeof( _supp_list ) + _supp_list.size() * sizeof( idx_t );
    size += sizeof( _indices ) + _indices.size() * sizeof( T3Point );

    return size;
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// TConstFnSpace
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
//
// constructor and destructor
//

template < typename value_t >
TConstFnSpace< value_t >::TConstFnSpace ( const TGrid * agrid )
        : TFnSpace< value_t >( agrid )
{
    construct();
}

template < typename value_t >
TConstFnSpace< value_t >::~TConstFnSpace ()
{
}
    
////////////////////////////////////////////////////////
//
// private functions
//

//
// construct function space by building indices
// and their support
//
template < typename value_t >
void
TConstFnSpace< value_t >::construct ()
{
    //
    // indices are built at the midpoints of each triangle
    //

    const size_t  n_triangles = this->_grid->n_triangles();

    this->_indices.resize( n_triangles );

    for ( size_t  i = 0; i < n_triangles; i++ )
    {
        const idx_t  vid[3] = { this->_grid->triangle( idx_t(i) ).vtx[0],
                                this->_grid->triangle( idx_t(i) ).vtx[1],
                                this->_grid->triangle( idx_t(i) ).vtx[2] };
            
        this->_indices[i]  = T3Point( this->_grid->vertex( vid[0] ) +
                                      this->_grid->vertex( vid[1] ) +
                                      this->_grid->vertex( vid[2] ) );
        this->_indices[i] *= 1.0 / 3.0;
    }// for

    //
    // support is defined by the single triangle enclosing index i
    // and since they are numbered in the same order, the mapping
    // is just the identity
    //

    this->_supp_list_ptr.resize( n_triangles+1, 0 );
    this->_supp_list.resize( n_triangles, 0 );

    for ( uint i = 0; i < n_triangles; i++ )
    {
        this->_supp_list_ptr[i] = i;
        this->_supp_list[i]     = i;
    }// for
        
    this->_supp_list_ptr[n_triangles] = uint( n_triangles );

    //
    // build mapping of triangle to triangle local indices
    // - one index per triangle
    //

    this->_tri_idx_ptr.resize( n_triangles+1, 0 );
    this->_tri_idx.resize( n_triangles+1, 0 );

    for ( uint i = 0; i < n_triangles; i++ )
    {
        this->_tri_idx_ptr[i] = i;
        this->_tri_idx[i]     = i;
    }// for
        
    this->_tri_idx_ptr[n_triangles] = uint( n_triangles );
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// TLinFnSpace
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
//
// constructor and destructor
//

template < typename value_t >
TLinearFnSpace< value_t >::TLinearFnSpace ( const TGrid * agrid )
        : TFnSpace< value_t >( agrid )
{
    construct();
}

template < typename value_t >
TLinearFnSpace< value_t >::~TLinearFnSpace ()
{
}
    
////////////////////////////////////////////////////////
//
// private functions
//

//
// construct function space by building indices
// and their support
//
template < typename value_t >
void
TLinearFnSpace< value_t >::construct ()
{
    //
    // just copy position of vertices to index-array
    //

    const size_t  n_vertices = this->_grid->n_vertices();
        
    this->_indices.resize( n_vertices );

    for ( idx_t  i = 0; i < idx_t( n_vertices ); ++i )
    {
        this->_indices[i] = this->_grid->vertex( i );
    }// for
        
    //
    // count number of triangles in support per index and
    // built start-end-array of positions in the support-list array
    //

    const size_t         n_triangles = this->_grid->n_triangles();
    std::vector< uint >  count( n_vertices, 0 );
    
    for ( idx_t  i = 0; i < idx_t( n_triangles ); ++i )
    {
        count[ this->_grid->triangle(i).vtx[0] ]++;
        count[ this->_grid->triangle(i).vtx[1] ]++;
        count[ this->_grid->triangle(i).vtx[2] ]++;
    }// for

    idx_t  pos = 0;

    this->_supp_list_ptr.resize( n_vertices+1, 0 );

    for ( size_t  i = 0; i < n_vertices; ++i )
    {
        this->_supp_list_ptr[i] = pos;
        pos                    += idx_t( count[i] );
    }// for

    this->_supp_list_ptr[n_vertices] = pos;

    //
    // fill support-list per index with the corresponding triangle indices
    // - pos now holds the total number of references between indices
    //   and triangles
    //

    this->_supp_list.resize( pos, 0 );
    
    for ( size_t  i = 0; i < n_vertices; ++i )
        count[i] = 0;

    for ( size_t  i = 0; i < n_triangles; ++i )
    {
        for ( uint  j = 0; j < 3; j++ )
        {
            const idx_t  vid = this->_grid->triangle(  idx_t(i) ).vtx[j];
            
            this->_supp_list[ this->_supp_list_ptr[ vid ] + count[vid] ] = idx_t( i );
            count[vid]++;
        }// for
    }// for

    //
    // build mapping of triangle to triangle local indices
    // - three indices per triangle
    //

    this->_tri_idx_ptr.resize( n_triangles+1, 0 );
    this->_tri_idx.resize( 3 * n_triangles, 0 );

    pos = 0;
    
    for ( size_t  i = 0; i < n_triangles; ++i )
    {
        this->_tri_idx_ptr[ i ] = pos;
        
        for ( uint j = 0; j < 3; j++ )
        {
            this->_tri_idx[ pos ] = this->_grid->triangle( idx_t(i) ).vtx[j];
            ++pos;
        }// for
    }// for

    this->_tri_idx_ptr[ n_triangles ] = pos;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

#define INST_ALL( type )                        \
    template class TFnSpace< type >;            \
    template class TConstFnSpace< type >;       \
    template class TLinearFnSpace< type >;

INST_ALL( float )
INST_ALL( double )

template class TFnSpace< T2Point >;

}// namespace Hpro
