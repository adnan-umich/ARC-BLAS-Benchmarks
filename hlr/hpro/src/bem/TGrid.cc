//
// Project     : HLIBpro
// File        : TGrid.hh
// Description : contains information about a grid
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/bem/TGrid.hh"

namespace Hpro
{

namespace B = BLAS;

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// TGrid
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
//
// constructor and destructor
//
///////////////////////////////////////////////////////////////////////

TGrid::TGrid ( const std::vector< T3Point > &     vertex_arr,
               const std::vector< triangle_t > &  triangle_arr )
{
    _vertices  = vertex_arr;
    _triangles = triangle_arr;
    
    comp_triangle_size_normal();
}

TGrid::TGrid ( const std::vector< T3Point > &     vertex_arr,
               const std::vector< triangle_t > &  triangle_arr,
               const std::vector< T3Point > &     vertex_normal )
{
    _vertices  = vertex_arr;
    _triangles = triangle_arr;

    add_vtx_normal( vertex_normal );
            
    comp_triangle_size_normal();
}


TGrid::~TGrid ()
{}

///////////////////////////////////////////////////////////////////////
//
// modify grid
//
///////////////////////////////////////////////////////////////////////

//
// translate grid by given vector
//
void
TGrid::translate ( const T3Point & t )
{
    //
    // change point positions and leave rest as it is
    //

    const size_t  nvertices = n_vertices();
    
    for ( size_t  i = 0; i < nvertices; i++ )
    {
        _vertices[i] += t;
    }// for
}

//
// scale grid by x, y and z components in given vector
//
void
TGrid::scale ( const T3Point & s )
{
    //
    // modify point positions and update normals, triangles sizes
    //

    const size_t  nvertices = n_vertices();
    
    for ( size_t  i = 0; i < nvertices; i++ )
    {
        _vertices[i] *= s;
    }// for

    comp_triangle_size_normal();
}

//
// rotate grid around vector <v> by angle alpha
//
void
TGrid::rotate ( const T3Point & v, const double alpha )
{
    //
    // form rotation matrix by alpha around v = (x,y,z) is given by
    //
    //   | x^2+(1-x^2) cos a        xy(1-cos a) + z sin a     xz(1-cos a) - y sin a |
    //   | xy(1-cos a) - z sin a    y^2 + (1-y^2) cos a       yz(1-cos a) + x sin a |
    //   | xz(1-cos a) + y sin a    yz(1-cos a) - x sin a     z^2 + (1-z^2) cos a   |
    //
    
    const double         c   = Math::cos( alpha );
    const double         s   = Math::sin( alpha );
    const double         x   = v[0];
    const double         y   = v[1];
    const double         z   = v[2];
    const double         xx  = x * x;
    const double         yy  = y * y;
    const double         zz  = z * z;
    const double         omc = 1.0 - c;
    B::Matrix< double >  M( 3, 3 );
    
    M(0,0) = xx + c * ( 1.0 - xx );
    M(1,0) = x * y * omc - z * s;
    M(2,0) = z * x * omc + y * s;

    M(0,1) = x * y * omc + z * s;
    M(1,1) = yy + c * ( 1.0 - yy );
    M(2,1) = y * z * omc - x * s;
    
    M(0,2) = z * x * omc - y * s;
    M(1,2) = y * z * omc + x * s;
    M(2,2) = zz + c * ( 1.0 - zz );

    //
    // rotate each vertex in grid
    //
    
    const size_t  nvertices = n_vertices();
    
    for ( size_t  i = 0; i < nvertices; i++ )
    {
        B::Vector< double >  coord( 3 ), rot( 3 );

        coord(0) = _vertices[i][0];
        coord(1) = _vertices[i][1];
        coord(2) = _vertices[i][2];

        B::mulvec( 1.0, M, coord, 0.0, rot );

        _vertices[i][0] = rot(0);
        _vertices[i][1] = rot(1);
        _vertices[i][2] = rot(2);
    }// for

    //
    // recompute normal directions and triangle sizes
    //
    
    comp_triangle_size_normal();
}

//
// switch triangle orientation, e.g. normal direction
//
void
TGrid::switch_tri_orient ()
{
    //
    // switch vertices of each triangle
    //

    for ( std::vector< triangle_t >::iterator  iter = _triangles.begin();
          iter != _triangles.end();
          ++iter )
    {
        std::swap( (*iter).vtx[1], (*iter).vtx[2] );
    }// for

    //
    // and recompute triangle normals
    //

    comp_triangle_size_normal();

    //
    // if vertex normals are present, switch them too
    //

    for ( std::vector< T3Point >::iterator  iter = _vtx_normal.begin();
          iter != _vtx_normal.end();
          ++iter )
    {
        (*iter) *= double(-1);
    }// for
}

//
// add given vertex normals to grid
//
void
TGrid::add_vtx_normal ( const std::vector< T3Point > &  vertex_normal )
{
    _vtx_normal = vertex_normal;

    if ( _vtx_normal.size() != _vertices.size() )
        HERROR( ERR_ARG, "(TGrid) add_vtx_normal",
                "size of vertex normal array is different from vertex array size" );
}

///////////////////////////////////////////////////////////////////////
//
// misc.
//
///////////////////////////////////////////////////////////////////////

//
// return size of grid in bytes
//
size_t
TGrid::byte_size () const
{
    size_t  size = 0;
    
    size += sizeof(_vertices)   + _vertices.size()   * sizeof(T3Point);
    size += sizeof(_triangles)  + _triangles.size()  * sizeof(triangle_t);
    size += sizeof(_tri_size)   + _tri_size.size()   * sizeof(double);
    size += sizeof(_tri_normal) + _tri_normal.size() * sizeof(T3Point);
    size += sizeof(_vtx_normal) + _vtx_normal.size() * sizeof(T3Point);

    return size;
}

///////////////////////////////////////////////////////////////////////
//
// internal routines
//
///////////////////////////////////////////////////////////////////////

//
// compute size and normal direction of each triangle
//
void
TGrid::comp_triangle_size_normal ()
{
    const size_t  ntriangles = n_triangles();
    
    _tri_size.resize( ntriangles );
    _tri_normal.resize( ntriangles );

    for ( size_t  i = 0; i < ntriangles; i++ )
    {
        idx_t   vid[3];
        T3Point X1, X2, v;

        vid[0] = _triangles[i].vtx[0];
        vid[1] = _triangles[i].vtx[1];
        vid[2] = _triangles[i].vtx[2];
        
        X1 = _vertices[vid[1]] - _vertices[vid[0]];
        X2 = _vertices[vid[2]] - _vertices[vid[0]];
        
        v = cross( X1, X2 );
        
        _tri_size[i] = v.norm2();
        
        v.normalise2();
        _tri_normal[i] = v;
    }// for
}

}// namespace
