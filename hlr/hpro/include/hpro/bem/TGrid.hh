#ifndef __HPRO_TGRID_HH
#define __HPRO_TGRID_HH
//
// Project     : HLIBpro
// File        : TGrid.hh
// Description : contains information about a grid
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/base/types.hh"
#include "hpro/base/config.hh"
#include "hpro/base/TPoint.hh"
#include "hpro/vector/TScalarVector.hh"
#include "hpro/bem/TBEMFunction.hh"

namespace Hpro
{
    
///////////////////////////////////////////////////////////////
//
// TGrid: stores coordinates and triangle information
//
///////////////////////////////////////////////////////////////

class TGrid
{
public:
    //!
    //! type for storing triangles
    //!
    struct triangle_t
    {
        idx_t  vtx[3];
    };
    
protected:
    //! coordinates of the vertices in the grid
    std::vector< T3Point >     _vertices;

    //! list of triangles defined by 3 vertex-indices per triangle
    std::vector< triangle_t >  _triangles;
    
    //! contains size of each triangle
    std::vector< double >      _tri_size;

    //! contains normal direction of each triangle
    std::vector< T3Point >     _tri_normal;
    
    //! contains normal direction for each vertex
    std::vector< T3Point >     _vtx_normal;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct grid with vertices \a vertex_arr and triangles
    //! in \a triangle_arr; only constant normal direction per
    //! triangle supported
    TGrid ( const std::vector< T3Point > &     vertex_arr,
            const std::vector< triangle_t > &  triangle_arr );

    //! construct grid with vertices \a vertex_arr and triangles
    //! in \a triangle_arr; the normal direction of a triangle
    //! is interpolated by normal direction of vertices
    //! (if tri_normal(idx,s,t) is used)
    TGrid ( const std::vector< T3Point > &     vertex_arr,
            const std::vector< triangle_t > &  triangle_arr,
            const std::vector< T3Point > &     vertex_normal );

    ~TGrid ();

    //////////////////////////////////////
    //
    // give access to local data
    //

    //! return number of vertices
    size_t                             n_vertices  ()                 const { return _vertices.size(); }

    //! return full vertex array
    const std::vector< T3Point > &     vertices    ()                 const  { return _vertices; }

    //! return vertex \a i
    const T3Point &                    vertex      ( const idx_t  i ) const { return _vertices[i]; }

    //! return number of triangles
    size_t                             n_triangles ()                 const { return _triangles.size(); }
    
    //! return full triangle array
    const std::vector< triangle_t > &  triangles   ()                 const { return _triangles; }
    
    //! return triangle \a i
    triangle_t                         triangle    ( const idx_t  i ) const { return _triangles[i]; }
    
    //////////////////////////////////////
    //
    // modify grid
    //

    //! translate grid by vector \a t
    void  translate  ( const T3Point &  t );

    //! scale grid by (x, y, z) defined by vector \a s
    void  scale      ( const T3Point &  s );

    //! rotate grid around vector \a v by angle \a alpha
    void  rotate     ( const T3Point &  v,
                       const double     alpha );

    //! switch triangle orientation, e.g. normal direction
    void  switch_tri_orient ();
    
    //////////////////////////////////////
    //
    // provide additional data
    //

    //! return size of triangle \a tri
    double   tri_size   ( const idx_t         tri ) const { return _tri_size[ tri ]; }

    //! return normal direction of triangle \a tri
    T3Point  tri_normal ( const idx_t         tri ) const { return _tri_normal[tri]; }
    
    //! return normal direction at local position (\a s, \a t)
    //! of triangle \a triidx; (\a s, \a t) is defined with respect to \a ptri,
    //! e.g. permuted, and not according to original vertex ordering of \a triidx
    T3Point  tri_normal ( const idx_t         triidx,
                          const triangle_t &  tri,
                          const double        s,
                          const double        t ) const
    {
        if ( has_vtx_normal() )
        {
            return ( (1.0-s-t) * _vtx_normal[tri.vtx[0]] +
                     s         * _vtx_normal[tri.vtx[1]] +
                     t         * _vtx_normal[tri.vtx[2]] ).normalise2();
        }// if
        else
        {
            // use standard normal direction
            return _tri_normal[triidx];
        }// else
    }

    //! return true, if grid has vertex normal data instead of triangle normals
    bool  has_vtx_normal () const { return ( _vtx_normal.size() > 0 ); }

    //! add given vertex normals to grid
    void  add_vtx_normal ( const std::vector< T3Point > &  vertex_normal );

    //////////////////////////////////////
    //
    // function evaluation
    //

    // evaluate function <fn> at index positions on grid
    // and build corresponding vector
    template < typename value_t >
    TScalarVector< value_t > *  eval ( const TBEMFunction< value_t > *  fn )  const;
    
    //////////////////////////////////////
    //
    // misc.
    //

    // return size of grid in bytes
    size_t byte_size () const;

protected:
    //////////////////////////////////////
    //
    // internal routines
    //

    // compute size and normal direction of each triangle
    void comp_triangle_size_normal ();
};

}// namespace

#endif  // __HPRO_TGRID_HH
