#ifndef __HPRO_TREFINABLEGRID_HH
#define __HPRO_TREFINABLEGRID_HH
//
// Project     : HLIBpro
// File        : TRefinableGrid.hh
// Description : defines a refinable grid class
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <functional>

#include "hpro/bem/TGrid.hh"

namespace Hpro
{
    
///////////////////////////////////////////////////////////////
//
// TRefinableGrid: extends normal grid with refinability
//
///////////////////////////////////////////////////////////////

//
// A TRefinableGrid stores the grid using
// - vertex array (as in TGrid)
// - edge array containing vertex indices for the end points
//   plus triangle indices for the adjoining triangles
// - triangle array (as in TGrid)
// - triangle-edge array containing edge indices for the local edges
//   and flags whether the edge direction should be swapped
//
// Indices into triangle and triangle-edge arrays should be consistent, i.e.,
// always referring to the same triangle.
//
// Refinement is always uniform.
//
class TRefinableGrid : public TGrid
{
public:

    struct edge_t {
        idx_t  v0, v1;  // end vertices
        idx_t  t0, t1;  // triangles on each side
    };

    struct triangle_edges_t {
        idx_t  e0, e1, e2;  // edges of triangle
        bool   s0, s1, s2;  // swap status of corresponding edge
    };

    using  edge_refine_func_t = std::function< T3Point ( const T3Point, const T3Point ) >;
    
protected:
    //! edges in triangle grid
    std::vector< edge_t >            _edges;
    
    //! mapping of triangles to corresponding edges
    std::vector< triangle_edges_t >  _triangle_edges;

    //! function object for vertex parametrization
    const edge_refine_func_t         _refine_func;

public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    //! standard constructor with basic data for refinement
    TRefinableGrid ( const std::vector< T3Point > &           vertices,
                     const std::vector< edge_t > &            edges,
                     const std::vector< triangle_t > &        triangles,
                     const std::vector< triangle_edges_t > &  triangle_edges );

    //! constructor with refinement data and parametrization
    TRefinableGrid ( const std::vector< T3Point > &           vertices,
                     const std::vector< edge_t > &            edges,
                     const std::vector< triangle_t > &        triangles,
                     const std::vector< triangle_edges_t > &  triangle_edges,
                     const edge_refine_func_t                 refine_func );

    ~TRefinableGrid ();

    //////////////////////////////////////
    //
    // access to internal data
    //

    // give access to edge array
    auto edges () const -> const std::vector< edge_t > &
    {
        return _edges;
    }

    // give access to triangle edge data
    auto triangle_edges () const -> const std::vector< triangle_edges_t > &
    {
        return _triangle_edges;
    }

    //////////////////////////////////////
    //
    // grid refinement
    //

    //! regular refine grid
    std::unique_ptr< TRefinableGrid >  refine  () const;
    
    //////////////////////////////////////
    //
    // misc.
    //

    // return size of grid in bytes
    size_t  byte_size  () const;
};

//
// construct grid for a sphere
//
auto
make_sphere   () -> std::unique_ptr< TRefinableGrid >;

auto
make_sphere2  () -> std::unique_ptr< TRefinableGrid >;

//
// construct grid for a cube
//
auto
make_cube     () -> std::unique_ptr< TRefinableGrid >;

//
// construct grid for a square
//
auto
make_square   () -> std::unique_ptr< TRefinableGrid >;

//
// either construct grid using internal grid construction
// as <gridname-lvl> or read grid from file <gridname>
//
auto
make_grid     ( const std::string &  gridname ) -> std::unique_ptr< TGrid >;

}// namespace Hpro

#endif  // __HPRO_TGRID_HH
