//
// Project     : HLIBpro
// File        : TRefinableGrid.cc
// Description : defines a refinable grid class
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/bem/TRefinableGrid.hh"
#include "hpro/io/TGridIO.hh"

namespace Hpro
{
    
using std::unique_ptr;
using std::make_unique;

//
// standard constructor with basic data for refinement
//
TRefinableGrid::TRefinableGrid ( const std::vector< T3Point > &           vertices,
                                 const std::vector< edge_t > &            edges,
                                 const std::vector< triangle_t > &        triangles,
                                 const std::vector< triangle_edges_t > &  triangle_edges )
        : TGrid( vertices, triangles )
        , _edges( edges )
        , _triangle_edges( triangle_edges )
{
    // TODO: test for consistency (size, indices, etc.)
}

//
// standard constructor with basic data for refinement
//
TRefinableGrid::TRefinableGrid ( const std::vector< T3Point > &           vertices,
                                 const std::vector< edge_t > &            edges,
                                 const std::vector< triangle_t > &        triangles,
                                 const std::vector< triangle_edges_t > &  triangle_edges,
                                 const edge_refine_func_t                 refine_func )
        : TGrid( vertices, triangles )
        , _edges( edges )
        , _triangle_edges( triangle_edges )
        , _refine_func( refine_func )
{
    // TODO: test for consistency (size, indices, etc.)
}



TRefinableGrid::~TRefinableGrid ()
{}

//////////////////////////////////////
//
// grid refinement
//

//
// regular refine grid
//
std::unique_ptr< TRefinableGrid >
TRefinableGrid::refine () const
{
    struct  childedges_t {
        idx_t  e0, e1;

        idx_t  child0 ( const bool swap ) { return ( swap ? e1 : e0 ); }
        idx_t  child1 ( const bool swap ) { return ( swap ? e0 : e1 ); }
    };

    const idx_t  UNSET = -1;
    
    //
    // first, loop through all edges, build midpoint and child edges
    //

    const size_t                 nedges     = _edges.size();
    const size_t                 nvertices  = _vertices.size();
    const size_t                 ntriangles = _triangles.size();
    std::vector< T3Point >       child_vertices( nvertices + nedges );        // old vertices + midpoint on each edge
    std::vector< edge_t >        child_edges( 2 * nedges + 3 * ntriangles );  // child edges plus inner triangle edges
    std::vector< childedges_t >  edge_to_children( nedges );

    // copy old vertices
    idx_t  vpos = 0;
    
    while ( vpos < idx_t(nvertices) )
    {
        child_vertices[vpos] = _vertices[vpos];
        vpos++;
    }// while

    // refine edges
    idx_t  epos = 0;
    
    for ( size_t  i = 0; i < nedges; ++i )
    {
        const auto  edge = _edges[ i ];
        const auto  v0   = _vertices[ edge.v0 ];
        const auto  v1   = _vertices[ edge.v1 ];
        T3Point     mid;

        // compute midpoint of edge 
        if ( _refine_func ) mid = _refine_func( v0, v1 );
        else                mid = 0.5 * ( v0 + v1 );

        // edge index -> index of midpoint = nvertices + edge index
        child_vertices[ vpos ] = mid;

        const auto  ce0 = edge_t{ edge.v0, vpos,  UNSET, UNSET };
        const auto  ce1 = edge_t{ vpos, edge.v1,  UNSET, UNSET };

        child_edges[ epos+0 ] = ce0;
        child_edges[ epos+1 ] = ce1;

        edge_to_children[ i ] = childedges_t{ epos, epos+1 };

        vpos++;
        epos += 2;
    }// for

    //
    // now loop over triangles and create child triangles
    // - also create new edges defining inner triangle
    //

    std::vector< triangle_t >        child_triangles( ntriangles * 4 );
    std::vector< triangle_edges_t >  child_triangle_edges( ntriangles * 4 );
    idx_t                            tpos = 0;

    for ( size_t  triangle = 0; triangle < ntriangles; ++triangle )
    {
        const auto  e0 = _triangle_edges[ triangle ].e0;
        const auto  e1 = _triangle_edges[ triangle ].e1;
        const auto  e2 = _triangle_edges[ triangle ].e2;
        const auto  s0 = _triangle_edges[ triangle ].s0;
        const auto  s1 = _triangle_edges[ triangle ].s1;
        const auto  s2 = _triangle_edges[ triangle ].s2;

        //
        // edges for inner triangle
        //

        const idx_t   vm0 = nvertices + e0;
        const idx_t   vm1 = nvertices + e1;
        const idx_t   vm2 = nvertices + e2;
        const edge_t  em02{ vm0, vm2,  tpos+0, tpos+3 };
        const edge_t  em10{ vm1, vm0,  tpos+1, tpos+3 };
        const edge_t  em21{ vm2, vm1,  tpos+2, tpos+3 };

        child_edges[ epos+0 ] = em02;
        child_edges[ epos+1 ] = em10;
        child_edges[ epos+2 ] = em21;
        
        //
        // child triangles
        //

        const idx_t  ce00 = edge_to_children[ e0 ].child0( s0 );
        const idx_t  ce01 = edge_to_children[ e0 ].child1( s0 );
        const idx_t  ce10 = edge_to_children[ e1 ].child0( s1 );
        const idx_t  ce11 = edge_to_children[ e1 ].child1( s1 );
        const idx_t  ce20 = edge_to_children[ e2 ].child0( s2 );
        const idx_t  ce21 = edge_to_children[ e2 ].child1( s2 );
        
        const triangle_edges_t  cte0{ ce21, ce00, epos+0,  s2, s0, false };
        const triangle_edges_t  cte1{ ce01, ce10, epos+1,  s0, s1, false };
        const triangle_edges_t  cte2{ ce11, ce20, epos+2,  s1, s2, false };
        const triangle_edges_t  cte3{ epos, epos+1, epos+2,  true, true, true };
        
        child_triangle_edges[ tpos+0 ] = cte0;
        child_triangle_edges[ tpos+1 ] = cte1;
        child_triangle_edges[ tpos+2 ] = cte2;
        child_triangle_edges[ tpos+3 ] = cte3;

        
        const triangle_t  ct0{ vm2, ( s0 ? _edges[ e0 ].v1 : _edges[ e0 ].v0 ), vm0 };
        const triangle_t  ct1{ vm0, ( s1 ? _edges[ e1 ].v1 : _edges[ e1 ].v0 ), vm1 };
        const triangle_t  ct2{ vm1, ( s2 ? _edges[ e2 ].v1 : _edges[ e2 ].v0 ), vm2 };
        const triangle_t  ct3{ vm0, vm2, vm1 };

        child_triangles[ tpos+0 ] = ct0;
        child_triangles[ tpos+1 ] = ct1;
        child_triangles[ tpos+2 ] = ct2;
        child_triangles[ tpos+3 ] = ct3;

        //
        // update adjoining triangle indices in outer child edges
        //

        if      ( child_edges[ ce00 ].t0 == UNSET ) child_edges[ ce00 ].t0 = tpos+0;
        else if ( child_edges[ ce00 ].t1 == UNSET ) child_edges[ ce00 ].t1 = tpos+0;
        else
            HERROR( ERR_CONSISTENCY, "", "all edge triangles set" );
        
        if      ( child_edges[ ce01 ].t0 == UNSET ) child_edges[ ce01 ].t0 = tpos+1;
        else if ( child_edges[ ce01 ].t1 == UNSET ) child_edges[ ce01 ].t1 = tpos+1;
        else
            HERROR( ERR_CONSISTENCY, "", "all edge triangles set" );
        
        if      ( child_edges[ ce10 ].t0 == UNSET ) child_edges[ ce10 ].t0 = tpos+1;
        else if ( child_edges[ ce10 ].t1 == UNSET ) child_edges[ ce10 ].t1 = tpos+1;
        else
            HERROR( ERR_CONSISTENCY, "", "all edge triangles set" );
        
        if      ( child_edges[ ce11 ].t0 == UNSET ) child_edges[ ce11 ].t0 = tpos+2;
        else if ( child_edges[ ce11 ].t1 == UNSET ) child_edges[ ce11 ].t1 = tpos+2;
        else
            HERROR( ERR_CONSISTENCY, "", "all edge triangles set" );
        
        if      ( child_edges[ ce20 ].t0 == UNSET ) child_edges[ ce20 ].t0 = tpos+2;
        else if ( child_edges[ ce20 ].t1 == UNSET ) child_edges[ ce20 ].t1 = tpos+2;
        else
            HERROR( ERR_CONSISTENCY, "", "all edge triangles set" );
        
        if      ( child_edges[ ce21 ].t0 == UNSET ) child_edges[ ce21 ].t0 = tpos+0;
        else if ( child_edges[ ce21 ].t1 == UNSET ) child_edges[ ce21 ].t1 = tpos+0;
        else
            HERROR( ERR_CONSISTENCY, "", "all edge triangles set" );

        //
        // increase index counters for next triangle
        //
        
        tpos += 4;
        epos += 3;
    }// for

    if ( _refine_func )
        return make_unique< TRefinableGrid >( child_vertices, child_edges, child_triangles, child_triangle_edges, _refine_func );
    else
        return make_unique< TRefinableGrid >( child_vertices, child_edges, child_triangles, child_triangle_edges );
}
    
//////////////////////////////////////
//
// misc.
//

//
// return size of grid in bytes
//
size_t
TRefinableGrid::byte_size () const
{
    return ( TGrid::byte_size() +
             sizeof( edge_t ) * _edges.size() +
             sizeof( triangle_edges_t ) * _triangle_edges.size() );
}

/////////////////////////////////////////////////////////////////////////////////
//
//
//

namespace
{

//
// assuming unit sphere around 0, we compute midpoint and normalize
// to move it onto the surface
//
T3Point
refine_sphere ( T3Point  v0,
                T3Point  v1 )
{
    T3Point  mid = 0.5 * ( v0 + v1 );

    mid.normalise2();

    return mid;
}

}// namespace anonymous

//
// construct grid for a sphere
//
auto
make_sphere () -> std::unique_ptr< TRefinableGrid >
{
    using edge_t           = TRefinableGrid::edge_t;
    using triangle_t       = TGrid::triangle_t;
    using triangle_edges_t = TRefinableGrid::triangle_edges_t;
    
    //
    // construct 6 vertices (base square plus two vertices at top and bottom)
    //
    
    std::vector< T3Point >  vertices( 6 );

    vertices[0] = T3Point( -1, -1,  0 );
    vertices[1] = T3Point(  1, -1,  0 );
    vertices[2] = T3Point(  1,  1,  0 );
    vertices[3] = T3Point( -1,  1,  0 );
    vertices[4] = T3Point(  0,  0, -1 );
    vertices[5] = T3Point(  0,  0,  1 );

    for ( auto &  v : vertices )
        v.normalise2();
    
    //
    // define the 12 edges 
    //

    std::vector< edge_t >  edges( 12 );
    
    edges[0]  = edge_t{ 0, 1,  0, 4 };
    edges[1]  = edge_t{ 2, 1,  1, 5 };
    edges[2]  = edge_t{ 2, 3,  2, 6 };
    edges[3]  = edge_t{ 0, 3,  3, 7 };

    edges[4]  = edge_t{ 1, 4,  0, 1 };
    edges[5]  = edge_t{ 4, 0,  0, 3 };
    edges[6]  = edge_t{ 4, 2,  1, 2 };
    edges[7]  = edge_t{ 3, 4,  2, 3 };

    edges[8]  = edge_t{ 1, 5,  4, 5 };
    edges[9]  = edge_t{ 5, 0,  4, 7 };
    edges[10] = edge_t{ 3, 5,  6, 7 };
    edges[11] = edge_t{ 5, 2,  5, 6 };

    //
    // and finally the 8 triangles
    //

    std::vector< triangle_t >  triangles( 8 );

    triangles[0] = triangle_t{ 0, 4, 1 };
    triangles[1] = triangle_t{ 1, 4, 2 };
    triangles[2] = triangle_t{ 2, 4, 3 };
    triangles[3] = triangle_t{ 3, 4, 0 };

    triangles[4] = triangle_t{ 0, 1, 5 };
    triangles[5] = triangle_t{ 1, 2, 5 };
    triangles[6] = triangle_t{ 2, 3, 5 };
    triangles[7] = triangle_t{ 3, 0, 5 };

    std::vector< triangle_edges_t >   triangle_edges( 8 );
    
    triangle_edges[0] = triangle_edges_t{ 0,  5,  4, true,  true,  true  };
    triangle_edges[1] = triangle_edges_t{ 1,  4,  6, false, false, false };
    triangle_edges[2] = triangle_edges_t{ 2,  6,  7, true,  true,  true  };
    triangle_edges[3] = triangle_edges_t{ 3,  7,  5, false, false, false };
                                                                                     
    triangle_edges[4] = triangle_edges_t{ 0,  8,  9, false, false, false };
    triangle_edges[5] = triangle_edges_t{ 1, 11,  8, true,  true,  true  };
    triangle_edges[6] = triangle_edges_t{ 2, 10, 11, false, false, false };
    triangle_edges[7] = triangle_edges_t{ 3,  9, 10, true,  true,  true  };

    return make_unique< TRefinableGrid >( vertices, edges, triangles, triangle_edges, refine_sphere );
}

//
// construct grid for a sphere (second version)
//
auto
make_sphere2 () -> std::unique_ptr< TRefinableGrid >
{
    using edge_t           = TRefinableGrid::edge_t;
    using triangle_t       = TGrid::triangle_t;
    using triangle_edges_t = TRefinableGrid::triangle_edges_t;
    
    //
    // construct 10 vertices
    //
    
    std::vector< T3Point >  vertices( 10 );

    vertices[0] = T3Point( -1, -1, -0.63 );
    vertices[1] = T3Point(  1, -1, -0.63 );
    vertices[2] = T3Point(  1, -1,  0.63 );
    vertices[3] = T3Point( -1, -1,  0.63 );
    vertices[4] = T3Point(  1,  1, -0.63 );
    vertices[5] = T3Point(  1,  1,  0.63 );
    vertices[6] = T3Point( -1,  1, -0.63 );
    vertices[7] = T3Point( -1,  1,  0.63 );
    
    vertices[8] = T3Point(  0,  0, -1 );
    vertices[9] = T3Point(  0,  0,  1 );

    for ( auto &  v : vertices )
        v.normalise2();
    
    //
    // and 24 edges
    //
    
    std::vector< edge_t >  edges( 24 );
    
    edges[0]  = edge_t{ 0, 1,  -1, -1 };
    edges[1]  = edge_t{ 1, 3,  -1, -1 };
    edges[2]  = edge_t{ 3, 0,  -1, -1 };
    edges[3]  = edge_t{ 2, 1,  -1, -1 };
    edges[4]  = edge_t{ 2, 3,  -1, -1 };
    edges[5]  = edge_t{ 1, 4,  -1, -1 };
    edges[6]  = edge_t{ 4, 2,  -1, -1 };
    edges[7]  = edge_t{ 5, 4,  -1, -1 };
    edges[8]  = edge_t{ 5, 2,  -1, -1 };
    edges[9]  = edge_t{ 4, 6,  -1, -1 };
    edges[10] = edge_t{ 6, 5,  -1, -1 };
    edges[11] = edge_t{ 7, 6,  -1, -1 };
    edges[12] = edge_t{ 7, 5,  -1, -1 };
    edges[13] = edge_t{ 6, 0,  -1, -1 };
    edges[14] = edge_t{ 0, 7,  -1, -1 };
    edges[15] = edge_t{ 3, 7,  -1, -1 };
    edges[16] = edge_t{ 0, 8,  -1, -1 };
    edges[17] = edge_t{ 1, 8,  -1, -1 };
    edges[18] = edge_t{ 4, 8,  -1, -1 };
    edges[19] = edge_t{ 6, 8,  -1, -1 };
    edges[20] = edge_t{ 3, 9,  -1, -1 };
    edges[21] = edge_t{ 2, 9,  -1, -1 };
    edges[22] = edge_t{ 5, 9,  -1, -1 };
    edges[23] = edge_t{ 7, 9,  -1, -1 };

    //
    // and finally the 16 triangles
    //
    
    std::vector< triangle_t >  triangles( 16 );

    triangles[0]  = triangle_t{ 0, 1, 3 };
    triangles[1]  = triangle_t{ 2, 3, 1 };
    triangles[2]  = triangle_t{ 1, 4, 2 };
    triangles[3]  = triangle_t{ 5, 2, 4 };
    triangles[4]  = triangle_t{ 4, 6, 5 };
    triangles[5]  = triangle_t{ 7, 5, 6 };
    triangles[6]  = triangle_t{ 6, 0, 7 };
    triangles[7]  = triangle_t{ 3, 7, 0 };
    
    triangles[8]  = triangle_t{ 8, 1, 0 };
    triangles[9]  = triangle_t{ 8, 4, 1 };
    triangles[10] = triangle_t{ 8, 6, 4 };
    triangles[11] = triangle_t{ 8, 0, 6 };
    triangles[12] = triangle_t{ 9, 3, 2 };
    triangles[13] = triangle_t{ 9, 2, 5 };
    triangles[14] = triangle_t{ 9, 5, 7 };
    triangles[15] = triangle_t{ 9, 7, 3 };
    
    std::vector< triangle_edges_t >   triangle_edges( 16 );

    triangle_edges[0]  = triangle_edges_t{  0,  1,  2, false, false, false };
    triangle_edges[1]  = triangle_edges_t{  4,  1,  3, false, true,  true  };
    triangle_edges[2]  = triangle_edges_t{  5,  6,  3, false, false, false };
    triangle_edges[3]  = triangle_edges_t{  8,  6,  7, false, true,  true  };
    triangle_edges[4]  = triangle_edges_t{  9, 10,  7, false, false, false };
    triangle_edges[5]  = triangle_edges_t{ 12, 10, 11, false, true,  true  };
    triangle_edges[6]  = triangle_edges_t{ 13, 14, 11, false, false, false };
    triangle_edges[7]  = triangle_edges_t{ 15, 14,  2, false, true,  true  };

    triangle_edges[8]  = triangle_edges_t{ 17,  0, 16, true,  true,  false };
    triangle_edges[9]  = triangle_edges_t{ 18,  5, 17, true,  true,  false };
    triangle_edges[10] = triangle_edges_t{ 19,  9, 18, true,  true,  false };
    triangle_edges[11] = triangle_edges_t{ 16, 13, 19, true,  true,  false };
    triangle_edges[12] = triangle_edges_t{ 20,  4, 21, true,  true,  false };
    triangle_edges[13] = triangle_edges_t{ 21,  8, 22, true,  true,  false };
    triangle_edges[14] = triangle_edges_t{ 22, 12, 23, true,  true,  false };
    triangle_edges[15] = triangle_edges_t{ 23, 15, 20, true,  true,  false };
    
    return make_unique< TRefinableGrid >( vertices, edges, triangles, triangle_edges, refine_sphere );
}

//
// construct grid for a cube
//
auto
make_cube   () -> std::unique_ptr< TRefinableGrid >
{
    using edge_t           = TRefinableGrid::edge_t;
    using triangle_t       = TGrid::triangle_t;
    using triangle_edges_t = TRefinableGrid::triangle_edges_t;
    
    //
    // construct 8 vertices
    //
    
    std::vector< T3Point >  vertices( 8 );

    vertices[0] = T3Point( -1, -1, -1 );
    vertices[1] = T3Point(  1, -1, -1 );
    vertices[2] = T3Point(  1,  1, -1 );
    vertices[3] = T3Point( -1,  1, -1 );
    vertices[4] = T3Point( -1, -1,  1 );
    vertices[5] = T3Point(  1, -1,  1 );
    vertices[6] = T3Point(  1,  1,  1 );
    vertices[7] = T3Point( -1,  1,  1 );

    //
    // and 18 edges
    //
    
    std::vector< edge_t >  edges( 18 );
    
    edges[0]  = edge_t{ 0, 1,  -1, -1 };
    edges[1]  = edge_t{ 1, 5,  -1, -1 };
    edges[2]  = edge_t{ 5, 0,  -1, -1 };

    edges[3]  = edge_t{ 1, 2,  -1, -1 };
    edges[4]  = edge_t{ 2, 6,  -1, -1 };
    edges[5]  = edge_t{ 6, 1,  -1, -1 };

    edges[6]  = edge_t{ 2, 3,  -1, -1 };
    edges[7]  = edge_t{ 3, 7,  -1, -1 };
    edges[8]  = edge_t{ 7, 2,  -1, -1 };

    edges[9]  = edge_t{ 3, 0,  -1, -1 };
    edges[10] = edge_t{ 0, 4,  -1, -1 };
    edges[11] = edge_t{ 4, 3,  -1, -1 };

    edges[12] = edge_t{ 5, 4,  -1, -1 };
    edges[13] = edge_t{ 6, 5,  -1, -1 };
    edges[14] = edge_t{ 7, 6,  -1, -1 };
    edges[15] = edge_t{ 4, 7,  -1, -1 };
    edges[16] = edge_t{ 6, 4,  -1, -1 };
    edges[17] = edge_t{ 0, 2,  -1, -1 };

    //
    // and finally the 12 triangles
    //
    
    std::vector< triangle_t >  triangles( 12 );

    triangles[0]  = triangle_t{ 0, 1, 5 };
    triangles[1]  = triangle_t{ 0, 5, 4 };
    
    triangles[2]  = triangle_t{ 1, 2, 6 };
    triangles[3]  = triangle_t{ 1, 6, 5 };
    
    triangles[4]  = triangle_t{ 2, 3, 7 };
    triangles[5]  = triangle_t{ 2, 7, 6 };
    
    triangles[6]  = triangle_t{ 3, 0, 4 };
    triangles[7]  = triangle_t{ 3, 4, 7 };
    
    triangles[8]  = triangle_t{ 4, 5, 6 };
    triangles[9]  = triangle_t{ 4, 6, 7 };
    
    triangles[10] = triangle_t{ 2, 1, 0 };
    triangles[11] = triangle_t{ 2, 0, 3 };
    
    std::vector< triangle_edges_t >   triangle_edges( 12 );

    triangle_edges[0]  = triangle_edges_t{  0,  1,  2, false, false, false };
    triangle_edges[1]  = triangle_edges_t{  2, 12, 10, true,  false, true  };
                                                                    
    triangle_edges[2]  = triangle_edges_t{  3,  4,  5, false, false, false };
    triangle_edges[3]  = triangle_edges_t{  5, 13,  1, true,  false, true  };
                                                                    
    triangle_edges[4]  = triangle_edges_t{  6,  7,  8, false, false, false };
    triangle_edges[5]  = triangle_edges_t{  8, 14,  4, true,  false, true  };
                                                                    
    triangle_edges[6]  = triangle_edges_t{  9, 10, 11, false, false, false };
    triangle_edges[7]  = triangle_edges_t{ 11, 15,  7, true,  false, true  };
                                                                    
    triangle_edges[8]  = triangle_edges_t{ 12, 13, 16, true,  true,  false };
    triangle_edges[9]  = triangle_edges_t{ 16, 14, 15, true,  true,  true  };
                                                                    
    triangle_edges[10] = triangle_edges_t{  3,  0, 17, true,  true,  false };
    triangle_edges[11] = triangle_edges_t{ 17,  9,  6, true,  true,  true  };
    
    return make_unique< TRefinableGrid >( vertices, edges, triangles, triangle_edges );
}

namespace
{

bool
is_near ( const double  f1,
          const double  f2 )
{
    if ( std::abs( f1 - f2 ) < 1e-14 )
        return true;
    else
        return false;
}

//
// assuming cylinder grid, normalise vertex in (x,y)-plane
// and make midpoint same distance from 0 on ende plates
//
T3Point
refine_cylinder ( T3Point  v0,
                  T3Point  v1 )
{
    T3Point  mid = 0.5 * ( v0 + v1 );

    // central axis is not modified
    if (( std::abs( mid.x() ) < 1e-12 ) &&
        ( std::abs( mid.y() ) < 1e-12 ))
        return mid;
        
    const T2Point  c0( v0.x(), v0.y() );
    const T2Point  c1( v1.x(), v1.y() );
    
    if ( is_near( std::abs( v0.z() ), 1.0 ) &&
         is_near( std::abs( v1.z() ), 1.0 ) &&
         is_near( v0.z(), v1.z() ))
    {
        //
        // convert to 2d polar coordinates
        //

        const T2Point  p0( c0.norm2(), std::atan2( c0.y(), c0.x() ) );
        const T2Point  p1( c1.norm2(), std::atan2( c1.y(), c1.x() ) );
        const T2Point  pmid = 0.5 * ( p0 + p1 );
        T2Point        cmid = spherical( pmid.y(), pmid.x() );

        // take care of possible sign change at 0
        if (( Math::sign( c0.x() ) == Math::sign( c1.x()   ) ) &&
            ( Math::sign( c0.x() ) != Math::sign( cmid.x() ) ))
            cmid[0] = -cmid.x();
            
        if (( Math::sign( c0.y() ) == Math::sign( c1.y()   ) ) &&
            ( Math::sign( c0.y() ) != Math::sign( cmid.y() ) ))
            cmid[1] = -cmid.y();
            
        return T3Point( cmid.x(), cmid.y(), mid.z() );
    }// if
    else
    {
        //
        // along cylinder: move vertex to surface by normalisation in xy-plane
        //
        
        T2Point  t2( mid.x(), mid.y() );

        t2.normalise2();
        
        return T3Point( t2.x(), t2.y(), mid.z() );
    }// else
}

}// namespace anonymous


//
// construct grid for a cylinder
//
auto
make_cylinder () -> std::unique_ptr< TRefinableGrid >
{
    auto                    cube  = make_cube();
    std::vector< T3Point >  vertices( cube->vertices() );
    const double            sqrt2 = Math::sqrt( 2.0 );
    
    for ( auto &  v : vertices )
    {
        v[0] = v[0] / sqrt2;
        v[1] = v[1] / sqrt2;
    }// for

    return make_unique< TRefinableGrid >( vertices,
                                          cube->edges(),
                                          cube->triangles(),
                                          cube->triangle_edges(),
                                          refine_cylinder );
}

//
// construct grid for a square
//
auto
make_square  () -> std::unique_ptr< TRefinableGrid >
{
    using edge_t           = TRefinableGrid::edge_t;
    using triangle_t       = TGrid::triangle_t;
    using triangle_edges_t = TRefinableGrid::triangle_edges_t;

    //
    // construct 4 vertices of square
    //
    
    std::vector< T3Point >  vertices( 4 );

    vertices[0] = T3Point( 0, 0, 0 );
    vertices[1] = T3Point( 1, 0, 0 );
    vertices[2] = T3Point( 1, 1, 0 );
    vertices[3] = T3Point( 0, 1, 0 );

    //
    // define the 5 edges 
    //

    std::vector< edge_t >  edges( 5 );
    
    edges[0]  = edge_t{ 0, 1,  0, -1 };
    edges[1]  = edge_t{ 1, 2,  0, -1 };
    edges[2]  = edge_t{ 2, 3,  1, -1 };
    edges[3]  = edge_t{ 3, 0,  1, -1 };
    edges[4]  = edge_t{ 2, 0,  0,  1 };

    //
    // and finally the 2 triangles
    //

    std::vector< triangle_t >  triangles( 2 );

    triangles[0] = triangle_t{ 0, 1, 2 };
    triangles[1] = triangle_t{ 0, 2, 3 };

    std::vector< triangle_edges_t >   triangle_edges( 2 );
    
    triangle_edges[0] = triangle_edges_t{ 0,  1,  4, false, false, false };
    triangle_edges[1] = triangle_edges_t{ 4,  2,  3, true,  false, false };

    return make_unique< TRefinableGrid >( vertices, edges, triangles, triangle_edges );
}

namespace
{

std::unique_ptr< TRefinableGrid >
make_grid ( const std::string &  name,
            const uint           lvl )
{
    if ( ! ( ( name == "sphere"  ) ||
             ( name == "sphere2" ) ||
             ( name == "cube"    ) ||
             ( name == "square"  ) ) )
        HERROR( ERR_ARG, "make_grid", "unknown grid name: " + name );
         
    std::unique_ptr< TRefinableGrid >  grid;
    
    if      ( name == "sphere"  ) grid = make_sphere();
    else if ( name == "sphere2" ) grid = make_sphere2();
    else if ( name == "cube"    ) grid = make_cube();
    else if ( name == "square"  ) grid = make_square();

    for ( uint  i = 0; i < lvl; ++i )
    {
        auto  rgrid = grid->refine();
        
        grid = std::move( rgrid );
    }// while

    return grid;
}

}// namespace anonymous

//
// construct grid for a square
//
std::unique_ptr< TGrid >
make_grid ( const std::string &  gridname )
{
    const auto  dashpos = gridname.find( '-' );

    if ( dashpos != std::string::npos )
    {
        const auto  basename = gridname.substr( 0, dashpos );
        const auto  lvl      = gridname.substr( dashpos+1, gridname.length() );

        if (( basename == "sphere" ) || ( basename == "sphere2" ) || ( basename == "cube" ) || ( basename == "square" ))
            return make_grid( basename, atoi( lvl.c_str() ) );
        else
            return read_grid( gridname );
    }// if
    else
    {
        if (( gridname == "sphere" ) || ( gridname == "sphere2" ) ||
            ( gridname == "cube"   ) || ( gridname == "square"  ))
            return make_grid( gridname, 0 );
        else
            return read_grid( gridname );
    }// else
}

}// namespace Hpro
