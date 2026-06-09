//
// Project     : HLIBpro
// File        : TConstEdgeFnSpace.cc
// Description : implements function space for constant normal linear tangential (CN/LT) edge elements
// Author      : Ronald Kriemann, Jonas Ballani
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/packed.hh"

#include "hpro/bem/TConstEdgeFnSpace.hh"

namespace Hpro
{

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// TConstEdgeFnSpace
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
//
// constructor and destructor
//
TConstEdgeFnSpace::TConstEdgeFnSpace ( const TGrid * agrid )
        : TFnSpace( agrid )
{
    construct();
}

TConstEdgeFnSpace::~TConstEdgeFnSpace ()
{
}

//
// evaluate basis function i on triangle defined by vtxidxs
// for all quadrature points given on the reference triangle,
// note that the value still has to be scaled by +/-0.5 * |e|/|tri|
//
bool
TConstEdgeFnSpace::eval_basis ( const idx_t                     i,
                                const TGrid::triangle_t &       tri,
                                const std::vector< T2Point > &  ref_points,
                                std::vector< T3Point > &        values ) const
{
    // find vertex opposite to edge i
    const edge_t  edge_i   = edge(i);
    idx_t         free_vtx = idx_t( _grid->n_vertices() );
    uint          ncommon  = 0;
    
    for ( uint j = 0; j < 3; j++ )
    {
        if ( ( edge_i.vtx[0] == tri.vtx[j] ) ||
             ( edge_i.vtx[1] == tri.vtx[j] ) )
            ncommon++;
        else
            free_vtx = tri.vtx[j];
    }// for

    if ( ncommon < 2 )
    {
        //
        // triangle does not lie in support
        //

        return false;
    }// if
    else
    {
        const size_t  npts = ref_points.size();

        if ( values.size() < npts )
            values.resize( npts );

        const T3Point  w(  _grid->vertex( free_vtx ) );
        const T3Point  v0( _grid->vertex( tri.vtx[0] ) );
        const T3Point  v1( _grid->vertex( tri.vtx[1] ) );
        const T3Point  v2( _grid->vertex( tri.vtx[2] ) );
        size_t         k = 0;

#if defined(__AVX__)

        using  packed_t = packed< double, ISA_AVX >;
        
        const packed_t  vw[3]  = { w[0],  w[1],  w[2]  };
        const packed_t  vv0[3] = { v0[0], v0[1], v0[2] };
        const packed_t  vv1[3] = { v1[0], v1[1], v1[2] };
        const packed_t  vv2[3] = { v2[0], v2[1], v2[2] };
        const packed_t  vone( double(1) );
        packed_t        vx[3];
        double          t[4];
        
        for ( ; k+3 < npts; k += 4 )
        {
            const packed_t  a1( ref_points[k][0],
                                ref_points[k+1][0],
                                ref_points[k+2][0],
                                ref_points[k+3][0] );
            const packed_t  a2( ref_points[k][1],
                                ref_points[k+1][1],
                                ref_points[k+2][1],
                                ref_points[k+3][1] );
            const packed_t  a0 = sub( sub( vone, a1 ), a2 );
            
            vx[0] = sub( add( add( mul( a0, vv0[0] ),
                                   mul( a1, vv1[0] ) ),
                              mul( a2, vv2[0] ) ),
                         vw[0] );
            vx[1] = sub( add( add( mul( a0, vv0[1] ),
                                   mul( a1, vv1[1] ) ),
                              mul( a2, vv2[1] ) ),
                         vw[1] );
            vx[2] = sub( add( add( mul( a0, vv0[2] ),
                                   mul( a1, vv1[2] ) ),
                              mul( a2, vv2[2] ) ),
                         vw[2] );

            store( vx[0], t );

            values[k][0]   = t[0];
            values[k+1][0] = t[1];
            values[k+2][0] = t[2];
            values[k+3][0] = t[3];

            store( vx[1], t );

            values[k][1]   = t[0];
            values[k+1][1] = t[1];
            values[k+2][1] = t[2];
            values[k+3][1] = t[3];

            store( vx[2], t );

            values[k][2]   = t[0];
            values[k+1][2] = t[1];
            values[k+2][2] = t[2];
            values[k+3][2] = t[3];
        }// for
        
#elif defined(__SSE2__)
        
        using  packed_t = packed< double, ISA_SSE2 >;

        const packed_t  vw[3]  = {  w[0],  w[1],  w[2] };
        const packed_t  vv0[3] = { v0[0], v0[1], v0[2] };
        const packed_t  vv1[3] = { v1[0], v1[1], v1[2] };
        const packed_t  vv2[3] = { v2[0], v2[1], v2[2] };
        const packed_t  vone( double(1) );
        packed_t        vx[3];
        double          t[6];
        
        for ( ; k+1 < npts; k += 2 )
        {
            const packed_t  a1( ref_points[k][0], ref_points[k+1][0] );
            const packed_t  a2( ref_points[k][1], ref_points[k+1][1] );;
            const packed_t  a0( sub( sub( vone, a1 ), a2 ) );

            vx[0] = sub( add( add( mul( a0, vv0[0] ),
                                   mul( a1, vv1[0] ) ),
                              mul( a2, vv2[0] ) ),
                         vw[0] );
            vx[1] = sub( add( add( mul( a0, vv0[1] ),
                                   mul( a1, vv1[1] ) ),
                              mul( a2, vv2[1] ) ),
                         vw[1] );
            vx[2] = sub( add( add( mul( a0, vv0[2] ),
                                   mul( a1, vv1[2] ) ),
                              mul( a2, vv2[2] ) ),
                         vw[2] );

            store( vx[0], t );
            store( vx[1], t+2 );
            store( vx[2], t+4 );

            values[k]   = T3Point( t[0], t[2], t[4] );
            values[k+1] = T3Point( t[1], t[3], t[5] );
        }// for
        
#endif

        //
        // fallback without vectorisation or handle remaining
        //
        
        for ( ; k < npts; k++ )
        {
            const double   a1 = ref_points[k][0];
            const double   a2 = ref_points[k][1];
            const double   a0 = 1.0 - a1 - a2;

            // transform coordinates
            const T3Point  x  = a0 * v0  +  a1 * v1  +  a2 * v2;

            values[k] = x - w;
        }// for
    }// else

    return true;
}

//
// evaluate basis function i on triangle defined by vtxidxs
// for all quadrature points which have been already transformed
// to the given triangle,
// note that the value still has to be scaled by +/-0.5 * |e|/|tri|
//
bool
TConstEdgeFnSpace::eval_basis ( const idx_t                     i,
                                const TGrid::triangle_t &       tri,
                                const std::vector< T3Point > &  quad_points,
                                std::vector< T3Point > &        values ) const
{
    const size_t  npts = quad_points.size();

    if ( values.size() < npts )
        values.resize( npts );

    // find vertex opposite to edge i
    edge_t  edge_i   = edge(i);
    idx_t   free_vtx = idx_t( _grid->n_vertices() );
    uint    ncommon  = 0;

    for ( uint  j = 0; j < 3; j++ )
    {
        if ( ( edge_i.vtx[0] == tri.vtx[j] ) || ( edge_i.vtx[1] == tri.vtx[j] ) )
            ncommon++;
        else
            free_vtx = tri.vtx[j];
    }// for

    // triangle does not lie in support
    if ( ncommon < 2 )
        return false;

    const T3Point  w( _grid->vertex( free_vtx ) );

    for ( size_t  k = 0; k < npts; k++ )
        values[k] = quad_points[k] - w;

    return true;
}

//
// scaling factor for basis function i on triangle tri = +/-0.5 * |e|/|tri|
//
double
TConstEdgeFnSpace::scaling_factor_basis ( const idx_t  i,
                                          const idx_t  tri ) const
{
    const double  edge_length = norm2(   _grid->vertex( edge(i).vtx[0] )
                                     - _grid->vertex( edge(i).vtx[1] ) );
    const double  tri_size    = 0.5 * _grid->tri_size( tri );

    if ( _supp_list[ _supp_list_ptr[i] ] == tri )          // positive triangle
        return double( 0.5 * edge_length / tri_size );
    else if ( _supp_list[ _supp_list_ptr[i] + 1 ] == tri ) // negative triangle
        return double( -0.5 * edge_length / tri_size );

    return 0;
}

//
// transform quadrature points on reference triangle to general triangle
// defined by its vertices
//
void
TConstEdgeFnSpace::ref_points2triangle ( const TGrid::triangle_t &       tri,
                                         const std::vector< T2Point > &  ref_points,
                                         std::vector< T3Point > &        quad_points ) const
{
    const size_t  npts = ref_points.size();

    if ( quad_points.size() < npts )
        quad_points.resize( npts );

    const T3Point  v0( _grid->vertex( tri.vtx[0] ) );
    const T3Point  v1( _grid->vertex( tri.vtx[1] ) );
    const T3Point  v2( _grid->vertex( tri.vtx[2] ) );
    size_t         k = 0;

#if defined(__AVX__)
        
    using  packed_t = packed< double, ISA_AVX >;
    
    const packed_t  vv0[3] = { v0[0], v0[1], v0[2] };
    const packed_t  vv1[3] = { v1[0], v1[1], v1[2] };
    const packed_t  vv2[3] = { v2[0], v2[1], v2[2] };
    const packed_t  vone( double(1) );
    double          t[3][4];
        
    for ( ; k+3 < npts; k += 4 )
    {
        const packed_t  a1( ref_points[k][0],
                            ref_points[k+1][0],
                            ref_points[k+2][0],
                            ref_points[k+3][0] );
        const packed_t  a2( ref_points[k][1],
                            ref_points[k+1][1],
                            ref_points[k+2][1],
                            ref_points[k+3][1] );
        const packed_t  a0( sub( sub( vone, a1 ), a2 ) );
        
        store( add( add( mul( a0, vv0[0] ),
                         mul( a1, vv1[0] ) ),
                    mul( a2, vv2[0] ) ),
               t[0] );
        store( add( add( mul( a0, vv0[1] ),
                         mul( a1, vv1[1] ) ),
                    mul( a2, vv2[1] ) ),
               t[1] );
        store( add( add( mul( a0, vv0[2] ),
                         mul( a1, vv1[2] ) ),
                    mul( a2, vv2[2] ) ),
               t[2] );
        
        quad_points[k]   = T3Point( t[0][0], t[1][0], t[2][0] );
        quad_points[k+1] = T3Point( t[0][1], t[1][1], t[2][1] );
        quad_points[k+2] = T3Point( t[0][2], t[1][2], t[2][2] );
        quad_points[k+3] = T3Point( t[0][3], t[1][3], t[2][3] );
    }// for

#elif defined(__SSE2__)
        
    using  packed_t = packed< double, ISA_SSE2 >;
    
    const packed_t  vv0[3] = { v0[0], v0[1], v0[2] };
    const packed_t  vv1[3] = { v1[0], v1[1], v1[2] };
    const packed_t  vv2[3] = { v2[0], v2[1], v2[2] };
    const packed_t  vone( double(1) );
    double          t[3][2];
        
    for ( ; k+2 < npts; k += 2 )
    {
        const packed_t  a1( ref_points[k][0],
                            ref_points[k+1][0] );
        const packed_t  a2( ref_points[k][1],
                            ref_points[k+1][1] );
        const packed_t  a0( sub( sub( vone, a1 ), a2 ) );
        
        store( add( add( mul( a0, vv0[0] ),
                         mul( a1, vv1[0] ) ),
                    mul( a2, vv2[0] ) ),
               t[0] );
        store( add( add( mul( a0, vv0[1] ),
                         mul( a1, vv1[1] ) ),
                    mul( a2, vv2[1] ) ),
               t[1] );
        store( add( add( mul( a0, vv0[2] ),
                         mul( a1, vv1[2] ) ),
                    mul( a2, vv2[2] ) ),
               t[2] );
        
        quad_points[k]   = T3Point( t[0][0], t[1][0], t[2][0] );
        quad_points[k+1] = T3Point( t[0][1], t[1][1], t[2][1] );
    }// for

#endif

    //
    // fallback and handle remaining
    //
    
    for ( ; k < npts; k++ )
    {
        const double  a1 = ref_points[k][0];
        const double  a2 = ref_points[k][1];
        const double  a0 = 1.0 - a1 - a2;

        // transform coordinates
        quad_points[k] = a0 * v0 + a1 * v1 + a2 * v2;
    }// for
}


////////////////////////////////////////////////////////
//
// private functions
//

//
// construct function space by building indices
// and their support
//
void
TConstEdgeFnSpace::construct ()
{
    const size_t  n_triangles = _grid->n_triangles();
    const size_t  n_vertices  = _grid->n_vertices();
    size_t        n_edges     = 0;

    std::vector< vtx_edge_list_t >  vtx_edge_list( n_vertices );

    vtx_edge_t  vtx_edge;
    idx_t       vid[3];
    idx_t       pos = 0;

    _tri_idx_ptr.resize( n_triangles+1, 0 );
    _tri_idx.resize( 3 * n_triangles, 0 );

    // build mapping of triangle to triangle local indices
    // - three indices per triangle which correspond
    // to the three edges of each triangle
    for ( uint i = 0; i < n_triangles; i++ )
    {
        _tri_idx_ptr[ i ] = pos;

        vid[0] = _grid->triangle(i).vtx[0];
        vid[1] = _grid->triangle(i).vtx[1];
        vid[2] = _grid->triangle(i).vtx[2];

        // iterate along edges of current triangle
        for ( uint j = 0; j < 3; j++ )
        {
            bool   is_in_edge_list = false;
            idx_t  edge_id;

            // check if current edge already possesses an identifier
            for ( const auto  ve : vtx_edge_list[vid[(j+1)%3]] )
            {
                if ( ve.vtx == vid[(j+2)%3] )
                {
                    is_in_edge_list = true;
                    edge_id = ve.edge;
                    break;
                }// if
            }// for

            if ( ! is_in_edge_list )
            {
                // append element (vid[(j+1)%3],n_edges) to vtx_edge_list[ vid[(j+2)%3] ]
                vtx_edge.vtx  = vid[(j+1)%3];
                vtx_edge.edge = uint( n_edges );
                vtx_edge_list[ vid[(j+2)%3] ].push_back( vtx_edge );

                // append element (vid[(j+2) % 3],n_edges) to vtx_edge_list[ vid[(j+1)%3] ]
                vtx_edge.vtx  = vid[(j+2)%3];
                vtx_edge.edge = uint( n_edges );
                vtx_edge_list[ vid[(j+1)%3] ].push_back( vtx_edge );

                // current edge gets new identifier and is added
                // to indices belonging to current triangle
                _tri_idx[ pos ] = idx_t( n_edges );
                n_edges++;
            }//if
            else
            {
                // current edge possesses already an identifier and
                // is added to indices belonging to current triangle
                _tri_idx[ pos ] = edge_id;
            }//else
            
            pos++;
        }//for
    }//for

    _tri_idx_ptr[ n_triangles ] = pos;

    // the support of each edge basis function consists
    // of the two triangles adjacent to that edge
    _supp_list.resize( 2 * n_edges, 0 );
    _supp_list_ptr.resize( n_edges+1, 0 );

    std::vector< bool >  flag( n_edges, false );

    pos = 0;

    for ( size_t  i = 0; i < n_triangles; i++ )
    {
        // iterate along edges of current triangle
        for ( uint j = 0; j < 3; j++ )
        {
            const idx_t  eid = _tri_idx[ pos ];
            
            // check if first triangle has already been assigned
            // to support of current edge
            if ( ! flag[eid] )
            {
                _supp_list[ 2 * eid ] =  idx_t(i);
                flag[eid] = true;
            }// if
            else
            {
                _supp_list[ 2 * eid + 1 ] =  idx_t(i);
            }//else
            pos++;
        }// for
    }// for

    for ( size_t  i = 0; i <= n_edges; i++ )
        _supp_list_ptr[i] = 2 *  idx_t(i);

    // indices are built at the midpoints of edges
    _indices.resize( n_edges );
    _edge_list.resize( n_edges );
    
    idx_t  vtx0[3], vtx1[3];

    for ( uint i = 0; i < n_edges; i++ )
    {
        const idx_t  tri0 = _supp_list[ 2 * i ];
        const idx_t  tri1 = _supp_list[ 2 * i + 1];

        for ( uint j = 0; j < 3; j++ )
        {
            vtx0[j] = _grid->triangle(tri0).vtx[j];
            vtx1[j] = _grid->triangle(tri1).vtx[j];
        }// for

        uint  ncommon = 0;

        for ( uint j = 0; j < 3; j++ )
        for ( uint l = 0; l < 3; l++ )
        {
            if ( vtx0[j] == vtx1[l] )
            {
                std::swap( vtx0[ncommon], vtx0[j] );
                std::swap( vtx1[ncommon], vtx1[l] );
                ncommon++;
                break;
            }// if
        }// for

        _edge_list[i].vtx[0] = vtx0[0];
        _edge_list[i].vtx[1] = vtx0[1];
        _indices[i] = 0.5 * ( _grid->vertex( vtx0[0] ) + _grid->vertex( vtx0[1] ) );
    }// for
}

}// namespace Hpro
