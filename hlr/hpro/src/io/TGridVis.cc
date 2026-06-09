//
// Project     : HLIBpro
// File        : TGridVis.cc
// Description : grid visualisation classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>
#include <algorithm>
#include <fstream>
#include <memory>

#include <boost/filesystem.hpp>

#include "hpro/blas/Algebra.hh"
#include "hpro/cluster/TBBox.hh"
#include "hpro/bem/TConstEdgeFnSpace.hh"

#include "baseio.hh"
#include "TPSPrinter.hh"
#include "TColourMap.hh"

#include "hpro/io/TGridVis.hh"

namespace Hpro
{

using std::unique_ptr;
using std::string;
using std::ofstream;

namespace fs = boost::filesystem;
namespace B  = BLAS;

namespace
{

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// local functions (forward decl.)
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// forward decl. needed by draw_rect
class TIsoProjection;

//
// assign minimal/maximal values to bounding box
//
void
minmax ( const T3Point & u, TBBox & bbox );

//
// print 3D rectangle
//
void
draw_rect ( T2DPrinter & prn, const TIsoProjection & P,
            const T3Point & n,
            const T3Point & u0, const T3Point & u1,
            const T3Point & u2, const T3Point & u3 );

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// local types
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//
// isometric projection for 3d->2d conversion
//
class TIsoProjection
{
private:
    // projection matrix (3x3)
    B::Matrix< double >  _mat;

public:
    /////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TIsoProjection ()
            : _mat( 3, 3 )
    {
        T3Point  dir( 1, 0, 0 );

        set_proj( dir );
    }
    
    TIsoProjection ( const T3Point & dir )
            : _mat( 3, 3 )
    {
        set_proj( dir );
    }

    /////////////////////////////////////////////
    //
    // projection methods
    //

    // define projection by viewing direction
    void set_proj ( const T3Point & dir )
    {
        T3Point  n( dir ), u, v;
        
        // n is normalised direction of view
        n.normalise2();

        // u = (0,1,0)^T x n
        u = cross( T3Point( 0, 1, 0 ), n ).normalise2();

        // v = n x u
        v = cross( n, u ).normalise2();

        _mat(0,0) = u[0]; _mat(0,1) = v[0]; _mat(0,2) = n[0];
        _mat(1,0) = u[1]; _mat(1,1) = v[1]; _mat(1,2) = n[1];
        _mat(2,0) = u[2]; _mat(2,1) = v[2]; _mat(2,2) = n[2];
        
        B::invert( _mat );
    }

    // project vector u to camera space and write result to v
    void project ( const T3Point & u, T3Point & v ) const
    {
        v[0] = _mat(0,0)*u[0] + _mat(0,1)*u[1] + _mat(0,2)*u[2];
        v[1] = _mat(1,0)*u[0] + _mat(1,1)*u[1] + _mat(1,2)*u[2];
        v[2] = _mat(2,0)*u[0] + _mat(2,1)*u[1] + _mat(2,2)*u[2];
    }

    // project bounding box <bbox1> to camera space and write
    // result to <bbox2>
    void project ( const TBBox & bbox1, TBBox & bbox2 ) const
    {
        T3Point  u[8], v[8];
        
        u[0] = T3Point( bbox1.min()[0], bbox1.min()[1], bbox1.min()[2] );
        u[1] = T3Point( bbox1.max()[0], bbox1.min()[1], bbox1.min()[2] );
        u[2] = T3Point( bbox1.max()[0], bbox1.min()[1], bbox1.max()[2] );
        u[3] = T3Point( bbox1.min()[0], bbox1.min()[1], bbox1.max()[2] );
        u[4] = T3Point( bbox1.min()[0], bbox1.max()[1], bbox1.min()[2] );
        u[5] = T3Point( bbox1.max()[0], bbox1.max()[1], bbox1.min()[2] );
        u[6] = T3Point( bbox1.max()[0], bbox1.max()[1], bbox1.max()[2] );
        u[7] = T3Point( bbox1.min()[0], bbox1.max()[1], bbox1.max()[2] );

        for ( uint i = 0; i < 8; i++ )
            project( u[i], v[i] );

        bbox2.min() = TPoint( v[0][0], v[0][1], v[0][2] );
        bbox2.max() = TPoint( v[0][0], v[0][1], v[0][2] );

        for ( uint i = 1; i < 8; i++ )
            minmax( v[i], bbox2 );
    }
};

//
// compare triangles according to depth (wrt. given projection) of midpoint
// - the triangles are referenced by index to triangle-array
//
class TTriCompare
{
protected:
    // data needed to compute triangle midpoints
    const TGrid *           _grid;
    const TIsoProjection *  _proj;

public:
    // constructor
    TTriCompare ( const TGrid *           grid,
                  const TIsoProjection *  proj )
    {
        _grid = grid;
        _proj = proj;
    }
    
    // return true if depth(t1) < depth(t2)
    bool operator () ( const uint  t1,
                       const uint  t2 ) const
    {
        //
        // compare projected depth (in z-coord.) of mid point of triangles
        //

        TGrid::triangle_t  tri1 = _grid->triangle( t1 );
        TGrid::triangle_t  tri2 = _grid->triangle( t2 );
        T3Point            u0, u1;
        T3Point            v0, v1;

        // compute midpoint of triangle 1
        u0  = _grid->vertex( tri1.vtx[0] ) + _grid->vertex( tri1.vtx[1] ) + _grid->vertex( tri1.vtx[2] );
        u0 *= 1.0 / 3.0;
        
        // compute midpoint of triangle 2
        u1  = _grid->vertex( tri2.vtx[0] ) + _grid->vertex( tri2.vtx[1] ) + _grid->vertex( tri2.vtx[2] );
        u1 *= 1.0 / 3.0;
        
        _proj->project( u0, v0 );
        _proj->project( u1, v1 );

        // compare depth, e.g. z-coordinate
        return v0[2] < v1[2];
    }
};

}// namespace anonymous

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// PostScript based grid visualisation
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

TGridVisBase::TGridVisBase ()
        : _min_fn_val(0)
        , _max_fn_val(0)
        , _cmap( "default" )
{}

void
TGridVisBase::set_func_value_interval ( const double  minval,
                                        const double  maxval )
{
    _min_fn_val = minval;
    _max_fn_val = std::max( maxval, minval ); // ensure maxval >= minval
}

void
TGridVisBase::set_colourmap ( const std::string &  cmap )
{
    _cmap = cmap;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// PostScript based grid visualisation
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//
// constructor and destructor
//

T2DGridVis::T2DGridVis ()
        : _view( T3Point( 1, 0, 0 ) )
        , _lighting( false )
        , _draw_bbox( false )
        , _draw_axis( false )
        , _draw_contour( true )
{}

T2DGridVis::T2DGridVis ( const T3Point aview_dir )
        : _view( aview_dir )
        , _lighting( false )
        , _draw_bbox( false )
        , _draw_axis( false )
        , _draw_contour( true )
{
    if ( _view.norm2() < Limits::epsilon<double>() )
        _view = T3Point( 1, 0, 0 );
    else
        _view.normalise2();
}

//////////////////////////////////////
//
// options
//

//
// set viewing direction
//
T2DGridVis &
T2DGridVis::view_dir  ( const T3Point &  view )
{
    _view = view;

    return *this;
}

//
// turn on/off lighting
//
T2DGridVis &
T2DGridVis::lighting ( const bool  b )
{
    _lighting = b;

    return *this;
}
    
//
// print bounding box
//
T2DGridVis &
T2DGridVis::draw_bbox ( const bool  b )
{
    _draw_bbox = b;

    return *this;
}
    
//
// print coord. axis
//
T2DGridVis &
T2DGridVis::draw_axis ( const bool  b )
{
    _draw_axis = b;

    return *this;
}
    
//
// print contour of triangle
//
T2DGridVis &
T2DGridVis::draw_contour ( const bool  b )
{
    _draw_contour = b;

    return *this;
}

///////////////////////////////////////////////////////////////////////////
//
// general printing method for 2D printers
//

namespace
{

//
// print grid with colours according to values in given vector
//
template < typename printer_t,
           typename fnspace_t,
           typename value_t >
void
print ( const T2DGridVis &          vis,
        const TGrid *               grid,
        const fnspace_t *           fnspace,
        const TVector< value_t > *  vec,
        const string &              filename,
        const T3Point               view_dir,
        const bool                  do_lighting,
        const bool                  draw_bbox,
        const bool                  draw_axis,
        const bool                  draw_contour )
{
    static_assert( std::is_floating_point< value_t >::value );
    
    using  real_t = value_t;
    
    //
    // first check if there is anything to print
    //
    
    if (( grid == nullptr ) || ( grid->n_vertices() == 0 ) || ( grid->n_triangles() == 0 ))
        return;

    if (( vec != nullptr ) && (  fnspace == nullptr ))
        HERROR( ERR_ARG, "(T2DGridVis) print", "function space is nullptr" );
    
    const bool  use_vec_value = (( vec != nullptr ) && ( fnspace != nullptr ));

    //
    // determine minimal and maximal values
    //

    real_t  minval = vis.min_func_value();
    real_t  maxval = vis.max_func_value();

    if ( use_vec_value && ( minval != maxval ))
    {
        minval = vec->entry( 0 );
        maxval = minval;
        
        for ( uint i = 1; i < vec->size(); i++ )
        {
            const auto  val = vec->entry( i );
            
            minval = std::min( minval, val );
            maxval = std::max( maxval, val );
        }// for
        
        HDEBUG( to_string( "(T2DGridVis) prints :  min/max = %.3e / %.3e", minval, maxval ) );
    }// if
    
    //
    // define projection:
    //  - isometric projection
    //  - only x and y in new coordinate system are used
    //  - n is the vector from where to view
    //

    // TPoint          n( 1, -1, 1 );
    TIsoProjection  P( view_dir ); 
    T3Point         u, v;

    //
    // determine bounding box
    //

    TBBox  bbox, proj_bbox;

    for ( uint i = 0; i < grid->n_vertices(); i++ )
    {
        u = grid->vertex(i);

        if ( i == 0 )
        {
            bbox.min().assign( 1.0, TPoint( u[0], u[1], u[2] ) );
            bbox.max().assign( 1.0, TPoint( u[0], u[1], u[2] ) );
        }// if
        else
        {
            for ( uint j = 0; j < 3; j++ )
            {
                bbox.min()[j] = std::min( bbox.min()[j], u[j] );
                bbox.max()[j] = std::max( bbox.max()[j], u[j] );
            }// for
        }// else
    }// while

    P.project( bbox, proj_bbox );

    //
    // determine magnification factor such that output
    // roughly fits in 500x500 window
    //
    
    double  len_x = proj_bbox.max()[0] - proj_bbox.min()[0];
    double  len_y = proj_bbox.max()[1] - proj_bbox.min()[1];
    double  mag   = 1.0;

    if ( std::max( len_x, len_y ) * mag > 500 )
    {
        while ( std::max( len_x, len_y ) * mag > 500 )
            mag /= 2.0;
    }// if
    else
    {
        while ( std::max( len_x, len_y ) * mag < 500 )
            mag *= 2.0;
        mag /= 2.0;
    }// else
    
    //
    // create printer with calculated bounding box
    //

    auto  prn = printer_t( uint(mag*len_x), uint(mag*len_y), filename );

    prn.begin();
    prn.scale( 1.0, -1.0 );
    prn.scale( mag, mag );
    prn.translate( -proj_bbox.min()[0], -proj_bbox.max()[1] );

    // prn.set_line_width( std::min( len_x, len_y ) / 250.0 );
    prn.set_line_width( 1.0 / (100.0 * mag) );

    //
    // print back of bounding box
    //

    if ( draw_bbox )
    {
        const T3Point  mn( -view_dir[0], -view_dir[1], -view_dir[2] );
        const T3Point  u0( bbox.min()[0], bbox.min()[1], bbox.min()[2] );
        const T3Point  u1( bbox.max()[0], bbox.min()[1], bbox.min()[2] );
        const T3Point  u2( bbox.max()[0], bbox.min()[1], bbox.max()[2] );
        const T3Point  u3( bbox.min()[0], bbox.min()[1], bbox.max()[2] );
        const T3Point  u4( bbox.min()[0], bbox.max()[1], bbox.min()[2] );
        const T3Point  u5( bbox.max()[0], bbox.max()[1], bbox.min()[2] );
        const T3Point  u6( bbox.max()[0], bbox.max()[1], bbox.max()[2] );
        const T3Point  u7( bbox.min()[0], bbox.max()[1], bbox.max()[2] );

        prn.set_gray( 128 );

        draw_rect( prn, P, mn, u0, u1, u2, u3 ); // front
        draw_rect( prn, P, mn, u1, u5, u6, u2 ); // right
        draw_rect( prn, P, mn, u5, u4, u7, u6 ); // back
        draw_rect( prn, P, mn, u4, u0, u3, u7 ); // left
        draw_rect( prn, P, mn, u3, u2, u6, u7 ); // up
        draw_rect( prn, P, mn, u0, u4, u5, u1 ); // down
    }// if               

    //
    // print coordinate system
    //

    if ( draw_axis )
    {
        double  mind = bbox.max()[0] - bbox.min()[0];

        for ( uint i = 1; i < 3; i++ )
            mind = std::min( mind, bbox.max()[i]-bbox.min()[i] );

        prn.set_font( "Helvetica-Bold", std::min( len_x, len_y ) / 25.0 );
        // prn.set_line_width( std::min( len_x, len_y ) / 100.0 );

        prn.set_rgb( 255, 0, 0 );
        u = T3Point( bbox.min() ) + T3Point( mind, 0, 0 );
        P.project( u, v );
        P.project( T3Point( bbox.min() ), u );
        prn.draw_line( u[0], u[1], v[0], v[1] );

        prn.set_rgb( 0, 192, 0 );
        u = T3Point( bbox.min() ) + T3Point( 0, mind, 0 );
        P.project( u, v );
        P.project( T3Point( bbox.min() ), u );
        prn.draw_line( u[0], u[1], v[0], v[1] );

        prn.set_rgb( 0, 0, 255 );
        u = T3Point( bbox.min() ) + T3Point( 0, 0, mind );
        P.project( u, v );
        P.project( T3Point( bbox.min() ), u );
        prn.draw_line( u[0], u[1], v[0], v[1] );
    }// if
    
    //
    // determine visible triangles and sort them back to front
    //

    const bool    do_backface_culling = false;
    const size_t  n_triangles         = grid->n_triangles();
    auto          vis_triangles       = std::vector< uint >( n_triangles );
    uint          n_vis               = 0;

    if ( do_backface_culling )
    {
        for ( size_t  i = 0; i < n_triangles; i++ )
        {
            const TGrid::triangle_t  tri = grid->triangle( idx_t( i ) );
            T3Point                  t0, t1, t2;
            T3Point                  normal;

            t0 = grid->vertex( tri.vtx[0] );
            t1 = grid->vertex( tri.vtx[1] );
            t2 = grid->vertex( tri.vtx[2] );
        
            //
            // check visibility (backface culling)
            //
        
            u = t1 - t0;
            v = t2 - t0;
        
            normal = cross( u, v ).normalise2();
        
            // if ( dot( normal, _view ) > 0 )
            //     continue;

            vis_triangles[n_vis++] = uint(i);
        }// for
    
        vis_triangles.resize( n_vis );
    }// if
    else
    {
        for ( size_t  i = 0; i < n_triangles; i++ )
            vis_triangles[i] = uint(i);
        n_vis = uint(n_triangles);
    }// else

    TTriCompare  tricmp( grid, & P );
    
    sort( vis_triangles.begin(), vis_triangles.end(), tricmp );
    
    //
    // print coloured triangles
    //

    unique_ptr< TColourMap >  cmap( Hpro::colourmap( vis.colourmap(), 1000 ) );
    
    T3Point        t0, t1, t2;
    T3Point        v0, v1, v2;
    T3Point        normal;
    colour_t       colour   = rgb( 180, 200, 220 );  // base colour of triangles
    colour_t       ambient  = rgb( 20, 20, 20 );     // ambient light colour
    double         diff_fac = 0.95;                  // diffusion coefficient of triangles
    T3Point        light( view_dir );                // position of light (viewing point)

    light *= -10.0;
        
    // prn.set_line_width( std::min( len_x, len_y ) / 100000.0 );
    prn.set_rgb( 0, 0, 0 );

    for ( uint i = 0; i < vis_triangles.size(); i++ )
    {
        const uint               tri_idx = vis_triangles[i];
        const TGrid::triangle_t  tri     = grid->triangle( tri_idx );
        colour_t                 tri_col;

        t0 = grid->vertex( tri.vtx[0] );
        t1 = grid->vertex( tri.vtx[1] );
        t2 = grid->vertex( tri.vtx[2] );

        if ( use_vec_value )
        {
            //
            // determine colour of triangle by evaluating basis functions at
            // triangle center and taking average  (TO DO: for now only rough 
            // approximation due to limited functionality in TFnSpace)
            //
        
            double  val         = 0;
            auto    tri_indices = fnspace->triangle_indices( tri_idx );

            for ( auto  idx : tri_indices )
                val += vec->entry( idx );

            if ( tri_indices.size() > 0 )
                val *= double(1) / double(tri_indices.size());

            tri_col = cmap->fentry( (val - minval) / (maxval - minval) );
        }// if
        else
            tri_col = colour;
        
        //
        // do lighting
        //
            
        if ( do_lighting )
        {
            double    angle = 0;
            colour_t  col1;
            colour_t  col2;
            
            u = t1 - t0;
            v = t2 - t0;

            normal = cross( u, v ).normalise2();

            u = (1.0 / 3.0) * ( t0 + t1 + t2 );
            v = ( light - u ).normalise2();
            
            angle = dot( v, normal );

            col1.red   = uint( std::max( 0.0, (diff_fac * angle * double(tri_col.red))   ) ) + ambient.red;
            col1.green = uint( std::max( 0.0, (diff_fac * angle * double(tri_col.green)) ) ) + ambient.green;
            col1.blue  = uint( std::max( 0.0, (diff_fac * angle * double(tri_col.blue))  ) ) + ambient.blue;

            // do same for other side of two-sided triangles
            u = t2 - t0;
            v = t1 - t0;

            normal = cross( u, v ).normalise2();

            u = (1.0 / 3.0) * ( t0 + t1 + t2 );
            v = ( light - u ).normalise2();
            
            angle = dot( v, normal );

            col2.red   = uint( std::max( 0.0, (diff_fac * angle * double(tri_col.red))   ) ) + ambient.red;
            col2.green = uint( std::max( 0.0, (diff_fac * angle * double(tri_col.green)) ) ) + ambient.green;
            col2.blue  = uint( std::max( 0.0, (diff_fac * angle * double(tri_col.blue))  ) ) + ambient.blue;

            tri_col = rgb( std::max( col1.red,   col2.red   ),
                           std::max( col1.green, col2.green ),
                           std::max( col1.blue,  col2.blue  ) );
        }// if

        P.project( t0, v0 );
        P.project( t1, v1 );
        P.project( t2, v2 );

        prn.set_colour( tri_col );
        prn.fill_triangle( v0[0], v0[1], v1[0], v1[1], v2[0], v2[1] );

        if ( draw_contour )
        {
            prn.set_gray( 0 );
            prn.draw_triangle( v0[0], v0[1], v1[0], v1[1], v2[0], v2[1] );
        }// if
    }// for

    //
    // print front of bounding box
    //

    if ( draw_bbox )
    {
        const T3Point  mn( view_dir[0], view_dir[1], view_dir[2] );
        const T3Point  u0( bbox.min()[0], bbox.min()[1], bbox.min()[2] );
        const T3Point  u1( bbox.max()[0], bbox.min()[1], bbox.min()[2] );
        const T3Point  u2( bbox.max()[0], bbox.min()[1], bbox.max()[2] );
        const T3Point  u3( bbox.min()[0], bbox.min()[1], bbox.max()[2] );
        const T3Point  u4( bbox.min()[0], bbox.max()[1], bbox.min()[2] );
        const T3Point  u5( bbox.max()[0], bbox.max()[1], bbox.min()[2] );
        const T3Point  u6( bbox.max()[0], bbox.max()[1], bbox.max()[2] );
        const T3Point  u7( bbox.min()[0], bbox.max()[1], bbox.max()[2] );

        prn.set_gray( 128 );
        // prn.set_line_width( std::min( len_x, len_y ) / 250.0 );

        draw_rect( prn, P, mn, u0, u1, u2, u3 ); // front
        draw_rect( prn, P, mn, u1, u5, u6, u2 ); // right
        draw_rect( prn, P, mn, u5, u4, u7, u6 ); // back
        draw_rect( prn, P, mn, u4, u0, u3, u7 ); // left
        draw_rect( prn, P, mn, u3, u2, u6, u7 ); // up
        draw_rect( prn, P, mn, u0, u4, u5, u1 ); // down
    }// if               

    prn.end();
}

}// namespace anonymous

//
// Postscript Format
//
void
TPSGridVis::print ( const TGrid *        grid,
                    const std::string &  filename ) const
{
    Hpro::print< TPSPrinter, TConstFnSpace< float >, float >( *this, grid, nullptr, nullptr, filename, _view, _lighting, _draw_bbox, _draw_axis, _draw_contour );
}

template < typename fnspace_t,
           typename value_t >
void TPSGridVis::print ( const TGrid *              grid,
                         const fnspace_t *          fnspace,
                         const TVector< value_t > * vec,
                         const std::string &        filename ) const
{
    Hpro::print< TPSPrinter, fnspace_t, value_t >( *this, grid, fnspace, vec, filename, _view, _lighting, _draw_bbox, _draw_axis, _draw_contour );
}

#define INST_PS( type1, type2 )                                         \
    template                                                            \
    void TPSGridVis::print< type1, type2 > ( const TGrid *            , \
                                             const type1 * ,            \
                                             const TVector< type2 > *  , \
                                             const std::string &      ) const;

INST_PS( TFnSpace< float >, float )
INST_PS( TFnSpace< double >, double )
INST_PS( TConstFnSpace< float >, float )
INST_PS( TConstFnSpace< double >, double )
INST_PS( TLinearFnSpace< float >, float )
INST_PS( TLinearFnSpace< double >, double )
INST_PS( TConstEdgeFnSpace, double )

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// VTK based grid visualisation
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//
// ctor.
//
TVTKGridVis::TVTKGridVis ()
        : _draw_normal( false )
{
}

//
// print normal direction
//
TVTKGridVis &
TVTKGridVis::draw_normal ( const bool  b )
{
    _draw_normal = b;

    return *this;
}
    
//
// print <grid> to file <filename>
//
void
TVTKGridVis::print ( const TGrid *   grid,
                     const string &  filename ) const
{
    //
    // first check if there is anything to print
    //
    
    if (( grid == nullptr ) || ( grid->n_vertices() == 0 ) || ( grid->n_triangles() == 0 ))
        return;

    //
    // header
    //

    unique_ptr< std::ostream >  out_ptr( open_write( add_extension( filename, "vtk" ) ) );
    std::ostream &              out = * out_ptr.get();

    out << "# vtk DataFile Version 2.0" << std::endl
        << "HLIBpro grid" << std::endl
        << "ASCII" << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;

    //
    // print vertices
    //

    out << "POINTS " << grid->n_vertices() << " FLOAT" << std::endl;
        
    for ( uint  i = 0; i < grid->n_vertices(); ++i )
    {
        T3Point  vtx( grid->vertex( i ) );
            
        out << vtx.x() << " " << vtx.y() << " " << vtx.z() << std::endl;
    }// for

    //
    // print triangles
    //

    out << "CELLS " << grid->n_triangles() << " " << 4 * grid->n_triangles() << std::endl;

    for ( uint i = 0; i < grid->n_triangles(); i++ )
    {
        const TGrid::triangle_t  tri = grid->triangle(i);

        out << "3 " << tri.vtx[0] << " " << tri.vtx[1] << " " << tri.vtx[2] << std::endl;
    }// for

    //
    // cell types
    //
    
    out << "CELL_TYPES " << grid->n_triangles() << std::endl;
        
    for ( uint  i = 0; i < grid->n_triangles(); ++i )
    {
        out << "5 ";
    }// for
    out << std::endl;
}

//
// print grid with colours according to values in given vector
//
template < typename fnspace_t,
           typename value_t >
void
TVTKGridVis::print ( const TGrid *              grid,
                     const fnspace_t *          fnspace,
                     const TVector< value_t > * vec,
                     const string &             filename ) const
{
    static_assert( std::is_floating_point< value_t >::value );
    
    using  real_t = value_t;

    //
    // first check if there is anything to print
    //
    
    if (( grid == nullptr ) || ( grid->n_vertices() == 0 ) || ( grid->n_triangles() == 0 ))
        return;

    if (( vec != nullptr ) && (  fnspace == nullptr ))
        HERROR( ERR_ARG, "(TVTKGridVis) print", "function space is nullptr" );
    
    const bool  use_vec_value = (( vec != nullptr ) && ( fnspace != nullptr ));
    
    //
    // determine minimal and maximal values
    //

    real_t  minval = min_func_value();
    real_t  maxval = max_func_value();

    if ( use_vec_value && ( minval == maxval ))
    {
        minval = vec->entry( 0 );
        maxval = minval;
        
        for ( uint i = 1; i < vec->size(); i++ )
        {
            const auto  val = vec->entry( i );
            
            minval = std::min( minval, val );
            maxval = std::max( maxval, val );
        }// for
        
        HDEBUG( to_string( "(TVTKGridVis) prints :  min/max = %.3e / %.3e", minval, maxval ) );
    }// if

    //
    // header
    //

    unique_ptr< std::ostream >  out_ptr( open_write( add_extension( filename, "vtk" ) ) );
    std::ostream &              out = * out_ptr.get();

    out << "# vtk DataFile Version 2.0" << std::endl
        << "HLIBpro grid" << std::endl
        << "ASCII" << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;

    //
    // print vertices
    //

    out << "POINTS " << grid->n_vertices() << " FLOAT" << std::endl;
        
    for ( uint  i = 0; i < grid->n_vertices(); ++i )
    {
        T3Point  vtx( grid->vertex( i ) );
            
        out << vtx.x() << " " << vtx.y() << " " << vtx.z() << std::endl;
    }// for

    //
    // print triangles
    //

    out << "CELLS " << grid->n_triangles() << " " << 4 * grid->n_triangles() << std::endl;

    for ( uint i = 0; i < grid->n_triangles(); i++ )
    {
        const TGrid::triangle_t  tri = grid->triangle(i);

        out << "3 " << tri.vtx[0] << " " << tri.vtx[1] << " " << tri.vtx[2] << std::endl;
    }// for


    //
    // cell types
    //
    
    out << "CELL_TYPES " << grid->n_triangles() << std::endl;
        
    for ( uint  i = 0; i < grid->n_triangles(); ++i )
    {
        out << "5 ";
    }// for
    out << std::endl;

    //
    // cell colour
    //

    out << "CELL_DATA " << grid->n_triangles() << std::endl
        << "COLOR_SCALARS cellcolour 3" << std::endl;
        
    unique_ptr< TColourMap >  cmap( Hpro::colourmap( this->colourmap(), 1000 ) );
    
    for ( uint  i = 0; i < grid->n_triangles(); ++i )
    {
        //
        // determine colour of triangle by evaluating basis functions at
        // triangle center and taking average  (TO DO: for now only rough 
        // approximation due to limited functionality in TFnSpace)
        //
        
        double  val         = 0;
        auto    tri_indices = fnspace->triangle_indices( i );

        for ( auto  idx : tri_indices )
            val += vec->entry( idx );

        if ( tri_indices.size() > 0 )
            val *= double(1) / double(tri_indices.size());

        // clip value
        val = std::max< double >( minval, std::min< double >( maxval, val ) );
        
        const colour_t  col( cmap->fentry( (val - minval) / (maxval - minval) ) );
        
        out << double( col.red   ) / 255.0 << " "
            << double( col.green ) / 255.0 << " "
            << double( col.blue  ) / 255.0 << " "
            << std::endl;
    }// if
}

#define INST_VTK( type1, type2 )                                        \
    template                                                            \
    void TVTKGridVis::print< type1 > ( const TGrid *            ,       \
                                       const type1 * ,                  \
                                       const TVector< type2 > *  ,      \
                                       const string &            ) const;

INST_VTK( TFnSpace< float >, float )
INST_VTK( TFnSpace< double >, double )
INST_VTK( TFnSpace< T2Point >, double )
INST_VTK( TConstFnSpace< float >, float )
INST_VTK( TConstFnSpace< double >, double )
INST_VTK( TLinearFnSpace< float >, float )
INST_VTK( TLinearFnSpace< double >, double )
INST_VTK( TConstEdgeFnSpace, double )

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// local functions (implementation)
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

namespace
{

//
// assign minimal/maximal values to bounding box
//
void
minmax ( const T3Point & u, TBBox & bbox )
{
    for ( uint j = 0; j < 2; j++ )
    {
        bbox.min()[j] = std::min( bbox.min()[j], u[j] );
        bbox.max()[j] = std::max( bbox.max()[j], u[j] );
    }// for
}

//
// print 3D rectangle
//
void
draw_rect ( T2DPrinter &           prn,
            const TIsoProjection & P,
            const T3Point &        n,
            const T3Point &        u0,
            const T3Point &        u1,
            const T3Point &        u2,
            const T3Point &        u3 )
{
    T3Point  w1( u0 ), w2( u0 ); 

    // backface culling
    {
        T3Point normal;
        
        w1 = u1 - u0;
        w2 = u2 - u0;
    
        normal = cross( w1, w2 ).normalise2();
    
        if ( dot( normal, n ) > 0 )
            return;
    }
    
    P.project( u0, w1 );
    P.project( u1, w2 );
    prn.draw_line( w1[0], w1[1], w2[0], w2[1] );

    w1 = w2;
    P.project( u2, w2 );
    prn.draw_line( w1[0], w1[1], w2[0], w2[1] );
    
    w1 = w2;
    P.project( u3, w2 );
    prn.draw_line( w1[0], w1[1], w2[0], w2[1] );
    
    w1 = w2;
    P.project( u0, w2 );
    prn.draw_line( w1[0], w1[1], w2[0], w2[1] );
}

}// namespace anonymous

//
// functional versions
//
void
print_ps ( const TGrid *        grid,
           const std::string &  filename )
{
    TPSGridVis  vis;

    vis.print( grid, filename );
}

void
print_vtk ( const TGrid *        grid,
            const std::string &  filename )
{
    TVTKGridVis  vis;

    vis.print( grid, filename );
}

}// namespace Hpro
