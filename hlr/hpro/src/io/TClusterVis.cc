//
// Project     : HLIBpro
// File        : TClusterVis.cc
// Description : classes for cluster tree and block cluster tree visualisation
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <fstream>
#include <list>
#include <deque>
#include <sstream>
#include <memory>
#include <vector>

#include "treealg.hh"
#include "baseio.hh"
#include "TPSPrinter.hh"
#include "colour.hh"
#include "TColourMap.hh"

#include "hpro/io/TCoordVis.hh"

#include "hpro/io/TClusterVis.hh"

namespace Hpro
{

using std::unique_ptr;
using std::string;
using std::list;
using std::ofstream;
using std::ostringstream;

namespace
{

//
// colours for PostScript output
//

// border colour
const colour_t  BORDER_COLOUR      = Tango::Black;

// background of block
const colour_t  ADM_BG_COLOUR      = mix( Tango::Chameleon1,  0.25, Tango::White );
const colour_t  INADM_BG_COLOUR    = mix( Tango::ScarletRed1, 0.75, Tango::White );
const colour_t  OMIT_BG_COLOUR     = mix( Tango::SkyBlue1,    0.25, Tango::White );
const colour_t  BLOCKED_BG_COLOUR  = mix( Tango::Black,       0.1,  Tango::White );

}// namespace anonymous

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// class T2DClusterVis::option_t
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// default ctor
//
T2DClusterVis::option_t::option_t ()
        : tree( true )
        , node_procs( false )
        , id( true )
{}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// T2DClusterVis
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// turn on/off drawing of tree
//
T2DClusterVis &
T2DClusterVis::tree    ( const bool  b )
{
    _opt.tree = b;

    return *this;
}
    
//
// turn on/off drawing of per-node processors
//
T2DClusterVis &
T2DClusterVis::node_procs  ( const bool  b )
{
    _opt.node_procs = b;

    return *this;
}

//
// turn on/off drawing of node IDs
//
T2DClusterVis &
T2DClusterVis::id  ( const bool  b )
{
    _opt.id = b;

    return *this;
} 

//
// output of cluster trees
//
void
T2DClusterVis::print ( const TCluster *  c,
                       const string &    name ) const
{
    #define POS(node) (double(node->last() + node->first()) / 2.0)

    if ( c == nullptr )
        HERROR( ERR_ARG, "(T2DClusterVis) write_ct", "argument is nullptr" );
    
    //
    // determine width and depth of tree
    //
    
    list< const TCluster * >  nodes, new_nodes;
    uint                      depth = uint(Tree::depth( c ));
    uint                      width = uint(c->size());
    const double              lvl_shrink = 0.75;
    const double              lvl_fac    = double(width) / double(depth);
    double                    lvl_pos    = 0.0;
    double                    lvl_height = lvl_fac;
    

    for ( uint i = 0; i < depth; i++ )
    {
        lvl_pos    += lvl_height;
        lvl_height *= lvl_shrink;
    }// while

    lvl_pos += 0.35 * lvl_height;
    
    //
    // adjust window
    //
    
    unique_ptr< T2DPrinter >  prn( get_printer( 500, 250, name ) );
    
    prn->begin();
    prn->scale( 500.0 / width, 250.0 / (lvl_pos) );
    prn->set_line_width( 0.5 );

    //
    // draw clustertree
    //

    const TCluster *  node;
    const TCluster *  son;
    int               level = 0;

    //
    // draw rectangles for processors first
    //

    if ( _opt.node_procs )
    {
        new_nodes.push_back( c );

        while ( ! new_nodes.empty() )
        {
            nodes = new_nodes;
            new_nodes.clear();
        
            while ( ! nodes.empty() )
            {
                node = behead( nodes );

                //
                // print coloured rectangle for node
                //

                prn->set_colour( gray( 0 ) );
                prn->draw_rect( double( node->first() ), double( level ) - 0.5,
                               double( node->last() ),  double( level ) + 0.5 );

                //
                // draw connection to sons
                //
            
                for ( uint i = 0; i < node->nsons(); i++ )
                {
                    if ((son = node->son( i )) != nullptr)
                        new_nodes.push_back( son );
                }// for

            }// while

            level++;
        }// while
    }// if
    
    //
    // now the cluster tree
    //

    prn->set_font( "Helvetica", 8 );
    
    lvl_pos    = 0.0;
    lvl_height = lvl_fac;
    
    if ( _opt.tree )
    {
        level = 0;
        new_nodes.push_back( c );

        while ( ! new_nodes.empty() )
        {
            nodes = new_nodes;
            new_nodes.clear();

            while ( ! nodes.empty() )
            {
                node = behead( nodes );

                if ( _opt.id )
                {
                    prn->draw_text( POS(node), lvl_pos, to_string( "%d", node->id() ), JUST_LEFT );
                }// if
                
                //
                // draw connection to sons
                //
            
                for ( uint i = 0; i < node->nsons(); i++ )
                {
                    if ((son = node->son( i )) != nullptr)
                    {
                        prn->draw_line( POS(node), lvl_pos, POS(son), lvl_pos + lvl_height );
                                   
                        new_nodes.push_back( son );
                    }// if
                }// for
            }// while

            level++;

            lvl_pos    += lvl_height;
            lvl_height *= lvl_shrink;
        }// while
    }// if

    //
    // draw indices
    //

    prn->save();
    prn->set_font( "Helvetica", 4 );
    
    uint  mod;

    if      ( c->size() <= 100  ) mod = 10;
    else if ( c->size() <= 1000 ) mod = 25;
    else if ( c->size() <= 5000 ) mod = 50;
    else                          mod = 100;
    
    lvl_pos -= 0.8 * lvl_height;
    
    for ( idx_t  i = c->first(); i <= c->last(); i++ )
    {
        if ( i % mod == 0 )
            prn->draw_text( double(i), lvl_pos, to_string( "%d", i ), JUST_CENTER );
        else
            prn->fill_circle( double(i), lvl_pos - 0.1, 0.1 );
    }// for

    if (( c->last() % mod != 0 ) && ( c->last() % mod > (mod/2) ))
        prn->draw_text( double(c->last()), lvl_pos, to_string( "%d", c->last() ), JUST_CENTER );
        
    prn->restore();
    
    //
    // and finish
    //
    
    prn->end();
}

//
// specialisations
//
T2DPrinter *
TPSClusterVis::get_printer ( const double         width,
                             const double         height,
                             const std::string &  filename ) const
{
    return new TPSPrinter( uint(width), uint(height), add_extension( filename, "eps" ) );
}

/////////////////////////////////////////////////////////////////
//
// VTK visualisation
//
/////////////////////////////////////////////////////////////////

void
TVTKClusterVis::print ( const TCluster *,
                        const std::string & ) const
{
    HERROR( ERR_NOT_IMPL, "", "use coordinate version" );
}

void
TVTKClusterVis::print ( const TCluster *      root,
                        const TCoordinate &   coord,
                        const TPermutation &  perm_i2e,
                        const std::string &   filename ) const
{
    //
    // colorise coordinate of indices per cluster per level
    //

    const size_t              ncoord  = coord.ncoord();
    list< const TCluster * >  clusters;
    uint                      lvl = 0;
    TVTKCoordVis              coord_vis;

    clusters.push_back( root );
    
    while ( ! clusters.empty() )
    {
        uint                      color = 1;
        std::vector< uint >       label( ncoord, 0 );
        list< const TCluster * >  sons;

        for ( auto  cl : clusters )
        {
            for ( auto  idx : *cl )
                label[ perm_i2e( idx ) ] = color;

            for ( uint  i = 0; i < cl->nsons(); ++i )
            {
                if ( cl->son( i ) != nullptr )
                    sons.push_back( cl->son( i ) );
            }// for

            color++;

            if ( color > 7 )
                color = 1; // maximal number of colors!
        }// for
            
        const auto  lvl_name = filename + to_string( "_%d", lvl );

        coord_vis.print( & coord, label, lvl_name );

        clusters.clear();
        clusters = sons;
        lvl++;
    }// while
}

void
print_vtk ( const TCluster *      cl,
            const TCoordinate &   coord,
            const TPermutation &  perm_i2e,
            const std::string &   filename )
{
    TVTKClusterVis  cl_vis;

    cl_vis.print( cl, coord, perm_i2e, filename );
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// class T2DBlockClusterVis::option_t
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// default ctor
//
T2DBlockClusterVis::option_t::option_t ()
        : border( true )
        , background( true )
        , all_nodes( true )
        , loc_is( false )
        , id( false )
        , inner_id( false )
        , procs( false )
        , single_proc( false )
        , pid( 0 )
        , legend( false )
        , max_size_ratio( 10000 )
        , adaptive_lw( false )
        , base_line_width( 1.0 )
{}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// T2DBlockClusterVis
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// options
//
T2DBlockClusterVis &
T2DBlockClusterVis::border ( const bool  b )
{
    _opt.border = b;
    
    return *this;
}
        
T2DBlockClusterVis &
T2DBlockClusterVis::background ( const bool  b )
{
    _opt.background = b;
    
    return *this;
}
        
T2DBlockClusterVis &
T2DBlockClusterVis::all_nodes ( const bool  b )
{
    _opt.all_nodes = b;
    
    return *this;
}
        
T2DBlockClusterVis &
T2DBlockClusterVis::loc_is ( const bool  b )
{
    _opt.loc_is = b;
    
    return *this;
}

T2DBlockClusterVis &
T2DBlockClusterVis::id ( const bool  b )
{
    _opt.id = b;
    
    return *this;
}
        
T2DBlockClusterVis &
T2DBlockClusterVis::inner_id ( const bool  b )
{
    _opt.inner_id = b;
    
    return *this;
}
        
T2DBlockClusterVis &
T2DBlockClusterVis::procs ( const bool  b )
{
    _opt.procs = b;
    
    return *this;
}
        
T2DBlockClusterVis &
T2DBlockClusterVis::single_proc ( const bool  b )
{
    _opt.single_proc = b;
    
    return *this;
}
        
T2DBlockClusterVis &
T2DBlockClusterVis::pid ( const uint  p )
{
    _opt.pid         = p;
    _opt.single_proc = true;
    
    return *this;
}

T2DBlockClusterVis &
T2DBlockClusterVis::legend ( const bool  b )
{
    _opt.legend = b;
    
    return *this;
}

T2DBlockClusterVis &
T2DBlockClusterVis::max_size_ratio ( const double  r )
{
    _opt.max_size_ratio = std::max( r, 0.0 );
    
    return *this;
}

T2DBlockClusterVis &
T2DBlockClusterVis::adaptive_lw  ( const bool  b )
{
    _opt.adaptive_lw = b;
    
    return *this;
}

T2DBlockClusterVis &
T2DBlockClusterVis::base_line_width  ( const double  lw )
{
    _opt.base_line_width = std::max( 0.0, lw );
    
    return *this;
}
    

//
// output of blockcluster trees
//
namespace
{

void
print_rec ( const TBlockCluster *            cluster,
            T2DPrinter *                     prn,
            const uint                       lvl,
            const uint                       max_level,
            const size_t                     max_size,
            const T2DBlockClusterVis::option_t &  opts,
            const TColourMap &               cmap,
            const bool                       show_legend,
            const double                     fn_size )
{
    if ( cluster == nullptr )
        HERROR( ERR_ARG, "print", "argument is nullptr" );
    
    //
    // first recurse
    //
    
    const TProcSet  cl_procs( cluster->procs() );
    const size_t    max_cl_size( std::max( cluster->rowcl()->size(), cluster->colcl()->size() ) );
    bool            is_leaf = cluster->is_leaf();
    const bool      omit    = ( double(max_size) / double(max_cl_size) > opts.max_size_ratio );
    
    // collect sons if ...
    if ( ! is_leaf && ! omit &&                                                             // still large enough and
         ( opts.all_nodes ||                                                                // should print all nodes or
           ( opts.single_proc && cl_procs.is_in( opts.pid ) && ( cl_procs.size() > 1 )) ||  // should print single proc and is local or
           ( opts.procs && (( cl_procs.size() > 1 ) || ( cl_procs == PROCSET_INVALID )))))  // should print proc sets and shared
    {
        for ( uint  i = 0; i < cluster->nsons(); ++i )
        {
            if ( cluster->son(i) != nullptr )
                print_rec( cluster->son(i), prn, lvl+1, max_level, max_size, opts, cmap, show_legend, fn_size );
        }// for
    }// if
    else
    {
        is_leaf = true;
    }// if

    //
    // now print this cluster
    //

    if ( opts.adaptive_lw )
        prn->set_line_width( opts.base_line_width * ( max_level - lvl + 1 ) );
    else
        prn->set_line_width( opts.base_line_width );
        
    
    //
    // print node if leaf or single processor detected (and not all nodes explicitly wanted);
    // if exceeding size ratio, draw at least border
    //
    
    const double  min0 = double(cluster->rowcl()->first());
    const double  max0 = double(cluster->rowcl()->last() + 1);
    const double  min1 = double(cluster->colcl()->first());
    const double  max1 = double(cluster->colcl()->last() + 1);
        
    if ( ! cluster->is_leaf() && omit )
    {
        if ( opts.background )
        {
            prn->set_colour( OMIT_BG_COLOUR );
            prn->fill_rect( min1, min0, max1, max0 );
        }// if
        prn->set_colour( BORDER_COLOUR );
        prn->draw_rect( min1, min0, max1, max0 );
    }// if
    else if (( cluster->is_leaf() || cluster->is_adm() ) || ( ! opts.all_nodes && ( cl_procs.size() == 1 ) ))
    {
        colour_t  bg_col = gray( 255 );
        colour_t  fg_col = gray( 0 );
            
        if ( opts.single_proc )
        {
            if ( cl_procs.is_in( opts.pid ) )
            {
                bg_col = rgb( 192, 0, 0 );
                prn->set_colour( bg_col );
                prn->fill_rect( min1, min0, max1, max0 );
            }// if
        }// if
        else if ( opts.procs )
        {
            if (( cl_procs != PROCSET_INVALID ) && ( cl_procs.size() == 1 ))
                bg_col = cmap.ientry( cl_procs.first() );
                
            if ( bg_col.luminance() < 112 ) fg_col = gray( 255 );
            else                            fg_col = gray( 0 );
                
            prn->set_colour( bg_col );
            prn->fill_rect( min1, min0, max1, max0 );
        }// if
        else
        {
            if ( cluster->is_leaf() && opts.background )
            {
                if ( cluster->is_adm() )
                    bg_col = ADM_BG_COLOUR;
                else
                    bg_col = INADM_BG_COLOUR;
                    
                prn->set_colour( bg_col );
                prn->fill_rect( min1, min0, max1, max0 );
            }// if
        }// else
            
        if ( opts.border )
        {
            prn->set_colour( BORDER_COLOUR );
            prn->draw_rect( min1, min0, max1, max0 );
        }// if

        prn->set_colour( fg_col );
            
        // print processor number
        if ( ! show_legend && ! opts.single_proc && opts.procs && ( cl_procs != PROCSET_INVALID ) && ( cl_procs.size() == 1 ))
        {
            const double  font_size = std::min( std::max( 2.0, double( std::min( cluster->rowcl()->size(),
                                                                                 cluster->colcl()->size() ) ) / 4.0 ),
                                                200.0 );
            prn->save();
            prn->set_font( "Courier-Bold", font_size );
            prn->draw_text( ( max1 + min1 ) / 2.0, ( max0 + min0 + font_size ) / 2.0,
                            to_string( "%d", cl_procs.first() ), JUST_CENTER );
            prn->restore();
        }// if

        // print block index set
        if ( opts.loc_is )
        {
            prn->draw_text( double(max1 + min1) / 2.0,
                            double(max0 + min0) / 2.0 - fn_size,
                            cluster->rowcl()->to_string(), JUST_CENTER );
            prn->draw_text( double(max1 + min1) / 2.0,
                            double(max0 + min0) / 2.0,
                            "x", JUST_CENTER );
            prn->draw_text( double(max1 + min1) / 2.0,
                            double(max0 + min0) / 2.0 + fn_size,
                            cluster->colcl()->to_string(), JUST_CENTER );
        }// if
    }// if
    else // if ( is_leaf )
    {
        //
        // considered as leaf: draw at least border if whished
        //
            
        if ( opts.border )
        {
            prn->set_colour( BORDER_COLOUR );
            prn->draw_rect( min1, min0, max1, max0 );
        }// if
    }// if

    if ( opts.id && ( opts.inner_id || is_leaf ))
    {
        const double  font_size = std::min( std::max( 2.0, double( std::min( cluster->rowcl()->size(),
                                                                             cluster->colcl()->size() ) ) / 4.0 ),
                                            200.0 );

        prn->save();
        prn->set_font( "Helvetica-Bold", font_size );
        prn->draw_text( double(max1 - (font_size / 4.0) ),
                        double(min0) + font_size,
                        to_string( "%d", cluster->id() ),
                        JUST_RIGHT );
        prn->restore();
    }// if
}

}// namespace anonymous

void
T2DBlockClusterVis::print ( const TBlockCluster *  bct,
                            const string &         name ) const
{
    if ( bct == nullptr )
        HERROR( ERR_ARG, "(T2DBlockClusterVis) print", "argument is nullptr" );
    
    uint  max_p = 0;

    // determine maximal processor number
    if ( bct->procs() != PROCSET_INVALID )
        max_p = bct->procs().size();
    
    //
    // go over cluster and draw nodes
    //

    const size_t  nrows        = size_t( bct->rowcl()->size() );
    const size_t  ncols        = size_t( bct->colcl()->size() );
    const size_t  max_size     = std::max( std::max( nrows, ncols ), size_t(1) );
    const size_t  min_size     = std::max( std::min( nrows, ncols ), size_t(1) );
    const double  ratio        = double(min_size) / double(max_size);
    const double  legend_ofs   = 5.0;
    const double  legend_width = 50.0;
    const bool    show_legend  = _opt.legend && ( max_p > 1 ) && ! _opt.single_proc;
    const double  extra_width  = ( show_legend ? legend_ofs + legend_width : 0 );
    const size_t  size_x       = ( ncols == max_size ? 500 : size_t(500.0 * ratio) );
    const size_t  size_y       = ( nrows == max_size ? 500 : size_t(500.0 * ratio) );
    auto          prn          = unique_ptr< T2DPrinter >( get_printer( uint(size_x + extra_width), uint(size_y), name ) );

    //
    // set up viewport so we can use indices as coord
    //

    const double  fn_size = 3.0;
    
    prn->begin();
    prn->save();
    prn->scale( double(size_x) / double(bct->colcl()->size()),
               double(size_y) / double(bct->rowcl()->size()) );
    prn->translate( double(bct->colcl()->first()),
                   double(bct->rowcl()->first()) );

    prn->set_line_width( 0.01 );
    prn->set_font( "Helvetica", fn_size );

    //
    // print block cluster tree
    //

    // TRainbowColourMap  base_cmap( max_p );
    // TShuffleColourMap  cmap( base_cmap, max_p );
    TTangoColourMap    cmap( max_p );
    const uint         max_level = bct->depth();

    print_rec( bct, prn.get(), 1, max_level, max_size, _opt, cmap, show_legend, fn_size );

    prn->restore();
    
    //
    // print legend
    //

    if ( show_legend )
    {
        const double  box_width  = legend_width;
        const double  box_height = double(size_y) / double(max_p);
        const double  lfn_size   = std::min( box_height / 3.0, 20.0 );

        prn->save();
        
        prn->set_font( "Helvetica-Bold", lfn_size );
        prn->set_line_width( 0.1 );
        
        for ( size_t  p = 0; p < max_p; ++p )
        {
            const double  x0 = size_x + legend_ofs;
            const double  y0 = p * box_height;
            const double  x1 = x0 + box_width;
            const double  y1 = y0 + box_height;

            colour_t  bg_col = cmap.ientry( p );
            colour_t  fg_col;
            
            if ( bg_col.luminance() < 96 ) fg_col = gray( 255 );
            else                           fg_col = gray( 0 );
            
            prn->set_colour( bg_col );
            prn->fill_rect( x0, y0, x1, y1 );
            prn->set_colour( gray( 0 ) );
            prn->draw_rect( x0, y0, x1, y1 );

            prn->set_colour( fg_col );
            prn->draw_text( (x0 + x1) / 2.0,
                           (y0 + y1) / 2.0 + lfn_size / 2.0,
                           to_string( "%d", p ), JUST_CENTER );
        }// for

        prn->restore();
    }// if
    
    prn->end();
}


//
// specialisations
//
T2DPrinter *
TPSBlockClusterVis::get_printer ( const double         width,
                                  const double         height,
                                  const std::string &  filename ) const
{
    return new TPSPrinter( uint(width), uint(height), add_extension( filename, "eps" ) );
}

/////////////////////////////////////////////////////////////////
//
// VTK visualisation
//
/////////////////////////////////////////////////////////////////

void
TVTKBlockClusterVis::print ( const TBlockCluster *  bcl,
                             const std::string &    filename ) const
{
    //
    // first check if there is anything to print
    //
    
    if ( bcl == nullptr )
        return;

    //
    // header
    //

    unique_ptr< std::ostream >  out_ptr( open_write( add_extension( filename, "vtk" ) ) );
    std::ostream &              out = * out_ptr.get();

    out << "# vtk DataFile Version 2.0" << std::endl
        << "HLIBpro blockclustertree" << std::endl
        << "ASCII" << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;

    std::list< const TBlockCluster * >  nodes;
    
    //
    // print vertices
    //

    std::deque< T3Point >  vertices;
    const auto             dist  = -1.0 * std::max( bcl->rowcl()->size(), bcl->colcl()->size() ) / 3.0;
    uint                   level = 0;

    nodes.push_back( bcl );

    while ( ! nodes.empty() )
    {
        std::list< const TBlockCluster * >  sons;
        
        for ( auto  node : nodes )
        {
            T3Point  vtx0( node->colcl()->first(), -node->rowcl()->first(), level * dist );
            T3Point  vtx1( node->colcl()->last(),  -node->rowcl()->first(), level * dist );
            T3Point  vtx2( node->colcl()->last(),  -node->rowcl()->last(),  level * dist );
            T3Point  vtx3( node->colcl()->first(), -node->rowcl()->last(),  level * dist );

            vertices.push_back( vtx0 );
            vertices.push_back( vtx1 );
            vertices.push_back( vtx2 );
            vertices.push_back( vtx3 );

            for ( uint  i = 0; i < node->nsons(); ++i )
            {
                if ( node->son( i ) != nullptr )
                    sons.push_back( node->son( i ) );
            }// for
        }// for

        nodes.clear();
        nodes = std::move( sons );
        ++level;
    }// while

    // collect_vertices( bcl, vertices, 0,  );
        
    out << "POINTS " << vertices.size() << " FLOAT" << std::endl;

    for ( auto  vtx : vertices )
    {
        out << vtx.x() << " " << vtx.y() << " " << vtx.z() << std::endl;
    }// for

    //
    // print cells
    //

    const auto  n_nodes = Tree::nnodes( bcl );
    
    out << "CELLS " << n_nodes << " " << 5 * n_nodes << std::endl;

    size_t  vtx_idx = 0;

    nodes.push_back( bcl );

    while ( ! nodes.empty() )
    {
        std::list< const TBlockCluster * >  sons;
        
        for ( auto  node : nodes )
        {
            out << "4 " << vtx_idx << " " << vtx_idx+1 << " " << vtx_idx+2 << " " << vtx_idx+3 << std::endl;

            vtx_idx += 4;

            for ( uint  i = 0; i < node->nsons(); ++i )
            {
                if ( node->son( i ) != nullptr )
                    sons.push_back( node->son( i ) );
            }// for
        }// for

        nodes.clear();
        nodes = std::move( sons );
    }// while
    
    //
    // cell types
    //
    
    out << "CELL_TYPES " << n_nodes << std::endl;
        
    for ( uint  i = 0; i < n_nodes; ++i )
    {
        out << "9 ";
    }// for
    out << std::endl;

    //
    // cell colour
    //

    out << "CELL_DATA " << n_nodes << std::endl
        << "COLOR_SCALARS cellcolour 3" << std::endl;
        
    nodes.push_back( bcl );

    while ( ! nodes.empty() )
    {
        std::list< const TBlockCluster * >  sons;
        
        for ( auto  node : nodes )
        {
            colour_t  col;
            
            if ( node->is_leaf() )
            {
                if ( node->is_adm() ) col = Tango::ScarletRed1;
                else                  col = Tango::Chameleon1;
            }// if
            else
                col = Tango::Aluminium2;

            out << double( col.red   ) / 255.0 << " "
                << double( col.green ) / 255.0 << " "
                << double( col.blue  ) / 255.0 << " "
                << std::endl;

            for ( uint  i = 0; i < node->nsons(); ++i )
            {
                if ( node->son( i ) != nullptr )
                    sons.push_back( node->son( i ) );
            }// for
        }// for

        nodes.clear();
        nodes = std::move( sons );
    }// while
}

void
print_vtk ( const TBlockCluster *  cl,
            const std::string &    filename )
{
    TVTKBlockClusterVis  vis;

    vis.print( cl, filename );
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// GraphViz visualisation
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void
TGVClusterVis::print ( const TCluster * c,
                       const string &   filename ) const
{
    ofstream  os( add_extension( filename, "dot" ).c_str() );
    
    //
    // print block cluster tree in GraphViz format
    //

    os << "digraph G {"  << std::endl
       << "ranksep = 1;" << std::endl;

    os << "node [ label = \"\", style = \"filled\" ];" << std::endl;
    os << "edge [ arrowhead = \"vee\" ];" << std::endl;

    //
    // go through tree and write edge relation
    //

    double                    diam = 1.0;
    size_t                    lvl  = 0;
    list< const TCluster * >  nodes;

    nodes.push_back( c );
    while ( ! nodes.empty() )
    {
        list< const TCluster * >  sons;
        
        while ( ! nodes.empty() )
        {
            const TCluster *  cluster = nodes.front();
            ostringstream     cname;

            cname << cluster << '-' << lvl;

            nodes.pop_front();

            os << "\"" << cname.str()
               << "\" [ width = " << diam << ", height = " << diam << "];" << std::endl;

            for ( uint  j = 0; j < cluster->nsons(); ++j )
            {
                const TCluster *  son = cluster->son( j );
                
                if ( son != nullptr )
                {
                    ostringstream  sname;

                    sname << *son << '-' << lvl+1;
                    
                    os << "\"" << cname.str() << "\"" << " -> " << "\"" << sname.str() << "\"" << ";" << std::endl;
                    
                    sons.push_back( son );
                }// if
            }// for
        }// while

        nodes = sons;
        diam  = std::max( 0.75 * diam, 0.01 );
        lvl++;
    }// while
    
    // finish
    os << "}" << std::endl;
}

void
TGVBlockClusterVis::print ( const TBlockCluster * c,
                            const string        & filename ) const
{
    ofstream  os( add_extension( filename, "dot" ).c_str() );
    
    //
    // print block cluster tree in GraphViz format
    //

    os << "digraph G {"  << std::endl
       << "ranksep = 3;" << std::endl;

    os << "node [ label = \"\", style = \"filled\" ];" << std::endl;
    os << "edge [ arrowhead = \"vee\" ];" << std::endl;
    
    // << "size = \"8,8\";" << std::endl;

    //
    // go through tree and write edge relation
    //

    size_t                         lvl       = 0;
    const size_t                   size_root = std::min( c->rowcl()->size(), c->colcl()->size() );
    list< const TBlockCluster * >  nodes;

    nodes.push_back( c );
    while ( ! nodes.empty() )
    {
        list< const TBlockCluster * >  sons;
        
        while ( ! nodes.empty() )
        {
            const TBlockCluster *  cluster = nodes.front();
            ostringstream          cname;
            string                 name_opts;

            cname << cluster->is() << '-' << lvl;
            
            nodes.pop_front();

            // set additional options for node, e.g. colour
            if ( lvl == 0 )
                name_opts += " root = \"true\", ";
            
            if ( _pid == -1 )
            {
                // set colour according to node type
                if ( cluster->is_leaf() )
                {
                    if ( cluster->is_adm() )
                        name_opts += "color = green, ";
                    else
                        name_opts += "color = red, ";
                }// if
                else
                    name_opts += "color = gray, ";
            }// if
            else
            {
                // set colour according to processor locality
                if ( cluster->procs().is_in( _pid ) )
                    name_opts += "color = blue, ";
                else
                    name_opts += "color = gray, ";
            }// else
            
            // set node size according to ratio between local and global cluster dimensions
            const size_t  size_cl = std::min( cluster->rowcl()->size(), cluster->colcl()->size() );
            const double  diam    = std::min( 2.5, 4.0 * double( size_cl ) / double( size_root ) );
                
            os << "\"" << cname.str()
               << "\" [ " + name_opts + " width = " << diam << ", height = " << diam << " ];"
               << std::endl;

            // collect sons and draw edges
            for ( uint  j = 0; j < cluster->nsons(); ++j )
            {
                const TBlockCluster *  son = cluster->son( j );
                
                if ( son != nullptr )
                {
                    ostringstream  sname;
                    string         ename_opts = " ";

                    sname << son->is() << '-' << lvl+1;
                    
                    if ( _pid != -1 )
                    {
                        if ( cluster->procs().is_in( _pid ) && son->procs().is_in( _pid ) )
                            ename_opts = " [ color = blue ]";
                    }// if
                    
                    os << "\"" << cname.str() << "\"" << " -> " << "\"" << sname.str() << "\"" << ename_opts << ";"
                       << std::endl;
                    
                    sons.push_back( son );
                }// if
            }// for
        }// while

        nodes = sons;
        lvl++;
    }// while
    
    // finish
    os << "}" << std::endl;
}

//
// functional versions
//
void
print_ps ( const TCluster *     cl,
           const std::string &  filename )
{
    TPSClusterVis  vis;

    vis.print( cl, filename );
}

void
print_ps ( const TBlockCluster *  cl,
           const std::string &    filename )
{
    TPSBlockClusterVis  vis;

    vis.print( cl, filename );
}

void
print_gv ( const TCluster *      cl,
           const std::string &   filename )
{
    TGVClusterVis  vis;

    vis.print( cl, filename );
}

void
print_gv ( const TBlockCluster *  cl,
           const std::string &    filename )
{
    TGVBlockClusterVis  vis;

    vis.print( cl, filename );
}


}// namespace Hpro
