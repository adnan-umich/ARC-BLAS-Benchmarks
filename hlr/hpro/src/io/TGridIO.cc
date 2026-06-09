//
// Project     : HLIBpro
// File        : TGridIO.cc
// Description : contains grid input/output classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <cstdio>
#include <cctype>
#include <vector>
#include <list>
#include <fstream>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "list.hh"
#include "baseio.hh"

#include "hpro/io/TGridIO.hh"

// use preferred Windows variant of sscanf
#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
#define sscanf  sscanf_s
#endif

namespace Hpro
{

namespace fs = boost::filesystem;

using std::string;
using std::list;
using std::unique_ptr;
using std::make_unique;


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// local types
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

namespace
{

struct edge_t { uint vtx[2]; };

}// namespace anonymous

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// TAutoGridIO
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

unique_ptr< TGrid >
TAutoGridIO::read  ( const string & filename ) const
{
    if ( filename == "" )
        HERROR( ERR_ARG, "(TAutoGridIO) read", "empty filename" );

    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TAutoGridIO) read", filename );
    
    unique_ptr< TGridIO >  gio;
        
    switch ( guess_format( filename ) )
    {
        case FMT_HPRO_GRID : gio = make_unique< THproGridIO >(); break;
        case FMT_PLY       : gio = make_unique< TPlyGridIO >(); break;
        case FMT_SURFMESH  : gio = make_unique< TSurfMeshGridIO >(); break;
        case FMT_GMSH      : gio = make_unique< TGMSHGridIO >(); break;

        case FMT_UNKNOWN   :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoGridIO) read", "" );
    }// switch

    return gio->read( filename );
}

void
TAutoGridIO::write  ( const TGrid *   grid,
                      const string &  filename ) const
{
    if ( filename == "" )
        HERROR( ERR_ARG, "(TAutoGridIO) write", "empty filename" );

    unique_ptr< TGridIO >  gio;
        
    switch ( guess_format( filename ) )
    {
        case FMT_HPRO_GRID : gio = make_unique< THproGridIO >(); break;
        case FMT_PLY       : gio = make_unique< TPlyGridIO >(); break;
        case FMT_SURFMESH  : gio = make_unique< TSurfMeshGridIO >(); break;
        case FMT_GMSH      : gio = make_unique< TGMSHGridIO >(); break;

        case FMT_UNKNOWN   :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoGridIO) write", "" );
    }// switch

    gio->write( grid, filename );
}

unique_ptr< TGrid >
read_grid  ( const std::string &  filename )
{
    TAutoGridIO  gio;

    return gio.read( filename );
}

void
write_grid  ( const TGrid *        grid,
              const std::string &  filename )
{
    TAutoGridIO  gio;

    gio.write( grid, filename );
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// THproGridIO
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//
// read in a grid from a file
//

unique_ptr< TGrid >
THproGridIO::read ( const string & filename ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(THproGridIO) read", filename );
    
    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();

    /////////////////////////////////////////////////////////////////
    //
    // read number of vertices, edges and faces
    //

    string  line;
    auto    parts      = std::vector< string >();
    uint    n_vertices = 0;
    uint    n_edges    = 0;
    uint    n_faces    = 0;

    while ( in.good() )
    {
        std::getline( in, line );
        
        // remove comments
        line = line.substr( 0, line.find( "#" ) );
        split( line, " \t", parts );

        if ( parts.size() < 2 )
            continue;

        if ( parts[0] == "nv" )
        {
            if ( n_vertices > 0 )
                HERROR( ERR_GRID_FORMAT, "(THproGridIO) read",
                       "multiple definition of vertex number" );

            n_vertices = str_to_int( parts[1] );
        }// if
        else if ( parts[0] == "ne" )
        {
            if ( n_edges > 0 )
                HERROR( ERR_GRID_FORMAT, "(THproGridIO) read",
                       "multiple definition of edge number" );

            n_edges = str_to_int( parts[1] );
        }// if
        else if ( parts[0] == "nt" )
        {
            if ( n_faces > 0 )
                HERROR( ERR_GRID_FORMAT, "(THproGridIO) read",
                       "multiple definition of triangle number" );

            n_faces = str_to_int( parts[1] );
        }// if
        else
        {
            if (( parts[0] == "v" ) || ( parts[0] == "e" ) || ( parts[0] == "t" ))
                break;
        }// else
    }// while

    if ( n_vertices == 0 )
        HERROR( ERR_GRID_FORMAT, "(THproGridIO) read", "no vertices in grid" );

    // if ( n_edges == 0 )
    //     HERROR( ERR_GRID_FORMAT, "(THproGridIO) read", "no edges in grid" );

    if ( n_faces == 0 )
        HERROR( ERR_GRID_FORMAT, "(THproGridIO) read", "no faces in grid" );
    
    /////////////////////////////////////////////////////////////////
    //
    // read vertices
    //

    std::vector< T3Point >  vertices( n_vertices );

    for ( uint i = 0; i < n_vertices; i++ )
    {
        // remove comments and leading whitespaces
        line = line.substr( 0, line.find( "#" ) );
        while (( line[0] != '\0' ) && isspace( line[0] ) )
            line = line.substr( 1 );

        if ( line.length() == 0 )
        {
            std::getline( in, line );

            --i;
            continue;
        }// if

        if ( line[0] != 'v' )
        {
            if (( line[0] == 'e' ) || ( line[0] == 't' ))
                break;

            HERROR( ERR_GRID_FORMAT, "(THproGridIO) read", "expected vertex" );
        }// if

        //
        // read vertex coordinates
        //

        double  x, y, z;
        int     id;

        sscanf( line.c_str(), "v %d x %lf y %lf z %lf", & id, & x, & y, & z );

        //
        // append vertex
        //

        vertices[id] = T3Point( x, y, z );
        
        std::getline( in, line );
    }// for

    /////////////////////////////////////////////////////////////////
    //
    // (over-) read edges
    //

    for ( uint i = 0; i < n_edges; i++ )
    {
        // remove comments and leading whitespaces
        line = line.substr( 0, line.find( "#" ) );
        while (( line[0] != '\0' ) && isspace( line[0] ) )
            line = line.substr( 1 );

        if ( line.length() == 0 )
        {
            std::getline( in, line );

            --i;
            continue;
        }// if

        if ( line[0] != 'e' )
        {
            if ( line[0] == 't' )
                break;

            HERROR( ERR_GRID_FORMAT, "(THproGridIO) read", "expected edge" );
        }// if

        std::getline( in, line );
    }// for

    /////////////////////////////////////////////////////////////////
    //
    // read triangles
    //

    std::vector< TGrid::triangle_t >  triangles( n_faces );
    
    for ( uint i = 0; i < n_faces; i++ )
    {
        //
        // remove comments and leading whitespaces
        //

        line = line.substr( 0, line.find( "#" ) );
        while (( line[0] != '\0' ) && isspace( line[0] ) )
            line = line.substr( 1 );

        if ( line.length() == 0 )
        {
            std::getline( in, line );

            --i;
            continue;
        }// if

        if ( line[0] != 't' )
            HERROR( ERR_GRID_FORMAT, "(THproGridIO) read", "expected triangle" );

        //
        // read cell depending on dimension
        //

        uint  id;
        uint  vid[3];

        sscanf( line.c_str(), "t %u v %u v %u v %u",
                & id, & vid[0], & vid[1], & vid[2] );

        //
        // check vertices and edges
        //

        for ( uint j = 0; j < 3; j++ )
            if ( vid[j] > n_vertices-1 )
                HERROR( ERR_GRID_DATA, "(THproGridIO) read",
                       "invalid vertex id " + to_string( "%d", vid[j] ) );

        triangles[id].vtx[0] = vid[0];
        triangles[id].vtx[1] = vid[1];
        triangles[id].vtx[2] = vid[2];

        if ( i < n_faces-1 )
            std::getline( in, line );
    }// for
    
    /////////////////////////////////////////////////////////////////
    //
    // finally create grid
    //

    return make_unique< TGrid >( vertices, triangles );
}

void
THproGridIO::write  ( const TGrid *   grid,
                      const string &  filename ) const
{
    if ( grid == nullptr )
        return;

    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();

    //
    // header
    //

    out << "# HLIBpro grid file" << std::endl
        << "# Library version: " << CFG::version() << std::endl
        << "#" << std::endl
        << "# number of vertices" << std::endl
        << "nv " << grid->n_vertices() << std::endl
        << "# number of edges" << std::endl
        << "ne " << 0 << std::endl
        << "# number of triangles" << std::endl
        << "nt " << grid->n_triangles() << std::endl;

    //
    // vertices
    //
    
    out << "#" << std::endl
        << "# vertex list" << std::endl
        << "#" << std::endl;

    for ( idx_t  i = 0; i < idx_t(grid->n_vertices()); ++i )
    {
        const T3Point  vtx( grid->vertex( i ) );

        out << "v " << i
            << " x " << vtx.x()
            << " y " << vtx.y()
            << " z " << vtx.z()
            << std::endl;
    }// for

    //
    // triangles
    //
    
    out << "#" << std::endl
        << "# triangle list" << std::endl
        << "#" << std::endl;

    for ( idx_t  i = 0; i < idx_t(grid->n_triangles()); ++i )
    {
        const TGrid::triangle_t  tri( grid->triangle( i ) );

        out << "t " << i
            << " v " << tri.vtx[0]
            << " v " << tri.vtx[1]
            << " v " << tri.vtx[2]
            << std::endl;
    }// for
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// TPlyGridIO
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

namespace
{

// different ply file formats
enum plyfmt_t  { PLY_ASCII, PLY_BIN_LITTLE, PLY_BIN_BIG };

// different ply element types
enum plyelem_t { PLY_VTX, PLY_FACE, PLY_EDGE };

// different ply property types
enum plyprop_t { PLY_VTX_X, PLY_VTX_Y, PLY_VTX_Z,
                 PLY_VTX_NX, PLY_VTX_NY, PLY_VTX_NZ,
                 PLY_FACE_LIST, PLY_PROP_OTHER } ;

// differenty ply property types
enum plytype_t { PLY_NONE, PLY_INT, PLY_UINT, PLY_FLOAT };

//
// parse given string for ply data type
//
plytype_t
ply_str_to_type ( const string & str )
{
    if (( str == "char" ) || ( str == "int" ) || ( str == "int8" ) ||
        ( str == "int16" ) || ( str == "int32" ))
        return PLY_INT;
    else if (( str == "uchar" ) || ( str == "uint" ) || ( str == "uint8" ) ||
             ( str == "uint16" ) || ( str == "uint32" ))
        return PLY_UINT;
    else if (( str == "float" ) || ( str == "float32" ) || ( str == "float64" ))
        return PLY_FLOAT;
    else
        HERROR( ERR_FMT_PLY, "str_to_plytype", "unknown ply data type \"" + str + "\"" );
}

//
// return coordinate value in string
//
double
ply_coord_prop ( const string & str, const plytype_t type )
{
    if      ( type == PLY_INT   ) return str_to_dbl( str );
    else if ( type == PLY_UINT  ) return str_to_dbl( str );
    else if ( type == PLY_FLOAT ) return str_to_dbl( str );

    return 0.0;
}

}// namespace anonymous

//
// read in a grid from a file
//
unique_ptr< TGrid >
TPlyGridIO::read ( const string & filename ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TPlyGridIO) read", filename );

    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();

    /////////////////////////////////////////////////////////////////
    //
    // parse header of ply file to get information about content
    //

    string             line;
    auto               parts        = std::vector< string >();
    plyfmt_t           format       = PLY_ASCII;
    uint               nvertices    = 0;
    uint               nfaces       = 0;
    uint               nedges       = 0;
    bool               end_header   = false;
    bool               elem_vtx     = false;
    bool               elem_face    = false;
    plytype_t          vtx_coord[3] = { PLY_NONE, PLY_NONE, PLY_NONE };
    plytype_t          vtx_norm[3]  = { PLY_NONE, PLY_NONE, PLY_NONE };
    plytype_t          vtx_count    = PLY_NONE;
    plytype_t          vtx_name     = PLY_NONE;
    list< plyelem_t >  elements;
    list< plyprop_t >  vtx_properties, face_properties;
    
    //
    // first line should be "ply"
    //
    
    std::getline( in, line );
    boost::trim_right( line );
    
    if ( line != "ply" )
        HERROR( ERR_FMT_PLY, "(TPlyGridIO) read", "not a ply file" );

    //
    // parse until end of header was signaled
    //
    
    while ( ! end_header && in.good() )
    {
        std::getline( in, line );
        split( line, " \t", parts );

        if (( parts.size() == 0 ) || ( parts[0] == "comment" ))
            continue;

        // reset element status
        if ( parts[0] != "property" )
        {
            elem_vtx  = false;
            elem_face = false;
        }// if

        //
        // decide upon entry type
        //

        if ( parts[0] == "end_header" )
        {
            end_header = true;
        }// if
        else if ( parts[0] == "format" )
        {
            if      ( parts[1] == "ascii" )
                format = PLY_ASCII;
            else if ( parts[1] == "binary_little_endian" )
                format = PLY_BIN_LITTLE;
            else if ( parts[1] == "binary_big_endian" )
                format = PLY_BIN_BIG;
            else
                HERROR( ERR_FMT_PLY, "(TPlyGridIO) read",
                        "unsupported ply file format \"" + parts[1] + "\"" );
        }// if
        else if ( parts[0] == "element" )
        {
            if ( parts.size() <= 2 )
                HERROR( ERR_FMT_PLY, "(TPlyGridIO) read", "missing cardinality in element entry" );
            
            if ( parts[1] == "vertex" )
            {
                nvertices = str_to_int( parts[2] );
                elem_vtx  = true;
                elements.push_back( PLY_VTX );
            }// if
            else if ( parts[1] == "face" )
            {
                nfaces    = str_to_int( parts[2] );
                elem_face = true;
                elements.push_back( PLY_FACE );
            }// if
            else if ( parts[1] == "edge" )
            {
                nedges = str_to_int( parts[2] );
                elements.push_back( PLY_EDGE );
            }// if
        }// if
        else if ( parts[0] == "property" )
        {
            if ( elem_vtx )
            {
                if ( parts.size() <= 2 )
                    HERROR( ERR_FMT_PLY, "(TPlyGridIO) read",
                            "invalid property format, expected \"property <type> <name>\"" );
                
                if      ( parts[2] == "x"  )
                {
                    vtx_coord[0] = ply_str_to_type( parts[1] );
                    vtx_properties.push_back( PLY_VTX_X );
                }// if
                else if ( parts[2] == "y"  )
                {
                    vtx_coord[1] = ply_str_to_type( parts[1] );
                    vtx_properties.push_back( PLY_VTX_Y );
                }// if
                else if ( parts[2] == "z"  )
                {
                    vtx_coord[2] = ply_str_to_type( parts[1] );
                    vtx_properties.push_back( PLY_VTX_Z );
                }// if
                else if ( parts[2] == "nx" )
                {
                    vtx_norm[0]  = ply_str_to_type( parts[1] );
                    vtx_properties.push_back( PLY_VTX_NX );
                }// if
                else if ( parts[2] == "ny" )
                {
                    vtx_norm[1]  = ply_str_to_type( parts[1] );
                    vtx_properties.push_back( PLY_VTX_NY );
                }// if
                else if ( parts[2] == "nz" )
                {
                    vtx_norm[2]  = ply_str_to_type( parts[1] );
                    vtx_properties.push_back( PLY_VTX_NZ );
                }// if
                else
                    vtx_properties.push_back( PLY_PROP_OTHER );
            }// if
            else if ( elem_face )
            {
                if ( parts.size() <= 4 )
                    HERROR( ERR_FMT_PLY, "(TPlyGridIO) read",
                            "invalid property format, expected \"property list <type1> <type2> vertex_index\"" );
                
                if ( parts[1] == "list"  )
                {
                    vtx_count = ply_str_to_type( parts[2] );

                    if (( vtx_count != PLY_INT ) && ( vtx_count != PLY_UINT ))
                        HERROR( ERR_FMT_PLY, "(TPlyGridIO) read",
                                "unsupported vertex count type; only integer supported" );
                        
                    vtx_name  = ply_str_to_type( parts[3] );

                    if (( vtx_name != PLY_INT ) && ( vtx_name != PLY_UINT ))
                        HERROR( ERR_FMT_PLY, "(TPlyGridIO) read",
                                "unsupported vertex name type; only integer supported" );

                    face_properties.push_back( PLY_FACE_LIST );
                }// if
                else
                    face_properties.push_back( PLY_PROP_OTHER );
            }// if
            else
            {
                // IGNORE
            }// else
        }// if
        else
        {
            // IGNORE
        }// else
    }// while

    if ( format != PLY_ASCII )
        HERROR( ERR_FMT_PLY, "(TPlyGridIO) read", "binary ply file format not yet supported" );

    /////////////////////////////////////////////////////////////////
    //
    // now read actual data according to order defined in header
    //

    std::vector< T3Point >            vertices( nvertices );
    std::vector< T3Point >            vtx_normal( nvertices );
    std::vector< TGrid::triangle_t >  triangles( nfaces );

    while ( elements.size() > 0 )
    {
        const plyelem_t  elem = behead( elements );
        
        if ( elem == PLY_VTX )
        {
            for ( uint i = 0; i < nvertices; i++ )
            {
                // read current line
                std::getline( in, line );
                split( line, " \t", parts );
                
                //
                // read properties per vertex
                //

                uint  j = 0;

                for ( list< plyprop_t >::const_iterator  iter = vtx_properties.begin();
                      iter != vtx_properties.end(); ++iter, ++j )
                {
                    const plyprop_t  prop = *iter;
                    
                    if      ( prop == PLY_VTX_X  ) vertices[i][0]   = ply_coord_prop( parts[j], vtx_coord[0] );
                    else if ( prop == PLY_VTX_Y  ) vertices[i][1]   = ply_coord_prop( parts[j], vtx_coord[1] );
                    else if ( prop == PLY_VTX_Z  ) vertices[i][2]   = ply_coord_prop( parts[j], vtx_coord[2] );
                    else if ( prop == PLY_VTX_NX ) vtx_normal[i][0] = ply_coord_prop( parts[j], vtx_norm[0] );
                    else if ( prop == PLY_VTX_NY ) vtx_normal[i][1] = ply_coord_prop( parts[j], vtx_norm[1] );
                    else if ( prop == PLY_VTX_NZ ) vtx_normal[i][2] = ply_coord_prop( parts[j], vtx_norm[2] );
                }// for
            }// for
        }// if
        else if ( elem == PLY_FACE )
        {
            for ( uint i = 0; i < nfaces; i++ )
            {
                // read current line
                std::getline( in, line );
                split( line, " \t", parts );

                //
                // read properties per face
                //

                uint  j = 0;

                for ( list< plyprop_t >::const_iterator  iter = face_properties.begin();
                      iter != face_properties.end();
                      ++iter, ++j )
                {
                    const plyprop_t  prop = *iter;
                    
                    if ( prop == PLY_FACE_LIST )
                    {
                        if ( int(parts.size()) - int(j) <= 0 )
                            HERROR( ERR_FMT_PLY, "(TPlyGridIO) read", "missing data in vertex list" );
                                    
                        const uint nvtx = str_to_int( parts[j++] );

                        if ( nvtx != 3 )
                            HERROR( ERR_FMT_PLY, "(TPlyGridIO) read", "only triangles are supported" );

                        if ( int(parts.size()) - int(j) < 3 )
                            HERROR( ERR_FMT_PLY, "(TPlyGridIO) read", "missing data in vertex list" );
                        
                        for ( uint k = 0; k < nvtx; k++ )
                            triangles[i].vtx[k] = str_to_int( parts[j++] );

                        // swap triangle to get outward normal direction
                        //swap( (*triangles.get())[i].vtx[1], (*triangles.get())[i].vtx[2] );
                        
                        // decrease since "++j" is in for-header
                        --j;
                    }// if
                }// for
            }// for
        }// if
        else if ( elem == PLY_EDGE )
        {
            HINFO( to_string( "(TPlyGridIO) read : reading %d edges", nedges ) );
        }// if
    }// while

    /////////////////////////////////////////////////////////////////
    //
    // finally create grid
    //

    if (( vtx_norm[0] != PLY_NONE ) &&
        ( vtx_norm[1] != PLY_NONE ) &&
        ( vtx_norm[2] != PLY_NONE ))
        return make_unique< TGrid >( vertices, triangles, vtx_normal );
    else
        return make_unique< TGrid >( vertices, triangles );
}

void
TPlyGridIO::write  ( const TGrid *   grid,
                     const string &  filename ) const
{
    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();

    //
    // write header
    //
    
    out << "ply" << std::endl
        << "format ascii 1.0" << std::endl
        << "comment HLIBpro generated" << std::endl
        << "element vertex " << grid->n_vertices() << std::endl
        << "property float x" << std::endl
        << "property float y" << std::endl
        << "property float z" << std::endl
        << "element face " << grid->n_triangles() << std::endl
        << "property list uchar int vertex_indices" << std::endl
        << "end_header" << std::endl;

    //
    // write vertices
    //

    for ( idx_t  i = 0; i < idx_t(grid->n_vertices()); ++i )
    {
        const T3Point  v( grid->vertex( i ) );

        out << v.x() << " "
            << v.y() << " "
            << v.z() << std::endl;
    }// for

    //
    // write triangles
    //

    for ( idx_t  i = 0; i < idx_t(grid->n_triangles()); ++i )
    {
        const TGrid::triangle_t  tri( grid->triangle( i ) );

        out << "3 "
            << tri.vtx[0] << " "
            << tri.vtx[1] << " "
            << tri.vtx[2] << std::endl;
    }// for
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// TSurMeshGridIO
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//
// read in a grid from a file
//

unique_ptr< TGrid >
TSurfMeshGridIO::read ( const string & filename ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TSurfMeshGridIO) read", filename );
    
    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();

    /////////////////////////////////////////////////////////////////
    //
    // read vertices
    //

    string   line( size_t(1024), '\0' );

    std::getline( in, line );
    boost::trim_right( line );

    if ( line != "surfacemesh" )
        HERROR( ERR_GRID_FORMAT, "(TSurfMeshGridIO) read", "invalid format identifier" );
        
    /////////////////////////////////////////////////////////////////
    //
    // read vertices
    //

    std::vector< string >   parts;
    uint                    n_vertices = 0;
    std::vector< T3Point >  vertices;

    std::getline( in, line );
    n_vertices = str_to_int( line );
    vertices.resize( n_vertices );
    
    for ( uint i = 0; i < n_vertices; i++ )
    {
        std::getline( in, line );
        split( line, " \t", parts );
        
        if ( parts.size() < 3 )
            HERROR( ERR_GRID_FORMAT, "(TSurfMeshGridIO) read", "expected 3 coordinates per vertex" );

        //
        // read vertex coordinates
        //

        double  x, y, z;

        sscanf( line.c_str(), "%lf %lf %lf", & x, & y, & z );

        //
        // append vertex
        //
        
        vertices[i] = T3Point( x, y, z );
    }// for

    /////////////////////////////////////////////////////////////////
    //
    // read faces
    //

    std::vector< TGrid::triangle_t >  triangles;
    uint                              n_faces = 0;
    
    std::getline( in, line );
    n_faces = str_to_int( line );
    triangles.resize( n_faces );
    
    for ( uint i = 0; i < n_faces; i++ )
    {
        std::getline( in, line );
        split( line, " \t", parts );
        
        if ( parts.size() < 3 )
            HERROR( ERR_GRID_FORMAT, "(TSurfMeshGridIO) read", "expected 3 vertex ids per face" );

        //
        // read vertex ids
        //

        TGrid::triangle_t  tri;

        tri.vtx[0] = str_to_int( parts[0] );
        tri.vtx[1] = str_to_int( parts[1] );
        tri.vtx[2] = str_to_int( parts[2] );
        
        for ( uint j = 0; j < 3; j++ )
        {
            if ( tri.vtx[j] == 0 )
                HERROR( ERR_GRID_DATA, "(TSurfMeshGridIO) read", "invalid vertex id 0" );
            if ( tri.vtx[j] > idx_t(n_vertices) )
                HERROR( ERR_GRID_DATA, "(TSurfMeshGridIO) read",
                        "invalid vertex id " + to_string( "%d", tri.vtx[j] ) );

            tri.vtx[j]--;
        }// for
        
        //
        // setup triangle
        //

        triangles[i] = tri;
    }// for
    
    /////////////////////////////////////////////////////////////////
    //
    // finally create grid
    //

    return make_unique< TGrid >( vertices, triangles );
}

void
TSurfMeshGridIO::write  ( const TGrid *,
                          const string & ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// TGMSHGridIO
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//
// read in a grid from a file
//

unique_ptr< TGrid >
TGMSHGridIO::read ( const string & filename ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TGMSHGridIO) read", filename );
    
    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();

    /////////////////////////////////////////////////////////////////
    //
    // file format
    //

    string                 line( size_t(1024), '\0' );
    std::vector< string >  parts;

    std::getline( in, line );
    boost::trim_right( line );

    if ( line != "$MeshFormat" )
        HERROR( ERR_GRID_FORMAT, "(TGMSHGridIO) read", "file not in GMSH format" );
        
    std::getline( in, line );
    split( line, " \t", parts );

    if ( parts.size() != 3 )
        HERROR( ERR_GRID_FORMAT, "(TGMSHGridIO) read", "file not in GMSH format" );
    
    if ( int( str_to_dbl( parts[0] ) ) != 2 )
        HERROR( ERR_GRID_FORMAT, "(TGMSHGridIO) read", "unsupported format version" );

    std::getline( in, line );
    boost::trim_right( line );

    if ( line != "$EndMeshFormat" )
        HERROR( ERR_GRID_FORMAT, "(TGMSHGridIO) read", "file not in GMSH format" );

    /////////////////////////////////////////////////////////////////
    //
    // read vertices
    //

    size_t                  nvertices = 0;
    std::vector< T3Point >  vertices;
    
    std::getline( in, line );
    boost::trim_right( line );

    if ( line != "$Nodes" )
        HERROR( ERR_GRID_FORMAT, "(TGMSHGridIO) read", "file not in GMSH format" );

    std::getline( in, line );

    nvertices = str_to_int( line );
    vertices.resize( nvertices );

    for ( size_t  i = 0; i < nvertices; ++i )
    {
        std::getline( in, line );
        split( line, " \t", parts );
        
        if ( parts.size() < 4 )
            HERROR( ERR_GRID_FORMAT, "(TGMSHGridIO) read", "wrong vertex format" );

        //
        // read vertex coordinates
        //

        int     id;
        double  x, y, z;

        sscanf( line.c_str(), "%d %lf %lf %lf", & id, & x, & y, & z );

        id--;
        if ( id >= int(nvertices) )
            HERROR( ERR_GRID_FORMAT, "(TGMSHGridIO) read", "vertex ID > no. of vertices" );
            
        //
        // append vertex
        //
        
        vertices[ id ] = T3Point( x, y, z );
    }// for

    std::getline( in, line );
    boost::trim_right( line );

    if ( line != "$EndNodes" )
        HERROR( ERR_GRID_FORMAT, "(TGMSHGridIO) read", "file not in GMSH format" );

    /////////////////////////////////////////////////////////////////
    //
    // read triangles
    //

    size_t                            nfaces     = 0;
    size_t                            ntriangles = 0;
    std::vector< TGrid::triangle_t >  faces;
    std::vector< bool >               vtx_used( nvertices, false );
    
    std::getline( in, line );
    boost::trim_right( line );

    if ( line != "$Elements" )
        HERROR( ERR_GRID_FORMAT, "(TGMSHGridIO) read", "file not in GMSH format" );

    std::getline( in, line );

    nfaces = str_to_int( line );
    faces.resize( nfaces );

    for ( size_t  i = 0; i < nfaces; ++i )
    {
        std::getline( in, line );
        split( line, " \t", parts );
        
        // only reading triangle elements: type == 2
        if ( str_to_int( parts[1] ) != 2 )
            continue;

        // triangle has 8 entries
        if ( parts.size() < 8 )
            HERROR( ERR_GRID_FORMAT, "(TGMSHGridIO) read", "wrong element format" );

        // read vertex IDs
        TGrid::triangle_t  tri;

        tri.vtx[0] = str_to_int( parts[5] ) - 1;
        tri.vtx[1] = str_to_int( parts[6] ) - 1;
        tri.vtx[2] = str_to_int( parts[7] ) - 1;
        
        if ( ( tri.vtx[0] >= idx_t( nvertices ) ) ||
             ( tri.vtx[1] >= idx_t( nvertices ) ) ||
             ( tri.vtx[2] >= idx_t( nvertices ) ) )
            HERROR( ERR_GRID_FORMAT, "(TGMSHGridIO) read", "vertex ID > no. of vertices" );

        vtx_used[ tri.vtx[0] ] = true;
        vtx_used[ tri.vtx[1] ] = true;
        vtx_used[ tri.vtx[2] ] = true;
        
        //
        // append triangle
        //
        
        faces[ ntriangles++ ] = tri;
    }// for

    faces.resize( ntriangles );
    
    std::getline( in, line );
    boost::trim_right( line );

    if ( line != "$EndElements" )
        HERROR( ERR_GRID_FORMAT, "(TGMSHGridIO) read", "file not in GMSH format" );

    //
    // check if unused vertices are in grid
    //

    bool  have_unused_vertices = false;
    
    for ( size_t i = 0; i < nvertices; ++i )
    {
        if ( ! vtx_used[i] )
        {
            have_unused_vertices = true;
            break;
        }// if
    }// for

    if ( have_unused_vertices )
    {
        //
        // remove unused vertices
        //
        
        size_t                n_used_vertices = 0;
        std::vector< idx_t >  mapping( nvertices );

        // remove unused vertices
        for ( size_t  i = 0; i < nvertices; ++i )
        {
            if ( vtx_used[i] )
            {
                mapping[i] = idx_t(n_used_vertices);

                if ( i != n_used_vertices )
                    vertices[ n_used_vertices ] = vertices[i];
                
                n_used_vertices++;
            }// if
        }// for

        vertices.resize( n_used_vertices );

        // update indices of triangles
        for ( auto & tri : faces )
        {
            tri.vtx[0] = mapping[ tri.vtx[0] ];
            tri.vtx[1] = mapping[ tri.vtx[1] ];
            tri.vtx[2] = mapping[ tri.vtx[2] ];
        }// for
    }// if
    
    /////////////////////////////////////////////////////////////////
    //
    // finally create grid
    //

    return make_unique< TGrid >( vertices, faces );
}

void
TGMSHGridIO::write  ( const TGrid *   grid,
                      const string &  filename ) const
{
    if ( grid == nullptr )
        return;

    auto    out_ptr = open_write( filename );
    auto &  out     = *out_ptr;

    //
    // header
    //

    out << "$MeshFormat" << std::endl
        << "2.0 0 8" << std::endl
        << "$EndMeshFormat" << std::endl;

    //
    // nodes (id coord0 coord1 coord2)
    //

    out << "$Nodes" << std::endl
        << grid->n_vertices() << std::endl;
    
    for ( uint  i = 0; i < grid->n_vertices(); ++i )
    {
        const auto  coord = grid->vertex( i );
        
        out << i+1 << " " << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
    }// for

    out << "$EndNodes" << std::endl;

    //
    // triangles (type=2 2 2 1 1 vtx0 vtx1 vtx2)
    //

    out << "$Elements" << std::endl
        << grid->n_triangles() << std::endl;
    
    for ( uint  i = 0; i < grid->n_triangles(); ++i )
    {
        const auto  tri = grid->triangle( i );

        out << "2 2 2 1 1 " << tri.vtx[0]+1 << " " << tri.vtx[1]+1 << " " << tri.vtx[2]+1 << std::endl;
    }// for

    out << "$EndElements" << std::endl;
}

}// namespace
