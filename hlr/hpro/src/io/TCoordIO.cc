//
// Project     : HLIBpro
// File        : TCoordIO.cc
// Description : classes for coordinate input/output
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "baseio.hh"

#include "hpro/io/TCoordIO.hh"

namespace Hpro
{

using std::unique_ptr;
using std::make_unique;
using std::string;
using std::vector;

using boost::spirit::qi::long_;
using boost::spirit::qi::double_;
using boost::spirit::qi::_1;
using boost::spirit::qi::phrase_parse;
using boost::spirit::ascii::space;
using boost::phoenix::ref;

namespace fs = boost::filesystem;

///////////////////////////////////////////////////
// 
// TCoordIO
//

void
TCoordIO::write ( const TCoordinate *, const string & ) const
{
    HERROR( ERR_NOT_IMPL, "(TCoordIO) write", "" );
}

unique_ptr< TCoordinate >
TCoordIO::read ( const string & ) const
{
    HERROR( ERR_NOT_IMPL, "(TCoordIO) read", "" );
}

///////////////////////////////////////////////////
// 
// TAutoCoordIO
//

void
TAutoCoordIO::write ( const TCoordinate * coo, const string & filename ) const
{
    unique_ptr< TCoordIO >  cio;
        
    switch ( guess_format( filename ) )
    {
        case FMT_SAMG   : cio = make_unique< TSAMGCoordIO >();   break;
        case FMT_HPRO   : cio = make_unique< THLibCoordIO >();   break;
        case FMT_MATLAB : cio = make_unique< TMatlabCoordIO >(); break;
        case FMT_MTX    : cio = make_unique< TMMCoordIO >();     break;
        case FMT_GMSH   : cio = make_unique< TGMSHCoordIO >();   break;

        case FMT_UNKNOWN :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoCoordIO) read", "" );
            return;
    }// switch

    if ( cio.get() == nullptr )
        return;

    cio->write( coo, filename );
}

unique_ptr< TCoordinate >
TAutoCoordIO::read  ( const string & filename ) const
{
    if ( filename == "" )
        HERROR( ERR_ARG, "(TAutoCoordIO) read", "empty filename" );

    if ( ! boost::filesystem::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TAutoCoordIO) read", filename );
    
    unique_ptr< TCoordIO >  cio;
        
    switch ( guess_format( filename ) )
    {
        case FMT_SAMG   : cio = make_unique< TSAMGCoordIO >();   break;
        case FMT_HPRO   : cio = make_unique< THLibCoordIO >();   break;
        case FMT_MATLAB : cio = make_unique< TMatlabCoordIO >(); break;
        case FMT_MTX    : cio = make_unique< TMMCoordIO >();     break;
        case FMT_GMSH   : cio = make_unique< TGMSHCoordIO >();   break;

        case FMT_UNKNOWN :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoCoordIO) read", "" );
            return nullptr;
    }// switch

    if ( cio.get() == nullptr )
        return nullptr;

    return cio->read( filename );
}

unique_ptr< TCoordinate >
read_coord  ( const char *  filename )
{
    TAutoCoordIO  cio;

    return cio.read( filename );
}

void
write_coord  ( const TCoordinate *  coord,
               const char *         filename )
{
    TAutoCoordIO  cio;

    cio.write( coord, filename );
}

///////////////////////////////////////////////////
// 
// TSAMGCoordIO
//

void
TSAMGCoordIO::write ( const TCoordinate * coo, const string & filename ) const
{
    if ( coo == nullptr )
        HERROR( ERR_ARG, "(TSAMGCoordIO) write", "coordinates are nullptr" );
    
    //
    // read coord file
    //

    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();

    const uint    dim    = coo->dim();
    const size_t  nnodes = coo->ncoord();
    string        line;

    out << nnodes << " " << dim << std::endl;

    for ( size_t  i = 0; i < nnodes; i++ )
    {
        const double * v = coo->coord( idx_t(i) );

        for ( uint  j = 0; j < dim; ++j )
            out << boost::format( "%g" ) % v[j];
        out << std::endl;
    }// for
}

unique_ptr< TCoordinate >
TSAMGCoordIO::read  ( const string & filename ) const
{
    //
    // read coord file
    //

    if ( ! boost::filesystem::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TSAMGCoordIO) read", filename );
    
    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();

    uint                dim    = 0;
    ulong               nnodes = 0;
    vector< double * >  vertices;
    string              line;
    vector< string >    parts;

    std::getline( in, line );
    split( line, " \t\n", parts );

    if ( parts.size() < 2 )
        HERROR( ERR_FMT_SAMG, "(TSAMGCoordIO) read", "expected <nvertices> <dim>" );
    
    nnodes = str_to_int( parts[0] );
    dim    = str_to_int( parts[1] );

    vertices.resize( nnodes );

    for ( uint i = 0; i < nnodes; i++ )
    {
        double * v = new double[ dim ];
        
        std::getline( in, line );
        split( line, " \t\n", parts );

        if ( parts.size() < dim )
            HERROR( ERR_FMT_SAMG, "(TSAMGCoordIO) read",
                    to_string( "expected %d floats per coordinate", dim ) );
        
        for ( uint  j = 0; j < dim; ++j )
        {
            v[j] = str_to_dbl( parts[j] );
        }// for

        vertices[i] = v;
    }// for

    return make_unique< TCoordinate >( vertices, dim );
}

///////////////////////////////////////////////////
// 
// TMMCoordIO
//

namespace
{

enum mtxformat_t { MTX_ARRAY, MTX_COORD };
enum mtxfield_t  { MTX_REAL, MTX_INTEGER, MTX_COMPLEX, MTX_PATTERN };
enum mtxsym_t    { MTX_GENERAL, MTX_SYM, MTX_SKEWSYM, MTX_HERM };

}// namespace anonymous

void
TMMCoordIO::write ( const TCoordinate *, const string & ) const
{
    HERROR( ERR_NOT_IMPL, "(TMMCoordIO) write", "" );
}

unique_ptr< TCoordinate >
TMMCoordIO::read ( const string & filename ) const
{
    //
    // read coordinate file
    //

    if ( ! boost::filesystem::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TMMCoordIO) read", filename );
    
    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();
    
    ///////////////////////////////////////////////////
    //
    // read header and determine format of file
    //

    string            line;
    vector< string >  parts;
    mtxformat_t       mat_format = MTX_ARRAY;
    mtxfield_t        mat_field  = MTX_REAL;
    mtxsym_t          mat_sym    = MTX_GENERAL;

    while ( in.good() )
    {
        std::getline( in, line );
        split( line, " \t\r\n", parts );

        if ( parts[0] == "%%MatrixMarket" )
        {
            if ( parts.size() < 3 )
                HERROR( ERR_FMT_MTX, "(TMMCoordIO) read",
                        "missing header data, expected <object> <format>" );
            
            if ( parts[1] != "matrix" )
                HERROR( ERR_FMT_MTX, "(TMMCoordIO) read", "file does not contain matrix" );
            
            if      ( parts[2] == "coordinate" ) mat_format = MTX_COORD;
            else if ( parts[2] == "array" )      mat_format = MTX_ARRAY;
            else
                HERROR( ERR_FMT_MTX, "(TMMCoordIO) read", "unknown matrix format \"" + parts[2] + "\"" );

            // look for qualifiers
            for ( uint i = 3; i < parts.size(); i++ )
            {
                if      ( parts[i] == "real"           ) mat_field = MTX_REAL;
                else if ( parts[i] == "integer"        ) mat_field = MTX_INTEGER;
                else if ( parts[i] == "complex"        ) mat_field = MTX_COMPLEX;
                else if ( parts[i] == "pattern"        ) mat_field = MTX_PATTERN;
                else if ( parts[i] == "general"        ) mat_sym   = MTX_GENERAL;
                else if ( parts[i] == "symmetric"      ) mat_sym   = MTX_SYM;
                else if ( parts[i] == "skew-symmetric" ) mat_sym   = MTX_SKEWSYM;
                else if ( parts[i] == "hermitian"      ) mat_sym   = MTX_HERM;
                else
                    HERROR( ERR_FMT_MTX, "(TMMCoordIO) read", "unknown qualifier \"" + parts[i] + "\"" );
            }// for
            
            break;
        }// if
    }// while

    if ( mat_sym != MTX_GENERAL )
        HERROR( ERR_NOT_IMPL, "(TMMCoordIO) read", "only unsymmetric matrices supported" );

    if ( mat_format == MTX_COORD )
        HERROR( ERR_NOT_IMPL, "(TMMCoordIO) read", "only \"array\" format supported" );

    if ( mat_field != MTX_REAL )
        HERROR( ERR_NOT_IMPL, "(TMMCoordIO) read", "only \"real\" matrices supported" );
        
    ///////////////////////////////////////////////////
    //
    // read coordinate matrix
    //
    
    if ( mat_format == MTX_ARRAY )
    {
        //
        // read size of vector
        //

        uint  nrows = 0, ncols = 0;

        while ( in.good() )
        {
            std::getline( in, line );

            if ( line[0] == '%' )
                continue;

            const bool  r = phrase_parse( line.begin(), line.end(),
                                          ( long_[ ref(nrows) = _1 ] >>
                                            long_[ ref(ncols) = _1 ] ),
                                          space );

            if ( ! r )
                HERROR( ERR_FMT_MTX, "(TMMCoordIO) read", "missing data, expected <nrows> <ncols>" );

            break;
        }// while

        HINFO( to_string( "(TMMCoordIO) read : reading matrix of size %dx%d", nrows, ncols ) );

        const uint          dim    = ncols;
        const ulong         nnodes = nrows;
        vector< double * >  vertices( nnodes );

        if ( dim == 3 )
        {
            for ( ulong  i = 0; i < nnodes; ++i )
            {
                std::getline( in, line );
                
                double      x, y, z;
                const bool  r   = phrase_parse( line.begin(), line.end(),
                                                ( double_[ ref(x) = _1 ] >>
                                                  double_[ ref(y) = _1 ] >>
                                                  double_[ ref(z) = _1 ] ),
                                                space );
                
                if ( ! r )
                    HERROR( ERR_FMT_MTX, "(TMMVectorIO) read",
                            "missing data, expected <float> <float> <float>" );

                double * v  = new double[ 3 ];

                v[0] = x;
                v[1] = y;
                v[2] = z;
                
                vertices[i] = v;
            }// for
        }// if
        else
            HERROR( ERR_NOT_IMPL, "", "" );

        return make_unique< TCoordinate >( vertices, dim );
    }// if
    else
        HERROR( ERR_NOT_IMPL, "", "" );

    return nullptr;
}

///////////////////////////////////////////////////
// 
// TPLTMGCoordIO
//

void
TPLTMGCoordIO::write ( const TCoordinate *, const string & ) const
{
    HERROR( ERR_NOT_IMPL, "(TPLTMGCoordIO) write", "" );
}

unique_ptr< TCoordinate >
TPLTMGCoordIO::read ( const string & filename ) const
{
    //
    // read coordinate file
    //

    if ( ! boost::filesystem::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TPLTMGCoordIO) read", filename );
    
    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();
    
    uint                dim    = 0;
    ulong               nnodes = 0;
    vector< double * >  vertices;
    string              line;
    vector< string >    parts;
    
    std::getline( in, line );
    split( line, " \t\n", parts );

    if ( parts.size() < 2 )
        HERROR( ERR_FMT_SAMG, "(TPLTMGCoordIO) read", "expected <nvertices> <dim>" );
    
    nnodes = str_to_int( parts[0] );
    dim    = str_to_int( parts[1] );

    vertices.resize( nnodes );

    for ( uint i = 0; i < nnodes; i++ )
    {
        double * v  = new double[ dim ];
        uint     id = 0;
        
        std::getline( in, line );
        split( line, " \t\n", parts );

        if ( parts.size() < dim+1 )
            HERROR( ERR_FMT_SAMG, "(TPLTMGCoordIO) read",
                    to_string( "expected <id> and %d floats per coordinate", dim ) );

        id = str_to_int( parts[0] );
        
        for ( uint  j = 0; j < dim; ++j )
        {
            v[j] = str_to_dbl( parts[j+1] );
        }// for

        id--;
        
        vertices[id] = v;
    }// for

    return make_unique< TCoordinate >( vertices, dim );
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// TGMSHCoordIO
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//
// read in a grid from a file
//

unique_ptr< TCoordinate >
TGMSHCoordIO::read ( const string & filename ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TGMSHCoordIO) read", filename );
    
    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();

    /////////////////////////////////////////////////////////////////
    //
    // file format
    //

    string            line( size_t(1024), '\0' );
    vector< string >  parts;

    std::getline( in, line );
    boost::trim_right( line );

    if ( line != "$MeshFormat" )
        HERROR( ERR_GRID_FORMAT, "(TGMSHCoordIO) read", "file not in GMSH format" );
        
    std::getline( in, line );
    split( line, " \t", parts );

    if ( parts.size() != 3 )
        HERROR( ERR_GRID_FORMAT, "(TGMSHCoordIO) read", "file not in GMSH format" );
    
    if ( int( str_to_dbl( parts[0] ) ) != 2 )
        HERROR( ERR_GRID_FORMAT, "(TGMSHCoordIO) read", "unsupported format version" );

    std::getline( in, line );
    boost::trim_right( line );

    if ( line != "$EndMeshFormat" )
        HERROR( ERR_GRID_FORMAT, "(TGMSHCoordIO) read", "file not in GMSH format" );

    /////////////////////////////////////////////////////////////////
    //
    // read vertices
    //

    size_t             nvertices = 0;
    vector< T3Point >  vertices;
    
    std::getline( in, line );
    boost::trim_right( line );

    if ( line != "$Nodes" )
        HERROR( ERR_GRID_FORMAT, "(TGMSHCoordIO) read", "file not in GMSH format" );

    std::getline( in, line );

    nvertices = str_to_int( line );
    vertices.resize( nvertices );

    for ( size_t  i = 0; i < nvertices; ++i )
    {
        std::getline( in, line );
        split( line, " \t", parts );
        
        if ( parts.size() < 4 )
            HERROR( ERR_GRID_FORMAT, "(TGMSHCoordIO) read", "wrong vertex format" );

        //
        // read vertex coordinates
        //

        int     id;
        double  x, y, z;

        sscanf( line.c_str(), "%d %lf %lf %lf", & id, & x, & y, & z );

        id--;
        if ( id >= int(nvertices) )
            HERROR( ERR_GRID_FORMAT, "(TGMSHCoordIO) read", "vertex ID > no. of vertices" );
            
        //
        // append vertex
        //
        
        vertices[ id ] = T3Point( x, y, z );
    }// for

    std::getline( in, line );
    boost::trim_right( line );

    if ( line != "$EndNodes" )
        HERROR( ERR_GRID_FORMAT, "(TGMSHCoordIO) read", "file not in GMSH format" );

    return make_unique< TCoordinate >( vertices );
}

void
TGMSHCoordIO::write  ( const TCoordinate *,
                       const string & ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

}// namespace Hpro
