//
// Project     : HLIBpro
// File        : baseio.cc
// Description : basic IO related functions and classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/config.h"

#include <fstream>

#include <boost/filesystem.hpp>

#if HPRO_HAS_BOOST_IOSTREAMS == 1
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#endif

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

#include "hpro/base/error.hh"
#include "hpro/base/System.hh"

#include "baseio.hh"

namespace Hpro
{

namespace fs = boost::filesystem;

#if HPRO_HAS_BOOST_IOSTREAMS == 1
namespace io = boost::iostreams;
#endif

using std::string;
using std::vector;
using std::unique_ptr;

using boost::spirit::qi::long_;
using boost::spirit::qi::double_;
using boost::spirit::qi::_1;
using boost::spirit::qi::phrase_parse;
using boost::spirit::ascii::space;
using boost::phoenix::ref;
    
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// 
// byte swapping (little <-> big endian)
//
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

template <> void
swap_bytes<bool> ( bool  & )
{}

template <> void
swap_bytes<int8_t> ( int8_t & )
{}

template <> void
swap_bytes<uint8_t> ( uint8_t & )
{}

template <> void
swap_bytes<int16_t> ( int16_t  & d )
{
    int8_t * s = reinterpret_cast< int8_t * >( & d );
    std::swap( s[0], s[1] );
}

template <> void
swap_bytes<unsigned short> ( uint16_t  & d )
{
    uint8_t * s = reinterpret_cast< uint8_t * >( & d );
    std::swap( s[0], s[1] );
}

template <> void
swap_bytes<int32_t> ( int32_t  & d )
{
    int8_t * s = reinterpret_cast< int8_t * >( & d );
    std::swap( s[0], s[3] );
    std::swap( s[1], s[2] );
}

template <> void
swap_bytes<uint32_t> ( uint32_t  & d )
{
    uint8_t * s = reinterpret_cast< uint8_t * >( & d );
    std::swap( s[0], s[3] );
    std::swap( s[1], s[2] );
}

template <> void
swap_bytes<int64_t> ( int64_t  & d )
{
    int8_t * s = reinterpret_cast< int8_t * >( & d );

    std::swap( s[0], s[7] );
    std::swap( s[1], s[6] );
    std::swap( s[2], s[5] );
    std::swap( s[3], s[4] );
}

template <> void
swap_bytes<uint64_t> ( uint64_t & d )
{
    uint8_t * s = reinterpret_cast< uint8_t * >( & d );

    std::swap( s[0], s[7] );
    std::swap( s[1], s[6] );
    std::swap( s[2], s[5] );
    std::swap( s[3], s[4] );
}

template <> void
swap_bytes<float> ( float & d )
{
    char * s = reinterpret_cast< char * >( & d );
    std::swap( s[0], s[3] );
    std::swap( s[1], s[2] );
}

template <> void
swap_bytes<double> ( double & d )
{
    char * s = reinterpret_cast< char * >( & d );
    std::swap( s[0], s[7] );
    std::swap( s[1], s[6] );
    std::swap( s[2], s[5] );
    std::swap( s[3], s[4] );
}

template <> void
swap_bytes< std::complex< float > > ( std::complex< float > & d )
{
    float  r = std::real( d );
    float  i = std::imag( d );

    swap_bytes< float >( r );
    swap_bytes< float >( i );

    d = std::complex< float >( r, i );
}

template <> void
swap_bytes< std::complex< double > > ( std::complex< double > & d )
{
    double  r = std::real( d );
    double  i = std::imag( d );

    swap_bytes< double >( r );
    swap_bytes< double >( i );

    d = std::complex< double >( r, i );
}

///////////////////////////////////////////////////
// 
// open file with suitable stream, e.g. for
// compression
//

std::unique_ptr< std::istream >
open_read ( const std::string &  filename )
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "open_read", filename );
    
#if HPRO_HAS_BOOST_IOSTREAMS == 1
    
    auto      in  = std::make_unique< io::filtering_istream >();
    fs::path  filepath( filename );
    string    ext = filepath.extension().string();

    boost::to_lower( ext );
    
    if      ( ext == ".gz"  ) in->push( io::gzip_decompressor()  );
    else if ( ext == ".bz2" ) in->push( io::bzip2_decompressor() );
    else if ( ext == ".zip" ) in->push( io::zlib_decompressor()  );
    
    in->push( io::file_descriptor_source( filepath, std::ios::in | std::ios::binary ) );
    
    return in;

#else

    return std::make_unique< std::ifstream >( filename.c_str(), std::ios::in | std::ios::binary );
    
#endif
}

std::unique_ptr< std::ostream >
open_write ( const std::string &  filename )
{
#if HPRO_HAS_BOOST_IOSTREAMS == 1
    
    auto      out = std::make_unique< io::filtering_ostream >();
    fs::path  filepath( filename );
    string    ext = filepath.extension().string();
    
    boost::to_lower( ext );

    if      ( ext == ".gz"  ) out->push( io::gzip_compressor()  );
    else if ( ext == ".bz2" ) out->push( io::bzip2_compressor() );
    else if ( ext == ".zip" ) out->push( io::zlib_compressor()  );
    
    out->push( io::file_descriptor_sink( filepath, std::ios::out | std::ios::binary ) );
    
    return out;

#else

    return std::make_unique< std::ofstream >( filename.c_str(), std::ios::out | std::ios::binary );
    
#endif
}

///////////////////////////////////////////////////
// 
// return extension without .gz etc.
//

string
extension ( const string &  filename )
{
    const fs::path  filepath( filename );
    fs::path        ext( filepath.extension() );

    if ( ext == ".gz" )
    {
        // look further
        const fs::path  gzpath = filepath.stem();
        const fs::path  gzext  = gzpath.extension();

        if ( ! gzext.empty() )
            ext = gzext;
    }// if

    string  sext( ext.string() );

    boost::trim_left_if( sext, boost::is_any_of( "." ) );
    boost::to_lower( sext );

    return sext;
}

///////////////////////////////////////////////////
// 
// file format detection
//

fmt_t
guess_format ( const string &  filename )
{
    /////////////////////////////////////////////////////////
    //
    // read content of file and try to guess format
    //

    if ( boost::filesystem::exists( filename ) )
    {
        try
        {
            //
            // HLIBpro format
            //

            {
                auto            in_ptr = open_read( filename );
                std::istream &  in     = * in_ptr.get();
        
                uint16_t  endian;
                char      text[126];
        
                in.read( text, 126 );
                in.read( reinterpret_cast< char * >( & endian ), sizeof(endian) );

                if (( endian == 0x484d ) || ( endian == 0x4d48 ))
                {
                    bool  byteswap = false;
            
                    if      ( endian == 0x484d ) byteswap = false;
                    else if ( endian == 0x4d48 ) byteswap = true;

                    uint16_t type;

                    in.read( reinterpret_cast< char * >( & type ), sizeof(type) );

                    if ( byteswap ) swap_bytes( type );

                    if ( type <= 8 )
                    {
                        // up to here, HLIBpro-type ("HM") and matrix-type is correct
                        // so we assume, it is a HLIBpro file
                        return FMT_HPRO;
                    }// if
                }// if
            }

            //
            // Matlab format
            //

            {
                auto            in_ptr = open_read( filename );
                std::istream &  in     = * in_ptr.get();
        
                int16_t  version, endian;
                char     text[124];
        
                in.read( text, 124 );
                in.read( reinterpret_cast< char * >( & version ), sizeof(version)  );
                in.read( reinterpret_cast< char * >( & endian ),  sizeof(endian)  );

                if (( endian == 0x4d49 ) || ( endian == 0x494d ))
                {
                    bool  byteswap = false;
            
                    if      ( endian == 0x4d49 ) byteswap = false;
                    else if ( endian == 0x494d ) byteswap = true;

                    int32_t ltype, type, size;

                    in.read( reinterpret_cast< char * >( & ltype ), sizeof(ltype) );

                    if ( byteswap ) swap_bytes( type );

                    if ( (( ltype >> 16 ) & 0xff) != 0 )
                    {
                        // small data element
                        type = ltype & 0xff;
                        size = ( ltype >> 16 ) & 0xff;
                    }// if
                    else
                    {
                        // normal data element
                        type = ltype;
                
                        in.read( reinterpret_cast< char * >( & size ), sizeof(size) );
                        if ( byteswap ) swap_bytes( size );
                    }// else
            
                    if (( type == 14 ) || ( type == 15 ))
                    {
                        // up to here, Matlab-type ("MI") and matrix-/compression-type is correct
                        // so we assume, it is a Matlab file
                        return FMT_MATLAB;
                    }// if
                }// if
            }

            //
            // HDF5 format
            //
            
            {
                const uchar     format_signature[] = { 0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a };
                const uint      sign_size = sizeof( format_signature ) / sizeof( uchar );
                auto            in_ptr    = open_read( filename );
                std::istream &  in        = * in_ptr.get();
                uchar           buf[ sign_size ];
                bool            is_hdf5 = true;

                in.read( reinterpret_cast< char * >( buf ), sign_size );

                for ( uint  i = 0; i < sign_size; ++i )
                    if ( buf[i] != format_signature[i] )
                        is_hdf5 = false;

                if ( is_hdf5 )
                    return FMT_HDF5;
            }
            
            //
            // Harwell/Boeing format
            //

            {
                auto            in_ptr = open_read( filename );
                std::istream &  in     = * in_ptr.get();
        
                string  line, mxtype;
            
                std::getline( in, line );
                std::getline( in, line );
                std::getline( in, line );
        
                mxtype = line.substr( 0, 3 );
        
                if (( mxtype[0] == 'R' ) || ( mxtype[0] == 'C' ) ||
                    ( mxtype[0] == 'P' ) || ( mxtype[0] == 'I' ))
                {
                    if (( mxtype[1] == 'U' ) || ( mxtype[1] == 'S' ) || ( mxtype[1] == 'H' ))
                    {
                        if (( mxtype[2] == 'A' ) || ( mxtype[2] == 'E' ))
                        {
                            // matrix type is valid, so assume Harwell/Boeing format
                            return FMT_HB;
                        }// if
                    }// if
                }// if
            }

            //
            // Matrix-Market format
            //

            {
                auto            in_ptr = open_read( filename );
                std::istream &  in     = * in_ptr.get();
        
                string            line;
                vector< string >  parts;
                bool              found_header = false;
                bool              is_mtx       = true;
            
                for ( uint l = 0; is_mtx && ( l < 100 ); l++ )
                {
                    std::getline( in, line );
                    split( line, " \t\r\n", parts );

                    if ( parts[0] == "%%MatrixMarket" )
                    {
                        found_header = true;
                    
                        if ( parts.size() < 3 )
                        {
                            is_mtx = false;
                            continue;
                        }// if
            
                        if ( parts[1] != "matrix" )
                        {
                            is_mtx = false;
                            continue;
                        }// if
            
                        if (( parts[2] != "coordinate" ) && ( parts[2] != "array" ))
                        {
                            is_mtx = false;
                            continue;
                        }// if

                        // look for qualifiers
                        for ( uint i = 3; i < parts.size(); i++ )
                        {
                            if ( ! (( parts[i] == "real" ) ||
                                    ( parts[i] == "integer" ) ||
                                    ( parts[i] == "complex" ) ||
                                    ( parts[i] == "pattern" ) ||
                                    ( parts[i] == "general" ) ||
                                    ( parts[i] == "symmetric" ) ||
                                    ( parts[i] == "skew-symmetric" ) ||
                                    ( parts[i] == "hermitian" )) )
                            {
                                is_mtx = false;
                                break;
                            }// if
                        }// for
            
                        break;
                    }// if
                }// for

                if ( found_header && is_mtx )
                    return FMT_MTX;
            }

            //
            // grid formats
            //

            {
                auto            in_ptr = open_read( filename );
                std::istream &  in     = * in_ptr.get();
        
                string  line;

                for ( uint i = 0; i < 10; i++ )
                {
                    std::getline( in, line );
                    boost::trim_right( line );
            
                    if ( line == "surfacemesh" )
                    {
                        return FMT_SURFMESH;
                    }// if
                    else if ( line == "ply" )
                    {
                        return FMT_PLY;
                    }// if
                    else if ( line == "$MeshFormat" )
                    {
                        return FMT_GMSH;
                    }// if
                    else if (( line.find( "nv " ) != string::npos ) ||
                             ( line.find( "ne " ) != string::npos ) ||
                             ( line.find( "nt " ) != string::npos ))
                    {
                        return FMT_HPRO_GRID;
                    }// if
                }// for
            }
        }// try
        catch ( Error & )
        {
            // some kind of error occured, fall back to filename based detection
        }// catch
    }// if

    //
    // try to guess format by looking at extension
    //

    const string  ext = extension( filename );

    if (( ext == "amg" ) || ( ext == "rhs" ) || ( ext == "sol" ) || (  ext == "coo" ))
        return FMT_SAMG;
    else if (( ext == "hb"  ) || ( ext == "rb" ) ||
             ( ext == "rua" ) || ( ext == "rsa" ) || ( ext == "psa" ))
        return FMT_HB;
    else if (( ext == "mat" ) || ( ext == "m" ))
        return FMT_MATLAB;
    else if ( ext == "hm" )
        return FMT_HPRO;
    else if (( ext == "hdf" ) || ( ext == "h5" ) || ( ext == "hdf5" ))
        return FMT_HDF5;
    else if (( ext == "mtx" ) || ( ext == "mm" ))
        return FMT_MTX;
    else if ( ext == "tri"  )
        return FMT_HPRO_GRID;
    else if ( ext == "ply"  )
        return FMT_PLY;
    else if ( ext == "sur"  )
        return FMT_SURFMESH;
    else if (( ext == "msh" ) || ( ext == "gmsh" ) || ( ext == "gmsh2" ))
        return FMT_GMSH;
    else
        return FMT_UNKNOWN;
}

// external functions for corresponding formats
variant_id_t hb_guess_value_type     ( const string &  filename );
variant_id_t matlab_guess_value_type ( const string &  filename, const string &  name );
variant_id_t hpro_guess_value_type   ( const string &  filename );
variant_id_t mtx_guess_value_type    ( const string &  filename );

variant_id_t
guess_value_type ( const string &  filename,
                   const string &  name )
{
    auto  fmt = guess_format( filename );

    switch ( fmt )
    {
        case FMT_SAMG   : return REAL_FP64;
        case FMT_HB     : return hb_guess_value_type( filename );
        case FMT_MATLAB : return matlab_guess_value_type( filename, name );
        case FMT_HPRO   : return hpro_guess_value_type( filename );
        case FMT_MTX    : return mtx_guess_value_type( filename );
        default         :
            HERROR( ERR_FMT_UNKNOWN, "guess_value_type", filename );
    }// switch
}

///////////////////////////////////////////////////
// 
// lexical cast with default value
//

long
str_to_int ( const string &  str,
             const long      def_val )
{
    long  i = 0;
    bool  r = phrase_parse( str.begin(), str.end(), long_[ ref(i) = _1 ], space );

    if ( ! r )
        return def_val;

    return i;
}

long
str_to_int ( const char *  begin,
             const char *  end,
             const long    def_val )
{
    long  i = 0;
    bool  r = phrase_parse( begin, end, long_[ ref(i) = _1 ], space );

    if ( ! r )
        return def_val;

    return i;
}

double
str_to_dbl ( const string &  str,
             const double    def_val )
{
    double  d = 0;
    bool    r = phrase_parse( str.begin(), str.end(), double_[ ref(d) = _1 ], space );

    if ( ! r )
        return def_val;

    return d;
}

double
str_to_dbl ( const char *  begin,
             const char *  end,
             const double  def_val )
{
    double  d = 0;
    bool    r = phrase_parse( begin, end, double_[ ref(d) = _1 ], space );

    if ( ! r )
        return def_val;

    return d;
}

///////////////////////////////////////////////////
//
// append file suffix if not already present
//
std::string
add_extension ( const std::string &  filename,
                const std::string &  suffix )
{
    //
    // check, if a suffix is present, e.g. ".???"
    //

    const fs::path  filepath( filename );
    
    if ( filepath.has_extension() )
    {
        //
        // return original filename
        //
        
        return filename;
    }// if
    else
    {
        //
        // append suffix
        //

        std::string  new_name;

        new_name = filename + "." + suffix;

        return new_name;
    }// else
}

}// namespace Hpro
