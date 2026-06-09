//
// Project     : HLIBpro
// File        : matlab_io.cc
// Description : Matlab related IO functions and classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <cstring>
#include <vector>
#include <sstream>
#include <memory>

#include "hpro/config.h"

#if HPRO_HAS_BOOST_IOSTREAMS == 1
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#endif

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include "hpro/matrix/TSparseMatrix.hh"
#include "hpro/matrix/THMatrix.hh"
#include "hpro/matrix/structure.hh"

#include "baseio.hh"

#include "hpro/io/TMatrixIO.hh"
#include "hpro/io/TVectorIO.hh"
#include "hpro/io/TCoordIO.hh"

namespace Hpro
{

namespace fs = boost::filesystem;

#if HPRO_HAS_BOOST_IOSTREAMS == 1
namespace io = boost::iostreams;
#endif

using std::unique_ptr;
using std::make_unique;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// local defines
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// padding for file acces in Matlab
#define PAD( s )             ((s) <= 4 ? (4 - (s)) : (((s)%8 == 0 ? 0 : 8 - ((s)%8))))
#define PAD_SEEKG( s, file ) { if (( (s) != 0 ) && ( PAD(s) != 0 )) { (file).seekg( PAD(s), std::ios_base::cur ); } }
#define PAD_SEEKP( s, file ) { if (( (s) != 0 ) && ( PAD(s) != 0 )) { (file).seekp( PAD(s), std::ios_base::cur ); } }

// for automatic swapping of bytes
#define BYTE_SWAP( t )  { if ( byteswap ) swap_bytes( t ); }

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// local types
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

namespace
{

//
// Matlab element and matrix types with their names
//
enum matmi_t
{
    MI_INT8 = 1,
    MI_UINT8,
    MI_INT16,
    MI_UINT16,
    MI_INT32,
    MI_UINT32,
    MI_SINGLE, 
    MI_RESERVED8,
    MI_DOUBLE,
    MI_RESERVED10,
    MI_RESERVED11,
    MI_INT64,
    MI_UINT64,
    MI_MATRIX,
    MI_COMPRESSED,
    MI_UTF8,
    MI_UTF16,
    MI_UTF32
};

enum matmx_t
{
    MX_CELL = 1,
    MX_STRUCTURE,
    MX_OBJECT,
    MX_CHAR_ARRAY,
    MX_SPARSE_ARRAY,
    MX_DBL_ARRAY,
    MX_SGL_ARRAY,
    MX_INT8_ARRAY,
    MX_UINT8_ARRAY,
    MX_INT16_ARRAY,
    MX_UINT16_ARRAY,
    MX_INT32_ARRAY,
    MX_UINT32_ARRAY
};

}// namespace anonymous

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// local functions
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

namespace
{

//
// template functions for reading and writing data
//
template < typename value_t >
value_t
mat5_read ( std::istream &  in,
            const bool      byteswap )
{
    value_t  val;
    
    in.read( reinterpret_cast< char * >( & val ), sizeof(value_t) );
    BYTE_SWAP( val );

    return val;
}

template < typename value_t >
void
mat5_read ( std::istream &  in,
            value_t *       data,
            const uint      nelem,
            const bool      byteswap )
{
    in.read( reinterpret_cast< char * >(  data ), sizeof(value_t)*nelem );
    
    if ( byteswap && ( sizeof(value_t) > 1 ))
    {
        for ( uint i = 0; i < nelem; i++ )
            BYTE_SWAP( data[i] );
    }// if
}

template < typename value_t >
void
mat5_write ( std::ostream &  out,
             value_t         data )
{
    out.write( reinterpret_cast< const char * >( & data ), sizeof(value_t) );
}

template < typename value_t >
void
mat5_write ( std::ostream &   out,
             const value_t *  data,
             const size_t     nelem )
{
    out.write( reinterpret_cast< const char * >(  data ), sizeof(value_t)*nelem );
}

void
mat5_write ( std::ostream &       out,
             const std::string &  str )
{
    BOOST_FOREACH( int8_t  ch, str )
    {
        mat5_write< int8_t >( out, ch );
    }// for
}

//
// read element tag with type and size
//
void
mat5_readtag ( std::istream & in,
               matmi_t &      type,
               int32_t &      size,
               const bool     byteswap )
{
    int32_t  ltype = 0;

    //
    // read first part of tag with type information
    //
    
    in.read( reinterpret_cast< char * >( & ltype ), sizeof(ltype) );
    BYTE_SWAP( ltype );

    //
    // type might also contain size of "small elements"
    //
    
    if ( (( ltype >> 16 ) & 0xff) != 0 )
    {
        // small data element
        type = matmi_t( ltype & 0xff );
        size = ( ltype >> 16 ) & 0xff;
    }// if
    else
    {
        // normal data element
        type = matmi_t( ltype );
        
        in.read( reinterpret_cast< char * >( & size ), sizeof(size) );
        BYTE_SWAP( size );
    }// else
}

//
// write element tag with type and size
//
void
mat5_writetag ( std::ostream & out,
                const int32_t  atype,
                const int32_t  size )
{
    int32_t  type = atype;
    
    if ( size <= 4 )
    {
        // write small element tag
        type |= (size & 0xff) << 16;
        out.write( reinterpret_cast< const char * >( & type ), sizeof(type) );
    }// if
    else
    {
        // write standard element tag
        out.write( reinterpret_cast< const char * >( & type ), sizeof(type) );
        out.write( reinterpret_cast< const char * >( & size ), sizeof(size) );
    }// else
}

//
// read coefficient of array and return as double
//
template < typename value_t >
value_t
mat5_read_coeff ( std::istream &  in,
                  const bool      byteswap,
                  const int32_t   type )
{
    value_t  coeff;

    switch ( type )
    {
        case MI_INT8   : coeff = value_t(mat5_read<int8_t>(   in, byteswap )); break;
        case MI_UINT8  : coeff = value_t(mat5_read<uint8_t>(  in, byteswap )); break;
        case MI_INT16  : coeff = value_t(mat5_read<int16_t>(  in, byteswap )); break;
        case MI_UINT16 : coeff = value_t(mat5_read<uint16_t>( in, byteswap )); break;
        case MI_INT32  : coeff = value_t(mat5_read<int32_t>(  in, byteswap )); break;
        case MI_UINT32 : coeff = value_t(mat5_read<uint32_t>( in, byteswap )); break;
        case MI_SINGLE : coeff = value_t(mat5_read<float>(    in, byteswap )); break;
        case MI_DOUBLE : coeff = value_t(mat5_read<double>(   in, byteswap )); break;
                        
    default:
        HERROR( ERR_NOT_IMPL, "mat5_read_coeff",
                to_string( "unsupported Matlab type %d", type ) );
    }// switch

    return coeff;
}

//
// read sparse double array
//
template < typename value_t >
void
mat5_read_sparse ( std::istream &               in,
                   const bool                   byteswap,
                   const std::vector< uint > &  dims,
                   const uint                   nonzero,
                   TMatrix< value_t > **        M,
                   TVector< value_t > **        v )
{
    using  real_t = real_type_t< value_t >;
    
    std::unique_ptr< TSparseMatrix< value_t > >  S;
    std::unique_ptr< TScalarVector< value_t > >  x;

    if ( v != nullptr )
    {
        int  len = (dims[0] != 1 ? dims[0] : dims[1]);

        x.reset( new TScalarVector< value_t >( len, 0 ) );
    }// if
    else
    {
        S.reset( new TSparseMatrix< value_t >(  dims[0], dims[1] ) );
    }// else
    
    //
    // read row indices
    //

    matmi_t                 type = matmi_t( 0 );
    int32_t                 size = 0;
    size_t                  nnz  = nonzero;
    std::vector< int32_t >  rowptr( nnz );

    mat5_readtag( in, type, size, byteswap );

    if ( type != MI_INT32 )
        HERROR( ERR_FMT_MATLAB, "mat5_read_sparse", "expected INT32" );

    // adjust number of non-zero
    if ( size != int32_t( sizeof(int32_t) * nnz ) )
        nnz = size / sizeof(int32_t);
    
    for ( size_t  i = 0; i < nnz; i++ )
    {
        const int32_t row = mat5_read<int32_t>( in, byteswap );

        if ( row >= int32_t(dims[0]) )
            HERROR( ERR_NOT_IMPL, "mat5_read_sparse", "rectangular sparse matrices" );

        rowptr[i] = row;
    }// for
    
    PAD_SEEKG( size, in );
    
    //
    // read column indices
    //

    std::vector< int32_t >  colind( dims[1]+1 );

    mat5_readtag( in, type, size, byteswap );

    if ( type != MI_INT32 )
        HERROR( ERR_FMT_MATLAB, "mat5_read_sparse", "expected INT32" );

    if ( size != int32_t( sizeof(int32_t) * (dims[1]+1) ) )
        HERROR( ERR_SIZE, "mat5_read_sparse", "while reading column entries" );
    
    for ( uint i = 0; i <= dims[1]; i++ )
    {
        const int32_t  col = mat5_read<int32_t>( in, byteswap );

        if ( col > int32_t( nnz ) )
            HERROR( ERR_CONSISTENCY, "mat5_read_sparse", "column index too large" );

        colind[i] = col;
    }// for
    
    PAD_SEEKG( size, in );

    // again, adjust number of non-zero with last number in colind
    if ( int32_t( nnz ) != colind[dims[1]] )
        nnz = colind[dims[1]];
    
    // build correct row indices for vectors in case of
    // column-vectors
    if (( x.get() != nullptr ) && ( dims[0] == 1 ))
    {
        uint  pos = 0;
        
        for ( uint i = 0; i < dims[1]; i++ )
            if ( colind[i] != colind[i+1] )
                rowptr[pos++] = i;
    }// if

    // build row-indices in CRS data of sparse matrices
    if ( S.get() != nullptr )
    {
        S->init( nnz );

        //
        // build row-wise numbering for sparse matrix
        //
        
        std::vector< int32_t >  tmp( dims[0], 0 );
        
        for ( uint j = 0; j < dims[1]; j++ )
        {
            const uint lb = colind[j];
            const uint ub = colind[j+1];

            for ( uint i = lb; i < ub; i++ )
                tmp[rowptr[i]]++;
        }// for

        //
        // build row indices in sparse matrix
        //

        uint  pos = 0;
        
        for ( uint i = 0; i < dims[0]; i++ )
        {
            S->rowptr(i) = pos;

            pos += tmp[i];
        }// for
        
        S->rowptr(dims[0]) = idx_t(nnz);
    }// if        
    
    //
    // read real parts
    //

    mat5_readtag( in, type, size, byteswap );

    if ( x.get() != nullptr )
    {
        for ( size_t  i = 0; i < nnz; i++ )
        {
            const auto  coeff = mat5_read_coeff< real_t >( in, byteswap, type );

            x->set_entry( rowptr[i], coeff );
        }// for
    }// if
    else
    {
        std::vector< int32_t >  tmp( dims[0], 0 );
        
        for ( uint col = 0; col < dims[1]; col++ )
        {
            const uint lb = colind[col];
            const uint ub = colind[col+1];

            for ( uint i = lb; i < ub; i++ )
            {
                const idx_t  row   = rowptr[i];
                const idx_t  idx   = S->rowptr(row) + tmp[row];
                const auto   coeff = mat5_read_coeff< real_t >( in, byteswap, type );
                
                S->colind(idx) = col;
                S->coeff(idx) = coeff;

                tmp[row]++;
            }// for
        }// for
    }// else
    
    PAD_SEEKG( size, in );

    //
    // read imaginary parts
    //

    if ( is_complex_type< value_t >::value )
    {
        mat5_readtag( in, type, size, byteswap );
        
        if ( x.get() != nullptr )
        {
            for ( size_t  i = 0; i < nnz; i++ )
            {
                const auto  coeff = mat5_read_coeff< real_t >( in, byteswap, type );

                x->set_entry( rowptr[i], get_value< value_t >::compose( std::real( x->entry( rowptr[i] ) ), coeff ) );
            }// for
        }// if
        else
        {
            std::vector< int32_t >  tmp( dims[0], 0 );
        
            for ( uint col = 0; col < dims[1]; col++ )
            {
                const uint lb = colind[col];
                const uint ub = colind[col+1];

                for ( uint i = lb; i < ub; i++ )
                {
                    const idx_t  row   = rowptr[i];
                    const idx_t  idx   = S->rowptr(row) + tmp[row];
                    const auto   coeff = mat5_read_coeff< real_t >( in, byteswap, type );
            
                    S->coeff( idx ) = get_value< value_t >::compose( std::real( S->coeff(idx) ), coeff );

                    tmp[row]++;
                }// for
            }// for
        }// else

        PAD_SEEKG( size, in );
    }// if

    //
    // finally sort entries in sparse matrices
    //

    if ( S.get() != nullptr )
    {
        S->sort_entries();
        S->test_symmetry();
    }// if

    //
    // return local objects or delete them
    //
    
    if ( x.get() != nullptr )
    {
        if ( *v == nullptr )
            *v = x.release();
    }// if
    else if ( S.get() != nullptr )
    {
        if ( *M == nullptr )
            *M = S.release();
    }// if
}

//
// read full double matrix
//
template < typename value_t >
void
mat5_read_dblarray ( std::istream &               in,
                     const bool                   byteswap,
                     const std::vector< uint > &  dims,
                     TMatrix< value_t > **        M,
                     TVector< value_t > **        v )
{
    using  real_t = real_type_t< value_t >;
    
    std::unique_ptr< TDenseMatrix< value_t > >   D;
    std::unique_ptr< TScalarVector< value_t > >  x;

    if ( v != nullptr )
    {
        const uint  len = (dims[0] != 1 ? dims[0] : dims[1]);

        x = std::make_unique< TScalarVector< value_t > >( len, 0 );
    }// if
    else
    {
        D = std::make_unique< TDenseMatrix< value_t > >();

        if ( D == nullptr )
            HERROR( ERR_MEM, "mat5_read_dblarray", "" );
        
        D->set_size( dims[0], dims[1] );
    }// else
    
    //
    // read real parts
    //

    matmi_t  type = matmi_t( 0 );
    int32_t  size = 0;

    mat5_readtag( in, type, size, byteswap );

    if ( x.get() != nullptr )
    {
        const uint  len = (dims[0] != 1 ? dims[0] : dims[1]);

        for ( uint i = 0; i < len; i++ )
        {
            const auto  coeff = mat5_read_coeff< real_t >( in, byteswap, type );
            
            x->set_entry( i, coeff );
        }// for
    }// if
    else
    {
        for ( uint j = 0; j < dims[1]; j++ )
            for ( uint i = 0; i < dims[0]; i++ )
            {
                const auto  coeff = mat5_read_coeff< real_t >( in, byteswap, type );
            
                D->set_entry( i, j, coeff );
            }// for
    }// else
    
    PAD_SEEKG( size, in );

    //
    // read imaginary parts
    //

    if ( is_complex_type< value_t >::value )
    {
        mat5_readtag( in, type, size, byteswap );
        
        if ( x.get() != nullptr )
        {
            const uint  len = (dims[0] != 1 ? dims[0] : dims[1]);
            
            for ( uint i = 0; i < len; i++ )
            {
                const auto  coeff = mat5_read_coeff< real_t >( in, byteswap, type );
            
                x->set_entry( i, get_value< value_t >::compose( std::real( x->entry( i ) ), coeff ) );
            }// for
        }// if
        else
        {
            for ( uint j = 0; j < dims[1]; j++ )
                for ( uint i = 0; i < dims[0]; i++ )
                {
                    const auto  coeff = mat5_read_coeff< real_t >( in, byteswap, type );
                    
                    D->set_entry( i, j, get_value< value_t >::compose( std::real( D->entry( i, j ) ), coeff ) );
                }// for
        }// else
        
        PAD_SEEKG( size, in );
    }// if

    //
    // test for symmetry or hermitian
    //

    if ( D.get() != nullptr )
    {
        if ( is_complex_type< value_t >::value )
        {
            bool  is_sym  = true;
            bool  is_herm = true;
        
            for ( uint i = 0; i < dims[0]; i++ )
                for ( uint j = i+1; j < dims[1]; j++ )
                {
                    if ( D->entry( i, j ) != D->entry( j, i ) )
                        is_sym = false;

                    if ( D->entry( i, j ) != Math::conj( D->entry( j, i ) ) )
                        is_herm = false;

                    if ( ! is_sym && ! is_herm )
                        break;
                }// for

            if ( is_sym )
                D->set_symmetric();
            else if ( is_herm )
                D->set_hermitian();
        }// if
        else
        {
            bool  is_sym = true;
        
            for ( uint i = 0; i < dims[0]; i++ )
                for ( uint j = i+1; j < dims[1]; j++ )
                    if ( D->entry( i, j ) != D->entry( j, i ) )
                    {
                        is_sym = false;
                        break;
                    }// if

            if ( is_sym )
                D->set_symmetric();
        }// else
    }// if

    //
    // return local objects or delete them
    //
    
    if ( x.get() != nullptr )
    {
        if ( *v == nullptr )
            *v = x.release();
    }// if
    else if ( D.get() != nullptr )
    {
        if ( *M == nullptr )
            *M = D.release();
    }// if
}

//
// read a single element from file
//
template < typename value_t >
void
mat5_read_element ( std::istream &         in,
                    const bool             byteswap,
                    const std::string &    name,
                    const bool             is_subfield,
                    const std::string &    fieldname,
                    TMatrix< value_t > **  M,
                    TVector< value_t > **  v )
{
    int32_t  element_size = 0;
    matmi_t  type = matmi_t( 0 );
        
    mat5_readtag( in, type, element_size, byteswap );

    int32_t  size = element_size;
    auto     pos  = in.tellg();
        
    if ( type == MI_COMPRESSED )
    {
        #if HPRO_HAS_BOOST_IOSTREAMS == 1
        
        std::string  inbuf( size, '\0' );

        in.read( const_cast< char * >( inbuf.data() ), size );
        
        //
        // read first 8 bytes to obtain size
        //
        
        uint32_t  tmp[2] = { 0, 0 };

        {
            io::filtering_istream  zin;
            std::istringstream     strin( inbuf );

            zin.push( io::zlib_decompressor()  );
            zin.push( strin );

            zin.read( reinterpret_cast< char * >( tmp ), sizeof(tmp) );
        }

        //
        // uncompress all and continue reading matlab format
        // - can not uncompress on the fly using zlib_decompressor as it
        //   does not support tell/seek
        // - "+8" gives correct size
        //
        
        std::string  outbuf( tmp[1] + 8, '\0' );

        {
            io::filtering_istream  zin;
            std::istringstream     strin( inbuf );

            zin.push( io::zlib_decompressor()  );
            zin.push( strin );

            zin.read( const_cast< char * >( outbuf.data() ), tmp[1] + 8 );
        }
        
        std::istringstream  strout( outbuf );
        
        mat5_read_element( strout, byteswap, name, is_subfield, fieldname, M, v );

        // go to point after datastructure in case we did not
        // read everything
        in.seekg( pos );
        in.seekg( element_size, std::ios_base::cur );
        
        PAD_SEEKG( element_size, in );
                
        return;

        #else
        
        HERROR( ERR_NOT_IMPL, "mat5_read_element", "compression not supported" );
        
        #endif
    }// if
    else if ( type != MI_MATRIX )
        HERROR( ERR_NOT_IMPL, "mat5_read_element", "only miMatrix is supported" );

    //
    // read array flags
    //

    mat5_readtag( in, type, size, byteswap );

    if ( type != MI_UINT32 )
        HERROR( ERR_FMT_MATLAB, "mat5_read_element", "expected UINT32" );

    if ( size != 8 )
        HERROR( ERR_SIZE, "mat5_read_element", "expected 2 * UINT32" );

    const uint32_t flags  = mat5_read<uint32_t>( in, byteswap );
    const bool     imag   = ((flags & 0x0800) != 0);
    const matmx_t  aclass = matmx_t( flags & 0xff );

    if ( imag && ! is_complex_type< value_t >::value )
        HERROR( ERR_REAL_CMPLX, "mat5_read_element", "found complex valued data but real valued requested" );
    
    PAD_SEEKG( size, in ); 
    
    //
    // read number of non-zeroes in case of sparse matrix
    // or junk in case of a dense matrix
    //
    
    const uint32_t  nnz = mat5_read<uint32_t>( in, byteswap );
    PAD_SEEKG( size, in );
    
    //
    // read dimensions 
    //
    
    mat5_readtag( in, type, size, byteswap );

    if ( type != MI_INT32 )
        HERROR( ERR_FMT_MATLAB, "mat5_read_dblarray", "expected INT32" );

    const uint           ndims    = size / 4;
    std::vector< uint >  dims( ndims );
    uint                 nentries = 1;
    bool                 skip     = false;  // indicates skipping of data element

    for ( uint i = 0; i < ndims; i++ )
    {
        dims[i]   = mat5_read<int32_t>( in, byteswap );
        nentries *= dims[i];
    }// for
    
    PAD_SEEKG( size, in ); 

    if ( ndims != 2 )
        HERROR( ERR_NOT_IMPL, "mat5_read_dblarray", "only 2-dim. matrices supported" );

    if ( aclass != MX_STRUCTURE )
    {
        // check if vector is sought but matrix found
        if (( M == nullptr ) && ( v != nullptr ) && ( dims[0] > 1 ) && ( dims[1] > 1 ))
            skip = true;
        
        // check if matrix is sought but vector found
        if (( M != nullptr ) && ( v == nullptr ) && (( dims[0] == 1 ) || ( dims[1] == 1 )))
            skip = true;
    }// if
    
    //
    // read name of array
    //
    
    mat5_readtag( in, type, size, byteswap );

    if ( type != MI_INT8 )
        HERROR( ERR_FMT_MATLAB, "mat5_read_dblarray", "expected INT8" );

    std::string  element_name( size_t(size), '\0' );

    // element_name[size] = '\0';
    
    in.read( const_cast< char * >( element_name.data() ), size );
    PAD_SEEKG( size, in );

    // in case of a subfield: adjust element name
    if (( element_name == "" ) && is_subfield )
        element_name = fieldname;

    // adjust "skipping" if requested name equals element_name
    if (( name != "" ) && ( name != element_name ))
        skip = true;
        
    //
    // read fields based on class value
    //

    if ( ! skip )
    {
        if ( aclass == MX_STRUCTURE )
        {
            mat5_readtag( in, type, size, byteswap );

            if ( type != MI_INT32 )
                HERROR( ERR_FMT_MATLAB, "mat5_read_element", "expected INT32" );

            const int32_t  fnlen = mat5_read<int32_t>( in, byteswap );
        
            mat5_readtag( in, type, size, byteswap );

            if ( type != MI_INT8 )
                HERROR( ERR_FMT_MATLAB, "mat5_read_element", "expected INT8" );

            const uint                  fields = size / fnlen;
            std::vector< std::string >  fieldnames( fields );

            for ( uint i = 0; i < fields; i++ )
            {
                fieldnames[i].resize( fnlen+1 );
                fieldnames[i][fnlen] = '\0';
                in.read( const_cast< char * >( fieldnames[i].data() ), fnlen );
            }// for

            PAD_SEEKG( size, in );

            for ( uint i = 0; i < fields; i++ )
            {
                mat5_read_element( in, byteswap, name, true, fieldnames[i], M, v );

                if (( M != nullptr ) && ( *M != nullptr ))
                    break;
                else if (( v != nullptr ) && ( *v != nullptr ))
                    break;
            }// for
        }// if
        else if ( aclass == MX_SPARSE_ARRAY )
        {
            if (( name == "" ) || ( name == element_name ))
                mat5_read_sparse< value_t >( in, byteswap, dims, nnz, M, v );
        }// if
        else if ( aclass == MX_DBL_ARRAY )
        {
            if (( name == "" ) || ( name == element_name ))
                mat5_read_dblarray< value_t >( in, byteswap, dims, M, v );
        }// if
        else if ( aclass == MX_SGL_ARRAY )
        {
            if (( name == "" ) || ( name == element_name ))
                mat5_read_dblarray< value_t >( in, byteswap, dims, M, v );
        }// if
        // else
        //     HERROR( ERR_TYPE, "mat5_read_element", "only dense and sparse matrices supported" );
    }// if

    // go to point after datastructure in case we did not
    // read everything
    in.seekg( pos );
    in.seekg( element_size, std::ios_base::cur );
        
    PAD_SEEKG( element_size, in ); 
}

bool
mat5_guess_type ( std::istream &         in,
                  const bool             byteswap,
                  const std::string &    name,
                  const bool             is_subfield,
                  const std::string &    fieldname,
                  variant_id_t &         value_type )
{
    bool     found        = false;
    int32_t  element_size = 0;
    matmi_t  type         = matmi_t( 0 );
        
    mat5_readtag( in, type, element_size, byteswap );

    int32_t  size = element_size;
    auto     pos  = in.tellg();
        
    if ( type == MI_COMPRESSED )
    {
        #if HPRO_HAS_BOOST_IOSTREAMS == 1
        
        std::string  inbuf( size, '\0' );

        in.read( const_cast< char * >( inbuf.data() ), size );
        
        //
        // read first 8 bytes to obtain size
        //
        
        uint32_t  tmp[2] = { 0, 0 };

        {
            io::filtering_istream  zin;
            std::istringstream     strin( inbuf );

            zin.push( io::zlib_decompressor()  );
            zin.push( strin );

            zin.read( reinterpret_cast< char * >( tmp ), sizeof(tmp) );
        }

        //
        // uncompress all and continue reading matlab format
        // - can not uncompress on the fly using zlib_decompressor as it
        //   does not support tell/seek
        // - "+8" gives correct size
        //
        
        std::string  outbuf( tmp[1] + 8, '\0' );

        {
            io::filtering_istream  zin;
            std::istringstream     strin( inbuf );

            zin.push( io::zlib_decompressor()  );
            zin.push( strin );

            zin.read( const_cast< char * >( outbuf.data() ), tmp[1] + 8 );
        }
        
        std::istringstream  strout( outbuf );

        while ( strout.good() )
        {
            if ( mat5_guess_type( strout, byteswap, name, is_subfield, fieldname, value_type ) )
            {
                found = true;
                break;
            }// if
        }// while

        // go to point after datastructure in case we did not
        // read everything
        in.seekg( pos );
        in.seekg( element_size, std::ios_base::cur );
        
        PAD_SEEKG( element_size, in );
                
        return found;
#else
        
        HERROR( ERR_NOT_IMPL, "mat5_guess_type", "compression not supported" );
        
#endif
    }// if
    else if ( type != MI_MATRIX )
        HERROR( ERR_NOT_IMPL, "mat5_guess_type", "only miMatrix is supported" );

    //
    // read array flags
    //

    mat5_readtag( in, type, size, byteswap );

    if ( type != MI_UINT32 )
        HERROR( ERR_FMT_MATLAB, "mat5_guess_type", "expected UINT32" );

    if ( size != 8 )
        HERROR( ERR_SIZE, "mat5_guess_type", "expected 2 * UINT32" );

    const uint32_t flags  = mat5_read<uint32_t>( in, byteswap );
    const bool     imag   = ((flags & 0x0800) != 0);
    const matmx_t  aclass = matmx_t( flags & 0xff );

    PAD_SEEKG( size, in ); 
    
    //
    // read number of non-zeroes in case of sparse matrix
    // or junk in case of a dense matrix
    //
    
    mat5_read<uint32_t>( in, byteswap );
    PAD_SEEKG( size, in );
    
    //
    // read dimensions 
    //
    
    mat5_readtag( in, type, size, byteswap );

    if ( type != MI_INT32 )
        HERROR( ERR_FMT_MATLAB, "mat5_guess_type", "expected INT32" );

    const uint           ndims    = size / 4;
    std::vector< uint >  dims( ndims );
    uint                 nentries = 1;
    bool                 skip     = false;  // indicates skipping of data element

    for ( uint i = 0; i < ndims; i++ )
    {
        dims[i]   = mat5_read<int32_t>( in, byteswap );
        nentries *= dims[i];
    }// for
    
    PAD_SEEKG( size, in ); 

    //
    // read name of array
    //
    
    mat5_readtag( in, type, size, byteswap );

    if ( type != MI_INT8 )
        HERROR( ERR_FMT_MATLAB, "mat5_guess_type", "expected INT8" );

    std::string  element_name( size_t(size), '\0' );

    // element_name[size] = '\0';
    
    in.read( const_cast< char * >( element_name.data() ), size );
    PAD_SEEKG( size, in );

    // in case of a subfield: adjust element name
    if (( element_name == "" ) && is_subfield )
        element_name = fieldname;

    // adjust "skipping" if requested name equals element_name
    if (( name != "" ) && ( name != element_name ))
        skip = true;
        
    //
    // read fields based on class value
    //

    if ( ! skip )
    {
        if ( aclass == MX_STRUCTURE )
        {
            mat5_readtag( in, type, size, byteswap );

            if ( type != MI_INT32 )
                HERROR( ERR_FMT_MATLAB, "mat5_guess_type", "expected INT32" );

            const int32_t  fnlen = mat5_read<int32_t>( in, byteswap );
        
            mat5_readtag( in, type, size, byteswap );

            if ( type != MI_INT8 )
                HERROR( ERR_FMT_MATLAB, "mat5_guess_type", "expected INT8" );

            const uint                  fields = size / fnlen;
            std::vector< std::string >  fieldnames( fields );

            for ( uint i = 0; i < fields; i++ )
            {
                fieldnames[i].resize( fnlen+1 );
                fieldnames[i][fnlen] = '\0';
                in.read( const_cast< char * >( fieldnames[i].data() ), fnlen );
            }// for

            PAD_SEEKG( size, in );

            for ( uint i = 0; i < fields; i++ )
            {
                if ( mat5_guess_type( in, byteswap, name, true, fieldnames[i], value_type ) )
                {
                    found = true;
                    break;
                }// if
            }// for
        }// if
        else
        {
            if ( imag )
                value_type = COMPLEX_FP64;
            else
                value_type = REAL_FP64;

            found = true;
        }// else
    }// if

    // go to point after datastructure in case we did not
    // read everything
    in.seekg( pos );
    in.seekg( element_size, std::ios_base::cur );
        
    PAD_SEEKG( element_size, in );

    return found;
}

//
// special write methods
//
template < typename value_t >
void
write_dense ( std::ostream &                   out,
              const BLAS::Matrix< value_t > &  D,
              const std::string &              matname )
{
    using  real_t = real_type_t< value_t >;
    
    // remember current position for size info
    const auto  pos = out.tellp();
        
    // size = 5 ensures standard element tag
    mat5_writetag( out, MI_MATRIX, 5 );

    //
    // write array flags
    // (if complex valued)
    //

    uint32_t  flags = ( is_complex_type< value_t >::value ? 0x0800 : 0 );
    uint32_t  tmp   = 0;

    flags |= is_single_prec< value_t >::value ? MX_SGL_ARRAY : MX_DBL_ARRAY;
    
    mat5_writetag( out, MI_UINT32, 8 );
    mat5_write<uint32_t>( out, flags );
    mat5_write<uint32_t>( out, tmp   );

    //
    // write dimensions
    //

    int32_t  rows = int32_t( D.nrows() );
    int32_t  cols = int32_t( D.ncols() );

    mat5_writetag( out, MI_INT32, 8 );
    mat5_write<int32_t>( out, rows );
    mat5_write<int32_t>( out, cols );

    //
    // write name
    //

    int32_t  size = int32_t( matname.length()*sizeof(int8_t) );
    
    mat5_writetag( out, MI_INT8, size );
    mat5_write( out, matname );

    PAD_SEEKP( size, out );
    
    //
    // write real part of array
    // - no padding neccessary since sizeof(double) == 8
    //

    matmi_t  mi_type = is_single_prec< value_t >::value ? MI_SINGLE : MI_DOUBLE;
    
    mat5_writetag( out, mi_type, int32_t( sizeof(real_t)*rows*cols ) );

    for ( int j = 0; j < cols; j++ )
        for ( int i = 0; i < rows; i++ )
        {
            mat5_write< real_t >( out, std::real( D(i,j) ) );
        }// for
    
    //
    // write imaginary part of array
    //

    if ( is_complex_type< value_t >::value )
    {
        mat5_writetag( out, mi_type, int32_t( sizeof(real_t)*rows*cols ) );

        for ( int j = 0; j < cols; j++ )
            for ( int i = 0; i < rows; i++ )
            {
                mat5_write< real_t >( out, std::imag( D(i,j) ) );
            }// for
    }// else

    //
    // compute size of object and rewrite array tag
    // (remove size of tag)
    //

    size = int32_t( out.tellp() - pos - 8 );

    out.seekp( pos );
        
    mat5_writetag( out, MI_MATRIX, size );
}

template < typename value_t >
void
write_rank ( std::ostream &                out,
             const TRkMatrix< value_t > *  R,
             const std::string &           matname )
{
    //
    // split low-rank matrix in A and B and write matname_A and matname_B
    // first : matrix A
    //
    
    // remember current position for size info
    const auto  pos = out.tellp();
        
    // size = 5 ensures standard element tag
    mat5_writetag( out, MI_MATRIX, 5 );

    //
    // write array flags
    // (if complex valued)
    //

    uint32_t  flags = ( R->is_complex() ? 0x0800 : 0 );
    uint32_t  tmp   = 0;

    flags |= MX_DBL_ARRAY;
    
    mat5_writetag( out, MI_UINT32, 8 );
    mat5_write<uint32_t>( out, flags );
    mat5_write<uint32_t>( out, tmp   );

    //
    // write dimensions
    //

    int32_t  rows = int32_t(R->rows());
    int32_t  rank = int32_t(R->rank());

    mat5_writetag( out, MI_INT32, 8 );
    mat5_write<int32_t>( out, rows );
    mat5_write<int32_t>( out, rank );

    //
    // write name
    //

    std::string   tname = matname + "_A";
    int32_t       size  = int(tname.length()*sizeof(int8_t));
    
    mat5_writetag( out, MI_INT8, size );
    mat5_write( out, tname );

    PAD_SEEKP( size, out );
    
    //
    // write real part of array
    // - no padding neccessary since sizeof(double) == 8
    //

    mat5_writetag( out, MI_DOUBLE, int32_t(sizeof(double)*rows*rank) );

    for ( int k = 0; k < rank; k++ )
        for ( int i = 0; i < rows; i++ )
        {
            mat5_write<double>( out, std::real( R->blas_mat_A()( i, k ) ) );
        }// for
    
    //
    // write imaginary part of array
    //

    if ( is_complex_type< value_t >::value )
    {
        mat5_writetag( out, MI_DOUBLE, int32_t(sizeof(double)*rows*rank) );

        for ( int k = 0; k < rank; k++ )
            for ( int i = 0; i < rows; i++ )
            {
                mat5_write<double>( out, std::imag( R->blas_mat_A()( i, k ) ) );
            }// for
    }// if

    // compute size of object and rewrite array tag
    // (remove size of tag)
    auto  fpos = out.tellp();

    size = int32_t( fpos - pos - 8 );

    out.seekp( pos );
    mat5_writetag( out, MI_MATRIX, size );
    out.seekp( fpos );
    
    //
    // now write matrix B
    //

    // remember current position for size info
    fpos = out.tellp();
        
    // size = 5 ensures standard element tag
    mat5_writetag( out, MI_MATRIX, 5 );

    //
    // write array flags
    // (if complex valued)
    //

    tmp = 0;

    mat5_writetag( out, MI_UINT32, 8 );
    mat5_write<uint32_t>( out, flags );
    mat5_write<uint32_t>( out, tmp   );

    //
    // write dimensions
    //

    int32_t  cols = int32_t(R->cols());

    mat5_writetag( out, MI_INT32, 8 );
    mat5_write<int32_t>( out, cols );
    mat5_write<int32_t>( out, rank );

    //
    // write name
    //

    tname = matname + "_B";
    size  = int32_t(tname.length()*sizeof(int8_t));
    
    mat5_writetag( out, MI_INT8, size );
    mat5_write( out, tname );

    PAD_SEEKP( size, out );
    
    //
    // write real part of array
    // - no padding neccessary since sizeof(double) == 8
    //

    mat5_writetag( out, MI_DOUBLE, int32_t(sizeof(double)*cols*rank) );

    for ( int k = 0; k < rank; k++ )
        for ( int i = 0; i < cols; i++ )
        {
            mat5_write<double>( out, std::real( R->blas_mat_B()( i, k ) ) );
        }// for
    
    //
    // write imaginary part of array
    //

    if ( is_complex_type< value_t >::value )
    {
        mat5_writetag( out, MI_DOUBLE, int32_t(sizeof(double)*cols*rank) );

        for ( int k = 0; k < rank; k++ )
            for ( int i = 0; i < cols; i++ )
            {
                mat5_write<double>( out, std::imag( R->blas_mat_B()( i, k ) ) );
            }// for
    }// if

    // compute size of object and rewrite array tag
    // (remove size of tag)
    size = int32_t( out.tellp() - fpos - 8);

    out.seekp( fpos );
    mat5_writetag( out, MI_MATRIX, size );
}

template < typename value_t >
void
write_sparse ( std::ostream &                    out,
               const TSparseMatrix< value_t > *  S,
               const std::string &               matname )
{
    // remember current position for size info
    auto  fpos = out.tellp();
        
    // size = 5 ensures standard element tag
    mat5_writetag( out, MI_MATRIX, 5 );
    
    //
    // write array flags
    // (if complex valued)
    //

    uint32_t  flags = ( S->is_complex() ? 0x0800 : 0 );
    uint32_t  nnz   = uint32_t(S->n_non_zero());

    flags |= MX_SPARSE_ARRAY;
    
    mat5_writetag( out, MI_UINT32, 8 );
    mat5_write<uint32_t>( out, flags );
    mat5_write<uint32_t>( out, nnz   );

    //
    // write dimensions
    //

    int32_t  rows = int32_t(S->rows());
    int32_t  cols = int32_t(S->cols());

    mat5_writetag( out, MI_INT32, 8 );
    mat5_write<int32_t>( out, rows );
    mat5_write<int32_t>( out, cols );

    //
    // write name
    //

    int32_t  size = int32_t(matname.length()*sizeof(int8_t));
    
    mat5_writetag( out, MI_INT8, size );
    mat5_write( out, matname );

    PAD_SEEKP( size, out );
    
    //
    // build column-wise numbering for sparse matrix
    //

    std::vector< int32_t >  tmp( S->cols(), 0 );
    
    for ( int i = 0; i < rows; i++ )
    {
        const idx_t  lb = S->rowptr(i);
        const idx_t  ub = S->rowptr(i+1);
        
        for ( idx_t  j = lb; j < ub; j++ )
            tmp[ S->colind(j) ]++;
    }// for

    //
    // build column indices in CCS format
    //

    std::vector< int32_t >  ccs_colind( cols+1 );
    int                     pos = 0;
    
    for ( int i = 0; i < cols; i++ )
    {
        ccs_colind[i] = pos;
        pos           += tmp[i];
    }// for
    
    ccs_colind[cols] = nnz;
        
    //
    // build row indices in CCS format
    //

    std::vector< int32_t >  ccs_rowptr( nnz );

    for ( int i = 0; i < cols; i++ )
        tmp[i] = 0;

    for ( int i = 0; i < rows; i++ )
    {
        const idx_t  lb = S->rowptr(i);
        const idx_t  ub = S->rowptr(i+1);

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col = S->colind(j);
            const idx_t  idx = ccs_colind[col] + tmp[col];
                
            ccs_rowptr[ idx ] = i;
            tmp[col]++;
        }// for
    }// for
    
    //
    // build new coeff array in CCS format
    //

    std::vector< double >  ccs_coeff( nnz );

    for ( int i = 0; i < cols; i++ )
        tmp[i] = 0;
    
    for ( int i = 0; i < rows; i++ )
    {
        const idx_t  lb = S->rowptr(i);
        const idx_t  ub = S->rowptr(i+1);

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col = S->colind(j);
            const idx_t  idx = ccs_colind[col] + tmp[col];

            ccs_coeff[ idx ] = std::real( S->coeff(j) );
            
            tmp[col]++;
        }// for
    }// for
    
    //
    // write row indices
    //

    mat5_writetag( out, MI_INT32, int32_t( sizeof(int32_t)*nnz ) );
    mat5_write<int32_t>( out, ccs_rowptr.data(), nnz );
    PAD_SEEKP( sizeof(int32_t)*nnz, out );
    
    //
    // write column indices
    //

    mat5_writetag( out, MI_INT32, int32_t( sizeof(int32_t)*(cols+1) ) );
    mat5_write<int32_t>( out, ccs_colind.data(), cols+1 );
    PAD_SEEKP( sizeof(int32_t)*(cols+1), out );

    //
    // write real parts
    //

    mat5_writetag( out, MI_DOUBLE, int32_t( sizeof(double)*nnz ) );

    for ( uint i = 0; i < nnz; i++ )
        mat5_write<double>( out, ccs_coeff[i] );

    //
    // read imaginary parts
    //

    if ( is_complex_type< value_t >::value )
    {
        for ( int i = 0; i < cols; i++ )
            tmp[i] = 0;
        
        for ( uint i = 0; i < nnz; i++ )
            ccs_coeff[i] = 0.0;
        
        for ( int i = 0; i < rows; i++ )
        {
            const idx_t  lb = S->rowptr(i);
            const idx_t  ub = S->rowptr(i+1);
            
            for ( idx_t  j = lb; j < ub; j++ )
            {
                const idx_t  col = S->colind(j);
                const idx_t  idx = ccs_colind[col] + tmp[col];
                
                ccs_coeff[ idx ] = std::imag( S->coeff(j) );
                
                tmp[col]++;
            }// for
        }// for
    
        mat5_writetag( out, MI_DOUBLE, int32_t( sizeof(double)*nnz ) );
        mat5_write<double>( out, ccs_coeff.data(), nnz );
    }// if
    
    //
    // compute size of object and rewrite array tag
    // (remove size of tag)
    //

    size = int32_t( long(out.tellp()) - fpos - 8 );

    out.seekp( fpos );
        
    mat5_writetag( out, MI_MATRIX, size );
}

template < typename value_t >
void
write_scalar ( std::ostream &                   out,
               const BLAS::Vector< value_t > &  v,
               const std::string &              vecname )
{
    // remember current position for size info
    auto  pos = out.tellp();
        
    // size = 5 ensures standard element tag
    mat5_writetag( out, MI_MATRIX, 5 );

    //
    // write array flags
    // (if complex valued)
    //

    uint32_t  flags = ( is_complex_type< value_t >::value ? 0x0800 : 0 );
    uint32_t  tmp   = 0;

    flags |= MX_DBL_ARRAY;
    
    mat5_writetag( out, MI_UINT32, 8 );
    mat5_write<uint32_t>( out, flags );
    mat5_write<uint32_t>( out, tmp   );

    //
    // write dimensions
    //

    int32_t  dim = int32_t( v.length() );

    mat5_writetag( out, MI_INT32, 8 );
    mat5_write<int32_t>( out, dim );
    mat5_write<int32_t>( out, 1 );

    //
    // write name
    //

    int32_t  size = int32_t( vecname.length()*sizeof(int8_t) );
    
    mat5_writetag( out, MI_INT8, size );
    mat5_write( out, vecname );

    PAD_SEEKP( size, out );
    
    //
    // write real part of array
    // - no padding neccessary since sizeof(double) == 8
    //

    mat5_writetag( out, MI_DOUBLE, int32_t( sizeof(double)*dim ) );

    for ( int i = 0; i < dim; i++ )
        mat5_write<double>( out, std::real( v(i) ) );
    
    //
    // write imaginary part of array
    //

    if ( is_complex_type< value_t >::value )
    {
        mat5_writetag( out, MI_DOUBLE, int32_t( sizeof(double)*dim ) );

        for ( int i = 0; i < dim; i++ )
            mat5_write<double>( out, std::imag( v(i) ) );
    }// else

    //
    // compute size of object and rewrite array tag
    // (remove size of tag)
    //

    size = int32_t( out.tellp() - pos - 8 );

    out.seekp( pos );
        
    mat5_writetag( out, MI_MATRIX, size );
}

}// namespace anonymous

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// class TMatlabMatrixIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// write matrix to file
//
template < typename value_t >
void
TMatlabMatrixIO::write ( const TMatrix< value_t > *  A,
                         const std::string &         filename ) const
{
    return write( A, filename, "M" );
}

//
// write matrix with given name (options are obsolete here)
//
template < typename value_t >
void
TMatlabMatrixIO::write ( const TMatrix< value_t > *  A,
                         const std::string &         filename,
                         const std::string &         matname ) const
{
    if ( A == nullptr )
        HERROR( ERR_ARG, "(TMatlabMatrixIO) write", "matrix is NULL" );

    std::unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &                   out = * out_ptr.get();

    int16_t  version   = 0x0100;
    int16_t  endian    = 0x4d49; // = MI
    char     buf[116]  = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char     subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    if ( is_dense( A ) )
    {
        write_dense( out, cptrcast( A, TDenseMatrix< value_t > )->blas_mat(), matname );
    }// if
    else if ( is_lowrank( A ) )
    {
        write_rank( out, cptrcast( A, TRkMatrix< value_t > ), matname );
    }// if
    else if ( is_sparse( A ) )
    {
        write_sparse( out, cptrcast( A, TSparseMatrix< value_t > ), matname );
    }// if
    else
        HERROR( ERR_MAT_TYPE, "(TMatlabMatrixIO) write", "unsupported matrix type " + A->typestr() );
}

//
// write linear operator \a A with name \a mname to file \a fname
//
template < typename value_t >
void
TMatlabMatrixIO::write ( const TLinearOperator< value_t > *  A,
                         const std::string &                 /* fname */,
                         const std::string &                 /* mname */ ) const
{
    HERROR( ERR_MAT_TYPE, "(TMatlabMatrixIO) write", A->typestr() );
}

//
// write BLAS matrix with given name
//
template < typename value_t >
void
TMatlabMatrixIO::write ( const BLAS::Matrix< value_t > &  A,
                         const std::string &              filename,
                         const std::string &              matname ) const
{
    std::unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &                   out = * out_ptr.get();

    int16_t  version   = 0x0100;
    int16_t  endian    = 0x4d49; // = MI
    char     buf[116]  = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char     subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    write_dense( out, A, matname );
}

//
// read matrix from file (assuming only one entry available)
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TMatlabMatrixIO::read ( const std::string & filename ) const
{
    return read< value_t >( filename, "" );
}

//
// read matrix with name <matname> from file
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TMatlabMatrixIO::read  ( const std::string & filename,
                         const std::string & matname ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TMatlabMatrixIO) read", filename );
    
    std::unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &                   in = * in_ptr.get();

    int16_t               version, endian;
    bool                  byteswap = false;
    TMatrix< value_t > *  M        = nullptr;

    in.seekg( 124 );
    in.read( reinterpret_cast< char * >( & version ), sizeof(version) );
    in.read( reinterpret_cast< char * >( & endian ),  sizeof(endian)  );
    
    if      ( endian == 0x4d49 )  byteswap = false; // MI
    else if ( endian == 0x494d )  byteswap = true;  // IM
    else
        HERROR( ERR_FMT_MATLAB, "(TMatlabMatrixIO) read",
                "endianess indicator wrong" );

    BYTE_SWAP( version );
        
    while (( M == nullptr ) && in.good() )
        mat5_read_element< value_t >( in, byteswap, matname, false, "", & M, nullptr );
    
    return std::unique_ptr< TMatrix< value_t > >( M );
}

#define INST_MATRIXIO( type )                               \
    template void TMatlabMatrixIO::write< type > ( const TMatrix< type > *, const std::string & ) const; \
    template void TMatlabMatrixIO::write< type > ( const TMatrix< type > *, const std::string &, const std::string & ) const; \
    template void TMatlabMatrixIO::write< type > ( const TLinearOperator< type > *, const std::string & , const std::string & ) const; \
    template void TMatlabMatrixIO::write< type > ( const BLAS::Matrix< type > &, const std::string &, const std::string & ) const; \
    template std::unique_ptr< TMatrix< type > > TMatlabMatrixIO::read< type > ( const std::string & ) const; \
    template std::unique_ptr< TMatrix< type > > TMatlabMatrixIO::read< type > ( const std::string &, const std::string & ) const;

INST_MATRIXIO( float )
INST_MATRIXIO( double )
INST_MATRIXIO( std::complex< float > )
INST_MATRIXIO( std::complex< double > )
    
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// class TMatlabVectorIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// write vector with name "v" to file <filename>
//
template < typename value_t >
void
TMatlabVectorIO::write ( const TVector< value_t > *  v,
                         const std::string &         filename ) const
{
    write( v, filename, "v" );
}

//
// write vector with name <vecname> to file <filename>
//
template < typename value_t >
void
TMatlabVectorIO::write ( const TVector< value_t > *  v,
                         const std::string &         filename,
                         const std::string &         vecname ) const
{
    if ( v == nullptr )
        HERROR( ERR_ARG, "(TMatlabVectorIO) write", "vector is NULL" );

    if ( ! IS_TYPE( v, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TMatlabVectorIO) write", v->typestr() );

    std::unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &                   out = * out_ptr.get();
    
    const auto  x        = cptrcast( v, TScalarVector< value_t > );
    int16_t     version  = 0x0100;
    int16_t     endian   = 0x4d49; // = MI
    char        buf[116] = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char        subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    write_scalar( out, x->blas_vec(), vecname );
}

template < typename value_t >
void
TMatlabVectorIO::write ( const BLAS::Vector< value_t > &  v,
                         const std::string &              filename,
                         const std::string &              vecname ) const
{
    std::unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &                   out = * out_ptr.get();
    
    int16_t  version  = 0x0100;
    int16_t  endian   = 0x4d49; // = MI
    char     buf[116] = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char     subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    write_scalar( out, v, vecname );
}

//
// read first vector from file <filename>
//
template < typename value_t >
std::unique_ptr< TVector< value_t > >
TMatlabVectorIO::read ( const std::string &  filename ) const
{
    return read< value_t >( filename, "" );
}

//
// read vector named <vecname> from file (vecname = "" means first vector avail.)
//
template < typename value_t >
std::unique_ptr< TVector< value_t > >
TMatlabVectorIO::read ( const std::string &  filename,
                        const std::string &  vecname ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TMatlabVectorIO) read", filename );
    
    std::unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();

    int16_t               version, endian;
    bool                  byteswap = false;
    TVector< value_t > *  v = nullptr;

    in.seekg( 124 );
    in.read( reinterpret_cast< char * >( & version ), sizeof(version) );
    in.read( reinterpret_cast< char * >( & endian ),  sizeof(endian)  );
    
    if      ( endian == 0x4d49 )  byteswap = false; // MI
    else if ( endian == 0x494d )  byteswap = true;  // IM
    else
        HERROR( ERR_FMT_MATLAB, "(TMatlabVectorIO) read", "endianess indicator wrong" );

    BYTE_SWAP( version );
        
    while (( v == nullptr ) && in.good() )
        mat5_read_element< value_t >( in, byteswap, vecname, false, "", nullptr, & v );
    
    return std::unique_ptr< TVector< value_t > >( v );
}

#define INST_VECTORIO( type )                                       \
    template void TMatlabVectorIO::write< type > ( const TVector< type > *, const std::string & ) const; \
    template void TMatlabVectorIO::write< type > ( const TVector< type > *, const std::string &, const std::string & ) const; \
    template void TMatlabVectorIO::write< type > ( const BLAS::Vector< type > &, const std::string &, const std::string & ) const; \
    template std::unique_ptr< TVector< type > > TMatlabVectorIO::read< type > ( const std::string &, const std::string & ) const; \
    template std::unique_ptr< TVector< type > > TMatlabVectorIO::read< type > ( const std::string & ) const;

INST_VECTORIO( float )
INST_VECTORIO( double )
INST_VECTORIO( std::complex< float > )
INST_VECTORIO( std::complex< double > )
    
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// class TMatlabCoordIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void
TMatlabCoordIO::write ( const TCoordinate *  coord,
                        const std::string &  filename ) const
{
    write( coord, filename, "coord" );
}

void
TMatlabCoordIO::write ( const TCoordinate *  coord,
                        const std::string &  filename,
                        const std::string &  coordname ) const
{
    if ( coord == nullptr )
        HERROR( ERR_ARG, "(TMatlabCoordIO) write", "NULL coordinates" );
    
    //
    // convert coordinates to dense matrix and write that
    //

    const auto  ncoord = idx_t(coord->ncoord());
    const auto  ndim   = idx_t(coord->dim());
    
    BLAS::Matrix< double >  D( ncoord, ndim );

    for ( idx_t  i = 0; i < ncoord; ++i )
    {
        auto  p = coord->coord( i );
        
        for ( idx_t  j = 0; j < ndim; ++j )
        {
            D( i, j ) = p[j];
        }// for
    }// for

    TMatlabMatrixIO  mio;

    mio.write( D, filename, coordname );
}

//
// read first vertices in file <filename>
//
std::unique_ptr< TCoordinate >
TMatlabCoordIO::read ( const std::string &  filename ) const
{
    return read( filename, "" );
}

//
// read vertices named <cooname> from file <filename>
//
std::unique_ptr< TCoordinate >
TMatlabCoordIO::read ( const std::string &  filename,
                       const std::string &  cooname ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TMatlabCoordIO) read", filename );
    
    std::unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();

    int16_t  version, endian;
    bool     byteswap = false;
    auto     M        = std::unique_ptr< TMatrix< double > >();

    in.seekg( 124 );
    in.read( reinterpret_cast< char * >( & version ), sizeof(version) );
    in.read( reinterpret_cast< char * >( & endian ),  sizeof(endian)  );
    
    if      ( endian == 0x4d49 )  byteswap = false; // MI
    else if ( endian == 0x494d )  byteswap = true;  // IM
    else
        HERROR( ERR_FMT_MATLAB, "(TMatlabCoordIO) read",
                "endianess indicator wrong" );

    BYTE_SWAP( version );
        
    while (( M.get() == nullptr ) && in.good() )
    {
        TMatrix< double > *  MT = nullptr;
        
        mat5_read_element< double >( in, byteswap, cooname, false, "", & MT, nullptr );

        M = std::unique_ptr< TMatrix< double > >( MT );
        
        if ( ! is_dense( *M ) )
        {
            M.reset( nullptr );
        }// if
    }// while

    if ( M.get() != nullptr )
    {
        //
        // convert dense matrix into coordinates
        //

        auto                     D      = ptrcast( M.get(), TDenseMatrix< double > );
        const idx_t              ncoord = idx_t( D->rows() );
        const idx_t              dim    = idx_t( D->cols() );
        std::vector< double * >  vertices( ncoord );

        for ( idx_t  i = 0; i < ncoord; i++ )
        {
            vertices[i] = new double[ dim ];

            for ( idx_t  j = 0; j < dim; j++ )
                vertices[i][j] = D->entry( i, j );
        }// for

        return std::make_unique< TCoordinate >( vertices, uint( dim ) );
    }// if

    return nullptr;
}

variant_id_t
matlab_guess_value_type ( const std::string &  filename,
                          const std::string &  name )
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TMatlabMatrixIO) read", filename );
    
    std::unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &                   in = * in_ptr.get();

    int16_t  version, endian;
    bool     byteswap = false;

    in.seekg( 124 );
    in.read( reinterpret_cast< char * >( & version ), sizeof(version) );
    in.read( reinterpret_cast< char * >( & endian ),  sizeof(endian)  );
    
    if      ( endian == 0x4d49 )  byteswap = false; // MI
    else if ( endian == 0x494d )  byteswap = true;  // IM
    else
        HERROR( ERR_FMT_MATLAB, "(TMatlabMatrixIO) read",
                "endianess indicator wrong" );

    BYTE_SWAP( version );

    variant_id_t  vtype = REAL_FP64;
    
    while ( in.good() )
    {
        if ( mat5_guess_type( in, byteswap, name, false, "", vtype ) )
            break;
    }// while
    
    return vtype;
}

}// namespace Hpro
