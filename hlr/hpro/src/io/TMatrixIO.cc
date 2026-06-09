//
// Project     : HLIBpro
// File        : TMatrixIO.cc
// Description : classes for matrix input/output
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <string>
#include <sstream>
#include <vector>
#include <memory>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

#include "hpro/config.h"

#if HPRO_USE_HDF5 == 1

#  if (defined(__GNUC__) || defined(__clang__)) && !defined(__ICC)
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wold-style-cast"
#  endif

#define H5_USE_110_API // define to use classical API

#include "H5Cpp.h"

#  if (defined(__GNUC__) || defined(__clang__)) && !defined(__ICC)
#    pragma GCC diagnostic pop
#  endif

#endif

#if HPRO_USE_NETCDF == 1
#  include <netcdf.h>
#endif

#include "hpro/matrix/TSparseMatrix.hh"
#include "hpro/matrix/THMatrix.hh"
#include "hpro/matrix/structure.hh"

#include "baseio.hh"

#include "hpro/io/TMatrixIO.hh"

namespace Hpro
{

namespace fs = boost::filesystem;

using boost::to_lower_copy;
using boost::spirit::qi::long_;
using boost::spirit::qi::double_;
using boost::spirit::qi::_1;
using boost::spirit::qi::phrase_parse;
using boost::spirit::ascii::space;
using boost::phoenix::ref;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// class TAutoMatrixIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

template < typename value_t >
void
TAutoMatrixIO::write ( const TMatrix< value_t > *  A,
                       const std::string &   filename ) const
{
    write( A, filename, "" );
}

template < typename value_t >
void
TAutoMatrixIO::write ( const TMatrix< value_t > *  A,
                       const std::string &              filename,
                       const std::string &              matname ) const
{
    switch ( guess_format( filename ) )
    {
        case FMT_HPRO    :
        {
            auto  mio = std::make_unique< THLibMatrixIO >();

            mio->write( A, filename );
        }
        break;
            
        case FMT_MATLAB  :
        {
            auto  mio = std::make_unique< TMatlabMatrixIO >();

            mio->write( A, filename, matname );
        }
        break;
            
        case FMT_SAMG    :
        {
            auto  mio = std::make_unique< TSAMGMatrixIO >();

            mio->write( A, filename );
        }
        break;

        case FMT_HB      :
        {
            auto  mio = std::make_unique< THBMatrixIO >();

            mio->write( A, filename );
        }
        break;

        case FMT_MTX     :
        {
            auto  mio = std::make_unique< TMMMatrixIO >();

            mio->write( A, filename );
        }
        break;

        case FMT_HDF5    :
        {
            auto  mio = std::make_unique< THDF5MatrixIO >();

            mio->write( A, filename, matname );
        }
        break;

        case FMT_UNKNOWN :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoMatrixIO) write", "" );
            return;
    }// switch
}

template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TAutoMatrixIO::read  ( const std::string &  filename ) const
{
    return read< value_t >( filename, "" );
}

template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TAutoMatrixIO::read  ( const std::string &  filename,
                       const std::string &  matname ) const
{
    if ( filename == "" )
        HERROR( ERR_ARG, "(TAutoMatrixIO) read", "empty filename" );

    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TAutoMatrixIO) read", filename );
    
    switch ( guess_format( filename ) )
    {
        case FMT_HPRO    :
        {
            auto  mio = std::make_unique< THLibMatrixIO >();

            return mio->read< value_t >( filename );
        }
        break;
            
        case FMT_MATLAB  :
        {
            auto  mio = std::make_unique< TMatlabMatrixIO >();

            return mio->read< value_t >( filename, matname );
        }
        break;
            
        case FMT_SAMG    :
        {
            auto  mio = std::make_unique< TSAMGMatrixIO >();

            return mio->read< value_t >( filename );
        }
        break;

        case FMT_HB      :
        {
            auto  mio = std::make_unique< THBMatrixIO >();

            return mio->read< value_t >( filename );
        }
        break;

        case FMT_MTX     :
        {
            auto  mio = std::make_unique< TMMMatrixIO >();

            return mio->read< value_t >( filename );
        }
        break;

        case FMT_UNKNOWN :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoMatrixIO) write", "" );
            return std::unique_ptr< TMatrix< value_t > >();
    }// switch

    return std::unique_ptr< TMatrix< value_t > >();
}

template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
read_matrix  ( const std::string &  filename )
{
    return read_matrix< value_t >( filename, "" );
}

template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
read_matrix  ( const std::string &  filename,
               const std::string &  matname )
{
    TAutoMatrixIO  mio;

    return mio.read< value_t >( filename, matname );
}

template < typename value_t >
void
write_matrix  ( const TMatrix< value_t > *  A,
                const std::string &         filename )
{
    write_matrix( A, filename, "" );
}

template < typename value_t >
void
write_matrix  ( const TMatrix< value_t > *  A,
                const std::string &         filename,
                const std::string &         matname )
{
    TAutoMatrixIO  mio;

    mio.write( A, filename, matname );
}

template < typename value_t >
std::unique_ptr< TLinearOperator< value_t > >
read_linop  ( const std::string &  filename )
{
    THLibMatrixIO  mio;

    return mio.read_linop< value_t >( filename );
}

template < typename value_t >
void
write_linop ( const TLinearOperator< value_t > *  A,
              const std::string &                 filename )
{
    THLibMatrixIO  mio;

    mio.write_linop< value_t >( A, filename );
}

#define INST_AUTO( type ) \
    template void TAutoMatrixIO::write< type >                              ( const TMatrix< type > *, const std::string & ) const; \
    template void TAutoMatrixIO::write< type >                              ( const TMatrix< type > *, const std::string &, const std::string & ) const; \
    template std::unique_ptr< TMatrix< type > > TAutoMatrixIO::read< type > ( const std::string & ) const; \
    template std::unique_ptr< TMatrix< type > > TAutoMatrixIO::read< type > ( const std::string &, const std::string & ) const; \
    template std::unique_ptr< TMatrix< type > > read_matrix< type >         ( const std::string & ); \
    template void write_matrix< type >                                      ( const TMatrix< type > *, const std::string & ); \
    template std::unique_ptr< TLinearOperator< type > > read_linop< type >  ( const std::string & ); \
    template void write_linop< type >                                       ( const TLinearOperator< type > *, const std::string & );

INST_AUTO( float )
INST_AUTO( double )
INST_AUTO( std::complex< float > )
INST_AUTO( std::complex< double > )
    
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// class TOctaveMatrixIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

template < typename value_t >
void
TOctaveMatrixIO::write ( const TMatrix< value_t > *  A,
                         const std::string &         name ) const
{
    if ( A == nullptr )
        HERROR( ERR_ARG, "(TOctaveMatrixIO) write", "argument is nullptr" );
    
    auto            out_ptr = open_write( name );
    std::ostream &  out     = * out_ptr.get();
    
    if ( ! is_dense( A ) )
        HERROR( ERR_MAT_TYPE, "(TOctaveMatrixIO) write", "only dense matrices supported" );

    //
    // convert matrix to dense and then write dense matrix
    //

    auto  D = cptrcast( A, TDenseMatrix< value_t > );

    //
    // write header
    //

    const size_t  n = D->rows();
    const size_t  m = D->cols();

    out << "# Created by HLIBpro" << std::endl
        << "# name: A" << std::endl;
    if ( D->is_complex() ) out << "# type: complex matrix" << std::endl;
    else                   out << "# type: matrix" << std::endl;
    out << "# rows: " << n << std::endl
        << "# columns: " << m << std::endl;

    //
    // write matrix
    //

    if ( is_complex_type< value_t >::value )
    {
        for ( idx_t  i = 0; i < idx_t(n); i++ )
        {
            for ( idx_t  j = 0; j < idx_t(m); j++ )
            {
                const auto  f = D->entry(i,j);
                
                out << boost::format( "(%.16e,%.16e) " ) % std::real(f) % std::imag(f);
            }// for
            
            out << std::endl;
        }// for
    }// if
    else
    {
        for ( idx_t  i = 0; i < idx_t(n); i++ )
        {
            for ( idx_t  j = 0; j < idx_t(m); j++ )
                out << boost::format( "%.16e " ) % std::real( D->entry(i,j) );
            
            out << std::endl;
        }// for
    }// else
}

#define INST_OCTAVE( type ) \
    template void TOctaveMatrixIO::write< type > ( const TMatrix< type > *, const std::string & ) const;

INST_OCTAVE( float )
INST_OCTAVE( double )
INST_OCTAVE( std::complex< float > )
INST_OCTAVE( std::complex< double > )
    

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// class TSAMGMatrixIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// write matrix to file
//
template < typename value_t >
void
TSAMGMatrixIO::write ( const TMatrix< value_t > *  A,
                       const std::string &         filename ) const
{
    using  real_t  = real_type_t< value_t >;
    
    if ( A == nullptr )
        HERROR( ERR_ARG, "(TSAMGMatrixIO) write", "matrix is nullptr" );

    // can not write complex matrices
    if ( is_complex_type< value_t >::value )
        HERROR( ERR_NREAL, "(TSAMGMatrixIO) write", "SAMG only supports real valued data" );

    //
    // set up name for format file
    //

    fs::path  filepath( filename );
    fs::path  ext = filepath.extension();
    std::string    compress_ext;
    std::string    frm_filename;

    while (( to_lower_copy( ext.string() ) == ".gz"    ) ||
           ( to_lower_copy( ext.string() ) == ".bzip2" ) ||
           ( to_lower_copy( ext.string() ) == ".zlib"  ))
    {
        compress_ext = compress_ext + ext.string();
        filepath.replace_extension( "" );
        ext = filepath.extension();
    }// while

    if ( to_lower_copy( ext.string() ) == ".amg" )
        filepath.replace_extension( "" );
    
    frm_filename = filepath.string() + ".frm" + compress_ext;
            
    //
    // write files
    //
    
    if ( IS_TYPE( A, TSparseMatrix ) )
    {
        auto  S = cptrcast( A, TSparseMatrix< value_t > );

        //
        // open format file to write basic information
        //
        
        const size_t           nnz = S->n_non_zero();
        const size_t           n   = S->rows();
        ulong                  type;
        bool                   row_sum_zero = true;
        std::vector< real_t >  diag_values( S->rows(), 0.0 );
        
        // determine, if row-sum is zero in matrix
        for ( idx_t  row = 0; row < idx_t(S->rows()); row++ )
        {
            const idx_t  lb = S->rowptr(row);
            const idx_t  ub = S->rowptr(row+1);
            real_t       f  = 0.0;
            bool         found_diag = false;
            
            for ( idx_t j = lb; j < ub; j++ )
            {
                if ( S->colind(j) == idx_t(row) )
                {
                    diag_values[row] = std::real( S->coeff(j) );
                    found_diag       = true;
                }// if
                
                f += std::real( S->coeff(j) );
            }// for

            if ( ! found_diag )
                HERROR( ERR_NOT_IMPL, "", "add zero diagonal elements" );
            
            if ( Math::abs( f ) > real_t(1e-32) )
                row_sum_zero = false;
        }// for
        
        if ( row_sum_zero ) type = (S->is_symmetric() ? 11 : 21);
        else                type = (S->is_symmetric() ? 12 : 22);

        
        //
        // write format to frm file
        //
        
        auto            out_ptr = open_write( frm_filename );
        std::ostream &  frm_out = * out_ptr.get();
        
        frm_out << "f  4" << std::endl
                << nnz << " " << n << " " << type << " 0 0" << std::endl;
    
        //
        // write matrix to amg file
        //
        
        out_ptr = open_write( filename );

        std::ostream &  amg_out = * out_ptr.get();
        
        // first the row-index
        for ( uint row = 0; row <= n; row++ )
            amg_out << S->rowptr(row)+1 << std::endl;
        
        // then the column-index, but make sure that the diagonal entry comes first
        for ( uint row = 0; row < n; row++ )
        {
            const idx_t  lb = S->rowptr(row);
            const idx_t  ub = S->rowptr(row+1);
            
            // first the diagonal entry
            amg_out << row+1 << std::endl;
            
            // then the rest
            for ( idx_t  j = lb; j < ub; j++ )
            {
                if ( S->colind(j) != idx_t(row) )
                    amg_out << S->colind(j)+1 << std::endl;
            }// for
        }// for
        
        // finally the coefficients, again with the diag-entry first
        for ( uint row = 0; row < n; row++ )
        {
            const idx_t  lb = S->rowptr(row);
            const idx_t  ub = S->rowptr(row+1);
            
            amg_out << boost::format( "%.16E" ) % double(diag_values[row]) << std::endl;
            
            for ( idx_t  j = lb; j < ub; j++ )
            {
                if ( S->colind(j) != idx_t(row) )
                    amg_out << boost::format( "%.16E" ) % double(std::real(S->coeff(j))) << std::endl;
            }// for
        }// for
    }// if
    else if ( is_dense( A ) )
    {
        //
        // the matrix is converted to a dense matrix and then written
        // in sparse format
        //

        auto                 D            = cptrcast( A, TDenseMatrix< value_t > );
        ulong                nnz;
        const size_t         n            = D->rows();
        const size_t         m            = D->cols();
        ulong                type;
        bool                 row_sum_zero = true;
        std::vector< uint >  row_ptr( n+1, 0 );
        uint                 pos          = 0;
        
        // determine, if row-sum is zero in matrix
        for ( idx_t  i = 0; i < idx_t(n); i++ )
        {
            real_t  f = 0.0;

            row_ptr[i] = pos;
            
            for ( idx_t  j = 0; j < idx_t(m); j++ )
            {
                const auto  val = std::real( D->entry( i, j ) );

                if (( i == j ) || ( val != 0.0 ))
                    pos++;
                
                f += val;
            }// for
            
            if ( Math::abs( f ) > real_t(1e-32) )
                row_sum_zero = false;
        }// for

        row_ptr[n] = pos;
        nnz        = pos;
        
        if ( row_sum_zero ) type = (D->is_symmetric() ? 11 : 21);
        else                type = (D->is_symmetric() ? 12 : 22);
        
        //
        // write format to frm file
        //
        
        auto            out_ptr = open_write( frm_filename );
        std::ostream &  frm_out = * out_ptr.get();
        
        frm_out << "f  4" << std::endl
                << nnz << " " << n << " " << type << " 0 0" << std::endl;
    
        //
        // write matrix to amg file
        //
        
        out_ptr = open_write( filename );

        std::ostream &  amg_out = * out_ptr.get();
        
        // first the row-index
        for ( uint i = 0; i <= n; i++ )
            amg_out << row_ptr[i] + 1 << std::endl;
        
        // then the column-index, but make sure that the diagonal entry comes first
        for ( uint i = 0; i < n; i++ )
        {
            // first the diagonal entry
            amg_out << i+1 << std::endl;
            
            // then the rest
            for ( uint j = 0; j < m; j++ )
            {
                if (( j != i ) && ( std::real(D->entry( i, j )) != real_t(0) ))
                    amg_out << j+1 << std::endl;
            }// for
        }// for
        
        // finally the coefficients, again with the diag-entry first
        for ( uint i = 0; i < n; i++ )
        {
            amg_out << boost::format( "%.16E" ) % double(std::real(D->entry( i, i ))) << std::endl;
            
            for ( uint j = 0; j < m; j++ )
            {
                const double  val = std::real( D->entry( i, j ) );
                
                if (( j != i ) && ( val != 0.0 ))
                    amg_out << boost::format( "%.16E" ) % val << std::endl;
            }// for
        }// for
    }// if
    else
        HERROR( ERR_MAT_TYPE, "(TSAMGMatrixIO) write", "unsupported matrix type " + A->typestr() );
}

//
// read matrix from file
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TSAMGMatrixIO::read ( const std::string &  filename ) const
{
    //
    // open format file to read further information
    //

    fs::path     amgpath( filename );
    fs::path     ext = amgpath.extension();
    std::string  compress_ext;
    fs::path     frmpath;

    while (( to_lower_copy( ext.string() ) == ".gz"    ) ||
           ( to_lower_copy( ext.string() ) == ".bzip2" ) ||
           ( to_lower_copy( ext.string() ) == ".zlib"  ))
    {
        compress_ext = compress_ext + ext.string();
        amgpath.replace_extension( "" );
        ext = amgpath.extension();
    }// while
    
    if ( ext == ".amg" )
        frmpath = amgpath.replace_extension( ".frm" );
    else
        frmpath = amgpath.string() + ".frm";

    std::unique_ptr< std::istream >  in_ptr;
    
    if ( fs::exists( frmpath ) )
        in_ptr = open_read( frmpath.string() );
    else
    {
        if ( fs::exists( frmpath.string() + compress_ext ) )
            in_ptr = open_read( frmpath.string() + compress_ext );
        else
            HERROR( ERR_FNEXISTS, "(TSAMGMatrixIO) read", frmpath.string() + "(" + compress_ext + ")" );
    }// else

    std::istream &  frm_in = * in_ptr.get();
    
    ulong                       nnz, n, type;
    std::string                 line;
    std::vector< std::string >  parts;

    std::getline( frm_in, line );
    std::getline( frm_in, line );
    split( line, " \t\r\n", parts );

    if ( parts.size() < 3 )
        HERROR( ERR_FMT_SAMG, "(TSAMGMatrixIO) read", "expected <nnz> <nrows> <type> in format file" );
    
    nnz  = str_to_int( parts[0] );
    n    = str_to_int( parts[1] );
    type = str_to_int( parts[2] );

    //
    // read matrix
    //

    in_ptr = open_read( filename );

    std::istream &  amg_in = * in_ptr.get();

    auto  S = std::make_unique< TSparseMatrix< value_t > >();

    S->set_size( n, n );
    S->init( nnz );

    // first fill row-indices
    for ( uint i = 0; i <= n; i++ )
    {
        ulong  idx;
            
        std::getline( amg_in, line );
        idx = str_to_int( line )-1;

        if ( idx > nnz )
            HERROR( ERR_FMT_SAMG, "(TSAMGMatrixIO) read", "row pointer index > number of non-zeroes" );
        
        S->rowptr(i) = idx;
    }// for

    // now fill column indices
    for ( uint i = 0; i < nnz; i++ )
    {
        ulong  idx;
            
        std::getline( amg_in, line );
        idx = str_to_int( line )-1;

        if ( idx >= n )
            HERROR( ERR_FMT_SAMG, "(TSAMGMatrixIO) read", "column index > number of columns" );
        
        S->colind(i) = idx;
    }// for

    // and finally the coefficients
    for ( uint i = 0; i < nnz; i++ )
    {
        double  val;
            
        std::getline( amg_in, line );
        val = str_to_dbl( line );

        S->coeff(i) = value_t(val);
    }// for

    S->sort_entries();

    if (( type == 11 ) || ( type == 12 )) S->set_symmetric();
    else                                  S->set_nonsym();

    return std::unique_ptr< TMatrix< value_t > >( S.release() );
}

#define INST_SAMG( type ) \
    template void TSAMGMatrixIO::write< type > ( const TMatrix< type > *, const std::string & ) const; \
    template std::unique_ptr< TMatrix< type > > TSAMGMatrixIO::read< type > ( const std::string & ) const;

INST_SAMG( float )
INST_SAMG( double )
INST_SAMG( std::complex< float > )
INST_SAMG( std::complex< double > )

///////////////////////////////////////////////////
// 
// input and output in PLTMG format
//

//
// write matrix with name
//
template < typename value_t >
void
TPLTMGMatrixIO::write ( const TMatrix< value_t > *,
                        const std::string & ) const
{
    HERROR( ERR_NOT_IMPL, "(TPLTMGMatrixIO) write", "" );
}

//
// read matrix from file
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TPLTMGMatrixIO::read  ( const std::string &  filename ) const
{
    //
    // read matrix file
    //

    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TPLTMGMatrixIO) read", filename );
    
    auto            in_ptr = open_read( filename );
    std::istream &  in     = * in_ptr.get();
    std::string     line;
    uint            n;

    std::getline( in, line );
    n = str_to_int( line );
    std::getline( in, line );
    std::getline( in, line );

    //
    // read rhs
    //

    for ( uint i = 0; i < n; i++ )
        std::getline( in, line );

    //
    // read sparse matrix
    //

    auto  pos = in.tellg();
    uint  nnz = 0;
    
    std::vector< uint >        rows( n, 0 );
    std::vector< std::string >  parts;

    // first count number of entries per row to set up CRS structure
    while ( in.good() )
    {
        uint  i, j;
        
        std::getline( in, line );
        split( line, " \t\n", parts );
        if ( parts.size() < 2 )
            HERROR( ERR_FMT_PLTMG, "(TPLTMGMatrixIO) read", "expecting <row> <col>" );

        i = str_to_int( parts[0] );
        j = str_to_int( parts[1] );

        --i; --j;

        rows[i]++;
        nnz++;
    }// while

    // now reread file and fill CRS data
    auto   S   = std::make_unique< TSparseMatrix< value_t > >( n, n );
    idx_t  ofs = 0;

    S->init( nnz );

    in.seekg( pos );

    for ( uint i = 0; i < n; i++ )
    {
        S->rowptr(i) = ofs;
        ofs += rows[i];
        rows[i] = 0; // reset for later usage
    }// for

    S->rowptr(n) = nnz;

    while ( in.good() )
    {
        uint    i, j;
        double  entry;
        
        std::getline( in, line );
        split( line, " \t\n", parts );
        if ( parts.size() < 3 )
            HERROR( ERR_FMT_PLTMG, "(TPLTMGMatrixIO) read", "expecting <row> <col> <coeff>" );

        i     = str_to_int( parts[0] );
        j     = str_to_int( parts[1] );
        entry = str_to_dbl( parts[2] );

        --i; --j;

        ofs = S->rowptr(i) + rows[i];

        S->colind(ofs) = j;
        S->coeff(ofs)  = value_t(entry);
        
        rows[i]++;
    }// while
    
    return std::unique_ptr< TMatrix< value_t > >( S.release() );
}

#define INST_PLTMG( type ) \
    template void TPLTMGMatrixIO::write< type > ( const TMatrix< type > *, const std::string & ) const; \
    template std::unique_ptr< TMatrix< type > > TPLTMGMatrixIO::read< type > ( const std::string & ) const;

INST_PLTMG( float )
INST_PLTMG( double )
INST_PLTMG( std::complex< float > )
INST_PLTMG( std::complex< double > )

///////////////////////////////////////////////////
// 
// input and output in MatrixMarket format
//

namespace
{

enum mtxformat_t { MTX_ARRAY, MTX_COORD };
enum mtxfield_t  { MTX_REAL, MTX_INTEGER, MTX_COMPLEX, MTX_PATTERN };
enum mtxsym_t    { MTX_GENERAL, MTX_SYM, MTX_SKEWSYM, MTX_HERM };

}// namespace anonymous

//
// write matrix to <fname>
//
template < typename value_t >
void
TMMMatrixIO::write ( const TMatrix< value_t > *, const std::string & ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// read matrix from <fname>
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TMMMatrixIO::read  ( const std::string &  filename ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TMMMatrixIO) read", filename );
    
    auto            in_ptr = open_read( filename );
    std::istream &  in     = * in_ptr.get();

    ///////////////////////////////////////////////////
    //
    // read header and determine format of file
    //

    std::string   line;
    auto          parts = std::vector< std::string >();
    mtxformat_t   mat_format = MTX_ARRAY;
    mtxfield_t    mat_field  = MTX_REAL;
    mtxsym_t      mat_sym    = MTX_GENERAL;

    while ( in.good() )
    {
        std::getline( in, line );
        split( line, " \t\r\n", parts );

        if (( parts[0] == "%%MatrixMarket" ) || ( parts[0] == "%MatrixMarket" ))
        {
            if ( parts.size() < 3 )
                HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read",
                        "missing header data, expected <object> <format>" );
            
            if ( parts[1] != "matrix" )
                HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read", "file does not contain matrix" );
            
            if      ( parts[2] == "coordinate" ) mat_format = MTX_COORD;
            else if ( parts[2] == "array" )      mat_format = MTX_ARRAY;
            else
                HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read", "unknown matrix format \"" + parts[2] + "\"" );

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
                else if ( parts[i] != ""               )
                    HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read", "unknown qualifier \"" + parts[i] + "\"" );
            }// for
            
            break;
        }// if
    }// while

    if ( mat_sym == MTX_SKEWSYM )
        HERROR( ERR_NOT_IMPL, "(TMMMatrixIO) read", "skewsym matrices" );

    if (( mat_field == MTX_COMPLEX ) && ! is_complex_type< value_t >::value )
        HERROR( ERR_REAL_CMPLX, "(TMMMatrixIO) read", "complex data found but real data requested" );
    
    ///////////////////////////////////////////////////
    //
    // read content depending on format
    //

    if ( mat_format == MTX_COORD )
    {
        //
        // read number of rows, columns and non-zero entries
        //

        long  nrows = 0, ncols = 0, nnz = 0;

        while ( in.good() )
        {
            std::getline( in, line );

            if ( line[0] == '%' )
                continue;

            const bool  r = phrase_parse( line.begin(), line.end(),
                                          ( long_[ ref(nrows) = _1 ] >>
                                            long_[ ref(ncols) = _1 ] >>
                                            long_[ ref(nnz) = _1 ] ),
                                          space );

            if ( ! r )
                HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read", "missing data, expected <nrows> <ncols> <nnz>" );
            
            break;
        }// while

        HINFO( to_string( "(TMMMatrixIO) read : reading sparse %dx%d matrix with %d non-zero entries",
                          nrows, ncols, nnz ) );

        //
        // now, first go through file and count entries per row
        //

        const bool           is_sym   = (( mat_sym == MTX_SYM ) || ( mat_sym == MTX_HERM ));
        auto                 file_pos = in.tellg();
        std::vector< uint >  row_count( nrows, 0 );
        size_t               real_nnz = 0;
        
        for ( long  i = 0; i < nnz; i++ )
        {
            std::getline( in, line );

            uint        row, col;
            const bool  r = phrase_parse( line.begin(), line.end(),
                                          ( long_[ ref(row) = _1 ] >>
                                            long_[ ref(col) = _1 ] ),
                                          space );

            if ( ! r )
                HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read", "missing data, expected <nrows> <ncols>" );
            
            row--;
            col--;
            
            if ( is_sym && ( row != col ))
            {
                row_count[col]++;
                real_nnz++;
            }// if
            
            row_count[row]++;
            real_nnz++;
        }// for

        //
        // finally, reread the file and construct CRS data for the
        // sparse matrix with the previously collected data
        //

        auto  S   = std::make_unique< TSparseMatrix< value_t > >( nrows, ncols );
        uint  pos = 0;

        S->init( real_nnz );

        // first the row pointers
        for ( long  i = 0; i < nrows; i++ )
        {
            S->rowptr(i) = pos;
            pos         += row_count[i];
            row_count[i] = 0; // reset for later usage
        }// for
        S->rowptr(nrows) = pos;

        if ( pos != real_nnz )
            HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read",
                    "number of non-zero entries and counted entries differs" );
        
        in.seekg( file_pos );
        
        for ( long  i = 0; i < nnz; i++ )
        {
            long  row = 0, col = 0;
            
            std::getline( in, line );

            if ( mat_field == MTX_COMPLEX )
            {
                double      re = 0, im = 0;
                const bool  r  = phrase_parse( line.begin(), line.end(),
                                               ( long_[ ref(row) = _1 ] >>
                                                 long_[ ref(col) = _1 ] >>
                                                 double_[ ref(re) = _1 ] >>
                                                 double_[ ref(im) = _1 ] ),
                                               space );
                
                if ( ! r )
                    HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read",
                            "missing data, expected <row> <cols> <real> <imag>" );

                // adjust since starts counting at 1
                row--;
                col--;
                
                {
                    const idx_t    idx = S->rowptr(row) + row_count[row];
                
                    S->colind(idx) = col;
                    S->coeff(idx)  = get_value< value_t >::compose( re, im );
                }

                if ( is_sym && ( row != col ))
                {
                    const idx_t  idx = S->rowptr(col) + row_count[col];
                
                    S->colind(idx) = row;
                    S->coeff(idx)  = get_value< value_t >::compose( re, im );
                }// if
            }// if
            else if ( mat_field == MTX_REAL )
            {
                double      val = 0;
                const bool  r   = phrase_parse( line.begin(), line.end(),
                                                ( long_[ ref(row) = _1 ] >>
                                                  long_[ ref(col) = _1 ] >>
                                                  double_[ ref(val) = _1 ] ),
                                                space );
                
                if ( ! r )
                    HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read",
                            "missing data, expected <row> <cols> <coeff>" );
                
                // adjust since starts counting at 1
                row--;
                col--;

                {
                    const idx_t  idx = S->rowptr(row) + row_count[row];
                
                    S->colind(idx) = col;
                    S->coeff(idx)  = value_t(val);
                }

                if ( is_sym && ( row != col ))
                {
                    const idx_t  idx = S->rowptr(col) + row_count[col];
                
                    S->colind(idx) = row;
                    S->coeff(idx)  = value_t(val);
                }// if
            }// if
            else if ( mat_field == MTX_PATTERN )
            {
                const double  val = 1;
                const bool    r   = phrase_parse( line.begin(), line.end(),
                                                  ( long_[ ref(row) = _1 ] >>
                                                    long_[ ref(col) = _1 ] ),
                                                  space );
                
                if ( ! r )
                    HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read",
                            "missing data, expected <row> <cols> <coeff>" );
                
                // adjust since starts counting at 1
                row--;
                col--;

                {
                    const idx_t  idx = S->rowptr(row) + row_count[row];
                
                    S->colind(idx) = col;
                    S->coeff(idx)  = value_t(val);
                }
                
                if ( is_sym && ( row != col ))
                {
                    const idx_t  idx = S->rowptr(col) + row_count[col];
                
                    S->colind(idx) = row;
                    S->coeff(idx)  = value_t(val);
                }// if
            }// else
            
            row_count[row]++;
        }// for

        return std::unique_ptr< TMatrix< value_t > >( S.release() );
    }// if
    else
        HERROR( ERR_NOT_IMPL, "(TMMMatrixIO) read", "array format" );
}

#define INST_MM( type ) \
    template void TMMMatrixIO::write< type > ( const TMatrix< type > *, const std::string & ) const; \
    template std::unique_ptr< TMatrix< type > > TMMMatrixIO::read< type > ( const std::string & ) const;

INST_MM( float )
INST_MM( double )
INST_MM( std::complex< float > )
INST_MM( std::complex< double > )

variant_id_t
mtx_guess_value_type ( const std::string &  filename )
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "mtx_guess_value_type", filename );
    
    auto            in_ptr = open_read( filename );
    std::istream &  in     = * in_ptr.get();

    ///////////////////////////////////////////////////
    //
    // read header and determine format of file
    //

    std::string   line;
    auto          parts = std::vector< std::string >();
    mtxfield_t    mat_field  = MTX_REAL;

    while ( in.good() )
    {
        std::getline( in, line );
        split( line, " \t\r\n", parts );

        if (( parts[0] == "%%MatrixMarket" ) || ( parts[0] == "%MatrixMarket" ))
        {
            // look for qualifiers
            for ( uint i = 3; i < parts.size(); i++ )
            {
                if      ( parts[i] == "real"    ) mat_field = MTX_REAL;
                else if ( parts[i] == "integer" ) mat_field = MTX_INTEGER;
                else if ( parts[i] == "complex" ) mat_field = MTX_COMPLEX;
                else if ( parts[i] == "pattern" ) mat_field = MTX_PATTERN;
            }// for
            
            break;
        }// if
    }// while

    if ( mat_field == MTX_COMPLEX )
        return COMPLEX_FP64;
    else
        return REAL_FP64;
}

///////////////////////////////////////////////////
// 
// input and output in HDF5 format
//

#if HPRO_USE_HDF5 == 1
    
using namespace H5;

namespace
{

template < typename value_t >
void
h5_write_dense ( H5File *                         file,
                 const std::string &              gname,
                 const std::string &              mname,
                 const BLAS::Matrix< value_t > &  A )
{
    if ( is_complex_type< value_t >::value )
        HERROR( ERR_CONSISTENCY, "", "" );
    
    hsize_t    dims[2] = { A.ncols(), A.nrows() };
    DataSpace  dataspace( 2, dims );

    //
    // define datatype
    //

    std::unique_ptr< FloatType >  datatype;
    
    if ( is_single_prec< value_t >::value )
        datatype = std::make_unique< FloatType >( PredType::NATIVE_FLOAT );
    else
        datatype = std::make_unique< FloatType >( PredType::NATIVE_DOUBLE );

    // // little endian
    // datatype->setOrder( H5T_ORDER_LE );
            
    //
    // Create dataset for matrix data
    //
            
    DataSet dataset = file->createDataSet( gname + "/" + mname, *datatype, dataspace );
            
    //
    // write the data to the dataset using default memory space, file
    // space, and transfer properties.
    //

    if ( is_single_prec< value_t >::value )
        dataset.write( A.data(), PredType::NATIVE_FLOAT );
    else
        dataset.write( A.data(), PredType::NATIVE_DOUBLE );
}

template < typename value_t >
void
h5_write_dense ( H5File *                                         file,
                 const std::string &                              gname,
                 const std::string &                              mname,
                 const BLAS::Matrix< std::complex< value_t > > &  A )
{
    if ( is_complex_type< value_t >::value )
        HERROR( ERR_CONSISTENCY, "", "" );
    
    hsize_t    dims[2] = { A.ncols(), A.nrows() };
    DataSpace  dataspace( 2, dims );

    //
    // define datatype
    //

    CompType  datatype( sizeof( std::complex< value_t > ) );
    
    if ( is_single_prec< value_t >::value )
    {
        datatype.insertMember( "real", 0,               PredType::NATIVE_FLOAT );
        datatype.insertMember( "imag", sizeof(value_t), PredType::NATIVE_FLOAT );
    }// if
    else
    {
        datatype.insertMember( "real", 0,               PredType::NATIVE_DOUBLE );
        datatype.insertMember( "imag", sizeof(value_t), PredType::NATIVE_DOUBLE );
    }// else
    
    //
    // Create dataset for matrix data
    //
            
    DataSet dataset = file->createDataSet( gname + "/" + mname, datatype, dataspace );
            
    //
    // write the data to the dataset using default memory space, file
    // space, and transfer properties.
    //

    if ( is_single_prec< value_t >::value )
        dataset.write( A.data(), datatype );
    else
        dataset.write( A.data(), datatype );
}

//
// write structural part of matrix (type, indexsets, etc.)
//
template < typename value_t >
void
h5_write_struct ( H5File *                    file,
                  const std::string &         gname,
                  const TMatrix< value_t > *  M,
                  const std::string &         type )
{
    // local type for storing data
    struct  structure_t
    {
        char *  type;
        long    id;
        long    row_first, row_last;
        long    col_first, col_last;
    };

    structure_t  data;
    TIndexSet    row_is = M->row_is();
    TIndexSet    col_is = M->col_is();

    data.type      = const_cast< char * >( type.c_str() );
    data.id        = M->id();
    data.row_first = row_is.first();
    data.row_last  = row_is.last();
    data.col_first = col_is.first();
    data.col_last  = col_is.last();

    hid_t  str_type = H5Tcopy( H5T_C_S1 );
    
    H5Tset_size( str_type, H5T_VARIABLE );

    hsize_t    dims[] = { 1 };
    DataSpace  data_space( 1, dims );
    CompType   data_type( sizeof(structure_t) );
        
    data_type.insertMember( "type",      HOFFSET(structure_t, type),      str_type );
    data_type.insertMember( "id",        HOFFSET(structure_t, id),        PredType::NATIVE_LONG );
    data_type.insertMember( "row_first", HOFFSET(structure_t, row_first), PredType::NATIVE_LONG );
    data_type.insertMember( "row_last",  HOFFSET(structure_t, row_last),  PredType::NATIVE_LONG );
    data_type.insertMember( "col_first", HOFFSET(structure_t, col_first), PredType::NATIVE_LONG );
    data_type.insertMember( "col_last",  HOFFSET(structure_t, col_last),  PredType::NATIVE_LONG );

    auto  data_set = file->createDataSet( gname + "/structure", data_type, data_space );

    data_set.write( & data, data_type );

    data_space.close();
    data_type.close();
    data_set.close();
}

template < typename value_t >
void
h5_write ( H5File *                    file,
           const std::string &         gname,
           const TMatrix< value_t > *  M );

template < typename value_t >
void
h5_write ( H5File *                         file,
           const std::string &              gname,
           const TDenseMatrix< value_t > *  M )
{
    h5_write_struct( file, gname, M, "dense" );
    h5_write_dense( file, gname, "D", M->blas_mat() );
}

template < typename value_t >
void
h5_write ( H5File *                      file,
           const std::string &           gname,
           const TRkMatrix< value_t > *  M )
{
    h5_write_struct( file, gname, M, "lowrank" );

    h5_write_dense( file, gname, "U", M->blas_mat_A() );
    h5_write_dense( file, gname, "V", M->blas_mat_B() );
}

template < typename value_t >
void
h5_write ( H5File *                         file,
           const std::string &              gname,
           const TBlockMatrix< value_t > *  M )
{
    h5_write_struct( file, gname, M, "structured" );

    //
    // write sub block structure and sub block ids
    //

    struct  structure_t
    {
        long    nblock_rows, nblock_cols;
        char *  ids;
    };

    std::ostringstream  ids;

    for ( uint i = 0; i < M->nblock_rows(); ++i )
    {
        for ( uint j = 0; j < M->nblock_cols(); ++j )
        {
            ids << M->block(i,j)->id() << ",";
        }// for
    }// for
                
    structure_t  data;

    data.nblock_rows = M->nblock_rows();
    data.nblock_cols = M->nblock_rows();
    data.ids         = const_cast< char * >( ids.str().c_str() );

    hid_t  str_type = H5Tcopy( H5T_C_S1 );
    
    H5Tset_size( str_type, H5T_VARIABLE );

    hsize_t    dims[] = { 1 };
    DataSpace  data_space( 1, dims );
    CompType   data_type( sizeof(structure_t) );
        
    data_type.insertMember( "nblock_rows", HOFFSET(structure_t, nblock_rows), PredType::NATIVE_LONG );
    data_type.insertMember( "nblock_cols", HOFFSET(structure_t, nblock_cols), PredType::NATIVE_LONG );
    data_type.insertMember( "ids",         HOFFSET(structure_t, ids),         str_type );

    auto  data_set = file->createDataSet( gname + "/blockstruct", data_type, data_space );

    data_set.write( & data, data_type );

    data_space.close();
    data_type.close();
    data_set.close();
    
    //
    // write sub blocks
    //
    
    auto  group = std::make_unique< Group >( file->createGroup( gname + "/subblocks" ) );

    // std::cout << gname + "/subblocks" << std::endl;
    
    for ( uint i = 0; i < M->nblock_rows(); ++i )
    {
        for ( uint j = 0; j < M->nblock_cols(); ++j )
        {
            auto  sgname = gname + "/subblocks/" + to_string( "%d", M->block(i,j)->id() );

            h5_write( file, sgname, M->block(i,j) );
        }// for
    }// for
}

template < typename value_t >
void
h5_write ( H5File *                    file,
           const std::string &         gname,
           const TMatrix< value_t > *  M )
{
    // std::cout << gname << std::endl;
    
    auto  group = std::make_unique< Group >( file->createGroup( gname ) );
    
    if ( is_dense( M ) )
    {
        h5_write( file, gname, cptrcast( M, TDenseMatrix< value_t > ) );
    }// if
    else if ( is_lowrank( M ) )
    {
        h5_write( file, gname, cptrcast( M, TRkMatrix< value_t > ) );
    }// if
    else if ( is_blocked( M ) )
    {
        h5_write( file, gname, cptrcast( M, TBlockMatrix< value_t > ) );
    }// if
    else
        HERROR( ERR_MAT_TYPE, "(THDF5MatrixIO) h5_write", M->typestr() );
}

// template < typename value_t >
// const BLAS::Matrix< value_t >
// h5_read_dense ( H5File *             file,
//                 const std::string &  mname )
// {
//     if ( is_complex_type< value_t >::value )
//         HERROR( ERR_CONSISTENCY, "", "" );

//     auto  data_name = std::string( "" );
//     auto  iter_op   = [] ( H5::H5Object &      loc,
//                            const std::string   attr_name,
//                            const H5O_info_t *  oinfo, 
//                            void *              operator_data ) -> int
//     {
//         std::string *  dname = static_cast< std::string * >( operator_data );
        
//         // if ( mname != "" )
//         // {
//         //     // use given matrix name if in file
//         //     if ( mname == attr_name )
//         //         *dname = attr_name;
//         // }// if
//         // else
//         if ( attr_name != "." )
//         {
//             // use first name encountered
//             if ( *dname == "" )
//                 *dname = attr_name;
//             else
//             {
//                 // reset if not expected value
//                 if (( attr_name != *dname + "/type" ) &&
//                     ( attr_name != *dname + "/value" ))
//                     *dname = "";
//             }// else
//         }// if
            
//         // std::cout << attr_name << std::endl;
        
//         return 0;
//     };
    
//     file->visit( H5_INDEX_NAME, H5_ITER_INC, iter_op, & data_name, 0 );

//     // std::cout << "data : " << data_name << std::endl;
    
//     auto  dataset    = file->openDataSet( data_name + "/value" );
//     auto  type_class = dataset.getTypeClass();

//     // if ( type_class == H5T_FLOAT )
//     //     std::cout << "double" << std::endl;
//     // else
//     //     std::cout << type_class << std::endl;

//     auto  dataspace = dataset.getSpace();
//     auto  ndims     = dataspace.getSimpleExtentNdims();
//     auto  dims      = std::vector< hsize_t >( ndims );
    
//     dataspace.getSimpleExtentDims( dims.data() );

//     if ( ndims != 2 )
//         std::cout << "not a matrix" << std::endl;
                                                
//     // std::cout << dims[0] << " × " << dims[1] << std::endl;

//     auto  M = BLAS::Matrix< value_t >( dims[0], dims[1] );

//     dataset.read( M.data(), PredType::NATIVE_DOUBLE );
    
//     return  M;
// }

herr_t
visit_func ( hid_t               /* loc_id */,
             const char *        name,
             const H5O_info_t *  info,
             void *              operator_data )
{
    std::string *  dname = static_cast< std::string * >( operator_data );
    std::string    oname = name;

    if ( oname[0] != '.' )
    {
        if ( info->type == H5O_TYPE_GROUP )
        {
            // use first name encountered
            if ( *dname == "" )
                *dname = oname;
        }// if
        else if ( info->type == H5O_TYPE_DATASET )
        {
            if ( *dname != "" )
            {
                if ( oname == *dname + "/value" )     // actual dataset
                    *dname = *dname + "/value";
                else if ( oname != *dname + "/type" ) // just type info
                    *dname = "";
            }// if
            else
            {
                *dname = oname;                       // directly use dataset
            }// else
        }// if
        
        // switch ( info->type )
        // {
        //     case H5O_TYPE_GROUP:
        //         std::cout << oname << "  (Group)" << std::endl;
        //         break;
        //     case H5O_TYPE_DATASET:
        //         std::cout << oname << "  (Dataset)" << std::endl;
        //         break;
        //     case H5O_TYPE_NAMED_DATATYPE:
        //         std::cout << oname << "  (Datatype)" << std::endl;
        //         break;
        //     default:
        //         std::cout << oname << "  (Unknown)" << std::endl;
        // }// switch
    }// if
    
    return 0;
}

template < typename value_t >
const BLAS::Matrix< value_t >
h5_read_dense_c ( hid_t                file,
                  const std::string &  /* mname */ )
{
    if ( is_complex_type< value_t >::value )
        HERROR( ERR_CONSISTENCY, "", "" );

    auto  data_name = std::string( "" );
    
    auto  status    = H5Ovisit( file, H5_INDEX_NAME, H5_ITER_INC, visit_func, & data_name );

    if ( status != 0 )
        return BLAS::Matrix< value_t >( 0, 0 );

    // std::cout << "found: " << data_name << std::endl;
    
    // check if any compatible type found
    if ( data_name == "" )
        return BLAS::Matrix< value_t >( 0, 0 );

    // data_name = data_name + "/value";
    
    auto  dataset   = H5Dopen( file, data_name.c_str(), H5P_DEFAULT );
    auto  dataspace = H5Dget_space( dataset );
    auto  ndims     = H5Sget_simple_extent_ndims( dataspace );
    auto  dims      = std::vector< hsize_t >( ndims );

    if ( ndims != 2 )
        std::cout << "not a matrix" << std::endl;
                                                
    H5Sget_simple_extent_dims( dataspace, dims.data(), nullptr );
    
    // std::cout << dims[0] << " × " << dims[1] << std::endl;

    auto  M = BLAS::Matrix< value_t >( dims[1], dims[0] );

    status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, M.data() );
    
    if ( status != 0 )
        return BLAS::Matrix< value_t >( 0, 0 );
    
    return  M;
}

}// namespace anonymous

#endif

//
// write matrix \a A to file \a fname
//
template < typename value_t >
void
THDF5MatrixIO::write ( const TMatrix< value_t > *  A,
                       const std::string &         fname ) const
{
    write( A, fname, "M" );
}

//
// write matrix \a A with name \a mname to file \a fname
//
#if HPRO_USE_HDF5 == 1
template < typename value_t >
void
THDF5MatrixIO::write ( const TMatrix< value_t > *  A,
                       const std::string &         fname,
                       const std::string &         mname ) const
{
    try
    {
        H5File  file( fname, H5F_ACC_TRUNC );

        if ( A == nullptr )
            return;

        h5_write( & file, "/" + mname, A );
    }// try
    catch( FileIException &      error ) { error.printErrorStack(); }
    catch( DataSetIException &   error ) { error.printErrorStack(); }
    catch( DataSpaceIException & error ) { error.printErrorStack(); }
    catch( DataTypeIException &  error ) { error.printErrorStack(); }
}
#else
template < typename value_t >
void
THDF5MatrixIO::write ( const TMatrix< value_t > * ,
                       const std::string &  ,
                       const std::string &   ) const
{
    HERROR( ERR_NOHDF5, "(THDF5MatrixIO) write", "" );
}
#endif


#if HPRO_USE_HDF5 == 1
template < typename value_t >
void
THDF5MatrixIO::write ( const BLAS::Matrix< value_t > &  A,
                       const std::string &              fname,
                       const std::string &              mname ) const
{
    try
    {
        H5File  file( fname, H5F_ACC_TRUNC );

        h5_write_dense( & file, "/", mname, A );
    }// try
    catch( FileIException &      error ) { error.printErrorStack(); }
    catch( DataSetIException &   error ) { error.printErrorStack(); }
    catch( DataSpaceIException & error ) { error.printErrorStack(); }
    catch( DataTypeIException &  error ) { error.printErrorStack(); }
}
#else
template < typename value_t >
void
THDF5MatrixIO::write ( const BLAS::Matrix< value_t > &,
                       const std::string &            ,
                       const std::string &            ) const
{
    HERROR( ERR_NOHDF5, "(THDF5MatrixIO) write", "" );
}
#endif

//
// read and return matrix from file \a fname
//
#if HPRO_USE_HDF5 == 1
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
THDF5MatrixIO::read  ( const std::string &  filename ) const
{
    try
    {
        auto  file = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
        auto  M    = h5_read_dense_c< value_t >( file, "" );

        return std::make_unique< TDenseMatrix< value_t > >( is( 0, M.nrows()-1 ), is( 0, M.ncols()-1 ), std::move( M ) );
    }// try
    catch( FileIException &      error ) { error.printErrorStack(); }
    catch( DataSetIException &   error ) { error.printErrorStack(); }
    catch( DataSpaceIException & error ) { error.printErrorStack(); }
    catch( DataTypeIException &  error ) { error.printErrorStack(); }

    return std::unique_ptr< TMatrix< value_t > >();
}
#else
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
THDF5MatrixIO::read  ( const std::string &  ) const
{
    HERROR( ERR_NOHDF5, "(THDF5MatrixIO) read", "" );

    return std::unique_ptr< TMatrix< value_t > >();
}
#endif

#define INST_HDF5( type ) \
    template void THDF5MatrixIO::write< type > ( const TMatrix< type > *, const std::string & ) const; \
    template void THDF5MatrixIO::write< type > ( const BLAS::Matrix< type > &, const std::string &, const std::string & ) const; \
    template std::unique_ptr< TMatrix< type > > THDF5MatrixIO::read< type > ( const std::string & ) const;

INST_HDF5( float )
INST_HDF5( double )
INST_HDF5( std::complex< float > )
INST_HDF5( std::complex< double > )

///////////////////////////////////////////////////
// 
// input and output in H2Lib format
//

//
// write matrix \a A to file \a fname
//
template < typename value_t >
void
TH2LibMatrixIO::write ( const TMatrix< value_t > *  /* A */,
                        const std::string &         /* filename */ ) const
{
    #if HPRO_USE_NETCDF == 0
    
    HERROR( ERR_NONETCDF, "(TH2LibMatrixIO) write", "" );

    #else

    HERROR( ERR_NOT_IMPL, "(TH2LibMatrixIO) write", "" );

    #endif
}

#if HPRO_USE_NETCDF == 1

#define NETCDF_CHECK( call, args, func )                                \
    {                                                                   \
        auto  ncres = call args ;                                       \
        if ( ncres != NC_NOERR )                                        \
            HERROR( ERR_NETCDF, func, to_string( "error %d during %s", ncres, #call ) ); \
    }

namespace
{

std::unique_ptr< TCluster >
netcdf_read_cl ( int       nc_file,
                 int       nc_sons,
                 int       nc_size,
                 int       nc_coeff,
                 int       dim,
                 idx_t     nidx,
                 size_t &  nclusters,
                 size_t &  ncoeffs )
{
    size_t     start  = nclusters;
    size_t     count  = 1;
    ptrdiff_t  stride = 1;

    // number of sons
    uint  nsons = 0;
    
    NETCDF_CHECK( nc_get_vars, ( nc_file, nc_sons, &start, &count, &stride, &nsons ), "netcdf_read_cl" );

    // size of cluster
    uint  size = 0;
    
    NETCDF_CHECK( nc_get_vars, ( nc_file, nc_size, &start, &count, &stride, &size ), "netcdf_read_cl" );

    // create new cluster
    auto  cl = std::make_unique< TCluster >( nidx, nidx + size - 1 );

    cl->set_nsons( nsons );
    
    nclusters++;

    // read sons
    if ( nsons > 0 )
    {
        auto  son_idx = nidx;
        
        for ( uint  i = 0; i < nsons; i++ )
        {
            auto  son = netcdf_read_cl( nc_file, nc_sons, nc_size, nc_coeff, dim, son_idx, nclusters, ncoeffs );

            son_idx += son->size();
            cl->set_son( i, son.release() );
        }// for
    }// if

    // read bounding box
    auto  bb = std::vector< double >( dim );
    
    start = ncoeffs;
    count = dim;
    ncoeffs += count;
    NETCDF_CHECK( nc_get_vars, ( nc_file, nc_coeff, &start, &count, &stride, bb.data() ), "netcdf_read_cl" );

    start += dim;
    ncoeffs += count;
    NETCDF_CHECK( nc_get_vars, ( nc_file, nc_coeff, &start, &count, &stride, bb.data() ), "netcdf_read_cl" );

    return  cl;
}

std::unique_ptr< TCluster >
netcdf_read_cl ( int                  nc_file,
                 const std::string &  prefix )
{
    auto    name = std::string();
    char    dimname[ NC_MAX_NAME + 1 ];

    // clusters dimension
    size_t  clusters    = 0;
    int     nc_clusters = 0;
    
    name = prefix + "_clusters";
    NETCDF_CHECK( nc_inq_dimid, ( nc_file, name.c_str(), &nc_clusters ), "netcdf_read_cl" );
    NETCDF_CHECK( nc_inq_dim,   ( nc_file, nc_clusters, dimname, & clusters ), "netcdf_read_cl" );

    // coeffs dimension
    size_t  coeffs    = 0;
    int     nc_coeffs = 0;
    
    name = prefix + "_coeffs";
    NETCDF_CHECK( nc_inq_dimid, ( nc_file, name.c_str(), &nc_coeffs ), "netcdf_read_cl" );
    NETCDF_CHECK( nc_inq_dim,   ( nc_file, nc_coeffs, dimname, &coeffs ), "netcdf_read_cl" );

    // totalsize dimension
    size_t  totalsize    = 0;
    int     nc_totalsize = 0;
    
    name = prefix + "_totalsize";
    NETCDF_CHECK( nc_inq_dimid, ( nc_file, name.c_str(), &nc_totalsize ), "netcdf_read_cl" );
    NETCDF_CHECK( nc_inq_dim,   ( nc_file, nc_totalsize, dimname, &totalsize ), "netcdf_read_cl" );

    // dim dimension
    size_t  dim    = 0;
    int     nc_dim = 0;
    
    name = prefix + "_dim";
    NETCDF_CHECK( nc_inq_dimid, ( nc_file, name.c_str(), &nc_dim ), "netcdf_read_cl" );
    NETCDF_CHECK( nc_inq_dim,   ( nc_file, nc_dim, dimname, &dim ), "netcdf_read_cl" );

    // sons variable
    int  nc_sons = 0;
    
    name = prefix + "_sons";
    NETCDF_CHECK( nc_inq_varid, ( nc_file, name.c_str(), &nc_sons ), "netcdf_read_cl" );

    // size variable
    int  nc_size = 0;
    
    name = prefix + "_size";
    NETCDF_CHECK( nc_inq_varid, ( nc_file, name.c_str(), &nc_size ), "netcdf_read_cl" );

    // idx variable
    int  nc_idx = 0;
    
    name = prefix + "_idx";
    NETCDF_CHECK( nc_inq_varid, ( nc_file, name.c_str(), &nc_idx ), "netcdf_read_cl" );

    // coeff variable
    int  nc_coeff = 0;
    
    name = prefix + "_coeff";
    NETCDF_CHECK( nc_inq_varid, ( nc_file, name.c_str(), &nc_coeff ), "netcdf_read_cl" );

    // read index
    auto  idx = std::vector< uint >( totalsize );
    
    nc_get_var( nc_file, nc_idx, idx.data() );

    // read coefficients from NetCDF variables
    idx_t   nidx      = 0;
    size_t  nclusters = 0;
    size_t  ncoeffs   = 0;
    
    auto  cl = netcdf_read_cl( nc_file, nc_sons, nc_size, nc_coeff, dim, nidx, nclusters, ncoeffs );
    
    if ( nclusters != clusters )
        HERROR( ERR_NETCDF, "netcdf_read_cl", "inconsistent number of cluster indices" );
    
    if ( ncoeffs != coeffs )
        HERROR( ERR_NETCDF, "netcdf_read_cl", "inconsistent number of cluster coefficients" );

    return  cl;
}

template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
netcdf_read_mat ( int               nc_file,
                  int               nc_type,
                  int               nc_rows,
                  int               nc_cols,
                  int               nc_rank,
                  int               nc_rsons,
                  int               nc_csons,
                  int               nc_coeffs,
                  size_t &          blockidx,
                  size_t &          coeffidx,
                  const TCluster *  rowcl,
                  const TCluster *  colcl )
{
    size_t     start  = blockidx;
    size_t     count  = 1;
    ptrdiff_t  stride = 1;
    int  val = 0;
    
    // number of rows
    NETCDF_CHECK( nc_get_vars, ( nc_file, nc_rows, &start, &count, &stride, &val ), "netcdf_read_mat" );

    const uint  nrows = val;

    // number of columns
    NETCDF_CHECK( nc_get_vars, ( nc_file, nc_cols, &start, &count, &stride, &val ), "netcdf_read_mat" );

    const uint  ncols = val;

    // type of matrix
    NETCDF_CHECK( nc_get_vars, ( nc_file, nc_type, &start, &count, &stride, &val ), "netcdf_read_mat" );

    const uint  type = val;

    //
    // create matrix based on given type
    //

    auto  M = std::unique_ptr< TMatrix< value_t > >();
    
    if ( type == 1 )
    {
        //
        // dense matrix
        //

        auto  D = BLAS::Matrix< value_t >( nrows, ncols );

        start     = coeffidx;
        count     = nrows * ncols;
        coeffidx += count;
        NETCDF_CHECK( nc_get_vars, ( nc_file, nc_coeffs, &start, &count, &stride, D.data() ), "netcdf_read_mat" );

        // next block
        blockidx++;

        M = std::make_unique< TDenseMatrix< value_t > >( *rowcl, *colcl, std::move( D ) );
    }// if
    else if ( type == 2 )
    {
        //
        // lowrank matrix
        //

        // rank
        NETCDF_CHECK( nc_get_vars, ( nc_file, nc_rank, &start, &count, &stride, &val ), "netcdf_read_mat" );

        const uint  rank = val;

        // read U
        auto  U = BLAS::Matrix< value_t >( nrows, rank );
        
        start     = coeffidx;
        count     = nrows * rank;
        coeffidx += count;
        NETCDF_CHECK( nc_get_vars, ( nc_file, nc_coeffs, &start, &count, &stride, U.data() ), "netcdf_read_mat" );

        // read V
        auto  V = BLAS::Matrix< value_t >( ncols, rank );
        
        start     = coeffidx;
        count     = ncols * rank;
        coeffidx += count;
        NETCDF_CHECK( nc_get_vars, ( nc_file, nc_coeffs, &start, &count, &stride, V.data() ), "netcdf_read_mat" );

        // next block
        blockidx++;

        M = std::make_unique< TRkMatrix< value_t > >( *rowcl, *colcl, std::move( U ), std::move( V ) );
    }// if
    else if ( type == 3 )
    {
        //
        // structured matrix
        //

        // number of row sons
        NETCDF_CHECK( nc_get_vars, ( nc_file, nc_rsons, &start, &count, &stride, &val ), "netcdf_read_mat" );

        const uint  nrsons = val;

        if ( nrsons != rowcl->nsons() )
            HERROR( ERR_MAT_STRUCT, "netcdf_read_mat", "matrix size incompatible with cluster" );
        
        // number of column sons
        NETCDF_CHECK( nc_get_vars, ( nc_file, nc_csons, &start, &count, &stride, &val ), "netcdf_read_mat" );

        const uint  ncsons = val;

        if ( ncsons != colcl->nsons() )
            HERROR( ERR_MAT_STRUCT, "netcdf_read_mat", "matrix size incompatible with cluster" );

        //
        // create matrix and recursively read sub blocks
        //

        auto  B = std::make_unique< TBlockMatrix< value_t > >( *rowcl, *colcl );

        B->set_block_struct( nrsons, ncsons );

        // next block
        blockidx++;
        
        for ( uint  j = 0; j < ncsons; ++j )
        {
            for ( uint  i = 0; i < nrsons; ++i )
            {
                auto  B_ij = netcdf_read_mat< value_t >( nc_file, nc_type, nc_rows, nc_cols, nc_rank,
                                                         nc_rsons, nc_csons, nc_coeffs,
                                                         blockidx, coeffidx, rowcl->son(i), colcl->son(j) );

                B->set_block( i, j, B_ij.release() );
            }// for
        }// for

        M = std::move( B );
    }// if
    else
        HERROR( ERR_NETCDF, "netcdf_read_mat", to_string( "unknown matrix type %d", type ) );

    M->set_id( blockidx-1 );
    
    return M;
}

template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
netcdf_read_mat ( int                  nc_file,
                  const std::string &  prefix,
                  const TCluster *     rowcl,
                  const TCluster *     colcl )
{
    auto    name = std::string();
    char    dimname[ NC_MAX_NAME + 1 ];

    // blocks
    size_t  blocks    = 0;
    int     nc_blocks = 0;
    
    name = prefix + "_blocks";
    NETCDF_CHECK( nc_inq_dimid, ( nc_file, name.c_str(), &nc_blocks ), "netcdf_read_mat" );
    NETCDF_CHECK( nc_inq_dim,   ( nc_file, nc_blocks, dimname, & blocks ), "netcdf_read_mat" );

    // coeffs
    size_t  coeffs    = 0;
    int     nc_coeffs = 0;
    
    name = prefix + "_coeffs";
    NETCDF_CHECK( nc_inq_dimid, ( nc_file, name.c_str(), &nc_coeffs ), "netcdf_read_mat" );
    NETCDF_CHECK( nc_inq_dim,   ( nc_file, nc_coeffs, dimname, &coeffs ), "netcdf_read_mat" );

    // type
    int  nc_type = 0;
    
    name = prefix + "_type";
    NETCDF_CHECK( nc_inq_varid, ( nc_file, name.c_str(), &nc_type ), "netcdf_read_mat" );

    // rows
    int  nc_rows = 0;
    
    name = prefix + "_rows";
    NETCDF_CHECK( nc_inq_varid, ( nc_file, name.c_str(), &nc_rows ), "netcdf_read_mat" );

    // cols
    int  nc_cols = 0;
    
    name = prefix + "_cols";
    NETCDF_CHECK( nc_inq_varid, ( nc_file, name.c_str(), &nc_cols ), "netcdf_read_mat" );

    // rank
    int  nc_rank = 0;
    
    name = prefix + "_rank";
    NETCDF_CHECK( nc_inq_varid, ( nc_file, name.c_str(), &nc_rank ), "netcdf_read_mat" );

    // rsons
    int  nc_rsons = 0;
    
    name = prefix + "_rsons";
    NETCDF_CHECK( nc_inq_varid, ( nc_file, name.c_str(), &nc_rsons ), "netcdf_read_mat" );

    // csons
    int  nc_csons = 0;
    
    name = prefix + "_csons";
    NETCDF_CHECK( nc_inq_varid, ( nc_file, name.c_str(), &nc_csons ), "netcdf_read_mat" );

    // coeff variable
    int  nc_coeff = 0;
    
    name = prefix + "_coeff";
    NETCDF_CHECK( nc_inq_varid, ( nc_file, name.c_str(), &nc_coeff ), "netcdf_read_mat" );

    //
    // read actual matrix
    //

    size_t  blockidx = 0;
    size_t  coeffidx = 0;
    
    auto  M = netcdf_read_mat< value_t >( nc_file, nc_type, nc_rows, nc_cols, nc_rank,
                                          nc_rsons, nc_csons, nc_coeff,
                                          blockidx, coeffidx, rowcl, colcl );
    
    if ( blocks != blockidx )
        HERROR( ERR_NETCDF, "netcdf_read_mat", "inconsistent number of blocks" );
    
    if ( coeffs != coeffidx )
        HERROR( ERR_NETCDF, "netcdf_read_mat", "inconsistent number of coefficients" );

    return  M;
}

}// namespace anonymous

#endif

template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TH2LibMatrixIO::read  ( const std::string &  filename ) const
{
    #if HPRO_USE_NETCDF == 0
    
    HERROR( ERR_NONETCDF, "(TH2LibMatrixIO) read", "" );

    return std::unique_ptr< TMatrix< value_t > >();

    #else

    int  nc_file = 0;
    int  nc_sons = 0;
    
    NETCDF_CHECK( nc_open, ( filename.c_str(), NC_NOWRITE, &nc_file ), "(TH2LibMatrixIO) read" );

    auto  rct = netcdf_read_cl( nc_file, "rc" );
    auto  cct = decltype( rct )();
    
    if ( nc_inq_varid( nc_file, "cc_sons", &nc_sons ) == NC_NOERR )
        cct = netcdf_read_cl( nc_file, "cc" );
    else
        cct = std::unique_ptr< TCluster >( rct->copy() );

    auto  M = netcdf_read_mat< value_t >( nc_file, "ma", rct.get(), cct.get() );

    return  M;
    
    #endif
}

#define INST_H2Lib( type ) \
    template void TH2LibMatrixIO::write< type > ( const TMatrix< type > *, const std::string & ) const; \
    template std::unique_ptr< TMatrix< type > > TH2LibMatrixIO::read< type > ( const std::string & ) const;

INST_H2Lib( float )
INST_H2Lib( double )
INST_H2Lib( std::complex< float > )
INST_H2Lib( std::complex< double > )

}// namespace Hpro
