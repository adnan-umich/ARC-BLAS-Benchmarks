//
// Project     : HLIBpro
// File        : HBIO.cc
// Description : classes for input/output in Harwll/Boeing format
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>
#include <string>
#include <memory>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

#include "baseio.hh"

#include "hpro/io/TVectorIO.hh"
#include "hpro/io/TMatrixIO.hh"

#include "hpro/matrix/TSparseMatrix.hh"
#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/structure.hh"

#include "hpro/vector/TScalarVector.hh"

namespace Hpro
{

namespace fs = boost::filesystem;

using boost::algorithm::to_upper;

///////////////////////////////////////////////////
// 
// input and output in Harwell-Boeing and the
// Rutherford-Boeing format
//

namespace
{

using boost::spirit::qi::char_;
using boost::spirit::qi::uint_;
using boost::spirit::qi::_1;
using boost::spirit::qi::phrase_parse;
using boost::spirit::ascii::space;
using boost::phoenix::ref;
    
//
// parse Fortran integer and float formats to determine
// number and lenght of corresponding fields
//
void
parse_int_fmt ( const std::string &  fmt,
                uint &               elem_len,
                uint &               n_elem )
{
    uint        n = 0, l = 0;
    const bool  r = phrase_parse( fmt.begin(), fmt.end(),
                                  ( -char_('(') >>
                                    uint_[ ref(n) = _1 ] >>
                                    // use | operator because of problems with string syntax and optimization
                                    ( char_('i') | char_('I') ) >>
                                    uint_[ ref(l) = _1 ] >>
                                    -char_(')') ),
                                  space );

    if ( ! r )
        HERROR( ERR_FMT_HB, "parse_int_fmt", "invalid integer format" );

    n_elem   = n;
    elem_len = l;
}

void
parse_flt_fmt ( const std::string &  fmt,
                uint &               elem_len,
                uint &               n_elem )
{
    uint        n = 0, l = 0;
    const bool  r = phrase_parse( fmt.begin(), fmt.end(),
                                  ( -char_('(') >>
                                    uint_[ ref(n) = _1 ] >>
                                    // use | operator because of problems with string syntax and optimization
                                    ( char_('d') | char_('e') | char_('f') | char_('D') | char_('E') | char_('F') ) >>
                                    uint_[ ref(l) = _1 ] >>
                                    char_('.') >>
                                    uint_ >>
                                    -char_(')') ),
                                  space );

    if ( ! r )
        HERROR( ERR_FMT_HB, "parse_flt_fmt", "invalid float format" );

    n_elem   = n;
    elem_len = l;
}

//
// special routines to write ints and floats with equal
// width in file
//
void
hb_write_int ( std::ostream &  out, const uint  val )
{
    out << boost::format( " %-9d" ) % val;
}

void
hb_write_flt ( std::ostream &  out, const double  val )
{
    if ( val < 0.0 )
        out << boost::format( " %-.16E" ) % val;
    else
        out << boost::format( "  %-.16E" ) % val;
}

//
// return number of non-zeroes in lower left part plus diagonal
//
template < typename value_t >
size_t
lower_left_nonzeroes ( const TSparseMatrix< value_t > *  S )
{
    if ( S == nullptr )
        return 0;

    size_t  llnnz = 0;

    for ( idx_t  row = 0; row < idx_t(S->rows()); row++ )
    {
        const idx_t  lb = S->rowptr(row);
        const idx_t  ub = S->rowptr(row+1);

        for ( idx_t  j = lb; j < ub; ++j )
            if ( S->colind(j) <= row )
                llnnz++;
    }// for

    return llnnz;
}

}// namespace anonymous

//
// write matrix with name
//
template < typename value_t >
void
THBMatrixIO::write ( const TMatrix< value_t > *  A,
                     const std::string &         fname ) const
{
    if ( IS_TYPE( A, TSparseMatrix ) )
    {
        //
        // convert CRS format into CCS format
        //

        auto           M      = cptrcast( A, TSparseMatrix< value_t > );
        const size_t   nrows  = M->rows();
        const size_t   ncols  = M->cols();
        const bool     is_sym = M->is_symmetric() || M->is_hermitian();
        const size_t   nnz    = ( is_sym ? lower_left_nonzeroes( M ) : M->n_non_zero() );
        auto           colptr = std::vector< idx_t >( ncols + 1, 0 );

        // count number of entries per column
        for ( idx_t  row = 0; row < idx_t(nrows); row++ )
        {
            const idx_t  lb = M->rowptr(idx_t(row));
            const idx_t  ub = M->rowptr(idx_t(row)+1);

            if ( is_sym )
            {
                for ( idx_t  j = lb; j < ub; j++ )
                    if ( M->colind(j) <= row ) // count only lower left part
                        colptr[ M->colind(j) ]++;
            }// if
            else
            {
                for ( idx_t  j = lb; j < ub; j++ )
                    colptr[ M->colind(j) ]++;
            }// else
        }// for
        
        // compute real colptr array
        idx_t  pos = 0;
    
        for ( size_t  j = 0; j < ncols; j++ )
        {
            const idx_t  n = colptr[j];
        
            colptr[j] = pos;
            pos      += n;
        }// for
        colptr[ncols] = pos;

        // set up rowind and coeff arrays
        std::vector< idx_t >    rowind( nnz, 0 );
        std::vector< value_t >  coeffs( nnz, value_t(0) );
        std::vector< idx_t >    colpos( ncols, 0 );

        for ( idx_t  row = 0; row < idx_t(nrows); row++ )
        {
            const idx_t  lb = M->rowptr(row);
            const idx_t  ub = M->rowptr(row+1);

            if ( is_sym )
            {
                for ( idx_t  j = lb; j < ub; j++ )
                {
                    const idx_t  col = M->colind(j);
                    const idx_t  idx = colptr[col] + colpos[col];

                    if ( col > row ) // only lower left part
                        continue;
                    
                    rowind[ idx ] = row;
                    coeffs[ idx ] = M->coeff( j );
                    
                    colpos[col]++;
                }// for
            }// if
            else
            {
                for ( idx_t  j = lb; j < ub; j++ )
                {
                    const idx_t  col = M->colind(j);
                    const idx_t  idx = colptr[col] + colpos[col];
            
                    rowind[ idx ] = row;
                    coeffs[ idx ] = M->coeff( j );
                    
                    colpos[col]++;
                }// for
            }// else
        }// for
    
        //
        // finally, write matrix to file
        //
    
        auto            out_ptr = open_write( fname );
        std::ostream &  out     = * out_ptr.get();

        // title and key
        out << to_string( "%-72s%-8s", "Created by HLIBpro", "" ) << std::endl;

        // line with number of lines for each type
        const uint  ptrcrd = uint( (ncols+1) / 7 + ( (ncols+1) % 7 == 0 ? 0 : 1 ) ); // 7 colptr per line plus rest
        const uint  indcrd = uint( nnz / 7 + ( nnz % 7 == 0 ? 0 : 1 ) );             // 7 rowind per line plus rest
        const uint  nvals  = uint( is_complex_type< value_t >::value ? 2*nnz : nnz );
        const uint  valcrd = nvals / 3 + ( nvals % 3 == 0 ? 0 : 1 );         // 3 values per line plus rest
        const uint  rhscrd = 0;                                              // no right hand sides
        const uint  totcrd = ptrcrd + indcrd + valcrd + rhscrd;

        out << to_string( "%-14d%-14d%-14d%-14d", totcrd, ptrcrd, indcrd, valcrd, rhscrd ) << std::endl;

        // matrix type, dimension of matrix, etc.
        char  mxtype[3];

        if ( is_complex_type< value_t >::value )
            mxtype[0] = 'C';
        else
            mxtype[0] = 'R';

        if      ( M->is_symmetric() ) mxtype[1] = 'S';
        else if ( M->is_hermitian() ) mxtype[1] = 'H';
        else                          mxtype[1] = 'U';

        mxtype[2] = 'A';

        out << mxtype[0] << mxtype[1] << mxtype[2]
            << boost::format( "           %-14d%-14d%-14d%-14d" ) % nrows % ncols % nnz % 0 << std::endl;

        // pointer, index and value format
        out << "(7i10)          (7i10)          (3e24.16)           (0e0.0)" << std::endl;

        // write column pointers
        pos = 0;
        for ( uint i = 0; i <= ncols; i++ )
        {
            hb_write_int( out, uint(colptr[i])+1 );
            if ( (++pos) % 7 == 0 )
                out << std::endl;
        }// for
    
        if ( pos % 7 != 0 )
            out << std::endl;

        // write row indices
        pos = 0;
        for ( uint i = 0; i < nnz; i++ )
        {
            hb_write_int( out, uint(rowind[i])+1 );
            if ( (++pos) % 7 == 0 )
                out << std::endl;
        }// for
    
        if ( pos % 7 != 0 )
            out << std::endl;

        // write coeffs
        if ( is_complex_type< value_t >::value )
        {
            pos = 0;
        
            for ( uint i = 0; i < nnz; i++ )
            {
                hb_write_flt( out, std::real( coeffs[i] ) );
                if ( (++pos) % 3 == 0 ) out << std::endl;
                hb_write_flt( out, std::imag( coeffs[i] ) );
                if ( (++pos) % 3 == 0 ) out << std::endl;
            }// for
        }// if
        else
        {
            pos = 0;
        
            for ( uint i = 0; i < nnz; i++ )
            {
                hb_write_flt( out, std::real( coeffs[i] ) );
                if ( (++pos) % 3 == 0 ) out << std::endl;
            }// for
        }// if
    }// if
    else if ( is_dense( A ) )
    {
        //
        // convert to dense and write all entries
        // (even if zero)
        //

        auto          M        = cptrcast( A, TDenseMatrix< value_t > );
        const size_t  nrows    = M->rows();
        const size_t  ncols    = M->cols();
        const bool    is_sym   = M->is_symmetric() || M->is_hermitian();
        const size_t  nnz_diag = std::min( nrows, ncols );
        const size_t  nnz      = ( is_sym ?
                                   ((nrows * ncols) - nnz_diag) / 2 + nnz_diag :
                                   nrows * ncols );
        
        auto            out_ptr = open_write( fname );
        std::ostream &  out     = * out_ptr.get();
        
        // title and key
        out << to_string( "%-72s%-8s", "Created by HLIBpro", "" ) << std::endl;

        // line with number of lines for each type
        const uint  ptrcrd = uint( (ncols+1) / 7 + ( (ncols+1) % 7 == 0 ? 0 : 1 ) ); // 7 colptr per line plus rest
        const uint  indcrd = uint( nrows / 7 + ( nrows % 7 == 0 ? 0 : 1 ) );         // 7 rowind per line plus rest
        const uint  nvals  = uint( is_complex_type< value_t >::value ? 2*nnz : nnz );
        const uint  valcrd = nvals / 3 + ( nvals % 3 == 0 ? 0 : 1 ); // 3 values per line plus rest
        const uint  rhscrd = 0;                                      // no right hand sides
        const uint  totcrd = ptrcrd + indcrd + valcrd + rhscrd;

        out << to_string( "%-14d%-14d%-14d%-14d", totcrd, ptrcrd, indcrd, valcrd, rhscrd ) << std::endl;

        // matrix type, dimension of matrix, etc.
        std::string  mxtype( 4, '\0' );

        if ( is_complex_type< value_t >::value )
            mxtype[0] = 'C';
        else
            mxtype[0] = 'R';

        if      ( M->is_symmetric() ) mxtype[1] = 'S';
        else if ( M->is_hermitian() ) mxtype[1] = 'H';
        else                          mxtype[1] = 'U';

        mxtype[2] = 'A';

        out << mxtype << boost::format( "            %-14d%-14d%-14d%-14d" ) % nrows % ncols % nnz % 0 << std::endl;

        // pointer, index and value format
        out << "(7i10)          (7i10)          (3e24.16)           (0e0.0)" << std::endl;

        // write column pointers
        uint  pos = 0;
        uint  cnt = 1;

        for ( uint i = 0; i < ncols; i++ )
        {
            hb_write_int( out, cnt );
            cnt += uint( is_sym ? nrows - i : nrows );
            if ( (++pos) % 7 == 0 )
                out << std::endl;
        }// for
        out << boost::format( " %-9d" ) % nnz << std::endl;

        // write row indices
        pos = 0;
        for ( uint j = 0; j < ncols; j++ )
            for ( uint i = (is_sym ? j : 0); i < nrows; i++ )
            {
                hb_write_int( out, i );
                if ( (++pos) % 7 == 0 )
                    out << std::endl;
            }// for
    
        if ( pos % 7 != 0 )
            out << std::endl;

        // write coeffs
        if ( is_complex_type< value_t >::value )
        {
            pos = 0;
        
            for ( uint j = 0; j < ncols; j++ )
                for ( uint i = (is_sym ? j : 0); i < nrows; i++ )
                {
                    const auto val = M->entry( i, j );
                    
                    hb_write_flt( out, std::real( val ) );
                    if ( (++pos) % 3 == 0 ) out << std::endl;
                    hb_write_flt( out, std::imag( val ) );
                    if ( (++pos) % 3 == 0 ) out << std::endl;
                }// for
        }// if
        else
        {
            pos = 0;
        
            for ( uint j = 0; j < ncols; j++ )
                for ( uint i = (is_sym ? j : 0); i < nrows; i++ )
                {
                    hb_write_flt( out, std::real( M->entry( i, j ) ) );
                    if ( (++pos) % 3 == 0 ) out << std::endl;
                }// for
        }// if
    }// else
}

//
// read matrix from file
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
THBMatrixIO::read  ( const std::string &  fname ) const
{
    using  real_t = real_type_t< value_t >;
    
    auto            in_ptr = open_read( fname );
    std::istream &  in     = * in_ptr.get();
    
    // get title and key
    std::string  line;
    std::string  title( 73, ' ' ), key( 9, ' ' );
    
    std::getline( in, line );
    title = line.substr( 0, 72 );
    key   = line.substr( 72, 8 );
    
    // retrieve total number of lines, number of lines for column pointers,
    // number of lines for row indices, number of lines for coefficients,
    // number of lines for right-hand sides
    // UNUSED : uint  totcrd;
    uint  ptrcrd, indcrd, valcrd, rhscrd;

    std::getline( in, line );
    // UNUSED : totcrd = str_to_int( line.substr(  0, 14 ) );
    std::string  tmp;

    ptrcrd = str_to_int( line.substr( 14, 14 ) );
    indcrd = str_to_int( line.substr( 28, 14 ) );
    valcrd = str_to_int( line.substr( 42, 14 ) );
    rhscrd = str_to_int( line.substr( 56, 14 ) );
    
    // retrieve matrix type, number of rows, number of columns,
    // number of non-zeroes, number of elemental matrix entries,
    std::string  mxtype;
    uint         nrows, ncols, nnzero;
    // UNUSED : uint    neltvl;
    
    std::getline( in, line );
    mxtype = line.substr( 0, 3 );
    nrows  = str_to_int( line.substr( 14, 14 ) );
    ncols  = str_to_int( line.substr( 28, 14 ) );
    nnzero = str_to_int( line.substr( 42, 14 ) );
    // UNUSED : neltvl = str_to_int( line.substr( 56, 14 ) );

    // retrieve pointer format, row index format, coefficient format,
    // right-hand side format
    std::string  ptrfmt, indfmt, valfmt, rhsfmt;
    bool         no_val = false;
    
    std::getline( in, line );
    ptrfmt = line.substr(  0, 16 ); to_upper( ptrfmt );
    indfmt = line.substr( 16, 16 ); to_upper( indfmt );

    if (( line.size() > 32 ) && ( line.size() - 32 >= 20 ))
    {
        valfmt = line.substr( 32, 20 ); to_upper( valfmt );
    }// if
    else
        no_val = true;

    if (( line.size() > 52 ) && ( line.size() - 52 >= 20 ))
    {
        rhsfmt = line.substr( 52, 20 ); to_upper( rhsfmt );
    }// if

    // optional data of right-hand sides
    // UNUSED : std::string  rhstype;
    // UNUSED : uint    nrhs = 0, nrhsix = 0;

    if ( rhscrd != 0 )
    {
        // UNUSED : std::getline( in, line );
        // UNUSED : rhstype = line.substr( 0, 3 );
        // UNUSED : nrhs    = str_to_int( line.substr(  14, 14 ) );
        // UNUSED : nrhsix  = str_to_int( line.substr(  28, 14 ) );
    }// if

    //
    // eval fields above
    //

    bool  is_complex = false;
    bool  is_sym     = false;
    bool  is_herm    = false;
    bool  is_pattern = false;
    bool  is_integer = false;

    mxtype[0] = char( toupper( mxtype[0] ) );
    mxtype[1] = char( toupper( mxtype[1] ) );
    mxtype[2] = char( toupper( mxtype[2] ) );
    
    if      ( mxtype[0] == 'R' ) is_complex = false;
    else if ( mxtype[0] == 'C' ) is_complex = true;
    else if ( mxtype[0] == 'P' ) is_pattern = true;
    else if ( mxtype[0] == 'I' ) is_integer = true;

    if      ( mxtype[1] == 'U' ) is_sym  = false;
    else if ( mxtype[1] == 'S' ) is_sym  = true;
    else if ( mxtype[1] == 'H' ) is_herm = true;

    if ( is_complex && ! is_complex_type< value_t >::value )
        HERROR( ERR_REAL_CMPLX, "(THBMatrixIO) read", "found complex valued data but real valued matrix request" );

    if ( mxtype[2] == 'E' )
        HERROR( ERR_NOT_IMPL, "(THBMatrixIO) read", "elemental matrices not supported" );

    if ( no_val && ! is_pattern )
        HERROR( ERR_NOT_IMPL, "(THBMatrixIO) read", "no value format but not pattern" );
    
    //
    // read column pointers (indices to row index array)
    //

    auto  colptr        = std::vector< uint >( ncols+1, 0 );
    uint  pos           = 0;
    uint  elem_len      = 0;
    uint  elem_per_line = 0;

    parse_int_fmt( ptrfmt, elem_len, elem_per_line );
    
    for ( uint i = 0; i < ptrcrd; i++ )
    {
        std::getline( in, line );
        
        for ( uint j = 0; (j < elem_per_line) && (pos < ncols+1); j++, pos++ )
        {
            const uint  idx = str_to_int( & line[j*elem_len], & line[(j+1)*elem_len] );

            if ( idx-1 > nnzero )
                HERROR( ERR_FMT_HB, "(THBMatio) read", "column pointer > no. of nonzeros" );
                
            colptr[pos] = idx-1;
        }// for
    }// for
    
    //
    // read row indices
    //

    std::vector< uint >  rowind( nnzero, 0 );

    parse_int_fmt( indfmt, elem_len, elem_per_line );

    pos = 0;
    
    for ( uint i = 0; i < indcrd; i++ )
    {
        std::getline( in, line );
        
        for ( uint j = 0; (j < elem_per_line) && (pos < nnzero); j++, pos++ )
        {
            const uint  idx = str_to_int( & line[j*elem_len], & line[(j+1)*elem_len] );

            if ( idx > nrows )
                HERROR( ERR_FMT_HB, "(THBMatio) read", "row index > number of rows" );
            
            rowind[pos] = idx-1;
        }// for
    }// for

    // test for coefficients in upper half
    if ( is_sym || is_herm )
    {
        for ( uint col = 0; col < ncols; col++ )
        {
            const uint lb = colptr[col];
            const uint ub = colptr[col+1];
            
            for ( uint j = lb; j < ub; j++ )
            {
                const uint row = rowind[j];

                if ( row < col )
                    HERROR( ERR_FMT_HB, "(THBMatio) read", "element in upper right part of symmetric matrix" );
            }// for
        }// for
    }// if
            
    //
    // read coefficients
    //

    std::vector< value_t > coeffs( nnzero, value_t( 0 ) );

    if ( is_pattern )
    {
        // just fill with ones
        for ( uint i = 0; i < nnzero; i++ )
            coeffs[i] = value_t(1);
    }// if
    else if ( is_integer )
    {
        parse_int_fmt( valfmt, elem_len, elem_per_line );

        pos = 0;
        
        for ( uint i = 0; i < valcrd; i++ )
        {
            std::getline( in, line );
            
            for ( uint j = 0; (j < elem_per_line) && (pos < nnzero); j++, pos++ )
            {
                const double val = str_to_dbl( & line[j*elem_len], & line[(j+1)*elem_len] );
                
                coeffs[pos] = value_t(val);
            }// for
        }// for
    }// if
    else
    {
        parse_flt_fmt( valfmt, elem_len, elem_per_line );

        pos = 0;

        if ( is_complex )
        {
            auto  rval = real_t(0);
            auto  ival = real_t(0);
            bool  real_part = true;
            
            for ( uint i = 0; i < valcrd; i++ )
            {
                std::getline( in, line );
                
                for ( uint j = 0; (j < elem_per_line) && (pos < nnzero); j++ )
                {
                    const double val = str_to_dbl( & line[j*elem_len], & line[(j+1)*elem_len] );

                    if ( real_part )
                        rval = real_t(val);
                    else
                    {
                        ival        = real_t(val);
                        coeffs[pos] = get_value< value_t >::compose( rval, ival );

                        // increase pos only if complete coeff was read
                        pos++;
                    }// else

                    real_part = ! real_part;
                }// for
            }// for
        }// if
        else
        {
            for ( uint i = 0; i < valcrd; i++ )
            {
                std::getline( in, line );
                
                for ( uint j = 0; (j < elem_per_line) && (pos < nnzero); j++, pos++ )
                {
                    const double  val = str_to_dbl( & line[j*elem_len], & line[(j+1)*elem_len] );
                
                    coeffs[pos] = value_t(val);
                }// for
            }// for
        }// else
    }// if

    in_ptr.reset( nullptr );
    
    //
    // convert compressed column storage to compressed row storage
    // and build matrix
    // - if symmetric or hermitian, copy the entries to the upper
    //   part
    //

    auto    counter     = std::vector< uint >( nrows, 0 );
    auto    S           = std::make_unique< TSparseMatrix< value_t > >( nrows, ncols );
    size_t  full_nnzero = 0;

    // count elements per row
    for ( uint col = 0; col < ncols; col++ )
    {
        const uint lb = colptr[col];
        const uint ub = colptr[col+1];

        for ( uint j = lb; j < ub; j++ )
        {
            const uint row = rowind[j];
            
            counter[ row ]++;
            full_nnzero++;
            
            // add another entry for upper, right part
            if (( row != col ) && ( is_sym || is_herm ))
            {
                counter[ col ]++;
                full_nnzero++;
            }// if
        }// for
    }// for

    S->init( full_nnzero );

    // create row pointers (indices to column index array)
    pos = 0;
    
    for ( uint row = 0; row < nrows; row++ )
    {
        S->rowptr(row) = pos;
        pos           += counter[row];
    }// for

    S->rowptr(nrows) = idx_t(full_nnzero);

    // fill column indices and coefficients
    for ( uint row = 0; row < nrows; row++ ) counter[row] = 0;
    
    for ( idx_t  col = 0; col < idx_t(ncols); col++ )
    {
        const idx_t  lb = colptr[col];
        const idx_t  ub = colptr[col+1];

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  row = rowind[j];
            
            // standard entry
            {
                const idx_t  idx = S->rowptr(row) + counter[row];
            
                S->colind( idx ) = col;
                S->coeff(  idx ) = coeffs[j];

                counter[row]++;
            }
            
            // entry in upper, right part in symmetric case
            if (( row != col ) && ( is_sym || is_herm ))
            {
                const idx_t  idx = S->rowptr(col) + counter[col];
                
                S->colind( idx ) = row;
                S->coeff(  idx ) = coeffs[j];

                counter[col]++;
            }// if
        }// for
    }// for

    if ( is_sym  ) S->set_symmetric();
    if ( is_herm ) S->set_hermitian();

    return std::unique_ptr< TMatrix< value_t > >( S.release() );
}

//
// write vector to file <fname>
//
template < typename value_t >
void
THBVectorIO::write ( const TVector< value_t > *  x,
                     const std::string &         filename ) const
{
    const size_t  nrows  = x->size();
    const size_t  ncols  = 1;
    const size_t  nnz    = nrows;
        
    auto            out_ptr = open_write( filename );
    std::ostream &  out     = * out_ptr.get();
        
    // title and key
    out << to_string( "%-72s%-8s", "Created by HLIBpro", "" ) << std::endl;

    // line with number of lines for each type
    const uint  ptrcrd = uint( (ncols+1) / 7 + ( (ncols+1) % 7 == 0 ? 0 : 1 ) ); // 7 colptr per line plus rest
    const uint  indcrd = uint( nrows / 7 + ( nrows % 7 == 0 ? 0 : 1 ) );         // 7 rowind per line plus rest
    const uint  nvals  = uint( x->is_complex() ? 2*nnz : nnz );
    const uint  valcrd = nvals / 3 + ( nvals % 3 == 0 ? 0 : 1 ); // 3 values per line plus rest
    const uint  rhscrd = 0;                                      // no right hand sides
    const uint  totcrd = ptrcrd + indcrd + valcrd + rhscrd;

    out << to_string( "%-14d%-14d%-14d%-14d", totcrd, ptrcrd, indcrd, valcrd, rhscrd ) << std::endl;

    // matrix type, dimension of matrix, etc.
    std::string  mxtype( 4, '\0' );

    if ( is_complex_type< value_t >::value )
        mxtype[0] = 'C';
    else
        mxtype[0] = 'R';

    mxtype[1] = 'U';
    mxtype[2] = 'A';

    out << mxtype << boost::format( "            %-14d%-14d%-14d%-14d" ) % nrows % ncols % nnz % 0 << std::endl;

    // pointer, index and value format
    out << "(7i10)          (7i10)          (3e24.16)           (0e0.0)" << std::endl;

    // write column pointers
    out << boost::format( " %-9d" ) % 1 
        << boost::format( " %-9d" ) % nnz << std::endl;

    // write row indices
    uint  pos = 0;

    for ( uint i = 0; i < nrows; i++ )
    {
        hb_write_int( out, i+1 );
        if ( (++pos) % 7 == 0 )
            out << std::endl;
    }// for
    
    if ( pos % 7 != 0 )
        out << std::endl;

    // write coeffs
    if ( is_complex_type< value_t >::value )
    {
        pos = 0;
        
        for ( uint i = 0; i < nrows; i++ )
        {
            const auto  val = x->entry( i );
            
            hb_write_flt( out, std::real( val ) );
            if ( (++pos) % 3 == 0 ) out << std::endl;
            hb_write_flt( out, std::imag( val ) );
            if ( (++pos) % 3 == 0 ) out << std::endl;
        }// for
    }// if
    else
    {
        pos = 0;
        
        for ( uint i = 0; i < nrows; i++ )
        {
            hb_write_flt( out, std::real( x->entry( i ) ) );
            if ( (++pos) % 3 == 0 ) out << std::endl;
        }// for
    }// if
}

//
// read vector from file <fname>
//
template < typename value_t >
std::unique_ptr< TVector< value_t > >
THBVectorIO::read ( const std::string & fname ) const
{
    using  real_t = real_type_t< value_t >;
    
    auto            in_ptr = open_read( fname );
    std::istream &  in     = * in_ptr.get();
   
    // get title and key
    std::string  line;
    std::string  title( 73, ' ' ), key( 9, ' ' );
    
    std::getline( in, line );
    title = line.substr( 0, 72 );
    key   = line.substr( 72, 8 );
    
    // retrieve total number of lines, number of lines for column pointers,
    // number of lines for row indices, number of lines for coefficients,
    // number of lines for right-hand sides
    // UNUSED : uint  totcrd;
    uint  ptrcrd, indcrd, valcrd, rhscrd;

    std::getline( in, line );
    // UNUSED : totcrd = def_lexical_cast< uint >( line.substr(  0, 14 ) );
    ptrcrd = str_to_int( line.substr( 14, 14 ) );
    indcrd = str_to_int( line.substr( 28, 14 ) );
    valcrd = str_to_int( line.substr( 42, 14 ) );
    rhscrd = str_to_int( line.substr( 56, 14 ) );
    
    // retrieve matrix type, number of rows, number of columns,
    // number of non-zeroes, number of elemental matrix entries,
    std::string  mxtype;
    uint         nrows, ncols, nnzero;
    // UNUSED : uint    neltvl;
    
    std::getline( in, line );
    mxtype = line.substr( 0, 3 );
    nrows  = str_to_int( line.substr( 14, 14 ) );
    ncols  = str_to_int( line.substr( 28, 14 ) );
    nnzero = str_to_int( line.substr( 42, 14 ) );
    // UNUSED : neltvl = str_to_int( line.substr( 56, 14 ) );

    if ( ncols != 1 )
        HERROR( ERR_FMT_HB, "(THBVectorIO) read_vector", "only 1 column per vector supported" );
    
    // retrieve pointer format, row index format, coefficient format,
    // right-hand side format
    std::string  ptrfmt, indfmt, valfmt, rhsfmt;
    bool         no_val = false;
    
    std::getline( in, line );
    ptrfmt = line.substr(  0, 16 ); to_upper( ptrfmt );
    indfmt = line.substr( 16, 16 ); to_upper( indfmt );

    if (( line.size() > 32 ) && ( line.size() - 32 >= 20 ))
    {
        valfmt = line.substr( 32, 20 ); to_upper( valfmt );
    }// if
    else
        no_val = true;
    
    if (( line.size() > 52 ) && ( line.size() - 52 >= 20 ))
    {
        rhsfmt = line.substr( 52, 20 ); to_upper( rhsfmt );
    }// if
    
    // optional data of right-hand sides
    // UNUSED : std::string  rhstype;
    // UNUSED : uint    nrhs = 0, nrhsix = 0;

    if ( rhscrd != 0 )
    {
        // UNUSED : std::getline( in, line );
        // UNUSED : rhstype = line.substr( 0, 3 );
        // UNUSED : nrhs    = str_to_int( line.substr(  14, 14 ) );
        // UNUSED : nrhsix  = str_to_int( line.substr(  28, 14 ) );
    }// if

    //
    // eval fields above
    //

    bool is_complex = false;
    bool is_pattern = false;
    bool is_integer = false;

    if      ( mxtype[0] == 'R' ) is_complex = false;
    else if ( mxtype[0] == 'C' ) is_complex = true;
    else if ( mxtype[0] == 'P' ) is_pattern = true;
    else if ( mxtype[0] == 'I' ) is_integer = true;

    if ( is_complex && ! is_complex_type< value_t >::value )
        HERROR( ERR_REAL_CMPLX, "(THBMatrixIO) read", "found complex valued data but real valued vector request" );

    if ( mxtype[2] == 'E' )
        HERROR( ERR_NOT_IMPL, "(THBVectorIO) read", "elemental matrices not supported" );

    if ( no_val && ! is_pattern )
        HERROR( ERR_NOT_IMPL, "(THBMatrixIO) read", "no value format but not pattern" );
    
    //
    // read column pointers (indices to row index array)
    //

    auto  colptr        = std::vector< uint >( ncols+1, 0 );
    uint  pos           = 0;
    uint  elem_len      = 0;
    uint  elem_per_line = 0;

    parse_int_fmt( ptrfmt, elem_len, elem_per_line );
    
    for ( uint i = 0; i < ptrcrd; i++ )
    {
        std::getline( in, line );
        
        for ( uint j = 0; (j < elem_per_line) && (pos < ncols+1); j++, pos++ )
        {
            const uint  idx = str_to_int( & line[j*elem_len], & line[(j+1)*elem_len] );

            colptr[pos] = idx-1;
        }// for
    }// for
    
    //
    // read row indices
    //

    std::vector< uint >  rowind( nnzero, 0 );

    parse_int_fmt( indfmt, elem_len, elem_per_line );

    pos = 0;
    
    for ( uint i = 0; i < indcrd; i++ )
    {
        std::getline( in, line );
        
        for ( uint j = 0; (j < elem_per_line) && (pos < nnzero); j++, pos++ )
        {
            const uint  idx = str_to_int( & line[j*elem_len], & line[(j+1)*elem_len] );

            rowind[pos] = idx-1;
        }// for
    }// for
    
    //
    // read coefficients and stored them directly in vector
    //

    auto  v = std::make_unique< TScalarVector< value_t > >( nrows, 0 );

    if ( is_pattern )
    {
        // just fill with ones
        for ( uint i = 0; i < nnzero; i++ )
            v->set_entry( rowind[i], 1.0 );
    }// if
    else if ( is_integer )
    {
        parse_int_fmt( valfmt, elem_len, elem_per_line );

        pos = 0;
        
        for ( uint i = 0; i < valcrd; i++ )
        {
            std::getline( in, line );
            
            for ( uint j = 0; (j < elem_per_line) && (pos < nnzero); j++, pos++ )
            {
                const double val = str_to_dbl( & line[j*elem_len], & line[(j+1)*elem_len] );
                
                v->set_entry( rowind[pos], value_t( val ) );
            }// for
        }// for
    }// if
    else
    {
        parse_flt_fmt( valfmt, elem_len, elem_per_line );

        pos = 0;

        if ( is_complex_type< value_t >::value )
        {
            auto  rval = real_t(0);
            auto  ival = real_t(0);
            bool  real_part = true;
            
            for ( uint i = 0; i < valcrd; i++ )
            {
                std::getline( in, line );
            
                for ( uint j = 0; (j < elem_per_line) && (pos < nnzero); j++ )
                {
                    const double  val = str_to_dbl( & line[j*elem_len], & line[(j+1)*elem_len] );

                    if ( real_part )
                        rval = real_t( val );
                    else
                    {
                        ival = real_t( val );
                        v->set_entry( rowind[pos], get_value< value_t >::compose( rval, ival ) );

                        // increase pos only if complete coeff was read
                        pos++; 
                    }// else

                    real_part = ! real_part;
                }// for
            }// for
        }// if
        else
        {
            for ( uint i = 0; i < valcrd; i++ )
            {
                std::getline( in, line );
            
                for ( uint j = 0; (j < elem_per_line) && (pos < nnzero); j++, pos++ )
                {
                    const auto  val = str_to_dbl( & line[j*elem_len], & line[(j+1)*elem_len] );
                
                    v->set_entry( rowind[pos], real_t( val ) );
                }// for
            }// for
        }// else
    }// if

    return std::unique_ptr< TVector< value_t > >( v.release() );
}

variant_id_t
hb_guess_value_type ( const std::string &  filename )
{
    auto            in_ptr = open_read( filename );
    std::istream &  in     = * in_ptr.get();
   
    // get title and key
    std::string  line;
    std::string  mxtype;
    
    std::getline( in, line );
    std::getline( in, line );
    std::getline( in, line );
    mxtype = line.substr( 0, 3 );
    
    if ( mxtype[0] == 'C' )
        return COMPLEX_FP64;
    else
        return REAL_FP64;
}

//
// explicit template instantiation
//
template void THBMatrixIO::write< float >                  ( const TMatrix< float > *,                  const std::string & ) const;
template void THBMatrixIO::write< double >                 ( const TMatrix< double > *,                 const std::string & ) const;
template void THBMatrixIO::write< std::complex< float > >  ( const TMatrix< std::complex< float > > *,  const std::string & ) const;
template void THBMatrixIO::write< std::complex< double > > ( const TMatrix< std::complex< double > > *, const std::string & ) const;

template std::unique_ptr< TMatrix< float > >                  THBMatrixIO::read< float >                  ( const std::string & fname ) const;
template std::unique_ptr< TMatrix< double > >                 THBMatrixIO::read< double >                 ( const std::string & fname ) const;
template std::unique_ptr< TMatrix< std::complex< float > > >  THBMatrixIO::read< std::complex< float > >  ( const std::string & fname ) const;
template std::unique_ptr< TMatrix< std::complex< double > > > THBMatrixIO::read< std::complex< double > > ( const std::string & fname ) const;

template void THBVectorIO::write< float >                  ( const TVector< float > *,                  const std::string & ) const;
template void THBVectorIO::write< double >                 ( const TVector< double > *,                 const std::string & ) const;
template void THBVectorIO::write< std::complex< float > >  ( const TVector< std::complex< float > > *,  const std::string & ) const;
template void THBVectorIO::write< std::complex< double > > ( const TVector< std::complex< double > > *, const std::string & ) const;

template std::unique_ptr< TVector< float > >                  THBVectorIO::read< float >                  ( const std::string & fname ) const;
template std::unique_ptr< TVector< double > >                 THBVectorIO::read< double >                 ( const std::string & fname ) const;
template std::unique_ptr< TVector< std::complex< float > > >  THBVectorIO::read< std::complex< float > >  ( const std::string & fname ) const;
template std::unique_ptr< TVector< std::complex< double > > > THBVectorIO::read< std::complex< double > > ( const std::string & fname ) const;

}// namespace Hpro
