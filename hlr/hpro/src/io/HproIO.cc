//
// Project     : HLIBpro
// File        : HLibIO.cc
// Description : HLIB-format IO classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <cstring>
#include <vector>
#include <memory>
#include <mutex>

#include <boost/filesystem.hpp>

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TBlockMatrix.hh"
#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"
#include "hpro/matrix/TSparseMatrix.hh"
#include "hpro/matrix/THMatrix.hh"
#include "hpro/matrix/TZeroMatrix.hh"

#include "hpro/vector/TScalarVector.hh"
#include "hpro/vector/TBlockVector.hh"

#include "baseio.hh"

#include "hpro/io/TVectorIO.hh"
#include "hpro/io/TMatrixIO.hh"
#include "hpro/io/TCoordIO.hh"

namespace Hpro
{

namespace B = BLAS;

using std::unique_ptr;
using std::make_unique;

namespace
{

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// local defines
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// for automatic swapping of bytes
#define BYTE_SWAP( t )  { if ( byteswap ) swap_bytes( t ); }

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// local types
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// identifiers for types
//
enum { HMFILE_NULL        = 0,

       // matrices
       HMFILE_MAT_BLOCK     = 1,
       HMFILE_MAT_DENSE     = 2,
       HMFILE_MAT_RANK      = 3,
       HMFILE_MAT_SPARSE    = 4,
       HMFILE_MAT_H         = 5,
       HMFILE_MAT_PERM      = 6,
       HMFILE_MAT_FAC       = 7,
       HMFILE_MAT_FACINV    = 8,
       HMFILE_MAT_H2        = 9,
       HMFILE_MAT_UNIFORM   = 10,
       HMFILE_MAT_ZERO      = 11,

       // factorisation subtypes
       HMFILE_FAC_LU        = 50,
       HMFILE_FAC_LDL       = 51,
       HMFILE_FAC_LL        = 52,
       HMFILE_FAC_LDU       = 53,
       HMFILE_FAC_NONE      = 60,  // dummy value 
       
       // vectors
       HMFILE_VEC_BLOCK     = 100,
       HMFILE_VEC_SCALAR    = 101,

       // coordinates
       HMFILE_COORD_VTX     = 200,  // only vertex data (single coordinate)
       HMFILE_COORD_BBOX    = 201,  // only bounding box (two coordinates: min + max)
       HMFILE_COORD_VTXBBOX = 202   // vertex + bounding box (three coord: vtx + min + max)
};

const uint16_t ENDIANESS_INDICATOR  = 0x484d;              // "HM" in ASCII
const uint16_t SAME_ENDIANESS       = ENDIANESS_INDICATOR; // same endianess as on file host
const uint16_t REV_ENDIANESS        = 0x4d48;              // reversed endianess as on file host

const uint16_t HLIBPRO_FILE_VERSION = 0x03;                // file format version

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// local functions
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// return globally unique ID (for Version1 matrices)
//
int
get_id ()
{
    static std::mutex  id_mutex;
    static int         id_counter = 0;
    int                ret_val    = 0;

    {
        std::scoped_lock  lock( id_mutex );

        ret_val = id_counter++;
    }

    return ret_val;
}

//
// template functions for reading and writing data
//
template < typename value_t >
value_t
hlib_read ( std::istream &  in,
            const bool      byteswap )
{
    value_t  val;
    
    in.read( reinterpret_cast< char * >( & val ), sizeof(value_t) );
    BYTE_SWAP( val );

    return val;
}

template < typename value_t >
void
hlib_read ( std::istream &  in,
            value_t *       data,
            const uint      nelem,
            const bool      byteswap )
{
    if ( nelem > 0 )
    {
        in.read( reinterpret_cast< char * >( data ), sizeof(value_t)*nelem );
    
        if ( byteswap && ( sizeof(value_t) > 1 ))
        {
            for ( uint i = 0; i < nelem; i++ )
                BYTE_SWAP( data[i] );
        }// if
    }// if
}

// template <>
// void
// hlib_read ( std::istream &  in, complex * data, const uint nelem, const bool byteswap )
// {
//     if ( nelem > 0 )
//     {
//         in.read( reinterpret_cast< char * >( data ), sizeof(complex)*nelem );

//         if ( byteswap )
//         {
//             for ( uint i = 0; i < nelem; i++ )
//             {
//                 BYTE_SWAP( std::real(data[i]) );
//                 BYTE_SWAP( std::imag(data[i]) );
//             }// for
//         }// if
//     }// if
// }

template < typename value_t >
void
hlib_read ( std::istream &          in,
            const bool              byteswap,
            B::Matrix< value_t > &  mat )
{
    for ( idx_t  j = 0; j < idx_t( mat.ncols() ); ++j )
    {
        for ( idx_t  i = 0; i < idx_t( mat.nrows() ); ++i )
        {
            mat( i, j ) = hlib_read< value_t >( in, byteswap );
        }// for
    }// for
}

template < typename value_t >
void
hlib_write ( std::ostream &  out,
             const value_t   data )
{
    out.write( reinterpret_cast< const char * >( & data ), sizeof(value_t) );
}

template <>
void
hlib_write ( std::ostream &  out, const float  data )
{
    double  t = double( data );
    
    out.write( reinterpret_cast< const char * >( & t ), sizeof(t) );
}

template <>
void
hlib_write ( std::ostream &  out, const std::complex< float >  data )
{
    double  t[2] = { double( std::real( data ) ), double( std::imag( data ) ) };
    
    out.write( reinterpret_cast< const char * >( & t ), sizeof(t) );
}

template <>
void
hlib_write ( std::ostream &  out, const std::complex< double >  data )
{
    double  t[2] = { double( std::real( data ) ), double( std::imag( data ) ) };
    
    out.write( reinterpret_cast< const char * >( & t ), sizeof(t) );
}

template < typename value_t >
void
hlib_write_mat ( std::ostream &                out,
                 const B::Matrix< value_t > &  mat )
{
    for ( idx_t  j = 0; j < idx_t( mat.ncols() ); ++j )
    {
        for ( idx_t  i = 0; i < idx_t( mat.nrows() ); ++i )
        {
            hlib_write( out, mat( i, j ) );
        }// for
    }// for
}

//
// write header in HLIBpro format
//
void
write_header ( std::ostream &  out )
{
    char      text[126] = "HMATRIX v2.0, Created by HLibPro";
    uint16_t  endian   = ENDIANESS_INDICATOR;
    uint16_t  version  = HLIBPRO_FILE_VERSION;

    for ( size_t i = strlen(text); i < 126; i++ ) text[i] = ' ';
    
    out.write( text, 126 );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
}

//
// read header in HLIBpro format
//
void
read_header ( std::istream &  in,
              bool &          byteswap,
              uint16_t &      version )
{
    char      text[126];
    uint16_t  endian  = ENDIANESS_INDICATOR;
    
    in.read( text, sizeof(text) );
    in.read( reinterpret_cast< char * >( & endian ),  sizeof(endian)  );
    in.read( reinterpret_cast< char * >( & version ), sizeof(version) );

    if      ( endian == SAME_ENDIANESS ) byteswap = false;
    else if ( endian == REV_ENDIANESS  ) byteswap = true;
    else
        HERROR( ERR_FMT_HFORMAT, "read_header", "endianess indicator invalid" );
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// HLIBpro write methods
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// forward decl.
//
template < typename value_t >
void
write_mat  ( std::ostream &             out,
             const TMatrix< value_t > * A );

//
// write header for matrices
//
template < typename value_t >
void
write_mat_header ( std::ostream &             out,
                   const TMatrix< value_t > * A,
                   uint                       type )
{
    uint16_t  ltype      = uint16_t(type);
    int32_t   id         = int32_t(A->id());
    uint32_t  row_ofs    = uint32_t(A->row_ofs());
    uint32_t  col_ofs    = uint32_t(A->col_ofs());
    uint32_t  proc_first = uint32_t(A->procs().first());
    uint32_t  proc_last  = uint32_t(A->procs().last());
    uint8_t   form       = uint8_t(A->form());
    uint8_t   is_complex = uint8_t(is_complex_type< value_t >::value);
    uint32_t  rows       = uint32_t(A->rows());
    uint32_t  cols       = uint32_t(A->cols());

    out.write( reinterpret_cast< const char * >( & ltype ),      sizeof(ltype) );
    out.write( reinterpret_cast< const char * >( & id ),         sizeof(id) );
    out.write( reinterpret_cast< const char * >( & row_ofs ),    sizeof(row_ofs) );
    out.write( reinterpret_cast< const char * >( & col_ofs ),    sizeof(col_ofs) );
    out.write( reinterpret_cast< const char * >( & proc_first ), sizeof(proc_first) );
    out.write( reinterpret_cast< const char * >( & proc_last ),  sizeof(proc_last) );
    out.write( reinterpret_cast< const char * >( & form ),       sizeof(form) );
    out.write( reinterpret_cast< const char * >( & is_complex ), sizeof(is_complex) );
    out.write( reinterpret_cast< const char * >( & rows ),       sizeof(rows) );
    out.write( reinterpret_cast< const char * >( & cols ),       sizeof(cols) );
}

//
// special write methods
//
template < typename value_t >
void
write_mat_block  ( std::ostream &                   out,
                   const TBlockMatrix< value_t > *  B )
{
    const bool      sym        = B->is_symmetric() || B->is_hermitian();
    const uint32_t  block_rows = B->block_rows();
    const uint32_t  block_cols = B->block_cols();

    
    write_mat_header( out, B, HMFILE_MAT_BLOCK );
    
    out.write( reinterpret_cast< const char * >( & block_rows ), sizeof(uint32_t) );
    out.write( reinterpret_cast< const char * >( & block_cols ), sizeof(uint32_t) );

    if ( sym )
    {
        for ( uint i = 0; i < block_rows; i++ )
            for ( uint j = 0; j <= i; j++ )
            {
                write_mat( out, B->block(i,j) );
            }// for
    }// if
    else
    {
        for ( uint i = 0; i < block_rows; i++ )
            for ( uint j = 0; j < block_cols; j++ )
            {
                write_mat( out, B->block(i,j) );
            }// for
    }// else
}

template < typename value_t >
void
write_mat_dense  ( std::ostream &                   out,
                   const TDenseMatrix< value_t > *  D )
{
    using  real_t = real_type_t< value_t >;
    
    const size_t  rows       = D->rows();
    const size_t  cols       = D->cols();

    write_mat_header( out, D, HMFILE_MAT_DENSE );

    if ( sizeof(real_t) == sizeof(double) )
        out.write( reinterpret_cast< const char * >( D->blas_mat().data() ), sizeof(value_t) * rows * cols );
    else
    {
        if ( is_complex_type< value_t >::value )
        {
            for ( idx_t  j = 0; j < idx_t(cols); ++j )
                for ( idx_t  i = 0; i < idx_t(rows); ++i )
                {
                    const double z[2] = { std::real( D->blas_mat()(i,j) ), std::imag( D->blas_mat()(i,j) ) }; 
                    
                    out.write( reinterpret_cast< const char * >( z ), 2 * sizeof(double) );
                }// for
        }// if
        else
        {
            for ( idx_t  j = 0; j < idx_t(cols); ++j )
                for ( idx_t  i = 0; i < idx_t(rows); ++i )
                {
                    const double  f = std::real( D->blas_mat()(i,j) );
                
                    out.write( reinterpret_cast< const char * >( & f ), sizeof(double) );
                }// for
        }// else
    }// else
}

template < typename value_t >
void
write_mat_rank ( std::ostream &                out,
                 const TRkMatrix< value_t > *  R )
{
    using  real_t = real_type_t< value_t >;
    
    const uint32_t  rank = uint32_t(R->rank());

    write_mat_header( out, R, HMFILE_MAT_RANK );

    out.write( reinterpret_cast< const char * >( & rank ), sizeof(uint32_t) );

    if ( sizeof(real_t) == sizeof(double) )
    {
        out.write( reinterpret_cast< const char * >( R->blas_mat_A().data() ), sizeof(value_t) * R->rows() * rank );
        out.write( reinterpret_cast< const char * >( R->blas_mat_B().data() ), sizeof(value_t) * R->cols() * rank );
    }// if
    else
    {
        if ( is_complex_type< value_t >::value )
        {
            for ( idx_t  k = 0; k < idx_t(R->rank()); ++k )
                for ( idx_t  i = 0; i < idx_t(R->rows()); ++i )
                {
                    const double z[2] = { std::real( R->blas_mat_A()(i,k) ), std::imag( R->blas_mat_A()(i,k) ) };
                    
                    out.write( reinterpret_cast< const char * >( z ), 2 * sizeof(double) );
                }// for

            for ( idx_t  k = 0; k < idx_t(R->rank()); ++k )
                for ( idx_t  i = 0; i < idx_t(R->cols()); ++i )
                {
                    const double z[2] = { std::real( R->blas_mat_B()(i,k) ), std::imag( R->blas_mat_B()(i,k) ) };
                    
                    out.write( reinterpret_cast< const char * >( z ), 2 * sizeof(double) );
                }// for
        }// if
        else
        {
            for ( idx_t  k = 0; k < idx_t(R->rank()); ++k )
                for ( idx_t  i = 0; i < idx_t(R->rows()); ++i )
                {
                    const double  f = std::real( R->blas_mat_A()(i,k) );
                    
                    out.write( reinterpret_cast< const char * >( & f ), sizeof(double) );
                }// for

            for ( idx_t  k = 0; k < idx_t(R->rank()); ++k )
                for ( idx_t  i = 0; i < idx_t(R->cols()); ++i )
                {
                    const double  f = std::real( R->blas_mat_B()(i,k) );
                    
                    out.write( reinterpret_cast< const char * >( & f ), sizeof(double) );
                }// for
        }// else
    }// else
}

template < typename value_t >
void
write_mat_sparse ( std::ostream &                    out,
                   const TSparseMatrix< value_t > *  S )
{
    using  real_t = real_type_t< value_t >;
    
    const uint32_t  nnz = uint32_t(S->n_non_zero());

    write_mat_header( out, S, HMFILE_MAT_SPARSE );
               
    out.write( reinterpret_cast< const char * >( & nnz ), sizeof(nnz) );

    if ( sizeof(idx_t) == sizeof(int32_t) )
    {
        out.write( reinterpret_cast< const char * >( & const_cast< TSparseMatrix< value_t > * >( S )->rowptr(0) ),
                   sizeof(int32_t) * (S->rows() + 1) );
        out.write( reinterpret_cast< const char * >( & const_cast< TSparseMatrix< value_t > * >( S )->colind(0) ), sizeof(int32_t) * nnz );
    }// if
    else
    {
        const size_t   nrows = S->rows();
        
        for ( idx_t i = 0; i <= idx_t(nrows); i++ )
        {
            hlib_write< int32_t >( out, S->rowptr( i ) );
        }// for
        
        for ( idx_t  i = 0; i < idx_t(nnz); i++ )
        {
            hlib_write< int32_t >( out, S->colind( i ) );
        }// for
    }// else
    
    if ( sizeof(real_t) == sizeof(double) )
        out.write( reinterpret_cast< const char * >( & const_cast< TSparseMatrix< value_t > * >( S )->coeff(0) ),
                   sizeof(value_t) * nnz );
    else
    {
        for ( idx_t  i = 0; i < idx_t(nnz); i++ )
        {
            hlib_write< value_t >( out, S->coeff(i) );
        }// for
    }// else
}

void
write_mat_perm ( std::ostream &        out,
                 const TPermutation *  P )
{
    const uint32_t  n = uint32_t(P->size());
    
    out.write( reinterpret_cast< const char * >( & n ), sizeof(uint32_t) );

    for ( uint32_t  i = 0; i < n; i++ )
    {
        const int32_t  v = int32_t((*P)[i]);
        
        out.write( reinterpret_cast< const char * >( & v ), sizeof(int32_t) );
    }// for
}

template < typename value_t >
void
write_mat_h ( std::ostream &               out,
              const THMatrix< value_t > *  H )
{
    const bool      sym        = H->is_symmetric() || H->is_hermitian();
    const uint32_t  block_rows = H->block_rows();
    const uint32_t  block_cols = H->block_cols();
    const uint8_t   has_perm   = H->has_perm();
    
    write_mat_header( out, H, HMFILE_MAT_H );
    
    out.write( reinterpret_cast< const char * >( & block_rows ), sizeof(uint32_t) );
    out.write( reinterpret_cast< const char * >( & block_cols ), sizeof(uint32_t) );

    // write permutations if present
    out.write( reinterpret_cast< const char * >( & has_perm ), sizeof(uint8_t) );

    if ( H->has_perm() )
    {
        write_mat_perm( out, & H->row_perm_e2i() );
        write_mat_perm( out, & H->row_perm_i2e() );
        write_mat_perm( out, & H->col_perm_e2i() );
        write_mat_perm( out, & H->col_perm_i2e() );
    }// if

    // now write subblocks
    if ( sym )
    {
        for ( uint i = 0; i < block_rows; i++ )
            for ( uint j = 0; j <= i; j++ )
            {
                write_mat( out, H->block(i,j) );
            }// for
    }// if
    else
    {
        for ( uint i = 0; i < block_rows; i++ )
            for ( uint j = 0; j < block_cols; j++ )
            {
                write_mat( out, H->block(i,j) );
            }// for
    }// else
}

template < typename value_t >
void
write_mat_zero ( std::ostream &                  out,
                 const TZeroMatrix< value_t > *  M )
{
    write_mat_header( out, M, HMFILE_MAT_ZERO );
}

template < typename value_t >
void
write_mat ( std::ostream &              out,
            const TMatrix< value_t > *  A )
{
    if ( A == nullptr )
    {
        const uint16_t  type = HMFILE_NULL;

        out.write( reinterpret_cast< const char * >( & type ), sizeof(type) );
    }// if
    else if ( IS_TYPE( A, THMatrix       ) ) write_mat_h(       out, cptrcast( A, THMatrix< value_t >       ) );
    else if ( IS_TYPE( A, TBlockMatrix   ) ) write_mat_block(   out, cptrcast( A, TBlockMatrix< value_t >   ) );
    else if ( IS_TYPE( A, TDenseMatrix   ) ) write_mat_dense(   out, cptrcast( A, TDenseMatrix< value_t >   ) );
    else if ( IS_TYPE( A, TRkMatrix      ) ) write_mat_rank(    out, cptrcast( A, TRkMatrix< value_t >      ) );
    else if ( IS_TYPE( A, TSparseMatrix  ) ) write_mat_sparse(  out, cptrcast( A, TSparseMatrix< value_t >  ) );
    else if ( IS_TYPE( A, TZeroMatrix    ) ) write_mat_zero(    out, cptrcast( A, TZeroMatrix< value_t >    ) );
    else
        HERROR( ERR_MAT_TYPE, "(THproMatrixIO) write", A->typestr() );
}

namespace LinOp
{

template < typename value_t >
void
write ( std::ostream &                      out,
        const TLinearOperator< value_t > *  A );

template < typename value_t >
void
write ( std::ostream &                      out,
        const TLinearOperator< value_t > *  A )
{
    if ( A == nullptr )
    {
        const uint16_t  type = HMFILE_NULL;

        out.write( reinterpret_cast< const char * >( & type ), sizeof(type) );
    }// if
    else if ( dynamic_cast< const TMatrix< value_t > * >( A ) != nullptr )
    {
        const TMatrix< value_t > *  M = cptrcast( A, TMatrix< value_t > );
        
        if      ( IS_TYPE( M, THMatrix       ) ) write_mat_h(       out, cptrcast( M, THMatrix< value_t >       ) );
        else if ( IS_TYPE( M, TBlockMatrix   ) ) write_mat_block(   out, cptrcast( M, TBlockMatrix< value_t >   ) );
        else if ( IS_TYPE( M, TDenseMatrix   ) ) write_mat_dense(   out, cptrcast( M, TDenseMatrix< value_t >   ) );
        else if ( IS_TYPE( M, TRkMatrix      ) ) write_mat_rank(    out, cptrcast( M, TRkMatrix< value_t >      ) );
        else if ( IS_TYPE( M, TSparseMatrix  ) ) write_mat_sparse(  out, cptrcast( M, TSparseMatrix< value_t >  ) );
        else if ( IS_TYPE( M, TZeroMatrix    ) ) write_mat_zero(    out, cptrcast( M, TZeroMatrix< value_t >    ) );
        else
            HERROR( ERR_MAT_TYPE, "(THproMatrixIO) write", M->typestr() );
    }// if
    else
        HERROR( ERR_MAT_TYPE, "(THproMatrixIO) write", "" );
}

}// namespace LinOp

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// HLIBpro read methods
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// forward decl.
//
template < typename value_t >
TMatrix< value_t > *
read_mat  ( std::istream &  in,
            bool            byteswap,
            const int       version );

//
// read matrix header
//
template < typename value_t >
void
read_mat_header ( std::istream &        in,
                  TMatrix< value_t > *  A,
                  bool                  byteswap,
                  const int             version )
{
    uint32_t   row_ofs    = 0;
    uint32_t   col_ofs    = 0;
    matform_t  form       = MATFORM_NONSYM;
    bool       is_complex = false;

    if ( version >= 2 )
    {
        int32_t    id         = -1;
        uint32_t   proc_first = 0;
        uint32_t   proc_last  = 0;
        
        id         = hlib_read<int32_t>(  in, byteswap );
        row_ofs    = hlib_read<uint32_t>( in, byteswap );
        col_ofs    = hlib_read<uint32_t>( in, byteswap );
        proc_first = hlib_read<uint32_t>( in, byteswap );
        proc_last  = hlib_read<uint32_t>( in, byteswap );
        form       = matform_t( hlib_read<uint8_t>( in, byteswap ) );
        is_complex = hlib_read<uint8_t>(  in, byteswap ) != 0;

        A->set_id( id );
        A->set_ofs( row_ofs, col_ofs );
        A->set_procs( ps( proc_first, proc_last ) );
        A->set_form( form );
    }// if
    else
    {
        row_ofs    = hlib_read<uint32_t>( in, byteswap );
        col_ofs    = hlib_read<uint32_t>( in, byteswap );
        form       = matform_t( hlib_read<uint8_t>( in, byteswap ) );
        is_complex = hlib_read<uint8_t>(  in, byteswap ) != 0;

        A->set_id( get_id() );
        A->set_ofs( row_ofs, col_ofs );
        A->set_procs( ps_single( NET::pid() ) );
        A->set_form( form );
    }// else

    if ( is_complex && ! is_complex_type< value_t >::value )
        HERROR( ERR_REAL_CMPLX, "(HLibIO) read_mat_header", "complex valued data but real value requested" );
}

//
// special read methods
//
template < typename value_t >
TMatrix< value_t > *
read_mat_block  ( std::istream &  in,
                  bool            byteswap,
                  const int       version )
{
    auto      A = make_unique< TBlockMatrix< value_t > >();
    uint32_t  rows, cols;
    uint32_t  block_rows, block_cols;

    read_mat_header( in, A.get(), byteswap, version );
    rows       = hlib_read<uint32_t>( in, byteswap );
    cols       = hlib_read<uint32_t>( in, byteswap );
    block_rows = hlib_read<uint32_t>( in, byteswap );
    block_cols = hlib_read<uint32_t>( in, byteswap );

    A->set_size( rows, cols );
    A->set_block_struct( block_rows, block_cols );

    if ( A->is_nonsym() )
    {
        for ( uint i = 0; i < block_rows; i++ )
            for ( uint j = 0; j < block_cols; j++ )
                A->set_block( i, j, read_mat< value_t >( in, byteswap, version ) );
    }// if
    else
    {
        for ( uint i = 0; i < block_rows; i++ )
            for ( uint j = 0; j <= i; j++ )
                A->set_block( i, j, read_mat< value_t >( in, byteswap, version ) );
    }// else
    
    return A.release();
}

template < typename value_t >
TMatrix< value_t > *
read_mat_dense  ( std::istream &  in,
                  bool            byteswap,
                  const int       version )
{
    using  real_t = real_type_t< value_t >;

    auto      A = make_unique< TDenseMatrix< value_t > >();
    uint32_t  rows, cols;

    read_mat_header( in, A.get(), byteswap, version );
    rows = hlib_read<uint32_t>( in, byteswap );
    cols = hlib_read<uint32_t>( in, byteswap );

    A->set_size( rows, cols );

    if ( sizeof(real_t) == sizeof(double) )
        hlib_read<value_t>( in, A->blas_mat().data(), rows*cols, byteswap );
    else
    {
        if ( is_complex_type< value_t >::value )
        {
            for ( idx_t  j = 0; j < idx_t(cols); j++ )
                for ( idx_t  i = 0; i < idx_t(rows); i++ )
                {
                    double z[2];

                    hlib_read<double>( in, z, 2, byteswap );
                    
                    A->blas_mat()(i,j) = get_value< value_t >::compose( z[0], z[1] );
                }// for
        }// if
        else
        {
            for ( idx_t  j = 0; j < idx_t(cols); j++ )
                for ( idx_t  i = 0; i < idx_t(rows); i++ )
                    A->blas_mat()(i,j) = value_t( hlib_read<double>( in, byteswap ) );
        }// else
    }// else
    
    return A.release();
}

template < typename value_t >
TMatrix< value_t > *
read_mat_rank ( std::istream &  in,
                bool            byteswap,
                const int       version )
{
    using  real_t = real_type_t< value_t >;

    auto      A = make_unique< TRkMatrix< value_t > >();
    uint32_t  rows, cols, rank;

    read_mat_header( in, A.get(), byteswap, version );
    rows = hlib_read<uint32_t>( in, byteswap );
    cols = hlib_read<uint32_t>( in, byteswap );
    rank = hlib_read<uint32_t>( in, byteswap );

    A->set_size( rows, cols, rank );

    if ( sizeof(real_t) == sizeof(double) )
    {
        hlib_read< value_t >( in, A->blas_mat_A().data(), rank*rows, byteswap );
        hlib_read< value_t >( in, A->blas_mat_B().data(), rank*cols, byteswap );
    }// if
    else
    {
        if ( is_complex_type< value_t >::value )
        {
            for ( idx_t  k = 0; k < idx_t(rank); k++ )
                for ( idx_t  i = 0; i < idx_t(rows); i++ )
                {
                    double z[2];

                    hlib_read<double>( in, z, 2, byteswap );
                    A->blas_mat_A()(i,k) = get_value< value_t >::compose( z[0], z[1] );
                }// for

            for ( idx_t  k = 0; k < idx_t(rank); k++ )
                for ( idx_t  i = 0; i < idx_t(cols); i++ )
                {
                    double z[2];

                    hlib_read<double>( in, z, 2, byteswap );
                    A->blas_mat_B()(i,k) = get_value< value_t >::compose( z[0], z[1] );
                }// for
        }// if
        else
        {
            for ( idx_t  k = 0; k < idx_t(rank); k++ )
                for ( idx_t  i = 0; i < idx_t(rows); i++ )
                    A->blas_mat_A()(i,k) = value_t( hlib_read<double>( in, byteswap ) );

            for ( idx_t  k = 0; k < idx_t(rank); k++ )
                for ( idx_t  i = 0; i < idx_t(cols); i++ )
                    A->blas_mat_B()(i,k) = value_t( hlib_read<double>( in, byteswap ) );
        }// else
    }// if

    return A.release();
}

template < typename value_t >
TMatrix< value_t > *
read_mat_sparse ( std::istream &  in,
                  bool            byteswap,
                  const int       version )
{
    using  real_t = real_type_t< value_t >;

    auto      A = make_unique< TSparseMatrix< value_t > >();
    uint32_t  rows, cols, nnz;

    read_mat_header( in, A.get(), byteswap, version );
    rows = hlib_read<uint32_t>( in, byteswap );
    cols = hlib_read<uint32_t>( in, byteswap );
    nnz  = hlib_read<uint32_t>( in, byteswap );

    A->set_size( rows, cols );
    A->init( nnz );

    for ( uint32_t  i = 0; i <= rows; i++ )
        A->rowptr(i) = hlib_read<int32_t>( in, byteswap );
    
    for ( uint32_t  i = 0; i < nnz; i++ )
        A->colind(i) = hlib_read<int32_t>( in, byteswap );

    if ( sizeof(real_t) == sizeof(double) )
        hlib_read< value_t >( in, & A->coeff(0), nnz, byteswap );
    else
    {
        if ( is_complex_type< value_t >::value )
        {
            for ( idx_t  i = 0; i < idx_t( nnz ); i++ )
            {
                double z[2];

                hlib_read<double>( in, z, 2, byteswap );
                A->coeff(i) = get_value< value_t >::compose( z[0], z[1] );
            }// for
        }// if
        else
        {
            for ( idx_t  i = 0; i < idx_t( nnz ); i++ )
                A->coeff(i) = value_t( hlib_read<double>( in, byteswap ) );
        }// else
    }// else
        
    return A.release();
}

TPermutation *
read_mat_perm  ( std::istream &  in, bool byteswap, const int )
{
    auto      P = make_unique< TPermutation >();
    uint32_t  size;

    size = hlib_read<uint32_t>( in, byteswap );

    P->resize( size );

    for ( uint32_t  i = 0; i < size; i++ )
        (*P)[i] = hlib_read<int32_t>( in, byteswap );
    
    return P.release();
}

template < typename value_t >
TMatrix< value_t > *
read_mat_h ( std::istream &  in,
             bool            byteswap,
             const int       version )
{
    auto      H = make_unique< THMatrix< value_t > >();
    uint32_t  rows, cols;
    uint32_t  block_rows, block_cols;
    bool      has_perm;

    read_mat_header( in, H.get(), byteswap, version );
    rows       = hlib_read<uint32_t>( in, byteswap );
    cols       = hlib_read<uint32_t>( in, byteswap );
    block_rows = hlib_read<uint32_t>( in, byteswap );
    block_cols = hlib_read<uint32_t>( in, byteswap );
    has_perm   = hlib_read<uint8_t>(  in, byteswap ) != 0;

    H->set_size( rows, cols );
    H->set_block_struct( block_rows, block_cols );

    if ( has_perm )
    {
        // TODO: check correct type (permutation)
        {
            auto  P_e2i = std::unique_ptr< TPermutation >( read_mat_perm( in, byteswap, version ) );
            auto  P_i2e = std::unique_ptr< TPermutation >( read_mat_perm( in, byteswap, version ) );
            
            H->set_row_perm( * P_e2i.get(), * P_i2e.get() );
        }

        {
            auto  P_e2i = std::unique_ptr< TPermutation >( read_mat_perm( in, byteswap, version ) );
            auto  P_i2e = std::unique_ptr< TPermutation >( read_mat_perm( in, byteswap, version ) );
            
            H->set_col_perm( * P_e2i.get(), * P_i2e.get() );
        }
    }// if
    
    if ( H->is_nonsym() )
    {
        for ( uint i = 0; i < block_rows; i++ )
            for ( uint j = 0; j < block_cols; j++ )
                H->set_block( i, j, read_mat< value_t >( in, byteswap, version ) );
    }// if
    else
    {
        for ( uint i = 0; i < block_rows; i++ )
            for ( uint j = 0; j <= i; j++ )
                H->set_block( i, j, read_mat< value_t >( in, byteswap, version ) );
    }// else
    
    H->comp_min_max_idx();
    
    return H.release();
}

template < typename value_t >
TZeroMatrix< value_t > *
read_mat_zero  ( std::istream &  in,
                 bool            byteswap,
                 const int       version )
{
    auto      M = make_unique< TZeroMatrix< value_t > >();
    uint32_t  rows, cols;

    read_mat_header( in, M.get(), byteswap, version );
    
    rows = hlib_read<uint32_t>( in, byteswap );
    cols = hlib_read<uint32_t>( in, byteswap );
    M->set_size( rows, cols );
    
    return M.release();
}

template < typename value_t >
TMatrix< value_t > *
read_mat  ( std::istream &  in,
            bool            byteswap,
            const int       version )
{
    uint16_t              type = HMFILE_NULL;
    TMatrix< value_t > *  M    = nullptr;
    
    type = hlib_read<uint16_t>( in, byteswap );

    if      ( type == HMFILE_NULL        ) M = nullptr;
    else if ( type == HMFILE_MAT_H       ) M = read_mat_h< value_t >(       in, byteswap, version );
    else if ( type == HMFILE_MAT_BLOCK   ) M = read_mat_block< value_t >(   in, byteswap, version );
    else if ( type == HMFILE_MAT_DENSE   ) M = read_mat_dense< value_t >(   in, byteswap, version );
    else if ( type == HMFILE_MAT_RANK    ) M = read_mat_rank< value_t >(    in, byteswap, version );
    else if ( type == HMFILE_MAT_SPARSE  ) M = read_mat_sparse< value_t >(  in, byteswap, version );
    else if ( version >= 2 )
    {
        if ( type == HMFILE_MAT_ZERO    ) M = read_mat_zero< value_t >(    in, byteswap, version );
    }// if
    else if ( version == 1 )
    {
        if      ( type == HMFILE_MAT_PERM    ) HERROR( ERR_MAT_TYPE, "(THproMatrixIO) read", "unsupported in v2.0" );
        else if ( type == HMFILE_MAT_FAC     ) HERROR( ERR_MAT_TYPE, "(THproMatrixIO) read", "unsupported in v2.0" );
        else if ( type == HMFILE_MAT_FACINV  ) HERROR( ERR_MAT_TYPE, "(THproMatrixIO) read", "unsupported in v2.0" );
    }// if
    else
        HERROR( ERR_MAT_TYPE, "(THproMatrixIO) read", "" );

    return M;
}

namespace LinOp
{

template < typename value_t >
TLinearOperator< value_t > *
read  ( std::istream &  in,
        bool            byteswap,
        const int       version );

template < typename value_t >
TLinearOperator< value_t > *
read  ( std::istream &  in,
        bool            byteswap,
        const int       version )
{
    uint16_t type = HMFILE_NULL;
    
    type = hlib_read<uint16_t>( in, byteswap );

    if      ( type == HMFILE_NULL        ) return nullptr;
    // Matrices
    else if ( type == HMFILE_MAT_H       ) return read_mat_h< value_t >(        in, byteswap, version );
    else if ( type == HMFILE_MAT_BLOCK   ) return read_mat_block< value_t >(    in, byteswap, version );
    else if ( type == HMFILE_MAT_DENSE   ) return read_mat_dense< value_t >(    in, byteswap, version );
    else if ( type == HMFILE_MAT_RANK    ) return read_mat_rank< value_t >(     in, byteswap, version );
    else if ( type == HMFILE_MAT_SPARSE  ) return read_mat_sparse< value_t >(   in, byteswap, version );
    // linear Operators
    // Rest
    else if ( version >= 2 )
    {
        if      ( type == HMFILE_MAT_ZERO ) return read_mat_zero< value_t >(    in, byteswap, version );
    }// if
    else if ( version == 1 )
    {
        if ( type == HMFILE_MAT_PERM ) HERROR( ERR_MAT_TYPE, "(THproMatrixIO) read", "old permutation unsupported in v2.0" );
    }// if
    else
        HERROR( ERR_MAT_TYPE, "(THproMatrixIO) read", "" );
    
    return nullptr;
}

}// namespace LinOp

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// special write methods for vectors
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

template < typename value_t >
void
write_vec ( std::ostream &              out,
            const TVector< value_t > *  v );

template < typename value_t >
void
write_vec_header ( std::ostream &              out,
                   const TVector< value_t > *  v,
                   const uint                  type )
{
    uint16_t  ltype      = uint16_t(type);
    uint32_t  ofs        = uint32_t(v->ofs());
    uint8_t   is_complex = uint8_t(is_complex_type< value_t >::value);
    uint32_t  size       = uint32_t(v->size());

    out.write( reinterpret_cast< const char * >( & ltype ),      sizeof(uint16_t) );
    out.write( reinterpret_cast< const char * >( & ofs ),        sizeof(uint32_t) );
    out.write( reinterpret_cast< const char * >( & is_complex ), sizeof(uint8_t) );
    out.write( reinterpret_cast< const char * >( & size ),       sizeof(uint32_t) );
}

template < typename value_t >
void
write_vec_block ( std::ostream &                   out,
                  const TBlockVector< value_t > *  v )
{
    uint32_t  blocks = v->n_blocks();
    
    write_vec_header( out, v, HMFILE_VEC_BLOCK );
    
    out.write( reinterpret_cast< const char * >( & blocks ), sizeof(uint32_t) );

    for ( uint i = 0; i < blocks; i++ )
        write_vec( out, v->block(i) );
}

template < typename value_t >
void
write_vec_scalar ( std::ostream &                    out,
                   const TScalarVector< value_t > *  v )
{
    using  real_t = real_type_t< value_t >;

    uint32_t  size = uint32_t(v->size());

    write_vec_header( out, v, HMFILE_VEC_SCALAR );

    if ( sizeof(real_t) == sizeof(double) )
        out.write( reinterpret_cast< const char * >( v->blas_vec().data() ), sizeof(value_t) * size );
    else
    {
        if ( is_complex_type< value_t >::value )
        {
            for ( uint i = 0; i < size; i++ )
            {
                const double z[2] = { std::real( v->entry(i) ), std::imag( v->entry(i) ) };
                
                out.write( reinterpret_cast< const char * >( z ), sizeof(double) * 2 );
            }// for
        }// if
        else
        {
            for ( uint i = 0; i < size; i++ )
            {
                const double  f = std::real( v->entry(i) );
                
                out.write( reinterpret_cast< const char * >( & f ), sizeof(double) );
            }// for
        }// else
    }// else
}

template < typename value_t >
void
write_vec ( std::ostream &              out,
            const TVector< value_t > *  v )
{
    if ( v == nullptr )
    {
        uint16_t  type = HMFILE_NULL;

        out.write( reinterpret_cast< const char * >( & type ), sizeof(uint16_t) );
    }// if
    else if ( IS_TYPE( v, TBlockVector  ) ) write_vec_block(  out, cptrcast( v, TBlockVector< value_t >  ) );
    else if ( IS_TYPE( v, TScalarVector ) ) write_vec_scalar( out, cptrcast( v, TScalarVector< value_t > ) );
    else
        HERROR( ERR_VEC_TYPE, "(THproVectorIO) write", v->typestr() );
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// special read methods for vectors
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

template < typename value_t >
TVector< value_t > *
read_vec ( std::istream &  in,
           const bool      byteswap );

template < typename value_t >
void
read_vec_header ( std::istream &        in,
                  TVector< value_t > *  v,
                  const bool            byteswap )
{
    uint32_t  ofs        = 0;
    bool      is_complex = false;

    ofs        = hlib_read<uint32_t>( in, byteswap );
    is_complex = hlib_read<uint8_t>(  in, byteswap ) != 0;

    v->set_ofs( ofs );

    if ( is_complex && ! is_complex_type< value_t >::value )
        HERROR( ERR_REAL_CMPLX, "(HLibIO) read_vec_header", "complex valued data but real value requested" );
}

template < typename value_t >
TVector< value_t > *
read_vec_block ( std::istream &  in,
                 const bool      byteswap )
{
    auto      v = new TBlockVector< value_t >();
    uint32_t  blocks;

    read_vec_header( in, v, byteswap );
    hlib_read<uint32_t>( in, byteswap ); // overread size
    blocks = hlib_read<uint32_t>( in, byteswap );

    v->set_block_struct( blocks );

    for ( uint i = 0; i < blocks; i++ )
        v->set_block( i, read_vec< value_t >( in, byteswap ) );
    
    return v;
}

template < typename value_t >
TVector< value_t > *
read_vec_scalar ( std::istream &  in,
                  const bool      byteswap )
{
    using  real_t = real_type_t< value_t >;

    auto      v = new TScalarVector< value_t >();
    uint32_t  size;

    read_vec_header( in, v, byteswap );
    size = hlib_read<uint32_t>( in, byteswap );

    v->set_size( size );

    if ( sizeof(real_t) == sizeof(double) )
        hlib_read< value_t >( in, v->blas_vec().data(), size, byteswap );
    else
    {
        if ( is_complex_type< value_t >::value )
        {
            for ( idx_t  i = 0; i < idx_t(size); i++ )
            {
                double  z[2];

                hlib_read<double>( in, z, 2, byteswap );
                v->blas_vec()(i) = get_value< value_t >::compose( z[0], z[1] );
            }// for
        }// if
        else
        {
            for ( idx_t  i = 0; i < idx_t(size); i++ )
                v->blas_vec()(i) = value_t( hlib_read<double>( in, byteswap ) );
        }// else
    }// else
    
    return v;
}

template < typename value_t >
TVector< value_t > *
read_vec ( std::istream &  in,
           const bool      byteswap )
{
    uint16_t type = HMFILE_NULL;
    
    type = hlib_read<uint16_t>( in, byteswap );

    if      ( type == HMFILE_NULL       ) return nullptr;
    else if ( type == HMFILE_VEC_BLOCK  ) return read_vec_block< value_t >(  in, byteswap );
    else if ( type == HMFILE_VEC_SCALAR ) return read_vec_scalar< value_t >( in, byteswap );
    else
        HERROR( ERR_VEC_TYPE, "(THproVectorIO) read", "" );
    
    return nullptr;
}

}// namespace anonymous

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// class THproMatrixIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// write matrix with name
//
template < typename value_t >
void
THproMatrixIO::write ( const TMatrix< value_t > *  A,
                       const std::string &         filename ) const
{
    if ( A == nullptr )
        HERROR( ERR_ARG, "(THproMatrixIO) write", "argument is nullptr" );

    auto            out_ptr = open_write( filename );
    std::ostream &  out     = * out_ptr.get();

    write_header( out );

    //
    // write matrix
    //
    
    Hpro::write_mat( out, A );
}

//
// read matrix from file
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
THproMatrixIO::read ( const std::string &  filename ) const
{
    if ( ! boost::filesystem::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(THproMatrixIO) read", filename );
    
    auto            in_ptr = open_read( filename );
    std::istream &  in     = * in_ptr.get();
    
    //
    // read (or skip) header
    //
    
    bool      byteswap = false;
    uint16_t  version;

    read_header( in, byteswap, version );
    
    return std::unique_ptr< TMatrix< value_t > >( Hpro::read_mat< value_t >( in, byteswap, version ) );
}

//
// write matrix with name
//
template < typename value_t >
void
THproMatrixIO::write_linop ( const TLinearOperator< value_t > *  A,
                             const std::string &                 filename ) const
{
    if ( A == nullptr )
        HERROR( ERR_ARG, "(THproMatrixIO) write_linop", "argument is nullptr" );
    
    auto            out_ptr = open_write( filename );
    std::ostream &  out     = * out_ptr.get();

    write_header( out );

    //
    // write matrix
    //
    
    LinOp::write( out, A );
}

//
// read matrix from file
//
template < typename value_t >
std::unique_ptr< TLinearOperator< value_t > >
THproMatrixIO::read_linop ( const std::string &  filename ) const
{
    auto            in_ptr = open_read( filename );
    std::istream &  in     = * in_ptr.get();
    
    //
    // read (or skip) header
    //
    
    bool      byteswap = false;
    uint16_t  version;

    read_header( in, byteswap, version );
    
    return std::unique_ptr< TLinearOperator< value_t > >( LinOp::read< value_t >( in, byteswap, version ) );
}

#define INST_MATRIXIO( type )                   \
    template void THproMatrixIO::write< type > ( const TMatrix< type > *, const std::string & ) const; \
    template std::unique_ptr< TMatrix< type > > THproMatrixIO::read< type > ( const std::string & ) const; \
    template void THproMatrixIO::write_linop< type > ( const TLinearOperator< type > *, const std::string & ) const; \
    template std::unique_ptr< TLinearOperator< type > > THproMatrixIO::read_linop< type > ( const std::string & ) const;

INST_MATRIXIO( float )
INST_MATRIXIO( double )
INST_MATRIXIO( std::complex< float > )
INST_MATRIXIO( std::complex< double > )

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// class THproVectorIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// write vector to file <filename>
//
template < typename value_t >
void
THproVectorIO::write ( const TVector< value_t > *  v,
                       const std::string &         filename ) const
{
    if ( v == nullptr )
        HERROR( ERR_ARG, "(THproVectorIO) write_vector", "vector is nullptr" );
    
    auto            out_ptr = open_write( filename );
    std::ostream &  out     = * out_ptr.get();

    write_header( out );

    //
    // write vector
    //

    Hpro::write_vec( out, v );
}

//
// read vector from file <filename>
//
template < typename value_t >
std::unique_ptr< TVector< value_t > >
THproVectorIO::read ( const std::string & filename ) const
{
    if ( ! boost::filesystem::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(THproVectorIO) read", filename );
    
    auto            in_ptr = open_read( filename );
    std::istream &  in     = * in_ptr.get();

    //
    // read (or skip) header
    //
    
    bool      byteswap = false;
    uint16_t  version;

    read_header( in, byteswap, version );
    
    return std::unique_ptr< TVector< value_t > >( Hpro::read_vec< value_t >( in, byteswap ) );
}

#define INST_VECTORIO( type )                   \
    template void THproVectorIO::write< type > ( const TVector< type > *, const std::string & ) const; \
    template std::unique_ptr< TVector< type > > THproVectorIO::read< type > ( const std::string & ) const;

INST_VECTORIO( float )
INST_VECTORIO( double )
INST_VECTORIO( std::complex< float > )
INST_VECTORIO( std::complex< double > )

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// class THproCoordIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// write vertices to <filename>
//
void
THproCoordIO::write ( const TCoordinate *  coo,
                      const std::string &  filename ) const
{
    if ( coo == nullptr )
        HERROR( ERR_ARG, "(THproCoordIO) write", "coordinates are nullptr" );
    
    auto            out_ptr = open_write( filename );
    std::ostream &  out     = * out_ptr.get();

    write_header( out );

    //
    // write <nvertices>, <dim> and <type>
    //

    const uint16_t type      = HMFILE_COORD_VTX;
    const uint32_t nvertices = uint32_t(coo->ncoord());
    const uint16_t tdim      = uint16_t(coo->dim());
    
    out.write( reinterpret_cast< const char * >( & type ), sizeof(type) );
    out.write( reinterpret_cast< const char * >( & nvertices ), sizeof(nvertices) );
    out.write( reinterpret_cast< const char * >( & tdim ), sizeof(tdim) );

    //
    // write coordinates
    //

    for ( uint i = 0; i < nvertices; i++ )
        out.write( reinterpret_cast< const char * >( coo->coord(i) ), coo->dim() * sizeof(double) );
}

//
// read vertices from <filename>
//
std::unique_ptr< TCoordinate >
THproCoordIO::read ( const std::string & filename ) const
{
    if ( ! boost::filesystem::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(THproCoordIO) read", filename );
    
    auto            in_ptr = open_read( filename );
    std::istream &  in     = * in_ptr.get();

    //
    // read (or skip) header
    //
    
    bool      byteswap = false;
    uint16_t  version;
    
    read_header( in, byteswap, version );
    
    //
    // read number of coordinates, spatial dimension and information
    // about stored data
    //

    const uint16_t  type      = hlib_read<uint16_t>( in, byteswap );
    const uint32_t  nvertices = hlib_read<uint32_t>( in, byteswap );
    const uint16_t  dim       = hlib_read<uint16_t>( in, byteswap );

    if ( type != HMFILE_COORD_VTX )
        HERROR( ERR_NOT_IMPL, "(THproCoordIO) read", "only supporting pure vertex data" );
    
    //
    // now read coordinates according to <type> info
    //

    std::vector< double * >  vertices( nvertices, nullptr );

    for ( uint i = 0; i < nvertices; i++ )
    {
        if ( type == HMFILE_COORD_VTX )
        {
            double * coord = new double[dim];

            hlib_read<double>( in, coord, dim, byteswap );

            vertices[i] = coord;
        }// if
    }// for

    return make_unique< TCoordinate >( vertices, dim );
}

/////////////////////////////////////////////////////////////////
//
// explicit template instantiation
//
/////////////////////////////////////////////////////////////////

variant_id_t
hpro_guess_value_type ( const std::string &  filename )
{
    if ( ! boost::filesystem::exists( filename ) )
        HERROR( ERR_FNEXISTS, "hpro_guess_value_type", filename );
    
    auto            in_ptr = open_read( filename );
    std::istream &  in     = * in_ptr.get();
    
    //
    // read (or skip) header
    //
    
    bool      byteswap = false;
    uint16_t  version;

    read_header( in, byteswap, version );

    auto  type = hlib_read<uint16_t>( in, byteswap );
    
    bool  is_complex = false;

    if (( type == 100 ) || ( type == 101 ))
    {
        //
        // vectors
        //
        
        hlib_read<uint32_t>( in, byteswap );
        is_complex = hlib_read<uint8_t>(  in, byteswap ) != 0;
    }// if
    else if (( type > 0 ) && ( type < 100 ))
    {
        //
        // matrices
        //
        
        if ( version >= 2 )
        {
            hlib_read<int32_t>(  in, byteswap );
            hlib_read<uint32_t>( in, byteswap );
            hlib_read<uint32_t>( in, byteswap );
            hlib_read<uint32_t>( in, byteswap );
            hlib_read<uint32_t>( in, byteswap );
            hlib_read<uint8_t>(  in, byteswap );
            is_complex = hlib_read<uint8_t>( in, byteswap ) != 0;
        }// if
        else
        {
            hlib_read<uint32_t>( in, byteswap );
            hlib_read<uint32_t>( in, byteswap );
            hlib_read<uint8_t>(  in, byteswap );
            is_complex = hlib_read<uint8_t>( in, byteswap ) != 0;
        }// else
    }// else

    if ( is_complex )
        return COMPLEX_FP64;
    else
        return REAL_FP64;
}

}// namespace Hpro
