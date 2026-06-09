//
// Project     : HLIBpro
// File        : Matrix.cc
// Description : implements dense matrix class for BLAS operations
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

#include "hpro/base/config.hh"
#include "hpro/base/System.hh"

#include "TRNG.hh"

#include "hpro/io/TMatrixIO.hh"
#include "hpro/io/TVectorIO.hh"

#include "hpro/blas/Vector.hh"

#include "hpro/blas/Matrix.hh"

namespace Hpro
{

namespace BLAS
{

// special Range type for "whole" indexset
Range  Range::all( -1, -1 );

//
// test data for invalid values, e.g. INF and NAN
//
template < typename value_t >
void
Matrix< value_t >::check_data  () const
{
    for ( idx_t  j = 0; j < idx_t( ncols() ); ++j )
    {
        for ( idx_t  i = 0; i < idx_t( nrows() ); ++i )
        {
            const value_t  val = (*this)(i,j);
            
            if ( Math::is_inf( val ) )
                HERROR( ERR_INF, "(Matrix) check_data", "" );
            
            if ( Math::is_nan( val ) )
                HERROR( ERR_NAN, "(Matrix) check_data", "" );
        }// for
    }// for
}

//
// stream output
//
template < typename T >
std::ostream &
operator << ( std::ostream &           os,
              const MatrixBase< T > &  M )
{
    if ( CFG::IO::use_matlab_syntax )
        os << '[';
    
    for ( idx_t  i = 0; i < idx_t( M.nrows() ); ++i )
    {
        for ( idx_t  j = 0; j < idx_t( M.ncols() )-1; ++j )
        {
            os << M(i,j) << ", ";
        }// for

        if ( M.ncols() > 0 )
        {
            if ( CFG::IO::use_matlab_syntax )
                os << M(i,idx_t(M.ncols())-1) << ';';
            os << std::endl;
        }// if
    }// for
        
    if ( CFG::IO::use_matlab_syntax )
        os << ']';
    
    return os;
}

template < typename T >
std::ostream &
operator << ( std::ostream &           os,
              const VectorBase< T > &  v )
{
    if ( CFG::IO::use_matlab_syntax )
        os << '[';
    
    for ( idx_t  i = 0; i < idx_t(v.length())-1; i++ )
        os << v(i) << ", ";
    
    if ( v.length() > 0 )
        os << v( idx_t(v.length())-1 );
    
    if ( CFG::IO::use_matlab_syntax )
        os << "]'";
    
    return os;
}

//
// create identity matrix
//
template < typename T >
Matrix< T >
identity ( const size_t  n )
{
    Matrix< T >  M( n, n );

    for ( idx_t  i = 0; i < idx_t( n ); i++ )
        M( i, i ) = T(1);

    return M;
}

//
// create random matrix
//
template < typename T >
Matrix< T >
random ( const size_t  nrows,
         const size_t  ncols )
{
    TRNG         rng;
    Matrix< T >  M( nrows, ncols );

    for ( idx_t  j = 0; j < idx_t( ncols ); j++ )
        for ( idx_t  i = 0; i < idx_t( nrows ); i++ )
            M( i, j ) = rng.rand< T >( typename real_type< T >::type_t( 1 ) );

    return M;
}

//
// create random vector
//
template < typename T >
Vector< T >
random ( const size_t  length )
{
    TRNG         rng;
    Vector< T >  v( length );

    for ( idx_t  i = 0; i < idx_t( length ); ++i )
        v( i ) = rng.rand< T >( typename real_type< T >::type_t( 1 ) );

    return v;
}

//
// explicit template instantiation
//

template class Matrix< float >;
template class Matrix< double >;
template class Matrix< std::complex< float > >;
template class Matrix< std::complex< double > >;


template std::ostream &
operator << < Matrix< float > >             ( std::ostream &,
                                              const MatrixBase< Matrix< float > > & );
template std::ostream &
operator << < Matrix< double > >            ( std::ostream &,
                                              const MatrixBase< Matrix< double > > & );

template std::ostream &
operator << < Matrix< std::complex< float > > >  ( std::ostream &,
                                                   const MatrixBase< Matrix< std::complex< float > > > & );
template std::ostream &
operator << < Matrix< std::complex< double > > > ( std::ostream &,
                                                   const MatrixBase< Matrix< std::complex< double > > > & );


template std::ostream &
operator << < Vector< float > >             ( std::ostream &,
                                              const VectorBase< Vector< float > > & );
template std::ostream &
operator << < Vector< double > >            ( std::ostream &,
                                              const VectorBase< Vector< double > > & );

template std::ostream &
operator << < Vector< std::complex< float > > >  ( std::ostream &,
                                                   const VectorBase< Vector< std::complex< float > > > & );
template std::ostream &
operator << < Vector< std::complex< double > > > ( std::ostream &,
                                                   const VectorBase< Vector< std::complex< double > > > & );


template Matrix< float >
identity< float >             ( const size_t  n );
    
template Matrix< double >
identity< double >            ( const size_t  n );

template Matrix< std::complex< float > >
identity< std::complex< float > >  ( const size_t  n );
    
template Matrix< std::complex< double > >
identity< std::complex< double > > ( const size_t  n );


template Matrix< float >
random< float >             ( const size_t  nrows,
                              const size_t  ncols );
    
template Matrix< double >
random< double >            ( const size_t  nrows,
                              const size_t  ncols );
    
template Matrix< std::complex< float > >
random< std::complex< float > >  ( const size_t  nrows,
                                   const size_t  ncols );
    
template Matrix< std::complex< double > >
random< std::complex< double > > ( const size_t  nrows,
                                   const size_t  ncols );
    
template Vector< float >
random< float >             ( const size_t  length );
    
template Vector< double >
random< double >            ( const size_t  length );
    
template Vector< std::complex< float > >
random< std::complex< float > >  ( const size_t  length );
    
template Vector< std::complex< double > >
random< std::complex< double > > ( const size_t  length );
    
}// namespace BLAS

///////////////////////////////////////////////////////////////////////////
//
// matrix permutation functions
//
///////////////////////////////////////////////////////////////////////////

namespace BLAS
{

namespace
{

//
// copy row/column <j> to row/column <i> in matrix M
//
template <typename value_t>
void
copy ( Matrix< value_t > &  M,
       const bool           sort_rows,
       const idx_t          i,
       const idx_t          j )
{
    if ( sort_rows )
    {
        auto  row_i = M.row( i );
        auto  row_j = M.row( j );
            
        copy( row_j, row_i );
    }// if
    else
    {
        auto  col_i = M.column( i );
        auto  col_j = M.column( j );
            
        copy( col_j, col_i );
    }// else
}

//
// copy row/column <i> to vector <tmp>
//
template <typename T>
void
copy ( Matrix< T > &  M,
       const bool     sort_rows,
       const idx_t    i,
       Vector< T > &  tmp )
{
    if ( sort_rows )
    {
        auto  row_i = M.row( i );
            
        copy( row_i, tmp );
    }// if
    else
    {
        auto  col_i = M.column( i );
            
        copy( col_i, tmp );
    }// else
}

//
// copy vector <tmp> to row/column <i>
//
template <typename T>
void
copy ( Vector< T > &  tmp,
       const bool     sort_rows,
       const idx_t    i,
       Matrix< T > &  M )
{
    if ( sort_rows )
    {
        auto  row_i = M.row( i );
            
        copy( tmp, row_i );
    }// if
    else
    {
        auto  col_i = M.column( i );
            
        copy( tmp, col_i );
    }// else
}

//
// swap the rows/columns <i> and <j> in matrix M
//
template <typename T>
void
swap ( Matrix< T > &  M,
       const bool     sort_rows,
       const idx_t    i,
       const idx_t    j,
       Vector< T > &  tmp )
{
    if ( sort_rows )
    {
        auto  row_i = M.row( i );
        auto  row_j = M.row( j );

        copy( row_i, tmp );
        copy( row_j, row_i );
        copy( tmp,   row_j );
    }// if
    else
    {
        auto  col_i = M.column( i );
        auto  col_j = M.column( j );

        copy( col_i, tmp );
        copy( col_j, col_i );
        copy( tmp,   col_j );
    }// else
}

//
// sort rows/columns of matrix according to given permutation
//
template <typename T>
void
matrix_sort( Matrix< T > &   M,
             const bool      sort_rows,
             TPermutation &  perm,
             const idx_t     lb,
             const idx_t     ub,
             Vector< T > &   tmp )
{
    if ( lb >= ub ) return;

    if ( (ub - lb) < 20 )
    {
        //
        // apply insertion sort for small ranges
        //

        for ( idx_t  i = lb+1; i <= ub; i++ )
        {
            const idx_t  v = perm[i];
            idx_t        j = i-1;

            copy( M, sort_rows, i, tmp );
            
            while (( j >= 0 ) && ( perm[j] > v ))
            {
                copy( M, sort_rows, j+1, j );
                perm[j+1] = perm[j];
                j--;
            }// if

            copy( tmp, sort_rows, j+1, M );
            perm[j+1] = v;
        }// for
    }// if
    else
    {
        //
        // apply quick sort for standard ranges
        //

        idx_t        i         = lb;
        idx_t        j         = ub;
        const idx_t  mid       = (lb + ub) / 2;
        idx_t        choice[3] = { perm[lb], perm[mid], perm[ub] };
        idx_t        pivot;

        // choose pivot (median-of-three)
        if ( choice[0] > choice[1] ) std::swap( choice[0], choice[1] );
        if ( choice[0] > choice[2] ) std::swap( choice[0], choice[2] );
        if ( choice[1] > choice[2] ) std::swap( choice[1], choice[2] );
        pivot = choice[1];

        // partition
        while ( i < j )
        {
            while ( perm[i] < pivot   ) i++;
            while ( pivot   < perm[j] ) j--;

            if ( i < j )
            {
                swap( M, sort_rows, i, j, tmp );
                std::swap( perm[i], perm[j] );
            }// if
        }// while

        // recursion
        matrix_sort( M, sort_rows, perm, lb, i-1, tmp );
        matrix_sort( M, sort_rows, perm, i+1, ub, tmp );
    }// else
}

}// namespace anonymous

//
// permute rows/columns in matrix
//
template < typename value_t >
void
permute ( Matrix< value_t > &   M,
          const TPermutation &  row_perm,
          const TPermutation &  col_perm )
{
#if 1

    //
    // apply permutations by sorting the rows/columns
    // simultaneously to sorting permutations, thereby
    // reverting order
    //
    
    TPermutation       perm;
    Vector< value_t >  tmp;

    {
        tmp  = Vector< value_t >( M.ncols() );
        perm = row_perm;
        
        matrix_sort( M, true, perm, 0, idx_t(M.nrows())-1, tmp );
    }// if
    
    {
        tmp  = Vector< value_t >( M.nrows() );
        perm = col_perm;
        
        matrix_sort( M, false, perm, 0, idx_t(M.ncols())-1, tmp );
    }// if
        
#else
    
    auto        T = copy( M );
    const uint  n = T.nrows();
    const uint  m = T.ncols();
    
    //
    // apply row permutation first by copying the old rows
    // to the new position defined by the permutation
    //

    {
        for ( uint i = 0; i < n; i++ )
        {
            const uint  pi = row_perm.permute( i );
            
            auto  M_i  = M.row( i );
            auto  T_pi = T.row( pi );
            
            copy( M_i, T_pi );
        }// for

        copy( T, M );
    }// if

    //
    // now apply column permutation in the same way
    //

    {
        for ( uint j = 0; j < m; j++ )
        {
            const uint  pj = col_perm.permute( j );

            auto  M_j  = M.row( j );
            auto  T_pj = T.row( pj );
            
            copy( M_j, T_pj );
        }// for

        copy( T, M );
    }// if
    
#endif
}

template < typename value_t >
void
permute ( Vector< value_t > &   v,
          const TPermutation &  perm )
{
    auto  t = v.copy();

    for ( size_t  i = 0; i < v.length(); ++i )
        v( perm[i] ) = t(i);
}

template < typename value_t >
void
permute_inv ( Vector< value_t > &   v,
              const TPermutation &  perm )
{
    auto  t = v.copy();

    for ( size_t  i = 0; i < v.length(); ++i )
        v( i ) = t( perm[i] );
}

template void permute< float >                       ( BLAS::Matrix< float > &,                  const TPermutation &,  const TPermutation & );
template void permute< double >                      ( BLAS::Matrix< double > &,                 const TPermutation &,  const TPermutation & );
template void permute< std::complex< float > >       ( BLAS::Matrix< std::complex< float > > &,  const TPermutation &,  const TPermutation & );
template void permute< std::complex< double > >      ( BLAS::Matrix< std::complex< double > > &, const TPermutation &,  const TPermutation & );

template void permute< float >                       ( BLAS::Vector< float > &,                  const TPermutation & );
template void permute< double >                      ( BLAS::Vector< double > &,                 const TPermutation & );
template void permute< std::complex< float > >       ( BLAS::Vector< std::complex< float > > &,  const TPermutation & );
template void permute< std::complex< double > >      ( BLAS::Vector< std::complex< double > > &, const TPermutation & );

template void permute_inv< float >                   ( BLAS::Vector< float > &,                  const TPermutation & );
template void permute_inv< double >                  ( BLAS::Vector< double > &,                 const TPermutation & );
template void permute_inv< std::complex< float > >   ( BLAS::Vector< std::complex< float > > &,  const TPermutation & );
template void permute_inv< std::complex< double > >  ( BLAS::Vector< std::complex< double > > &, const TPermutation & );

}// namespace BLAS

//////////////////////////////////////////////////////////
//
// debug helpers
//

namespace DBG
{

//!
//! write matrix to file
//!
template < typename T >
void
write ( const BLAS::Matrix< T > & M,
        const std::string &       filename,
        const std::string &       matname )
{
    TMatlabMatrixIO  mio;
    
    mio.write( M, filename, matname );
}

// instantiate
template void write< float >   ( const BLAS::Matrix< float > &,  const std::string &, const std::string & );
template void write< double >  ( const BLAS::Matrix< double > &, const std::string &, const std::string & );

template void write< std::complex< float > >   ( const BLAS::Matrix< std::complex< float > > &,  const std::string &, const std::string & );
template void write< std::complex< double > >  ( const BLAS::Matrix< std::complex< double > > &, const std::string &, const std::string & );

//!
//! write vector to file
//!
template < typename T >
void
write ( const BLAS::Vector< T > & v,
        const std::string &       filename,
        const std::string &       matname )
{
    TMatlabVectorIO  vio;
    
    vio.write( v, filename, matname );
}

// instantiate
template void write< float >    ( const BLAS::Vector< float > &,   const std::string &, const std::string & );
template void write< double >   ( const BLAS::Vector< double > &,  const std::string &, const std::string & );

template void write< std::complex< float > >   ( const BLAS::Vector< std::complex< float > > &,  const std::string &, const std::string & );
template void write< std::complex< double > >  ( const BLAS::Vector< std::complex< double > > &, const std::string &, const std::string & );

}// namespace DBG

}// namespace Hpro
