#ifndef __HPRO_STRUCTURE_HH
#define __HPRO_STRUCTURE_HH
//
// Project     : HLIBpro
// File        : structure.hh
// Description : test functions for various matrix structures
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"
#include "hpro/matrix/TSparseMatrix.hh"
#include "hpro/matrix/THMatrix.hh"
#include "hpro/matrix/TGhostMatrix.hh"

namespace Hpro
{

//!
//! return true, if given op is a matrix
//!
template < typename value_t > bool is_matrix       ( const TLinearOperator< value_t > *  M ) noexcept { return IS_TYPE( M, TMatrix ); }
template < typename value_t > bool is_matrix       ( const TLinearOperator< value_t > &  M ) noexcept { return is_matrix( &M ); }

//!
//! return true, if matrix is an H-matrix
//!
template < typename value_t > bool is_hmat         ( const TMatrix< value_t > *  A ) noexcept { return IS_TYPE(  A, THMatrix ); }
template < typename value_t > bool is_hmat         ( const TMatrix< value_t > &  A ) noexcept { return is_hmat( &A ); }

//!
//! return true, if matrix is a block matrix
//!
template < typename value_t > bool is_blocked      ( const TMatrix< value_t > *  A ) noexcept { return (A != nullptr) && A->is_blocked(); }
template < typename value_t > bool is_blocked      ( const TMatrix< value_t > &  A ) noexcept { return A.is_blocked(); }

//!
//! return true, if all matrices are blocked
//!
template < typename value_t > bool is_blocked_all  ( const TMatrix< value_t > *  A ) noexcept { return is_blocked( A ); }
template < typename value_t > bool is_blocked_all  ( const TMatrix< value_t > &  A ) noexcept { return is_blocked( A ); }

template < typename value_t, typename... T >
bool is_blocked_all  ( const TMatrix< value_t > *  A, T...  mats ) { return is_blocked( A ) && is_blocked_all( mats... ); }
template < typename value_t, typename... T >
bool is_blocked_all  ( const TMatrix< value_t > &  A, T...  mats ) { return is_blocked( A ) && is_blocked_all( mats... ); }


//!
//! return true, if matrix is a dense matrix
//!
template < typename value_t >
bool is_dense  ( const TMatrix< value_t > *  A ) noexcept { return IS_TYPE( A, TDenseMatrix ); }

template < typename value_t >
bool is_dense  ( const TMatrix< value_t > &  A ) noexcept { return is_dense( &A ); }

//!
//! return true, if matrix is a lowrank matrix
//!
template < typename value_t >
bool is_lowrank  ( const TMatrix< value_t > *  A ) noexcept { return IS_TYPE( A, TRkMatrix ); }

template < typename value_t >
bool is_lowrank  ( const TMatrix< value_t > &  A ) noexcept { return is_lowrank( &A ); }

//!
//! return true, if matrix is a sparse matrix
//!
template < typename value_t >
bool is_sparse  ( const TMatrix< value_t > *  A ) noexcept { return IS_TYPE( A, TSparseMatrix ); }

template < typename value_t >
bool is_sparse  ( const TMatrix< value_t > &  A ) noexcept { return is_sparse( &A ); }

//!
//! return true, if matrix is a ghost matrix
//!
template < typename value_t >
bool is_ghost ( const TMatrix< value_t > *  A ) noexcept { return IS_TYPE( A, TGhostMatrix ); }

template < typename value_t >
bool is_ghost ( const TMatrix< value_t > &  A ) noexcept { return is_ghost( &A ); }

//!
//! return true if \a A is block diagonal
//!
template < typename value_t >
bool is_diag  ( const TMatrix< value_t > *  A );

template < typename value_t >
bool is_diag  ( const TMatrix< value_t > &  A ) { return is_diag( &A ); }

//!
//! return true if \a A has domain-decomposition format, e.g.
//! only block diagonal and non-zero last block row/column
//!
template < typename value_t >
bool is_dd  ( const TMatrix< value_t > *  A );

template < typename value_t >
bool is_dd  ( const TMatrix< value_t > &  A ) { return is_dd( &A ); }

//!
//! return true if \a A is block matrix with flat hierarchy
//!
template < typename value_t >
bool is_flat  ( const TMatrix< value_t > *  A );

template < typename value_t >
bool is_flat  ( const TMatrix< value_t > &  A ) { return is_flat( &A ); }

//!
//! return true if \a A is a diagonal block
//!
template < typename value_t >
bool is_on_diag  ( const TMatrix< value_t > *  A ) noexcept { return ( A != nullptr ) && ( A->row_is() == A->col_is() ); }

template < typename value_t >
bool is_on_diag  ( const TMatrix< value_t > &  A ) noexcept { return A.row_is() == A.col_is(); }

//!
//! return true if \a A is in lower left part of matrix
//!
template < typename value_t >
bool is_in_lower_left ( const TMatrix< value_t > *  A ) noexcept { return ( A != nullptr ) && A->row_is().is_right_or_equal_to( A->col_is() ); }

template < typename value_t >
bool is_in_lower_left ( const TMatrix< value_t > &  A ) noexcept { return is_in_lower_lef( &A ); }

//!
//! return true if \a A is in lower left part of matrix
//!
template < typename value_t >
bool is_strictly_in_lower_left ( const TMatrix< value_t > *  A ) noexcept { return ( A != nullptr ) && A->row_is().is_strictly_right_of( A->col_is() ); }

template < typename value_t >
bool is_strictly_in_lower_left ( const TMatrix< value_t > &  A ) noexcept { return is_strictly_in_lower_left( &A ); }

//!
//! return true if \a A is in upper right part of matrix
//!
template < typename value_t >
bool is_in_upper_right ( const TMatrix< value_t > *  A ) noexcept { return ( A != nullptr ) && A->col_is().is_right_or_equal_to( A->row_is() ); }

template < typename value_t >
bool is_in_upper_right ( const TMatrix< value_t > &  A ) noexcept { return is_in_upper_right( &A ); }

//!
//! return true if \a A is in upper right part of matrix
//!
template < typename value_t >
bool is_strictly_in_upper_right ( const TMatrix< value_t > *  A ) noexcept { return ( A != nullptr ) && A->col_is().is_strictly_right_of( A->row_is() ); }

//!
//! return true if \a A is lower left block triangular matrix
//!
template < typename value_t >
bool is_lower_left  ( const TMatrix< value_t > *  A );

template < typename value_t >
bool is_lower_left  ( const TMatrix< value_t > &  A ) { return is_lower_left( &A ); }

//!
//! return true if \a A is upper right block triangular matrix
//!
template < typename value_t >
bool is_upper_right  ( const TMatrix< value_t > *  A );

template < typename value_t >
bool is_upper_right  ( const TMatrix< value_t > &  A ) { return is_upper_right( &A ); }

//
// return true if \a A is corresponding to leaf block
//
template < typename value_t >
bool is_leaf  ( const TMatrix< value_t > *  A ) noexcept { return ( A != nullptr ) && ! A->is_blocked(); }

template < typename value_t >
bool is_leaf  ( const TMatrix< value_t > &  A ) noexcept { return ! A.is_blocked(); }

//!
//! return true, if matrix is considered too small for parallel treatement
//!
template < typename value_t >
bool is_small  ( const TMatrix< value_t > *  A ) noexcept { return ( A != nullptr ) && ( std::min( A->rows(), A->cols() ) <= CFG::Arith::max_seq_size ); }

template < typename value_t >
bool is_small  ( const TMatrix< value_t > &  A ) noexcept { return std::min( A.rows(), A.cols() ) <= CFG::Arith::max_seq_size; }

//!
//! return true, if any/all matrix is considered too small for parallel treatement
//!
template < typename value_t >
bool is_small_any  ( const TMatrix< value_t > *  A ) noexcept { return is_small( A ); }

template < typename value_t >
bool is_small_any  ( const TMatrix< value_t > &  A ) noexcept { return is_small( A ); }

template < typename value_t, typename... T >
bool is_small_any  ( const TMatrix< value_t > *  A, T...  mats ) { return is_small( A ) || is_small_any( mats... ); }

template < typename value_t, typename... T >
bool is_small_any  ( const TMatrix< value_t > &  A, T...  mats ) { return is_small( A ) || is_small_any( mats... ); }

//!
//! return true if matrix is zero, e.g. all entries zero
//!
template < typename value_t >
bool is_zero  ( const TMatrix< value_t > *  A ); 

template < typename value_t >
bool is_zero  ( const TMatrix< value_t > &  A ) { return is_zero( &A ); }

//!
//! return maximal ID in all submatrices
//!
template < typename value_t >
int max_id  ( const TMatrix< value_t > *  A );

template < typename value_t >
int max_id  ( const TMatrix< value_t > &  A ) { return max_id( &A ); }

//
// return first level on diagonal on which leaf blocks are found
//
template < typename value_t >
uint get_first_diag_leaf_level ( const TMatrix< value_t > *  A );

template < typename value_t >
uint get_first_diag_leaf_level ( const TMatrix< value_t > &  A ) { return get_first_diag_leaf_level( &A ); }

//!
//! return matrix with rowis( A_row ) × colis( A_col )
//!
template < typename value_t >
TMatrix< value_t > * 
product_block  ( const TMatrix< value_t > *  A_row,
                 const TMatrix< value_t > *  A_col );

template < typename value_t >
TMatrix< value_t > * 
product_block  ( const TMatrix< value_t > &  A_row,
                 const TMatrix< value_t > &  A_col )
{
    return product_block( &A_row, &A_col );
}

//!
//! return matrix with rowis( op( A_row ) ) × colis( op( A_col ) )
//!
template < typename value_t >
TMatrix< value_t > *
product_block  ( const matop_t               op_row,
                 const TMatrix< value_t > *  A_row,
                 const matop_t               op_col,
                 const TMatrix< value_t > *  A_col );

template < typename value_t >
TMatrix< value_t > *
product_block  ( const matop_t               op_row,
                 const TMatrix< value_t > &  A_row,
                 const matop_t               op_col,
                 const TMatrix< value_t > &  A_col )
{
    return product_block( op_row, &A_row, op_col, &A_col );
}

//!
//! return sub block of M with given id
//!
template < typename value_t >
const TMatrix< value_t > *
get_block  ( const TMatrix< value_t > *  M,
             const int                   id );

template < typename value_t >
const TMatrix< value_t > *
get_block  ( const TMatrix< value_t > &  M,
             const int                   id )
{
    return get_block( &M, id );
}

//!
//! return sub block of M with given block index set
//!
template < typename value_t >
const TMatrix< value_t > *
get_block  ( const TMatrix< value_t > *  M,
             const TBlockIndexSet &      bis );

template < typename value_t >
const TMatrix< value_t > *
get_block  ( const TMatrix< value_t > &  M,
             const TBlockIndexSet &      bis )
{
    return get_block( &M, bis );
}

//!
//! return number of sub blocks in matrix (including inner blocks)
//!
template < typename value_t >
size_t get_nblocks  ( const TMatrix< value_t > *  M );

template < typename value_t >
size_t get_nblocks  ( const TMatrix< value_t > &  M ) { return get_nblocks( &M ); }

//!
//! return number of leaf blocks in matrix
//!
template < typename value_t >
size_t get_nblocks_leaf  ( const TMatrix< value_t > *  M );

template < typename value_t >
size_t get_nblocks_leaf  ( const TMatrix< value_t > &  M ) { return get_nblocks_leaf( &M ); }

//!
//! return number of low-rank blocks in matrix
//!
template < typename value_t >
size_t get_nblocks_lowrank  ( const TMatrix< value_t > *  M );

template < typename value_t >
size_t get_nblocks_lowrank  ( const TMatrix< value_t > &  M ) { return get_nblocks_lowrank( &M ); }

//!
//! return number of dense blocks in matrix
//!
template < typename value_t >
size_t get_nblocks_dense  ( const TMatrix< value_t > *  M );

template < typename value_t >
size_t get_nblocks_dense  ( const TMatrix< value_t > &  M ) { return get_nblocks_dense( &M ); }

//!
//! return number of level in matrix
//!
template < typename value_t >
size_t get_nlevel         ( const TMatrix< value_t > &  M );

template < typename value_t >
size_t get_nlevel         ( const TMatrix< value_t > *  M )
{
    if ( M == nullptr )
        return 0;
    
    return get_nlevel( *M );
}

}// namespace Hpro

#endif  // __HPRO_STRUCTURE_HH

