//
// Project     : HLIBpro
// File        : structure.cc
// Description : test functions for various matrix structures
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/error.hh"

#include "hpro/matrix/TBlockMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"

#include "hpro/parallel/NET.hh"

#include "hpro/matrix/structure.hh"

namespace Hpro
{

namespace
{

//
// return true if <A> is block diagonal
//
template < typename value_t >
bool
is_diag_local ( const TMatrix< value_t > * A )
{
    if ( A == nullptr )
        return false;
    
    if ( ! is_blocked( A ) )
        return false;

    //
    // check if only diagonal blocks and last row/column
    // are non-empty
    //
    
    const auto  BA  = cptrcast( A, TBlockMatrix< value_t > );
    const uint  br  = BA->block_rows();
    const uint  bc  = BA->block_cols();

    // must have at least 2 diagonal blocks
    if ( br == 1 )
    {
        return false;
    }// if
    else
    {
        for ( uint i = 0; i < br; i++ )
        {
            for ( uint j = 0; j < bc; j++ )
            {
                if ( i == j )
                    continue;
            
                const auto  A_ij = BA->block( i, j );

                // if matrix is not existent or is low-rank with 0 rank, continue
                if (( A_ij == nullptr ) || A_ij->is_zero() )
                    continue;

                // matrix exists and has content, so no DD format
                return false;
            }// for
        }// for
    }// else

    return true;
}

}// namespace anonymous

template < typename value_t >
bool
is_diag ( const TMatrix< value_t > *  A )
{
    bool  res = is_diag_local( A );
    
    if ( A->procs().size() <= 1 )
        return res;

    //
    // exchange local result to form global result
    //

    // size_t  cres = ( res ? 1 : 0 );
    
    // cres = NET::reduce_all( A->procs(), cres, NET::OP_SUM );

    // if ( uint( cres ) != A->procs().size() )
    //     res = false;
    
    return res;
}

namespace
{

//
// return true if <A> has domain-decomposition format, e.g.
// only diagonal blocks and non-zero last row/column
//
template < typename value_t >
bool
is_dd_local ( const TMatrix< value_t > * A )
{
    if ( A == nullptr )
        return false;
    
    if ( ! is_blocked( A ) )
        return false;

    //
    // check if only diagonal blocks and last row/column
    // are non-empty
    //
    
    const auto  BA  = cptrcast( A, TBlockMatrix< value_t > );
    const uint  br  = BA->block_rows();
    const uint  bc  = BA->block_cols();

    // check for at least two subdomains
    if (( br < 3 ) || ( bc < 3 ))
    {
        return false;
    }// if
    else
    {
        for ( uint i = 0; i < br; i++ )
        {
            //
            // test upper right part
            //

            for ( uint j = i+1; j < bc-1; j++ )
            {
                const auto  A_ij = BA->block( i, j );

                // if matrix is not existent or is low-rank with 0 rank, continue
                if (( A_ij == nullptr ) || A_ij->is_zero() )
                    continue;

                // matrix exists and has content, so no DD format
                return false;
            }// for

            //
            // test lower left part
            //

            for ( uint j = i+1; j < br-1; j++ )
            {
                const auto  A_ij = BA->block( j, i );

                // if matrix is not existent or is low-rank with 0 rank, continue
                if (( A_ij == nullptr ) || A_ij->is_zero() )
                    // ( is_lowrank( A_ij ) && ( cptrcast( A_ij, TRkMatrix< value_t > )->rank() == 0 )))
                    continue;

                // matrix exists and has content, so no DD format
                return false;
            }// for
        }// for
    }// else

    return true;
}

}// namespace anonymous

template < typename value_t >
bool
is_dd ( const TMatrix< value_t > *  A )
{
    bool  res = is_dd_local( A );
    
    if ( A->procs().size() <= 1 )
        return res;

    //
    // exchange local result to form global result
    //

    // size_t  cres = ( res ? 1 : 0 );
    
    // cres = NET::reduce_all( A->procs(), cres, NET::OP_SUM );

    // if ( uint( cres ) != A->procs().size() )
    //     res = false;
    
    return res;
}

//
// return true, if A is block matrix with flat hierarchy
//
template < typename value_t >
bool
is_flat ( const TMatrix< value_t > *  A )
{
    if ( A == nullptr )
        return false;
    
    if ( ! is_blocked( A ) )
        return false;

    if ( A->parent() != nullptr )
        return false;
    
    auto  B = cptrcast( A, TBlockMatrix< value_t > );

    // // minimal number of blocks are assumed
    // if ( std::min( B->block_rows(), B->block_cols() ) > 2 )
    //     return false;
    
    for ( uint i = 0; i < B->block_rows(); ++i )
        for ( uint j = 0; j < B->block_cols(); ++j )
            if ( is_blocked( B->block( i, j ) ) )
                return false;

    return true;
}

//
// return true if <A> is lower left block matrix
//
template < typename value_t >
bool
is_lower_left ( const TMatrix< value_t > * A )
{
    if ( A == nullptr )
        return false;
    
    if ( ! is_blocked( A ) )
        return false;

    //
    // check if all blocks in upper right part are empty
    //
    
    const auto  BA = cptrcast( A, TBlockMatrix< value_t > );
    const uint  br = BA->block_rows();
    const uint  bc = BA->block_cols();

    for ( uint i = 0; i < br; i++ )
    {
        //
        // test upper right part
        //

        for ( uint j = i+1; j < bc; j++ )
        {
            const auto  A_ij = BA->block( i, j );

            // if matrix is not existent or is low-rank with 0 rank, continue
            if (( A_ij == nullptr ) || A_ij->is_zero() )
                continue;

            // matrix exists and has content, so not lower triangular
            return false;
        }// for
    }// for

    return true;
}

//
// return true if <A> is upper right block matrix
//
template < typename value_t >
bool
is_upper_right ( const TMatrix< value_t > * A )
{
    if ( A == nullptr )
        return false;
    
    if ( ! is_blocked( A ) )
        return false;

    //
    // check if all blocks in upper right part are empty
    //
    
    const auto  BA = cptrcast( A, TBlockMatrix< value_t > );
    const uint  br = BA->block_rows();
    const uint  bc = BA->block_cols();

    for ( uint i = 0; i < bc; i++ )
    {
        //
        // test lower left part
        //

        for ( uint j = i+1; j < br; j++ )
        {
            const auto  A_ij = BA->block( j, i );

            // if matrix is not existent or is low-rank with 0 rank, continue
            if (( A_ij == nullptr ) || A_ij->is_zero() )
                continue;

            // matrix exists and has content, so not upper triangular
            return false;
        }// for
    }// for

    return true;
}

//!
//! return true if matrix is zero, e.g. all entries zero
//!
template < typename value_t >
bool
is_zero ( const TMatrix< value_t > *  A )
{
    // NULL matrices are considered zero!
    if ( A == nullptr )
        return true;
    else
        return A->is_zero();
}

//
// return maximal ID in all submatrices
//
template < typename value_t >
int
max_id ( const TMatrix< value_t > *  A )
{
    int  id = A->id();
        
    // recurse
    if ( A->is_blocked() )
    {
        const auto  B = cptrcast( A, TBlockMatrix< value_t > );
        
        for ( uint j = 0; j < B->block_cols(); ++j )
        {
            for ( uint i = 0; i < B->block_rows(); ++i )
            {
                const auto  B_ij = B->block( i, j );
                
                if ( B_ij == nullptr )
                    continue;

                id = std::max( id, max_id( B_ij ) );
            }// for
        }// for
    }// if

    return id;
}

//////////////////////////////////////////////////////////
//
// return vector containing diagonal coefficients of A
//

namespace
{

template < typename value_t >
void
diagonal ( const TMatrix< value_t > *  A,
           TScalarVector< value_t > *  diag )
{
    if ( A == nullptr )
        HERROR( ERR_ARG, "diagonal", "matrix is nullptr" );

    if ( A->row_is() != A->col_is() )
        return;
    
    if ( is_blocked( A ) )
    {
        auto  B = cptrcast( A, TBlockMatrix< value_t > );
        
        for ( uint  i = 0; i < std::min( B->block_rows(), B->block_cols() ); ++i )
        {
            diagonal( B->block( i, i ), diag );
        }// for
    }// if
    else
    {
        const idx_t  ofs = A->row_is().first();

        for ( auto  i : A->row_is() )
            diag->set_entry( i - diag->ofs(), A->entry( i - ofs, i - ofs ) );
    }// else
}

}// namespace anonymous

template < typename value_t >
std::unique_ptr< TVector< value_t > >
diagonal ( const TMatrix< value_t > *  A )
{
    if ( A == nullptr )
        HERROR( ERR_ARG, "diagonal", "matrix is nullptr" );
    
    auto  diag = std::make_unique< TScalarVector< value_t > >( is( 0, std::min< idx_t >( A->rows(), A->cols() )-1 ) );

    diagonal( A, diag.get() );

    return std::unique_ptr< TVector< value_t > >( diag.release() );
}

//////////////////////////////////////////////////////////
//
// return copy of all diagonal blocks of A with full hierarchy
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
copy_diag   ( const TMatrix< value_t > * A )
{
    if ( A == nullptr )
        return nullptr;
    
    if ( is_blocked( A ) )
    {
        auto  BA = cptrcast( A, TBlockMatrix< value_t > );
        auto  T  = std::make_unique< TBlockMatrix< value_t > >();

        T->copy_struct_from( A );
        T->set_block_struct( BA->block_rows(), BA->block_cols() );

        for ( uint i = 0; i < std::min( BA->block_rows(), BA->block_cols() ); i++ )
            T->set_block( i, i, copy_diag( BA->block(i,i) ).release() );

        return T;
    }// if
    else
        return A->copy();
}

//
// return first level on diagonal on which leaf blocks are found
//
namespace
{

template < typename value_t >
uint
get_first_diag_leaf_level ( const TMatrix< value_t > *  A,
                            const uint                  lvl )
{
    uint  max_lvl = lvl;

    if ( is_leaf( A ) || is_small( A ) )
        return lvl;
    
    if ( ! is_leaf( A ) && ! is_small( A ) && is_blocked( A ) )
    {
        const auto  B = cptrcast( A, TBlockMatrix< value_t > );

        max_lvl = 0;

        for ( uint  i = 0; i < B->block_rows(); ++i )
        {
            if ( B->block( i, i ) != nullptr )
            {
                auto  lvl_ii = get_first_diag_leaf_level( B->block( i, i ), lvl+1 );
                    
                if ( max_lvl == 0 ) max_lvl = lvl_ii;
                else                max_lvl = std::min( max_lvl, lvl_ii );
            }// if
        }// for
    }// if

    return max_lvl;
}

}// namespace anonymous

template < typename value_t >
uint
get_first_diag_leaf_level ( const TMatrix< value_t > *  A )
{
    return get_first_diag_leaf_level( A, 0 );
}

//////////////////////////////////////////////////////////
//
// return matrix with rowis( A_row ) × colis( A_col )
//
template < typename value_t >
TMatrix< value_t > *
product_block ( const TMatrix< value_t > *  A_row,
                const TMatrix< value_t > *  A_col )
{
    //
    // iterate over both lists in case of symmetric matrices (upper half missing)
    //
    
    for ( auto  M_row = A_row->next_in_block_row(); M_row != nullptr; M_row = M_row->next_in_block_row() )
    {
        if ( M_row->col_is() == A_col->col_is() )
            return M_row;
    }// for

    for ( auto  M_col = A_col->next_in_block_col(); M_col != nullptr; M_col = M_col->next_in_block_col() )
    {
        if ( M_col->row_is() == A_row->row_is() )
            return M_col;
    }// for

    return nullptr;
}

template < typename value_t >
TMatrix< value_t > *
product_block ( const matop_t               op_row,
                const TMatrix< value_t > *  A_row,
                const matop_t               op_col,
                const TMatrix< value_t > *  A_col )
{
    //
    // iterate over both lists in case of symmetric matrices (upper half missing)
    //

    for ( auto  M_row = A_row->next_in_block_row( op_row ); M_row != nullptr; M_row = M_row->next_in_block_row( op_row ) )
    {
        if ( M_row->col_is( op_row ) == A_col->col_is( op_col ) )
            return M_row;
    }// for

    for ( auto  M_col = A_col->next_in_block_col( op_col ); M_col != nullptr; M_col = M_col->next_in_block_col( op_col ) )
    {
        if ( M_col->row_is( op_col ) == A_row->row_is( op_row ) )
            return M_col;
    }// for

    return nullptr;
}

//
// return sub block of M with given id
//
template < typename value_t >
const TMatrix< value_t > *
get_block ( const TMatrix< value_t > *  M,
            const int                   id )
{
    if ( M == nullptr )
        return nullptr;
    
    if ( M->id() == id )
        return M;

    if ( is_blocked( M ) )
    {
        auto  B = cptrcast( M, TBlockMatrix< value_t > );

        for ( uint  i = 0; i < B->block_rows(); ++i )
            for ( uint  j = 0; j < B->block_cols(); ++j )
            {
                auto  M_ij = get_block( B->block( i, j ), id );

                if ( M_ij != nullptr )
                    return M_ij;
            }// for
    }// if

    return nullptr;
}

//
// return matrix block in A corresponding to block index set bis
//
template < typename value_t >
const TMatrix< value_t > *
get_block ( const TMatrix< value_t > *  A,
            const TBlockIndexSet &      bis )
{
    if ( A == nullptr )
        return nullptr;
    
    if ( bis.is_subset_of( A->block_is() ) )
    {
        if ( bis == A->block_is() )
            return A;

        if ( is_blocked( A ) )
        {
            const auto  B = cptrcast( A, TBlockMatrix< value_t > );

            for ( uint i = 0; i < B->block_rows(); ++i )
            {
                for ( uint j = 0; j < B->block_rows(); ++j )
                {
                    auto  B_ij = B->block( i, j );
                    
                    if (( B_ij != nullptr ) && bis.is_subset_of( B_ij->block_is() ) )
                    {
                        return get_block( B_ij, bis );
                    }// if
                }// for
            }// for
        }// if

        return nullptr;
    }// if
    else
    {
        // try to look in parent matrix
        return get_block( A->parent(), bis );
    }// else
}

//!
//! return number of sub blocks in matrix (including inner blocks)
//!
template < typename value_t >
size_t
get_nblocks ( const TMatrix< value_t > *   M )
{
    if ( M == nullptr )
        return 0;
    
    if ( is_blocked( M ) )
    {
        auto    B = cptrcast( M, TBlockMatrix< value_t > );
        size_t  n = 1;  // this block

        for ( uint  i = 0; i < B->nblock_rows(); ++i )
            for ( uint  j = 0; j < B->nblock_cols(); ++j )
                n += get_nblocks( B->block( i, j ) );

        return n;
    }// if
    else
        return 1;
}

//!
//! return number of leaf blocks in matrix
//!
template < typename value_t >
size_t
get_nblocks_leaf ( const TMatrix< value_t > *   M )
{
    if ( M == nullptr )
        return 0;
    
    if ( is_blocked( M ) )
    {
        auto    B = cptrcast( M, TBlockMatrix< value_t > );
        size_t  n = 0;

        for ( uint  i = 0; i < B->nblock_rows(); ++i )
            for ( uint  j = 0; j < B->nblock_cols(); ++j )
                n += get_nblocks_leaf( B->block( i, j ) );

        return n;
    }// if
    else
        return 1;
}   

//!
//! return number of low-rank sub blocks in matrix
//!
template < typename value_t >
size_t
get_nblocks_lowrank ( const TMatrix< value_t > *   M )
{
    if ( M == nullptr )
        return 0;
    
    if ( is_blocked( M ) )
    {
        auto    B = cptrcast( M, TBlockMatrix< value_t > );
        size_t  n = 0;

        for ( uint  i = 0; i < B->nblock_rows(); ++i )
            for ( uint  j = 0; j < B->nblock_cols(); ++j )
                n += get_nblocks_lowrank( B->block( i, j ) );

        return n;
    }// if
    else if ( is_lowrank( M ) )
        return 1;
    else
        return 0;
}   

//!
//! return number of dense sub blocks in matrix
//!
template < typename value_t >
size_t
get_nblocks_dense  ( const TMatrix< value_t > *   M )
{
    if ( M == nullptr )
        return 0;
    
    if ( is_blocked( M ) )
    {
        auto    B = cptrcast( M, TBlockMatrix< value_t > );
        size_t  n = 0;

        for ( uint  i = 0; i < B->nblock_rows(); ++i )
            for ( uint  j = 0; j < B->nblock_cols(); ++j )
                n += get_nblocks_dense( B->block( i, j ) );

        return n;
    }// if
    else if ( is_dense( M ) )
        return 1;
    else
        return 0;
}   

//!
//! return number of level in matrix
//!
template < typename value_t >
size_t
get_nlevel ( const TMatrix< value_t > &  M )
{
    if ( is_blocked( M ) )
    {
        auto    B    = cptrcast( &M, TBlockMatrix< value_t > );
        size_t  nlvl = 0;

        for ( uint  i = 0; i < B->nblock_rows(); ++i )
            for ( uint  j = 0; j < B->nblock_cols(); ++j )
            {
                if ( B->block( i, j ) != nullptr )
                    nlvl = std::max( nlvl, get_nlevel( *(B->block( i, j )) ) );
            }// for

        return nlvl + 1; // plus this level
    }// if
    else
        return 1;
}

//////////////////////////////////////////////////////////
//
// extract real/imaginary parts
//

//
// return real part of matrix \a M
//
template < typename value_t >
std::unique_ptr< TMatrix< real_type_t< value_t > > >
restrict_re ( const TMatrix< value_t > *  M,
              const TTruncAcc &           acc )
{
    using  real_t = real_type_t< value_t >;
    
    if ( M == nullptr )
        HERROR( ERR_ARG, "restrict_im", "matrix is null" );

    // return copy in case of real valued matrix
    if constexpr ( ! is_complex_type< value_t >::value )
        return M->copy();
    
    if ( is_hmat( M ) )
    {
        auto  HM = cptrcast( M, THMatrix< value_t > );
        auto  C  = std::make_unique< THMatrix< real_t > >();
        auto  HC = ptrcast( C.get(), THMatrix< real_t > );

        HC->copy_struct_from_all( HM );
        HC->comp_min_max_idx();
        HC->set_row_perm( HM->row_perm_e2i(), HM->row_perm_i2e() );
        HC->set_col_perm( HM->col_perm_e2i(), HM->col_perm_i2e() );
        
        for ( uint  i = 0; i < HC->nblock_rows(); ++i )
        {
            for ( uint  j = 0; j < HC->nblock_cols(); ++j )
            {
                if ( HM->block( i, j ) != nullptr )
                {
                    auto  H_ij = restrict_re( HM->block( i, j ), acc );
                    
                    H_ij->set_parent( HC );
                    HC->set_block( i, j, H_ij.release() );
                }// if
            }// for
        }// for

        C->set_id( M->id() );
        
        return C;
    }// if
    else if ( is_blocked( M ) )
    {
        auto  BM = cptrcast( M, TBlockMatrix< value_t > );
        auto  C  = std::make_unique< TBlockMatrix< real_t > >();
        auto  BC = ptrcast( C.get(), TBlockMatrix< real_t > );

        BC->copy_struct_from_all( BM );
        
        for ( uint  i = 0; i < BC->nblock_rows(); ++i )
        {
            for ( uint  j = 0; j < BC->nblock_cols(); ++j )
            {
                if ( BM->block( i, j ) != nullptr )
                {
                    auto  B_ij = restrict_re( BM->block( i, j ), acc );
                    
                    B_ij->set_parent( BC );
                    BC->set_block( i, j, B_ij.release() );
                }// if
            }// for
        }// for

        C->set_id( M->id() );
        
        return C;
    }// if
    else if ( is_dense( M ) )
    {
        //
        // restrict to imaginary part
        //
        
        auto  DM = cptrcast( M, TDenseMatrix< value_t > )->blas_mat();
        auto  C  = std::make_unique< TDenseMatrix< real_t > >( M->row_is(), M->col_is() );
        auto  DC = C->blas_mat();
        
        for ( size_t  j = 0; j < DC.ncols(); ++j )
            for ( size_t  i = 0; i < DC.nrows(); ++i )
                DC(i,j) = std::real( DM(i,j) );
        
        C->set_id( M->id() );
        C->set_form( M->form() );
        C->set_procs( M->procs() );
        
        return C;
    }// if
    else if ( is_lowrank( M ) )
    {
        //
        // approximate restriction to real part
        //
            
        auto  RM   = cptrcast( M, TRkMatrix< value_t > );
        auto  A    = RM->blas_mat_A();
        auto  B    = RM->blas_mat_B();
        auto  lacc = acc( RM->row_is(), RM->col_is() );

        //
        // U × V' = (U_r + i·U_i) × (V_r + i·V_i)'
        //        = (U_r + i·U_i) × (V_r' - i·V_i')
        //        = (U_r × V_r' + U_i × V_i' ) + i·( U_i × V_r' - U_r × V_i' )
        //        =       [U_r,  U_i] × [V_r, V_i]'
        //          + i · [U_i, -U_r] × [V_r, V_i]'
        //

        const auto  nrows = RM->nrows();
        const auto  ncols = RM->ncols();
        const auto  k     = RM->rank();
        auto        Ure   = BLAS::Matrix< real_t >( nrows, 2*k );
        auto        Vre   = BLAS::Matrix< real_t >( ncols, 2*k );

        // copy [U_r, U_i] to Ure
        for ( uint  l = 0; l < k; ++l )
        {
            for ( uint  i = 0; i < nrows; ++i )
            {
                Ure(i,l)   = std::real( A(i,l) );
                Ure(i,k+l) = std::imag( A(i,l) );
            }// for
        }// for

        // copy [V_r, V_i] to Vre
        for ( uint  l = 0; l < k; ++l )
        {
            for ( uint  i = 0; i < ncols; ++i )
            {
                Vre(i,l)   = std::real( B(i,l) );
                Vre(i,k+l) = std::imag( B(i,l) );
            }// for
        }// for

        // reduce rank
        BLAS::truncate_svd( Ure, Vre, lacc );
            
        auto  R = std::make_unique< TRkMatrix< real_t > >( M->row_is(), M->col_is(), std::move( Ure ), std::move( Vre ) );

        R->set_id( M->id() );
        R->set_form( M->form() );
        R->set_procs( M->procs() );

        return R;
    }// if
    else
        HERROR( ERR_MAT_TYPE, "restrict_re", "unsupported matrix type " + M->typestr() );
}

//
// return imaginary part of matrix \a M
//
template < typename value_t >
std::unique_ptr< TMatrix< real_type_t< value_t > > >
restrict_im ( const TMatrix< value_t > *  M,
              const TTruncAcc &           acc )
{
    using  real_t = real_type_t< value_t >;
    
    constexpr bool  real_valued = ! is_complex_type< value_t >::value;
    
    if ( M == nullptr )
        HERROR( ERR_ARG, "restrict_im", "matrix is null" );
    
    if ( is_hmat( M ) )
    {
        auto  HM = cptrcast( M, THMatrix< value_t > );
        auto  C  = std::make_unique< THMatrix< real_t > >();
        auto  HC = ptrcast( C.get(), THMatrix< real_t > );

        HC->copy_struct_from_all( HM );
        HC->comp_min_max_idx();
        HC->set_row_perm( HM->row_perm_e2i(), HM->row_perm_i2e() );
        HC->set_col_perm( HM->col_perm_e2i(), HM->col_perm_i2e() );
        
        for ( uint  i = 0; i < HC->nblock_rows(); ++i )
        {
            for ( uint  j = 0; j < HC->nblock_cols(); ++j )
            {
                if ( HM->block( i, j ) != nullptr )
                {
                    auto  H_ij = restrict_im( HM->block( i, j ), acc );
                    
                    H_ij->set_parent( HC );
                    HC->set_block( i, j, H_ij.release() );
                }// if
            }// for
        }// for

        C->set_id( M->id() );
        
        return C;
    }// if
    else if ( is_blocked( M ) )
    {
        auto  BM = cptrcast( M, TBlockMatrix< value_t > );
        auto  C  = std::make_unique< TBlockMatrix< real_t > >();
        auto  BC = ptrcast( C.get(), TBlockMatrix< real_t > );

        BC->copy_struct_from_all( BM );
        
        for ( uint  i = 0; i < BC->nblock_rows(); ++i )
        {
            for ( uint  j = 0; j < BC->nblock_cols(); ++j )
            {
                if ( BM->block( i, j ) != nullptr )
                {
                    auto  B_ij = restrict_im( BM->block( i, j ), acc );
                    
                    B_ij->set_parent( BC );
                    BC->set_block( i, j, B_ij.release() );
                }// if
            }// for
        }// for

        C->set_id( M->id() );
        
        return C;
    }// if
    else if ( is_dense( M ) )
    {
        //
        // restrict to imaginary part
        //
            
        auto  C  = std::make_unique< TDenseMatrix< real_t > >( M->row_is(), M->col_is() );

        if ( ! real_valued )
        {
            auto  DM = cptrcast( M, TDenseMatrix< value_t > )->blas_mat();
            auto  DC = C->blas_mat();
            
            for ( size_t  j = 0; j < DC.ncols(); ++j )
                for ( size_t  i = 0; i < DC.nrows(); ++i )
                    DC(i,j) = std::imag( DM(i,j) );
        }// if
        
        C->set_id( M->id() );
        C->set_form( M->form() );
        C->set_procs( M->procs() );
        
        return C;
    }// if
    else if ( is_lowrank( M ) )
    {
        //
        // return zero matrix in case of real value
        //
        
        if ( real_valued )
        {
            auto  R = std::make_unique< TRkMatrix< real_t > >( M->row_is(), M->col_is() );
            
            R->set_id( M->id() );
            R->set_form( M->form() );
            R->set_procs( M->procs() );

            return R;
        }// if
        
        //
        // approximate restriction to imaginary part
        //
            
        auto  RM   = cptrcast( M, TRkMatrix< value_t > );
        auto  A    = RM->blas_mat_A();
        auto  B    = RM->blas_mat_B();
        auto  lacc = acc( RM->row_is(), RM->col_is() );

        //
        // U × V' = (U_r + i·U_i) × (V_r + i·V_i)'
        //        = (U_r + i·U_i) × (V_r' - i·V_i')
        //        = (U_r × V_r' + U_i × V_i' ) + i·( U_i × V_r' - U_r × V_i' )
        //        =       [U_r,  U_i] × [V_r, V_i]'
        //          + i · [U_i, -U_r] × [V_r, V_i]'
        //

        const auto  nrows = RM->nrows();
        const auto  ncols = RM->ncols();
        const auto  k     = RM->rank();
        auto        Uim   = BLAS::Matrix< real_t >( nrows, 2*k );
        auto        Vim   = BLAS::Matrix< real_t >( ncols, 2*k );

        // copy [U_i, -U_r] to Uim
        for ( uint  l = 0; l < k; ++l )
        {
            for ( uint  i = 0; i < nrows; ++i )
            {
                Uim(i,l)   =  std::imag( A(i,l) );
                Uim(i,k+l) = -std::real( A(i,l) );
            }// for
        }// for

        // copy [V_r, V_i] to Vim
        for ( uint  l = 0; l < k; ++l )
        {
            for ( uint  i = 0; i < ncols; ++i )
            {
                Vim(i,l)   = std::real( B(i,l) );
                Vim(i,k+l) = std::imag( B(i,l) );
            }// for
        }// for

        // reduce rank
        BLAS::truncate_svd( Uim, Vim, lacc );
            
        auto  R = std::make_unique< TRkMatrix< real_t > >( M->row_is(), M->col_is(), std::move( Uim ), std::move( Vim ) );

        R->set_id( M->id() );
        R->set_form( M->form() );
        R->set_procs( M->procs() );

        return R;
    }// if
    else
        HERROR( ERR_MAT_TYPE, "restrict_im", "unsupported matrix type " + M->typestr() );
    
}

//
// explicit template instantiations
//

#define INST_ALL( type )                                                \
    template bool is_diag ( const TMatrix< type > * );                  \
    template bool is_dd ( const TMatrix< type > * );                    \
    template bool is_flat ( const TMatrix< type > * );                  \
    template bool is_lower_left ( const TMatrix< type > * A );          \
    template bool is_upper_right ( const TMatrix< type > * A );         \
    template bool is_zero ( const TMatrix< type > * );                  \
    template int max_id ( const TMatrix< type > * );                    \
    template std::unique_ptr< TVector< type > > diagonal ( const TMatrix< type > * ); \
    template std::unique_ptr< TMatrix< type > > copy_diag   ( const TMatrix< type > * A ); \
    template uint get_first_diag_leaf_level ( const TMatrix< type > * ); \
    template TMatrix< type > * product_block ( const TMatrix< type > *, const TMatrix< type > * ); \
    template TMatrix< type > * product_block ( const matop_t, const TMatrix< type > *, const matop_t, const TMatrix< type > * ); \
    template const TMatrix< type > * get_block ( const TMatrix< type > *, const int ); \
    template const TMatrix< type > * get_block ( const TMatrix< type > *, const TBlockIndexSet & ); \
    template size_t get_nblocks ( const TMatrix< type > * );            \
    template size_t get_nblocks_leaf ( const TMatrix< type > * );       \
    template size_t get_nblocks_lowrank ( const TMatrix< type > * );    \
    template size_t get_nblocks_dense  ( const TMatrix< type > * );     \
    template size_t get_nlevel ( const TMatrix< type > & );             \
    template std::unique_ptr< TMatrix< real_type_t< type > > > restrict_re ( const TMatrix< type > *, const TTruncAcc & ); \
    template std::unique_ptr< TMatrix< real_type_t< type > > > restrict_im ( const TMatrix< type > *, const TTruncAcc & );

INST_ALL( float )
INST_ALL( double )
INST_ALL( std::complex< float > )
INST_ALL( std::complex< double > )

}// namespace Hpro
