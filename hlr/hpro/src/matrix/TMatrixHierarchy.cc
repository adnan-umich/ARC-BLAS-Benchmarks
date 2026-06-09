//
// Project     : HLIBpro
// File        : TMatrixHierarchy.cc
// Description : represents a level-wise hierarchy of matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <iostream>
#include <list>
#include <deque>
#include <algorithm>

#include "hpro/matrix/TBlockMatrix.hh"
#include "hpro/matrix/structure.hh"

#include "hpro/matrix/TMatrixHierarchy.hh"

namespace Hpro
{

//////////////////////////////////////////////////////////////
//
// local functions
//
//////////////////////////////////////////////////////////////

namespace
{

//
// comparison for block_list_t::sort
//
// bool
// rowis_compare ( TMatrix * A, TMatrix * B )
// {
//     return A->row_is().is_strictly_left_of( B->row_is() );
// }

// bool
// colis_compare ( TMatrix * A, TMatrix * B )
// {
//     return A->col_is().is_strictly_left_of( B->col_is() );
// }

//
// comparison for mat_storage_t::sort
//
template < typename value_t >
bool
list_compare ( typename TSparseBlockMatrix< value_t >::block_list_t *  list1,
               typename TSparseBlockMatrix< value_t >::block_list_t *  list2 )
{
    return list1->front()->row_is().is_strictly_left_of( list2->front()->row_is() );
}
               
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
// TSparseBlockMatrix
//
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//
// destructor
//
template < typename value_t >
TSparseBlockMatrix< value_t >::~TSparseBlockMatrix ()
{
    // NO DELETION, BECAUSE THIS MATRIX TYPE IS (FOR NOW) ONLY
    // USED TO HOLD POINTERS TO OTHERWISE REFERENCED MATRICES
    
    // // delete submatrices
    // for ( mat_storage_t::iterator  iter_row = _blocks.begin();
    //       iter_row != _blocks.end();
    //       ++iter_row )
    // {
    //     for ( block_list_t::iterator  iter_col = (*iter_row)->begin();
    //           iter_col != (*iter_row)->end();
    //           ++iter_col )
    //     {
    //         delete *iter_col;
    //     }
    // }

    //
    // delete block lists
    //
    
    for ( auto  iter = _block_rows.begin(); iter != _block_rows.end(); ++iter )
    {
        delete iter->second;
    }// for
    
    for ( auto  iter = _block_cols.begin(); iter != _block_cols.end(); ++iter )
    {
        delete iter->second;
    }// for
}

//
// access individual submatrix addressed by is0 × is1
//
template < typename value_t >
TMatrix< value_t > *
TSparseBlockMatrix< value_t >::block ( const TIndexSet &  is0,
                                       const TIndexSet &  is1 )
{
    block_list_t *  row = _block_rows[ is0 ];

    if ( row != nullptr )
    {
        for ( auto  M : * row )
        {
            if ( M->col_is() == is1 )
                return M;
        }// for
    }// if

    return nullptr;
}

//
// return submatrix t×s with is0 × is1 ⊆ t×s
template < typename value_t >
TMatrix< value_t > *
TSparseBlockMatrix< value_t >::block_containing ( const TIndexSet &  is0,
                                       const TIndexSet &  is1 )
{
    for ( auto  entry : _block_rows )
    {
        if ( is0.is_subset_of( entry.first ) )
        {
            for ( auto  M : * entry.second )
            {
                if ( is1.is_subset_of( M->col_is() ) )
                    return M;
            }// for
        }// if
    }// for
    
    return nullptr;
}
        

//
// insert matrix A into sparse block matrix
//
template < typename value_t >
void
TSparseBlockMatrix< value_t >::insert_block ( TMatrix< value_t > *  A )
{
    TIndexSet  is0( A->row_is() );
    TIndexSet  is1( A->col_is() );

    //
    // insert into block rows
    //
    
    auto  bl_row = _block_rows[ is0 ];

    if ( bl_row == nullptr )
    {
        //
        // insert new list into block storage
        //

        bl_row = new block_list_t;
        bl_row->push_front( A );
        _block_rows[ is0 ] = bl_row;

        // use row blocklist also in sparse block matrix storage
        _blocks.push_back( bl_row );
    }// if
    else
    {
        bl_row->push_back( A );
    }// else

    //
    // insert into block columns
    //
    
    auto  bl_col = _block_cols[ is1 ];

    if ( bl_col == nullptr )
    {
        bl_col = new block_list_t;
        bl_col->push_front( A );
        _block_cols[ is1 ] = bl_col;
    }// if
    else
    {
        bl_col->push_back( A );
    }// else
}

//
// sort data for efficient (and correct) access
//
template < typename value_t >
void
TSparseBlockMatrix< value_t >::sort ()
{
    _blocks.sort( list_compare< value_t > );
    
    // for ( auto & entry : _block_rows )
    // {
    //     entry.second->sort( rowis_compare );
    // }
    
    // for ( auto & entry : _block_cols )
    // {
    //     entry.second->sort( colis_compare );
    // }
    
    // // sort list of lists
    // std::sort( _blocks.begin(), _blocks.end(), list_compare );

    // for ( auto & entry : _block_rows )
    // {
    //     std::sort( entry.second->begin(), entry.second->end(), rowis_compare );
    // }

    // for ( auto & entry : _block_cols )
    // {
    //     std::sort( entry.second->begin(), entry.second->end(), colis_compare );
    // }
}

//
// return size in bytes used by this object
//
template < typename value_t >
size_t
TSparseBlockMatrix< value_t >::byte_size () const
{
    size_t  size = 0;
    
    for ( auto  iter_row = _blocks.cbegin(); iter_row != _blocks.cend(); ++iter_row )
    {
        for ( auto  iter_col = (*iter_row)->cbegin(); iter_col != (*iter_row)->cend(); ++iter_col )
        {
            size += (*iter_col)->byte_size();
        }// for
    }// for

    return size;
}

//
// print content
//
template < typename value_t >
void
TSparseBlockMatrix< value_t >::print ( const uint  ofs ) const
{
    uint  row = 0;
    
    for ( auto  iter_row = _blocks.cbegin(); iter_row != _blocks.cend(); ++iter_row )
    {
        for ( uint  i = 0; i < ofs; ++i )
            std::cout << ' ';

        std::cout << "block row " << row << std::endl;
        
        for ( auto  iter_col = (*iter_row)->cbegin(); iter_col != (*iter_row)->cend(); ++iter_col )
        {
            for ( uint  i = 0; i < ofs + 2; ++i )
                std::cout << ' ';
            
            std::cout << (*iter_col)->typestr() << ' '
                      << (*iter_col)->block_is() << std::endl;
        }// for

        ++row;
    }// for
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
// TMatrixHierarchy
//
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//
// ctors and dtor
//

//
// construct an empty hierarchy
//
template < typename value_t >
TMatrixHierarchy< value_t >::TMatrixHierarchy ()
{
}

//
// construct a hierarchy based on given block matrix
//
template < typename value_t >
TMatrixHierarchy< value_t >::TMatrixHierarchy ( TMatrix< value_t > *  A,
                                                const bool            with_blocked )
{
    //
    // for now, we need as a 
    //
    
    //
    // perform BFS through blockmatrix and create sparse block matrices
    // for each level (with at least one leaf)
    //

    size_t                                        lvl = 0;            // counter for current level
    std::deque< TMatrix< value_t > * >            matrices;           // matrices on current level
    std::list< TSparseBlockMatrix< value_t > * >  hierarchy_list;     // temporary storage for matrices in hierarchy

    matrices.push_front( A );
    
    while ( ! matrices.empty() )
    {
        //
        // test, if current level holds a leaf matrix
        //
        
        bool  found_leaf = false;
        
        for ( auto  M : matrices )
        {
            if ( ! M->is_blocked() )
            {
                found_leaf = true;
                break;
            }// if
        }// for

        // 
        // if leaf is in current level set or if also inner nodes are
        // requested, build sparse block matrix
        //

        if ( found_leaf || with_blocked )
        {
            auto  SB = std::make_unique< TSparseBlockMatrix< value_t > >();

            for ( auto  M : matrices )
            {
                SB->insert_block( M );
            }// for

            SB->sort();
            _hierarchy.push_back( SB.release() );
        }// if
        else
            _hierarchy.push_back( ptrcast( nullptr, TSparseBlockMatrix< value_t > ) );

        //
        // collect list of son matrices of current level
        //

        std::deque< TMatrix< value_t > * >  sons;
        
        for ( auto  M : matrices )
        {
            if ( is_blocked( M ) )
            {
                auto  BM = ptrcast( M, TBlockMatrix< value_t > );
                                             
                for ( uint  bi = 0; bi < BM->block_rows(); bi++ )
                {
                    for ( uint  bj = 0; bj < BM->block_cols(); bj++ )
                    {
                        if ( BM->block( bi, bj ) != nullptr )
                            sons.push_back( BM->block( bi, bj ) );
                    }// for
                }// for
            }// if
        }// for

        matrices = std::move( sons );
        lvl++;
    }// while
}

//
// destructor
//
template < typename value_t >
TMatrixHierarchy< value_t >::~TMatrixHierarchy ()
{
    for ( auto *  SB : _hierarchy )
    {
        if ( SB != nullptr )
            delete SB;
    }// for
}

/////////////////////////////////////////////////
//
// misc.
//

//
// return size in bytes used by this object
//
template < typename value_t >
size_t
TMatrixHierarchy< value_t >::byte_size () const
{
    size_t  size = 0;

    for ( size_t  lvl = 0; lvl < _hierarchy.size(); lvl++ )
    {
        if ( _hierarchy[ lvl ] != nullptr )
            size += _hierarchy[ lvl ]->byte_size();
    }// for

    return size;
}

//
// print content
//
template < typename value_t >
void
TMatrixHierarchy< value_t >::print ( const uint  ofs ) const
{
    for ( size_t  lvl = 0; lvl < _hierarchy.size(); lvl++ )
    {
        if ( _hierarchy[ lvl ] != nullptr )
        {
            for ( uint  i = 0; i < ofs; ++i )
                std::cout << ' ';

            std::cout << "Level " << lvl << std::endl;

            _hierarchy[ lvl ]->print( ofs + 2 );
        }// if
    }// for
}

template class TMatrixHierarchy< float >;
template class TMatrixHierarchy< double >;
template class TMatrixHierarchy< std::complex< float > >;
template class TMatrixHierarchy< std::complex< double > >;
    
}// namespace
