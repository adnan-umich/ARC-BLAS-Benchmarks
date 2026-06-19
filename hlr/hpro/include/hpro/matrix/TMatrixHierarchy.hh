#ifndef __HPRO_TMATRIXHIERARCHY_HH
#define __HPRO_TMATRIXHIERARCHY_HH
//
// Project     : HLIBpro
// File        : TMatrixHierarchy.hh
// Description : represents a level-wise hierarchy of matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>
#include <list>
#include <deque>
#include <map>

#include "hpro/matrix/TMatrix.hh"

namespace Hpro
{

//!
//! \class TSparseBlockMatrix
//! \brief Represents a n×m block matrix with only a small number of non-null
//!        sub matrices stored in an efficient way
//!
template < typename value_t >
class TSparseBlockMatrix
{
public:
    //! \typedef block_list_t
    //! list for block rows and columns
    using  block_list_t  = std::list< TMatrix< value_t > * >;
        
    //! \typedef mat_storage_t
    //! storage type for sparse block matrix
    using  mat_storage_t = std::list< block_list_t * >;
        
    //! \typedef list_map_t
    //! mapping of indexsets to block lists
    using  list_map_t    = std::map< TIndexSet, block_list_t *, TIndexSet::map_cmp_t >;
        
private:
    //! sparse matrix of submatrices
    mat_storage_t  _blocks;
        
    //! mapping of indexsets to block rows and block columns
    list_map_t     _block_rows, _block_cols;
        
public:
    //////////////////////////////////////////////////////////////
    //
    // ctors and dtor
    //
        
    //! construct block matrix without any submatrices
    TSparseBlockMatrix () {}
        
    //! destruct sparse block matrix
    ~TSparseBlockMatrix ();

    // sort data for efficient (and correct) access
    void  sort ();
        
    //////////////////////////////////////////////////////////////
    //
    // submatrix management functions
    //
        
    //! access individual submatrix addressed by is0 × is1
    TMatrix< value_t > *
    block ( const TIndexSet &  is0,
            const TIndexSet &  is1 );
        
    //! return submatrix t×s with is0 × is1 ⊆ t×s
    TMatrix< value_t > *
    block_containing ( const TIndexSet &  is0,
                       const TIndexSet &  is1 );
        
    //! insert matrix A into sparse block matrix
    void
    insert_block ( TMatrix< value_t > *  A );
        
    //! return sparse block matrix
    const mat_storage_t *
    blocks () const
    {
        return & _blocks;
    }
        
    //! return block row corresponding to indexset \a is
    block_list_t *
    block_row ( const TIndexSet &  is )
    {
        return _block_rows[ is ];
    }
        
    //! return block column corresponding to indexset \a is
    block_list_t *
    block_col ( const TIndexSet &  is )
    {
        return _block_cols[ is ];
    }

    /////////////////////////////////////////////////
    //
    // misc.
    //
        
    //! return size in bytes used by this object
    size_t byte_size () const;

    //! print content of block matrix
    void   print ( const uint  ofs = 0 ) const;
};

//!
//! \class   TMatrixHierarchy
//! \brief   Represents a level-wise hierarchy of matrices.
//! \details Stores a hierarchy of block matrices, such that each level matrix
//!          holds all matrices (or leaves) of that level of a given matrix.
//!
template < typename value_t >
class TMatrixHierarchy
{
private:
    //! the matrix hierarchy
    std::deque< TSparseBlockMatrix< value_t > * >  _hierarchy;

public:
    //////////////////////////////////////////////////////////////
    //
    // ctors and dtor
    //

    //! construct an empty hierarchy
    TMatrixHierarchy ();
    
    //! construct a hierarchy based on given matrix A
    //! \param  A             matrix to construct hierarchy with
    //! \param  with_blocked  if true, also block matrices will be stored in hierarchy
    TMatrixHierarchy ( TMatrix< value_t > *  A,
                       const bool            with_blocked );

    //! destructor
    ~TMatrixHierarchy ();
    
    //////////////////////////////////////////////////////////////
    //
    // manage hierarchy
    //

    //! return number of levels in hierarchy
    size_t  n_levels () const { return _hierarchy.size(); }
    
    //! return blockmatrix for given level
    TSparseBlockMatrix< value_t > *
    matrix ( const size_t  lvl )
    {
        if ( lvl < _hierarchy.size() )
            return _hierarchy[ lvl ];
        else
            return nullptr;
    }

    //! return matrix defined by is0×is1 on given level
    TMatrix< value_t > *
    block ( const size_t       lvl,
            const TIndexSet &  row_is,
            const TIndexSet &  col_is )
    {
        TMatrix< value_t > *  M = nullptr;
        
        if ( lvl < _hierarchy.size() )
        {
            TSparseBlockMatrix< value_t > *  B = _hierarchy[ lvl ];

            if ( B != nullptr )
                M = B->block( row_is, col_is );
        }// if

        return M;
    }

    //! return matrix t×s with is0×is1 ⊆ t×s on given level
    TMatrix< value_t > *
    block_containing ( const size_t       lvl,
                       const TIndexSet &  row_is,
                       const TIndexSet &  col_is )
    {
        TMatrix< value_t > *  M = nullptr;
        
        if ( lvl < _hierarchy.size() )
        {
            TSparseBlockMatrix< value_t > *  B = _hierarchy[ lvl ];

            if ( B != nullptr )
                M = B->block_containing( row_is, col_is );
        }// if

        return M;
    }

    /////////////////////////////////////////////////
    //
    // misc.
    //

    //! return size in bytes used by this object
    size_t byte_size () const;

    //! print content of matrix hierarchy
    void   print ( const uint  ofs = 0 ) const;
};

}// namespace

#endif  // __HPRO_TMATRIXHIERARCHY_HH
