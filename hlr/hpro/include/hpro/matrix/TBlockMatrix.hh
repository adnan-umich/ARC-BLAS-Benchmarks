#ifndef __HPRO_TBLOCKMATRIX_HH
#define __HPRO_TBLOCKMATRIX_HH
//
// Project     : HLIBpro
// File        : TBlockMatrix.hh
// Description : class for a matrix consisting of submatrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <list>
#include <vector>

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"

namespace Hpro
{

// local matrix type
DECLARE_TYPE( TBlockMatrix );

//!
//! \ingroup  Matrix_Module
//! \class    TBlockMatrix
//! \brief    Class for a n×m block matrix of TMatrix sub matrices.
//!
template < typename T_value >
class TBlockMatrix : public TMatrix< T_value >
{
public:
    using  value_t = T_value;
    using  real_t  = typename real_type< value_t >::type_t;

private:
    //! @cond

    //! number of rows of matrix
    size_t                               _rows;

    //! number of columns of matrix
    size_t                               _cols;
    
    //! number of block rows of matrix
    uint                                 _block_rows;

    //! number of block columns of matrix
    uint                                 _block_cols;
    
    //! array of sub-blocks of this matrix (column-wise storage !!!)
    std::vector< TMatrix< value_t > * >  _blocks;

    //! @endcond
    
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct block matrix with size and block structure defined by \a bct
    TBlockMatrix ( const TBlockCluster * bct = nullptr )
            : TMatrix< value_t >( bct )
            , _block_rows(0)
            , _block_cols(0)
    {
        if ( bct != nullptr ) set_cluster( bct );
        else                  set_block_struct( 1, 1 );
    }

    //! construct block matrix with over block index set \a row_is × \a col_is
    TBlockMatrix ( const TIndexSet &  row_is,
                   const TIndexSet &  col_is )
            : TMatrix< value_t >()
            , _block_rows(0)
            , _block_cols(0)
    {
        this->set_ofs(  row_is.first(), col_is.first() );
        set_size( row_is.size(),  col_is.size() );
    }

    //! construct block matrix with over block index set \a bis
    TBlockMatrix ( const TBlockIndexSet &  bis )
            : TMatrix< value_t >()
            , _block_rows(0)
            , _block_cols(0)
    {
        this->set_ofs(  bis.row_is().first(), bis.col_is().first() );
        set_size( bis.row_is().size(),  bis.col_is().size() );
    }

    //! dtor
    virtual ~TBlockMatrix ();

    //! set cluster of matrix
    virtual void          set_cluster ( const TBlockCluster * c );
    
    ///////////////////////////////////////////
    //
    // block-handling
    //

    //! return number of rows
    virtual size_t        rows () const { return _rows; }

    //! return number of columns
    virtual size_t        cols () const { return _cols; }

    //! set size of matrix
    virtual void          set_size ( const size_t  r, const size_t  c ) { _rows = r; _cols = c; }
    
    //! return number of block-rows
    uint                  block_rows   ()                    const { return _block_rows; }
    uint                  block_rows   ( const matop_t  op ) const { return ( op == apply_normal
                                                                              ? block_rows()
                                                                              : block_cols() ); }
    uint                  nblock_rows  ()                    const { return block_rows(); }
    uint                  nblock_rows  ( const matop_t  op ) const { return block_rows( op ); }
    
    //! return number of block-columns
    uint                  block_cols   ()                    const { return _block_cols; }
    uint                  block_cols   ( const matop_t  op ) const { return ( op == apply_normal
                                                                              ? block_cols()
                                                                              : block_rows() ); }
    uint                  nblock_cols  ()                    const { return block_cols(); }
    uint                  nblock_cols  ( const matop_t  op ) const { return block_cols( op ); }

    //! return number of sub blocks of matrix
    uint                  no_of_blocks () const { return _block_rows*_block_cols; }

    //! set block structure
    void                  set_block_struct ( const uint n, const uint m );

    //! set matrix format
    virtual void          set_form         ( const matform_t  f );
    
    //! set matrix format
    virtual void          set_form         ( const matform_t         f,
                                             const recursion_type_t  rec_type );
    
    //! return matrix coefficient (\a i, \a j)
    virtual value_t       entry  ( const idx_t i, const idx_t j ) const;

    //! return block at block index (\a i, \a j)
    TMatrix< value_t > *        block ( const uint i, const uint j )       { return _blocks[(j*_block_rows)+i]; }

    //! return block at block index (\a i, \a j)
    const TMatrix< value_t > *  block ( const uint i, const uint j ) const { return _blocks[(j*_block_rows)+i]; }

    //! return block at block index (\a i, \a j)
    TMatrix< value_t > *        block ( const uint     i,
                                        const uint     j,
                                        const matop_t  op )                { return ( op == apply_normal
                                                                                      ? _blocks[(j*_block_rows)+i]
                                                                                      : _blocks[(i*_block_rows)+j] ); }

    //! return block at block index (\a i, \a j)
    const TMatrix< value_t > *  block ( const uint     i,
                                        const uint     j,
                                        const matop_t  op ) const          { return ( op == apply_normal
                                                                                      ? _blocks[(j*_block_rows)+i]
                                                                                      : _blocks[(i*_block_rows)+j] ); }

    //! set matrix block at block index (\a i,\a j) to matrix \a A
    void                        set_block ( const uint             i,
                                            const uint             j,
                                            TMatrix< value_t > *   A )
    {
        _blocks[(j*_block_rows)+i] = A;

        if ( A != nullptr )
            A->set_parent( this );
    }

    //! replace matrix block \a A by matrix \a B (don't delete A !)
    void                  replace_block ( TMatrix< value_t > * A,
                                          TMatrix< value_t > * B );
    
    //! delete block (i,j)
    void                  delete_block  ( const uint i, const uint j )
    {
        delete _blocks[(j*_block_rows)+i];
        _blocks[(j*_block_rows)+i] = nullptr;
    }
    
    //! return subblock of matrix corresponding to block cluster \a t
    TMatrix< value_t > *  bc_block ( const TBlockCluster * t ) const;

    //! return subblock of matrix corresponding to block cluster (\a tau, \a sigma)
    TMatrix< value_t > *  bc_block ( const TCluster * tau, const TCluster * sigma ) const;

    //! clear pointers to all sub blocks
    void                  clear_blocks ();
    
    //! return true, if matrix is blocked
    virtual bool          is_blocked () const { return true; }

    //! set processor set of matrix
    virtual void          set_procs  ( const TProcSet &        ps,
                                       const recursion_type_t  rec_type = nonrecursive );

    ///////////////////////////////////////////
    //
    // management of update accumulator
    //

    //! apply stored updates U to local matrix M, e.g., M = M + U,
    //! with accuracy \a acc
    virtual void  apply_updates ( const TTruncAcc &       acc,
                                  const recursion_type_t  recursion );
    
    //! return true, if matrix has updates not yet applied
    virtual bool  has_updates   ( const recursion_type_t  recursion ) const;

    /////////////////////////////////////////////////
    //
    // BLAS-routines
    //

    //! transpose matrix
    virtual void transpose ();
    
    //! conjugate matrix coefficients
    virtual void conjugate ();
    
    //! compute this ≔ α·this
    virtual void scale ( const value_t alpha );
    
    //! compute this ≔ this + α · matrix
    virtual void add ( const value_t alpha, const TMatrix< value_t > * matrix );
    
    //! compute y ≔ β·y + α·op(M)·x, with M = this
    virtual void mul_vec ( const value_t               alpha,
                           const TVector< value_t > *  x,
                           const value_t               beta,
                           TVector< value_t > *        y,
                           const matop_t               op = apply_normal ) const;
    using TMatrix< value_t >::mul_vec;
        
    ///////////////////////////////////////////////////////////
    //
    // linear operator mapping
    //

    //! same as above but only the dimension of the vector spaces is tested,
    //! not the corresponding index sets
    virtual void  apply_add   ( const value_t                    alpha,
                                const BLAS::Vector< value_t > &  x,
                                BLAS::Vector< value_t > &        y,
                                const matop_t                    op = apply_normal ) const;

    virtual void  apply_add   ( const value_t                    alpha,
                                const BLAS::Matrix< value_t > &  X,
                                BLAS::Matrix< value_t > &        Y,
                                const matop_t                    op = apply_normal ) const;

    using TMatrix< value_t >::apply_add;

    /////////////////////////////////////////////////
    //
    // misc.
    //

    //! collect matrix blocks corresponding to leaves in list \a leaf_list
    template < typename T_list >
    void collect_leaves ( T_list &  leaf_list ) const
    {
        for ( uint i = 0; i < block_rows(); i++ )
            for ( uint j = 0; j < block_cols(); j++ )
            {
                const TMatrix< value_t >  * A_ij = block(i,j);
            
                if ( A_ij == nullptr )
                    continue;
            
                if ( ! IS_TYPE( A_ij, TBlockMatrix ) )
                    leaf_list.push_back( const_cast< TMatrix< value_t > * >( A_ij ) );
                else
                    cptrcast( A_ij, TBlockMatrix )->collect_leaves( leaf_list );
            }// for
    }

    //! truncate all rank-blocks in matrix to accuracy \a acc
    virtual void truncate ( const TTruncAcc & acc );
    
    //! print matrix to stdout
    virtual void print ( const uint ofs = 0 ) const;
    
    //! return matrix of same class (but no content)
    virtual auto create () const -> std::unique_ptr< TMatrix< value_t > >
    {
        return std::make_unique< TBlockMatrix< value_t > >();
    }

    //! return copy of matrix
    virtual auto copy         () const -> std::unique_ptr< TMatrix< value_t > >;

    //! copy matrix wrt. accuracy \a acc and optional coarsening
    virtual auto copy         ( const TTruncAcc & acc,
                                const bool        coarsen = false ) const -> std::unique_ptr< TMatrix< value_t > >;

    //! return structural copy of matrix
    virtual auto copy_struct  () const -> std::unique_ptr< TMatrix< value_t > >;

    //! copy matrix into \a A
    virtual void copy_to      ( TMatrix< value_t > *  A ) const;

    //! copy matrix into \a A with accuracy \a acc and optional coarsening
    virtual void copy_to      ( TMatrix< value_t > *  A,
                                const TTruncAcc &     acc,
                                const bool            coarsen = false ) const;

    //! return size in bytes used by this object
    virtual size_t byte_size () const;

    //! return size of (floating point) data in bytes handled by this object
    virtual size_t data_byte_size () const;

    //! copy complete structural information from given matrix
    virtual void copy_struct_from ( const TMatrix< value_t > *  M );
    
    //! copy complete structural information from given matrix with different type
    template < typename T_value_M >
    void copy_struct_from_all ( const TBlockMatrix< T_value_M > * M )
    {
        TMatrix< value_t >::copy_struct_from_all( M );
        set_block_struct( M->block_rows(), M->block_cols() );
    }
        
    //
    // RTTI
    //

    HPRO_RTTI_DERIVED( TBlockMatrix, TMatrix< value_t > )

    //
    // serialisation
    //

    //! read data from stream \a s and copy to matrix
    virtual void read  ( TByteStream &  s );

    //! use data from stream \a s to build matrix
    virtual void build ( TByteStream &  s );

    //! write data to stream \a s
    virtual void write ( TByteStream &  s ) const;

    //! returns size of object in bytestream
    virtual size_t bs_size () const;

    //
    // tests
    //

    //! test data for invalid values, e.g. INF and NAN
    virtual void check_data () const;
};

}// namespace Hpro

#endif  // __HPRO_TBLOCKMATRIX_HH
