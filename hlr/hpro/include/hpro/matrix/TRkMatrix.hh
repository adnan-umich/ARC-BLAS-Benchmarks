#ifndef __HPRO_TRKMATRIX_HH
#define __HPRO_TRKMATRIX_HH
//
// Project     : HLIBpro
// File        : TRkMatrix.hh
// Description : class for rank-k-matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/TTruncAcc.hh"
#include "hpro/blas/Matrix.hh"
#include "hpro/matrix/TMatrix.hh"
#include "hpro/vector/TScalarVector.hh"

namespace Hpro
{

// forward decl.
template < typename value_t >
class TDenseMatrix;

// local matrix type
DECLARE_TYPE( TRkMatrix );

//!
//! \ingroup Matrix_Module
//! \class   TRkMatrix
//! \brief   Represents low rank matrices in factored form: \f$ M = A B^H \f$.
//!
template < typename T_value >
class TRkMatrix : public TMatrix< T_value >
{
public:
    using  value_t = T_value;
    using  real_t  = typename real_type< value_t >::type_t;

protected:
    //! @cond
    
    // factors of low-rank representation
    BLAS::Matrix< value_t >  _mat_A, _mat_B;

    // size of the vectors in A and B
    size_t                   _rows, _cols;
    
    // current rank of the matrix
    size_t                   _rank;
    
    //! @endcond
    
public:
    /////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct zero sized low-rank matrix
    TRkMatrix ();

    //! construct low-rank matrix of size \a rows × \a cols
    TRkMatrix ( const size_t                     rows,
                const size_t                     cols );

    //! construct low-rank matrix of size defined by block index set
    TRkMatrix ( const TBlockIndexSet &           block_is );
    
    //! construct low-rank matrix of size defined by block index set
    TRkMatrix ( const TIndexSet &                arow_is,
                const TIndexSet &                acol_is );
    
    //! construct low-rank matrix of size defined by block index set
    //! and real factors \a A and \a B
    TRkMatrix ( const TIndexSet &                arow_is,
                const TIndexSet &                acol_is,
                const BLAS::Matrix< value_t > &  A,
                const BLAS::Matrix< value_t > &  B );

    //! construct low-rank matrix of size defined by block index set
    //! and real factors \a A and \a B (move version)
    TRkMatrix ( const TIndexSet &                arow_is,
                const TIndexSet &                acol_is,
                BLAS::Matrix< value_t > &&       A,
                BLAS::Matrix< value_t > &&       B );

    //! construct low-rank matrix of size defined by block cluster
    TRkMatrix ( const TBlockCluster *            bc );

    //! copy constructor
    TRkMatrix ( const TRkMatrix< value_t > &     A );

    //! destructor
    ~TRkMatrix ()
    {}

    /////////////////////////////////////////////////
    //
    // access local variables
    //

    //! return true, if matrix is zero
    virtual bool    is_zero   () const { return (( _rank == 0 ) && ! this->accumulator().has_updates()); }
    
    //
    // access actual data of matrix factors as BLAS::Matrix
    //
    
    //! return pointer to internal matrix data of matrix A
    BLAS::Matrix< value_t > &           blas_mat_A  ()       { return _mat_A; }
    
    //! return const pointer to internal matrix data of matrix A
    const BLAS::Matrix< value_t > &     blas_mat_A  () const { return _mat_A; }

    //! return pointer to internal matrix data of matrix B
    BLAS::Matrix< value_t > &           blas_mat_B  ()       { return _mat_B; }
    
    //! return const pointer to internal matrix data of matrix B
    const BLAS::Matrix< value_t > &     blas_mat_B  () const { return _mat_B; }

    //
    // access individual vectors in A and B as BLAS vectors
    //
    
    //! return vector A_i
    const BLAS::Vector< value_t >   blas_vec_A  ( const idx_t  i ) const { return _mat_A.column( i ); }

    //! return vector B_i
    const BLAS::Vector< value_t >   blas_vec_B  ( const idx_t  i ) const { return _mat_B.column( i ); }
    
    //
    // access individual vectors in A and B as H vectors
    //
    
    //! return vector A_i
    const TScalarVector< value_t >  vec_A   ( const idx_t  i ) const
    {
        return TScalarVector< value_t >( this->row_is(), _mat_A.column( i ) );
    }

    //! return vector B_i
    const TScalarVector< value_t >  vec_B   ( const idx_t  i ) const
    {
        return TScalarVector< value_t >( this->col_is(), _mat_B.column( i ) );
    }
    
    //
    // manage rank of matrix
    //
    
    //! return rank of matrix
    size_t          rank         () const { return _rank; }

    //! set rank of matrix without truncating data
    void            set_rank     ( const size_t  k );

    //! compute actual rank of matrix (remove zero vectors)
    void            comp_rank    ();

    //
    // access size information
    //
    
    //! return number of rows 
    virtual size_t  rows         () const { return _rows; }

    //! return number of columns 
    virtual size_t  cols         () const { return _cols; }

    //
    // change size
    //
    
    //! set new cluster over which matrix is defined and change size accordingly
    virtual void    set_cluster  ( const TBlockCluster * c );

    //! set size and rank of matrix (if zero == true fill new memory with zeros)
    virtual void    set_size     ( const size_t  n,
                                   const size_t  m,
                                   const size_t  k );
    
    //! set new size but keep current rank of matrix
    virtual void    set_size     ( const size_t  n,
                                   const size_t  m )
    {
        return set_size( n, m, _rank );
    }

    //
    // access size information
    //
    
    //! return coefficent (\a i, \a j) of matrix
    value_t            entry        ( const idx_t  i, const idx_t j ) const;

    ///////////////////////////////////////////
    //
    // management of update accumulator
    //

    //! apply stored updates U to local matrix M, e.g., M = M + U,
    //! with accuracy \a acc
    virtual void    apply_updates ( const TTruncAcc &       acc,
                                    const recursion_type_t  recursion );
    
    /////////////////////////////////////////////////
    //
    // manage Rk-matrices
    //

    //! transpose matrix
    virtual void    transpose    ();
    
    //! conjugate matrix coefficients
    virtual void    conjugate    ();
    
    //! truncate matrix  w.r.t. accuracy \a acc
    virtual void    truncate     ( const TTruncAcc & acc );

    //! set this ≔ A·B^H
    virtual void    set_lrmat    ( const BLAS::Matrix< value_t > &  A,
                                   const BLAS::Matrix< value_t > &  B );
    
    virtual void    set_lrmat    ( BLAS::Matrix< value_t > &&  A,
                                   BLAS::Matrix< value_t > &&  B )
    {
        if (( A.nrows() != _mat_A.nrows() ) ||
            ( B.nrows() != _mat_B.nrows() ) ||
            ( A.ncols() != B.ncols()))
            HERROR( ERR_ARG, "(TRkMatrix) set_lrmat", "input matrices have invalid dimension" );
        
        _rank  = A.ncols();
        _mat_A = std::move( A );
        _mat_B = std::move( B );
    }
    
    //! compute this ≔ this + α·A·B^H and truncate result w.r.t. \a acc (real valued variant)
    void            add_rank     ( const value_t                    alpha,
                                   const BLAS::Matrix< value_t > &  A,
                                   const BLAS::Matrix< value_t > &  B,
                                   const TTruncAcc &                acc );
    
    //! add a dense matrix and truncate w.r.t. \a acc (real valued variant)
    void            add_dense    ( const value_t                    alpha,
                                   const BLAS::Matrix< value_t > &  D,
                                   const TTruncAcc &                acc );
    
    //! copy a densematrix (nxm) as low-rank matrix (rank = min{n,m})
    void            copy_dense   ( const TDenseMatrix< value_t > *  A,
                                   const TTruncAcc &                acc );
    
    /////////////////////////////////////////////////
    //
    // BLAS-routines
    //

    //! scale matrix by constant factor
    virtual void         scale      ( const value_t       f );
    
    //! compute y ≔ α·M·x + β·y, where M is either this, this^T or this^H depending on \a op
    virtual void         mul_vec    ( const value_t               alpha,
                                      const TVector< value_t > *  x,
                                      const value_t               beta,
                                      TVector< value_t > *        y,
                                      const matop_t               op = apply_normal ) const;
    using TMatrix< value_t >::mul_vec;

    //! compute this = this + α·A without truncation
    virtual void         add        ( const value_t               alpha,
                                      const TMatrix< value_t > *  A );
        
    //! compute matrix product α·op_A(this)·op_B(B) 
    virtual TRkMatrix *  mul_right  ( const value_t               alpha,
                                      const TMatrix< value_t > *  B,
                                      const matop_t               op_A,
                                      const matop_t               op_B ) const;

    //! compute matrix product α·op_A(A)·op_B(this) 
    virtual TRkMatrix *  mul_left   ( const value_t               alpha,
                                      const TMatrix< value_t > *  A,
                                      const matop_t               op_A,
                                      const matop_t               op_B ) const;

    ///////////////////////////////////////////////////////////
    //
    // linear operator mapping
    //

    //! same as above but only the dimension of the vector spaces is tested,
    //! not the corresponding index sets
    virtual void  apply_add   ( const value_t                     alpha,
                                const BLAS::Vector< value_t > &   x,
                                BLAS::Vector< value_t > &         y,
                                const matop_t                     op = apply_normal ) const;

    using TMatrix< value_t >::apply_add;
    
    /////////////////////////////////////////////////
    //
    // misc methods
    //

    //
    // virtual constructors
    //

    //! return object of same type
    virtual auto   create       () const -> std::unique_ptr< TMatrix< value_t > >
    {
        return std::make_unique< TRkMatrix< value_t > >();
    }

    //! return copy of matrix
    virtual auto   copy         () const -> std::unique_ptr< TMatrix< value_t > >;

    //! return copy matrix wrt. given accuracy; if \a do_coarsen is set, perform coarsening
    virtual auto   copy         ( const TTruncAcc &  acc,
                                  const bool         do_coarsen = false ) const -> std::unique_ptr< TMatrix< value_t > >;

    //! return structural copy of matrix
    virtual auto   copy_struct  () const -> std::unique_ptr< TMatrix< value_t > >;

    // copy matrix data to \a A
    virtual void   copy_to      ( TMatrix< value_t > *          A ) const;

    // copy matrix data to \a A and truncate w.r.t. \acc with optional coarsening
    virtual void   copy_to      ( TMatrix< value_t > *          A,
                                  const TTruncAcc &  acc,
                                  const bool         do_coarsen = false ) const;
    
    //! return size in bytes used by this object
    virtual size_t byte_size    () const;

    //! return size of (floating point) data in bytes handled by this object
    virtual size_t data_byte_size () const { return _mat_A.data_byte_size() + _mat_B.data_byte_size(); }
    
    //
    // RTTI
    //

    HPRO_RTTI_DERIVED( TRkMatrix, TMatrix< value_t > )

    //
    // serialisation
    //

    //! read matrix from byte stream
    virtual void read  ( TByteStream & s );

    //! construct matrix from byte stream
    virtual void build ( TByteStream & s );
    
    //! write matrix into byte stream
    virtual void write ( TByteStream & s ) const;

    //! return size of object in a bytestream
    virtual size_t bs_size () const;

    //
    // tests
    //

    //! test data for invalid values, e.g. INF and NAN
    virtual void check_data () const;
};

//
// wrapper for matrix_A/B and cmatrix_A/B to avoid
// if clauses ( if ( is_complex() ) )
//

template < typename value_t >  BLAS::Matrix< value_t > &        blas_mat_A  ( TRkMatrix< value_t > *        A ) { return A->blas_mat_A(); }
template < typename value_t >  const BLAS::Matrix< value_t > &  blas_mat_A  ( const TRkMatrix< value_t > *  A ) { return A->blas_mat_A(); }
template < typename value_t >  BLAS::Matrix< value_t > &        blas_mat_B  ( TRkMatrix< value_t > *        A ) { return A->blas_mat_B(); }
template < typename value_t >  const BLAS::Matrix< value_t > &  blas_mat_B  ( const TRkMatrix< value_t > *  A ) { return A->blas_mat_B(); }

template < typename value_t >  BLAS::Matrix< value_t > &        blas_mat_A  ( TRkMatrix< value_t > &        A ) { return A.blas_mat_A(); }
template < typename value_t >  const BLAS::Matrix< value_t > &  blas_mat_A  ( const TRkMatrix< value_t > &  A ) { return A.blas_mat_A(); }
template < typename value_t >  BLAS::Matrix< value_t > &        blas_mat_B  ( TRkMatrix< value_t > &        A ) { return A.blas_mat_B(); }
template < typename value_t >  const BLAS::Matrix< value_t > &  blas_mat_B  ( const TRkMatrix< value_t > &  A ) { return A.blas_mat_B(); }

template < typename value_t >  BLAS::Matrix< value_t > &        blas_mat_A  ( std::unique_ptr< TRkMatrix< value_t > > & M ) { return M->blas_mat_A(); }
template < typename value_t >  BLAS::Matrix< value_t > &        blas_mat_B  ( std::unique_ptr< TRkMatrix< value_t > > & M ) { return M->blas_mat_B(); }

}// namespace Hpro

#endif  // __HPRO_TRKMATRIX_HH
