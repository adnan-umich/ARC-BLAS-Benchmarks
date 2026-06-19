#ifndef __HPRO_TDENSEMATRIX_HH
#define __HPRO_TDENSEMATRIX_HH
//
// Project     : HLIBpro
// File        : TDenseMatrix.hh
// Description : class for dense matrices of arbitrary (small) size
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/types.hh"
#include "hpro/blas/Matrix.hh"

#include "hpro/cluster/TPermutation.hh"

#include "hpro/matrix/TMatrix.hh"

#include "hpro/vector/TScalarVector.hh"

namespace Hpro
{

// local matrix type
DECLARE_TYPE( TDenseMatrix );

//!
//! \ingroup Matrix_Module
//! \class   TDenseMatrix
//! \brief   Represent a dense matrix
//!
template < typename T_value >
class TDenseMatrix : public TMatrix< T_value >
{
public:
    using  value_t = T_value;
    using  real_t  = typename real_type< value_t >::type_t;

protected:
    //! @cond
    
    //! number of rows 
    size_t                   _rows;

    //! number of columns rows 
    size_t                   _cols;

    //! matrix data
    BLAS::Matrix< value_t >  _mat;

    //! @endcond
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct zero sized matrix
    TDenseMatrix ()
            : _rows(0), _cols(0)
    {}

    //! construct matrix of size \a n × \a m
    TDenseMatrix ( const size_t  n,
                   const size_t  m )
            : TMatrix< value_t >()
            , _rows( 0 )
            , _cols( 0 )
    {
        set_size( n, m );
    }
    
    //! construct matrix with size defined by \a arow_is × \a acol_is
    TDenseMatrix ( const TIndexSet &   arow_is,
                   const TIndexSet &   acol_is )
            : TMatrix< value_t >()
            , _rows( 0 )
            , _cols( 0 )
    {
        this->set_block_is( TBlockIndexSet( arow_is, acol_is ) );
    }

    //! construct matrix with size defined by \a arow_is × \a acol_is
    //! and data given by \a M
    TDenseMatrix ( const TIndexSet &                arow_is,
                   const TIndexSet &                acol_is,
                   const BLAS::Matrix< value_t > &  M )
            : TMatrix< value_t >()
            , _rows( arow_is.size() )
            , _cols( acol_is.size() )
            , _mat( M )
    {
        // do not call set_block_is to avoid size initialisation
        this->set_ofs( arow_is.first(), acol_is.first() );
    }

    //! construct matrix with size defined by \a arow_is × \a acol_is
    //! and move data from \a M
    TDenseMatrix ( const TIndexSet &           arow_is,
                   const TIndexSet &           acol_is,
                   BLAS::Matrix< value_t > &&  M )
            : TMatrix< value_t >()
            , _rows( arow_is.size() )
            , _cols( acol_is.size() )
            , _mat( std::move( M ) )
    {
        // do not call set_block_is to avoid size initialisation
        this->set_ofs( arow_is.first(), acol_is.first() );
    }
    
    //! copy constructor
    TDenseMatrix ( const TDenseMatrix< value_t > &  mat )
            : TMatrix< value_t >()
            , _rows(0)
            , _cols(0)
    {
        mat.copy_to( this );
    }

    //! construct matrix with size defined by block cluster
    TDenseMatrix ( const TBlockCluster *  bct )
            : TMatrix< value_t >( bct )
            , _rows(0)
            , _cols(0)
    {
        set_cluster( bct );
    }

    //! destructor
    virtual ~TDenseMatrix ()
    {}

    ////////////////////////////////////////////////
    //
    // manage internal data
    //

    //! return number of rows
    virtual size_t  rows () const { return _rows; }

    //! return number of columns
    virtual size_t  cols () const { return _cols; }

    //! set size of matrix to \a n × \a m
    virtual void    set_size ( const size_t  n,
                               const size_t  m );
    
    //! set size as defined by block cluster \a c
    virtual void    set_cluster ( const TBlockCluster * c );

    //! return true, if matrix is dense
    virtual bool    is_dense () const { return true; }
    
    //
    // return internal BLAS matrices
    //

    //! return BLAS matrix
    BLAS::Matrix< value_t > &        blas_mat  ()       { return _mat; }

    //! return constant BLAS matrix
    const BLAS::Matrix< value_t > &  blas_mat  () const { return _mat; }

    //
    // return rows and columns
    //

    //! return row \a i as vector object
    const TScalarVector< value_t >
    row  ( const idx_t  i ) const
    {
        return TScalarVector< value_t >( this->col_is(),
                                         BLAS::Vector< value_t >( _mat, i, BLAS::Range( 0, idx_t(_cols)-1 ) ) );
    }

    //! return column \a i as vector object
    const TScalarVector< value_t >
    column  ( const idx_t  i ) const
    {
        return TScalarVector< value_t >( this->row_is(),
                                         BLAS::Vector< value_t >( _mat, BLAS::Range( 0, idx_t(_rows)-1 ), i ) );
    }

    //! return row \a i as BLAS vector
    const BLAS::Vector< value_t >
    blas_row ( const idx_t  i ) const
    {
        return BLAS::Vector< value_t >( _mat, i, BLAS::Range( 0, idx_t(_cols)-1 ) );
    }

    //! return column \a i as BLAS vector
    const BLAS::Vector< value_t >
    blas_col ( const idx_t  i ) const 
    {
        return BLAS::Vector< value_t >( _mat, BLAS::Range( 0, idx_t(_rows)-1 ), i );
    }

    //
    // and to pointers to matrix entries
    //

    //! return pointer to data starting at coefficient (\a i, \a j)
    value_t *  entry_ptr ( const idx_t i, const idx_t j )
    {
        HASSERT( (i < idx_t(_rows)) && ( j < idx_t(_cols)), ERR_ARR_BOUND, "(TDenseMatrix) entry_ptr", "" );
        return & _mat( i, j );
    }

    ///////////////////////////////////////////////
    //
    // access coeff.
    //

    //! return coefficient (\a i, \a j)
    value_t  entry ( const idx_t i, const idx_t j ) const
    {
        HASSERT( (i < idx_t(_rows)) && ( j < idx_t(_cols)), ERR_ARR_BOUND, "(TDenseMatrix) entry", "" );
        return _mat( i, j );
    }

    //! set coefficient (\a i, \a j) to \a f
    void set_entry  ( const idx_t i, const idx_t j, const value_t  f )
    {
        HASSERT( (i < idx_t(_rows)) && ( j < idx_t(_cols)), ERR_ARR_BOUND, "(TDenseMatrix) set_entry", "" );
        _mat( i, j ) = f;
    }

    //! add \a f to coefficient (\a i, \a j)
    void add_entry  ( const idx_t i, const idx_t j, const value_t f )
    {
        HASSERT( (i < idx_t(_rows)) && ( j < idx_t(_cols)), ERR_ARR_BOUND, "(TDenseMatrix) add_entry", "" );
        _mat( i, j ) += f;
    }

    /////////////////////////////////////////////////
    //
    // BLAS-routines
    //

    // scale matrix by constant factor \a f
    virtual void                   scale      ( const value_t            f );
    
    //! compute y ≔ α·op(this)·x + β·y
    virtual void                   mul_vec    ( const value_t               alpha,
                                                const TVector< value_t > *  x,
                                                const value_t               beta,
                                                TVector< value_t > *        y,
                                                const matop_t               op = apply_normal ) const;
    using TMatrix< value_t >::mul_vec;

    //! compute this = this + α·A without truncation
    virtual void                  add        ( const value_t                alpha,
                                               const TMatrix< value_t > *   A );

    //! compute matrix product α·op_A(this)·op_B(B) 
    virtual TMatrix< value_t > *  mul_right  ( const value_t                alpha,
                                               const TMatrix< value_t > *   B,
                                               const matop_t                op_A,
                                               const matop_t                op_B ) const;
    
    //! compute matrix product α·op_A(A)·op_B(this) 
    virtual TMatrix< value_t > *  mul_left   ( const value_t                alpha,
                                               const TMatrix< value_t > *   A,
                                               const matop_t                op_A,
                                               const matop_t                op_B ) const;
    
    //! compute α·this + β·op(M), where op(M) is subblock of this or this a subblock of op(M))
    void                          add_block  ( const value_t                    alpha,
                                               const value_t                    beta,
                                               const TDenseMatrix< value_t > *  M,
                                               const matop_t                    op = apply_normal );
    
    ///////////////////////////////////////////////////////////
    //
    // linear operator mapping
    //

    //! same as above but only the dimension of the vector spaces is tested,
    //! not the corresponding index sets
    virtual void  apply_add   ( const value_t                       alpha,
                                const BLAS::Vector< value_t > &     x,
                                BLAS::Vector< value_t > &           y,
                                const matop_t                       op = apply_normal ) const;

    using TMatrix< value_t >::apply_add;

    /////////////////////////////////////////////////
    //
    // misc. methods
    //

    //! transpose matrix
    virtual void transpose ();
    
    //! conjugate matrix coefficients
    virtual void conjugate ();
    
    //! truncate matrix to given accuracy: nothing to be done
    virtual void truncate ( const TTruncAcc & ) {}

    //! copy operator
    TDenseMatrix< value_t > &  operator = ( const TDenseMatrix< value_t > &  mat );

    //
    // virtual constructors
    //

    //! return matrix object of same class as this
    virtual auto  create       () const -> std::unique_ptr< TMatrix< value_t > >
    {
        return std::make_unique< TDenseMatrix< value_t > >();
    }

    //! return copy of matrix
    virtual auto  copy         () const -> std::unique_ptr< TMatrix< value_t > >;
    using TMatrix< value_t >::copy;

    //! return structural copy of matrix
    virtual auto  copy_struct  () const -> std::unique_ptr< TMatrix< value_t > >;

    //! copy matrix into A
    virtual void  copy_to      ( TMatrix< value_t > * A ) const;
    using TMatrix< value_t >::copy_to;
    
    //! permute rows/columns in matrix according to \a row_pern and \a col_perm
    //! - nullptr is treated as identity
    void permute ( const TPermutation *  row_perm,
                   const TPermutation *  col_perm );

    //! return size in bytes used by this object
    virtual size_t byte_size () const;

    //! return size of (floating point) data in bytes handled by this object
    virtual size_t data_byte_size () const { return _mat.data_byte_size(); }

    //
    // RTTI
    //

    HPRO_RTTI_DERIVED( TDenseMatrix, TMatrix< value_t > )

    //
    // serialisation
    //

    //! read matrix from byte stream
    virtual void    read     ( TByteStream & s );

    //! construct matrix based on data in byte stream
    virtual void    build    ( TByteStream & s );

    //! write matrix to byte stream
    virtual void    write    ( TByteStream & s ) const;

    //! returns size of object in bytestream
    virtual size_t  bs_size  () const;

    //
    // tests
    //

    //! test data for invalid values, e.g. INF and NAN
    virtual void check_data () const;
};

//
// wrapper for entry_ptr/centry_ptr (relative to global indexset)
//
template < typename value_t >
value_t *
entry_ptr       ( const TDenseMatrix< value_t > *  A,
                  const idx_t                      i,
                  const idx_t                      j )
{
    return const_cast< TDenseMatrix< value_t > * >( A )->entry_ptr( i - A->row_ofs(), j - A->col_ofs() );
}

//
// wrapper for entry_ptr/centry_ptr (relative to local indexset)
//
template < typename value_t >
value_t *
entry_ptr_loc  ( const TDenseMatrix< value_t > *  A,
                 const idx_t                      i,
                 const idx_t                      j )

{
    return const_cast< TDenseMatrix< value_t > * >( A )->entry_ptr( i, j );
}

//
// wrapper for set_entry/set_centry
//
template <typename value_t>
void
set_entry ( TDenseMatrix< value_t > *  A,
            const idx_t                i,
            const idx_t                j,
            const value_t              val )
{
    return A->set_entry( i, j, val );
}

//
// wrapper for BLAS matrix
//
template < typename value_t > const BLAS::Matrix< value_t > &  blas_mat  ( const TDenseMatrix< value_t > *  A ) { return A->blas_mat(); }
template < typename value_t >       BLAS::Matrix< value_t > &  blas_mat  (       TDenseMatrix< value_t > *  A ) { return A->blas_mat(); }

template < typename value_t > const BLAS::Matrix< value_t > &  blas_mat  ( const TDenseMatrix< value_t > &  A ) { return A.blas_mat(); }
template < typename value_t >       BLAS::Matrix< value_t > &  blas_mat  (       TDenseMatrix< value_t > &  A ) { return A.blas_mat(); }

template < typename value_t >       BLAS::Matrix< value_t > &  blas_mat  ( std::unique_ptr< TDenseMatrix< value_t > > &  M ) { return M->blas_mat(); }

//
// wrapper for row/column (relative to local indexset)
//
template < typename  value_t > const BLAS::Vector< value_t >  blas_row  ( const TDenseMatrix< value_t > *  A, const idx_t  i ) { return A->blas_row( i ); }
template < typename  value_t > const BLAS::Vector< value_t >  blas_col  ( const TDenseMatrix< value_t > *  A, const idx_t  i ) { return A->blas_col( i ); }

template < typename  value_t > const BLAS::Vector< value_t >  blas_row  ( const TDenseMatrix< value_t > &  A, const idx_t  i ) { return A.blas_row( i ); }
template < typename  value_t > const BLAS::Vector< value_t >  blas_col  ( const TDenseMatrix< value_t > &  A, const idx_t  i ) { return A.blas_col( i ); }

}// namespace Hpro

#endif  // __HPRO_TDENSEMATRIX_HH
