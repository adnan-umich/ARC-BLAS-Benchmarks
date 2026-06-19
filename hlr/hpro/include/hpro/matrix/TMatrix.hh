#ifndef __HPRO_TMATRIX_HH
#define __HPRO_TMATRIX_HH
//
// Project     : HLIBpro
// File        : TMatrix.hh
// Description : baseclass for all matrix-classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/TStreamable.hh"
#include "hpro/base/TTruncAcc.hh"
#include "hpro/blas/types.hh"

#include "hpro/cluster/TBlockCluster.hh"

#include "hpro/vector/TScalarVector.hh"

#include "hpro/parallel/TMutex.hh"

#include "hpro/matrix/TLinearOperator.hh"
#include "hpro/matrix/TUpdateAccumulator.hh"

namespace Hpro
{

//
// describes form of matrix
//
enum matform_t
{
    unsymmetric    = 0,      // matrix is unsymmetric
    symmetric      = 1,      // matrix is symmetric: A = A^T
    hermitian      = 2,      // matrix is hermitian: A = A^H
    
    MATFORM_NONSYM = 0,      // matrix is non-symmetric
    MATFORM_SYM    = 1,      // matrix is symmetric: A = A^T
    MATFORM_HERM   = 2       // matrix is hermitian: A = A^H
};
    
// local matrix type
DECLARE_TYPE( TMatrix );

// forward decl.
template < typename value_t >
class TBlockMatrix;

//!
//! \ingroup Matrix_Module
//! \class   TMatrix
//! \brief   Base class for all matrices, defining basic properties, e.g.
//!          underlying block index and processor set
//!
template < typename T_value >
class TMatrix :
        public TLinearOperator< T_value >,
        public TStreamable,
        public TLockable
{
public:
    using  value_t = T_value;
    using  real_t  = typename real_type< value_t >::type_t;
    
private:
    //! @cond
    
    //! globally unique id
    int                        _id;

    //! corresponding block cluster
    const TBlockCluster *      _cluster;

    //! row offset of the block index set
    idx_t                      _row_ofs;

    //! column offset of the block index set
    idx_t                      _col_ofs;

    //! processors this matrix is on
    TProcSet                   _procs;
    
    //! flag to indicate form of matrix
    matform_t                  _matform;
    
    //! hierarchy data
    TBlockMatrix< value_t > *  _parent;
    TMatrix< value_t > *       _prev_in_block_row;
    TMatrix< value_t > *       _next_in_block_row;
    TMatrix< value_t > *       _prev_in_block_col;
    TMatrix< value_t > *       _next_in_block_col;
    TMatrix< value_t > *       _row_diag;
    TMatrix< value_t > *       _col_diag;

    // accumulator
    TUpdateAccumulator< value_t >  _accumulator;
    
    //! @endcond
    
 public:
    ////////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct zero sized matrix
    TMatrix ();

    //! construct matrix of size defined by block cluster \a bcl
    TMatrix ( const TBlockCluster *   bcl );

    //! construct matrix of size defined by block index set \a rowis × \a colis
    TMatrix ( const TIndexSet &       rowis,
              const TIndexSet &       colis );

    //! construct matrix of size defined by block index set \a bis
    TMatrix ( const TBlockIndexSet &  bis );

    //! copy constructor
    TMatrix ( const TMatrix< value_t > &  A );

    //! dtor
    virtual ~TMatrix () {}

    ////////////////////////////////////////////////////////
    //
    // handle structure of matrix
    //

    //! return ID
    int               id             () const { return _id; }

    //! set ID
    void              set_id         ( const int  aid ) { _id = aid; }

    //! return number of rows
    virtual size_t    rows           () const = 0;
    virtual size_t    nrows          () const { return rows(); }

    //! return number of columns
    virtual size_t    cols           () const = 0;
    virtual size_t    ncols          () const { return cols(); }
    
    //! return number of rows of op(M)
    virtual size_t    nrows          ( const matop_t  op ) const { return ( op == apply_normal ? nrows() : ncols() ); }

    //! return number of columns of op(M)
    virtual size_t    ncols          ( const matop_t  op ) const { return ( op == apply_normal ? ncols() : nrows() ); }
    
    //! return row index set
    TIndexSet         row_is         () const { return TIndexSet( row_ofs(), idx_t(row_ofs() + rows()) - 1 ); }

    //! return column index set
    TIndexSet         col_is         () const { return TIndexSet( col_ofs(), idx_t(col_ofs() + cols()) - 1 ); }

    //! return block index set
    TBlockIndexSet    block_is       () const { return TBlockIndexSet( row_is(), col_is() ); }
    
    //! return row index set w.r.t. given matrix operation
    TIndexSet         row_is         ( const matop_t  op ) const { return op == apply_normal ? row_is() : col_is(); }

    //! return row index set w.r.t. given matrix operation
    TIndexSet         col_is         ( const matop_t  op ) const { return op == apply_normal ? col_is() : row_is(); }

    //! return row index set w.r.t. given matrix operation
    TBlockIndexSet    block_is       ( const matop_t  op ) const { return ( op == apply_normal ?
                                                                            TBlockIndexSet( row_is(), col_is() ) :
                                                                            TBlockIndexSet( col_is(), row_is() ) ); }

    //! return first index (number) in row
    virtual idx_t     row_ofs        () const                     { return _row_ofs; }
    virtual idx_t     row_ofs        ( const matop_t  op ) const  { return ( op == apply_normal ? _row_ofs : _col_ofs ); }

    //! return first index (number) in column
    virtual idx_t     col_ofs        () const                     { return _col_ofs; }
    virtual idx_t     col_ofs        ( const matop_t  op ) const  { return ( op == apply_normal ? _col_ofs : _row_ofs ); }

    //! set size of matrix
    virtual void      set_size       ( const size_t n, const size_t m ) = 0;

    //! set index set offsets
    virtual void      set_ofs        ( const idx_t  r, const idx_t  c ) { _row_ofs = r; _col_ofs = c; }

    //! set block index set of matrix
    virtual void      set_block_is   ( const TBlockIndexSet &  is )
    {
        set_size( is.row_is().size(), is.col_is().size() );
        set_ofs( is.row_is().first(), is.col_is().first() );
    }

    //
    // various properties, e.g. matrix format (symmetry ...)
    //
    
    //! return true if matrix is unsymmetric
    bool              is_nonsym       () const              { return is_unsymmetric(); }
    bool              is_unsymmetric  () const              { return _matform == unsymmetric; }

    //! return true if matrix is symmetric
    bool              is_symmetric    () const              { return (( _matform == symmetric ) ||
                                                                     ( this->is_real() && ( _matform == hermitian ) )); }

    //! return true if matrix is hermitian
    bool              is_hermitian    () const              { return (( _matform == hermitian ) ||
                                                                     ( this->is_real() && ( _matform == symmetric ) )); }

    //! return matrix format
    matform_t         form            () const              { return _matform; }

    //! set matrix to be unsymmetric
    void              set_nonsym      ()                    { set_unsymmetric(); }
    void              set_unsymmetric ()                    { set_form( unsymmetric ); }

    //! set matrix to be symmetric
    void              set_symmetric   ()                    { set_form( symmetric ); }

    //! set matrix to be hermitian
    void              set_hermitian   ()                    { set_form( hermitian ); }

    //! set matrix format
    virtual void      set_form        ( const matform_t f ) { _matform = f; }

    //! return true, if matrix is zero
    virtual bool      is_zero         () const              { return false; }

    //! return true, if matrix is blocked
    virtual bool      is_blocked      () const              { return false; }

    //! return true, if matrix is dense
    virtual bool      is_dense        () const              { return false; }
    
    //
    // linear operator properties
    //
    
    //! return true, if operator is self adjoint
    virtual bool      is_self_adjoint () const             { return is_hermitian(); }

    //
    // handle matrix processor
    //

    //! return matrix processor set
    const TProcSet &  procs          () const                { return _procs; }

    //! return number of processors in local set
    uint              nprocs         () const                { return _procs.size(); }

    //! set processor set of matrix
    virtual void      set_procs      ( const TProcSet &        ps,
                                       const recursion_type_t  rec_type = nonrecursive );

    //! return true if matrix is distributed
    bool              is_distributed () const                { return nprocs() > 1; }
    
    ///////////////////////////////////////////
    //!
    //! copy complete structural information from given matrix,
    //! e.g. size, offsets and field type
    //!

    template < typename T_value_M >
    void copy_struct_from_all ( const TMatrix< T_value_M > * M )
    {
        set_id( M->id() );
        set_form( M->form() );
        set_ofs( M->row_ofs(), M->col_ofs() );
        set_size( M->rows(), M->cols() );
        set_procs( M->procs() );
    }
    
    virtual void copy_struct_from ( const TMatrix * M )
    {
        copy_struct_from_all( M );
    }

    ///////////////////////////////////////////
    //
    // access to hierarchy data
    //

    TBlockMatrix< value_t > *  parent            () const { return _parent; }
    
    TMatrix< value_t > *       next_in_block_row () const { return _next_in_block_row; }
    TMatrix< value_t > *       prev_in_block_row () const { return _prev_in_block_row; }
    TMatrix< value_t > *       next_in_block_col () const { return _next_in_block_col; }
    TMatrix< value_t > *       prev_in_block_col () const { return _prev_in_block_col; }

    TMatrix< value_t > *       next_in_block_row ( const matop_t  op ) const { return op == apply_normal ? _next_in_block_row : _next_in_block_col; }
    TMatrix< value_t > *       prev_in_block_row ( const matop_t  op ) const { return op == apply_normal ? _prev_in_block_row : _prev_in_block_col; }
    TMatrix< value_t > *       next_in_block_col ( const matop_t  op ) const { return op == apply_normal ? _next_in_block_col : _next_in_block_row; }
    TMatrix< value_t > *       prev_in_block_col ( const matop_t  op ) const { return op == apply_normal ? _prev_in_block_col : _prev_in_block_row; }

    TMatrix< value_t > *       row_diag () const { return _row_diag; }
    TMatrix< value_t > *       col_diag () const { return _col_diag; }
    
    void  set_parent            ( TBlockMatrix< value_t > *  M )
    {
        _parent = M;
    }
    
    void  set_next_in_block_row ( TMatrix< value_t > * M ) { _next_in_block_row = M; }
    void  set_prev_in_block_row ( TMatrix< value_t > * M ) { _prev_in_block_row = M; }
    void  set_next_in_block_col ( TMatrix< value_t > * M ) { _next_in_block_col = M; }
    void  set_prev_in_block_col ( TMatrix< value_t > * M ) { _prev_in_block_col = M; }

    void  set_block_row_neighbours ( TMatrix< value_t > *  prev,
                                     TMatrix< value_t > *  next )
    {
        _prev_in_block_row = prev;
        _next_in_block_row = next;
    }
    
    void  set_block_col_neighbours ( TMatrix< value_t > *  prev,
                                     TMatrix< value_t > *  next )
    {
        _prev_in_block_col = prev;
        _next_in_block_col = next;
    }

    void  set_row_diag             ( TMatrix< value_t > *  diag ) { _row_diag = diag; }
    void  set_col_diag             ( TMatrix< value_t > *  diag ) { _col_diag = diag; }

    // set hierarchy data automatically
    void  set_hierarchy_data       ();
    
    ///////////////////////////////////////////
    //
    // management of update accumulator
    //

    //! access accumulator object
    TUpdateAccumulator< value_t > &        accumulator  ()       { return _accumulator; }
    const TUpdateAccumulator< value_t > &  accumulator  () const { return _accumulator; }

    //! add update matrix
    void          add_update             ( const TMatrix< value_t > *  M,
                                           const TTruncAcc &           acc );
    
    //! add update U to set of recursive pending updates
    void          add_pending_direct     ( TDirectMatrixUpdate< value_t > *     U );

    //! add update U to set of recursive pending updates
    void          add_pending_recursive  ( TRecursiveMatrixUpdate< value_t > *  U );
        
    //! apply stored updates U to local matrix M, e.g., M = M + U,
    //! with accuracy \a acc
    virtual void  apply_updates          ( const TTruncAcc &       acc,
                                           const recursion_type_t  rec_type );

    //! return true, if matrix has updates not yet applied
    virtual bool  has_updates            ( const recursion_type_t  recursion ) const;

    //! return true, if parent matrix has updates not yet applied
    virtual bool  has_parent_updates     ( const recursion_type_t  recursion ) const;

    ///////////////////////////////////////////
    //
    // single coefficient access
    //

    //! return index (\a i,\a j) of the matrix (real valued)
    //! - \a i and \a j are global indices, e.g. not w.r.t. local index set
    virtual value_t        entry  ( const idx_t  i, const idx_t  j ) const;

    ///////////////////////////////////////////
    //
    // access cluster
    //

    //! return corresponding block cluster of matrix
    const TBlockCluster * cluster () const { return _cluster; }

    //! set block cluster of matrix
    virtual void set_cluster ( const TBlockCluster * c );

    //! set block cluster of matrix (with forced setting of cluster variable)
    virtual void set_cluster_force ( const TBlockCluster * c )
    {
        set_cluster( c );

        _cluster = c;
    }

    /////////////////////////////////////////////////
    //
    // linear operator algebra
    //

    //!
    //! mapping function of linear operator \f$A\f$, e.g. \f$ y := A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply      ( const TVector< value_t > *  x,
                               TVector< value_t > *        y,
                               const matop_t               op = apply_normal ) const
    {
        mul_vec( value_t(1), x, value_t(0), y, op );
    }

    //!
    //! mapping function with update: \f$ y := y + \alpha A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply_add  ( const value_t               alpha,
                               const TVector< value_t > *  x,
                               TVector< value_t > *        y,
                               const matop_t               op = apply_normal ) const
    {
        mul_vec( alpha, x, value_t(1), y, op );
    }

    virtual void  apply_add   ( const value_t               alpha,
                                const TMatrix< value_t > *  X,
                                TMatrix< value_t > *        Y,
                                const matop_t               op = apply_normal ) const;

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
    
    ///////////////////////////////////////////////////////////
    //
    // access to vector space elements
    //

    //! return dimension of domain
    virtual size_t  domain_dim () const { return ncols(); }
    
    //! return dimension of range
    virtual size_t  range_dim  () const { return nrows(); }
    
    //! return vector in domain space
    virtual auto    domain_vector  () const -> std::unique_ptr< TVector< value_t > > { return col_vector(); }

    //! return vector in range space
    virtual auto    range_vector   () const -> std::unique_ptr< TVector< value_t > > { return row_vector(); }
    
    /////////////////////////////////////////////////
    //
    // algebra routines
    //

    //! transpose matrix
    virtual void transpose ();
    
    //! conjugate matrix coefficients
    virtual void conjugate ();
    
    //! truncate matrix to accuracy \a acc
    virtual void truncate  ( const TTruncAcc & acc ) = 0;

    //! compute this ≔ α·this
    virtual void scale     ( const value_t  alpha );
    
    //! compute this ≔ this + α · matrix
    virtual void add       ( const value_t               alpha,
                             const TMatrix< value_t > *  matrix );

    //! compute y ≔ β·y + α·op(M)·x, with M = this
    virtual void mul_vec   ( const value_t               alpha,
                             const TVector< value_t > *  x,
                             const value_t               beta,
                             TVector< value_t >       *  y,
                             const matop_t               op = apply_normal ) const;
    
    //! compute α·op(A)·op(B), with A = this
    virtual TMatrix< value_t > *  mul_right ( const value_t               alpha,
                                              const TMatrix< value_t > *  B,
                                              const matop_t               op_A,
                                              const matop_t               op_B ) const;
    
    //! compute α·op(A)·op(B), with B = this
    virtual TMatrix< value_t > *  mul_left  ( const value_t               alpha,
                                              const TMatrix< value_t > *  A,
                                              const matop_t               op_A,
                                              const matop_t               op_B ) const;

    /////////////////////////////////////////////////
    //
    // misc.
    //

    //
    // memory consumption
    //
    
    //! return size in bytes used by this object
    virtual size_t byte_size () const;

    //! return size of (floating point) data in bytes handled by this object
    virtual size_t data_byte_size () const { return 0; }

    //! return size in bytes used by this distributed object,
    //! i.e. of all distributed sub matrices
    virtual size_t global_byte_size () const;

    //
    // RTTI
    //

    HPRO_RTTI_DERIVED( TMatrix, TLinearOperator< value_t > )

    //
    // virtual constructor
    //

    //! return matrix of same class (but no content)
    virtual auto  create       () const -> std::unique_ptr< TMatrix< value_t > > = 0;

    //! return copy of matrix
    virtual auto  copy         () const -> std::unique_ptr< TMatrix< value_t > >;

    //! return copy of matrix with accuracy \a acc and optional coarsening
    virtual auto  copy         ( const TTruncAcc & acc,
                                 const bool        coarsen = false ) const -> std::unique_ptr< TMatrix< value_t > >;

    //! return structural copy of matrix (with zeroed/empty data)
    virtual auto  copy_struct  () const -> std::unique_ptr< TMatrix< value_t > >;

    //! copy data from matrix \a A
    virtual void  copy_from    ( const TMatrix< value_t > * A );
    
    //! copy matrix into matrix \a A
    virtual void  copy_to      ( TMatrix< value_t > *         A ) const;
    
    //! copy matrix into matrix \a A with accuracy \a acc and optional coarsening
    virtual void  copy_to      ( TMatrix< value_t > *         A,
                                 const TTruncAcc &            acc,
                                 const bool                   coarsen = false ) const;

    //! return appropriate row vector object for matrix
    virtual auto row_vector    () const -> std::unique_ptr< TVector< value_t > >  
    {
        return std::make_unique< TScalarVector< value_t > >( rows(), row_ofs() );
    }
    
    //! return appropriate column vector object for matrix
    virtual auto col_vector    () const -> std::unique_ptr< TVector< value_t > >
    {
        return std::make_unique< TScalarVector< value_t > >( cols(), col_ofs() );
    }
    
    //
    // serialisation
    //

    //! read data from stream \a s and copy to matrix
    virtual void read  ( TByteStream & s );

    //! use data from stream \a s to build matrix
    virtual void build ( TByteStream & s );

    //! write data to stream \a s
    virtual void write ( TByteStream & s ) const;

    //! returns size of object in bytestream
    virtual size_t bs_size () const;

    //
    // parallel methods
    //

    //! sum up \a nparts parallel copies of matrix
    //! - if bs != NULL it will be used
    virtual void sum ( const TProcSet  & p,
                       const uint        pid,
                       const uint        nparts,
                       TByteStream     * bs,
                       const TTruncAcc & acc );

    //
    // tests
    //

    //! test data for invalid values, e.g. INF and NAN
    virtual void check_data () const;
    
    //
    // output
    //
    
    //! print basic info about matrix to stdout
    virtual void print ( const uint ofs = 0 ) const;

};

//////////////////////////////////////////////////////////
//
// functional versions
//

template < typename value_t >
value_t
entry ( const TMatrix< value_t > *  M,
        const idx_t                 i,
        const idx_t                 j )
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "entry", "M is null" );
    
    return M->entry( i, j );
}

template < typename value_t >
value_t
entry ( const TMatrix< value_t > &  M,
        const idx_t                 i,
        const idx_t                 j )
{
    return M.entry( i, j );
}

template < typename value_t,
           typename alpha_t >
void
scale ( const alpha_t         alpha,
        TMatrix< value_t > *  M )
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "scale", "M is null" );
    
    M->scale( value_t( alpha ) );
}

template < typename value_t,
           typename alpha_t >
void
scale ( const alpha_t         alpha,
        TMatrix< value_t > &  M )
{
    M.scale( value_t( alpha ) );
}

template < typename value_t >
void
mul_vec ( const value_t               alpha,
          const TMatrix< value_t > *  A,
          const TVector< value_t > *  x,
          const value_t               beta,
          TVector< value_t > *        y,
          const matop_t               op )
{
    if ( A == nullptr )
        HERROR( ERR_ARG, "mul_vec", "A is null" );
    
    A->mul_vec( alpha, x, beta, y, op );
}

template < typename value_t,
           typename alpha_t,
           typename beta_t >
void
mul_vec ( const alpha_t               alpha,
          const TMatrix< value_t > &  A,
          const TVector< value_t > &  x,
          const beta_t                beta,
          TVector< value_t > &        y,
          const matop_t               op )
{
    A.mul_vec( value_t(alpha), & x, value_t(beta), & y, op );
}

//
// copy matrices
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
copy ( const TMatrix< value_t > *  M )
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "copy", "matrix is null" );
    
    return M->copy();
}

template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
copy ( const TMatrix< value_t > &  M )
{
    return M.copy();
}

//////////////////////////////////////////////////////////
//
// extracts various parts of matrix
//

//! return vector containing diagonal coefficients of A
template < typename value_t >
std::unique_ptr< TVector< value_t > >
diagonal ( const TMatrix< value_t > *  A );

template < typename value_t >
std::unique_ptr< TVector< value_t > >
diagonal ( const TMatrix< value_t > &  A )
{
    return diagonal( &A );
}

//!
//! return copy of all diagonal blocks of A with full hierarchy
//!
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
copy_diag   ( const TMatrix< value_t > * A );

template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
copy_diag   ( const TMatrix< value_t > & A )
{
    return copy_diag( &A );
}

//!
//! return real part of matrix \a A
//!
template < typename value_t >
std::unique_ptr< TMatrix< real_type_t< value_t > > >
restrict_re  ( const TMatrix< value_t > *  A,
               const TTruncAcc &           acc );

//!
//! return imaginary part of matrix \a A
//!
template < typename value_t >
std::unique_ptr< TMatrix< real_type_t< value_t > > >
restrict_im  ( const TMatrix< value_t > *  A,
               const TTruncAcc &           acc );

//////////////////////////////////////////////////////////
//
// debug helpers
//

namespace DBG
{

//!
//! write Matrix to file
//!
template < typename value_t >
void
write ( const TMatrix< value_t > *  M,
        const std::string &         filename,
        const std::string &         matname );

template < typename value_t >
void
write ( const TMatrix< value_t > &  M,
        const std::string &         filename,
        const std::string &         matname );

}// namespace DBG

}// namespace Hpro

#endif  // __HPRO_TMATRIX_HH
