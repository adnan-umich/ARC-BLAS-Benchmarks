//
// Project     : HLIBpro
// File        : TRkMatrix.cc
// Description : class for rank-k-matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/blas/Algebra.hh"

#include "hpro/vector/TBlockVector.hh"
#include "hpro/vector/convert.hh"

#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/structure.hh"

#include "hpro/matrix/TRkMatrix.hh"

namespace Hpro
{

namespace B = BLAS;

using std::unique_ptr;
using std::make_unique;

namespace
{

template < typename value_t >
void
task_mulvec ( const value_t                 alpha,
              const matop_t                 op,
              const B::Matrix< value_t > &  A,
              const B::Matrix< value_t > &  B,
              const B::Vector< value_t > &  x,
              B::Vector< value_t > &        y )
{
    const bool  test_dense = false;

    B::Vector< value_t >  res_exact;
    
    if ( test_dense )
    {
        const auto  D = B::prod( value_t(1), A, B::adjoint( B ) );

        res_exact = B::mulvec( alpha, B::mat_view( op, D ), x );
        B::add( value_t(1), y, res_exact );
    }// if

    // B::Matrix< value_t >  CA( A, copy_value );
    // B::Matrix< value_t >  CB( B, copy_value );
    // B::Vector< value_t >  cx( x, copy_value );
    // B::Vector< value_t >  cy( y, copy_value );
    
    const auto            rank = A.ncols();
    B::Vector< value_t >  v( rank );
                
    if ( op == apply_normal )
    {
        B::mulvec( alpha, B::adjoint( B ), x, value_t(0), v );
        B::mulvec( value_t(1), A, v, value_t(1), y );
    }// if
    else if ( is_complex_type< value_t >::value && ( op == apply_transposed ))
    {
        //
        // (AB^H)^T x = conj(B) A^T x
        //   1. v = A^T x
        //   2. t = B conj(v)
        //   3. y = y + conj(t)
        //
        
        B::mulvec( alpha, B::transposed( A ), x, value_t(0), v );
        B::conj( v );

        B::Vector< value_t >  t( y.length() );
            
        B::mulvec( value_t(1), B, v, value_t(0), t );
        B::conj( t );
        B::add( value_t(1), t, y );
    }// else
    else // if ( op == apply_adjoint ) or value_t == real and op == transposed/adjoint
    {
        B::mulvec( alpha, B::adjoint( A ), x, value_t(0), v );
        B::mulvec( value_t(1), B, v, value_t(1), y );
    }// else

    if ( test_dense )
    {
        B::add( value_t(-1), y, res_exact );

        const auto  f = B::norm2( res_exact );
        
        if ( f > 1e-10 )
        {
            LOG::printf( "%.8e", f );
        }// if
    }// if
}

}// namespace anonymous

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// TRkMatrix
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////
//
// constructor and destructor
//

template < typename value_t >
TRkMatrix< value_t >::TRkMatrix ()
        : _rows( 0 )
        , _cols( 0 )
        , _rank( 0 )
{}

template < typename value_t >
TRkMatrix< value_t >::TRkMatrix ( const size_t  r,
                                  const size_t  c )
        : _rows( r )
        , _cols( c )
        , _rank( 0 )
{}

template < typename value_t >
TRkMatrix< value_t >::TRkMatrix ( const TIndexSet &   arow_is,
                                  const TIndexSet &   acol_is )
        : TMatrix< value_t >()
        , _rows( 0 )
        , _cols( 0 )
        , _rank( 0 )
{
    this->set_block_is( TBlockIndexSet( arow_is, acol_is ) );
}

template < typename value_t >
TRkMatrix< value_t >::TRkMatrix ( const TIndexSet &             arow_is,
                                  const TIndexSet &             acol_is,
                                  const B::Matrix< value_t > &  A,
                                  const B::Matrix< value_t > &  B )
        : TMatrix< value_t >()
        , _mat_A( A )
        , _mat_B( B )
        , _rows( arow_is.size() )
        , _cols( acol_is.size() )
        , _rank( _mat_A.ncols() )
{
    if ( A.ncols() != B.ncols() )
        HERROR( ERR_ARG, "(TRkMatrix) ctor", "rank of A and B differs" );

    // do not call set_block_is to avoid size initialisation
    this->set_ofs( arow_is.first(), acol_is.first() );
}

template < typename value_t >
TRkMatrix< value_t >::TRkMatrix ( const TIndexSet &        arow_is,
                                  const TIndexSet &        acol_is,
                                  B::Matrix< value_t > &&  A,
                                  B::Matrix< value_t > &&  B )
        : TMatrix< value_t >()
        , _mat_A( std::move( A ) )
        , _mat_B( std::move( B ) )
        , _rows( arow_is.size() )
        , _cols( acol_is.size() )
        , _rank( _mat_A.ncols() )
{
    if ( A.ncols() != B.ncols() )
        HERROR( ERR_ARG, "(TRkMatrix) ctor", "rank of A and B differs" );

    // do not call set_block_is to avoid size initialisation
    this->set_ofs( arow_is.first(), acol_is.first() );
}

template < typename value_t >
TRkMatrix< value_t >::TRkMatrix ( const TBlockIndexSet &  ablock_is )
        : TMatrix< value_t >()
        , _rows( 0 )
        , _cols( 0 )
        , _rank( 0 )
{
    this->set_block_is( ablock_is );
}

template < typename value_t >
TRkMatrix< value_t >::TRkMatrix ( const TBlockCluster * bcl )
        : TMatrix< value_t >( bcl )
        , _rows( 0 )
        , _cols( 0 )
        , _rank( 0 )
{
    if ( bcl != nullptr )
    {
        set_size( bcl->rowcl()->size(), bcl->colcl()->size() );
    }// if
}

template < typename value_t >
TRkMatrix< value_t >::TRkMatrix ( const TRkMatrix &  A )
        : TMatrix< value_t >()
        , _rows( 0 )
        , _cols( 0 )
        , _rank( 0 )
{
    set_cluster( A.cluster() );
    set_size( A.rows(), A.cols(), A.rank() );
    this->set_ofs( A.row_ofs(), A.col_ofs() );
    
    B::copy( A._mat_A, _mat_A );
    B::copy( A._mat_B, _mat_B );
}

//
// set rank of matrix
//

template < typename value_t >
void
TRkMatrix< value_t >::set_rank ( const size_t  k )
{
    if ( k == _rank )
        return;

    // allocate memory for new arrays
    set_size( _rows, _cols, k );
}

//
// usual matrix access (read-only)
//
template < typename value_t >
value_t
TRkMatrix< value_t >::entry ( const idx_t i, const idx_t j ) const
{
    value_t  f = value_t(0);
        
    for ( uint k = 0; k < _rank; k++ )
        f += _mat_A( i, k ) * Math::conj( _mat_B( j, k ) );
    
    return f;
}

//
// update size of matrices if cluster changed
//
template < typename value_t >
void
TRkMatrix< value_t >::set_cluster ( const TBlockCluster * c )
{
    TMatrix< value_t >::set_cluster( c );

    if ( c != nullptr )
        set_size( c->rowcl()->size(), c->colcl()->size(), _rank );
}

/////////////////////////////////////////////////
//
// management of update accumulator
//

//
// apply stored updates U to local matrix M, e.g., M = M + U,
// with accuracy \a acc
//
template < typename value_t >
void
TRkMatrix< value_t >::apply_updates ( const TTruncAcc &       /* acc */,
                                      const recursion_type_t )
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

/////////////////////////////////////////////////
//
// truncate the rank via bestapproximation in frobenius-norm
//

template < typename value_t >
void
TRkMatrix< value_t >::truncate ( const TTruncAcc & acc )
{
    if ( _rank == 0 )
        return;
        
    TScopedLock  mlock( *this );

    _rank = B::truncate( _mat_A, _mat_B, acc( this ) );

    if (( _mat_A.ncols() != _rank ) || ( _mat_B.ncols() != _rank ))
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) truncate", "" );
}

//
// copy given dense matrix into local rank-matrix
//
template < typename value_t >
void
TRkMatrix< value_t >::copy_dense ( const TDenseMatrix< value_t > *  A,
                                   const TTruncAcc &                acc )
{
    if ((A->rows() != _rows) || (A->cols() != _cols))
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) copy_dense", "argument has wrong dimensions" );

    TScopedLock   mlock( *this );
    const size_t  n = A->rows();
    const size_t  m = A->cols();

    if ( acc.is_exact() )
    {
        //
        // copy dense matrix directly
        //

        if ( A->rows() > A->cols() )
        {
            const size_t          rk = A->cols();
            B::Matrix< value_t >  TA( rows(), rk );
            B::Matrix< value_t >  TB( cols(), rk );
                
            // copy A to R_A
            B::copy( A->blas_mat(), TA );
                
            // and set R_B to e_1 ... e_rk
            B::fill( value_t(0), TB );
            for ( uint i = 0; i < rk; i++ )
                TB( i, i ) = value_t(1);

            _mat_A = std::move( TA );
            _mat_B = std::move( TB );
            _rank   = rk;
        }// if
        else
        {
            const size_t          rk = A->rows();
            B::Matrix< value_t >  TA( rows(), rk );
            B::Matrix< value_t >  TB( cols(), rk );
                
            // copy A^T to R_B; conjugate since B^H == A
            B::copy( B::adjoint( A->blas_mat() ), TB );
                
            // and set R_A to e_1 ... e_rk
            B::fill( value_t(0), TA );
            for ( uint i = 0; i < rk; i++ )
                TA( i, i ) = value_t(1);

            _mat_A = std::move( TA );
            _mat_B = std::move( TB );
            _rank   = rk;
        }// else
    }// if
    else
    {
        //
        // compute low-rank approximation of dense matrix
        //
            
        B::Matrix< value_t > D( n, m );
        
        B::copy( A->blas_mat(), D );
        _rank = B::approx( D, acc( this ), _mat_A, _mat_B );
    }// else

    if (( _mat_A.ncols() != _rank ) || ( _mat_B.ncols() != _rank ))
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) copy_dense", "" );
}

//
// set this ≔ A·B^H
//
template < typename value_t >
void
TRkMatrix< value_t >::set_lrmat ( const BLAS::Matrix< value_t > &  A,
                                  const BLAS::Matrix< value_t > &  B )
{
    if (( A.nrows() != _mat_A.nrows() ) ||
        ( B.nrows() != _mat_B.nrows() ) ||
        ( A.ncols() != B.ncols()))
        HERROR( ERR_ARG, "(TRkMatrix) set_lrmat", "input matrices have invalid dimension" );

    if ( _rank == A.ncols() )
    {
        B::copy( A, _mat_A );
        B::copy( B, _mat_B );
    }// if
    else
    {
        _mat_A = std::move( B::Matrix< value_t >( A, copy_value ) );
        _mat_B = std::move( B::Matrix< value_t >( B, copy_value ) );
        _rank  = A.ncols();
    }// else
}

//
// add a rank-k-matrix and truncate result
//
template < typename value_t >
void
TRkMatrix< value_t >::add_rank ( const value_t                 alpha,
                                 const B::Matrix< value_t > &  A,
                                 const B::Matrix< value_t > &  B,
                                 const TTruncAcc &             acc )
{
    if ( A.ncols() != B.ncols() )
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) add_rank", "given matrices have different ncols" );

    const size_t  k = A.ncols();
    
    // check trivial case
    if (( k == 0 ) || ( alpha == value_t(0) ))
        return;
    
    //
    // build new matrix holding the sum
    //

    TScopedLock           mlock( *this );
    B::Matrix< value_t >  TA( _rows, _rank + k );
    B::Matrix< value_t >  TB( _cols, _rank + k );

    // copy this
    if ( _rank > 0 )
    {
        const B::Range        old_vecs( 0, idx_t(_rank)-1 );
        B::Matrix< value_t >  partA( TA, B::Range( 0, idx_t(_rows)-1 ), old_vecs );
        B::Matrix< value_t >  partB( TB, B::Range( 0, idx_t(_cols)-1 ), old_vecs );

        B::copy( _mat_A, partA );
        B::copy( _mat_B, partB );
    }// if

    // copy given matrix
    if ( k > 0 )
    {
        const B::Range        new_vecs( idx_t(_rank), idx_t(_rank + k) - 1 );
        B::Matrix< value_t >  partA( TA, B::Range( 0, idx_t(_rows)-1 ), new_vecs );
        B::Matrix< value_t >  partB( TB, B::Range( 0, idx_t(_cols)-1 ), new_vecs );

        B::copy( A, partA );
        B::copy( B, partB );

        if ( alpha != value_t(1) )
        {
            if ( A.nrows() < B.nrows() ) B::scale( alpha,             partA );
            else                         B::scale( Math::conj(alpha), partB );
        }// if
    }// if

    //
    // truncate result
    //

    if ( CFG::Arith::zero_sum_trunc || ( _rank > 0 ))
        _rank = B::truncate( TA, TB, acc( this ) );
    else
        _rank = k;

    _mat_A = std::move( TA );
    _mat_B = std::move( TB );
}

//
// add a dense matrix and truncate
//
template < typename value_t >
void
TRkMatrix< value_t >::add_dense ( const value_t                 alpha,
                                  const B::Matrix< value_t > &  D,
                                  const TTruncAcc &             acc )
{
    if ( alpha == value_t(0) )
        return;
    
    TScopedLock           mlock( *this );
    B::Matrix< value_t >  TD( _rows, _cols );
    
    B::copy( D, TD );

    if ( alpha != value_t(1) )
        B::scale( alpha, TD );
    
    // add this to D
    if ( _rank > 0 )
        B::prod( value_t(1), _mat_A, adjoint(_mat_B), value_t(1), TD );
    
    // truncate and store result in this
    _rank = B::approx( TD, acc( this ), _mat_A, _mat_B );
}

/////////////////////////////////////////////////
//
// BLAS-routines (real valued)
//

//
// scale matrix by constant factor
//
template < typename value_t >
void
TRkMatrix< value_t >::scale ( const value_t  f )
{
    if (( f == value_t(1) ) || ( rank() == 0 ))
        return;

    TScopedLock  mlock( *this );

    if ( f == value_t(0) )
    {
        B::fill( value_t(0), blas_mat_A() );
        B::fill( value_t(0), blas_mat_B() );
    }// if
    else
    {
        //
        // scale matrix with lower norm
        //

        if ( B::normF( blas_mat_A() ) < B::normF( blas_mat_B() ) )
            B::scale( f, blas_mat_A() );
        else
            B::scale( f, blas_mat_B() );
    }// else
}
    
//
// matrix-vector-multiplication
//
template < typename value_t >
void
TRkMatrix< value_t >::mul_vec ( const value_t               alpha,
                                const TVector< value_t > *  x,
                                const value_t               beta,
                                TVector< value_t > *        y,
                                const matop_t               op ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TRkMatrix) mul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(TRkMatrix) mul_vec", "y = nullptr" );
    
    if ( op == apply_normal )
    {
        if (( this->row_is() != y->is() ) || ( this->col_is() != x->is() ))
            HERROR( ERR_INDEXSET, "(TRkMatrix) mul_vec", "incompatible vector index set" );
    }// if
    else
    {
        if (( this->col_is() != y->is() ) || ( this->row_is() != x->is() ))
            HERROR( ERR_INDEXSET, "(TRkMatrix) mul_vec", "incompatible vector index set" );
    }// if
    
    //
    // check if scalar vector
    //
    
    if ( IS_TYPE( x, TScalarVector ) && IS_TYPE( y, TScalarVector ) )
    {
        auto *  sy =  ptrcast( y, TScalarVector< value_t > );
        auto *  sx = cptrcast( x, TScalarVector< value_t > );

        //
        // we assume, that the index-sets are coherent
        // and we can access by local indices
        //

        // scale left-side
        if      ( beta == value_t(0) ) sy->fill( value_t(0) );
        else if ( beta != value_t(1) ) sy->scale( beta );

        if (( alpha == value_t(0) ) || ( _rank == 0 ))
            return;

        task_mulvec< value_t >( alpha, op, blas_mat_A(), blas_mat_B(), sx->blas_vec(), sy->blas_vec() );
    }// if
    else
        HERROR( ERR_VEC_TYPE, "(TRkMatrix) cmul_vec", y->typestr() + " += TRkMatrix * " + x->typestr() );
}

//
// compute this ≔ this + α · M
//
template < typename value_t >
void
TRkMatrix< value_t >::add ( const value_t               alpha,
                            const TMatrix< value_t > *  M )
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "(TRkMatrix) add", "nullptr argument" );
    
    if ( ! IS_TYPE( M, TRkMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TRkMatrix) add", M->typestr() );

    auto  R = cptrcast( M, TRkMatrix );

    //
    // some trivial checks
    //

    if (( alpha == value_t(0) ) || ( R->rank() == 0 ))
        return;

    if ((rows() != R->rows()) || (cols() != R->cols()))
        HERROR( ERR_ARG, "(TRkMatrix) cadd", "argument has wrong dimensions" );

    //
    // In general, when adding two rank-k-matrices, we get a rank-2k-matrix.
    // we do not truncate the rank here (just in case, you want something
    // different), so you have to do this in some other place.
    //

    const B::Range  new_vecs( idx_t(_rank), idx_t(_rank + R->rank()) - 1 );

    set_rank( _rank + R->rank() );

    B::Matrix< value_t >  A_new( _mat_A, B::Range( 0, idx_t(_rows)-1 ), new_vecs );
    B::Matrix< value_t >  B_new( _mat_B, B::Range( 0, idx_t(_cols)-1 ), new_vecs );
        
    // copy new vectors
    B::copy( R->blas_mat_A(), A_new );
    B::copy( R->blas_mat_B(), B_new );
        
    // scale vectors (choose smaller ones)
    if ( alpha != value_t(1) )
    {
        if (_rows < _cols) B::scale( alpha, A_new );
        else               B::scale( alpha, B_new );
    }// if
}

//
// matrix-matrix-multiplication
//
template < typename value_t >
TRkMatrix< value_t > *
TRkMatrix< value_t >::mul_right ( const value_t               alpha,
                                  const TMatrix< value_t > *  B,
                                  const matop_t               op_A,
                                  const matop_t               op_B ) const
{
    const bool  trans_A = (op_A != apply_normal);
    const bool  trans_B = (op_B != apply_normal);
    
    if ( B == nullptr )
        HERROR( ERR_ARG, "(TRkMatrix) mul_right", "argument is nullptr" );

    if ( (trans_A ? rows() : cols()) != (trans_B ? B->cols() : B->rows()) )
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) mul_right", "incompatible dimensions" );
    
    ////////////////////////////////////////////////////////////////
    //
    // build temporary matrix holding the result
    //

    const TRkMatrix *  A = this;
    const size_t       n = ( trans_A ? A->cols() : A->rows() );
    const size_t       m = ( trans_B ? B->rows() : B->cols() );
    const size_t       k = A->rank();
    auto               C = make_unique< TRkMatrix >();

    C->set_size( n, m, k );
    C->set_ofs( (trans_A ? A->col_ofs() : A->row_ofs()),
                (trans_B ? B->row_ofs() : B->col_ofs()) );

    if ( k == 0 )
        return C.release();
    
    ////////////////////////////////////////////////////////////////
    //
    // multiply
    //

    //
    // depending on op_A, copy either A- or B-vectors to A(C)
    //
    
    if ( op_A == apply_normal )
        B::copy( A->blas_mat_A(), C->blas_mat_A() );
    else if ( op_A == apply_transposed )
    {
        B::copy( A->blas_mat_B(), C->blas_mat_A() );
        B::conj( C->blas_mat_A() );
    }// if
    else
        B::copy( A->blas_mat_B(), C->blas_mat_A() );

    //
    // compute B-vectors of C
    //
    
    const size_t    ma   = ( trans_A ? A->rows()    : A->cols()    );
    const idx_t     xofs = ( trans_B ? B->col_ofs() : B->row_ofs() );
    
    if ( op_A == apply_normal )
    {
        if ( op_B == apply_normal )
        {
            // B(C)^H = B(A)^H * B => B(C) = B^H * B(A)
            for ( uint i = 0; i < k; i++ )
            {
                auto  vx = A->vec_B( i );
                auto  vy = C->vec_B( i );
                    
                B->mul_vec( value_t(1), & vx, value_t(0), & vy, apply_adjoint );
            }// for
        }// if
        else if ( op_B == apply_transposed )
        {
            TScalarVector< value_t >  sx( ma, xofs );
            
            // B(C)^H = B(A)^H * B^T => B(C) = conj(B) * B(A) = conj( B * conj(B(A)) )
            for ( uint i = 0; i < k; i++ )
            {
                auto  B_i = A->blas_mat_B().column( i );
                auto  vy  = C->vec_B( i );
                    
                B::copy( B_i, sx.blas_vec() );
                B::conj( sx.blas_vec() );
                B->mul_vec( value_t(1), & sx, value_t(0), & vy, apply_normal );
                B::conj( vy.blas_vec() );
            }// for
        }// if
        else if ( op_B == apply_adjoint )
        {
            // B(C)^H = B(A)^H * B^H => B(C) = B * B(A)
            for ( uint i = 0; i < k; i++ )
            {
                auto  vx = A->vec_B( i );
                auto  vy = C->vec_B( i );
                    
                B->mul_vec( value_t(1), & vx, value_t(0), & vy, apply_normal );
            }// for
        }// if
    }// if
    else if ( op_A == apply_transposed )
    {
        if ( op_B == apply_normal )
        {
            // B(C)^H = A(A)^T * B => B(C) = B^H * conj(A(A)) = conj(B^T * A)
            for ( uint i = 0; i < k; i++ )
            {
                auto  vx = A->vec_A( i );
                auto  vy = C->vec_B( i );
                    
                B->mul_vec( value_t(1), & vx, value_t(0), & vy, apply_transposed );
            }// for

            B::conj( C->blas_mat_B() );
        }// if
        else if ( op_B == apply_transposed )
        {
            // B(C)^H = A(A)^T * B^T => B(C) = conj(B) * conj(A(A)) = conj(B * A(A))
            for ( uint i = 0; i < k; i++ )
            {
                auto  vx = A->vec_A( i );
                auto  vy = C->vec_B( i );
                    
                B->mul_vec( value_t(1), & vx, value_t(0), & vy, apply_normal );
            }// for

            B::conj( C->blas_mat_B() );
        }// if
        else if ( op_B == apply_adjoint )
        {
            TScalarVector< value_t >  sx( ma, xofs );
                
            // B(C)^H = A(A)^T * B^H => B(C) = B * conj(A(A))
            for ( uint i = 0; i < k; i++ )
            {
                auto  A_i = A->blas_mat_A().column( i );
                auto  vy  = C->vec_B( i );
                    
                B::copy( A_i, sx.blas_vec() );
                B::conj( sx.blas_vec() );
                B->mul_vec( value_t(1), & sx, value_t(0), & vy, apply_normal );
            }// for
        }// if
    }// if
    else if ( op_A == apply_adjoint )
    {
        if ( op_B == apply_normal )
        {
            // B(C)^H = A(A)^H * B => B(C) = B^H * A(A)
            for ( uint i = 0; i < k; i++ )
            {
                auto  vx = A->vec_A( i );
                auto  vy = C->vec_B( i );
                    
                B->mul_vec( value_t(1), & vx, value_t(0), & vy, apply_adjoint );
            }// for
        }// if
        else if ( op_B == apply_transposed )
        {
            TScalarVector< value_t >  sx( ma, xofs );
                
            // B(C)^H = A(A)^H * B^T => B(C) = conj(B) * A(A) = conj(B * conj(A(A)))
            for ( uint i = 0; i < k; i++ )
            {
                auto  A_i = A->blas_mat_A().column( i );
                auto  vy  = C->vec_B( i );
                    
                B::copy( A_i, sx.blas_vec() );
                B::conj( sx.blas_vec() );
                B->mul_vec( value_t(1), & sx, value_t(0), & vy, apply_normal );
                B::conj( vy.blas_vec() );
            }// for
        }// if
        else if ( op_B == apply_adjoint )
        {
            // B(C)^H = A(A)^H * B^H => B(C) = B * A(A)
            for ( uint i = 0; i < k; i++ )
            {
                auto  vx = A->vec_A( i );
                auto  vy = C->vec_B( i );
                    
                B->mul_vec( value_t(1), & vx, value_t(0), & vy, apply_normal );
            }// for
        }// if
    }// if

    //
    // finally, scale C by alpha
    //
    
    if ( alpha != value_t(1) )
    {
        if ( n > m ) B::scale( Math::conj(alpha), C->blas_mat_B() );
        else         B::scale( alpha,             C->blas_mat_A() );
    }// if

    return C.release();
}

template < typename value_t >
TRkMatrix< value_t > *
TRkMatrix< value_t >::mul_left ( const value_t               alpha,
                                 const TMatrix< value_t > *  A,
                                 const matop_t               op_A,
                                 const matop_t               op_B ) const
{
    const bool  trans_A = (op_A != apply_normal);
    const bool  trans_B = (op_B != apply_normal);
    
    if ( A == nullptr )
        HERROR( ERR_ARG, "(TRkMatrix) cmul_left", "argument is nullptr" );

    if ( (trans_A ? A->rows() : A->cols()) != (trans_B ? cols() : rows()) )
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) cmul_left", "incompatible dimensions" );

    ////////////////////////////////////////////////////////////////
    //
    // build temporary matrix holding the result
    //

    const auto    B = this;
    const size_t  n = ( trans_A ? A->cols() : A->rows() );
    const size_t  m = ( trans_B ? B->rows() : B->cols() );
    const size_t  k = B->rank();
    auto          C = make_unique< TRkMatrix >();

    C->set_size( n, m, k );
    C->set_ofs( (trans_A ? A->col_ofs() : A->row_ofs()),
                (trans_B ? B->row_ofs() : B->col_ofs()) );
    
    if ( k == 0 )
        return C.release();
    
    ////////////////////////////////////////////////////////////////
    //
    // multiply
    //

    //
    // copy B-/A-vectors
    //
    
    if ( op_B == apply_normal )
        B::copy( B->blas_mat_B(), C->blas_mat_B() );
    else if ( op_B == apply_transposed )
    {
        B::copy( B->blas_mat_A(), C->blas_mat_B() );
        B::conj( C->blas_mat_B() );
    }// if
    else if ( op_B == apply_adjoint )
        B::copy( B->blas_mat_A(), C->blas_mat_B() );

    //
    // compute A- or B-vectors of C
    //
    
    const size_t  nb   = ( trans_B ? B->cols()    : B->rows()    );
    const idx_t   xofs = ( trans_A ? A->row_ofs() : A->col_ofs() );

    if ( op_B == apply_normal )
    {
        // multiply A * A(B)
        for ( uint i = 0; i < k; i++ )
        {
            auto  vx = B->vec_A( i );
            auto  vy = C->vec_A( i );
                
            A->mul_vec( value_t(1), & vx, value_t(0), & vy, op_A );
        }// for
    }// if
    else if ( op_B == apply_transposed )
    {
        TScalarVector< value_t >  sx( nb, xofs );
            
        // multiply A * conj(B(B))
        for ( uint i = 0; i < k; i++ )
        {
            auto  B_i = B->blas_mat_B().column( i );
            auto  vy  = C->vec_A( i );
                
            B::copy( B_i, sx.blas_vec() );
            B::conj( sx.blas_vec() );
            A->mul_vec( value_t(1), & sx, value_t(0), & vy, op_A );
        }// for
    }// if
    else if ( op_B == apply_adjoint )
    {
        // multiply A * B(B)
        for ( uint i = 0; i < k; i++ )
        {
            auto  vx = B->vec_B( i );
            auto  vy = C->vec_A( i );
                
            A->mul_vec( value_t(1), & vx, value_t(0), & vy, op_A );
        }// for
    }// if

    //
    // finally, scale C with alpha
    //
    
    if ( alpha != value_t(1) )
    {
        if ( n > m ) B::scale( Math::conj(alpha), C->blas_mat_B() );
        else         B::scale( alpha,             C->blas_mat_A() );
    }// if

    return C.release();
}
    
///////////////////////////////////////////////////////////
//
// linear operator mapping
//

template < typename value_t >
void
TRkMatrix< value_t >::apply_add   ( const value_t                    alpha,
                                    const BLAS::Vector< value_t > &  x,
                                    BLAS::Vector< value_t > &        y,
                                    const matop_t                    op ) const
{
    HASSERT( ( x.length() == ncols( op ) ) && ( y.length() == nrows( op ) ),
             ERR_ARG, "(TRkMatrix) apply_add", "incompatible vector dimensions" );

    switch ( op )
    {
        case apply_normal :
        {
            const auto  t = BLAS::mulvec( alpha, BLAS::adjoint( blas_mat_B() ), x );
                
            BLAS::mulvec( value_t(1), blas_mat_A(), t, value_t(1), y );
        }
        break;

        case apply_transposed :
        {
            const auto  t = BLAS::mulvec( value_t(1), BLAS::transposed( blas_mat_A() ), x );
            auto        s = BLAS::mulvec( value_t(1), blas_mat_B(), t );
                
            BLAS::conj( s );
            BLAS::add( alpha, s, y );
        }
        break;

        case apply_adjoint :
        {
            const auto  t = BLAS::mulvec( alpha, BLAS::adjoint( blas_mat_A() ), x );
                
            BLAS::mulvec( value_t(1), blas_mat_B(), t, value_t(1), y );
        }
        break;

        default:  // apply_conjugate
        {
            const auto  t = BLAS::mulvec( value_t(1), BLAS::transposed( blas_mat_B() ), x );
            auto        s = BLAS::mulvec( value_t(1), blas_mat_A(), t );
                
            BLAS::conj( s );
            BLAS::add( alpha, s, y );
        }
    }// switch
}

/////////////////////////////////////////////////
//
// misc methods
//

//
// transpose matrix
//
template < typename value_t >
void
TRkMatrix< value_t >::transpose ()
{
    //
    // (A·B^H)^T = B^H^T A^T = conj(B) conj(A)^H
    //
    
    conjugate();

    {
        TScopedLock  mlock( *this );

        std::swap( _mat_A, _mat_B );
        std::swap( _rows, _cols );
    }

    TMatrix< value_t >::transpose();
}
    
//
// conjugate matrix coefficients
//
template < typename value_t >
void
TRkMatrix< value_t >::conjugate ()
{
    if ( this->is_complex() && ! this->is_hermitian() )
    {
        TScopedLock  mlock( *this );
        
        B::conj( blas_mat_A() );
        B::conj( blas_mat_B() );
    }// if
}
    
//
// allocate memory (and update data)
//
template < typename value_t >
void
TRkMatrix< value_t >::set_size ( const size_t  n,
                                 const size_t  m,
                                 const size_t  new_rank )
{
    //
    // first check if we really have to do something
    //
    if (( _rows == n ) && ( _cols == m ) && ( new_rank == _rank ))
        return;

    TScopedLock  mlock( *this );

    //
    // a rank 0 means -> remove all data
    //

    if ( new_rank == 0 )
    {
        // delete old data
        _mat_A = B::Matrix< value_t >( n, 0 );
        _mat_B = B::Matrix< value_t >( m, 0 );

        _rows = n;
        _cols = m;
        _rank = 0;
    }// if
    else
    {
        const  B::Range  loc_row_is( 0, idx_t(n)-1 );
        const  B::Range  loc_col_is( 0, idx_t(m)-1 );
        const  B::Range  old_vec_is( 0, idx_t( std::min( _rank, new_rank ) ) - 1 );
            
        //
        // handle matrix A
        //

        if ( _rows == n )
        {
            // only copy if different rank
            if ( _rank != new_rank )
            {
                B::Matrix< value_t >  new_A( n, new_rank );
                B::Matrix< value_t >  sub_old_A( _mat_A, loc_row_is, old_vec_is );
                B::Matrix< value_t >  sub_new_A( new_A, loc_row_is, old_vec_is );

                B::copy( sub_old_A, sub_new_A );

                _mat_A = std::move( new_A );
            }// if
        }// if
        else
        {
            // simply create (nullified) array
            _rows   = n;
            _mat_A = B::Matrix< value_t >( n, new_rank );
        }// else

        //
        // handle matrix B
        //

        if ( _cols == m )
        {
            // only copy if different rank
            if ( _rank != new_rank )
            {
                B::Matrix< value_t >  new_B( m, new_rank );
                B::Matrix< value_t >  sub_old_B( _mat_B, loc_col_is, old_vec_is );
                B::Matrix< value_t >  sub_new_B( new_B, loc_col_is, old_vec_is );

                B::copy( sub_old_B, sub_new_B );

                _mat_B = std::move( new_B );
            }// if
        }// if
        else
        {
            // simply create (nullified) arrays
            _cols   = m;
            _mat_B = B::Matrix< value_t >( m, new_rank );
        }// else

        // adjust rank
        _rank = new_rank;
    }// else

    if ( rows() != _mat_A.nrows() )
        HERROR( ERR_CONSISTENCY, "", "" );

    if ( cols() != _mat_B.nrows() )
        HERROR( ERR_CONSISTENCY, "", "" );
}

//
// virtual constructor
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TRkMatrix< value_t >::copy () const
{
    auto  M = TMatrix< value_t >::copy();
    auto  R = ptrcast( M.get(), TRkMatrix );

    if ( this->cluster() != nullptr )
        R->set_cluster( this->cluster() );

    R->set_rank( _rank );

    if ( _rank > 0 )
    {
        B::copy( blas_mat_A(), R->blas_mat_A() );
        B::copy( blas_mat_B(), R->blas_mat_B() );
    }// if

    return M;
}

//
// copy matrix wrt. given accuracy
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TRkMatrix< value_t >::copy ( const TTruncAcc & acc, const bool ) const
{
    auto  M = copy();
    auto  R = ptrcast( M.get(), TRkMatrix );

    R->truncate( acc( this ) );
    
    return M;
}

//
// return structural copy
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TRkMatrix< value_t >::copy_struct () const
{
    auto  M = TMatrix< value_t >::copy_struct();
    auto  R = ptrcast( M.get(), TRkMatrix );

    if ( this->cluster() != nullptr )
        R->set_cluster( this->cluster() );

    R->set_rank( 0 );

    return M;
}

//
// copy matrix into A
//
template < typename value_t >
void
TRkMatrix< value_t >::copy_to ( TMatrix< value_t > * A ) const
{
    TMatrix< value_t >::copy_to( A );
    
    if ( ! IS_TYPE( A, TRkMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TRkMatrix) copy_to", A->typestr() );
    
    TRkMatrix * R = ptrcast( A, TRkMatrix );

    if (( R->rows() != rows() ) || ( R->cols() != cols() ))
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) copy_to", "matrix has wrong dimension" );

    R->set_size( rows(), cols(), _rank );
    R->set_ofs( this->row_ofs(), this->col_ofs() );

    if ( _rank > 0 )
    {
        B::copy( blas_mat_A(), R->blas_mat_A() );
        B::copy( blas_mat_B(), R->blas_mat_B() );
    }// if
}

template < typename value_t >
void
TRkMatrix< value_t >::copy_to ( TMatrix< value_t > * A, const TTruncAcc & acc, const bool ) const
{
    TMatrix< value_t >::copy_to( A );
    
    if ( ! IS_TYPE( A, TRkMatrix ) )
    {
        HERROR( ERR_NOT_IMPL, "", "" );
        // convert( this, A, acc );
        // return;
    }// if
    
    TRkMatrix * R = ptrcast( A, TRkMatrix );

    if (( R->rows() != rows() ) || ( R->cols() != cols() ))
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) copy_to", "matrix has wrong dimension" );

    R->set_size( rows(), cols(), _rank );
    R->set_ofs( this->row_ofs(), this->col_ofs() );

    if ( _rank > 0 )
    {
        B::copy( blas_mat_A(), R->blas_mat_A() );
        B::copy( blas_mat_B(), R->blas_mat_B() );
        
        R->truncate( acc( this ) );
    }// if
}

//
// return size in bytes used by this object
//
template < typename value_t >
size_t
TRkMatrix< value_t >::byte_size () const
{
    size_t  s = TMatrix< value_t >::byte_size() + 3 * sizeof(size_t);

    s += _mat_A.byte_size();
    s += _mat_B.byte_size();

    return s;
}

//
// serialisation
//

template < typename value_t >
void
TRkMatrix< value_t >::read  ( TByteStream & s )
{
    TMatrix< value_t >::read( s );

    size_t  n, m, k;
    
    s.get( n );
    s.get( m );
    s.get( k );

    set_size( n, m, k );

    s.get( _mat_A.data(), sizeof( value_t ) * _rows * _rank );
    s.get( _mat_B.data(), sizeof( value_t ) * _cols * _rank );
}

template < typename value_t >
void
TRkMatrix< value_t >::build  ( TByteStream & s )
{
    TMatrix< value_t >::build( s );

    size_t  n, m, k;
    
    s.get( n );
    s.get( m );
    s.get( k );

    set_size( n, m, k );

    s.get( _mat_A.data(), sizeof( value_t ) * _rows * _rank );
    s.get( _mat_B.data(), sizeof( value_t ) * _cols * _rank );
}

template < typename value_t >
void
TRkMatrix< value_t >::write ( TByteStream & s ) const
{
    TMatrix< value_t >::write( s );

    s.put( _rows );
    s.put( _cols );
    s.put( _rank );

    s.put( _mat_A.data(), sizeof( value_t ) * _rows * _rank );
    s.put( _mat_B.data(), sizeof( value_t ) * _cols * _rank );
}

//
// returns size of object in bytestream
//
template < typename value_t >
size_t
TRkMatrix< value_t >::bs_size () const
{
    return (TMatrix< value_t >::bs_size() + sizeof(_rows) + sizeof(_cols) + sizeof(_rank) +
            (_rank * (_rows + _cols) * sizeof(value_t)));
}

//
// test data for invalid values, e.g. INF and NAN
//
template < typename value_t >
void
TRkMatrix< value_t >::check_data () const
{
    _mat_A.check_data();
    _mat_B.check_data();
}

//
// explicit instantiation
//

template class TRkMatrix< float >;
template class TRkMatrix< double >;
template class TRkMatrix< std::complex< float > >;
template class TRkMatrix< std::complex< double > >;

}// namespace Hpro
