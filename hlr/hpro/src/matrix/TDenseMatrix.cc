//
// Project     : HLIBpro
// File        : TDenseMatrix.cc
// Description : class for dense matrices of arbitrary (small) size
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/error.hh"

#include "hpro/vector/TScalarVector.hh"

#include "hpro/matrix/TRkMatrix.hh"
#include "hpro/matrix/structure.hh"

#include "hpro/matrix/TDenseMatrix.hh"

namespace Hpro
{

// namespace abbr.
namespace B = BLAS;

using std::unique_ptr;
using std::make_unique;

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// local defines, functions
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

namespace
{

// defines maximal ratio of rows vs. columns before considered
// as low rank matrix
const size_t  MAX_DIM_RATIO = 30;

//
// dense MVM with optional parallelisation
//
template < typename value_t >
void
task_mulvec ( const value_t                 alpha,
              const matop_t                 op,
              const B::Matrix< value_t > &  M,
              const B::Vector< value_t > &  x,
              const value_t                 beta,
              B::Vector< value_t > &        y )
{
    B::mulvec( alpha, B::mat_view( op, M ), x, beta, y );
}

}// namespace anonymous

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// TDenseMatrix
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

//
// set size of matrix
//
template < typename value_t >
void
TDenseMatrix< value_t >::set_size ( const size_t n, const size_t m )
{
    TScopedLock  mlock( *this );
    
    if (( n != _rows ) || (m != _cols))
    {
        if ( n * m > 0 )
        {
            _mat  = std::move( B::Matrix< value_t >( n, m ) );
            _rows = n;
            _cols = m;
        }// if
        else
        {            
            _mat  = std::move( B::Matrix< value_t >() );
            _rows = 0;
            _cols = 0;
        }// else
    }// if
}

///////////////////////////////////////////////
//
// block-handling
//

//
// do α·this + β·op(M)
// (either M is sub block of this or this a sub block of M)
//
template < typename value_t >
void
TDenseMatrix< value_t >::add_block ( const value_t                    alpha,
                                     const value_t                    beta,
                                     const TDenseMatrix< value_t > *  M,
                                     const matop_t                    op )
{
    TScopedLock  mlock( *this );
    
    if ( op == apply_normal )
    {
        //
        // check if bigger or smaller
        //

        if ( this->block_is().is_subset( M->block_is() ) )
        {
            //
            // M is sub matrix of this
            //
            
            B::Matrix< value_t >  Rthis( _mat,
                                         intersect( this->row_is(), M->row_is() ) - this->row_ofs(),
                                         intersect( this->col_is(), M->col_is() ) - this->col_ofs() );

            // this ≔ α·this
            if ( alpha != value_t(1) )
            {
                if ( alpha == value_t(0) )
                    B::fill( value_t(0), Rthis );
                else
                    B::scale( alpha, Rthis );
            }// if
            
            // this ≔ this + β·M
            if ( beta != value_t(0) )
            {
                B::add( beta, M->blas_mat(), Rthis );
            }// else
        }// if
        else if ( M->block_is().is_subset( this->block_is() ) )
        {
            //
            // this is sub matrix of M
            //
            
            // this ≔ α·this
            if ( alpha != value_t(1) )
            {
                if ( alpha == value_t(0) )
                    B::fill( value_t(0), _mat );
                else
                    B::scale( alpha, _mat );
            }// if

            // this ≔ this + β·M
            if ( beta != value_t(0) )
            {
                B::Matrix< value_t >  RM( M->blas_mat(),
                                          intersect( this->row_is(), M->row_is() ) - M->row_ofs(),
                                          intersect( this->col_is(), M->col_is() ) - M->col_ofs() );

                B::add( beta, RM, _mat );
            }// if
        }// if
        else
            HERROR( ERR_INDEXSET, "(TDenseMatrix) add_block",
                   "given matrix is neither a subblock nor a superblock" );
    }// if
    else if ( op == apply_transposed )
    {
        TBlockIndexSet  MT_bs( Hpro::transpose( M->block_is() ) );

        if ( this->block_is().is_subset( MT_bs ) )
        {
            //
            // M^T is sub matrix of this
            //
            
            B::Matrix< value_t >  Rthis( _mat,
                                         intersect( this->row_is(), MT_bs.row_is() ) - this->row_ofs(),
                                         intersect( this->col_is(), MT_bs.col_is() ) - this->col_ofs() );

            // this ≔ α·this
            if ( alpha != value_t(1) )
            {
                if ( alpha == value_t(0) )
                    B::fill( value_t(0), Rthis );
                else
                    B::scale( alpha, Rthis );
            }// if

            // this ≔ this + β·M
            if ( beta != value_t(0) )
            {
                B::add( beta, B::transposed( M->blas_mat() ), Rthis );
            }// else
        }// if
        else
            HERROR( ERR_INDEXSET, "(TDenseMatrix) add_block",
                   "given matrix is neither a subblock nor a superblock" );
    }// if
    else // if ( op == apply_adjoint )
    {
        TBlockIndexSet  MT_bs( Hpro::transpose( M->block_is() ) );

        if ( this->block_is().is_subset( MT_bs ) )
        {
            //
            // M^T is sub matrix of this
            //
            
            B::Matrix< value_t >  Rthis( _mat,
                                         intersect( this->row_is(), MT_bs.row_is() ) - this->row_ofs(),
                                         intersect( this->col_is(), MT_bs.col_is() ) - this->col_ofs() );

            // this ≔ α·this
            if ( alpha != value_t(1) )
            {
                if ( alpha == value_t(0) )
                    B::fill( value_t(0), Rthis );
                else
                    B::scale( alpha, Rthis );
            }// if

            // this ≔ this + β·M
            if ( beta != value_t(0) )
            {
                B::add( beta, B::adjoint( M->blas_mat() ), Rthis );
            }// else
        }// if
        else
            HERROR( ERR_INDEXSET, "(TDenseMatrix) add_block",
                   "given matrix is neither a subblock nor a superblock" );
    }// else
}

///////////////////////////////////////////////
//
// handle size (via cluster)
//

template < typename value_t >
void
TDenseMatrix< value_t >::set_cluster ( const TBlockCluster * c )
{
    TMatrix< value_t >::set_cluster( c );

    // set size to the size of the cluster
    if ( c != nullptr )
        set_size( c->rowcl()->size(), c->colcl()->size() );
}

/////////////////////////////////////////////////
//
// BLAS-routines
//

//
// scale matrix by constant factor
//
template < typename value_t >
void
TDenseMatrix< value_t >::scale ( const value_t  f )
{
    if ( f == value_t(1) )
        return;

    TScopedLock  mlock( *this );

    if ( f == value_t(0) )
        B::fill( value_t(0), _mat );
    else
        B::scale( f, _mat );
}

//
// matrix-vector-multiplication
//
template < typename value_t >
void
TDenseMatrix< value_t >::mul_vec ( const value_t               alpha,
                                   const TVector< value_t > *  x,
                                   const value_t               beta,
                                   TVector< value_t > *        y,
                                   const matop_t               op ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TDenseMatrix) cmul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(TDenseMatrix) cmul_vec", "y = nullptr" );
    
    if ( op == apply_normal )
    {
        if (( this->row_is() != y->is() ) || ( this->col_is() != x->is() ))
            HERROR( ERR_INDEXSET, "(TDenseMatrix) cmul_vec", "incompatible vector index set" );
    }// if
    else
    {
        if (( this->col_is() != y->is() ) || ( this->row_is() != x->is() ))
            HERROR( ERR_INDEXSET, "(TDenseMatrix) cmul_vec", "incompatible vector index set" );
    }// if
    
    if ( ! IS_TYPE( x, TScalarVector ) || ! IS_TYPE( y, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TDenseMatrix) cmul_vec",
                y->typestr() + " += TDenseMatrix * " + x->typestr() );
        
    //
    // multiply
    //

    auto  s_y = ptrcast(  y, TScalarVector< value_t > );
    auto  s_x = cptrcast( x, TScalarVector< value_t > );
                
    task_mulvec( alpha, op, blas_mat(), s_x->blas_vec(), beta, s_y->blas_vec() );
}

//
// compute this ≔ this + a·M
//
template < typename value_t >
void
TDenseMatrix< value_t >::add ( const value_t               a,
                               const TMatrix< value_t > *  M )
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "(TDenseMatrix) add", "argument is nullptr" );
    
    if ( ! IS_TYPE( M, TDenseMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TDenseMatrix) add", M->typestr() );

    const TDenseMatrix *  D = cptrcast( M, TDenseMatrix );

    //
    // some checks
    //
    
    if ( a == value_t(0) )
        return;

    if ( this->block_is() != D->block_is() )
        HERROR( ERR_MAT_SIZE, "(TDenseMatrix) cadd", "" );

    //
    // add matrix
    //

    TScopedLock  mlock( *this );
    
    B::add( a, D->blas_mat(), blas_mat() );
}

//
// matrix-matrix-multiplication: alpha op_A(this) * op_B(B)
//
template < typename value_t >
TMatrix< value_t > *
TDenseMatrix< value_t >::mul_right ( const value_t               alpha,
                                     const TMatrix< value_t > *  B,
                                     const matop_t               op_A,
                                     const matop_t               op_B ) const
{
    const bool  trans_A = (op_A != apply_normal);
    const bool  trans_B = (op_B != apply_normal);
    
    if ( B == nullptr )
        HERROR( ERR_ARG, "(TDenseMatrix) cmul_right", "argument is nullptr" );

    if ( (trans_A ? rows() : cols()) != (trans_B ? B->cols() : B->rows()) )
        HERROR( ERR_MAT_SIZE, "(TDenseMatrix) cmul_right", "" );
    
    ////////////////////////////////////////////////////////////////
    //
    // check if number of coefficients in destination matrix
    // is much larger than in source matrices, e.g. multiplication
    // of low-rank matrices as dense matrices
    //

    const auto    A = this;
    const size_t  n = ( trans_A ? A->cols() : A->rows() );
    const size_t  m = ( trans_B ? B->rows() : B->cols() );
    
    ////////////////////////////////////////////////////////////////
    //
    // build temporary matrix holding the result
    //

    auto  C = make_unique< TDenseMatrix< value_t > >();

    C->set_size( n, m );
    C->set_ofs( (trans_A ? A->col_ofs() : A->row_ofs()),
                (trans_B ? B->row_ofs() : B->col_ofs()) );
    C->scale( value_t(0) );

    ////////////////////////////////////////////////////////////////
    //
    // multiply
    //

    // check if B is also dense and speed up things
    if ( Hpro::is_dense( B ) )
    {
        const auto dB = cptrcast( B, TDenseMatrix< value_t > );

        B::prod( alpha,
                 mat_view( op_A, A->blas_mat() ),
                 mat_view( op_B, dB->blas_mat() ),
                 value_t(1), C->blas_mat() );
    }// if
    else
    {
        //
        // compute A * B for each row of A
        //

        const size_t  ma   = ( trans_A ? A->rows()    : A->cols()    );
        const idx_t   xofs = ( trans_B ? B->col_ofs() : B->row_ofs() );
        const idx_t   yofs = ( trans_B ? B->row_ofs() : B->col_ofs() );
        auto          vy   = TScalarVector< value_t >( m, yofs );

        if ( op_A == apply_normal )
        {
            auto  sx = TScalarVector< value_t >( ma, xofs );

            if ( op_B == apply_normal )
            {
                // compute C^T = B^T * A^T
                for ( uint i = 0; i < n; i++ )
                {
                    auto  A_i = A->blas_mat().row( i );
                    auto  C_i = C->blas_mat().row( i );
                        
                    B::copy( A_i, sx.blas_vec() );
                    B->mul_vec( value_t(1), & sx, value_t(0), & vy, apply_transposed );
                    B::copy( vy.blas_vec(), C_i );
                }// for
            }// if
            else if ( op_B == apply_transposed )
            {
                // compute C^T = B * A^T
                for ( uint i = 0; i < n; i++ )
                {
                    auto  A_i = A->blas_mat().row( i );
                    auto  C_i = C->blas_mat().row( i );

                    B::copy( A_i, sx.blas_vec() );
                    B->mul_vec( value_t(1), & sx, value_t(0), & vy, apply_normal );
                    B::copy( vy.blas_vec(), C_i );
                }// for
            }// if
            else if ( op_B == apply_adjoint )
            {
                // compute C^H = B * A^H
                for ( uint i = 0; i < n; i++ )
                {
                    auto  A_i = A->blas_mat().row( i );
                    auto  C_i = C->blas_mat().row( i );

                    B::copy( A_i, sx.blas_vec() );
                    B::conj( sx.blas_vec() );
                    B->mul_vec( value_t(1), & sx, value_t(0), & vy, apply_normal );
                    B::copy( vy.blas_vec(), C_i );
                }// for

                B::conj( C->blas_mat() );
            }// if
        }// if
        else if ( op_A == apply_transposed )
        {
            if ( op_B == apply_normal )
            {
                // compute C^T = B^T * A
                for ( uint i = 0; i < n; i++ )
                {
                    auto  C_i = C->blas_mat().row( i );
                    auto  vx  = A->column( i );
                        
                    B->mul_vec( value_t(1), & vx, value_t(0), & vy, apply_transposed );
                    B::copy( vy.blas_vec(), C_i );
                }// for
            }// if
            else if ( op_B == apply_transposed )
            {
                // compute C^T = B * A
                for ( uint i = 0; i < n; i++ )
                {
                    auto  C_i = C->blas_mat().row( i );
                    auto  vx  = A->column( i );
                        
                    B->mul_vec( value_t(1), & vx, value_t(0), & vy, apply_normal );
                    B::copy( vy.blas_vec(), C_i );
                }// for
            }// if
            else if ( op_B == apply_adjoint )
            {
                TScalarVector< value_t >  sx( ma, xofs );

                // compute C^H = B * conj(A)
                for ( uint i = 0; i < n; i++ )
                {
                    auto  A_i = A->blas_mat().row( i );
                    auto  C_i = C->blas_mat().row( i );

                    B::copy( A_i, sx.blas_vec() );
                    B->mul_vec( value_t(1), & sx, value_t(0), & vy, apply_normal );
                    B::copy( vy.blas_vec(), C_i );
                }// for

                B::conj( C->blas_mat() );
            }// if
        }// if
        else if ( op_A == apply_adjoint )
        {
            if ( op_B == apply_normal )
            {
                // compute C^H = B^H * A
                for ( uint i = 0; i < n; i++ )
                {
                    auto  C_i = C->blas_mat().row( i );
                    auto  vx = A->column( i );
                        
                    B->mul_vec( value_t(1), & vx, value_t(0), & vy, apply_adjoint );
                    B::copy( vy.blas_vec(), C_i );
                }// for

                B::conj( C->blas_mat() );
            }// if
            else if ( op_B == apply_transposed )
            {
                TScalarVector< value_t >  sx( ma, xofs );

                // compute C^T = B * conj(A)
                for ( uint i = 0; i < n; i++ )
                {
                    auto  A_i = A->blas_mat().row( i );
                    auto  C_i = C->blas_mat().row( i );

                    B::copy( A_i, sx.blas_vec() );
                    B->mul_vec( value_t(1), & sx, value_t(0), & vy, apply_normal );
                    B::copy( vy.blas_vec(), C_i );
                }// for
            }// if
            else if ( op_B == apply_adjoint )
            {
                // compute C^H = B * A
                for ( uint i = 0; i < n; i++ )
                {
                    auto  C_i = C->blas_mat().row( i );
                    auto  vx  = A->column( i );
                        
                    B->mul_vec( value_t(1), & vx, value_t(0), & vy, apply_normal );
                    B::copy( vy.blas_vec(), C_i );
                }// for

                B::conj( C->blas_mat() );
            }// if
        }// if

        //
        // finally, scale C by alpha
        //

        if ( alpha != value_t(1) )
            B::scale( alpha, C->blas_mat() );
    }// else

    return C.release();
}

//
// matrix-matrix-multiplication: alpha op_A(A)·op_B(this)
//
template < typename value_t >
TMatrix< value_t > *
TDenseMatrix< value_t >::mul_left ( const value_t               alpha,
                                    const TMatrix< value_t > *  A,
                                    const matop_t               op_A,
                                    const matop_t               op_B ) const
{
    const bool  trans_A = (op_A != apply_normal);
    const bool  trans_B = (op_B != apply_normal);
    
    if ( A == nullptr )
        HERROR( ERR_ARG, "(TDenseMatrix) cmul_left", "argument is nullptr" );

    if ( (trans_A ? A->rows() : A->cols()) != (trans_B ? cols() : rows()) )
        HERROR( ERR_MAT_SIZE, "(TDenseMatrix) cmul_left", "" );

    ////////////////////////////////////////////////////////////////
    //
    // check if number of coefficients in destination matrix
    // is much larger than in source matrices, e.g. multiplication
    // of low-rank matrices as dense matrices
    //

    const auto    B = this;
    const size_t  n = ( trans_A ? A->cols() : A->rows() );
    const size_t  m = ( trans_B ? B->rows() : B->cols() );
    
    ////////////////////////////////////////////////////////////////
    //
    // build temporary matrix holding the result
    //

    auto  C = make_unique< TDenseMatrix< value_t > >();

    C->set_size( n, m );
    C->set_ofs( (trans_A ? A->col_ofs() : A->row_ofs()),
                (trans_B ? B->row_ofs() : B->col_ofs()) );
    C->scale( value_t(0) );
    
    ////////////////////////////////////////////////////////////////
    //
    // multiply
    //

    // check if A is also dense and speed up things
    if ( IS_TYPE( A, TDenseMatrix ) )
    {
        const auto  dA = cptrcast( A, TDenseMatrix< value_t > );

        B::prod( alpha,
                 mat_view( op_A, dA->blas_mat() ),
                 mat_view( op_B, B->blas_mat() ),
                 value_t(1), C->blas_mat() );
    }// if
    else
    {
        //
        // multiply each column of op(B) with A and store result in C
        //

        const size_t   nb   = ( trans_B ? B->cols()    : B->rows()    );
        const idx_t    xofs = ( trans_A ? A->row_ofs() : A->col_ofs() );

        if ( op_B == apply_normal )
        {
            // compute C = op(A) * B
            for ( uint i = 0; i < m; i++ )
            {
                auto  vx = B->column( i );
                auto  vy = C->column( i );
                    
                A->mul_vec( value_t(1), & vx, value_t(0), & vy, op_A );
            }// for
        }// if
        else if ( op_B == apply_transposed )
        {
            TScalarVector< value_t > sx( nb, xofs );
                
            // compute C = op(A) * B^T
            for ( uint i = 0; i < m; i++ )
            {
                auto  B_i = B->blas_mat().row( i );
                auto  vy  = C->column( i );
                   
                B::copy( B_i, sx.blas_vec() );
                A->mul_vec( value_t(1), & sx, value_t(0), & vy, op_A );
            }// for
        }// if
        else if ( op_B == apply_adjoint )
        {
            TScalarVector< value_t > sx( nb, xofs );
                
            // compute C = op(A) * B^H
            for ( uint i = 0; i < m; i++ )
            {
                auto  B_i = B->blas_mat().row( i );
                auto  vy  = C->column( i );
                    
                B::copy( B_i, sx.blas_vec() );
                B::conj( sx.blas_vec() );
                A->mul_vec( value_t(1), & sx, value_t(0), & vy, op_A );
            }// for
        }// if

        //
        // finally, scale C by alpha
        //

        if ( alpha != value_t(1) )
        {
            B::scale( alpha, C->blas_mat() );
        }// if
    }// else

    return C.release();
}

///////////////////////////////////////////////////////////
//
// linear operator mapping
//

//! same as above but only the dimension of the vector spaces is tested,
//! not the corresponding index sets
template < typename value_t >
void
TDenseMatrix< value_t >::apply_add   ( const value_t                    alpha,
                                       const BLAS::Vector< value_t > &  x,
                                       BLAS::Vector< value_t > &        y,
                                       const matop_t                    op ) const
{
    HASSERT( ( x.length() == ncols( op ) ) && ( y.length() == nrows( op ) ),
             ERR_ARG, "(TDenseMatrix) apply_add", "incompatible vector dimensions" );

    BLAS::mulvec( alpha, mat_view( op, Hpro::blas_mat( this ) ), x, value_t(1), y );
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
TDenseMatrix< value_t >::transpose ()
{
    {
        TScopedLock  mlock( *this );
    
        B::Matrix< value_t >  T( cols(), rows() );
            
        B::copy( B::transposed( blas_mat() ), T );
        blas_mat() = std::move( T );
        
        std::swap( _rows, _cols );
    }

    TMatrix< value_t >::transpose();
}
    
//
// conjugate matrix coefficients
//
template < typename value_t >
void
TDenseMatrix< value_t >::conjugate ()
{
    if ( is_complex_type< value_t >::value )
    {
        TScopedLock  mlock( *this );
        
        B::conj( blas_mat() );
    }// if
}

//
// virtual constructor
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TDenseMatrix< value_t >::copy () const
{
    auto  M = TMatrix< value_t >::copy();
    auto  D = ptrcast( M.get(), TDenseMatrix< value_t > );

    B::copy( blas_mat(), D->blas_mat() );

    return M;
}

//
// copy matrix into A
//
template < typename value_t >
void
TDenseMatrix< value_t >::copy_to ( TMatrix< value_t > * A ) const
{
    TMatrix< value_t >::copy_to( A );
    
    if ( ! IS_TYPE( A, TDenseMatrix ) )
    {
        HERROR( ERR_NOT_IMPL, "", "" );
        // convert( this, A, acc_exact );
        // return;
    }// if
    
    auto  D = ptrcast( A, TDenseMatrix< value_t > );

    D->set_size( rows(), cols() );

    B::copy( blas_mat(), D->blas_mat() );
}

//
// return structural copy
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TDenseMatrix< value_t >::copy_struct () const
{
    return TMatrix< value_t >::copy_struct();
}

//
// return size in bytes used by this object
//
template < typename value_t >
size_t
TDenseMatrix< value_t >::byte_size () const
{
    size_t  size = TMatrix< value_t >::byte_size() + 2*sizeof(size_t) + sizeof(B::Matrix< value_t >);

    size += _mat.byte_size();

    return size;
}

//
// serialisation
//
template < typename value_t >
void
TDenseMatrix< value_t >::read  ( TByteStream & s )
{
    TMatrix< value_t >::read( s );

    size_t  n, m;
    
    s.get( n );
    s.get( m );

    set_size( n, m );

    s.get( blas_mat().data(), sizeof( value_t ) * _rows * _cols );
}

template < typename value_t >
void
TDenseMatrix< value_t >::build ( TByteStream & s )
{
    TMatrix< value_t >::build( s );

    size_t  n, m;
    
    s.get( n );
    s.get( m );

    set_size( n, m );

    s.get( blas_mat().data(), sizeof( value_t ) * _rows * _cols );
}

template < typename value_t >
void
TDenseMatrix< value_t >::write ( TByteStream & s ) const
{
    TMatrix< value_t >::write( s );

    s.put( _rows );
    s.put( _cols );

    s.put( blas_mat().data(), sizeof( value_t ) * _rows * _cols );
}

//
// returns size of object in bytestream
//
template < typename value_t >
size_t
TDenseMatrix< value_t >::bs_size () const
{
    return (TMatrix< value_t >::bs_size() + sizeof(_rows) + sizeof(_cols) +
            (_rows * _cols * sizeof(value_t) ) );
}

///////////////////////////////////////////////////////////////////////////
//
// matrix permutation functions
//
///////////////////////////////////////////////////////////////////////////

namespace
{

//
// copy row/column <j> to row/column <i> in matrix M
//
template <typename T>
void
copy ( B::Matrix< T > & M,
       const bool       sort_rows,
       const idx_t      i,
       const idx_t      j )
{
    if ( sort_rows )
    {
        B::Vector< T >  row_i( M.row( i ) );
        B::Vector< T >  row_j( M.row( j ) );
            
        B::copy( row_j, row_i );
    }// if
    else
    {
        B::Vector< T >  col_i( M.column( i ) );
        B::Vector< T >  col_j( M.column( j ) );
            
        B::copy( col_j, col_i );
    }// else
}

//
// copy row/column <i> to vector <tmp>
//
template <typename T>
void
copy ( B::Matrix< T > &  M,
       const bool        sort_rows,
       const idx_t       i,
       B::Vector< T > &  tmp )
{
    if ( sort_rows )
    {
        B::Vector< T >  row_i( M.row( i ) );
            
        B::copy( row_i, tmp );
    }// if
    else
    {
        B::Vector< T >  col_i( M.column( i ) );
            
        B::copy( col_i, tmp );
    }// else
}

//
// copy vector <tmp> to row/column <i>
//
template <typename T>
void
copy ( B::Vector< T > &  tmp,
       const bool        sort_rows,
       const idx_t       i,
       B::Matrix< T > &  M )
{
    if ( sort_rows )
    {
        B::Vector< T >  row_i( M.row( i ) );
            
        B::copy( tmp, row_i );
    }// if
    else
    {
        B::Vector< T >  col_i( M.column( i ) );
            
        B::copy( tmp, col_i );
    }// else
}

//
// swap the rows/columns <i> and <j> in matrix M
//
template <typename T>
void
swap ( B::Matrix< T > &  M,
       const bool        sort_rows,
       const idx_t       i,
       const idx_t       j,
       B::Vector< T > &  tmp )
{
    if ( sort_rows )
    {
        B::Vector< T >  row_i( M.row( i ) );
        B::Vector< T >  row_j( M.row( j ) );

        B::copy( row_i, tmp );
        B::copy( row_j, row_i );
        B::copy( tmp,   row_j );
    }// if
    else
    {
        B::Vector< T >  col_i( M.column( i ) );
        B::Vector< T >  col_j( M.column( j ) );

        B::copy( col_i, tmp );
        B::copy( col_j, col_i );
        B::copy( tmp,   col_j );
    }// else
}

//
// sort rows/columns of matrix according to given permutation
//
template <typename T>
void
matrix_sort( B::Matrix< T > &  M,
             const bool        sort_rows,
             TPermutation &    perm,
             const idx_t       lb,
             const idx_t       ub,
             B::Vector< T > &  tmp )
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
TDenseMatrix< value_t >::permute ( const TPermutation *  row_perm,
                                   const TPermutation *  col_perm )
{
    BLAS::permute( blas_mat(), *row_perm, *col_perm );
    
// #if 1

//     //
//     // apply permutations by sorting the rows/columns
//     // simultaneously to sorting permutations, thereby
//     // reverting order
//     //
    
//     TPermutation   perm;

//     B::Vector< value_t >  tmp;

//     if ( row_perm != nullptr )
//     {
//         tmp  = B::Vector< value_t >( cols() );
//         perm = * row_perm;
//         matrix_sort( blas_mat(), true, perm, 0, idx_t(rows())-1, tmp );
//     }// if
    
//     if ( col_perm != nullptr )
//     {
//         tmp  = B::Vector< value_t >( rows() );
//         perm = * col_perm;
//         matrix_sort( blas_mat(), false, perm, 0, idx_t(cols())-1, tmp );
//     }// if
        
// #else
    
//     unique_ptr< TDenseMatrix >  T( ptrcast( copy().release(), TDenseMatrix ) );
//     const uint                  n = rows();
//     const uint                  m = cols();
    
//     //
//     // apply row permutation first by copying the old rows
//     // to the new position defined by the permutation
//     //

//     if ( row_perm != nullptr )
//     {
//         for ( uint i = 0; i < n; i++ )
//         {
//             const uint pi = row_perm->permute( i );
            
//             B::Vector< value_t >  this_i( blas_mat().row( i ) );
//             B::Vector< value_t >  T_pi( T->blas_mat().row( pi ) );
            
//             B::copy( this_i, T_pi );
//         }// for

//         B::copy( T->blas_mat(), blas_mat() );
//     }// if

//     //
//     // now apply column permutation in the same way
//     //

//     if ( col_perm != nullptr )
//     {
//         for ( uint j = 0; j < m; j++ )
//         {
//             const uint pj = col_perm->permute( j );

//             B::Vector< value_t >  this_j( blas_mat().column( j ) );
//             B::Vector< value_t >  T_pj( T->blas_mat().column( pj ) );
            
//             B::copy( this_j, T_pj );
//         }// for

//         B::copy( T->blas_mat(), blas_mat() );
//     }// if
    
// #endif
}

//
// test data for invalid values, e.g. INF and NAN
//
template < typename value_t >
void
TDenseMatrix< value_t >::check_data () const
{
    _mat.check_data();
}

template class TDenseMatrix< float >;
template class TDenseMatrix< double >;
template class TDenseMatrix< std::complex< float > >;
template class TDenseMatrix< std::complex< double > >;

}// namespace
