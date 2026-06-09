//
// Project     : HLIBpro
// File        : TBlockMatrix.cc
// Description : class for a matrix consisting of submatrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <list>
#include <iostream>

#include "hpro/vector/TBlockVector.hh"

#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"
#include "hpro/matrix/TBSHMBuilder.hh"
#include "hpro/matrix/structure.hh"

#include "hpro/matrix/TBlockMatrix.hh"

namespace Hpro
{

// namespace abbr.
namespace B = BLAS;

///////////////////////////////////////////
//
// constructor and destructor
//
template < typename value_t >
TBlockMatrix< value_t >::~TBlockMatrix ()
{
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
            if ( block(i,j) != nullptr )
                delete block(i,j);

    _blocks.resize( 0 );
}

//
// set cluster of matrix
//
template < typename value_t >
void
TBlockMatrix< value_t >::set_cluster ( const TBlockCluster * c )
{
    TMatrix< value_t >::set_cluster( c );

    if ( c != nullptr )
    {
        set_block_struct( c->nrows(), c->ncols() );
        set_size( c->rowcl()->size(), c->colcl()->size() );
    }// if
}

//
// set block-structure
//
template < typename value_t >
void
TBlockMatrix< value_t >::set_block_struct ( const uint bn, const uint bm )
{
    // adjust values
    uint n = bn;
    uint m = bm;
    
    if (( n == 0 ) && ( m != 0 )) n = 1;
    if (( n != 0 ) && ( m == 0 )) m = 1;
    
    std::vector< TMatrix< value_t > * >  tmp( n * m );
    const uint                           mrc = std::min(_block_rows*_block_cols, n*m);

    for ( uint i = 0; i < mrc; i++ )
        tmp[i] = _blocks[i];

    for ( uint i = mrc; i < _block_rows*_block_cols; i++ )
        delete _blocks[i];
    
    for ( uint i = mrc; i < n*m; i++ )
        tmp[i] = nullptr;
    
    _blocks     = std::move( tmp );
    _block_rows = n;
    _block_cols = m;
}

//
// set matrix format
//
template < typename value_t >
void
TBlockMatrix< value_t >::set_form ( const matform_t  f )
{
    TMatrix< value_t >::set_form( f );

    // also change for diagonal sub blocks
    for ( uint i = 0; i < std::min( block_rows(), block_cols() ); ++i )
    {
        if ( block( i, i ) != nullptr )
            block( i, i )->set_form( f );
    }// for
}
    
template < typename value_t >
void
TBlockMatrix< value_t >::set_form ( const matform_t         f,
                                    const recursion_type_t  t )
{
    if ( t == recursive )
    {
        set_form( f );
    }// if
    else
    {
        // only change locally
        TMatrix< value_t >::set_form( f );
    }// else
}
    
//
// replace matrix block <A> by matrix <B> (don't delete A !)
//
template < typename value_t >
void
TBlockMatrix< value_t >::replace_block ( TMatrix< value_t > * A, TMatrix< value_t > * B )
{
    if ( A == nullptr )
        HWARNING( "in (TBlockMatrix) replace_block : replacing nullptr block" );

    TScopedLock  slock( * this );
    
    for ( uint i = 0; i < block_rows(); i++ )
    {
        for ( uint j = 0; j < block_cols(); j++ )
        {
            if ( block( i, j ) == A )
            {
                set_block( i, j, B );
                return;
            }// if
        }// for
    }// for
    
    HERROR( ERR_MAT_STRUCT, "(TBlockMatrix) replace_block", "given matrix not a sub block" );
}

//
// clear pointers to all subblocks
//
template < typename value_t >
void
TBlockMatrix< value_t >::clear_blocks ()
{
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
            set_block( i, j, nullptr );
}
    
//
// access single matrix coefficient
//
template < typename value_t >
value_t
TBlockMatrix< value_t >::entry ( const idx_t  row,
                                 const idx_t  col ) const
{
    // translate index pair to lower half in case of symmetry
    idx_t  srow     = row;
    idx_t  scol     = col;
    bool   swap_idx = false;

    if ( ! this->is_nonsym() && ( srow < scol ))
    {
        std::swap( srow, scol );
        swap_idx = true;
    }// if

    const idx_t  rofs = this->row_ofs();
    const idx_t  cofs = this->col_ofs();
    
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
        {
            const TMatrix< value_t > * A_ij = block(i,j);

            if ( A_ij == nullptr )
                continue;
            
            if ( A_ij->block_is().is_in( rofs + srow, cofs + scol ) )
            {
                const auto  val = A_ij->entry( srow - A_ij->row_ofs() + rofs,
                                               scol - A_ij->col_ofs() + cofs );
                
                if ( swap_idx && this->is_hermitian() )
                    return Math::conj( val );
                else
                    return val;
            }// if
        }// for

    return value_t(0);
}

//
// set processor set of matrix
//
template < typename value_t >
void
TBlockMatrix< value_t >::set_procs ( const TProcSet &        ps,
                                     const recursion_type_t  t )
{
    TMatrix< value_t >::set_procs( ps, nonrecursive );

    if ( t == recursive )
    {
        for ( uint i = 0; i < block_rows(); i++ )
            for ( uint j = 0; j < block_cols(); j++ )
            {
                TMatrix< value_t > * A_ij = block(i,j);

                if ( A_ij != nullptr )
                    A_ij->set_procs( ps, t );
            }// for
    }// if
}

/////////////////////////////////////////////////
//
// update accumulator
//
template < typename value_t >
void
TBlockMatrix< value_t >::apply_updates ( const TTruncAcc &       /* acc */,
                                         const recursion_type_t  /* recursion */ )
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// return true, if matrix has updates not yet applied
//
template < typename value_t >
bool
TBlockMatrix< value_t >::has_updates ( const recursion_type_t  recursion ) const
{
    if ( TMatrix< value_t >::has_updates( recursion ) )
        return true;

    if ( recursion )
    {
        for ( uint  i = 0; i < block_rows(); ++i )
            for ( uint  j = 0; j < block_cols(); ++j )
            {
                if (( block( i, j ) != nullptr ) && block( i, j )->has_updates( recursion ) )
                    return true;
            }// for
    }// if

    return false;
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
TBlockMatrix< value_t >::scale ( const value_t  f )
{
    if ( f == value_t(1) )
        return;
    
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
        {
            TMatrix< value_t >  * A_ij = block(i,j);
            
            if ( A_ij != nullptr )
                A_ij->scale( f );
        }// for
}
    
//
// matrix-vector-mult.
//
template < typename value_t >
void
TBlockMatrix< value_t >::mul_vec ( const value_t               alpha,
                                   const TVector< value_t > *  x,
                                   const value_t               beta,
                                   TVector< value_t > *        y,
                                   const matop_t               op ) const
{
    // if ( ! is_small( this ) || ( procs().size() > 1 ))
    // {
    //     Hpro::mul_vec( procs(), alpha, this, x, beta, y, op );
    //     return;
    // }// if
    
    if ( x == nullptr ) HERROR( ERR_ARG, "(TBlockMatrix) mul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(TBlockMatrix) mul_vec", "y = nullptr" );

    if (( this->row_is( op ) != y->is() ) || ( this->col_is( op ) != x->is() ))
        HERROR( ERR_INDEXSET, "(TBlockMatrix) mul_vec", "incompatible vector index set" );
    
    //
    // decide upon type
    //

    if ( IS_TYPE( x, TScalarVector ) && IS_TYPE( y, TScalarVector ) )
    {
        auto  sx = cptrcast( x, TScalarVector< value_t > );
        auto  sy =  ptrcast( y, TScalarVector< value_t > );

        //
        // get pointer to actual data by comparing local indices and indices
        // in blockmatrix
        //
        
        if ( beta != value_t(1) )
            y->scale( beta );
        
        if ( alpha == value_t(0) )
            return;
        
        if ( op == apply_normal )
        {
            const matop_t  sym_op = this->is_symmetric() ? apply_transposed : apply_adjoint;
            
            for ( uint i = 0; i < block_rows(); i++ )
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    const auto  A_ij = block(i,j);
                    
                    if ( A_ij == nullptr )
                        continue;

                    {
                        // multiply with sub block
                        auto  tx = sub_vector( sx, A_ij->col_is() );
                        auto  ty = sub_vector( sy, A_ij->row_is() );
                
                        A_ij->mul_vec( alpha, & tx, value_t(1), & ty, op );
                    }

                    // in case of symmetry, the same again but transposed
                    if ( ! this->is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        auto  tx = sub_vector( sx, A_ij->row_is() );
                        auto  ty = sub_vector( sy, A_ij->col_is() );

                        A_ij->mul_vec( alpha, & tx, value_t(1), & ty, sym_op );
                    }// if
                }// for
        }// if
        else if ( op == apply_transposed )
        {
            for ( uint i = 0; i < block_rows(); i++ )
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    const auto  A_ij = block(i,j);
                    
                    if ( A_ij == nullptr )
                        continue;

                    {
                        // multiply with sub block
                        auto  tx = sub_vector( sx, A_ij->row_is() );
                        auto  ty = sub_vector( sy, A_ij->col_is() );
                    
                        A_ij->mul_vec( alpha, & tx, value_t(1), & ty, op );
                    }

                    // in case of symmetry, the same again but transposed
                    if ( ! this->is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        auto  tx = sub_vector( sx, A_ij->col_is() );
                        auto  ty = sub_vector( sy, A_ij->row_is() );
                        
                        if ( A_ij->is_complex() && this->is_hermitian() )
                        {
                            // multiply with conjugated but not transposed matrix B
                            // via: y = conj(B) x  ⇔  conj(y) = B conj(x)
                            TScalarVector< value_t >  bx( A_ij->col_is() );
                            TScalarVector< value_t >  by( A_ij->row_is() );

                            B::copy( tx.blas_vec(), bx.blas_vec() );
                            B::conj( bx.blas_vec() );
                                            
                            A_ij->mul_vec( alpha, & bx, value_t(1), & by, apply_normal );

                            B::conj( by.blas_vec() );
                            B::add( value_t(1), by.blas_vec(), ty.blas_vec() );
                        }// if
                        else
                        {
                            A_ij->mul_vec( alpha, & tx, value_t(1), & ty, apply_normal );
                        }// else
                    }// if
                }// for
        }// if
        else if ( op == apply_adjoint )
        {
            for ( uint i = 0; i < block_rows(); i++ )
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    const auto  A_ij = block(i,j);
                    
                    if ( A_ij == nullptr )
                        continue;

                    {
                        // multiply with sub block
                        auto  tx = sub_vector( sx, A_ij->row_is() );
                        auto  ty = sub_vector( sy, A_ij->col_is() );
                    
                        A_ij->mul_vec( alpha, & tx, value_t(1), & ty, op );
                    }

                    // in case of symmetry, the same again but transposed
                    if ( ! this->is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        auto  tx = sub_vector( sx, A_ij->col_is() );
                        auto  ty = sub_vector( sy, A_ij->row_is() );
                        
                        if ( this->is_symmetric() && A_ij->is_complex() )
                        {
                            // multiply with conjugated but not transposed matrix B
                            // via: y = conj(B) x <=> conj(y) = B conj(x)
                            TScalarVector< value_t >  bx( A_ij->col_is() );
                            TScalarVector< value_t >  by( A_ij->row_is() );

                            B::copy( tx.blas_vec(), bx.blas_vec() );
                            B::conj( bx.blas_vec() );
                                            
                            A_ij->mul_vec( alpha, & bx, value_t(1), & by, apply_normal );

                            B::conj( by.blas_vec() );
                            B::add( value_t(1), by.blas_vec(), ty.blas_vec() );
                        }// if
                        else
                        {
                            A_ij->mul_vec( alpha, & tx, value_t(1), & ty, apply_normal );
                        }// else
                    }// if
                }// for
        }// if
    }// if
    else if ( IS_TYPE( x, TBlockVector ) && IS_TYPE( y, TBlockVector ) )
    {
        auto bx = cptrcast( x, TBlockVector< value_t > );
        auto by =  ptrcast( y, TBlockVector< value_t > );

        //
        // get pointer to actual data by comparing local indices and indices
        // in blockmatrix
        //
        
        if ( beta != value_t(1) )
            y->scale( beta );
        
        if ( alpha == value_t(0) )
            return;
        
        if ( op == apply_normal )
        {
            const matop_t  sym_op = this->is_symmetric() ? apply_transposed : apply_adjoint;
            
            for ( uint i = 0; i < block_rows(); i++ )
            {
                auto  y_i = by->block( i );
                
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    auto  A_ij = block(i,j);
                    auto  x_j  = bx->block( j );
                    
                    if ( A_ij == nullptr )
                        continue;
                    
                    A_ij->mul_vec( alpha, x_j, value_t(1), y_i, op );

                    // in case of symmetry, the same again but transposed
                    if ( ! this->is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        A_ij->mul_vec( alpha, bx->block( i ), value_t(1), by->block( j ), sym_op );
                    }// if
                }// for
            }// for
        }// if
        else if ( op == apply_transposed )
        {
            for ( uint i = 0; i < block_rows(); i++ )
            {
                auto  x_i = bx->block( i );
                
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    auto  A_ij = block(i,j);
                    auto  y_j  = by->block( j );
                    
                    if ( A_ij == nullptr )
                        continue;
                    
                    A_ij->mul_vec( alpha, x_i, value_t(1), y_j, op );

                    // in case of symmetry, the same again but transposed
                    if ( ! this->is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        if ( A_ij->is_complex() && this->is_hermitian() )
                        {
                            // multiply with conjugated but not transposed matrix B
                            // via: y = conj(B) x <=> conj(y) = B conj(x)
                            HERROR( ERR_NOT_IMPL, "", "" );
                        }// if
                        else
                        {
                            A_ij->mul_vec( alpha, bx->block(j), value_t(1), by->block(i), apply_normal );
                        }// else
                    }// if
                }// for
            }// for
        }// if
        else if ( op == apply_adjoint )
        {
            for ( uint i = 0; i < block_rows(); i++ )
            {
                auto  x_i = bx->block( i );
                
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    auto  A_ij = block(i,j);
                    auto  y_j  = by->block( j );
                    
                    if ( A_ij == nullptr )
                        continue;
                    
                    A_ij->mul_vec( alpha, x_i, value_t(1), y_j, op );

                    // in case of symmetry, the same again but transposed
                    if ( ! this->is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        if ( this->is_symmetric() && A_ij->is_complex() )
                        {
                            // multiply with conjugated but not transposed matrix B
                            // via: y = conj(B) x <=> conj(y) = B conj(x)
                            HERROR( ERR_NOT_IMPL, "", "" );
                        }// if
                        else
                        {
                            A_ij->mul_vec( alpha, bx->block(j), value_t(1), by->block(i), apply_normal );
                        }// else
                    }// if
                }// for
            }// for
        }// if
    }// if
    else
        HERROR( ERR_VEC_TYPE, "(TBlockMatrix) mul_vec",
                y->typestr() + " = TBlockMatrix * " + x->typestr() );
}

//
// compute this = this + alpha * matrix
//
template < typename value_t >
void
TBlockMatrix< value_t >::add ( const value_t               /* alpha */,
                               const TMatrix< value_t > *  /* B */ )
{
    HERROR( ERR_NOT_IMPL, "", "" );
    // return Hpro::add( alpha, B, value_t(1), this, TTruncAcc( Limits::epsilon<real_t>() ) );
}

///////////////////////////////////////////////////////////
//
// linear operator mapping
//

//
// same as above but only the dimension of the vector spaces is tested,
// not the corresponding index sets
//
template < typename value_t >
void
TBlockMatrix< value_t >::apply_add  ( const value_t                    alpha,
                                      const BLAS::Vector< value_t > &  x,
                                      BLAS::Vector< value_t > &        y,
                                      const matop_t                    op ) const
{
    if ( alpha == value_t(0) )
        return;
        
    if ( op == apply_normal )
    {
        const matop_t  sym_op = this->is_symmetric() ? apply_transposed : apply_adjoint;
            
        for ( uint i = 0; i < block_rows(); i++ )
        {
            for ( uint j = 0; j < block_cols(); j++ )
            {
                const auto  A_ij = block(i,j);
                    
                if ( A_ij == nullptr )
                    continue;

                {
                    // multiply with sub block
                    const B::Vector< value_t >  tx( x, A_ij->col_is() - this->col_ofs() );
                    B::Vector< value_t >        ty( y, A_ij->row_is() - this->row_ofs() );
                
                    A_ij->apply_add( alpha, tx, ty, op );
                }

                // in case of symmetry, the same again but transposed
                if ( ! this->is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                {
                    const B::Vector< value_t >  tx( x, A_ij->row_is() - this->row_ofs() );
                    B::Vector< value_t >        ty( y, A_ij->col_is() - this->col_ofs() );

                    A_ij->apply_add( alpha, tx, ty, sym_op );
                }// if
            }// for
        }// for
    }// if
    else if ( op == apply_transposed )
    {
        for ( uint i = 0; i < block_rows(); i++ )
        {
            for ( uint j = 0; j < block_cols(); j++ )
            {
                const auto  A_ij = block(i,j);
                    
                if ( A_ij == nullptr )
                    continue;

                {
                    // multiply with sub block
                    const B::Vector< value_t >  tx( x, A_ij->row_is() - this->row_ofs() );
                    B::Vector< value_t >        ty( y, A_ij->col_is() - this->col_ofs() );
                    
                    A_ij->apply_add( alpha, tx, ty, op );
                }

                // in case of symmetry, the same again but transposed
                if ( ! this->is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                {
                    const B::Vector< value_t >  tx( x, A_ij->col_is() - this->col_ofs() );
                    B::Vector< value_t >        ty( y, A_ij->row_is() - this->row_ofs() );
                        
                    if ( this->is_hermitian() )
                    {
                        // multiply with conjugated but not transposed matrix B
                        // via: y = conj(B) x  ⇔  conj(y) = B conj(x)
                        HERROR( ERR_NOT_IMPL, "", "" );
                    }// if
                    else
                    {
                        A_ij->apply_add( alpha, tx, ty, apply_normal );
                    }// else
                }// if
            }// for
        }// for
    }// if
    else if ( op == apply_adjoint )
    {
        for ( uint i = 0; i < block_rows(); i++ )
            for ( uint j = 0; j < block_cols(); j++ )
            {
                const auto  A_ij = block(i,j);
                    
                if ( A_ij == nullptr )
                    continue;

                {
                    // multiply with sub block
                    const B::Vector< value_t >  tx( x, A_ij->row_is() - this->row_ofs() );
                    B::Vector< value_t >        ty( y, A_ij->col_is() - this->col_ofs() );
                    
                    A_ij->apply_add( alpha, tx, ty, op );
                }

                // in case of symmetry, the same again but transposed
                if ( ! this->is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr ))
                {
                    const B::Vector< value_t >  tx( x, A_ij->col_is() - this->col_ofs() );
                    B::Vector< value_t >        ty( y, A_ij->row_is() - this->row_ofs() );
                        
                    if ( this->is_symmetric() )
                    {
                        // multiply with conjugated but not transposed matrix B
                        // via: y = conj(B) x <=> conj(y) = B conj(x)
                        HERROR( ERR_NOT_IMPL, "", "" );
                    }// if
                    else
                    {
                        A_ij->apply_add( alpha, tx, ty, apply_normal );
                    }// else
                }// if
            }// for
    }// if
}

template < typename value_t >
void
TBlockMatrix< value_t >::apply_add  ( const value_t                    alpha,
                                      const BLAS::Matrix< value_t > &  X,
                                      BLAS::Matrix< value_t > &        Y,
                                      const matop_t                    op ) const
{
    if ( alpha == value_t(0) )
        return;
        
    if ( op == apply_normal )
    {
        const matop_t  sym_op = this->is_symmetric() ? apply_transposed : apply_adjoint;
            
        for ( uint i = 0; i < nblock_rows(); i++ )
        {
            for ( uint j = 0; j < nblock_cols(); j++ )
            {
                const auto  A_ij = block(i,j);
                    
                if ( A_ij == nullptr )
                    continue;

                {
                    // multiply with sub block
                    auto  TX = B::Matrix< value_t >( X, A_ij->col_is() - this->col_ofs(), B::Range::all );
                    auto  TY = B::Matrix< value_t >( Y, A_ij->row_is() - this->row_ofs(), B::Range::all );
                
                    A_ij->apply_add( alpha, TX, TY, op );
                }

                // in case of symmetry, the same again but transposed
                if ( ! this->is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                {
                    auto  TX = B::Matrix< value_t >( X, A_ij->row_is() - this->row_ofs(), B::Range::all );
                    auto  TY = B::Matrix< value_t >( Y, A_ij->col_is() - this->col_ofs(), B::Range::all );

                    A_ij->apply_add( alpha, TX, TY, sym_op );
                }// if
            }// for
        }// for
    }// if
    else if ( op == apply_transposed )
    {
        for ( uint i = 0; i < nblock_rows(); i++ )
        {
            for ( uint j = 0; j < nblock_cols(); j++ )
            {
                const auto  A_ij = block(i,j);
                    
                if ( A_ij == nullptr )
                    continue;

                {
                    // multiply with sub block
                    auto  TX = B::Matrix< value_t >( X, A_ij->row_is() - this->row_ofs(), B::Range::all );
                    auto  TY = B::Matrix< value_t >( Y, A_ij->col_is() - this->col_ofs(), B::Range::all );
                    
                    A_ij->apply_add( alpha, TX, TY, op );
                }

                // in case of symmetry, the same again but transposed
                if ( ! this->is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                {
                    auto  TX = B::Matrix< value_t >( X, A_ij->col_is() - this->col_ofs(), B::Range::all );
                    auto  TY = B::Matrix< value_t >( Y, A_ij->row_is() - this->row_ofs(), B::Range::all );
                        
                    if ( this->is_hermitian() )
                    {
                        // multiply with conjugated but not transposed matrix B
                        // via: y = conj(B) x  ⇔  conj(y) = B conj(x)
                        HERROR( ERR_NOT_IMPL, "", "" );
                    }// if
                    else
                    {
                        A_ij->apply_add( alpha, TX, TY, apply_normal );
                    }// else
                }// if
            }// for
        }// for
    }// if
    else if ( op == apply_adjoint )
    {
        for ( uint i = 0; i < block_rows(); i++ )
            for ( uint j = 0; j < block_cols(); j++ )
            {
                const auto  A_ij = block(i,j);
                    
                if ( A_ij == nullptr )
                    continue;

                {
                    // multiply with sub block
                    auto  TX = B::Matrix< value_t >( X, A_ij->row_is() - this->row_ofs(), B::Range::all );
                    auto  TY = B::Matrix< value_t >( Y, A_ij->col_is() - this->col_ofs(), B::Range::all );
                    
                    A_ij->apply_add( alpha, TX, TY, op );
                }

                // in case of symmetry, the same again but transposed
                if ( ! this->is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr ))
                {
                    auto  TX = B::Matrix< value_t >( X, A_ij->col_is() - this->col_ofs(), B::Range::all );
                    auto  TY = B::Matrix< value_t >( Y, A_ij->row_is() - this->row_ofs(), B::Range::all );
                        
                    if ( this->is_symmetric() )
                    {
                        // multiply with conjugated but not transposed matrix B
                        // via: y = conj(B) x <=> conj(y) = B conj(x)
                        HERROR( ERR_NOT_IMPL, "", "" );
                    }// if
                    else
                    {
                        A_ij->apply_add( alpha, TX, TY, apply_normal );
                    }// else
                }// if
            }// for
    }// if
}

///////////////////////////////////////////
//
// block-handling
//

//
// return subblock of matrix corresponding to t
//
template < typename value_t >
TMatrix< value_t > *
TBlockMatrix< value_t >::bc_block ( const TBlockCluster * t ) const
{
    if ( t == nullptr )
        HERROR( ERR_ARG, "(TBlockMatrix) bc_block", "given block cluster is nullptr" );
    
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
        {
            const TMatrix< value_t >  * A_ij = block(i,j);
            
            if (( A_ij != nullptr ) && ( A_ij->cluster() == t ))
                return const_cast< TMatrix< value_t > * >( A_ij );
        }// for
    
    return nullptr;
}

template < typename value_t >
TMatrix< value_t > *
TBlockMatrix< value_t >::bc_block ( const TCluster * tau, const TCluster * sigma ) const
{
    if (( tau == nullptr ) || ( sigma == nullptr ))
        HERROR( ERR_ARG, "(TBlockMatrix) bc_block", "given clusters are nullptr" );
    
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
        {
            const TMatrix< value_t >  * A_ij = block(i,j);
            
            if ( ( A_ij != nullptr ) && ( A_ij->cluster() != nullptr ) &&
                 ( A_ij->cluster()->rowcl() == tau) &&
                 ( A_ij->cluster()->colcl() == sigma ) )
                return const_cast< TMatrix< value_t > * >( A_ij );
        }// for
    
    return nullptr;
}

//
// truncate all rank-blocks in matrix to given accuracy
//
template < typename value_t >
void
TBlockMatrix< value_t >::truncate ( const TTruncAcc & acc )
{
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
        {
            TMatrix< value_t >  * A_ij = block(i,j);
            
            if ( A_ij == nullptr )
                continue;

            A_ij->truncate( acc );
        }// for
}
    
//
// transpose matrix
//
template < typename value_t >
void
TBlockMatrix< value_t >::transpose ()
{
    if ( this->is_symmetric() )
    {
        // nothing to do, because transposed matrix is the same
        return;
    }// if
    else if ( this->is_hermitian() )
    {
        // A = A^H ⇒ A^T = conj(A)
        conjugate();
    }// if
    else
    {
        //
        // transpose sub matrices and transpose internal
        // pointer matrix
        //

        std::vector< TMatrix< value_t > * >  T( block_cols() * block_rows() );

        for ( uint  j = 0; j < block_cols(); ++j )
            for ( uint  i = 0; i < block_rows(); ++i )
            {
                if ( block( i, j ) != nullptr )
                    block( i, j )->transpose();

                T[ i * block_cols() + j ] = block( i, j );
            }// for
        
        for ( uint  i = 0; i < block_cols() * block_rows(); ++i )
            _blocks[ i ] = T[ i ];

        std::swap( _block_rows, _block_cols );
        std::swap( _rows, _cols );
    }// else

    TMatrix< value_t >::transpose();
}

//
// conjugate matrix coefficients
//
template < typename value_t >
void
TBlockMatrix< value_t >::conjugate ()
{
    if ( this->is_hermitian() )
    {
        return;
    }// if
    else
    {
        for ( uint  j = 0; j < block_cols(); ++j )
            for ( uint  i = 0; i < block_rows(); ++i )
            {
                if ( block( i, j ) != nullptr )
                    block( i, j )->conjugate();
            }// for
    }// else
}

//
// output of matrix
//
template < typename value_t >
void
TBlockMatrix< value_t >::print ( const uint ofs ) const
{
    for ( uint i = 0; i < ofs; i++ )
        std::cout << ' ';

    std::cout << this->typestr()
              << " ( " << rows() << " x " << cols()
              << ", +" << this->row_ofs() << "+" << this->col_ofs()
              << ", " << block_rows() << " x " << block_cols()
              << " )"
              << std::endl;
    
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
            if ( block(i,j) != nullptr )
            {
                block(i,j)->print( ofs + 4 );
                if (( i < block_rows()-1 ) || ( j < block_cols()-1 ))
                    std::cout << std::endl;
            }// if
}

//
// virtual constructors
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TBlockMatrix< value_t >::copy () const
{
    auto  M = TMatrix< value_t >::copy();
    auto  B = ptrcast( M.get(), TBlockMatrix< value_t > );

    B->set_block_struct( block_rows(), block_cols() );

    //
    // copy each of the blocks and add to new matrix
    //

    for ( uint  i = 0; i < block_rows(); ++i )
    {
        for ( uint  j = 0; j < block_cols(); ++j )
        {
            const TMatrix< value_t > *  B_ij = block(i,j);
            
            if ( B_ij != nullptr )
            {
                auto  T = B_ij->copy();
                
                T->set_parent( B );
                B->set_block( i, j, T.release() );
            }// if
        }// for
    }// for
        
    return M;
}

//
// copy matrix wrt. given accuracy and coarsening
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TBlockMatrix< value_t >::copy ( const TTruncAcc & acc, const bool coarsen ) const
{
    auto  M = TMatrix< value_t >::copy();
    auto  B = ptrcast( M.get(), TBlockMatrix< value_t > );

    B->set_block_struct( block_rows(), block_cols() );

    //
    // copy each of the blocks and add to new matrix
    //

    for ( uint  i = 0; i < block_rows(); ++i )
    {
        for ( uint  j = 0; j < block_cols(); ++j )
        {
            const TMatrix< value_t > *  M_ij = block(i,j);

            if ( M_ij == nullptr )
            {
                B->set_block( i, j, nullptr );
            }// if
            else
            {
                std::unique_ptr< TMatrix< value_t > >  T_ij( M_ij->copy( acc( M_ij ), coarsen ) );

                // if ( coarsen )
                // {
                //     TCoarsen  mcoarsen;
                            
                //     T_ij = unique_ptr< TMatrix< value_t > >( mcoarsen.coarsen( T_ij.release(), acc( M_ij ) ) );
                // }// if
                    
                B->set_block( i, j, T_ij.release() );
            }// else
        }// for
    }// for
    
    return M;
}

//
// return structural copy
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TBlockMatrix< value_t >::copy_struct () const
{
    auto  M = TMatrix< value_t >::copy_struct();
    auto  B = ptrcast( M.get(), TBlockMatrix< value_t > );

    B->set_block_struct( block_rows(), block_cols() );

    //
    // copy each of the blocks and add to new matrix
    //

    for ( uint  i = 0; i < block_rows(); ++i )
    {
        for ( uint  j = 0; j < block_cols(); ++j )
        {
            const TMatrix< value_t > *  B_ij = block(i,j);

            if ( B_ij != nullptr )
            {
                auto  T = B_ij->copy_struct();

                T->set_parent( B );
                B->set_block( i, j, T.release() );
            }// if
        }// for
    }// for
        
    return M;
}

//
// copy matrix into A
//
template < typename value_t >
void
TBlockMatrix< value_t >::copy_to ( TMatrix< value_t > * A ) const
{
    TMatrix< value_t >::copy_to( A );
    
    if ( ! IS_TYPE( A, TBlockMatrix ) )
    {
        HERROR( ERR_NOT_IMPL, "", "" );
        // convert( this, A, acc_exact );
        // return;
    }// if

    TBlockMatrix  * b_A = ptrcast( A, TBlockMatrix< value_t > );
    
    //
    // go over all blocks, get assoc. blocks in A and
    // copy them one by one
    //
    
    for ( uint  i = 0; i < block_rows(); ++i )
    {
        for ( uint  j = 0; j < block_cols(); ++j )
        {
            if ( block(i,j) == nullptr )
                continue;
            
            const TMatrix< value_t > * B = block(i,j);
            TMatrix< value_t >       * C = b_A->block(i,j);
            
            if ( C != nullptr )
                B->copy_to( C );
            else
                HERROR( ERR_MAT_STRUCT, "(TBlockMatrix) copy_to",
                        "could not find corresponding block in A" );
        }// for
    }// for
}

template < typename value_t >
void
TBlockMatrix< value_t >::copy_to ( TMatrix< value_t > *  A,
                                   const TTruncAcc &     acc,
                                   const bool            coarsen ) const
{
    TMatrix< value_t >::copy_to( A );
    
    if ( ! IS_TYPE( A, TBlockMatrix ) )
    {
        HERROR( ERR_NOT_IMPL, "", "" );
        // convert( this, A, acc );
        // return;
    }// if

    TBlockMatrix  * b_A = ptrcast( A, TBlockMatrix< value_t > );
    
    //
    // go over all blocks, get assoc. blocks in A and
    // copy them one by one
    //

    for ( uint  i = 0; i < block_rows(); ++i )
    {
        for ( uint  j = 0; j < block_cols(); ++j )
        {
            if ( block(i,j) == nullptr )
                continue;
            
            const TMatrix< value_t > * B_ij = block(i,j);
            TMatrix< value_t >       * C_ij = b_A->block(i,j);
            
            if ( C_ij == nullptr )
                HERROR( ERR_MAT_STRUCT, "(TBlockMatrix) copy_to",
                        "could not find corresponding block in A" );
            
            if ( IS_TYPE( B_ij, TBlockMatrix ) )
                cptrcast( B_ij, TBlockMatrix< value_t > )->copy_to( C_ij, acc( B_ij ), coarsen );
            else if ( IS_TYPE( B_ij, TRkMatrix ) )
                cptrcast( B_ij, TRkMatrix< value_t > )->copy_to( C_ij, acc( B_ij ) );
            else
                B_ij->copy_to( C_ij );

            // if ( coarsen && ( acc.rank() == 0 ))
            // {
            //     TCoarsen  mcoarsen;

            //     C_ij = mcoarsen.coarsen( C_ij, acc( B_ij ) );
            // }// if
            
            b_A->set_block( i, j, C_ij );
        }// for
    }// for
}

//
// copy complete structural information from given matrix
//
template < typename value_t >
void
TBlockMatrix< value_t >::copy_struct_from ( const TMatrix< value_t > * M )
{
    TMatrix< value_t >::copy_struct_from( M );

    if ( IS_TYPE( M, TBlockMatrix ) )
    {
        const TBlockMatrix * B = cptrcast( M, TBlockMatrix< value_t > );

        set_block_struct( B->block_rows(), B->block_cols() );
    }// if
}
    
//
// return size in bytes used by this object
//
template < typename value_t >
size_t
TBlockMatrix< value_t >::byte_size () const
{
    size_t  size = 0;
    
    size += TMatrix< value_t >::byte_size();
    size += 2*sizeof(size_t) + 2*sizeof(uint) + sizeof(_blocks) + _block_rows*_block_cols*sizeof(TMatrix< value_t >*);

    //
    // compute memory consumption of sub blocks
    //

    for ( uint i = 0; i < block_rows(); i++ )
    {
        for ( uint j = 0; j < block_cols(); j++ )
        {
            auto  A_ij = block(i,j);
            
            if ( A_ij != nullptr )
            {
                size += ulong(A_ij->byte_size());
            }// if
        }// for
    }// for

    return size;
}

//
// return size of (floating point) data in bytes handled by this object
//
template < typename value_t >
size_t
TBlockMatrix< value_t >::data_byte_size () const
{
    size_t  size = 0;
    
    //
    // compute memory consumption of sub blocks
    //

    for ( uint i = 0; i < block_rows(); i++ )
    {
        for ( uint j = 0; j < block_cols(); j++ )
        {
            auto  A_ij = block(i,j);
            
            if ( A_ij != nullptr )
                size += A_ij->data_byte_size();
        }// for
    }// for

    return size;
}

//
// serialisation
//

template < typename value_t >
void
TBlockMatrix< value_t >::read ( TByteStream & s )
{
    TMatrix< value_t >::read( s );

    size_t  n, m;

    s.get( n );
    s.get( m );

    set_size( n, m );
    
    uint  bn, bm;

    s.get( bn );
    s.get( bm );

    // compare with current setting
    if (( _block_rows != bn ) || ( _block_cols != bm ))
        HERROR( ERR_MAT_STRUCT, "(TBlockMatrix) read", "wrong block-structure in stream" );

    std::vector< char >  present( _block_rows * _block_cols );
    
    s.get( present.data(), _block_rows * _block_cols * sizeof(char) );
    
    // iterate over submatrices and read stream
    for ( uint i = 0; i < _block_rows*_block_cols; i++ )
    {
        if ( present[i] != 0 )
        {
            if ( _blocks[i] != nullptr )
                _blocks[i]->read( s );
            else
                HERROR( ERR_CONSISTENCY, "(TBlockMatrix) read", "nullptr submatrix but matrix in stream" );
        }// if
        else
        {
            if ( _blocks[i] != nullptr )
                HERROR( ERR_CONSISTENCY, "(TBlockMatrix) read", "non-nullptr submatrix but no matrix in stream" );
        }// else
    }// for
}

template < typename value_t >
void
TBlockMatrix< value_t >::build ( TByteStream & s )
{
    TMatrix< value_t >::build( s );

    size_t  n, m;

    s.get( n );
    s.get( m );

    set_size( n, m );
    
    uint  bn, bm;
    
    s.get( bn );
    s.get( bm );

    set_block_struct( bn, bm );

    //
    // build submatrices
    //

    std::vector< char >  present( _block_rows * _block_cols );
    TBSHMBuilder         builder;
    
    s.get( present.data(), _block_rows * _block_cols * sizeof(char) );
    
    for ( uint i = 0; i < _block_rows*_block_cols; i++ )
    {
        if ( present[i] != 0 )
        {
            auto  B_ij = builder.build< value_t >( s );
            
            if ( B_ij.get() == nullptr )
                HWARNING( "in (TBlockMatrix) build : submatrix is nullptr" );

            _blocks[i] = B_ij.release();
        }// if
        else
            _blocks[i] = nullptr;
    }// for
}

template < typename value_t >
void
TBlockMatrix< value_t >::write ( TByteStream & s ) const
{
    TMatrix< value_t >::write( s );

    s.put( _rows );
    s.put( _cols );

    s.put( _block_rows );
    s.put( _block_cols );

    // write boolean array indicating non-nullptr submatrices
    std::vector< char >  present( _block_rows * _block_cols );
    
    for ( uint i = 0; i < _block_rows*_block_cols; i++ )
    {
        if ( _blocks[i] != nullptr )
            present[i] = 1;
        else
            present[i] = 0;
    }// for

    s.put( present.data(), _block_rows*_block_cols*sizeof(char) );
    
    // iterate over submatrices and read stream
    for ( uint i = 0; i < _block_rows*_block_cols; i++ )
    {
        if ( _blocks[i] != nullptr )
            _blocks[i]->write( s );
    }// for
}

//
// returns size of object in bytestream
//
template < typename value_t >
size_t
TBlockMatrix< value_t >::bs_size () const
{
    size_t  size = ( TMatrix< value_t >::bs_size() +
                     sizeof(_rows) + sizeof(_cols) +
                     sizeof(_block_rows) + sizeof(_block_cols) +
                     _block_rows * _block_cols * sizeof(char) );
        
    for ( uint i = 0; i < _block_rows*_block_cols; i++ )
    {
        if ( _blocks[i] != nullptr )
            size += _blocks[i]->bs_size();
    }// for

    return size;
}

//
// test data for invalid values, e.g. INF and NAN
//
template < typename value_t >
void
TBlockMatrix< value_t >::check_data () const
{
    for ( uint  i = 0; i < block_rows(); ++i )
    {
        for ( uint  j = 0; j < block_cols(); ++j )
        {
            if ( block(i,j) != nullptr )
                block(i,j)->check_data();
        }// for
    }// for
}

//
// explicit instantiation
//

template class TBlockMatrix< float >;
template class TBlockMatrix< double >;
template class TBlockMatrix< std::complex< float > >;
template class TBlockMatrix< std::complex< double > >;

}// namespace Hpro
