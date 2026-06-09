//
// Project     : HLIBpro
// File        : TMatrix.cc
// Description : baseclass for all matrix-classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <iostream>

#include "hpro/base/error.hh"

#include "hpro/parallel/NET.hh"

#include "hpro/io/TMatrixIO.hh"

#include "hpro/matrix/structure.hh"
#include "hpro/matrix/TMatrixHierarchy.hh"
#include "hpro/matrix/TBlockMatrix.hh"

#include "hpro/matrix/TMatrix.hh"

namespace Hpro
{

//
// construct zero sized matrix
//
template < typename value_t >
TMatrix< value_t >::TMatrix ()
        : _id(-1)
        , _cluster(nullptr)
        , _row_ofs(0)
        , _col_ofs(0)
        , _procs(PROCSET_INVALID)
        , _matform(MATFORM_NONSYM)
        , _parent( nullptr )
        , _prev_in_block_row( nullptr )
        , _next_in_block_row( nullptr )
        , _prev_in_block_col( nullptr )
        , _next_in_block_col( nullptr )
        , _row_diag( nullptr )
        , _col_diag( nullptr )
{}

//
// construct matrix of size defined by block cluster \a c
//
template < typename value_t >
TMatrix< value_t >::TMatrix ( const TBlockCluster *  bc )
        : _id(-1)
        , _cluster(nullptr)
        , _row_ofs(0)
        , _col_ofs(0)
        , _procs(PROCSET_INVALID)
        , _matform(MATFORM_NONSYM)
        , _parent( nullptr )
        , _prev_in_block_row( nullptr )
        , _next_in_block_row( nullptr )
        , _prev_in_block_col( nullptr )
        , _next_in_block_col( nullptr )
        , _row_diag( nullptr )
        , _col_diag( nullptr )
{
    set_cluster( bc );
}

//
// construct matrix of size defined by block cluster \a c
//
template < typename value_t >
TMatrix< value_t >::TMatrix ( const TIndexSet &  arowis,
                              const TIndexSet &  acolis )
        : _id(-1)
        , _cluster(nullptr)
        , _row_ofs(0)
        , _col_ofs(0)
        , _procs(PROCSET_INVALID)
        , _matform(MATFORM_NONSYM)
        , _parent( nullptr )
        , _prev_in_block_row( nullptr )
        , _next_in_block_row( nullptr )
        , _prev_in_block_col( nullptr )
        , _next_in_block_col( nullptr )
        , _row_diag( nullptr )
        , _col_diag( nullptr )
{
    set_block_is( bis( arowis, acolis ) );
}

//
// construct matrix of size defined by block cluster \a c
//
template < typename value_t >
TMatrix< value_t >::TMatrix ( const TBlockIndexSet &  bis )
        : _id(-1)
        , _cluster(nullptr)
        , _row_ofs(0)
        , _col_ofs(0)
        , _procs(PROCSET_INVALID)
        , _matform(MATFORM_NONSYM)
        , _parent( nullptr )
        , _prev_in_block_row( nullptr )
        , _next_in_block_row( nullptr )
        , _prev_in_block_col( nullptr )
        , _next_in_block_col( nullptr )
        , _row_diag( nullptr )
        , _col_diag( nullptr )
{
    set_block_is( bis );
}

//
// copy constructor
//
template < typename value_t >
TMatrix< value_t >::TMatrix ( const TMatrix< value_t > &  A )
        : _id(-1)
        , _cluster(nullptr)
        , _row_ofs(A.row_ofs())
        , _col_ofs(A.col_ofs())
        , _procs(PROCSET_INVALID)
        , _matform(MATFORM_NONSYM)
        , _parent( nullptr )
        , _prev_in_block_row( nullptr )
        , _next_in_block_row( nullptr )
        , _prev_in_block_col( nullptr )
        , _next_in_block_col( nullptr )
        , _row_diag( nullptr )
        , _col_diag( nullptr )
{
    _id      = A._id;
    // _cluster = A._cluster;
    _procs   = A._procs;
    _matform = A._matform;
}

//
// access cluster
//
template < typename value_t >
void
TMatrix< value_t >::set_cluster ( const TBlockCluster * c )
{
    if ( c != nullptr )
    {
        _id      = c->id();
        _row_ofs = c->rowcl()->first();
        _col_ofs = c->colcl()->first();
        set_procs( c->procs() );
    }// if
}

namespace
{

//
// set hierarchy information in all matrix blocks
//
template < typename value_t >
void
rec_set_hierarchy ( TMatrix< value_t > *           M,
                    TMatrixHierarchy< value_t > *  H,
                    uint                           level )
{
    {
        auto                  block_row = H->matrix( level )->block_row( M->row_is() );
        TMatrix< value_t > *  prev      = nullptr;
        TMatrix< value_t > *  curr      = nullptr;
        TMatrix< value_t > *  next      = nullptr;
        TMatrix< value_t > *  diag      = nullptr;

        for ( auto  mat : * block_row )
        {
            if ( is_on_diag( mat ) )
            {
                diag = mat;
                break;
            }// if
        }// for

        for ( auto  mat : * block_row )
        {
            prev = curr;
            curr = next;
            next = mat;

            if ( curr != nullptr )
            {
                curr->set_block_row_neighbours( prev, next );
                curr->set_row_diag( diag );
            }// if

            // sets data on last element in list
            if ( next != nullptr )
            {
                next->set_block_row_neighbours( curr, nullptr );
                next->set_row_diag( diag );
            }// if
        }// for
    }

    {
        auto                  block_col = H->matrix( level )->block_col( M->col_is() );
        TMatrix< value_t > *  prev      = nullptr;
        TMatrix< value_t > *  curr      = nullptr;
        TMatrix< value_t > *  next      = nullptr;
        TMatrix< value_t > *  diag      = nullptr;

        for ( auto  mat : * block_col )
        {
            if ( is_on_diag( mat ) )
            {
                diag = mat;
                break;
            }// if
        }// for
        
        for ( auto  mat : * block_col )
        {
            prev = curr;
            curr = next;
            next = mat;

            if ( curr != nullptr )
            {
                curr->set_block_col_neighbours( prev, next );
                curr->set_col_diag( diag );
            }// if

            // sets data on last element in list
            if ( next != nullptr )
            {
                next->set_block_col_neighbours( curr, nullptr );
                next->set_col_diag( diag );
            }// if
        }// for
    }

    if ( is_blocked( M ) )
    {
        auto  B = ptrcast( M, TBlockMatrix< value_t > );

        for ( uint i = 0; i < B->nblock_rows(); ++i )
        {
            for ( uint j = 0; j < B->nblock_cols(); ++j )
            {
                if ( B->block( i, j ) != nullptr )
                    rec_set_hierarchy( B->block( i, j ), H, level+1 );
            }// for
        }// for
    }// if
}

//
// set hierarchy information in all matrix blocks
//
template < typename value_t >
void
flat_set_hierarchy ( TBlockMatrix< value_t > *      M,
                     TMatrixHierarchy< value_t > *  H,
                     uint                           level )
{
    for ( uint i = 0; i < M->nblock_rows(); ++i )
    {
        TMatrix< value_t > *  M_i = nullptr;

        // look for non-null block in current block row
        for ( uint j = 0; j < M->nblock_cols(); ++j )
        {
            if ( M->block( i, j ) != nullptr )
            {
                M_i = M->block( i, j );
                break;
            }// if
        }// for

        if ( M_i == nullptr )
            continue;

        //
        // set pointers in current block row
        //
        
        auto                  block_row = H->matrix( level+1 )->block_row( M_i->row_is() );
        TMatrix< value_t > *  prev      = nullptr;
        TMatrix< value_t > *  curr      = nullptr;
        TMatrix< value_t > *  next      = nullptr;
        TMatrix< value_t > *  diag      = nullptr;

        for ( auto  mat : * block_row )
        {
            if ( is_on_diag( mat ) )
            {
                diag = mat;
                break;
            }// if
        }// for

        for ( auto  mat : * block_row )
        {
            prev = curr;
            curr = next;
            next = mat;

            if ( curr != nullptr )
            {
                curr->set_block_row_neighbours( prev, next );
                curr->set_row_diag( diag );
            }// if

            // sets data on last element in list
            if ( next != nullptr )
            {
                next->set_block_row_neighbours( curr, nullptr );
                next->set_row_diag( diag );
            }// if
        }// for
    }// for

    for ( uint j = 0; j < M->nblock_cols(); ++j )
    {
        TMatrix< value_t > *  M_i = nullptr;

        // look for non-null block in current block column
        for ( uint i = 0; i < M->nblock_rows(); ++i )
        {
            if ( M->block( i, j ) != nullptr )
            {
                M_i = M->block( i, j );
                break;
            }// if
        }// for

        if ( M_i == nullptr )
            continue;

        //
        // set pointers in current block column
        //
        
        auto                  block_col = H->matrix( level+1 )->block_col( M_i->col_is() );
        TMatrix< value_t > *  prev      = nullptr;
        TMatrix< value_t > *  curr      = nullptr;
        TMatrix< value_t > *  next      = nullptr;
        TMatrix< value_t > *  diag      = nullptr;

        for ( auto  mat : * block_col )
        {
            if ( is_on_diag( mat ) )
            {
                diag = mat;
                break;
            }// if
        }// for
        
        for ( auto  mat : * block_col )
        {
            prev = curr;
            curr = next;
            next = mat;

            if ( curr != nullptr )
            {
                curr->set_block_col_neighbours( prev, next );
                curr->set_col_diag( diag );
            }// if

            // sets data on last element in list
            if ( next != nullptr )
            {
                next->set_block_col_neighbours( curr, nullptr );
                next->set_col_diag( diag );
            }// if
        }// for
    }

    //
    // proceed to subblocks if they are not leaves
    //
    
    for ( uint i = 0; i < M->nblock_rows(); ++i )
    {
        for ( uint j = 0; j < M->nblock_cols(); ++j )
        {
            if (( M->block( i, j ) != nullptr ) && is_blocked( M->block( i, j ) ) )
                rec_set_hierarchy( M->block( i, j ), H, level+1 );
        }// for
    }// for
}

}// namespace anonymous

//
// set hierarchy data automatically
//
template < typename value_t >
void
TMatrix< value_t >::set_hierarchy_data ()
{
    //
    // go up hierarchy as much as possible
    //

    auto  M        = this;
    auto  M_parent = M->parent();

    while ( M_parent != nullptr )
    {
        M        = M_parent;
        M_parent = M->parent();
    }// while

    //
    // adjust hierarchy pointers
    //

    auto  H = std::make_unique< TMatrixHierarchy< value_t > >( M, true );

    if ( is_flat( M ) )
        flat_set_hierarchy( ptrcast( M, TBlockMatrix< value_t > ), H.get(), 0 );
    else
        rec_set_hierarchy( M, H.get(), 0 );
}

namespace
{

//
// recursively apply M to all dense subblocks of M
// - return true if all subblocks of M are dense
//
template < typename value_t >
bool
apply_only_dense ( TMatrix< value_t > *        M,
                   const TMatrix< value_t > *  Upd )
{
    bool  all_dense = true;
    
    if (( M == nullptr ) || ( Upd == nullptr ))
        return all_dense;
    
    if ( M->is_blocked() )
    {
        auto    B = ptrcast( M, TBlockMatrix< value_t > );

        for ( uint  i = 0; i < B->nblock_rows(); ++i )
            for ( uint  j = 0; j < B->nblock_cols(); ++j )
            {
                auto  res = apply_only_dense( B->block( i, j ), Upd );

                all_dense = all_dense && res;
            }// for
    }// if
    else if ( M->is_dense() )
    {
        auto  D = ptrcast( M, TDenseMatrix< value_t > );
        
        if ( Upd->is_dense() )
        {
            D->add_block( value_t(1), value_t(1), cptrcast( Upd, TDenseMatrix< value_t > ), apply_normal );
        }// if
        else if ( is_lowrank( Upd ) )
        {
            auto         RUpd = cptrcast( Upd, TRkMatrix< value_t > );
            TScopedLock  lock( * D );
                
            // copy A(A)·B(A)^H to dense block
            BLAS::Matrix< value_t >  matA( blas_mat_A< value_t >( RUpd ),
                                           intersect( D->row_is(), Upd->row_is() ) - Upd->row_ofs(),
                                           BLAS::Range::all );
            BLAS::Matrix< value_t >  matB( blas_mat_B< value_t >( RUpd ),
                                           intersect( D->col_is(), Upd->col_is() ) - Upd->col_ofs(),
                                           BLAS::Range::all );
                
            BLAS::prod( value_t(1), matA, adjoint(matB), value_t(1), blas_mat< value_t >( D ) );
        }// if
        else
            HERROR( ERR_CONSISTENCY, "apply_only_dense", Upd->typestr() );

        all_dense = true;
    }// if
    else
    {
        //
        // non-dense block detected
        //
        
        all_dense = false;
    }// else

    return all_dense;
}

}// namespace anonymous

//
// add update matrix
//
template < typename value_t >
void
TMatrix< value_t >::add_update ( const TMatrix< value_t > *  M,
                                 const TTruncAcc &           acc )
{
    if ( M == nullptr )
        return;

    const auto  bis  = block_is();
    const auto  bisM = M->block_is();

    if ( bis == bisM )
    {
        if ( ! ( CFG::Arith::dense_accu || CFG::Arith::lazy_eval ) )
        {
            //
            // apply the update immediately to all dense blocks and
            // don't store update if only dense subblocks exist
            //
            
            const auto  all_dense = apply_only_dense( this, M );
            
            if ( all_dense )
                return;
        }
    
        accumulator().add_update( M, acc, this );
    }// if
    else if ( ( bis != bisM ) && bis.is_subset_of( bisM ) )
    {
        //
        // as this is a parent update, we don't need to apply this to
        // dense subblocks as that was done when the update was initially applied
        //
        
        if ( ! ( CFG::Arith::dense_accu || CFG::Arith::lazy_eval ) && is_dense() )
            return;
    
        accumulator().init( this );
        accumulator().add_parent_update( M, acc );
    }// if
    else
        HERROR( ERR_CONSISTENCY, "", "update neither local nor from parent" );
}
    
//
// add update U to set of recursive pending updates
//
template < typename value_t >
void
TMatrix< value_t >::add_pending_direct ( TDirectMatrixUpdate< value_t > *  U )
{
    accumulator().add_pending_direct( U );
}

//
// add update U to set of recursive pending updates
//
template < typename value_t >
void
TMatrix< value_t >::add_pending_recursive ( TRecursiveMatrixUpdate< value_t > *  U )
{
    accumulator().add_pending_recursive( U );
}
        
//
// apply locally stored updates with given accuracy
//
template < typename value_t >
void
TMatrix< value_t >::apply_updates ( const TTruncAcc &       /* acc */,
                                    const recursion_type_t )
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// return true, if matrix has updates not yet applied
//
template < typename value_t >
bool
TMatrix< value_t >::has_updates ( const recursion_type_t ) const
{
    return false;
}

//
// return true, if matrix has updates not yet applied
//
template < typename value_t >
bool
TMatrix< value_t >::has_parent_updates ( const recursion_type_t  t ) const
{
    bool  updates = has_updates( nonrecursive );
    
    if (( _parent != nullptr ) && ( t == recursive ))
        updates = updates | _parent->has_parent_updates( t );

    return updates;
}

//
// set processor set of matrix
//
template < typename value_t >
void
TMatrix< value_t >::set_procs ( const TProcSet &  ps,
                                const recursion_type_t ) // recursive is ignored
{
    _procs = ps;
}

//
// output of matrix
//
template < typename value_t >
void
TMatrix< value_t >::print ( const uint  ofs ) const
{
    for ( uint i = 0; i < ofs; i++ )
        std::cout << ' ';

    std::cout << this->typestr() << block_is() << std::endl;
}

///////////////////////////////////////////
//
// single coefficient access
//

//
// return index (i,j) of the matrix
//
template < typename value_t >
value_t
TMatrix< value_t >::entry  ( const idx_t, const idx_t ) const
{
    HERROR( ERR_NOT_IMPL, "(TMatrix) entry", "" );
    return value_t(0);
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
TMatrix< value_t >::scale ( const value_t )
{
    HERROR( ERR_NOT_IMPL, "(TMatrix) scale", "" );
}

//
// this = this + a * matrix
//
template < typename value_t >
void
TMatrix< value_t >::add ( const value_t, const TMatrix< value_t > * )
{
    HERROR( ERR_NOT_IMPL, "(TMatrix) add", "" );
}

//
// matrix-vector-multiplication
//
template < typename value_t >
void
TMatrix< value_t >::mul_vec ( const value_t,
                              const TVector< value_t > *,
                              const value_t,
                              TVector< value_t > *,
                              const matop_t ) const
{
    HERROR( ERR_NOT_IMPL, "(TMatrix) mul_vec", "" );
}

//
// return alpha * this * B
//
template < typename value_t >
TMatrix< value_t > *
TMatrix< value_t >::mul_right ( const value_t, const TMatrix< value_t > *,
                                const matop_t, const matop_t ) const
{
    HERROR( ERR_NOT_IMPL, "(TMatrix) mul_right", "" );
    return nullptr;
}

//
// return alpha * A * this
//
template < typename value_t >
TMatrix< value_t > *
TMatrix< value_t >::mul_left ( const value_t, const TMatrix< value_t > *,
                               const matop_t, const matop_t ) const
{
    HERROR( ERR_NOT_IMPL, "(TMatrix) mul_left", "" );
    return nullptr;
}

///////////////////////////////////////////////////////////
//
// linear operator mapping
//

template < typename value_t >
void
TMatrix< value_t >::apply_add   ( const value_t               /* alpha */,
                                  const TMatrix< value_t > *  /* X */,
                                  TMatrix< value_t > *        /* Y */,
                                  const matop_t               /* op */ ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
    // multiply( alpha, op, this, apply_normal, X, value_t(1), Y, acc_exact );
}

//
// same as above but only the dimension of the vector spaces is tested,
// not the corresponding index sets
//
template < typename value_t >
void
TMatrix< value_t >::apply_add   ( const value_t                    alpha,
                                  const BLAS::Vector< value_t > &  x,
                                  BLAS::Vector< value_t > &        y,
                                  const matop_t                    op ) const
{
    TScalarVector< value_t >  sx( col_is( op ), x );
    TScalarVector< value_t >  sy( row_is( op ), y );

    apply_add( alpha, & sx, & sy, op );
}

template < typename value_t >
void
TMatrix< value_t >::apply_add   ( const value_t                    alpha,
                                  const BLAS::Matrix< value_t > &  X,
                                  BLAS::Matrix< value_t > &        Y,
                                  const matop_t                    op ) const
{
    if ( std::imag( alpha ) != real_type_t< value_t >(0) )
        HERROR( ERR_NOT_IMPL, "(TMatrix) apply_add", "only real valued alpha implemented" );
    
    TDenseMatrix< value_t >  DX( col_is( op ), is( 0, X.ncols()-1 ), X );
    TDenseMatrix< value_t >  DY( row_is( op ), is( 0, Y.ncols()-1 ), Y );

    apply_add( alpha, & DX, & DY, op );
}

//
// virtual constructor
//

//
// return copy of matrix
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TMatrix< value_t >::copy () const
{
    auto  M = create();

    M->set_id( id() );
    M->set_form( form() );
    M->set_ofs( row_ofs(), col_ofs() );
    M->set_size( rows(), cols() );
    
    M->_procs = _procs;

    // if ( _cluster != nullptr )
    //     M->set_cluster( _cluster );
    
    return M;
}

//
// return structural copy of matrix
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TMatrix< value_t >::copy_struct () const
{
    // no data here anyway
    return TMatrix< value_t >::copy();
}

//
// return copy of matrix wrt. given accuracy
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TMatrix< value_t >::copy ( const TTruncAcc &,
                           const bool ) const
{
    // by default return exact copy
    return copy();
}

//
// copy data from matrix \a A
//
template < typename value_t >
void
TMatrix< value_t >::copy_from ( const TMatrix< value_t > * A )
{
    if ( A == nullptr )
        HERROR( ERR_ARG, "(TMatrix) copy_from", "A is NULL" );

    A->copy_to( this );
}

//
// copy matrix into given matrix
//
template < typename value_t >
void
TMatrix< value_t >::copy_to ( TMatrix< value_t > * A ) const
{
    if ( A == nullptr )
        HERROR( ERR_ARG, "(TMatrix) copy_to", "argument is nullptr" );

    // if ( _cluster != nullptr )
    //     A->set_cluster( _cluster );

    A->_id      = _id;
    A->_row_ofs = _row_ofs;
    A->_col_ofs = _col_ofs;
    A->_matform = _matform;
    A->_procs   = _procs;
}

//
// copy matrix into given matrix wrt. given accuracy
//
template < typename value_t >
void
TMatrix< value_t >::copy_to ( TMatrix< value_t > *  A,
                              const TTruncAcc & ,
                              const bool ) const
{
    // by default copy exactly
    copy_to( A );
}

//
// transpose matrix
//
template < typename value_t >
void
TMatrix< value_t >::transpose ()
{
    // only change index set offsets and leave size change to
    // derived classes for efficiency reasons
    set_ofs( col_ofs(), row_ofs() );
}
    
//
// conjugate matrix coefficients
//
template < typename value_t >
void
TMatrix< value_t >::conjugate ()
{
}

//
// return size in bytes used by this object
//
template < typename value_t >
size_t
TMatrix< value_t >::byte_size () const
{
    return ( sizeof(_id) + sizeof(_cluster) + sizeof(_row_ofs) + sizeof(_col_ofs) +
             _procs.byte_size() + sizeof(_matform) + TLockable::byte_size() );
    // 7*sizeof(TMatrix<value_t>*) + sizeof(_accumulator) + TLockable::byte_size() );
}

//
// return size in bytes used by this distributed object,
// i.e. of all distributed sub matrices
//
template < typename value_t >
size_t
TMatrix< value_t >::global_byte_size () const
{
    // sum up over all nodes in local processor set
    size_t  local_size = byte_size();
    size_t  global_size;

    if ( procs().size() > 1 )
        NET::reduce_all( procs(), & local_size, & global_size, 1, NET::OP_SUM );
    else
        global_size = local_size;

    return global_size;
}

//
// serialisation
//

template < typename value_t >
void
TMatrix< value_t >::read ( TByteStream & s )
{
    typeid_t  t;

    TStreamable::read( s );
    
    s.get( t );

    if ( t != type() )
        HERROR( ERR_BS_TYPE, "(TMatrix) read", "at " + block_is().to_string() + " found \"" + RTTI::id_to_type( t ) + "\","
               + " expected \"" + this->typestr() + "\"" );

    s.get( _id );
    s.get( _row_ofs );
    s.get( _col_ofs );

    _procs.read( s );
    
    s.get( _matform );
}

template < typename value_t >
void
TMatrix< value_t >::build ( TByteStream & s )
{
    s.get( _id );
    s.get( _row_ofs );
    s.get( _col_ofs );

    _procs.read( s );
    
    s.get( _matform );
}

template < typename value_t >
void
TMatrix< value_t >::write ( TByteStream & s ) const
{
    const typeid_t  t = type();

    TStreamable::write( s );
    
    s.put( t );
    s.put( _id );
    s.put( _row_ofs );
    s.put( _col_ofs );

    _procs.write( s );
    
    s.put( _matform );
}

//
// returns size of object in bytestream
//
template < typename value_t >
size_t
TMatrix< value_t >::bs_size () const
{
    return ( TStreamable::bs_size() +
             sizeof(typeid_t) + sizeof(_id) +
             sizeof(_row_ofs) + sizeof(_col_ofs) +
             _procs.bs_size() + sizeof(_matform) );
}

//
// parallel methods
//

//
// sum up nparts parallel copies
// (if bs != nullptr it will be used)
//
template < typename value_t >
void
TMatrix< value_t >::sum ( const TProcSet &   /* ps */,
                          const uint         /* pid */,
                          const uint         /* parts */,
                          TByteStream *      /* bs */,
                          const TTruncAcc &  /* acc */ )
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// tests
//

//
// test data for invalid values, e.g. INF and NAN
//
template < typename value_t >
void
TMatrix< value_t >::check_data () const
{
}
    
//////////////////////////////////////////////////////////
//
// debug helpers
//

namespace DBG
{

//
// write vector to file
//
template < typename value_t >
void
write ( const TMatrix< value_t > *  M,
        const std::string &         filename,
        const std::string &         matname )
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "write", "matrix is null" );
    
    write( *M, filename, matname );
}

template < typename value_t >
void
write ( const TMatrix< value_t > &  M,
        const std::string &         filename,
        const std::string &         matname )
{
    TMatlabMatrixIO  mio( false );

    mio.write( & M, filename, matname );
}

//
// write vector to file
//
template < typename value_t >
void
write ( const TLinearOperator< value_t > *  M,
        const std::string &                 filename,
        const std::string &                 matname )
{
    TMatlabMatrixIO  mio( false );

    mio.write( M, filename, matname );
}

#define INST_WRITE( type ) \
    template void write< type > ( const TMatrix< type > *, const std::string &, const std::string & ); \
    template void write< type > ( const TMatrix< type > &, const std::string &, const std::string & ); \
    template void write< type > ( const TLinearOperator< type > &, const std::string &, const std::string & );

INST_WRITE( float )
INST_WRITE( double )
INST_WRITE( std::complex< float > )
INST_WRITE( std::complex< double > )

}// namespace DBG

template class TMatrix< float >;
template class TMatrix< double >;
template class TMatrix< std::complex< float > >;
template class TMatrix< std::complex< double > >;

}// namespace Hpro
