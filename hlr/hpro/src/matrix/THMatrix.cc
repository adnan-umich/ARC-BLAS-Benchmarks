//
// Project     : HLIBpro
// File        : THMatrix.cc
// Description : class for H-matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/parallel/NET.hh"

#include "hpro/matrix/structure.hh"
#include "hpro/matrix/TMatrixHierarchy.hh"

#include "hpro/matrix/THMatrix.hh"

namespace Hpro
{

using std::unique_ptr;

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//
// model independent H-Matrix implementation
//
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////
//
// constructor
//
template < typename value_t >
THMatrix< value_t >::THMatrix ( const TBlockCluster * bct )
        : TBlockMatrix< value_t >( bct )
        , _max_id( 0 )
{}

template < typename value_t >
THMatrix< value_t >::THMatrix ( const TBlockIndexSet &  bis )
        : TBlockMatrix< value_t >( bis )
        , _max_id( 0 )
{
}

template < typename value_t >
THMatrix< value_t >::THMatrix ( const TIndexSet &  arowis,
                                const TIndexSet &  acolis )
        : TBlockMatrix< value_t >( arowis, acolis )
        , _max_id( 0 )
{
}

/////////////////////////////////////////////////
//
// access internal data
//

//
// access single matrix coefficient
//
template < typename value_t >
value_t
THMatrix< value_t >::entry ( const idx_t row, const idx_t col ) const
{
    idx_t  srow = row;
    idx_t  scol = col;
    bool   swap_idx = false;
    
    // translate index pair to lower half in case of symmetry
    if ( ! this->is_nonsym() && ( srow < scol ))
    {
        std::swap( srow, scol );
        swap_idx = true;
    }// if
    
    const idx_t rofs = this->row_ofs();
    const idx_t cofs = this->col_ofs();
    
    for ( uint i = 0; i < this->block_rows(); i++ )
        for ( uint j = 0; j < this->block_cols(); j++ )
        {
            const TMatrix< value_t > * A_ij = this->block(i,j);

            if ( A_ij == nullptr )
                continue;
            
            if ( A_ij->block_is().is_in( rofs + srow, cofs + scol ) )
            {
                const auto  val = A_ij->entry( srow - A_ij->row_ofs() + rofs,
                                               scol - A_ij->col_ofs() + cofs );
                
                if ( swap_idx && this->is_hermitian() )
                    return Math::conj(val);
                else
                    return val;
            }// if
        }// for


    return value_t(0);
}

/////////////////////////////////////////////////
//
// BLAS-routines
//

//
// matrix-vector-mult.
//
template < typename value_t >
void
THMatrix< value_t >::mul_vec ( const value_t               alpha,
                               const TVector< value_t > *  x,
                               const value_t               beta,
                               TVector< value_t > *        y,
                               const matop_t               op ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(THMatrix) mul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(THMatrix) mul_vec", "y = nullptr" );

    //
    // multiply
    //

    TBlockMatrix< value_t >::mul_vec( alpha, x, beta, y, op );
}

/////////////////////////////////////////////////
//
// misc.
//

//
// return appropriate vector-types for matrix
//
template < typename value_t >
auto
THMatrix< value_t >::row_vector () const -> unique_ptr< TVector< value_t > >
{
    return std::make_unique< TScalarVector< value_t > >( this->rows(), this->row_ofs() );
}

template < typename value_t >
auto
THMatrix< value_t >::col_vector () const -> unique_ptr< TVector< value_t > >
{
    return std::make_unique< TScalarVector< value_t > >( this->cols(), this->col_ofs() );
}
    
//
// transpose matrix
//
template < typename value_t >
void
THMatrix< value_t >::transpose ()
{
    TBlockMatrix< value_t >::transpose();

    swap( _row_perm_e2i, _col_perm_e2i );
    swap( _row_perm_i2e, _col_perm_i2e );

    comp_min_max_idx();
}
    
//
// return size in bytes used by this object
//
template < typename value_t >
size_t
THMatrix< value_t >::byte_size () const
{
    size_t  size = TBlockMatrix< value_t >::byte_size();

    size += _row_perm_e2i.byte_size();
    size += _col_perm_e2i.byte_size();
    size += _row_perm_i2e.byte_size();
    size += _col_perm_i2e.byte_size();
    
#if 0
    size += _leaves.byte_size();
#endif

    return size;
}

//
// serialisation
//

template < typename value_t >
void
THMatrix< value_t >::read  ( TByteStream & s )
{
    TBlockMatrix< value_t >::read( s );
}

template < typename value_t >
void
THMatrix< value_t >::build ( TByteStream & s )
{
    TBlockMatrix< value_t >::build( s );
}

template < typename value_t >
void
THMatrix< value_t >::write ( TByteStream & s ) const
{
    TBlockMatrix< value_t >::write( s );
}

//
// returns size of object in bytestream
//
template < typename value_t >
size_t
THMatrix< value_t >::bs_size () const
{
    return (TBlockMatrix< value_t >::bs_size() +
            _row_perm_e2i.byte_size() +
            _col_perm_e2i.byte_size() +
            _row_perm_i2e.byte_size() +
            _col_perm_i2e.byte_size() );
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//
// model dependent H-Matrix implementation (BSP)
//
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

//
// compute min/max tau/sigma indices
//
template < typename value_t >
void
THMatrix< value_t >::comp_min_max_idx ()
{
    //
    // set up matrix hierarchy
    //

    this->set_hierarchy_data();
    _max_id = Hpro::max_id( this );
}

/////////////////////////////////////////////////
//
// misc.
//

//
// virtual constructor
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
THMatrix< value_t >::copy () const
{
    auto        M = TBlockMatrix< value_t >::copy();
    THMatrix *  H = ptrcast( M.get(), THMatrix< value_t > );

    H->comp_min_max_idx();
    H->set_row_perm( _row_perm_e2i, _row_perm_i2e );
    H->set_col_perm( _col_perm_e2i, _col_perm_i2e );
    
    return M;
}

//
// copy matrix wrt. given accuracy and coarsening
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
THMatrix< value_t >::copy ( const TTruncAcc & acc, const bool coarsen ) const
{
    auto        M = TBlockMatrix< value_t >::copy( acc, coarsen );
    THMatrix *  H = ptrcast( M.get(), THMatrix< value_t > );

    H->comp_min_max_idx();
    H->set_row_perm( _row_perm_e2i, _row_perm_i2e );
    H->set_col_perm( _col_perm_e2i, _col_perm_i2e );
    
    return M;
}

//
// return structural copy
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
THMatrix< value_t >::copy_struct () const
{
    auto        M = TBlockMatrix< value_t >::copy_struct();
    THMatrix *  H = ptrcast( M.get(), THMatrix< value_t > );

    H->comp_min_max_idx();
    H->set_row_perm( _row_perm_e2i, _row_perm_i2e );
    H->set_col_perm( _col_perm_e2i, _col_perm_i2e );
    
    return M;
}

//
// copy matrix into A
//
template < typename value_t >
void
THMatrix< value_t >::copy_to ( TMatrix< value_t > * A ) const
{
    TBlockMatrix< value_t >::copy_to( A );

    if ( ! IS_TYPE( A, THMatrix ) )
        HERROR( ERR_MAT_TYPE, "(THMatrix) copy_to", A->typestr() );

    THMatrix * H = ptrcast( A, THMatrix );
    
    H->comp_min_max_idx();
    H->set_row_perm( _row_perm_e2i, _row_perm_i2e );
    H->set_col_perm( _col_perm_e2i, _col_perm_i2e );
}

template < typename value_t >
void
THMatrix< value_t >::copy_to ( TMatrix< value_t > * A, const TTruncAcc & acc, const bool coarsen ) const
{
    TBlockMatrix< value_t >::copy_to( A, acc, coarsen );

    if ( ! IS_TYPE( A, THMatrix ) )
        HERROR( ERR_MAT_TYPE, "(THMatrix) copy_to", A->typestr() );

    THMatrix * H = ptrcast( A, THMatrix );
    
    H->comp_min_max_idx();
    H->set_row_perm( _row_perm_e2i, _row_perm_i2e );
    H->set_col_perm( _col_perm_e2i, _col_perm_i2e );
}

//
// copy complete structural information from given matrix
//
template < typename value_t >
void
THMatrix< value_t >::copy_struct_from ( const TMatrix< value_t > * M )
{
    TBlockMatrix< value_t >::copy_struct_from( M );

    if ( IS_TYPE( M, THMatrix ) )
    {
        const THMatrix * H = cptrcast( M, THMatrix );

        comp_min_max_idx();
        set_row_perm( H->row_perm_e2i(), H->row_perm_i2e() );
        set_col_perm( H->col_perm_e2i(), H->col_perm_i2e() );
    }// if
}
    
//
// collect leaves
//
template < typename value_t >
bool
size_sort ( TMatrix< value_t > * M1, TMatrix< value_t > * M2 )
{
    return std::max( M1->rows(), M1->cols() ) > std::max( M2->rows(), M2->cols() );
}

// void
// THMatrix< value_t >::build_leaf_list ()
// {
//     std::list< TMatrix< value_t > * >  tmp;
    
//     TBlockMatrix< value_t >::collect_leaves( tmp );

//     _leaves.resize( tmp.size() );

//     size_t  i = 0;
//     for ( auto M : tmp ) _leaves[i++] = M;

//     sort( _leaves.begin(), _leaves.end(), size_sort );
// }

//
// explicit instantiation
//

template class THMatrix< float >;
template class THMatrix< double >;
template class THMatrix< std::complex< float > >;
template class THMatrix< std::complex< double > >;

}// namespace Hpro
