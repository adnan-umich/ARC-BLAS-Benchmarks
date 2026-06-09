//
// Project     : HLIBpro
// File        : TSparseMatrix.cc
// Description : class for a sparse matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <iostream>

#include "hpro/matrix/TSparseMatrix.hh"

namespace Hpro
{

using std::unique_ptr;
using std::make_unique;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// local functions
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

namespace
{

//
// return number of non-zeroes in lower left part plus diagonal
//
template < typename value_t >
size_t
lower_left_nonzeroes ( const TSparseMatrix< value_t > *  S )
{
    if ( S == nullptr )
        return 0;

    size_t  llnnz = 0;

    for ( idx_t  row = 0; row < idx_t(S->rows()); row++ )
    {
        const idx_t  lb = S->rowptr(row);
        const idx_t  ub = S->rowptr(row+1);

        for ( idx_t  j = lb; j < ub; ++j )
            if ( S->colind(j) <= row )
                llnnz++;
    }// for

    return llnnz;
}

}// namespace anonymous

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TSparseMatrix
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

///////////////////////////////////////////
//
// constructor and destructor
//

template < typename value_t >
TSparseMatrix< value_t >::~TSparseMatrix ()
{
    _rowptr.resize( 0 );
    _colind.resize( 0 );
    _coeff.resize( 0 );
}

/////////////////////////////////////////////////
//
// access data
//

//
// set cluster of matrix
//
template < typename value_t >
void
TSparseMatrix< value_t >::set_cluster ( const TBlockCluster * c )
{
    if ( _colind.size() > 0 )
        HERROR( ERR_CONSISTENCY, "(TSparseMatrix) set_cluster",
                "changing sparse matrix in CRS format" );
    
    if ( c != nullptr )
    {
        if (( c->rowcl() != nullptr ) && ( c->colcl() != nullptr ))
            set_size( c->rowcl()->size(), c->colcl()->size() );
    }// if
}

//
// directly set dimension of matrix
//
template < typename value_t >
void
TSparseMatrix< value_t >::set_size ( const size_t  n,
                                     const size_t  m )
{
    if (( _rows != n ) || ( _cols != m ))
    {
        if ( _colind.size() > 0 )
            HERROR( ERR_CONSISTENCY, "(TSparseMatrix) set_size",
                    "can not change non-zero sparse matrix" );
        
        _rows = n;
        _cols = m;

        // initialise with no non zeroes
        init( 0 );
    }// if
}

//
// return a_ij
//
template < typename value_t >
value_t
TSparseMatrix< value_t >::entry ( const idx_t  ai,
                                  const idx_t  aj ) const
{
    if ( ! ( this->row_is().is_in( ai ) && this->col_is().is_in( aj ) ) )
        HERROR( ERR_MAT_SIZE, "(TSparseMatrix) entry",
                to_string( "(%d,%d) not in indexset", ai, aj ) );
        
    if ( this->is_complex() )
        HERROR( ERR_COMPLEX, "(TSparseMatrix) entry", "matrix is complex valued" );

    const idx_t  lb = _rowptr[ai];
    const idx_t  ub = _rowptr[ai+1];
    
    for ( idx_t  pos = lb; pos < ub; pos++ )
    {
        if ( _colind[pos] == aj )
            return _coeff[pos];
    }// for

    return value_t(0);
}

//
// set entry a_ij (add if not present)
//
template < typename value_t >
void
TSparseMatrix< value_t >::set_entry ( const idx_t    ai,
                                      const idx_t    aj,
                                      const value_t  c )
{
    if ( ! ( this->row_is().is_in( ai ) && this->col_is().is_in( aj ) ) )
        HERROR( ERR_MAT_SIZE, "(TSparseMatrix) set_entry",
                to_string( "(%d,%d) not in indexset", ai, aj ) );
    
    const idx_t  lb = _rowptr[ai];
    const idx_t  ub = _rowptr[ai+1];

    for ( idx_t  pos = lb; pos < ub; pos++ )
    {
        if ( _colind[pos] == aj )
        {
            _coeff[pos] = c;
            
            return;
        }// if
    }// for

    HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) set_entry",
            to_string( "(%d,%d) not in sparsity pattern", ai, aj ) );
}

//
// add to entry a_ij if present
//
template < typename value_t >
void
TSparseMatrix< value_t >::add_entry ( const idx_t    ai,
                                      const idx_t    aj,
                                      const value_t  c )
{
    if ( ! ( this->row_is().is_in( ai ) && this->col_is().is_in( aj ) ) )
        HERROR( ERR_MAT_SIZE, "(TSparseMatrix) add_entry",
                to_string( "(%d,%d) not in indexset", ai, aj ) );
    
    const idx_t  lb = _rowptr[ai];
    const idx_t  ub = _rowptr[ai+1];

    for ( idx_t  pos = lb; pos < ub; pos++ )
    {
        if ( _colind[pos] == aj )
        {
            _coeff[pos] += c;
            
            return;
        }// if
    }// for

    HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) add_entry",
            to_string( "(%d,%d) not in sparsity pattern", ai, aj ) );
}

//
// add to entry a_ij (insert if not present)
//
template < typename value_t >
bool
TSparseMatrix< value_t >::has_entry ( const idx_t  ai,
                                      const idx_t  aj ) const
{
    if ( ! ( this->row_is().is_in( ai ) && this->col_is().is_in( aj ) ) )
        HERROR( ERR_MAT_SIZE, "(TSparseMatrix) has_entry",
                to_string( "(%d,%d) not in indexset", ai, aj ) );
    
    const idx_t  lb = _rowptr[ai];
    const idx_t  ub = _rowptr[ai+1];

    for ( idx_t  pos = lb; pos < ub; pos++ )
    {
        if ( _colind[pos] == aj )
            return true;
    }// for
    
    return false;
}

//
// sort list per row wrt column
//
template < typename value_t >
void
TSparseMatrix< value_t >::sort_entries ()
{
    //
    // sort CRS data structure
    //
    
    for ( size_t i = 0; i < rows(); i++ )
    {
        const idx_t   lb = _rowptr[i];
        const idx_t   ub = _rowptr[i+1];
        
        for ( idx_t  j = lb; j < ub; j++ )
        {
            bool changed = false;
            
            for ( idx_t  l = lb; l < (ub-j+lb-1); l++ )
                if ( _colind[l] > _colind[l+1] )
                {
                    std::swap( _colind[l], _colind[l+1] );
                    std::swap( _coeff[l],  _coeff[l+1] );
                    
                    changed = true;
                }// if
            
            if ( ! changed )
                break;
        }// for
    }// for
}

//
// initialise CRS data for n entries
//
template < typename value_t >
void
TSparseMatrix< value_t >::init ( const size_t nnz )
{
    _coeff.resize( nnz );
    _colind.resize( nnz );
    _rowptr.resize( rows() + 1 );

    // nullify CRS data
    for ( size_t  i = 0; i <= rows(); i++ )
        _rowptr[i] = 0;
}

//
// import CRS data
//
template < typename value_t >
template < typename T_idx >
void
TSparseMatrix< value_t >::import_crs ( const size_t     nrows,
                                       const size_t     ncols,
                                       const size_t     nnz,
                                       const T_idx *    arowptr,
                                       const T_idx *    acolind,
                                       const value_t *  acoeffs )
{
    init( 0 );

    _rows = nrows;
    _cols = ncols;

    _rowptr.resize( nrows+1 );
    _colind.resize( nnz );

    for ( size_t i = 0; i <= nrows; ++i ) _rowptr[i] = arowptr[i];
    for ( size_t i = 0; i <  nnz;   ++i ) _colind[i] = acolind[i];

    _coeff.resize( nnz );

    for ( size_t i = 0; i < nnz; ++i )
        _coeff[i]  = acoeffs[i];
}

//
// import CCS data
//
template < typename value_t >
template < typename T_idx >
void
TSparseMatrix< value_t >::import_ccs ( const size_t     nrows,
                                       const size_t     ncols,
                                       const size_t     nnz,
                                       const T_idx *    acolptr,
                                       const T_idx *    arowind,
                                       const value_t *  acoeffs )
{
    init( 0 );

    _rows = nrows;
    _cols = ncols;

    _rowptr.resize( nrows+1 );
    _colind.resize( nnz );

    _coeff.resize( nnz );

    // count number of entries per column
    for ( idx_t  col = 0; col < idx_t(ncols); ++col )
    {
        const idx_t  lb = acolptr[ col ];
        const idx_t  ub = acolptr[ col+1 ];

        for ( idx_t  j = lb; j < ub; j++ )
            _rowptr[ arowind[j] ]++;
    }// for

    // compute real rowptr array
    idx_t  pos = 0;

    for ( size_t  row = 0; row < nrows; ++row )
    {
        const idx_t  n = _rowptr[row];

        _rowptr[row] = pos;
        pos         += n;
    }// for
    _rowptr[nrows] = pos;

    if ( pos != idx_t(nnz) )
        HERROR( ERR_CONSISTENCY, "(TSparseMatrix) import_ccs",
                "given number of non-zeroes differs from counted number" );

    // set up rowind and coeff arrays
    std::vector< idx_t >   rowpos( nrows );

    for ( idx_t  col = 0; col < idx_t(ncols); ++col )
    {
        const idx_t  lb = acolptr[ col ];
        const idx_t  ub = acolptr[ col+1 ];

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  row = arowind[j];
            const idx_t  idx = _rowptr[row] + rowpos[row];

            _colind[ idx ] = col;
            _coeff[  idx ] = acoeffs[ j ];

            rowpos[row]++;
        }// for
    }// for
}

//
// export internal data to CCS format (real valued)
//
template < typename value_t >
template < typename T_idx >
void
TSparseMatrix< value_t >::export_ccs  ( std::vector< T_idx > &    colptr,
                                        std::vector< T_idx > &    rowind,
                                        std::vector< value_t > &  acoeffs,
                                        const bool                use_sym ) const
{
    const size_t  nrows  = rows();
    const size_t  ncols  = cols();
    const bool    is_sym = use_sym && ( this->is_symmetric() || this->is_hermitian() );
    const size_t  nnz    = ( is_sym ? lower_left_nonzeroes( this ) : n_non_zero() );

    colptr.resize( ncols+1 );

    for ( idx_t  col = 0; col <= idx_t(ncols); ++col )
        colptr[ col ] = 0;

    // count number of entries per column
    for ( idx_t  row = 0; row < idx_t(nrows); ++row )
    {
        const idx_t  lb = rowptr(idx_t(row));
        const idx_t  ub = rowptr(idx_t(row)+1);

        if ( is_sym )
        {
            for ( idx_t  j = lb; j < ub; j++ )
                if ( colind(j) <= row ) // count only lower left part
                    colptr[ colind(j) ]++;
        }// if
        else
        {
            for ( idx_t  j = lb; j < ub; j++ )
                colptr[ colind(j) ]++;
        }// else
    }// for

    // compute real colptr array
    idx_t  pos = 0;

    for ( size_t  j = 0; j < ncols; j++ )
    {
        const idx_t  n = colptr[j];

        colptr[j] = T_idx(pos);
        pos      += n;
    }// for
    colptr[ncols] = T_idx(pos);

    // set up rowind and coeff arrays
    std::vector< idx_t >   colpos( ncols );

    rowind.resize( nnz );
    acoeffs.resize( nnz );

    for ( idx_t  row = 0; row < idx_t(nrows); row++ )
    {
        const idx_t  lb = rowptr(row);
        const idx_t  ub = rowptr(row+1);

        if ( is_sym )
        {
            for ( idx_t  j = lb; j < ub; j++ )
            {
                const idx_t  col = colind(j);
                const idx_t  idx = colptr[col] + colpos[col];

                if ( col > row ) // only lower left part
                    continue;

                rowind[ idx ]  = T_idx( row );
                acoeffs[ idx ] = coeff(j);

                colpos[col]++;
            }// for
        }// if
        else
        {
            for ( idx_t  j = lb; j < ub; j++ )
            {
                const idx_t  col = colind(j);
                const idx_t  idx = colptr[col] + colpos[col];

                rowind[ idx ]  = T_idx( row );
                acoeffs[ idx ] = coeff(j);

                colpos[col]++;
            }// for
        }// else
    }// for
}
        

//
// permute entries in sparse matrix
//
template < typename value_t >
void
TSparseMatrix< value_t >::permute ( const TPermutation &  rowperm,
                                    const TPermutation &  colperm )
{
    TPermutation  inv_rowperm( rowperm );

    inv_rowperm.invert();

    std::vector< idx_t >    new_rowptr;
    std::vector< idx_t >    new_colind;
    std::vector< value_t >  new_coeff;

    new_coeff.resize(  _coeff.size() );
    new_colind.resize( _colind.size() );
    new_rowptr.resize( _rowptr.size() );

    //
    // copy data from old to new position
    //
    
    idx_t  cpos = 0;
    idx_t  rpos = 0;
    
    for ( size_t i = 0; i < rows(); i++ )
    {
        const idx_t  old_row = inv_rowperm.permute( idx_t(i) );
        const idx_t  lb      = _rowptr[old_row];
        const idx_t  ub      = _rowptr[old_row+1];

        new_rowptr[ i ] = rpos;
        
        for ( idx_t  j = lb; j < ub; j++, rpos++, cpos++ )
        {
            const idx_t   new_col = colperm.permute( _colind[j] );

            new_coeff[cpos]  = _coeff[j];
            new_colind[cpos] = new_col;
        }// for
    }// for

    new_rowptr[ rows() ] = _rowptr[ rows() ];

    //
    // copy to local data
    //
    
    for ( size_t i = 0; i < _coeff.size();  i++ ) _coeff[i]  = new_coeff[i];
    for ( size_t i = 0; i < _colind.size(); i++ ) _colind[i] = new_colind[i];
    for ( size_t i = 0; i < _rowptr.size(); i++ ) _rowptr[i] = new_rowptr[i];
}

//
// check symmetry of matrix, return true if symmetric
//
template < typename value_t >
bool
TSparseMatrix< value_t >::test_symmetry () 
{
    //
    // go through sparsity pattern and compare entries
    // with transposed entries
    //

    const idx_t  n       = idx_t(rows());
    const idx_t  m       = idx_t(cols());
    bool         is_sym  = true;
    bool         is_herm = true;

    this->set_nonsym();
    
    for ( idx_t  i = 0; i < n; i++ )
    {
        // stop as soon as not symmetric and not hermitian
        if ( ! ( is_sym || is_herm ) )
            break;
        
        const idx_t  lb = _rowptr[i];
        const idx_t  ub = _rowptr[i+1];
        
        for ( idx_t j = lb; j < ub; j++ )
        {
            const idx_t  col = _colind[j];
            
            if (( col >= n ) || ( i >= m ))
                is_herm = is_sym = false;
            
            const value_t  f = entry( col, i );

            if ( _coeff[j] != f               ) is_sym  = false;
            if ( _coeff[j] != Math::conj( f ) ) is_herm = false;
        }// for
    }// for

    if      ( is_sym  ) this->set_symmetric();
    else if ( is_herm ) this->set_hermitian();
    
    return is_sym || is_herm;
}

//
// compute average/maximal number of entries per row
//
template < typename value_t >
size_t
TSparseMatrix< value_t >::avg_entries_per_row () const
{
    return _colind.size() / rows();
}
    
template < typename value_t >
size_t
TSparseMatrix< value_t >::max_entries_per_row () const
{
    size_t max_entries = 0;

    for ( size_t i = 0; i < rows(); i++ )
        max_entries = std::max( size_t( _rowptr[i+1] - _rowptr[i] ), max_entries );

    return max_entries;
}

//
// return true if matrix contains zero or elements < ε on diagonal
//
template < typename value_t >
bool
TSparseMatrix< value_t >::has_diag_zero ( const real_t  eps )
{
    const size_t  nrows = rows();
    const size_t  ncols = cols();

    for ( idx_t  i = 0; i < std::min( idx_t( nrows ), idx_t( ncols ) ); ++i )
    {
        if ( ! has_entry( i, i ) || ( Math::abs( entry( i, i ) ) < eps ))
            return true;
    }// for

    return false;
}

//
// return true if S is (weakly) diagonally dominant
//

namespace
{

//
// return ∑_j |S_ij|, j≠i
//
template < typename value_t >
real_type_t< value_t >
comp_off_diag ( const idx_t            i,
                const TSparseMatrix< value_t > *  S )
{
    using  real_t = typename TSparseMatrix< value_t >::real_t;
    
    real_t       retval = real_t(0);
    const idx_t  lb     = S->rowptr( i );
    const idx_t  ub     = S->rowptr( i+1 );

    for ( idx_t  j = lb; j < ub; j++ )
    {
        if ( S->colind( j ) != i )
            retval += Math::abs( S->coeff( j ) );
    }// for

    return retval;
}

}// namespace anonymous

template < typename value_t >
bool
TSparseMatrix< value_t >::is_diag_dom ( const bool  weak )
{
    if ( weak )
    {
        for ( auto  i : this->row_is() )
        {
            const real_t  off_diag = comp_off_diag( i, this );
            const real_t  diag     = Math::abs( entry( i, i ) );
            
            if ( diag < off_diag )
                return false;
        }// for
    }// if
    else
    {
        for ( auto  i : this->row_is() )
        {
            const real_t  off_diag = comp_off_diag( i, this );
            const real_t  diag     = Math::abs( entry( i, i ) );
            
            if ( diag <= off_diag )
                return false;
        }// for
    }// else

    return true;
}
    
//
// return factor α, such that S + α·I is diagonally dominant (this = S)
//
template < typename value_t >
typename TSparseMatrix< value_t >::real_t
TSparseMatrix< value_t >::diag_dom_factor ()
{
    real_t  alpha = real_t(0);
    
    for ( auto  i : this->row_is() )
    {
        const real_t  off_diag = comp_off_diag( i, this );
        const real_t  diag     = Math::abs( entry( i, i ) );

        if ( diag < off_diag )
            alpha = std::max( alpha, off_diag - diag );
    }// for

    return alpha;
}

//
// some tests on the matrix
//
template < typename value_t >
void
TSparseMatrix< value_t >::check_matrix () const
{
    //
    // print number of entries per row
    //

    size_t  min_entries = 0;
    size_t  max_entries = 0;
    size_t  min_pos     = 0;
    size_t  max_pos     = 0;
    
    for ( size_t i = 0; i < rows(); i++ )
    {
        size_t  nentries = _rowptr[i+1] - _rowptr[i];
        
        if (( min_entries == 0 ) || ( nentries < min_entries ))
        {
            min_pos     = i;
            min_entries = nentries;
        }// if
        
        if ( nentries > max_entries )
        {
            max_pos     = i;
            max_entries = nentries;
        }// if
    }// for

    std::cout << "(TSparseMatrix) check_matrix : "
              << "minimal/maximal/average degree per row = "
              << min_entries << "(at " << min_pos << ") / "
              << max_entries << "(at " << max_pos << ") / "
              << avg_entries_per_row()
              << std::endl;

    //
    // look for zero rows
    //

    size_t  n_zero_row = 0;
    
    for ( size_t row = 0; row < rows(); row++ )
    {
        const idx_t  lb      = _colind[row];
        const idx_t  ub      = _colind[row+1];
        value_t      row_sum = value_t(0);

        for ( idx_t j = lb; j < ub; j++ )
        {
            row_sum += _coeff[j];
        }// for

        if ( Math::abs( row_sum ) < 1e-8 )
        {
            n_zero_row++;
        }// if
    }// for

    std::cout << "(TSparseMatrix) check_matrix : "
              << n_zero_row << " rows with (nearly) zero row-sum"
              << std::endl;
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
TSparseMatrix< value_t >::scale ( const value_t  f )
{
    if ( f == value_t(1) )
        return;

    for ( size_t i = 0; i < _coeff.size(); i++ )
        _coeff[i] *= f;
}
    
//
// matrix-vector-multiplication
//
template < typename value_t >
void
TSparseMatrix< value_t >::mul_vec ( const value_t               alpha,
                                    const TVector< value_t > *  x,
                                    const value_t               beta,
                                    TVector< value_t > *        y,
                                    const matop_t               op ) const
{
    //
    // check for nullptr, index set and type
    //
    
    if ( x == nullptr ) HERROR( ERR_ARG, "(TSparseMatrix) mul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(TSparseMatrix) mul_vec", "y = nullptr" );

    if ( op == MATOP_NORM )
    {
        if (( this->row_is() != y->is() ) || ( this->col_is() != x->is() ))
            HERROR( ERR_VEC_STRUCT, "(TSparseMatrix) mul_vec", "incompatible vector index set" );
    }// if
    else
    {
        if (( this->col_is() != y->is() ) || ( this->row_is() != x->is() ))
            HERROR( ERR_VEC_STRUCT, "(TSparseMatrix) mul_vec", "incompatible vector index set" );
    }// if
    
    if ( ! IS_TYPE( x, TScalarVector ) || ! IS_TYPE( y, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TSparseMatrix) mul_vec",
                y->typestr() + " += TSparseMatrix * " + x->typestr() );

    //
    // apply β
    //

    if      ( beta == value_t(0) ) y->fill( value_t(0) );
    else if ( beta != value_t(1) ) y->scale( beta );

    //
    // stop, if nothing to do
    //
    
    if (( alpha == value_t(0) ) || ( n_non_zero() == 0 ))
        return;

    //
    // entry-wise multiplication
    //

    auto        sx = cptrcast( x, TScalarVector< value_t > );
    auto        sy =  ptrcast( y, TScalarVector< value_t > );
    const auto  n  = rows();
            
    if ( op == apply_normal )
    {
        //
        // normal matrix
        //
                
        for ( idx_t i = 0; i < idx_t(n); i++ )
        {
            const idx_t  lb = _rowptr[i];
            const idx_t  ub = _rowptr[i+1];
            value_t      s  = value_t(0);
                    
            for ( idx_t j = lb; j < ub; j++ )
                s += _coeff[j] * sx->blas_vec()( _colind[j] );

            sy->blas_vec()( i ) += alpha * s;
        }// for
    }// if
    else if ( op == apply_transposed )
    {
        //
        // transposed matrix
        //
                
        for ( idx_t  i = 0; i < idx_t(n); i++ )
        {
            const idx_t    lb = _rowptr[i];
            const idx_t    ub = _rowptr[i+1];
            const value_t  a  = alpha * sx->blas_vec()( i );
                
            for ( idx_t j = lb; j < ub; j++ )
                sy->blas_vec()( _colind[j] ) += _coeff[j] * a;
        }// for
    }// if
    else if ( op == MATOP_ADJ )
    {
        //
        // adjoint matrix
        //
                
        for ( idx_t  i = 0; i < idx_t(n); i++ )
        {
            const idx_t    lb = _rowptr[i];
            const idx_t    ub = _rowptr[i+1];
            const value_t  a  = alpha * sx->blas_vec()( i );

            for ( idx_t j = lb; j < ub; j++ )
                sy->blas_vec()( _colind[j] ) += Math::conj( _coeff[j] ) * a;
        }// for
    }// if
}

//
// compute this = this + a * M
//
template < typename value_t >
void
TSparseMatrix< value_t >::add ( const value_t               a,
                                const TMatrix< value_t > *  M )
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "(TSparseMatrix) add", "matrix is nullptr" );

    if ( ! IS_TYPE( M, TSparseMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TSparseMatrix) add", M->typestr() );

    if ( a == value_t(0) )
        return;
    
    auto  S = cptrcast( M, TSparseMatrix< value_t > );

    if (( _coeff.size()  != S->_coeff.size() ) ||
        ( _colind.size() != S->_colind.size() ))
        HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) add", "argument has different CRS format" );
    
    //
    // check for same structure
    //
    
    for ( size_t i = 0; i < _rowptr.size(); i++ )
        if ( _rowptr[i] != S->_rowptr[i] )
            HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) add", "argument has different CRS format" );
    
    for ( size_t i = 0; i < _colind.size(); i++ )
        if ( _colind[i] != S->_colind[i] )
            HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) add", "argument has different CRS format" );
    
    // finally: same format, just add entries
    // (assuming same ordering !!!)
    for ( size_t  i = 0; i < _coeff.size(); i++ )
        _coeff[i] += a * S->_coeff[i];
}

/////////////////////////////////////////////////
//
// input/output related
//

//
// output of matrix
//
template < typename value_t >
void
TSparseMatrix< value_t >::print ( const uint ofs ) const
{
    TMatrix< value_t >::print( ofs );
}

//
// serialisation
//

template < typename value_t >
void
TSparseMatrix< value_t >::read  ( TByteStream & s )
{
    return build( s );
}

template < typename value_t >
void
TSparseMatrix< value_t >::build ( TByteStream & s )
{
    TMatrix< value_t >::read( s );

    s.get( & _rows, sizeof(size_t) );
    s.get( & _cols, sizeof(size_t) );

    //
    // read CRS data
    //
    
    size_t  nnz;
    
    s.get( & nnz, sizeof(size_t) );
    
    init( nnz );
    
    s.get( _rowptr.data(), _rows+1 * sizeof(idx_t) );
    s.get( _colind.data(), nnz * sizeof(idx_t) );
    
    s.get( _coeff.data(), nnz * sizeof(value_t) );
}

template < typename value_t >
void
TSparseMatrix< value_t >::write ( TByteStream & s ) const
{
    TMatrix< value_t >::write( s );

    s.put( & _rows, sizeof(size_t) );
    s.put( & _cols, sizeof(size_t) );

    const size_t  nnz  = _colind.size();
        
    s.put( & nnz, sizeof(nnz) );
    s.put( const_cast< std::vector< idx_t > & >(_rowptr).data(), sizeof(size_t) * _rowptr.size() );
    s.put( const_cast< std::vector< idx_t > & >(_colind).data(), sizeof(size_t) * _colind.size() );
    s.put( const_cast< std::vector< value_t > & >(_coeff).data(), sizeof(value_t) * _coeff.size() );
}

//
// returns size of object in bytestream
//
template < typename value_t >
size_t
TSparseMatrix< value_t >::bs_size () const
{
    size_t  size = TMatrix< value_t >::bs_size() + 2*sizeof(size_t);

    size += sizeof(size_t);
    size += _rowptr.size() * sizeof(size_t);
    size += _colind.size() * sizeof(size_t);
    size += _coeff.size() * sizeof(value_t);

    return size;
}

/////////////////////////////////////////////////
//
// misc.
//

//
// transpose matrix
//
template < typename value_t >
void
TSparseMatrix< value_t >::transpose ()
{
    // to be done
    HERROR( ERR_NOT_IMPL, "(TSparseMatrix) transpose", "" );

    TMatrix< value_t >::transpose();
}
    
//
// conjugate matrix coefficients
//
template < typename value_t >
void
TSparseMatrix< value_t >::conjugate ()
{
    if ( this->is_complex() && ! this->is_hermitian() )
    {
        for ( size_t  i = 0; i < _coeff.size(); ++i )
            _coeff[i] = Math::conj( _coeff[i] );
    }// if
}

//
// restrict sparse matrix to given block index set
//
template < typename value_t >
unique_ptr< TSparseMatrix< value_t > >
TSparseMatrix< value_t >::restrict ( const TIndexSet &  rowis,
                                     const TIndexSet &  colis ) const
{
    if ( ! this->block_is().is_sub( TBlockIndexSet( rowis, colis ) ) )
        HERROR( ERR_MAT_STRUCT, "TSparseMatrix) restrict", "given block index set is not a local sub set" );
                
    auto  T = make_unique< TSparseMatrix< value_t > >( rowis.size(), colis.size() );

    T->set_block_is( TBlockIndexSet( rowis, colis ) );
    
    //
    // determine total number of nonzeroes and number per row in subblock
    //

    size_t                sub_nnzero = 0;
    std::vector< idx_t >  sub_rowptr( rowis.size()+1 );

    for ( auto  row : rowis )
    {
        const idx_t  lb = _rowptr[ row - this->row_ofs() ];
        const idx_t  ub = _rowptr[ row - this->row_ofs() + 1 ];

        sub_rowptr[ row - rowis.first() ] = idx_t(sub_nnzero);
        
        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col = _colind[j] + this->col_ofs();
            
            if ( colis.is_in( col ) )
                sub_nnzero++;
        }// for
    }// for
    sub_rowptr[ rowis.size() ] = idx_t(sub_nnzero);

    T->init( sub_nnzero );

    //
    // set up CRS data in sub matrix
    //

    for ( idx_t  i = 0; i < idx_t(sub_rowptr.size()); ++i )
        T->rowptr( i ) = sub_rowptr[i];

    idx_t  pos = 0;
    
    for ( auto  row : rowis )
    {
        const idx_t  lb = _rowptr[ row - this->row_ofs() ];
        const idx_t  ub = _rowptr[ row - this->row_ofs() + 1 ];

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col = _colind[j] + this->col_ofs();
            
            if ( colis.is_in( col ) )
            {
                T->colind( pos ) = col - colis.first();
                T->coeff( pos )  = _coeff[j];
                
                pos++;
            }// if
        }// for
    }// for

    return T;
}

//
// restrict sparse matrix to given block index set
// which is given in a different ordering
//
template < typename value_t >
unique_ptr< TSparseMatrix< value_t > >
TSparseMatrix< value_t >::restrict ( const TIndexSet &     rowis,
                                     const TPermutation *  rowperm,
                                     const TIndexSet &     colis,
                                     const TPermutation *  colperm ) const
{
    if ( ! this->block_is().is_sub( TBlockIndexSet( rowis, colis ) ) )
        HERROR( ERR_MAT_STRUCT, "TSparseMatrix) restrict", "given block index set is not a local sub set" );
                
    auto  T = make_unique< TSparseMatrix< value_t > >( rowis.size(), colis.size() );

    T->set_block_is( TBlockIndexSet( rowis, colis ) );
    
    //
    // determine total number of nonzeroes and number per row in subblock
    //

    size_t                sub_nnzero = 0;
    std::vector< idx_t >  sub_rowptr( rowis.size()+1 );

    for ( auto  row : rowis )
    {
        const idx_t  prow = ( rowperm != nullptr ? rowperm->permute( row ) : row );
        const idx_t  lb   = _rowptr[ prow - this->row_ofs() ];
        const idx_t  ub   = _rowptr[ prow - this->row_ofs() + 1 ];

        sub_rowptr[ row - rowis.first() ] = idx_t(sub_nnzero);
        
        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col  = _colind[j] + this->col_ofs();
            const idx_t  pcol = ( colperm != nullptr ? colperm->permute( col ) : col );
            
            if ( colis.is_in( pcol ) )
                sub_nnzero++;
        }// for
    }// for
    sub_rowptr[ rowis.size() ] = idx_t(sub_nnzero);

    T->init( sub_nnzero );

    //
    // set up CRS data in sub matrix
    //

    for ( idx_t  i = 0; i < idx_t(sub_rowptr.size()); ++i )
        T->rowptr( i ) = sub_rowptr[i];

    idx_t  pos = 0;
    
    for ( auto  row : rowis )
    {
        const idx_t  prow = ( rowperm != nullptr ? rowperm->permute( row ) : row );
        const idx_t  lb   = _rowptr[ prow - this->row_ofs() ];
        const idx_t  ub   = _rowptr[ prow - this->row_ofs() + 1 ];

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col  = _colind[j] + this->col_ofs();
            const idx_t  pcol = ( colperm != nullptr ? colperm->permute( col ) : col );
            
            if ( colis.is_in( pcol ) )
            {
                T->colind( pos ) = pcol - colis.first();
                T->coeff( pos )  = _coeff[j];

                pos++;
            }// if
        }// for
    }// for

    return T;
}

//
// return number of coefficients in sub block index set \a rowis × \a colis
//
template < typename value_t >
size_t
TSparseMatrix< value_t >::restrict_nonzeroes ( const TIndexSet &  rowis,
                                               const TIndexSet &  colis ) const
{
    if ( ! this->block_is().is_sub( TBlockIndexSet( rowis, colis ) ) )
        HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) restrict_nonzeroes", "given block index set is not a local sub set" );
                
    //
    // determine total number of nonzeroes and number per row in subblock
    //

    size_t  sub_nnzero = 0;

    for ( auto  row : rowis )
    {
        const idx_t  lb = _rowptr[ row - this->row_ofs() ];
        const idx_t  ub = _rowptr[ row - this->row_ofs() + 1 ];

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col = _colind[j] + this->col_ofs();
            
            if ( colis.is_in( col ) )
                sub_nnzero++;
        }// for
    }// for

    return sub_nnzero;
}

template < typename value_t >
size_t
TSparseMatrix< value_t >::restrict_nonzeroes ( const TIndexSet &     rowis,
                                               const TPermutation *  rowperm,
                                               const TIndexSet &     colis,
                                               const TPermutation *  colperm ) const
{
    if ( ! this->block_is().is_sub( TBlockIndexSet( rowis, colis ) ) )
        HERROR( ERR_MAT_STRUCT, "TSparseMatrix) restrict", "given block index set is not a local sub set" );
                
    //
    // determine total number of nonzeroes and number per row in subblock
    //

    size_t  sub_nnzero = 0;

    for ( auto  row : rowis )
    {
        const idx_t  prow = ( rowperm != nullptr ? rowperm->permute( row ) : row );
        const idx_t  lb   = _rowptr[ prow - this->row_ofs() ];
        const idx_t  ub   = _rowptr[ prow - this->row_ofs() + 1 ];

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col  = _colind[j] + this->col_ofs();
            const idx_t  pcol = ( colperm != nullptr ? colperm->permute( col ) : col );
            
            if ( colis.is_in( pcol ) )
                sub_nnzero++;
        }// for
    }// for

    return sub_nnzero;
}

//
// virtual constructors
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TSparseMatrix< value_t >::copy () const
{
    auto  M = TMatrix< value_t >::copy();
    auto  S = ptrcast( M.get(), TSparseMatrix< value_t > );

    if ( this->cluster() != nullptr )
        S->set_cluster( this->cluster() );
    else
        S->set_block_is( this->block_is() );

    const size_t  nnz = _colind.size();

    S->init( nnz );

    for ( size_t i = 0; i <= rows(); i++ )
        S->_rowptr[i] = _rowptr[i];

    for ( size_t i = 0; i < nnz; i++ )
    {
        S->_colind[i] = _colind[i];
        S->_coeff[i]  = _coeff[i];
    }// for

    return M;
}

//
// copy matrix into A
//
template < typename value_t >
void
TSparseMatrix< value_t >::copy_to ( TMatrix< value_t > * A ) const
{
    TMatrix< value_t >::copy_to( A );
    
    if ( ! IS_TYPE( A, TSparseMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TSparseMatrix) copy_to", A->typestr() );
    
    auto  M = ptrcast( A, TSparseMatrix< value_t > );

    if ( this->cluster() != nullptr )
        M->set_cluster( this->cluster() );
    else
        M->set_block_is( this->block_is() );

    const size_t  nnz = _colind.size();
        
    M->init( nnz );

    for ( size_t i = 0; i <= rows(); i++ )
        M->_rowptr[i] = _rowptr[i];
    
    for ( size_t i = 0; i < nnz; i++ )
    {
        M->_colind[i] = _colind[i];
        M->_coeff[i]  = _coeff[i];
    }// for
}

//
// return structural copy
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TSparseMatrix< value_t >::copy_struct () const
{
    auto  M = TMatrix< value_t >::copy_struct();
    auto  S = ptrcast( M.get(), TSparseMatrix< value_t > );

    if ( this->cluster() != nullptr )
        S->set_cluster( this->cluster() );
    else
        S->set_block_is( this->block_is() );

    S->init( 0 );

    return M;
}

//
// return size in bytes used by this object
//
template < typename value_t >
size_t
TSparseMatrix< value_t >::byte_size () const
{
    size_t  size = TMatrix< value_t >::byte_size() + 2*sizeof(size_t);

    size += _coeff.size()  * sizeof(value_t) + sizeof(_coeff);
    size += _colind.size() * sizeof(idx_t)   + sizeof(_colind);
    size += _rowptr.size() * sizeof(idx_t)   + sizeof(_rowptr);
    
    return size;
}

//
// test data for invalid values, e.g. INF and NAN
//
template < typename value_t >
void
TSparseMatrix< value_t >::check_data () const
{
    for ( auto  val : _coeff )
    {
        if ( Math::is_inf( val ) )
            HERROR( ERR_INF, "(TSparseMatrix) check_data", "" );
        
        if ( Math::is_nan( val ) )
            HERROR( ERR_NAN, "(TSparseMatrix) check_data", "" );
    }// for
}

//
// print histogram for entries per row in GNUPLOT format
//
template < typename value_t >
void
TSparseMatrix< value_t >::print_pattern_hist ( std::ostream & os ) const
{
    //
    // compute histogram
    //
    
    const size_t           max_deg = max_entries_per_row();
    std::vector< size_t >  hist( max_deg + 1 );

    for ( idx_t  i = 0; i < idx_t( rows() ); ++i )
    {
        hist[ rowptr(i+1) - rowptr(i) ]++;
    }// for

    //
    // print to stream
    //
    
    os << "# See GNUPLOT manual for a description of the format" << std::endl
       << std::endl
       << "# <number of coeffs per row>  <number of rows in matrix>" << std::endl;

    for ( size_t  i = 0; i <= max_deg; ++i )
        os << i << " " << hist[i] << std::endl;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

template class TSparseMatrix< float >;
template class TSparseMatrix< double >;
template class TSparseMatrix< std::complex< float > >;
template class TSparseMatrix< std::complex< double > >;

#define INST_IMPORT_CRS( T_idx, T_val )                                 \
    template void                                                       \
    TSparseMatrix< T_val >::import_crs< T_idx > ( const size_t   nrows, \
                                                  const size_t   ncols, \
                                                  const size_t   nnonzero, \
                                                  const T_idx *  rowptr, \
                                                  const T_idx *  colind, \
                                                  const T_val *  coeffs )

INST_IMPORT_CRS( int,   float );
INST_IMPORT_CRS( int,   double );
INST_IMPORT_CRS( int,   std::complex< float > );
INST_IMPORT_CRS( int,   std::complex< double > );
INST_IMPORT_CRS( int64_t, float );
INST_IMPORT_CRS( int64_t, double );
INST_IMPORT_CRS( int64_t, std::complex< float > );
INST_IMPORT_CRS( int64_t, std::complex< double > );
INST_IMPORT_CRS( long long int, float );
INST_IMPORT_CRS( long long int, double );
INST_IMPORT_CRS( long long int, std::complex< float > );
INST_IMPORT_CRS( long long int, std::complex< double > );

#define INST_IMPORT_CCS( T_idx, T_val )                                 \
    template void                                                       \
    TSparseMatrix< T_val >::import_ccs< T_idx > ( const size_t   nrows, \
                                                  const size_t   ncols, \
                                                  const size_t   nnonzero, \
                                                  const T_idx *  colptr, \
                                                  const T_idx *  rowind, \
                                                  const T_val *  coeffs )

INST_IMPORT_CCS( int,   float );
INST_IMPORT_CCS( int,   double );
INST_IMPORT_CCS( int,   std::complex< float > );
INST_IMPORT_CCS( int,   std::complex< double > );
INST_IMPORT_CCS( int64_t, float );
INST_IMPORT_CCS( int64_t, double );
INST_IMPORT_CCS( int64_t, std::complex< float > );
INST_IMPORT_CCS( int64_t, std::complex< double > );
INST_IMPORT_CCS( long long int, float );
INST_IMPORT_CCS( long long int, double );
INST_IMPORT_CCS( long long int, std::complex< float > );
INST_IMPORT_CCS( long long int, std::complex< double > );

#define INST_EXPORT_CCS( T_idx, T_val )                                 \
    template void                                                       \
    TSparseMatrix< T_val >::export_ccs< T_idx > ( std::vector< T_idx > &  colptr, \
                                                  std::vector< T_idx > &  rowind, \
                                                  std::vector< T_val > &  coeffs, \
                                                  const bool              use_sym ) const

INST_EXPORT_CCS( int,   float );
INST_EXPORT_CCS( int,   double );
INST_EXPORT_CCS( int,   std::complex< float > );
INST_EXPORT_CCS( int,   std::complex< double > );
INST_EXPORT_CCS( int64_t, float );
INST_EXPORT_CCS( int64_t, double );
INST_EXPORT_CCS( int64_t, std::complex< float > );
INST_EXPORT_CCS( int64_t, std::complex< double > );
INST_EXPORT_CCS( long long int, float );
INST_EXPORT_CCS( long long int, double );
INST_EXPORT_CCS( long long int, std::complex< float > );
INST_EXPORT_CCS( long long int, std::complex< double > );

}// namespace Hpro
