#ifndef __HPRO_BLAS_MATRIXBASE_HH
#define __HPRO_BLAS_MATRIXBASE_HH
//
// Project     : HLIBpro
// File        : MatrixBase.hh
// Description : base class for all matrices and matrix views
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/blas/Range.hh"

namespace Hpro
{

namespace BLAS
{

//
// trait for giving access to matrix properties
//
template < typename T_matrix > struct matrix_trait;

//
// signals, that T is of matrix type
//
template < typename T >
struct is_matrix
{
    static const bool  value = false;
};

#if __cplusplus >= 201703L
template <typename T> inline constexpr bool is_matrix_v = is_matrix< T >::value;
#endif

//!
//! different matrix views in BLAS
//!
enum blasview_t
{
    BLAS_NORMAL     = 'N',
    BLAS_TRANSPOSED = 'T',
    BLAS_ADJOINT    = 'C'
};

//!
//! \ingroup  BLAS_Module
//! \class    MatrixBase
//! \brief    defines basic interface for matrices
//!
template < typename T_derived >
class MatrixBase
{
public:
    //! scalar value type of matrix
    using  value_t = typename matrix_trait< T_derived >::value_t;

public:
    //
    // data access
    //

    //! return number of rows of matrix
    size_t       nrows        () const noexcept { return derived().nrows(); }
    
    //! return number of columns of matrix
    size_t       ncols        () const noexcept { return derived().ncols(); }

    //! return number of rows of matrix
    Range        row_range    () const noexcept { return Range( 0, idx_t(nrows())-1 ); }
    
    //! return number of columns of matrix
    Range        col_range    () const noexcept { return Range( 0, idx_t(ncols())-1 ); }

    //! return coefficient (i,j)
    value_t      operator ()  ( const idx_t i, const idx_t j ) const noexcept
    {
        return derived()(i,j);
    }

    //! return reference to coefficient (i,j)
    value_t &    operator ()  ( const idx_t i, const idx_t j ) noexcept
    {
        return derived()(i,j);
    }

    //! return pointer to internal data
    value_t *    data         () const noexcept { return derived().data(); }

    //! return stride w.r.t. row index set
    size_t       row_stride   () const noexcept { return derived().row_stride(); }

    //! return stride w.r.t. column index set
    size_t       col_stride   () const noexcept { return derived().col_stride(); }

    //! return BLAS matrix view of matrix object
    blasview_t   blas_view    () const noexcept { return derived().blas_view(); }

    //! return number of rows of actual BLAS matrix
    size_t       blas_nrows   () const noexcept { return derived().blas_nrows(); }

    //! return number of columns of actual BLAS matrix
    size_t       blas_ncols   () const noexcept { return derived().blas_ncols(); }

private:
    //! convert to derived type
    T_derived &        derived  ()       noexcept { return * static_cast<       T_derived * >( this ); }
    const T_derived &  derived  () const noexcept { return * static_cast< const T_derived * >( this ); }
};

//
// signals, that T is of matrix type
//
template < typename T >
struct is_matrix< MatrixBase< T > >
{
    static const bool  value = is_matrix< T >::value;
};


//! stream output for matrices
template < typename T >
std::ostream & operator << ( std::ostream & os, const MatrixBase< T > & M );

}// namespace BLAS

}// namespace Hpro

#endif  // __HPRO_BLAS_MATRIXBASE_HH
