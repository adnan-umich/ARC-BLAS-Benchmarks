#ifndef __HPRO_BLAS_CUDA_HH
#define __HPRO_BLAS_CUDA_HH
//
// Project     : HLIBpro
// File        : blas/cuda.hh
// Description : CUDA related functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2024. All Rights Reserved.
//

#include <hpro/blas/types.hh>
#include <hpro/blas/Vector.hh>
#include <hpro/blas/Matrix.hh>

namespace Hpro { namespace CUDA {

////////////////////////////////////////////////////////////////////////////////
//
// management functions
//
////////////////////////////////////////////////////////////////////////////////

//! (index) type for cuda handles
using cuda_handle_t = int;

//! signals invalid cuda handle
constexpr int  INVALID_HANDLE = -1;

//!
//! \ingroup  BLAS_Module
//! \brief    initialize CUDA related data structures
//!
void
init ();

//!
//! \ingroup  BLAS_Module
//! \brief    finish CUDA related data structures
//!
void
done ();

//!
//! \ingroup  BLAS_Module
//! \brief    get free handle for CUDA (stream/cuBlas/cuSolver)
//!           - if no free handle exists, INVALID_HANDLE is returned
//!
cuda_handle_t
request_handle ();

//!
//! \ingroup  BLAS_Module
//! \brief    release previously aquired handle
//!
void
release_handle ( const cuda_handle_t  handle );

////////////////////////////////////////////////////////////////////////////////
//
// algebra functions
//
////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup  BLAS_Module
//! \brief   Compute QR factorisation of the n×m matrix \a A.
//!          On exit \a A contains Q
//! \return  \c true if QR was computed successfully and \c false otherwise
//!
template < typename value_t >
bool
qr  ( const cuda_handle_t        handle,
      BLAS::Matrix< value_t > &  A,
      BLAS::Matrix< value_t > &  R );

//!
//! \ingroup  BLAS_Module
//! \brief    compute SVD decomposition \f$ A = U·S·V^H \f$ of the nxm matrix \a A with
//!           n×min(n,m) matrix U, min(n,m)×min(n,m) matrix S (diagonal)
//!           and m×min(n,m) matrix V; \a A will be overwritten with U upon exit
//! \return   \c true if SVD was computed successfully and \c false otherwise
//!
template < typename value_t >
bool
svd  ( const cuda_handle_t                       handle,
       BLAS::Matrix< value_t > &                 A,
       BLAS::Vector< real_type_t< value_t > > &  S,
       BLAS::Matrix< value_t > &                 V );

}}// namespace Hpro::CUDA

#endif // __HPRO_BLAS_CUDA_HH
