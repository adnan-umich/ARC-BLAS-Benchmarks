#ifndef __HLIB_BLAS_HH
#define __HLIB_BLAS_HH
//
// Project     : HLib
// File        : hlib-blas.hh
// Description : HLIBpro include file containing basic linear algebra related headers
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//
//! \defgroup BLAS_Module BLAS
//!
//! This modules provides most low level algebra functions, e.g.
//! vector dot products, matrix multiplication, factorisation and
//! singular value decomposition. See also \ref TUTORBLAS for an introduction.
//!
//! To include all BLAS algebra functions and classes add
//! \code
//! #include <hlib-blas.hh>
//! \endcode
//! to your source files.
//!

#include <hpro/blas/Vector.hh>
#include <hpro/blas/Matrix.hh>
#include <hpro/blas/Algebra.hh>
#include <hpro/blas/Range.hh>
#include <hpro/blas/gpu.hh>
#include <hpro/blas/matrix_view.hh>
#include <hpro/blas/test.hh>

#endif // __HLIB_BLAS_HH
