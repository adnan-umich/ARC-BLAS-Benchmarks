#ifndef __HPRO_BLAS_LAPACK_HH
#define __HPRO_BLAS_LAPACK_HH
//
// Project     : HLIBpro
// File        : lapack.hh
// Description : definition of LAPACK functions in c-format
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <atomic>
#include <cstdint>

#include "hpro/config.h"

#if HPRO_USE_MKL == 1
#include <mkl_cblas.h>
#endif

#include "hpro/base/traits.hh"
#include "hpro/base/types.hh"

#include "hpro/blas/flops.hh"

namespace Hpro
{

namespace BLAS
{

//
// FLOP counting
//

extern std::atomic< std::uint64_t >  FLOPS;

#if HPRO_COUNT_FLOPS == 1
#  define ADD_FLOPS( n )  Hpro::BLAS::FLOPS += (n)
#else
#  define ADD_FLOPS( n )
#endif

//
// return flop numbers
//
inline
uint64_t
get_flops ()
{
    return FLOPS.load();
}

//
// reset flop numbers
//
inline
void
reset_flops ()
{
    FLOPS = 0;
}

}// namespace blas

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// BLAS integer type
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

// define 32-bit integers vs. 64-bit integers
#if HPRO_USE_ILP64 == 1
using  blas_int_t = int64_t;   // ILP64
#else
using  blas_int_t = int32_t;   // LP64
#endif

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// declaration of external BLAS functions
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

//
// Win32 DLL export declaration
///
#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
// #  define  CFUNCDECL __declspec(dllimport)
#  define  CFUNCDECL
#else
#  define  CFUNCDECL
#endif

extern "C" {

///////////////////////////////////////////////////////////////////
//
// real-valued functions (single precision)
//
///////////////////////////////////////////////////////////////////

// set (part of) A to alpha/beta
CFUNCDECL
void
slaset_  ( const char *         UPLO,
           const blas_int_t *   M,
           const blas_int_t *   N,
           const float *        ALPHA,
           const float *        BETA,
           float *              A,
           const blas_int_t *   LDA   );

// copy (part of) A to B
CFUNCDECL
void
slacpy_  ( const char *         UPLO,
           const blas_int_t *   M,
           const blas_int_t *   N,
           const float *        A,
           const blas_int_t *   LDA,
           float *              B,
           const blas_int_t *   LDB  );

// compute dotproduct between x and y
CFUNCDECL
float
sdot_ ( const blas_int_t *      n,
        const float *           dx,
        const blas_int_t *      incx,
        const float *           dy,
        const blas_int_t *      incy);

// copy vector x into vector y
CFUNCDECL
void
scopy_ ( const blas_int_t *     n,
         const float *          dx,
         const blas_int_t *     incx,
         float       *          dy,
         const blas_int_t *     incy);

// compute y = y + a * x
CFUNCDECL
void
saxpy_ ( const blas_int_t *     n,
         const float *          da,
         const float *          dx,
         const blas_int_t *     incx,
         float *                dy,
         const blas_int_t *     incy );

// compute sum of absolut values of x
CFUNCDECL
float
sasum_  ( const blas_int_t *    n,
          const float *         dx,
          const blas_int_t *    incx );

// return euclidean norm of x
CFUNCDECL
float
snrm2_  ( const blas_int_t *    n,
          const float *         dx,
          const blas_int_t *    incx );

// scale x by a
CFUNCDECL
void
sscal_ ( const blas_int_t *     n,
         const float *          da,
         float *                dx,
         const blas_int_t *     incx );

// interchange x and y
CFUNCDECL
void
sswap_ ( const blas_int_t *     n,
         float *                dx,
         const blas_int_t *     incx,
         float *                dy,
         const blas_int_t *     incy );

// finds index with element having maximal absolut value
CFUNCDECL
blas_int_t
isamax_ ( const blas_int_t *    n,
          const float   *,
          const blas_int_t * );

// y = alpha A x + beta y
CFUNCDECL
void
sgemv_ ( const char *           trans,
         const blas_int_t *     M,
         const blas_int_t *     N,
         const float *          alpha,
         const float *          A,
         const blas_int_t *     lda,
         const float *          dx,
         const blas_int_t *     incx,
         const float *          beta,
         float *                dy,
         const blas_int_t *     incy );

// compute c = alpha * a * b + beta * c
CFUNCDECL
void
sgemm_ ( const char *           transa,
         const char *           transb,
         const blas_int_t *     n,
         const blas_int_t *     m,
         const blas_int_t *     k,
         const float *          alpha,
         const float *          a,
         const blas_int_t *     lda,
         const float *          b,
         const blas_int_t *     ldb,
         const float *          beta,
         float *                c,
         const blas_int_t *     ldc );

// computes y = op(A)*x for upper/lower triangular A (y overwrites x)
blas_int_t
strmv_ ( const char *           uplo,
         const char *           trans,
         const char *           diag,
         const blas_int_t *     n,
         const float *          A,
         const blas_int_t *     ldA,
         float *                x,
         const blas_int_t *     incx );

// computes B = alpha * op(A) * B or B = alpha * B * op(A)
CFUNCDECL
void
strmm_ ( const char *           side,
         const char *           uplo,
         const char *           transa,
         const char *           diag,
         const blas_int_t *     n,
         const blas_int_t *     m,
         const float *          alpha,
         const float *          A,
         const blas_int_t *     lda,
         float *                B,
         const blas_int_t *     ldb );

// solves A x = b or A^T x = b with triangular A
CFUNCDECL
void
strsv_ ( const char *           uplo,
         const char *           trans,
         const char *           diag,
         const blas_int_t *     n,
         const float *          A,
         const blas_int_t *     lda,
         float *                x,
         const blas_int_t *     incx );

// solves A X = B or A^T X = B with triangular A
CFUNCDECL
void
strsm_ ( const char *           side,
         const char *           uplo,
         const char *           trans,
         const char *           diag,
         const blas_int_t *     n,
         const blas_int_t *     m,
         const float *          alpha,
         const float *          A,
         const blas_int_t *     lda,
         float *                B,
         const blas_int_t *     ldb );

// A = alpha * x * y^T + A, A \in \R^{m x n}
CFUNCDECL
void
sger_ ( const blas_int_t *      m,
        const blas_int_t *      n,
        const float *           alpha,
        const float *           x,
        const blas_int_t *      incx,
        const float *           y,
        const blas_int_t *      incy,
        float *                 A,
        const blas_int_t *      lda );

///////////////////////////////////////////////////////////////////
//
// real-valued functions (double precision)
//
///////////////////////////////////////////////////////////////////

// set (part of) A to alpha/beta
CFUNCDECL
void
dlaset_  ( const char *         UPLO,
           const blas_int_t *   M,
           const blas_int_t *   N,
           const double *       ALPHA,
           const double *       BETA,
           double *             A,
           const blas_int_t *   LDA   );

// copy (part of) A to B
CFUNCDECL
void
dlacpy_  ( const char *         UPLO,
           const blas_int_t *   M,
           const blas_int_t *   N,
           const double *       A,
           const blas_int_t *   LDA,
           double *             B,
           const blas_int_t *   LDB  );

// compute dotproduct between x and y
CFUNCDECL
double
ddot_ ( const blas_int_t *      n,
        const double *          dx,
        const blas_int_t *      incx,
        const double *          dy,
        const blas_int_t *      incy);
    
// copy vector x into vector y
CFUNCDECL
void
dcopy_ ( const blas_int_t *     n,
         const double *         dx,
         const blas_int_t *     incx,
         double       *         dy,
         const blas_int_t *     incy);

// compute y = y + a * x
CFUNCDECL
void
daxpy_ ( const blas_int_t *     n,
         const double *         da,
         const double *         dx,
         const blas_int_t *     incx,
         double *               dy,
         const blas_int_t *     incy );

// compute sum of absolut values of x
CFUNCDECL
double
dasum_  ( const blas_int_t *    n,
          const double *        dx,
          const blas_int_t *    incx );

// return euclidean norm of x
CFUNCDECL
double
dnrm2_  ( const blas_int_t *    n,
          const double *        dx,
          const blas_int_t *    incx );

// scale x by a
CFUNCDECL
void
dscal_ ( const blas_int_t *     n,
         const double *         da,
         double *               dx,
         const blas_int_t *     incx );

// interchange x and y
CFUNCDECL
void
dswap_ ( const blas_int_t *     n,
         double *               dx,
         const blas_int_t *     incx,
         double *               dy,
         const blas_int_t *     incy );

// finds index with element having maximal absolut value
CFUNCDECL
blas_int_t
idamax_ ( const blas_int_t *    n,
          const double   *,
          const blas_int_t * );
    
// y = alpha A x + beta y
CFUNCDECL
void
dgemv_ ( const char *           trans,
         const blas_int_t *     M,
         const blas_int_t *     N,
         const double *         alpha,
         const double *         A,
         const blas_int_t *     lda,
         const double *         dx,
         const blas_int_t *     incx,
         const double *         beta,
         double *               dy,
         const blas_int_t *     incy );

// compute c = alpha * a * b + beta * c
CFUNCDECL
void
dgemm_ ( const char *           transa,
         const char *           transb,
         const blas_int_t *     n,
         const blas_int_t *     m,
         const blas_int_t *     k,
         const double *         alpha,
         const double *         a,
         const blas_int_t *     lda,
         const double *         b,
         const blas_int_t *     ldb,
         const double *         beta,
         double *               c,
         const blas_int_t *     ldc );

// computes y = op(A)*x for upper/lower triangular A (y overwrites x)
void
dtrmv_ ( const char *           uplo,
         const char *           trans,
         const char *           diag,
         const blas_int_t *     n,
         const double *         A,
         const blas_int_t *     ldA,
         double *               x,
         const blas_int_t *     incx );

// computes B = alpha * op(A) * B or B = alpha * B * op(A)
CFUNCDECL
void
dtrmm_ ( const char *           side,
         const char *           uplo,
         const char *           transa,
         const char *           diag,
         const blas_int_t *     n,
         const blas_int_t *     m,
         const double *         alpha,
         const double *         A,
         const blas_int_t *     lda,
         double *               B,
         const blas_int_t *     ldb );

// solves A x = b or A^T x = b with triangular A
CFUNCDECL
void
dtrsv_ ( const char *           uplo,
         const char *           trans,
         const char *           diag,
         const blas_int_t *     n,
         const double *         A,
         const blas_int_t *     lda,
         double *               x,
         const blas_int_t *     incx );
    
// solves A X = B or A^T X = B with triangular A
CFUNCDECL
void
dtrsm_ ( const char *           side,
         const char *           uplo,
         const char *           trans,
         const char *           diag,
         const blas_int_t *     n,
         const blas_int_t *     m,
         const double *         alpha,
         const double *         A,
         const blas_int_t *     lda,
         double *               B,
         const blas_int_t *     ldb );

// A = alpha * x * y^T + A, A \in \R^{m x n}
CFUNCDECL
void
dger_ ( const blas_int_t *      m,
        const blas_int_t *      n,
        const double *          alpha,
        const double *          x,
        const blas_int_t *      incx,
        const double *          y,
        const blas_int_t *      incy,
        double *                A,
        const blas_int_t *      lda );

///////////////////////////////////////////////////////////////////
//
// complex-valued functions (single precision)
//
///////////////////////////////////////////////////////////////////

// set (part of) A to alpha/beta
CFUNCDECL
void
claset_  ( const char *                 UPLO,
           const blas_int_t *           M,
           const blas_int_t *           N,
           const std::complex<float> *       ALPHA,
           const std::complex<float> *       BETA,
           std::complex<float> *             A,
           const blas_int_t *           LDA   );

// copy (part of) A to B
CFUNCDECL
void
clacpy_  ( const char *                 UPLO,
           const blas_int_t *           M,
           const blas_int_t *           N,
           const std::complex<float> *       A,
           const blas_int_t *           LDA,
           std::complex<float> *             B,
           const blas_int_t *           LDB  );

// compute dotproduct between x and y
CFUNCDECL
void
xcdotu_ ( const blas_int_t *            n,
          const std::complex<float> *        dx,
          const blas_int_t *            incx,
          const std::complex<float> *        dy,
          const blas_int_t *            incy,
          std::complex<float> *              retval );
CFUNCDECL
void
xcdotc_ ( const blas_int_t *            n,
          const std::complex<float> *        dx,
          const blas_int_t *            incx,
          const std::complex<float> *        dy,
          const blas_int_t *            incy,
          std::complex<float> *              retval );

// copy vector x into vector y
CFUNCDECL
void
ccopy_ ( const blas_int_t *             n,
         const std::complex<float> *         dx,
         const blas_int_t *             incx,
         std::complex<float> *               dy,
         const blas_int_t *             incy);

// compute y = y + a * x
CFUNCDECL
void
caxpy_ ( const blas_int_t *             n,
         const std::complex<float> *         da,
         const std::complex<float> *         dx,
         const blas_int_t *             incx,
         std::complex<float> *               dy,
         const blas_int_t *             incy );

// compute euclidean norm
CFUNCDECL
float
scnrm2_ ( const blas_int_t *            n,
          const std::complex<float> *        x,
          const blas_int_t *            incx );

// scale x by a
CFUNCDECL
void
cscal_  ( const blas_int_t *            n,
          const std::complex<float> *        da,
          std::complex<float> *              dx,
          const blas_int_t *            incx );
CFUNCDECL
void
csscal_ ( const blas_int_t *            n,
          const float *                 da,
          std::complex<float> *              dx,
          const blas_int_t *            incx );

// interchange x and y
CFUNCDECL
void
cswap_ ( const blas_int_t *     n,
         std::complex<float> *       dx,
         const blas_int_t *     incx,
         std::complex<float> *       dy,
         const blas_int_t *     incy );

// finds index with element having maximal absolut value
CFUNCDECL
blas_int_t
icamax_ ( const blas_int_t *    n,
          const std::complex<float> *,
          const blas_int_t * );

// y = alpha A x + beta y
CFUNCDECL
void
cgemv_ ( const char *                   trans,
         const blas_int_t *             M,
         const blas_int_t *             N,
         const std::complex<float> *         alpha,
         const std::complex<float> *         A,
         const blas_int_t *             lda,
         const std::complex<float> *         dx,
         const blas_int_t *             incx,
         const std::complex<float> *         beta,
         std::complex<float> *               dy,
         const blas_int_t *             incy );

// compute c = alpha * a * b + beta * c
CFUNCDECL
void
cgemm_ ( const char *                   transa,
         const char *                   transb,
         const blas_int_t *             n,
         const blas_int_t *             m,
         const blas_int_t *             k,
         const std::complex<float> *         alpha,
         const std::complex<float> *         a,
         const blas_int_t *             lda,
         const std::complex<float> *         b,
         const blas_int_t *             ldb,
         const std::complex<float> *         beta,
         std::complex<float> *               c,
         const blas_int_t *             ldc );

// computes y = op(A)*x for upper/lower triangular A (y overwrites x)
void
ctrmv_ ( const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const blas_int_t *             n,
         const std::complex<float> *         A,
         const blas_int_t *             ldA,
         std::complex<float> *               x,
         const blas_int_t *             incx );

// computes B = alpha * op(A) * B or B = alpha * B * op(A)
CFUNCDECL
void
ctrmm_ ( const char *                   side,
         const char *                   uplo,
         const char *                   transa,
         const char *                   diag,
         const blas_int_t *             m,
         const blas_int_t *             n,
         const std::complex<float> *         alpha,
         const std::complex<float> *         a,
         const blas_int_t *             lda,
         std::complex<float> *               b,
         const blas_int_t *             ldb );

// solves A x = b or A^T x = b with triangular A
CFUNCDECL
void
ctrsv_ ( const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const blas_int_t *             n,
         const std::complex<float> *         a,
         const blas_int_t *             lda,
         std::complex<float> *               x,
         const blas_int_t *             incx );

// solves A X = B or A^T X = B with triangular A
CFUNCDECL
void
ctrsm_ ( const char *                   side,
         const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const blas_int_t *             n,
         const blas_int_t *             m,
         const std::complex< float > *       alpha,
         const std::complex< float > *       A,
         const blas_int_t *             lda,
         std::complex< float > *             B,
         const blas_int_t *             ldb );

// rank-1 update: A = alpha * x * y^T + A
CFUNCDECL
void
cgeru_ ( const blas_int_t *             m,
         const blas_int_t *             n,
         const std::complex<float> *         alpha,
         const std::complex<float> *         x,
         const blas_int_t *             incx,
         const std::complex<float> *         y,
         const blas_int_t *             incy,
         std::complex<float> *               A,
         const blas_int_t *             lda );

// rank-1 update: A = alpha * x * conj(y^T) + A
CFUNCDECL
void
cgerc_ ( const blas_int_t *             m,
         const blas_int_t *             n,
         const std::complex<float> *         alpha,
         const std::complex<float> *         x,
         const blas_int_t *             incx,
         const std::complex<float> *         y,
         const blas_int_t *             incy,
         std::complex<float> *               A,
         const blas_int_t *             lda );

///////////////////////////////////////////////////////////////////
//
// complex-valued functions (double precision)
//
///////////////////////////////////////////////////////////////////

// set (part of) A to alpha/beta
CFUNCDECL
void
zlaset_  ( const char *                 UPLO,
           const blas_int_t *           M,
           const blas_int_t *           N,
           const std::complex<double> *      ALPHA,
           const std::complex<double> *      BETA,
           std::complex<double> *            A,
           const blas_int_t *           LDA   );

// copy (part of) A to B
CFUNCDECL
void
zlacpy_  ( const char *                 UPLO,
           const blas_int_t *           M,
           const blas_int_t *           N,
           const std::complex<double> *      A,
           const blas_int_t *           LDA,
           std::complex<double> *            B,
           const blas_int_t *           LDB  );

// compute dotproduct between x and y
CFUNCDECL
void
xzdotu_ ( const blas_int_t *            n,
          const std::complex<double> *       dx,
          const blas_int_t *            incx,
          const std::complex<double> *       dy,
          const blas_int_t *            incy,
          std::complex<double> *             retval );
CFUNCDECL
void
xzdotc_ ( const blas_int_t *            n,
          const std::complex<double> *       dx,
          const blas_int_t *            incx,
          const std::complex<double> *       dy,
          const blas_int_t *            incy,
          std::complex<double> *             retval );

// copy vector x into vector y
CFUNCDECL
void
zcopy_ ( const blas_int_t *             n,
         const std::complex<double> *        dx,
         const blas_int_t *             incx,
         std::complex<double>       *        dy,
         const blas_int_t *             incy);

// compute y = y + a * x
CFUNCDECL
void
zaxpy_ ( const blas_int_t *             n,
         const std::complex<double> *        da,
         const std::complex<double> *        dx,
         const blas_int_t *             incx,
         std::complex<double> *              dy,
         const blas_int_t *             incy );

// compute euclidean norm
CFUNCDECL
double
dznrm2_ ( const blas_int_t *            n,
          const std::complex<double> *       x,
          const blas_int_t *            incx );

// scale x by a
CFUNCDECL
void
zscal_  ( const blas_int_t *            n,
          const std::complex<double> *       da,
          std::complex<double> *             dx,
          const blas_int_t *            incx );
CFUNCDECL
void
zdscal_ ( const blas_int_t *            n,
          const double *                da,
          std::complex<double> *             dx,
          const blas_int_t *            incx );

// interchange x and y
CFUNCDECL
void
zswap_ ( const blas_int_t *             n,
         std::complex<double> *              dx,
         const blas_int_t *             incx,
         std::complex<double> *              dy,
         const blas_int_t *             incy );

// finds index with element having maximal absolut value
CFUNCDECL
blas_int_t
izamax_ ( const blas_int_t *            n,
          const std::complex<double> *,
          const blas_int_t * );
    
// y = alpha A x + beta y
CFUNCDECL
void
zgemv_ ( const char *                   trans,
         const blas_int_t *             M,
         const blas_int_t *             N,
         const std::complex<double> *        alpha,
         const std::complex<double> *        A,
         const blas_int_t *             lda,
         const std::complex<double> *        dx,
         const blas_int_t *             incx,
         const std::complex<double> *        beta,
         std::complex<double> *              dy,
         const blas_int_t *             incy );

// compute c = alpha * a * b + beta * c
CFUNCDECL
void
zgemm_ ( const char *                   transa,
         const char *                   transb,
         const blas_int_t *             n,
         const blas_int_t *             m,
         const blas_int_t *             k,
         const std::complex<double> *        alpha,
         const std::complex<double> *        a,
         const blas_int_t *             lda,
         const std::complex<double> *        b,
         const blas_int_t *             ldb,
         const std::complex<double> *        beta,
         std::complex<double> *              c,
         const blas_int_t *             ldc );

// computes y = op(A)*x for upper/lower triangular A (y overwrites x)
void
ztrmv_ ( const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const blas_int_t *             n,
         const std::complex<double> *        A,
         const blas_int_t *             ldA,
         std::complex<double> *              x,
         const blas_int_t *             incx );

// computes B = alpha * op(A) * B or B = alpha * B * op(A)
CFUNCDECL
void
ztrmm_ ( const char *                   side,
         const char *                   uplo,
         const char *                   transa,
         const char *                   diag,
         const blas_int_t *             m,
         const blas_int_t *             n,
         const std::complex<double> *        alpha,
         const std::complex<double> *        a,
         const blas_int_t *             lda,
         std::complex<double> *              b,
         const blas_int_t *             ldb );

// solves A x = b or A^T x = b with triangular A
CFUNCDECL
void
ztrsv_ ( const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const blas_int_t *             n,
         const std::complex<double> *        a,
         const blas_int_t *             lda,
         std::complex<double> *              x,
         const blas_int_t *             incx );

// solves A X = B or A^T X = B with triangular A
CFUNCDECL
void
ztrsm_ ( const char *                   side,
         const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const blas_int_t *             n,
         const blas_int_t *             m,
         const std::complex< double > *      alpha,
         const std::complex< double > *      A,
         const blas_int_t *             lda,
         std::complex< double > *            B,
         const blas_int_t *             ldb );

// rank-1 update: A = alpha * x * y^T + A
CFUNCDECL
void
zgeru_ ( const blas_int_t *             m,
         const blas_int_t *             n,
         const std::complex<double> *        alpha,
         const std::complex<double> *        x,
         const blas_int_t *             incx,
         const std::complex<double> *        y,
         const blas_int_t *             incy,
         std::complex<double> *              A,
         const blas_int_t *             lda );

// rank-1 update: A = alpha * x * conj(y^T) + A
CFUNCDECL
void
zgerc_ ( const blas_int_t *             m,
         const blas_int_t *             n,
         const std::complex<double> *        alpha,
         const std::complex<double> *        x,
         const blas_int_t *             incx,
         const std::complex<double> *        y,
         const blas_int_t *             incy,
         std::complex<double> *              A,
         const blas_int_t *             lda );
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// definition of external Lapack functions
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

extern "C" {
    
//////////////////////////////////////////////////////////////
//
// real-valued functions (single precision)
//
//////////////////////////////////////////////////////////////

// solve linear system of equations
CFUNCDECL
void
sgesv_   ( const blas_int_t *   n,
           const blas_int_t *   nrhs,
           float *              A,
           const blas_int_t *   lda,
           blas_int_t *         ipiv,
           float *              B,
           const blas_int_t *   ldb,
           blas_int_t *         info );

// compute inverse of triangular system
CFUNCDECL
void
strtri_  ( const char *         uplo,
           const char *         diag,
           const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           blas_int_t *         info );

// compute eigenvalues and eigenvectors of a tridiagonal, symmetric matrix
CFUNCDECL
void
sstev_   ( const char *         jobz,
           const blas_int_t *   n,
           float *              D,
           float *              E,
           float *              Z,
           const blas_int_t *   ldz,
           float *              work,
           blas_int_t *         info );

// compute eigenvalues and eigenvectors of a symmetric matrix
CFUNCDECL
void
ssyev_   ( const char *         jobz,
           const char *         uplo,
           const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           float *              w,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );
    
CFUNCDECL
void
ssyevd_  ( const char *         jobz,
           const char *         uplo,
           const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           float *              w,
           float *              work,
           const blas_int_t *   lwork,
           const blas_int_t *   iwork,
           const blas_int_t *   liwork,
           blas_int_t *         info );
    
// compute selected eigenvalues and eigenvectors of a symmetric matrix
CFUNCDECL
void ssyevx_ ( const char * jobz, const char * range, const char * uplo,
               const blas_int_t * n, float * A, const blas_int_t * ldA,
               const float * vl, const float * vu, const blas_int_t * il, const blas_int_t * iu,
               const float * abstol, blas_int_t * m, float * W, float * Z, const blas_int_t * ldZ,
               float * work, const blas_int_t * lwork, blas_int_t * iwork, blas_int_t * ifail,
               blas_int_t * info );
    
// compute singular-value-decomposition
CFUNCDECL
void
sgesvd_  ( const char *         jobu,
           const char *         jobv,
           const blas_int_t *   n,
           const blas_int_t *   m,
           float *              A,
           const blas_int_t *   lda,
           float *              S,
           float *              U,
           const blas_int_t *   ldu,
           float *              V,
           const blas_int_t *   ldv,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
sgesdd_  ( const char *         job,
           const blas_int_t *   n,
           const blas_int_t *   m,
           float *              A,
           const blas_int_t *   lda,
           float *              S,
           float *              U,
           const blas_int_t *   ldu,
           float *              VT,
           const blas_int_t *   ldvt,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         iwork,
           blas_int_t *         info );

CFUNCDECL
void
sgesvj_  ( const char *         joba,
           const char *         jobu,
           const char *         jobv,
           const blas_int_t *   m,
           const blas_int_t *   n,
           float *              a,
           const blas_int_t *   lda,
           float *              sva,
           const blas_int_t *   mv,
           float *              v,
           const blas_int_t *   ldv,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
sgejsv_  ( const char *         joba,
           const char *         jobu,
           const char *         jobv,
           const char *         jobr,
           const char *         jobt,
           const char *         jobp,
           const blas_int_t *   m,
           const blas_int_t *   n,
           float *              a,
           const blas_int_t *   lda,
           float *              sva,
           float *              u,
           const blas_int_t *   ldu,
           float *              v,
           const blas_int_t *   ldv,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         iwork,
           blas_int_t *         info );

// compute QR-factorisation
CFUNCDECL
void
sgeqrf_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           float *              tau,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
sgeqr2_ ( const blas_int_t *  nrows,
          const blas_int_t *  ncols,
          float *             A,
          const blas_int_t *  lda,
          float *             tau,
          float *             work,
          const blas_int_t *  info );

CFUNCDECL
void
sorgqr_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           const blas_int_t *   k,
           float *              A,
           const blas_int_t *   lda,
           const float *        tau,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
sorg2r_ ( const blas_int_t *        nrows,
          const blas_int_t *        ncols,
          const blas_int_t *        nref,
          float *                   A,
          const blas_int_t *        lda,
          const float *             tau,
          float *                   work,
          const blas_int_t *        info );

CFUNCDECL
void
sgeqp3_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           blas_int_t *         jpvt,
           float *              tau,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

#if HPRO_HAS_GEQP3_TRUNC == 1

CFUNCDECL
void
sgeqp3trunc_  ( const blas_int_t *   m,
                const blas_int_t *   n,
                float *              A,
                const blas_int_t *   lda,
                blas_int_t *         jpvt,
                float *              tau,
                blas_int_t *         ntrunc,
                float *              atrunc,
                float *              rtrunc,
                float *              work,
                const blas_int_t *   lwork,
                blas_int_t *         info );

#endif

// compute LQ factorisation of given matrix A
CFUNCDECL
void
sgelqf_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           float *              tau,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
sgelq2_ ( const blas_int_t *  nrows,
          const blas_int_t *  ncols,
          float *             A,
          const blas_int_t *  lda,
          float *             tau,
          float *             work,
          const blas_int_t *  info );

CFUNCDECL
void
sorglq_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           const blas_int_t *   k,
           float *              A,
           const blas_int_t *   lda,
           const float *        tau,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
sorg2l_ ( const blas_int_t *        nrows,
          const blas_int_t *        ncols,
          const blas_int_t *        nref,
          float *                   A,
          const blas_int_t *        lda,
          const float *             tau,
          float *                   work,
          const blas_int_t *        info );


// compute LU factorisation of given matrix A
CFUNCDECL
void
sgetrf_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           blas_int_t *         ipiv,
           blas_int_t *         info );

// compute inverse of A (using result from getrf)
CFUNCDECL
void
sgetri_  ( const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           blas_int_t *         ipiv,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

// determine machine parameters
CFUNCDECL
float
slamch_  ( char *               cmach );

// generate plane-rotation
CFUNCDECL
blas_int_t
slartg_  ( float *              f,
           float *              g,
           float *              cs,
           float *              sn,
           float *              r );

// compute householder reflection
CFUNCDECL
void
slarfg_ ( const blas_int_t *      n,
          const float *           alpha,
          const float *           x,
          const blas_int_t *      incx,
          const float *           tau );

// apply householder reflection
CFUNCDECL
void
slarf_  ( const char *            side,
          const blas_int_t *      n,
          const blas_int_t *      m,
          const float *           V,
          const blas_int_t *      incv,
          const float *           tau,
          float *                 C,
          const blas_int_t *      ldc,
          const float *           work );

//////////////////////////////////////////////////////////////
//
// real-valued functions (double precision)
//
//////////////////////////////////////////////////////////////

// solve linear system of equations
CFUNCDECL
void
dgesv_   ( const blas_int_t *   n,
           const blas_int_t *   nrhs,
           double *             A,
           const blas_int_t *   lda,
           blas_int_t *         ipiv,
           double *             B,
           const blas_int_t *   ldb,
           blas_int_t *         info );

// compute inverse of triangular system
CFUNCDECL
void
dtrtri_  ( const char *         uplo,
           const char *         diag,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           blas_int_t *         info );
    
// compute eigenvalues and eigenvectors of a tridiagonal, symmetric matrix
CFUNCDECL
void
dstev_   ( const char *         jobz,
           const blas_int_t *   n,
           double *             D,
           double *             E,
           double *             Z,
           const blas_int_t *   ldz,
           double *             work,
           blas_int_t *         info );
    
// compute eigenvalues and eigenvectors of a symmetric matrix
CFUNCDECL
void
dsyev_   ( const char *         jobz,
           const char *         uplo,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           double *             w,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );
    
CFUNCDECL
void
dsyevd_  ( const char *         jobz,
           const char *         uplo,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           double *             w,
           double *             work,
           const blas_int_t *   lwork,
           const blas_int_t *   iwork,
           const blas_int_t *   liwork,
           blas_int_t *         info );
    
// compute selected eigenvalues and eigenvectors of a symmetric matrix
CFUNCDECL
void
dsyevx_  ( const char *         jobz,
           const char *         range,
           const char *         uplo,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   ldA,
           const double *       vl,
           const double *       vu,
           const blas_int_t *   il,
           const blas_int_t *   iu,
           const double *       abstol,
           blas_int_t *         m,
           double *             W,
           double *             Z,
           const blas_int_t *   ldZ,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         iwork,
           blas_int_t *         ifail,
           blas_int_t *         info );
    
// compute singular-value-decomposition
CFUNCDECL
void
dgesvd_  ( const char *         jobu,
           const char *         jobv,
           const blas_int_t *   n,
           const blas_int_t *   m,
           double *             A,
           const blas_int_t *   lda,
           double *             S,
           double *             U,
           const blas_int_t *   ldu,
           double *             V,
           const blas_int_t *   ldv,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
dgesdd_  ( const char *         jobz,
           const blas_int_t *   n,
           const blas_int_t *   m,
           double *             A,
           const blas_int_t *   lda,
           double *             S,
           double *             U,
           const blas_int_t *   ldu,
           double *             VT,
           const blas_int_t *   ldvt,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         iwork,
           blas_int_t *         info );

CFUNCDECL
void
dgejsv_  ( const char *         joba,
           const char *         jobu,
           const char *         jobv,
           const char *         jobr,
           const char *         jobt,
           const char *         jobp,
           const blas_int_t *   m,
           const blas_int_t *   n,
           double *             a,
           const blas_int_t *   lda,
           double *             sva,
           double *             u,
           const blas_int_t *   ldu,
           double *             v,
           const blas_int_t *   ldv,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         iwork,
           blas_int_t *         info );

CFUNCDECL
void
dgesvj_  ( const char *         joba,
           const char *         jobu,
           const char *         jobv,
           const blas_int_t *   m,
           const blas_int_t *   n,
           double *             a,
           const blas_int_t *   lda,
           double *             sva,
           const blas_int_t *   mv,
           double *             v,
           const blas_int_t *   ldv,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

// compute QR-factorisation
CFUNCDECL
void
dgeqrf_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           double *             tau,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
dgeqr2_ ( const blas_int_t *  nrows,
          const blas_int_t *  ncols,
          double *            A,
          const blas_int_t *  lda,
          double *            tau,
          double *            work,
          const blas_int_t *  info );

CFUNCDECL
void
dorgqr_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           const blas_int_t *   k,
           double *             A,
           const blas_int_t *   lda,
           const double *       tau,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
dorg2r_ ( const blas_int_t *        nrows,
          const blas_int_t *        ncols,
          const blas_int_t *        nref,
          double *                  A,
          const blas_int_t *        lda,
          const double *            tau,
          double *                  work,
          const blas_int_t *        info );

CFUNCDECL
void
dgeqp3_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           blas_int_t *         jpvt,
           double *             tau,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

#if HPRO_HAS_GEQP3_TRUNC == 1

CFUNCDECL
void
dgeqp3trunc_  ( const blas_int_t *   m,
                const blas_int_t *   n,
                double *             A,
                const blas_int_t *   lda,
                blas_int_t *         jpvt,
                double *             tau,
                blas_int_t *         ntrunc,
                double *             atrunc,
                double *             rtrunc,
                double *             work,
                const blas_int_t *   lwork,
                blas_int_t *         info );

#endif

// compute LQ-factorisation
CFUNCDECL
void
dgelqf_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           double *             tau,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
dgelq2_ ( const blas_int_t *  nrows,
          const blas_int_t *  ncols,
          double *            A,
          const blas_int_t *  lda,
          double *            tau,
          double *            work,
          const blas_int_t *  info );

CFUNCDECL
void
dorglq_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           const blas_int_t *   k,
           double *             A,
           const blas_int_t *   lda,
           const double *       tau,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
dorg2l_ ( const blas_int_t *        nrows,
          const blas_int_t *        ncols,
          const blas_int_t *        nref,
          double *                  A,
          const blas_int_t *        lda,
          const double *            tau,
          double *                  work,
          const blas_int_t *        info );


// compute LU factorisation of given matrix A
CFUNCDECL
void
dgetrf_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           blas_int_t *         ipiv,
           blas_int_t *         info );
    
// compute inverse of A (using result from getrf)
CFUNCDECL
void
dgetri_  ( const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           blas_int_t *         ipiv,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );
    
// determine machine parameters
CFUNCDECL
double
dlamch_  ( char *               cmach );

// generate plane-rotation
CFUNCDECL
blas_int_t
dlartg_  ( double *             f,
           double *             g,
           double *             cs,
           double *             sn,
           double *             r );

// compute householder reflection
CFUNCDECL
void
dlarfg_ ( const blas_int_t *      n,
          const double *          alpha,
          const double *          x,
          const blas_int_t *      incx,
          const double *          tau );

// apply householder reflection
CFUNCDECL
void
dlarf_  ( const char *            side,
          const blas_int_t *      n,
          const blas_int_t *      m,
          const double *          V,
          const blas_int_t *      incv,
          const double *          tau,
          double *                C,
          const blas_int_t *      ldc,
          const double *          work );

//////////////////////////////////////////////////////////////
//
// complex-valued functions (single precision)
//
//////////////////////////////////////////////////////////////

// solve linear system of equations
CFUNCDECL
void
cgesv_   ( const blas_int_t *           n,
           const blas_int_t *           nrhs,
           std::complex<float> *        a,
           const blas_int_t *           lda,
           blas_int_t *                 ipiv,
           std::complex<float> *        b,
           const blas_int_t *           ldb,
           blas_int_t *                 info );

// compute inverse of triangular system
CFUNCDECL
void
ctrtri_  ( const char *                 uplo,
           const char *                 diag,
           const blas_int_t *           n,
           std::complex< float > *      A,
           const blas_int_t *           lda,
           blas_int_t *                 info );
    
// compute eigenvalues and eigenvectors of a symmetric matrix
CFUNCDECL
void
cheev_   ( const char *                 jobz,
           const char *                 uplo,
           const blas_int_t *           n,
           std::complex<float> *        A,
           const blas_int_t *           lda,
           float *                      w,
           std::complex<float> *        work,
           const blas_int_t *           lwork,
           float *                      rwork,
           blas_int_t *                 info );
    
CFUNCDECL
void
cheevd_  ( const char *                 jobz,
           const char *                 uplo,
           const blas_int_t *           n,
           std::complex<float> *        A,
           const blas_int_t *           lda,
           float *                      w,
           std::complex<float> *        work,
           const blas_int_t *           lwork,
           float *                      rwork,
           const blas_int_t *           lrwork,
           const blas_int_t *           iwork,
           const blas_int_t *           liwork,
           blas_int_t *                 info );
    
// compute singular-value-decomposition
CFUNCDECL
void
cgesvd_  ( const char *                 jobu,
           const char *                 jobv,
           const blas_int_t *           n,
           const blas_int_t *           m,
           std::complex<float> *        A,
           const blas_int_t *           lda,
           float *                      S,
           std::complex<float> *        U,
           const blas_int_t *           ldu,
           std::complex<float> *        V,
           const blas_int_t *           ldv,
           std::complex<float> *        work,
           const blas_int_t *           lwork,
           float *                      rwork,
           blas_int_t *                 info );

CFUNCDECL
void
cgesdd_  ( const char *                 job,
           const blas_int_t *           n,
           const blas_int_t *           m,
           std::complex<float> *        A,
           const blas_int_t *           lda,
           float *                      S,
           std::complex<float> *        U,
           const blas_int_t *           ldu,
           std::complex<float> *        VT,
           const blas_int_t *           ldvt,
           std::complex<float> *        work,
           const blas_int_t *           lwork,
           float *                      rwork,
           const blas_int_t *           iwork,
           blas_int_t *                 info );

CFUNCDECL
void
cgesvj_  ( const char *                 joba,
           const char *                 jobu,
           const char *                 jobv,
           const blas_int_t *           n,
           const blas_int_t *           m,
           std::complex<float> *        A,
           const blas_int_t *           lda,
           float *                      S,
           const blas_int_t *           mv,
           std::complex<float> *        V,
           const blas_int_t *           ldv,
           std::complex<float> *        cwork,
           const blas_int_t *           lwork,
           float *                      rwork,
           const blas_int_t *           lrwork,
           blas_int_t *                 info );

// compute QR-factorisation
CFUNCDECL
void
cgeqrf_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           std::complex<float> *        A,
           const blas_int_t *           lda,
           std::complex<float> *        tau,
           std::complex<float> *        work,
           const blas_int_t *           lwork,
           blas_int_t *                 info );

CFUNCDECL
void
cgeqr2_ ( const blas_int_t *        nrows,
          const blas_int_t *        ncols,
          std::complex< float > *   A,
          const blas_int_t *        lda,
          std::complex< float > *   tau,
          std::complex< float > *   work,
          const blas_int_t *        info );

CFUNCDECL
void
cungqr_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           const blas_int_t *           k,
           std::complex<float> *        a,
           const blas_int_t *           lda,
           const std::complex<float> *  tau,
           std::complex<float> *        work,
           const blas_int_t *           lwork,
           blas_int_t *                 info );

CFUNCDECL
void
cung2r_ ( const blas_int_t *            nrows,
          const blas_int_t *            ncols,
          const blas_int_t *            nref,
          std::complex< float > *       A,
          const blas_int_t *            lda,
          const std::complex< float > * tau,
          std::complex< float > *       work,
          const blas_int_t *            info );

CFUNCDECL
void
cgeqp3_  ( const blas_int_t *     m,
           const blas_int_t *     n,
           std::complex<float> *  A,
           const blas_int_t *     lda,
           blas_int_t *           jpvt,
           std::complex<float> *  tau,
           std::complex<float> *  work,
           const blas_int_t *     lwork,
           float *                rwork,
           blas_int_t *           info );

#if HPRO_HAS_GEQP3_TRUNC == 1

CFUNCDECL
void
cgeqp3trunc_  ( const blas_int_t *     m,
                const blas_int_t *     n,
                std::complex<float> *  A,
                const blas_int_t *     lda,
                blas_int_t *           jpvt,
                std::complex<float> *  tau,
                blas_int_t *           ntrunc,
                float *                atrunc,
                float *                rtrunc,
                std::complex<float> *  work,
                const blas_int_t *     lwork,
                float *                rwork,
                blas_int_t *           info );

#endif

// compute LQ-factorisation
CFUNCDECL
void
cgelqf_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           std::complex<float> *        A,
           const blas_int_t *           lda,
           std::complex<float> *        tau,
           std::complex<float> *        work,
           const blas_int_t *           lwork,
           blas_int_t *                 info );

CFUNCDECL
void
cgelq2_ ( const blas_int_t *        nrows,
          const blas_int_t *        ncols,
          std::complex< float > *   A,
          const blas_int_t *        lda,
          std::complex< float > *   tau,
          std::complex< float > *   work,
          const blas_int_t *        info );

CFUNCDECL
void
cunglq_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           const blas_int_t *           k,
           std::complex<float> *        a,
           const blas_int_t *           lda,
           const std::complex<float> *  tau,
           std::complex<float> *        work,
           const blas_int_t *           lwork,
           blas_int_t *                 info );

CFUNCDECL
void
cung2l_ ( const blas_int_t *            nrows,
          const blas_int_t *            ncols,
          const blas_int_t *            nref,
          std::complex< float > *       A,
          const blas_int_t *            lda,
          const std::complex< float > * tau,
          std::complex< float > *       work,
          const blas_int_t *            info );

// compute LU factorisation of given matrix A
CFUNCDECL
void
cgetrf_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           std::complex<float> *        a,
           const blas_int_t *           lda,
           blas_int_t *                 ipiv,
           blas_int_t *                 info);

// compute inverse of A (using result from getrf)
CFUNCDECL
void
cgetri_  ( const blas_int_t *           n,
           std::complex<float> *        a,
           const blas_int_t *           lda,
           blas_int_t *                 ipiv,
           std::complex<float> *        work,
           const blas_int_t *           lwork,
           blas_int_t *                 info);

// compute householder reflection
CFUNCDECL
void
clarfg_ ( const blas_int_t *            n,
          const std::complex< float > * alpha,
          const std::complex< float > * x,
          const blas_int_t *            incx,
          const std::complex< float > * tau );

// apply householder reflection
CFUNCDECL
void
clarf_  ( const char *                  side,
          const blas_int_t *            n,
          const blas_int_t *            m,
          const std::complex< float > * V,
          const blas_int_t *            incv,
          const std::complex< float > * tau,
          std::complex< float > *       C,
          const blas_int_t *            ldc,
          const std::complex< float > * work );

//////////////////////////////////////////////////////////////
//
// complex-valued functions (double precision)
//
//////////////////////////////////////////////////////////////

// solve linear system of equations
CFUNCDECL
void
zgesv_   ( const blas_int_t *           n,
           const blas_int_t *           nrhs,
           std::complex<double> *       a,
           const blas_int_t *           lda,
           blas_int_t *                 ipiv,
           std::complex<double> *       b,
           const blas_int_t *           ldb,
           blas_int_t *                 info );

// compute inverse of triangular system
CFUNCDECL
void
ztrtri_  ( const char *                 uplo,
           const char *                 diag,
           const blas_int_t *           n,
           std::complex< double > *     A,
           const blas_int_t *           lda,
           blas_int_t *                 info );
    
// compute eigenvalues and eigenvectors of a symmetric matrix
CFUNCDECL
void
zheev_   ( const char *                 jobz,
           const char *                 uplo,
           const blas_int_t *           n,
           std::complex<double> *       A,
           const blas_int_t *           lda,
           double *                     w,
           std::complex<double> *       work,
           const blas_int_t *           lwork,
           double *                     rwork,
           blas_int_t *                 info );
    
CFUNCDECL
void
zheevd_  ( const char *                 jobz,
           const char *                 uplo,
           const blas_int_t *           n,
           std::complex<double> *       A,
           const blas_int_t *           lda,
           double *                     w,
           std::complex<double> *       work,
           const blas_int_t *           lwork,
           double *                     rwork,
           const blas_int_t *           lrwork,
           const blas_int_t *           iwork,
           const blas_int_t *           liwork,
           blas_int_t *                 info );
    
// compute singular-value-decomposition
CFUNCDECL
void
zgesvd_  ( const char *                 jobu,
           const char *                 jobv,
           const blas_int_t *           n,
           const blas_int_t *           m,
           std::complex<double> *       A,
           const blas_int_t *           lda,
           double *                     S,
           std::complex<double> *       U,
           const blas_int_t *           ldu,
           std::complex<double> *       V,
           const blas_int_t *           ldv,
           std::complex<double> *       work,
           const blas_int_t *           lwork,
           double *                     rwork,
           blas_int_t *                 info );

CFUNCDECL
void
zgesdd_  ( const char *                 job,
           const blas_int_t *           n,
           const blas_int_t *           m,
           std::complex<double> *       A,
           const blas_int_t *           lda,
           double *                     S,
           std::complex<double> *       U,
           const blas_int_t *           ldu,
           std::complex<double> *       VT,
           const blas_int_t *           ldvt,
           std::complex<double> *       work,
           const blas_int_t *           lwork,
           double *                     rwork,
           const blas_int_t *           iwork,
           blas_int_t *                 info );

CFUNCDECL
void
zgesvj_  ( const char *                 joba,
           const char *                 jobu,
           const char *                 jobv,
           const blas_int_t *           n,
           const blas_int_t *           m,
           std::complex<double> *       A,
           const blas_int_t *           lda,
           double *                     S,
           const blas_int_t *           mv,
           std::complex<double> *       V,
           const blas_int_t *           ldv,
           std::complex<double> *       cwork,
           const blas_int_t *           lwork,
           double *                     rwork,
           const blas_int_t *           lrwork,
           blas_int_t *                 info );

// compute QR-factorisation
CFUNCDECL
void
zgeqrf_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           std::complex<double> *       A,
           const blas_int_t *           lda,
           std::complex<double> *       tau,
           std::complex<double> *       work,
           const blas_int_t *           lwork,
           blas_int_t *                 info );

CFUNCDECL
void
zgeqr2_ ( const blas_int_t *        nrows,
          const blas_int_t *        ncols,
          std::complex< double > *  A,
          const blas_int_t *        lda,
          std::complex< double > *  tau,
          std::complex< double > *  work,
          const blas_int_t *        info );

CFUNCDECL
void
zungqr_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           const blas_int_t *           k,
           std::complex<double> *       a,
           const blas_int_t *           lda,
           const std::complex<double> * tau,
           std::complex<double> *       work,
           const blas_int_t *           lwork,
           blas_int_t *                 info );

CFUNCDECL
void
zung2r_ ( const blas_int_t *              nrows,
          const blas_int_t *              ncols,
          const blas_int_t *              nref,
          std::complex< double > *        A,
          const blas_int_t *              lda,
          const std::complex< double > *  tau,
          std::complex< double > *        work,
          const blas_int_t *              info );

CFUNCDECL
void
zgeqp3_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           std::complex<double> *       A,
           const blas_int_t *           lda,
           blas_int_t *                 jpvt,
           std::complex<double> *       tau,
           std::complex<double> *       work,
           const blas_int_t *           lwork,
           double *                     rwork,
           blas_int_t *                 info );

#if HPRO_HAS_GEQP3_TRUNC == 1

CFUNCDECL
void
zgeqp3trunc_  ( const blas_int_t *      m,
                const blas_int_t *      n,
                std::complex<double> *  A,
                const blas_int_t *      lda,
                blas_int_t *            jpvt,
                std::complex<double> *  tau,
                blas_int_t *            ntrunc,
                double *                atrunc,
                double *                rtrunc,
                std::complex<double> *  work,
                const blas_int_t *      lwork,
                double *                rwork,
                blas_int_t *            info );

#endif

// compute LQ-factorisation
CFUNCDECL
void
zgelqf_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           std::complex<double> *       A,
           const blas_int_t *           lda,
           std::complex<double> *       tau,
           std::complex<double> *       work,
           const blas_int_t *           lwork,
           blas_int_t *                 info );

CFUNCDECL
void
zgelq2_ ( const blas_int_t *        nrows,
          const blas_int_t *        ncols,
          std::complex< double > *  A,
          const blas_int_t *        lda,
          std::complex< double > *  tau,
          std::complex< double > *  work,
          const blas_int_t *        info );

CFUNCDECL
void
zunglq_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           const blas_int_t *           k,
           std::complex<double> *       a,
           const blas_int_t *           lda,
           const std::complex<double> * tau,
           std::complex<double> *       work,
           const blas_int_t *           lwork,
           blas_int_t *                 info );

CFUNCDECL
void
zung2l_ ( const blas_int_t *              nrows,
          const blas_int_t *              ncols,
          const blas_int_t *              nref,
          std::complex< double > *        A,
          const blas_int_t *              lda,
          const std::complex< double > *  tau,
          std::complex< double > *        work,
          const blas_int_t *              info );

// compute LU factorisation of given matrix A
CFUNCDECL
void
zgetrf_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           std::complex<double> *       a,
           const blas_int_t *           lda,
           blas_int_t *                 ipiv,
           blas_int_t *                 info);

// compute inverse of A (using result from getrf)
CFUNCDECL
void
zgetri_  ( const blas_int_t *           n,
           std::complex<double> *       a,
           const blas_int_t *           lda,
           blas_int_t *                 ipiv,
           std::complex<double> *       work,
           const blas_int_t *           lwork,
           blas_int_t *                 info);

// compute householder reflection
CFUNCDECL
void
zlarfg_ ( const blas_int_t *              n,
          const std::complex< double > *  alpha,
          const std::complex< double > *  x,
          const blas_int_t *              incx,
          const std::complex< double > *  tau );

// apply householder reflection
CFUNCDECL
void
zlarf_  ( const char *                    side,
          const blas_int_t *              n,
          const blas_int_t *              m,
          const std::complex< double > *  V,
          const blas_int_t *              incv,
          const std::complex< double > *  tau,
          std::complex< double > *        C,
          const blas_int_t *              ldc,
          const std::complex< double > *  work );

//////////////////////////////////////////////////////////////
//
// misc. helpers
//
//////////////////////////////////////////////////////////////

// return problem dependent parameters
blas_int_t
ilaenv_ ( const blas_int_t *    ispec,
          const char *          name,
          const char *          opts,
          const blas_int_t *    n1,
          const blas_int_t *    n2,
          const blas_int_t *    n3,
          const blas_int_t *    n4 );

}// extern "C"



namespace BLAS
{

//! @cond

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// wrappers for BLAS functions
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// *laset
//
#define HPRO_LASET_FUNC( type, func )                                  \
    inline void laset ( const char       uplo,                         \
                        const blas_int_t m,                            \
                        const blas_int_t n,                            \
                        const type       alpha,                        \
                        const type       beta,                         \
                        type *           A,                            \
                        const blas_int_t lda ) {                       \
        func( & uplo, & m, & n, & alpha, & beta, A, & lda ); }

HPRO_LASET_FUNC( float,                  slaset_ )
HPRO_LASET_FUNC( double,                 dlaset_ )
HPRO_LASET_FUNC( std::complex< float >,  claset_ )
HPRO_LASET_FUNC( std::complex< double >, zlaset_ )

#undef HPRO_LASET_FUNC

//
// *lacpy
//
#define HPRO_LACPY_FUNC( type, func )                       \
    inline void lacpy ( const char       uplo,              \
                        const blas_int_t m,                 \
                        const blas_int_t n,                 \
                        const type *     A,                 \
                        const blas_int_t lda,               \
                        type *           B,                 \
                        const blas_int_t ldb ) {            \
        func( & uplo, & m, & n, A, & lda, B, & ldb ); }

HPRO_LACPY_FUNC( float,             slacpy_ )
HPRO_LACPY_FUNC( double,            dlacpy_ )
HPRO_LACPY_FUNC( std::complex< float >,  clacpy_ )
HPRO_LACPY_FUNC( std::complex< double >, zlacpy_ )

#undef HPRO_LACPY_FUNC

//
// *scal
//
#define HPRO_SCAL_FUNC( type, func )                                    \
    inline void scal ( const blas_int_t n,                              \
                       const type       alpha,                          \
                       type *           x,                              \
                       const blas_int_t incx )                          \
        { func( & n, & alpha, x, & incx ); }

HPRO_SCAL_FUNC( float,                  sscal_ )
HPRO_SCAL_FUNC( double,                 dscal_ )
HPRO_SCAL_FUNC( std::complex< float >,  cscal_ )
HPRO_SCAL_FUNC( std::complex< double >, zscal_ )

//
// *copy
//
#define HPRO_COPY_FUNC( type, func )                                    \
    inline void copy ( const blas_int_t n, const type * x,              \
                       const blas_int_t incx,                           \
                       type * y, const blas_int_t incy )                \
    { func( & n, x, & incx, y, & incy ); }
                                                       
HPRO_COPY_FUNC( float,                  scopy_ )
HPRO_COPY_FUNC( double,                 dcopy_ )
HPRO_COPY_FUNC( std::complex< float >,  ccopy_ )
HPRO_COPY_FUNC( std::complex< double >, zcopy_ )

//
// *swap
//
#define HPRO_SWAP_FUNC( type, func )                                    \
    inline  void swap ( const blas_int_t n, type * x,                   \
                        const blas_int_t incx,                          \
                        type * y, const blas_int_t incy )               \
        { func( & n, x, & incx, y, & incy ); }
                                                       
HPRO_SWAP_FUNC( float,                  sswap_ )
HPRO_SWAP_FUNC( double,                 dswap_ )
HPRO_SWAP_FUNC( std::complex< float >,  cswap_ )
HPRO_SWAP_FUNC( std::complex< double >, zswap_ )

//
// i*amax
//
#define HPRO_MAX_IDX_FUNC( type, func )                                 \
    inline blas_int_t max_idx ( const blas_int_t n,                     \
                                type * x, const blas_int_t incx )       \
    { return func( & n, x, & incx ); }
                                                       
HPRO_MAX_IDX_FUNC( float,                  isamax_ )
HPRO_MAX_IDX_FUNC( double,                 idamax_ )
HPRO_MAX_IDX_FUNC( std::complex< float >,  icamax_ )
HPRO_MAX_IDX_FUNC( std::complex< double >, izamax_ )

//
// *axpy
//
#define HPRO_AXPY_FUNC( type, func, flops )                             \
    inline void axpy ( const blas_int_t n,                              \
                       const type alpha, const type * x,                \
                       const blas_int_t incx,                           \
                       type * y, const blas_int_t incy )                \
    {                                                                   \
        ADD_FLOPS( FLOPS_##flops( n ) );                                \
        func( & n, & alpha, x, & incx, y, & incy );                     \
    }

HPRO_AXPY_FUNC( float,                  saxpy_, SAXPY )
HPRO_AXPY_FUNC( double,                 daxpy_, DAXPY )
HPRO_AXPY_FUNC( std::complex< float >,  caxpy_, CAXPY )
HPRO_AXPY_FUNC( std::complex< double >, zaxpy_, ZAXPY )

//
// *dot(c)
//
inline
float
dot ( const blas_int_t  n,
      const float *  x, const blas_int_t  incx,
      const float *  y, const blas_int_t  incy )
{
    ADD_FLOPS( FLOPS_SDOT( n ) );
    return sdot_( & n, x, & incx, y, & incy );
}

inline
double
dot ( const blas_int_t  n,
      const double *  x, const blas_int_t  incx,
      const double *  y, const blas_int_t  incy )
{
    ADD_FLOPS( FLOPS_DDOT( n ) );
    return ddot_( & n, x, & incx, y, & incy );
}

inline
std::complex< float >
dot ( const blas_int_t  n,
      const std::complex< float > *  x, const blas_int_t  incx,
      const std::complex< float > *  y, const blas_int_t  incy )
{
    ADD_FLOPS( FLOPS_CDOT( n ) );

    std::complex< float >  res;
    
    #if HPRO_USE_MKL == 1

    cblas_cdotc_sub( n, x, incx, y, incy, & res );
    
    #else
    
    xcdotc_( & n, x, & incx, y, & incy, & res );
    
    #endif
    
    return res;
}

inline
std::complex< double >
dot ( const blas_int_t  n,
      const std::complex< double > *  x, const blas_int_t  incx,
      const std::complex< double > *  y, const blas_int_t  incy )
{
    ADD_FLOPS( FLOPS_ZDOT( n ) );

    std::complex< double >  res;
    
    #if HPRO_USE_MKL == 1

    cblas_zdotc_sub( n, x, incx, y, incy, & res );
    
    #else
    
    xzdotc_( & n, x, & incx, y, & incy, & res );
    
    #endif
    
    return res;
}

//
// *dotu
//
inline
float
dotu ( const blas_int_t n,
       const float * x, const blas_int_t incx,
       const float * y, const blas_int_t incy )
{
    ADD_FLOPS( FLOPS_DDOT( n ) );
    
    return dot( n, x, incx, y, incy );
}

inline
double
dotu ( const blas_int_t n,
       const double * x, const blas_int_t incx,
       const double * y, const blas_int_t incy )
{
    ADD_FLOPS( FLOPS_DDOT( n ) );
    
    return dot( n, x, incx, y, incy );
}

inline
std::complex< float >
dotu ( const blas_int_t  n,
       const std::complex< float > *  x, const blas_int_t  incx,
       const std::complex< float > *  y, const blas_int_t  incy )
{
    ADD_FLOPS( FLOPS_CDOT( n ) );
    
    std::complex< float >  res;
    
    #if HPRO_USE_MKL == 1

    cblas_cdotu_sub( n, x, incx, y, incy, & res );
    
    #else
    
    xcdotu_( & n, x, & incx, y, & incy, & res );
    
    #endif
    
    return res;
}

inline
std::complex< double >
dotu ( const blas_int_t  n,
       const std::complex< double > *  x, const blas_int_t  incx,
       const std::complex< double > *  y, const blas_int_t  incy )
{
    ADD_FLOPS( FLOPS_ZDOT( n ) );
    
    std::complex< double >  res;
    
    #if HPRO_USE_MKL == 1

    cblas_zdotu_sub( n, x, incx, y, incy, & res );
    
    #else
    
    xzdotu_( & n, x, & incx, y, & incy, & res );
    
    #endif
    
    return res;
}

//
// *norm2
//
#define HPRO_NORM2_FUNC( type, func, flops )                            \
    inline typename real_type< type >::type_t                           \
    norm2 ( const blas_int_t n, const type * x, const blas_int_t incx ) \
    {                                                                   \
        ADD_FLOPS( FLOPS_##flops( n ) );                                \
        return func( & n, x, & incx );                                  \
    }

HPRO_NORM2_FUNC( float,                  snrm2_,  SDOT )
HPRO_NORM2_FUNC( double,                 dnrm2_,  DDOT )
HPRO_NORM2_FUNC( std::complex< float >,  scnrm2_, CDOT )
HPRO_NORM2_FUNC( std::complex< double >, dznrm2_, ZDOT )

//
// *ger
//
#define HPRO_GER_FUNC( type, func, flops )                              \
    inline                                                              \
    void ger ( const blas_int_t n, const blas_int_t m, const type alpha, \
               const type * x, const blas_int_t incx,                   \
               const type * y, const blas_int_t incy,                   \
               type * A, const blas_int_t ldA )                         \
    {                                                                   \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        func( & n, & m, & alpha, x, & incx, y, & incy, A, & ldA );      \
    }

HPRO_GER_FUNC( float,                  sger_,  SGER )
HPRO_GER_FUNC( double,                 dger_ , DGER )
HPRO_GER_FUNC( std::complex< float >,  cgerc_, CGER )
HPRO_GER_FUNC( std::complex< double >, zgerc_, ZGER )

//
// *geru
//
#define HPRO_GERU_FUNC( type, func, flops )                             \
    inline                                                              \
    void geru ( const blas_int_t n, const blas_int_t m, const type alpha, \
               const type * x, const blas_int_t incx,                   \
               const type * y, const blas_int_t incy,                   \
               type * A, const blas_int_t ldA )                         \
    {                                                                   \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        func( & n, & m, & alpha, x, & incx, y, & incy, A, & ldA );      \
    }

HPRO_GERU_FUNC( float,                  sger_,  SGER )
HPRO_GERU_FUNC( double,                 dger_,  DGER )
HPRO_GERU_FUNC( std::complex< float >,  cgeru_, CGER )
HPRO_GERU_FUNC( std::complex< double >, zgeru_, ZGER )

//
// *gemv
//
#define HPRO_GEMV_FUNC( type, func, flops )                             \
    inline                                                              \
    void gemv ( const char trans, const blas_int_t n, const blas_int_t m, const type alpha, \
                const type * A, const blas_int_t ldA,                   \
                const type * x, const blas_int_t incx,                  \
                const type beta, type * y, const blas_int_t incy ) {    \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        char  ltrans = trans;                                           \
        if ( ! is_complex_type< type >::value && ( trans == 'C' ) )     \
            ltrans = 'T';                                               \
        func( & ltrans, & n, & m, & alpha, A, & ldA, x, & incx, & beta, y, & incy ); }

HPRO_GEMV_FUNC( float,                  sgemv_, SGEMV )
HPRO_GEMV_FUNC( double,                 dgemv_, DGEMV )
HPRO_GEMV_FUNC( std::complex< float >,  cgemv_, CGEMV )
HPRO_GEMV_FUNC( std::complex< double >, zgemv_, ZGEMV )

//
// *trmv
//
#define HPRO_TRMV_FUNC( type, func, flops )                             \
    inline                                                              \
    void trmv ( const char uplo, const char trans, const char diag,     \
                const blas_int_t n, const type * A, const blas_int_t ldA, \
                type * x, const blas_int_t incx ) {                     \
        ADD_FLOPS( FLOPS_##flops( n ) );                                \
        char  ltrans = trans;                                           \
        if ( ! is_complex_type< type >::value && ( trans == 'C' ) )     \
            ltrans = 'T';                                               \
        func( & uplo, & ltrans, & diag, & n, A, & ldA, x, & incx ); }

HPRO_TRMV_FUNC( float,                  strmv_, STRMV )
HPRO_TRMV_FUNC( double,                 dtrmv_, DTRMV )
HPRO_TRMV_FUNC( std::complex< float >,  ctrmv_, CTRMV )
HPRO_TRMV_FUNC( std::complex< double >, ztrmv_, ZTRMV )

//
// *trsv
//
#define HPRO_TRSV_FUNC( type, func, flops )                             \
    inline                                                              \
    void trsv ( const char uplo, const char trans, const char diag,     \
                const blas_int_t n, const type * A, const blas_int_t ldA, \
                type * b, const blas_int_t incb ) {                     \
        ADD_FLOPS( FLOPS_##flops( n ) );                                \
        char  ltrans = trans;                                           \
        if ( ! is_complex_type< type >::value && ( trans == 'C' ) )     \
            ltrans = 'T';                                               \
        func( & uplo, & ltrans, & diag, & n, A, & ldA, b, & incb ); }

HPRO_TRSV_FUNC( float,                  strsv_, STRSV )
HPRO_TRSV_FUNC( double,                 dtrsv_, DTRSV )
HPRO_TRSV_FUNC( std::complex< float >,  ctrsv_, CTRSV )
HPRO_TRSV_FUNC( std::complex< double >, ztrsv_, ZTRSV )

//
// *gemm
//
#define HPRO_GEMM_FUNC( type, func, flops )                             \
    inline                                                              \
    void gemm ( const char transA, const char transB,                   \
                const blas_int_t n, const blas_int_t m, const blas_int_t k, \
                const type alpha,                                       \
                const type * A, const blas_int_t ldA,                   \
                const type * B, const blas_int_t ldB,                   \
                const type beta, type * C, const blas_int_t ldC ) {     \
        ADD_FLOPS( FLOPS_##flops( n, m, k ) );                          \
        char  ltransA = transA;                                         \
        char  ltransB = transB;                                         \
        if ( ! is_complex_type< type >::value && ( transA == 'C' ) )    \
            ltransA = 'T';                                              \
        if ( ! is_complex_type< type >::value && ( transB == 'C' ) )    \
            ltransB = 'T';                                              \
        func( & ltransA, & ltransB, & n, & m, & k,                      \
              & alpha, A, & ldA, B, & ldB,                              \
              & beta, C, & ldC ); }

HPRO_GEMM_FUNC( float,                  sgemm_, SGEMM )
HPRO_GEMM_FUNC( double,                 dgemm_, DGEMM )
HPRO_GEMM_FUNC( std::complex< float >,  cgemm_, CGEMM )
HPRO_GEMM_FUNC( std::complex< double >, zgemm_, ZGEMM )

//
// *trmm
//
#define HPRO_TRMM_FUNC( type, func, flops )                             \
    inline                                                              \
    void trmm ( const char side, const char uplo, const char trans, const char diag, \
                const blas_int_t n, const blas_int_t m, const type alpha, const type * A, \
                const blas_int_t ldA, type * B, const blas_int_t ldB ) { \
        ADD_FLOPS( FLOPS_##flops( side, n, m ) );                       \
        char  ltrans = trans;                                           \
        if ( ! is_complex_type< type >::value && ( trans == 'C' ) )     \
            ltrans = 'T';                                               \
        func( & side, & uplo, & ltrans, & diag, & n, & m, & alpha, A, & ldA, B, & ldB ); }

HPRO_TRMM_FUNC( float,                  strmm_, STRMM )
HPRO_TRMM_FUNC( double,                 dtrmm_, CTRMM )
HPRO_TRMM_FUNC( std::complex< float >,  ctrmm_, DTRMM )
HPRO_TRMM_FUNC( std::complex< double >, ztrmm_, ZTRMM )

//
// *trsm
//
#define HPRO_TRSM_FUNC( type, func, flops )                             \
    inline                                                              \
    void trsm ( const char side, const char uplo, const char trans, const char diag, \
                const blas_int_t n, const blas_int_t m, const type alpha, const type * A, \
                const blas_int_t ldA, type * B, const blas_int_t ldB ) { \
        ADD_FLOPS( FLOPS_##flops( side, n, m ) );                       \
        char  ltrans = trans;                                           \
        if ( ! is_complex_type< type >::value && ( trans == 'C' ) )     \
            ltrans = 'T';                                               \
        func( & side, & uplo, & ltrans, & diag, & n, & m, & alpha, A, & ldA, B, & ldB ); }

HPRO_TRSM_FUNC( float,                  strsm_, STRSM )
HPRO_TRSM_FUNC( double,                 dtrsm_, CTRSM )
HPRO_TRSM_FUNC( std::complex< float >,  ctrsm_, DTRSM )
HPRO_TRSM_FUNC( std::complex< double >, ztrsm_, ZTRSM )



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// wrappers for LAPACK functions
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// *gesv
//
#define HPRO_GESV_FUNC( type, func )                                \
    inline void gesv ( const blas_int_t  n,                         \
                       const blas_int_t  nrhs,                      \
                       type *            A,                         \
                       const blas_int_t  ldA,                       \
                       blas_int_t *      ipiv,                      \
                       type *            B,                         \
                       const blas_int_t  ldB,                       \
                       blas_int_t &      info ) {                   \
        info = 0;                                                   \
        func( & n, & nrhs, A, & ldA, ipiv, B, & ldB, & info ); }

HPRO_GESV_FUNC( float,                  sgesv_ )
HPRO_GESV_FUNC( double,                 dgesv_ )
HPRO_GESV_FUNC( std::complex< float >,  cgesv_ )
HPRO_GESV_FUNC( std::complex< double >, zgesv_ )

#undef HPRO_GESV_FUNC

//
// *trtri
//
#define HPRO_TRTRI_FUNC( type, func )                                   \
    inline  void trtri ( const char        uplo,                        \
                         const char        diag,                        \
                         const blas_int_t  n,                           \
                         type *            A,                           \
                         const blas_int_t  ldA,                         \
                         blas_int_t &      info ) {                     \
        info = 0;                                                       \
        func( & uplo, & diag, & n, A, & ldA, & info ); }

HPRO_TRTRI_FUNC( float,                  strtri_ )
HPRO_TRTRI_FUNC( double,                 dtrtri_ )
HPRO_TRTRI_FUNC( std::complex< float >,  ctrtri_ )
HPRO_TRTRI_FUNC( std::complex< double >, ztrtri_ )

#undef HPRO_TRTRI_FUNC

//
// *getrf
//
#define HPRO_GETRF_FUNC( type, func, flops )                            \
    inline void getrf ( const blas_int_t  n,                            \
                        const blas_int_t  m,                            \
                        type *            A,                            \
                        const blas_int_t  ldA,                          \
                        blas_int_t *      ipiv,                         \
                        blas_int_t &      info ) {                      \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, ipiv, & info ); }

HPRO_GETRF_FUNC( float,                  sgetrf_, SGETRF )
HPRO_GETRF_FUNC( double,                 dgetrf_, CGETRF )
HPRO_GETRF_FUNC( std::complex< float >,  cgetrf_, DGETRF )
HPRO_GETRF_FUNC( std::complex< double >, zgetrf_, ZGETRF )

#undef HPRO_GETRF_FUNC

//
// *getri
//
#define HPRO_GETRI_FUNC( type, func, flops )                           \
    inline void getri ( const blas_int_t  n,                           \
                        type *            A,                           \
                        const blas_int_t  ldA,                         \
                        blas_int_t *      ipiv,                        \
                        type *            work,                        \
                        blas_int_t        lwork,                       \
                        blas_int_t &      info ) {                     \
        ADD_FLOPS( FLOPS_##flops( n ) );                               \
        info = 0;                                                      \
        func( & n, A, & ldA, ipiv, work, & lwork, & info ); }

HPRO_GETRI_FUNC( float,                  sgetri_, SGETRI )
HPRO_GETRI_FUNC( double,                 dgetri_, CGETRI )
HPRO_GETRI_FUNC( std::complex< float >,  cgetri_, DGETRI )
HPRO_GETRI_FUNC( std::complex< double >, zgetri_, ZGETRI )

#undef HPRO_GETRI_FUNC

//
// *geqrf
//
#define HPRO_GEQRF_FUNC( type, func, flops )                            \
    inline void geqrf ( const blas_int_t  n,                            \
                        const blas_int_t  m,                            \
                        type *            A,                            \
                        const blas_int_t  ldA,                          \
                        type *            tau,                          \
                        type *            work,                         \
                        blas_int_t        lwork,                        \
                        blas_int_t &      info ) {                      \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, tau, work, & lwork, & info ); }

HPRO_GEQRF_FUNC( float,                  sgeqrf_, SGEQRF )
HPRO_GEQRF_FUNC( double,                 dgeqrf_, CGEQRF )
HPRO_GEQRF_FUNC( std::complex< float >,  cgeqrf_, DGEQRF )
HPRO_GEQRF_FUNC( std::complex< double >, zgeqrf_, ZGEQRF )

#undef HPRO_GEQRF_FUNC

//
// *geqr2
//
#define HPRO_GEQR2_FUNC( type, func, flops )                            \
    inline void geqr2 ( const blas_int_t  n,                            \
                        const blas_int_t  m,                            \
                        type *            A,                            \
                        const blas_int_t  ldA,                          \
                        type *            tau,                          \
                        type *            work,                         \
                        blas_int_t &      info ) {                      \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, tau, work, & info ); }

HPRO_GEQR2_FUNC( float,                  sgeqr2_, SGEQRF )
HPRO_GEQR2_FUNC( double,                 dgeqr2_, CGEQRF )
HPRO_GEQR2_FUNC( std::complex< float >,  cgeqr2_, DGEQRF )
HPRO_GEQR2_FUNC( std::complex< double >, zgeqr2_, ZGEQRF )

#undef HPRO_GEQR2_FUNC

//
// *orgqr
//
#define HPRO_ORGQR_FUNC( type, func, flops )                            \
    inline void orgqr ( const blas_int_t  n,                            \
                        const blas_int_t  m,                            \
                        const blas_int_t  k,                            \
                        type *            A,                            \
                        const blas_int_t  ldA,                          \
                        const type *      tau,                          \
                        type *            work,                         \
                        blas_int_t        lwork,                        \
                        blas_int_t &      info ) {                      \
        ADD_FLOPS( FLOPS_##flops( n, m, k ) );                          \
        info = 0;                                                       \
        func( & n, & m, & k, A, & ldA, tau, work, & lwork, & info ); }

HPRO_ORGQR_FUNC( float,             sorgqr_, SORGQR )
HPRO_ORGQR_FUNC( double,            dorgqr_, DORGQR )
HPRO_ORGQR_FUNC( std::complex< float >,  cungqr_, CUNGQR )
HPRO_ORGQR_FUNC( std::complex< double >, zungqr_, ZUNGQR )

#undef HPRO_ORGQR_FUNC

//
// *org2r
//
#define HPRO_ORG2R_FUNC( type, func, flops )                            \
    inline void org2r ( const blas_int_t  n,                            \
                        const blas_int_t  m,                            \
                        const blas_int_t  k,                            \
                        type *            A,                            \
                        const blas_int_t  ldA,                          \
                        const type *      tau,                          \
                        type *            work,                         \
                        blas_int_t &      info ) {                      \
        ADD_FLOPS( FLOPS_##flops( n, m, k ) );                          \
        info = 0;                                                       \
        func( & n, & m, & k, A, & ldA, tau, work, & info ); }

HPRO_ORG2R_FUNC( float,                  sorg2r_, SORG2R )
HPRO_ORG2R_FUNC( double,                 dorg2r_, DORG2R )
HPRO_ORG2R_FUNC( std::complex< float >,  cung2r_, CUNG2R )
HPRO_ORG2R_FUNC( std::complex< double >, zung2r_, ZUNG2R )

#undef HPRO_ORG2R_FUNC

//
// *geqp3
//
#define HPRO_GEQP3_FUNC( type, func )                                   \
    inline void geqp3 ( const blas_int_t  n,                            \
                        const blas_int_t  m,                            \
                        type *            A,                            \
                        const blas_int_t  ldA,                          \
                        blas_int_t *      jpvt,                         \
                        type *            tau,                          \
                        type *            work,                         \
                        blas_int_t        lwork,                        \
                        typename real_type< type >::type_t *,           \
                        blas_int_t &      info ) {                      \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, jpvt, tau, work, & lwork, & info ); }

HPRO_GEQP3_FUNC( float,  sgeqp3_ )
HPRO_GEQP3_FUNC( double, dgeqp3_ )

#undef HPRO_GEQP3_FUNC

#define HPRO_GEQP3_FUNC( type, func )                                   \
    inline void geqp3 ( const blas_int_t  n,                            \
                        const blas_int_t  m,                            \
                        type *            A,                            \
                        const blas_int_t  ldA,                          \
                        blas_int_t *      jpvt,                         \
                        type *            tau,                          \
                        type *            work,                         \
                        blas_int_t        lwork,                        \
                        typename real_type< type >::type_t *  rwork,    \
                        blas_int_t &      info ) {                      \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, jpvt, tau, work, & lwork, rwork, & info ); }

HPRO_GEQP3_FUNC( std::complex< float >,  cgeqp3_ )
HPRO_GEQP3_FUNC( std::complex< double >, zgeqp3_ )

#undef HPRO_GEQP3_FUNC

#if HPRO_HAS_GEQP3_TRUNC == 1

//
// *geqp3trunc
//
#define HPRO_GEQP3TRUNC_FUNC( type, func )                              \
    inline void geqp3trunc ( const blas_int_t  n,                       \
                             const blas_int_t  m,                       \
                             type *            A,                       \
                             const blas_int_t  ldA,                     \
                             blas_int_t *      jpvt,                    \
                             type *            tau,                     \
                             blas_int_t &      ntrunc,                  \
                             typename real_type< type >::type_t  atrunc, \
                             typename real_type< type >::type_t  rtrunc, \
                             type *            work,                    \
                             blas_int_t        lwork,                   \
                             typename real_type< type >::type_t *,      \
                             blas_int_t &      info ) {                 \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, jpvt, tau, & ntrunc, & atrunc, & rtrunc, work, & lwork, & info ); \
    }

HPRO_GEQP3TRUNC_FUNC( float,  sgeqp3trunc_ )
HPRO_GEQP3TRUNC_FUNC( double, dgeqp3trunc_ )

#undef HPRO_GEQP3TRUNC_FUNC

#define HPRO_GEQP3TRUNC_FUNC( type, func )                              \
    inline void geqp3trunc ( const blas_int_t  n,                       \
                             const blas_int_t  m,                       \
                             type *            A,                       \
                             const blas_int_t  ldA,                     \
                             blas_int_t *      jpvt,                    \
                             type *            tau,                     \
                             blas_int_t &      ntrunc,                  \
                             typename real_type< type >::type_t  atrunc, \
                             typename real_type< type >::type_t  rtrunc, \
                             type *            work,                    \
                             blas_int_t        lwork,                   \
                             typename real_type< type >::type_t *  rwork, \
                             blas_int_t &      info ) {                 \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, jpvt, tau, & ntrunc, & atrunc, & rtrunc, work, & lwork, rwork, & info ); }

HPRO_GEQP3TRUNC_FUNC( std::complex< float >,  cgeqp3trunc_ )
HPRO_GEQP3TRUNC_FUNC( std::complex< double >, zgeqp3trunc_ )

#undef HPRO_GEQP3TRUNC_FUNC

#endif

//
// *gelqf
//
#define HPRO_GELQF_FUNC( type, func, flops )                            \
    inline void gelqf ( const blas_int_t  n,                            \
                        const blas_int_t  m,                            \
                        type *            A,                            \
                        const blas_int_t  ldA,                          \
                        type *            tau,                          \
                        type *            work,                         \
                        blas_int_t        lwork,                        \
                        blas_int_t &      info ) {                      \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, tau, work, & lwork, & info ); }

HPRO_GELQF_FUNC( float,                  sgelqf_, SGELQF )
HPRO_GELQF_FUNC( double,                 dgelqf_, CGELQF )
HPRO_GELQF_FUNC( std::complex< float >,  cgelqf_, DGELQF )
HPRO_GELQF_FUNC( std::complex< double >, zgelqf_, ZGELQF )

#undef HPRO_GELQF_FUNC

//
// *gelq2
//
#define HPRO_GELQ2_FUNC( type, func, flops )                            \
    inline void gelq2 ( const blas_int_t  n,                            \
                        const blas_int_t  m,                            \
                        type *            A,                            \
                        const blas_int_t  ldA,                          \
                        type *            tau,                          \
                        type *            work,                         \
                        blas_int_t &      info ) {                      \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, tau, work, & info ); }

HPRO_GELQ2_FUNC( float,                  sgelq2_, SGELQF )
HPRO_GELQ2_FUNC( double,                 dgelq2_, CGELQF )
HPRO_GELQ2_FUNC( std::complex< float >,  cgelq2_, DGELQF )
HPRO_GELQ2_FUNC( std::complex< double >, zgelq2_, ZGELQF )

#undef HPRO_GELQ2_FUNC

//
// *orglq
//
#define HPRO_ORGLQ_FUNC( type, func, flops )                            \
    inline void orglq ( const blas_int_t  n,                            \
                        const blas_int_t  m,                            \
                        const blas_int_t  k,                            \
                        type *            A,                            \
                        const blas_int_t  ldA,                          \
                        const type *      tau,                          \
                        type *            work,                         \
                        blas_int_t        lwork,                        \
                        blas_int_t &      info ) {                      \
        ADD_FLOPS( FLOPS_##flops( n, m, k ) );                          \
        info = 0;                                                       \
        func( & n, & m, & k, A, & ldA, tau, work, & lwork, & info ); }

HPRO_ORGLQ_FUNC( float,             sorglq_, SORGLQ )
HPRO_ORGLQ_FUNC( double,            dorglq_, DORGLQ )
HPRO_ORGLQ_FUNC( std::complex< float >,  cunglq_, CUNGLQ )
HPRO_ORGLQ_FUNC( std::complex< double >, zunglq_, ZUNGLQ )

#undef HPRO_ORGLQ_FUNC

//
// *org2l
//
#define HPRO_ORG2L_FUNC( type, func, flops )                            \
    inline void org2l ( const blas_int_t  n,                            \
                        const blas_int_t  m,                            \
                        const blas_int_t  k,                            \
                        type *            A,                            \
                        const blas_int_t  ldA,                          \
                        const type *      tau,                          \
                        type *            work,                         \
                        blas_int_t &      info ) {                      \
        ADD_FLOPS( FLOPS_##flops( n, m, k ) );                          \
        info = 0;                                                       \
        func( & n, & m, & k, A, & ldA, tau, work, & info ); }

HPRO_ORG2L_FUNC( float,                  sorg2l_, SORGQL )
HPRO_ORG2L_FUNC( double,                 dorg2l_, DORGQL )
HPRO_ORG2L_FUNC( std::complex< float >,  cung2l_, CUNGQL )
HPRO_ORG2L_FUNC( std::complex< double >, zung2l_, ZUNGQL )

#undef HPRO_ORG2R_FUNC

//
// *syev/*heev
//
#define HPRO_HEEV_FUNC( type, func )                                    \
    inline void heev ( const char            jobz,                      \
                       const char            uplo,                      \
                       const blas_int_t      n,                         \
                       type *                A,                         \
                       const blas_int_t      ldA,                       \
                       real_type_t< type > * W,                         \
                       type *                work,                      \
                       const blas_int_t      lwork,                     \
                       real_type_t< type > *,                           \
                       blas_int_t &          info ) {                   \
        info = 0;                                                       \
        func( & jobz, & uplo, & n, A, & ldA, W, work, & lwork, & info ); }

HPRO_HEEV_FUNC( float,  ssyev_ )
HPRO_HEEV_FUNC( double, dsyev_ )

#undef HPRO_HEEV_FUNC

#define HPRO_HEEV_FUNC( type, func )                                    \
    inline void heev ( const char            jobz,                      \
                       const char            uplo,                      \
                       const blas_int_t      n,                         \
                       type *                A,                         \
                       const blas_int_t      ldA,                       \
                       real_type_t< type > * W,                         \
                       type *                work,                      \
                       const blas_int_t      lwork,                     \
                       real_type_t< type > * rwork,                     \
                       blas_int_t & info ) {                            \
        info = 0;                                                       \
        func( & jobz, & uplo, & n, A, & ldA, W, work, & lwork, rwork, & info ); }

HPRO_HEEV_FUNC( std::complex< float >,  cheev_ )
HPRO_HEEV_FUNC( std::complex< double >, zheev_ )

#undef HPRO_HEEV_FUNC

//
// *syevd/*heevd
//
#define HPRO_HEEVD_FUNC( type, func )                                   \
    inline void heevd ( const char            jobz,                     \
                        const char            uplo,                     \
                        const blas_int_t      n,                        \
                        type *                A,                        \
                        const blas_int_t      ldA,                      \
                        real_type_t< type > * W,                        \
                        type *                work,                     \
                        const blas_int_t      lwork,                    \
                        real_type_t< type > *,                          \
                        const blas_int_t      ,                         \
                        const blas_int_t *    iwork,                    \
                        const blas_int_t      liwork,                   \
                        blas_int_t &          info ) {                  \
        info = 0;                                                       \
        func( & jobz, & uplo, & n, A, & ldA, W, work, & lwork, iwork, & liwork, & info ); }

HPRO_HEEVD_FUNC( float,  ssyevd_ )
HPRO_HEEVD_FUNC( double, dsyevd_ )

#undef HPRO_HEEVD_FUNC

#define HPRO_HEEVD_FUNC( type, func )                                   \
    inline void heevd ( const char            jobz,                     \
                        const char            uplo,                     \
                        const blas_int_t      n,                        \
                        type *                A,                        \
                        const blas_int_t      ldA,                      \
                        real_type_t< type > * W,                        \
                        type *                work,                     \
                        const blas_int_t      lwork,                    \
                        real_type_t< type > * rwork,                    \
                        const blas_int_t      lrwork,                   \
                        const blas_int_t *    iwork,                    \
                        const blas_int_t      liwork,                   \
                        blas_int_t & info ) {                           \
        info = 0;                                                       \
        func( & jobz, & uplo, & n, A, & ldA, W, work, & lwork, rwork, & lrwork, iwork, & liwork, & info ); }

HPRO_HEEVD_FUNC( std::complex< float >,  cheevd_ )
HPRO_HEEVD_FUNC( std::complex< double >, zheevd_ )

#undef HPRO_HEEVD_FUNC

//
// *syevx/*heevx
//
#define HPRO_HEEVX_FUNC( type, func )                                   \
    inline void heevx ( const char                  jobz,               \
                        const char                  uplo,               \
                        const blas_int_t            n,                  \
                        type *                      A,                  \
                        const blas_int_t            ldA,                \
                        const blas_int_t            il,                 \
                        const blas_int_t            iu,                 \
                        blas_int_t &                m,                  \
                        real_type< type >::type_t *  W,                 \
                        type *                      Z,                  \
                        const blas_int_t            ldZ,                \
                        type *                      work,               \
                        const blas_int_t            lwork,              \
                        blas_int_t *                iwork,              \
                        blas_int_t &                info ) {            \
        char                      range = 'I';                          \
        real_type< type >::type_t  vl, vu;                              \
        real_type< type >::type_t  abstol = 0;                          \
        blas_int_t                       ifail  = 0;                    \
        info = 0;                                                       \
        func( & jobz, & range, & uplo, & n, A, & ldA, & vl, & vu, & il, & iu, \
              & abstol, & m, W, Z, & ldZ, work, & lwork, iwork, & ifail, & info ); }

HPRO_HEEVX_FUNC( float,  ssyevx_ )
HPRO_HEEVX_FUNC( double, dsyevx_ )

#undef HPRO_HEEVX_FUNC

//
// sstev/dstev
//
#define HPRO_STEV_FUNC( type, func )                                  \
    inline void stev ( const char        jobz,                        \
                       const blas_int_t  n,                           \
                       type *            D,                           \
                       type *            E,                           \
                       type *            Z,                           \
                       const blas_int_t  ldZ,                         \
                       type *            work,                        \
                       blas_int_t &      info ) {                     \
        info = 0;                                                     \
        func( & jobz, & n, D, E, Z, & ldZ, work, & info ); }

HPRO_STEV_FUNC( float,  sstev_ )
HPRO_STEV_FUNC( double, dstev_ )

#undef HPRO_STEV_FUNC

//
// *gesvd
//
#define HPRO_GESVD_FUNC( type, func, flops )                            \
    inline void gesvd ( const char                  jobu,               \
                        const char                  jobv,               \
                        const blas_int_t            n,                  \
                        const blas_int_t            m,                  \
                        type *                      A,                  \
                        const blas_int_t            ldA,                \
                        real_type< type >::type_t *  S,                 \
                        type *                      U,                  \
                        const blas_int_t            ldU,                \
                        type *                      VT,                 \
                        const blas_int_t            ldVT,               \
                        type *                      work,               \
                        const blas_int_t            lwork,              \
                        real_type< type >::type_t *,                    \
                        blas_int_t &                info ) {            \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & jobu, & jobv, & n, & m, A, & ldA, S, U, & ldU, VT, & ldVT, work, & lwork, & info ); }

HPRO_GESVD_FUNC( float,  sgesvd_, SGESVD )
HPRO_GESVD_FUNC( double, dgesvd_, DGESVD )

#undef  HPRO_GESVD_FUNC

#define HPRO_GESVD_FUNC( type, func, flops )                            \
    inline void gesvd ( const char                  jobu,               \
                        const char                  jobv,               \
                        const blas_int_t            n,                  \
                        const blas_int_t            m,                  \
                        type *                      A,                  \
                        const blas_int_t            ldA,                \
                        real_type< type >::type_t *  S,                 \
                        type *                      U,                  \
                        const blas_int_t            ldU,                \
                        type *                      VT,                 \
                        const blas_int_t            ldVT,               \
                        type *                      work,               \
                        const blas_int_t            lwork,              \
                        real_type< type >::type_t *  rwork,             \
                        blas_int_t &                info ) {            \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & jobu, & jobv, & n, & m, A, & ldA, S, U, & ldU, VT, & ldVT, work, & lwork, rwork, & info ); }

HPRO_GESVD_FUNC( std::complex< float >,  cgesvd_, CGESVD )
HPRO_GESVD_FUNC( std::complex< double >, zgesvd_, CGESVD )

#undef  HPRO_GESVD_FUNC

//
// *gesdd
//
#define HPRO_GESDD_FUNC( type, func, flops )                            \
    inline void gesdd ( const char                  jobz,               \
                        const blas_int_t            n,                  \
                        const blas_int_t            m,                  \
                        type *                      A,                  \
                        const blas_int_t            ldA,                \
                        real_type< type >::type_t *  S,                 \
                        type *                      U,                  \
                        const blas_int_t            ldU,                \
                        type *                      VT,                 \
                        const blas_int_t            ldVT,               \
                        type *                      work,               \
                        const blas_int_t            lwork,              \
                        real_type< type >::type_t *,                    \
                        blas_int_t *                iwork,              \
                        blas_int_t &                info ) {            \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & jobz, & n, & m, A, & ldA, S, U, & ldU, VT, & ldVT, work, & lwork, iwork, & info ); }

HPRO_GESDD_FUNC( float,  sgesdd_, SGESDD )
HPRO_GESDD_FUNC( double, dgesdd_, DGESDD )

#undef  HPRO_GESDD_FUNC

#define HPRO_GESDD_FUNC( type, func, flops )                            \
    inline void gesdd ( const char                  jobz,               \
                        const blas_int_t            n,                  \
                        const blas_int_t            m,                  \
                        type *                      A,                  \
                        const blas_int_t            ldA,                \
                        real_type< type >::type_t *  S,                 \
                        type *                      U,                  \
                        const blas_int_t            ldU,                \
                        type *                      VT,                 \
                        const blas_int_t            ldVT,               \
                        type *                      work,               \
                        const blas_int_t            lwork,              \
                        real_type< type >::type_t *  rwork,             \
                        blas_int_t *                iwork,              \
                        blas_int_t &                info ) {            \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & jobz, & n, & m, A, & ldA, S, U, & ldU, VT, & ldVT, work, & lwork, rwork, iwork, & info ); }

HPRO_GESDD_FUNC( std::complex< float >,  cgesdd_, CGESDD )
HPRO_GESDD_FUNC( std::complex< double >, zgesdd_, ZGESDD )

#undef  HPRO_GESDD_FUNC

#define HPRO_GESVJ_FUNC( type, func )                                   \
    inline void gesvj ( const char                  joba,               \
                        const char                  jobu,               \
                        const char                  jobv,               \
                        const blas_int_t            n,                  \
                        const blas_int_t            m,                  \
                        type *                      A,                  \
                        const blas_int_t            ldA,                \
                        real_type< type >::type_t * S,                  \
                        const blas_int_t            mv,                 \
                        type *                      V,                  \
                        const blas_int_t            ldV,                \
                        type *                      cwork,              \
                        const blas_int_t            lwork,              \
                        real_type< type >::type_t *,                    \
                        const blas_int_t,                               \
                        blas_int_t &                info ) {            \
        info = 0;                                                       \
        func( & joba, & jobu, & jobv, & n, & m, A, & ldA, S, & mv, V, & ldV, cwork, & lwork, & info ); }

HPRO_GESVJ_FUNC( float,  sgesvj_ )
HPRO_GESVJ_FUNC( double, dgesvj_ )

#undef  HPRO_GESVJ_FUNC

#define HPRO_GESVJ_FUNC( type, func )                                   \
    inline void gesvj ( const char                  joba,               \
                        const char                  jobu,               \
                        const char                  jobv,               \
                        const blas_int_t            n,                  \
                        const blas_int_t            m,                  \
                        type *                      A,                  \
                        const blas_int_t            ldA,                \
                        real_type< type >::type_t * S,                  \
                        const blas_int_t            mv,                 \
                        type *                      V,                  \
                        const blas_int_t            ldV,                \
                        type *                      cwork,              \
                        const blas_int_t            lwork,              \
                        real_type< type >::type_t * rwork,              \
                        const blas_int_t            lrwork,             \
                        blas_int_t &                info ) {            \
        info = 0;                                                       \
        func( & joba, & jobu, & jobv, & n, & m, A, & ldA, S, & mv, V, & ldV, cwork, & lwork, rwork, & lrwork, & info ); }

HPRO_GESVJ_FUNC( std::complex< float >,  cgesvj_ )
HPRO_GESVJ_FUNC( std::complex< double >, zgesvj_ )

#undef  HPRO_GESVJ_FUNC

//
// *larfg
//
#define HPRO_LARFG_FUNC( type, func )                                   \
    inline void larfg ( const blas_int_t  n,                            \
                        type &            alpha,                        \
                        type *            x,                            \
                        const blas_int_t  incx,                         \
                        type &            tau ) {                       \
        func( & n, & alpha, x, & incx, & tau ); }

HPRO_LARFG_FUNC( float,                  slarfg_ )
HPRO_LARFG_FUNC( double,                 dlarfg_ )
HPRO_LARFG_FUNC( std::complex< float >,  clarfg_ )
HPRO_LARFG_FUNC( std::complex< double >, zlarfg_ )

#undef HPRO_LARFG_FUNC

//
// *larf
//
// apply householder reflection
#define HPRO_LARF_FUNC( type, func )                                    \
    inline void larf ( const char        side,                          \
                       const blas_int_t  n,                             \
                       const blas_int_t  m,                             \
                       const type *      V,                             \
                       const blas_int_t  incv,                          \
                       const type        tau,                           \
                       type *            C,                             \
                       const blas_int_t  ldc,                           \
                       type *            work ) {                       \
        func( & side, & n, & m, V, & incv, & tau, C, & ldc, work ); }

HPRO_LARF_FUNC( float,                  slarf_ )
HPRO_LARF_FUNC( double,                 dlarf_ )
HPRO_LARF_FUNC( std::complex< float >,  clarf_ )
HPRO_LARF_FUNC( std::complex< double >, zlarf_ )

#undef HPRO_LARF_FUNC

//! @endcond

}// namespace BLAS

}// namespace Hpro

#endif // __HPRO_BLAS_LAPACK_HH
