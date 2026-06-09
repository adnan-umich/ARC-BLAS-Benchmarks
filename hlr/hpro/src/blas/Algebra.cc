//
// Project     : HLIBpro
// File        : Algebra.cc
// Description : provide linear algebra functions for BLAS matrices and vectors
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>
#include <list>
#include <atomic>
#include <random>
#include <tuple>
#include <mutex>

#include <hpro/base/config.hh>
#include <hpro/blas/gpu.hh>

#include "TRNG.hh"

#include <hpro/blas/Algebra.hh>

namespace Hpro
{

namespace BLAS
{

// global flop counter
std::atomic< std::uint64_t >  FLOPS( 0 );

namespace
{

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// local defines and constants
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// set to 1 to replace very small values with zero
#define CHECK_SMALL        0

// set to 1 to enable ACA based SVD
#define USE_ACASVD         0

// set to 1 to check ACA SVD for approx. errors
#define CHECK_ACASVD       0

// set to 1 to check LAPACK SVD for approx. errors
#define CHECK_GESVD        0

// set to 1 to check low-rank truncation for approx. errors
#define CHECK_TRUNCATE     0

// set to 1 to activate statistics of algebra functions
#define FUNC_STAT          0

// LAPACK constant for a workspace query
const blas_int_t           LAPACK_WS_QUERY = blas_int_t(-1);

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// timing counters for various LAPACK functions
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

namespace LOGGING
{

#if FUNC_STAT == 1

#define  HAS_LOGGING

using  mutex_t = std::mutex;
using  lock_t  = std::scoped_lock;

#  define  DEF_COUNTER( counter ) \
    std::atomic< size_t >  n##counter( 0 ); \
    mutex_t                m##counter; \
    double                 t##counter = 0

#  define LOG_COUNTER( counter )  { LOGGING::n##counter += 1; }
#  define LOG_TIC( timer )        auto  tic_##timer = Time::Wall::now();
#  define LOG_TOC( timer )        {                                     \
        auto  toc_##timer = Time::Wall::now();                          \
        LOGGING::lock_t lock_##timer( LOGGING::m##timer );              \
        LOGGING::t##timer += (toc_##timer - tic_##timer); }

DEF_COUNTER( approx );
DEF_COUNTER( trunc );
DEF_COUNTER( gesvd );
DEF_COUNTER( gesdd );
DEF_COUNTER( gesvj );
DEF_COUNTER( qr );
DEF_COUNTER( gemm );

#else

#  define LOG_COUNTER( counter )
#  define LOG_TIC( timer )      
#  define LOG_TOC( timer )      

#endif

}// namespace LOGGING

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// local functions
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// replace small values by zero
//
#if CHECK_SMALL == 1
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
check_small ( T1 &  A )
{
    using  real_t = typename real_type< T >::type_t;

    const real_t  eps = Math::square( Math::square( Limits::epsilon< real_t >() ) );
    
    for ( idx_t  j = 0; j < idx_t( A.ncols() ); ++j )
        for ( idx_t  i = 0; i < idx_t( A.nrows() ); ++i )
            if ( Math::abs( A[i] ) < eps )
                A(i,j) = 0.0;
}
#endif

//
// check for invalid data
//
//#if CHECK_INF_NAN == 1
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
check_inf_nan ( const T1 &    A,
                const char *  fnname,
                const char *  msg )
{
    for ( idx_t  j = 0; j < idx_t( A.ncols() ); ++j )
        for ( idx_t  i = 0; i < idx_t( A.nrows() ); ++i )
        {
            if ( Math::is_inf( A(i,j) ) )
                HERROR( ERR_INF, fnname, msg );
            
            if ( Math::is_nan( A(i,j) ) )
                HERROR( ERR_NAN, fnname, msg );
        }// for
}
//#endif

//
// defines for the above functions
//

#if CHECK_SMALL == 0
#define DO_CHECK_SMALL( A )
#else
#define DO_CHECK_SMALL( A )  check_small( (A) );
#endif

// #if CHECK_INF_NAN == 0
// #define DO_CHECK_INF_NAN( A, fn, msg )
// #else
#define DO_CHECK_INF_NAN( A, fn, msg )  if ( CFG::BLAS::check_inf_nan ) { check_inf_nan( (A), (fn), (msg) ); }
//#endif

// local random number generator
std::random_device                rd{};
std::mt19937                      generator{ rd() };
std::normal_distribution<>        normal_distr{ 0, 1 };
std::uniform_real_distribution<>  uniform_distr{ 0, 1 };

//
// return random number in range [-max,max]
// (range applies also to imaginary part)
//
template < typename T >
T
rand ( const typename real_type< T >::type_t  max )
{
    // return 2.0 * distr( mt ) * max - max;
    return 2.0 * T( uniform_distr( generator ) ) * max - max;
}

template <>
std::complex< float >
rand< std::complex< float > > ( const float  max )
{
    return std::complex< float >( rand< float >( max ), rand< float >( max ) );
}

}// namespace anonymous

//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// Basic Matrix algebra
//
////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

template < typename T1 >
void
fill_rand ( Matrix< T1 > &  M )
{
    const idx_t  n = idx_t(M.nrows());
    const idx_t  m = idx_t(M.ncols());
    
    for ( idx_t j = 0; j < m; ++j )
        for ( idx_t i = 0; i < n; ++i )
            M(i,j) = rand< T1 >( typename real_type< T1 >::type_t( 1 ) );
}

template < typename T1 >
void
fill_rand ( Matrix< T1 > &  M,
            TRNG &          rng )
{
    const idx_t  n = idx_t(M.nrows());
    const idx_t  m = idx_t(M.ncols());
    
    for ( idx_t j = 0; j < m; ++j )
        for ( idx_t i = 0; i < n; ++i )
            M(i,j) = rng.rand< T1 >( typename real_type< T1 >::type_t( 1 ) );
}

// local random number generator
template < typename value_t >
void
fill_rand_normal ( Vector< value_t > &  v )
{
    const idx_t  n = idx_t(v.length());

    std::normal_distribution<>  distr{ 0, 1 };
    
    for ( idx_t i = 0; i < n; ++i )
        v(i) = value_t( distr( generator ) );
}

template < typename value_t >
void
fill_rand_normal ( Matrix< value_t > &  M )
{
    const idx_t  n = idx_t(M.nrows());
    const idx_t  m = idx_t(M.ncols());

    std::normal_distribution<>  distr{ 0, 1 };
    
    for ( idx_t j = 0; j < m; ++j )
        for ( idx_t i = 0; i < n; ++i )
            M(i,j) = value_t( distr( generator ) );
}

//
// return spectral norm of M
//
template < typename value_t >
real_type_t< value_t >
norm2 ( const Matrix< value_t > &  M )
{
    using  real_t = real_type_t< value_t >;

    const size_t  minnm = std::min( M.nrows(), M.ncols() );

    if ( minnm == 0 )
        return real_t(0);
    
    Matrix< value_t >  A( M.nrows(), M.ncols() );
    Vector< real_t >   S( minnm );

    copy( M, A );
    sv( A, S );

    if ( S.length() == 0 )
        return real_t(0);
    
    real_t  sv_max = S(0);

    for ( idx_t  i = 1; i < idx_t( S.length() ); ++i )
    {
        if ( S(i) == real_t(0) )
            continue;

        sv_max = std::max( sv_max, S(i) );
    }// for

    return sv_max;
}

//
// return condition of M
//
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value,
                  typename real_type< typename T1::value_t >::type_t >
cond ( const T1 &  M )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    const size_t  minnm = std::min( M.nrows(), M.ncols() );

    if ( minnm == 0 )
        return real_t(0);
    
    Matrix< value_t >  A( M.nrows(), M.ncols() );
    Vector< real_t >   S( minnm );

    copy( M, A );
    sv( A, S );

    if ( S.length() == 0 )
        return real_t(0);
    
    real_t  sv_max = S(0);
    real_t  sv_min = S(0);

    for ( idx_t  i = 1; i < idx_t( S.length() ); ++i )
    {
        if ( S(i) == real_t(0) )
            continue;

        sv_max = std::max( sv_max, S(i) );
        sv_min = std::min( sv_min, S(i) );
    }// for

    if ( sv_min == real_t(0) )
        return real_t(0);
    
    return sv_max / sv_min;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// Advanced Matrix algebra
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//!
//! invert matrix \a A; \a A will be overwritten with A^-1 upon exit
//!
template < typename value_t >
void
invert ( Matrix< value_t > &  A )
{
    // using  value_t = typename T1::value_t;

    if ( A.nrows() != A.ncols() )
        HERROR( ERR_MAT_SIZE, "(BLAS) invert", "matrix is not square" );
    
    blas_int_t                  info = 0;
    Hpro::vector< blas_int_t >  ipiv( A.nrows() );
    
    DO_CHECK_INF_NAN( A, "(BLAS) invert", "in input matrix" );

    MKL_SEQ_START;
    
    getrf( blas_int_t( A.nrows() ), blas_int_t( A.nrows() ), A.data(),
           blas_int_t( A.col_stride() ), ipiv.data(), info );
    
    if      ( info < 0 ) HERROR( ERR_ARG, "(BLAS) invert",
                                 to_string( "argument %d LAPACK::getrf", -info ) );
    else if ( info > 0 ) HERROR( ERR_MAT_SINGULAR, "(BLAS) invert",
                                 to_string( "in LAPACK::getrf (info = %d)", info ) );
    
    DO_CHECK_INF_NAN( A, "(BLAS) invert", "after getrf" );

    //
    // workspace query
    //
    
    value_t  work_query = value_t(0);

    getri( blas_int_t( A.nrows() ), A.data(), blas_int_t( A.col_stride() ),
           ipiv.data(), & work_query, LAPACK_WS_QUERY, info );

    if ( info < 0 ) HERROR( ERR_ARG, "(BLAS) invert",
                            to_string( "argument %d LAPACK::getri", -info ) );
    
    const blas_int_t         lwork = blas_int_t( std::real( work_query ) );
    Hpro::vector< value_t >  work( lwork, value_t(0) );
    
    getri( blas_int_t( A.nrows() ), A.data(), blas_int_t( A.col_stride() ),
           ipiv.data(), work.data(), lwork, info );

    MKL_SEQ_END;
    
    if      ( info < 0 ) HERROR( ERR_ARG, "(BLAS) invert",
                                 to_string( "argument %d LAPACK::getri", -info ) );
    else if ( info > 0 ) HERROR( ERR_MAT_SINGULAR, "(BLAS) invert",
                                 to_string( "in LAPACK::getri (info = %d)", info ) );
    
    DO_CHECK_INF_NAN( A, "(BLAS) invert", "in output matrix" );
}

//!
//! invert lower or upper triangular matrix \a A with unit or non-unit diagonal;
//! \a A will be overwritten with A^-1 upon exit
//!
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
invert ( T1 &               A,
         const tri_type_t   tri_type,
         const diag_type_t  diag_type )
{
    if ( A.nrows() != A.ncols() )
        HERROR( ERR_MAT_SIZE, "(BLAS) invert", "matrix is not square" );
    
    blas_int_t  info = 0;
    
    DO_CHECK_INF_NAN( A, "(BLAS) invert", "in input matrix" );
    
    MKL_SEQ_START;
    
    trtri( char(tri_type), char(diag_type), blas_int_t( A.nrows() ), A.data(),
           blas_int_t( A.col_stride() ), info );

    MKL_SEQ_END;
    
    if      ( info < 0 ) HERROR( ERR_ARG, "(BLAS) invert",
                                 to_string( "argument %d LAPACK::trtri", -info ) );
    else if ( info > 0 ) HERROR( ERR_MAT_SINGULAR, "(BLAS) invert",
                                 to_string( "in LAPACK::trtri (info = %d)", info ) );
    
    DO_CHECK_INF_NAN( A, "(BLAS) invert", "in output matrix" );
}

//!
//! \ingroup  BLAS_Module
//! \brief    compute pseudo inverse of matrix \a A with precision \a acc
//!
//! \detail   Compute pseudo inverse B of matrix \a A up to precision \a acc,
//!           e.g. \f$\|A-B\|\le \epsilon\f$ with \f$\epsilon\f$ defined by \a acc.
//!
template <typename T>
void
pseudo_invert ( Matrix< T > &  A )
{
    using  value_t = T;
    using  real_t  = typename real_type< value_t >::type_t;

    DO_CHECK_INF_NAN( A, "(BLAS) pseudo_invert", "in input matrix" );

    const auto         n = idx_t(A.nrows());
    const auto         m = idx_t(A.ncols());
    const auto         min_nm = std::min( n, m );
    Matrix< value_t >  U( n, m );
    Matrix< value_t >  V( m, min_nm );
    Vector< real_t >   S( min_nm );
            
    copy( A, U );
    svd( U, S, V );

    for ( idx_t  i = 0; i < min_nm; ++i )
    {
        if ( S(i) != real_t(0) )
            S(i) = real_t(1) / S(i);
    }// for
    
    prod_diag( V, S, min_nm );

    if ( n != m )
        A = std::move( Matrix< value_t >( m, n ) );
    
    prod( value_t(1), V, adjoint(U), value_t(0), A );

    DO_CHECK_INF_NAN( A, "(BLAS) pseudo_invert", "in output matrix" );
}

//!
//! \ingroup  BLAS_Module
//! \brief    compute pseudo inverse of matrix \a A with precision \a acc
//!
//! \detail   Compute pseudo inverse B of matrix \a A up to precision \a acc,
//!           e.g. \f$\|A-B\|\le \epsilon\f$ with \f$\epsilon\f$ defined by \a acc.
//!
template <typename T>
void
pseudo_invert ( Matrix< T > &      A,
                const TTruncAcc &  acc )
{
    using  value_t = T;
    using  real_t  = typename real_type< value_t >::type_t;

    DO_CHECK_INF_NAN( A, "(BLAS) pseudo_invert", "in input matrix" );

    const auto         n = idx_t(A.nrows());
    const auto         m = idx_t(A.ncols());
    const auto         min_nm = std::min( n, m );
    Matrix< value_t >  U( n, m );
    Matrix< value_t >  V( m, min_nm );
    Vector< real_t >   S( min_nm );
            
    copy( A, U );
    svd( U, S, V );

    const auto  k = idx_t(acc.trunc_rank( S ));

    for ( idx_t  i = 0; i < k; ++i )
    {
        if ( S(i) != real_t(0) )
            S(i) = real_t(1) / S(i);
    }// for
    
    prod_diag( V, S, k );

    Matrix< value_t >  Uk( U, Range( 0, n-1 ), Range( 0, k-1 ) );
    Matrix< value_t >  Vk( V, Range( 0, m-1 ), Range( 0, k-1 ) );

    if ( n != m )
        A = std::move( Matrix< value_t >( m, n ) );
    
    prod( value_t(1), Vk, adjoint(Uk), value_t(0), A );

    DO_CHECK_INF_NAN( A, "(BLAS) pseudo_invert", "in output matrix" );
}

//!
//! compute LU factorisation of the n×m matrix \a A with
//! n×min(n,m) matrix L and min(n,m)xm matrix U; \a A will be
//! overwritten with L upon exit
//!
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
lu     ( T1 &  A )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    DO_CHECK_SMALL(   A );
    DO_CHECK_INF_NAN( A, "(BLAS) lu", "in input matrix" );
        
    const idx_t   nm1 = idx_t( A.nrows() )-1;
    const idx_t   mm1 = idx_t( A.ncols() )-1;
    const idx_t   mrc = std::min( nm1, mm1 );
    const real_t  eps = Math::square( Math::square( Limits::epsilon< real_t >() ) );
    
    for ( idx_t j = 0; j <= mrc; ++j )
    {
        // test for singularity
        if ( Math::abs( A(j,j) ) < eps )
            HERROR( ERR_DIV_ZERO, "(BLAS) lu", to_string( "a_(%d,%d) = 0", j, j ) );
        
        // compute elements j+1:n of j-th column
        if ( j < nm1 )
        {
            Vector< value_t >  col( A, Range( j+1, nm1 ), j );

            BLAS::scale( value_t(1.0) / A(j,j), col );
        }// if
        
        // update trailing submatrix
        if ( j < mrc )
        {
            Vector< value_t >  col( A, Range( j+1, nm1 ), j );
            Vector< value_t >  row( A, j, Range( j+1, mm1 ) );
            Matrix< value_t >  sub( A, Range( j+1, nm1 ), Range( j+1, mm1 ) );
            
            add_r1u( value_t(-1.0), col, row, sub );
        }// if
    }// for

    DO_CHECK_SMALL(   A );
    DO_CHECK_INF_NAN( A, "(BLAS) lu", "in output matrix" );
}

//!
//! compute L·L^T factorisation of given symmetric n×n matrix \a A
//!
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
llt    ( T1 &  A )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    if ( A.nrows() != A.ncols() )
        HERROR( ERR_MAT_SIZE, "(BLAS) llt", "matrix is not square" );
    
    DO_CHECK_SMALL(   A );
    DO_CHECK_INF_NAN( A, "(BLAS) llt", "in input matrix" );

    const idx_t  n = idx_t( A.nrows() );
    
    for ( idx_t  j = 0; j < n; ++j )
    {
        //
        // compute L_jj and test if <= 0
        //

        Vector< value_t >  rowj( A, j, Range( 0, j-1 ) );
        value_t                 ajj = A(j,j) - value_t( dotu( rowj, rowj ) );
            
        if ( std::real( ajj ) <= real_t(0) )
            HERROR( ERR_MAT_NPOSDEF, "(BLAS) llt", to_string( "A_%d,%d = %.8e", j, j, std::real(ajj) ) );
            
        ajj = Math::sqrt( ajj );
            
        A(j,j) = ajj;

        //
        // compute elements j+1:n of column j
        //
            
        if ( j < n-1 )
        {
            Vector< value_t >  colj( A, Range( j+1, n-1 ), j );
            
            // apply previous rows
            if ( j > 0 )
            {
                Matrix< value_t >  Asub( A, Range( j+1, n-1 ), Range( 0, j-1 ) );

                mulvec( value_t(-1), Asub, rowj, value_t(1), colj );
            }// if
            
            scale( value_t(1) / ajj, colj );
        }// if
    }// for
        
    // reset strictly upper part of L
    for ( idx_t  i = 1; i < n; ++i )
    {
        Vector< value_t >  Ai( A, Range( 0, i-1 ), i );
        
        fill( value_t(0), Ai );
    }// for
        
    DO_CHECK_SMALL(   A );
    DO_CHECK_INF_NAN( A, "(BLAS) llt", "in output matrix" );
}

//!
//! compute L·L^H factorisation of given hermitian n×n matrix \a A
//!
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
llh    ( T1 &  A )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    if ( A.nrows() != A.ncols() )
        HERROR( ERR_MAT_SIZE, "(BLAS) llh", "matrix is not square" );
    
    DO_CHECK_SMALL(   A );
    DO_CHECK_INF_NAN( A, "(BLAS) llh", "in input matrix" );

    const real_t  eps = Math::square( Math::square( Limits::epsilon< real_t >() ) );
    const idx_t   n   = idx_t( A.nrows() );
    
    for ( idx_t  j = 0; j < n; ++j )
    {
        if ( Math::abs( std::imag( A(j,j) ) ) > eps )
            HERROR( ERR_MAT_NHERM, "(BLAS) llh", to_string( "A_%d,%d = ", j, j ) + to_string( A(j,j) ) );
        
        //
        // compute L_jj and test if <= 0
        //

        Vector< value_t >  Aj( A, j, Range( 0, j-1 ) );
        value_t                 ajj = A(j,j) - value_t( dot( Aj, Aj ) );
            
        if ( std::real( ajj ) <= 0.0 )
        {
            A(j,j) = ajj;
            HERROR( ERR_MAT_NPOSDEF, "(BLAS) llh", to_string( "A_%d,%d = %.8e", j, j, std::real(ajj) ) );
        }// if
            
        ajj = Math::sqrt( ajj );
            
        A(j,j) = ajj;

        //
        // compute elements j+1:n of column j
        //
            
        if ( j < n-1 )
        {
            Vector< value_t >  colj( A, Range( j+1, n-1 ), j );

            // apply previous rows
            if ( j > 0 )
            {
                Matrix< value_t >  Asub( A, Range( j+1, n-1 ), Range( 0, j-1 ) );
                Vector< value_t >  rowj( A, j, Range( 0, j-1 ) );

                conj( rowj );
                mulvec( value_t(-1), Asub, rowj, value_t(1), colj );
                conj( rowj );
            }// if
            
            scale( value_t(1) / ajj, Aj );
        }// if
    }// for
        
    // reset strictly upper part of L
    for ( idx_t  i = 1; i < n; ++i )
    {
        Vector< value_t >  Ai( A, Range( 0, i-1 ), i );
        
        fill( value_t(0), Ai );
    }// for
        
    DO_CHECK_SMALL(   A );
    DO_CHECK_INF_NAN( A, "(BLAS) llh", "in output matrix" );
}

//!
//! compute L·D·L^T factorisation of given symmetric n×n matrix \a A
//!
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
ldlt   ( T1 &  A )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    if ( A.nrows() != A.ncols() )
        HERROR( ERR_MAT_SIZE, "(BLAS) ldlt", "matrix is not square" );

    DO_CHECK_SMALL(   A );
    DO_CHECK_INF_NAN( A, "(BLAS) ldlt", "in input matrix" );

    const idx_t   n   = idx_t( A.nrows() );
    const real_t  eps = Math::square( Math::square( Limits::epsilon< real_t >() ) );

    for ( idx_t  j = 0; j < n; ++j )
    {
        const value_t  d_j = A( j, j );
            
        // test for singularity
        if ( Math::abs( d_j ) < eps )
            HERROR( ERR_DIV_ZERO, "(BLAS) ldlt", to_string( "a_(%d,%d) = 0", j, j ) );
        
        // compute L(j+1:n,j)
        if ( j < n-1 )
        {
            Vector< value_t >  colj( A, Range( j+1, n-1 ), j );
            
            scale( value_t(1) / d_j, colj );
        }// if
        
        // update A(j+1:n,j+1:n)
        if ( j < n-1 )
        {
            for ( idx_t  jj = j+1; jj < n; ++jj )
            {
                // for ( idx_t  ii = jj; ii < n; ++ii )
                //     A( ii, jj ) -= A( ii, j ) * d_j * A( jj, j );
                const value_t           f = d_j * A( jj, j );
                Vector< value_t >  col_j(  A, Range( jj, n-1 ), j );
                Vector< value_t >  col_jj( A, Range( jj, n-1 ), jj );

                add( -f, col_j, col_jj );
            }// for
        }// if
    }// for
    
    // reset strictly upper part of L
    for ( idx_t  i = 1; i < n; ++i )
    {
        Vector< value_t >  Ai( A, Range( 0, i-1 ), i );
        
        fill( value_t(0), Ai );
    }// for
        
    DO_CHECK_SMALL(   A );
    DO_CHECK_INF_NAN( A, "(BLAS) ldlt", "in output matrix" );
}

//!
//! compute L·D·L^H factorisation of given hermitian n×n matrix \a A
//!
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
ldlh   ( T1 &  A )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    if ( A.nrows() != A.ncols() )
        HERROR( ERR_MAT_SIZE, "(BLAS) ldlh", "matrix is not square" );

    DO_CHECK_SMALL(   A );
    DO_CHECK_INF_NAN( A, "(BLAS) ldlh", "in input matrix" );

    const idx_t   n   = idx_t( A.nrows() );
    const real_t  eps = Math::square( Math::square( Limits::epsilon< real_t >() ) );

    for ( idx_t  j = 0; j < n; ++j )
    {
        const value_t  d_j = A( j, j );
            
        // test for singularity
        if ( Math::abs( d_j ) < eps )
            HERROR( ERR_DIV_ZERO, "(BLAS) ldlh", to_string( "a_(%d,%d) = 0", j, j ) );
        
        // compute L(j+1:n,j)
        if ( j < n-1 )
        {
            Vector< value_t >  colj( A, Range( j+1, n-1 ), j );
            
            scale( value_t(1) / d_j, colj );
        }// if
        
        // update A(j+1:n,j+1:n)
        if ( j < n-1 )
        {
            for ( idx_t  jj = j+1; jj < n; ++jj )
            {
                // for ( idx_t  ii = jj; ii < n; ++ii )
                //     A( ii, jj ) -= A( ii, j ) * d_j * Hpro::conj( A( jj, j ) );
                
                const value_t           f = d_j * Math::conj( A( jj, j ) );
                Vector< value_t >  col_j(  A, Range( jj, n-1 ), j );
                Vector< value_t >  col_jj( A, Range( jj, n-1 ), jj );

                add( -f, col_j, col_jj );
            }// for
        }// if
    }// for
    
    // reset strictly upper part of L
    for ( idx_t  i = 1; i < n; ++i )
    {
        Vector< value_t >  Ai( A, Range( 0, i-1 ), i );
        
        fill( value_t(0), Ai );
    }// for
        
    DO_CHECK_SMALL(   A );
    DO_CHECK_INF_NAN( A, "(BLAS) ldlh", "in output matrix" );
}

namespace
{

// template version of GPU::*_qr_size
template < typename value_t >  size_t  gpu_qr_size ();
template <>                    size_t  gpu_qr_size< float >                  () { return CFG::GPU::fp32_qr_size; }
template <>                    size_t  gpu_qr_size< double >                 () { return CFG::GPU::fp64_qr_size; }
template <>                    size_t  gpu_qr_size< std::complex< float > >  () { return CFG::GPU::fp32_qr_size; }
template <>                    size_t  gpu_qr_size< std::complex< double > > () { return CFG::GPU::fp64_qr_size; }

}// namespace anonymous

//!
//! compute LQ factorisation of the n×m matrix \a A with
//! m×m matrix Q and nxm matrix L (n >= m); \a A will be
//! overwritten with Q upon exit
//!
template < typename value_t >
void
lq ( Matrix< value_t > &  A,
     Matrix< value_t > &  L )
{
    LOG_COUNTER( qr );
    LOG_TIC( qr );

    if ( false && ( std::min( A.nrows(), A.ncols() ) <= 32 ))
    {
        MKL_SEQ_START;

        const auto               nrows = A.nrows();
        const auto               ncols = A.ncols();
        blas_int_t               info  = 0;
        Hpro::vector< value_t >  tau( ncols );
        Hpro::vector< value_t >  work( ncols );
    
        gelq2( nrows, ncols, A.data(), A.col_stride(), tau.data(), work.data(), info );

        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) lq", to_string( "argument %d to LAPACK::gelq2", -info ) );
    
        if (( L.nrows() != nrows ) || ( L.ncols() != ncols ))
            L = std::move( Matrix< value_t >( nrows, ncols ) );
    
        for ( size_t  i = 0; i < nrows; ++i )
            for ( size_t  j = 0; j <= std::min( i, ncols-1 ); ++j )
                L(i,j) = A(i,j);

        org2l( nrows, ncols, ncols, A.data(), A.col_stride(), tau.data(), work.data(), info );

        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) lq", to_string( "argument %d to LAPACK::org2l", -info ) );
    
        MKL_SEQ_END;
    }// if
    else
    {
        const auto  nrows = blas_int_t( A.nrows() );
        const auto  ncols = blas_int_t( A.ncols() );
        blas_int_t  info  = 0;

        //
        // workspace query
        //

        value_t  dummy      = value_t(0); // non-NULL workspace for latest Intel MKL
        value_t  work_query = value_t(0);

        gelqf( nrows, ncols, A.data(), blas_int_t( A.col_stride() ), & dummy, & work_query, LAPACK_WS_QUERY, info );

        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) lq", to_string( "argument %d to LAPACK::gelqf", -info ) );
    
        auto   lwork = blas_int_t( std::real( work_query ) );

        orglq( nrows, ncols, nrows, A.data(), blas_int_t( A.col_stride() ), & dummy, & work_query, LAPACK_WS_QUERY, info );
    
        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) lq", to_string( "argument %d to LAPACK::orglq", -info ) );

        // adjust work space size
        lwork = std::max( lwork, blas_int_t( std::real( work_query ) ) );
    
        Hpro::vector< value_t >  tmp_space( lwork + ncols );
        value_t *                work = tmp_space.data();
        value_t *                tau  = work + lwork;

        //
        // compute Householder vectors and R
        //

        MKL_SEQ_START;

        gelqf( nrows, ncols, A.data(), blas_int_t( A.col_stride() ), tau, work, lwork, info );
    
        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) lq", to_string( "argument %d to LAPACK::gelqf", -info ) );

        //
        // copy lower triangular matrix to L
        //

        if (( blas_int_t( L.nrows() ) != nrows ) || ( blas_int_t( L.ncols() ) != ncols ))
            L = std::move( Matrix< value_t >( nrows, ncols ) );
        else
            fill( value_t(0), L );
    
        for ( blas_int_t  i = 0; i < nrows; i++ )
        {
            Vector< value_t >  colA( A, i, Range( 0, std::min( i, ncols-1 ) ) );
            Vector< value_t >  colL( L, i, Range( 0, std::min( i, ncols-1 ) ) );

            copy( colA, colL );
        }// for

        //
        // compute Q
        //
    
        orglq( nrows, ncols, nrows, A.data(), blas_int_t( A.col_stride() ), tau, work, lwork, info );

        MKL_SEQ_END;

        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) lq", to_string( "argument %d to LAPACK::orglq", -info ) );
    }
    
    LOG_TOC( qr );
}

//!
//! compute QR factorisation of the n×m matrix \a A with
//! n×m matrix Q and mxm matrix R (n >= m); \a A will be
//! overwritten with Q upon exit
//!
template < typename value_t >
void
qr ( Matrix< value_t > &  A,
     Matrix< value_t > &  R )
{
    LOG_COUNTER( qr );
    LOG_TIC( qr );

    //
    // use LQ on A' if nrows(A) < ncols(A)
    //
    
    if ( A.nrows() < A.ncols() )
    {
        auto  AT = Matrix< value_t >( A.ncols(), A.nrows() );
        auto  L  = Matrix< value_t >();

        // DBG::write( A, "A.mat", "A" );
        
        copy( adjoint(A), AT );
        // DBG::write( AT, "AT.mat", "AT" );

        lq( AT, L );

        // DBG::write( AT, "Q.mat", "Q" );
        // DBG::write( L, "L.mat", "L" );

        if (( A.nrows() != AT.ncols() ) || ( A.ncols() != AT.nrows() ))
            A = std::move( Matrix< value_t >( AT.ncols(), AT.nrows() ) );
                           
        copy( adjoint(AT), A );
        
        if (( R.nrows() != L.ncols() ) || ( R.ncols() != L.nrows() ))
            R = std::move( Matrix< value_t >( L.ncols(), L.nrows() ) );
                           
        copy( adjoint(L), R );
        
        return;
    }// if
    
    //
    // try GPU if large enough
    //
    
    bool  use_cpu = true;

    if ( std::min( A.nrows(), A.ncols() ) >= gpu_qr_size< value_t >() )
    {
        LOG_COUNTER( gpu );
        LOG_TIC( gpu );
        
        const auto  handle = GPU::request_handle();

        if ( handle != GPU::INVALID_HANDLE )
        {
            if ( GPU::qr( handle, A, R ) )
                use_cpu = false;

        //     const auto  nrows = A.nrows();
        //     const auto  ncols = A.ncols();
        //     const auto  minrc = std::min( nrows, ncols );
        //     auto        tau   = BLAS::Vector< value_t >( minrc );
                
        //     if ( GPU::qr( handle, A, tau ) )
        //     {
        //         //
        //         // copy upper triangular matrix to R
        //         //

        //         if (( blas_int_t( R.nrows() ) != ncols ) || ( blas_int_t( R.ncols() ) != ncols ))
        //             R = std::move( Matrix< value_t >( ncols, ncols ) );
        //         else
        //             fill( value_t(0), R );
    
        //         for ( blas_int_t  i = 0; i < ncols; i++ )
        //         {
        //             Vector< value_t >  colA( A, Range( 0, i ), i );
        //             Vector< value_t >  colR( R, Range( 0, i ), i );

        //             copy( colA, colR );
        //         }// for

        //         //
        //         // compute Q
        //         //

        //         value_t  dummy      = value_t(0); // non-NULL workspace for Intel MKL
        //         value_t  work_query = value_t(0);
        //         auto     info       = blas_int_t(0);

        //         orgqr( nrows, ncols, ncols, A.data(), blas_int_t( A.col_stride() ), & dummy, & work_query, LAPACK_WS_QUERY, info );
    
        //         if ( info < 0 )
        //             HERROR( ERR_ARG, "(BLAS) qr", to_string( "argument %d to LAPACK::orgqr", -info ) );

        //         // adjust work space size
        //         auto  lwork = blas_int_t( std::real( work_query ) );
        //         auto  work  = Hpro::vector< value_t >( lwork );

        //         MKL_SEQ_START;
                
        //         orgqr( nrows, ncols, ncols, A.data(), blas_int_t( A.col_stride() ), tau.data(), work.data(), lwork, info );

        //         MKL_SEQ_END;
                
        //         if ( info < 0 )
        //             HERROR( ERR_ARG, "(BLAS) qr", to_string( "argument %d to LAPACK::orgqr", -info ) );

        //         // do not fall back to CPU based computation
        //         GPU::release_handle( handle );
        //         LOG_TOC( gpu );
        //         LOG_TOC( qr );
        //         return;
        //     }// if
        }// if

        GPU::release_handle( handle );
        LOG_TOC( gpu );
    }// if

    //
    // else use CPU
    //
    
    if ( std::min( A.nrows(), A.ncols() ) <= 32 )
    {

        MKL_SEQ_START;

        const auto               nrows = A.nrows();
        const auto               ncols = A.ncols();
        const auto               minrc = std::min( nrows, ncols );
        blas_int_t               info  = 0;
        Hpro::vector< value_t >  tau( ncols );
        Hpro::vector< value_t >  work( ncols );
    
        geqr2( nrows, ncols, A.data(), nrows, tau.data(), work.data(), info );

        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) qr", to_string( "argument %d to LAPACK::geqr2", -info ) );
    
        if (( R.nrows() != minrc ) || ( R.ncols() != ncols ))
            R = std::move( Matrix< value_t >( minrc, ncols ) );
    
        for ( size_t  i = 0; i < ncols; ++i )
            for ( size_t  j = 0; j <= std::min( i, minrc ); ++j )
                R(j,i) = A(j,i);

        org2r( nrows, ncols, ncols, A.data(), nrows, tau.data(), work.data(), info );

        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) qr", to_string( "argument %d to LAPACK::org2r", -info ) );
    
        MKL_SEQ_END;
    
    }// if
    else
    {
        const blas_int_t  n     = blas_int_t( A.nrows() );
        const blas_int_t  m     = blas_int_t( A.ncols() );
        blas_int_t        info  = 0;

        //
        // workspace query
        //

        value_t  dummy      = value_t(0); // non-NULL workspace for latest Intel MKL
        value_t  work_query = value_t(0);

        geqrf( n, m, A.data(), blas_int_t( A.col_stride() ), & dummy, & work_query, LAPACK_WS_QUERY, info );

        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) qr", to_string( "argument %d to LAPACK::geqrf", -info ) );
    
        blas_int_t   lwork = blas_int_t( std::real( work_query ) );

        orgqr( n, m, m, A.data(), blas_int_t( A.col_stride() ), & dummy, & work_query, LAPACK_WS_QUERY, info );
    
        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) qr", to_string( "argument %d to LAPACK::orgqr", -info ) );

        // adjust work space size
        lwork = std::max( lwork, blas_int_t( std::real( work_query ) ) );
    
        Hpro::vector< value_t >  tmp_space( lwork + m );
        value_t *                work  = tmp_space.data();
        value_t *                tau   = work + lwork;

        //
        // compute Householder vectors and R
        //

        MKL_SEQ_START;

        geqrf( n, m, A.data(), blas_int_t( A.col_stride() ), tau, work, lwork, info );
    
        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) qr", to_string( "argument %d to LAPACK::geqrf", -info ) );

        //
        // copy upper triangular matrix to R
        //

        if (( blas_int_t( R.nrows() ) != m ) || ( blas_int_t( R.ncols() ) != m ))
            R = std::move( Matrix< value_t >( m, m ) );
        else
            fill( value_t(0), R );
    
        for ( blas_int_t  i = 0; i < m; i++ )
        {
            Vector< value_t >  colA( A, Range( 0, i ), i );
            Vector< value_t >  colR( R, Range( 0, i ), i );

            copy( colA, colR );
        }// for

        //
        // compute Q
        //
    
        orgqr( n, m, m, A.data(), blas_int_t( A.col_stride() ), tau, work, lwork, info );

        MKL_SEQ_END;

        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) qr", to_string( "argument %d to LAPACK::orgqr", -info ) );

    }
    
    LOG_TOC( qr );
}

//
// compute QR factorisation of the tall-and-skinny n×m matrix \a A, m ≪ n,
// with n×m matrix Q and mxm matrix R (n >= m); \a A will be overwritten
// with Q upon exit
//
template < typename value_t >
void
tsqr  ( Matrix< value_t > &  A,
        Matrix< value_t > &  R,
        const size_t         ntile2 )
{
    const size_t  nrows   = A.nrows();
    const size_t  ncols   = A.ncols();
    const bool    use_tbb = false;
    const size_t  ntile   = std::max( ntile2 > 0
                                      ? ntile2                // use user provided value
                                      : (use_tbb
                                         ? ( ncols >= 16
                                             ? nrows / 4      // two sub-divides 
                                             : ( ncols >= 8
                                                 ? nrows / 2  // one sub-divide
                                                 : nrows ) )
                                         : nrows ),
                                      ncols );                // ensure TS property

    if ( ncols > nrows )
        HERROR( ERR_ARG, "(BLAS) tsqr", "#cols > #rows" );
    
    if ( nrows > ntile )
    {
        auto  mid   = nrows / 2;
        auto  rows0 = Range( 0, mid-1 );
        auto  rows1 = Range( mid, nrows-1 );
        auto  Q0    = Matrix< value_t >( A, rows0, Range::all, copy_value );
        auto  Q1    = Matrix< value_t >( A, rows1, Range::all, copy_value );
        auto  R0    = Matrix< value_t >( ncols, ncols );
        auto  R1    = Matrix< value_t >( ncols, ncols );

        //
        // A = | Q0 R0 | = | Q0   | | R0 | = | Q0   | Q2 R
        //     | Q1 R1 |   |   Q1 | | R1 |   |   Q1 | 
        //
        
        tsqr( Q0, R0, ntile );
        tsqr( Q1, R1, ntile );

        auto  Q2  = Matrix< value_t >( 2*ncols, ncols );
        auto  Q20 = Matrix< value_t >( Q2, Range( 0,     ncols-1   ), Range::all );
        auto  Q21 = Matrix< value_t >( Q2, Range( ncols, 2*ncols-1 ), Range::all );

        copy( R0, Q20 );
        copy( R1, Q21 );

        qr( Q2, R );

        //
        // Q = | Q0    | Q    (overwrite A)
        //     |    Q1 |
        //
        
        auto  Q_0  = Matrix< value_t >( A, rows0, Range::all );
        auto  Q_1  = Matrix< value_t >( A, rows1, Range::all );

        prod( value_t(1), Q0, Q20, value_t(0), Q_0 );
        prod( value_t(1), Q1, Q21, value_t(0), Q_1 );
    }// if
    else
    {
        qr( A, R );
    }// else
}

//
// Compute QR factorisation with column pivoting of the matrix \a A.
//
template < typename value_t >
void
qrp ( Matrix< value_t > &          A,
      Matrix< value_t > &          R,
      std::vector< blas_int_t > &  P )
{
    LOG_COUNTER( qr );
    LOG_TIC( qr );
    
    using  real_t  = typename real_type< value_t >::type_t;
    
    const blas_int_t  n     = blas_int_t( A.nrows() );
    const blas_int_t  m     = blas_int_t( A.ncols() );
    const blas_int_t  minnm = std::min( n, m );
    blas_int_t        info  = 0;

    if ( blas_int_t( P.size() ) != m )
        P.resize( m );
    
    //
    // workspace query
    //

    value_t                 dummy      = value_t(0); // non-NULL workspace for latest Intel MKL
    value_t                 work_query = value_t(0);
    Hpro::vector< real_t >  rwork( 2*m );

    geqp3( n, m, A.data(), blas_int_t( A.col_stride() ), P.data(),
           & dummy, & work_query, LAPACK_WS_QUERY, rwork.data(), info );

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) qrp", to_string( "argument %d to LAPACK::geqp3", -info ) );
    
    blas_int_t   lwork = blas_int_t( std::real( work_query ) );

    orgqr( n, minnm, minnm, A.data(), blas_int_t( A.col_stride() ), & dummy,
           & work_query, LAPACK_WS_QUERY, info );
    
    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) qrp", to_string( "argument %d to LAPACK::orgqr", -info ) );

    // adjust work space size
    lwork = std::max( lwork, blas_int_t( std::real( work_query ) ) );
    
    Hpro::vector< value_t >  tmp_space( lwork + m );
    value_t *                work  = tmp_space.data();
    value_t *                tau   = work + lwork;

    //
    // compute Householder vectors and R
    //
    
    MKL_SEQ_START;
    
    geqp3( n, m, A.data(), blas_int_t( A.col_stride() ), P.data(), tau, work, lwork, rwork.data(), info );
    
    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) qrp", to_string( "argument %d to LAPACK::geqp3", -info ) );

    //
    // copy upper triangular matrix to R
    //

    if (( blas_int_t( R.nrows() ) != m ) || ( blas_int_t( R.ncols() ) != m ))
        R = std::move( Matrix< value_t >( m, m ) );
    else
        fill( value_t(0), R );
    
    for ( blas_int_t  i = 0; i < m; i++ )
    {
        const auto         irange = Range( 0, std::min( i, minnm-1 ) );
        Vector< value_t >  colA( A, irange, i );
        Vector< value_t >  colR( R, irange, i );

        copy( colA, colR );
    }// for
        
    //
    // compute Q
    //
    
    orgqr( n, minnm, minnm, A.data(), blas_int_t( A.col_stride() ), tau, work, lwork, info );

    if ( n < m )
    {
        //
        // realloc Q
        //

        auto  Q  = Matrix< value_t >( n, n );
        auto  Ai = Matrix< value_t >( A, Range::all, Range( 0, n-1 ) );

        copy( Ai, Q );

        A = std::move( Q );
    }// if
    
    MKL_SEQ_END;

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) qr", to_string( "argument %d to LAPACK::orgqr", -info ) );

    //
    // adjount indices in P (1-counted to 0-counted)
    //

    for ( blas_int_t  i = 0; i < m; i++ )
        --P[i];
    
    LOG_TOC( qr );
}

//
// Compute QR factorisation with column pivoting of the matrix \a A.
//
template < typename T >
std::enable_if_t< is_matrix< T >::value, idx_t >
qrp_trunc ( T &                              A,
            Matrix< typename T::value_t > &  R,
            std::vector< blas_int_t > &      P,
            const TTruncAcc &                acc )
{
    using  value_t = typename T::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    #if HPRO_HAS_GEQP3_TRUNC == 1

    //
    // do QRP and determine ranks on-the-fly, stopping if accuracy was reached
    //
    
    LOG_COUNTER( qr );
    LOG_TIC( qr );
    
    const blas_int_t  n     = blas_int_t( A.nrows() );
    const blas_int_t  m     = blas_int_t( A.ncols() );
    blas_int_t        info  = 0;

    if ( blas_int_t( P.size() ) != m )
        P.resize( m );
    
    //
    // workspace query
    //

    blas_int_t              ntrunc     = ( acc.is_fixed_rank() ? acc.rank() : std::min( n, m ) );
    real_t                  atrunc     = acc.abs_eps();
    real_t                  rtrunc     = ( acc.is_fixed_rank() ? real_t(0)  : acc.rel_eps() );
    value_t                 dummy      = value_t(0); // non-NULL workspace for latest Intel MKL
    value_t                 work_query = value_t(0);
    Hpro::vector< real_t >  rwork( 2*m );

    geqp3trunc( n, m, A.data(), blas_int_t( A.col_stride() ), P.data(), & dummy,
                ntrunc, atrunc, rtrunc,
                & work_query, LAPACK_WS_QUERY, rwork.data(), info );

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) qrp", to_string( "argument %d to LAPACK::geqp3", -info ) );
    
    blas_int_t   lwork = blas_int_t( std::real( work_query ) );

    orgqr( n, m, m, A.data(), blas_int_t( A.col_stride() ), & dummy,
           & work_query, LAPACK_WS_QUERY, info );
    
    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) qrp", to_string( "argument %d to LAPACK::orgqr", -info ) );

    // adjust work space size
    lwork = std::max( lwork, blas_int_t( std::real( work_query ) ) );
    
    Hpro::vector< value_t >  tmp_space( lwork + m );
    value_t *                work  = tmp_space.data();
    value_t *                tau   = work + lwork;

    //
    // compute Householder vectors and R
    //

    MKL_SEQ_START;

    geqp3trunc( n, m, A.data(), blas_int_t( A.col_stride() ), P.data(), tau,
                ntrunc, atrunc, rtrunc,
                work, lwork, rwork.data(), info );
    
    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) qrp", to_string( "argument %d to LAPACK::geqp3", -info ) );

    // adjust number of columns in A
    auto  m_trunc = ntrunc;
    
    //
    // copy upper triangular matrix to R
    //

    if (( blas_int_t( R.nrows() ) != m_trunc ) || ( blas_int_t( R.ncols() ) != m ))
        R = std::move( Matrix< value_t >( m_trunc, m ) );
    else
        fill( value_t(0), R );
    
    for ( blas_int_t  i = 0; i < m_trunc; i++ )
    {
        Vector< value_t >  colA( A, Range( 0, i ), i );
        Vector< value_t >  colR( R, Range( 0, i ), i );

        copy( colA, colR );
    }// for

    if ( m > m_trunc )
    {
        Matrix< value_t >  A_rest( A, Range( 0, m_trunc-1 ), Range( m_trunc, m-1 ) );
        Matrix< value_t >  R_rest( R, Range( 0, m_trunc-1 ), Range( m_trunc, m-1 ) );

        copy( A_rest, R_rest );
    }// if

    //
    // compute Q
    //
    
    orgqr( n, m, m_trunc, A.data(), blas_int_t( A.col_stride() ), tau, work, lwork, info );

    MKL_SEQ_END;

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) qr", to_string( "argument %d to LAPACK::orgqr", -info ) );

    //
    // adjust indices in P (1-counted to 0-counted)
    //

    for ( blas_int_t  i = 0; i < m; i++ )
        --P[i];
    
    LOG_TOC( qr );

    DO_CHECK_INF_NAN( A, "(BLAS) qrp_trunc", "in output matrix A" );
    DO_CHECK_INF_NAN( R, "(BLAS) qrp_trunc", "in output matrix R" );

    return R.nrows();
    
    #else

    //
    // do full QRP and determine rank later
    //
    
    const idx_t  k_old = idx_t( A.ncols() );
        
    qrp( A, R, P );
            
    //
    // determine norms of submatrices R(0:i,0:i) and R(i+1:k,i+1:k)
    //
    
    Vector< real_t >  S( k_old );
    
    for ( int i = 0; i < k_old; ++i )
    {
        Matrix< value_t >  R2( R, Range( i, k_old-1 ), Range( i, k_old-1 ) );
        
        S( i ) = normF( R2 );
    }// for
    
    // determine truncated rank based on norms
    return idx_t( acc.trunc_rank( S ) );
    
    #endif
}

//
// compute eigenvalues and eigenvectors of matrix \a M
//
template < typename value_t >
void
eigen ( Matrix< value_t > &  M,
        Vector< value_t > &  eig_val,
        Matrix< value_t > &  eig_vec )
{
    using  real_t  = typename real_type< value_t >::type_t;
    
    //
    // first check, if M is symmetric (or hermitian
    //

    bool  is_herm = false;

    const size_t  n = M.nrows();
    const size_t  m = M.ncols();

    if ( n == m )
    {
        is_herm = true;

        for ( idx_t  j = 1; j < idx_t(m); ++j )
        {
            for ( idx_t  i = 0; i < j; ++i )
            {
                if ( M(i,j) != Math::conj( M(j,i) ) ) is_herm = false;
            }// for

            if ( ! is_herm )
                break;
        }// for
    }// if

    blas_int_t  info = 0;
    
    if ( is_herm )
    {
        value_t  work_query = value_t(0);
        real_t   rdummy     = 0;

        if (( eig_vec.nrows() != M.nrows() ) || ( eig_vec.ncols() != M.ncols() ))
            eig_vec= std::move( Matrix< value_t >( M.nrows(), M.ncols() ) );
        
        // work space query
        heev( 'V', 'L', blas_int_t(n), eig_vec.data(), blas_int_t(eig_vec.col_stride()),
              & rdummy, & work_query, LAPACK_WS_QUERY, & rdummy, info );

        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) eigen", to_string( "argument %d to LAPACK::*(sy|he)ev", -info ) );
        
        const blas_int_t        lwork = blas_int_t( std::real( work_query ) );
        Vector< real_t >        seig_val( n );
        vector< value_t >       work( lwork );
        Hpro::vector< real_t >  rwork( is_complex_type< value_t >::value ? 3*n-2 : 0 );

        copy( M, eig_vec );

        MKL_SEQ_START;
        
        heev( 'V', 'L', blas_int_t(n), eig_vec.data(), blas_int_t(eig_vec.col_stride()),
              seig_val.data(), work.data(), lwork, rwork.data(), info );

        MKL_SEQ_END;
        
        if ( eig_val.length() != n )
            eig_val = std::move( Vector< value_t >( n ) );
        
        for ( idx_t  i = 0; i < idx_t(n); ++i )
            eig_val(i) = seig_val(i);
        
        if ( info < 0 )
            HERROR( ERR_ARG, "(BLAS) eigen", to_string( "argument %d to LAPACK::*(sy|he)ev", -info ) );
    }// if
    else // unsymmetric
    {
        HERROR( ERR_NOT_IMPL, "(BLAS) eigen", "unsymmetric case" );
    }// else
}

//
// compute eigenvalues and eigenvectors of the hermitian matrix \a M
//
template < typename value_t >
void
eigen_herm ( Matrix< value_t > &                 M,
             Vector< real_type_t< value_t > > &  eig_val,
             Matrix< value_t > &                 eig_vec )
{
    using  real_t  = typename real_type< value_t >::type_t;

    #if 0

    const size_t  n          = M.nrows();
    blas_int_t    info       = 0;
    value_t       work_query = value_t(0);
    real_t        rdummy     = 0;

    if (( eig_vec.nrows() != M.nrows() ) || ( eig_vec.ncols() != M.ncols() ))
        eig_vec= std::move( Matrix< value_t >( M.nrows(), M.ncols() ) );
        
    // work space query
    heev( 'V', 'L', blas_int_t(n), eig_vec.data(), blas_int_t(eig_vec.col_stride()),
          & rdummy, & work_query, LAPACK_WS_QUERY, & rdummy, info );

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) eigen_herm", to_string( "argument %d to LAPACK::*(sy|he)ev", -info ) );
        
    if ( eig_val.length() != n )
        eig_val = std::move( Vector< real_t >( n ) );
    
    const blas_int_t         lwork = blas_int_t( std::real( work_query ) );
    Hpro::vector< value_t >  work( lwork );
    Hpro::vector< real_t >   rwork( is_complex_type< value_t >::value ? 3*n-2 : 0 );
    
    copy( M, eig_vec );
    
    MKL_SEQ_START;
    
    heev( 'V', 'L', blas_int_t(n), eig_vec.data(), blas_int_t(eig_vec.col_stride()),
          eig_val.data(), work.data(), lwork, rwork.data(), info );
    
    MKL_SEQ_END;
    
    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) eigen_herm", to_string( "argument %d to LAPACK::*(sy|he)ev", -info ) );

    #else

    const size_t  n           = M.nrows();
    blas_int_t    info        = 0;
    value_t       work_query  = value_t(0);
    real_t        rwork_query = real_t(0);
    blas_int_t    iwork_query = 0;
    real_t        rdummy      = 0;

    if (( eig_vec.nrows() != M.nrows() ) || ( eig_vec.ncols() != M.ncols() ))
        eig_vec= std::move( Matrix< value_t >( M.nrows(), M.ncols() ) );
        
    // work space query
    heevd( 'V', 'L', blas_int_t(n), eig_vec.data(), blas_int_t(eig_vec.col_stride()),
           & rdummy,
           & work_query,  LAPACK_WS_QUERY,
           & rwork_query, LAPACK_WS_QUERY,
           & iwork_query, LAPACK_WS_QUERY,
           info );

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) eigen", to_string( "argument %d to LAPACK::*(sy|he)evd", -info ) );
        
    if ( eig_val.length() != n )
        eig_val = std::move( Vector< real_t >( n ) );
    
    const blas_int_t            lwork  = blas_int_t( std::real( work_query ) );
    const blas_int_t            lrwork = blas_int_t( rwork_query );
    const blas_int_t            liwork = iwork_query;
    Hpro::vector< value_t >     work( lwork );
    Hpro::vector< real_t >      rwork( lrwork );
    Hpro::vector< blas_int_t >  iwork( liwork );
    
    copy( M, eig_vec );
    
    MKL_SEQ_START;
    
    heevd( 'V', 'L', blas_int_t(n), eig_vec.data(), blas_int_t(eig_vec.col_stride()),
           eig_val.data(),
           work.data(),  lwork,
           rwork.data(), lrwork,
           iwork.data(), liwork,
           info );
    
    MKL_SEQ_END;
    
    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) eigen", to_string( "argument %d to LAPACK::*(sy|he)evd", -info ) );

    #endif
}

//
// compute eigenvalues and eigenvectors of matrix \a M
//
template < typename value_t >
void
eigen_herm ( Matrix< value_t > &  M,
             const Range &        eig_range,
             Vector< value_t > &  eig_val,
             Matrix< value_t > &  eig_vec )
{
    const size_t  n          = M.nrows();
    blas_int_t    m          = 0;
    blas_int_t    info       = 0;
    value_t       work_query = value_t(0);

    if (( eig_vec.nrows() != M.nrows() ) || ( eig_vec.ncols() != eig_range.size() ))
        eig_vec = std::move( Matrix< value_t >( M.nrows(), eig_range.size() ) );

    if ( eig_val.length() != eig_range.size() )
        eig_val = std::move( Vector< value_t >( eig_range.size() ) );

    Vector< value_t >  tmp_eig_val( M.ncols() );
    blas_int_t         idummy     = 0;
        
    // work space query
    heevx( 'V', 'L', blas_int_t(n), M.data(), blas_int_t(M.col_stride()),
           blas_int_t(eig_range.first())+1, blas_int_t(eig_range.last())+1,
           m, tmp_eig_val.data(), eig_vec.data(), blas_int_t(eig_vec.col_stride()),
           & work_query, LAPACK_WS_QUERY, & idummy, info );

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) eigen", to_string( "argument %d to LAPACK::*(sy|he)evx", -info ) );
        
    const blas_int_t            lwork = blas_int_t( std::real( work_query ) );
    Hpro::vector< value_t >     work( lwork );
    Hpro::vector< blas_int_t >  iwork( 5*n );
    
    MKL_SEQ_START;
        
    heevx( 'V', 'L', blas_int_t(n), M.data(), blas_int_t(M.col_stride()),
           blas_int_t(eig_range.first())+1, blas_int_t(eig_range.last())+1,
           m, tmp_eig_val.data(), eig_vec.data(), blas_int_t(eig_vec.col_stride()),
           work.data(), lwork, iwork.data(), info );

    MKL_SEQ_END;

    if ( m != blas_int_t(eig_range.size()) )
    {
        HERROR( ERR_NOT_IMPL, "", "" );
    }// if

    for ( blas_int_t  i = 0; i < m; ++i )
        eig_val(i) = tmp_eig_val(i);
    
    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) eigen", to_string( "argument %d to LAPACK::*(sy|he)evx", -info ) );
}

//!
//! compute eigenvalues and eigenvectors 
//! of the \b symmetric, \b tridiagonal matrix defines by diagonal
//! coefficients in \a diag and off-diagonal coefficients \a subdiag
//!
template < typename value_t >
void
eigen_tri ( Vector< value_t > &  diag,
            Vector< value_t > &  subdiag,
            Vector< value_t > &  eig_val,
            Matrix< value_t > &  eig_vec )
{
    const size_t             n = diag.length();
    Hpro::vector< value_t >  work( std::max<size_t>( 1, 2*n-2 ) );
    blas_int_t               info = 0;

    if ( subdiag.length() != n-1 )
        HERROR( ERR_VEC_SIZE, "(BLAS) eigen_tri", "sub diagonal has wrong dimension" );

    if (( eig_vec.nrows() != n ) || ( eig_vec.ncols() != n ))
        eig_vec = std::move( Matrix< value_t >( n, n ) );
    
    if ( eig_val.length() != n )
        eig_val = std::move( Vector< value_t >( n ) );

    copy( diag, eig_val );
    
    MKL_SEQ_START;
    
    stev( 'V', blas_int_t(n), eig_val.data(), subdiag.data(),
          eig_vec.data(), blas_int_t(eig_vec.col_stride()),
          work.data(), info );

    MKL_SEQ_END;

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) eigen_tri", to_string( "argument %d to LAPACK::*stev", -info ) );
}

//
// compute SVD of matrix A using LAPACK
//
template < typename T1,
           typename T2,
           typename T3 >
std::enable_if_t< is_matrix< T1 >::value &&
                  is_vector< T2 >::value &&
                  is_matrix< T3 >::value &&
                  std::is_same< typename T1::value_t, typename T3::value_t >::value &&
                  std::is_same< typename real_type< typename T1::value_t >::type_t,
                                typename T2::value_t >::value, void >
gesvd    ( T1 &   A,
           T2 &   S,
           T3 &   V )
{
    LOG_COUNTER( gesvd );
    LOG_TIC( gesvd );
    
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    const blas_int_t    n      = blas_int_t( A.nrows() );
    const blas_int_t    m      = blas_int_t( A.ncols() );
    const blas_int_t    min_nm = std::min( n, m );
    blas_int_t          info   = 0;
    value_t             work_query = value_t(0);
    Matrix< value_t >   VT( min_nm, A.ncols() );
    value_t             vdummy = 0;
    real_t              rdummy = 0;

    if ( S.length() != size_t( min_nm ) )
        S = std::move( Vector< real_t >( min_nm ) );
    
    // work space query
    gesvd( 'O', 'S',
           blas_int_t( A.nrows() ),
           blas_int_t( A.ncols() ),
           A.data(),
           blas_int_t( A.col_stride() ),
           S.data(),
           & vdummy,
           blas_int_t( A.col_stride() ),
           & vdummy,
           blas_int_t( VT.col_stride() ),
           & work_query,
           LAPACK_WS_QUERY,
           & rdummy,
           info );

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) gesvd", to_string( "argument %d to LAPACK::gesvd", -info ) );
        
    const blas_int_t         lwork = blas_int_t( std::real( work_query ) );
    Hpro::vector< value_t >  work( lwork );
    Hpro::vector< real_t >   rwork( 5 * min_nm );
    
    /////////////////////////////////////////////////////////////////
    #if CHECK_GESVD == 1
    Matrix< value_t >  TA( n, m );

    copy( A, TA );
    #endif
    /////////////////////////////////////////////////////////////////

    MKL_SEQ_START;

    gesvd( 'O', 'S',
           blas_int_t( A.nrows() ),
           blas_int_t( A.ncols() ),
           A.data(),
           blas_int_t( A.col_stride() ),
           S.data(),
           & vdummy,
           blas_int_t( A.col_stride() ),
           VT.data(),
           blas_int_t( VT.col_stride() ),
           work.data(),
           lwork,
           rwork.data(),
           info );

    MKL_SEQ_END;

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) gesvd", to_string( "argument %d to LAPACK::gesvd", -info ) );
    else if ( info > 0 )
        HERROR( ERR_NCONVERGED, "(BLAS) gesvd", to_string( "in LAPACK::gesvd (info = %d)", info ) );

    if (( V.nrows() != VT.ncols() ) || ( V.ncols() != size_t(min_nm) ))
        V = std::move( Matrix< value_t >( VT.ncols(), min_nm ) );
    
    copy( adjoint( VT ), V );

    /////////////////////////////////////////////////////////////////
    #if CHECK_GESVD == 1
    {
        Matrix< value_t >       TU( n, min_nm );
        Matrix< value_t >  Ak( A, Range( 0, n-1 ), Range( 0, min_nm-1 ) );
        Matrix< value_t >  Vk( V, Range( 0, m-1 ), Range( 0, min_nm-1 ) );
        
        copy( Ak, TU );
        prod_diag( TU, S, min_nm );
        prod( value_t(-1), TU, adjoint( Vk ), value_t(1), TA );
        
        std::cout << "(BLAS) svd : |M-svd(M)| = " << normF( TA ) << std::endl;
    }
    #endif
    /////////////////////////////////////////////////////////////////

    LOG_TOC( gesvd );
}

//
// compute SVD of matrix A using LAPACK (only left/right sing. vectors)
//
template < typename T1,
           typename T2 >
std::enable_if_t< is_matrix< T1 >::value &&
                  is_vector< T2 >::value &&
                  std::is_same< typename real_type< typename T1::value_t >::type_t,
                                typename T2::value_t >::value, void >
gesvd ( T1 &        A,
        T2 &        S,
        const bool  left )
{
    LOG_COUNTER( gesvd );
    LOG_TIC( gesvd );
    
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    const blas_int_t    n      = blas_int_t( A.nrows() );
    const blas_int_t    m      = blas_int_t( A.ncols() );
    const blas_int_t    min_nm = std::min( n, m );
    blas_int_t          info   = 0;
    value_t             work_query = value_t(0);
    value_t             vdummy = 0;
    real_t              rdummy = 0;
    
    if ( S.length() != size_t( min_nm ) )
        S = std::move( Vector< real_t >( min_nm ) );

    // work space query
    gesvd( left ? 'O' : 'N',
           left ? 'N' : 'O',
           blas_int_t( A.nrows() ),
           blas_int_t( A.ncols() ),
           A.data(),
           blas_int_t( A.col_stride() ),
           S.data(),
           & vdummy,
           blas_int_t( A.col_stride() ),
           & vdummy, 1,
           & work_query,
           LAPACK_WS_QUERY,
           & rdummy,
           info );

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) gesvd", to_string( "argument %d to LAPACK::gesvd", -info ) );
        
    const blas_int_t         lwork = blas_int_t( std::real( work_query ) );
    Hpro::vector< value_t >  work( lwork );
    Hpro::vector< real_t >   rwork( 5 * min_nm );

    MKL_SEQ_START;
    
    gesvd( left ? 'O' : 'N',
           left ? 'N' : 'O',
           blas_int_t( A.nrows() ),
           blas_int_t( A.ncols() ),
           A.data(),
           blas_int_t( A.col_stride() ),
           S.data(),
           & vdummy,
           blas_int_t( A.col_stride() ),
           & vdummy, 1,
           work.data(),
           lwork,
           rwork.data(),
           info );

    MKL_SEQ_END;

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) gesvd", to_string( "argument %d to LAPACK::gesvd", -info ) );
    else if ( info > 0 )
        HERROR( ERR_NCONVERGED, "(BLAS) gesvd", to_string( "in LAPACK::gesvd (info = %d)", info ) );

    LOG_TOC( gesvd );
}

//
// compute singular values of matrix A using LAPACK 
//
template < typename T1,
           typename T2 >
std::enable_if_t< is_matrix< T1 >::value &&
                  is_vector< T2 >::value &&
                  std::is_same< typename real_type< typename T1::value_t >::type_t,
                                typename T2::value_t >::value, void >
gesvd ( T1 &        A,
        T2 &        S )
{
    LOG_COUNTER( gesvd );
    LOG_TIC( gesvd );
    
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    const blas_int_t  n      = blas_int_t( A.nrows() );
    const blas_int_t  m      = blas_int_t( A.ncols() );
    const blas_int_t  min_nm = std::min( n, m );
    blas_int_t        info   = 0;
    value_t           work_query = value_t(0);
    value_t           vdummy = 0;
    real_t            rdummy = 0;

    if ( S.length() != size_t( min_nm ) )
        S = std::move( Vector< real_t >( min_nm ) );

    // work space query
    gesvd( 'N', 'N',
           blas_int_t( A.nrows() ), blas_int_t( A.ncols() ),
           A.data(),
           blas_int_t( A.col_stride() ),
           S.data(),
           & vdummy,
           blas_int_t( A.col_stride() ),
           & vdummy, 1,
           & work_query,
           LAPACK_WS_QUERY,
           & rdummy,
           info );

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) gesvd", to_string( "argument %d to LAPACK::gesvd", -info ) );
        
    const blas_int_t         lwork = blas_int_t( std::real( work_query ) );
    Hpro::vector< value_t >  work( lwork );
    Hpro::vector< real_t >   rwork( 5 * min_nm );

    MKL_SEQ_START;
    
    gesvd( 'N', 'N',
           blas_int_t( A.nrows() ),
           blas_int_t( A.ncols() ), A.data(),
           blas_int_t( A.col_stride() ),
           S.data(),
           & vdummy,
           blas_int_t( A.col_stride() ),
           & vdummy, 1,
           work.data(),
           lwork,
           rwork.data(),
           info );

    MKL_SEQ_END;

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) gesvd", to_string( "argument %d to LAPACK::gesvd", -info ) );
    else if ( info > 0 )
        HERROR( ERR_NCONVERGED, "(BLAS) gesvd", to_string( "in LAPACK::gesvd (info = %d)", info ) );

    LOG_TOC( gesvd );
}

//
// compute SVD of matrix A using LAPACK
//
template < typename T1,
           typename T2,
           typename T3 >
std::enable_if_t< is_matrix< T1 >::value &&
                  is_vector< T2 >::value &&
                  is_matrix< T3 >::value &&
                  std::is_same< typename T1::value_t, typename T3::value_t >::value &&
                  std::is_same< typename real_type< typename T1::value_t >::type_t,
                                typename T2::value_t >::value, void >
gesdd    ( T1 &  A,
           T2 &  S,
           T3 &  V )
{
    LOG_COUNTER( gesdd );
    LOG_TIC( gesdd );
    
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    const blas_int_t   n          = blas_int_t( A.nrows() );
    const blas_int_t   m          = blas_int_t( A.ncols() );
    const blas_int_t   min_nm     = std::min( n, m );
    blas_int_t         info       = 0;
    value_t            work_query = value_t(0);
    value_t            vdummy     = value_t(0);
    real_t             rdummy     = real_t(0);
    blas_int_t         idummy     = blas_int_t(0);
    Matrix< value_t >  VT( min_nm, A.ncols() );

    if ( S.length() != size_t( min_nm ) )
        S = std::move( Vector< real_t >( min_nm ) );

    // work space query
    gesdd( 'S',
           blas_int_t( A.nrows() ),
           blas_int_t( A.ncols() ),
           A.data(),
           blas_int_t( A.col_stride() ),
           S.data(),
           & vdummy,
           blas_int_t( A.col_stride() ),
           & vdummy,
           blas_int_t( VT.col_stride() ),
           & work_query,
           LAPACK_WS_QUERY,
           & rdummy,
           & idummy,
           info );

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) gesdd", to_string( "argument %d to LAPACK::gesdd", -info ) );
        
    const blas_int_t            lwork = blas_int_t( std::real( work_query ) );
    Hpro::vector< value_t >     work( lwork );
    Hpro::vector< blas_int_t >  iwork( 8 * min_nm );
    Hpro::vector< real_t >      rwork( min_nm*std::max(5*min_nm+7,2*std::max(n,m)+2*min_nm+1) );
    Matrix< value_t >           U(  A.nrows(), min_nm );
    
    /////////////////////////////////////////////////////////////////
    #if CHECK_GESVD == 1
    Matrix< value_t >  TA( n, m );

    copy( A, TA );
    #endif
    /////////////////////////////////////////////////////////////////

    MKL_SEQ_START;

    gesdd( 'S',
           blas_int_t( A.nrows() ),
           blas_int_t( A.ncols() ),
           A.data(),
           blas_int_t( A.col_stride() ),
           S.data(),
           U.data(),
           blas_int_t( U.col_stride() ),
           VT.data(),
           blas_int_t( VT.col_stride() ),
           work.data(),
           lwork,
           rwork.data(),
           iwork.data(),
           info );

    MKL_SEQ_END;

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) gesdd", to_string( "argument %d to LAPACK::gesdd", -info ) );
    else if ( info > 0 )
        HERROR( ERR_NCONVERGED, "(BLAS) gesdd", to_string( "in LAPACK::gesdd (info = %d)", info ) );

    if (( V.nrows() != VT.ncols() ) || ( V.ncols() != size_t(min_nm) ))
        V = std::move( Matrix< value_t >( VT.ncols(), min_nm ) );

    copy( adjoint( VT ), V );
    copy( U, A );

    /////////////////////////////////////////////////////////////////
    #if CHECK_GESVD == 1
    {
        Matrix< value_t >       TU( n, min_nm );
        Matrix< value_t >  Ak( A, Range( 0, n-1 ), Range( 0, min_nm-1 ) );
        Matrix< value_t >  Vk( V, Range( 0, m-1 ), Range( 0, min_nm-1 ) );
        
        copy( Ak, TU );
        prod_diag( TU, S, min_nm );
        prod( value_t(-1), TU, adjoint( Vk ), value_t(1), TA );
        
        std::cout << "(BLAS) svd : |M-svd(M)| = " << normF( TA ) << std::endl;
    }
    #endif
    /////////////////////////////////////////////////////////////////

    LOG_TOC( gesdd );
}

template < typename T1,
           typename T2,
           typename T3 >
std::enable_if_t< is_matrix< T1 >::value &&
                  is_vector< T2 >::value &&
                  is_matrix< T3 >::value &&
                  std::is_same< typename T1::value_t, typename T3::value_t >::value &&
                  std::is_same< typename real_type< typename T1::value_t >::type_t,
                                typename T2::value_t >::value, void >
gesvj    ( T1 &   A,
           T2 &   S,
           T3 &   V )
{
    #if HPRO_HAS_GESVJ == 1
    
    LOG_COUNTER( gesvj );
    LOG_TIC( gesvj );
    
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    const blas_int_t    n      = blas_int_t( A.nrows() );
    const blas_int_t    m      = blas_int_t( A.ncols() );
    const blas_int_t    min_nm = std::min( n, m );
    blas_int_t          info   = 0;
    blas_int_t          mv     = 0;

    if ( S.length() != size_t( min_nm ) )
        S = std::move( Vector< real_t >( min_nm ) );

    if (( V.nrows() != A.ncols() ) || ( blas_int_t(V.ncols()) != min_nm ))
        V = std::move( Matrix< value_t >( A.ncols(), min_nm ) );
    
    const blas_int_t         lwork = std::max< blas_int_t >( 6, n + m );
    Hpro::vector< value_t >  cwork( lwork );
    Hpro::vector< real_t >   rwork( lwork );

    // TODO: accuracy based on TTruncAcc
    cwork[0] = 1e8;
    rwork[0] = 1e8;

    MKL_SEQ_START;
    
    gesvj( 'G', 'U', 'V',
           blas_int_t( A.nrows() ),
           blas_int_t( A.ncols() ),
           A.data(),
           blas_int_t( A.col_stride() ),
           S.data(),
           mv,
           V.data(),
           blas_int_t( V.col_stride() ),
           cwork.data(),
           lwork,
           rwork.data(),
           lwork,
           info );

    MKL_SEQ_END;

    if ( info < 0 )
        HERROR( ERR_ARG, "(BLAS) gesvj", to_string( "argument %d to LAPACK::gesvd", -info ) );
    else if ( info > 0 )
        HERROR( ERR_NCONVERGED, "(BLAS) gesvj", "did not converge" );

    LOG_TOC( gesvj );

    #else

    HERROR( ERR_NOT_IMPL, "gesvj", "no support for gesvj available" );

    #endif
}

//
// compute SVD of matrix A (given in \a U) using cross aproximation with
// full pivoting, followed by low rank truncation
//
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
acasvd ( T1 &                                                            U,
         Vector< typename real_type< typename T1::value_t >::type_t > &  S,
         Matrix< typename T1::value_t > &                                V )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    ///////////////////////////////////////////////////////////////////
    //
    // alternative SVD with ACA-Full as first step and recompression
    // using SVD afterwards
    //
    ///////////////////////////////////////////////////////////////////

    const idx_t   n      = idx_t( U.nrows() );
    const idx_t   m      = idx_t( U.ncols() );
    const uint    min_nm = std::min( n, m );

    if ( S.length() != min_nm )
        S = std::move( Vector< real_t >( min_nm ) );

    if (( V.nrows() != m ) || ( V.ncols() != min_nm ))
        V = std::move( Matrix< value_t >( m, min_nm ) );
        
    /////////////////////////////////////////////////////////////////
    #if CHECK_ACASVD == 1
    Matrix< value_t >  TA( n, m );

    copy( U, TA );
    #endif
    /////////////////////////////////////////////////////////////////
    
    //
    // do ACA full on M and put vector pairs into A and B (U and VT)
    // - M will be 
    //

    Matrix< value_t >  M( n, m );
    Matrix< value_t >  A( U );
    Matrix< value_t >  B( V );
    const Range        row_is( 0, n-1 );
    const Range        col_is( 0, m-1 );
    const real_t       eps  = Limits::epsilon<real_t>();
    real_t             apr  = eps;
    const uint         maxk = min_nm;
    idx_t              k    = 0;
    
    BLAS::copy( U, M );
    
    while ( k < maxk )
    {
        //
        // look for maximal element
        //

        idx_t    pivot_i, pivot_j;
        value_t  pivot_val;
        
        max_idx( M, pivot_i, pivot_j );
        pivot_val = M( pivot_i, pivot_j );

        // stop if maximal element is 0
        if ( Math::abs( pivot_val ) == 0.0 )
            break;

        Vector< value_t >  pivot_col( M, row_is, pivot_j );
        Vector< value_t >  pivot_row( M, pivot_i, col_is );
        
        //
        // copy row and column into A/B and update D
        //

        Vector< value_t >  col( A, row_is, k );
        Vector< value_t >  row( B, col_is, k );

        copy( pivot_col, col );
        copy( pivot_row, row );
        conj( row );  // row is stored in B in transposed/adjoint form
        scale( value_t(1) / Math::conj(pivot_val), row );

        ++k;
        
        //
        // look at norm of residual
        //
            
        const real_t  norm = norm2( col ) * norm2( row );

        if      ( k == 1 )     apr *= norm;
        else if ( norm < apr ) { k--; break; }

        //
        // update dense matrix
        //
        
        add_r1( value_t(-1), col, row, M );
    }// while

    if ( k == 0 )
    {
        //
        // set singular values to 0
        //

        fill( real_t(0), S );

        return;
    }// if

    const Range  rank_is( 0, k-1 );

    /////////////////////////////////////////////////////////////////
    #if CHECK_ACASVD == 1
    {
        Matrix< value_t >  TT( n, m );
        Matrix< value_t >  Ak( A, row_is, rank_is );
        Matrix< value_t >  Bk( B, col_is, rank_is );

        prod( value_t(1), Ak, adjoint(Bk), value_t(0), TT );
        add( value_t(-1), TA, TT );
        std::cout << "(BLAS) acasvd : |M-aca(M)|/|M| = " << 
                  << normF( TT ) / normF( TA ) << std::endl;
    }
    #endif
    /////////////////////////////////////////////////////////////////
    
    //
    // first k vectors in A and B are our low rank factorisation
    // now use A·B^H = Q_A·R_A·(Q_B·R_B)^H
    //               = Q_A·R_A·R_B^H·Q_B^H
    // and compute SVD for R_A · R_B^H to compute (Q_A·U)·S·(V^H·Q_B^H)
    // as SVD for given matrix
    //        

    Matrix< value_t >  Ak( A, row_is, rank_is );
    Matrix< value_t >  Bk( B, col_is, rank_is );
    Matrix< value_t >  QA( n, k );
    Matrix< value_t >  QB( m, k );
    Matrix< value_t >  RA( k, k );
    Matrix< value_t >  RB( k, k );

    //
    // do QR-factorisation of A and B
    //

    copy( Ak, QA );
    copy( Bk, QB );
    
    qr( QA, RA );
    qr( QB, RB );
    
    // RA = RA·RB^T
    Matrix< value_t >  R( k, k );

    prod( value_t(1), RA, adjoint(RB), value_t(0), R );
    
    //
    // SVD(R) = U S V^T
    //

    gesvd( R, S, RB );

    //
    // compute final SVD by U = QA · U_R and V^T = V_R^T · Q_B^T
    //

    Matrix< value_t >  Uk( U, row_is, rank_is );
    Matrix< value_t >  Vk( V, col_is, rank_is );

    prod( value_t(1), QA, R,  value_t(0), Uk );
    prod( value_t(1), QB, RB, value_t(0), Vk );
    
    // BLAS::gemm( MATOP_NORM, MATOP_NORM, n, k, k, 1.0, QA,  n, U_R, k, 0.0, U,  n );
    // BLAS::gemm( MATOP_NORM, MATOP_ADJ, k, m, k, 1.0, V_R, k, QB,  m, 0.0, VT, min_nm );

    /////////////////////////////////////////////////////////////////
    #if CHECK_ACASVD == 1
    {
        Matrix< value_t > TT( n, m );
        Matrix< value_t > TU( n, k );
        
        BLAS::copy( Uk, TU );
        prod_diag( TU, S, k );
        prod( value_t(1), TU, adjoint(Vk), value_t(0), TT );
        add( value_t(-1), TA, TT );
         
        std::cout << "(BLAS) acasvd : |M-svd(M)|/|M| = " << 
                  << normF( TT ) / normF( TA ) << std::endl;
    }
    #endif // CHECK_ACASVD
    /////////////////////////////////////////////////////////////////
}

//
// compute SVD of A (given in \a U) using LAPACK gesvd functions        
//
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
lasvd ( T1 &                                                            U,
        Vector< typename real_type< typename T1::value_t >::type_t > &  S,
        Matrix< typename T1::value_t > &                                V )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    const idx_t            n  = idx_t( U.nrows() );
    const idx_t            m  = idx_t( U.ncols() );
    idx_t                  ln = n;
    idx_t                  lm = m;
    Hpro::vector< idx_t >  row_pos;
    Hpro::vector< idx_t >  col_pos;

    if ( CFG::BLAS::check_zeroes )
    {
        ln = 0;
        lm = 0;
    
        //
        // check if we have zero columns or rows in the matrix
        // and reduce matrix to non-zero rows/columns
        //

        row_pos.resize( n );
        col_pos.resize( m );

        // check columns
        {
            const Range  row_is( 0, n-1 );
        
            for ( idx_t  j = 0; j < m; j++ )
            {
                Vector< value_t >  col_j( U, row_is, j );
        
                if ( norm2( col_j ) != real_t(0) )
                {
                    if ( lm != j )
                    {
                        Vector< value_t >  col_new( U, row_is, lm );

                        copy( col_j, col_new );
                    }// if

                    col_pos[j] = lm;
                    lm++;
                }// if
                else
                    col_pos[j] = -1;
            }// for
        }
    
        // check rows
        {
            const Range  col_is( 0, lm-1 );
        
            for ( idx_t  i = 0; i < n; i++ )
            {
                Vector< value_t >  row_i( U, i, col_is );

                if ( norm2( row_i ) != real_t(0) )
                {
                    if ( ln != i )
                    {
                        Vector< value_t >  row_new( U, ln, col_is );

                        copy( row_i, row_new );
                    }// if
                
                    row_pos[i] = ln;
                    ln++;
                }// if
                else
                    row_pos[i] = -1;
            }// for
        }
    }// if
    
    if (( ln == 0 ) || ( lm == 0 ))
    {
        // matrices are empty
        fill( real_t(0),  S );
        fill( value_t(0), V );
    }// if
    else if (( ln == n ) && ( lm == m ))
    {
        #if HPRO_HAS_GESVJ == 1
        if ( CFG::BLAS::use_gesvj && ( U.nrows() >= U.ncols() ))
        {
            gesvj( U, S, V );
        }// if
        else
        #endif
        {
            if ( std::min(n,m) < CFG::BLAS::gesvd_limit )
                gesvd( U, S, V );
            else
                gesdd( U, S, V );
        }// else
    }// if
    else
    {
        if ( CFG::BLAS::check_zeroes )
        {
            //
            // perform restricted SVD
            //
        
            const idx_t        lmin = std::min( ln, lm );
            Matrix< value_t >  TU( U, Range( 0, ln-1 ), Range( 0, lm-1 ) );
            Matrix< value_t >  TV( V, Range( 0, lm-1 ), Range( 0, lmin-1 ) );
            Vector< real_t>    TS( S, Range( 0, lmin-1 ) );

            #if HPRO_HAS_GESVJ == 1
            if ( CFG::BLAS::use_gesvj && ( U.nrows() >= U.ncols() ))
            {
                gesvj( TU, TS, TV );
            }// if
            else
            #endif
            {
                if ( lmin < CFG::BLAS::gesvd_limit )
                    gesvd( TU, TS, TV );
                else
                    gesdd( TU, TS, TV );
            }// else

            //
            // refill all previously removed rows/columns with zero
            // and put rows/column into original position
            //

            for ( idx_t  i = n-1; i >= 0; --i )
            {
                Vector< value_t >  u_i( U.row( i ) );
                const idx_t        i_temp = row_pos[i];
             
                if ( i_temp == -1 )
                {
                    fill( value_t(0), u_i );
                }// if
                else if ( i_temp != i )
                {
                    Vector< value_t >  u_p( U.row( i_temp ) );

                    copy( u_p, u_i );
                }// else
            }// for
        
            for ( idx_t  i = m-1; i >= 0; --i )
            {
                Vector< value_t >  v_i( V.row( i ) );
                const idx_t        i_temp = col_pos[i];
            
                if ( i_temp == -1 )
                {
                    fill( value_t(0), v_i );
                }// if
                else if ( i_temp != i )
                {
                    Vector< value_t >  v_p( V.row( i_temp ) );

                    copy( v_p, v_i );
                }// else
            }// for
        }// if
    }// else
}

//
// compute SVD of A (given in \a U) using LAPACK gesvd functions
// but return only left/right singular vectors
//
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
lasvd ( T1 &                                                            U,
        Vector< typename real_type< typename T1::value_t >::type_t > &  S,
        const bool                                                      left )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    const idx_t            n  = idx_t( U.nrows() );
    const idx_t            m  = idx_t( U.ncols() );
    idx_t                  ln = n;
    idx_t                  lm = m;
    Hpro::vector< idx_t >  row_pos;
    Hpro::vector< idx_t >  col_pos;

    if ( CFG::BLAS::check_zeroes )
    {
        ln = 0;
        lm = 0;
    
        //
        // check if we have zero columns or rows in the matrix
        // and reduce matrix to non-zero rows/columns
        //

        row_pos.resize( n );
        col_pos.resize( m );

        // check columns
        {
            const Range  row_is( 0, n-1 );
        
            for ( idx_t  j = 0; j < m; j++ )
            {
                Vector< value_t >  col_j( U, row_is, j );
        
                if ( norm2( col_j ) != real_t(0) )
                {
                    if ( lm != j )
                    {
                        Vector< value_t >  col_new( U, row_is, lm );

                        copy( col_j, col_new );
                    }// if

                    col_pos[j] = lm;
                    lm++;
                }// if
                else
                    col_pos[j] = -1;
            }// for
        }
    
        // check rows
        {
            const Range  col_is( 0, lm-1 );
        
            for ( idx_t  i = 0; i < n; i++ )
            {
                Vector< value_t >  row_i( U, i, col_is );

                if ( norm2( row_i ) != real_t(0) )
                {
                    if ( ln != i )
                    {
                        Vector< value_t >  row_new( U, ln, col_is );

                        copy( row_i, row_new );
                    }// if
                
                    row_pos[i] = ln;
                    ln++;
                }// if
                else
                    row_pos[i] = -1;
            }// for
        }
    }// if
    
    if (( ln == 0 ) || ( lm == 0 ))
    {
        // matrices are empty
        S = std::move( Vector< real_t >() );
    }// if
    else if (( ln == n ) && ( lm == m ))
    {
        gesvd( U, S, left );
    }// if
    else
    {
        if ( CFG::BLAS::check_zeroes )
        {
            //
            // perform restricted SVD
            //
        
            const idx_t        lmin = std::min( ln, lm );
            Matrix< value_t >  TU( U, Range( 0, ln-1 ), Range( 0, lm-1 ) );
            Vector< real_t>    TS( S, Range( 0, lmin-1 ) );

            gesvd( TU, TS, left );

            //
            // refill all previously removed rows/columns with zero
            // and put rows/column into original position
            //

            for ( idx_t  i = n-1; i >= 0; --i )
            {
                Vector< value_t >  u_i( U.row( i ) );
                const idx_t        i_temp = row_pos[i];
            
                if ( i_temp == -1 )
                {
                    fill( value_t(0), u_i );
                }// if
                else if ( i_temp != i )
                {
                    Vector< value_t >  u_p( U.row( i_temp ) );

                    copy( u_p, u_i );
                }// else
            }// for
        }// if
    }// else
}

//
// compute SVD of A (given in \a U) using LAPACK gesvd functions
// but return only singular values
//
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
lasvd ( T1 &                                                            U,
        Vector< typename real_type< typename T1::value_t >::type_t > &  S )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    const idx_t            n  = idx_t( U.nrows() );
    const idx_t            m  = idx_t( U.ncols() );
    idx_t                  ln = n;
    idx_t                  lm = m;
    Hpro::vector< idx_t >  row_pos;
    Hpro::vector< idx_t >  col_pos;

    if ( CFG::BLAS::check_zeroes )
    {
        ln = 0;
        lm = 0;
    
        //
        // check if we have zero columns or rows in the matrix
        // and reduce matrix to non-zero rows/columns
        //

        row_pos.resize( n );
        col_pos.resize( m );

        // check columns
        {
            const Range  row_is( 0, n-1 );
        
            for ( idx_t  j = 0; j < m; j++ )
            {
                Vector< value_t >  col_j( U, row_is, j );
        
                if ( norm2( col_j ) != real_t(0) )
                {
                    if ( lm != j )
                    {
                        Vector< value_t >  col_new( U, row_is, lm );

                        copy( col_j, col_new );
                    }// if

                    col_pos[j] = lm;
                    lm++;
                }// if
                else
                    col_pos[j] = -1;
            }// for
        }
    
        // check rows
        {
            const Range  col_is( 0, lm-1 );
        
            for ( idx_t  i = 0; i < n; i++ )
            {
                Vector< value_t >  row_i( U, i, col_is );

                if ( norm2( row_i ) != real_t(0) )
                {
                    if ( ln != i )
                    {
                        Vector< value_t >  row_new( U, ln, col_is );

                        copy( row_i, row_new );
                    }// if
                
                    row_pos[i] = ln;
                    ln++;
                }// if
                else
                    row_pos[i] = -1;
            }// for
        }
    }// if
    
    if (( ln == 0 ) || ( lm == 0 ))
    {
        // matrices are empty
        fill( real_t(0),  S );
    }// if
    else if (( ln == n ) && ( lm == m ))
    {
        gesvd( U, S );
    }// if
    else
    {
        if ( CFG::BLAS::check_zeroes )
        {
            //
            // perform restricted SVD
            //
        
            const idx_t        lmin = std::min( ln, lm );
            Matrix< value_t >  TU( U, Range( 0, ln-1 ), Range( 0, lm-1 ) );
            Vector< real_t>    TS( S, Range( 0, lmin-1 ) );

            gesvd( TU, TS );
        }// if
    }// else
}

//!
//! compute SVD decomposition A = U·S·V^H of the nxm matrix \a A with
//! n×min(n,m) matrix U, min(n,m)×min(n,m) matrix S (diagonal)
//! and m×min(n,m) matrix V; \a A will be overwritten with U upon exit
//!
template <typename T>
void
svd_double ( Matrix< T > &                               U,
             Vector< typename real_type< T >::type_t > & S,
             Matrix< T > &                               V )
{
    #if USE_ACASVD == 1
    acasvd( U, S, V );
    #else
    lasvd( U, S, V );
    #endif
}

namespace
{

// template version of GPU::*_svd_size
template < typename value_t >  size_t  gpu_svd_size ();
template <>                    size_t  gpu_svd_size< float >                  () { return CFG::GPU::fp32_svd_size; }
template <>                    size_t  gpu_svd_size< double >                 () { return CFG::GPU::fp64_svd_size; }
template <>                    size_t  gpu_svd_size< std::complex< float > >  () { return CFG::GPU::fp32_svd_size; }
template <>                    size_t  gpu_svd_size< std::complex< double > > () { return CFG::GPU::fp64_svd_size; }

}// namespace anonymous


template < typename value_t >
void
svd    ( Matrix< value_t > &                                U,
         Vector< typename real_type< value_t >::type_t > &  S,
         Matrix< value_t > &                                V )
{
    DO_CHECK_SMALL(   U );
    DO_CHECK_INF_NAN( U, "(BLAS) svd", "in input matrix" );

    //
    // use double precision version to increase accuracy
    //
    
    if ( CFG::BLAS::use_double_prec && is_single_prec< value_t >::value )
    {
        return svd_double( U, S, V );
    }// if
    
    //
    // try GPU
    //

    bool  use_cpu = true;
        
    if ( std::min( U.nrows(), U.ncols() ) >= gpu_svd_size< value_t >() )
    {
        const auto  handle = GPU::request_handle();

        if ( handle != GPU::INVALID_HANDLE )
        {
            if ( U.nrows() < U.ncols() )
            {
                // only support nrows >= ncols, so transpose
                if (( V.nrows() != U.ncols() ) || ( V.ncols() != U.nrows() ))
                    V = Matrix< value_t >( U.ncols(), U.nrows() );

                copy( adjoint( U ), V );

                // LOG::print( "using GPU (adjoint)" );
                
                if ( GPU::svd( handle, V, S, U ) )
                    use_cpu = false;
            }// if
            else
            {
                // LOG::print( "using GPU" );

                if ( GPU::svd( handle, U, S, V ) )
                    use_cpu = false;
            }// else
        }// if

        GPU::release_handle( handle );
    }// if

    //
    // use (or fall back to) CPU
    //
    
    if ( use_cpu )
    {
        #if USE_ACASVD == 1

        acasvd( U, S, V );
        
        #else  // USE_ACASVD == 0
        
        lasvd( U, S, V );
        
        #endif
    }// if

    DO_CHECK_INF_NAN( U, "(BLAS) svd", "in output matrix" );
    DO_CHECK_INF_NAN( V, "(BLAS) svd", "in output matrix" );
}

//
// uses double precision to compute SVD for single precision data
//
template <>
void
svd_double< float > ( Matrix< float > &  U,
                      Vector< float > &  S,
                      Matrix< float > &  V )
{
    Matrix< double >  Ud( U.nrows(), U.ncols() );
    Vector< double >  Sd;
    Matrix< double >  Vd;

    for ( idx_t  j = 0; j < idx_t(U.ncols()); ++j )
        for ( idx_t  i = 0; i < idx_t(U.nrows()); ++i )
            Ud( i, j ) = U( i, j );
    
    svd( Ud, Sd, Vd );

    if ( S.length() != Sd.length() )
        S = Vector< float >( Sd.length() );
    
    if ( V.nrows() * V.ncols() != Vd.nrows() * Vd.ncols() )
        V = Matrix< float >( Vd.nrows(), Vd.ncols() );
    
    for ( idx_t  j = 0; j < idx_t(Ud.ncols()); ++j )
        for ( idx_t  i = 0; i < idx_t(Ud.nrows()); ++i )
            U( i, j ) = float( Ud( i, j ) );
    
    for ( idx_t  i = 0; i < idx_t(Sd.length()); ++i )
        S( i ) = float( Sd( i ) );
    
    for ( idx_t  j = 0; j < idx_t(Vd.ncols()); ++j )
        for ( idx_t  i = 0; i < idx_t(Vd.nrows()); ++i )
            V( i, j ) = float( Vd( i, j ) );
}

template <>
void
svd_double< std::complex< float > > ( Matrix< std::complex< float > > &  U,
                                      Vector< float > &             S,
                                      Matrix< std::complex< float > > &  V )
{
    Matrix< std::complex< double > >  Ud( U.nrows(), U.ncols() );
    Vector< double >             Sd;
    Matrix< std::complex< double > >  Vd;

    for ( idx_t  j = 0; j < idx_t(U.ncols()); ++j )
        for ( idx_t  i = 0; i < idx_t(U.nrows()); ++i )
            Ud( i, j ) = U( i, j );
    
    svd( Ud, Sd, Vd );

    if ( S.length() != Sd.length() )
        S = Vector< float >( Sd.length() );
    
    if ( V.nrows() * V.ncols() != Vd.nrows() * Vd.ncols() )
        V = Matrix< std::complex< float > >( Vd.nrows(), Vd.ncols() );
    
    for ( idx_t  j = 0; j < idx_t(Ud.ncols()); ++j )
        for ( idx_t  i = 0; i < idx_t(Ud.nrows()); ++i )
            U( i, j ) = std::complex< float >( float(std::real( Ud( i, j ) )),
                                               float(std::imag( Ud( i, j ) )) );
    
    for ( idx_t  i = 0; i < idx_t(Sd.length()); ++i )
        S( i ) = float( Sd( i ) );
    
    for ( idx_t  j = 0; j < idx_t(Vd.ncols()); ++j )
        for ( idx_t  i = 0; i < idx_t(Vd.nrows()); ++i )
            V( i, j ) = std::complex< float >( float(std::real( Vd( i, j ) )),
                                               float(std::imag( Vd( i, j ) )) );
}


//!
//! compute SVD decomposition A = U·S·V^H of the nxm matrix \a A
//! but return only the singular values S ∈ ℝ^min(n,m);
//! \a A will be overwritten with U upon exit
//!
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
svd    ( T1 &                                                            A,
         Vector< typename real_type< typename T1::value_t >::type_t > &  S,
         const bool                                                      left )
{
    DO_CHECK_SMALL(   A );
    DO_CHECK_INF_NAN( A, "(BLAS) svd", "in input matrix" );

    lasvd( A, S, left );
}

//!
//! compute SVD decomposition A = U·S·V^H of the nxm matrix \a A
//! but return only the singular values S ∈ ℝ^min(n,m);
//! \a A will be overwritten with U upon exit
//!
template < typename value_t >
void
sv     ( Matrix< value_t > &                 A,
         Vector< real_type_t< value_t > > &  S )
{
    DO_CHECK_SMALL(   A );
    DO_CHECK_INF_NAN( A, "(BLAS) sv", "in input matrix" );

    lasvd( A, S );
}

//!
//! compute SVD decomposition \f$ M = A·B^H = U·S·V^H \f$ of the nxm 
//! low-rank matrix \a M but return only the singular values S ∈ ℝ^min(n,m);
//! \a A and \a B will be overwritten upon exit
//!
template < typename value_t >
void
sv ( Matrix< value_t > &                 A,
     Matrix< value_t > &                 B,
     Vector< real_type_t< value_t > > &  S )
{
    // using  value_t = typename T1::value_t;

    DO_CHECK_SMALL(   A );
    DO_CHECK_SMALL(   B );
    DO_CHECK_INF_NAN( A, "(BLAS) sv", "in input matrix A" );
    DO_CHECK_INF_NAN( B, "(BLAS) sv", "in input matrix B" );

    if ( A.ncols() != B.ncols() )
        HERROR( ERR_MAT_SIZE, "(BLAS) sv", "rank in A and B differs" );

    const idx_t   n   = idx_t( A.nrows() );
    const idx_t   m   = idx_t( B.nrows() );
    const idx_t   k   = idx_t( A.ncols() );
    const idx_t   mrc = std::min(n, m);

    if ( k >= mrc )
    {
        //
        // A·B^H = U·S·V^H
        //
            
        Matrix< value_t >  M( n, m );

        prod( value_t(1), A, adjoint(B), value_t(0), M );

        sv( M, S );
    }// if
    else
    {
        Matrix< value_t >  RA( k, k );
        Matrix< value_t >  RB( k, k );
        Matrix< value_t >  R( k, k );

        //
        // A·B^H = Q_A·R_A·(Q_B·R_B)^H = Q_A·(R_A·R_B^H)·Q_B^H = Q_A·(U·S·V^H)·Q_B^H
        //
        
        qr( A, RA );
        qr( B, RB );
        
        prod( value_t(1), RA, adjoint(RB), value_t(0), R );
            
        sv( R, S );
    }// else
}

//!
//! approximate given dense matrix \a M by low rank matrix
//! according to accuracy \a acc. The low rank matrix will
//! be stored in \a A and \a B
//!
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, size_t >
approx_svd ( T1 &                              M,
             const TTruncAcc &                 acc,
             Matrix< typename T1::value_t > &  A,
             Matrix< typename T1::value_t > &  B )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    DO_CHECK_SMALL(   M );
    DO_CHECK_INF_NAN( M, "(BLAS) approx_svd", "in input matrix" );

    //
    // perform SVD of M
    //

    const idx_t        n   = idx_t( M.nrows() );
    const idx_t        m   = idx_t( M.ncols() );
    const idx_t        mrc = std::min(n,m);
    Vector< real_t >   S( mrc );
    Matrix< value_t >  V( m, mrc );

    svd( M, S, V );
        
    // determine truncated rank based on singular values
    const idx_t  k = idx_t( acc.trunc_rank( S ) );

    //
    // u_i -> a_i, v_i -> b_i
    // scale smaller one of A and B by S
    //

    const Range        row_is( 0, n-1 );
    const Range        col_is( 0, m-1 );
    Matrix< value_t >  Uk( M, row_is, Range( 0, k-1 ) );
    Matrix< value_t >  Vk( V, col_is, Range( 0, k-1 ) );
    
    A = std::move( Matrix< value_t >( n, k ) );
    B = std::move( Matrix< value_t >( m, k ) );

    copy( Uk, A );
    copy( Vk, B );

    if ( n < m ) prod_diag( A, S, k );
    else         prod_diag( B, S, k );

    DO_CHECK_SMALL(   A );
    DO_CHECK_SMALL(   B );
    DO_CHECK_INF_NAN( A, "(BLAS) approx_svd", "in output matrix" );
    DO_CHECK_INF_NAN( B, "(BLAS) approx_svd", "in output matrix" );

    return k;
}

template < typename T >
std::pair< Matrix< T >, Matrix< T > >
approx_svd ( Matrix< T > &      M,
             const TTruncAcc &  acc )
{
    using  value_t = T;
    using  real_t  = typename real_type< value_t >::type_t;

    //
    // perform SVD of M
    //

    const idx_t        n   = idx_t( M.nrows() );
    const idx_t        m   = idx_t( M.ncols() );
    const idx_t        mrc = std::min(n,m);
    Vector< real_t >   S( mrc );
    Matrix< value_t >  V( m, mrc );

    svd( M, S, V );
        
    // determine truncated rank based on singular values
    const idx_t  k = idx_t( acc.trunc_rank( S ) );

    //
    // u_i -> a_i, v_i -> b_i
    // scale smaller one of A and B by S
    //

    const Range        row_is( 0, n-1 );
    const Range        col_is( 0, m-1 );
    Matrix< value_t >  Uk( M, row_is, Range( 0, k-1 ) );
    Matrix< value_t >  Vk( V, col_is, Range( 0, k-1 ) );
    
    Matrix< value_t >  A( n, k );
    Matrix< value_t >  B( m, k );

    copy( Uk, A );
    copy( Vk, B );

    if ( n < m ) prod_diag( A, S, k );
    else         prod_diag( B, S, k );

    return { std::move( A ), std::move( B ) };
}

//!
//! approximate given dense matrix \a M by low rank matrix
//! according to accuracy \a acc. The low rank matrix will
//! be stored in \a A and \a B
//!
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, size_t >
approx_rrqr ( T1 &                              M,
              const TTruncAcc &                 acc,
              Matrix< typename T1::value_t > &  A,
              Matrix< typename T1::value_t > &  B )
{
    using  value_t = typename T1::value_t;

    DO_CHECK_SMALL(   M );
    DO_CHECK_INF_NAN( M, "(BLAS) approx_rrqr", "in input matrix" );

    //
    // check dimensions since QR only works with n >= m,
    // hence if m > n, transpose and call again with exchanged A/B
    //

    const idx_t  n = idx_t( M.nrows() );
    const idx_t  m = idx_t( M.ncols() );

    if ( m > n )
    {
        Matrix< value_t >  MH( m, n );

        copy( adjoint( M ), MH );

        const auto  k = approx_rrqr( MH, acc, B, A );

        return k;
    }// if
    
    //
    // perform column pivoted QR of M
    //

    const idx_t                mrc = std::min(n,m);
    Matrix< value_t >          R( mrc, m );
    std::vector< blas_int_t >  P( m, 0 );
    const auto                 k = qrp_trunc( M, R, P, acc );

    //
    // restrict first k columns
    //

    Matrix< value_t >  Qk( M, Range::all, Range( 0, k-1 ) );

    // A = Q_k
    A = std::move( Matrix< value_t >( n, k ) );
    copy( Qk, A );

    // copy first k columns of R^T to B, i.e., first k rows of R
    Matrix< value_t >  Rk( R, Range( 0, k-1 ), Range::all );
    Matrix< value_t >  TB( m, k );
    
    copy( adjoint( Rk ), TB );
    
    // then permute rows of TB (do P·R^T) and copy to B
    B = std::move( Matrix< value_t >( m, k ) );
    
    for ( int i = 0; i < m; ++i )
    {
        auto  j    = P[i];
        auto  TB_i = TB.row( i );
        auto  B_j  = B.row( j );

        copy( TB_i, B_j );
    }// for
    
    DO_CHECK_SMALL(   A );
    DO_CHECK_SMALL(   B );
    DO_CHECK_INF_NAN( A, "(BLAS) approx_rrqr", "in output matrix" );
    DO_CHECK_INF_NAN( B, "(BLAS) approx_rrqr", "in output matrix" );

    return k;
}

//
// compute basis for column space (range) of M
//
template < typename value_t >
Matrix< value_t >
column_basis ( const Matrix< value_t > &  M,
               const TTruncAcc &          acc,
               const uint                 power_steps,
               const uint                 oversampling )
{
    const idx_t  n = idx_t( M.nrows() );
    const idx_t  m = idx_t( M.ncols() );

    if ( acc.is_fixed_rank() )
    {
        const auto  k = idx_t(acc.rank());
        auto        T = random< value_t >( m, k + oversampling );
        auto        Y = prod( value_t(1), M, T );

        //
        // power iteration
        //
            
        Matrix< value_t >  MtQ( m, k + oversampling );
        Matrix< value_t >  R( k + oversampling, k + oversampling );
        
        for ( uint  j = 0; j < power_steps; ++j )
        {
            qr( Y, R );
            prod( value_t(1), adjoint(M), Y, value_t(0), MtQ );

            qr( MtQ, R );
            prod( value_t(1), M, MtQ, value_t(0), Y );
        }// for

        qr( Y, R );

        return Y;
    }// if
    else
    {
        Matrix< value_t >  A( M, copy_value );
        auto               norm_M  = normF( M );
        const auto         rel_eps = acc.rel_eps();
        const auto         abs_eps = acc.abs_eps();
        const uint         bsize   = std::min< uint >( std::max< uint >( CFG::BLAS::sample_size, 1 ),
                                                       std::min< uint >( n, m ) );
        const uint         nblocks = std::min< uint >( n, m ) / bsize;
        auto               Qs = std::list< Matrix< value_t > >();
        Matrix< value_t >  T_i( m, bsize );
        TRNG               rng;

        #if 1
        using  real_t  = typename real_type< value_t >::type_t;
        
        for ( uint  i = 0; i < nblocks; ++i )
        {
            fill_rand_normal( T_i );
            
            auto  Q_i  = prod( value_t(1), M, T_i ); // Y_i
            auto  TQ_i = copy( Q_i );
            
            for ( auto  Q_j : Qs )
            {
                const auto  QhQi = prod( value_t(1), adjoint( Q_j ), TQ_i );

                prod( value_t(-1), Q_j, QhQi, value_t(1), Q_i );
            }// for

            real_t  norm_Qi = real_t(0);

            for ( uint  j = 0; j < bsize; ++j )
            {
                const auto  Qi_j = Q_i.column( j );

                norm_Qi = std::max( norm_Qi, norm2( Qi_j ) );
            }// for

            // use first approximation also as approximation of norm of M
            if ( i == 0 )
                norm_M = norm_Qi;

            //
            // power iteration
            //
            
            Matrix< value_t >  R( bsize, bsize );
            Matrix< value_t >  AtQ( m, bsize );
            
            for ( uint  j = 0; j < power_steps; ++j )
            {
                // copy( Y_i, Q_i );
                qr( Q_i, R );
                prod( value_t(1), adjoint(A), Q_i, value_t(0), AtQ );
                
                qr( AtQ, R );
                prod( value_t(1), A, AtQ, value_t(0), Q_i );  // Q_i = Y_i
            }// for
            
            // copy( Y_i, Q_i );
            qr( Q_i, R );
            
            //
            // project Q_i away from previous Q_j
            //
            //    Q_i = Q_i - [ Q_0 .. Q_i-1 ] [ Q_0 .. Q_i-1 ]^H Q_i = Q_i - Σ_j=0^i-1 Q_j Q_j^H Q_i
            //
                
            if ( i > 0 )
            {
                Matrix< value_t >  C_i( Q_i, copy_value );
                Matrix< value_t >  QjtQi( bsize, bsize );
                
                for ( const auto &  Q_j : Qs )
                {
                    prod( value_t(1), adjoint(Q_j), C_i, value_t(0), QjtQi );
                    prod( value_t(-1), Q_j, QjtQi, value_t(1), Q_i );
                }// for
                
                qr( Q_i, R );
            }// if
            
            //
            // A = A - Q_i Q_i^t A
            //

            // auto  B_i = prod( value_t(1), adjoint(Q_i), A );
            
            // prod( value_t(-1), Q_i, B_i, value_t(1), A );
            
            // const auto  norm_A = normF( A );

            // LOG::printf( "|A| = %.6e, %6e    %6e, %6e", norm_A, abs_eps, norm_A / norm_M, rel_eps );

            Qs.push_back( std::move( Q_i ) );
            
            if (( norm_Qi <= abs_eps ) || (( norm_Qi ) <= rel_eps * norm_M ))
                break;
        }// for
        #else
        for ( uint  i = 0; i < nblocks; ++i )
        {
            fill_rand_normal( T_i );
            
            auto  Q_i = prod( value_t(1), A, T_i ); // Y_i
            
            //
            // power iteration
            //
            
            Matrix< value_t >  R( bsize, bsize );
            Matrix< value_t >  AtQ( m, bsize );
            
            for ( uint  j = 0; j < q; ++j )
            {
                // copy( Y_i, Q_i );
                qr( Q_i, R );
                prod( value_t(1), adjoint(A), Q_i, value_t(0), AtQ );
                
                qr( AtQ, R );
                prod( value_t(1), A, AtQ, value_t(0), Q_i );  // Q_i = Y_i
            }// for
            
            // copy( Y_i, Q_i );
            qr( Q_i, R );
            
            //
            // project Q_i away from previous Q_j
            //
            //    Q_i = Q_i - [ Q_0 .. Q_i-1 ] [ Q_0 .. Q_i-1 ]^H Q_i = Q_i - Σ_j=0^i-1 Q_j Q_j^H Q_i
            //
                
            if ( i > 0 )
            {
                Matrix< value_t >  C_i( Q_i, copy_value );
                Matrix< value_t >  QjtQi( bsize, bsize );
                
                for ( const auto &  Q_j : Qs )
                {
                    prod( value_t(1), adjoint(Q_j), C_i, value_t(0), QjtQi );
                    prod( value_t(-1), Q_j, QjtQi, value_t(1), Q_i );
                }// for
                
                qr( Q_i, R );
            }// if
            
            //
            // A = A - Q_i Q_i^t A
            //

            auto  B_i = prod( value_t(1), adjoint(Q_i), A );
            
            prod( value_t(-1), Q_i, B_i, value_t(1), A );
            
            const auto  norm_A = normF( A );

            // LOG::printf( "|A| = %.6e, %6e    %6e, %6e", norm_A, abs_eps, norm_A / norm_M, rel_eps );

            Qs.push_back( std::move( Q_i ) );
            
            if (( norm_A < abs_eps ) || (( norm_A / norm_M ) < rel_eps ))
                break;
        }// for
        #endif
        
        //
        // collect Q_i's into final result
        //

        Matrix< value_t >  Q( n, Qs.size() * bsize );
        idx_t              pos = 0;

        for ( const auto &  Q_i : Qs )
        {
            Matrix< value_t >  Q_sub( Q, Range::all, Range( pos * bsize, (pos+1)*bsize - 1 ) );

            copy( Q_i, Q_sub );
            ++pos;
        }// for

        return Q;
    }// else
}

//
// compute basis for column space (range) of M = A·B^H
//
template < typename value_t >
Matrix< value_t >
column_basis ( const Matrix< value_t > &  IA,
               const Matrix< value_t > &  IB,
               const TTruncAcc &          acc,
               const uint                 power_steps,
               const uint                 oversampling )
{
    const idx_t  n    = idx_t( IA.nrows() );
    const idx_t  m    = idx_t( IB.nrows() );
    const idx_t  rank = idx_t( IA.ncols() );
    
    if ( acc.is_fixed_rank() )
    {
        auto         A = copy( IA );
        auto         B = copy( IB );
        const idx_t  k   = idx_t(acc.rank());
        auto         T   = random< value_t >( n, k + oversampling );
        auto         BtT = prod( value_t(1), adjoint(B), T );
        auto         Y   = prod( value_t(1), A, BtT );

        //
        // power iteration
        //
            
        Matrix< value_t >  AtQ( rank, k + oversampling );
        Matrix< value_t >  BAtQ( m, k + oversampling );
        Matrix< value_t >  R( k + oversampling, k + oversampling );
        
        for ( uint  j = 0; j < power_steps; ++j )
        {
            // [Y,R] = qr(Y); MtQ = M^H·Y = B·A^H·Y
            qr( Y, R );
            prod( value_t(1), adjoint(A), Y, value_t(0), AtQ );
            prod( value_t(1), B, AtQ, value_t(0), BAtQ );

            // [Q,R] = qr(B·A^H·Y); Y = A·B^H·Q
            qr( BAtQ, R );
            prod( value_t(1), adjoint(B), BAtQ, value_t(0), AtQ );
            prod( value_t(1), A, AtQ, value_t(0), Y );
        }// for

        qr( Y, R );

        return Y;
    }// if
    else
    {
        Matrix< value_t >  A( IA, copy_value );
        Matrix< value_t >  B( IB, copy_value );
        const auto         norm_0  = lr_normF( A, B );
        const auto         rel_eps = acc.rel_eps();
        const auto         abs_eps = acc.abs_eps();
        const uint         bsize   = std::min< uint >( std::max< uint >( CFG::BLAS::sample_size, 1 ),
                                                       std::min< uint >( n, m ) );
        const uint         nblocks = std::min( n, m ) / bsize;
        auto               Qs = std::list< Matrix< value_t > >();
        Matrix< value_t >  T_i( m, bsize );
        Matrix< value_t >  BAtQ( m, bsize );
        TRNG               rng;

        // auto  test_ortho = [&Qs,bsize] ()
        // {
        //     const uint         nQ = Qs.size();
        //     Matrix< value_t >  QtQ( nQ * bsize, nQ * bsize );
        //     uint               i = 0;

        //     for ( const auto & Qi : Qs )
        //     {
        //         uint  j = 0;
                
        //         for ( const auto & Qj : Qs )
        //         {
        //             Matrix< value_t >  QtQ_ij( QtQ,
        //                                        Range( i*bsize, (i+1)*bsize-1 ),
        //                                        Range( j*bsize, (j+1)*bsize-1 ) );

        //             prod( value_t(1), adjoint(Qi), Qj, value_t(1), QtQ_ij );
        //             ++j;
        //         }// for

        //         ++i;
        //     }// for

        //     const auto  f = norm2( QtQ );

        //     if ( f - 1.0 > 1e-12 )
        //     {
        //         DBG::printf( "|Q^H Q| - 1 = %.6e", f - 1 );
        //         DBG::write( QtQ, "QtQ.mat", "QtQ" );
        //         DBG::breakpoint();

        //         i = 0;
                
        //         for ( const auto & Qi : Qs )
        //         {
        //             auto  mname = to_string( "Q%02d", i );
        //             auto  fname = to_string( "Q%02d.mat", i );

        //             DBG::write( Qi, fname.c_str(), mname.c_str() );
        //             ++i;
        //         }// for
        //     }// if
        // };
        
        for ( uint  i = 0; i < nblocks; ++i )
        {
            fill_rand_normal( T_i );

            auto  BtT = prod( value_t(1), adjoint(B), T_i );
            auto  Q_i = prod( value_t(1), A, BtT ); // Y_i

            //
            // power iteration
            //
            
            Matrix< value_t >  R( bsize, bsize );
            Matrix< value_t >  AtQ( rank, bsize );
            
            for ( uint  j = 0; j < power_steps; ++j )
            {
                // copy( Y_i, Q_i );
                qr( Q_i, R );
                prod( value_t(1), adjoint(A), Q_i, value_t(0), AtQ );
                prod( value_t(1), B, AtQ, value_t(0), BAtQ );
                
                qr( BAtQ, R );
                prod( value_t(1), adjoint(B), BAtQ, value_t(0), AtQ );
                prod( value_t(1), A, AtQ, value_t(0), Q_i );  // Q_i = Y_i
            }// for
            
            // copy( Y_i, Q_i );
            qr( Q_i, R );
            
            //
            // project Q_i away from previous Q_j
            //
                
            if ( i > 0 )
            {
                Matrix< value_t >  C_i( Q_i, copy_value );
                Matrix< value_t >  QjtQi( bsize, bsize );
                
                for ( const auto &  Q_j : Qs )
                {
                    prod( value_t(1), adjoint(Q_j), C_i, value_t(0), QjtQi );
                    prod( value_t(-1), Q_j, QjtQi, value_t(1), Q_i );
                }// for
                
                qr( Q_i, R );
            }// if

            //
            // M = M - Q_i Q_i^T M = A·B^H - Q_i Q_i^T A·B^H = (A - Q_i Q_i^T A) B^H
            //

            auto  QtA = prod( value_t(1), adjoint(Q_i), A );

            prod( value_t(-1), Q_i, QtA, value_t(1), A );
            
            const auto  norm_i = lr_normF( A, B );

            // LOG::printf( "|A| = %.6e, %6e    %6e, %6e", norm_A, abs_eps, norm_A / norm_M, rel_eps );

            Qs.push_back( std::move( Q_i ) );

            // test_ortho();
            
            if (( norm_i < abs_eps ) || (( norm_i / norm_0 ) < rel_eps ))
                break;
        }// for

        //
        // collect Q_i's into final result
        //

        Matrix< value_t >  Q( n, Qs.size() * bsize );
        idx_t              pos = 0;

        for ( const auto &  Q_i : Qs )
        {
            Matrix< value_t >  Q_sub( Q, Range::all, Range( pos * bsize, (pos+1)*bsize - 1 ) );

            copy( Q_i, Q_sub );
            ++pos;
        }// for

        return Q;
    }// else
}


//
// approximate by randomized SVD
//
// compute column basis Q for A
// assumption:  A = Q·Q^H·A
//                = Q·B             [ Q_B, R_B ] = qr(B^H)
//                = Q·R_B^H·Q_B^H   [ U, S, V^H ]  = svd( R_B )
//                = Q·V·S·U^H·B^H
//
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, size_t >
approx_randsvd ( T1 &                              M,
                 const TTruncAcc &                 acc,
                 Matrix< typename T1::value_t > &  A,
                 Matrix< typename T1::value_t > &  B,
                 const uint                        power_steps,
                 const uint                        oversampling )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    const idx_t  n = idx_t( M.nrows() );
    const idx_t  m = idx_t( M.ncols() );

    //
    // compute column basis
    //
    
    auto  Q = column_basis( M, acc, power_steps, oversampling );
    auto  k = idx_t(Q.ncols());

    // B = Q^H · M  or B^H = M^H · Q
    auto  BT = prod( value_t(1), adjoint(M), Q );

    Matrix< value_t >  R_B( k, k );
    Matrix< value_t >  V( k, k );
    Vector< real_t >   S;

    // B^T = Q_B R_B  (Q_B overwrites B)
    qr( BT, R_B );

    // R_B = U·S·V^H
    svd( R_B, S, V );

    // determine truncated rank based on singular values
    k = idx_t( acc.trunc_rank( S ) );

    // A = Y · V_k, B = B^T · U_k
    Matrix< value_t >  Uk( R_B, Range::all, Range( 0, k-1 ) );
    Matrix< value_t >  Vk( V,   Range::all, Range( 0, k-1 ) );
    
    A = prod( value_t(1), Q,  Vk );
    B = prod( value_t(1), BT, Uk );

    if ( n < m ) prod_diag( A, S, k );
    else         prod_diag( B, S, k );

    return k;
}

//
// wrapper around different approximation functions
//
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, size_t >
approx ( T1 &                              M,
         const TTruncAcc &                 acc,
         Matrix< typename T1::value_t > &  A,
         Matrix< typename T1::value_t > &  B )
{
    LOG_COUNTER( approx );

    switch ( CFG::BLAS::approx_method )
    {
        default :
        case use_svd     : return approx_svd(  M, acc, A, B );
        case use_rrqr    : return approx_rrqr( M, acc, A, B );
        case use_randsvd : return approx_randsvd( M, acc, A, B );
    }// switch
}



//!
//! truncate given \a IA · \a IB^H low rank matrix (\a IA being n×k
//! and \a IB being m×k) with respect to given accuracy \a acc;
//! store truncated matrix in IA and IB
//!
template <typename T>
size_t
truncate_svd ( Matrix< T > &      IA,
               Matrix< T > &      IB,
               const TTruncAcc &  acc )
{
    using  value_t = T;
    using  real_t  = typename real_type< value_t >::type_t;

    if ( IA.ncols() != IB.ncols() )
        HERROR( ERR_MAT_SIZE, "(BLAS) truncate_svd", "rank in A and B differs" );

    DO_CHECK_INF_NAN( IA, "(BLAS) truncate_svd", "in input matrix A" );
    DO_CHECK_INF_NAN( IB, "(BLAS) truncate_svd", "in input matrix B" );

    const idx_t  n = idx_t( IA.nrows() );
    const idx_t  m = idx_t( IB.nrows() );

    /////////////////////////////////////////////////////////////////
    #if CHECK_TRUNCATE == 1
    Matrix< value_t >   TM( n, m );

    prod( value_t(1), IA, adjoint(IB), value_t(0), TM );
    #endif
    /////////////////////////////////////////////////////////////////

    //
    // remove zero vector pairs in A/B
    //

    const idx_t   irank = idx_t( IA.ncols() );
    idx_t         k     = 0;
    const Range   row_is( 0, n-1 );
    const Range   col_is( 0, m-1 );

    for ( idx_t  i = 0; i < irank; i++ )
    {
        Vector< value_t >  a_i( IA, row_is, i );
        Vector< value_t >  b_i( IB, col_is, i );

        if (( norm2( a_i ) > real_t(0) ) && ( norm2( b_i ) > real_t(0) ))
        {
            // copy if vectors have been removed
            if ( k < i )
            {
                Vector< value_t >  a_k( IA, row_is, k );
                Vector< value_t >  b_k( IB, col_is, k );
                
                copy( a_i, a_k );
                copy( b_i, b_k );
            }// if
                
            k++;
        }// if
    }// for

    //
    // we don't increase rank
    //

    const idx_t  acc_rank = idx_t( acc.rank() );
    
    if ( k == 0 )
    {
        // reset matrices
        IA = std::move( Matrix< value_t >( n, 0 ) );
        IB = std::move( Matrix< value_t >( m, 0 ) );
        
        return 0;
    }// if

    if ( k <= acc_rank )
    {
        // copy restricted vector pairs
        Matrix< value_t >  A( n, k );
        Matrix< value_t >  B( m, k );
        Matrix< value_t >  RA( IA, row_is, Range( 0, k-1 ) );
        Matrix< value_t >  RB( IB, col_is, Range( 0, k-1 ) );

        copy( RA, A );
        copy( RB, B );

        IA = std::move( A );
        IB = std::move( B );
        
        return k;
    }// if

    //
    // if k is bigger than the possible rank,
    // we create a dense-matrix and do truncation
    // via full SVD
    //

    const Range        irank_is( 0, k-1 );
    Matrix< value_t >  A( IA, row_is, irank_is );
    Matrix< value_t >  B( IB, col_is, irank_is );
    const idx_t        mrc   = std::min(n, m);
    idx_t              orank = 0;
        
    if ( k >= mrc )
    {
        //
        // build U = A*B^T
        //
            
        Matrix< value_t >  M( n, m );

        prod( value_t(1), A, adjoint(B), value_t(0), M );
            
        //
        // truncate to rank-k
        //

        TTruncAcc  lacc( acc );

        lacc.set_max_rank( k );

        orank = idx_t( approx( M, lacc, IA, IB ) ); // approx _reallocates_ IA/IB
    }// if
    else
    {
        //
        // do QR-factorisation of A and B
        //

        Matrix< value_t >  QA, QB, RA, RB;

        QA = std::move( Matrix< value_t >( A.nrows(), k ) );
        RA = std::move( Matrix< value_t >( k, k ) );
        
        copy( A, QA );
        qr( QA, RA );
        
        QB = std::move( Matrix< value_t >( B.nrows(), k ) );
        RB = std::move( Matrix< value_t >( k, k ) );
        
        copy( B, QB );
        qr( QB, RB );

        //
        // R = R_A · upper_triangular(QB)^H = R_B^H
        //
        
        Matrix< value_t >  R( k, k );

        prod( value_t(1), RA, adjoint(RB), value_t(0), R );
        // BLAS::trmm( 'R', 'U', 'C', 'N', k, k, 1.0, RT, k, R, k );
        
        //
        // SVD(R) = U S V^H
        //
            
        Vector< real_t >   S( k );
        Matrix< value_t >  U( std::move( R ) );
        Matrix< value_t >  V( std::move( RB ) );
            
        svd( U, S, V );
        
        // determine truncated rank based on singular values
        orank = idx_t( acc.trunc_rank( S ) );

        //
        // only build new vectors, if rank is decreased
        //
        
        if ( orank < irank )
        {
            //
            // build new matrices A and B
            //

            const Range  orank_is( 0, orank-1 );

            LOG_COUNTER( gemm );
            LOG_TIC( gemm );
                
            // A := Q_A · U
            Matrix< value_t >  Urank( U, irank_is, orank_is );
            
            // U := U·S
            prod_diag( Urank, S, orank );
            
            IA = std::move( Matrix< value_t >( n, orank ) );
            prod( value_t(1), QA, Urank, value_t(0), IA );

            
            // B := Q_B · conj(V)
            Matrix< value_t >  Vrank( V, irank_is, orank_is );

            // conj( Vrank );
            IB = std::move( Matrix< value_t >( m, orank ) );
            prod( value_t(1), QB, Vrank, value_t(0), IB );

            LOG_TOC( gemm );
        }// if
    }// else
            
    /////////////////////////////////////////////////////////////////
    #if CHECK_TRUNCATE == 1
    const real_t  norm_M = normF( TM );
        
    prod( value_t(-1), IA, adjoint(IB), value_t(1), TM );
    DBG::printf( "(BLAS) truncate_svd : |M-AB^H|/|M| = %.6e", normF( TM ) / norm_M );
    #endif
    /////////////////////////////////////////////////////////////////
    
    DO_CHECK_INF_NAN( IA, "(BLAS) truncate_svd", "in output matrix A" );
    DO_CHECK_INF_NAN( IB, "(BLAS) truncate_svd", "in output matrix B" );

    return orank;
}

template <typename T>
std::pair< Matrix< T >, Matrix< T > >
truncate2_svd  ( const Matrix< T > &  A,
                 const Matrix< T > &  B,
                 const TTruncAcc &    acc )
{
    using  value_t = T;
    using  real_t  = typename real_type< value_t >::type_t;

    if ( A.ncols() != B.ncols() )
        HERROR( ERR_MAT_SIZE, "(BLAS) truncate_svd", "rank in A and B differs" );

    const idx_t  n     = idx_t( A.nrows() );
    const idx_t  m     = idx_t( B.nrows() );
    const idx_t  irank = idx_t( A.ncols() );

    //
    // don't increase rank
    //

    const idx_t  acc_rank = idx_t( acc.rank() );

    Matrix< T >  OA, OB;
    
    if ( irank == 0 )
    {
        // reset matrices
        OA = std::move( Matrix< value_t >( n, 0 ) );
        OB = std::move( Matrix< value_t >( m, 0 ) );
    }// if

    if ( irank <= acc_rank )
    {
        OA = std::move( Matrix< value_t >( A, copy_value ) );
        OB = std::move( Matrix< value_t >( B, copy_value ) );
    }// if

    //
    // if k is bigger than the possible rank,
    // we create a dense-matrix and do truncation
    // via full SVD
    //

    const idx_t  mrc   = std::min(n, m);
    idx_t        orank = 0;
        
    if ( acc_rank >= mrc )
    {
        //
        // build U = A*B^T
        //
            
        Matrix< value_t >  M( n, m );

        prod( value_t(1), A, adjoint(B), value_t(0), M );
            
        //
        // truncate to rank-k
        //

        TTruncAcc  lacc( acc );

        lacc.set_max_rank( acc_rank );

        std::tie( OA, OB ) = approx_svd( M, lacc );
    }// if
    else
    {
        //
        // do QR-factorisation of A and B
        //

        Matrix< value_t >  QA, QB, RA, RB;

        QA = std::move( Matrix< value_t >( A.nrows(), irank ) );
        RA = std::move( Matrix< value_t >( irank, irank ) );
        
        copy( A, QA );
        qr( QA, RA );
        
        QB = std::move( Matrix< value_t >( B.nrows(), irank ) );
        RB = std::move( Matrix< value_t >( irank, irank ) );
        
        copy( B, QB );
        qr( QB, RB );

        //
        // R = R_A · upper_triangular(QB)^H = R_B^H
        //
        
        Matrix< value_t >  R( irank, irank );

        prod( value_t(1), RA, adjoint(RB), value_t(0), R );
        
        //
        // SVD(R) = U S V^H
        //
            
        Vector< real_t >   S( irank );
        Matrix< value_t >  U( std::move( R ) );
        Matrix< value_t >  V( std::move( RB ) );
            
        svd( U, S, V );
        
        // determine truncated rank based on singular values
        orank = idx_t( acc.trunc_rank( S ) );

        //
        // only build new vectors, if rank is decreased
        //
        
        if ( orank < irank )
        {
            //
            // build new matrices A and B
            //

            const Range  irank_is( 0, irank-1 );
            const Range  orank_is( 0, orank-1 );

            LOG_COUNTER( gemm );
            LOG_TIC( gemm );
                
            // A := Q_A · U
            Matrix< value_t >  Urank( U, irank_is, orank_is );
            
            // U := U·S
            prod_diag( Urank, S, orank );
            OA = prod( value_t(1), QA, Urank );
            
            // B := Q_B · conj(V)
            Matrix< value_t >  Vrank( V, irank_is, orank_is );

            OB = prod( value_t(1), QB, Vrank );

            LOG_TOC( gemm );
        }// if
        else
        {
            OA = std::move( Matrix< value_t >( A, copy_value ) );
            OB = std::move( Matrix< value_t >( B, copy_value ) );
        }// else
    }// else

    return { std::move( OA ), std::move( OB ) };
}

namespace
{

//
// copy row/column <j> to row/column <i> in matrix M
//
template <typename T>
void
copy ( Matrix< T > & M,
       const idx_t   i,
       const idx_t   j )
{
    Vector< T >  col_i( M.column( i ) );
    Vector< T >  col_j( M.column( j ) );
            
    copy( col_j, col_i );
}

//
// copy row/column <i> to vector <tmp>
//
template <typename T>
void
copy ( Matrix< T > &  M,
       const idx_t    i,
       Vector< T > &  tmp )
{
    Vector< T >  col_i( M.column( i ) );
            
    copy( col_i, tmp );
}

//
// copy vector <tmp> to row/column <i>
//
template <typename T>
void
copy ( Vector< T > &  tmp,
       const idx_t    i,
       Matrix< T > &  M )
{
    Vector< T >  col_i( M.column( i ) );
            
    copy( tmp, col_i );
}

//
// swap the rows/columns <i> and <j> in matrix M
//
template <typename T>
void
swap ( Matrix< T > &  M,
       const idx_t    i,
       const idx_t    j,
       Vector< T > &  tmp )
{
    Vector< T >  col_i( M.column( i ) );
    Vector< T >  col_j( M.column( j ) );
    
    copy( col_i, tmp );
    copy( col_j, col_i );
    copy( tmp,   col_j );
}

//
// sort rows/columns of matrix according to given permutation
//
template <typename T>
void
matrix_sort( Matrix< T > &                M,
             std::vector< blas_int_t > &  P,
             const idx_t                  lb,
             const idx_t                  ub,
             Vector< T > &                tmp )
{
    if ( lb >= ub ) return;

    if ( (ub - lb) < 20 )
    {
        //
        // apply insertion sort for small ranges
        //

        for ( idx_t  i = lb+1; i <= ub; i++ )
        {
            const idx_t  v = P[i];
            idx_t        j = i-1;

            copy( M, i, tmp );
            
            while (( j >= 0 ) && ( P[j] > v ))
            {
                copy( M, j+1, j );
                P[j+1] = P[j];
                j--;
            }// if

            copy( tmp, j+1, M );
            P[j+1] = v;
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
        idx_t        choice[3] = { P[lb], P[mid], P[ub] };
        idx_t        pivot;

        // choose pivot (median-of-three)
        if ( choice[0] > choice[1] ) std::swap( choice[0], choice[1] );
        if ( choice[0] > choice[2] ) std::swap( choice[0], choice[2] );
        if ( choice[1] > choice[2] ) std::swap( choice[1], choice[2] );
        pivot = choice[1];

        // partition
        while ( i < j )
        {
            while ( P[i] < pivot   ) i++;
            while ( pivot   < P[j] ) j--;

            if ( i < j )
            {
                swap( M, i, j, tmp );
                std::swap( P[i], P[j] );
            }// if
        }// while

        // recursion
        matrix_sort( M, P, lb, i-1, tmp );
        matrix_sort( M, P, i+1, ub, tmp );
    }// else
}

}// namespace anonymous

//
// permute matrix columns
//
template <typename T>
void
permute_col ( Matrix< T > &                M,
              std::vector< blas_int_t > &  P )
{
    Vector< T >                 tmp( M.nrows() );
    Hpro::vector< blas_int_t >  Pt( P.size() );
        
    for ( size_t  i = 0; i < P.size(); ++i )
        Pt[ P[i] ] = i;
            
    matrix_sort( M, Pt, 0, M.ncols()-1, tmp );
}

//!
//! truncate given \a IA · \a IB^H low rank matrix (\a IA being n×k
//! and \a IB being m×k) with respect to given accuracy \a acc;
//! store truncated matrix in IA and IB
//!
template <typename T>
size_t
truncate_rrqr ( Matrix< T > &      A,
                Matrix< T > &      B,
                const TTruncAcc &  acc )
{
    using  value_t = T;

    if ( A.ncols() != B.ncols() )
        HERROR( ERR_MAT_SIZE, "(BLAS) truncate_rrqr", "rank in A and B differs" );

    DO_CHECK_INF_NAN( A, "(BLAS) truncate_rrqr", "in input matrix A" );
    DO_CHECK_INF_NAN( B, "(BLAS) truncate_rrqr", "in input matrix B" );

    const idx_t  n = idx_t( A.nrows() );
    const idx_t  m = idx_t( B.nrows() );
    
    /////////////////////////////////////////////////////////////////
    #if CHECK_TRUNCATE == 1
    Matrix< value_t >   TM( n, m );

    prod( value_t(1), A, adjoint(B), value_t(0), TM );
    #endif
    /////////////////////////////////////////////////////////////////

    //
    // we don't increase rank
    //

    const idx_t  k_old    = idx_t( A.ncols() );
    const idx_t  acc_rank = idx_t( acc.rank() );
    
    if ( k_old == 0 )
    {
        // reset matrices
        A = std::move( Matrix< value_t >( n, 0 ) );
        B = std::move( Matrix< value_t >( m, 0 ) );
        
        return 0;
    }// if

    if ( k_old <= acc_rank )
        return k_old;

    //
    // if k is bigger than the possible rank,
    // we create a dense-matrix and do truncation
    // via full SVD
    //

    const idx_t  mrc   = std::min(n, m);
    idx_t        k_new = 0;
        
    if ( k_old >= mrc )
    {
        //
        // build U = A·B^T
        //
            
        Matrix< value_t >  M( n, m );

        prod( value_t(1), A, adjoint(B), value_t(0), M );
            
        //
        // truncate to rank-k
        //

        TTruncAcc  lacc( acc );

        lacc.set_max_rank( k_old );

        k_new = idx_t( approx( M, lacc, A, B ) ); // "approx" _reallocates_ A/B
    }// if
    else
    {
        //
        // QR-factorisation of B
        //

        Matrix< value_t >  QB( B );
        Matrix< value_t >  RB( k_old, k_old );

        qr( QB, RB );

        //
        // QA = A·RB'
        //

        Matrix< value_t >  QA( A.nrows(), k_old );

        prod( T(1), A, adjoint(RB), T(0), QA );
        
        //
        // compute column-pivoted QR of A
        //

        Matrix< value_t >          RA( k_old, k_old );
        std::vector< blas_int_t >  P( k_old, 0 );

        //
        // do partial QRP based on given accuracy
        //
        
        k_new = qrp_trunc( QA, RA, P, acc );
        
        //
        // restrict first k_new columns
        //

        Matrix< value_t >  Qk( QA, Range::all,          Range( 0, k_new-1 ) );
        Matrix< value_t >  Rk( RA, Range( 0, k_new-1 ), Range( 0, k_old-1 ) );

        // A = QA · RA
        A = std::move( Matrix< value_t >( n, k_new ) );
        copy( Qk, A );

        // B = QB · P  (B' = P' · QB')

        #if 1

        // permute QB by copy (more memory but thread safe)
        Matrix< value_t >  QB_P( m, k_old );
        
        for ( int i = 0; i < k_old; ++i )
        {
            auto  j      = P[i];
            auto  QB_P_i = QB_P.column( i );
            auto  Q_j    = QB.column( j );

            copy( Q_j, QB_P_i );
        }// for

        B = std::move( Matrix< value_t >( m, k_new ) );

        prod( value_t(1), QB_P, adjoint( Rk ), value_t(1), B );
        
        #else

        // permute QB by sorting (less memory but not thread safe for unknown reasons)
        permute_col( QB, P );

        B = std::move( Matrix< value_t >( m, k_new ) );

        prod( value_t(1), QB, adjoint( Rk ), value_t(1), B );

        #endif
        
    }// else
    
    /////////////////////////////////////////////////////////////////
    #if CHECK_TRUNCATE == 1
    using real_t  = typename real_type< value_t >::type_t;
    const real_t  norm_M = normF( TM );
        
    prod( value_t(-1), A, adjoint(B), value_t(1), TM );
    DBG::printf( "(BLAS) truncate_rrqr : |M-AB^H|/|M| = %.6e", normF( TM ) / norm_M );

    // if ( normF( TM ) / norm_M > 1e-3 )
    //     std::exit( 0 );
    #endif
    /////////////////////////////////////////////////////////////////

    
    DO_CHECK_INF_NAN( A, "(BLAS) truncate_rrqr", "in output matrix A" );
    DO_CHECK_INF_NAN( B, "(BLAS) truncate_rrqr", "in output matrix B" );

    return k_new;
}

//!
//! truncate given \a A · \a B^H low rank matrix (\a A being n×k
//! and \a B being m×k) with respect to given accuracy \a acc;
//! store truncated matrix in A and B
//!
template <typename T>
size_t
truncate_rand ( Matrix< T > &      A,
                Matrix< T > &      B,
                const TTruncAcc &  acc )
{
    using  value_t = T;
    using  real_t  = typename real_type< value_t >::type_t;

    if ( A.ncols() != B.ncols() )
        HERROR( ERR_MAT_SIZE, "(BLAS) truncate_rand", "rank in A and B differs" );

    DO_CHECK_INF_NAN( A, "(BLAS) truncate_rand", "in input matrix A" );
    DO_CHECK_INF_NAN( B, "(BLAS) truncate_rand", "in input matrix B" );

    const idx_t  n = idx_t( A.nrows() );
    const idx_t  m = idx_t( B.nrows() );

    /////////////////////////////////////////////////////////////////
    #if CHECK_TRUNCATE == 1
    Matrix< value_t >   TM( n, m );
    Matrix< value_t >   TA( A, copy_value );
    Matrix< value_t >   TB( B, copy_value );

    prod( value_t(1), A, adjoint(B), value_t(0), TM );
    #endif
    /////////////////////////////////////////////////////////////////

    //
    // we don't increase rank
    //

    const idx_t  k_old    = idx_t( A.ncols() );
    const idx_t  acc_rank = idx_t( acc.rank() );
    
    if ( k_old == 0 )
    {
        // reset matrices
        A = std::move( Matrix< value_t >( n, 0 ) );
        B = std::move( Matrix< value_t >( m, 0 ) );
        
        return 0;
    }// if

    if ( k_old <= acc_rank )
        return k_old;

    //
    // if k is bigger than the possible rank,
    // we create a dense-matrix and do truncation
    // via full SVD
    //

    const idx_t  mrc   = std::min(n, m);
    idx_t        k_new = 0;
        
    if ( k_old >= mrc )
    {
        //
        // build U = A·B^T
        //
            
        Matrix< value_t >  M( n, m );

        prod( value_t(1), A, adjoint(B), value_t(0), M );
            
        //
        // truncate to rank-k
        //

        TTruncAcc  lacc( acc );

        lacc.set_max_rank( k_old );

        k_new = idx_t( approx( M, lacc, A, B ) ); // "approx" _reallocates_ A/B
    }// if
    else
    {
        //
        // compute column basis
        //

        auto  Q      = column_basis( A, B, acc, CFG::BLAS::power_steps, CFG::BLAS::oversampling );
        auto  k_base = idx_t(Q.ncols());

        // Q^H · A · B^H  = (B·A^H·Q)^H
        auto  AtQ  = prod( value_t(1), adjoint(A), Q );
        auto  BAtQ = prod( value_t(1), B, AtQ );

        Matrix< value_t >  R( k_base, k_base );
        Matrix< value_t >  V( k_base, k_base );
        Vector< real_t >   S;

        // (B·A^H·Q)^H = Q_B R
        qr( BAtQ, R );
        
        // R_B = U·S·V^H
        svd( R, S, V );
        
        // determine truncated rank based on singular values
        k_new = idx_t( acc.trunc_rank( S ) );

        // A = Y · V_k, B = B^T · U_k
        Matrix< value_t >  Uk( R, Range::all, Range( 0, k_new-1 ) );
        Matrix< value_t >  Vk( V, Range::all, Range( 0, k_new-1 ) );

        A = prod( value_t(1), Q,    Vk );
        B = prod( value_t(1), BAtQ, Uk );

        if ( n < m ) prod_diag( A, S, k_new );
        else         prod_diag( B, S, k_new );
    }// else
    
    /////////////////////////////////////////////////////////////////
    #if CHECK_TRUNCATE == 1
    Matrix< value_t >   TMC( TM, copy_value );
    const real_t  norm_M = normF( TM );
        
    prod( value_t(-1), A, adjoint(B), value_t(1), TM );
    DBG::printf( "(BLAS) truncate_rand : |M-AB^H|/|M| = %.6e", normF( TM ) / norm_M );
    #endif
    /////////////////////////////////////////////////////////////////

    
    DO_CHECK_INF_NAN( A, "(BLAS) truncate_rand", "in output matrix A" );
    DO_CHECK_INF_NAN( B, "(BLAS) truncate_rand", "in output matrix B" );

    return k_new;
}

template <typename T>
size_t
truncate ( Matrix< T > &      IA,
           Matrix< T > &      IB,
           const TTruncAcc &  acc )
{
    LOG_COUNTER( trunc );
    
    switch ( CFG::BLAS::trunc_method )
    {
        default:
        case use_svd  : return truncate_svd(  IA, IB, acc );
        case use_rrqr : return truncate_rrqr( IA, IB, acc );
        case use_rand : return truncate_rand( IA, IB, acc );
    }// switch
}

//
// construct factorisation A = Q·R of \a A, with orthonormal Q
// - A is overwritten with Q upon exit
//
template < typename T1 >
std::enable_if_t< is_matrix< T1 >::value, void >
factorise_ortho ( T1 &                              A,
                  Matrix< typename T1::value_t > &  R )
{
    using  value_t = typename T1::value_t;

    if ( std::min( A.nrows(), A.ncols() ) == 0 )
    {
        R = std::move( Matrix< value_t >() );
        return;
    }// if

    //
    // compute QR factorisation of A
    //
    
    qr( A, R );
}

//
// construct approximate factorisation A = Q·R of \a A, with orthonormal Q
// - approximation quality is definded by \a acc
// - A is overwritten with Q upon exit
//
template < typename T >
void
factorise_ortho ( Matrix< T > &      A,
                  Matrix< T > &      R,
                  const TTruncAcc &  acc )
{
    using  value_t = T;
    using  real_t  = typename real_type< value_t >::type_t;
    
    const size_t       n = A.nrows();
    const size_t       m = A.ncols();
    const size_t       q = std::min( n, m );
    Matrix< value_t >  Q_hat;
    Vector< real_t >   S( q );

    if ( q == 0 )
    {
        R = std::move( Matrix< value_t >() );
        return;
    }// if

    #if 1
    
    //
    // 1. compute QR factorisation of A
    // 2. compute SVD of R and restrict Q according
    //    to singular values and given accuracy
    // 3. R ≔ Q^H A 
    //
    
    if ( n >= m )
    {
        Matrix< value_t >  V( m, m );
        Matrix< value_t >  Q( n, m );

        copy( A, Q );
        
        // compute QR of A
        qr( Q, R );

        // compute SVD of R
        svd( R, S, V );

        // use Q·U, e.g. Q·R, since X = QUSV^H
        Q_hat = std::move( Matrix< value_t >( n, m ) );
        prod( value_t(1), Q, R, value_t(0), Q_hat );
    }// if
    else
    {
        Matrix< value_t >  V( n, n );
        
        // compute QR of A^H
        Q_hat = std::move( Matrix< value_t >( m, n ) );
        copy( adjoint( A ), Q_hat );
        qr( Q_hat, R );

        // compute SVD of R^H
        conj_transpose( R );
        svd( R, S, V );

        // use U, e.g. R, since X = USV^HQ^H
        Q_hat = std::move( R );
    }// else

    // compute new rank and assign final matrices
    size_t  rank = acc.trunc_rank( S );

    Matrix< value_t >  Q( n, rank );
    Matrix< value_t >  Q_hat_l( Q_hat, Range( 0, idx_t(n)-1 ), Range( 0, idx_t(rank)-1 ) );

    copy( Q_hat_l, Q );

    #else

    //
    // 1. compute SVD of A = U·S·V^H, restricted to left singular vectors
    // 2. assign first vectors (according to accuracy ) of U to Q
    // 3. R ≔ Q^H A 
    //

    // compute left singular vectors and values
    Matrix< value_t >  U( n, m );

    copy( A, U );
    svd( U, S );
    
    // compute new rank
    size_t  rank = acc.trunc_rank( S );
    
    Matrix< value_t >  Q( n, rank );
    Matrix< value_t >  U_l( U, Range( 0, n-1 ), Range( 0, rank-1 ) );
    
    copy( U_l, Q );
    
    #endif
    
    R = std::move( Matrix< value_t >( rank, m ) );
    prod( value_t(1), adjoint( Q ), A, value_t(0), R );
    A = std::move( Q );
}

//
// print statistics for Algebra functions
//
void
print_statistics ()
{
    #ifdef HAS_LOGGING
    
    size_t  napprox = LOGGING::napprox;
    size_t  ntrunc  = LOGGING::ntrunc;
    size_t  ngesvd  = LOGGING::ngesvd;
    size_t  ngesdd  = LOGGING::ngesdd;
    size_t  ngesvj  = LOGGING::ngesvj;
    size_t  nqr     = LOGGING::nqr;
    size_t  ngemm   = LOGGING::ngemm;

    if ( napprox > 0 ) LOG::printf( "approx: #%d ",           napprox );
    if ( ntrunc  > 0 ) LOG::printf( "trunc:  #%d ",           ntrunc );
    if ( ngesvd  > 0 ) LOG::printf( "gesvd:  #%d / t = %.2f", ngesvd, LOGGING::tgesvd );
    if ( ngesdd  > 0 ) LOG::printf( "gesdd:  #%d / t = %.2f", ngesdd, LOGGING::tgesdd );
    if ( ngesvj  > 0 ) LOG::printf( "gesvj:  #%d / t = %.2f", ngesvj, LOGGING::tgesvj );
    if ( nqr     > 0 ) LOG::printf( "qr:     #%d / t = %.2f", nqr,    LOGGING::tqr );
    if ( ngemm   > 0 ) LOG::printf( "gemm:   #%d / t = %.2f", ngemm,  LOGGING::tgemm );

    #endif
}

//
// reset statistics for Algebra functions
//
void
reset_statistics ()
{
    #ifdef HAS_LOGGING
    
    LOGGING::napprox = 0;
    LOGGING::ntrunc  = 0;
    LOGGING::ngesvd  = 0;
    LOGGING::ngesdd  = 0;
    LOGGING::ngesvj  = 0;
    LOGGING::nqr     = 0;
    LOGGING::ngemm   = 0;

    LOGGING::tapprox = 0;
    LOGGING::ttrunc  = 0;
    LOGGING::tgesvd  = 0;
    LOGGING::tgesdd  = 0;
    LOGGING::tgesvj  = 0;
    LOGGING::tqr     = 0;
    LOGGING::tgemm   = 0;

    #endif
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//
// include template instantiations
//
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

#include "Algebra_Instantiations.icc"

}// namespace BLAS

}// namespace Hpro
