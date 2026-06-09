//
// Project     : HLIBpro
// File        : blas/cuda.cc
// Description : CUDA related functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2024. All Rights Reserved.
//

#include <array>
#include <mutex>

#include <hpro/config.h>

#if HPRO_USE_CUDA == 1
#  include <cuda_runtime.h>
#  include <cublas_v2.h>
#  include <cusolverDn.h>
#endif

#include <hpro/blas/Algebra.hh>
#include <hpro/blas/cuda.hh>

namespace Hpro { namespace CUDA {

////////////////////////////////////////////////////////////////////////////////
//
// management functions
//
////////////////////////////////////////////////////////////////////////////////

#if HPRO_USE_CUDA == 1

// wrapper for cuda, cuBlas and cuSolver functions
#define HPRO_CUDA_CHECK( func, args )                                   \
    {                                                                   \
        auto  result = func args ;                                      \
        if ( result != cudaSuccess )                                   \
        {                                                               \
            HERROR( ERR_CUDA, #func, Hpro::to_string( "CUDA code: %d", result ) ); \
        }                                                               \
    }

#define HPRO_CUBLAS_CHECK( func, args )                                 \
    {                                                                   \
        auto  result = func args ;                                      \
        if ( result != CUBLAS_STATUS_SUCCESS )                          \
        {                                                               \
            HERROR( ERR_CUDA, #func, Hpro::to_string( "CUBLAS code: %d", result ) ); \
        }                                                               \
    }

#define HPRO_CUSOLVER_CHECK( func, args )                               \
    {                                                                   \
        auto  result = func args ;                                      \
        if ( result != CUSOLVER_STATUS_SUCCESS )                        \
        {                                                               \
            HERROR( ERR_CUDA, #func, Hpro::to_string( "CUSOLVER code: %d", result ) ); \
        }                                                               \
    }

//
// mapping of default types to cuBLAS types
//
template < typename T > struct cuda_traits
{
    using  cuda_type = T;
    using  ptr_type  = T *;
    
    static constexpr cudaDataType  blas_type = CUDA_R_32F;
};

template <>
struct cuda_traits< float >
{
    using  cuda_type = float;
    using  ptr_type  = float *;
    
    static constexpr cudaDataType  blas_type = CUDA_R_32F;
};

template <>
struct cuda_traits< double >
{
    using  cuda_type = double;
    using  ptr_type  = double *;
    
    static constexpr cudaDataType  blas_type = CUDA_R_64F;
};

template <>
struct cuda_traits< std::complex< float > >
{
    using  cuda_type = cuFloatComplex;
    using  ptr_type  = cuFloatComplex *;
    
    static constexpr cudaDataType  blas_type = CUDA_C_32F;
};

template <>
struct cuda_traits< std::complex< double > >
{
    using  cuda_type = cuDoubleComplex;
    using  ptr_type  = cuDoubleComplex *;
    
    static constexpr cudaDataType  blas_type = CUDA_C_64F;
};

//
// joined handle for cuBLAS and cuSolverDn
//
struct handle_t
{
    cudaStream_t        stream;
    cublasHandle_t      blas;
    cusolverDnHandle_t  solver;
};

#endif

namespace
{

#if HPRO_USE_CUDA == 1

//
// signals availability of CUDA
//
bool  has_cuda = false;

//
// concurrent streams
//
constexpr int                        max_handles = 16;
std::array< handle_t, max_handles >  handles;
std::array< bool, max_handles >      handle_used;
std::mutex                           handle_mtx;

#endif

}// namespace anonymous

//
// initialize CUDA related data structures
//
void
init ()
{
    #if HPRO_USE_CUDA == 1

    //
    // first check number of devices
    //
    
    int  ndev = 0;
    
    if ( cudaGetDeviceCount( & ndev ) != cudaSuccess )
        return;

    HINFO( to_string( "(CUDA) init : found %d CUDA device(s)", ndev ) );
    
    if ( ndev == 0 )
        return;

    //
    // now try to set up streams
    //
    
    bool  found_error = false;
    uint  error_index = 0;
    
    for ( uint  i = 0; i < max_handles; ++i )
    {
        if ( cudaStreamCreateWithFlags( & handles[i].stream, cudaStreamNonBlocking ) != cudaSuccess )
        {
            HWARNING( "(CUDA) init : error in cudaStreamCreateWithFlags; disabling CUDA" );
            found_error = true;
            error_index = i;
            break;
        }// if

        // after this, assumption is that CUDA is present and working!!!
        
        HPRO_CUBLAS_CHECK(   cublasCreate,        ( & handles[i].blas ) );
        HPRO_CUBLAS_CHECK(   cublasSetStream,     (   handles[i].blas, handles[i].stream ) );
        
        HPRO_CUSOLVER_CHECK( cusolverDnCreate,    ( & handles[i].solver ) );
        HPRO_CUSOLVER_CHECK( cusolverDnSetStream, (   handles[i].solver, handles[i].stream ) );

        handle_used[i] = false;
    }// for

    if ( found_error )
    {
        // clean up all initialized data
        for ( uint  i = 0; i < error_index; ++i )
        {
            HPRO_CUSOLVER_CHECK( cusolverDnDestroy,  ( handles[i].solver ) );
            HPRO_CUBLAS_CHECK(   cublasDestroy,      ( handles[i].blas ) );
            HPRO_CUDA_CHECK(     cudaStreamDestroy,  ( handles[i].stream ) );
        }// for
    }// if
    else
    {
        has_cuda = true;
    }// else

    #endif
}

//
// finish CUDA related data structures
//
void
done ()
{
    #if HPRO_USE_CUDA == 1

    if ( ! has_cuda )
        return;

    has_cuda = false;
    
    for ( uint  i = 0; i < max_handles; ++i )
    {
        if ( handle_used[i] )
            HERROR( ERR_CUDA, "(CUDA) done", "handle still in use" );
        
        HPRO_CUSOLVER_CHECK( cusolverDnDestroy,  ( handles[i].solver ) );
        HPRO_CUBLAS_CHECK(   cublasDestroy,      ( handles[i].blas ) );
        HPRO_CUDA_CHECK(     cudaStreamDestroy,  ( handles[i].stream ) );
    }// for

    #endif
}

//
// get free handle for CUDA (stream/cuBlas/cuSolver)
//
cuda_handle_t
request_handle ()
{
    #if HPRO_USE_CUDA == 1
    
    if ( ! has_cuda )
        return INVALID_HANDLE;
    
    auto  lock = std::scoped_lock( handle_mtx );

    for ( uint  i = 0; i < max_handles; ++i )
    {
        if ( ! handle_used[i] )
        {
            handle_used[i] = true;
            return i;
        }// if
    }// for

    #endif
    
    return INVALID_HANDLE;
}

//
// release previously aquired handle
//
#if HPRO_USE_CUDA == 1
void
release_handle ( const cuda_handle_t  handle )
{
    if ( handle == INVALID_HANDLE )
        return;

    auto  lock = std::scoped_lock( handle_mtx );

    if (( handle < 0 ) || ( handle >= max_handles ))
        HERROR( ERR_CUDA, "(CUDA) release_handle", "unknown handle given" );

    if ( ! handle_used[ handle ] )
        HERROR( ERR_CUDA, "(CUDA) release_handle", "given handle not in use" );

    handle_used[ handle ] = false;

}
#else
void release_handle ( const cuda_handle_t  ) {}
#endif

////////////////////////////////////////////////////////////////////////////////
//
// device handling
//
////////////////////////////////////////////////////////////////////////////////

#if HPRO_USE_CUDA == 1

namespace
{

//
// allocate device memory
//
template < typename value_t >
value_t *
device_alloc ( const size_t  n )
{
    value_t *  ptr    = nullptr;
    auto       retval = cudaMalloc( & ptr, n * sizeof(value_t) );

    if ( retval != cudaSuccess )
    {
        HWARNING( "(CUDA) device_alloc : error in \"cudaMalloc\" : " + std::string( cudaGetErrorString( retval ) ) );
        return nullptr;
    }// if

    return  ptr;
}

//
// deallocate device memory
//
template < typename value_t >
void
device_free ( value_t *  ptr )
{
    const auto  retval = cudaFree( ptr );

    if ( retval != cudaSuccess )
        HWARNING( "(CUDA) device_alloc : error in \"cudaFree\" : " + std::string( cudaGetErrorString( retval ) ) );
}

//
// copy data to device
//
template < typename value_t >
bool
to_device ( handle_t                                   handle,
            const BLAS::Matrix< value_t > &            M_host,
            typename cuda_traits< value_t >::ptr_type  M_dev )
{
    const auto  retval = cudaMemcpyAsync( M_dev, M_host.data(), sizeof(value_t) * M_host.nrows() * M_host.ncols(), cudaMemcpyHostToDevice, handle.stream );
    
    if ( retval != cudaSuccess )
    {
        HWARNING( to_string( "(CUDA) from_device : error in \"cublasSetMatrixAsync\" (%d)", retval ) );
        return false;
    }// if
    else
        return true;
}

//
// copy data from device
//
template < typename value_t >
bool
from_device ( handle_t                                   handle,
              typename cuda_traits< value_t >::ptr_type  M_dev,
              BLAS::Matrix< value_t > &                  M_host )
{
    const auto  retval = cudaMemcpyAsync( M_host.data(), M_dev, sizeof(value_t) * M_host.nrows() * M_host.ncols(), cudaMemcpyDeviceToHost, handle.stream );

    if ( retval != cudaSuccess )
    {
        HWARNING( to_string( "(CUDA) from_device : error in \"cublasGetMatrixAsync\" (%d)", retval ) );
        return false;
    }// if
    else
        return true;
}

template < typename value_t >
bool
from_device ( handle_t                                   handle,
              typename cuda_traits< value_t >::ptr_type  v_dev,
              BLAS::Vector< value_t > &                  v_host )
{
    const auto  retval = cudaMemcpyAsync( v_host.data(), v_dev, sizeof(value_t) * v_host.length(), cudaMemcpyDeviceToHost, handle.stream );

    if ( retval != cudaSuccess )
    {
        HWARNING( to_string( "(CUDA) from_device : error in \"cublasGetVectorAsync\" (%d)", retval ) );
        return false;
    }// if
    else
        return true;
}

template < typename value_t >
value_t
from_device ( handle_t                                   handle,
              typename cuda_traits< value_t >::ptr_type  dev_data )
{
    value_t  data;

    const auto  retval = cudaMemcpyAsync( & data, dev_data, sizeof(typename cuda_traits< value_t >::cuda_type), cudaMemcpyDeviceToHost, handle.stream );

    if ( retval != cudaSuccess )
        HWARNING( to_string( "(CUDA) from_device : error in \"cudaMemcpAsync\" (%d)", retval ) );
    
    return data;
}

}// namespace anonymous

#endif

////////////////////////////////////////////////////////////////////////////////
//
// QR related functions
//
////////////////////////////////////////////////////////////////////////////////

#if HPRO_USE_CUDA == 1

namespace
{

#define HPRO_CUDA_ORGQR_BUFSIZE( vtype, func )                      \
    cusolverStatus_t                                                \
    orgqr_buffersize (                                              \
        cusolverDnHandle_t  handle,                                 \
        int                 nrows,                                  \
        int                 ncols,                                  \
        int                 nelem,                                  \
        const vtype *       A,                                      \
        int                 ldA,                                    \
        const vtype *       tau,                                    \
        int *               lwork)                                  \
{ return func( handle, nrows, ncols, nelem, A, ldA, tau, lwork ); }

HPRO_CUDA_ORGQR_BUFSIZE( float,           cusolverDnSorgqr_bufferSize )
HPRO_CUDA_ORGQR_BUFSIZE( double,          cusolverDnDorgqr_bufferSize )
HPRO_CUDA_ORGQR_BUFSIZE( cuComplex,       cusolverDnCungqr_bufferSize )
HPRO_CUDA_ORGQR_BUFSIZE( cuDoubleComplex, cusolverDnZungqr_bufferSize )

#undef HPRO_CUDA_ORGQR_BUFSIZE

#define HPRO_CUDA_ORGQR( vtype, func )                                  \
    cusolverStatus_t                                                    \
    orgqr ( cusolverDnHandle_t  handle,                                 \
            int                 nrows,                                  \
            int                 ncols,                                  \
            int                 nelem,                                  \
            vtype *             A,                                      \
            int                 ldA,                                    \
            const vtype *       tau,                                    \
            vtype *             work,                                   \
            int                 lwork,                                  \
            int *               devInfo)                                \
    { return func( handle, nrows, ncols, nelem, A, ldA, tau, work, lwork, devInfo ); }

HPRO_CUDA_ORGQR( float,           cusolverDnSorgqr )
HPRO_CUDA_ORGQR( double,          cusolverDnDorgqr )
HPRO_CUDA_ORGQR( cuComplex,       cusolverDnCungqr )
HPRO_CUDA_ORGQR( cuDoubleComplex, cusolverDnZungqr )

#undef HPRO_CUDA_ORGQR

}// namespace anonymous

template < typename value_t >
bool
qr  ( const cuda_handle_t        handle_idx,
      BLAS::Matrix< value_t > &  M,
      BLAS::Matrix< value_t > &  R )
{
    if ( ! has_cuda )
        return false;

    if ( handle_idx == INVALID_HANDLE )
        return false;

    using  cuda_value_t = typename cuda_traits< value_t >::cuda_type;

    auto        handle = handles[ handle_idx ];
    const auto  nrows  = M.nrows();
    const auto  ncols  = M.ncols();
    const auto  minnm  = std::min( nrows, ncols );

    const auto  type_M = cuda_traits< value_t >::blas_type;

    auto  dev_M    = device_alloc< cuda_value_t >( nrows * ncols );
    auto  dev_tau  = device_alloc< cuda_value_t >( minnm );
    auto  dev_info = device_alloc< int >( 1 );

    if (! (( dev_M    != nullptr ) &&
           ( dev_tau  != nullptr ) &&
           ( dev_info != nullptr )) )
    {
        HWARNING( "(CUDA) qr : could not allocate memory on device" );
        
        device_free( dev_info );
        device_free( dev_tau );
        device_free( dev_M );

        return false;
    }// if
    
    to_device( handle, M, dev_M );

    cusolverDnParams_t  params;

    HPRO_CUSOLVER_CHECK( cusolverDnCreateParams, ( & params ) ); // TODO: make it global?

    //
    // allocate work space
    //
    
    size_t  lwork_dev = 0;
    size_t  lwork_hst = 0;

    {
        auto  retval = cusolverDnXgeqrf_bufferSize( handle.solver, params,
                                                    nrows, ncols,
                                                    type_M, dev_M, nrows,
                                                    type_M, dev_tau,
                                                    type_M,
                                                    & lwork_dev, & lwork_hst );
        
        if ( retval != CUSOLVER_STATUS_SUCCESS )
        {
            HWARNING( to_string( "(CUDA) qr : error in \"cusolverDnXgeqrf\" (code: %d)", retval ) );
        
            cusolverDnDestroyParams( params );
            device_free( dev_info );
            device_free( dev_tau );
            device_free( dev_M );

            return false;
        }// if
    }

    {
        int   lwork  = 0;
        auto  retval = orgqr_buffersize( handle.solver,
                                         nrows, ncols, minnm,
                                         dev_M, nrows,
                                         dev_tau,
                                         & lwork );
        
        if ( retval != CUSOLVER_STATUS_SUCCESS )
        {
            HWARNING( to_string( "(CUDA) qr : error in \"cusolverDnXorgqr_bufferSize\" (code: %d)", retval ) );
        
            cusolverDnDestroyParams( params );
            device_free( dev_info );
            device_free( dev_tau );
            device_free( dev_M );

            return false;
        }// if

        lwork_dev = std::max< size_t >( lwork_dev, lwork );
    }

    auto  hst_work = std::vector< value_t >( lwork_hst );
    auto  dev_work = device_alloc< cuda_value_t >( lwork_dev );

    {
        auto  retval = cusolverDnXgeqrf( handle.solver, params,
                                         nrows, ncols,
                                         type_M, dev_M, nrows,
                                         type_M, dev_tau,
                                         type_M,
                                         dev_work, lwork_dev,
                                         hst_work.data(), lwork_hst,
                                         dev_info );
        
        if ( retval != CUSOLVER_STATUS_SUCCESS )
        {
            HWARNING( to_string( "(CUDA) qr : error in \"cusolverDnXgeqrf\" (code: %d)", retval ) );
        
            device_free( dev_work );
            cusolverDnDestroyParams( params );
            device_free( dev_info );
            device_free( dev_tau );
            device_free( dev_M );

            return false;
        }// if
    }

    auto  info = from_device< int >( handle, dev_info );

    if ( info != 0 )
    {
        if ( info < 0 )
        { HWARNING( "(CUDA) qr : " + to_string( "error in argument %d", info ) ); }

        device_free( dev_work );
        cusolverDnDestroyParams( params );
        device_free( dev_info );
        device_free( dev_tau );
        device_free( dev_M );

        return false;
    }// if

    //
    // get R
    //

    if ( ! from_device( handle, dev_M,   M   ) )
    {
        device_free( dev_work );
        cusolverDnDestroyParams( params );
        device_free( dev_info );
        device_free( dev_tau );
        device_free( dev_M );

        return false;
    }// if

    if (( blas_int_t( R.nrows() ) != ncols ) || ( blas_int_t( R.ncols() ) != ncols ))
        R = std::move( BLAS::Matrix< value_t >( ncols, ncols ) );
    else
        fill( value_t(0), R );
        
    for ( blas_int_t  i = 0; i < ncols; i++ )
    {
        BLAS::Vector< value_t >  colM( M, BLAS::Range( 0, i ), i );
        BLAS::Vector< value_t >  colR( R, BLAS::Range( 0, i ), i );
            
        copy( colM, colR );
    }// for

    //
    // compute Q
    //

    {
        auto  retval = orgqr( handle.solver,
                              nrows, ncols, minnm,
                              dev_M, nrows,
                              dev_tau,
                              dev_work, lwork_dev,
                              dev_info );
        
        if ( retval != CUSOLVER_STATUS_SUCCESS )
        {
            HWARNING( to_string( "(CUDA) qr : error in \"cusolverDnXorgqr\" (code: %d)", retval ) );
        
            device_free( dev_work );
            cusolverDnDestroyParams( params );
            device_free( dev_info );
            device_free( dev_tau );
            device_free( dev_M );

            return false;
        }// if
    }

    bool  retval = true;
    
    if ( ! from_device( handle, dev_M,   M   ) ) retval = false;

    HPRO_CUSOLVER_CHECK( cusolverDnDestroyParams, ( params ) );
    device_free( dev_work );
    device_free( dev_info );
    device_free( dev_tau );
    device_free( dev_M );

    return retval;
}

#else

template < typename value_t >
bool
qr  ( const cuda_handle_t        ,
      BLAS::Matrix< value_t > &  ,
      BLAS::Matrix< value_t > &   )
{
    return false;
}

#endif

#define INST_QR( T )                            \
    template bool                               \
    qr< T >  ( const cuda_handle_t  handle_idx, \
               BLAS::Matrix< T > &  M,          \
               BLAS::Matrix< T > &  tau )

INST_QR( float );
INST_QR( double );
INST_QR( std::complex< float > );
INST_QR( std::complex< double > );

#undef INST_QR

////////////////////////////////////////////////////////////////////////////////
//
// SVD related functions
//
////////////////////////////////////////////////////////////////////////////////

//
// compute SVD decomposition \f$ A = U·S·V^H \f$ of the nxm matrix \a A with
// n×min(n,m) matrix U, min(n,m)×min(n,m) matrix S (diagonal)
// and m×min(n,m) matrix V; \a A will be overwritten with U upon exit
//
#if HPRO_USE_CUDA == 1
    
template < typename value_t >
bool
svd  ( const cuda_handle_t                       handle_idx,
       BLAS::Matrix< value_t > &                 A,
       BLAS::Vector< real_type_t< value_t > > &  S,
       BLAS::Matrix< value_t > &                 V )
{
    if ( ! has_cuda )
        return false;

    if ( handle_idx == INVALID_HANDLE )
        return false;

    // only support for nrows >= ncols
    if ( A.nrows() < A.ncols() )
        return false;

    using  real_t       = real_type_t< value_t >;
    using  cuda_value_t = typename cuda_traits< value_t >::cuda_type;
    using  cuda_real_t  = typename cuda_traits< real_t >::cuda_type;
    
    auto        handle = handles[ handle_idx ];
    const auto  nrows  = A.nrows();
    const auto  ncols  = A.ncols();
    const auto  minnm  = std::min( nrows, ncols );

    //
    // allocate data memory on device
    //

    const auto  type_A = cuda_traits< value_t >::blas_type;
    const auto  type_S = cuda_traits< real_t >::blas_type;
    
    auto  dev_A    = device_alloc< cuda_value_t >( nrows * ncols );
    auto  dev_U    = device_alloc< cuda_value_t >( nrows * nrows );
    auto  dev_S    = device_alloc< cuda_real_t >( minnm );
    auto  dev_VT   = device_alloc< cuda_value_t >( ncols * ncols );
    auto  dev_info = device_alloc< int >( 1 );

    if (! (( dev_A    != nullptr ) &&
           ( dev_U    != nullptr ) &&
           ( dev_S    != nullptr ) &&
           ( dev_VT   != nullptr ) &&
           ( dev_info != nullptr )) )
    {
        HWARNING( "(CUDA) svd : could not allocate memory on device" );
        
        device_free( dev_info );
        device_free( dev_VT );
        device_free( dev_S );
        device_free( dev_U );
        device_free( dev_A );

        return false;
    }// if
        
    to_device( handle, A, dev_A );

    cusolverDnParams_t  params;

    HPRO_CUSOLVER_CHECK( cusolverDnCreateParams, ( & params ) ); // TODO: make it global?
    
    size_t      lwork_dev = 0;
    size_t      lwork_hst = 0;
    const char  jobU      = 'S';
    const char  jobV      = 'S';
    
    // double  pertubation_error;
    
    {
        // auto  retval = cusolverDnXgesvdp_bufferSize( handle.solver, params,
        //                                              CUSOLVER_EIG_MODE_VECTOR,
        //                                              1, // economy mode
        //                                              nrows, ncols,
        //                                              type_A, dev_A, nrows,
        //                                              type_S, dev_S,
        //                                              type_A, dev_U, nrows,
        //                                              type_A, dev_VT, minnm,
        //                                              type_A,
        //                                              & lwork_dev, & lwork_hst );
        
        auto  retval = cusolverDnXgesvd_bufferSize( handle.solver, params,
                                                    jobU, jobV,
                                                    nrows, ncols,
                                                    type_A, dev_A, nrows,
                                                    type_S, dev_S,
                                                    type_A, dev_U, nrows,
                                                    type_A, dev_VT, minnm,
                                                    type_A,
                                                    & lwork_dev, & lwork_hst );

        if ( retval != CUSOLVER_STATUS_SUCCESS )
        {
            HWARNING( to_string( "(CUDA) svd : error in \"cusolverDnXgesvd*_bufferSize\" (code: %d)", retval ) );
            
            cusolverDnDestroyParams( params );
            device_free( dev_info );
            device_free( dev_VT );
            device_free( dev_S );
            device_free( dev_U );
            device_free( dev_A );
            return false;
        }// if
    }

    auto  hst_work = std::vector< value_t >( lwork_hst );
    auto  dev_work = device_alloc< value_t >( lwork_dev );
        
    {
        // auto  retval = cusolverDnXgesvdp( handle.solver, params,
        //                                   CUSOLVER_EIG_MODE_VECTOR,
        //                                   1, // economy mode
        //                                   nrows, ncols,
        //                                   type_A, dev_A, nrows,
        //                                   type_S, dev_S,
        //                                   type_A, dev_U, nrows,
        //                                   type_A, dev_VT, minnm,
        //                                   type_A,
        //                                   dev_work, lwork_dev,
        //                                   hst_work.data(), lwork_hst,
        //                                   dev_info, & pertubation_error );

        auto  retval = cusolverDnXgesvd( handle.solver,
                                         params,
                                         jobU, jobV,
                                         nrows, ncols,
                                         type_A, dev_A, nrows,
                                         type_S, dev_S,
                                         type_A, dev_U, nrows,
                                         type_A, dev_VT, minnm,
                                         type_A,
                                         dev_work, lwork_dev,
                                         hst_work.data(), lwork_hst,
                                         dev_info );

        if ( retval != CUSOLVER_STATUS_SUCCESS )
        {
            HWARNING( to_string( "(CUDA) svd : error in \"cusolverDnXgesvd*\" (code: %d)", retval ) );
            
            cusolverDnDestroyParams( params );
            device_free( dev_work );
            device_free( dev_info );
            device_free( dev_VT );
            device_free( dev_S );
            device_free( dev_U );
            device_free( dev_A );
            return false;
        }// if
    }

    // if ( cudaStreamSynchronize( handle.stream ) != cudaSuccess )
    //     HWARNING( "(CUDA) svd : error in cuStreamSynchronize %d" );
        
    auto  info = from_device< int >( handle, dev_info );

    if ( info != 0 )
    {
        if ( info < 0 ) { HWARNING( "(CUDA) svd : " + to_string( "error in argument %d", info ) ); }
        else            { HWARNING( "(CUDA) svd : " + to_string( "no convergence during SVD: %d", info ) ); }

        device_free( dev_work );
        cusolverDnDestroyParams( params );
        device_free( dev_info );
        device_free( dev_VT );
        device_free( dev_S );
        device_free( dev_U );
        device_free( dev_A );
        return false;
    }// if
    
    //
    // get data from device
    //

    bool  retval = true;
    auto  VT     = BLAS::Matrix< value_t >( minnm, ncols );
    
    if (( A.nrows() != nrows ) || ( A.ncols() != minnm ))
        A = std::move( BLAS::Matrix< value_t >( nrows, minnm ) );
    
    if ( S.length() != minnm )
        S = std::move( BLAS::Vector< real_t >( minnm ) );

    if ( ! from_device( handle, dev_U,  A  ) ) retval = false;
    if ( ! from_device( handle, dev_S,  S  ) ) retval = false;
    if ( ! from_device( handle, dev_VT, VT ) ) retval = false;

    if (( V.nrows() != ncols ) || ( V.ncols() != minnm ))
        V = std::move( BLAS::Matrix< value_t >( ncols, minnm ) );
    
    BLAS::copy( BLAS::adjoint( VT ), V );
    
    HPRO_CUSOLVER_CHECK( cusolverDnDestroyParams, ( params ) ); // TODO: make it global?
    device_free( dev_work );
    device_free( dev_info );
    device_free( dev_VT );
    device_free( dev_S );
    device_free( dev_U );
    device_free( dev_A );

    return retval;
}

#else

template < typename value_t >
bool
svd  ( const cuda_handle_t                       ,
       BLAS::Matrix< value_t > &                 ,
       BLAS::Vector< real_type_t< value_t > > &  ,
       BLAS::Matrix< value_t > &                  )
{
    return false;
}

#endif

#define INST_SVD( T )                                           \
    template bool                                               \
    svd< T >  ( const cuda_handle_t                 handle_idx, \
                BLAS::Matrix< T > &                 U,          \
                BLAS::Vector< real_type_t< T > > &  S,          \
                BLAS::Matrix< T > &                 V )

INST_SVD( float );
INST_SVD( double );
INST_SVD( std::complex< float > );
INST_SVD( std::complex< double > );

#undef INST_SVD

}}// namespace Hpro::CUDA
