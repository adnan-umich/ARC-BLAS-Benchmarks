//
// Project     : HLIBpro
// File        : blas/hip.cc
// Description : HIP related functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2024. All Rights Reserved.
//

#include <array>
#include <mutex>

#include <hpro/config.h>

#if HPRO_USE_HIP == 1
#  ifdef __HIP_PLATFORM_NVIDIA__
#    include <cuda_runtime.h>
#  endif

#  include <hipblas/hipblas.h>
#  include <hipsolver/hipsolver.h>
#endif

#include <hpro/blas/Algebra.hh>
#include <hpro/blas/hip.hh>

namespace Hpro { namespace HIP {

////////////////////////////////////////////////////////////////////////////////
//
// management functions
//
////////////////////////////////////////////////////////////////////////////////

#if HPRO_USE_HIP == 1

// wrapper for hip, cuBlas and hipsolver functions
#define HPRO_HIP_CHECK( func, args )                                    \
    {                                                                   \
        auto  result = func args ;                                      \
        if ( result != hipSuccess )                                     \
        {                                                               \
            HERROR( ERR_HIP, #func, Hpro::to_string( "HIP code: %d", result ) ); \
        }                                                               \
    }

#define HPRO_HIPSOLVER_CHECK( func, args )                              \
    {                                                                   \
        auto  result = func args ;                                      \
        if ( result != HIPSOLVER_STATUS_SUCCESS )                       \
        {                                                               \
            HERROR( ERR_HIP, #func, Hpro::to_string( "HIPSOLVER code: %d", result ) ); \
        }                                                               \
    }

//
// mapping of default types to cuBLAS types
//
template < typename T > struct hip_traits
{
    using  hip_type  = T;
    using  ptr_type  = T *;
    using  real_type = T;
};

template <>
struct hip_traits< float >
{
    using  hip_type  = float;
    using  ptr_type  = hip_type *;
    using  real_type = float;
};

template <>
struct hip_traits< double >
{
    using  hip_type  = double;
    using  ptr_type  = hip_type *;
    using  real_type = double;
};

template <>
struct hip_traits< std::complex< float > >
{
    using  hip_type  = hipFloatComplex;
    using  ptr_type  = hip_type *;
    using  real_type = float;
};

template <>
struct hip_traits< std::complex< double > >
{
    using  hip_type  = hipDoubleComplex;
    using  ptr_type  = hip_type *;
    using  real_type = double;
};

template <>
struct hip_traits< hipFloatComplex >
{
    using  hip_type  = hipFloatComplex;
    using  ptr_type  = hip_type *;
    using  real_type = float;
};

template <>
struct hip_traits< hipDoubleComplex >
{
    using  hip_type  = hipDoubleComplex;
    using  ptr_type  = hip_type *;
    using  real_type = double;
};

//
// joined handle for cuBLAS and hipsolverDn
//
struct handle_t
{
    hipStream_t        stream;
    hipsolverHandle_t  solver;
};

#endif

namespace
{

#if HPRO_USE_HIP == 1

//
// signals availability of HIP
//
bool  has_hip = false;

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
// initialize HIP related data structures
//
void
init ()
{
    #if HPRO_USE_HIP == 1

    //
    // first check number of devices
    //
    
    int  ndev = 0;
    
    if ( hipGetDeviceCount( & ndev ) != hipSuccess )
        return;

    HINFO( to_string( "(HIP) init : found %d HIP device(s)", ndev ) );
    
    if ( ndev == 0 )
        return;

    // get number of compute units
    int   dev_id   = 0;
    auto  dev_prop = hipDeviceProp_t();
    
    if ( hipGetDeviceProperties( & dev_prop, dev_id ) != hipSuccess )
        return;

    const auto  cu_count = dev_prop.multiProcessorCount;
    
    //
    // now try to set up streams
    // - 4 CUs per stream (try to keep 32 divisible by it)
    //

    uint32_t  cu_num      = 4;
    uint32_t  cu_pattern  = 0b1111; // four CUs
    uint32_t  cu_mask[8]  = { 0, 0, 0, 0, 0, 0, 0, 0 }; // for 8 * 32 = 256 compute units
    uint32_t  cu_idx      = 0;      // start position for active CUs in mask
    bool      found_error = false;
    uint      error_index = 0;
    
    for ( uint  i = 0; i < max_handles; ++i )
    {
        const auto  arr_ofs = cu_idx / 32;
        const auto  bit_ofs = ( cu_idx - arr_ofs * 32 );

        cu_mask[ arr_ofs ] = cu_pattern << bit_ofs;
        cu_idx            += cu_num;
        
        if ( hipExtStreamCreateWithCUMask( & handles[i].stream, 8, cu_mask ) != hipSuccess )
        {
            HWARNING( "(HIP) init : error in hipExtStreamCreateWithCUMask; disabling HIP" );
            found_error = true;
            error_index = i;
            break;
        }// if

        // after this, assumption is that HIP is present and working!!!
        
        HPRO_HIPSOLVER_CHECK( hipsolverCreate,    ( & handles[i].solver ) );
        HPRO_HIPSOLVER_CHECK( hipsolverSetStream, (   handles[i].solver, handles[i].stream ) );

        handle_used[i] = false;
    }// for

    if ( found_error )
    {
        // clean up all initialized data
        for ( uint  i = 0; i < error_index; ++i )
        {
            HPRO_HIPSOLVER_CHECK( hipsolverDestroy,  ( handles[i].solver ) );
            HPRO_HIP_CHECK(       hipStreamDestroy,  ( handles[i].stream ) );
        }// for
    }// if
    else
    {
        has_hip = true;
    }// else

    #endif
}

//
// finish HIP related data structures
//
void
done ()
{
    #if HPRO_USE_HIP == 1

    if ( ! has_hip )
        return;

    has_hip = false;
    
    for ( uint  i = 0; i < max_handles; ++i )
    {
        if ( handle_used[i] )
            HERROR( ERR_HIP, "(HIP) done", "handle still in use" );
        
        HPRO_HIPSOLVER_CHECK( hipsolverDestroy,  ( handles[i].solver ) );
        HPRO_HIP_CHECK(       hipStreamDestroy,  ( handles[i].stream ) );
    }// for

    #endif
}

//
// get free handle for HIP (stream/cuBlas/hipsolver)
//
hip_handle_t
request_handle ()
{
    #if HPRO_USE_HIP == 1
    
    if ( ! has_hip )
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
#if HPRO_USE_HIP == 1
void
release_handle ( const hip_handle_t  handle )
{
    if ( handle == INVALID_HANDLE )
        return;

    auto  lock = std::scoped_lock( handle_mtx );

    if (( handle < 0 ) || ( handle >= max_handles ))
        HERROR( ERR_HIP, "(HIP) release_handle", "unknown handle given" );

    if ( ! handle_used[ handle ] )
        HERROR( ERR_HIP, "(HIP) release_handle", "given handle not in use" );

    handle_used[ handle ] = false;

}
#else
void release_handle ( const hip_handle_t  ) {}
#endif

////////////////////////////////////////////////////////////////////////////////
//
// device handling
//
////////////////////////////////////////////////////////////////////////////////

#if HPRO_USE_HIP == 1

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
    auto       retval = hipMalloc( & ptr, n * sizeof(value_t) );

    if ( retval != hipSuccess )
    {
        HWARNING( "(HIP) device_alloc : error in \"hipMalloc\" : " + std::string( hipGetErrorString( retval ) ) );
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
    const auto  retval = hipFree( ptr );

    if ( retval != hipSuccess )
        HWARNING( "(HIP) device_alloc : error in \"hipFree\" : " + std::string( hipGetErrorString( retval ) ) );
}

//
// copy data to device
//
template < typename value_t >
bool
to_device ( handle_t                                  handle,
            const BLAS::Matrix< value_t > &           M_host,
            typename hip_traits< value_t >::ptr_type  M_dev )
{
    const auto  retval = hipMemcpyWithStream( M_dev, M_host.data(), sizeof(value_t) * M_host.nrows() * M_host.ncols(), hipMemcpyHostToDevice, handle.stream );
    
    if ( retval != hipSuccess )
    {
        HWARNING( to_string( "(HIP) to_device : error in \"hipMemcpyWithStream\" (%d)", retval ) );
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
from_device ( handle_t                                  handle,
              typename hip_traits< value_t >::ptr_type  M_dev,
              BLAS::Matrix< value_t > &                 M_host )
{
    const auto  retval = hipMemcpyWithStream( M_host.data(), M_dev, sizeof(value_t) * M_host.nrows() * M_host.ncols(), hipMemcpyDeviceToHost, handle.stream );

    if ( retval != hipSuccess )
    {
        HWARNING( to_string( "(HIP) from_device : error in \"hipMemcpyWithStream\" (%d)", retval ) );
        return false;
    }// if
    else
        return true;
}

template < typename value_t >
bool
from_device ( handle_t                                  handle,
              typename hip_traits< value_t >::ptr_type  v_dev,
              BLAS::Vector< value_t > &                 v_host )
{
    const auto  retval = hipMemcpyWithStream( v_host.data(), v_dev, sizeof(value_t) * v_host.length(), hipMemcpyDeviceToHost, handle.stream );

    if ( retval != hipSuccess )
    {
        HWARNING( to_string( "(HIP) from_device : error in \"hipMemcpyWithStream\" (%d)", retval ) );
        return false;
    }// if
    else
        return true;
}

template < typename value_t >
value_t
from_device ( handle_t                                  handle,
              typename hip_traits< value_t >::ptr_type  dev_data )
{
    value_t  data;

    const auto  retval = hipMemcpyWithStream( & data, dev_data, sizeof(typename hip_traits< value_t >::hip_type), hipMemcpyDeviceToHost, handle.stream );

    if ( retval != hipSuccess )
        HWARNING( to_string( "(HIP) from_device : error in \"hipMemcpyWithStream\" (%d)", retval ) );
    
    return data;
}

}// namespace anonymous

#endif

////////////////////////////////////////////////////////////////////////////////
//
// QR related functions
//
////////////////////////////////////////////////////////////////////////////////

#if HPRO_USE_HIP == 1

namespace
{

//
// geqrf
//

#define HPRO_HIP_GEQRF_BUFSIZE( vtype, func )       \
    hipsolverStatus_t                               \
    geqrf_bufferSize ( hipsolverHandle_t  handle,   \
                       int                nrows,    \
                       int                ncols,    \
                       vtype *            A,        \
                       int                ldA,      \
                       int *              lwork )   \
    { return func( handle, nrows, ncols, A, ldA, lwork ); }

HPRO_HIP_GEQRF_BUFSIZE( float,            hipsolverSgeqrf_bufferSize )
HPRO_HIP_GEQRF_BUFSIZE( double,           hipsolverDgeqrf_bufferSize )
HPRO_HIP_GEQRF_BUFSIZE( hipFloatComplex,  hipsolverCgeqrf_bufferSize )
HPRO_HIP_GEQRF_BUFSIZE( hipDoubleComplex, hipsolverZgeqrf_bufferSize )

#undef HPRO_HIP_GEQRF_BUFSIZE

#define HPRO_HIP_GEQRF( vtype, func )                                   \
    hipsolverStatus_t                                                   \
    geqrf ( hipsolverHandle_t  handle,                                  \
            int                nrows,                                   \
            int                ncols,                                   \
            vtype *            A,                                       \
            int                ldA,                                     \
            vtype *            tau,                                     \
            vtype *            work,                                    \
            int                lwork,                                   \
            int *              devInfo )                                \
    { return func( handle, nrows, ncols, A, ldA, tau, work, lwork, devInfo ); }

HPRO_HIP_GEQRF( float,            hipsolverSgeqrf )
HPRO_HIP_GEQRF( double,           hipsolverDgeqrf )
HPRO_HIP_GEQRF( hipFloatComplex,  hipsolverCgeqrf )
HPRO_HIP_GEQRF( hipDoubleComplex, hipsolverZgeqrf )

#undef HPRO_HIP_GEQRF

//
// orgqr
//

#define HPRO_HIP_ORGQR_BUFSIZE( vtype, func )       \
    hipsolverStatus_t                               \
    orgqr_bufferSize ( hipsolverHandle_t  handle,   \
                       int                nrows,    \
                       int                ncols,    \
                       int                nelem,    \
                       vtype *            A,        \
                       int                ldA,      \
                       vtype *            tau,      \
                       int *              lwork )   \
    { return func( handle, nrows, ncols, nelem, A, ldA, tau, lwork ); }

HPRO_HIP_ORGQR_BUFSIZE( float,            hipsolverSorgqr_bufferSize )
HPRO_HIP_ORGQR_BUFSIZE( double,           hipsolverDorgqr_bufferSize )
HPRO_HIP_ORGQR_BUFSIZE( hipFloatComplex,  hipsolverCungqr_bufferSize )
HPRO_HIP_ORGQR_BUFSIZE( hipDoubleComplex, hipsolverZungqr_bufferSize )

#undef HPRO_HIP_ORGQR_BUFSIZE

#define HPRO_HIP_ORGQR( vtype, func )                                   \
    hipsolverStatus_t                                                   \
    orgqr ( hipsolverHandle_t  handle,                                  \
            int                nrows,                                   \
            int                ncols,                                   \
            int                nelem,                                   \
            vtype *            A,                                       \
            int                ldA,                                     \
            vtype *            tau,                                     \
            vtype *            work,                                    \
            int                lwork,                                   \
            int *              devInfo )                                \
    { return func( handle, nrows, ncols, nelem, A, ldA, tau, work, lwork, devInfo ); }

HPRO_HIP_ORGQR( float,            hipsolverSorgqr )
HPRO_HIP_ORGQR( double,           hipsolverDorgqr )
HPRO_HIP_ORGQR( hipFloatComplex,  hipsolverCungqr )
HPRO_HIP_ORGQR( hipDoubleComplex, hipsolverZungqr )

#undef HPRO_HIP_ORGQR


}// namespace anonymous

template < typename value_t >
bool
qr ( const hip_handle_t         handle_idx,
     BLAS::Matrix< value_t > &  A,
     BLAS::Matrix< value_t > &  R )
{
    if ( ! has_hip )
        return false;

    if ( handle_idx == INVALID_HANDLE )
        return false;

    // only support for nrows >= ncols
    if ( A.nrows() < A.ncols() )
        return false;

    using  real_t      = real_type_t< value_t >;
    using  hip_value_t = typename hip_traits< value_t >::hip_type;
    
    auto        handle = handles[ handle_idx ];
    const auto  nrows  = A.nrows();
    const auto  ncols  = A.ncols();
    const auto  mindim = std::min( nrows, ncols );

    //
    // allocate memory on device and copy data
    //

    auto  dev_A    = device_alloc< hip_value_t >( nrows * ncols );
    auto  dev_tau  = device_alloc< hip_value_t >( mindim );
    auto  dev_info = device_alloc< int >( 1 );

    if (! (( dev_A    != nullptr ) &&
           ( dev_tau  != nullptr ) &&
           ( dev_info != nullptr )) )
    {
        HWARNING( "(HIP) qr : could not allocate memory on device" );
        
        device_free( dev_info );
        device_free( dev_tau );
        device_free( dev_A );

        return false;
    }// if
        
    if ( ! to_device( handle, A, dev_A ) )
    {
        device_free( dev_info );
        device_free( dev_tau );
        device_free( dev_A );

        return false;
    }// if

    //
    // determine size of and allocate work space
    //
    
    int  lwork_dev = 0;
    
    {
        int   work_geqrf = 0;
        int   work_orgqr = 0;

        {
            auto  retval = geqrf_bufferSize( handle.solver, nrows, ncols, dev_A, nrows, & work_geqrf );

            if ( retval != HIPSOLVER_STATUS_SUCCESS )
            {
                HWARNING( to_string( "(HIP) qr : error in \"hipsolverXgeqrf_bufferSize\" (code: %d)", retval ) );
                
                device_free( dev_info );
                device_free( dev_tau );
                device_free( dev_A );
                return false;
            }// if
        }

        {
            auto  retval = orgqr_bufferSize( handle.solver, nrows, ncols, mindim, dev_A, nrows, dev_tau, & work_orgqr );

            if ( retval != HIPSOLVER_STATUS_SUCCESS )
            {
                HWARNING( to_string( "(HIP) qr : error in \"hipsolverXorgqr_bufferSize\" (code: %d)", retval ) );
                
                device_free( dev_info );
                device_free( dev_tau );
                device_free( dev_A );
                return false;
            }// if
        }

        lwork_dev = std::max( work_geqrf, work_orgqr );
    }

    auto  dev_work = device_alloc< hip_value_t >( lwork_dev );

    if ( dev_work == nullptr )
    {
        HWARNING( "(HIP) qr : could not allocate memory on device" );
        
        device_free( dev_info );
        device_free( dev_tau );
        device_free( dev_A );

        return false;
    }// if

    //
    // perform QR
    //
    
    {
        auto  retval = geqrf( handle.solver,
                              nrows, ncols,
                              dev_A, nrows,
                              dev_tau,
                              dev_work, lwork_dev,
                              dev_info );

        if ( retval != HIPSOLVER_STATUS_SUCCESS )
        {
            HWARNING( to_string( "(HIP) qr : error in \"hipsolverXgeqrf\" (code: %d)", retval ) );
            
            device_free( dev_work );
            device_free( dev_info );
            device_free( dev_tau );
            device_free( dev_A );
            return false;
        }// if
    }

    auto  info = from_device< int >( handle, dev_info );

    if ( info != 0 )
    {
        if ( info < 0 ) { HWARNING( "(HIP) qr : " + to_string( "error in argument %d", info ) ); }
        else            { HWARNING( "(HIP) qr : " + to_string( "no convergence during SVD: %d", info ) ); }

        device_free( dev_work );
        device_free( dev_info );
        device_free( dev_tau );
        device_free( dev_A );
        return false;
    }// if
    
    //
    // copy upper triangular matrix to R
    //

    if ( ! from_device( handle, dev_A,  A ) )
    {
        HWARNING( "(HIP) qr : error in getting R from device" );
            
        device_free( dev_work );
        device_free( dev_info );
        device_free( dev_tau );
        device_free( dev_A );
        return false;
    }// if
    
    
    if (( blas_int_t( R.nrows() ) != ncols ) || ( blas_int_t( R.ncols() ) != ncols ))
        R = std::move( BLAS::Matrix< value_t >( ncols, ncols ) );
    else
        fill( value_t(0), R );
        
    for ( blas_int_t  i = 0; i < ncols; i++ )
    {
        BLAS::Vector< value_t >  colA( A, BLAS::Range( 0, i ), i );
        BLAS::Vector< value_t >  colR( R, BLAS::Range( 0, i ), i );
            
        copy( colA, colR );
    }// for

    //
    // compute Q
    //

    {
        auto  retval = orgqr( handle.solver,
                              nrows, ncols, mindim,
                              dev_A, nrows,
                              dev_tau,
                              dev_work, lwork_dev,
                              dev_info );

        if ( retval != HIPSOLVER_STATUS_SUCCESS )
        {
            HWARNING( to_string( "(HIP) qr : error in \"hipsolverXorgqr\" (code: %d)", retval ) );
            
            device_free( dev_work );
            device_free( dev_info );
            device_free( dev_tau );
            device_free( dev_A );
            return false;
        }// if
    }
    
    if ( ! from_device( handle, dev_A,  A ) )
    {
        HWARNING( "(HIP) qr : error in getting Q from device" );
            
        device_free( dev_work );
        device_free( dev_info );
        device_free( dev_tau );
        device_free( dev_A );
        return false;
    }// if

    device_free( dev_work );
    device_free( dev_info );
    device_free( dev_tau );
    device_free( dev_A );

    return true;
}

#else  // HPRO_USE_HIP

template < typename value_t >
bool
qr ( const hip_handle_t         handle_idx,
     BLAS::Matrix< value_t > &  A,
     BLAS::Matrix< value_t > &  R )
{
    return false;
}

#endif // HPRO_USE_HIP

#define HPRO_HIP_INST_QR( T )                   \
    template bool                               \
    qr< T >  ( const hip_handle_t   handle_idx, \
               BLAS::Matrix< T > &  A,          \
               BLAS::Matrix< T > &  R )

HPRO_HIP_INST_QR( float );
HPRO_HIP_INST_QR( double );
HPRO_HIP_INST_QR( std::complex< float > );
HPRO_HIP_INST_QR( std::complex< double > );

#undef HPRO_HIP_INST_QR

////////////////////////////////////////////////////////////////////////////////
//
// SVD related functions
//
////////////////////////////////////////////////////////////////////////////////

#if HPRO_USE_HIP == 1

namespace
{

//
// wrapper for gesvd_buffersize and gesvd
//
template < typename value_t >
hipsolverStatus_t
gesvd_bufferSize ( handle_t     handle,
                   signed char  jobu,
                   signed char  jobv,
                   int          nrows,
                   int          ncols,
                   int *        lwork );

#define HPRO_HIP_GESVD_BUFSIZE( vtype, func )                           \
    template <>                                                         \
    hipsolverStatus_t                                                   \
    gesvd_bufferSize< vtype > ( handle_t      handle,                   \
                                signed char   jobu,                     \
                                signed char   jobv,                     \
                                int           nrows,                    \
                                int           ncols,                    \
                                int *         lwork )                   \
    { return func( handle.solver, jobu, jobv, nrows, ncols, lwork ); }

HPRO_HIP_GESVD_BUFSIZE( float,            hipsolverSgesvd_bufferSize )
HPRO_HIP_GESVD_BUFSIZE( double,           hipsolverDgesvd_bufferSize )
HPRO_HIP_GESVD_BUFSIZE( hipFloatComplex,  hipsolverCgesvd_bufferSize )
HPRO_HIP_GESVD_BUFSIZE( hipDoubleComplex, hipsolverZgesvd_bufferSize )

#undef HPRO_HIP_GESVD_BUFSIZE

#define HPRO_HIP_GESVD( vtype, rtype, func )    \
    hipsolverStatus_t                           \
    gesvd ( handle_t     handle,                \
            signed char  jobu,                  \
            signed char  jobv,                  \
            int          nrows,                 \
            int          ncols,                 \
            vtype *      A, int lda,            \
            rtype *      S,                     \
            vtype *      U, int ldu,            \
            vtype *      V, int ldv,            \
            vtype *      work, int lwork,       \
            rtype *      rwork,                 \
            int *        info )                 \
{ return func( handle.solver, jobu, jobv, nrows, ncols, A, lda, S, U, ldu, V, ldv, work, lwork, rwork, info ); }

HPRO_HIP_GESVD(  float,  float, hipsolverSgesvd )
HPRO_HIP_GESVD( double, double, hipsolverDgesvd )
HPRO_HIP_GESVD(  hipFloatComplex,  float, hipsolverCgesvd )
HPRO_HIP_GESVD( hipDoubleComplex, double, hipsolverZgesvd )

#undef HPRO_HIP_GESVD

}// namespace anonymous

#endif

//
// compute SVD decomposition \f$ A = U·S·V^H \f$ of the nxm matrix \a A with
// n×min(n,m) matrix U, min(n,m)×min(n,m) matrix S (diagonal)
// and m×min(n,m) matrix V; \a A will be overwritten with U upon exit
//
#if HPRO_USE_HIP == 1
    
template < typename value_t >
bool
svd  ( const hip_handle_t                        handle_idx,
       BLAS::Matrix< value_t > &                 A,
       BLAS::Vector< real_type_t< value_t > > &  S,
       BLAS::Matrix< value_t > &                 V )
{
    if ( ! has_hip )
        return false;

    if ( handle_idx == INVALID_HANDLE )
        return false;

    // only support for nrows >= ncols
    if ( A.nrows() < A.ncols() )
        return false;

    using  real_t      = real_type_t< value_t >;
    using  hip_value_t = typename hip_traits< value_t >::hip_type;
    using  hip_real_t  = typename hip_traits< real_t >::real_type;
    
    auto        handle = handles[ handle_idx ];
    const auto  nrows  = A.nrows();
    const auto  ncols  = A.ncols();
    const auto  minnm  = std::min( nrows, ncols );

    //
    // allocate data memory on device
    //

    auto          dev_A     = device_alloc< hip_value_t >( nrows * ncols );
    auto          dev_U     = device_alloc< hip_value_t >( nrows * nrows );
    auto          dev_S     = device_alloc< hip_real_t >( minnm );
    auto          dev_VT    = device_alloc< hip_value_t >( ncols * ncols );
    auto          dev_info  = device_alloc< int >( 1 );
    hip_real_t *  dev_rwork = nullptr;

    if (! (( dev_A    != nullptr ) &&
           ( dev_U    != nullptr ) &&
           ( dev_S    != nullptr ) &&
           ( dev_VT   != nullptr ) &&
           ( dev_info != nullptr )) )
    {
        HWARNING( "(HIP) svd : could not allocate memory on device" );
        
        device_free( dev_info );
        device_free( dev_VT );
        device_free( dev_S );
        device_free( dev_U );
        device_free( dev_A );

        return false;
    }// if
        
    to_device( handle, A, dev_A );

    int         lwork_dev = 0;
    const char  jobU      = 'S';
    const char  jobV      = 'S';
    
    {
        auto  retval = gesvd_bufferSize< hip_value_t >( handle, jobU, jobV, nrows, ncols, & lwork_dev );

        if ( retval != HIPSOLVER_STATUS_SUCCESS )
        {
            HWARNING( to_string( "(HIP) svd : error in \"hipsolverXgesvd_bufferSize\" (code: %d)", retval ) );
            
            device_free( dev_info );
            device_free( dev_VT );
            device_free( dev_S );
            device_free( dev_U );
            device_free( dev_A );
            return false;
        }// if
    }

    auto  dev_work = device_alloc< hip_value_t >( lwork_dev );
        
    {
        auto  retval = gesvd( handle,
                              jobU, jobV,
                              nrows, ncols,
                              dev_A, nrows,
                              dev_S,
                              dev_U, nrows,
                              dev_VT, minnm,
                              dev_work, lwork_dev,
                              dev_rwork,
                              dev_info );

        if ( retval != HIPSOLVER_STATUS_SUCCESS )
        {
            HWARNING( to_string( "(HIP) svd : error in \"hipsolverXgesvd\" (code: %d)", retval ) );
            
            device_free( dev_work );
            device_free( dev_info );
            device_free( dev_VT );
            device_free( dev_S );
            device_free( dev_U );
            device_free( dev_A );
            return false;
        }// if
    }

    auto  info = from_device< int >( handle, dev_info );

    if ( info != 0 )
    {
        if ( info < 0 ) { HWARNING( "(HIP) svd : " + to_string( "error in argument %d", info ) ); }
        else            { HWARNING( "(HIP) svd : " + to_string( "no convergence during SVD: %d", info ) ); }

        device_free( dev_work );
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
svd  ( const hip_handle_t                       ,
       BLAS::Matrix< value_t > &                 ,
       BLAS::Vector< real_type_t< value_t > > &  ,
       BLAS::Matrix< value_t > &                  )
{
    return false;
}

#endif

#define HPRO_HIP_INST_SVD( T )                                  \
    template bool                                               \
    svd< T >  ( const hip_handle_t                  handle_idx, \
                BLAS::Matrix< T > &                 U,          \
                BLAS::Vector< real_type_t< T > > &  S,          \
                BLAS::Matrix< T > &                 V )

HPRO_HIP_INST_SVD( float );
HPRO_HIP_INST_SVD( double );
HPRO_HIP_INST_SVD( std::complex< float > );
HPRO_HIP_INST_SVD( std::complex< double > );

#undef HPRO_HIP_INST_SVD

}}// namespace Hpro::HIP
