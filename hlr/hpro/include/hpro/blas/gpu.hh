#ifndef __HPRO_BLAS_GPU_HH
#define __HPRO_BLAS_GPU_HH
//
// Project     : HLIBpro
// File        : blas/gpu.hh
// Description : wrapper for GPU related functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2024. All Rights Reserved.
//

#include <hpro/blas/hip.hh>
#include <hpro/blas/cuda.hh>

namespace Hpro { namespace GPU {

#if HPRO_USE_CUDA == 1

using gpu_handle_t = CUDA::cuda_handle_t;

using CUDA::INVALID_HANDLE;
using CUDA::init;
using CUDA::done;
using CUDA::request_handle;
using CUDA::release_handle;
using CUDA::qr;
using CUDA::svd;

#elif HPRO_USE_HIP == 1

using gpu_handle_t = HIP::hip_handle_t;

using HIP::INVALID_HANDLE;
using HIP::init;
using HIP::done;
using HIP::request_handle;
using HIP::release_handle;
using HIP::qr;
using HIP::svd;

#else

constexpr int  INVALID_HANDLE = -1;

using gpu_handle_t = int;

inline void init () {}
inline void done () {}

inline gpu_handle_t request_handle () { return INVALID_HANDLE; }
inline void         release_handle ( const gpu_handle_t ) {}

template < typename value_t >
bool
qr  ( const gpu_handle_t         ,
      BLAS::Matrix< value_t > &  ,
      BLAS::Matrix< value_t > &   )
{
    return false;
}

template < typename value_t >
bool
svd  ( const gpu_handle_t                        ,
       BLAS::Matrix< value_t > &                 ,
       BLAS::Vector< real_type_t< value_t > > &  ,
       BLAS::Matrix< value_t > &                  )
{
    return false;
}

#endif

}}// namespace Hpro::GPU

#endif // __HPRO_BLAS_GPU_HH
