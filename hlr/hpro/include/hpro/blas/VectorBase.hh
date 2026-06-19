#ifndef __HPRO_BLAS_VECTORBASE_HH
#define __HPRO_BLAS_VECTORBASE_HH
//
// Project     : HLIBpro
// File        : Vector.hh
// Description : defines basic interface for vectors
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <ostream>

#include "hpro/blas/Range.hh"

namespace Hpro
{

namespace BLAS
{

//
// trait for giving access to vector properties
//
template < typename T_vector >
struct vector_trait;

//
// signals, that T is of vector type
//
template < typename T >
struct is_vector
{
    static const bool  value = false;
};

#if __cplusplus >= 201703L
template <typename T> inline constexpr bool is_vector_v = is_vector< T >::value;
#endif

//!
//! \ingroup  BLAS_Module
//! \class    Vector
//! \brief    defines basic interface for vectors
//!
template < typename T_derived >
class VectorBase
{
public:
    //! scalar value type of vector
    using  value_t = typename vector_trait< T_derived >::value_t;

public:

    //! return length of vector
    size_t    length      ()                  const noexcept { return derived().length(); }

    //! return stride of index set
    size_t    stride      ()                  const noexcept { return derived().length(); }

    //! return coefficient at position \a i
    value_t   operator () ( const idx_t   i ) const noexcept { return derived()( i ); }
    
    //! return reference to coefficient at position \a i
    value_t & operator () ( const idx_t   i )       noexcept { return derived()( i ); }
    
    // //! return reference to sub vector defined by \a r
    // Vector    operator () ( const Range & r ) const { return derived()( r ); }

    //! give access to internal data
    value_t * data        ()                  const noexcept { return derived().data(); }

private:
    //! convert to derived type
    T_derived &        derived  ()       noexcept { return * static_cast<       T_derived * >( this ); }
    const T_derived &  derived  () const noexcept { return * static_cast< const T_derived * >( this ); }
};

//
// stream output
//
template < typename T >
std::ostream &
operator << ( std::ostream & os, const VectorBase< T > & v );

}// namespace BLAS

}// namespace Hpro

#endif  // __HPRO_BLAS_VECTORBASE_HH
