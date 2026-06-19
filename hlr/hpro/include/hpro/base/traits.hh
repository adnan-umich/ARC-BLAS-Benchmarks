#ifndef __HPRO_BLAS_TRAITS_HH
#define __HPRO_BLAS_TRAITS_HH
//
// Project     : HLIBpro
// File        : traits.hh
// Description : provide type traits for BLAS functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/types.hh"

namespace Hpro
{

//!
//! signals complex valued types
//!
template <typename T> struct is_complex_type                            { static constexpr bool value = false; };
template <>           struct is_complex_type< std::complex< float > >   { static constexpr bool value = true; };
template <>           struct is_complex_type< std::complex< double > >  { static constexpr bool value = true; };

#if __cplusplus >= 201703L
template <typename T> inline constexpr bool is_complex_type_v = is_complex_type< T >::value;
#endif

//!
//! signals single/double precision
//!
template <typename T> struct is_single_prec                            { static constexpr bool value = false; };
template <>           struct is_single_prec< float >                   { static constexpr bool value = true; };
template <>           struct is_single_prec< std::complex< float > >   { static constexpr bool value = true; };

template <typename T> struct is_double_prec                            { static constexpr bool value = false; };
template <>           struct is_double_prec< double >                  { static constexpr bool value = true; };
template <>           struct is_double_prec< std::complex< double > >  { static constexpr bool value = true; };

#if __cplusplus >= 201703L
template <typename T> inline constexpr bool is_single_prec_v = is_single_prec< T >::value;
template <typename T> inline constexpr bool is_double_prec_v = is_double_prec< T >::value;
#endif

//!
//! provide real valued type forming base of T
//!
template <typename T> struct real_type                       { using  type_t = T; };
template <typename T> struct real_type< std::complex< T > >  { using  type_t = T; };

template <typename T> using real_type_t = typename real_type< T >::type_t;

//!
//! result type for type combinations
//!
template < typename T1, typename T2 > struct promote_type { using type_t = T1; };

template <> struct promote_type< float,  float  > { using type_t = float;  };
template <> struct promote_type< float,  double > { using type_t = double; };
template <> struct promote_type< double, float  > { using type_t = double; };
template <> struct promote_type< double, double > { using type_t = double; };
    
template <> struct promote_type< std::complex< float >,  std::complex< float >  > { using type_t = std::complex< float >;  };
template <> struct promote_type< std::complex< float >,  std::complex< double > > { using type_t = std::complex< double >; };
template <> struct promote_type< std::complex< double >, std::complex< float >  > { using type_t = std::complex< double >; };
template <> struct promote_type< std::complex< double >, std::complex< double > > { using type_t = std::complex< double >; };
    
template < typename T1, typename T2 > using promote_type_t = typename promote_type< T1, T2 >::type_t;
    
//!
//! access value_t of classes
//!
template <typename T> struct value_type { using  type_t = typename T::value_t; };
template <typename T> using value_type_t = typename value_type< T >::type_t;


//
// test if T1 and T2 are of the same type
//
template < typename T1, typename T2 > struct is_same_type           { static constexpr bool  value = false; };
template < typename T1 >              struct is_same_type< T1, T1 > { static constexpr bool  value = true;  };

#if __cplusplus >= 201703L
template < typename T1, typename T2 > inline constexpr bool is_same_type_v = is_same_type< T1, T2 >::value;
#endif

}// namespace Hpro

#endif  // __HPRO_BLAS_TRAITS_HH
