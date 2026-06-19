#ifndef __HPRO_TYPES_HH
#define __HPRO_TYPES_HH
//
// Project     : HLIBpro
// File        : types.hh
// Description : type definitions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif

#  include <afx.h>
#endif

#include <algorithm>

#include <memory>

#if __cplusplus >= 201103L
#include <utility>
#endif

#include <stdlib.h>

#include "hpro/base/basetypes.hh"
#include "hpro/config.h"

namespace Hpro
{

////////////////////////////////////////////////////////
//
// for handling of variants with respect to value type
//

enum variant_id_t
{
    REAL_FP32    = 0,
    REAL_FP64    = 1,
    COMPLEX_FP32 = 2,
    COMPLEX_FP64 = 3
};

////////////////////////////////////////////////////////
//
// various types for typical function arguments
//

//
// recursion handling
//
enum recursion_type_t : bool
{
    nonrecursive = false,
    recursive    = true
};

////////////////////////////////////////////////////////
//
// terminal output types
//

//!
//! \enum   term_charset_t
//! \brief  character set in terminal output
//!
enum term_charset_t
{
    ascii_charset,    //!< use ascii character set
    unicode_charset,  //!< use unicode character set
    auto_charset      //!< autodetect terminal capabilities
};

//!
//! \enum   progress_charset_t
//! \brief  enumeration of supported character sets in TStreamProgressBar
//!
enum term_color_t
{
    no_color,         //!< do not use color
    use_color,        //!< use color if available
    auto_color        //!< autodetect terminal capabilities
};

////////////////////////////////////////////////////////
//
// type checking and conversion in C-style
//

#define TYPE_ID( type )        TYPE_##type
#define DECLARE_TYPE( type )   const uint TYPE_ID( type ) = Hpro::RTTI::register_type( #type )
#define IS_TYPE( obj, type )   Hpro::__internal_is_type( (obj), TYPE_ID( type ) )

template <class T>
bool
__internal_is_type  ( const T *       obj,
                      const typeid_t  type )
{
    if ( obj == nullptr ) return false;
    else                  return obj->is_type( type );
}

template <class T>
bool
is_type  ( const T *       obj,
           const typeid_t  type )
{
    return __internal_is_type( obj, type );
}

// cast pointer type with minimal C++ type-checking
#define ptrcast(  ptr, type )  static_cast< type * >( ptr )
#define cptrcast( ptr, type )  static_cast< const type * >( ptr )

#define HPRO_RTTI_BASE( cname )                                         \
    virtual Hpro::typeid_t  type () const                               \
    { return TYPE_ID( cname ); }

#define HPRO_RTTI_DERIVED( cname, parent )                              \
    virtual Hpro::typeid_t  type    () const                            \
    { return TYPE_ID( cname ); }                                        \
    virtual bool      is_type ( const Hpro::typeid_t  t ) const         \
    { return (t == TYPE_ID( cname )) || parent::is_type( t ); }

}// namespace Hpro

////////////////////////////////////////////////////////
//
// disable copy operator or constructor to prevent
// warnings
//

#define DISABLE_COPY_OP( T_class ) \
    void operator = ( const T_class & ) = delete

////////////////////////////////////////////////////////
//
// adjustments for safe pointers
//

namespace Hpro
{

#if __cplusplus >= 201103L || ( defined(_MSC_VER) && _MSC_VER < 1914 )

//
// extend type checking to auto pointers
//
template <class T>
bool
__internal_is_type  ( const std::unique_ptr< T > &  obj,
                      const typeid_t                type )
{
    if ( obj.get() == NULL ) return false;
    else                     return obj->is_type( type );
}

#endif

}// namespace Hpro

#endif // __HPRO_TYPES_HH
