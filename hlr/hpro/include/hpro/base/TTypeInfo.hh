#ifndef __HPRO_TTYPEINFO_HH
#define __HPRO_TTYPEINFO_HH
//
// Project     : HLIBpro
// File        : TTypeInfo.hh
// Description : baseclass for RTTI aware classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <string>

#include "hpro/base/types.hh"
#include "hpro/base/System.hh"

namespace Hpro
{

//!
//! \class  TTypeInfo
//! \brief  Provides basic interface and methods for RTTI.
//!
//!         In \HLIBpro an independent RTTI system is used to simplify
//!         type checks and provide readable output, e.g. in error messages.
//!
//!         All classes that should be handled by this RTTI have to be derived
//!         from TTypeInfo and must implement the functions "type" and "is_type"
//!         (for hierarchical type checks, e.g. test against base classes).
//!
class TTypeInfo
{
public:
    //! return type ID of object
    virtual typeid_t     type     ()                   const = 0;

    //! return true if local object is of given type ID \a t 
    virtual bool         is_type  ( const typeid_t t ) const { return t == type(); }

    //! return string representation of type
    virtual std::string  typestr  ()                   const { return RTTI::id_to_type( type() ); }
};

//!
//! functional version of TTypeInfo::typestr (also for non-TTypeInfo objects)
//!
template < typename T >
std::string
typestr ( const T &  obj )
{
    return obj.typestr();
}

template < typename T >
std::string
typestr ( const T *  obj )
{
    return obj->typestr();
}

}// namespace Hpro

#endif  // __HPRO_TTYPEINFO_HH
