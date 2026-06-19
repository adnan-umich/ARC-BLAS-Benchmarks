#ifndef __HPRO_BASETYPES_HH
#define __HPRO_BASETYPES_HH
//
// Project     : HLIBpro
// File        : basetypes.hh
// Description : basic type definitions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <cstdlib>
#include <complex>

#include "hpro/config.h"

namespace Hpro
{

//
// abbr. for standard types
//
using  uchar = unsigned char;
using  uint  = unsigned int;
using  ulong = unsigned long;

//!
//! \typedef  idx_t
//! \brief    type for indices
//!
using  idx_t = long;

//!
//! \typedef  id_t
//! \brief    type for identifiers
//!
using  id_t  = ulong;

//!
//! type for type IDs
//!
using  typeid_t = uint;

}// namespace Hpro

// for compatibility
namespace HLIB = Hpro;

#endif  // __HPRO_BASETYPES_HH
