#ifndef __HPRO_LIST_HH
#define __HPRO_LIST_HH
//
// Project     : HLIBpro
// File        : list.hh
// Description : helper functions for lists
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <algorithm>
#include <list>

namespace Hpro
{
    
//!
//! return first item and remove it from list
//!
template <typename T_container >
typename T_container::value_type
behead ( T_container &  acont )
{
    typename T_container::value_type  t = acont.front();

    acont.pop_front();
    
    return t;
}
         
//!
//! return true if \a adata is in \a alist
//!
template <typename T_container >
bool
is_in ( const T_container &                       acont,
        const typename T_container::value_type &  adata )
{
    return ( std::find( acont.begin(), acont.end(), adata ) != acont.end() );
}
         
}// namespace Hpro

#endif   // __HPRO_LIST_HH
