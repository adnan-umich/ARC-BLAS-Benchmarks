#ifndef __HPRO_UNORDERED_MAP_HH
#define __HPRO_UNORDERED_MAP_HH
//
// Project     : HLIBpro
// File        : unordered_map.hh
// Description : provide map class (unordered_map preferred)
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <tuple>
#include <unordered_map>

#include "hpro/cluster/TIndexSet.hh"

/////////////////////////////////////////////////////////////////////////////
//
// index set hash functions
//
/////////////////////////////////////////////////////////////////////////////

namespace std
{

template <>
struct hash< std::pair< int, int > >
{
    size_t operator () ( const std::pair< int, int > &  arg ) const
    {
        const auto  a = arg.first;
        const auto  b = arg.second;
        
        return  ( a >= b
                  ? a * a + a + b
                  : a + b * b );
    }
};

template <>
struct hash< Hpro::TIndexSet >
{
    size_t operator () ( const Hpro::TIndexSet &  is ) const
    {
        return ( std::hash< Hpro::idx_t >()( is.first() ) +
                 std::hash< Hpro::idx_t >()( is.last()  ) );
    }
};

template <>
struct hash< Hpro::TBlockIndexSet >
{
    size_t operator () ( const Hpro::TBlockIndexSet &  is ) const
    {
        return ( std::hash< Hpro::TIndexSet >()( is.row_is() ) +
                 std::hash< Hpro::TIndexSet >()( is.col_is() ) +
                 std::hash< Hpro::idx_t >()( is.row_is().first() ) ); // ← to make it unique w.r.t. transposition
    }
};

}// namespace std

#endif  // __HPRO_UNORDERED_MAP_HH
