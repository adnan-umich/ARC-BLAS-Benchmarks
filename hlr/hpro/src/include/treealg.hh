#ifndef __HPRO_TREEALG_HH
#define __HPRO_TREEALG_HH
//
// Project     : HLIBpro
// File        : treealg.hh
// Description : provides tree algorithms
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/basetypes.hh"

#include "list.hh"

namespace Hpro
{

namespace Tree
{

//
// return depth of tree by doing BFS
//
template <class T>
size_t
depth ( const T * node )
{
    if ( node == nullptr )
        return 0;
    
    size_t                  depth = 0;
    std::list< const T * >  nodes;

    nodes.push_back( node );

    while ( ! nodes.empty() )
    {
        std::list< const T * >  sons;
        
        while ( ! nodes.empty() )
        {
            const T * n = behead( nodes );

            for ( uint  i = 0; i < n->nsons(); i++ )
                if ( n->son(i) != nullptr )
                    sons.push_back( n->son(i) );
        }// while

        if ( ! sons.empty() )
        {
            nodes = sons;
            sons.clear();
            
            ++depth;
        }// if
    }// while

    return depth;
}

//
// return number of nodes in tree
//
template <class T>
size_t
nnodes ( const T * node )
{
    if ( node == nullptr )
        return 0;
    
    size_t                  nnodes = 0;
    std::list< const T * >  nodes;

    nodes.push_back( node );

    while ( ! nodes.empty() )
    {
        const T * n = behead( nodes );

        ++nnodes;
        
        for ( uint  i = 0; i < n->nsons(); i++ )
            if ( n->son(i) != nullptr )
                nodes.push_back( n->son(i) );
    }// while

    return nnodes;
}

//
// return number of leaves in tree
//
template <class T>
size_t
nleaves ( const T * node )
{
    if ( node == nullptr )
        return 0;
    
    size_t                  nleaves = 0;
    std::list< const T * >  nodes;

    nodes.push_back( node );

    while ( ! nodes.empty() )
    {
        const T * n = behead( nodes );

        if ( n->nsons() == 0 )
            ++nleaves;
        else
        {
            for ( uint  i = 0; i < n->nsons(); i++ )
                if ( n->son(i) != nullptr )
                    nodes.push_back( n->son(i) );
        }// else
    }// while

    return nleaves;
}

}// namespace Tree

}// namespace Hpro

#endif  // __HPRO_TREEALG_HH
