//
// Project     : HLIBpro
// File        : TScotchAlgPartStrat.cc
// Description : class for algebraic clustertree construction with Scoth
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <cstdio>
#include <vector>

#include "hpro/config.h"

#if HPRO_USE_SCOTCH == 1
extern "C" {
#include "scotch/scotch.h"
}
#endif

#include "hpro/cluster/TAlgPartStrat.hh"

namespace Hpro
{

//////////////////////////////////////////////
//
// Scotch declaration and wrapper
//

#if HPRO_USE_SCOTCH == 1

namespace
{

//
// SCOTCH is not multi-thread safe, therefore guard with mutex
//
TMutex  SCOTCH_mutex;

//
// C++ call to Scotch
//
void
partition_graph ( const TGraph &               graph,
                  std::vector< SCOTCH_Num > &  part )
{
    TScopedLock  lock( SCOTCH_mutex );
    
    SCOTCH_Graph               sgraph;
    int                        retval;
    const size_t               nnodes = graph.nnodes();
    const size_t               nedges = graph.nedges();
    std::vector< SCOTCH_Num >  adj_ptrs( nnodes+1 );
    std::vector< SCOTCH_Num >  adj_nodes( nedges );

    for ( size_t  i = 0; i < nedges; ++i )
        adj_nodes[i] = SCOTCH_Num( graph.adj_nodes()[i] );
    
    for ( size_t  i = 0; i <= nnodes; ++i )
        adj_ptrs[i] = SCOTCH_Num( graph.adj_list_ptr()[i] );
    
    //
    // define graph
    //

    if (( retval = SCOTCH_graphInit( & sgraph ) ) != 0 )
        HERROR( ERR_SCOTCH, "(partition_graph) SCOTCH_graphInit", "" );
    
    retval = SCOTCH_graphBuild( & sgraph,
                                0,                   // base val (C: 0, Fortran: 1)
                                SCOTCH_Num(nnodes),  // number of vertices
                                adj_ptrs.data(),     // vertex pointers
                                nullptr,             // vertex end pointers
                                nullptr,             // vertex weights
                                nullptr,             // vertex labels
                                SCOTCH_Num(nedges),  // number of edges
                                adj_nodes.data(),    // edges
                                nullptr              // edge weights
                                );

    if ( retval != 0 )
        HERROR( ERR_SCOTCH, "(partition_graph) SCOTCH_graphBuild", "invalid graph" );

    if ( SCOTCH_graphCheck( & sgraph ) )
        HERROR( ERR_SCOTCH, "(partition_graph) SCOTCH_graphCheck", "graph data not consistent" );

    //
    // partition graph
    //
    
    SCOTCH_Strat  strat;

    retval = SCOTCH_stratInit( & strat );
    
    if ( retval != 0 )
        HERROR( ERR_SCOTCH, "(partition_graph) SCOTCH_stratInit", "" );
    
    retval = SCOTCH_graphPart( & sgraph, 2, & strat, part.data() );

    if ( retval != 0 )
        HERROR( ERR_SCOTCH, "(partition_graph) SCOTCH_graphPart", "error during graph partitioning" );

    SCOTCH_stratExit( & strat );
    SCOTCH_graphExit( & sgraph );
}

}// namespace anonymous

#endif

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TScotchAlgPartStrat (implementation)
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//
// partition graph into two sets
//
#if HPRO_USE_SCOTCH == 1
void
TScotchAlgPartStrat::partition ( const TGraph &  graph,
                                 TNodeSet &      left,
                                 TNodeSet &      right ) const
{
    //
    // use Scotch to partition graph
    //

    const size_t               nnodes   = graph.nnodes();
    size_t                     nleft    = 0;
    size_t                     nright   = 0;
    std::vector< SCOTCH_Num >  part( nnodes, 0 );

    partition_graph( graph, part );

    for ( size_t  i = 0; i < nnodes; ++i )
    {
        if ( part[i] == 0 ) nleft++;
        else                nright++;
    }// for

    //
    // build left and right sets
    //

    left.resize( nleft );
    right.resize( nright );

    for ( node_t  node = 0; node < node_t(nnodes); ++node )
    {
        if ( part[node] == 0 ) left.append( node );
        else                   right.append( node );
    }// for
}

#else

void
TScotchAlgPartStrat::partition ( const TGraph &  ,
                                 TNodeSet &      ,
                                 TNodeSet &       ) const
{
    HERROR( ERR_NOSCOTCH, "(TScotchAlgPartStratAlgo) partition", "" );
}

#endif

}// namespace Hpro
