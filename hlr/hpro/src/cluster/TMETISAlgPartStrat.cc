//
// Project     : HLIBpro
// File        : TMETISAlgPartStrat.cc
// Description : class for algebraic clustertree construction
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <time.h>

#include <vector>

#include "hpro/config.h"

#if HPRO_USE_METIS == 1
extern "C" {
#include <metis.h>
}
#endif

#if defined(METIS_VER_MAJOR) && METIS_VER_MAJOR >= 5
// rename METIS type to prevent conflicts with Hpro::idx_t
using  metis_idx_t = idx_t;
#else
using  metis_idx_t = int;
#endif

#include "hpro/cluster/TAlgPartStrat.hh"

namespace Hpro
{

//////////////////////////////////////////////
//
// METIS declaration and wrapper
//

#if HPRO_USE_METIS == 1

namespace
{

//
// metis options
//

#if defined(METIS_VER_MAJOR) && METIS_VER_MAJOR >= 5
metis_idx_t  METIS_OPTIONS[ METIS_NOPTIONS ];
#else
int          METIS_OPTIONS[5] = { 0, 0, 0, 0, 0 };
#endif

//
// C++ call to METIS
//
void
partition_graph ( const TGraph &                graph,
                  std::vector< metis_idx_t > &  part )
{
    // local versions of parameters to METIS
    metis_idx_t                 nnodes     = metis_idx_t(graph.nnodes());
    metis_idx_t                 nedges     = metis_idx_t(graph.nedges());
    metis_idx_t                 npart      = 2;
    metis_idx_t                 edgecut    = 0;
    std::vector< metis_idx_t >  adj_ptrs( nnodes+1 );
    std::vector< metis_idx_t >  adj_nodes( nedges );

    
    for ( metis_idx_t  i = 0; i <= nnodes; ++i )
        adj_ptrs[i] = metis_idx_t( graph.adj_list_ptr()[i] );
    
    for ( metis_idx_t  i = 0; i < nedges; ++i )
        adj_nodes[i] = metis_idx_t( graph.adj_nodes()[i] );
    
    
    if ( graph.has_edge_weights() )
    {
        // transform edge weights of type 'weight_t' to integer values for METIS
        std::vector< metis_idx_t >  edge_weights( nedges );
        
        for ( metis_idx_t  i = 0; i < nedges; ++i )
        {
            edge_weights[ i ] = metis_idx_t( graph.edge_weight( i ) );
        }// for
    
        #if defined(METIS_VER_MAJOR) && METIS_VER_MAJOR >= 5
        
        metis_idx_t  nconstraints = 1;
        
        METIS_PartGraphRecursive( & nnodes,
                                  & nconstraints,
                                  adj_ptrs.data(),
                                  adj_nodes.data(),
                                  nullptr,
                                  nullptr,
                                  edge_weights.data(),
                                  & npart,
                                  nullptr,
                                  nullptr,
                                  METIS_OPTIONS,
                                  & edgecut,
                                  part.data() );
        
        #else

        metis_idx_t  numflag    = 0;
        metis_idx_t  wgtflag    = 1; // with edge weights
        
        METIS_PartGraphRecursive( & nnodes,          // may result in SEGFAULTS! Use "METIS_PartGraphKway" instead
                                  adj_ptrs.data(),
                                  adj_nodes.data(),
                                  nullptr,
                                  edge_weights.data(),
                                  & wgtflag,
                                  & numflag,
                                  & npart,
                                  METIS_OPTIONS,
                                  & edgecut,
                                  part.data() );

        #endif
    }// if
    else
    {
        #if defined(METIS_VER_MAJOR) && METIS_VER_MAJOR >= 5
        
        metis_idx_t  nconstraints = 1;
        
        METIS_PartGraphRecursive( & nnodes,
                                  & nconstraints,
                                  adj_ptrs.data(),
                                  adj_nodes.data(),
                                  nullptr,
                                  nullptr,
                                  nullptr,
                                  & npart,
                                  nullptr,
                                  nullptr,
                                  METIS_OPTIONS,
                                  & edgecut,
                                  part.data() );
        
        #else
        
        metis_idx_t  numflag    = 0;
        metis_idx_t  wgtflag    = 0; // no weights
        
        METIS_PartGraphRecursive( & nnodes,
                                  adj_ptrs.data(),
                                  adj_nodes.data(),
                                  nullptr,
                                  nullptr,
                                  & wgtflag,
                                  & numflag,
                                  & npart,
                                  METIS_OPTIONS,
                                  & edgecut,
                                  part.data() );
        
        #endif
    }// else
}

}// namespace anonymous

#endif

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TMETISAlgPartStrat (implementation)
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//
// ctor
//
TMETISAlgPartStrat::TMETISAlgPartStrat ( const bool use_random )
{
    #if HPRO_USE_METIS == 1 && defined(METIS_VER_MAJOR) && METIS_VER_MAJOR >= 5

    static TMutex  mutex;

    { // guard option array
        TScopedLock  lock( mutex );
        
        METIS_SetDefaultOptions( METIS_OPTIONS );

        // define C-style numbering
        METIS_OPTIONS[ METIS_OPTION_NUMBERING ] = 0;

        // seed for randomness in METIS
        METIS_OPTIONS[ METIS_OPTION_SEED ] = ( use_random ? ::time( nullptr ) : 0 );
    }
        
    #endif
}

//
// partition graph into two sets
//
void
TMETISAlgPartStrat::partition ( const TGraph &  graph,
                                TNodeSet &      left,
                                TNodeSet &      right ) const
{
    //
    // use METIS to partition graph
    //

    const size_t                nnodes   = graph.nnodes();
    size_t                      nleft    = 0;
    size_t                      nright   = 0;
    std::vector< metis_idx_t >  part( nnodes, 0 );

#if HPRO_USE_METIS == 1
    partition_graph( graph, part );
#else
    HERROR( ERR_NOMETIS, "(TMETISAlgPartStrat) partition", "" );
#endif

    for ( size_t  i = 0; i < nnodes; i++ )
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

}// namespace
