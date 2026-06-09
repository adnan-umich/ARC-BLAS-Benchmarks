//
// Project     : HLIBpro
// File        : TMongooseAlgPartStrat.cc
// Description : class for algebraic clustertree construction
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <time.h>

#include <vector>

#include "hpro/config.h"

#if HPRO_USE_MONGOOSE == 1
#include <Mongoose.hpp>
#endif

#include "hpro/cluster/TAlgPartStrat.hh"

namespace Hpro
{

//////////////////////////////////////////////
//
// Mongoose declaration and wrapper
//

#if HPRO_USE_MONGOOSE == 1

namespace
{

//
// C++ call to Mongoose
//
void
partition_graph ( const TGraph &        graph,
                  std::vector< int > &  part )
{
    // local versions of parameters to Mongoose
    Mongoose::Int                 nnodes = int(graph.nnodes());
    Mongoose::Int                 nedges = int(graph.nedges());
    std::vector< Mongoose::Int >  adj_ptrs( nnodes+1 );
    std::vector< Mongoose::Int >  adj_nodes( nedges );

    for ( int  i = 0; i <= nnodes; ++i )
        adj_ptrs[i] = Mongoose::Int( graph.adj_list_ptr()[i] );
    
    for ( int  i = 0; i < nedges; ++i )
        adj_nodes[i] = Mongoose::Int( graph.adj_nodes()[i] );
    
    if ( graph.has_edge_weights() )
    {
        HERROR( ERR_NOT_IMPL, "", "" );
    }// if
    else
    {
        auto  options = Mongoose::EdgeCut_Options::create();
        auto  mgraph  = Mongoose::Graph::create( nnodes, nedges, adj_ptrs.data(), adj_nodes.data() );
        auto  edgecut = Mongoose::edge_cut( mgraph, options );

        part.resize( nnodes );
        
        for ( size_t  i = 0; i < size_t(nnodes); ++i )
        {
            if ( edgecut->partition[i] ) part[i] = 1;
            else                         part[i] = 0;
        }// for

        edgecut->~EdgeCut();
        mgraph->~Graph();
        options->~EdgeCut_Options();
    }// else
}

}// namespace anonymous

#endif

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TMongooseAlgPartStrat (implementation)
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//
// ctor
//
TMongooseAlgPartStrat::TMongooseAlgPartStrat ()
{
}

//
// partition graph into two sets
//
void
TMongooseAlgPartStrat::partition ( const TGraph &  graph,
                                TNodeSet &      left,
                                TNodeSet &      right ) const
{
    //
    // use Mongoose to partition graph
    //

    const size_t        nnodes   = graph.nnodes();
    size_t              nleft    = 0;
    size_t              nright   = 0;
    std::vector< int >  part( nnodes, 0 );

#if HPRO_USE_MONGOOSE == 1
    partition_graph( graph, part );
#else
    HERROR( ERR_NOMONGOOSE, "(TMongooseAlgPartStrat) partition", "" );
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

}// namespace Hpro
