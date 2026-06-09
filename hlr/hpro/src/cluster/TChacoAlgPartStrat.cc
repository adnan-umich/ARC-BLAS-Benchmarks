//
// Project     : HLIBpro
// File        : TChacoAlgPartStrat.cc
// Description : class for algebraic clustertree construction with Chaco
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include <stdio.h>

#include "hpro/config.h"
#include "hpro/cluster/TAlgPartStrat.hh"

namespace Hpro
{

//////////////////////////////////////////////
//
// Chaco declaration and wrapper
//

#if HPRO_USE_CHACO == 1

// external call to CHACO routine
extern "C" {
    int interface ( int       nvtxs,                /* number of vertices in full graph */
                    int     * start,                /* start of edge list for each vertex */
                    int     * adjacency,            /* edge list data */
                    int     * vwgts,                /* weights for all vertices */
                    float   * ewgts,                /* weights for all edges */
                    float   * x,                    /* coordinates for inertial method */
                    float   * y,
                    float   * z,
                    char    * outassignname,        /* name of assignment output file */
                    char    * outfilename,          /* output file name */
                    short   * assignment,           /* set number of each vtx (length n) */
                    int       architecture,         /* 0 => hypercube, d => d-dimensional mesh */
                    int       ndims_tot,            /* total number of cube dimensions to divide */
                    int       mesh_dims[3],         /* dimensions of mesh of processors */
                    double  * goal,                 /* desired set sizes for each set */
                    int       global_method,        /* global partitioning algorithm */
                    int       local_method,         /* local partitioning algorithm */
                    int       rqi_flag,             /* should I use RQI/Symmlq eigensolver? */
                    int       vmax,                 /* how many vertices to coarsen down to? */
                    int       ndims,                /* number of eigenvectors (2^d sets) */
                    double    eigtol,               /* tolerance on eigenvectors */
                    long      seed                  /* for random graph mutations */
                    );
}

namespace
{

//
// C++ call to Scotch
//
void
partition_graph ( const uint              n,
                  const vector< idx_t > & xadj,
                  const vector< idx_t > & adjncy,
                  vector< idx_t > &       part )
{
    const uint       nnz = adjncy.size();
    vector< int >    adjacency( nnz, 0 );
    vector< short >  tpart( n, 0 );
    static int       mesh_dims[3] = { 0, 0, 0 };
    
    for ( uint i = 0; i <  nnz; i++ ) adjacency[i] = adjncy[i]+1;
    
    interface( n, const_cast< int * >( xadj.carray() ), const_cast< int * >( adjacency.carray() ),
               nullptr, nullptr,
               nullptr, nullptr, nullptr,
               nullptr, nullptr,
               tpart.carray(),
               0, 1, mesh_dims, nullptr,
               1, 1,  // global,local method
               0,     // 0: Lanczos, 1: RQI/Symmlq
               20,    // vmax for small graphs
               1,     // bisection
               1e-2,  // eigen solver tolerance
               1      // seed for random numbers
               );  

    for ( uint i = 0; i < n; i++ )
        part[i] = tpart[i];
}

}// namespace anonymous

#endif

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TChacoAlgPartStrat (implementation)
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//
// partition graph into two sets
//
void
TChacoAlgPartStrat::partition ( const TGraph &  graph,
                                TNodeSet &      left,
                                TNodeSet &      right ) const
{
    //
    // use Chaco to partition graph
    //

    const size_t          nnodes   = graph.nnodes();
    size_t                nleft    = 0;
    size_t                nright   = 0;
    std::vector< idx_t >  part( nnodes, 0 );

#if HPRO_USE_CHACO == 1
    partition_graph( nnodes, graph.adj_list_ptr(), graph.adj_nodes(), part );
#else
    HERROR( ERR_NOCHACO, "(TChacoAlgPartStratAlgo) partition", "" );
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
