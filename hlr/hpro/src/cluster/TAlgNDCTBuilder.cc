//
// Project     : HLIBpro
// File        : TNDAlgCTBuilder.cc
// Description : class for alg. ct. construction with nested dissection
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "treealg.hh"
#include "scheduler.hh"
#include "tensor.hh"

#include "hpro/cluster/TAlgCTBuilder.hh"

namespace Hpro
{

// namespace abbr.
namespace B = BLAS;

namespace
{

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// local defines
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// iterator over nodesets/graphs
// #define FOREACH( i, set )  for ( size_t  i = 0; i < (set).n_nodes(); i++ )

//
// (un)visit nodes in nodesets
//
template < typename T >
void
visit ( T &                    nodes,
        std::vector< bool > &  visited )
{
    for ( auto  node : nodes )
        visited[ node ] = true;
}

template < typename T >
void
unvisit ( T &                    nodes,
          std::vector< bool > &  visited )
{
    for ( auto  node : nodes )
        visited[ node ] = false;
}

// enable/disable debugging
#define PRINT     if ( false )

// enable/disable for computing SCCs during vertex_sep partitioning
#define COMP_SCC  1

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// local types
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// different node labels (allows or'ing)
enum { NONE       = 0x0,
       LEFT       = 0x1,
       RIGHT      = 0x2,
       VERTEX_SEP = 0x4,
       LOCAL      = 0x8 };

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// local variables
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// infinity depth
const uint    DEPTH_INF    = Limits::max<uint>();

// minimal size for parallel calls
const size_t  MIN_PAR_SIZE = 250;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// local functions
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// construct a connected graph for the vertex spepator using the connectivity information of 'graph'
std::unique_ptr< TGraph >
build_vtxsep_graph  ( const TGraph &    graph,
                      const TNodeSet &  vtxsep );

// compute avergage depth of domain sub trees
double
avg_dom_depth       ( const TCluster *  node );

}// namespace anonymous

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// TAlgNDCTBuilder
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////
//
// constructor and destructor
//

TAlgNDCTBuilder::TAlgNDCTBuilder ( const TAlgCTBuilder *  alg_ct_builder,
                                   const uint             n_min,
                                   const uint             min_leaf_lvl )
        : TAlgCTBuilder( nullptr, n_min, min_leaf_lvl )
        , _alg_ct_builder( alg_ct_builder )
        , _sync_interface_depth( CFG::Cluster::sync_interface_depth )
{
    if ( _alg_ct_builder == nullptr )
        HERROR( ERR_ARG, "(TAlgNDCTBuilder)", "ct builder is NULL" );
}

TAlgNDCTBuilder::~TAlgNDCTBuilder ()
{
}

//////////////////////////////////////////////
//
// clustertree creation
//

//
// divide for cluster tree with given partition <part>
//
std::unique_ptr< TCluster >
TAlgNDCTBuilder::divide ( const TGraph &             graph,
                          const uint                 lvl,
                          TPermutation &             perm,
                          const idx_t                idx_ofs,
                          const uint                 n_min,
                          any_const_sparse_matrix_t  S,
                          const uint                 max_lvl,
                          std::atomic< int > &       id ) const
{
    PRINT graph.print( "graph.dot" );

    //
    // stop recursion if number of DoF is small enough
    //

    if (( lvl >= this->_min_leaf_lvl ) && ( graph.nnodes() <= n_min ))
        return build_leaf( graph, idx_ofs, perm, id );

    if ( lvl > max_lvl )
    {
        HWARNING( to_string( "in (TAlgNDCTBuilder) divide : maximal tree depth reached; depth = %d", lvl ) );
        return build_leaf( graph, idx_ofs, perm, id );
    }// if

    ////////////////////////////////////////////////////////////////////
    //
    // partition graph into two disjoint nodesets
    //

    TNodeSet  left, right;

    if ( CFG::Cluster::build_scc )
        this->scc_partition( graph, left, right );
    else
        partition( graph, left, right );

    PRINT
    {
        std::vector< uint > label( graph.nnodes() );

        for ( auto  node : left  ) label[ node ] = LEFT;
        for ( auto  node : right ) label[ node ] = RIGHT;

        graph.print( "partition.dot", label );
    }

    // if one of the subsets is empty we would apply partitioning to
    // the same set one level below, therefore we stop recursion here
    if ( left.nnodes()  == 0 ) return build_leaf( graph, idx_ofs, perm, id );
    if ( right.nnodes() == 0 ) return build_leaf( graph, idx_ofs, perm, id );

    ////////////////////////////////////////////////////////////////////
    //
    // sum up connections between individual sets and order accordingly
    //

    check_flow( graph, left, right, S );

    ////////////////////////////////////////////////////////////////////
    //
    // build vertex separator
    //

    TNodeSet  vtxsep;

    build_vtx_sep( graph, left, right, vtxsep );

    PRINT
    {
        std::vector< uint > label( graph.nnodes() );

        for ( auto  node : left  )  label[ node ] = LEFT;
        for ( auto  node : right )  label[ node ] = RIGHT;
        for ( auto  node : vtxsep ) label[ node ] = VERTEX_SEP;

        graph.print( "vtxsep.dot", label );
    }

    // consistency check
    if ( left.nnodes() + right.nnodes() + vtxsep.nnodes() != graph.nnodes() )
        HERROR( ERR_CONSISTENCY, "(TAlgNDCTBuilder) divide", "lost nodes during divide" );

    //
    // recheck for empty sons
    //

    if ( left.nnodes()  == 0 ) return build_leaf( graph, idx_ofs, perm, id );
    if ( right.nnodes() == 0 ) return build_leaf( graph, idx_ofs, perm, id );

    ////////////////////////////////////////////////////////////////////
    //
    // restrict <graph> to left and right subset
    //

    std::unique_ptr< TGraph >  son_graph[2];

    son_graph[0] = graph.restrict( left  );
    son_graph[1] = graph.restrict( right );

    ////////////////////////////////////////////////////////////////////
    //
    // recursive call for both sublists (and interface)
    //

    std::unique_ptr< TCluster >  son[2], if_son;
    const idx_t                  left_ofs  = idx_ofs;
    const idx_t                  right_ofs = left_ofs + idx_t( son_graph[0]->nnodes() );

    auto  build_soncl = 
        [this,lvl,n_min,S,max_lvl,&perm,&id] ( const TGraph &  sgraph,
                                               const idx_t     sofs ) -> std::unique_ptr< TCluster >
        {
            // build first son
            auto  cl = divide( sgraph, lvl+1, perm, sofs, n_min, S, max_lvl, id );

            if ( cl.get() == nullptr )
                HERROR( ERR_NULL, "(TAlgNDCTBuilder) divide", "son cluster is NULL" );

            if ( cl->size() != sgraph.nnodes() )
                HERROR( ERR_CONSISTENCY, "(TAlgNDCTBuilder) divide", "" );

            cl->set_domain( true );

            return cl;
        };

    son[0] = build_soncl( * son_graph[0], left_ofs  );
    son[1] = build_soncl( * son_graph[1], right_ofs );
    
    // build cluster for vertex separator
    if ( vtxsep.nnodes() > 0 )
    {
        const auto    if_idx_ofs  = idx_t( right_ofs + son_graph[1]->nnodes() );
        double        reduction   = 0.0;
        // const size_t  dom_max_lvl = max( max( son[0]->depth(), son[1]->depth() ), size_t( 1 ) );
        const double  dom_avg_lvl = std::max( std::max( avg_dom_depth( son[0].get() ), avg_dom_depth( son[1].get() ) ), 1.0 );
        const uint    nd_max_lvl  = std::min( uint(lvl + 1 + dom_avg_lvl), max_lvl );

        // compute optimal reduction of clustersize to achieve same depth as domain cluster trees
        //
        //                         1
        //                       ─────
        //              ⎛ n_min ⎞ lvl
        // reduction =  ⎜───────⎟
        //              ⎝ |if|  ⎠
        //

        if ( dom_avg_lvl > 0 )
            reduction = Math::pow( double(n_min) / double(vtxsep.nnodes()),
                                   double(1)     / double(dom_avg_lvl) );

#if 0

        TNodeSet  surrounding( graph.nnodes() );

        for ( auto  i : graph )
            surrounding.append( node_t(i) );

        if_son = divide_if( graph, surrounding, vtxsep, lvl+1, if_idx_ofs, perm, nd_max_lvl, n_min,
                            TOptClusterSize( vtxsep.nnodes(), reduction ), id );

#else

        // construct a connected graph for the vertex spepator using the connectivity information of 'graph'
        auto  vtxsep_graph = build_vtxsep_graph( graph, vtxsep );

        if_son = divide_if( * vtxsep_graph.get(), lvl+1, perm, if_idx_ofs, n_min, S, nd_max_lvl,
                            TOptClusterSize( vtxsep.nnodes(), reduction ), id );

#endif

        if ( if_son.get() == nullptr )
            HERROR( ERR_NULL, "(TAlgNDCTBuilder) divide", "interface cluster is NULL" );

        if ( if_son->size() != vtxsep.nnodes() )
            HERROR( ERR_CONSISTENCY, "(TAlgNDCTBuilder) divide", "" );
    }// if

    ////////////////////////////////////////////////////////////////////
    //
    // finally create new cluster
    //

    auto  cluster = std::make_unique< TCluster >();

    cluster->set_id( id++ );
    
    if ( if_son.get() != nullptr )
    {
        cluster->set_first_last( std::min( std::min( son[0]->first(), son[1]->first() ), if_son->first() ),
                                 std::max( std::max( son[0]->last(),  son[1]->last() ),  if_son->last()  ) );
        cluster->set_nsons( 3 );
        cluster->set_son( 0, son[0].release() );
        cluster->set_son( 1, son[1].release() );
        cluster->set_son( 2, if_son.release() );
    }// if
    else
    {
        cluster->set_first_last( std::min( son[0]->first(), son[1]->first() ),
                                 std::max( son[0]->last(),  son[1]->last()  ) );
        cluster->set_nsons( 2 );
        cluster->set_son( 0, son[0].release() );
        cluster->set_son( 1, son[1].release() );
    }// else

    return cluster;
}

//
// divide a given graph of interface indices
//
std::unique_ptr< TCluster >
TAlgNDCTBuilder::divide_if ( const TGraph &           graph,
                             const TNodeSet &         surrounding,
                             const TNodeSet &         nodes,
                             const uint               lvl,
                             const idx_t              idx_ofs,
                             TPermutation &           perm,
                             const uint               max_lvl,
                             const uint               n_min,
                             const TOptClusterSize &  csize,
                             std::atomic< int > &     id ) const
{
    //
    // check if number of indices in this cluster
    // is too small to divide
    //

    if ((( lvl >= this->_min_leaf_lvl ) && ( nodes.nnodes() <= n_min )) || ( lvl >= max_lvl ))
        return build_leaf( graph, nodes, idx_ofs, perm, id );

    //
    // if size of cluster is less then optimal size, do not split
    //

    if ( _sync_interface_depth && ! csize.is_optimal( nodes.nnodes() ) )
    {
        auto  cluster = std::make_unique< TCluster >();

        cluster->set_id( id++ );
        cluster->set_nsons( 1 );

        auto  son = divide_if( graph, surrounding, nodes, lvl+1, idx_ofs, perm,
                               max_lvl, n_min, csize.recurse(), id );

        if ( son == nullptr )
            HERROR( ERR_NULL, "(TAlgNDCTBuilder) divide_if", "son cluster is NULL" );

        cluster->set_first_last( son->first(), son->last() );
        cluster->set_son( 0, son.release() );

        return cluster;
    }// if

    ////////////////////////////////////////////////////////////////////
    //
    // before partitioning, compute SCCs in graph defined by
    // <surrounding> and test for connectivity or restrict
    // local graph
    //

    const size_t  nnodes = surrounding.nnodes();
    TNodeSet      left( nnodes ), right( nnodes );
    TNodeSet      loc_sur;

#if COMP_SCC == 1
    {
        std::list< TNodeSet >  sccs;

        build_scc( graph, surrounding, sccs );

        if ( sccs.size() > 1 )
        {
            //
            // first sort out thise SCCs without interface nodes
            //

            std::vector< char >      label( graph.nnodes(), NONE );
            std::list< TNodeSet * >  nonempty_scc;

            for ( auto  node : nodes ) label[ node ] = LOCAL;

            for ( auto &  scc : sccs )
            {
                bool  has_if_node = false;

                for ( auto  node : scc )
                {
                    if ( label[ node ] == LOCAL )
                    {
                        has_if_node = true;
                        break;
                    }// if
                }// for

                if ( has_if_node )
                    nonempty_scc.push_back( & scc );
            }// for

            if ( nonempty_scc.size() > 1 )
            {
                //
                // first, put all components with just one node together
                //

                size_t  nsingle = 0;

                for ( auto  scc : nonempty_scc )
                    if ( scc->nnodes() == 1 )
                        nsingle++;

                TNodeSet  singlenodes( nsingle ); // declared outside since needed later

                if ( nsingle > 0 )
                {
                    std::list< TNodeSet * >  tmplist;

                    while ( ! nonempty_scc.empty() )
                    {
                        TNodeSet *  set = behead( nonempty_scc );

                        if ( set->nnodes() == 1 )
                            singlenodes.append( (*set)[0] );
                        else
                            tmplist.push_back( set );
                    }// while

                    nonempty_scc = tmplist;
                    nonempty_scc.push_back( & singlenodes );
                }// if

                //
                // partition SCCs directly
                //

                TMFitSched             sched;
                std::vector< int >     part;
                std::vector< double >  costs( nonempty_scc.size(), 0.0 );
                size_t                 i = 0;

                for ( auto  scc : nonempty_scc )
                    costs[i++] = double( scc->nnodes() );

                sched.schedule( 2, part, costs );

                i = 0;
                for ( auto  scc : nonempty_scc )
                {
                    if ( part[i++] == 0 )
                    {
                        for ( auto  node : *scc )
                            left.append( node );
                    }// if
                    else
                    {
                        for ( auto  node : *scc )
                            right.append( node );
                    }// if
                }// for
            }// if
            else
            {
                //
                // restrict surrounding to surrounding of remaining SCC
                //

                partition( graph, *(nonempty_scc.front()), nodes, left, right, loc_sur );
            }// else
        }// if
        else
        {
            partition( graph, surrounding, nodes, left, right, loc_sur );
        }// else
    }

#else

    partition( graph, surrounding, nodes, left, right, loc_sur );

#endif

    if ( left.nnodes() + right.nnodes() != nodes.nnodes() )
        HERROR( ERR_CONSISTENCY, "(TAlgNDCTBuilder) divide_if", "lost nodes in partitioning" );

    PRINT
    {
        std::vector< uint > label( graph.nnodes(), NONE );

        for ( auto  node : graph.nodes() ) label[ node ] = LOCAL;
        for ( auto  node : left          ) label[ node ] = LEFT;
        for ( auto  node : right         ) label[ node ] = RIGHT;

        graph.print( "partition.dot", label );
    }

    //
    // if one of the subsets is empty => stop recursion
    //

    if ( left.nnodes()  == 0 ) return build_leaf( graph, nodes, idx_ofs, perm, id );
    if ( right.nnodes() == 0 ) return build_leaf( graph, nodes, idx_ofs, perm, id );

    ////////////////////////////////////////////////////////////////////
    //
    // recursive call for both sublists (and interface)
    //

    std::unique_ptr< TCluster >  son[2];
    const auto                   left_ofs  = idx_ofs;
    const auto                   right_ofs = idx_t(left_ofs + left.nnodes());

    auto  build_soncl =
        [this,lvl,max_lvl,n_min,&graph,&loc_sur,&perm,&csize,&id] ( const TNodeSet &  snodes,
                                                                    const idx_t       sofs ) -> std::unique_ptr< TCluster >
        {
            auto  cl = divide_if( graph, loc_sur, snodes, lvl+1, sofs, perm, max_lvl, n_min, csize.recurse(), id );
            
            if ( cl.get() == nullptr )
                HERROR( ERR_NULL, "(TAlgNDCTBuilder) divide_if", "son cluster is NULL" );
            
            if ( cl->size() != snodes.nnodes() )
                HERROR( ERR_CONSISTENCY, "(TAlgNDCTBuilder) divide_if", "" );

            return cl;
        };

    son[0] = build_soncl( left,  left_ofs  );
    son[1] = build_soncl( right, right_ofs );
    
    ////////////////////////////////////////////////////////////////////
    //
    // finally create new cluster
    //

    auto  cluster = std::make_unique< TCluster >();

    cluster->set_id( id++ );
    cluster->set_first_last( std::min( son[0]->first(), son[1]->first() ),
                             std::max( son[0]->last(),  son[1]->last()  ) );
    cluster->set_nsons( 2 );
    cluster->set_son( 0, son[0].release() );
    cluster->set_son( 1, son[1].release() );

    return cluster;
}

//
// build cluster tree for vertex separator \a graph
//
std::unique_ptr< TCluster >
TAlgNDCTBuilder::divide_if ( const TGraph &             graph,
                             const uint                 lvl,
                             TPermutation &             perm,
                             const idx_t                idx_ofs,
                             const uint                 n_min,
                             any_const_sparse_matrix_t  S,
                             const uint                 max_lvl,
                             const TOptClusterSize &    csize,
                             std::atomic< int > &       id ) const
{
    PRINT graph.print( "graph.dot" );

    //
    // check if number of indices in this cluster
    // is too small to divide
    //

    if ((( lvl >= this->_min_leaf_lvl ) && ( graph.nnodes() <= n_min )) || ( lvl >= max_lvl ))
        return build_leaf( graph, idx_ofs, perm, id );

    //
    // if size of cluster is less then optimal size, do not split
    //

    if ( _sync_interface_depth && ! csize.is_optimal( graph.nnodes() ) )
    {
        auto  cluster = std::make_unique< TCluster >();

        cluster->set_nsons( 1 );

        auto  son = divide_if( graph, lvl+1, perm, idx_ofs, n_min, S, max_lvl, csize.recurse(), id );

        if ( son.get() == nullptr )
            HERROR( ERR_NULL, "(TAlgNDCTBuilder) divide_if", "son cluster is NULL" );

        cluster->set_first_last( son->first(), son->last() );
        cluster->set_son( 0, son.release() );

        return cluster;
    }// if


    ////////////////////////////////////////////////////////////////////
    //
    // partition graph into two disjoint nodesets
    //

    TNodeSet  left, right;

    if ( CFG::Cluster::build_scc )
        this->scc_partition( graph, left, right );
    else
        partition( graph, left, right );

    // consistency check
    if ( left.nnodes() + right.nnodes() != graph.nnodes() )
        HERROR( ERR_CONSISTENCY, "(TAlgNDCTBuilder) divide_if", "lost nodes during partition" );

    PRINT
    {
        std::vector< uint > label( graph.nnodes() );

        for ( auto  node : left  ) label[ node ] = LEFT;
        for ( auto  node : right ) label[ node ] = RIGHT;

        graph.print( "lrgraph.dot", label );
    }

    HDEBUG( to_string( "(TAlgNDCTBuilder) divide_if : left = %d, right = %d",
                       left.nnodes(), right.nnodes() ) );


    // if one of the subsets is empty we would apply partitioning to
    // the same set one level below, therefore we stop recursion here
    if ( left.nnodes()  == 0 ) return build_leaf( graph, idx_ofs, perm, id );
    if ( right.nnodes() == 0 ) return build_leaf( graph, idx_ofs, perm, id );


    ////////////////////////////////////////////////////////////////////
    //
    // restrict <graph> to left and right subset
    //

    std::unique_ptr< TGraph >  son_graph[2];

    son_graph[0] = graph.restrict( left  );
    son_graph[1] = graph.restrict( right );

    ////////////////////////////////////////////////////////////////////
    //
    // recursive call for both sublists
    //

    std::unique_ptr< TCluster >  son[2];
    const auto                   left_ofs  = idx_ofs;
    const auto                   right_ofs = idx_t(left_ofs + left.nnodes());

    auto  build_soncl =
        [this,lvl,n_min,max_lvl,S,&perm,&csize,&id] ( const TGraph &  sgraph,
                                                      const idx_t     sofs ) -> std::unique_ptr< TCluster >
        {
            auto  cl = divide_if( sgraph, lvl+1, perm, sofs, n_min, S, max_lvl, csize.recurse(), id );

            if ( cl.get() == nullptr )
                HERROR( ERR_NULL, "(TAlgNDCTBuilder) divide_if", "son cluster is NULL" );

            if ( cl->size() != sgraph.nnodes() )
                HERROR( ERR_CONSISTENCY, "(TAlgNDCTBuilder) divide_if", "" );

            return cl;
        };

    son[0] = build_soncl( * son_graph[0], left_ofs  );
    son[1] = build_soncl( * son_graph[1], right_ofs );
    
    ////////////////////////////////////////////////////////////////////
    //
    // finally create new cluster
    //

    auto  cluster = std::make_unique< TCluster >();

    cluster->set_id( id++ );
    cluster->set_first_last( std::min( son[0]->first(), son[1]->first() ),
                             std::max( son[0]->last(),  son[1]->last()  ) );

    cluster->set_nsons( 2 );
    cluster->set_son( 0, son[0].release() );
    cluster->set_son( 1, son[1].release() );

    return cluster;
}

//
// partitioning algorithm for vertex separators
//
void
TAlgNDCTBuilder::partition ( const TGraph &    graph,
                             const TNodeSet &  surrounding,
                             const TNodeSet &  nodes,
                             TNodeSet &        left,
                             TNodeSet &        right,
                             TNodeSet &        loc_sur ) const
{
    ////////////////////////////////////////////////////////////////////
    //
    // look for suitable startnode, e.g. node with maximal distance
    // to all other nodes in graph, for the connectivity, use the
    // surrounding component
    //

    const size_t         nnodes = surrounding.nnodes();
    std::vector< char >  label( graph.nnodes(), NONE );

    // mark all nodes in local subgraph
    for ( auto  node : surrounding ) label[ node ] = LOCAL;

    // and label all nodes in vertex separator
    for ( auto  node : nodes ) label[ node ] = VERTEX_SEP;

    PRINT
    {
        std::vector< uint > label2( graph.nnodes() );

        for ( auto  node : surrounding ) label2[ node ] = LOCAL;
        for ( auto  node : nodes )       label2[ node ] = VERTEX_SEP;

        graph.print( "graph.dot", label2 );
    }

    //
    // look for two sets in graph with maximal (or nearly maximal)
    // distance to start with
    //

    uint                 max_depth = DEPTH_INF;
    std::vector< bool >  visited( graph.nnodes(), true );  // marks _all_ nodes in graph visited

    right.append( nodes[0] );

    while ( true )
    {
        left = right;

        // unvisits only nodes in surrounding, rest remains visited
        unvisit( surrounding, visited );

        //
        // do BFS in matrix graph with given start node
        //

        const uint  depth = bfs_vtxsep( graph, left, right, visited, label, uint(nodes.nnodes()) );

        // if no larger depth could be found, stop iteration
        if (( max_depth == DEPTH_INF ) || ( depth > max_depth ))
        {
            max_depth = depth;

            // restrict left/right set to nodes in interface
            restrict_vtx( left,  label );
            restrict_vtx( right, label );

            if ( right.nnodes() == 0 )
                HERROR( ERR_CONSISTENCY, "blibli", "" );
        }// if
        else if ( depth <= max_depth )
            break;
    }// while

    // restrict right set because we did not do so in the loop
    restrict_vtx( right, label );

    // check for empty subsets
    if (( left.nnodes()  == 0 ) || ( right.nnodes() == 0 ))
        return;

    //
    // nodes which are visited in last BFS are enough for
    // connectivity, so reduce <surrounding>
    //

    size_t  nsur = 0;

    for ( auto  node : surrounding )
    {
        if ( visited[ node ] )
            nsur++;
    }// for

    loc_sur.resize( nsur );

    for ( auto  node : surrounding )
    {
        if ( visited[ node ] ) loc_sur.append( node );
        else                   label[ node ] = NONE;
    }// for

    PRINT
    {
        std::vector< uint >  label2( graph.nnodes(), NONE );

        for ( auto  node : loc_sur ) label2[ node ] = LOCAL;
        for ( auto  node : left )    label2[ node ] = LEFT;
        for ( auto  node : right )   label2[ node ] = RIGHT;

        graph.print( "graph.dot", label2 );
    }

    //
    // remove all but one in each list to avoid
    // not connected sets right from the start
    //

    const node_t  left_start  = left[0];
    const node_t  right_start = right[0];

    left.remove_all();
    left.append( left_start );
    
    right.remove_all();
    right.append( right_start );

    ////////////////////////////////////////////////////////////////////
    //
    // now that we have starting nodes for both sets, do a BFS
    // simultaneously from both sides until all nodes have been visited
    //

    TNodeSet  tleft(  nnodes );
    TNodeSet  tright( nnodes );
    TNodeSet  succ(   nnodes );
    bool      finished_left  = false;
    bool      finished_right = false;

    unvisit( loc_sur, visited );
    visited[ left_start  ] = true;
    visited[ right_start ] = true;

    tleft  = left;
    tright = right;

    while ( left.nnodes() + right.nnodes() < nodes.nnodes() )
    {
        //
        // decide, which set to update
        //

        bool  update_left  = false;
        bool  update_right = false;

        if      ( finished_left  )                   update_right = true;
        else if ( finished_right )                   update_left  = true;
        else if ( left.nnodes() < right.nnodes() ) update_left  = true;
        else                                         update_right = true;

        if ( update_left )
        {
            bfs_step( graph, tleft, succ, visited, label, LOCAL | VERTEX_SEP );

            // if no successor was found, the left set can not be enhanced
            if ( succ.nnodes() == 0 )
                finished_left = true;

            // update number of interface nodes per set
            for ( auto  node : succ  )
            {
                if ( label[ node ] == VERTEX_SEP )
                    left.append( node );
            }// for

            tleft = succ;
        }// if

        if ( update_right )
        {
            bfs_step( graph, tright, succ, visited, label, LOCAL | VERTEX_SEP );

            // if no successor was found, the right set can not be enhanced
            if ( succ.nnodes() == 0 )
                finished_right = true;

            // update number of interface nodes per set
            for ( auto  node : succ )
            {
                if ( label[ node ] == VERTEX_SEP )
                    right.append( node );
            }// for

            tright = succ;
        }// else
    }// while

}

//
// build leaf in a cluster tree defined by <nodes>
//
std::unique_ptr< TCluster >
TAlgNDCTBuilder::build_leaf ( const TGraph &        graph,
                              const TNodeSet &      nodes,
                              const idx_t           idx_ofs,
                              TPermutation &        perm,
                              std::atomic< int > &  id ) const
{
    if ( nodes.nnodes() == 0 )
        HERROR( ERR_ARG, "(TAlgNDCTBuilder) build_leaf", "empty indexset" );

    const idx_t  lb = idx_ofs;
    idx_t        ub = idx_ofs;

    for ( auto  node : nodes )
        perm[ graph.global_name()[ node ] ] = ub++;

    auto  cl = std::make_unique< TCluster >( lb, ub-1 );

    cl->set_id( id++ );

    return cl;
}

//
// build vertex separator between given sets
//
void
TAlgNDCTBuilder::build_vtx_sep ( const TGraph & graph,
                                 TNodeSet &     left,
                                 TNodeSet &     right,
                                 TNodeSet &     vtxsep ) const
{
    //
    // mark nodes in left and right set
    //

    const size_t         nnodes = graph.nnodes();
    std::vector< char >  label( nnodes, NONE );

    for ( auto  node : left  ) label[ node ] = LEFT;
    for ( auto  node : right ) label[ node ] = RIGHT;

    PRINT
    {
        std::vector< uint >  tlabel( nnodes, NONE );

        for ( size_t  i = 0; i < nnodes; i++ )
            tlabel[i] = label[i];

        graph.print( "partition.dot", tlabel );
    }

    //
    // detect interface by comparing all nodes in larger set with neighbours;
    // if a neighbour belongs to a different set, the node is on the vertex separator
    //

    TNodeSet * domain;

    if ( left.nnodes() > right.nnodes() ) domain = & left;
    else                                    domain = & right;

    vtxsep.resize( nnodes );

    for ( auto  node : (*domain) )
    {
        if ( label[ node ] == VERTEX_SEP )
            continue;

        //
        // go through neighbour list and look for nodes of the opposite set
        //

        for ( auto  neigh : graph.adj_nodes( node ) )
        {
            if (( label[neigh] != VERTEX_SEP ) && ( label[neigh] != label[node] ))
            {
                // put node into vertex separator
                label[ node ] = VERTEX_SEP;
                vtxsep.append( node );
                break;
            }// if
        }// for
    }// for

    if ( vtxsep.nnodes() == 0 )
        return;

    //
    // rebuild domain set
    //

    #if COMP_SCC == 1

    {
        //
        // put all from previous domain _not_ in vertex separator back into domain
        //
        
        char  domain_lbl = NONE;

        for ( auto  node : (*domain) )
        {
            if ( label[ node ] != VERTEX_SEP )
            {
                domain_lbl = label[ node ];
                break;
            }// if
        }// for

        (*domain).remove_all();

        if ( domain_lbl != NONE )
            for ( size_t  i = 0; i < nnodes; i++ )
                if ( label[ i ] == domain_lbl )
                    (*domain).append( node_t(i) );
    }

    #else

    //
    // do a BFS from one node of the left/right set to look for unreachable
    // nodes, which are now isolated by the interface
    //

    TNodeSet             nodes( domain->nnodes() );
    TNodeSet             succ( nnodes );
    std::vector< bool >  visited( nnodes, false );
    node_t               start      = node_t( nnodes );
    char                 domain_lbl = NONE;

    for ( auto  node : (*domain) )
    {
        if ( label[ node ] != VERTEX_SEP )
        {
            start      = node;
            domain_lbl = label[start];
            break;
        }// if
    }// for

    domain->remove_all();

    // only do this, if domain has a node left
    if ( start != nnodes )
    {
        (*domain).append( start );
        nodes.append( start );
        visited[ start ] = true;

        do
        {
            bfs_step( graph, nodes, succ, visited, label, domain_lbl );

            for ( auto  node : succ )
                (*domain).append( node );

            nodes = succ;
        } while ( succ.nnodes() > 0 );
    }// if

    //
    // all remaining, not connected nodes are now collected
    // by doing a BFS from the set of interface nodes;
    // since the surrounding graph is connected, they are
    // reachable by vertex spepator nodes
    //

    nodes = vtxsep;

    visit( nodes, visited );

    do
    {
        bfs_step( graph, nodes, succ, visited, label, domain_lbl );

        for ( auto  node : succ )
            vtxsep.append( node );

        nodes = succ;
    } while ( succ.nnodes() > 0 );

    #endif
}

//
// do a complete BFS starting at <start> in <graph> but finish
// if <max_nnodes> VERTEX_SEP nodes have been visited
// - visited == true is assumed for all nonlocal nodes
//
uint
TAlgNDCTBuilder::bfs_vtxsep ( const TGraph &              graph,
                              TNodeSet &                  start,
                              TNodeSet &                  last,
                              std::vector< bool > &       visited,
                              const std::vector< char > & label,
                              const uint                  max_nnodes ) const
{
    uint        depth     = 0;
    size_t      vis_nodes = 0;
    TNodeSet    succ(  graph.nnodes() );
    TNodeSet    nodes( graph.nnodes() );

    visit( start, visited );

    for ( auto  node : start )
        if ( label[ node ] == VERTEX_SEP )
            vis_nodes++;

    nodes = start;

    do
    {
        last = nodes;
        succ.remove_all();

        for ( auto  node : nodes )
        {
            //
            // visit successors of node
            //

            for ( auto  neigh : graph.adj_nodes( node ) )
            {
                if (( label[ neigh ] & (VERTEX_SEP | LOCAL) ) && ! visited[ neigh ] )
                {
                    visited[ neigh ] = true;

                    if ( label[ neigh ] == VERTEX_SEP )
                        vis_nodes++;

                    succ.append( neigh );
                }// if
            }// for
        }// for

        nodes = succ;

        depth++;

        if ( vis_nodes == max_nnodes )
        {
            // remember this last set
            last = succ;
            break;
        }// if
    } while ( succ.nnodes() > 0 );

    if ( vis_nodes != max_nnodes )
        HERROR( ERR_CONSISTENCY, "blublu", "" );

    return depth-1;
}

//
// collect successors to given nodes
// - search is restricted to nodes with label == <local>
//
void
TAlgNDCTBuilder::bfs_step ( const TGraph &               graph,
                            TNodeSet &                   nodes,
                            TNodeSet &                   succ,
                            std::vector< bool > &        visited,
                            const std::vector< char > &  label,
                            const char                   local  ) const
{
    succ.remove_all();

    for ( auto  node : nodes )
    {
        //
        // visit successors of node
        //

        for ( auto  neigh : graph.adj_nodes( node ) )
        {
            if (( label[ neigh ] & local ) && ! visited[ neigh ] )
            {
                visited[ neigh ] = true;
                succ.append( neigh );
            }// if
        }// for
    }// for
}

//
// restrict given node set to nodes in vertex separator
//
void
TAlgNDCTBuilder::restrict_vtx ( TNodeSet &                   nodes,
                                const std::vector< char > &  label ) const
{
    size_t  pos = 0;

    for ( auto  node : nodes )
        if ( label[ node ] == VERTEX_SEP )
            nodes[pos++] = node;

    nodes.set_nnodes( pos );
}

//
// build SCCs of graph restricted to <surrounding>
//
void
TAlgNDCTBuilder::build_scc ( const TGraph &           graph,
                             const TNodeSet &         surrounding,
                             std::list< TNodeSet > &  scc ) const
{
    //
    // do a BFS through all nodes and build a component
    // via reachability
    //

    const size_t         nnodes   = surrounding.nnodes();
    size_t               nvisited = 0;
    TNodeSet             nodes( nnodes );
    TNodeSet             succ(  nnodes );
    std::vector< bool >  visited( graph.nnodes(), false );
    std::vector< char >  label( graph.nnodes(), NONE );

    // mark all nodes in local subgraph
    for ( auto  node : surrounding ) label[ node ] = LOCAL;

    for ( auto  graph_node : graph.nodes() )
    {
        const node_t  start_node = graph_node;

        if ( ! visited[ start_node ] )
        {
            TNodeSet  component( nnodes - nvisited );

            nodes.remove_all();
            succ.remove_all();
            nodes.append( start_node );

            while ( nodes.nnodes() > 0 )
            {
                for ( auto  node : nodes )
                {
                    if ( label[ node ] != LOCAL )
                        continue;

                    component.append( node );
                    visited[ node ] = true;

                    for ( auto  neigh : graph.adj_nodes( node ) )
                    {
                        if ( ! visited[ neigh ] && ( label[ neigh ] == LOCAL ))
                        {
                            visited[neigh] = true;
                            succ.append( neigh );
                            nvisited++;
                        }// if
                    }// for
                }// for

                nodes = succ;
                succ.remove_all();
            }// while

            //
            // add SCC to return list
            //

            scc.push_back( component );
        }// if
    }// for
}

namespace
{

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// Routines to construct a connected graph for the vertex separator
// node set, while approximation the distance relation of the nodes
// of the vertex spepator.
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//
// check if 'edge_matrix' represents a connected graph
//
bool
matrix_connected ( const tensor2< char > &  edge_matrix )
{
    const size_t         nnodes                  = edge_matrix.dim0();
    size_t               number_of_visited_nodes = 0;
    std::vector< bool >  is_visited( nnodes );

    //
    // do a BFS search starting at one node and check
    // if all other nodes can be reached from it
    //

    // initialise first BFS step
    TNodeSet  bfs_new( nnodes );
    TNodeSet  bfs_old( nnodes );

    bfs_old.append( node_t(0) );
    is_visited[ 0 ] = true;
    number_of_visited_nodes++;

    //
    // do a BFS until no more new nodes are found
    //
    while( bfs_old.nnodes() > 0 )
    {
        for ( auto  node : bfs_old )
        {
            for ( idx_t  i = 0; i < idx_t( nnodes ); i ++ )
            {
                if ( ( edge_matrix( node, i ) ) && ( ! is_visited[ i ] ) )
                {
                    bfs_new.append( node_t( i ) );
                    is_visited[ i ] = true;
                    number_of_visited_nodes++;
                }// if
            }// for
        }// for

        // prepare next step of BFS
        bfs_old.remove_all();

        for ( auto  node : bfs_new )
        {
            bfs_old.append( node );
        }// for

        bfs_new.remove_all();

    }// while

    if( number_of_visited_nodes == nnodes )
        return true;
    else
        return false;
}

//
// construct a connected graph for the vertex spepator using the
// connectivity information of 'graph'
//
std::unique_ptr< TGraph >
build_vtxsep_graph ( const TGraph &    graph,
                     const TNodeSet &  vtxsep )
{
    std::unique_ptr< TGraph >  vtx_graph( graph.create() );
    const size_t               nnodes        = graph.nnodes();
    const size_t               nnodes_vtxsep = vtxsep.nnodes();

    if ( nnodes_vtxsep == 0 )
        HERROR( ERR_CONSISTENCY, "build_vtxsep_graph", "vertex separator empty" );

    // if vertex spepator consists of only one node return graph with only one node
    if ( nnodes_vtxsep == 1 )
    {
        const node_t  node = vtxsep[0];

        vtx_graph->init( 1, 0, true );

        vtx_graph->global_name()[ 0 ] = graph.global_name()[ node ];

        return vtx_graph;
    }// if

    ///////////////////////////////////////////////////////////////////////
    //
    // create matchings for the nodes of the vertex separator using the BFS
    // and collect the neighbourhood information of these matchings
    //
    ///////////////////////////////////////////////////////////////////////

    std::vector< node_t >  match_name( nnodes );
    std::vector< bool >    is_matched( nnodes );
    node_t                 id1 = 0;

    // initialise first matchings
    for ( auto  node : vtxsep )
    {
        match_name[node] = id1;
        is_matched[node] = true;

        id1++;
    }// for

    // use 'edge_matrix' to store the neighbourhood information of the matchings
    tensor2< char >  edge_matrix( nnodes_vtxsep, nnodes_vtxsep );

    for ( idx_t  i = 0; i < idx_t( nnodes_vtxsep ); ++i )
    {
        for ( idx_t  j = 0; j < idx_t( nnodes_vtxsep ); ++j )
        {
            edge_matrix( i, j ) = false;
        }// for
    }// for

    // initialise node sets of the BFS
    TNodeSet  bfs_new( nnodes );
    TNodeSet  bfs_old( nnodes );

    for ( auto  node : vtxsep )
        bfs_old.append( node );

    //
    // find matchings with the BFS as long the matchings are not connected
    //

    size_t  nedges = 0;

    while( true )
    {
        // number of new edges which are found in one BFS step
        size_t  nedges_new = 0;

        // do one BFS step, match nodes and look for connections between matchings
        for ( auto  node : bfs_old )
        {
            const node_t  node_name = match_name[ node ];

            // check the matchings of the adjacent nodes of 'node'
            for ( auto  neigh : graph.adj_nodes( node ) )
            {
                // check if neigh is already matched
                if ( is_matched[ neigh ] )
                {
                    const node_t  neigh_name = match_name[ neigh ];

                    //
                    // if a new connection is found between the different matchings 'neigh_name'
                    // and 'node_name' then add the corresponding edges to the 'edge_matrix'
                    //

                    if (( neigh_name != node_name ) && ( ! edge_matrix( neigh_name, node_name ) ))
                    {
                        edge_matrix( neigh_name, node_name ) = true;
                        edge_matrix( node_name, neigh_name ) = true;
                        nedges     += 2;
                        nedges_new += 2;
                    }// if
                }// if
                else
                {
                    //
                    // because 'neigh' is not matched, add him to the node set of
                    // the new visited nodes and match him with 'node'
                    //

                    bfs_new.append( neigh );
                    is_matched[ neigh ] = true;
                    match_name[ neigh ] = node_name;
                }// else
            }// for
        }// for

        // check if matchings represent a connected graph
        if ( nedges_new > 0 )
        {
            if ( matrix_connected( edge_matrix ) )
                break;
        }// if

        // prepare next step of BFS
        bfs_old.remove_all();

        for ( auto  node : bfs_new )
            bfs_old.append( node );

        bfs_new.remove_all();

        if ( bfs_old.nnodes() == 0 )
        {
            graph.print( "graph.dot" );
            HERROR( ERR_CONSISTENCY, "(TAlgNDCTBuilder) build_vtxsep_graph", "" );
        }// if
    }// while

    ///////////////////////////////////////////////////////////////////////
    //
    // construct a connected graph for the vertex separator with the
    // neighbourhood information of the matchings saved in 'edge_matrix'
    //
    ///////////////////////////////////////////////////////////////////////

    // build edge info for 'vtx_graph'
    vtx_graph->init( nnodes_vtxsep, nedges, true );

    size_t  pos = 0;

    for ( idx_t  i = 0; i < idx_t( nnodes_vtxsep ); i++ )
    {
        vtx_graph->adj_list_ptr()[ i ] = idx_t(pos);

        for ( idx_t  j = 0; j < idx_t( nnodes_vtxsep ); j++ )
        {
            if ( edge_matrix( i, j ) )
            {
                vtx_graph->adj_nodes()[ pos ] = node_t( j );
                pos++;
            }// if
        }// for
    }// for

    vtx_graph->adj_list_ptr()[ nnodes_vtxsep ] = idx_t( pos );

    // set global name of nodes of 'vtx_graph'
    node_t  id2 = 0;

    for ( auto  node : vtxsep )
    {
        vtx_graph->global_name()[ id2 ] = graph.global_name()[ node ];

        id2++;
    }// for

    return vtx_graph;
}

//
// compute avergage depth of domain sub trees
//
double
avg_dom_depth ( const TCluster *  node )
{
    using  cl_list_t = std::list< const TCluster * >;

    size_t     depth_sum = 0;  // sum of depth of all sub trees
    size_t     n_leafs   = 0;  // number of leafs
    cl_list_t  nodes;
    size_t     lvl = 0;

    nodes.push_back( node );

    while ( ! nodes.empty() )
    {
        cl_list_t  sons;

        while ( ! nodes.empty() )
        {
            const TCluster *  cl = behead( nodes );

            if ( cl->is_leaf() )
            {
                depth_sum += lvl;
                ++n_leafs;
            }// if
            else
            {
                for ( uint  i = 0; i < cl->nsons(); ++i )
                {
                    const TCluster *  son_i = cl->son( i );

                    if (( son_i == nullptr ) || ! son_i->is_domain() )
                        continue;

                    sons.push_back( son_i );
                }// for
            }// else
        }// while

        nodes = sons;
        ++lvl;
    }// while

    // average size is algebraic mean
    if ( n_leafs > 0 )
        return double( depth_sum ) / double( n_leafs );
    else
        return 0;
}

}// namespace anonymous

}// namespace Hpro
