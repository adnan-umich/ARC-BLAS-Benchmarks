//
// Project     : HLIBpro
// File        : TAlgCTBuilder.cc
// Description : class for algebraic clustertree construction
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <iostream>
#include <tuple>

#include "hpro/base/config.hh"

#include "list.hh"
#include "scheduler.hh"

#include "hpro/cluster/TAlgCTBuilder.hh"

namespace Hpro
{

namespace
{

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// local defines
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// for debugging
#define PRINT if ( false )

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// local types
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

enum { NONE = 0, LEFT = 1, RIGHT = 2, INTERFACE = 3 };

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// local variables
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// undefined permutation
const idx_t   UNDEF_PERM   = -1;

// minimal size for parallel calls
const size_t  MIN_PAR_SIZE = 250;

}// namespace anonymous


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// TAlgCTBuilder
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//
// constructor and destructor
//
TAlgCTBuilder::TAlgCTBuilder ( TAlgPartStrat *  part_strat,
                               const uint       n_min,
                               const uint       min_leaf_lvl )
        : _part_strat( part_strat ),
          _n_min( n_min ),
          _min_leaf_lvl(min_leaf_lvl)
{
    _high_deg_fac = 0;
    _edge_weights_mode = edge_weights_off;
}

//
// activate/deactivate (if \a fac = 0) high degree node separation
//
void
TAlgCTBuilder::set_high_deg_fac ( const uint  fac )
{
    _high_deg_fac = fac;
}

//
// activate/deactivate the use of edge weights for the graph bi-partitioning
//
void
TAlgCTBuilder::set_edge_weights_mode ( const edge_weights_mode_t  edge_weights_mode )
{
    _edge_weights_mode = edge_weights_mode;
}

//////////////////////////////////////////////
//
// clustertree creation
//

//
// build tree out of sparse matrix
//
std::unique_ptr< TClusterTree >
TAlgCTBuilder::build ( any_const_sparse_matrix_t  S,
                       const idx_t                idx_ofs ) const
{
    const auto    nrows = std::visit( [] ( auto && S_ptr ) { return S_ptr->nrows(); }, S );
    const size_t  maxid = nrows;
    TPermutation  perm( maxid );

    for ( size_t  i = 0; i < maxid; ++i )
        perm[i] = UNDEF_PERM;
    
    ////////////////////////////////////////////////////////////////////
    //
    // adjust n_min
    //

    const uint  n_min = ( _n_min == 0 ? adjust_n_min( S ) : _n_min );

    ////////////////////////////////////////////////////////////////////
    //
    // build graph out of sparse matrix
    // - filter out highly connected nodes as defined by _high_deg_fac
    //

    //
    // look for highly connected nodes and remove them
    //

    const size_t    nnodes = nrows;
    std::vector< bool >  node_mask;
    std::list< node_t >  high_deg_nodes;

    if ( _high_deg_fac > 0 )
    {
        const size_t  avg_degree = std::visit( [] ( auto && S_ptr ) { return S_ptr->avg_entries_per_row(); }, S );
        auto          rowptr     = std::visit( [] ( auto && S_ptr ) { return S_ptr->rowptr(); }, S );
        const size_t  max_degree = _high_deg_fac * avg_degree;
        size_t        nmasked    = 0;

        node_mask.resize( nnodes, false );
        
        for ( node_t  node = 0; node < node_t( nnodes ); ++node )
        {
            const size_t  degree = rowptr[node+1] - rowptr[node];
            
            if ( degree > max_degree )
            {
                HINFO( to_string( "(TAlgCTBuilder) build : high degree node : %d (degree = %d vs. %d average)",
                                  node, degree, avg_degree ) );
                node_mask[node] = true;
                high_deg_nodes.push_back( node );
                ++nmasked;
            }// if
        }// for

        HNOTICE( to_string( "(TAlgCTBuilder) build : found %d high degree nodes", nmasked ) );
    }// if

    //
    // build graph (without masked nodes)
    //
    
    std::unique_ptr< TGraph >  G( _edge_weights_mode == edge_weights_off
                                  ? std::make_unique< TGraph >( S, node_mask )
                                  : ( _edge_weights_mode == edge_weights_on
                                      ? std::make_unique< TEWGraph >( S, node_mask, false )
                                      : std::make_unique< TEWGraph >( S, node_mask, true ) )
                                  );

    ////////////////////////////////////////////////////////////////////
    //
    // compute clustertree
    //

    auto  id   = std::atomic< int >( 0 );
    auto  root = divide( *G, 0, perm, idx_ofs, n_min, S, uint(nnodes / 2), id );

    if ( root.get() == nullptr )
        HERROR( ERR_NULL, "(TAlgCTBuilder) build", "root cluster is NULL" );
    
    ////////////////////////////////////////////////////////////////////
    //
    // add highly connected nodes as another cluster and extend tree
    //

    if ( ! high_deg_nodes.empty() )
    {
        const idx_t  lb = idx_ofs + idx_t(root->size());
        idx_t        ub = lb;

        while ( ! high_deg_nodes.empty() )
            perm[ behead( high_deg_nodes ) ] = ub++;
        
        auto  high_cl  = std::make_unique< TCluster >( lb, ub-1 );
        auto  new_root = std::make_unique< TCluster >( idx_ofs, ub-1 );

        new_root->set_nsons( 2 );
        new_root->set_son( 0, root.release() );
        new_root->set_son( 1, high_cl.release() );

        // assign new root
        root = std::move( new_root );
    }// if
    
    ////////////////////////////////////////////////////////////////////
    //
    // consistency check: are all DoFs covered by clustertree
    //
    
    bool  all_handled = true;
    
    for ( size_t  i = 0; i < maxid; i++ )
    {
        if ( perm[ i ] == UNDEF_PERM )
        {
            perm[ i ]   = 0;
            all_handled = false;
            HERROR( ERR_CONSISTENCY, "(TAlgCTBuilder) build", to_string( "DoF %d is not handled", i ) );
        }// if
    }// for

    if ( ! all_handled )
    {
        root.reset( nullptr );
        return nullptr;
    }// if
    
    //
    // finally put all together in a tree
    //

    auto  perm_e2i = std::make_unique< TPermutation >( perm );
    auto  perm_i2e = std::make_unique< TPermutation >( perm );

    perm_i2e->invert();
    
    return std::make_unique< TClusterTree >( root.release(), perm_e2i.release(), perm_i2e.release() );
}

//
// divide a given graph and build corresponding cluster tree
//
std::unique_ptr< TCluster >
TAlgCTBuilder::divide ( const TGraph &             graph,
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

    if (( lvl >= _min_leaf_lvl ) && ( graph.nnodes() <= n_min ))
        return build_leaf( graph, idx_ofs, perm, id );

    if ( lvl > max_lvl )
    {
        HWARNING( to_string( "in (TAlgCTBuilder) divide : maximal tree depth reached; depth = %d", lvl ) );
        return build_leaf( graph, idx_ofs, perm, id );
    }// if
    
    ////////////////////////////////////////////////////////////////////
    //
    // partition graph into two disjoint nodesets
    //

    TNodeSet  left, right;

    if ( CFG::Cluster::build_scc )
        scc_partition( graph, left, right );
    else
        partition( graph, left, right );

    PRINT
    {
        std::vector< uint > label( graph.nnodes() );

        for ( auto  node : left  ) label[ node ] = LEFT;
        for ( auto  node : right ) label[ node ] = RIGHT;
        
        graph.print( "lrgraph.dot", label );
    }

    HDEBUG( to_string( "(TAlgCTBuilder) divide : left = %d, right = %d",
                       left.nnodes(), right.nnodes() ) );
    
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
    // restrict <graph> to left and right subset
    //

    auto  left_graph  = graph.restrict( left  );
    auto  right_graph = graph.restrict( right );
   
    ////////////////////////////////////////////////////////////////////
    //
    // recursive call for both sublists
    //

    std::unique_ptr< TCluster > son[2];
    const idx_t                 left_ofs  = idx_ofs;
    const idx_t                 right_ofs = left_ofs + idx_t(left_graph->nnodes());

    auto  build_soncl =
        [this,lvl,n_min,S,max_lvl,&perm,&id] ( const TGraph &  sgraph,
                                               const idx_t     sofs ) -> std::unique_ptr< TCluster >
        {
            auto  cl = divide( sgraph, lvl+1, perm, sofs, n_min, S, max_lvl, id );

            if ( cl.get() == nullptr )
                HERROR( ERR_NULL, "(TAlgCTBuilder) divide", "son cluster is NULL" );

            return cl;
        };

    son[0] = build_soncl( * left_graph.get(),  left_ofs  );
    son[1] = build_soncl( * right_graph.get(), right_ofs );
    
    ////////////////////////////////////////////////////////////////////
    //
    // finally create new cluster
    //
    
    auto  cluster = std::make_unique< TCluster >( std::min( son[0]->first(), son[1]->first() ),
                                                  std::max( son[0]->last(),  son[1]->last()  ) );

    cluster->set_id( id++ );
    cluster->set_nsons( 2 );
    cluster->set_son( 0, son[0].release() );
    cluster->set_son( 1, son[1].release() );
    
    return cluster;
}

//
// compute graph bi-partitioning of \a graph and store result in \a left and \a right
//
void
TAlgCTBuilder::partition ( const TGraph &  graph,
                           TNodeSet &      left,
                           TNodeSet &      right ) const
{
    if ( _part_strat != nullptr )
        _part_strat->partition( graph, left, right );
    else
        HERROR( ERR_NULL, "(TAlgCTBuilder) partition", "partitioning strategy is NULL" );
}

//
// check for SCCs before actual partitioning
//
void
TAlgCTBuilder::scc_partition ( const TGraph &  graph,
                               TNodeSet &      left,
                               TNodeSet &      right ) const
{
    if ( graph.nnodes() == 0 )
        return;

    //
    // first check, if the graph is connected; if not, we
    // perform SCC detection and simply partition those
    // SCCs as best as possible
    //
    
    const size_t      nnodes = graph.nnodes();
    std::list< TNodeSet >  sccs;

    left.resize( nnodes );
    right.resize( nnodes );
    
    graph.build_scc( sccs );

    if ( sccs.size() > 1 )
    {
        #if 1
        //
        // partition components directly
        //

        TMFitSched        sched;
        std::vector< int >     part;
        std::vector< double >  costs( sccs.size(), 0.0 );
        uint              i = 0;

        for ( const auto &  scc : sccs )
            costs[i++] = double( scc.nnodes() );

        sched.schedule( 2, part, costs );

        i = 0;
        for ( const auto &  scc : sccs )
        {
            if ( part[i++] == 0 )
            {
                for ( auto  node : scc )
                    left.append( node );
            }// if
            else
            {
                for ( auto  node : scc )
                    right.append( node );
            }// if
        }// for

        #else
        
        //
        // reorder SCCs by putting single node SCCs back
        //

        size_t                   nsingle = 0;
        std::list< TNodeSet * >  scc_ptr;
                
        for ( const auto &  scc : sccs )
        {
            if ( scc.nnodes() == 1 )
                nsingle++;
            else
                scc_ptr.push_back( & scc );
        }// for

        TNodeSet  singlenodes( nsingle );
                
        if ( nsingle > 0 )
        {
            for ( const auto &  scc : sccs )
            {
                if ( scc.nnodes() == 1 )
                    singlenodes.append( scc[0] );
            }// for

            scc_ptr.push_back( & singlenodes );
        }// if
                
        //
        // partition components directly
        //

        TMFitSched      sched;
        std::vector< int >   part;
        std::vector< real >  costs( scc.size(), 0.0 );
        uint            i = 0;

        for ( auto  scc : scc_ptr )
            costs[i++] = scc->nnodes();

        sched.schedule( 2, part, costs );

        i = 0;
        for ( auto  scc : scc_ptr )
        {
            if ( part[i++] == 0 )
            {
                for ( auto  node : * scc )
                    left.append( node );
            }// if
            else
            {
                for ( auto  node : * scc )
                    right.append( node );
            }// if
        }// for
        #endif
    }// if
    else
    {
        // if only a single SCC exists, use real partitioning algorithm
        partition( graph, left, right );
    }// else
}

//
// build leaf in a cluster tree
//
std::unique_ptr< TCluster >
TAlgCTBuilder::build_leaf ( const TGraph &        graph,
                            const idx_t           idx_ofs,
                            TPermutation &        perm,
                            std::atomic< int > &  id ) const
{
    if ( graph.nnodes() == 0 )
        HERROR( ERR_ARG, "(TAlgCTBuilder) build_leaf", "empty indexset" );

    const idx_t  lb = idx_ofs;
    idx_t        ub = idx_ofs;

    for ( size_t  i = 0; i < graph.nnodes(); i++ )
        perm[ graph.global_name()[i] ] = ub++;

    auto  cl = std::make_unique< TCluster >( lb, ub-1 );

    cl->set_id( id++ );

    return  cl;
}

//
// adjust n_min based on sparse matrix if default value of 0
// was given in constructor
// - called in "build" for each given <S>
//
uint
TAlgCTBuilder::adjust_n_min ( any_const_sparse_matrix_t  S ) const
{
    const size_t  max_per_row = std::visit( [] ( auto &&  S_ptr ) { return S_ptr->max_entries_per_row(); }, S );
    const size_t  avg_per_row = std::visit( [] ( auto &&  S_ptr ) { return S_ptr->avg_entries_per_row(); }, S );
    const size_t  n_min       = (max_per_row > 2 * avg_per_row ? avg_per_row : max_per_row);

    return uint( n_min < 20 ? 20 : n_min );
}

//
// analyze connections between subgraphs <left> and <right>
// and swap if necessary
//
namespace
{

template < typename value_t >
std::pair< size_t, real_type_t< value_t > >
count_conn ( const TNodeSet &                  nodes,
             const TSparseMatrix< value_t > *  S,
             const TGraph &                    graph,
             const std::vector< uint > &       label )
{
    size_t  edgecut      = 0;
    auto    connectivity = real_type_t< value_t >( 0 );
            
    for ( auto  node : nodes )
    {
        for ( auto  neigh : graph.adj_nodes( node ) )
        {
            if ( label[ node ] != label[ neigh ] )
            {
                connectivity += Math::abs( S->entry( graph.global_name()[node], graph.global_name()[neigh] ) );
                ++edgecut;
            }// if
        }// for
    }// for

    return { edgecut, connectivity };
}

}// namespace anonymous

void
TAlgCTBuilder::check_flow ( const TGraph &             graph,
                            TNodeSet &                 left,
                            TNodeSet &                 right,
                            any_const_sparse_matrix_t  S ) const
{
    using real_t = double;
    
    ////////////////////////////////////////////////////////////////////
    //
    // sum up connections between individual sets and order accordingly
    //

    size_t         n_edgecut1    = 0;
    size_t         n_edgecut2    = 0;
    real_t         left_to_right = 0.0;
    real_t         right_to_left = 0.0;
    std::vector< uint > label( graph.nnodes() );
    
    for ( auto  node : left  ) label[ node ] = LEFT;
    for ( auto  node : right ) label[ node ] = RIGHT;

    std::visit( [&] ( auto && S_ptr ) { std::tie( n_edgecut1, left_to_right ) = count_conn( left,  S_ptr, graph, label ); }, S );
    std::visit( [&] ( auto && S_ptr ) { std::tie( n_edgecut2, right_to_left ) = count_conn( right, S_ptr, graph, label ); }, S );

    PRINT
    {
        const size_t  n_edgecut = n_edgecut1 + n_edgecut2;
    
        std::cout << "#ec = "   << n_edgecut
                  << ", L → R " << left_to_right 
                  << ", R → L " << right_to_left
                  << std::endl;
    }// PRINT

    // apply consistent order for all clusters
    if ( left_to_right > right_to_left )
        std::swap( left, right );
}

}// namespace Hpro
