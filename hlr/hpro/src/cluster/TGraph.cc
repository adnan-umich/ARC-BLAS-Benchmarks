//
// Project     : HLIBpro
// File        : TGraph.cc
// Description : classes for representing graphs
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <set>
#include <map>
#include <fstream>

#include "hpro/cluster/TAlgCTBuilder.hh"

#include "hpro/cluster/TGraph.hh"

namespace Hpro
{

using std::list;
using std::set;
using std::map;
using std::string;

namespace
{

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// local variables
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

// for log_base(·) 
// (the smaller base, the more important are large weights)
// - computed value due to bug in clang compiler
const double  LOG_BASE = 0.405465108108164; // std::log( 1.5 );
    
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// local functions
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//
// linear mapping of coefficients to weights
//
template < typename T_weight,
           typename T_coeff >
T_weight
transform_lin ( const T_coeff  coeff,
                const T_coeff  base )
{
    if ( std::is_integral< T_weight >::value )
    {
        if ( std::is_floating_point< T_coeff >::value )
        {
            if ( base == T_coeff(0) )
                HERROR( ERR_CONSISTENCY, "transform_lin", "base parameter is zero" );
            
            return T_weight( 80000 * Math::ceil( Math::abs( coeff / base ) ) ); //WARNING: danger of overflow?
        }// if
        else
        {
            return T_weight( coeff );
        }// else
    }// if
    else
    {
        return T_weight( 80000 * coeff );
    }// else
}

//
// logarithmic mapping of coefficients to weights
//
template < typename T_weight,
           typename T_coeff >
T_weight
transform_log ( const T_coeff  coeff,
                const T_coeff  base )
{
    if ( std::is_integral< T_weight >::value )
    {
        if ( std::is_floating_point< T_coeff >::value )
        {
            if ( base == T_coeff(0) )
                HERROR( ERR_CONSISTENCY, "transform_log", "base parameter is zero" );

            T_coeff  c = Math::abs( coeff / base );
    
            if ( c < 1.0 ) c = 1.0; 
    
            return T_weight( Math::ceil( ( Math::log( c ) / LOG_BASE ) + 1 ) );
        }// if
        else
        {
            return T_weight( coeff );
        }// else
    }// if
    else
    {
        return T_weight( coeff );
    }// else
}

}// namespace anonymous

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TGraph
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//
// constructor
//
TGraph::TGraph ( any_const_sparse_matrix_t    S,
                 const std::vector< bool > &  node_mask )
{
    std::visit(
        [this,&node_mask] ( auto &&  S_ptr )
        {
            if ( S_ptr == nullptr )
                return;
    
            ////////////////////////////////////////////////////////////////////
            //
            // construct undirected, e.g. symmetric, graph out of sparse matrix
            // first copy structure from S to edge-data and add missing
            // symmetry counterparts, then copy data into real graph
            //

            //
            // build mapping from masked nodes to original names
            //

            const size_t           n_nodes   = S_ptr->nrows();
            const bool             has_mask  = ( node_mask.size() == n_nodes );
            std::vector< node_t >  name_map( n_nodes );
            node_t                 node_name = 0;
            size_t                 nexcluded = 0;

            for ( node_t  node = 0; node < node_t( n_nodes ); node++ )
            {
                if ( has_mask && node_mask[node] )
                    ++nexcluded;
                else
                    name_map[node] = node_name++;
            }// for
    
            //
            // compute edge information, with masked nodes removed
            //

            std::vector< size_t >  degrees( n_nodes );
            size_t                 n_edges = 0;

            for ( node_t  node = 0; node < node_t( n_nodes ); node++ )
            {
                if ( has_mask && node_mask[node] )
                    continue;
                
                const auto  lb = S_ptr->rowptr(node);
                const auto  ub = S_ptr->rowptr(node+1);
            
                for ( idx_t  j = lb; j < ub; j++ )
                {
                    const auto  neigh = S_ptr->colind(j);

                    // no masked nodes or single loops
                    if (( has_mask && node_mask[neigh] ) || ( node == neigh ))
                        continue;
                
                    ++degrees[node];
                    ++n_edges;

                    // also add edge from neighbour if missing
                    if ( ! S_ptr->has_entry( neigh, node ) )
                    {
                        ++degrees[neigh];
                        ++n_edges;
                    }// if
                }// for
            }// for
    
            //
            // build graph without masked nodes
            //
    
            idx_t  pos = 0;
    
            init( n_nodes - nexcluded, n_edges, true );

            // build adj_list_ptr-data first
            for ( size_t  i = 0; i < n_nodes; i++ )
            {
                if ( has_mask && node_mask[i] )
                    continue;
        
                const node_t  node = name_map[i];

                adj_list_ptr()[node] = pos;
                pos                 += idx_t( degrees[i] );
            }// for    

            adj_list_ptr()[ n_nodes - nexcluded ] = pos;

            // now copy the edges (use <degrees> as counter)
            for ( node_t  node = 0; node < node_t( n_nodes ); node++ )
                degrees[node] = 0;
    
            for ( node_t  node = 0; node < node_t( n_nodes ); node++ )
            {
                if ( has_mask && node_mask[node] )
                    continue;
        
                const node_t  lnode = name_map[node];
                const idx_t   lb    = S_ptr->rowptr(node);
                const idx_t   ub    = S_ptr->rowptr(node+1);
            
                for ( idx_t  j = lb; j < ub; j++ )
                {
                    const node_t  neigh = S_ptr->colind(j);

                    if (( has_mask && node_mask[ neigh ] ) || ( neigh == node ))
                        continue;

                    const node_t  lneigh = name_map[ neigh ];
                    const idx_t   idx    = adj_list_ptr()[lnode] + idx_t(degrees[lnode]);
            
                    adj_nodes()[ idx ] = lneigh;
                    degrees[lnode]++;
            
                    if ( ! S_ptr->has_entry( neigh, node ) )
                    {
                        const idx_t  nidx = adj_list_ptr()[lneigh] + idx_t(degrees[lneigh]);

                        adj_nodes()[ nidx ] = lnode;
                        degrees[lneigh]++;
                    }// if
                }// while
            }// for

            // add global names
            for ( size_t  i = 0; i < n_nodes; i++ )
            {
                if ( has_mask && node_mask[i] )
                    continue;
        
                global_name()[name_map[i]] = node_t(i);
            }// for
        },
        S );
}

//
// return node with minimal/maximal degree wrt. marked subgraph
//
node_t
TGraph::min_degree_node () const
{
    if ( nnodes() == 0 )
        HERROR( ERR_ARG, "(TGraph) min_degree_node", "empty graph" );
    
    node_t  min_node   = 0;
    size_t  min_degree = degree( min_node );
    
    for ( node_t  node = 1; node < node_t( nnodes() ); node++ )
    {
        const size_t  d = degree( node );

        if ( d < min_degree )
        {
            min_node   = node;
            min_degree = d;
        }// if
    }// for

    return min_node;
}

node_t
TGraph::max_degree_node () const
{
    if ( nnodes() == 0 )
        HERROR( ERR_ARG, "(TGraph) max_degree_node", "empty graph" );
    
    node_t  max_node   = 0;
    size_t  max_degree = degree( max_node );
    
    for ( node_t  node = 1; node < node_t( nnodes() ); node++ )
    {
        const size_t  d = degree( node );

        if ( d < max_degree )
        {
            max_node   = node;
            max_degree = d;
        }// if
    }// for

    return max_node;
}

//
// compute strongly connected components in undirected graphs
//
void
TGraph::build_scc ( list< TNodeSet > &  scc ) const
{
    const size_t         n_nodes  = nnodes();
    size_t               nvisited = 0;
    TNodeSet             curr_nodes( n_nodes );
    TNodeSet             succ_nodes(  n_nodes );
    std::vector< bool >  visited( n_nodes, false );
    
    for ( auto  start_node : nodes() )
    {
        if ( ! visited[ start_node ] )
        {
            TNodeSet  component( n_nodes - nvisited );

            curr_nodes.remove_all();
            succ_nodes.remove_all();
            curr_nodes.append( start_node );

            while ( curr_nodes.nnodes() > 0 )
            {
                for ( auto  node : curr_nodes  )
                {
                    component.append( node );
                    visited[node] = true;

                    for ( auto  neigh : adj_nodes( node ))
                    {
                        if ( ! visited[ neigh ] )
                        {
                            visited[neigh] = true;
                            succ_nodes.append( neigh );
                            nvisited++;
                        }// if
                    }// for
                }// for

                curr_nodes = succ_nodes;
                succ_nodes.remove_all();
            }// while

            component.resize( component.nnodes() );
            scc.push_back( component );
        }// if
    }// for
}

//
// compute strongly connected components in undirected graphs
// w.r.t. marked nodes
//
void
TGraph::build_scc ( list< TNodeSet > &           scc,
                    const std::vector< uint > &  label,
                    const uint                   mark ) const
{
    //
    // do a BFS through all nodes and build a component
    // via reachability
    //

    const size_t         n_nodes  = nnodes();
    size_t               nvisited = 0;
    TNodeSet             curr_nodes( n_nodes );
    TNodeSet             succ_nodes(  n_nodes );
    std::vector< bool >  visited( n_nodes, false );
    
    for ( auto  start_node : nodes() )
    {
        if ( ! visited[ start_node ] && ( label[start_node] == mark ) )
        {
            TNodeSet  component( n_nodes - nvisited );

            curr_nodes.remove_all();
            succ_nodes.remove_all();
            curr_nodes.append( start_node );

            while ( curr_nodes.nnodes() > 0 )
            {
                for ( auto  node : curr_nodes )
                {
                    component.append( node );
                    visited[node] = true;

                    for ( auto  neigh : adj_nodes( node ) )
                    {
                        if ( ! visited[ neigh ] && ( label[ neigh ] == mark ))
                        {
                            visited[neigh] = true;
                            succ_nodes.append( neigh );
                            nvisited++;
                        }// if
                    }// for
                }// for

                curr_nodes = succ_nodes;
                succ_nodes.remove_all();
            }// while

            component.resize( component.nnodes() );
            scc.push_back( component );
        }// if
    }// for
}

//
// restrict <graph> to nodes in <nodes> and return resulting subgraph
//
std::unique_ptr< TGraph >
TGraph::restrict ( const TNodeSet &  anodes ) const
{
    std::unique_ptr< TGraph >  subgraph( this->create() );

    {
        if ( has_edge_weights() != subgraph->has_edge_weights() )
            HERROR( ERR_CONSISTENCY, "(TGraph) restrict", "graph and subgraph are not of the same type" );
    
        //
        // mark local nodes
        //

        std::vector< bool >  local( nnodes(), false );

        for ( auto  node : anodes )
            local[ node ] = true;
    
        //
        // count edge number and map node names
        //

        size_t                 n_edges = 0;
        std::vector< node_t >  namemap( nnodes(), 0 );
        node_t                 id     = 0;
    
        for ( auto  node : anodes )
        {
            namemap[node] = id++;
        
            for ( auto  neigh : adj_nodes( node ) )
            {
                if ( local[neigh] )
                    n_edges++;
            }// for
        }// for

        //
        // init subgraph, copy global names and fill edge data
        //
    
        const bool  global = ( global_name().size() > 0 );

        if ( global )
            subgraph->init( anodes.nnodes(), n_edges, true );
        else
            subgraph->init( anodes.nnodes(), n_edges, false );

        id = 0;
    
        for ( auto  node : anodes )
        {
            const node_t  newnode = namemap[node];
        
            subgraph->adj_list_ptr()[ newnode ] = id;

            if ( global )
                subgraph->global_name()[ newnode ]  = global_name()[ node ];
        
            for ( auto  neigh_weight : adj_nodes_weights( node ) )
            {
                const node_t  neigh = neigh_weight.first;
                
                if ( local[ neigh ] )
                {
                    subgraph->adj_nodes()[ id ] = namemap[neigh];
                
                    if ( has_edge_weights() )
                        subgraph->set_edge_weight( id , neigh_weight.second );
                
                    id++;
                }// if
            }// for
        }// for

        subgraph->adj_list_ptr()[ anodes.nnodes() ] = id;
    }

    //
    // DEBUG
    //
    
    if ( false )
    {
        //
        // test new graph
        //

        if ( subgraph->nnodes() != anodes.nnodes() ) HERROR( ERR_CONSISTENCY, "", "" );

        if ( subgraph->_adj_list_ptr.size() != subgraph->nnodes() + 1 )
            HERROR( ERR_CONSISTENCY, "", "" );

        for ( auto  node : subgraph->nodes() )
        {
            for ( auto  neigh : subgraph->adj_nodes( node ) )
            {
                if ( neigh >= node_t(subgraph->nnodes()) )
                    HERROR( ERR_CONSISTENCY, "", "" );

                bool  found = false;
                
                for ( auto  neigh2 : subgraph->adj_nodes( neigh ) )
                {
                    if ( neigh2 == node )
                    {
                        found = true;
                        break;
                    }// if
                }// for

                if ( ! found )
                    HERROR( ERR_CONSISTENCY, "", "" );
            }// for
        }// for
    }// if
    
    return subgraph;
}

//
// compute vertex separator between given subgraphs
// (defined by <label>) and put nodes into interface
//
void
TGraph::vertex_separator ( std::vector< uint > &  label,
                           const TNodeSet &       left,
                           const TNodeSet &       right,
                           TNodeSet &             vertex_sep,
                           const uint             if_label ) const
{
#if 0
    //
    // use vertex cover approximation algorihm applied to
    // edge separator to compute vertex separator, e.g. find
    // maximal matching E' and add all endpoints of edges in
    // E' to vertex cover
    //
    // edge seperator is determined by comparing labels of
    // endpoints, if they differ, the edge is part of sep.
    //

    const uint           n_nodes = n_nodes();
    std::vector< uint >  edge_del( nedges(), false );

    vertex_sep.resize( n_nodes );
    
    for ( uint node = 0; node < n_nodes; node++ )
    {
        const uint lb = adj_list_ptr()[node];
        const uint ub = adj_list_ptr()[node+1];

        for ( uint j = lb; j < ub; j++ )
        {
            // check if edge was previously deleted from graph
            if ( edge_del[j] )
                continue;
            
            const uint neigh = adj_nodes()[j];

            if (( label[node]  != if_label ) &&
                ( label[neigh] != if_label ) &&
                ( label[node]  != label[neigh] ) )
            {
                //
                // found another edge in separator, now mark
                // both nodes and remove all edges incident
                // to <node> or <neigh>
                //

                label[node]  = if_label;
                label[neigh] = if_label;

                for ( uint i = lb; i < ub; i++ )
                    edge_del[i] = true;

                const uint nlb = adj_list_ptr()[neigh];
                const uint nub = adj_list_ptr()[neigh+1];

                for ( uint i = nlb; i < nub; i++ )
                    edge_del[i] = true;

                break;
            }// if
        }// for
    }// for

    //
    // finally collect marked nodes to given set
    //

    for ( uint i = 0; i < n_nodes; i++ )
        if ( label[i] == if_label )
            vertex_sep.append( i );
    
#else
    
    //
    // detect interface by comparing all nodes in graph with adjacent nodes;
    // if a neighbour belongs to a different set, the node is on the interface
    //

    const TNodeSet * domain;
        
    if ( left.nnodes() > right.nnodes() )
        domain = & left;
    else
        domain = & right;
    
    vertex_sep.resize( nnodes() );
    
    for ( auto  node : (*domain) )
    {
        if ( label[ node ] == if_label )
            continue;
        
        //
        // go through neighbour list and look for nodes of the opposite set
        //
        
        for ( auto  neigh : adj_nodes( node ) )
        {
            if (( label[neigh] != if_label ) && ( label[neigh] != label[node] ))
            {
                // put node into interface
                label[ node ] = if_label;
                vertex_sep.append( node );
                break;
            }// if
        }// for
    }// for

#endif
}
    
//
// print graph
//
void
TGraph::print ( const std::string &  filename,
                const bool           global ) const
{
    std::ofstream  file( filename.c_str() );

    //
    // print graph in graphviz format
    //

    file << "graph G { " << std::endl
         << "size = \"8,8\";" << std::endl;

    //
    // print edges in subgraphs and adjacent to subgraph
    //
    
    for ( auto  node : nodes() )
    {
        for ( auto  neigh : adj_nodes( node ) )
        {
            if ( neigh > node )
            {
                if ( global )
                    file << to_string( "%d -- %d;", global_name()[node], global_name()[neigh] ) << std::endl;
                else
                    file << to_string( "%d -- %d;", node, neigh ) << std::endl;
            }// if
        }// for
    }// for
    
    // finish
    file << "} " << std::endl;
}

//
// print graph with labels
//
void
TGraph::print ( const std::string &          filename,
                const std::vector< uint > &  label,
                const bool                   global ) const
{
    std::ofstream  file( filename.c_str() );

    //
    // print graph in graphviz format
    //

    file << "graph G { " << std::endl
         << "size = \"8,8\";" << std::endl;

    //
    // determine different labels and assign colours 
    //

    set< uint >        lbl_set;
    map< uint, uint >  col_map;
    uint               ncolor = 1;

    for ( size_t  i = 0; i < nnodes(); i++ )
    {
        // label = 0 means unlabeled
        if ( label[i] == 0 )
            continue;
        
        if ( lbl_set.find( label[i] ) == lbl_set.end() )
        {
            lbl_set.insert( label[i] );
            col_map[ label[i] ] = ncolor;
            ++ncolor;
        }// if
    }// for

    if ( ncolor > 7 )
        HERROR( ERR_ARG, "(TGraph) print", "too many different non-zero labels (max = 6)" );
    
    //
    // label nodes
    //

    std::vector< string >  colours( 7 );

    colours[1] = "red";
    colours[2] = "steelblue";
    colours[3] = "green";
    colours[4] = "yellow";
    colours[5] = "pink";
    colours[6] = "brown";

    for ( auto  lbl : lbl_set )
    {
        file << "node [style=filled,fillcolor=" << colours[ col_map[ lbl ] ] << "];" << std::endl;

        for ( size_t  i = 0; i < nnodes(); i++ )
        {
            if ( label[i] == lbl )
            {
                file << ( global ? global_name()[i] : i ) << "; ";
            }// if
        }// for
        
        file << std::endl;
    }// for
    
    //
    // remaining nodes (label == 0 ) are not filled
    //
    
    file << "node [style=solid];" << std::endl;

    //
    // print edges in subgraphs and adjacent to subgraph
    //
    
    for ( auto  node : nodes() )
    {
        for ( auto  neigh : adj_nodes( node ) )
        {
            if ( neigh > node )
            {
                if ( global )
                    file << to_string( "%d -- %d;", global_name()[node], global_name()[neigh] ) << std::endl;
                else
                    file << to_string( "%d -- %d;", node, neigh ) << std::endl;
            }// if
        }// for
    }// for
    
    // finish
    file << "} " << std::endl;
}

//
// write in Chaco/Jostle/Metis file format
//
void
TGraph::write ( const std::string &  filename ) const
{
    std::ofstream  file( filename.c_str() );

    // first line: number of nodes and edges (counted only in one direction)
    file << nnodes() << ' ' << nedges() / 2 << " 0" << std::endl;
    
    // now, <n_nodes> lines with adjacent nodes
    for ( auto  node : nodes() )
    {
        for ( auto  neigh : adj_nodes( node ) )
        {
            // node names numbered from 1..n_nodes
            file << ' ' << neigh+1 << ' ';
        }// for

        file << std::endl;
    }// for
}
    
//
// copy operator
//
TGraph &
TGraph::operator = ( const TGraph & graph )
{
    _adj_list_ptr = graph._adj_list_ptr;
    _adj_nodes    = graph._adj_nodes;
    _global_name  = graph._global_name;

    return *this;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TEWGraph
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//
// constructor
//
TEWGraph::TEWGraph ( any_const_sparse_matrix_t    S,
                     const std::vector< bool > &  node_mask,
                     const bool                   sym_edge_weights )
{
    std::visit(
        [this,&node_mask,sym_edge_weights] ( auto &&  S_ptr )
        {
            using  value_t = typename std::remove_pointer_t< std::decay_t< decltype(S_ptr) > >::value_t;
            using  real_t  = real_type_t< value_t >;
    
            if ( S_ptr == nullptr || ( S_ptr->n_non_zero() == 0 ) )
                return;

            ////////////////////////////////////////////////////////////////////
            //
            // construct undirected, e.g. symmetric, graph out of sparse matrix
            // first copy structure from S to edge-data and add missing
            // symmetry counterparts, then copy data into real graph
            //

            //
            // build mapping from masked nodes to original names
            //

            const size_t           n_nodes    = S_ptr->rows();
            const bool             has_mask   = ( node_mask.size() == n_nodes );
            std::vector< node_t >  name_map( n_nodes );
            node_t                 node_name  = 0;
            size_t                 nexcluded  = 0;

            for ( node_t  node = 0; node < node_t( n_nodes ); node++ )
            {
                if ( has_mask && node_mask[node] )
                    ++nexcluded;
                else
                    name_map[node] = node_name++;
            }// for
    
            //
            // compute edge information, with masked nodes removed;
            // also: determine minimal coefficient (wrt. absolute value)
            // for edge weights
            //

            std::vector< size_t >  degrees( n_nodes );
            size_t                 n_edges = 0;
            real_t                 min_val = real_t(0);
    
            for ( node_t  node = 0; node < node_t( n_nodes ); node++ )
            {
                if ( has_mask && node_mask[node] )
                    continue;
        
                const idx_t  lb = S_ptr->rowptr(node);
                const idx_t  ub = S_ptr->rowptr(node+1);
            
                for ( idx_t j = lb; j < ub; j++ )
                {
                    const node_t  neigh = S_ptr->colind(j);

                    // no masked nodes or single loops
                    if (( has_mask && node_mask[neigh] ) || ( node == neigh ))
                        continue;
                
                    ++degrees[node];
                    ++n_edges;

                    // update minimal coefficient
                    const auto  val = Math::abs( S_ptr->coeff( j ) );

                    if ( val > real_t(0) )
                    {
                        if ( min_val == 0 ) min_val = val;
                        else                min_val = std::min( min_val, val );
                    }// if
            
                    // also add edge from neighbour if missing
                    if ( ! S_ptr->has_entry( neigh, node ) )
                    {
                        ++degrees[neigh];
                        ++n_edges;
                    }// if
                }// for
            }// for
    
            //
            // build graph without masked nodes
            //
    
            idx_t  pos = 0;
    
            init( n_nodes - nexcluded, n_edges, true );

            // build adj_list_ptr-data first
            for ( size_t  i = 0; i < n_nodes; i++ )
            {
                if ( has_mask && node_mask[i] )
                    continue;
        
                const node_t  node = name_map[i];

                adj_list_ptr()[node] = pos;
                pos                 += idx_t( degrees[i] );
            }// for    

            adj_list_ptr()[ n_nodes - nexcluded ] = pos;

            // now copy the edges (use <degrees> as counter)
            for ( node_t  node = 0; node < node_t( n_nodes ); node++ )
                degrees[node] = 0;
    
            for ( node_t  node = 0; node < node_t( n_nodes ); node++ )
            {
                if ( has_mask && node_mask[node] )
                    continue;
        
                const node_t  lnode      = name_map[node];
                const idx_t   lb         = S_ptr->rowptr(node);
                const idx_t   ub         = S_ptr->rowptr(node+1);
                const auto    coeff_node = Math::abs( S_ptr->entry( node, node ) );
            
                for ( idx_t  j = lb; j < ub; j++ )
                {
                    const node_t  neigh = S_ptr->colind(j);

                    if (( has_mask && node_mask[ neigh ] ) || ( neigh == node ))
                        continue;

                    const auto      coeff_neigh = Math::abs( S_ptr->entry( neigh, neigh ) );
                    const node_t    lneigh = name_map[ neigh ];
                    const idx_t     idx    = adj_list_ptr()[lnode] + idx_t(degrees[lnode]);
                    const auto      coeff  = Math::abs( S_ptr->coeff( j ) );

                    // transform the real valued coefficients to values of type 'weight_t'
                    //           │    γ · |a_ij|   │
                    // weight =  │ ─────────────── │, with γ = 80.000 (see transform_lin)
                    //           │ |a_ii| + |a_jj| │
                    //
                    const weight_t  weight = transform_lin< weight_t, real_t >( Math::abs( coeff / ( coeff_node + coeff_neigh ) ), 1.0 );
                    // const weight_t  weight = transform_log< weight_t, real_t >( coeff, min_val );
            
                    adj_nodes()[ idx ] = lneigh;
                    set_edge_weight( idx, weight );
                    degrees[lnode]++;
            
                    if ( ! S_ptr->has_entry( neigh, node ) )
                    {
                        const idx_t  nidx = adj_list_ptr()[lneigh] + idx_t(degrees[lneigh]);

                        adj_nodes()[ nidx ] = lnode;
                        if ( sym_edge_weights )
                            set_edge_weight( nidx, weight );      // symmetric weight
                        else
                            set_edge_weight( nidx, weight_t(0) ); // unsymmetric weight
                        degrees[lneigh]++;
                    }// if
                    else if ( sym_edge_weights )
                    {
                        // use maximal weight of both edges
                        const auto      ncoeff  = Math::abs( S_ptr->entry( neigh, node ) );
                        const weight_t  nweight = transform_log< weight_t, real_t >( ncoeff, min_val );
                
                        set_edge_weight( idx, std::max( weight, nweight ) );
                    }// if
                }// while
            }// for

            // add global names
            for ( size_t  i = 0; i < n_nodes; i++ )
            {
                if ( has_mask && node_mask[i] )
                    continue;
        
                global_name()[name_map[i]] = node_t(i);
            }// for
        },
        S );
}

//
// copy operator
//
TEWGraph &
TEWGraph::operator = ( const TEWGraph & graph )
{
    _adj_list_ptr = graph._adj_list_ptr;
    _adj_nodes    = graph._adj_nodes;
    _global_name  = graph._global_name;
    _edge_weights = graph._edge_weights;

    return *this;
}


//
// return the minimal absolute value of the edge weights 
// of the graph which is bigger than zero
//
weight_t
TEWGraph::min_edge_weight ( ) const
{
     const size_t  n_edges = nedges();
  
     if (( n_edges < 1 ) || ! has_edge_weights() )
        HERROR( ERR_CONSISTENCY, "(TEWGraph) min_edge_weight", "no edge weight in TGraph" );
     
     // initialise min_weight with a correct value
     weight_t  min_weight = 0;
     
     for ( size_t  i = 0; i < n_edges; i++ )
     {
         const weight_t  weight = edge_weight( idx_t( i ) );
         
         if ( weight > 0 )
         {
             min_weight = weight;
             break;
         }// if
     }// for
     
     // compute min_weight
     for ( size_t  i = 0; i < n_edges; i++ )
     {
         const weight_t  weight = edge_weight( idx_t( i ) );

         if ( (weight < min_weight) && (weight > 0) )
             min_weight = weight;
     }// for
     
     if ( min_weight == 0 )
        HERROR( ERR_CONSISTENCY,  "(TEWGraph) min_edge_weight", "no minimal edge weight bigger than zero found" );
     
     return min_weight; 
}

//
// write in Chaco/Jostle/Metis file format
//
void
TEWGraph::write ( const std::string &  filename ) const
{
    std::ofstream  file( filename.c_str() );

    // first line: number of nodes and edges (counted only in one direction)
    file << nnodes() << ' ' << nedges() / 2 << " 1" << std::endl;
    
    // now, <n_nodes> lines with adjacent nodes
    for ( auto  node : nodes() )
    {
        for ( auto  neigh_weight : adj_nodes_weights( node ) )
        {
            const node_t  neigh = neigh_weight.first;
            
            // node names numbered from 1..n_nodes
            file << ' ' << neigh+1 << ' ' << neigh_weight.second << ' ';
        }// for

        file << std::endl;
    }// for
}

}// namespace Hpro
