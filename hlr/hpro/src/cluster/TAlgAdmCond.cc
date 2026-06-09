//
// Project     : HLIBpro
// File        : TAlgAdmCond.cc
// Description : algebraic admissibility condition for sparse matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <deque>
#include <unordered_map>

#include "list.hh"

#include "hpro/base/error.hh"

#include "hpro/cluster/TAlgAdmCond.hh"

namespace Hpro
{

//////////////////////////////////////////////
//
// local defines
//

// enables output
//#define PRINT  if ( false )

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
// TAlgAdmCond
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////
//
// constructor and destructor
//

TAlgAdmCond::TAlgAdmCond ( any_const_sparse_matrix_t  S,
                           const TPermutation *       perm_i2e )
{
    _mat          = S;
    _row_perm_i2e = perm_i2e;
    _col_perm_i2e = perm_i2e;

    std::visit( [] ( auto && S_ptr ) { if ( S_ptr == nullptr ) HERROR( ERR_ARG, "(TAlgAdmCond)", "sparse matrix is nullptr" ); }, _mat );

    if ( perm_i2e != nullptr )
    {
        _row_perm_e2i = new TPermutation( * _row_perm_i2e );
        _row_perm_e2i->invert();
    
        _col_perm_e2i = new TPermutation( * _col_perm_i2e );
        _col_perm_e2i->invert();
    }// if
    else
    {
        _row_perm_e2i = nullptr;
        _col_perm_e2i = nullptr;
    }// else
}

TAlgAdmCond::TAlgAdmCond ( any_const_sparse_matrix_t  S,
                           const TPermutation *       row_perm_i2e,
                           const TPermutation *       col_perm_i2e )
{
    _mat          = S;
    _row_perm_i2e = row_perm_i2e;
    _col_perm_i2e = col_perm_i2e;

    std::visit( [] ( auto && S_ptr ) { if ( S_ptr == nullptr ) HERROR( ERR_ARG, "(TAlgAdmCond)", "sparse matrix is nullptr" ); }, _mat );
    
    if ( row_perm_i2e != nullptr )
    {
        _row_perm_e2i = new TPermutation( * _row_perm_i2e );
        _row_perm_e2i->invert();
    }// if
    else
        _row_perm_e2i = nullptr;
        
    if ( col_perm_i2e != nullptr )
    {
        _col_perm_e2i = new TPermutation( * _col_perm_i2e );
        _col_perm_e2i->invert();
    }// if
    else
        _col_perm_e2i = nullptr;
}

TAlgAdmCond::~TAlgAdmCond ()
{
    delete _row_perm_e2i;
    delete _col_perm_e2i;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
// TStdAlgAdmCond
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////
//
// constructor and destructor
//

TStdAlgAdmCond::TStdAlgAdmCond ( const double               eta,
                                 any_const_sparse_matrix_t  S,
                                 const TPermutation *       perm_i2e )
        : TAlgAdmCond( S, perm_i2e )
        , _eta(eta)
{
    std::visit( [this] ( auto && S_ptr ) { _visited.resize( S_ptr->rows() ); }, S );
}

TStdAlgAdmCond::TStdAlgAdmCond ( const double               eta,
                                 any_const_sparse_matrix_t  S,
                                 const TPermutation *       row_perm_i2e,
                                 const TPermutation *       col_perm_i2e )
        : TAlgAdmCond( S, row_perm_i2e, col_perm_i2e )
        , _eta(eta)
{
    std::visit( [this] ( auto && S_ptr ) { _visited.resize( S_ptr->rows() ); }, S );
}

///////////////////////////////////////////
//
// check block-cluster if admissible
//

bool
TStdAlgAdmCond::is_adm ( const TBlockCluster * c ) const
{
    const TCluster  * rowcl, * colcl;

    rowcl = c->rowcl();
    colcl = c->colcl();

    if ( rowcl == colcl )
        return false;

    std::visit( [] ( auto && S ) { if ( S == nullptr ) HERROR( ERR_ARG, "(TStdAlgAdmCond)", "sparse matrix is nullptr" ); }, _mat );

    //
    // in case of nested dissection, offdiagonal, pure-domain clusters are admissible
    //

    if ( rowcl->is_domain() && colcl->is_domain() )
        return true;

    //
    // determine diameter of both clusters
    //

    const uint  diam_rowcl = diameter( rowcl, _row_perm_i2e, _row_perm_e2i );
    const uint  diam_colcl = diameter( colcl, _row_perm_i2e, _col_perm_e2i );

    //
    // now define minimal path-length between both sets
    // and compare distance between clusters with it
    //
    
    uint  min_len = uint( std::ceil( std::min( double( diam_rowcl ), double( diam_colcl ) ) / _eta ) );

    return cmp_dist( rowcl, colcl, std::max<uint>( min_len, 1 ) );
}

//
// determine diameter of a cluster
//
uint
TStdAlgAdmCond::diameter ( const TCluster *      c,
                           const TPermutation *  perm_i2e,
                           const TPermutation *  perm_e2i ) const
{
    //
    // look for two nodes in graph with maximal distance
    //

    uint            max_depth = 0;
    TNodeSet        left( c->size() ), right( c->size() );

    // start with first node in cluster
    if ( perm_i2e != nullptr ) right.append( perm_i2e->permute( c->first() ) );
    else                       right.append( c->first() );
    
    while ( true )
    {
        //
        // copy all local nodes from right to left set (as new start)
        //
        
        left.remove_all();
        
        for ( auto  node : right )
        {
            if ( is_local( c, node, perm_e2i ) )
                left.append( node );
        }// for

        //
        // do BFS in matrix graph with given start node
        //

        const uint depth = bfs( left, right, c, perm_i2e, perm_e2i );

        //
        // if no longer depth could be found, stop iteration
        //

        if      ( depth >  max_depth ) max_depth = depth;
        else if ( depth <= max_depth ) break;
    }// while

    return 2 * max_depth;
}

//
// perform a BFS from set <start> in matrix and store
// last visited nodes in <last>; stop BFS if all locally
// marked nodes have been visited; return the depth
// of the BFS iteration
//
uint
TStdAlgAdmCond::bfs ( TNodeSet &            start,
                      TNodeSet &            last,
                      const TCluster *      tau,
                      const TPermutation *  perm_i2e,
                      const TPermutation *  perm_e2i ) const
{
    uint      depth = 0;
    TNodeSet  succ(  tau->size() );
    TNodeSet  nodes( tau->size() );
    uint      n_loc_visited = 0;
    auto      rowptr = std::visit( [] ( auto && S ) { return S->rowptr(); }, _mat );
    auto      colind = std::visit( [] ( auto && S ) { return S->colind(); }, _mat );

    for ( auto  node : start )
    {
        _visited[ node ] = true;

        // count local nodes in start set
        if ( is_local( tau, node, perm_e2i ) )
            n_loc_visited++;
    }// for
    
    nodes = start;
    while (( n_loc_visited < tau->size() ) && ( nodes.nnodes() > 0 ))
    {
        last = nodes;

        // bfs_step( graph, nodes, succ, visited );
        {
            succ.remove_all();
    
            for ( auto  node : nodes )
            {
                //
                // visit successors of node
                //

                const idx_t  lb = rowptr[ node ];
                const idx_t  ub = rowptr[ node+1 ];
        
                for ( idx_t  j = lb; j < ub; j++ )
                {
                    const node_t  neigh  = colind[ j ];
                    const bool    is_loc = is_local( tau, neigh, _col_perm_e2i ); // WARNING: different row
                    // and column perm. may
                    // lead to problems

                    if (( neigh != node ) && is_loc && ! _visited[ neigh ] )
                    {
                        _visited[ neigh ] = true;
                        if ( is_loc )
                            n_loc_visited++;
                        succ.append( neigh );
                    }// if
                }// for
            }// for

            nodes.remove_all();
        }
        
        nodes = succ;
        succ.remove_all();
        depth++;
    }// while

    // unmark visited nodes
    for ( idx_t  i = tau->first(); i <= tau->last(); i++ )
    {
        const node_t  node = ( perm_i2e != nullptr ? perm_i2e->permute( i ) : i );

        _visited[ node ] = false;
    }// for
    
    return depth-1;
}

//
// return true, if distance between two clusters is
// bigger than min_dist
//
bool
TStdAlgAdmCond::cmp_dist ( const TCluster *  tau,
                           const TCluster *  sigma,
                           const uint        min_dist ) const
{
    std::list< node_t >  visnodes;
    std::list< node_t >  left;
    std::list< node_t >  succ;

    // collect nodes in <tau> in start set for BFS iteration
    for ( idx_t  i = tau->first(); i <= tau->last(); i++ )
    {
        const node_t  node = ( _row_perm_i2e != nullptr ? _row_perm_i2e->permute( i ) : i );

        left.push_back( node );
        _visited[ node ] = true;
        visnodes.push_back( node );
    }// for

    //
    // do a BFS from <tau> up to <min_dist> steps and check
    // if a node in <sigma> was reached
    //

    auto  rowptr = std::visit( [] ( auto && S ) { return S->rowptr(); }, _mat );
    auto  colind = std::visit( [] ( auto && S ) { return S->colind(); }, _mat );
    
    for ( uint i = 0; i < min_dist; i++ )
    {
        //
        // determine successors of <left> set
        //

        while ( ! left.empty() )
        {
            const node_t  node = behead( left );
            const idx_t   lb   = rowptr[ node ];
            const idx_t   ub   = rowptr[ node+1 ];

            for ( idx_t  l = lb; l < ub; l++ )
            {
                const node_t  neigh = colind[ l ];

                if ( ! _visited[ neigh ] )
                {
                    // return immediately, if neighbour is in <sigma>
                    if ( is_local( sigma, neigh, _col_perm_e2i ) )
                    {
                        // unmark visited nodes
                        while ( visnodes.size() > 0 )
                            _visited[ behead( visnodes ) ] = false;
                        
                        return false;
                    }// if

                    // otherwise remember neighbour for next iteration step
                    succ.push_back( neigh );
                    _visited[neigh] = true;
                }// if
            }// for
        }// for

        // next iteration with successor nodes
        left = succ;
        succ.clear();
    }// for

    // unmark visited nodes
    while ( ! visnodes.empty() )
        _visited[ behead( visnodes ) ] = false;
    
    return true;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
// TWeakAlgAdmCond
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////
//
// constructor and destructor
//

TWeakAlgAdmCond::TWeakAlgAdmCond ( any_const_sparse_matrix_t  S,
                                   const TPermutation *       perm_i2e,
                                   const uint                 distance,
                                   const uint                 connectivity )
        : TAlgAdmCond( S, perm_i2e )
        , _distance( std::max< uint >( distance, 1 ) )
        , _connectivity( connectivity )
{
}

TWeakAlgAdmCond::TWeakAlgAdmCond ( any_const_sparse_matrix_t  S,
                                   const TPermutation *       row_perm_i2e,
                                   const TPermutation *       col_perm_i2e,
                                   const uint                 distance,
                                   const uint                 connectivity )
        : TAlgAdmCond( S, row_perm_i2e, col_perm_i2e )
        , _distance( std::max< uint >( distance, 1 ) )
        , _connectivity( connectivity )
{
}

///////////////////////////////////////////
//
// check block-cluster if admissible
//

namespace
{

//
// check if there is a path of length <path_length> from <rowcl> to <colcl> in <S>
//
template < typename value_t >
bool
has_path ( const TCluster *                  rowcl,
           const TCluster *                  colcl,
           const TSparseMatrix< value_t > *  S,
           const uint                        path_length,
           const TPermutation *              row_perm,
           const TPermutation *              col_perm )
{
    std::deque< idx_t >                nodes;
    std::unordered_map< idx_t, bool >  visited;

    auto  was_visited  = [&visited] ( const idx_t  idx ) { return ( visited.find( idx ) != visited.end() ); };
    auto  row_perm_i2e = [row_perm] ( const idx_t  idx ) { return ( row_perm != nullptr ? row_perm->permute( idx ) : idx ); };
    auto  col_perm_e2i = [col_perm] ( const idx_t  idx ) { return ( col_perm != nullptr ? col_perm->permute( idx ) : idx ); };
    
    for ( auto  rowidx : * rowcl )
    {
        nodes.push_back( rowidx );
        visited[ rowidx ] = true;
    }// for
    
    for ( uint  length = 0; length < path_length; ++length )
    {
        std::deque< idx_t >  succ;
        
        for ( auto  rowidx : nodes )
        {
            const idx_t  row = row_perm_i2e( rowidx );
            const idx_t  lb  = S->rowptr( row );
            const idx_t  ub  = S->rowptr( row + 1 );
            
            for ( idx_t  j = lb; j < ub; j++ )
            {
                const idx_t  colidx = col_perm_e2i( S->colind( j ) );
                
                if ( colcl->is_in( colidx ) )
                    return true;

                if (( length != path_length-1 ) && ! was_visited( colidx ) )
                {
                    succ.push_back( colidx );
                    visited[ colidx ] = true;
                }// if
            }// for
        }// for

        if ( length != path_length-1 )
            nodes = succ;
    }// for

    return false;
}

//
// check if there is an edge connecting <rowcl> with <colcl> in <S>
//
template < typename value_t >
bool
has_edge ( const TCluster *                  rowcl,
           const TCluster *                  colcl,
           const TSparseMatrix< value_t > *  S,
           const TPermutation *              row_perm,
           const TPermutation *              col_perm )
{
    auto  row_perm_i2e = [row_perm] ( const idx_t  idx ) { return ( row_perm != nullptr ? row_perm->permute( idx ) : idx ); };
    auto  col_perm_e2i = [col_perm] ( const idx_t  idx ) { return ( col_perm != nullptr ? col_perm->permute( idx ) : idx ); };
    
    for ( auto  rowidx : * rowcl )
    {
        const idx_t  row = row_perm_i2e( rowidx );
        const idx_t  lb  = S->rowptr( row );
        const idx_t  ub  = S->rowptr( row + 1 );
        
        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  colidx = col_perm_e2i( S->colind( j ) );
            
            if ( colcl->is_in( colidx ) )
                return true;
        }// for
    }// for

    return false;
}

//
// count number of direct connections between rowcl and colcl
//
template < typename value_t >
uint
count_edges ( const TCluster *                  rowcl,
              const TCluster *                  colcl,
              const TSparseMatrix< value_t > *  S,
              const TPermutation *              row_perm,
              const TPermutation *              col_perm )
{
    auto  row_perm_i2e = [row_perm] ( const idx_t  idx ) { return ( row_perm != nullptr ? row_perm->permute( idx ) : idx ); };
    auto  col_perm_e2i = [col_perm] ( const idx_t  idx ) { return ( col_perm != nullptr ? col_perm->permute( idx ) : idx ); };
    uint  nedges       = 0;
    
    for ( auto  rowidx : * rowcl )
    {
        const idx_t  row = row_perm_i2e( rowidx );
        const idx_t  lb  = S->rowptr( row );
        const idx_t  ub  = S->rowptr( row + 1 );
        
        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  colidx = col_perm_e2i( S->colind( j ) );
            
            if ( colcl->is_in( colidx ) )
                nedges++;
        }// for
    }// for

    return nedges;
}

}// namespace anonymous

bool
TWeakAlgAdmCond::is_adm ( const TBlockCluster * c ) const
{
    const TCluster  * rowcl, * colcl;

    rowcl = c->rowcl();
    colcl = c->colcl();

    if ( rowcl == colcl )
        return false;

    std::visit( [] ( auto && S ) { if ( S == nullptr ) HERROR( ERR_ARG, "(TWeakAlgAdmCond)", "sparse matrix is nullptr" ); }, _mat );

    //
    // in case of nested dissection, offdiagonal, pure-domain clusters are admissible
    //

    if ( rowcl->is_domain() && colcl->is_domain() )
        return true;

    const bool  is_unsymmetric = std::visit( [] ( auto && S ) { return S->is_unsymmetric(); }, _mat );
    
    if (( _distance > 1 ) || ( _connectivity == 0 ))
    {
        //
        // test for path between rowcl and colcl
        // - check both directions in case of unsymmetric matrix
        //

        bool  row_to_col = false;
        bool  col_to_row = false;
    
        if ( _distance == 1 )
            row_to_col = std::visit( [=] ( auto && S ) { return has_edge( rowcl, colcl, S,            _row_perm_i2e, _col_perm_e2i ); }, _mat );
        else
            row_to_col = std::visit( [=] ( auto && S ) { return has_path( rowcl, colcl, S, _distance, _row_perm_i2e, _col_perm_e2i ); }, _mat );
        
        if ( row_to_col )
            return false;
        
        if ( is_unsymmetric )
        {
            if ( _distance == 1 )
                col_to_row = std::visit( [=] ( auto && S ) { return has_edge( colcl, rowcl, S,            _col_perm_i2e, _row_perm_e2i ); }, _mat );
            else
                col_to_row = std::visit( [=] ( auto && S ) { return has_path( colcl, rowcl, S, _distance, _col_perm_i2e, _row_perm_e2i ); }, _mat );
        }// if

        if ( row_to_col || col_to_row )
            return false;
    
        return true;
    }// if
    else
    {
        //
        // count number of edges between clusters (both directions)
        //

        uint  nrow2col = 0;
        uint  ncol2row = 0;

        nrow2col = std::visit( [=] ( auto && S ) { return count_edges( rowcl, colcl, S, _row_perm_i2e, _col_perm_e2i ); }, _mat );
        ncol2row = std::visit( [=] ( auto && S ) { return count_edges( colcl, rowcl, S, _col_perm_i2e, _row_perm_e2i ); }, _mat );

        if ( std::max( nrow2col, ncol2row ) > 2 )
            return false;
    
        return true;
    }// else
}

}// namespace Hpro
