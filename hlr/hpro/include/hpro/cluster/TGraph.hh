#ifndef __HPRO_TGRAPH_HH
#define __HPRO_TGRAPH_HH
//
// Project     : HLIBpro
// File        : TGraph.hh
// Description : classes for representing graphs
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>
#include <list>

#include "hpro/cluster/TNodeSet.hh"
#include "hpro/matrix/TSparseMatrix.hh"

namespace Hpro
{

//!
//! \ingroup  Cluster_Module
//! \typedef  weight_t
//! \brief    type for edge weights
//!
using  weight_t = int;   // only 'real' or 'int' should be used. 'unsigned int' will produce error 
    
    
//!
//! \ingroup  Cluster_Module
//! \class    TGraph
//! \brief    Class for a undirected graph stored as adjacency matrix in CRS representation.
//!
class TGraph
{
public:
    //!
    //! class representing node set in graph
    //!
    struct TNodes
    {
    public:
        struct iterator
        {
        private:
            const TGraph *  _graph;
            idx_t           _pos;

        public:
            iterator ( const TGraph *  graph,
                       const idx_t     pos )
                    : _graph( graph )
                    , _pos( pos )
            {}

            iterator  operator ++ ()    noexcept { return iterator( _graph, ++_pos ); }
            iterator  operator ++ (int) noexcept { return iterator( _graph, _pos++ ); }

            iterator  operator -- ()    noexcept { return iterator( _graph, --_pos ); }
            iterator  operator -- (int) noexcept { return iterator( _graph, _pos-- ); }

            node_t    operator *  () const noexcept { return _pos; }

            bool      operator == ( const iterator &  it ) const noexcept
            {
                return (_pos == it._pos) && (_graph == it._graph);
            }

            bool      operator != ( const iterator &  it ) const noexcept
            {
                return (_pos != it._pos) || (_graph != it._graph);
            }
        };
        
        using  const_iterator = iterator;
        
    private:
        const TGraph *  _graph;

    public:
        TNodes ( const TGraph *  graph )
                : _graph( graph )
        {
            if ( _graph == nullptr )
                HERROR( ERR_ARG, "(TNodes) ctor", "graph is nullptr" );
        }

        iterator        begin  () const noexcept { return iterator( _graph, 0 ); }
        iterator        end    () const noexcept { return iterator( _graph, idx_t(_graph->nnodes()) ); }
        
        const_iterator  cbegin () const noexcept { return const_iterator( _graph, 0 ); }
        const_iterator  cend   () const noexcept { return const_iterator( _graph, idx_t(_graph->nnodes()) ); }
    };
    
    //!
    //! represents adjacent nodes in a graph
    //!
    class TAdjNodes
    {
    public:
        struct iterator
        {
        private:
            const TGraph *  _graph;
            idx_t           _pos;

        public:
            iterator ( const TGraph *  graph,
                       const idx_t     pos )
                    : _graph( graph )
                    , _pos( pos )
            {}

            iterator  operator ++ ()    noexcept { return iterator( _graph, ++_pos ); }
            iterator  operator ++ (int) noexcept { return iterator( _graph, _pos++ ); }

            iterator  operator -- ()    noexcept { return iterator( _graph, --_pos ); }
            iterator  operator -- (int) noexcept { return iterator( _graph, _pos-- ); }

            node_t    operator *  () const noexcept { return _graph->_adj_nodes[ _pos ]; }

            bool      operator == ( const iterator &  it ) const noexcept
            {
                return (_pos == it._pos) && (_graph == it._graph);
            }

            bool      operator != ( const iterator &  it ) const noexcept
            {
                return (_pos != it._pos) || (_graph != it._graph);
            }

            idx_t     position    () const noexcept { return _pos; }
        };

        using  const_iterator = iterator;
        
    private:
        //! graph we address
        const TGraph *  _graph;

        // node to which other nodes are adjacent
        const node_t    _node;
        
    public:
        //!
        //! ctor: init. iterator
        //!
        TAdjNodes ( const TGraph *  graph,
                    const node_t    node )
                : _graph( graph )
                , _node( node )
        {
            if ( _graph == nullptr )
                HERROR( ERR_ARG, "(TAdjNodes) ctor", "graph is nullptr" );
        }

        iterator        begin  () const noexcept { return iterator( _graph, _graph->_adj_list_ptr[ _node ] ); }
        iterator        end    () const noexcept { return iterator( _graph, _graph->_adj_list_ptr[ _node+1 ] ); }
        
        const_iterator  cbegin () const noexcept { return const_iterator( _graph, _graph->_adj_list_ptr[ _node ] ); }
        const_iterator  cend   () const noexcept { return const_iterator( _graph, _graph->_adj_list_ptr[ _node+1 ] ); }
    };

    //!
    //! represents adjacent nodes with edge weights in a graph
    //!
    class TAdjNodesWeights
    {
    public:
        struct iterator
        {
        private:
            const TGraph *  _graph;
            idx_t           _pos;

        public:
            iterator ( const TGraph *  graph,
                       idx_t           pos )
                    : _graph( graph )
                    , _pos( pos )
            {}

            iterator  operator ++ ()    noexcept { return iterator( _graph, ++_pos ); }
            iterator  operator ++ (int) noexcept { return iterator( _graph, _pos++ ); }

            iterator  operator -- ()    noexcept { return iterator( _graph, --_pos ); }
            iterator  operator -- (int) noexcept { return iterator( _graph, _pos-- ); }

            auto      operator *  () const noexcept -> std::pair< node_t, weight_t >
            {
                return std::pair< node_t, weight_t >( _graph->_adj_nodes[ _pos ],
                                                      _graph->edge_weight( _pos ) );
            }

            bool      operator == ( const iterator &  it ) const noexcept
            {
                return (_pos == it._pos) && (_graph == it._graph);
            }

            bool      operator != ( const iterator &  it ) const noexcept
            {
                return (_pos != it._pos) || (_graph != it._graph);
            }

            idx_t     position    () const noexcept { return _pos; }
        };

        using  const_iterator = iterator;
        
    private:
        //! graph we address
        const TGraph *  _graph;

        // node to which other nodes are adjacent
        const node_t    _node;
        
    public:
        //!
        //! ctor: init. iterator
        //!
        TAdjNodesWeights ( const TGraph *  graph,
                           const node_t    node )
                : _graph( graph )
                , _node( node )
        {
            if ( _graph == nullptr )
                HERROR( ERR_ARG, "(TAdjNodesWeights) ctor", "graph is nullptr" );
        }

        iterator        begin  () const noexcept { return iterator( _graph, _graph->_adj_list_ptr[ _node ] ); }
        iterator        end    () const noexcept { return iterator( _graph, _graph->_adj_list_ptr[ _node+1 ] ); }
        
        const_iterator  cbegin () const noexcept { return const_iterator( _graph, _graph->_adj_list_ptr[ _node ] ); }
        const_iterator  cend   () const noexcept { return const_iterator( _graph, _graph->_adj_list_ptr[ _node+1 ] ); }
    };
    
protected:
    //! contains pointers into adjacency list for each node
    std::vector< idx_t >   _adj_list_ptr;
    
    //! contains list of adjacent nodes for each node
    std::vector< node_t >  _adj_nodes;

    //! global name of local nodes
    std::vector< node_t >  _global_name;

public:
    //////////////////////////////////////////////
    //
    // constructor
    //

    //! construct empty graph
    TGraph () {}

    //! construct graph using symmetrised sparsity pattern of matrix \a S;
    //! nodes marked in \a nodemask, i.e. nodemask[node] = true, will be excluded
    //! from graph
    TGraph ( any_const_sparse_matrix_t    S,
             const std::vector< bool > &  nodemask );

    //! dtor
    virtual ~TGraph () {}
    
    //////////////////////////////////////////////
    //
    // access graph data
    //

    //! return number of nodes in graph
    size_t     nnodes            () const noexcept { return _adj_list_ptr.size() > 0 ? _adj_list_ptr.size() - 1 : 0; }
    
    //! return number of edges in graph
    size_t     nedges            () const noexcept { return _adj_nodes.size(); }

    //! return true, if graph stores global names
    bool       has_global_names  () const noexcept { return _global_name.size() > 0; }

    
    //! return degree of \a node, e.g. number of incident edges
    size_t     degree            ( const node_t  node ) const noexcept
    {
        return _adj_list_ptr[node+1] - _adj_list_ptr[node];
    }

    //! return set of graph nodes
    TNodes     nodes             () const noexcept
    {
        return TNodes( this );
    }
    
    //! return set of adjacent nodes in graph
    TAdjNodes  adj_nodes         ( const node_t  node ) const noexcept
    {
        return TAdjNodes( this, node );
    }

    //! return set of adjacent nodes with edge weights in graph
    TAdjNodesWeights  adj_nodes_weights ( const node_t  node ) const noexcept
    {
        return TAdjNodesWeights( this, node );
    }

    
    //! initialise graph for \a nnnodes nodes and \a nedges edges
    //! - if \a global is true, global names will also be initialised
    virtual void init ( const size_t  annodes,
                        const size_t  anedges,
                        const bool    aglobal = false )
    {
        _adj_nodes.resize( anedges );
        _adj_list_ptr.resize( annodes + 1 );
        if ( aglobal )
            _global_name.resize( annodes );
    }
    
    

    //////////////////////////////////////////////
    //
    // give access to the internal data structures
    //

    // //! return adjacency data
    std::vector< node_t > &        adj_nodes    ()       noexcept { return _adj_nodes; }
    const std::vector< node_t > &  adj_nodes    () const noexcept { return _adj_nodes; }

    // //! return adjacency pointer data
    std::vector< idx_t > &         adj_list_ptr ()       noexcept { return _adj_list_ptr; }
    const std::vector< idx_t > &   adj_list_ptr () const noexcept { return _adj_list_ptr; }

    //! return global name data
    std::vector< node_t > &        global_name  ()       noexcept { return _global_name; }
    const std::vector< node_t > &  global_name  () const noexcept { return _global_name; }

    //////////////////////////////////////////////
    //
    // various graph algorithms
    //

    //! return node with minimal degree in given graph
    node_t  min_degree_node () const;

    //! return node with maximal degree in given graph
    node_t  max_degree_node () const;

    //! compute strongly connected components and build graph for each component
    // void build_scc ( std::list< TGraph > & scc ) const;

    //! compute strongly connected components but restrict computation to nodes
    //! marked by \a mark, where node marks are stored in \a label
    void build_scc ( std::list< TNodeSet > &      scc,
                     const std::vector< uint > &  label,
                     const uint                   mark ) const;

    //! compute strongly connected components and build nodeset for each component
    void build_scc ( std::list< TNodeSet > &      scc ) const;

    //! restrict graph to nodes in \a nodes and return resulting subgraph
    auto restrict  ( const TNodeSet &  nodes ) const -> std::unique_ptr< TGraph >;

    //! compute vertex separator between subgraphs \a left and \a right
    //! (also defined by \a label) and put nodes into \a vertex_sep
    void vertex_separator ( std::vector< uint > &  label,
                            const TNodeSet &       left,
                            const TNodeSet &       right,
                            TNodeSet &             vertex_sep,
                            const uint             if_label ) const;
    
    //////////////////////////////////////////////
    //
    // edge weight management
    //

    //! return false to show that Graph has no edge weights
    virtual bool      has_edge_weights  () const noexcept { return false; }
    
    //! return edge weight for \a i'th edge
    virtual weight_t  edge_weight       ( const idx_t ) const noexcept { return weight_t( 1 ); }
    
    //! set edge weight for \a i'th edge
    virtual void      set_edge_weight   ( const idx_t,
                                          const weight_t ) 
    {
        HERROR( ERR_CONSISTENCY, "(TGraph) set_edge_weight", "TGraph has no edge weights" );
    }    
    
    //! return the minimal absolute value of the edge weights of TEWGraph
    virtual weight_t  min_edge_weight   () const
    {
        HERROR( ERR_CONSISTENCY, "(TGraph) min_edge_weight", "TGraph has no edge weights" );
        return weight_t( 1 );
    }
    
    //////////////////////////////////////////////
    //
    // output
    //

    //! print graph in GraphViz format to \a filename
    void print ( const std::string &          filename,
                 const bool                   global = false ) const;

    //! print graph in GraphViz format to \a filename with labels in \a label
    void print ( const std::string &          filename ,
                 const std::vector< uint > &  label,
                 const bool                   global = false ) const;

    //! export graph in Chaco/Jostle/Metis format to \a filename
    virtual void write ( const std::string &  filename ) const;          
    
    //////////////////////////////////////////////
    //
    // misc.
    //

    //! copy operator
    TGraph & operator = ( const TGraph & graph );

    //! virtual constructor, e.g. return object of same type
    virtual TGraph *  create () const { return new TGraph; }
    
};

//!
//! \ingroup  Cluster_Module
//! \class    TEWGraph
//! \brief    Represents undirected graph with edge weights.
//!
class TEWGraph : public TGraph
{
private:
    //! edge weights
    std::vector< weight_t >   _edge_weights;

public:
    //////////////////////////////////////////////
    //
    // constructor
    //

    //! construct empty graph
    TEWGraph () {}
   
    //! construct graph using symmetrised sparsity pattern of matrix \a S;
    //! nodes marked in \a nodemask, i.e. nodemask[node] = true, will be excluded
    //! from graph; matrix coefficients will be used for edge weights
    TEWGraph ( any_const_sparse_matrix_t    S,
               const std::vector< bool > &  nodemask,
               const bool                   sym_edge_weights );
    
    //! construct an edge weighted graph using a normal graph 
    //! and initialising edge weights with 1
    TEWGraph ( const TGraph &  graph ) 
    {
        _adj_list_ptr = graph.adj_list_ptr();
        _adj_nodes    = graph.adj_nodes();
        _global_name  = graph.global_name();
        _edge_weights.resize( graph.nedges(), weight_t( 1 ) );
    }  

    //! dtor
    virtual ~TEWGraph () {}

    //////////////////////////////////////////////
    //
    // give access to the internal data structures
    //

    //! initialise graph for \a nnnodes nodes and \a nedges edges
    //! - if \a global is true, global names will also be initialised
    virtual void init ( const size_t  annodes,
                        const size_t  anedges,
                        const bool    aglobal = false )
    {
        TGraph::init( annodes, anedges, aglobal );

        // initialise edge weights with standard weight
        _edge_weights.resize( anedges, weight_t( 1 ) );
    }

    //////////////////////////////////////////////
    //
    // edge weight management
    //

    //! return edge weight for \a i'th edge
    virtual weight_t  edge_weight       ( const idx_t  i ) const noexcept
    {
        return std::abs( _edge_weights[i] );
    }
    
    //! set edge weight for \a i'th edge
    virtual void      set_edge_weight   ( const idx_t     i,
                                          const weight_t  w )
    {
        _edge_weights[i] = w;
    }

    //! return false to show that Graph has no edge weights
    virtual bool      has_edge_weights  () const noexcept { return true; }
    
    //! return the minimal absolute value of the edge weights of TEWGraph which is bigger than zero
    virtual weight_t  min_edge_weight   () const;

    //////////////////////////////////////////////
    //
    // output
    //

    //! export graph in Chaco/Jostle/Metis format to \a filename
    virtual void write ( const std::string &  filename ) const;          
    
    //////////////////////////////////////////////
    //
    // misc.
    //

    //! copy operator
    TEWGraph & operator = ( const TEWGraph & graph );
    
    //! virtual constructor, e.g. return object of same type
    virtual TGraph *  create () const { return new TEWGraph; }
};

}// namespace Hpro

#endif // __HPRO_TGRAPH_HH
