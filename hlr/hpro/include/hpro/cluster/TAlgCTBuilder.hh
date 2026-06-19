#ifndef __HPRO_TALGCTBUILDER_HH
#define __HPRO_TALGCTBUILDER_HH
//
// Project     : HLIBpro
// File        : TAlgCTBuilder.hh
// Description : class for algebraic clustertree construction
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>
#include <atomic>

#include "hpro/cluster/types.hh"
#include "hpro/cluster/TCluster.hh"
#include "hpro/cluster/TClusterTree.hh"
#include "hpro/cluster/TGraph.hh"
#include "hpro/cluster/TAlgPartStrat.hh"
#include "hpro/matrix/TSparseMatrix.hh"

namespace Hpro
{

//!
//! \ingroup  Cluster_Module
//! \class    TAlgCTBuilder
//! \brief    Base class for cluster tree construction algorithms based
//!           on graph partitioning with graph defined by a sparse matrix.
//!
class TAlgCTBuilder
{
protected:
    // graph partitioning strategy
    TAlgPartStrat *      _part_strat;
    
    // minimal size of a cluster, i.e. not smaller than this
    const uint           _n_min;

    // minimal level for leaves, i.e. no leaves on a level less than this
    const uint           _min_leaf_lvl;

    // factor for determining high degree nodes (≥ high_deg_fac · avg_degree);
    // if 0, no separation is performed
    uint                 _high_deg_fac;
    
    // defines usage of edge weights for graph partitioning 
    edge_weights_mode_t  _edge_weights_mode;
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct cluster tree builder with given partition strategy and tree parameters
    TAlgCTBuilder ( TAlgPartStrat *  part_strat,
                    const uint       n_min        = CFG::Cluster::nmin,
                    const uint       min_leaf_lvl = 0 );

    //! dtor
    virtual ~TAlgCTBuilder () {}

    //////////////////////////////////////////////
    //
    // options
    //

    //! activate/deactivate (if \a fac = 0) high degree node separation
    void  set_high_deg_fac      ( const uint  fac );
    
    //! set mode for edge weights
    void  set_edge_weights_mode ( const edge_weights_mode_t  edge_weights_mode );

    //////////////////////////////////////////////
    //
    // clustertree creation
    //

    //! build cluster tree using connectivity in sparse matrix \a S
    std::unique_ptr< TClusterTree >
    build   ( any_const_sparse_matrix_t  S,
              const idx_t                idx_ofs = 0 ) const;

    //! divide graph \a graph and build corresponding cluster tree
    virtual std::unique_ptr< TCluster >
    divide  ( const TGraph &             graph,
              const uint                 lvl,
              TPermutation &             perm,
              const idx_t                idx_ofs,
              const uint                 n_min,
              any_const_sparse_matrix_t  S,
              const uint                 max_lvl,
              std::atomic< int > &       id ) const;

    //! compute graph bi-partitioning of \a graph and store result in \a left and \a right
    virtual void           partition     ( const TGraph &        graph,
                                           TNodeSet &            left,
                                           TNodeSet &            right ) const;

protected:
    //! same as \see partition, but first check \a graph for connected components
    virtual void           scc_partition ( const TGraph &        graph,
                                           TNodeSet &            left,
                                           TNodeSet &            right ) const;

    //! build leaf node for indices in \a graph
    virtual std::unique_ptr< TCluster >  build_leaf  ( const TGraph &        graph,
                                                       const idx_t           idx_ofs,
                                                       TPermutation &        perm,
                                                       std::atomic< int > &  id ) const;
    
    //! adjust n_min based on sparse matrix if default value of 0 was given in constructor
    virtual uint           adjust_n_min  ( any_const_sparse_matrix_t S ) const;

    //! analyze connections between sub graphs \a left and \a right and swap if necessary
    virtual void           check_flow    ( const TGraph &             graph,
                                           TNodeSet &                 left,
                                           TNodeSet &                 right,
                                           any_const_sparse_matrix_t  S ) const;

    DISABLE_COPY_OP( TAlgCTBuilder );
};

//!
//! \ingroup  Cluster_Module
//! \class    TAlgNDCTBuilder
//! \brief    Enhances algebraic clustering by nested dissection.
//!
class TAlgNDCTBuilder : public TAlgCTBuilder
{
private:
    // base clustering algorithm to use for partitioning
    const TAlgCTBuilder *  _alg_ct_builder;

    // mode for handling interface cluster tree depth
    bool                              _sync_interface_depth;

protected:
    //!
    //! \class  TOptClusterSize
    //! \brief  Controls optimal cluster size per tree level.
    //!
    class TOptClusterSize
    {
    private:
        //! optimal size of cluster 
        const size_t   _opt_size;

        //! optimal reduction of indices per level
        const double   _reduction;

    public:
        //! construct size control
        TOptClusterSize ( const size_t  opt_size  = 0,
                          const double  reduction = 0.0 )
                : _opt_size(opt_size), _reduction( reduction )
        {}

        //! return true, if \a n is of optimal size
        bool  is_optimal ( const size_t  n ) const { return n >= _opt_size; }

        //! return size control for next level in tree
        TOptClusterSize  recurse () const
        {
            return TOptClusterSize( size_t( double(_opt_size) * _reduction ),
                                    _reduction );
        }
    };
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct cluster tree builder with \a alg_ct_builder as base algorithm
    TAlgNDCTBuilder ( const TAlgCTBuilder *  alg_ct_builder,
                      const uint             n_min        = CFG::Cluster::nmin,
                      const uint             min_leaf_lvl = 0 );

    //! dtor
    virtual ~TAlgNDCTBuilder ();

    //////////////////////////////////////////////
    //
    // options
    //

    //! set mode for handling interface cluster tree depth, e.g. synchronise with domain clusters
    void  sync_interface_depth ( const bool  b )
    {
        _sync_interface_depth = b;
    }
    
    //////////////////////////////////////////////
    //
    // clustertree creation
    //

    //! divide graph \a graph and build corresponding cluster tree
    //! (instead of bipartition, also build vertex separator son)
    virtual std::unique_ptr< TCluster >
    divide  ( const TGraph &             graph,
              const uint                 lvl,
              TPermutation &             perm,
              const idx_t                idx_ofs,
              const uint                 n_min,
              any_const_sparse_matrix_t  S,
              const uint                 max_lvl,
              std::atomic< int > &       id ) const;
    
    //! partition graph using base algorithm
    virtual void
    partition  ( const TGraph &        graph,
                 TNodeSet &            left,
                 TNodeSet &            right ) const
    {
        _alg_ct_builder->partition( graph, left, right );
    }

protected:
    
    //! divide vertex separator \a vtxsep using connectivity defined by \a graph
    std::unique_ptr< TCluster >
    divide_if ( const TGraph &           graph,
                const TNodeSet &         surrounding,
                const TNodeSet &         vtxsep,
                const uint               lvl,
                const idx_t              idx_ofs,
                TPermutation &           perm,
                const uint               max_lvl,
                const uint               n_min,
                const TOptClusterSize &  csize,
                std::atomic< int > &     id ) const;

    //! build cluster tree for vertex separator \a graph
    std::unique_ptr< TCluster >
    divide_if ( const TGraph &             graph,
                const uint                 lvl,
                TPermutation &             perm,
                const idx_t                idx_ofs,
                const uint                 n_min,
                any_const_sparse_matrix_t  S,
                const uint                 max_lvl,
                const TOptClusterSize &    csize,
                std::atomic< int > &       id ) const;
    
    //! build leaf node for indices in \a nodes
    virtual std::unique_ptr< TCluster >
    build_leaf ( const TGraph &        graph,
                 const TNodeSet &      nodes,
                 const idx_t           idx_ofs,
                 TPermutation &        perm,
                 std::atomic< int > &  id ) const;
    using TAlgCTBuilder::build_leaf;
    
    //! graph partitioning algorithm for vertex separators
    virtual void partition     ( const TGraph &              graph,
                                 const TNodeSet &            surrounding,
                                 const TNodeSet &            vtxsep,
                                 TNodeSet &                  left,
                                 TNodeSet &                  right,
                                 TNodeSet &                  loc_sur ) const;

    //! build vertex separator \a vtxsep between \a left and \a right
    virtual void build_vtx_sep ( const TGraph &              graph,
                                 TNodeSet &                  left,
                                 TNodeSet &                  right,
                                 TNodeSet &                  vtxsep ) const;

    //! do a complete BFS starting at \a start in \a graph but finish
    //! if \a max_nnodes "VERTEX_SEP" nodes have been visited
    //! - visited == true is assumed for all nonlocal nodes
    virtual uint bfs_vtxsep    ( const TGraph &              graph,
                                 TNodeSet &                  start,
                                 TNodeSet &                  last,
                                 std::vector< bool > &       visited,
                                 const std::vector< char > & label,
                                 const uint                  max_nnodes ) const;
    
    //! collect successors to nodes in \a nodes
    //! - search is restricted to nodes with label == \a local
    virtual void bfs_step      ( const TGraph &              graph,
                                 TNodeSet &                  nodes,
                                 TNodeSet &                  succ,
                                 std::vector< bool > &       visited,
                                 const std::vector< char > & label,
                                 const char                  local  ) const;

    //! restrict node set \a nodes to nodes in vertex separator, defined by \a label
    virtual void restrict_vtx  ( TNodeSet &                  nodes,
                                 const std::vector< char > & label ) const;

    //! build strongly connected components of \a graph restricted to \a surrounding
    virtual void build_scc     ( const TGraph &              graph,
                                 const TNodeSet &            surrounding,
                                 std::list< TNodeSet > &     scc ) const;
};

}// namespace Hpro

#endif  // __HPRO_TALGCTBUILDER_HH
