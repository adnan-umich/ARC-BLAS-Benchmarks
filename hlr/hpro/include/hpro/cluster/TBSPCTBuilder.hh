#ifndef __HPRO_TBSPCTBUILDER_HH
#define __HPRO_TBSPCTBUILDER_HH
//
// Project     : HLIBpro
// File        : TBSPCTBuilder.hh
// Description : build clustertrees for coordinate based indexsets via BSP
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>
#include <deque>

#include "hpro/base/config.hh"

#include "hpro/cluster/types.hh"
#include "hpro/cluster/TNodeSet.hh"
#include "hpro/cluster/TGeomCluster.hh"
#include "hpro/cluster/TClusterTree.hh"
#include "hpro/cluster/TBSPPartStrat.hh"
#include "hpro/cluster/TCoordinate.hh"
#include "hpro/cluster/TPermutation.hh"

#include "hpro/matrix/TSparseMatrix.hh"

namespace Hpro
{

//!
//! \ingroup  Cluster_Module
//! \class    TGeomCTBuilder
//! \brief    Base class for all cluster tree constructors based on geometry data.
//!
class TGeomCTBuilder
{
protected:
    //!
    //! \struct  data_t
    //! \brief   Datatype for internal argument transfer
    //!
    struct data_t
    {
        const TCoordinate *  coord;        // coordinate information of indices
        TPermutation *       perm;         // permutation to built
        const uint           nmin;         // minimale leaf size
        const uint           min_leaf_lvl; // minimal level where leaves might appear
        const uint           max_lvl;      // maximal level to reach in clustering
        std::atomic< int >   id;           // unique cluster id
    };

    //!
    //! \class  TOptClusterSize
    //! \brief  Controls optimal cluster size per tree level.
    //!
    class TOptClusterSize
    {
    private:
        //! optimal size of cluster 
        size_t   _opt_size;

        //! optimal reduction of indices per level
        double   _reduction;

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

protected:
    //! minimal size of a cluster, i.e. not smaller than this
    uint                   _n_min;
 
    //! minimal level on which leaves may occur
    uint                   _min_leaf_lvl;

    //! maximal level of tree (0: automatic choice)
    uint                   _max_lvl;

    //! flag for adjusting bounding volumes of nodes
    bool                   _adjust_bvol;

    //! flag for sorting sub clusters w.r.t. size
    bool                   _sort_wrt_size;
    

public:
    //! construct cluster tree builder
    TGeomCTBuilder ( const uint  n_min        = CFG::Cluster::nmin,
                     const uint  min_leaf_lvl = 0 );

    // dtor
    virtual ~TGeomCTBuilder () {}
    
    //!
    //! build cluster tree out of given coordinate set
    //! \param   coord   : geometry information for each index
    //! \param   idx_ofs : start renumbering indices from \a idx_ofs
    //!
    virtual
    std::unique_ptr< TClusterTree >
    build  ( const TCoordinate *      coord,
             const idx_t              idx_ofs = 0 ) const;

    //! recursively build cluster tree for indices in \a dofs
    virtual
    std::unique_ptr< TGeomCluster >
    divide ( const TNodeSet &         dofs,
             const uint               lvl,
             const TBoundingVolume &  bvol,
             const TOptClusterSize &  csize,
             const idx_t              index_ofs,
             data_t &                 data ) const = 0;

    //
    // give access to local parameters
    //

    uint  n_min          () const { return _n_min;         }
    uint  min_leaf_lvl   () const { return _min_leaf_lvl;  }
    uint  max_lvl        () const { return _max_lvl;       }
    bool  adjust_bvol      () const { return _adjust_bvol;     }
    bool  sort_wrt_size  () const { return _sort_wrt_size; }

    //! set maximal leaf level
    void  set_max_lvl    ( const uint  l ) { _max_lvl = l; }
    
    //! set flag for adjusting bounding volume
    void  adjust_bvol      ( const bool  b ) { _adjust_bvol = b; }

    //! set flag for sorting son cluster wrt. size
    void  sort_wrt_size  ( const bool  b ) { _sort_wrt_size = b; }
    
protected:
    //! create a leaf in a clustertree containing indices in \a dofs
    virtual
    std::unique_ptr< TGeomCluster >
    build_leaf ( const TNodeSet &         dofs,
                 const uint               lvl,
                 const idx_t              index_ofs,
                 const TBoundingVolume &  bvol,
                 data_t &                 data ) const;

    //! compute bounding box of index set defined by \a dofs
    virtual
    TBBox
    compute_bbox    ( const TNodeSet &  dofs,
                      const data_t &    data ) const;

    //! compute bounding sphere of index set defined by \a dofs
    virtual
    TBSphere
    compute_bsphere ( const TNodeSet &  dofs,
                      const data_t &    data ) const;

    //! compute bounding volume of index set defined by \a dofs
    virtual
    TBoundingVolume
    compute_bvol    ( const TNodeSet &  dofs,
                      const data_t &    data ) const;

    //! update bounding volume of index set defined by \a dofs
    virtual
    void
    update_bvol     ( const TNodeSet &   dofs,
                      TBoundingVolume &  bvol,
                      const data_t &     data ) const;

    //! check and update bvol in case of degenerate axis, e.g. very small length
    virtual
    void
    check_bvol      ( TBoundingVolume &  bvol ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TBSPCTBuilder
//! \brief    Base class for all cluster tree constructors based on BSP.
//!
class TBSPCTBuilder : public TGeomCTBuilder
{
protected:
    //! type of partitioning strategy
    const TBSPPartStrat *  _part_strat;
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct BSP cluster tree builder with partitioning strategy \a part_strat
    TBSPCTBuilder ( const TBSPPartStrat *  part_strat,
                    const uint             n_min        = CFG::Cluster::nmin,
                    const uint             min_leaf_lvl = 0 );

    //! dtor
    virtual ~TBSPCTBuilder ();

    //////////////////////////////////////////////
    //
    // clustertree creation
    //

    //! recursively build cluster tree for indices in \a dofs
    virtual
    std::unique_ptr< TGeomCluster >
    divide ( const TNodeSet &         dofs,
             const uint               lvl,
             const TBoundingVolume &  bvol,
             const TOptClusterSize &  csize,
             const idx_t              index_ofs,
             data_t &                 data ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TBSPNDCTBuilder
//! \brief    Combines binary space partitioning with nested dissection
//!           based on connectivity defined by a sparse matrix.
//!
class TBSPNDCTBuilder : public TBSPCTBuilder
{
private:
    //! sparse matrix for connectivity between indices
    any_const_sparse_matrix_t  _sparse_mat;

    // mode for handling interface cluster tree depth
    bool                       _sync_interface_depth;

public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct cluster tree with partition strategy defined by \a part_strat
    //! and connectivity defined by \a S
    TBSPNDCTBuilder ( any_const_sparse_matrix_t  S,
                      const TBSPPartStrat *      part_strat,
                      const uint                 n_min        = CFG::Cluster::nmin,
                      const uint                 min_leaf_lvl = 0 );

    //! construct cluster tree with partition strategy defined by \a part_strat
    //! (must use "build( coord, S )")
    TBSPNDCTBuilder ( const TBSPPartStrat *  part_strat,
                      const uint             n_min        = CFG::Cluster::nmin,
                      const uint             min_leaf_lvl = 0 );

    virtual ~TBSPNDCTBuilder ();

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

    //! build nested dissection cluster tree out of coordinate set \a coord and
    //! additional connectivity defined by local sparse matrix
    virtual
    std::unique_ptr< TClusterTree >
    build     ( const TCoordinate *      coord,
                const idx_t              idx_ofs = 0 ) const;

    //! build nested dissection cluster tree out of coordinate set \a coord and
    //! additional connectivity defined by sparse matrix \a S
    virtual
    std::unique_ptr< TClusterTree >
    build     ( const TCoordinate *        coord,
                any_const_sparse_matrix_t  S,
                const idx_t                idx_ofs = 0 ) const;
    
    //! recursively build cluster tree
    virtual
    std::unique_ptr< TGeomCluster >
    divide    ( const TNodeSet &         dofs,
                const uint               lvl,
                const TBoundingVolume &  bvol,
                const TOptClusterSize &  csize,
                const idx_t              index_ofs,
                data_t &                 data ) const;

    //! recursively build cluster tree for interfaces clusters
    virtual
    std::unique_ptr< TGeomCluster >
    divide_if ( const TNodeSet &         dofs,
                const uint               lvl,
                const uint               max_lvl,
                const TBoundingVolume &  bvol,
                const TOptClusterSize &  csize,
                const idx_t              index_ofs,
                data_t &                 data ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TMBLRBSPCTBuilder
//! \brief    implement MBLR clustering
//!
class TMBLRCTBuilder : public TGeomCTBuilder
{
protected:
    //! type of partitioning strategy
    const TBSPPartStrat *  _part_strat;

    //! number of levels in MBLR clustering
    const size_t           _nlevel;
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct BSP cluster tree builder with partitioning strategy \a part_strat
    TMBLRCTBuilder ( const size_t           nlevel,
                     const TBSPPartStrat *  part_strat,
                     const uint             n_min = CFG::Cluster::nmin );

    //! dtor
    virtual ~TMBLRCTBuilder ();

    //////////////////////////////////////////////
    //
    // clustertree creation
    //

    //! build cluster tree out of given coordinate set
    virtual
    std::unique_ptr< TClusterTree >
    build ( const TCoordinate *  coord,
            const idx_t          index_ofs = 0 ) const;
    
    //! recursively build cluster tree for indices in \a dofs
    virtual
    std::unique_ptr< TGeomCluster >
    divide ( const TNodeSet &         dofs,
             const uint               lvl,
             const TBoundingVolume &  bvol,
             const TOptClusterSize &  csize,
             const idx_t              index_ofs,
             data_t &                 data ) const;

    //! recursively build cluster tree for leaf-level partitioning in \a leaves
    virtual
    std::unique_ptr< TGeomCluster >
    divide ( const std::list< TNodeSet > &         leaves,
             const std::list< TBoundingVolume > &  leaf_bvol,
             const TBoundingVolume &               bvol,
             const uint                            lvl,
             const idx_t                           index_ofs,
             data_t &                              data ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TSFCCTBuilder
//! \brief    construct cluster tree using space filling curves
//!
class TSFCCTBuilder : public TGeomCTBuilder
{
public:
    //!
    //! various options
    //!
    enum cluster_type_t
    {
        binary,     // standard, binary partitioning
        blr         // single level partitioning (no hierarchy)
    };

    enum split_type_t
    {
        volume,     // split based on volume
        diameter,   // split based on diameter
        cardinality // split based on cardinality
    };

protected:
    //! partitioning type
    cluster_type_t  _cl_type;

    //! splitting type
    split_type_t    _split_type;
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor
    TSFCCTBuilder ( const cluster_type_t  cl_type      = binary,
                    const split_type_t    split_type   = cardinality,
                    const uint            n_min        = CFG::Cluster::nmin,
                    const uint            min_leaf_lvl = 0 );

    //! dtor
    virtual ~TSFCCTBuilder ();

    //////////////////////////////////////////////
    //
    // clustertree creation
    //

    //! build cluster tree out of given coordinate set \a coord
    virtual std::unique_ptr< TClusterTree >  build  ( const TCoordinate *      coord,
                                                      const idx_t              idx_ofs = 0 ) const;
    
    //! recursively build cluster tree for indices in \a dofs
    virtual std::unique_ptr< TGeomCluster >  divide ( const TNodeSet &         dofs,
                                                      const uint               lvl,
                                                      const TBoundingVolume &  bvol,
                                                      const TOptClusterSize &  csize,
                                                      const idx_t              index_ofs,
                                                      data_t &                 data ) const;

protected:
    //! split into <nsub> sub sets
    virtual
    std::pair< std::vector< TNodeSet >,
               std::vector< TBoundingVolume > >
    split ( const uint        nsub,
            const TNodeSet &  dofs,
            data_t &          data ) const;
};

}// namespace Hpro

#endif  // __HPRO_TBSPCTBUILDER_HH
