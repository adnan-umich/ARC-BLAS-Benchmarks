#ifndef __HPRO_TBSPPARTSTRAT_HH
#define __HPRO_TBSPPARTSTRAT_HH
//
// Project     : HLIBpro
// File        : TBSPPartStrat.hh
// Description : partitioning strategies for geometrical clustering
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/cluster/types.hh"
#include "hpro/cluster/TNodeSet.hh"
#include "hpro/cluster/TGeomCluster.hh"
#include "hpro/cluster/TCoordinate.hh"

#include "hpro/matrix/TSparseMatrix.hh"

namespace Hpro
{

//!
//! \ingroup  Cluster_Module
//! \class    TBSPPartStrat
//! \brief    Base class for partitioning strategies for geometrical BSP clustering.
//!
class TBSPPartStrat
{
public:
    // dtor
    virtual ~TBSPPartStrat () {}
    
    //!
    //! partition \a dofs into sub sets and store result in \a partition.
    //! \a bbox holds bounding box of \a dofs cluster. Bounding boxes of
    //! sons will be stored in \a son_bbox.
    //!
    virtual void partition ( const TCoordinate *     coord,
                             const TNodeSet  &       dofs,
                             TNodeSet &              left,
                             TNodeSet &              right,
                             const TBBox &           bbox,
                             std::vector< TBBox > &  son_bbox,
                             const uint              depth ) const = 0;
};

//!
//! \ingroup  Cluster_Module
//! \class    TGeomBSPPartStrat
//! \brief    Partition according to geometrical volume of index sets.
//!
class TGeomBSPPartStrat : public TBSPPartStrat
{
private:
    // mode for choose splitting axis
    const split_axis_mode_t  _split_axis_mode;
    
public:
    //!
    //! ctor 
    //!
    TGeomBSPPartStrat ( const split_axis_mode_t  split_axis_mode = CFG::Cluster::split_mode );
    
    //!
    //! partition \a dofs into sub sets
    //!
    virtual void partition ( const TCoordinate *     coord,
                             const TNodeSet &        dofs,
                             TNodeSet &              left,
                             TNodeSet &              right,
                             const TBBox &           bbox,
                             std::vector< TBBox > &  son_bbox,
                             const uint              depth ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TCardBSPPartStrat
//! \brief    Partition according to cardinality of index sets.
//!
class TCardBSPPartStrat : public TBSPPartStrat
{
private:
    // mode for choose splitting axis
    const split_axis_mode_t   _split_axis_mode;
    
public:
    //!
    //! ctor
    //!
    TCardBSPPartStrat ( const split_axis_mode_t  split_axis_mode = CFG::Cluster::split_mode );
    
    //!
    //! partition \a dofs into sub sets
    //!
    virtual void partition ( const TCoordinate *     coord,
                             const TNodeSet &        dofs,
                             TNodeSet &              left,
                             TNodeSet &              right,
                             const TBBox &           bbox,
                             std::vector< TBBox > &  son_bbox,
                             const uint              depth ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TPCABSPPartStrat
//! \brief    Partition according to principle component analysis.
//!
class TPCABSPPartStrat : public TBSPPartStrat
{
private:
    // use cardinality balanced clustering
    const bool  _use_card;
    
public:
    //!
    //! ctor
    //!
    TPCABSPPartStrat ( const bool  use_card = false );
    
    //!
    //! partition \a dofs into sub sets
    //!
    virtual void partition ( const TCoordinate *     coord,
                             const TNodeSet  &       dofs,
                             TNodeSet &              left,
                             TNodeSet &              right,
                             const TBBox &           bbox,
                             std::vector< TBBox > &  son_bbox,
                             const uint              depth ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TNDBSPPartStrat
//! \brief    Special partition strategy to optimized nested dissection clustering
//!
template < typename value_t >
class TNDBSPPartStrat : public TBSPPartStrat
{
private:
    // sparse matrix containing connectivity information
    const TSparseMatrix< value_t > *  _sparse_mat;
    
    // do/don't use weighted edge cut
    const bool                        _use_edgecut_weights;
    
public:
    //!
    //! ctor
    //!
    TNDBSPPartStrat ( const TSparseMatrix< value_t > *  S,
                      const edgecut_weights_mode_t      edgecut_weights_mode = edgecut_weights_off );
    
    //!
    //! partition \a dofs into sub sets
    //!
    virtual void partition ( const TCoordinate *     coord,
                             const TNodeSet  &       dofs,
                             TNodeSet &              left,
                             TNodeSet &              right,
                             const TBBox &           bbox,
                             std::vector< TBBox > &  son_bbox,
                             const uint              depth ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TAutoBSPPartStrat
//! \brief    Automatic choice of best partitioning strategy.
//!
class TAutoBSPPartStrat : public TBSPPartStrat
{
private:
    // auto is based on geometrical and cardinality partitioning
    TGeomBSPPartStrat  _geom;
    TCardBSPPartStrat  _card;
    
public:
    //!
    //! ctor
    //!
    TAutoBSPPartStrat ( const split_axis_mode_t  split_axis_mode = CFG::Cluster::split_mode );
    
    //!
    //! partition \a dofs into sub sets
    //!
    virtual void partition ( const TCoordinate *     coord,
                             const TNodeSet &        dofs,
                             TNodeSet &              left,
                             TNodeSet &              right,
                             const TBBox &           bbox,
                             std::vector< TBBox > &  son_bbox,
                             const uint              depth ) const;
};

}// namespace

#endif // __HPRO_TBSPPARTSTRAT_HH
