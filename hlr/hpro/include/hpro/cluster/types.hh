#ifndef __HPRO_CLUSTER_TYPES_HH
#define __HPRO_CLUSTER_TYPES_HH
//
// Project     : HLIBpro
// File        : types.hh
// Description : types used for clustering functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

namespace Hpro
{

//!
//! \ingroup  Cluster_Module
//! \enum     split_axis_t
//! \brief    enumeration for different modes of choosing splitting axis
//!
enum split_axis_mode_t
{
    regular_split_axis,   //!< alternate splitting axis regularly (…-x-y-z-x-y-z-…)
    adaptive_split_axis   //!< choose splitting axis based on box dimensions
};

//!
//! \ingroup  Cluster_Module
//! \enum     adjust_depth_mode_t
//! \brief    turns on/off depth adjustment of sub trees
//!
enum adjust_depth_mode_t : bool
{
    adjust_depth_off = false,   //!< do not adjust depth
    adjust_depth_on  = true     //!< adjust depth of sub trees
};
    
//!
//! \ingroup  Cluster_Module
//! \enum     edgecut_weights_mode_t
//! \brief    turn on/off usage of weights in edge cut
//!
enum edgecut_weights_mode_t : bool
{
    edgecut_weights_off = false,    //!< do not use weights in edge cut
    edgecut_weights_on  = true      //!< use weights in edge cut, e.g. based on matrix coefficients
};

//!
//! \ingroup  Cluster_Module
//! \enum     edge_weights_mode_t
//! \brief    turn on/off usage of edge weights
//!
enum edge_weights_mode_t
{
    edge_weights_off,     //!< do not use edge weights, e.g. constant weights
    edge_weights_on,      //!< use edge weights, e.g. based on matrix coefficients
    edge_weights_on_sym   //!< use symmetric edge weights
};

//!
//! \ingroup  Cluster_Module
//! \enum     cluster_level_mode_t
//! \brief    turn on/off enforcement of same cluster level in block cluster
//!
enum cluster_level_mode_t : bool
{
    cluster_level_any  = false,   //!< allow arbitrary cluster levels in block cluster
    cluster_level_same = true     //!< enforce same cluster level in block cluster
};

//!
//! \ingroup  Cluster_Module
//! \enum     diam_mode_t
//! \brief    choose between minimum and maximum of cluster diameters in std. adm.
//!
enum diam_mode_t
{
    use_min_diam,     //!< use minimum of cluster diameters
    use_max_diam      //!< use maximum of cluster diameters
};

}// namespace Hpro

#endif // __HPRO_CLUSTER_TYPES_HH
