#ifndef __HPRO_TBLOCKCLUSTERTREE_HH
#define __HPRO_TBLOCKCLUSTERTREE_HH
//
// Project     : HLIBpro
// File        : TBlockClusterTree.hh
// Description : class for a cluster tree
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/cluster/TClusterTree.hh"
#include "hpro/cluster/TBlockCluster.hh"

namespace Hpro
{

//!
//! \ingroup  Cluster_Module
//! \class    TBlockClusterTree
//! \brief    Represents a block cluster tree.
//!
class TBlockClusterTree
{
protected:
    //! root of the block cluster tree
    TBlockCluster *       _root;
    
    //! row cluster tree
    const TClusterTree *  _row_ct;

    //! column cluster tree
    const TClusterTree *  _col_ct;

public:
    ////////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct block cluster tree with \a aroot as root and
    //! \a arow_ct and \a acol_ct as corresponding cluster trees
    TBlockClusterTree ( TBlockCluster *       aroot,
                        const TClusterTree *  arow_ct,
                        const TClusterTree *  acol_ct )
            : _root(aroot)
            , _row_ct(arow_ct)
            , _col_ct(acol_ct)
    {
        if ( _root == NULL )
            HERROR( ERR_ARG, "(TBlockClusterTree)", "root is NULL" );
    }

    //! delete block cluster tree (but not the cluster trees)
    virtual ~TBlockClusterTree ()
    {
        delete _root;
//         delete _row_ct;
//         delete _col_ct;
    }

    ////////////////////////////////////////////////////////
    //
    // access local data
    //

    //! return root of block cluster tree
    const TBlockCluster *  root    () const { return _root; }
    TBlockCluster *        root    ()       { return _root; }

    //! return row cluster tree
    const TClusterTree *   row_ct  () const { return _row_ct; }

    //! return column cluster tree
    const TClusterTree *   col_ct  () const { return _col_ct; }

    ////////////////////////////////////////////////////////
    //
    // access tree data
    //

    //! return number of nodes in tree
    uint  nnodes        () const { return _root->nnodes(); }
    
    //! return depth of tree
    uint  depth         () const { return _root->depth(); }
    
    //! compute sparsity constant of tree
    uint  compute_c_sp  ()                    const { return _root->compute_c_sp(); }

    //! compute sharing constant of tree with distribution onto \a nprocs processors
    uint  compute_c_sh  ( const uint nprocs ) const { return _root->compute_c_sh( nprocs ); }
    
    ////////////////////////////////////////////////////////
    //
    // misc. methods
    //

    //! collect leaves or nodes with depth \a adepth in tree
    virtual void    collect_leaves  ( std::list< TBlockCluster * > &  leaves,
                                      const int                       adepth = -1,
                                      const int                       level = 0 ) const
    {
        _root->collect_leaves( leaves, adepth, level );
    }

    //! return size in bytes used by this object
    virtual size_t  byte_size       () const
    {
        return _root->byte_size() + sizeof(TBlockCluster*) + 2 * sizeof(TClusterTree*);
    }
};

}// namespace Hpro

#endif  // __HPRO_TBLOCKCLUSTERTREE_HH
