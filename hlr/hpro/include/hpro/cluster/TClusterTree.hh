#ifndef __HPRO_TCLUSTERTREE_HH
#define __HPRO_TCLUSTERTREE_HH
//
// Project     : HLIBpro
// File        : TClusterTree.hh
// Description : class for a cluster tree
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/cluster/TCluster.hh"
#include "hpro/cluster/TPermutation.hh"

namespace Hpro
{

//!
//! \ingroup  Cluster_Module
//! \class    TClusterTree
//! \brief    Represents a cluster tree with permutation of index sets
//!
class TClusterTree
{
private:
    //! root of the cluster tree
    TCluster *           _root;

    //! permutations from external to internal numbering
    const TPermutation * _perm_e2i;

    //! permutations from internal to external numbering
    const TPermutation * _perm_i2e;

public:
    ////////////////////////////////////////////////////////
    //
    // constructor
    //

    //! construct cluster tree with root as \a rootcl and permutations
    //! \a perme2i and \a permi2e
    TClusterTree ( TCluster *            rootcl,
                   const TPermutation *  perme2i,
                   const TPermutation *  permi2e )
            : _root(rootcl)
            , _perm_e2i(perme2i)
            , _perm_i2e(permi2e)
    {
        if ( _root == NULL )
            HERROR( ERR_ARG, "(TClusterTree)", "root is NULL" );
    }

    //! deconstruct cluster tree and permutation objects
    virtual ~TClusterTree ()
    {
        delete _root;
        delete _perm_e2i;
        delete _perm_i2e;
    }

    ////////////////////////////////////////////////////////
    //
    // access local data
    //

    //! return root of cluster tree
    TCluster *           root     ()       { return _root; }

    //! return root of cluster tree
    const TCluster *     root     () const { return _root; }

    //! return external to internal permutation
    const TPermutation * perm_e2i () const { return _perm_e2i; }

    //! return internal to external permutation
    const TPermutation * perm_i2e () const { return _perm_i2e; }

    //! release ownership of data (use with care!!!)
    void                 release_data ()
    {
        _root     = nullptr;
        _perm_e2i = nullptr;
        _perm_i2e = nullptr;
    }
    
    ////////////////////////////////////////////////
    //
    // manage tree
    //

    //! return no of nodes
    virtual uint  nnodes () const
    {
        return _root->nnodes();
    }

    //! depth of tree
    virtual uint  depth () const
    {
        return _root->depth();
    }

    ////////////////////////////////////////////////////////
    //
    // misc. methods
    //

    //! collect leaves or nodes with depth \a adepth in \a list
    virtual void collect_leaves ( std::list< TCluster * > &  leaves,
                                  const int                  adepth = -1,
                                  const int                  level = 0 ) const
    {
        _root->collect_leaves( leaves, adepth, level );
    }

    //! flatten hierarchy of cluster tree
    virtual void flatten ()
    {
        Hpro::flatten( _root );
    }    
    
    //! return size in bytes used by this object
    virtual size_t byte_size () const
    {
        size_t  size = _root->byte_size();

        if ( _perm_e2i != NULL ) size += _perm_e2i->byte_size();
        if ( _perm_i2e != NULL ) size += _perm_i2e->byte_size();

        return size;
    }
};

}// namespace Hpro

#endif  // __HPRO_TCLUSTERTREE_HH
