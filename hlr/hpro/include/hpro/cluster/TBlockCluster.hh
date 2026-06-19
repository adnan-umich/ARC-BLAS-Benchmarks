#ifndef __HPRO_TBLOCKCLUSTER_HH
#define __HPRO_TBLOCKCLUSTER_HH
//
// Project     : HLIBpro
// File        : TBlockCluster.hh
// Description : class for a block cluster tree
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>
#include <list>

#include "hpro/cluster/TCluster.hh"

#include "hpro/parallel/TProcSet.hh"

namespace Hpro
{

//!
//! \ingroup  Cluster_Module
//! \class    TBlockCluster
//! \brief    Representing a node in a block cluster tree as product of two clusters.
//!
class TBlockCluster
{
protected:
    //!@cond
    
    //! globally unique id
    int                             _id;

    //! parent cluster
    TBlockCluster *                 _parent;
    
    //! row/column clusters
    TCluster *                      _rowcl;
    TCluster *                      _colcl;

    //! block layout
    size_t                          _nrows;
    size_t                          _ncols;
    
    //! son-clusters (nrows × ncols)
    std::vector< TBlockCluster * >  _sons;

    //! flag for admissiblity
    bool                            _adm;
    
    //! local processor set
    TProcSet                        _procs;
    
    //!@endcond
    
public:
    ////////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct empty, leaf block cluster with \a parent as parent node
    TBlockCluster ( TBlockCluster *  parent );

    //! construct leaf block cluster defined by \a rowcl × \a colcl with \a parent as parent node
    TBlockCluster ( TBlockCluster *  parent,
                    TCluster *       rowcl,
                    TCluster *       colcl );

    //! dtor
    virtual ~TBlockCluster ();
    
    ////////////////////////////////////////////////////////
    //
    // access local variables and flags
    //

    //! return ID
    int                    id           () const { return _id; }

    //! set ID
    void                   set_id       ( const int  aid ) { _id = aid; }

    //! return parent cluster
    TBlockCluster *        parent       ()       { return _parent; }
    const TBlockCluster *  parent       () const { return _parent; }

    //! set parent cluster to \a bct
    void                   set_parent   ( TBlockCluster * bct ) { _parent = bct; }

    //! return processor set
    const TProcSet &       procs        () const { return _procs; }

    //! set local processor set and if \a recursive is true also son sets to \a ps
    void                   set_procs    ( const TProcSet &  ps,
                                          const bool        recursive = false );

    //! return row cluster
    TCluster *             rowcl        ()       { return _rowcl; }
    const TCluster *       rowcl        () const { return _rowcl; }
    TIndexSet              rowis        () const { if ( _rowcl != nullptr ) { return *_rowcl; } else { return TIndexSet(); } }

    //! return column cluster
    TCluster *             colcl        ()       { return _colcl; }
    const TCluster *       colcl        () const { return _colcl; }
    TIndexSet              colis        () const { if ( _colcl != nullptr ) { return *_colcl; } else { return TIndexSet(); } }

    //! set row cluster
    void                   set_rowcl    ( TCluster *  cl );

    //! set column cluster
    void                   set_colcl    ( TCluster *  cl );

    //! set row and column cluster
    void                   set_clusters ( TCluster *  row_cl,
                                          TCluster *  col_cl );

    //! return block index set of block cluster
    TBlockIndexSet         is           () const
    {
        if (( _rowcl != nullptr ) && ( _colcl != nullptr ))
            return TBlockIndexSet( * _rowcl, * _colcl );
        else
            return TBlockIndexSet();
    }
    
    //! return true if block cluster is admissible
    bool                   is_adm       () const         { return _adm; }

    //! set admissibility of block cluster to \a b
    void                   set_adm      ( const bool b ) { _adm = b; }

    ////////////////////////////////////////////////
    //
    // tree interface
    //

    //! return number of sons
    virtual size_t     nsons         () const { return _sons.size(); }

    // //! set number of sons
    // virtual void    set_nsons     ( const size_t  n );

    // //! adjust number of sons to number of non-nullptr sons in local list
    // virtual void    adjust_nsons  ();

    //
    // access sons in a linear fashion
    //

    //! return \a i'th son
    virtual TBlockCluster *        son         ( const size_t  i )       { return _sons[i]; }
    virtual const TBlockCluster *  son         ( const size_t  i ) const { return _sons[i]; }

    //! set \a i'th son to \a son
    virtual void                   set_son     ( const size_t     i,
                                                 TBlockCluster *  son,
                                                 const bool       del_son = true );

    //! add \a son to set of sons at first free slot
    virtual void                   add_son     ( TBlockCluster *  son );

    //
    // access sons block-wise
    //

    //! change block layout
    virtual void                   set_layout  ( const size_t  nrows,
                                                 const size_t  ncols );

    //! return number of rows
    virtual size_t                 nrows       () const { return _nrows; }
    
    //! return number of columns
    virtual size_t                 ncols       () const { return _ncols; }
    
    //! return son at position (\a i, \a j)
    virtual TBlockCluster *        son         ( const size_t  i,
                                                 const size_t  j );
    virtual const TBlockCluster *  son         ( const size_t  i,
                                                 const size_t  j ) const;

    //! set son at position (\a i, \a j) to \a son
    virtual void                   set_son     ( const size_t     i,
                                                 const size_t     j,
                                                 TBlockCluster *  son,
                                                 const bool       del_son = true );

    //
    // access sons wrt. clusters
    //

    //! return son corresponding to block cluster \a rowcl × \a colcl
    virtual TBlockCluster *       son_cl       ( const TCluster *  rowcl,
                                                 const TCluster *  colcl );
    virtual const TBlockCluster * son_cl       ( const TCluster *  rowcl,
                                                 const TCluster *  colcl ) const;

    //
    // tree-properties
    //
    
    //! return true of node is leaf
    virtual bool  is_leaf () const
    {
        if ( nsons() == 0 )
            return true;
        
        for ( size_t  i = 0; i < nsons(); i++ )
            if ( _sons[i] != nullptr )
                return false;
        
        return true;
    }

    //! make node a leaf
    virtual void    make_leaf ();

    //! return number of nodes in tree
    virtual size_t  nnodes () const;

    //! return depth of tree
    virtual size_t  depth () const;

    ////////////////////////////////////////////////
    //
    // procset managment
    //
    
    //! adjust processor set such that local set is union of son sets
    void assign_procs ();
    
    ////////////////////////////////////////////////////////
    //
    // misc. methods
    //

    //! return true if \a i'th son is present
    bool has_son ( const size_t  i ) const { return (i < nsons() ? _sons[i] != nullptr : false); }

    //! return true if given cluster is a subcluster of this
    bool is_sub_cluster ( const TBlockCluster * c ) const;
    
    //! return object of same type
    virtual TBlockCluster * create () const { return new TBlockCluster( nullptr ); }

    //! collect leaves or nodes with depth \a depth in tree
    void collect_leaves ( std::list< TBlockCluster * > &  leaves,
                          const int                       depth = -1,
                          const int                       level = 0 ) const;

    //! compute sparsity constant of tree
    //! - if \a leaves is true, only leaf blocks are counted
    //! - if \a adm is true, only admissible blocks (leaves) are counted
    uint compute_c_sp ( const bool  leaves = false,
                        const bool  adm    = false ) const;

    //! compute sharing constant of tree for index set partition with \a nprocs processors
    uint compute_c_sh ( const uint  nprocs ) const;
    
    //! return copy of node/subtree
    virtual TBlockCluster * copy  () const;

    //! return size in bytes used by this object
    virtual size_t byte_size () const;

    //! stream output
    void print ( const uint ofs = 0 ) const;

    //! return string representation
    std::string  to_string () const
    {
        if ( rowcl() != nullptr )
            if ( colcl() != nullptr ) return rowcl()->to_string() + " × " + colcl()->to_string();
            else                      return rowcl()->to_string() + " × ∅";
        else
            if ( colcl() != nullptr ) return "∅ × " + colcl()->to_string();
            else                      return "∅ × ∅";
    }
};

//
//! compute sparsity constant of tree
//
inline
uint
compute_c_sp ( const TBlockCluster &  bc,
               const bool             leaves = false,
               const bool             adm    = false  )
{
    return bc.compute_c_sp( leaves, adm );
}

//!
//! \ingroup  Cluster_Module
//! \brief    flatten hierarchy of block cluster tree to first level with leaves
//! \detail   The first levels in the block cluster tree with no leaves are eliminated
//!           such that the tree root will directly have sons on this level.
//!
void flatten_leaf ( TBlockCluster * root ); 

/////////////////////////////////////////////////////////////////////////////
//
// apply functor for block cluster trees
//

template < typename Functor >
void apply ( const TBlockCluster *  bc,
             const Functor &        f )
{
    if ( bc == nullptr )
        return;

    f( bc );

    if ( ! bc->is_leaf() )
    {
        for ( size_t  i = 0; i < bc->nsons(); ++i )
            apply( bc->son( i ), f );
    }// if
}

}// namespace Hpro

#endif  // __HPRO_TBLOCKCLUSTER_HH
