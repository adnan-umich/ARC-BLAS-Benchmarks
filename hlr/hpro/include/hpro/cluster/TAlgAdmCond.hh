#ifndef __HPRO_TALGADMCOND_HH
#define __HPRO_TALGADMCOND_HH
//
// Project     : HLIBpro
// File        : TAlgAdmCond.hh
// Description : algebraic admissibility condition for sparse matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/cluster/TPermutation.hh"
#include "hpro/cluster/TNodeSet.hh"
#include "hpro/cluster/TAdmCondition.hh"

#include "hpro/matrix/TSparseMatrix.hh"

namespace Hpro
{

//!
//! \ingroup  Cluster_Module
//! \class    TAlgAdmCond
//! \brief    base class for algebraic admissibility conditions
//!
class TAlgAdmCond : public TAdmCondition
{
protected:
    //! @cond
    
    // sparse matrix defining the matrix graph
    any_const_sparse_matrix_t  _mat;

    // mapping of index-names from external (in sparse matrix)
    // to internal numbering (in cluster tree)
    const TPermutation *       _row_perm_i2e, * _col_perm_i2e;
    TPermutation *             _row_perm_e2i, * _col_perm_e2i;
    
    //! @endcond
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor with graph defined by \a S and mapping of internal to external
    //! indices defined by \a perm_i2e (row and column mappings identical)
    TAlgAdmCond ( any_const_sparse_matrix_t  S,
                  const TPermutation *              perm_i2e = nullptr );

    //! ctor with graph defined by \a S and mapping of internal to external
    //! indices defined by \a row_perm_i2e and \a col_perm_i2e
    TAlgAdmCond ( any_const_sparse_matrix_t  S,
                  const TPermutation *              row_perm_i2e,
                  const TPermutation *              col_perm_i2e );

    //! dtor
    virtual ~TAlgAdmCond ();

protected:

    DISABLE_COPY_OP( TAlgAdmCond );
};

//!
//! \ingroup  Cluster_Module
//! \class    TStdAlgAdmCond
//! \brief    Standard admissibility condition based on matrix graph criteria.
//!
class TStdAlgAdmCond : public TAlgAdmCond
{
protected:
    //! @cond
    
    // admissibility parameter
    const double                 _eta;
    
    // mark visited nodes (not mt-safe !!!)
    mutable std::vector< bool >  _visited;
    
    //! @endcond
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor
    TStdAlgAdmCond ( const double               eta,
                     any_const_sparse_matrix_t  S,
                     const TPermutation *       perm_i2e = nullptr );

    //! ctor
    TStdAlgAdmCond ( const double               eta,
                     any_const_sparse_matrix_t  S,
                     const TPermutation *       row_perm_i2e,
                     const TPermutation *       col_perm_i2e );

    //! dtor
    virtual ~TStdAlgAdmCond () {}

    ///////////////////////////////////////////
    //
    // check block-cluster if admissible
    //

    //! return true if \a cl is admissible
    virtual bool  is_adm    ( const TBlockCluster * cl ) const;

protected:
    //! determine diameter of cluster \a cl
    virtual uint  diameter  ( const TCluster *      cl,
                              const TPermutation *  perm_i2e,
                              const TPermutation *  perm_e2i ) const;

    //! Perform a BFS from set \a start in matrix and store last visited nodes
    //! in \a last. Stop BFS if all nodes in \a tau have been visited. Return the depth
    //! of the BFS iteration.
    virtual uint  bfs       ( TNodeSet &            start,
                              TNodeSet &            last,
                              const TCluster *      tau,
                              const TPermutation *  perm_i2e,
                              const TPermutation *  perm_e2i ) const;

    //! return true, if distance between \a tau and \a sigma is bigger than \a min_dist
    virtual bool  cmp_dist  ( const TCluster *      tau,
                              const TCluster *      sigma,
                              const uint            min_dist ) const;

    //! return true if \a node is local to cluster tree \a cl
    bool          is_local  ( const TCluster *      cl,
                              const node_t          node,
                              const TPermutation *  perm_e2i ) const
    {
        idx_t  idx;
        
        if ( perm_e2i != NULL ) idx = perm_e2i->permute( node );
        else                    idx = node;
        
        return (( idx >= cl->first() ) && ( idx <= cl->last() ));
    }

    DISABLE_COPY_OP( TStdAlgAdmCond );
};

//!
//! \ingroup  Cluster_Module
//! \class    TStdAlgAdmCond
//! \brief    Weak admissibility condition based on matrix graph criteria.
//!
class TWeakAlgAdmCond : public TAlgAdmCond
{
private:
    // upper limit for distance (path length) between clusters (default: 1)
    const uint  _distance;
    
    // upper limit for number of direct edges between clusters (default: 0)
    const uint  _connectivity;
    
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct object for algebraic weak admissibility based
    //! on connectivity in \a S
    //! - \a distance is the upper limit for the path length between clusters
    //! - \a connectivity is the upper limit for direct edges between clusters
    TWeakAlgAdmCond ( any_const_sparse_matrix_t  S,
                      const TPermutation *       perm_i2e     = nullptr,
                      const uint                 distance     = 1,
                      const uint                 connectivity = 0 );

    //! construct object for algebraic weak admissibility with row and column
    //! permutations \a row_perm_i2e and \a col_perm_i2e (from internal to external
    //! ordering)
    TWeakAlgAdmCond ( any_const_sparse_matrix_t  S,
                      const TPermutation *       row_perm_i2e,
                      const TPermutation *       col_perm_i2e,
                      const uint                 distance     = 1,
                      const uint                 connectivity = 0 );

    //! dtor
    virtual ~TWeakAlgAdmCond () {}

    ///////////////////////////////////////////
    //
    // check block-cluster if admissible
    //

    //! return true if \a cl is weakly admissible
    virtual bool is_adm ( const TBlockCluster *  c ) const;
};

}// namespace

#endif  // __HPRO_TALGADMCOND_HH
