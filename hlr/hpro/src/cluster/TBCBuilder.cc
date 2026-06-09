//
// Project     : HLIBpro
// File        : TBCBuilder.cc
// Description : builds a block-cluster-tree
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <mutex>

#include "hpro/parallel/NET.hh"

#include "hpro/cluster/TBCBuilder.hh"

namespace Hpro
{

// namespace
// {

// //
// // return globally unique ID
// //
// int
// get_id ()
// {
//     static std::mutex  id_mutex;
//     static int         id_counter = 0;
//     int                ret_val    = 0;

//     {
//         std::lock_guard< std::mutex >  lock( id_mutex );

//         ret_val = id_counter++;
//     }

//     return ret_val;
// }

// }// namespace anonymous

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// TBCBuilder (implementation)
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//
// construct block cluster tree builder
//
TBCBuilder::TBCBuilder ( const uint                  min_lvl,
                         const cluster_level_mode_t  cluster_lvl_mode )
        : _min_lvl( min_lvl )
        , _same_cluster_level( cluster_lvl_mode == cluster_level_same )
{}

////////////////////////////////////////////////
//
// build a block-clustertree
//

std::unique_ptr< TBlockClusterTree >
TBCBuilder::build ( const TClusterTree *   rowct,
                    const TClusterTree *   colct,
                    const TAdmCondition *  ac ) const
{
    if ((rowct == nullptr) || (colct == nullptr))
        HERROR( ERR_ARG, "(TBCBuilder) build", "cluster tree arguments are nullptr" );

    const TCluster * rowcl = rowct->root();
    const TCluster * colcl = colct->root();
    
    //
    // call recursive procedure for building the tree
    //

    auto  root = build( rowcl, colcl, ac );

    return std::make_unique< TBlockClusterTree >( root.release(), rowct, colct );
}

std::unique_ptr< TBlockCluster >
TBCBuilder::build ( const TCluster *       rowcl,
                    const TCluster *       colcl,
                    const TAdmCondition *  ac ) const
{
    if ((rowcl == nullptr) || (colcl == nullptr))
        HERROR( ERR_ARG, "(TBCBuilder) build", "cluster arguments are nullptr" );

    //
    // call recursive procedure for building the tree
    //

    auto  id   = std::atomic< int >( 0 );
    auto  root = create_bc( nullptr,
                            const_cast< TCluster * >( rowcl ),
                            const_cast< TCluster * >( colcl ) );
    
    rec_build( root.get(), ac, 0, id );

    return root;
}
////////////////////////////////////////////////
//
// recusivly build a block-clustertree
//

void
TBCBuilder::rec_build ( TBlockCluster *        bc,
                        const TAdmCondition *  ac,
                        const uint             level,
                        std::atomic< int > &   id ) const
{
    auto  rowcl = bc->rowcl();
    auto  colcl = bc->colcl();

    if ((rowcl == nullptr) || (colcl == nullptr))
        HERROR( ERR_NULL, "(TBCBuilder) rec_build", "block has NULL clusters" );
    
    //
    // decide if blockcluster is a leaf or if to refine
    //

    if ( ac->is_adm( bc ) && ( level >= _min_lvl ))
    {
        // mark as admissible
        bc->set_adm( true );
        bc->make_leaf();
    }// if
    else if (( rowcl != colcl ) && rowcl->is_domain() && colcl->is_domain() )
    {
        // nested dissection case: offdiagonal domain-domain cluster
        bc->set_adm( true );
        bc->make_leaf();
    }// if
    else if ( _same_cluster_level && ( rowcl->is_leaf() || colcl->is_leaf() ))
    {
        // stop if one of them is leaf
        bc->make_leaf();
    }// if
    else if ( rowcl->is_leaf() && colcl->is_leaf() ) 
    {
        // stop if both of them are leaves
        bc->make_leaf();
    }// else
    else
    {
        //
        // refine block cluster
        //

        refine( bc );
    
        //
        // recurse
        //

        if ( ! bc->is_leaf() )
        {
            const auto  n_rowcl = std::max< size_t >( rowcl->nsons(), 1 );
            const auto  n_colcl = std::max< size_t >( colcl->nsons(), 1 );

            bc->set_layout( n_rowcl, n_colcl );
        
            for ( size_t  i = 0; i < n_rowcl; ++i )
            {
                for ( size_t  j = 0; j < n_colcl; ++j )
                {
                    auto *  son_ij = bc->son( i, j );
                
                    if ( son_ij == nullptr )
                        HERROR( ERR_NULL, "(TBCBuilder) rec_build", "son is nullptr" );
                
                    rec_build( son_ij, ac, level+1, id );
                }// for
            }// for
        }// if
    }// else

    //
    // finally set ID to ensure ID(parent) >= max(id(sons))
    //

    bc->set_id( id++ );
}
                           
//
// create sons for given block cluster
//
void
TBCBuilder::refine ( TBlockCluster *  bc ) const
{
    auto  rowcl = bc->rowcl();
    auto  colcl = bc->colcl();

    if (( rowcl == nullptr ) || ( colcl == nullptr ))
        HERROR( ERR_NULL, "(TBCBuilder) refine", "block has NULL clusters" );

    if ( rowcl->is_leaf() && colcl->is_leaf() )
    {
        // std::cout << "row/col is leaf" << std::endl;
        return;
    }// if
    
    const auto  n_rowcl = std::max< size_t >( rowcl->nsons(), 1 );
    const auto  n_colcl = std::max< size_t >( colcl->nsons(), 1 );

    bc->set_layout( n_rowcl, n_colcl );

    for ( size_t  i = 0; i < n_rowcl; ++i )
    {
        auto  rowcl_i = (rowcl->is_leaf() ? rowcl : rowcl->son(i));

        if ( rowcl_i == nullptr )
            HERROR( ERR_NULL, "(TBCBuilder) refine", "son of row-cluster is nullptr" );
            
        for ( size_t  j = 0; j < n_colcl; ++j )
        {
            auto  colcl_j = (colcl->is_leaf() ? colcl : colcl->son(j));
                    
            if ( colcl_j == nullptr )
                HERROR( ERR_NULL, "(TBCBuilder) refine", "son of column-cluster is nullptr" );
            
            auto  son = create_bc( bc, rowcl_i, colcl_j );

            bc->set_son( i, j, son.release() );
        }// for
    }// for
}

//
// create new object for a block-cluster
//
std::unique_ptr< TBlockCluster >
TBCBuilder::create_bc ( TBlockCluster * parent ) const
{
    return std::make_unique< TBlockCluster >( parent );
}

//
// create new object for a block-cluster
//
std::unique_ptr< TBlockCluster >
TBCBuilder::create_bc ( TBlockCluster *  parent,
                        TCluster *       rowcl,
                        TCluster *       colcl ) const
{
    auto  bc = create_bc( parent );

    bc->set_rowcl( rowcl );
    bc->set_colcl( colcl );

    // assign to all processors by default
    bc->set_procs( TProcSet( 0, NET::nprocs()-1 ) );

    return bc;
}

}// namespace Hpro
