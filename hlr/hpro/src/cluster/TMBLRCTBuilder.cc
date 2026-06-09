//
// Project     : HLIBpro
// File        : TMBLRCTBuilder.cc
// Description : build clustertrees via MBLR
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <iostream>
#include <list>
#include <deque>
#include <vector>

#include "list.hh"
#include "scheduler.hh"

#include "hpro/cluster/TBSPCTBuilder.hh"

namespace Hpro
{

using std::list;
using std::make_unique;
using std::unique_ptr;

namespace
{

// minimal size for parallel calls
const size_t  MIN_PAR_SIZE       = 2500;

//
// comparison of indices based on coordinates
//
struct idxcoord_cmp_t
{
    // array containing vertex positions
    const TCoordinate & coord;
    
    // dimension we compare
    uint                cmp_dim;

    // ctor
    idxcoord_cmp_t ( const TCoordinate &  acoord,
                     const int            acmp_dim )
            : coord( acoord )
            , cmp_dim( acmp_dim )
    {}
    
    // compare indices
    bool  operator ()  ( const idx_t &  t1,
                         const idx_t &  t2 ) const
    {
        return ( coord.coord( t1 )[cmp_dim] < coord.coord( t2 )[cmp_dim] );
    }
};

//
// split nodes in <dofs> into sets of size <ndest>
//
void
split ( size_t               ndest,
        const TCoordinate &  coord,
        const TNodeSet  &    dofs,
        const TBBox &        bbox,
        list< TNodeSet > &   parts,
        list< TBBox > &      parts_bbox )
{
    //
    // adjust destination size to try to split into similar sized clusters
    //

    const size_t  ndofs = dofs.nnodes();

    ndest = size_t( double(ndofs) / std::round( double(ndofs) / double(ndest) ) );
        
    //
    // for split axis choose largest dimension
    //
    
    uint        split_dim  = 0;
    const uint  dim      = bbox.min().dim();
    double      max_size = -1.0;

    for ( uint i = 0; i < dim; i++ )
    {
        if ( (bbox.max()[i] - bbox.min()[i]) > max_size )
        {
            max_size = bbox.max()[i] - bbox.min()[i];
            split_dim  = i;
        }// if
    }// for

    //
    // sort indices wrt maximal dimension
    //

    idxcoord_cmp_t         idxcmp( coord, split_dim );
    std::vector< idx_t >   arr_dof( dofs.nnodes() );
    uint                   no = 0;

    for ( auto  dof : dofs )
        arr_dof[no++] = dof;

    sort( arr_dof.begin(), arr_dof.end(), idxcmp );
    
    //
    // now partition indexset into 2 parts and update
    // refined bounding boxes
    //

    size_t  npart_dofs = 0;
    idx_t   start_ofs  = 0;
    double  split_pos  = bbox.min()[ split_dim ];
    
    for ( size_t  i = 0; i < ndofs; i++ )
    {
        const node_t  id  = arr_dof[i];
        
        npart_dofs++;
        
        // update middle coordinate with new vertex
        split_pos = std::max( split_pos, coord.coord( id )[split_dim] );
        
        if (( npart_dofs >= ndest ) || ( i == (ndofs-1) ))
        {
            // avoid small remainders
            if (( i < ndofs ) && (( ndofs - i ) < size_t( ndest * 0.1 )))
            {
                npart_dofs += ndofs-i;
                i           = ndofs-1;
            }// if
            
            TNodeSet  subset( npart_dofs );

            for ( size_t  j = start_ofs; j <= i; ++j )
                subset.append( arr_dof[j] );

            parts.push_back( std::move( subset ) );

            auto  sub_bbox = bbox;

            sub_bbox.min()[ split_dim ] = split_pos;
            sub_bbox.max()[ split_dim ] = split_pos;

            parts_bbox.push_back( sub_bbox );

            npart_dofs = 0;
            start_ofs  = i+1;
        }// if
    }// for
}

//
// recursively construct leaf level clusters of (max) size nmin
//
void
spatial_sort ( const size_t             nmin,
               const TBSPPartStrat &    part,
               const TCoordinate &      coord,
               const TNodeSet  &        dofs,
               const TBBox &            bbox,
               const uint               lvl,
               list< idx_t > &          sorted )
{
    if ( dofs.size() <= nmin )
    {
        for ( auto  nodes : dofs )
            sorted.push_back( nodes );

        return;
    }// if
    
    //
    // split nodes
    //

    TNodeSet              left, right;
    std::vector< TBBox >  son_bbox;
    
    part.partition( &coord, dofs, left, right, bbox, son_bbox, lvl );

    //
    // recurse
    //

    if ( left.nnodes() == 0 )
    {
        spatial_sort( nmin, part, coord, right, son_bbox[1], lvl+1, sorted );
    }// if
    else if ( right.nnodes() == 0 )
    {
        spatial_sort( nmin, part, coord, left,  son_bbox[0], lvl+1, sorted );
    }// if
    else
    {
        spatial_sort( nmin, part, coord, left,  son_bbox[0], lvl+1, sorted );
        spatial_sort( nmin, part, coord, right, son_bbox[1], lvl+1, sorted );
    }// else
}

//
// build node for inner cluster with given sons
//
// std::unique_ptr< TGeomCluster >
// build_cluster ( std::list< std::unique_ptr< TGeomCluster > > &  son_clusters )
// {
//     idx_t  min_idx = 0;
//     idx_t  max_idx = 0;
//     TBBox  par_bbox;
//     bool   init    = false;
                
//     for ( auto &  son_cl : son_clusters )
//     {
//         if ( init )
//         {
//             par_bbox.extend( son_cl->bbox() );
//             min_idx = std::min( min_idx, son_cl->first() );
//             max_idx = std::max( max_idx, son_cl->last() );
//         }// if
//         else
//         {
//             par_bbox = son_cl->bbox();
//             min_idx  = son_cl->first();
//             max_idx  = son_cl->last();
//             init     = true;
//         }// else
//     }// for

//     auto   parent = std::make_unique< TGeomCluster >( min_idx, max_idx, par_bbox );
//     uint   spos   = 0;
                
//     parent->set_nsons( son_clusters.size() );

//     for ( auto &  son_cl : son_clusters )
//         parent->set_son( spos++, son_cl.release() );

//     return parent;
// }

}// namespace anonymous

//
// construct BSP cluster tree builder with partitioning strategy \a part_strat
//
TMBLRCTBuilder::TMBLRCTBuilder ( const size_t           nlevel,
                                 const TBSPPartStrat *  part_strat,
                                 const uint             n_min )
        : TGeomCTBuilder( n_min )
        , _part_strat( part_strat )
        , _nlevel( nlevel )
{
}

//
// dtor
//
TMBLRCTBuilder::~TMBLRCTBuilder ()
{
}

//
// build tree out of point cloud
//
unique_ptr< TClusterTree >
TMBLRCTBuilder::build ( const TCoordinate *  coord,
                        const idx_t          index_ofs ) const
{
    if ( coord == nullptr )
        return nullptr;

    //
    // setup nodes and permutation
    //

    const size_t  max_dof = coord->ncoord();
    TNodeSet      dofs( max_dof );
    
    for ( uint i = 0; i < max_dof; i++ )
        dofs.append( i );

    TPermutation  perm_e2i( max_dof );
    const idx_t   UNDEF_IDX = -1;
    

    for ( size_t  i = 0; i < max_dof; ++i )
        perm_e2i[i] = UNDEF_IDX;

    //
    // sort indices based on binary space partitioning
    //
    
    data_t             data = { coord, & perm_e2i, _n_min, _min_leaf_lvl, ( _max_lvl == 0 ? uint(max_dof / 2) : _max_lvl ), 0 };
    // TCardBSPPartStrat  part( adaptive_split_axis );
    list< idx_t >      sorted;
    auto               root_bvol = compute_bvol( dofs, data );

    spatial_sort( 4, *_part_strat, *coord, dofs, root_bvol.bbox(), 0, sorted );

    //
    // leaf level clusters: split index set into requested number of sub sets
    // and construct clusters
    //
    
    const size_t  n        = dofs.nnodes();
    const auto    alpha    = std::pow( double(_n_min) / double(n), 1.0 / double(_nlevel) );
    auto          clusters = std::deque< std::unique_ptr< TGeomCluster > >();

    {
        //
        // combine all nodes in same sub-sequence into single parent cluster
        //

        TNodeSet  cl_idxs( 2*_n_min );
        idx_t     ofs     = index_ofs;
        int       sub_seq = 0;  // current sub sequence
        uint      pos     = 0;  // index of i'th cluster in <clusters>
        
        for ( auto  idx : sorted )
        {
            if (( pos == n-1 ) || ( cl_idxs.nnodes() == _n_min )) 
            {
                // last one not yet in set
                if ( pos == n-1 )
                    cl_idxs.append( idx );

                auto  leaf = build_leaf( cl_idxs, _nlevel, ofs, compute_bvol( cl_idxs, data ), data );

                ofs += leaf->size();
                cl_idxs.remove_all();
                sub_seq++;
                clusters.push_back( std::move( leaf ) );
            }// if

            cl_idxs.append( idx );
            pos++;
        }// for
    }

    //
    // build clusters for remaining hierarchy level with same strategy
    //
    
    for ( int  lvl = _nlevel-1; lvl > 0; --lvl )
    {
        const size_t  ndest = std::round( double( n ) * std::pow( alpha, lvl ) );

        //
        // combine all nodes in same sub-sequence into single parent cluster
        //

        auto    par_clusters = std::deque< std::unique_ptr< TGeomCluster > >();
        auto    son_clusters = std::deque< std::unique_ptr< TGeomCluster > >();
        int     sub_seq      = 0;  // current sub sequence
        uint    pos          = 0;  // index of i'th cluster in <clusters>
        size_t  son_size     = 0;
        
        while ( true )
        {
            // if ( clusters.empty() || ( cl_part[pos] != sub_seq ))
            if ( clusters.empty() || ( son_size + clusters.front()->size() > ndest ))
            {
                // count number of dofs
                size_t  ndof = 0;
                
                for ( auto &  son_cl : son_clusters )
                    ndof += son_cl->size();

                TNodeSet  parent_dofs( ndof );
                idx_t     min_idx, max_idx;
                bool      first = true;
                        
                for ( auto &  son_cl : son_clusters )
                {
                    for ( idx_t  idx = son_cl->first(); idx <= son_cl->last(); ++idx )
                    {
                        parent_dofs.append( idx );

                        if ( first )
                        {
                            min_idx = son_cl->first();
                            max_idx = son_cl->last();
                            first   = false;
                        }// if
                        else
                        {
                            min_idx = std::min( min_idx, son_cl->first() );
                            max_idx = std::max( max_idx, son_cl->last() );
                        }// else
                    }// for
                }// for

                auto  bvol   = compute_bvol( parent_dofs, data );
                auto  parent = std::make_unique< TGeomCluster >( min_idx, max_idx, bvol );

                parent->set_nsons( son_clusters.size() );

                for ( uint  i = 0; i < son_clusters.size(); ++i )
                    parent->set_son( i, son_clusters[i].release() );
                
                par_clusters.push_back( std::move( parent ) );

                son_clusters.clear();
                sub_seq++;
                son_size = 0;

                if ( clusters.empty() )
                    break;
            }// if
            else
            {
                auto  cl = std::move( clusters.front() );

                clusters.pop_front();

                son_size += cl->size();
                son_clusters.push_back( std::move( cl ) );
                pos++;
            }// else
        }// for

        //
        // move to next level
        //

        clusters.clear();
        clusters = std::move( par_clusters );
    }// for
    
    //
    // finally build the root cluster
    //

    if ( clusters.empty() )
        HERROR( ERR_CONSISTENCY, "(TMBLRCTBuilder) build", "no son clusters for root cluster" );
    
    auto  root = std::make_unique< TGeomCluster >( index_ofs, max_dof + index_ofs - 1, root_bvol );

    root->set_nsons( clusters.size() );

    for ( uint  i = 0; i < clusters.size(); ++i )
        root->set_son( i, clusters[i].release() );
    
    //
    // consistency check: are all DoFs covered by clustertree
    //
    
    bool  all_handled = true;
    
    for ( uint i = 0; i < max_dof; i++ )
        if ( perm_e2i[ i ] == UNDEF_IDX )
        {
            perm_e2i[ i ] = 0;
            all_handled   = false;
            HERROR( ERR_CONSISTENCY, "(TMBLRCTBuilder) build", to_string( "DoF %d is not handled", i ) );
        }// if

    if ( ! all_handled )
    {
        root.reset( nullptr );
        return nullptr;
    }// if
    
    //
    // finally put all together in a tree
    //

    auto  permutation_e2i = make_unique< TPermutation >( perm_e2i );
    auto  permutation_i2e = make_unique< TPermutation >( perm_e2i );

    permutation_i2e->invert();
    
    auto  ct = make_unique< TClusterTree >( root.release(), permutation_e2i.release(), permutation_i2e.release() );

    return ct;
}

//
// recursively build cluster tree for indices in \a dofs
//
std::unique_ptr< TGeomCluster >
TMBLRCTBuilder::divide ( const TNodeSet &         dofs,
                         const uint               lvl,
                         const TBoundingVolume &  bvol,
                         const TOptClusterSize &  csize,
                         const idx_t              index_ofs,
                         data_t &                 data ) const
{
    // //
    // // check if number of indices in this cluster is too small to divide
    // // (and minimal leaf level was reached)
    // //

    // if (( dofs.nnodes() <= data.nmin ) || ( lvl >= _nlevel ))
    //     return build_leaf( dofs, lvl, index_ofs, bvol, data );

    // if ( lvl > data.max_lvl )
    // {
    //     // show warning only for default choice of max_lvl
    //     if ( _max_lvl == 0 )
    //         HWARNING( to_string( "in (TMBLRCTBuilder) divide : maximal tree depth reached; depth = %d", lvl ) );
        
    //     return build_leaf( dofs, lvl, index_ofs, bvol, data );
    // }// if
    
    // //
    // // if size of cluster leaves is at most optimal size don't split cluster
    // //

    // if ( ! csize.is_optimal( dofs.nnodes() ) )
    // {
    //     auto  son = divide( dofs, lvl+1, bvol, csize.recurse(), index_ofs, data );

    //     if ( son.get() == nullptr )
    //         HERROR( ERR_NULL, "(TMBLRCTBuilder) divide", "son cluster" );

    //     auto  cluster = make_unique< TGeomCluster >( son->first(), son->last(), son->bvol() );
        
    //     cluster->set_nsons( 1 );
    //     cluster->set_son( 0, son.release() );

    //     return cluster;
    // }// if
    
    // //
    // // divide sons
    // //
    // // - target size for subdivision is determined by number of MBLR levels with 
    // //   target size of leaves = nmin, e.g.
    // //
    // //    - n · α^l = n_min  and therefore  n_i = n·α^i
    // //

    // auto  dof_bbox = bvol.bbox();
    
    // // if ( _adjust_bvol )
    // //     dof_bbox = compute_bbox( dofs, data );

    // const size_t      n       = dofs.nnodes();
    // // const double      alpha   = std::pow( std::log( data.nmin ) / std::log( n ), 1.0 / ( _nlevel - lvl ) );
    // // const size_t      ndest   = std::max< size_t >( std::pow( n, alpha ), data.nmin );
    // const size_t      nblocks = ( lvl != _nlevel ? std::pow( double( n ) / double( data.nmin ), 1.0 / ( _nlevel - lvl ) ) : 1 );
    // const size_t      ndest   = std::round( double( n ) / double( nblocks ) );
    // list< TNodeSet >  subsets;
    // list< TBBox >     subbboxs;

    // // DBG::printf( "size / dest / α : %d / %d / %.4f", n, ndest, alpha );
    // // DBG::printf( "size / dest / nblocks : %d / %d / %d", n, ndest, nblocks );

    // split( ndest, * data.coord, dofs, dof_bbox, subsets, subbboxs );

    // if ( subsets.empty() )
    //     HERROR( ERR_CONSISTENCY, "", "" );
    
    // size_t  nsons = subsets.size();
    
    // // if only single set remains, can't divide further
    // if ( nsons == 1 )
    //     return build_leaf( dofs, lvl, index_ofs, bvol, data );

    // // DBG::printf( "#sons = %d", nsons );
    
    // // for ( auto &  sset : subsets )
    // //     DBG::printf( " son size = %d", sset.nnodes() );
    
    // // // consistency check
    // // if ( son_dofs[0].nnodes() + son_dofs[1].nnodes() != dofs.nnodes() )
    // //     HERROR( ERR_CONSISTENCY, "(TBSPCTBuilder) divide", "lost nodes during divide" );

    // //
    // // copy to vector and set up son index offsets for parallel execution
    // //

    // vector< TNodeSet >  son_dofs( nsons );
    // vector< TBBox >     son_bbox( nsons );
    // vector< idx_t >     son_ofs( nsons );
    // idx_t               ofs = index_ofs;
    // idx_t               idx = 0;

    // while ( ! subsets.empty() )
    // {
    //     TNodeSet         sdofs = behead( subsets );
    //     TBBox  sbbox = behead( subbboxs );

    //     son_ofs[ idx ] = ofs;

    //     ofs += sdofs.nnodes();

    //     son_dofs[ idx ] = std::move( sdofs );
    //     son_bbox[ idx ] = std::move( sbbox );
    //     ++idx;
    // }// for
    
    // //
    // // recursive call for building clustertrees with sons
    // //

    // vector< unique_ptr< TGeomCluster > >  sons( nsons );

    // auto  son_divide = [&son_dofs,&son_bbox,&son_ofs,&csize,&data,lvl,&sons,this] ( const uint  i )
    // {
    //     // DBG::indent( 1 );
        
    //     auto  cl = divide( son_dofs[i], lvl+1, son_bbox[i], csize.recurse(), son_ofs[i], data );
                
    //     // DBG::indent( -1 );
        
    //     if ( cl.get() == nullptr )
    //         HERROR( ERR_NULL, "(TMBLRCTBuilder) divide", "son cluster" );
        
    //     sons[i] = std::move( cl );
    // };
        
    // for ( uint  i = 0; i < nsons; ++i )
    //     son_divide( i );
    
    // //
    // // set bounding box of cluster
    // //

    // if ( _adjust_bbox )
    // {
    //     //
    //     // as union of bb of sons
    //     //

    //     bool  init = false;

    //     for ( auto &  son : sons )
    //     {
    //         if ( init )
    //             dof_bbox.join( son->bbox() );
    //         else
    //         {
    //             dof_bbox = son->bbox();
    //             init     = true;
    //         }// else
    //     }// for
    // }// if
    // else
    // {
    //     //
    //     // as provided by argument and updated with son-bboxes
    //     //
        
    //     for ( auto &  son : sons )
    //         dof_bbox.join( son->bbox() );
    // }// else

    // check_bbox( dof_bbox );

    // //
    // // finally build cluster
    // //

    // idx_t   min_idx = 0, max_idx = 0;

    // {
    //     bool  init = false;
        
    //     for ( auto &  son : sons )
    //     {
    //         if ( init )
    //         {
    //             min_idx = std::min( min_idx, son->first() );
    //             max_idx = std::max( max_idx, son->last() );
    //         }// if
    //         else
    //         {
    //             min_idx = son->first();
    //             max_idx = son->last();
    //             init    = true;
    //         }// else
    //     }// for
    // }
    
    // auto  cluster = make_unique< TGeomCluster >( min_idx, max_idx, dof_bbox );

    // cluster->set_nsons( nsons );

    // for ( uint  i = 0; i < nsons; ++i )
    //     cluster->set_son( i, sons[i].release() );
    
    // return cluster;
}


//
// recursively build cluster tree based on set of clusters
//
std::unique_ptr< TGeomCluster >
TMBLRCTBuilder::divide ( const list< TNodeSet > &         leaves,
                         const list< TBoundingVolume > &  leaf_bbox,
                         const TBoundingVolume &          bbox,
                         const uint                       lvl,
                         const idx_t                      index_ofs,
                         data_t &                         data ) const
{
    // //
    // // divide set of leaves into subsets based on target size
    // //
    // // - target size for subdivision is determined by number of MBLR levels with 
    // //   target size of leaves = nmin, e.g.
    // //
    // //    - n_i^α^(l - i) = n_min  with  α due to n^α^l = n_min   or
    // //    - n_(i-1) / b   = n_i    with  b due to n_min = n / b^l
    // //

    // // count number of nodes
    // size_t  n = 0;

    // for ( auto &  cl : leaves )
    //     n += cl.nnodes();

    // const size_t  nblocks = ( lvl != _nlevel ? std::pow( double( n ) / double( data.nmin ), 1.0 / ( _nlevel - lvl ) ) : 1 );
    // const size_t  ndest   = std::round( double( n ) / double( nblocks ) );

    // //
    // // if single leaf remaining or maximal level reached: construct cluster
    // //

    // if ( leaves.size() == 1 )
    //     return build_leaf( leaves.front(), lvl, index_ofs, leaf_bbox.front(), data );

    // if ( lvl == _nlevel )
    // {
    //     size_t           ndofs   = 0;
    //     TBBox  cl_bbox = leaf_bbox.front();

    //     for ( auto &  cl : leaves )
    //         ndofs += cl.nnodes();

    //     for ( auto &  bb : leaf_bbox )
    //         cl_bbox.join( bb );
               
    //     TNodeSet  dofs( ndofs );
        
    //     for ( auto &  cl : leaves )
    //         for ( auto  i : cl )
    //             dofs.append( i );

    //     return build_leaf( dofs, lvl, index_ofs, cl_bbox, data );
    // }// if
         
    // //
    // // divide cluster
    // //
    
    // auto   sons      = list< std::unique_ptr< TGeomCluster > >();
    // auto   iter_cl   = leaves.begin();
    // auto   iter_bbox = leaf_bbox.begin();
    // idx_t  ofs       = index_ofs;

    // while ( iter_cl != leaves.end() )
    // {
    //     size_t            son_size = 0;
    //     list< TNodeSet >  son_leaves;
    //     list< TBBox >     son_bboxes;
    //     TBBox             son_bbox;
    //     bool              init = false;

    //     while (( son_size < ndest ) && ( iter_cl != leaves.end() ))
    //     {
    //         // at least one cluster but otherwise do not exceed ndest
    //         if ( ! son_leaves.empty() && (( son_size + (*iter_cl).nnodes() ) > ndest ))
    //             break;
            
    //         if ( ! init )
    //         {
    //             son_bbox = *iter_bbox;
    //             init     = true;
    //         }// if
    //         else
    //             son_bbox.join( *iter_bbox );

    //         son_size += (*iter_cl).nnodes();
    //         son_leaves.push_back( *iter_cl++ );
    //         son_bboxes.push_back( *iter_bbox++ );
    //     }// while

    //     auto  son = divide( son_leaves, son_bboxes, son_bbox, lvl+1, ofs, data );

    //     if ( son.get() == nullptr )
    //         HERROR( ERR_NULL, "(TMBLRCTBuilder) divide", "son cluster" );
        
    //     sons.push_back( std::move( son ) );

    //     ofs += son_size;
    // }// while
        

    // //
    // // set bounding box of cluster
    // //

    // auto  cl_bbox = bbox;
    
    // if ( _adjust_bbox )
    // {
    //     //
    //     // as union of bb of sons
    //     //

    //     bool  init = false;

    //     for ( auto &  son : sons )
    //     {
    //         if ( init )
    //             cl_bbox.join( son->bbox() );
    //         else
    //         {
    //             cl_bbox = son->bbox();
    //             init    = true;
    //         }// else
    //     }// for
    // }// if
    // else
    // {
    //     //
    //     // as provided by argument and updated with son-bboxes
    //     //
        
    //     for ( auto &  son : sons )
    //         cl_bbox.join( son->bbox() );
    // }// else

    // check_bbox( cl_bbox );

    // //
    // // finally build cluster
    // //

    // idx_t   min_idx = 0, max_idx = 0;

    // {
    //     bool  init = false;
        
    //     for ( auto &  son : sons )
    //     {
    //         if ( init )
    //         {
    //             min_idx = std::min( min_idx, son->first() );
    //             max_idx = std::max( max_idx, son->last() );
    //         }// if
    //         else
    //         {
    //             min_idx = son->first();
    //             max_idx = son->last();
    //             init    = true;
    //         }// else
    //     }// for
    // }
    
    // auto  cluster = std::make_unique< TGeomCluster >( min_idx, max_idx, cl_bbox );
    // uint  pos     = 0;

    // cluster->set_nsons( sons.size() );

    // for ( auto &  son : sons )
    //     cluster->set_son( pos++, son.release() );
    
    // return cluster;
}

}// namespace Hpro
