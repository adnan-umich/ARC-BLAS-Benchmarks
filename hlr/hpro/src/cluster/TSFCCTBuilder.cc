//
// Project     : HLIBpro
// File        : TSFCCTBuilder.cc
// Description : build clustertrees via space filling curves
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2023. All Rights Reserved.
//

#include <array>
#include <algorithm>

#include <hpro/config.h>

#if HPRO_USE_CGAL == 1
#  include <CGAL/hilbert_sort.h>
#endif

#include "hpro/cluster/TBSPCTBuilder.hh"

namespace Hpro
{

namespace
{

// undefined index (for permutation consistency check)
constexpr idx_t   UNDEF_IDX    = -1;
    
//
// represents point with id for sorting
//
template < int dim >
struct id_point_t
{
    double  coord[dim];
    idx_t   id;
    
    id_point_t ( double  ax,
                 idx_t   aid )
            : id(aid)
    {
        static_assert( dim == 1, "dimension mismatch" );

        coord[0] = ax;
    }
    
    id_point_t ( double  ax,
                 double  ay,
                 idx_t   aid )
            : id(aid)
    {
        static_assert( dim == 2, "dimension mismatch" );

        coord[0] = ax;
        coord[1] = ay;
    }
    
    id_point_t ( double  ax,
                 double  ay,
                 double  az,
                 idx_t   aid )
            : id(aid)
    {
        static_assert( dim == 3, "dimension mismatch" );

        coord[0] = ax;
        coord[1] = ay;
        coord[2] = az;
    }

    double  operator [] ( const uint  i ) const { return coord[i]; }
};

template < int dim >
struct less_x {
    bool operator () ( const id_point_t< dim > &  p, 
                       const id_point_t< dim > &  q ) const { return p[0] < q[0]; }
};

template < int dim >
struct less_y {
    bool operator () ( const id_point_t< dim > &  p, 
                       const id_point_t< dim > &  q ) const { return p[1] < q[1]; }
};

template < int dim >
struct less_z {
    bool operator () ( const id_point_t< dim > &  p, 
                       const id_point_t< dim > &  q ) const { return p[2] < q[2]; }
};

struct sorting_traits_1_t
{
    using  Point_1  = id_point_t< 1 >;
    using  Less_x_1 = less_x< 1 >;
    
    Less_x_1  less_x_1_object() const { return Less_x_1(); }
};

struct sorting_traits_2_t
{
    using  Point_2  = id_point_t< 2 >;
    using  Less_x_2 = less_x< 2 >;
    using  Less_y_2 = less_y< 2 >;

    Less_x_2  less_x_2_object() const { return Less_x_2(); }
    Less_y_2  less_y_2_object() const { return Less_y_2(); }
};

struct sorting_traits_3_t
{
    using  Point_3  = id_point_t< 3 >;
    using  Less_x_3 = less_x< 3 >;
    using  Less_y_3 = less_y< 3 >;
    using  Less_z_3 = less_z< 3 >;
    
    Less_x_3  less_x_3_object() const { return Less_x_3(); }
    Less_y_3  less_y_3_object() const { return Less_y_3(); }
    Less_z_3  less_z_3_object() const { return Less_z_3(); }
};

}// namespace anonymous

//
// ctor
//
TSFCCTBuilder::TSFCCTBuilder ( const cluster_type_t  cl_type,
                               const split_type_t    split_type,
                               const uint            n_min,
                               const uint            min_leaf_lvl )
        : TGeomCTBuilder( n_min, min_leaf_lvl )
        , _cl_type( cl_type )
        , _split_type( split_type )
{}

//
// dtor
//
TSFCCTBuilder::~TSFCCTBuilder ()
{}

//////////////////////////////////////////////
//
// clustertree creation
//

//
// build cluster tree out of given coordinate set <coord>
//
std::unique_ptr< TClusterTree >
TSFCCTBuilder::build  ( const TCoordinate *  coord,
                        const idx_t          idx_ofs ) const
{
    #if HPRO_USE_CGAL == 1
    
    if ( coord == nullptr )
        return nullptr;

    //
    // setup array for son assignment and permutation
    //

    const size_t    ncoord = coord->ncoord();
    TPermutation    perm_e2i( ncoord );

    for ( size_t  i = 0; i < ncoord; ++i )
        perm_e2i[i] = UNDEF_IDX;
    
    //
    // build root of clustertree (this node ;)
    // copy indices into local list and divide
    //

    auto              root = std::unique_ptr< TGeomCluster >();
    data_t            data = { coord, & perm_e2i, _n_min, _min_leaf_lvl, ( _max_lvl == 0 ? uint(ncoord / 2) : _max_lvl ) };
    TNodeSet          dofs( ncoord );
    TBoundingVolume   bvol;
    TOptClusterSize   csize;

    //
    // use CGAL to sort coordinates (and dofs) by SFC
    //

    if ( coord->dim() == 1 )
    {
        using  point_t = id_point_t< 1 >;

        auto  points = std::vector< point_t >();

        points.reserve( ncoord );

        for ( size_t  i = 0; i < ncoord; ++i )
        {
            auto  c = coord->coord( i );
            auto  p = point_t( c[0], i );

            points.push_back( p );
        }// for
        
        std::sort( points.begin(), points.end(), [] ( const auto & p, const auto & q ) { return p[0] < q[0]; }  );
            
        for ( auto &  p : points )
            dofs.append( p.id );
    }// if
    else if ( coord->dim() == 2 )
    {
        using  point_t = id_point_t< 2 >;

        auto  points = std::vector< point_t >();

        points.reserve( ncoord );

        for ( size_t  i = 0; i < ncoord; ++i )
        {
            auto  c = coord->coord( i );
            auto  p = point_t( c[0], c[1], i );

            points.push_back( p );
        }// for
        
        auto  sst = sorting_traits_2_t();

        CGAL::hilbert_sort( points.begin(), points.end(), sst );

        for ( auto &  p : points )
            dofs.append( p.id );
    }// if
    else if ( coord->dim() == 3 )
    {
        using  point_t = id_point_t< 3 >;

        auto  points = std::vector< point_t >();

        points.reserve( ncoord );

        for ( size_t  i = 0; i < ncoord; ++i )
        {
            auto  c = coord->coord( i );
            auto  p = point_t( c[0], c[1], c[2], i );

            points.push_back( p );
        }// for
        
        auto  sst = sorting_traits_3_t();

        CGAL::hilbert_sort( points.begin(), points.end(), sst );

        for ( auto &  p : points )
            dofs.append( p.id );
    }// if
    else
        HERROR( ERR_DIM, "(TSFCBSPCTBuilder) build", "unsupported coordinate dimension" );
    
    //
    // divide nodeset
    //
    
    // compute bounding volume of root cluster
    bvol = compute_bvol( dofs, data );
    
    // and start partitioning
    root = divide( dofs, 0, bvol, csize, idx_ofs, data );
    
    if ( root.get() == nullptr )
        HERROR( ERR_NULL, "(TGeomCTBuilder) build", "root cluster" );
    
    //
    // consistency check: are all DoFs covered by clustertree
    //
    
    bool  all_handled = true;
    
    for ( uint i = 0; i < ncoord; i++ )
        if ( perm_e2i[ i ] == UNDEF_IDX )
        {
            perm_e2i[ i ] = 0;
            all_handled   = false;
            HERROR( ERR_CONSISTENCY, "(TGeomCTBuilder) build",
                    to_string( "DoF %d is not handled", i ) );
        }// if

    if ( ! all_handled )
    {
        root.reset( nullptr );
        return nullptr;
    }// if
    
    //
    // finally put all together in a tree
    //

    auto  permutation_e2i = std::make_unique< TPermutation >( perm_e2i );
    auto  permutation_i2e = std::make_unique< TPermutation >( perm_e2i );

    permutation_i2e->invert();
    
    auto  ct = std::make_unique< TClusterTree >( root.release(), permutation_e2i.release(), permutation_i2e.release() );

    return ct;

    #else

    HERROR( ERR_NOCGAL, "(TSFCCTBuilder) build", "" );
    
    #endif
}

//
// recursively build cluster tree for indices in <dofs>
//
std::unique_ptr< TGeomCluster >
TSFCCTBuilder::divide ( const TNodeSet &         dofs,
                        const uint               lvl,
                        const TBoundingVolume &  bvol,
                        const TOptClusterSize &  csize,
                        const idx_t              index_ofs,
                        data_t &                 data ) const
{
    //
    // check if number of indices in this cluster is too small to divide
    // (and minimal leaf level was reached)
    //

    if (( lvl >= data.min_leaf_lvl ) && ( dofs.nnodes() <= data.nmin ))
        return build_leaf( dofs, lvl, index_ofs, bvol, data );

    if ( lvl > data.max_lvl )
    {
        // show warning only for default choice of max_lvl
        if ( _max_lvl == 0 )
            HWARNING( to_string( "in (TBSPCTBuilder) divide : maximal tree depth reached; depth = %d", lvl ) );
        
        return build_leaf( dofs, lvl, index_ofs, bvol, data );
    }// if

    //
    // decide upon partitioning type how to split
    // ASSUMPTION: all dofs are ordered according to SFC and therefore, just split sets in half
    //

    if ( _cl_type == binary )
    {
        //
        // if size of cluster leaves is at most optimal size don't split cluster
        //

        if ( ! csize.is_optimal( dofs.nnodes() ) )
        {
            auto  son = divide( dofs, lvl+1, bvol, csize.recurse(), index_ofs, data );

            if ( son.get() == nullptr )
                HERROR( ERR_NULL, "(TBSPCTBuilder) divide", "son cluster" );

            auto  cluster = std::make_unique< TGeomCluster >( son->first(), son->last(), son->bvol() );
        
            cluster->set_nsons( 1 );
            cluster->set_son( 0, son.release() );

            return cluster;
        }// if
    
        //
        // divide into two sons
        //

        auto  [ son_sets, son_bvol ] = split( 2, dofs, data );
        
        HDEBUG( to_string( "(TSFCCTBuilder) divide : left = %d, right = %d", son_sets[0].nnodes(), son_sets[1].nnodes() ) );

        // if one of the subsets is empty => continue with this cluster
        if      ( son_sets[0].nnodes() == 0 ) return divide( son_sets[1], lvl+1, son_bvol[1], csize.recurse(), index_ofs, data );
        else if ( son_sets[1].nnodes() == 0 ) return divide( son_sets[0], lvl+1, son_bvol[0], csize.recurse(), index_ofs, data );

        // consistency check
        if ( son_sets[0].nnodes() + son_sets[1].nnodes() != dofs.nnodes() )
            HERROR( ERR_CONSISTENCY, "(TSFCCTBuilder) divide", "lost nodes during divide" );
    
        //
        // recursive call for building clustertrees with sons
        //

        auto         sons       = std::array< std::unique_ptr< TGeomCluster >, 2 >();
        const idx_t  son_ofs[2] = { index_ofs, index_ofs + idx_t( son_sets[0].nnodes() ) };

        auto  build_soncl =
            [this,lvl,&csize,&data] ( TNodeSet &               sdofs,
                                      const TBoundingVolume &  sbvol,
                                      const idx_t              sofs ) -> std::unique_ptr< TGeomCluster >
            {
                auto  cl = divide( sdofs, lvl+1, sbvol, csize.recurse(), sofs, data );
                
                if ( cl.get() == nullptr )
                    HERROR( ERR_NULL, "(TSFCCTBuilder) divide", "son cluster" );
                
                sdofs.remove_all();
                sdofs.resize( 0 );

                return cl;
            };

        sons[0] = build_soncl( son_sets[0], son_bvol[0], son_ofs[0]  );
        sons[1] = build_soncl( son_sets[1], son_bvol[1], son_ofs[1] );

        // //
        // // set bounding box of cluster
        // //

        // if ( _adjust_bvol )
        // {
        //     //
        //     // as union of bb of sons
        //     //

        //     cl_bvol = sons[0]->bvol();
        //     cl_bvol.join( sons[1]->bvol() );
        // }// if
        // else
        // {
        //     //
        //     // as provided by argument and updated with son-bvoles
        //     //
        
        //     cl_bvol.join( sons[0]->bvol() );
        //     cl_bvol.join( sons[1]->bvol() );
        // }// else

        // check_bvol( cl_bvol );

        //
        // finally build cluster
        //
    
        auto  cluster = std::make_unique< TGeomCluster >( std::min( sons[0]->first(), sons[1]->first() ),
                                                          std::max( sons[0]->last(),  sons[1]->last()  ),
                                                          bvol );

        cluster->set_nsons( 2 );
        cluster->set_son( 0, sons[0].release() );
        cluster->set_son( 1, sons[1].release() );
    
        return cluster;
    }// if
    else if ( _cl_type == blr )
    {
        //
        // divide into n/nmin sons
        //

        // auto  cl_bvol = bvol;

        // if ( _adjust_bvol )
        //     cl_bvol = compute_bvol( dofs, data );

        const uint  nsons                  = dofs.nnodes() / data.nmin;
        auto        [ son_sets, son_bvol ] = split( nsons, dofs, data );
        auto        son_ofs                = std::vector< idx_t >( nsons );
        size_t      pos                    = 0;

        for ( uint  son_idx = 0; son_idx < nsons; ++son_idx )
        {
            son_ofs[ son_idx ] = index_ofs + pos;
            pos               += son_sets[ son_idx ].size();
        }// for
        
        //
        // build son clusters
        //

        auto  sons = std::vector< std::unique_ptr< TGeomCluster > >( nsons );

        for ( uint  i = 0; i < nsons; ++i )
            sons[i] = build_leaf( son_sets[i], lvl+1, son_ofs[i], son_bvol[i], data );
        
        // //
        // // set bounding box of cluster
        // //

        // if ( _adjust_bvol )
        // {
        //     //
        //     // as union of bb of sons
        //     //

        //     cl_bvol = son_bvol[0];

        //     for ( uint  i = 1; i < nsons; ++i )
        //         cl_bvol.join( son_bvol[i] );
        // }// if
        // else
        // {
        //     //
        //     // as provided by argument and updated with son-bvoles
        //     //
        
        //     for ( uint  i = 0; i < nsons; ++i )
        //         cl_bvol.join( son_bvol[i] );
        // }// else

        // check_bvol( cl_bvol );

        //
        // finally build cluster
        //

        auto  first   = sons[0]->first();
        auto  last    = sons[0]->last();

        for ( uint  i = 1; i < nsons; ++i )
        {
            first = std::min( first, sons[i]->first() );
            last  = std::max( last,  sons[i]->last() );
        }// for
        
        auto  cluster = std::make_unique< TGeomCluster >( first, last, bvol );

        cluster->set_nsons( nsons );

        for ( uint  i = 0; i < nsons; ++i )
            cluster->set_son( i, sons[i].release() );
    
        return cluster;
    }// if
    else
        HERROR( ERR_CONSISTENCY, "(TSFCTBuilder) divide", "unknown partition type" );
}

//
// split into <nsub> sub sets
//
std::pair< std::vector< TNodeSet >,
           std::vector< TBoundingVolume > >
TSFCCTBuilder::split ( const uint        nsub,
                       const TNodeSet &  dofs,
                       data_t &          data ) const
{
    auto  son_sets = std::vector< TNodeSet >( nsub );
    auto  son_bvol = std::vector< TBoundingVolume >( nsub, TBoundingVolume() );
              
    if ( _split_type == cardinality )
    {
        //
        // split first nsub-1 into parts of size <sub_size>
        //
        
        const uint  sub_size = dofs.nnodes() / nsub;
        size_t      pos      = 0;

        for ( uint  son_idx = 0; son_idx < nsub-1; ++son_idx )
        {
            son_sets[ son_idx ].resize( sub_size );
            
            for ( idx_t  dof_idx = 0; dof_idx < sub_size; ++dof_idx, ++pos )
                son_sets[ son_idx ].append( dofs[ pos ] );

            son_bvol[ son_idx ] = compute_bvol( son_sets[ son_idx ], data );
        }// for

        //
        // last son gets the rest (to handle case of nsub * sub_size != dofs.size()) 
        //

        son_sets[ nsub-1 ].resize( dofs.nnodes() - pos );
            
        while ( pos < dofs.nnodes() )
        {
            son_sets[ nsub-1 ].append( dofs[ pos++ ] );
        }// else

        son_bvol[ nsub-1 ] = compute_bvol( son_sets[ nsub-1 ], data );
    }// if
    else // if (( _split_type == volume ) || ( _split_type == diameter ))
    {
        //
        // FOR NOW: assume two sub sets
        //

        if (( nsub != 2 ) || ( dofs.size() < 2 ))
            HERROR( ERR_NOT_IMPL, "", "" );

        //
        // start with single dof per subset and increase smaller subset
        // per iteration step
        //

        size_t  nassigned = 0;
        auto    tsets     = std::vector< TNodeSet >( nsub );
        auto    sub_pos   = std::vector< idx_t >( nsub );

        tsets[0].resize( dofs.size() );
        tsets[1].resize( dofs.size() );
        
        sub_pos[0] = 0;
        sub_pos[1] = dofs.size()-1;
        
        while ( nassigned < dofs.size() )
        {
            auto  vol0 = son_bvol[0].volume();
            auto  vol1 = son_bvol[1].volume();

            if (( std::min( vol0, vol1 ) < 1e-14 ) || ( _split_type == diameter ))
            {
                //
                // in case of degenerate boxes or if splitting per diameter, revert to diameter
                //
                
                vol0 = son_bvol[0].diameter();
                vol1 = son_bvol[1].diameter();
            }// if
            
            if ( vol0 <= vol1 )
            {
                tsets[0].append( dofs[ sub_pos[0] ] );
                // son_bvol[0].extend( TPoint( data.coord->dim(), data.coord->coord( dofs[ sub_pos[0] ] ) ) );
                son_bvol[0] = compute_bvol( tsets[0], data ); // TODO: could be expensive!
                sub_pos[0]++;
            }// if
            else
            {
                tsets[1].append( dofs[ sub_pos[1] ] );
                // son_bvol[1].extend( TPoint( data.coord->dim(), data.coord->coord( dofs[ sub_pos[1] ] ) ) );
                son_bvol[1] = compute_bvol( tsets[1], data ); // TODO: could be expensive!
                sub_pos[1]--;
            }// else

            nassigned++;
        }// while
        
        //
        // copy to real sub sets
        //

        for ( uint  son_idx = 0; son_idx < nsub; ++son_idx )
        {
            son_sets[ son_idx ].resize( tsets[ son_idx ].size() );

            for ( auto  dof : tsets[ son_idx ] )
                son_sets[ son_idx].append( dof );
        }// for
    }// else

    return { std::move( son_sets ), std::move( son_bvol ) };
}

}// namespace Hpro
