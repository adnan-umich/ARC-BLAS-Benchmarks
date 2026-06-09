//
// Project     : HLIBpro
// File        : bsp_part_strat.cc
// Description : partitioning strategies for geometrical clustering
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <algorithm>
#include <unordered_map>
#include <set>

#include "hpro/blas/Algebra.hh"

#include "hpro/cluster/TBSPPartStrat.hh"

namespace Hpro
{

using std::set;
using std::pair;
using std::unordered_map;

namespace B = BLAS;

////////////////////////////////////////////////////////////
//
// partition according to geometrical volume
//
////////////////////////////////////////////////////////////

//
// ctor
//

TGeomBSPPartStrat::TGeomBSPPartStrat ( const split_axis_mode_t  split_axis_mode )
        : _split_axis_mode( split_axis_mode )
{}

//
// basic method to partition given indexset
//
void
TGeomBSPPartStrat::partition ( const TCoordinate *     coord,
                               const TNodeSet  &       dofs,
                               TNodeSet &              left,
                               TNodeSet &              right,
                               const TBBox &           bbox,
                               std::vector< TBBox > &  son_bbox,
                               const uint              depth ) const
{
    if ( coord == nullptr )
        HERROR( ERR_ARG, "(TGeomBSPPartStrat) partition", "undefined (nullptr) coordinate set" );
    
    //
    // determine dimension with maximum scale
    //

    uint  max_dim = 0;
    
    if ( _split_axis_mode == regular_split_axis )
    {
        max_dim = depth % coord->dim();
    }// if
    else if ( _split_axis_mode == adaptive_split_axis )
    {
        //
        // choose largest dimension
        //
        
        const uint  dim      = bbox.min().dim();
        double      max_size = -1.0;

        max_dim = 0;
            
        for ( uint i = 0; i < dim; i++ )
        {
            if ( (bbox.max()[i] - bbox.min()[i]) > max_size )
            {
                max_size = bbox.max()[i] - bbox.min()[i];
                max_dim  = i;
            }// if
        }// for
    }// if

    //
    // set refined bounding boxes
    //

    const double  mid = (bbox.max()[max_dim] + bbox.min()[max_dim]) / 2;
    
    son_bbox.resize( 2 );
    
    son_bbox[0] = bbox;
    son_bbox[0].max()[ max_dim ] = mid;
    
    son_bbox[1] = bbox;
    son_bbox[1].min()[ max_dim ] = mid;
    
    //
    // now partition indexset into two parts
    //

    size_t  n_left  = 0;
    size_t  n_right = 0;
    
    for ( auto  dof : dofs )
    {
        if ( coord->coord( dof )[ max_dim ] <= mid )
            ++n_left;
        else
            ++n_right;
    }// for

    left.remove_all();
    right.remove_all();
    
    left.resize( n_left );
    right.resize( n_right );

    for ( auto  dof : dofs )
    {
        if ( coord->coord( dof )[ max_dim ] <= mid )
            left.append( dof );
        else
            right.append( dof );
    }// for
    
}

////////////////////////////////////////////////////////////
//
// partition according to cardinality
//
////////////////////////////////////////////////////////////

namespace
{

template < typename T_compare >
void
seq_sort ( std::vector< node_t > &  dofs,
           const size_t             lb,
           const size_t             ub,
           const T_compare &        less )
{
    std::sort( & dofs[lb], & dofs[ub-1], less );
}

template < typename T_compare >
void
par_sort ( std::vector< node_t > &  dofs,
           const size_t             lb,
           const size_t             ub,
           const T_compare &        less )
{
    if ( ub - lb < 100000 )
    {
        seq_sort( dofs, lb, ub, less );
        return;
    }// if
    
    const auto  mid = ( ub + lb ) / 2;
    
    par_sort( dofs, lb, mid, less );
    par_sort( dofs, mid, ub, less );

    //
    // merge subsets
    //

    // return if [lb,mid-1][mid,ub-1] is already sorted
    if ( less( dofs[mid-1], dofs[mid] ) )
        return;

    std::vector< node_t >  sdofs( ub - lb );
    size_t                 pos = 0;
    auto                   i1  = lb;
    auto                   i2  = mid; 
  
    while (( i1 < mid ) && ( i2 < ub ))
    { 
        // If element 1 is in right place 
        if ( less( dofs[i1], dofs[i2] ) )
        { 
            sdofs[ pos++ ] = dofs[i1];
            ++i1;
        }// if
        else
        { 
            sdofs[ pos++ ] = dofs[i2];
            ++i2;
        }// else
    }// while

    for ( ; i1 < mid; ++i1 )
        sdofs[ pos++ ] = dofs[i1];

    for ( ; i2 < ub; ++i2 )
        sdofs[ pos++ ] = dofs[i2];

    for ( size_t  i = 0; i < sdofs.size(); ++i )
        dofs[ i+lb ] = sdofs[i];
}

}// namespace anonymous

//
// sorting of indices
//
class TIdxCompare
{
protected:
    // array containing vertex positions
    const TCoordinate * _coord;
    
    // dimension we compare
    int                 _cmp_dim;

public:
    // constructor
    TIdxCompare ( const TCoordinate *  coord,
                  const int            cmp_dim )
            : _coord( coord ), _cmp_dim( cmp_dim )
    {
        if ( _coord == nullptr )
            HERROR( ERR_ARG, "(TIdxSort) compare", "no coordinates available" );
    }
    
    // compare indices
    bool  operator ()  ( const node_t &  t1,
                         const node_t &  t2 ) const
    {
        return ( _coord->coord( t1 )[_cmp_dim] < _coord->coord( t2 )[_cmp_dim] );
    }
};

//
// ctor
//

TCardBSPPartStrat::TCardBSPPartStrat ( const split_axis_mode_t  split_axis_mode )
        : _split_axis_mode( split_axis_mode )
{
}

//
// basic method to partition given indexset
//
void
TCardBSPPartStrat::partition ( const TCoordinate *     coord,
                               const TNodeSet  &       dofs,
                               TNodeSet &              left,
                               TNodeSet &              right,
                               const TBBox &           bbox,
                               std::vector< TBBox > &  son_bbox,
                               const uint              depth ) const
{
    if ( coord == nullptr )
        HERROR( ERR_ARG, "(TCardBSPPartStrat) partition", "undefined (nullptr) coordinate set" );
    
    //
    // determine dimension with maximum scale
    //

    uint  max_dim = 0;
    
    if ( _split_axis_mode == regular_split_axis )
    {
        max_dim = depth % coord->dim();
    }// if
    else if ( _split_axis_mode == adaptive_split_axis )
    {
        //
        // choose largest dimension
        //
        
        const uint  dim      = bbox.dim();
        double      max_size = -1.0;

        max_dim = 0;
            
        for ( uint i = 0; i < dim; i++ )
        {
            if ( (bbox.max()[i] - bbox.min()[i]) > max_size )
            {
                max_size = bbox.max()[i] - bbox.min()[i];
                max_dim  = i;
            }// if
        }// for
    }// else

    //
    // sort indices wrt maximal dimension
    //

    TIdxCompare            idxcmp( coord, max_dim );
    std::vector< node_t >  arr_dof( dofs.nnodes() );
    uint                   no = 0;

    for ( auto  dof : dofs )
        arr_dof[no++] = dof;

    // sort( arr_dof.begin(), arr_dof.end(), idxcmp );
    par_sort( arr_dof, 0, no, idxcmp );

    //
    // now partition indexset into 2 parts and update
    // refined bounding boxes
    //

    const size_t  size    = dofs.nnodes();
    const size_t  half    = size / 2;
    double        mid     = bbox.min()[max_dim];
    size_t        n_left  = 0;
    size_t        n_right = 0;
    
    son_bbox.resize( 2 );
    
    son_bbox[0] = bbox;
    son_bbox[1] = bbox;

    for ( size_t  i = 0; i < size; i++ )
    {
        const node_t  id  = arr_dof[i];

        if ( i < half )
        {
            ++n_left;

            // update middle coordinate with new vertex
            mid = std::max( mid, coord->coord( id )[max_dim] );
        }// if
        else
        {
            ++n_right;
        }// else
    }// for

    left.remove_all();
    right.remove_all();
    
    left.resize( n_left );
    right.resize( n_right );
    
    for ( size_t  i = 0; i < size; i++ )
    {
        const node_t  id  = arr_dof[i];

        if ( i < half ) left.append(  id );
        else            right.append( id );
    }// for
    
    son_bbox[0].max()[max_dim] = mid;
    son_bbox[1].min()[max_dim] = mid;
}

////////////////////////////////////////////////////////////
//
// partitioning according to principle component analysis
//
////////////////////////////////////////////////////////////

//
// compare index-angle pairs according to angle
//
bool
angle_cmp  ( const pair< idx_t, double > &  t1,
             const pair< idx_t, double > &  t2 )
{
    return ( t1.second < t2.second );
}

//
// ctor
//
TPCABSPPartStrat::TPCABSPPartStrat ( const bool  use_card )
        : _use_card( use_card )
{}

//
// basic method to partition given indexset
//
void
TPCABSPPartStrat::partition ( const TCoordinate * coord,
                              const TNodeSet &    dofs,
                              TNodeSet &          left,
                              TNodeSet &          right,
                              const TBBox &,
                              std::vector< TBBox > &,
                              const uint ) const
{
    if ( coord == nullptr )
        HERROR( ERR_ARG, "(TPCABSPPartStrat) partition", "undefined (nullptr) coordinate set" );
    
    //
    // determine center of cluster
    //

    const uint    dim   = coord->dim();
    const size_t  ndofs = dofs.nnodes();
    TPoint        center( dim );

    if ( ndofs == 0 )
        return;
    
    for ( auto  dof : dofs )
    {
        const auto  coo = coord->coord( dof );
        
        for ( uint j = 0; j < dim; ++j )
            center[j] += coo[ j ];
    }// for

    center.scale( 1.0 / double( ndofs ) );
    
    //
    // determine principle component
    //
    // let X   = (coord_0,coord_1,…) and
    //     X_c = X - center          (for each column)
    //

    TPoint  dir( dim );
    
    {
        //
        // compute covariance matrix C = 1/n X_c^T·X_c
        //
        
        auto  C = B::Matrix< double >( dim, dim );
        auto  x = B::Vector< double >( dim );

        for ( auto  dof : dofs )
        {
            const auto  coo = coord->coord( dof );
            
            // add (x_i - x_center) · (x_i - x_center)^T
            for ( uint j = 0; j < dim; ++j )
                x(j) = coo[j] - center[j];
            
            B::add_r1( 1.0, x, x, C );
        }// for

        B::scale( 1.0 / double( ndofs ), C );
        
        #if 0

        //
        // compute main direction (first singular vector of C)
        //
        
        B::Vector< double >  sv( dim );
        
        B::svd( C, sv );

        for ( uint j = 0; j < dim; ++j )
            dir[j] = C(j,0);
        
        #else

        //
        // compute main direction (eigen vector to largest eigen value)
        //
        
        auto  eig_vec = B::Matrix< double >( dim, dim );
        auto  eig_val = B::Vector< double >( dim );
        
        B::eigen_herm( C, eig_val, eig_vec );

        const idx_t  max_idx = B::max_idx( eig_val );
        
        for ( uint j = 0; j < dim; ++j )
            dir[j] = eig_vec(j,max_idx);
        
        #endif
    }
    
    //
    // sort coordinates by splitting at angle of main direction and center
    //

    size_t  nleft  = 0;
    size_t  nright = 0;
     
    if ( ! _use_card )
    {
        double  mid_angle = dot( dir, center );

        for ( auto  dof : dofs )
        {
            const auto  coo = coord->coord( dof );
            TPoint      x( dim );

            for ( uint j = 0; j < dim; ++j )
                x[j] = coo[j] - center[j];
        
            if ( dot( dir, x ) <= mid_angle ) ++nleft;
            else                              ++nright;
        }// for

        if (( nleft == 0 ) || ( nright == 0 ))
        {
            //
            // sort coordinates by splitting at center of all angles
            //

            // determine center of angles
            double  cmin  = 0.0;
            double  cmax  = 0.0;
            bool    first = true;
            TPoint  x( dim );

            for ( auto  dof : dofs )
            {
                const auto  coo = coord->coord( dof );

                for ( uint j = 0; j < dim; ++j )
                    x[j] = coo[j] - center[j];
        
                const auto  d = dot( dir, x );

                if ( first )
                {
                    cmin  = d; 
                    cmax  = d;
                    first = false;
                }// if
                else
                {
                    cmin = std::min( cmin, d );
                    cmax = std::max( cmax, d );
                }// else
            }// for

            // sort w.r.t. center of angles
            mid_angle = (cmin + cmax) / 2.0;

            nleft = nright = 0;
        
            for ( auto  dof : dofs )
            {
                const auto  coo = coord->coord( dof );

                for ( uint j = 0; j < dim; ++j )
                    x[j] = coo[j] - center[j];
            
                if ( dot( dir, x ) <= mid_angle ) ++nleft;
                else                              ++nright;
            }// for

            left.remove_all();
            right.remove_all();
    
            left.resize( nleft );
            right.resize( nright );

            for ( auto  dof : dofs )
            {
                const auto  coo = coord->coord( dof );

                for ( uint j = 0; j < dim; ++j )
                    x[j] = coo[j] - center[j];
            
                if ( dot( dir, x ) <= mid_angle )
                    left.append( dof );
                else
                    right.append( dof );
            }// for
        }// if
        else
        {
            left.remove_all();
            right.remove_all();
    
            left.resize( nleft );
            right.resize( nright );
            
            for ( auto  dof : dofs )
            {
                const auto  coo = coord->coord( dof );
                TPoint      x( dim );

                for ( uint j = 0; j < dim; ++j )
                    x[j] = coo[j] - center[j];
        
                if ( dot( dir, x ) <= mid_angle )
                    left.append( dof );
                else
                    right.append( dof );
            }// for
        }// else
    }// if
    
    //
    // if cardinality was chosen or if still one of the subclusters is empty, 
    // sort according to angle coose by cardinality
    //
        
    if (( nleft == 0 ) || ( nright == 0 ))
    {
        auto    idx_angle = std::vector< pair< idx_t, double > >( ndofs );
        size_t  pos = 0;
   
        for ( auto  dof : dofs )
        {
            TPoint  coo( dim, coord->coord( dof ) );

            coo.add( -1.0, center );

            auto  data = std::pair< idx_t, double >( dof, dot( dir, coo ) );

            idx_angle[pos++] = data;
        }// for

        sort( idx_angle.begin(), idx_angle.end(), angle_cmp );
        
        //
        // put first half in <left> and second half in <right>
        //

        const size_t  ndofs_half = ndofs / 2;
    
        pos = 0;

        left.remove_all();
        right.remove_all();
    
        left.resize( ndofs_half );
        right.resize( ndofs - ndofs_half );
        
        for ( auto &  idx : idx_angle )
        {
            if ( pos < ndofs_half ) left.append( idx.first );
            else                    right.append( idx.first );

            ++pos;
        }// for
    }// if
}

////////////////////////////////////////////////////////////
//
// partition optimised for nested dissection
//
////////////////////////////////////////////////////////////

//
// ctor
//
template < typename value_t >
TNDBSPPartStrat< value_t >::TNDBSPPartStrat ( const TSparseMatrix< value_t > *  S,
                                              const edgecut_weights_mode_t      edgecut_weights_mode )
        : _sparse_mat( S ),
          _use_edgecut_weights( edgecut_weights_mode == edgecut_weights_on )
{
    if ( _sparse_mat == nullptr )
        HERROR( ERR_ARG, "(TNDBSPPartStrat) ctor", "sparse matrix is nullptr" );
}

//
// basic method to partition given indexset
//
template < typename value_t >
void
TNDBSPPartStrat< value_t >::partition ( const TCoordinate *     coord,
                                        const TNodeSet  &       dofs,
                                        TNodeSet &              left,
                                        TNodeSet &              right,
                                        const TBBox &           bbox,
                                        std::vector< TBBox > &  son_bbox,
                                        const uint ) const
{
    using  idxmap_t = unordered_map< idx_t, idx_t >;
    
    if ( coord == nullptr )
        HERROR( ERR_ARG, "(TNDBSPPartStrat) partition", "undefined (nullptr) coordinate set" );
    
    //
    // try each dimension and choose the one with minimal edgecut
    //

    const size_t           ndofs       = dofs.nnodes();
    const size_t           half        = ndofs / 2;
    uint                   min_ec_dim  = 0;
    double                 min_edgecut = Limits::max< double >();
    std::vector< node_t >  arr_dof( dofs.nnodes() );
    std::set< idx_t >      local;
    std::vector< idx_t >   left_nodes( dofs.nnodes() );
    std::vector< idx_t >   loc2glo( ndofs );
    idxmap_t               glo2loc;
    idx_t                  pos = 0;
    std::vector< uint >    part( ndofs );

    for ( auto  dof : dofs )
    {
        local.insert( dof );
        loc2glo[ pos ] = dof;
        glo2loc[ dof ] = pos;

        ++pos;
    }// for

    
    for ( uint  dim = 0; dim < coord->dim(); ++dim )
    {
        //
        // sort indices wrt maximal dimension
        //

        auto  idxcmp = TIdxCompare( coord, dim );
        uint  no     = 0;

        for ( auto  dof : dofs )
            arr_dof[no++] = dof;

        sort( arr_dof.begin(), arr_dof.end(), idxcmp );
    
        //
        // partition indexset and compute edgecut
        //

        for ( size_t  i = 0; i < ndofs; i++ )
        {
            if ( i < half ) part[ i ] = 0;
            else            part[ i ] = 1;
        }// for

        //
        // compute edgecut
        //

        double  edgecut = 0.0;

        for ( size_t  i = 0; i < ndofs; i++ )
        {
            if ( part[i] == 0 )
            {
                const node_t  dof = loc2glo[i];
                const idx_t   lb  = _sparse_mat->rowptr( dof );
                const idx_t   ub  = _sparse_mat->rowptr( dof+1 );

                for ( idx_t  l = lb; l < ub; ++l )
                {
                    const idx_t  neigh = _sparse_mat->colind( l );

                    if ( local.find( neigh ) == local.end() )
                        break;

                    if ( part[ glo2loc[ neigh ] ] == 1 )
                    {
                        if ( _use_edgecut_weights )
                            edgecut += Math::abs( _sparse_mat->coeff( l ) );
                        else
                            edgecut++;
                    }// if
                }// for
            }// if
        }// for

        if ( edgecut < min_edgecut )
        {
            min_ec_dim  = dim;
            min_edgecut = edgecut;
        }// if
    }// for

    //
    // now partition again with minimal edgecut dimension
    //

    {
        //
        // sort indices wrt maximal dimension
        //

        auto  idxcmp = TIdxCompare( coord, min_ec_dim );
        uint  no     = 0;

        for ( auto  dof : dofs )
            arr_dof[no++] = dof;

        sort( arr_dof.begin(), arr_dof.end(), idxcmp );
    
        //
        // now partition indexset into 2 parts and update
        // refined bounding boxes
        //

        double  mid    = bbox.min()[min_ec_dim];
        size_t  nleft  = 0;
        size_t  nright = 0;
    
        son_bbox.resize( 2 );
    
        son_bbox[0] = bbox;
        son_bbox[1] = bbox;

        for ( size_t  i = 0; i < ndofs; i++ )
        {
            const node_t  id  = arr_dof[i];

            if ( i < half )
            {
                ++nleft;

                // update middle coordinate with new vertex
                mid = std::max( mid, coord->coord( id )[min_ec_dim] );
            }// if
            else
            {
                ++nright;
            }// else
        }// for

        left.remove_all();
        right.remove_all();
    
        left.resize( nleft );
        right.resize( nright );
        
        for ( size_t  i = 0; i < ndofs; i++ )
        {
            const node_t  id  = arr_dof[i];

            if ( i < half ) left.append(  id );
            else            right.append( id );
        }// for
        
        son_bbox[0].max()[min_ec_dim] = mid;
        son_bbox[1].min()[min_ec_dim] = mid;
    }
}

//
// explicit template instantiation
//

template class TNDBSPPartStrat< float >;
template class TNDBSPPartStrat< double >;
template class TNDBSPPartStrat< std::complex< float > >;
template class TNDBSPPartStrat< std::complex< double > >;

////////////////////////////////////////////////////////////
//
// automatic choice of best partitioning strategy
//
////////////////////////////////////////////////////////////

//
// ctor
//

TAutoBSPPartStrat::TAutoBSPPartStrat ( const split_axis_mode_t  split_axis_mode )
        : _geom( split_axis_mode )
        , _card( split_axis_mode )
{
}

//
// basic method to partition given indexset
//
void
TAutoBSPPartStrat::partition ( const TCoordinate *     coord,
                               const TNodeSet &        dofs,
                               TNodeSet &              left,
                               TNodeSet &              right,
                               const TBBox &           bbox,
                               std::vector< TBBox > &  son_bbox,
                               const uint              depth ) const
{
    _geom.partition( coord, dofs, left, right, bbox, son_bbox, depth );

    //
    // check if clusters are really out of balance
    // and recluster w.r.t. cardinality
    //

    const auto  min_n = std::min( left.nnodes(), right.nnodes() );
    const auto  max_n = std::max( left.nnodes(), right.nnodes() );

    // if difference between max and min is more than 90% of
    // whole set we recluster according to cardinality
    if ( double(max_n - min_n) / double(dofs.size()) > 0.90 )
        _card.partition( coord, dofs, left, right, bbox, son_bbox, depth );
}

}// namespace Hpro
