//
// Project     : HLIBpro
// File        : TBBox.cc
// Description : class for an axis aligned bounding box
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <hpro/cluster/TBSphere.hh>

#include <hpro/cluster/TBBox.hh>

namespace Hpro
{

//
// ctor
//
TBBox::TBBox ( const TBSphere &  bsphere )
{
    auto  dim    = bsphere.dim();
    auto  center = bsphere.center();
    auto  radius = TPoint( dim );

    for ( uint  i = 0; i < dim; ++i )
        radius[i] = bsphere.radius();

    _bb_min = center - radius;
    _bb_max = center + radius;
}

///////////////////////////////////////////////
//
// bounding box properties
//

//
// return true if given point is inside box
//
bool
TBBox::is_inside ( const TPoint & x ) const
{
    if ( x.dim() != min().dim() )
        HERROR( ERR_DIM, "(TBBox) is_inside", "given point has wrong dimension" );

    const uint   mdim = min().dim();
    const double eps  = 10.0 * Limits::epsilon<double>();

    for ( uint i = 0; i < mdim; i++ )
        if (( x[i] + eps < min()[i] ) || ( max()[i] < x[i] - eps ))
            return false;

    return true;
}
    
//
// return diameter of box
//
double
TBBox::diameter () const
{
    TPoint  d;

    d.assign( 1.0, _bb_max, -1.0, _bb_min );

    return d.norm2();
}

//
// return volume to given box
//
double
TBBox::volume () const
{
    const auto  mdim = min().dim();

    if ( mdim == 0 )
        return 0;
    
    auto  vol = double(1);

    for ( uint i = 0; i < mdim; i++ )
    {
        const auto  min_i = min()[i];
        const auto  max_i = max()[i];

        vol *= ( max_i - min_i );
    }// for

    return vol;
}

//
// return distance to given box
//
double
TBBox::distance ( const TBBox & bbox ) const
{
    //
    // compare all (!) corners of bounding box and
    // choose minimal distance
    //

    TPoint      dist;
    const uint  mdim = min().dim();

    dist.set_dim( mdim );
    
    for ( uint i = 0; i < mdim; i++ )
    {
        const double min0 = min()[i];
        const double max0 = max()[i];
        const double min1 = bbox.min()[i];
        const double max1 = bbox.max()[i];

        if      (min0 > max1) dist[i] = min0 - max1;
        else if (max0 < min1) dist[i] = min1 - max0;
        else                  dist[i] = 0;
    }// for

    return dist.norm2();
}

//
// return distance to given box but coordinates have
// periodicity defined by <period>
//
double
TBBox::distance ( const TBBox & bbox, const TPoint & period ) const
{
    //
    // compare all (!) corners of bounding box, with and
    // without shift and choose minimal distance
    //

    const uint  mdim = min().dim();
    TPoint      dist( mdim );
    double      min0, max0, min1, max1;
    
    HASSERT( mdim == period.dim(), ERR_DIM, "(TBBox) distance",
             "periodicity has different dimension than coordinates" );
    
    for ( uint i = 0; i < mdim; i++ )
    {
        min0 = min()[i];
        max0 = max()[i];

        //
        // compare standard boxes
        //
        
        min1 = bbox.min()[i];
        max1 = bbox.max()[i];

        if      (min0 > max1) dist[i] = min0 - max1;
        else if (max0 < min1) dist[i] = min1 - max0;
        else                  dist[i] = 0;

        //
        // compare given box shifted by <period>
        //
        
        min1 = bbox.min()[i] + period[i];
        max1 = bbox.max()[i] + period[i];

        if      (min0 > max1) dist[i] = std::min( dist[i], min0 - max1 );
        else if (max0 < min1) dist[i] = std::min( dist[i], min1 - max0 );
        else                  dist[i] = 0;

        //
        // compare given box shifted by " - <period> "
        //
        
        min1 = bbox.min()[i] - period[i];
        max1 = bbox.max()[i] - period[i];

        if      (min0 > max1) dist[i] = std::min( dist[i], min0 - max1 );
        else if (max0 < min1) dist[i] = std::min( dist[i], min1 - max0 );
        else                  dist[i] = 0;
    }// for

    return dist.norm2();
}

//
// join local bbox with given bbox
//
void
TBBox::join ( const TBBox & bbox )
{
    // if not initialised, use given box
    if ( min().dim() == 0 )
    {
        *this = bbox;
        return;
    }// if

    if (( min().dim() != bbox.min().dim() ) ||
        ( max().dim() != bbox.max().dim() ))
        HERROR( ERR_ARG, "(TBBox) join", "argument has different spatial dimension" );

    const uint  mdim = min().dim(); // assuming that min and max have same dimension

    for ( uint i = 0; i < mdim; i++ )
    {
        min()[i] = std::min( min()[i], bbox.min()[i] );
        max()[i] = std::max( max()[i], bbox.max()[i] );
    }// for
}
    
//
// extend local bbox by given point
//
void
TBBox::extend ( const TPoint &  p )
{
    if ( min().dim() == 0 )
    {
        _bb_min = p;
        _bb_max = p;
    }// if
    else
    {
        if (( min().dim() != p.dim() ) ||
            ( max().dim() != p.dim() ))
            HERROR( ERR_ARG, "(TBBox) extent", "given point has different spatial dimension" );
        
        const uint  mdim = min().dim(); // assuming that min and max have same dimension
        
        for ( uint i = 0; i < mdim; i++ )
        {
            min()[i] = std::min( min()[i], p[i] );
            max()[i] = std::max( max()[i], p[i] );
        }// for
    }// else
}

//
// return dimension of intersection with \a bbox
uint
TBBox::overlap_dim ( const TBBox &  bbox ) const
{
    //
    // reduce dimension of overlap with every empty axis intersection
    //

    const uint  d      = dim();
    uint        n_over = d;
    
    if ( bbox.dim() != d )
        HERROR( ERR_ARG, "(TBBox) overlap_dim", "different dimension in given bbox" );

    for ( uint  i = 0; i < d; i++ )
    {
        const auto  tmin = min()[i];  // τ
        const auto  tmax = max()[i];
            
        const auto  smin = bbox.min()[i]; // σ
        const auto  smax = bbox.max()[i];

        if (( tmax <= smin ) ||   // ├── τ ──┤├── σ ──┤
            ( smax <= tmin ))     // ├── σ ──┤├── τ ──┤
        {
            n_over--;
        }// if
    }// for
        
    return  n_over;
}

//
// check volume and adjust if degenerate
//
void
TBBox::check ()
{
    const double  BVOL_EPS = 1e-5;
    const auto    diam     = diameter();
    
    for ( uint  i = 0; i < dim(); i++ )
    {
        const auto  diff = max()[i] - min()[i];
        
        if ( diff < BVOL_EPS * diam )
        {
            min()[i] -= 0.5 * ( BVOL_EPS * diam - diff );
            max()[i] += 0.5 * ( BVOL_EPS * diam - diff );
        }// if
    }// for
}

///////////////////////////////////////////////
//
// misc.
//

//
// copy operator
//
TBBox &
TBBox::operator = ( const TBBox & box )
{
    _bb_min = box._bb_min;
    _bb_max = box._bb_max;

    return *this;
}

//
// return string representation
//
std::string
TBBox::to_string () const
{
    return "( " + _bb_min.to_string() + " ... " + _bb_max.to_string() + " )";
}
    
}// namespace
