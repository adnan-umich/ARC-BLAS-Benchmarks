//
// Project     : HLIBpro
// File        : cluster/TBoundingVolume
// Description : class handling bounding volumes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2024. All Rights Reserved.
//

#include <hpro/config.h>

#include <hpro/cluster/TBoundingVolume.hh>

namespace Hpro
{

///////////////////////////////////////////////
//
// bounding box properties
//

//
// return true if point \a x is inside box
//
bool
TBoundingVolume::is_inside ( const TPoint &  x ) const
{
    return _bbox.is_inside( x );
}
    
//
// return diameter of box
//
double
TBoundingVolume::diameter () const
{
    #if HPRO_USE_CGAL == 0
    return _bbox.diameter();
    #else
    return std::min( _bbox.diameter(), _bsphere.diameter() );
    #endif
}

//
// return volume of box
//
double
TBoundingVolume::volume () const
{
    #if HPRO_USE_CGAL == 0
    return _bbox.volume();
    #else
    return std::min( _bbox.volume(), _bsphere.volume() );
    #endif
}

//
// return distance to \a box
//
double
TBoundingVolume::distance ( const TBoundingVolume &  bvol ) const
{
    #if HPRO_USE_CGAL == 0
    
    return _bbox.distance( bvol.bbox() );
    
    #else

    const auto  dist_box = _bbox.distance( bvol.bbox() );
    const auto  dist_sph = _bsphere.distance( bvol.bsphere() );

    return std::max( dist_box, dist_sph );

    #endif
}

//
// return distance to \a box but coordinates have
// periodicity defined by \a period
//
double
TBoundingVolume::distance ( const TBoundingVolume &  bvol,
                            const TPoint &           period ) const
{
    return _bbox.distance( bvol.bbox(), period );
}

//
// return dimension of intersection with \a bbox
//
uint
TBoundingVolume::overlap_dim ( const TBoundingVolume &  bvol ) const
{
    return _bbox.overlap_dim( bvol.bbox() );
}

//
// extend local bbox by given bounding volume
//
void
TBoundingVolume::extend ( const TBoundingVolume &  /* bvol */ )
{
    HERROR( ERR_NOT_IMPL, "", "" );
}
    
//
// extend local bbox by given bbox
//
void
TBoundingVolume::extend ( const TBBox &  /* bbox */ )
{
    HERROR( ERR_NOT_IMPL, "", "" );
}
    
//
// check volume and adjust if degenerate
//
void
TBoundingVolume::check ()
{
    _bbox.check();
    _bsphere.check();
}

//
// return string representation
//
std::string
TBoundingVolume::to_string () const
{
    return _bbox.to_string();
}

}// namespace Hpro
