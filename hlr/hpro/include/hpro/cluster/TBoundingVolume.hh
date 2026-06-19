#ifndef __HPRO_TBOUNDINGVOLUME_HH
#define __HPRO_TBOUNDINGVOLUME_HH
//
// Project     : HLIBpro
// File        : TBoundingVolume.hh
// Description : class handling bounding volumes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2024. All Rights Reserved.
//

#include <hpro/cluster/TBBox.hh>
#include <hpro/cluster/TBSphere.hh>

namespace Hpro
{

/////////////////////////////////////////////////////////////////
//
// using bounding box and bounding sphere together to optimise
// distance/diameter
//
class TBoundingVolume
{
private:
    TBBox     _bbox;
    TBSphere  _bsphere;

public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor for empty bbox
    TBoundingVolume ()
    {}
    
    //! ctor with given bbox and bsphere
    TBoundingVolume ( const TBBox &     bbox,
                      const TBSphere &  bsphere )
            : _bbox( bbox )
            , _bsphere( bsphere )
    {
        if ( bbox.dim() != bsphere.dim() )
            HERROR( ERR_ARG, "(TBoundingVolume) ctor", "different spatial dimension in bbox/bsphere" );
    }

    //! ctor for given bbox (bsphere based on bbox)
    TBoundingVolume ( const TBBox &     bbox )
            : _bbox( bbox )
    {
        _bsphere = TBSphere( ( _bbox.max() + _bbox.min() ) * 0.5,
                             _bbox.diameter() * 0.5 );
    }
    
    //! copy ctor
    TBoundingVolume ( const TBoundingVolume & box )
    {
        *this = box;
    }
    
    ///////////////////////////////////////////////
    //
    // access local variables
    //

    //! return bounding box
    const TBBox &      bbox     () const { return _bbox; }
    
    //! return bounding sphere
    const TBSphere &   bsphere  () const { return _bsphere; }
    
    //! return minimal coordinate of bbox
    const TPoint &     min      () const { return _bbox.min(); }

    //! return maximal coordinate of bbox
    const TPoint &     max      () const { return _bbox.max(); }

    //! return spatial dimension of bbox
    uint               dim      () const { return _bbox.dim(); }
    
    ///////////////////////////////////////////////
    //
    // bounding box properties
    //

    //! return true if point \a x is inside box
    bool    is_inside   ( const TPoint &  x ) const;
    
    //! return diameter of box
    double  diameter    () const;

    //! return volume of box
    double  volume      () const;

    //! return distance to \a box
    double  distance    ( const TBoundingVolume &  bvol ) const;

    //! return distance to \a box but coordinates have
    //! periodicity defined by \a period
    double  distance    ( const TBoundingVolume &  bvol,
                          const TPoint &           period ) const;

    //! return dimension of intersection with \a bbox
    uint    overlap_dim ( const TBoundingVolume &  bvol ) const;

    //! extend local bbox by given bounding volume
    void    extend      ( const TBoundingVolume &  bvol );
    
    //! extend local bbox by given bbox
    void    extend      ( const TBBox &  bbox );
    
    //! check volume and adjust if degenerate
    void    check       ();
    
    ///////////////////////////////////////////////
    //
    // misc.
    //

    //! return size in bytes used by this object
    size_t  byte_size () const { return _bbox.byte_size() + _bsphere.byte_size(); }

    //! copy operator
    TBoundingVolume & operator = ( const TBoundingVolume &  bvol )
    {
        _bbox    = bvol._bbox;
        _bsphere = bvol._bsphere;
        
        return *this;
    }

    //! return string representation
    std::string  to_string () const;
};

}// namespace Hpro

#endif  // __HPRO_TBOUNDINGVOLUME_HH
