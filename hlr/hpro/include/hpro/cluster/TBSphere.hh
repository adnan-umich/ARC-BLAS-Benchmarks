#ifndef __HPRO_TBSPHERE_HH
#define __HPRO_TBSPHERE_HH
//
// Project     : HLIBpro
// File        : TBSphere.hh
// Description : class for a axis aligned bounding box 
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2024. All Rights Reserved.
//

#include <hpro/base/types.hh>
#include <hpro/base/TPoint.hh>

namespace Hpro
{

// forward decl. for ctor
class TBBox;

/////////////////////////////////////////////////////////////////
//
// bounding sphere represented by center point and radius
//
class TBSphere
{
private:
    // center and radius
    TPoint  _center;
    double  _radius;

public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor for empty bsphere
    TBSphere ()
            : _radius(0)
    {}
    
    //! ctor for defining center and radius
    TBSphere ( const TPoint &  acenter,
               const double    aradius )
            : _center( acenter )
            , _radius( aradius )
    {
        if ( _radius < 0 )
            HERROR( ERR_ARG, "(TBSphere) ctor", "negative radius" );
    }

    //! ctor same as above but for 2D
    TBSphere ( const T2Point &  acenter,
               const double     aradius )
            : _center( 2, acenter.vector() )
            , _radius( aradius )
    {
        if ( _radius < 0 )
            HERROR( ERR_ARG, "(TBSphere) ctor", "negative radius" );
    }

    //! ctor same as above but for 3D
    TBSphere ( const T3Point &  acenter,
               const double     aradius )
            : _center( 3, acenter.vector() )
            , _radius( aradius )
    {
        if ( _radius < 0 )
            HERROR( ERR_ARG, "(TBSphere) ctor", "negative radius" );
    }

    //! copy ctor
    TBSphere ( const TBSphere &  bsphere )
    {
        *this = bsphere;
    }

    //! convert bounding box to bounding sphere
    TBSphere ( const TBBox &  bbox );
    
    ///////////////////////////////////////////////
    //
    // access local variables
    //

    //! return center of bsphere
    const TPoint &  center () const { return _center; }

    //! return radius
    double          radius () const { return _radius; }

    //! return spatial dimension of bsphere
    uint            dim    () const { return _center.dim(); }
    
    ///////////////////////////////////////////////
    //
    // bounding bsphere properties
    //

    //! return true if point \a x is inside bsphere
    bool    is_inside ( const TPoint & x ) const;
    
    //! return diameter of bsphere
    double  diameter  () const { return 2 * _radius; }

    //! return volume of bsphere
    double  volume    () const;

    //! return distance to \a bsphere
    double  distance  ( const TBSphere &  bsphere ) const;
    
    //! return distance to \a bsphere (not supported)
    double  distance  ( const TBSphere &  bsphere,
                        const TPoint &    /* period */ ) const { return distance( bsphere ); }
    
    //! extend local bbox by given point
    void    extend      ( const TPoint &  p );

    //! extend local bbox by given bbox
    void    extend      ( const TBSphere &  bsphere );

    //! return dimension of intersection with \a bsphere
    uint    overlap_dim ( const TBSphere &  bsphere ) const;

    //! check volume and adjust if degenerate
    void    check       ();
    
    ///////////////////////////////////////////////
    //
    // misc.
    //

    //! return size in bytes used by this object
    size_t byte_size () const { return _center.byte_size() + sizeof(_radius); }

    //! copy operator
    TBSphere & operator = ( const TBSphere &  bsphere );

    //! return string representation
    std::string  to_string () const;
};

}// namespace Hpro

#endif  // __HPRO_TBSPHERE_HH
