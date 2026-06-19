#ifndef __HPRO_TBBOX_HH
#define __HPRO_TBBOX_HH
//
// Project     : HLIBpro
// File        : TBBox.hh
// Description : class for a axis aligned bounding box 
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <hpro/base/types.hh>
#include <hpro/base/TPoint.hh>

namespace Hpro
{

// forward decl. for ctor
class TBSphere;

/////////////////////////////////////////////////////////////////
//
// type for a bounding box
//
class TBBox
{
private:
    // minimal and maximal point of box
    TPoint  _bb_min, _bb_max;

public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor for empty bbox
    TBBox ()
    {}
    
    //! ctor for bbox [\a bbmin, \a bbmax ]
    TBBox ( const TPoint &  bbmin,
            const TPoint &  bbmax )
            : _bb_min( bbmin )
            , _bb_max( bbmax )
    {
        if ( bbmin.dim() != bbmax.dim() )
            HERROR( ERR_ARG, "(TBBox) ctor", "different spatial dimension in bbox coordinates" );
    }

    //! ctor for bbox [\a bbmin, \a bbmax ] (special version for T2Point)
    TBBox ( const T2Point &  bbmin,
            const T2Point &  bbmax )
            : _bb_min( 2, bbmin.vector() )
            , _bb_max( 2, bbmax.vector() )
    {}

    //! ctor for bbox [\a bbmin, \a bbmax ] (special version for T3Point)
    TBBox ( const T3Point &  bbmin,
            const T3Point &  bbmax )
            : _bb_min( 3, bbmin.vector() )
            , _bb_max( 3, bbmax.vector() )
    {}

    //! copy ctor
    TBBox ( const TBBox & box )
    {
        *this = box;
    }
    
    //! convert bounding sphere to bounding box
    TBBox ( const TBSphere &  bsphere );
    
    ///////////////////////////////////////////////
    //
    // access local variables
    //

    //! return minimal coordinate of bbox
    TPoint &       min ()       { return _bb_min; }
    const TPoint & min () const { return _bb_min; }

    //! return maximal coordinate of bbox
    TPoint &       max ()       { return _bb_max; }
    const TPoint & max () const { return _bb_max; }

    //! return spatial dimension of bbox
    uint           dim () const { return _bb_min.dim(); }
    
    //! set spatial dimension of bbox
    void           set_dim ( const uint  d )
    {
        _bb_min.set_dim( d );
        _bb_max.set_dim( d );
    }
    
    ///////////////////////////////////////////////
    //
    // bounding box properties
    //

    //! return true if point \a x is inside box
    bool    is_inside   ( const TPoint & x ) const;
    
    //! return diameter of box
    double  diameter    () const;

    //! return volume of box
    double  volume      () const;

    //! return distance to \a box
    double  distance    ( const TBBox & box ) const;

    //! return distance to \a box but coordinates have
    //! periodicity defined by \a period
    double  distance    ( const TBBox &   box,
                          const TPoint &  period ) const;

    //! join local bbox with \a box
    void    join        ( const TBBox & box );
    
    //! extend local bbox by given point
    void    extend      ( const TPoint &  p );

    //! extend local bbox by given bbox
    void    extend      ( const TBBox &  bbox ) { join( bbox ); }

    //! return dimension of intersection with \a bbox
    uint    overlap_dim ( const TBBox &  bbox ) const;
    
    //! check volume and adjust if degenerate
    void    check       ();
    
    ///////////////////////////////////////////////
    //
    // misc.
    //

    //! return size in bytes used by this object
    size_t byte_size () const
    { return _bb_min.byte_size() + _bb_max.byte_size(); }

    //! copy operator
    TBBox & operator = ( const TBBox & box );

    //! return string representation
    std::string  to_string () const;
};

//!
//! return intersection of bboxes
//!
inline
TBBox
intersection ( const TBBox &  bbox1,
               const TBBox &  bbox2 )
{
    if ( bbox1.dim() != bbox2.dim() )
        HERROR( ERR_ARG, "intersection", "bboxes have different dimension" );
    
    TBBox  inter( bbox1 );

    for ( uint  i = 0; i < bbox1.dim(); ++i )
    {
        const  auto  r1 = bbox1.min()[i];
        const  auto  r2 = bbox1.max()[i];

        const  auto  c1 = bbox2.min()[i];
        const  auto  c2 = bbox2.max()[i];

        if      ( c2 <  r1 ) { inter.min()[i] = inter.max()[i] = 0; }
        else if ( c2 <= r2 ) { inter.min()[i] = std::max( r1, c1 ); inter.max()[i] = c2; }
        else if ( c1 <= r2 ) { inter.min()[i] = c1;                 inter.max()[i] = r2; }
        else if ( c1 >  r2 ) { inter.min()[i] = inter.max()[i] = 0; }
        else                 { inter.min()[i] = inter.max()[i] = 0; } // not reachable
    }// for

    return inter;
}

}// namespace Hpro

#endif  // __HPRO_TBBOX_HH
