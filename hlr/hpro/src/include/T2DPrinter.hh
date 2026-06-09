#ifndef __HPRO_T2DPRINTER_HH
#define __HPRO_T2DPRINTER_HH
//
// Project     : HLIBpro
// File        : T2DPrinter.hh
// Description : baseclass for all 2d printers
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/types.hh"

#include "colour.hh"

namespace Hpro
{

//
// text justification
//
enum justification_t
{
    JUST_LEFT,
    JUST_RIGHT,
    JUST_CENTER
};

//
// basicly implements world <-> window transformation
// and some colour-space conversion
//
class T2DPrinter
{
protected:
    // is printer active
    bool    _active;
    
    // old pos (for lineto, moveto)
    double  _old_pos[2];

public:
    //////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    T2DPrinter ();

    virtual ~T2DPrinter ();

    //////////////////////////////////////////////////
    //
    // methods for drawing
    //

    // begin and end drawing
    virtual void begin () = 0;
    virtual void end   () = 0;

    //
    // transformation methods
    //
    
    // scale output by ( \a x, \a y )
    virtual void scale     ( const double x, const double y );
    
    // translate output by ( \a x, \a y )
    virtual void translate ( const double x, const double y );
    
    // rotate output by \a angle ccw
    virtual void rotate    ( const double angle );
    
    // save and restore current state
    virtual void save    ();
    virtual void restore ();
    
    //
    // 2D - drawing
    //

    virtual void draw_point ( const double x,  const double y ) = 0;
    virtual void draw_line  ( const double x1, const double y1,
                              const double x2, const double y2 ) = 0;

    virtual void draw_triangle ( const double x0, const double y0,
                                 const double x1, const double y1,
                                 const double x2, const double y2 ) = 0;
    virtual void fill_triangle ( const double x0, const double y0,
                                 const double x1, const double y1,
                                 const double x2, const double y2 ) = 0;
    
    virtual void draw_rect ( const double x1, const double y1,
                             const double x2, const double y2 ) = 0;
    
    virtual void fill_rect ( const double x1, const double y1,
                             const double x2, const double y2 ) = 0;
    
    virtual void draw_circle ( const double x, const double y, const double radius ) = 0;
    virtual void fill_circle ( const double x, const double y, const double radius ) = 0;
    
    virtual void move_to ( const double x, const double y );
    virtual void line_to ( const double x, const double y );

    virtual void draw_text( const double x, const double y,
                            const std::string & text,
                            const justification_t just    = JUST_LEFT,
                            const bool            outline = false ) = 0;

    //
    // set styles
    //

    virtual void set_gray ( const int g ) { return set_rgb( g, g, g ); }
    virtual void set_rgb  ( const int r, const int g, const int b ) = 0;
    virtual void set_hsv  ( const int h, const int s, const int v );

    virtual void set_colour     ( const colour_t       c )     = 0;
    virtual void set_line_width ( const double         width ) = 0;
    virtual void set_font       ( const std::string &  font,
                                  const double         size )  = 0;

protected:
    //
    // misc.
    //

    // convertes HSV to RGB space
    void hsv_to_rgb ( const int h, const int s, const int v,
                      int & r, int & g, int & b ) const;
};

}// namespace

#endif  // __HPRO_T2DPRINTER_HH
