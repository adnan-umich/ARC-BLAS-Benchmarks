#ifndef __HPRO_TPSPRINTER_HH
#define __HPRO_TPSPRINTER_HH
//
// Project     : HLIBpro
// File        : TPSPrinter.hh
// Description : class for a 2d PostScipt printer
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <ostream>
#include <vector>

#include "hpro/base/TPoint.hh"

#include "T2DPrinter.hh"

namespace Hpro
{

//
// write output in encapsulated-postcript into file
//
class TPSPrinter : public T2DPrinter
{
protected:
    // size of window to print to
    uint            _width, _height;
    
    // the output stream
    std::ostream *  _outstream;

    // if true, stream was created externally
    bool            _ext_stream;
    
public:
    //////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TPSPrinter ( const uint           width,
                 const uint           height,
                 const std::string &  name );

    TPSPrinter ( const uint      width,
                 const uint      height,
                 std::ostream &  outstream );

    virtual ~TPSPrinter ();

    //////////////////////////////////////////////////
    //
    // methods for drawing
    //

    // begin and end drawing
    virtual void begin ();
    virtual void end   ();

    //
    // transformation methods
    //
    
    // scale output
    virtual void scale     ( const double x, const double y );
    
    // translate output
    virtual void translate ( const double x, const double y );
    
    // rotate output by \a angle ccw
    virtual void rotate    ( const double angle );
    
    // save and restore current state
    virtual void save    ();
    virtual void restore ();
    
    //
    // 2D - drawing
    //

    virtual void draw_point ( const double x, const double y );
    virtual void draw_line  ( const double x1, const double y1,
                              const double x2, const double y2 );

    virtual void draw_triangle ( const double x0, const double y0,
                                 const double x1, const double y1,
                                 const double x2, const double y2 );
    virtual void fill_triangle ( const double x0, const double y0,
                                 const double x1, const double y1,
                                 const double x2, const double y2 );
    
    virtual void draw_rect ( const double x1, const double y1,
                             const double x2, const double y2 );
    
    virtual void fill_rect ( const double x1, const double y1,
                             const double x2, const double y2 );

    virtual void draw_poly ( const std::vector< T2Point > & points );
    virtual void fill_poly ( const std::vector< T2Point > & points );
    
    virtual void draw_circle ( const double x, const double y, const double radius );
    virtual void fill_circle ( const double x, const double y, const double radius );
    
    virtual void draw_text( const double x, const double y,
                            const std::string & text,
                            const justification_t just    = JUST_LEFT,
                            const bool            outline = false );

    //
    // set styles
    //

    virtual void set_gray ( const int g );
    virtual void set_rgb  ( const int r, const int g, const int b );
    virtual void set_hsv  ( const int h, const int s, const int v );

    virtual void set_colour     ( const colour_t       c );
    virtual void set_pen        ( const uint           type );
    virtual void set_line_width ( const double         width );
    virtual void set_font       ( const std::string &  font,
                                  const double         size );

protected:
    // actually update style, i.e. write changes to file
    void update_style ();

    // write a floating point number to file
    void write_float ( const double x );
};

}// namespace

#endif  // __HPRO_TPSPRINTER_HH
