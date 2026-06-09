//
// Project     : HLIBpro
// File        : T2DPrinter.cc
// Description : baseclass for all 2d printers
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/error.hh"

#include "T2DPrinter.hh"

namespace Hpro
{

//////////////////////////////////////////////////
//
// constructor and destructor
//

T2DPrinter::T2DPrinter ()
{
    _old_pos[0] = _old_pos[1] = 0.0;
    _active = false;
}

T2DPrinter::~T2DPrinter ()
{
}

//////////////////////////////////////////////////
//
// methods for drawing
//

// begin and end drawing
void
T2DPrinter::begin ()
{
    if ( _active )
        HERROR( ERR_CONSISTENCY, "(T2DPrinter) begin", "printer is still active" );
}

void
T2DPrinter::end ()
{
    if ( ! _active )
        HERROR( ERR_INIT, "(T2DPrinter) end", "" );
}

//
// transformation methods
//

//
// scale output
//
void
T2DPrinter::scale ( const double, const double )
{
    HERROR( ERR_NOT_IMPL, "(T2DPrinter) scale", "" );
}

//
// translate output
//
void
T2DPrinter::translate ( const double, const double )
{
    HERROR( ERR_NOT_IMPL, "(T2DPrinter) translate", "" );
}

//
// rotate output
//
void
T2DPrinter::rotate ( const double )
{
    HERROR( ERR_NOT_IMPL, "(T2DPrinter) rotate", "" );
}

//
// save and restore current state
//
void
T2DPrinter::save ()
{
    HERROR( ERR_NOT_IMPL, "(T2DPrinter) save", "" );
}

void
T2DPrinter::restore ()
{
    HERROR( ERR_NOT_IMPL, "(T2DPrinter) restore", "" );
}
    
//
// 2D - drawing
//

void
T2DPrinter::move_to ( const double x, const double y )
{
    _old_pos[0] = x;
    _old_pos[1] = y;
}

void
T2DPrinter::line_to ( const double x, const double y )
{
    draw_line( _old_pos[0], _old_pos[1], x, y );
    _old_pos[0] = x;
    _old_pos[1] = y;
}

//
// set color
//
void
T2DPrinter::set_hsv ( const int h, const int s, const int v )
{
    int r, g, b;

    hsv_to_rgb( h, s, v, r, g, b );
    return set_rgb( r, g, b );
}

//
// convertes HSV to RGB space
//
void
T2DPrinter::hsv_to_rgb ( const int H, const int s, const int v,
                         int & r, int & g, int & b ) const
{
    r = 0;
    g = 0;
    b = 0;

    if (( s == 0 ) || ( H == -1 ))
    {                                    // achromatic case
        r = g = b = v;
    }// if
    else
    {                                    // chromatic case
        uint  h = H;
        uint  f, p, q, t;
        
        if ( uint( h ) >= 360 )
            h %= 360;
        
        f  = h % 60;
        h /= 60;
        p  = uint((2 * v * (255 - s) + 255) / 510);

        if ( h & 1 )
        {
            q = uint((2 * v * (15300 - s * f) + 15300) / 30600);
            
            switch( h )
            {
                case 1: r = int(q); g = int(v); b = int(p); break;
                case 3: r = int(p); g = int(q); b = int(v); break;
                case 5: r = int(v); g = int(p); b = int(q); break;
            }// switch
        }// if
        else
        {
            t = uint((2 * v * (15300 - (s * (60 - f))) + 15300) / 30600);

            switch( h )
            {
                case 0: r = int(v); g = int(t); b = int(p); break;
                case 2: r = int(p); g = int(v); b = int(t); break;
                case 4: r = int(t); g = int(p); b = int(v); break;
            }// switch
        }// else
    }// else
}

}// namespace
