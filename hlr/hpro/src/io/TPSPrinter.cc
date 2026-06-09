//
// Project     : HLIBpro
// File        : TPPrinter.cc
// Description : class for a 2d PostScipt printer
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <fstream>

#include "hpro/base/error.hh"
#include "hpro/base/System.hh"

#include "TPSPrinter.hh"

namespace Hpro
{

// different pen-types
static char pen_types[][64] = { { "" },               // solid line
                                { "2" },              // . . . . 
                                { "4" },              // --  --  --  --
                                { "8" },              // ---   ---   ---
                                { "4 3 2 3" },        // -- . -- . -- 
                                { "8 5 2 5" },        // ---  .  ---  .
                                { "4 3 2 3 2 3" },    // -- . . -- . . --
                                { "8 5 2 5 2 5" } };  // ---  .  .  ---  .

//////////////////////////////////////////////////
//
// constructor and destructor
//

TPSPrinter::TPSPrinter ( const uint           width,
                         const uint           height,
                         const std::string &  fname )
{
    _width  = width;
    _height = height;
    
    _old_pos[0] = _old_pos[1] = 0.0;
    
    _active = false;

    _outstream  = new std::ofstream( fname.c_str() );
    _ext_stream = false;
}

TPSPrinter::TPSPrinter ( const uint      width,
                         const uint      height,
                         std::ostream &  outstream )
{
    _width  = width;
    _height = height;
    
    _old_pos[0] = _old_pos[1] = 0.0;
    
    _active = false;

    _outstream  = & outstream;
    _ext_stream = true;

    if ( _outstream == nullptr )
        HERROR( ERR_ARG, "(TPSPrinter) ctor", "nullptr output stream" );
}

TPSPrinter::~TPSPrinter ()
{
    if ( _active )
        end();

    if ( ! _ext_stream && ( _outstream != nullptr ) )
        delete _outstream;
}

//////////////////////////////////////////////////
//
// methods for drawing
//

// begin and end drawing
void
TPSPrinter::begin ()
{
    if ( _active )
        HERROR( ERR_CONSISTENCY, "(TPSPrinter) begin", "printer is still active" );

    _active = true;

    //
    // write header of eps-file
    //

    const int border = 5;
    const int ofs    = 20;

    // print preamble
    (*_outstream) << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl
                  << "%Title: " << std::endl
                  << "%Creator: TPSPrinter" << std::endl
                  << "%Orientation: Portrait" << std::endl
                  << "%%BoundingBox: " << (-border+ofs) << " " << (-border+ofs) << " " << (_width+border+ofs) << " " << (_height+border+ofs) << std::endl
                  << "%Pages: 0" << std::endl
                  << "%BeginSetup" << std::endl
                  << "%EndSetup" << std::endl
                  << "%Magnification: 1.0000" << std::endl
                  << "%EndComments" << std::endl

                  << "/PSPrinterDict 200 dict def" << std::endl
                  << "PSPrinterDict begin" << std::endl
    
                  << "/BD   { bind def } bind def" << std::endl 
                  << "/X    { exch } bind def" << std::endl
                  << std::endl
    
                  << "/NP   { newpath } bind def" << std::endl 
                  << "/M    { moveto } bind def" << std::endl 
                  << "/L    { lineto } bind def" << std::endl 
                  << "/RL   { rlineto } bind def" << std::endl 
                  << "/DA   { NP arc closepath } bind def" << std::endl 
                  << "/CP   { closepath } bind def" << std::endl 
                  << "/S    { stroke } bind def" << std::endl 
                  << "/F    { fill } bind def" << std::endl
                  << "/R    { roll } bind def" << std::endl
                  << std::endl
    
                  << "/DL   { 4 2 roll" << std::endl
                  << "        newpath" << std::endl
                  << "        moveto lineto" << std::endl
                  << "        closepath" << std::endl
                  << "        stroke } bind def" << std::endl
                  << "/DR   { newpath" << std::endl
                  << "        4 2 roll" << std::endl
                  << "        moveto exch dup 0 rlineto" << std::endl
                  << "        exch 0 exch rlineto" << std::endl
                  << "        -1 mul 0 rlineto" << std::endl
                  << "        closepath stroke } bind def" << std::endl
                  << "/FR   { newpath" << std::endl
                  << "        4 2 roll" << std::endl
                  << "        moveto exch dup 0 rlineto" << std::endl
                  << "        exch 0 exch rlineto" << std::endl
                  << "        -1 mul 0 rlineto" << std::endl
                  << "        closepath fill } bind def" << std::endl
                  << "/DT   { newpath" << std::endl
                  << "        moveto lineto lineto" << std::endl
                  << "        closepath stroke } bind def" << std::endl
                  << "/FT   { newpath" << std::endl
                  << "        moveto lineto lineto" << std::endl
                  << "        closepath fill } bind def" << std::endl
                  << "/DC   { newpath" << std::endl
                  << "        0 360 arc" << std::endl
                  << "        closepath stroke } bind def" << std::endl
                  << "/FC   { newpath" << std::endl
                  << "        0 360 arc" << std::endl
                  << "        closepath fill } bind def" << std::endl
                  << std::endl
        
                  << "/SW   { setlinewidth } bind def" << std::endl 
                  << "/SD   { setdash } bind def" << std::endl 
                  << "/SRGB { setrgbcolor } bind def" << std::endl 
                  << "/SG   { setgray } bind def" << std::endl 
                  << "/FF   { findfont } bind def" << std::endl 
                  << "/SF   { setfont } bind def" << std::endl 
                  << "/SCF  { scalefont } bind def" << std::endl 
                  << "/TB   { 3 1 roll" << std::endl
                  << "        translate 1 -1 scale" << std::endl
                  << "      } bind def" << std::endl
                  << "/TL   { gsave" << std::endl
                  << "          TB newpath 0 0 moveto show" << std::endl
                  << "        grestore } bind def" << std::endl
                  << "/TR   { gsave" << std::endl
                  << "          TB" << std::endl
                  << "          dup stringwidth pop" << std::endl
                  << "          newpath neg 0 moveto show" << std::endl
                  << "        grestore } bind def" << std::endl
                  << "/TC   { gsave" << std::endl
                  << "          TB" << std::endl
                  << "          dup stringwidth pop 2 div" << std::endl
                  << "          newpath neg 0 moveto show" << std::endl
                  << "        grestore } bind def" << std::endl
                  << "/OTL  { gsave" << std::endl
                  << "          TB 0 0 moveto true charpath stroke" << std::endl
                  << "        grestore } bind def" << std::endl
                  << "/OTR  { gsave" << std::endl
                  << "          TB" << std::endl
                  << "          dup stringwidth pop" << std::endl
                  << "          neg 0 moveto true charpath stroke" << std::endl
                  << "        grestore } bind def" << std::endl
                  << "/OTC  { gsave" << std::endl
                  << "          TB" << std::endl
                  << "          dup stringwidth pop 2 div" << std::endl
                  << "          neg 0 moveto true charpath stroke" << std::endl
                  << "        grestore } bind def" << std::endl

                  << "/GS   { gsave } bind def" << std::endl 
                  << "/GR   { grestore } bind def" << std::endl 
                  << std::endl
        
                  << "end" << std::endl
                  << "%EndProlog" << std::endl
                  << "PSPrinterDict begin" << std::endl
                  << "/psprnsavedpage save def" << std::endl;

    //
    // setup additional details
    //
    
    save();
    set_line_width( 1 );
    set_gray( 0 );
    scale( 1, -1 );
    translate( 0, -int(_height) );
    translate( ofs, -int(ofs) );
}

void
TPSPrinter::end ()
{
    if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) end", "" );

    restore();
        
    (*_outstream) << "psprnsavedpage restore" << std::endl
                  << "end" << std::endl
                  << "showpage" << std::endl
                  << "%EOF" << std::endl
                  << std::flush;
    // _file.close();
    
    _active = false;
}

//
// transformation methods
//

//
// scale output
//
void
TPSPrinter::scale ( const double x, const double y )
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) scale", "" );

    (*_outstream) << x << " " << y << " scale" << std::endl;
}

//
// translate output
//
void
TPSPrinter::translate ( const double x, const double y )
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) translate", "" );
    
    write_float( x );
    write_float( y );
    (*_outstream) << "translate" << std::endl;
}

//
// rotate output by \a angle
//
void
TPSPrinter::rotate ( const double angle )
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) rotate", "" );
    
    write_float( angle );
    (*_outstream) << "rotate" << std::endl;
}

//
// save and restore current state
//
void
TPSPrinter::save ()
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) save", "" );

    (*_outstream) << "GS" << std::endl;
}

void
TPSPrinter::restore ()
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) restore", "" );

    (*_outstream) << "GR" << std::endl;
}
    
//
// 2D - drawing
//

void
TPSPrinter::draw_point ( const double x, const double y )
{
    // simulate point by small circle
    return fill_circle( x, y, 1.5 );
}

void
TPSPrinter::draw_line ( const double x1, const double y1,
                        const double x2, const double y2 )
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) draw_line", "" );

    write_float( x1 );
    write_float( y1 );
    write_float( x2 );
    write_float( y2 );
    (*_outstream) << "DL" << std::endl;
}

void
TPSPrinter::draw_triangle ( const double x0, const double y0,
                            const double x1, const double y1,
                            const double x2, const double y2 )
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) draw_triangle", "" );

    write_float( x0 );
    write_float( y0 );
    write_float( x1 );
    write_float( y1 );
    write_float( x2 );
    write_float( y2 );
    (*_outstream) << "DT" << std::endl;
}

void
TPSPrinter::fill_triangle ( const double x0, const double y0,
                            const double x1, const double y1,
                            const double x2, const double y2 )
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) fill_triangle", "" );

    write_float( x0 );
    write_float( y0 );
    write_float( x1 );
    write_float( y1 );
    write_float( x2 );
    write_float( y2 );
    (*_outstream) << "FT" << std::endl;
}
    
void
TPSPrinter::draw_rect ( const double x1, const double y1, const double x2, const double y2 )
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) draw_rect", "" );

    write_float( x1 );
    write_float( y1 );
    write_float( x2-x1 );
    write_float( y2-y1 );
    (*_outstream) << "DR" << std::endl;
}

void
TPSPrinter::fill_rect ( const double x1, const double y1, const double x2, const double y2 )
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) fill_rect", "" );

    write_float( x1 );
    write_float( y1 );
    write_float( x2-x1 );
    write_float( y2-y1 );
    (*_outstream) << "FR" << std::endl;
}

void
TPSPrinter::draw_poly ( const std::vector< T2Point > & points )
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) draw_poly", "" );

    if ( points.size() < 2 )
        return;
    
    (*_outstream) << "NP" << std::endl;
    write_float( points[0][0] ); write_float( points[0][1] );
    (*_outstream) << "M" << std::endl;

    for ( uint i = 1; i < points.size(); i++ )
    {
        write_float( points[i][0] ); write_float( points[i][1] );
        (*_outstream) << "L" << std::endl;
    }// for
    (*_outstream) << "CP S" << std::endl;
}

void
TPSPrinter::fill_poly ( const std::vector< T2Point > & points )
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) fill_poly", "" );

    if ( points.size() < 2 )
        return;
    
    (*_outstream) << "NP" << std::endl;
    write_float( points[0][0] ); write_float( points[0][1] );
    (*_outstream) << "M" << std::endl;

    for ( uint i = 1; i < points.size(); i++ )
    {
        write_float( points[i][0] ); write_float( points[i][1] );
        (*_outstream) << "L" << std::endl;
    }// for
    (*_outstream) << "CP F" << std::endl;
}

void
TPSPrinter::draw_text( const double           x,
                       const double           y,
                       const std::string &    text,
                       const justification_t  just,
                       const bool             outline )
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) draw_text", "" );

    write_float( x );
    write_float( y );

    (*_outstream) << "(" << text << ") ";
    
    if ( outline )
        (*_outstream) << "O";

    switch ( just )
    {
    case JUST_LEFT   : (*_outstream) << "TL" << std::endl; break;
    case JUST_RIGHT  : (*_outstream) << "TR" << std::endl; break;
    case JUST_CENTER : (*_outstream) << "TC" << std::endl; break;
    }// switch
}

void
TPSPrinter::draw_circle ( const double x, const double y, const double radius )
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) draw_circle", "" );

    write_float( x );
    write_float( y );
    write_float( radius );
    (*_outstream) << "DC" << std::endl;
}

void
TPSPrinter::fill_circle ( const double x, const double y, const double radius )
{
	if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) fill_circle", "" );

    write_float( x );
    write_float( y );
    write_float( radius );
    (*_outstream) << "FC" << std::endl;
}

//
// set color
//
void
TPSPrinter::set_gray ( const int g )
{
    if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) set_gray", "" );

    write_float( float( std::max( 0, std::min( 255, g ) ) ) / 256.0 );
    (*_outstream) << "SG" << std::endl;
}

void
TPSPrinter::set_rgb ( const int r, const int g, const int b )
{
    if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) set_rgb", "" );

    set_colour( rgb( r, g, b ) );
    // write_float( float( max( 0, min( 255, r ) ) ) / 256.0 );
    // write_float( float( max( 0, min( 255, g ) ) ) / 256.0 );
    // write_float( float( max( 0, min( 255, b ) ) ) / 256.0 );
    // (*_outstream) << "SRGB" << std::endl;
}

void
TPSPrinter::set_hsv ( const int h, const int s, const int v )
{
    if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) set_hsv", "" );


    set_colour( hsv( h, s, v ) );

    // int r, g, b;
    // hsv_to_rgb( h, s, v, r, g, b );
    // return set_rgb( r, g, b );
}
    
void
TPSPrinter::set_colour ( const colour_t  c )
{
    if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) set_colour", "" );

    write_float( float( std::min( 255u, c.red   ) ) / 255.0 );
    write_float( float( std::min( 255u, c.green ) ) / 255.0 );
    write_float( float( std::min( 255u, c.blue  ) ) / 255.0 );
    (*_outstream) << "SRGB" << std::endl;
}

void
TPSPrinter::set_pen ( const uint type )
{
    if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) set_pen", "" );

    if ( type >= sizeof(pen_types)/(64*sizeof(char)) )
        HERROR( ERR_ARG, "(TPSPrinter) set_pen", "unkown pen" );

    (*_outstream) << "[ " << pen_types[type] << " ] 0 SD" << std::endl;
}

void
TPSPrinter::set_line_width ( const double width )
{
    if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) set_line_width", "" );

    write_float( width );
    (*_outstream) << "SW" << std::endl;
}

void
TPSPrinter::set_font( const std::string & font, const double size )
{
    if ( ! _active )
        HERROR( ERR_INIT, "(TPSPrinter) set_font", "" );

    (*_outstream) << "/" + font + " FF ";
    write_float( std::max( 0.01, size ) );
    (*_outstream) << "SCF SF" << std::endl;
}

//
// write a floating point number to file
//
void
TPSPrinter::write_float ( const double f )
{
    // (*_outstream) << to_string( "%f ", f );
    (*_outstream) << f << ' ';
    
    // const double  acc = 1e-1;
    
    // if ( Math::abs( f ) >= 1.0 )
    // {
    //     if ( Math::abs( f - double(int(f)) ) < acc )
    //         (*_outstream) << to_string( "%d ", int(f) ) );
    //     else if ( Math::abs( 10.0 * f - double(int(10.0*f)) ) < acc )
    //         (*_outstream) << to_string( "%.1f ", f ) );
    //     else if ( Math::abs( 100.0 * f - double(int(100.0*f)) ) < acc )
    //         (*_outstream) << to_string( "%.2f ", f ) );
    //     else if ( Math::abs( 1000.0 * f - double(int(1000.0*f)) ) < acc )
    //         (*_outstream) << to_string( "%.3f ", f ) );
    //     else
    //         (*_outstream) << to_string( "%f ", f ) );
    // }// if
    // else
    // {
    //     const real mag = int( Math::log10( f ) ); 
    //     const real mf  = Math::pow( 10, mag ) + f;

    //     if ( Math::abs( mf - double(int(mf)) ) < acc )
    //         (*_outstream) << to_string( "%d ", int(f) ) );
    //     else if ( Math::abs( 10.0 * mf - double(int(10.0*mf)) ) < acc )
    //         (*_outstream) << to_string( "%.1f ", f ) );
    //     else if ( Math::abs( 100.0 * mf - double(int(100.0*mf)) ) < acc )
    //         (*_outstream) << to_string( "%.2f ", f ) );
    //     else if ( Math::abs( 1000.0 * mf - double(int(1000.0*mf)) ) < acc )
    //         (*_outstream) << to_string( "%.3f ", f ) );
    //     else
    //         (*_outstream) << to_string( "%f ", f ) );
        
    // }// else
}

}// namespace
