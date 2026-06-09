#ifndef __HPRO_COLOUR_HH
#define __HPRO_COLOUR_HH
//
// Project     : HLIBpro
// File        : colour.hh
// Description : colour managment functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/types.hh"
#include "hpro/base/String.hh"

namespace Hpro
{

//!
//! \struct  colour_t
//! \brief   RGB colour class
//!
struct colour_t
{
    //! red colour channel value
    uint  red;

    //! green colour channel value
    uint  green;

    //! blue colour channel value
    uint  blue;

    
    //! construct black colour
    colour_t ()
            : red(0), green(0), blue(0)
    {}

    //! construct colour with given colour value
    colour_t ( const uint r, const uint g, const uint b )
            : red(   std::min( 255u, r ) ),
              green( std::min( 255u, g ) ),
              blue(  std::min( 255u, b ) )
    {}

    //! construct colour from CSS code, e.g. 0xrrggbb (hex notation)
    colour_t ( const uint  css_code )
            : red(0), green(0), blue(0)
    {
        red   = ( css_code & 0xff0000 ) >> 16;
        green = ( css_code & 0x00ff00 ) >> 8;
        blue  = ( css_code & 0x0000ff );
    }

    //! mix given colour in as: scale · this + (1 - scale) · col
    void mix ( const double    scale,
               const colour_t  col )
    {
        const double  rest = 1.0 - scale;
        
        red   = std::min( 255u, uint( scale * double(red)   + rest * double(col.red)   ) );
        green = std::min( 255u, uint( scale * double(green) + rest * double(col.green) ) );
        blue  = std::min( 255u, uint( scale * double(blue)  + rest * double(col.blue)  ) );
    }

    //! return luminance ∈ [0,255]
    uint luminance () const
    {
        return uint( 255.0 * ( 0.299 * ( double(red)   / 255.0 ) +
                               0.587 * ( double(green) / 255.0 ) +
                               0.114 * ( double(blue)  / 255.0 ) ) );
    }

    //! return html color code
    std::string  to_string () const
    {
        return Hpro::to_string( "%02x%02x%02x", red, green, blue );
    }
};

//!
//! \brief          construct colour out of RGB values
//! \param   r      red value,   r ∈ [0,255]
//! \param   g      green value, g ∈ [0,255]
//! \param   b      blue value,  b ∈ [0,255]
//! \return         colour
//!
inline
colour_t
rgb ( const uint  r,
      const uint  g,
      const uint  b )
{
    return colour_t( r, g, b );
}

//!
//! \brief          construct grayscale colour
//! \param   g      gray value, g ∈ [0,255]
//! \return         colour
//!
inline
colour_t
gray ( const uint  g )
{
    return colour_t( g, g, g );
}

//!
//! \brief          construct colour out of HSV values
//! \param   H      hue,        H ∈ [0,360]
//! \param   s      saturation, s ∈ [0,100]
//! \param   v      value,      v ∈ [0,100]
//! \return         colour
//!
colour_t
hsv ( const int  H,
      const int  s,
      const int  v );

//!
//! \brief          functional version of colour_t::mix (\sa colour_t::mix)
//! \param   col1   first colour to mix in
//! \param   scale  scaling factor for first colour
//! \param   col2   second colour to mix in
//! \return         mixed colour defined by scale · col1 + (1 - scale) · col2
//!
inline
colour_t
mix ( const colour_t  col1,
      const double    scale,
      const colour_t  col2 )
{
    colour_t  res = col1;

    res.mix( scale, col2 );

    return res;
}

namespace Tango
{

//
// standard colours from the Tango project
//
const colour_t  White( 255, 255, 255 );
const colour_t  Black(   0,   0,   0 );

const colour_t  Butter1( 252,233,79 );
const colour_t  Butter2( 237,212,0 );
const colour_t  Butter3( 196,160,0 );

const colour_t  Chameleon1( 138,226,52 );
const colour_t  Chameleon2( 115,210,22 );
const colour_t  Chameleon3( 78,154,6 );

const colour_t  Orange1( 252,175,62 );
const colour_t  Orange2( 245,121,0 );
const colour_t  Orange3( 206,92,0 );

const colour_t  SkyBlue1( 114,159,207 );
const colour_t  SkyBlue2( 52,101,164 );
const colour_t  SkyBlue3( 32,74,135 );

const colour_t  Plum1( 173,127,168 );
const colour_t  Plum2( 117,80,123 );
const colour_t  Plum3( 92,53,102 );

const colour_t  Chocolate1( 233,185,110 );
const colour_t  Chocolate2( 193,125,17 );
const colour_t  Chocolate3( 143,89,2 );

const colour_t  ScarletRed1( 239,41,41 );
const colour_t  ScarletRed2( 204,0,0 );
const colour_t  ScarletRed3( 164,0,0 );

const colour_t  Aluminium1( 238,238,236 );
const colour_t  Aluminium2( 211,215,207 );
const colour_t  Aluminium3( 186,189,182 );
const colour_t  Aluminium4( 136,138,133 );
const colour_t  Aluminium5( 85,87,83 );
const colour_t  Aluminium6( 46,52,54 );

}// namespace Tango

}// namespace Hpro

#endif  // __HPRO_COLOUR_HH
