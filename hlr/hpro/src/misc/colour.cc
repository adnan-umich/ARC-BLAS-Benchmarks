//
// Project     : H-Matrix
// File        : colour.cc
// Description : colour managment functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "colour.hh"

namespace Hpro
{

//
// construct colour out of HSV values
//
colour_t
hsv ( const int  H,
      const int  S,
      const int  V )
{
    double  r = 0;
    double  g = 0;
    double  b = 0;
        
    if (( S == 0 ) || ( H == -1 ))
    {                                    // achromatic case
        r = g = b = double(V) / 100.0;
    }// if
    else
    {                                    // chromatic case
        double  h = double(H);
        double  s = double(S) / 100.0;
        double  v = double(V) / 100.0;
        
        int     hi = int( std::floor( h / 60.0 ) );
        double  f  = ( h / 60.0 ) - double( hi );

        double  p = v * ( 1.0 - s );
        double  q = v * ( 1.0 - s * f );
        double  t = v * ( 1.0 - s * ( 1.0 - f ) );

        switch( hi )
        {
            case 0:
            case 6: r = v; g = t; b = p; break;
            case 1: r = q; g = v; b = p; break;
            case 2: r = p; g = v; b = t; break;
            case 3: r = p; g = q; b = v; break;
            case 4: r = t; g = p; b = v; break;
            case 5: r = v; g = p; b = q; break;
        }// switch
    }// else

    return colour_t( uint(r * 255.0), uint(g * 255.0), uint(b * 255.0) );
}

}// namespace Hpro
