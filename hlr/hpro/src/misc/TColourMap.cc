//
// Project     : HLIBpro
// File        : TColourMap.cc
// Description : provides various colourmaps
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <algorithm>
#include <vector>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include "hpro/base/error.hh"
#include "hpro/base/System.hh"

#include "TColourMap.hh"

namespace Hpro
{

namespace
{

// maximal colour channel value
const double MAX_COL_VAL = 255.0;

}// namespace anonymous

//
// return colourmap by name, e.g. jet, hot, coolwarm, etc.
//
TColourMap *
colourmap ( const std::string &  name,
            const size_t         n )
{
    std::string  tname = name;

    boost::to_lower( tname );
    boost::trim( tname );

    // adjust name for default colourmap
    if (( tname == "" ) || ( tname == "default" ))
        tname = "coolwarm";
    
    if      ( tname == "jet"      ) return new TJetColourMap( n );
    else if ( tname == "hot"      ) return new THotColourMap( n );
    else if ( tname == "copper"   ) return new TCopperColourMap( n );
    else if ( tname == "hsv"      ) return new THSVColourMap( n );
    else if ( tname == "bone"     ) return new TBoneColourMap( n );
    else if ( tname == "rainbow"  ) return new TRainbowColourMap( n );
    else if ( tname == "coolwarm" ) return new TCoolWarmColourMap( n );
    else if ( tname == "magma"    ) return new TMagmaColourMap( n );
    else if ( tname == "inferno"  ) return new TInfernoColourMap( n );
    else if ( tname == "plasma"   ) return new TPlasmaColourMap( n );
    else if ( tname == "viridis"  ) return new TViridisColourMap( n );
    else if ( tname == "parula"   ) return new TParulaColourMap( n );
    else
    {
        HWARNING( "unknown colourmap \"" + name + "\"; using \"default\"" );
        return colourmap( "default" );
    }// else

    return nullptr;
}

////////////////////////////////////////////////////////////////
//
// Remark:
//
//  Without different remarks, all colourmaps are from Octave 
//  written by Kai Habel <kai.habel@gmx.de>
//
////////////////////////////////////////////////////////////////

//
// TangoColourMap
//
void
TTangoColourMap::init ( const size_t n )
{
    //
    // base colors
    //

    static colour_t  tango[] = { 
        colour_t(0xef2929),
        colour_t(0x729fcf),
        colour_t(0x8ae234),
        colour_t(0xad7fa8),
        colour_t(0xfcaf3e),
        colour_t(0xfce94f),
        colour_t(0xe9b96e),
        colour_t(0xbabdb6),

        colour_t(0xcc0000),
        colour_t(0x3465a4),
        colour_t(0x73d216),
        colour_t(0x75507b),
        colour_t(0xf57900),
        colour_t(0xedd400),
        colour_t(0xc17d11),
        colour_t(0x888a85),

        colour_t(0xa40000),
        colour_t(0x204a87),
        colour_t(0x4e9a06),
        colour_t(0x5c3566),
        colour_t(0xce5c00),
        colour_t(0xc4a000),
        colour_t(0x8f5902),
        colour_t(0x555753) };
    
    const size_t  ntango = sizeof(tango) / sizeof(colour_t);
        
    //
    // just copy for requested number of colors
    //
    
    _map.resize( n );

    for ( size_t  i = 0; i < n; ++i )
        _map[i] = tango[ i % ntango ];
}

//
// JetColourMap
//
void
TJetColourMap::init ( const size_t n )
{
    _map.resize( n );
    
    for ( size_t i = 0; i < n; i++ )
    {
        const double  val = double(i) / double(n);
        double        r   = 0.0;
        double        g   = 0.0;
        double        b   = 0.0;

        if      ( val >= 3.0/8.0 && val < 5.0/8.0 ) r = 4.0 * val - 3.0/2.0;
        else if ( val >= 5.0/8.0 && val < 7.0/8.0 ) r = 1.0;
        else if ( val >= 7.0/8.0 )                  r = -4.0 * val + 9.0/2.0;
        
        if      ( val >= 1.0/8.0 && val < 3.0/8.0 ) g = 4.0 * val - 1.0/2.0;
        else if ( val >= 3.0/8.0 && val < 5.0/8.0 ) g = 1.0;
        else if ( val >= 5.0/8.0 && val < 7.0/8.0 ) g = -4.0 * val + 7.0/2.0;
        
        if      ( val <  1.0/8.0 )                  b = 4.0 * val + 1.0/2.0;
        else if ( val >= 1.0/8.0 && val < 3.0/8.0 ) b = 1.0;
        else if ( val >= 3.0/8.0 && val < 5.0/8.0 ) b = -4.0 * val + 5.0/2.0;

        _map[i] = rgb( uint(MAX_COL_VAL * r), uint(MAX_COL_VAL * g), uint(MAX_COL_VAL * b) );
    }// for
}

//
// HotColourMap
//
void
THotColourMap::init ( const size_t n )
{
    _map.resize( n );
    
    for ( size_t i = 0; i < n; i++ )
    {
        const double  val = double(i) / double(n);
        double        r   = 0.0;
        double        g   = 0.0;
        double        b   = 0.0;

        if ( val < 2.0/5.0 ) r = 5.0/2.0 * val;
        else                 r = 1.0;
        
        if      ( val >= 2.0/5.0 && val < 4.0/5.0 ) g = 5.0/2.0 * val - 1.0;
        else if ( val >= 4.0/5.0 )                  g = 1.0;
        
        if ( val >= 4.0/5.0 ) b = 5.0 * val - 4.0;

        _map[i] = rgb( uint(MAX_COL_VAL * r), uint(MAX_COL_VAL * g), uint(MAX_COL_VAL * b) );
    }// for
}

//
// CopperColourMap
//
void
TCopperColourMap::init ( const size_t n )
{
    _map.resize( n );
    
    for ( size_t i = 0; i < n; i++ )
    {
        const double  val = double(i) / double(n);
        double        r   = 0.0;
        double        g   = 0.0;
        double        b   = 0.0;

        if      ( val <  4.0/5.0 ) r = 5.0/4.0 * val;
        else if ( val >= 4.0/5.0 ) r = 1.0;

        g = 4.0/5.0 * val;
        b = 1.0/2.0 * val;

        _map[i] = rgb( uint(MAX_COL_VAL * r), uint(MAX_COL_VAL * g), uint(MAX_COL_VAL * b) );
    }// for
}

//
// HSVColourMap
//
void
THSVColourMap::init ( const size_t n )
{
    _map.resize( n );
    
    for ( size_t i = 0; i < n; i++ )
        _map[i] = hsv( int( 360.0 * double(i) / double(n) ),
                       100,
                       100 );
}

//
// BoneColourMap
//
void
TBoneColourMap::init ( const size_t n )
{
    _map.resize( n );
    
    for ( size_t i = 0; i < n; i++ )
    {
        const double  val = double(i) / double(n);
        double        r   = 0.0;
        double        g   = 0.0;
        double        b   = 0.0;

        if      ( val <  3.0/4.0 ) r = 7.0/8.0 * val;
        else if ( val >= 3.0/4.0 ) r = 11.0/8.0 * val - 3.0/8.0;
        
        if      ( val <  3.0/8.0                  ) g = 7.0/8.0 * val;
        else if ( val >= 3.0/8.0 && val < 3.0/4.0 ) g = 29.0/24.0 * val - 1.0/8.0;
        else if ( val >= 3.0/4.0                  ) g = 7.0/8.0 * val + 1.0/8.0;
        
        if      ( val <  3.0/8.0 ) b = 29.0/24.0 * val;
        else if ( val >= 3.0/8.0 ) b = 7.0/8.0 * val + 1.0/8.0;

        _map[i] = rgb( uint(MAX_COL_VAL * r), uint(MAX_COL_VAL * g), uint(MAX_COL_VAL * b) );
    }// for
}

//
// RainbowColourMap
//
void
TRainbowColourMap::init ( const size_t n )
{
    _map.resize( n );
    
    for ( size_t i = 0; i < n; i++ )
    {
        const double  val = double(i) / double(n);
        double        r   = 0.0;
        double        g   = 0.0;
        double        b   = 0.0;

        if      ( val <  2.0/5.0 )                  r = 1.0;
        else if ( val >= 2.0/5.0 && val < 3.0/5.0 ) r = -5.0 * val + 3.0;
        else if ( val >= 4.0/5.0 )                  r = 10.0/3.0 * val - 8.0/3.0;
        
        if      ( val <  2.0/5.0 )                  g = 5.0/2.0 * val;
        else if ( val >= 2.0/5.0 && val < 3.0/5.0 ) g = 1.0;
        else if ( val >= 3.0/5.0 && val < 4.0/5.0 ) g = -5.0 * val + 4.0;
        
        if      ( val >= 3.0/5.0 && val < 4.0/5.0 ) b = 5.0 * val - 3.0;
        else if ( val >= 4.0/5.0 )                  b = 1.0;

        _map[i] = rgb( uint(MAX_COL_VAL * r), uint(MAX_COL_VAL * g), uint(MAX_COL_VAL * b) );
    }// for
}

//
// array based colour maps function
//
void
init_cmap_with_data ( const size_t               n_cmap,
                      std::vector< colour_t > &  cmap,
                      const colour_t *           data,
                      const size_t               n_data )
{
    cmap.resize( n_cmap );

    if ( n_cmap <= n_data )
    {
        //
        // directly choose colour from interval [0,n_data)
        //
        
        for ( size_t  i = 0; i < n_cmap; ++i )
        {
            const size_t  pos = size_t( std::floor( double(i) / double(n_cmap) * (n_data-1) ) );

            cmap[i] = data[pos];
        }// for
    }// if
    else
    {
        //
        // interpolate between colours in interval [0,n_data)
        //
        
        for ( size_t  i = 0; i < n_cmap; ++i )
        {
            const size_t    pos1 = size_t( std::floor( double(i)   / double(n_cmap) * (n_data-1) ) );
            const size_t    pos2 = size_t( std::floor( double(i+1) / double(n_cmap) * (n_data-1) ) );
            const colour_t  c1 = data[pos1];
            const colour_t  c2 = data[pos2];;
            const colour_t  c( (c1.red   + c2.red  ) / 2,
                               (c1.green + c2.green) / 2,
                               (c1.blue  + c2.blue ) / 2 );
                               
            cmap[i] = c;
        }// for
    }// else
}
    
//
// CoolWarmColourMap
//
namespace
{

colour_t  coolwarm257[] = {
    colour_t( 59,76,192 ),
    colour_t( 60,78,194 ),
    colour_t( 61,80,195 ),
    colour_t( 62,81,197 ),
    colour_t( 63,83,198 ),
    colour_t( 64,85,200 ),
    colour_t( 66,87,201 ),
    colour_t( 67,88,203 ),
    colour_t( 68,90,204 ),
    colour_t( 69,92,206 ),
    colour_t( 70,93,207 ),
    colour_t( 71,95,209 ),
    colour_t( 73,97,210 ),
    colour_t( 74,99,211 ),
    colour_t( 75,100,213 ),
    colour_t( 76,102,214 ),
    colour_t( 77,104,215 ),
    colour_t( 79,105,217 ),
    colour_t( 80,107,218 ),
    colour_t( 81,109,219 ),
    colour_t( 82,110,221 ),
    colour_t( 84,112,222 ),
    colour_t( 85,114,223 ),
    colour_t( 86,115,224 ),
    colour_t( 87,117,225 ),
    colour_t( 89,119,226 ),
    colour_t( 90,120,228 ),
    colour_t( 91,122,229 ),
    colour_t( 93,123,230 ),
    colour_t( 94,125,231 ),
    colour_t( 95,127,232 ),
    colour_t( 96,128,233 ),
    colour_t( 98,130,234 ),
    colour_t( 99,131,235 ),
    colour_t( 100,133,236 ),
    colour_t( 102,135,237 ),
    colour_t( 103,136,238 ),
    colour_t( 104,138,239 ),
    colour_t( 106,139,239 ),
    colour_t( 107,141,240 ),
    colour_t( 108,142,241 ),
    colour_t( 110,144,242 ),
    colour_t( 111,145,243 ),
    colour_t( 112,147,243 ),
    colour_t( 114,148,244 ),
    colour_t( 115,150,245 ),
    colour_t( 116,151,246 ),
    colour_t( 118,153,246 ),
    colour_t( 119,154,247 ),
    colour_t( 120,156,247 ),
    colour_t( 122,157,248 ),
    colour_t( 123,158,249 ),
    colour_t( 124,160,249 ),
    colour_t( 126,161,250 ),
    colour_t( 127,163,250 ),
    colour_t( 129,164,251 ),
    colour_t( 130,165,251 ),
    colour_t( 131,167,252 ),
    colour_t( 133,168,252 ),
    colour_t( 134,169,252 ),
    colour_t( 135,171,253 ),
    colour_t( 137,172,253 ),
    colour_t( 138,173,253 ),
    colour_t( 140,174,254 ),
    colour_t( 141,176,254 ),
    colour_t( 142,177,254 ),
    colour_t( 144,178,254 ),
    colour_t( 145,179,254 ),
    colour_t( 147,181,255 ),
    colour_t( 148,182,255 ),
    colour_t( 149,183,255 ),
    colour_t( 151,184,255 ),
    colour_t( 152,185,255 ),
    colour_t( 153,186,255 ),
    colour_t( 155,187,255 ),
    colour_t( 156,188,255 ),
    colour_t( 158,190,255 ),
    colour_t( 159,191,255 ),
    colour_t( 160,192,255 ),
    colour_t( 162,193,255 ),
    colour_t( 163,194,255 ),
    colour_t( 164,195,254 ),
    colour_t( 166,196,254 ),
    colour_t( 167,197,254 ),
    colour_t( 168,198,254 ),
    colour_t( 170,199,253 ),
    colour_t( 171,199,253 ),
    colour_t( 172,200,253 ),
    colour_t( 174,201,253 ),
    colour_t( 175,202,252 ),
    colour_t( 176,203,252 ),
    colour_t( 178,204,251 ),
    colour_t( 179,205,251 ),
    colour_t( 180,205,251 ),
    colour_t( 182,206,250 ),
    colour_t( 183,207,250 ),
    colour_t( 184,208,249 ),
    colour_t( 185,208,248 ),
    colour_t( 187,209,248 ),
    colour_t( 188,210,247 ),
    colour_t( 189,210,247 ),
    colour_t( 190,211,246 ),
    colour_t( 192,212,245 ),
    colour_t( 193,212,245 ),
    colour_t( 194,213,244 ),
    colour_t( 195,213,243 ),
    colour_t( 197,214,243 ),
    colour_t( 198,214,242 ),
    colour_t( 199,215,241 ),
    colour_t( 200,215,240 ),
    colour_t( 201,216,239 ),
    colour_t( 203,216,238 ),
    colour_t( 204,217,238 ),
    colour_t( 205,217,237 ),
    colour_t( 206,217,236 ),
    colour_t( 207,218,235 ),
    colour_t( 208,218,234 ),
    colour_t( 209,219,233 ),
    colour_t( 210,219,232 ),
    colour_t( 211,219,231 ),
    colour_t( 213,219,230 ),
    colour_t( 214,220,229 ),
    colour_t( 215,220,228 ),
    colour_t( 216,220,227 ),
    colour_t( 217,220,225 ),
    colour_t( 218,220,224 ),
    colour_t( 219,220,223 ),
    colour_t( 220,221,222 ),
    colour_t( 221,221,221 ),
    colour_t( 222,220,219 ),
    colour_t( 223,220,218 ),
    colour_t( 224,219,216 ),
    colour_t( 225,219,215 ),
    colour_t( 226,218,214 ),
    colour_t( 227,218,212 ),
    colour_t( 228,217,211 ),
    colour_t( 229,216,209 ),
    colour_t( 230,216,208 ),
    colour_t( 231,215,206 ),
    colour_t( 232,215,205 ),
    colour_t( 232,214,203 ),
    colour_t( 233,213,202 ),
    colour_t( 234,212,200 ),
    colour_t( 235,212,199 ),
    colour_t( 236,211,197 ),
    colour_t( 236,210,196 ),
    colour_t( 237,209,194 ),
    colour_t( 238,209,193 ),
    colour_t( 238,208,191 ),
    colour_t( 239,207,190 ),
    colour_t( 240,206,188 ),
    colour_t( 240,205,187 ),
    colour_t( 241,204,185 ),
    colour_t( 241,203,184 ),
    colour_t( 242,202,182 ),
    colour_t( 242,201,181 ),
    colour_t( 243,200,179 ),
    colour_t( 243,199,178 ),
    colour_t( 244,198,176 ),
    colour_t( 244,197,174 ),
    colour_t( 245,196,173 ),
    colour_t( 245,195,171 ),
    colour_t( 245,194,170 ),
    colour_t( 245,193,168 ),
    colour_t( 246,192,167 ),
    colour_t( 246,191,165 ),
    colour_t( 246,190,163 ),
    colour_t( 246,188,162 ),
    colour_t( 247,187,160 ),
    colour_t( 247,186,159 ),
    colour_t( 247,185,157 ),
    colour_t( 247,184,156 ),
    colour_t( 247,182,154 ),
    colour_t( 247,181,152 ),
    colour_t( 247,180,151 ),
    colour_t( 247,178,149 ),
    colour_t( 247,177,148 ),
    colour_t( 247,176,146 ),
    colour_t( 247,174,145 ),
    colour_t( 247,173,143 ),
    colour_t( 247,172,141 ),
    colour_t( 247,170,140 ),
    colour_t( 247,169,138 ),
    colour_t( 247,167,137 ),
    colour_t( 247,166,135 ),
    colour_t( 246,164,134 ),
    colour_t( 246,163,132 ),
    colour_t( 246,161,131 ),
    colour_t( 246,160,129 ),
    colour_t( 245,158,127 ),
    colour_t( 245,157,126 ),
    colour_t( 245,155,124 ),
    colour_t( 244,154,123 ),
    colour_t( 244,152,121 ),
    colour_t( 244,151,120 ),
    colour_t( 243,149,118 ),
    colour_t( 243,147,117 ),
    colour_t( 242,146,115 ),
    colour_t( 242,144,114 ),
    colour_t( 241,142,112 ),
    colour_t( 241,141,111 ),
    colour_t( 240,139,109 ),
    colour_t( 240,137,108 ),
    colour_t( 239,136,106 ),
    colour_t( 238,134,105 ),
    colour_t( 238,132,103 ),
    colour_t( 237,130,102 ),
    colour_t( 236,129,100 ),
    colour_t( 236,127,99 ),
    colour_t( 235,125,97 ),
    colour_t( 234,123,96 ),
    colour_t( 233,121,95 ),
    colour_t( 233,120,93 ),
    colour_t( 232,118,92 ),
    colour_t( 231,116,90 ),
    colour_t( 230,114,89 ),
    colour_t( 229,112,88 ),
    colour_t( 228,110,86 ),
    colour_t( 227,108,85 ),
    colour_t( 227,106,83 ),
    colour_t( 226,104,82 ),
    colour_t( 225,102,81 ),
    colour_t( 224,100,79 ),
    colour_t( 223,98,78 ),
    colour_t( 222,96,77 ),
    colour_t( 221,94,75 ),
    colour_t( 220,92,74 ),
    colour_t( 218,90,73 ),
    colour_t( 217,88,71 ),
    colour_t( 216,86,70 ),
    colour_t( 215,84,69 ),
    colour_t( 214,82,67 ),
    colour_t( 213,80,66 ),
    colour_t( 212,78,65 ),
    colour_t( 210,75,64 ),
    colour_t( 209,73,62 ),
    colour_t( 208,71,61 ),
    colour_t( 207,69,60 ),
    colour_t( 205,66,59 ),
    colour_t( 204,64,57 ),
    colour_t( 203,62,56 ),
    colour_t( 202,59,55 ),
    colour_t( 200,57,54 ),
    colour_t( 199,54,53 ),
    colour_t( 198,51,52 ),
    colour_t( 196,49,50 ),
    colour_t( 195,46,49 ),
    colour_t( 193,43,48 ),
    colour_t( 192,40,47 ),
    colour_t( 190,37,46 ),
    colour_t( 189,34,45 ),
    colour_t( 188,30,44 ),
    colour_t( 186,26,43 ),
    colour_t( 185,22,41 ),
    colour_t( 183,17,40 ),
    colour_t( 181,11,39 ),
    colour_t( 180,4,38 )
};

}// namespace anonymous

void
TCoolWarmColourMap::init ( const size_t n )
{
    init_cmap_with_data( n, _map, coolwarm257, sizeof(coolwarm257) / sizeof(colour_t) );
                         
    // _map.resize( n );

    // if ( n <= 257 )
    // {
    //     //
    //     // directly choose colour from interval [0,256]
    //     //
        
    //     for ( size_t  i = 0; i < n; ++i )
    //     {
    //         const size_t  pos = size_t( std::floor( double(i)/double(n) * 256.0 ) );

    //         _map[i] = coolwarm257[pos];
    //     }// for
    // }// if
    // else
    // {
    //     //
    //     // interpolate between colours in interval [0,256]
    //     //
        
    //     for ( size_t  i = 0; i < n; ++i )
    //     {
    //         const size_t    pos1 = size_t( std::floor( double(i)/double(n) * 256.0 ) );
    //         const size_t    pos2 = size_t( std::floor( double(i+1)/double(n) * 256.0 ) );
    //         const colour_t  c1 = coolwarm257[pos1];
    //         const colour_t  c2 = coolwarm257[pos2];;
    //         const colour_t  c( (c1.red   + c2.red  ) / 2.0,
    //                            (c1.green + c2.green) / 2.0,
    //                            (c1.blue  + c2.blue ) / 2.0 );
                               
    //         _map[i] = c;
    //     }// for
    // }// else
}

//
// TMagmaColourMap
//

namespace
{

colour_t  magma_data[] = {
    colour_t( 0, 0, 3 ), 
    colour_t( 0, 0, 4 ), 
    colour_t( 0, 0, 6 ), 
    colour_t( 1, 0, 7 ), 
    colour_t( 1, 1, 9 ), 
    colour_t( 1, 1, 11 ), 
    colour_t( 2, 2, 13 ), 
    colour_t( 2, 2, 15 ), 
    colour_t( 3, 3, 17 ), 
    colour_t( 4, 3, 19 ), 
    colour_t( 4, 4, 21 ), 
    colour_t( 5, 4, 23 ), 
    colour_t( 6, 5, 25 ), 
    colour_t( 7, 5, 27 ), 
    colour_t( 8, 6, 29 ), 
    colour_t( 9, 7, 31 ), 
    colour_t( 10, 7, 34 ), 
    colour_t( 11, 8, 36 ), 
    colour_t( 12, 9, 38 ), 
    colour_t( 13, 10, 40 ), 
    colour_t( 14, 10, 42 ), 
    colour_t( 15, 11, 44 ), 
    colour_t( 16, 12, 47 ), 
    colour_t( 17, 12, 49 ), 
    colour_t( 18, 13, 51 ), 
    colour_t( 20, 13, 53 ), 
    colour_t( 21, 14, 56 ), 
    colour_t( 22, 14, 58 ), 
    colour_t( 23, 15, 60 ), 
    colour_t( 24, 15, 63 ), 
    colour_t( 26, 16, 65 ), 
    colour_t( 27, 16, 68 ), 
    colour_t( 28, 16, 70 ), 
    colour_t( 30, 16, 73 ), 
    colour_t( 31, 17, 75 ), 
    colour_t( 32, 17, 77 ), 
    colour_t( 34, 17, 80 ), 
    colour_t( 35, 17, 82 ), 
    colour_t( 37, 17, 85 ), 
    colour_t( 38, 17, 87 ), 
    colour_t( 40, 17, 89 ), 
    colour_t( 42, 17, 92 ), 
    colour_t( 43, 17, 94 ), 
    colour_t( 45, 16, 96 ), 
    colour_t( 47, 16, 98 ), 
    colour_t( 48, 16, 101 ), 
    colour_t( 50, 16, 103 ), 
    colour_t( 52, 16, 104 ), 
    colour_t( 53, 15, 106 ), 
    colour_t( 55, 15, 108 ), 
    colour_t( 57, 15, 110 ), 
    colour_t( 59, 15, 111 ), 
    colour_t( 60, 15, 113 ), 
    colour_t( 62, 15, 114 ), 
    colour_t( 64, 15, 115 ), 
    colour_t( 66, 15, 116 ), 
    colour_t( 67, 15, 117 ), 
    colour_t( 69, 15, 118 ), 
    colour_t( 71, 15, 119 ), 
    colour_t( 72, 16, 120 ), 
    colour_t( 74, 16, 121 ), 
    colour_t( 75, 16, 121 ), 
    colour_t( 77, 17, 122 ), 
    colour_t( 79, 17, 123 ), 
    colour_t( 80, 18, 123 ), 
    colour_t( 82, 18, 124 ), 
    colour_t( 83, 19, 124 ), 
    colour_t( 85, 19, 125 ), 
    colour_t( 87, 20, 125 ), 
    colour_t( 88, 21, 126 ), 
    colour_t( 90, 21, 126 ), 
    colour_t( 91, 22, 126 ), 
    colour_t( 93, 23, 126 ), 
    colour_t( 94, 23, 127 ), 
    colour_t( 96, 24, 127 ), 
    colour_t( 97, 24, 127 ), 
    colour_t( 99, 25, 127 ), 
    colour_t( 101, 26, 128 ), 
    colour_t( 102, 26, 128 ), 
    colour_t( 104, 27, 128 ), 
    colour_t( 105, 28, 128 ), 
    colour_t( 107, 28, 128 ), 
    colour_t( 108, 29, 128 ), 
    colour_t( 110, 30, 129 ), 
    colour_t( 111, 30, 129 ), 
    colour_t( 113, 31, 129 ), 
    colour_t( 115, 31, 129 ), 
    colour_t( 116, 32, 129 ), 
    colour_t( 118, 33, 129 ), 
    colour_t( 119, 33, 129 ), 
    colour_t( 121, 34, 129 ), 
    colour_t( 122, 34, 129 ), 
    colour_t( 124, 35, 129 ), 
    colour_t( 126, 36, 129 ), 
    colour_t( 127, 36, 129 ), 
    colour_t( 129, 37, 129 ), 
    colour_t( 130, 37, 129 ), 
    colour_t( 132, 38, 129 ), 
    colour_t( 133, 38, 129 ), 
    colour_t( 135, 39, 129 ), 
    colour_t( 137, 40, 129 ), 
    colour_t( 138, 40, 129 ), 
    colour_t( 140, 41, 128 ), 
    colour_t( 141, 41, 128 ), 
    colour_t( 143, 42, 128 ), 
    colour_t( 145, 42, 128 ), 
    colour_t( 146, 43, 128 ), 
    colour_t( 148, 43, 128 ), 
    colour_t( 149, 44, 128 ), 
    colour_t( 151, 44, 127 ), 
    colour_t( 153, 45, 127 ), 
    colour_t( 154, 45, 127 ), 
    colour_t( 156, 46, 127 ), 
    colour_t( 158, 46, 126 ), 
    colour_t( 159, 47, 126 ), 
    colour_t( 161, 47, 126 ), 
    colour_t( 163, 48, 126 ), 
    colour_t( 164, 48, 125 ), 
    colour_t( 166, 49, 125 ), 
    colour_t( 167, 49, 125 ), 
    colour_t( 169, 50, 124 ), 
    colour_t( 171, 51, 124 ), 
    colour_t( 172, 51, 123 ), 
    colour_t( 174, 52, 123 ), 
    colour_t( 176, 52, 123 ), 
    colour_t( 177, 53, 122 ), 
    colour_t( 179, 53, 122 ), 
    colour_t( 181, 54, 121 ), 
    colour_t( 182, 54, 121 ), 
    colour_t( 184, 55, 120 ), 
    colour_t( 185, 55, 120 ), 
    colour_t( 187, 56, 119 ), 
    colour_t( 189, 57, 119 ), 
    colour_t( 190, 57, 118 ), 
    colour_t( 192, 58, 117 ), 
    colour_t( 194, 58, 117 ), 
    colour_t( 195, 59, 116 ), 
    colour_t( 197, 60, 116 ), 
    colour_t( 198, 60, 115 ), 
    colour_t( 200, 61, 114 ), 
    colour_t( 202, 62, 114 ), 
    colour_t( 203, 62, 113 ), 
    colour_t( 205, 63, 112 ), 
    colour_t( 206, 64, 112 ), 
    colour_t( 208, 65, 111 ), 
    colour_t( 209, 66, 110 ), 
    colour_t( 211, 66, 109 ), 
    colour_t( 212, 67, 109 ), 
    colour_t( 214, 68, 108 ), 
    colour_t( 215, 69, 107 ), 
    colour_t( 217, 70, 106 ), 
    colour_t( 218, 71, 105 ), 
    colour_t( 220, 72, 105 ), 
    colour_t( 221, 73, 104 ), 
    colour_t( 222, 74, 103 ), 
    colour_t( 224, 75, 102 ), 
    colour_t( 225, 76, 102 ), 
    colour_t( 226, 77, 101 ), 
    colour_t( 228, 78, 100 ), 
    colour_t( 229, 80, 99 ), 
    colour_t( 230, 81, 98 ), 
    colour_t( 231, 82, 98 ), 
    colour_t( 232, 84, 97 ), 
    colour_t( 234, 85, 96 ), 
    colour_t( 235, 86, 96 ), 
    colour_t( 236, 88, 95 ), 
    colour_t( 237, 89, 95 ), 
    colour_t( 238, 91, 94 ), 
    colour_t( 238, 93, 93 ), 
    colour_t( 239, 94, 93 ), 
    colour_t( 240, 96, 93 ), 
    colour_t( 241, 97, 92 ), 
    colour_t( 242, 99, 92 ), 
    colour_t( 243, 101, 92 ), 
    colour_t( 243, 103, 91 ), 
    colour_t( 244, 104, 91 ), 
    colour_t( 245, 106, 91 ), 
    colour_t( 245, 108, 91 ), 
    colour_t( 246, 110, 91 ), 
    colour_t( 246, 112, 91 ), 
    colour_t( 247, 113, 91 ), 
    colour_t( 247, 115, 92 ), 
    colour_t( 248, 117, 92 ), 
    colour_t( 248, 119, 92 ), 
    colour_t( 249, 121, 92 ), 
    colour_t( 249, 123, 93 ), 
    colour_t( 249, 125, 93 ), 
    colour_t( 250, 127, 94 ), 
    colour_t( 250, 128, 94 ), 
    colour_t( 250, 130, 95 ), 
    colour_t( 251, 132, 96 ), 
    colour_t( 251, 134, 96 ), 
    colour_t( 251, 136, 97 ), 
    colour_t( 251, 138, 98 ), 
    colour_t( 252, 140, 99 ), 
    colour_t( 252, 142, 99 ), 
    colour_t( 252, 144, 100 ), 
    colour_t( 252, 146, 101 ), 
    colour_t( 252, 147, 102 ), 
    colour_t( 253, 149, 103 ), 
    colour_t( 253, 151, 104 ), 
    colour_t( 253, 153, 105 ), 
    colour_t( 253, 155, 106 ), 
    colour_t( 253, 157, 107 ), 
    colour_t( 253, 159, 108 ), 
    colour_t( 253, 161, 110 ), 
    colour_t( 253, 162, 111 ), 
    colour_t( 253, 164, 112 ), 
    colour_t( 254, 166, 113 ), 
    colour_t( 254, 168, 115 ), 
    colour_t( 254, 170, 116 ), 
    colour_t( 254, 172, 117 ), 
    colour_t( 254, 174, 118 ), 
    colour_t( 254, 175, 120 ), 
    colour_t( 254, 177, 121 ), 
    colour_t( 254, 179, 123 ), 
    colour_t( 254, 181, 124 ), 
    colour_t( 254, 183, 125 ), 
    colour_t( 254, 185, 127 ), 
    colour_t( 254, 187, 128 ), 
    colour_t( 254, 188, 130 ), 
    colour_t( 254, 190, 131 ), 
    colour_t( 254, 192, 133 ), 
    colour_t( 254, 194, 134 ), 
    colour_t( 254, 196, 136 ), 
    colour_t( 254, 198, 137 ), 
    colour_t( 254, 199, 139 ), 
    colour_t( 254, 201, 141 ), 
    colour_t( 254, 203, 142 ), 
    colour_t( 253, 205, 144 ), 
    colour_t( 253, 207, 146 ), 
    colour_t( 253, 209, 147 ), 
    colour_t( 253, 210, 149 ), 
    colour_t( 253, 212, 151 ), 
    colour_t( 253, 214, 152 ), 
    colour_t( 253, 216, 154 ), 
    colour_t( 253, 218, 156 ), 
    colour_t( 253, 220, 157 ), 
    colour_t( 253, 221, 159 ), 
    colour_t( 253, 223, 161 ), 
    colour_t( 253, 225, 163 ), 
    colour_t( 252, 227, 165 ), 
    colour_t( 252, 229, 166 ), 
    colour_t( 252, 230, 168 ), 
    colour_t( 252, 232, 170 ), 
    colour_t( 252, 234, 172 ), 
    colour_t( 252, 236, 174 ), 
    colour_t( 252, 238, 176 ), 
    colour_t( 252, 240, 177 ), 
    colour_t( 252, 241, 179 ), 
    colour_t( 252, 243, 181 ), 
    colour_t( 252, 245, 183 ), 
    colour_t( 251, 247, 185 ), 
    colour_t( 251, 249, 187 ), 
    colour_t( 251, 250, 189 ), 
    colour_t( 251, 252, 191 )    
};

}// namespace anonymous

void
TMagmaColourMap::init ( const size_t n )
{
    init_cmap_with_data( n, _map, magma_data, sizeof(magma_data) / sizeof(colour_t) );
}

//
// TInfernoColourMap
//

namespace
{

colour_t  inferno_data[] = {
    colour_t( 0, 0, 3 ), 
    colour_t( 0, 0, 4 ), 
    colour_t( 0, 0, 6 ), 
    colour_t( 1, 0, 7 ), 
    colour_t( 1, 1, 9 ), 
    colour_t( 1, 1, 11 ), 
    colour_t( 2, 1, 14 ), 
    colour_t( 2, 2, 16 ), 
    colour_t( 3, 2, 18 ), 
    colour_t( 4, 3, 20 ), 
    colour_t( 4, 3, 22 ), 
    colour_t( 5, 4, 24 ), 
    colour_t( 6, 4, 27 ), 
    colour_t( 7, 5, 29 ), 
    colour_t( 8, 6, 31 ), 
    colour_t( 9, 6, 33 ), 
    colour_t( 10, 7, 35 ), 
    colour_t( 11, 7, 38 ), 
    colour_t( 13, 8, 40 ), 
    colour_t( 14, 8, 42 ), 
    colour_t( 15, 9, 45 ), 
    colour_t( 16, 9, 47 ), 
    colour_t( 18, 10, 50 ), 
    colour_t( 19, 10, 52 ), 
    colour_t( 20, 11, 54 ), 
    colour_t( 22, 11, 57 ), 
    colour_t( 23, 11, 59 ), 
    colour_t( 25, 11, 62 ), 
    colour_t( 26, 11, 64 ), 
    colour_t( 28, 12, 67 ), 
    colour_t( 29, 12, 69 ), 
    colour_t( 31, 12, 71 ), 
    colour_t( 32, 12, 74 ), 
    colour_t( 34, 11, 76 ), 
    colour_t( 36, 11, 78 ), 
    colour_t( 38, 11, 80 ), 
    colour_t( 39, 11, 82 ), 
    colour_t( 41, 11, 84 ), 
    colour_t( 43, 10, 86 ), 
    colour_t( 45, 10, 88 ), 
    colour_t( 46, 10, 90 ), 
    colour_t( 48, 10, 92 ), 
    colour_t( 50, 9, 93 ), 
    colour_t( 52, 9, 95 ), 
    colour_t( 53, 9, 96 ), 
    colour_t( 55, 9, 97 ), 
    colour_t( 57, 9, 98 ), 
    colour_t( 59, 9, 100 ), 
    colour_t( 60, 9, 101 ), 
    colour_t( 62, 9, 102 ), 
    colour_t( 64, 9, 102 ), 
    colour_t( 65, 9, 103 ), 
    colour_t( 67, 10, 104 ), 
    colour_t( 69, 10, 105 ), 
    colour_t( 70, 10, 105 ), 
    colour_t( 72, 11, 106 ), 
    colour_t( 74, 11, 106 ), 
    colour_t( 75, 12, 107 ), 
    colour_t( 77, 12, 107 ), 
    colour_t( 79, 13, 108 ), 
    colour_t( 80, 13, 108 ), 
    colour_t( 82, 14, 108 ), 
    colour_t( 83, 14, 109 ), 
    colour_t( 85, 15, 109 ), 
    colour_t( 87, 15, 109 ), 
    colour_t( 88, 16, 109 ), 
    colour_t( 90, 17, 109 ), 
    colour_t( 91, 17, 110 ), 
    colour_t( 93, 18, 110 ), 
    colour_t( 95, 18, 110 ), 
    colour_t( 96, 19, 110 ), 
    colour_t( 98, 20, 110 ), 
    colour_t( 99, 20, 110 ), 
    colour_t( 101, 21, 110 ), 
    colour_t( 102, 21, 110 ), 
    colour_t( 104, 22, 110 ), 
    colour_t( 106, 23, 110 ), 
    colour_t( 107, 23, 110 ), 
    colour_t( 109, 24, 110 ), 
    colour_t( 110, 24, 110 ), 
    colour_t( 112, 25, 110 ), 
    colour_t( 114, 25, 109 ), 
    colour_t( 115, 26, 109 ), 
    colour_t( 117, 27, 109 ), 
    colour_t( 118, 27, 109 ), 
    colour_t( 120, 28, 109 ), 
    colour_t( 122, 28, 109 ), 
    colour_t( 123, 29, 108 ), 
    colour_t( 125, 29, 108 ), 
    colour_t( 126, 30, 108 ), 
    colour_t( 128, 31, 107 ), 
    colour_t( 129, 31, 107 ), 
    colour_t( 131, 32, 107 ), 
    colour_t( 133, 32, 106 ), 
    colour_t( 134, 33, 106 ), 
    colour_t( 136, 33, 106 ), 
    colour_t( 137, 34, 105 ), 
    colour_t( 139, 34, 105 ), 
    colour_t( 141, 35, 105 ), 
    colour_t( 142, 36, 104 ), 
    colour_t( 144, 36, 104 ), 
    colour_t( 145, 37, 103 ), 
    colour_t( 147, 37, 103 ), 
    colour_t( 149, 38, 102 ), 
    colour_t( 150, 38, 102 ), 
    colour_t( 152, 39, 101 ), 
    colour_t( 153, 40, 100 ), 
    colour_t( 155, 40, 100 ), 
    colour_t( 156, 41, 99 ), 
    colour_t( 158, 41, 99 ), 
    colour_t( 160, 42, 98 ), 
    colour_t( 161, 43, 97 ), 
    colour_t( 163, 43, 97 ), 
    colour_t( 164, 44, 96 ), 
    colour_t( 166, 44, 95 ), 
    colour_t( 167, 45, 95 ), 
    colour_t( 169, 46, 94 ), 
    colour_t( 171, 46, 93 ), 
    colour_t( 172, 47, 92 ), 
    colour_t( 174, 48, 91 ), 
    colour_t( 175, 49, 91 ), 
    colour_t( 177, 49, 90 ), 
    colour_t( 178, 50, 89 ), 
    colour_t( 180, 51, 88 ), 
    colour_t( 181, 51, 87 ), 
    colour_t( 183, 52, 86 ), 
    colour_t( 184, 53, 86 ), 
    colour_t( 186, 54, 85 ), 
    colour_t( 187, 55, 84 ), 
    colour_t( 189, 55, 83 ), 
    colour_t( 190, 56, 82 ), 
    colour_t( 191, 57, 81 ), 
    colour_t( 193, 58, 80 ), 
    colour_t( 194, 59, 79 ), 
    colour_t( 196, 60, 78 ), 
    colour_t( 197, 61, 77 ), 
    colour_t( 199, 62, 76 ), 
    colour_t( 200, 62, 75 ), 
    colour_t( 201, 63, 74 ), 
    colour_t( 203, 64, 73 ), 
    colour_t( 204, 65, 72 ), 
    colour_t( 205, 66, 71 ), 
    colour_t( 207, 68, 70 ), 
    colour_t( 208, 69, 68 ), 
    colour_t( 209, 70, 67 ), 
    colour_t( 210, 71, 66 ), 
    colour_t( 212, 72, 65 ), 
    colour_t( 213, 73, 64 ), 
    colour_t( 214, 74, 63 ), 
    colour_t( 215, 75, 62 ), 
    colour_t( 217, 77, 61 ), 
    colour_t( 218, 78, 59 ), 
    colour_t( 219, 79, 58 ), 
    colour_t( 220, 80, 57 ), 
    colour_t( 221, 82, 56 ), 
    colour_t( 222, 83, 55 ), 
    colour_t( 223, 84, 54 ), 
    colour_t( 224, 86, 52 ), 
    colour_t( 226, 87, 51 ), 
    colour_t( 227, 88, 50 ), 
    colour_t( 228, 90, 49 ), 
    colour_t( 229, 91, 48 ), 
    colour_t( 230, 92, 46 ), 
    colour_t( 230, 94, 45 ), 
    colour_t( 231, 95, 44 ), 
    colour_t( 232, 97, 43 ), 
    colour_t( 233, 98, 42 ), 
    colour_t( 234, 100, 40 ), 
    colour_t( 235, 101, 39 ), 
    colour_t( 236, 103, 38 ), 
    colour_t( 237, 104, 37 ), 
    colour_t( 237, 106, 35 ), 
    colour_t( 238, 108, 34 ), 
    colour_t( 239, 109, 33 ), 
    colour_t( 240, 111, 31 ), 
    colour_t( 240, 112, 30 ), 
    colour_t( 241, 114, 29 ), 
    colour_t( 242, 116, 28 ), 
    colour_t( 242, 117, 26 ), 
    colour_t( 243, 119, 25 ), 
    colour_t( 243, 121, 24 ), 
    colour_t( 244, 122, 22 ), 
    colour_t( 245, 124, 21 ), 
    colour_t( 245, 126, 20 ), 
    colour_t( 246, 128, 18 ), 
    colour_t( 246, 129, 17 ), 
    colour_t( 247, 131, 16 ), 
    colour_t( 247, 133, 14 ), 
    colour_t( 248, 135, 13 ), 
    colour_t( 248, 136, 12 ), 
    colour_t( 248, 138, 11 ), 
    colour_t( 249, 140, 9 ), 
    colour_t( 249, 142, 8 ), 
    colour_t( 249, 144, 8 ), 
    colour_t( 250, 145, 7 ), 
    colour_t( 250, 147, 6 ), 
    colour_t( 250, 149, 6 ), 
    colour_t( 250, 151, 6 ), 
    colour_t( 251, 153, 6 ), 
    colour_t( 251, 155, 6 ), 
    colour_t( 251, 157, 6 ), 
    colour_t( 251, 158, 7 ), 
    colour_t( 251, 160, 7 ), 
    colour_t( 251, 162, 8 ), 
    colour_t( 251, 164, 10 ), 
    colour_t( 251, 166, 11 ), 
    colour_t( 251, 168, 13 ), 
    colour_t( 251, 170, 14 ), 
    colour_t( 251, 172, 16 ), 
    colour_t( 251, 174, 18 ), 
    colour_t( 251, 176, 20 ), 
    colour_t( 251, 177, 22 ), 
    colour_t( 251, 179, 24 ), 
    colour_t( 251, 181, 26 ), 
    colour_t( 251, 183, 28 ), 
    colour_t( 251, 185, 30 ), 
    colour_t( 250, 187, 33 ), 
    colour_t( 250, 189, 35 ), 
    colour_t( 250, 191, 37 ), 
    colour_t( 250, 193, 40 ), 
    colour_t( 249, 195, 42 ), 
    colour_t( 249, 197, 44 ), 
    colour_t( 249, 199, 47 ), 
    colour_t( 248, 201, 49 ), 
    colour_t( 248, 203, 52 ), 
    colour_t( 248, 205, 55 ), 
    colour_t( 247, 207, 58 ), 
    colour_t( 247, 209, 60 ), 
    colour_t( 246, 211, 63 ), 
    colour_t( 246, 213, 66 ), 
    colour_t( 245, 215, 69 ), 
    colour_t( 245, 217, 72 ), 
    colour_t( 244, 219, 75 ), 
    colour_t( 244, 220, 79 ), 
    colour_t( 243, 222, 82 ), 
    colour_t( 243, 224, 86 ), 
    colour_t( 243, 226, 89 ), 
    colour_t( 242, 228, 93 ), 
    colour_t( 242, 230, 96 ), 
    colour_t( 241, 232, 100 ), 
    colour_t( 241, 233, 104 ), 
    colour_t( 241, 235, 108 ), 
    colour_t( 241, 237, 112 ), 
    colour_t( 241, 238, 116 ), 
    colour_t( 241, 240, 121 ), 
    colour_t( 241, 242, 125 ), 
    colour_t( 242, 243, 129 ), 
    colour_t( 242, 244, 133 ), 
    colour_t( 243, 246, 137 ), 
    colour_t( 244, 247, 141 ), 
    colour_t( 245, 248, 145 ), 
    colour_t( 246, 250, 149 ), 
    colour_t( 247, 251, 153 ), 
    colour_t( 249, 252, 157 ), 
    colour_t( 250, 253, 160 ), 
    colour_t( 252, 254, 164 )
};

}// namespace anonymous

void
TInfernoColourMap::init ( const size_t n )
{
    init_cmap_with_data( n, _map, inferno_data, sizeof(inferno_data) / sizeof(colour_t) );
}

//
// TPlasmaColourMap
//

namespace
{

colour_t  plasma_data[] = {
    colour_t( 12, 7, 134 ), 
    colour_t( 16, 7, 135 ), 
    colour_t( 19, 6, 137 ), 
    colour_t( 21, 6, 138 ), 
    colour_t( 24, 6, 139 ), 
    colour_t( 27, 6, 140 ), 
    colour_t( 29, 6, 141 ), 
    colour_t( 31, 5, 142 ), 
    colour_t( 33, 5, 143 ), 
    colour_t( 35, 5, 144 ), 
    colour_t( 37, 5, 145 ), 
    colour_t( 39, 5, 146 ), 
    colour_t( 41, 5, 147 ), 
    colour_t( 43, 5, 148 ), 
    colour_t( 45, 4, 148 ), 
    colour_t( 47, 4, 149 ), 
    colour_t( 49, 4, 150 ), 
    colour_t( 51, 4, 151 ), 
    colour_t( 52, 4, 152 ), 
    colour_t( 54, 4, 152 ), 
    colour_t( 56, 4, 153 ), 
    colour_t( 58, 4, 154 ), 
    colour_t( 59, 3, 154 ), 
    colour_t( 61, 3, 155 ), 
    colour_t( 63, 3, 156 ), 
    colour_t( 64, 3, 156 ), 
    colour_t( 66, 3, 157 ), 
    colour_t( 68, 3, 158 ), 
    colour_t( 69, 3, 158 ), 
    colour_t( 71, 2, 159 ), 
    colour_t( 73, 2, 159 ), 
    colour_t( 74, 2, 160 ), 
    colour_t( 76, 2, 161 ), 
    colour_t( 78, 2, 161 ), 
    colour_t( 79, 2, 162 ), 
    colour_t( 81, 1, 162 ), 
    colour_t( 82, 1, 163 ), 
    colour_t( 84, 1, 163 ), 
    colour_t( 86, 1, 163 ), 
    colour_t( 87, 1, 164 ), 
    colour_t( 89, 1, 164 ), 
    colour_t( 90, 0, 165 ), 
    colour_t( 92, 0, 165 ), 
    colour_t( 94, 0, 165 ), 
    colour_t( 95, 0, 166 ), 
    colour_t( 97, 0, 166 ), 
    colour_t( 98, 0, 166 ), 
    colour_t( 100, 0, 167 ), 
    colour_t( 101, 0, 167 ), 
    colour_t( 103, 0, 167 ), 
    colour_t( 104, 0, 167 ), 
    colour_t( 106, 0, 167 ), 
    colour_t( 108, 0, 168 ), 
    colour_t( 109, 0, 168 ), 
    colour_t( 111, 0, 168 ), 
    colour_t( 112, 0, 168 ), 
    colour_t( 114, 0, 168 ), 
    colour_t( 115, 0, 168 ), 
    colour_t( 117, 0, 168 ), 
    colour_t( 118, 1, 168 ), 
    colour_t( 120, 1, 168 ), 
    colour_t( 121, 1, 168 ), 
    colour_t( 123, 2, 168 ), 
    colour_t( 124, 2, 167 ), 
    colour_t( 126, 3, 167 ), 
    colour_t( 127, 3, 167 ), 
    colour_t( 129, 4, 167 ), 
    colour_t( 130, 4, 167 ), 
    colour_t( 132, 5, 166 ), 
    colour_t( 133, 6, 166 ), 
    colour_t( 134, 7, 166 ), 
    colour_t( 136, 7, 165 ), 
    colour_t( 137, 8, 165 ), 
    colour_t( 139, 9, 164 ), 
    colour_t( 140, 10, 164 ), 
    colour_t( 142, 12, 164 ), 
    colour_t( 143, 13, 163 ), 
    colour_t( 144, 14, 163 ), 
    colour_t( 146, 15, 162 ), 
    colour_t( 147, 16, 161 ), 
    colour_t( 149, 17, 161 ), 
    colour_t( 150, 18, 160 ), 
    colour_t( 151, 19, 160 ), 
    colour_t( 153, 20, 159 ), 
    colour_t( 154, 21, 158 ), 
    colour_t( 155, 23, 158 ), 
    colour_t( 157, 24, 157 ), 
    colour_t( 158, 25, 156 ), 
    colour_t( 159, 26, 155 ), 
    colour_t( 160, 27, 155 ), 
    colour_t( 162, 28, 154 ), 
    colour_t( 163, 29, 153 ), 
    colour_t( 164, 30, 152 ), 
    colour_t( 165, 31, 151 ), 
    colour_t( 167, 33, 151 ), 
    colour_t( 168, 34, 150 ), 
    colour_t( 169, 35, 149 ), 
    colour_t( 170, 36, 148 ), 
    colour_t( 172, 37, 147 ), 
    colour_t( 173, 38, 146 ), 
    colour_t( 174, 39, 145 ), 
    colour_t( 175, 40, 144 ), 
    colour_t( 176, 42, 143 ), 
    colour_t( 177, 43, 143 ), 
    colour_t( 178, 44, 142 ), 
    colour_t( 180, 45, 141 ), 
    colour_t( 181, 46, 140 ), 
    colour_t( 182, 47, 139 ), 
    colour_t( 183, 48, 138 ), 
    colour_t( 184, 50, 137 ), 
    colour_t( 185, 51, 136 ), 
    colour_t( 186, 52, 135 ), 
    colour_t( 187, 53, 134 ), 
    colour_t( 188, 54, 133 ), 
    colour_t( 189, 55, 132 ), 
    colour_t( 190, 56, 131 ), 
    colour_t( 191, 57, 130 ), 
    colour_t( 192, 59, 129 ), 
    colour_t( 193, 60, 128 ), 
    colour_t( 194, 61, 128 ), 
    colour_t( 195, 62, 127 ), 
    colour_t( 196, 63, 126 ), 
    colour_t( 197, 64, 125 ), 
    colour_t( 198, 65, 124 ), 
    colour_t( 199, 66, 123 ), 
    colour_t( 200, 68, 122 ), 
    colour_t( 201, 69, 121 ), 
    colour_t( 202, 70, 120 ), 
    colour_t( 203, 71, 119 ), 
    colour_t( 204, 72, 118 ), 
    colour_t( 205, 73, 117 ), 
    colour_t( 206, 74, 117 ), 
    colour_t( 207, 75, 116 ), 
    colour_t( 208, 77, 115 ), 
    colour_t( 209, 78, 114 ), 
    colour_t( 209, 79, 113 ), 
    colour_t( 210, 80, 112 ), 
    colour_t( 211, 81, 111 ), 
    colour_t( 212, 82, 110 ), 
    colour_t( 213, 83, 109 ), 
    colour_t( 214, 85, 109 ), 
    colour_t( 215, 86, 108 ), 
    colour_t( 215, 87, 107 ), 
    colour_t( 216, 88, 106 ), 
    colour_t( 217, 89, 105 ), 
    colour_t( 218, 90, 104 ), 
    colour_t( 219, 91, 103 ), 
    colour_t( 220, 93, 102 ), 
    colour_t( 220, 94, 102 ), 
    colour_t( 221, 95, 101 ), 
    colour_t( 222, 96, 100 ), 
    colour_t( 223, 97, 99 ), 
    colour_t( 223, 98, 98 ), 
    colour_t( 224, 100, 97 ), 
    colour_t( 225, 101, 96 ), 
    colour_t( 226, 102, 96 ), 
    colour_t( 227, 103, 95 ), 
    colour_t( 227, 104, 94 ), 
    colour_t( 228, 106, 93 ), 
    colour_t( 229, 107, 92 ), 
    colour_t( 229, 108, 91 ), 
    colour_t( 230, 109, 90 ), 
    colour_t( 231, 110, 90 ), 
    colour_t( 232, 112, 89 ), 
    colour_t( 232, 113, 88 ), 
    colour_t( 233, 114, 87 ), 
    colour_t( 234, 115, 86 ), 
    colour_t( 234, 116, 85 ), 
    colour_t( 235, 118, 84 ), 
    colour_t( 236, 119, 84 ), 
    colour_t( 236, 120, 83 ), 
    colour_t( 237, 121, 82 ), 
    colour_t( 237, 123, 81 ), 
    colour_t( 238, 124, 80 ), 
    colour_t( 239, 125, 79 ), 
    colour_t( 239, 126, 78 ), 
    colour_t( 240, 128, 77 ), 
    colour_t( 240, 129, 77 ), 
    colour_t( 241, 130, 76 ), 
    colour_t( 242, 132, 75 ), 
    colour_t( 242, 133, 74 ), 
    colour_t( 243, 134, 73 ), 
    colour_t( 243, 135, 72 ), 
    colour_t( 244, 137, 71 ), 
    colour_t( 244, 138, 71 ), 
    colour_t( 245, 139, 70 ), 
    colour_t( 245, 141, 69 ), 
    colour_t( 246, 142, 68 ), 
    colour_t( 246, 143, 67 ), 
    colour_t( 246, 145, 66 ), 
    colour_t( 247, 146, 65 ), 
    colour_t( 247, 147, 65 ), 
    colour_t( 248, 149, 64 ), 
    colour_t( 248, 150, 63 ), 
    colour_t( 248, 152, 62 ), 
    colour_t( 249, 153, 61 ), 
    colour_t( 249, 154, 60 ), 
    colour_t( 250, 156, 59 ), 
    colour_t( 250, 157, 58 ), 
    colour_t( 250, 159, 58 ), 
    colour_t( 250, 160, 57 ), 
    colour_t( 251, 162, 56 ), 
    colour_t( 251, 163, 55 ), 
    colour_t( 251, 164, 54 ), 
    colour_t( 252, 166, 53 ), 
    colour_t( 252, 167, 53 ), 
    colour_t( 252, 169, 52 ), 
    colour_t( 252, 170, 51 ), 
    colour_t( 252, 172, 50 ), 
    colour_t( 252, 173, 49 ), 
    colour_t( 253, 175, 49 ), 
    colour_t( 253, 176, 48 ), 
    colour_t( 253, 178, 47 ), 
    colour_t( 253, 179, 46 ), 
    colour_t( 253, 181, 45 ), 
    colour_t( 253, 182, 45 ), 
    colour_t( 253, 184, 44 ), 
    colour_t( 253, 185, 43 ), 
    colour_t( 253, 187, 43 ), 
    colour_t( 253, 188, 42 ), 
    colour_t( 253, 190, 41 ), 
    colour_t( 253, 192, 41 ), 
    colour_t( 253, 193, 40 ), 
    colour_t( 253, 195, 40 ), 
    colour_t( 253, 196, 39 ), 
    colour_t( 253, 198, 38 ), 
    colour_t( 252, 199, 38 ), 
    colour_t( 252, 201, 38 ), 
    colour_t( 252, 203, 37 ), 
    colour_t( 252, 204, 37 ), 
    colour_t( 252, 206, 37 ), 
    colour_t( 251, 208, 36 ), 
    colour_t( 251, 209, 36 ), 
    colour_t( 251, 211, 36 ), 
    colour_t( 250, 213, 36 ), 
    colour_t( 250, 214, 36 ), 
    colour_t( 250, 216, 36 ), 
    colour_t( 249, 217, 36 ), 
    colour_t( 249, 219, 36 ), 
    colour_t( 248, 221, 36 ), 
    colour_t( 248, 223, 36 ), 
    colour_t( 247, 224, 36 ), 
    colour_t( 247, 226, 37 ), 
    colour_t( 246, 228, 37 ), 
    colour_t( 246, 229, 37 ), 
    colour_t( 245, 231, 38 ), 
    colour_t( 245, 233, 38 ), 
    colour_t( 244, 234, 38 ), 
    colour_t( 243, 236, 38 ), 
    colour_t( 243, 238, 38 ), 
    colour_t( 242, 240, 38 ), 
    colour_t( 242, 241, 38 ), 
    colour_t( 241, 243, 38 ), 
    colour_t( 240, 245, 37 ), 
    colour_t( 240, 246, 35 ), 
    colour_t( 239, 248, 33 )
};

}// namespace anonymous

void
TPlasmaColourMap::init ( const size_t n )
{
    init_cmap_with_data( n, _map, plasma_data, sizeof(plasma_data) / sizeof(colour_t) );
}

//
// TViridisColourMap
//

namespace
{

colour_t  viridis_data[] = {
    colour_t( 68, 1, 84 ), 
    colour_t( 68, 2, 85 ), 
    colour_t( 68, 3, 87 ), 
    colour_t( 69, 5, 88 ), 
    colour_t( 69, 6, 90 ), 
    colour_t( 69, 8, 91 ), 
    colour_t( 70, 9, 92 ), 
    colour_t( 70, 11, 94 ), 
    colour_t( 70, 12, 95 ), 
    colour_t( 70, 14, 97 ), 
    colour_t( 71, 15, 98 ), 
    colour_t( 71, 17, 99 ), 
    colour_t( 71, 18, 101 ), 
    colour_t( 71, 20, 102 ), 
    colour_t( 71, 21, 103 ), 
    colour_t( 71, 22, 105 ), 
    colour_t( 71, 24, 106 ), 
    colour_t( 72, 25, 107 ), 
    colour_t( 72, 26, 108 ), 
    colour_t( 72, 28, 110 ), 
    colour_t( 72, 29, 111 ), 
    colour_t( 72, 30, 112 ), 
    colour_t( 72, 32, 113 ), 
    colour_t( 72, 33, 114 ), 
    colour_t( 72, 34, 115 ), 
    colour_t( 72, 35, 116 ), 
    colour_t( 71, 37, 117 ), 
    colour_t( 71, 38, 118 ), 
    colour_t( 71, 39, 119 ), 
    colour_t( 71, 40, 120 ), 
    colour_t( 71, 42, 121 ), 
    colour_t( 71, 43, 122 ), 
    colour_t( 71, 44, 123 ), 
    colour_t( 70, 45, 124 ), 
    colour_t( 70, 47, 124 ), 
    colour_t( 70, 48, 125 ), 
    colour_t( 70, 49, 126 ), 
    colour_t( 69, 50, 127 ), 
    colour_t( 69, 52, 127 ), 
    colour_t( 69, 53, 128 ), 
    colour_t( 69, 54, 129 ), 
    colour_t( 68, 55, 129 ), 
    colour_t( 68, 57, 130 ), 
    colour_t( 67, 58, 131 ), 
    colour_t( 67, 59, 131 ), 
    colour_t( 67, 60, 132 ), 
    colour_t( 66, 61, 132 ), 
    colour_t( 66, 62, 133 ), 
    colour_t( 66, 64, 133 ), 
    colour_t( 65, 65, 134 ), 
    colour_t( 65, 66, 134 ), 
    colour_t( 64, 67, 135 ), 
    colour_t( 64, 68, 135 ), 
    colour_t( 63, 69, 135 ), 
    colour_t( 63, 71, 136 ), 
    colour_t( 62, 72, 136 ), 
    colour_t( 62, 73, 137 ), 
    colour_t( 61, 74, 137 ), 
    colour_t( 61, 75, 137 ), 
    colour_t( 61, 76, 137 ), 
    colour_t( 60, 77, 138 ), 
    colour_t( 60, 78, 138 ), 
    colour_t( 59, 80, 138 ), 
    colour_t( 59, 81, 138 ), 
    colour_t( 58, 82, 139 ), 
    colour_t( 58, 83, 139 ), 
    colour_t( 57, 84, 139 ), 
    colour_t( 57, 85, 139 ), 
    colour_t( 56, 86, 139 ), 
    colour_t( 56, 87, 140 ), 
    colour_t( 55, 88, 140 ), 
    colour_t( 55, 89, 140 ), 
    colour_t( 54, 90, 140 ), 
    colour_t( 54, 91, 140 ), 
    colour_t( 53, 92, 140 ), 
    colour_t( 53, 93, 140 ), 
    colour_t( 52, 94, 141 ), 
    colour_t( 52, 95, 141 ), 
    colour_t( 51, 96, 141 ), 
    colour_t( 51, 97, 141 ), 
    colour_t( 50, 98, 141 ), 
    colour_t( 50, 99, 141 ), 
    colour_t( 49, 100, 141 ), 
    colour_t( 49, 101, 141 ), 
    colour_t( 49, 102, 141 ), 
    colour_t( 48, 103, 141 ), 
    colour_t( 48, 104, 141 ), 
    colour_t( 47, 105, 141 ), 
    colour_t( 47, 106, 141 ), 
    colour_t( 46, 107, 142 ), 
    colour_t( 46, 108, 142 ), 
    colour_t( 46, 109, 142 ), 
    colour_t( 45, 110, 142 ), 
    colour_t( 45, 111, 142 ), 
    colour_t( 44, 112, 142 ), 
    colour_t( 44, 113, 142 ), 
    colour_t( 44, 114, 142 ), 
    colour_t( 43, 115, 142 ), 
    colour_t( 43, 116, 142 ), 
    colour_t( 42, 117, 142 ), 
    colour_t( 42, 118, 142 ), 
    colour_t( 42, 119, 142 ), 
    colour_t( 41, 120, 142 ), 
    colour_t( 41, 121, 142 ), 
    colour_t( 40, 122, 142 ), 
    colour_t( 40, 122, 142 ), 
    colour_t( 40, 123, 142 ), 
    colour_t( 39, 124, 142 ), 
    colour_t( 39, 125, 142 ), 
    colour_t( 39, 126, 142 ), 
    colour_t( 38, 127, 142 ), 
    colour_t( 38, 128, 142 ), 
    colour_t( 38, 129, 142 ), 
    colour_t( 37, 130, 142 ), 
    colour_t( 37, 131, 141 ), 
    colour_t( 36, 132, 141 ), 
    colour_t( 36, 133, 141 ), 
    colour_t( 36, 134, 141 ), 
    colour_t( 35, 135, 141 ), 
    colour_t( 35, 136, 141 ), 
    colour_t( 35, 137, 141 ), 
    colour_t( 34, 137, 141 ), 
    colour_t( 34, 138, 141 ), 
    colour_t( 34, 139, 141 ), 
    colour_t( 33, 140, 141 ), 
    colour_t( 33, 141, 140 ), 
    colour_t( 33, 142, 140 ), 
    colour_t( 32, 143, 140 ), 
    colour_t( 32, 144, 140 ), 
    colour_t( 32, 145, 140 ), 
    colour_t( 31, 146, 140 ), 
    colour_t( 31, 147, 139 ), 
    colour_t( 31, 148, 139 ), 
    colour_t( 31, 149, 139 ), 
    colour_t( 31, 150, 139 ), 
    colour_t( 30, 151, 138 ), 
    colour_t( 30, 152, 138 ), 
    colour_t( 30, 153, 138 ), 
    colour_t( 30, 153, 138 ), 
    colour_t( 30, 154, 137 ), 
    colour_t( 30, 155, 137 ), 
    colour_t( 30, 156, 137 ), 
    colour_t( 30, 157, 136 ), 
    colour_t( 30, 158, 136 ), 
    colour_t( 30, 159, 136 ), 
    colour_t( 30, 160, 135 ), 
    colour_t( 31, 161, 135 ), 
    colour_t( 31, 162, 134 ), 
    colour_t( 31, 163, 134 ), 
    colour_t( 32, 164, 133 ), 
    colour_t( 32, 165, 133 ), 
    colour_t( 33, 166, 133 ), 
    colour_t( 33, 167, 132 ), 
    colour_t( 34, 167, 132 ), 
    colour_t( 35, 168, 131 ), 
    colour_t( 35, 169, 130 ), 
    colour_t( 36, 170, 130 ), 
    colour_t( 37, 171, 129 ), 
    colour_t( 38, 172, 129 ), 
    colour_t( 39, 173, 128 ), 
    colour_t( 40, 174, 127 ), 
    colour_t( 41, 175, 127 ), 
    colour_t( 42, 176, 126 ), 
    colour_t( 43, 177, 125 ), 
    colour_t( 44, 177, 125 ), 
    colour_t( 46, 178, 124 ), 
    colour_t( 47, 179, 123 ), 
    colour_t( 48, 180, 122 ), 
    colour_t( 50, 181, 122 ), 
    colour_t( 51, 182, 121 ), 
    colour_t( 53, 183, 120 ), 
    colour_t( 54, 184, 119 ), 
    colour_t( 56, 185, 118 ), 
    colour_t( 57, 185, 118 ), 
    colour_t( 59, 186, 117 ), 
    colour_t( 61, 187, 116 ), 
    colour_t( 62, 188, 115 ), 
    colour_t( 64, 189, 114 ), 
    colour_t( 66, 190, 113 ), 
    colour_t( 68, 190, 112 ), 
    colour_t( 69, 191, 111 ), 
    colour_t( 71, 192, 110 ), 
    colour_t( 73, 193, 109 ), 
    colour_t( 75, 194, 108 ), 
    colour_t( 77, 194, 107 ), 
    colour_t( 79, 195, 105 ), 
    colour_t( 81, 196, 104 ), 
    colour_t( 83, 197, 103 ), 
    colour_t( 85, 198, 102 ), 
    colour_t( 87, 198, 101 ), 
    colour_t( 89, 199, 100 ), 
    colour_t( 91, 200, 98 ), 
    colour_t( 94, 201, 97 ), 
    colour_t( 96, 201, 96 ), 
    colour_t( 98, 202, 95 ), 
    colour_t( 100, 203, 93 ), 
    colour_t( 103, 204, 92 ), 
    colour_t( 105, 204, 91 ), 
    colour_t( 107, 205, 89 ), 
    colour_t( 109, 206, 88 ), 
    colour_t( 112, 206, 86 ), 
    colour_t( 114, 207, 85 ), 
    colour_t( 116, 208, 84 ), 
    colour_t( 119, 208, 82 ), 
    colour_t( 121, 209, 81 ), 
    colour_t( 124, 210, 79 ), 
    colour_t( 126, 210, 78 ), 
    colour_t( 129, 211, 76 ), 
    colour_t( 131, 211, 75 ), 
    colour_t( 134, 212, 73 ), 
    colour_t( 136, 213, 71 ), 
    colour_t( 139, 213, 70 ), 
    colour_t( 141, 214, 68 ), 
    colour_t( 144, 214, 67 ), 
    colour_t( 146, 215, 65 ), 
    colour_t( 149, 215, 63 ), 
    colour_t( 151, 216, 62 ), 
    colour_t( 154, 216, 60 ), 
    colour_t( 157, 217, 58 ), 
    colour_t( 159, 217, 56 ), 
    colour_t( 162, 218, 55 ), 
    colour_t( 165, 218, 53 ), 
    colour_t( 167, 219, 51 ), 
    colour_t( 170, 219, 50 ), 
    colour_t( 173, 220, 48 ), 
    colour_t( 175, 220, 46 ), 
    colour_t( 178, 221, 44 ), 
    colour_t( 181, 221, 43 ), 
    colour_t( 183, 221, 41 ), 
    colour_t( 186, 222, 39 ), 
    colour_t( 189, 222, 38 ), 
    colour_t( 191, 223, 36 ), 
    colour_t( 194, 223, 34 ), 
    colour_t( 197, 223, 33 ), 
    colour_t( 199, 224, 31 ), 
    colour_t( 202, 224, 30 ), 
    colour_t( 205, 224, 29 ), 
    colour_t( 207, 225, 28 ), 
    colour_t( 210, 225, 27 ), 
    colour_t( 212, 225, 26 ), 
    colour_t( 215, 226, 25 ), 
    colour_t( 218, 226, 24 ), 
    colour_t( 220, 226, 24 ), 
    colour_t( 223, 227, 24 ), 
    colour_t( 225, 227, 24 ), 
    colour_t( 228, 227, 24 ), 
    colour_t( 231, 228, 25 ), 
    colour_t( 233, 228, 25 ), 
    colour_t( 236, 228, 26 ), 
    colour_t( 238, 229, 27 ), 
    colour_t( 241, 229, 28 ), 
    colour_t( 243, 229, 30 ), 
    colour_t( 246, 230, 31 ), 
    colour_t( 248, 230, 33 ), 
    colour_t( 250, 230, 34 ), 
    colour_t( 253, 231, 36 )
};

}// namespace anonymous

void
TViridisColourMap::init ( const size_t n )
{
    init_cmap_with_data( n, _map, viridis_data, sizeof(viridis_data) / sizeof(colour_t) );
}

//
// TParulaColourMap
//

namespace
{

colour_t  parula_data[] = {
    colour_t( 53, 42, 134 ), 
    colour_t( 53, 43, 137 ), 
    colour_t( 53, 45, 141 ), 
    colour_t( 53, 46, 144 ), 
    colour_t( 53, 48, 147 ), 
    colour_t( 54, 49, 150 ), 
    colour_t( 54, 51, 153 ), 
    colour_t( 54, 52, 156 ), 
    colour_t( 54, 54, 159 ), 
    colour_t( 54, 55, 162 ), 
    colour_t( 53, 57, 165 ), 
    colour_t( 53, 59, 169 ), 
    colour_t( 53, 60, 172 ), 
    colour_t( 52, 62, 175 ), 
    colour_t( 51, 63, 178 ), 
    colour_t( 51, 65, 181 ), 
    colour_t( 50, 67, 185 ), 
    colour_t( 48, 68, 188 ), 
    colour_t( 47, 70, 191 ), 
    colour_t( 45, 72, 194 ), 
    colour_t( 44, 74, 197 ), 
    colour_t( 41, 75, 201 ), 
    colour_t( 39, 77, 204 ), 
    colour_t( 36, 79, 207 ), 
    colour_t( 33, 82, 210 ), 
    colour_t( 29, 84, 213 ), 
    colour_t( 25, 86, 216 ), 
    colour_t( 20, 88, 218 ), 
    colour_t( 16, 91, 220 ), 
    colour_t( 12, 93, 222 ), 
    colour_t( 8, 94, 223 ), 
    colour_t( 5, 96, 224 ), 
    colour_t( 3, 98, 224 ), 
    colour_t( 2, 99, 225 ), 
    colour_t( 1, 101, 225 ), 
    colour_t( 1, 102, 225 ), 
    colour_t( 1, 103, 225 ), 
    colour_t( 1, 104, 225 ), 
    colour_t( 2, 106, 224 ), 
    colour_t( 2, 107, 224 ), 
    colour_t( 3, 108, 224 ), 
    colour_t( 4, 109, 223 ), 
    colour_t( 5, 110, 223 ), 
    colour_t( 6, 111, 223 ), 
    colour_t( 7, 112, 222 ), 
    colour_t( 8, 113, 222 ), 
    colour_t( 10, 114, 221 ), 
    colour_t( 11, 115, 221 ), 
    colour_t( 12, 116, 220 ), 
    colour_t( 13, 117, 220 ), 
    colour_t( 13, 118, 219 ), 
    colour_t( 14, 119, 219 ), 
    colour_t( 15, 120, 218 ), 
    colour_t( 16, 121, 217 ), 
    colour_t( 16, 122, 217 ), 
    colour_t( 17, 123, 216 ), 
    colour_t( 18, 123, 216 ), 
    colour_t( 18, 124, 215 ), 
    colour_t( 19, 125, 215 ), 
    colour_t( 19, 126, 214 ), 
    colour_t( 19, 127, 214 ), 
    colour_t( 19, 128, 213 ), 
    colour_t( 20, 129, 213 ), 
    colour_t( 20, 130, 212 ), 
    colour_t( 20, 131, 212 ), 
    colour_t( 20, 132, 211 ), 
    colour_t( 20, 133, 211 ), 
    colour_t( 19, 135, 211 ), 
    colour_t( 19, 136, 210 ), 
    colour_t( 19, 137, 210 ), 
    colour_t( 18, 138, 210 ), 
    colour_t( 17, 139, 210 ), 
    colour_t( 17, 140, 210 ), 
    colour_t( 16, 142, 210 ), 
    colour_t( 15, 143, 210 ), 
    colour_t( 14, 144, 209 ), 
    colour_t( 13, 146, 209 ), 
    colour_t( 12, 147, 209 ), 
    colour_t( 11, 148, 209 ), 
    colour_t( 10, 149, 209 ), 
    colour_t( 9, 151, 209 ), 
    colour_t( 8, 152, 209 ), 
    colour_t( 8, 153, 208 ), 
    colour_t( 7, 154, 208 ), 
    colour_t( 7, 155, 207 ), 
    colour_t( 6, 156, 207 ), 
    colour_t( 6, 157, 206 ), 
    colour_t( 6, 158, 206 ), 
    colour_t( 6, 159, 205 ), 
    colour_t( 6, 160, 204 ), 
    colour_t( 6, 161, 204 ), 
    colour_t( 5, 161, 203 ), 
    colour_t( 5, 162, 202 ), 
    colour_t( 5, 163, 201 ), 
    colour_t( 5, 164, 200 ), 
    colour_t( 5, 165, 200 ), 
    colour_t( 5, 165, 199 ), 
    colour_t( 5, 166, 198 ), 
    colour_t( 5, 167, 197 ), 
    colour_t( 6, 167, 196 ), 
    colour_t( 6, 168, 195 ), 
    colour_t( 6, 169, 194 ), 
    colour_t( 7, 169, 193 ), 
    colour_t( 7, 170, 192 ), 
    colour_t( 8, 171, 190 ), 
    colour_t( 9, 171, 189 ), 
    colour_t( 10, 172, 188 ), 
    colour_t( 11, 172, 187 ), 
    colour_t( 13, 173, 186 ), 
    colour_t( 14, 174, 185 ), 
    colour_t( 16, 174, 184 ), 
    colour_t( 17, 175, 182 ), 
    colour_t( 19, 175, 181 ), 
    colour_t( 20, 176, 180 ), 
    colour_t( 22, 177, 179 ), 
    colour_t( 24, 177, 177 ), 
    colour_t( 26, 178, 176 ), 
    colour_t( 28, 178, 175 ), 
    colour_t( 30, 179, 174 ), 
    colour_t( 32, 179, 172 ), 
    colour_t( 34, 180, 171 ), 
    colour_t( 36, 180, 170 ), 
    colour_t( 38, 181, 168 ), 
    colour_t( 40, 181, 167 ), 
    colour_t( 42, 182, 165 ), 
    colour_t( 44, 182, 164 ), 
    colour_t( 47, 183, 163 ), 
    colour_t( 49, 183, 161 ), 
    colour_t( 51, 184, 160 ), 
    colour_t( 54, 184, 158 ), 
    colour_t( 56, 185, 157 ), 
    colour_t( 59, 185, 155 ), 
    colour_t( 61, 185, 154 ), 
    colour_t( 64, 186, 152 ), 
    colour_t( 67, 186, 151 ), 
    colour_t( 69, 187, 149 ), 
    colour_t( 72, 187, 148 ), 
    colour_t( 75, 187, 146 ), 
    colour_t( 78, 188, 145 ), 
    colour_t( 81, 188, 143 ), 
    colour_t( 83, 188, 142 ), 
    colour_t( 86, 189, 140 ), 
    colour_t( 89, 189, 139 ), 
    colour_t( 92, 189, 137 ), 
    colour_t( 95, 189, 136 ), 
    colour_t( 98, 190, 134 ), 
    colour_t( 101, 190, 133 ), 
    colour_t( 104, 190, 131 ), 
    colour_t( 107, 190, 130 ), 
    colour_t( 110, 190, 129 ), 
    colour_t( 113, 190, 128 ), 
    colour_t( 116, 190, 126 ), 
    colour_t( 119, 190, 125 ), 
    colour_t( 121, 190, 124 ), 
    colour_t( 124, 191, 123 ), 
    colour_t( 127, 191, 122 ), 
    colour_t( 130, 191, 120 ), 
    colour_t( 132, 191, 119 ), 
    colour_t( 135, 191, 118 ), 
    colour_t( 138, 190, 117 ), 
    colour_t( 140, 190, 116 ), 
    colour_t( 143, 190, 115 ), 
    colour_t( 145, 190, 114 ), 
    colour_t( 148, 190, 113 ), 
    colour_t( 150, 190, 112 ), 
    colour_t( 153, 190, 111 ), 
    colour_t( 155, 190, 110 ), 
    colour_t( 158, 190, 109 ), 
    colour_t( 160, 190, 108 ), 
    colour_t( 162, 190, 107 ), 
    colour_t( 165, 190, 106 ), 
    colour_t( 167, 190, 105 ), 
    colour_t( 169, 189, 104 ), 
    colour_t( 171, 189, 104 ), 
    colour_t( 174, 189, 103 ), 
    colour_t( 176, 189, 102 ), 
    colour_t( 178, 189, 101 ), 
    colour_t( 180, 189, 100 ), 
    colour_t( 182, 189, 99 ), 
    colour_t( 185, 188, 98 ), 
    colour_t( 187, 188, 97 ), 
    colour_t( 189, 188, 97 ), 
    colour_t( 191, 188, 96 ), 
    colour_t( 193, 188, 95 ), 
    colour_t( 195, 187, 94 ), 
    colour_t( 197, 187, 93 ), 
    colour_t( 199, 187, 92 ), 
    colour_t( 202, 187, 91 ), 
    colour_t( 204, 187, 91 ), 
    colour_t( 206, 187, 90 ), 
    colour_t( 208, 186, 89 ), 
    colour_t( 210, 186, 88 ), 
    colour_t( 212, 186, 87 ), 
    colour_t( 214, 186, 86 ), 
    colour_t( 216, 186, 85 ), 
    colour_t( 218, 185, 85 ), 
    colour_t( 220, 185, 84 ), 
    colour_t( 222, 185, 83 ), 
    colour_t( 224, 185, 82 ), 
    colour_t( 226, 185, 81 ), 
    colour_t( 228, 185, 80 ), 
    colour_t( 230, 185, 79 ), 
    colour_t( 232, 185, 78 ), 
    colour_t( 234, 185, 77 ), 
    colour_t( 236, 185, 76 ), 
    colour_t( 238, 185, 75 ), 
    colour_t( 240, 185, 74 ), 
    colour_t( 242, 185, 72 ), 
    colour_t( 243, 185, 71 ), 
    colour_t( 245, 185, 70 ), 
    colour_t( 247, 186, 68 ), 
    colour_t( 249, 186, 67 ), 
    colour_t( 250, 187, 65 ), 
    colour_t( 251, 188, 63 ), 
    colour_t( 253, 189, 62 ), 
    colour_t( 253, 190, 60 ), 
    colour_t( 254, 191, 58 ), 
    colour_t( 254, 193, 57 ), 
    colour_t( 254, 194, 55 ), 
    colour_t( 254, 195, 54 ), 
    colour_t( 254, 197, 53 ), 
    colour_t( 254, 198, 52 ), 
    colour_t( 254, 199, 50 ), 
    colour_t( 253, 200, 49 ), 
    colour_t( 253, 202, 48 ), 
    colour_t( 252, 203, 47 ), 
    colour_t( 252, 204, 46 ), 
    colour_t( 251, 206, 45 ), 
    colour_t( 251, 207, 44 ), 
    colour_t( 250, 208, 43 ), 
    colour_t( 250, 209, 42 ), 
    colour_t( 249, 211, 41 ), 
    colour_t( 248, 212, 40 ), 
    colour_t( 248, 213, 39 ), 
    colour_t( 247, 215, 38 ), 
    colour_t( 247, 216, 37 ), 
    colour_t( 246, 217, 36 ), 
    colour_t( 246, 219, 35 ), 
    colour_t( 245, 220, 34 ), 
    colour_t( 245, 222, 33 ), 
    colour_t( 245, 223, 32 ), 
    colour_t( 244, 225, 30 ), 
    colour_t( 244, 226, 29 ), 
    colour_t( 244, 228, 28 ), 
    colour_t( 244, 230, 27 ), 
    colour_t( 244, 231, 26 ), 
    colour_t( 244, 233, 25 ), 
    colour_t( 244, 235, 24 ), 
    colour_t( 245, 237, 22 ), 
    colour_t( 245, 238, 21 ), 
    colour_t( 245, 240, 20 ), 
    colour_t( 246, 242, 19 ), 
    colour_t( 246, 244, 17 ), 
    colour_t( 247, 246, 16 ), 
    colour_t( 248, 248, 15 ), 
    colour_t( 248, 250, 13 )
};

}// namespace anonymous

void
TParulaColourMap::init ( const size_t n )
{
    init_cmap_with_data( n, _map, parula_data, sizeof(parula_data) / sizeof(colour_t) );
}

//
// RandomColourMap
//
void
TRandomColourMap::init ( const size_t n )
{
    const size_t  primes[] = { 5, 7, 11, 13, 17, 19, 23, 29 }; // enough as long as n < prod(primes)
    size_t        prime    = 3;
    const size_t  base_n   = n;

    // look for small prime, not a factor of n
    for ( size_t  i = 0; i < sizeof(primes)/sizeof(size_t); ++i )
    {
        if ( (n / primes[i]) * primes[i] != n )
        {
            prime = primes[i];
            break;
        }// if
    }// for
    
    _base_cmap.init( base_n );

    _map.resize( n );

    for ( size_t  i = 0; i < n; ++i )
    {
        _map[i] = _base_cmap.ientry( (i * prime) % n );
    }// for

    // std::random_shuffle( _map.begin(), _map.end() );
}

//
// ShuffleColourMap
//
void
TShuffleColourMap::init ( const size_t n )
{
    const size_t  base_n = n;
    const size_t  cycle  = 4;

    _base_cmap.init( base_n );

    _map.resize( n );

    for ( size_t  i = 0; i < n; ++i )
    {
        const size_t  pos = ((i*cycle) % n) + (i*cycle)/n;

        _map[i] = _base_cmap.ientry( pos );
    }// for
}

}// namespace Hpro
