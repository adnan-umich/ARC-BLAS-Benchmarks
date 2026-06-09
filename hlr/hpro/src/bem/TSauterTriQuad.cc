//
// Project     : HLIBpro
// File        : TSauterTriQuad.cc
// Description : provides quadrature rules for 2 triangles
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/packed.hh"
#include "hpro/base/error.hh"
#include "hpro/bem/TGaussQuad.hh"
#include "hpro/parallel/TMutex.hh"

#include "hpro/bem/TSauterTriQuad.hh"

namespace Hpro
{

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// local data
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

//
// caches for quadrature points with mutices
// (with this, costly recomputation of quadrature
//  points is prevented)
//

namespace
{

std::vector< std::vector< T2Point > * >  cache_eq_pts1, cache_eq_pts2;
std::vector< std::vector< double > >     cache_eq_wghts;
std::vector< std::vector< T2Point > * >  cache_ed_pts1, cache_ed_pts2;
std::vector< std::vector< double > >     cache_ed_wghts;
std::vector< std::vector< T2Point > * >  cache_vt_pts1, cache_vt_pts2;
std::vector< std::vector< double > >     cache_vt_wghts;
std::vector< std::vector< T2Point > * >  cache_ne_pts1, cache_ne_pts2;
std::vector< std::vector< double > >     cache_ne_wghts;

TMutex  cache_eq_mutex;
TMutex  cache_ed_mutex;
TMutex  cache_vt_mutex;
TMutex  cache_ne_mutex;

}// namespace

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// Gauss quadrature rule for triangle pairs developed
// by Stefan Sauter
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////
//
// constructor and destructor
//

TSauterTriQuad::TSauterTriQuad ()
        : _order_xi( 0 )
{
}

TSauterTriQuad::~TSauterTriQuad ()
{
}

//
// construct quadrature points and weights
//
void
TSauterTriQuad::build ( const uint                ncommon,
                        const uint                order,
                        std::vector< T2Point > &  t1,
                        std::vector< T2Point > &  t2,
                        std::vector< double > &   w ) const
{
    switch ( ncommon )
    {
    case 3  : quad_equ(  order, t1, t2, w ); break;
    case 2  : quad_edge( order, t1, t2, w ); break;
    case 1  : quad_vtx(  order, t1, t2, w ); break;
    case 0  : quad_dist( order, t1, t2, w ); break;
    default :
        HERROR( ERR_ARG, "(TSauterTriQuad) build", "unkown triangle type" );
    }// switch
}

void
TSauterTriQuad::build ( const uint               ncommon,
                        const uint               order,
                        std::vector< rule_t > &  rule ) const
{
    std::vector< T2Point >  tri1pts, tri2pts;
    std::vector< double >   weights;

    build( ncommon, order, tri1pts, tri2pts, weights );

    rule.resize( tri1pts.size() );

    for ( uint i = 0; i < tri1pts.size(); i++ )
    {
        rule[i].pttri1 = tri1pts[i];
        rule[i].pttri2 = tri2pts[i];
        rule[i].wght   = weights[i];
    }// for
}

//
// return quadrature points for given triangle pair
//
uint
TSauterTriQuad::points ( const uint                       order,
                         const std::vector< double * > &  tri1,
                         const std::vector< double * > &  tri2,
                         std::vector< T3Point > &         tri1_pts,
                         std::vector< T3Point > &         tri2_pts,
                         std::vector< double > &          weights ) const
{
    //
    // count number of common vertices and
    // order common vertices to the same position
    //
    
    uint  common_vertices = 0;
    uint  perm1[3] = { 0, 0, 0 };  // describes ordered vertices of triangle1
    uint  perm2[3] = { 0, 0, 0 };  // and triangle2

    for ( uint i = 0; i < 3; i++ )
        for ( uint j = 0; j < 3; j++ )
            if ( tri1[i] == tri2[j] )
            {
                perm1[common_vertices] = i;
                perm2[common_vertices] = j;
                common_vertices++;
                break;
            }// if
    
    if ( common_vertices > 3 )
        HERROR( ERR_CONSISTENCY, "(TSauterTriQuad) points",
                "given triangles have more than 3 common vertices" );
    
    //
    // assign remaining vertices
    //
    
    uint  other_vertices = common_vertices;
    
    for ( uint i = 0; i < 3 && other_vertices < 3; i++ )
    {
        uint  j;
        
        for ( j = 0; j < other_vertices; j++ )
            if ( tri1[i] == tri1[perm1[j]] )
                break;
        
        if ( j == other_vertices )
            perm1[other_vertices++] = i;
    }// for
    
    other_vertices = common_vertices;
    
    for ( uint i = 0; i < 3 && other_vertices < 3; i++ )
    {
        uint  j;
        
        for ( j = 0; j < other_vertices; j++ )
            if ( tri2[i] == tri2[perm2[j]] )
                break;
        
        if ( j == other_vertices )
            perm2[other_vertices++] = i;
    }// for

    //
    // compute quadrature points based on detected situation
    //

    std::vector< T2Point > * quad1_pts;
    std::vector< T2Point > * quad2_pts;
    
    switch ( common_vertices )
    {
    case 3 :
        cache_eq_mutex.lock();
        if (( cache_eq_pts1.size() > order ) && ( cache_eq_pts1[order]->size() > 0 ))
        {
            quad1_pts = cache_eq_pts1[order];
            quad2_pts = cache_eq_pts2[order];
        }// if
        else
        {
            quad1_pts = new std::vector< T2Point >;
            quad2_pts = new std::vector< T2Point >;
                
            cache_eq_pts1.resize( order+1 );
            cache_eq_pts2.resize( order+1 );
            cache_eq_wghts.resize( order+1 );

            quad_equ( order, * quad1_pts, * quad2_pts, cache_eq_wghts[order] );

            cache_eq_pts1[order] = quad1_pts;
            cache_eq_pts2[order] = quad2_pts;
        }// else
        cache_eq_mutex.unlock();

        weights.resize( cache_eq_wghts[order].size() );
            
        for ( uint i = 0; i < cache_eq_wghts[order].size(); i++ )
            weights[i] = cache_eq_wghts[order][i];
        
        break;
        
    case 2 :
        cache_ed_mutex.lock();
        if (( cache_ed_pts1.size() > order ) && ( cache_ed_pts1[order]->size() > 0 ))
        {
            quad1_pts = cache_ed_pts1[order];
            quad2_pts = cache_ed_pts2[order];
        }// if
        else
        {
            quad1_pts = new std::vector< T2Point >;
            quad2_pts = new std::vector< T2Point >;
                
            cache_ed_pts1.resize( order+1 );
            cache_ed_pts2.resize( order+1 );
            cache_ed_wghts.resize( order+1 );

            quad_edge( order, * quad1_pts, * quad2_pts, cache_ed_wghts[order] );

            cache_ed_pts1[order] = quad1_pts;
            cache_ed_pts2[order] = quad2_pts;
        }// else
        cache_ed_mutex.unlock();

        weights.resize( cache_ed_wghts[order].size() );
        
        for ( uint i = 0; i < cache_ed_wghts[order].size(); i++ )
            weights[i] = cache_ed_wghts[order][i];
        
        break;
        
    case 1 :
        cache_vt_mutex.lock();
        if (( cache_vt_pts1.size() > order ) && ( cache_vt_pts1[order]->size() > 0 ))
        {
            quad1_pts = cache_vt_pts1[order];
            quad2_pts = cache_vt_pts2[order];
        }// if
        else
        {
            quad1_pts = new std::vector< T2Point >;
            quad2_pts = new std::vector< T2Point >;
                
            cache_vt_pts1.resize( order+1 );
            cache_vt_pts2.resize( order+1 );
            cache_vt_wghts.resize( order+1 );

            quad_vtx(  order, * quad1_pts, * quad2_pts, cache_vt_wghts[order] );

            cache_vt_pts1[order] = quad1_pts;
            cache_vt_pts2[order] = quad2_pts;
        }// else
        cache_vt_mutex.unlock();
        
        weights.resize( cache_vt_wghts[order].size() );
        
        for ( uint i = 0; i < cache_vt_wghts[order].size(); i++ )
            weights[i] = cache_vt_wghts[order][i];
        
        break;
        
    case 0 :
    default:
        cache_ne_mutex.lock();
        if (( cache_ne_pts1.size() > order ) && ( cache_ne_pts1[order]->size() > 0 ))
        {
            quad1_pts = cache_ne_pts1[order];
            quad2_pts = cache_ne_pts2[order];
        }// if
        else
        {
            quad1_pts = new std::vector< T2Point >;
            quad2_pts = new std::vector< T2Point >;
                
            cache_ne_pts1.resize( order+1 );
            cache_ne_pts2.resize( order+1 );
            cache_ne_wghts.resize( order+1 );

            quad_dist( order, * quad1_pts, * quad2_pts, cache_ne_wghts[order] );

            cache_ne_pts1[order] = quad1_pts;
            cache_ne_pts2[order] = quad2_pts;
        }// else
        cache_ne_mutex.unlock();

        weights.resize( cache_ne_wghts[order].size() );
            
        for ( uint i = 0; i < cache_ne_wghts[order].size(); i++ )
            weights[i] = cache_ne_wghts[order][i];
        
        break;
    }// else

    HASSERT( quad1_pts->size() == quad2_pts->size(),
             ERR_CONSISTENCY, "(TSauterTriQuad) points",
             "different number of quadrature points per triangle" );
    
    //
    // transform quadrature points to local triangle coordinates
    //

    const size_t  npts = quad1_pts->size();

    tri1_pts.resize( npts );
    tri2_pts.resize( npts );

#if defined(__SSE2__)

    using  packed_t = packed< double, ISA_SSE2 >;
    
    const packed_t  one( 1.0 );
    packed_t        tmp0, tmp1;
    
    // load triangle coordinates
    const packed_t  tri10[3] = { tri1[perm1[0]][0], tri1[perm1[0]][1], tri1[perm1[0]][2] };
    const packed_t  tri11[3] = { tri1[perm1[1]][0], tri1[perm1[1]][1], tri1[perm1[1]][2] };
    const packed_t  tri12[3] = { tri1[perm1[2]][0], tri1[perm1[2]][1], tri1[perm1[2]][2] };

    const packed_t  tri20[3] = { tri2[perm2[0]][0], tri2[perm2[0]][1], tri2[perm2[0]][2] };
    const packed_t  tri21[3] = { tri2[perm2[1]][0], tri2[perm2[1]][1], tri2[perm2[1]][2] };
    const packed_t  tri22[3] = { tri2[perm2[2]][0], tri2[perm2[2]][1], tri2[perm2[2]][2] };

    for ( size_t  i = 0; i+1 < npts; i +=2 )
    {
        double    c[2];
        
        // load quadrature points for first triangle
        tmp0 = load< packed_t >( (*quad1_pts)[i].vector() );
        tmp1 = load< packed_t >( (*quad1_pts)[i+1].vector() );

        const packed_t  a1 = unpacklo( tmp0, tmp1 );
        const packed_t  a2 = unpackhi( tmp0, tmp1 );
        const packed_t  a0 = sub( sub( one, a1 ),  a2 );

        // load quadrature points for second triangle
        tmp0 = load< packed_t >( (*quad2_pts)[i].vector() );
        tmp1 = load< packed_t >( (*quad2_pts)[i+1].vector() );

        const packed_t  b1 = unpacklo( tmp0, tmp1 );
        const packed_t  b2 = unpackhi( tmp0, tmp1 );
        const packed_t  b0 = sub( sub( one, b1 ), b2 );

        // transform to triangle coordinates
        store( add( add( mul( a0, tri10[0] ), mul( a1, tri11[0] ) ), mul( a2, tri12[0] ) ), c );
        tri1_pts[i][0]   = c[0];
        tri1_pts[i+1][0] = c[1];

        store( add( add( mul( b0, tri20[0] ), mul( b1, tri21[0] ) ), mul( b2, tri22[0] ) ), c );
        tri2_pts[i][0]   = c[0];
        tri2_pts[i+1][0] = c[1];

        store( add( add( mul( a0, tri10[1] ), mul( a1, tri11[1] ) ), mul( a2, tri12[1] ) ), c );
        tri1_pts[i][1]   = c[0];
        tri1_pts[i+1][1] = c[1];
        
        store( add( add( mul( b0, tri20[1] ), mul( b1, tri21[1] ) ), mul( b2, tri22[1] ) ), c );
        tri2_pts[i][1]   = c[0];
        tri2_pts[i+1][1] = c[1];

        store( add( add( mul( a0, tri10[2] ), mul( a1, tri11[2] ) ), mul( a2, tri12[2] ) ), c );
        tri1_pts[i][2]   = c[0];
        tri1_pts[i+1][2] = c[1];
        
        store( add( add( mul( b0, tri20[2] ), mul( b1, tri21[2] ) ), mul( b2, tri22[2] ) ), c );
        tri2_pts[i][2]   = c[0];
        tri2_pts[i+1][2] = c[1];
    }// for

    if ( npts % 2 == 1 )
    {
        const size_t  i  = npts-1;
        
        const double  a1 = (*quad1_pts)[i][0];
        const double  a2 = (*quad1_pts)[i][1];
        const double  a0 = 1.0 - a1 - a2;
        
        const double  b1 = (*quad2_pts)[i][0];
        const double  b2 = (*quad2_pts)[i][1];
        const double  b0 = 1.0 - b1 - b2;
        
        tri1_pts[i][0] = a0 * tri1[perm1[0]][0] + a1 * tri1[perm1[1]][0] + a2 * tri1[perm1[2]][0];
        tri1_pts[i][1] = a0 * tri1[perm1[0]][1] + a1 * tri1[perm1[1]][1] + a2 * tri1[perm1[2]][1];
        tri1_pts[i][2] = a0 * tri1[perm1[0]][2] + a1 * tri1[perm1[1]][2] + a2 * tri1[perm1[2]][2];

        tri2_pts[i][0] = b0 * tri2[perm2[0]][0] + b1 * tri2[perm2[1]][0] + b2 * tri2[perm2[2]][0];
        tri2_pts[i][1] = b0 * tri2[perm2[0]][1] + b1 * tri2[perm2[1]][1] + b2 * tri2[perm2[2]][1];
        tri2_pts[i][2] = b0 * tri2[perm2[0]][2] + b1 * tri2[perm2[1]][2] + b2 * tri2[perm2[2]][2];
    }// for
    
#else
    
    for ( size_t  i = 0; i < npts; i++ )
    {
        const double  a1 = (*quad1_pts)[i][0];
        const double  a2 = (*quad1_pts)[i][1];
        const double  a0 = 1.0 - a1 - a2;
        
        const double  b1 = (*quad2_pts)[i][0];
        const double  b2 = (*quad2_pts)[i][1];
        const double  b0 = 1.0 - b1 - b2;
        
        tri1_pts[i][0] = a0 * tri1[perm1[0]][0] + a1 * tri1[perm1[1]][0] + a2 * tri1[perm1[2]][0];
        tri1_pts[i][1] = a0 * tri1[perm1[0]][1] + a1 * tri1[perm1[1]][1] + a2 * tri1[perm1[2]][1];
        tri1_pts[i][2] = a0 * tri1[perm1[0]][2] + a1 * tri1[perm1[1]][2] + a2 * tri1[perm1[2]][2];

        tri2_pts[i][0] = b0 * tri2[perm2[0]][0] + b1 * tri2[perm2[1]][0] + b2 * tri2[perm2[2]][0];
        tri2_pts[i][1] = b0 * tri2[perm2[0]][1] + b1 * tri2[perm2[1]][1] + b2 * tri2[perm2[2]][1];
        tri2_pts[i][2] = b0 * tri2[perm2[0]][2] + b1 * tri2[perm2[1]][2] + b2 * tri2[perm2[2]][2];
    }// for

#endif // __SSE2__

    return uint( npts );
}
    
///////////////////////////////////////////////////////////////////
//
// quadrature rules for various cases (see Sauter/Schwab)
// - equ  : triangles are equal
// - edge : triangles have common edge
// - vtx  : triangles have common vertex
// - dist : triangles have positive distance
//
void
TSauterTriQuad::quad_equ ( const uint                order,
                           std::vector< T2Point > &  t1,
                           std::vector< T2Point > &  t2,
                           std::vector< double > &   w ) const
{
    uint                   nn, q = 0;
    const uint             loc_order_xi = ( _order_xi == 0 ? order : Math::min( _order_xi, order ) );
    std::vector< double >  x1( order ), w1( order );
    std::vector< double >  x1_xi( loc_order_xi ), w1_xi( loc_order_xi );
    TGaussQuad             quad;

    quad.build( order, x1, w1 );
    quad.build( loc_order_xi, x1_xi, w1_xi );

    nn = 6 * loc_order_xi * order * order * order;

    t1.resize( nn );
    t2.resize( nn );
    w.resize( nn );

    for ( uint xi = 0; xi < loc_order_xi; xi++ )
    {
        const double wxi = w1_xi[xi];
        const double xxi = x1_xi[xi];
        
        for ( uint eta3 = 0; eta3 < order; eta3++ )
        {
            const double weta3 = wxi * w1[eta3];
            const double xeta3 = x1[eta3];
            
            for ( uint eta2 = 0; eta2 < order; eta2++ )
            {
                const double weta2 = weta3 * w1[eta2];
                const double xeta2 = x1[eta2];
                
                for ( uint eta1 = 0; eta1 < order; eta1++ )
                {
                    const double weta1 = weta2 * w1[eta1];
                    const double xeta1 = x1[eta1];
                    const double lw    = weta1 * xxi * xxi * xxi * xeta1 * xeta1 * xeta2;

                    if ( q + 5 >= nn )
                        HERROR( ERR_SIZE, "(TSauterTriQuad) quad_equ", "array index too large" );

                    t1[q][0]  = xxi;
                    t1[q][1]  = xxi * (1.0 - xeta1 + xeta1 * xeta2);
                    t2[q][0]  = xxi * (1.0 - xeta1 * xeta2 * xeta3);
                    t2[q][1]  = xxi * (1.0 - xeta1);
                    w[q]      = lw;
                    t1[q][0] -= t1[q][1];            // transform to usual unit triangle
                    t2[q][0] -= t2[q][1];
                    q++;

                    t1[q][0]  = xxi * (1.0 - xeta1 * xeta2 * xeta3);
                    t1[q][1]  = xxi * (1.0 - xeta1);
                    t2[q][0]  = xxi;
                    t2[q][1]  = xxi * (1.0 - xeta1 + xeta1 * xeta2);
                    w[q]      = lw;
                    t1[q][0] -= t1[q][1];
                    t2[q][0] -= t2[q][1];
                    q++;

                    t1[q][0]  = xxi;
                    t1[q][1]  = xxi * (xeta1 * (1.0 - xeta2 + xeta2 * xeta3));
                    t2[q][0]  = xxi * (1.0 - xeta1 * xeta2);
                    t2[q][1]  = xxi * (xeta1 * (1.0 - xeta2));
                    w[q]      = lw;
                    t1[q][0] -= t1[q][1];
                    t2[q][0] -= t2[q][1];
                    q++;

                    t1[q][0]  = xxi * (1.0 - xeta1 * xeta2);
                    t1[q][1]  = xxi * (xeta1 * (1.0 - xeta2));
                    t2[q][0]  = xxi;
                    t2[q][1]  = xxi * (xeta1 * (1.0 - xeta2 + xeta2 * xeta3));
                    w[q]      = lw;
                    t1[q][0] -= t1[q][1];
                    t2[q][0] -= t2[q][1];
                    q++;
          
                    t1[q][0]  = xxi * (1.0 - xeta1 * xeta2 * xeta3);
                    t1[q][1]  = xxi * (xeta1 * (1.0 - xeta2 * xeta3));
                    t2[q][0]  = xxi;
                    t2[q][1]  = xxi * (xeta1 * (1.0 - xeta2));
                    w[q]      = lw;
                    t1[q][0] -= t1[q][1];
                    t2[q][0] -= t2[q][1];
                    q++;

                    t1[q][0] = xxi;
                    t1[q][1] = xxi * (xeta1 * (1.0 - xeta2));
                    t2[q][0] = xxi * (1.0 - xeta1 * xeta2 * xeta3);
                    t2[q][1] = xxi * (xeta1 * (1.0 - xeta2 * xeta3));
                    w[q]     = lw;
                    t1[q][0] -= t1[q][1];
                    t2[q][0] -= t2[q][1];
                    q++;
                }// for
            }// for
        }// for
    }// for
}

void
TSauterTriQuad::quad_edge ( const uint                order,
                            std::vector< T2Point > &  t1,
                            std::vector< T2Point > &  t2,
                            std::vector< double > &   w ) const
{
    uint                   nn, q = 0;
    const uint             loc_order_xi = ( _order_xi == 0 ? order : Math::min( _order_xi, order ) );
    std::vector< double >  x1( order ), w1( order );
    std::vector< double >  x1_xi( loc_order_xi ), w1_xi( loc_order_xi );
    TGaussQuad             quad;

    quad.build( order, x1, w1 );
    quad.build( loc_order_xi, x1_xi, w1_xi );

    nn = 5 * loc_order_xi * order * order * order;

    t1.resize( nn );
    t2.resize( nn );
    w.resize( nn );

    for ( uint xi = 0; xi < loc_order_xi; xi++ )
    {
        const double wxi = w1_xi[xi];
        const double xxi = x1_xi[xi];
        
        for ( uint eta3 = 0; eta3 < order; eta3++ )
        {
            const double weta3 = wxi * w1[eta3];
            const double xeta3 = x1[eta3];
            
            for ( uint eta2 = 0; eta2 < order; eta2++ )
            {
                const double weta2 = weta3 * w1[eta2];
                const double xeta2 = x1[eta2];
                
                for ( uint eta1 = 0; eta1 < order; eta1++ )
                {
                    const double weta1 = weta2 * w1[eta1];
                    const double xeta1 = x1[eta1];
                    const double lw    = weta1 * xxi * xxi * xxi * xeta1 * xeta1 * xeta2;
                    
                    if ( q + 4 >= nn )
                        HERROR( ERR_SIZE, "(TSauterTriQuad) quad_edge", "array index too large" );

                    t1[q][0]  = xxi;
                    t1[q][1]  = xxi * (xeta1 * xeta3);
                    t2[q][0]  = xxi * (1.0 - xeta1 * xeta2);
                    t2[q][1]  = xxi * (xeta1 * (1.0 - xeta2));
                    w[q]      = weta1 * xxi * xxi * xxi * xeta1 * xeta1;
                    t1[q][0] -= t1[q][1];
                    t2[q][0] -= t2[q][1];
                    q++;

                    t1[q][0]  = xxi;
                    t1[q][1]  = xxi * (xeta1);
                    t2[q][0]  = xxi * (1.0 - xeta1 * xeta2 * xeta3);
                    t2[q][1]  = xxi * (xeta1 * xeta2 * (1.0 - xeta3));
                    w[q]      = lw;
                    t1[q][0] -= t1[q][1];
                    t2[q][0] -= t2[q][1];
                    q++;

                    t1[q][0]  = xxi * (1.0 - xeta1 * xeta2);
                    t1[q][1]  = xxi * (xeta1 * (1.0 - xeta2));
                    t2[q][0]  = xxi;
                    t2[q][1]  = xxi * (xeta1 * xeta2 * xeta3);
                    w[q]      = lw;
                    t1[q][0] -= t1[q][1];
                    t2[q][0] -= t2[q][1];
                    q++;

                    t1[q][0]  = xxi * (1.0 - xeta1 * xeta2 * xeta3);
                    t1[q][1]  = xxi * (xeta1 * xeta2 * (1.0 - xeta3));
                    t2[q][0]  = xxi;
                    t2[q][1]  = xxi * (xeta1);
                    w[q]      = lw;
                    t1[q][0] -= t1[q][1];
                    t2[q][0] -= t2[q][1];
                    q++;

                    t1[q][0]  = xxi * (1.0 - xeta1 * xeta2 * xeta3);
                    t1[q][1]  = xxi * (xeta1 * (1.0 - xeta2 * xeta3));
                    t2[q][0]  = xxi;
                    t2[q][1]  = xxi * (xeta1 * xeta2);
                    w[q]      = lw;
                    t1[q][0] -= t1[q][1];
                    t2[q][0] -= t2[q][1];
                    q++;
                }// for
            }// for
        }// for
    }// for
}

void
TSauterTriQuad::quad_vtx ( const uint                order,
                           std::vector< T2Point > &  t1,
                           std::vector< T2Point > &  t2,
                           std::vector< double > &   w ) const
{
    uint                   nn, q = 0;
    const uint             loc_order_xi = ( _order_xi == 0 ? order : Math::min( _order_xi, order ) );
    std::vector< double >  x1( order ), w1( order );
    std::vector< double >  x1_xi( loc_order_xi ), w1_xi( loc_order_xi );
    TGaussQuad             quad;

    quad.build( order, x1, w1 );
    quad.build( loc_order_xi, x1_xi, w1_xi );

    nn = 2 * loc_order_xi * order * order * order;

    t1.resize( nn );
    t2.resize( nn );
    w.resize( nn );

    for ( uint xi = 0; xi < loc_order_xi; xi++ )
    {
        const double wxi = w1_xi[xi];
        const double xxi = x1_xi[xi];
        
        for ( uint eta3 = 0; eta3 < order; eta3++ )
        {
            const double weta3 = wxi * w1[eta3];
            const double xeta3 = x1[eta3];
            
            for ( uint eta2 = 0; eta2 < order; eta2++ )
            {
                const double weta2 = weta3 * w1[eta2];
                const double xeta2 = x1[eta2];
                
                for ( uint eta1 = 0; eta1 < order; eta1++ )
                {
                    const double weta1 = weta2 * w1[eta1];
                    const double xeta1 = x1[eta1];
                    const double lw    = weta1 * xxi * xxi * xxi * xeta2;
                    
                    if ( q + 1 >= nn )
                        HERROR( ERR_SIZE, "(TSauterTriQuad) quad_vtx", "array index too large" );

                    t1[q][0]  = xxi;
                    t1[q][1]  = xxi * xeta1;
                    t2[q][0]  = xxi * xeta2;
                    t2[q][1]  = xxi * xeta2 * xeta3;
                    w[q]      = lw;
                    t1[q][0] -= t1[q][1];
                    t2[q][0] -= t2[q][1];
                    q++;
                    
                    t1[q][0]  = xxi * xeta2;
                    t1[q][1]  = xxi * xeta2 * xeta3;
                    t2[q][0]  = xxi;
                    t2[q][1]  = xxi * xeta1;
                    w[q]      = lw;
                    t1[q][0] -= t1[q][1];
                    t2[q][0] -= t2[q][1];
                    q++;
                }// for
            }// for
        }// for
    }// for
}

void
TSauterTriQuad::quad_dist ( const uint                order,
                            std::vector< T2Point > &  t1,
                            std::vector< T2Point > &  t2,
                            std::vector< double > &   w ) const
{
    uint                   nn, q = 0;
    const uint             loc_order_xi = ( _order_xi == 0 ? order : Math::min( _order_xi, order ) );
    std::vector< double >  x1( order ), w1( order );
    std::vector< double >  x1_xi( loc_order_xi ), w1_xi( loc_order_xi );
    TGaussQuad             quad;

    quad.build( order, x1, w1 );
    quad.build( loc_order_xi, x1_xi, w1_xi );

    nn = loc_order_xi * order * order * order;

    t1.resize( nn );
    t2.resize( nn );
    w.resize( nn );

    for ( uint xi = 0; xi < loc_order_xi; xi++ )
    {
        const double wxi = w1_xi[xi];
        const double xxi = x1_xi[xi];
        
        for ( uint eta3 = 0; eta3 < order; eta3++ )
        {
            const double weta3 = wxi * w1[eta3];
            const double xeta3 = x1[eta3];
            
            for ( uint eta2 = 0; eta2 < order; eta2++ )
            {
                const double weta2 = weta3 * w1[eta2];
                const double xeta2 = x1[eta2];
                
                for ( uint eta1 = 0; eta1 < order; eta1++ )
                {
                    const double weta1 = weta2 * w1[eta1];
                    const double xeta1 = x1[eta1];
                    
                    if ( q >= nn )
                        HERROR( ERR_SIZE, "(TSauterTriQuad) quad_dist", "array index too large" );

                    t1[q][0]  = xxi;
                    t1[q][1]  = xxi * xeta1;
                    t2[q][0]  = xeta2;
                    t2[q][1]  = xeta2 * xeta3;
                    w[q]      = weta1 * xxi * xeta2;
                    t1[q][0] -= t1[q][1];
                    t2[q][0] -= t2[q][1];
                    q++;
                }// for
            }// for
        }// for
    }// for
}

}// namespace
