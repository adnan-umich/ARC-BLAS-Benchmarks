//
// Project     : HLIBpro
// File        : TGaussQuad.cc
// Description : provides Gauss quadrature rules
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/blas/Algebra.hh"
#include "hpro/parallel/TMutex.hh"
#include "hpro/base/packed.hh"

#include "hpro/bem/TGaussQuad.hh"

namespace Hpro
{

namespace B = BLAS;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TGaussQuad
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////
//
// construct quadrature points and weights
// in [0,1]
//
void
TGaussQuad::build ( const uint               order,
                    std::vector< double > &  pos,
                    std::vector< double > &  wght ) const
{
    pos.resize( order );
    wght.resize( order );
    
    if ( order == 0 )
        return;
    
    B::Vector< double >  d( order );
    B::Vector< double >  e( order-1 );
    B::Vector< double >  eval( order );
    B::Matrix< double >  evec( order, order );

    for( uint i = 0; i < order; i++ )
        d(i) = 0.0;
    
    for( uint i = 0; i < order-1; i++ )
        e(i) = (double(i) + 1.0) / Math::sqrt( (2.0 * double(i) + 1.0) * (2.0 * double(i) + 3.0) );

    B::eigen_tri( d, e, eval, evec );
    
    for( uint i = 0; i < order; i++ )
    {
        pos[i]  = 0.5 * (eval(i) + 1.0);
        wght[i] = Math::square( evec( 0, i ) );
    }// for
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TTriGaussQuad
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//
// cache for quadrature points and weights
//
namespace
{

std::vector< std::vector< T2Point > * >  cache_points;
std::vector< std::vector< double > >     cache_weights;
TMutex                                   cache_mutex;

}// namespace anonymous

//
// build quadrature points for 2D simplex
//
void
TTriGaussQuad::build ( const uint                order,
                       std::vector< T2Point > &  pts,
                       std::vector< double > &   wghts ) const
{
    TGaussQuad             gauss1d;
    std::vector< double >  x1d, w1d;

    gauss1d.build( order, x1d, w1d );
        
    pts.resize( order * order );
    wghts.resize( order * order );

    for ( uint i = 0; i < order; i++ )
    {
        for ( uint j = 0; j < order; j++ )
        {
            wghts[i+j*order]  = w1d[i] * w1d[j] * x1d[i];
            pts[i+j*order][0] = x1d[i] * ( 1.0 - x1d[j] );
            pts[i+j*order][1] = x1d[i] * x1d[j];
        }// for
    }// for
}

//
// construct quadrature points (<x>) and
// weights (<w>) for triangle <tri> of
// order <order>
//
void
TTriGaussQuad::build ( const uint                      order,
                       const std::vector< T3Point > &  tri,
                       std::vector< T3Point > &        x,
                       std::vector< double > &         w ) const
{
    //
    // compute or just copy quadrature points
    //

    TScopedLock               lock( cache_mutex );
    std::vector< T2Point > *  points;
    
    if (( cache_points.size() > order ) && ( cache_points[order]->size() > 0 ))
    {
        points = cache_points[order];
        w      = cache_weights[order];
    }// if
    else
    {
        points = new std::vector< T2Point >( order * order );
        w.resize( order * order );

        build( order, * points, w );

        cache_points.resize(  order + 1 );
        cache_weights.resize( order + 1 );

        cache_points[order]  = points;
        cache_weights[order] = w;
    }// else

    //
    // transform quadrature points to triangle coordinates
    //

    const size_t  npts = points->size();

    x.resize( npts );

#if defined(__SSE2__)

    using  packed_t = packed< double, ISA_SSE2 >;
    
    const packed_t  one( 1.0 );
    packed_t        tmp0, tmp1;
    double          c[2];
    
    // load triangle coordinates
    const packed_t  v0[3] = { tri[0][0], tri[0][1], tri[0][2] };
    const packed_t  v1[3] = { tri[1][0], tri[1][1], tri[1][2] };
    const packed_t  v2[3] = { tri[2][0], tri[2][1], tri[2][2] };

    for ( size_t  i = 0; i+1 < npts; i +=2 )
    {
        // load quadrature points for first triangle
        tmp0 = load< packed_t >( (*points)[i].vector() );
        tmp1 = load< packed_t >( (*points)[i+1].vector() );

        const packed_t  a1( unpacklo( tmp0, tmp1 ) );
        const packed_t  a2( unpackhi( tmp0, tmp1 ) );
        const packed_t  a0( sub( sub( one, a1 ),  a2 ) );

        // transform to triangle coordinates
        store( add( add( mul( a0, v0[0] ), mul( a1, v1[0] ) ), mul( a2, v2[0] ) ), c );
        x[i][0]   = c[0];
        x[i+1][0] = c[1];

        store( add( add( mul( a0, v0[1] ), mul( a1, v1[1] ) ), mul( a2, v2[1] ) ), c );
        x[i][1]   = c[0];
        x[i+1][1] = c[1];

        store( add( add( mul( a0, v0[2] ), mul( a1, v1[2] ) ), mul( a2, v2[2] ) ), c );
        x[i][2]   = c[0];
        x[i+1][2] = c[1];
    }// for

    if ( npts % 2 == 1 )
    {
        const size_t  i  = npts-1;
        
        const double  a1 = (*points)[i][0];
        const double  a2 = (*points)[i][1];
        const double  a0 = 1.0 - a1 - a2;
        
        x[i] = a0 * tri[0] + a1 * tri[1] + a2 * tri[2];
    }// for
    
#else

    for ( size_t  i = 0; i < npts; i++ )
    {
        const double  a1 = (*points)[i][0];
        const double  a2 = (*points)[i][1];
        const double  a0 = 1.0 - a1 - a2;
        
        x[i] = a0 * tri[0] + a1 * tri[1] + a2 * tri[2];
    }// for

#endif
}

}// namespace
