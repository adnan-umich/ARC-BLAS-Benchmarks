//
// Project     : HLIBpro
// File        : TMaternCovCoeffFn.cc
// Description : matrix coefficients for matern covariance kernel
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/config.h"

#if HPRO_USE_GSL == 1
#  include <gsl/gsl_sf_bessel.h>
#  include <gsl/gsl_sf_gamma.h>
#else
#  include <cmath>
#endif

#include "hpro/misc/TMaternCovCoeffFn.hh"

namespace Hpro
{

namespace
{

//
// some constants
//
const double  SQRT_3 = std::sqrt( 3 );
const double  SQRT_5 = std::sqrt( 5 );

//
// compute modified bessel function
//
double
compute_Bessel ( const double  d,
                 const double  rho,
                 const double  nu,
                 const double  sigmasq,
                 const double  scale_fac ) // constant scaling factor for optimization
{
    if ( d < Limits::epsilon< double >() )
        return sigmasq;
    else
    {
        const double temp =  d / rho;

        #if HPRO_USE_GSL == 1
        return scale_fac * gsl_sf_bessel_Knu( nu, temp ) * std::pow( temp, nu );
        #else
        return scale_fac * std::cyl_bessel_k( nu, temp ) * std::pow( temp, nu );
        #endif
    }// else
}

}// namespace anonymous

//
// constructor
//
template < typename point_t >
TMaternCovCoeffFn< point_t >::TMaternCovCoeffFn ( const value_t                   sigma,
                                                  const value_t                   length,
                                                  const value_t                   nu,
                                                  const std::vector< point_t > &  vertices )
    : TCoeffFn< value_t >()
    , _length( length )
    , _nu( nu )
    , _sigmasq( sigma*sigma )
    #if HPRO_USE_GSL == 1
    , _scale_fac( _sigmasq / ( std::pow(2.0, nu-1) * gsl_sf_gamma(nu) ) )
    #else
    , _scale_fac( _sigmasq / ( std::pow(2.0, nu-1) * std::tgamma( nu ) ) )
    #endif
    , _row_vertices( vertices )
    , _col_vertices( vertices )
    , _matern_type( general )
{}

template < typename point_t >
TMaternCovCoeffFn< point_t >::TMaternCovCoeffFn ( const value_t                   sigma,
                                                  const value_t                   length,
                                                  const value_t                   nu,
                                                  const std::vector< point_t > &  row_vertices,
                                                  const std::vector< point_t > &  col_vertices )
    : TCoeffFn< value_t >()
    , _length( length )
    , _nu( nu )
    , _sigmasq( sigma*sigma )
    #if HPRO_USE_GSL == 1
    , _scale_fac( _sigmasq / ( std::pow(2.0, nu-1) * gsl_sf_gamma(nu) ) )
    #else
    , _scale_fac( _sigmasq / ( std::pow(2.0, nu-1) * std::tgamma( nu ) ) )
    #endif
    , _row_vertices( row_vertices )
    , _col_vertices( col_vertices )
    , _matern_type( general )
{}

//
// coefficient evaluation
//
template < typename point_t >
void
TMaternCovCoeffFn< point_t >::eval  ( const std::vector< idx_t > &  rowidxs,
                                      const std::vector< idx_t > &  colidxs,
                                      value_t *                     matrix ) const
{
    const size_t  n = rowidxs.size();
    const size_t  m = colidxs.size();

    for ( size_t  j = 0; j < m; ++j )
    {
        const idx_t      idx1 = colidxs[ j ];
        const point_t &  y    = _col_vertices[ idx1 ];
            
        for ( size_t  i = 0; i < n; ++i )
        {
            const idx_t      idx0 = rowidxs[ i ];
            const point_t &  x    = _row_vertices[ idx0 ];
            const value_t    dist = norm2( x - y );

            switch ( _matern_type )
            {
                case one_half :
                {
                    matrix[ j*n + i ] = _sigmasq * std::exp( - dist / _length );
                }
                break;

                case three_half :
                {
                    const auto  sqrt3drho = SQRT_3 * dist / _length;
                    
                    matrix[ j*n + i ] = _sigmasq * ( 1.0 + sqrt3drho ) * std::exp( - sqrt3drho );
                }
                break;

                case five_half :
                {
                    const auto  sqrt5drho = SQRT_5 * dist / _length;
                    
                    matrix[ j*n + i ] = _sigmasq * ( 1.0 + sqrt5drho + ( 5 * dist * dist ) / ( 3 * _length * _length ) ) * std::exp( - sqrt5drho );
                }
                break;

                case general :
                default :
                {
                    matrix[ j*n + i ] = compute_Bessel( dist, _length, _nu, _sigmasq, _scale_fac );
                }
            }// switch
        }// for
    }// for
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template class TMaternCovCoeffFn< TPoint >;
template class TMaternCovCoeffFn< T2Point >;
template class TMaternCovCoeffFn< T3Point >;

}// namespace Hpro
