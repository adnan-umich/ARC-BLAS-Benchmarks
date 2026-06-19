#ifndef __HPRO_TSAUTERTRIQUAD_HH
#define __HPRO_TSAUTERTRIQUAD_HH
//
// Project     : HLIBpro
// File        : TSauterTriQuad.hh
// Description : provides quadrature rules for 2 triangles
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/base/types.hh"
#include "hpro/base/TPoint.hh"

namespace Hpro
{
    
//
// gauss quadrature for triangles in Galerkin method
// based on work from Stefan Sauter
//
class TSauterTriQuad
{
public:
    // special type for quadrature points and weights
    struct rule_t {
        T2Point  pttri1;
        T2Point  pttri2;
        double   wght;
    };

protected:
    // order for ξ
    const uint  _order_xi;
    
public:
    /////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TSauterTriQuad ();

    ~TSauterTriQuad ();

    /////////////////////////////////////////////
    //
    // build quadrature points and weights
    //

    //! construct 2d quadrature points of order \a order for 
    //! a pair of triangles with \a ncommon common vertices;
    //! the quadrature points will be stored in \a tri1_pts and
    //! \a tri2_pts and the quadrature weights in \a w 
    void build ( const uint                ncommon,
                 const uint                order,
                 std::vector< T2Point > &  tri1_pts,
                 std::vector< T2Point > &  tri2_pts,
                 std::vector< double > &   w ) const;

    void build ( const uint               ncommon,
                 const uint               order,
                 std::vector< rule_t > &  rule ) const;

    //! build quadrature points for given triangle pair
    //! - returns the actual number of points
    //! - tri1_pts, tri2_pts and weights have to be at least of size 6*order^4
    uint points ( const uint                       order,
                  const std::vector< double * > &  tri1,
                  const std::vector< double * > &  tri2,
                  std::vector< T3Point > &         tri1_pts,
                  std::vector< T3Point > &         tri2_pts,
                  std::vector< double >  &         weights ) const;
    
protected:
    //
    // quadrature rules for various cases
    //

    // equal triangles
    void quad_equ  ( const uint                order,
                     std::vector< T2Point > &  tri1_pos,
                     std::vector< T2Point > &  tri2_pos,
                     std::vector< double > &   w ) const;

    // triangles with a common edge
    void quad_edge ( const uint                order,
                     std::vector< T2Point > &  tri1_pos,
                     std::vector< T2Point > &  tri2_pos,
                     std::vector< double > &   w ) const;

    // triangles with a common vertex
    void quad_vtx  ( const uint                order,
                     std::vector< T2Point > &  tri1_pos,
                     std::vector< T2Point > &  tri2_pos,
                     std::vector< double > &   w ) const;

    // separated triangles
    void quad_dist ( const uint                order,
                     std::vector< T2Point > &  tri1_pos,
                     std::vector< T2Point > &  tri2_pos,
                     std::vector< double > &   w ) const;
};

}// namespace

#endif // __HPRO_TSAUTERTRIQUAD_HH
