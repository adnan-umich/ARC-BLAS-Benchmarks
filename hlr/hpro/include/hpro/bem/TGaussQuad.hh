#ifndef __HPRO_TGAUSSQUAD_HH
#define __HPRO_TGAUSSQUAD_HH
//
// Project     : HLIBpro
// File        : TGaussQuad.hh
// Description : provides Gauss quadrature rules
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/base/types.hh"
#include "hpro/base/TPoint.hh"

namespace Hpro
{
    
//
// implements 1D Gauss quadrature rules
//
class TGaussQuad
{
public:
    /////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TGaussQuad () {}

    /////////////////////////////////////////////
    //
    // build quadrature points and weights
    //

    // construct quadrature points \a pos and weights \a wghts in [0,1]
    void build ( const uint               order,
                 std::vector< double > &  pos,
                 std::vector< double > &  wghts ) const;
};

//
// implements 2D Gauss quadrature rules for triangles
//
class TTriGaussQuad
{
public:
    /////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TTriGaussQuad () {}

    /////////////////////////////////////////////
    //
    // build quadrature points and weights
    //

    //
    // build quadrature points for 2D simplex
    //
    void build ( const uint                      order,
                 std::vector< T2Point > &        pts,
                 std::vector< double > &         wghts ) const;
                 
    //
    // construct quadrature points (<x>) and
    // weights (<w>) for triangle <tri> of
    // order <order>
    //
    void build ( const uint                      order,
                 const std::vector< T3Point > &  tri,
                 std::vector< T3Point > &        pts,
                 std::vector< double > &         wghts ) const;
};

}// namespace

#endif // __HPRO_TGAUSSQUAD_HH
