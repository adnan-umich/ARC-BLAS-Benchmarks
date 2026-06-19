#ifndef __HLIB_BEM_HH
#define __HLIB_BEM_HH
//
// Project     : HLib
// File        : hlib-bem.hh
// Description : HLIBpro include file headers for BEM apps
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//
//! \defgroup BEM_Module BEM
//!
//! This modules provides all functions and class for boundary element
//! methods (BEM), e.g. triangular grids, bilinear forms and
//! quadrature rules.
//!
//! \code
//! #include <hlib-bem.hh>
//! \endcode
//

#include <hpro/bem/TGrid.hh>
#include <hpro/bem/TBEMBF.hh>
#include <hpro/bem/TBFCoeffFn.hh>
#include <hpro/bem/TMassBF.hh>
#include <hpro/bem/TLaplaceBF.hh>
#include <hpro/bem/THelmholtzBF.hh>
#include <hpro/bem/TAcousticBF.hh>
#include <hpro/bem/TExpBF.hh>
#include <hpro/bem/TBEMRHS.hh>
#include <hpro/bem/TGaussQuad.hh>
#include <hpro/bem/TRefinableGrid.hh>
#include <hpro/bem/TSauterTriQuad.hh>

#include <hpro/bem/TConstEdgeFnSpace.hh>
#include <hpro/bem/TMaxwellBF.hh>
#include <hpro/bem/TMaxwellRHS.hh>

#endif // __HLIB_BEM_HH
