#ifndef __HLIB_SOLVER_HH
#define __HLIB_SOLVER_HH
//
// Project     : HLib
// File        : hlib-solver.hh
// Description : HLIBpro include file containing iterative solver related headers
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//
//! \defgroup Solver_Module Solver
//!
//! This modules provides classes for iterative solvers, e.g. CG, BiCG-Stab, GMRES.
//!
//! \code
//! #include <hlib-solver.hh>
//! \endcode
//

#include <hpro/solver/TSolver.hh>
#include <hpro/solver/TLinearIteration.hh>
#include <hpro/solver/TBiCGStab.hh>
#include <hpro/solver/TCG.hh>
#include <hpro/solver/TCGS.hh>
#include <hpro/solver/TMINRES.hh>
#include <hpro/solver/TGMRES.hh>
#include <hpro/solver/TTFQMR.hh>
#include <hpro/solver/TAutoSolver.hh>

#endif  // __HLIB_SOLVER_HH
