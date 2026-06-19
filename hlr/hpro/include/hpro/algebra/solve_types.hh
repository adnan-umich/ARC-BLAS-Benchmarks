#ifndef __HPRO_SOLVE_TYPES_HH
#define __HPRO_SOLVE_TYPES_HH
//
// Project     : HLIBpro
// File        : solve_types.hh
// Description : special basic types for evaluation and solve algorithms
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/algebra/types.hh"

namespace Hpro
{

/////////////////////////////////////////////////////////////////////////////////////
//
// special types for solve algorithm
//
/////////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup Algebra_Module
//! \struct  solve_option_t
//! \brief   options for how to solve with given matrix
//!
struct solve_option_t
{
    //! is matrix to be evaluated point or block wise
    eval_type_t     eval;

    //! is diagonal unit or not
    diag_type_t     diag;

    //! do diagonal blocks hold inverse or not
    storage_type_t  storage;
    
    //! constructor
    solve_option_t ( const eval_type_t     aeval,
                     const diag_type_t     adiag,
                     const storage_type_t  astor = CFG::Arith::storage_type )
            : eval( aeval )
            , diag( adiag )
            , storage( astor )
    {}
};

/////////////////////////////////////////////////////////////////////////////////////
//
// special types for eval algorithm
//
/////////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup Algebra_Module
//! \struct  eval_option_t
//! \brief   options for how to evaluate given matrix
//!
struct eval_option_t
{
    //! is matrix to be evaluated point or block wise
    eval_type_t     eval;

    //! is diagonal unit or not
    diag_type_t     diag;

    //! do diagonal blocks hold inverse or not
    storage_type_t  storage;
    
    //! constructor
    eval_option_t ( const eval_type_t     aeval,
                    const diag_type_t     adiag,
                    const storage_type_t  astor = CFG::Arith::storage_type )
            : eval( aeval )
            , diag( adiag )
            , storage( astor )
    {}
};

}// namespace Hpro

#endif  // __HPRO_SOLVE_TYPES_HH
