#ifndef __HPRO_TKACZMARZ_HH
#define __HPRO_TKACZMARZ_HH
//
// Project     : HLIBpro
// File        : TKaczmarz.hh
// Description : iterative solver based on Kaczmarz algorithm
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/solver/TSolver.hh"

namespace Hpro
{

//! @cond

//!
//! \ingroup  Solver_Module
//! \class    TKaczmarz
//! \brief    Implements Kaczmarz iteration.
//!
class TKaczmarz : public TSolver
{
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct Kaczmarz solver object with corresponding stop criteria
    TKaczmarz ( const TStopCriterion &  stop_crit = TStopCriterion() );

    virtual ~TKaczmarz ();

    ////////////////////////////////////////////////
    //
    // solving the system
    //

    //! solve A·x = b with optional preconditioner \a W
    template < typename value_t,
               typename value_pre_t >
    void solve ( const TLinearOperator< value_t > *      A,
                 TVector< value_t > *                    x,
                 const TVector< value_t > *              b,
                 const TLinearOperator< value_pre_t > *  W    = nullptr,
                 TSolverInfo *                           info = nullptr ) const;

    //! generic implementation for "virtual" solve method
    virtual
    void solve ( any_const_operator_t  A,
                 any_vector_t          x,
                 any_const_vector_t    b,
                 any_const_operator_t  W,
                 TSolverInfo *         info = nullptr ) const;
    virtual
    void solve ( any_const_operator_t  A,
                 any_vector_t          x,
                 any_const_vector_t    b,
                 TSolverInfo *         info = nullptr ) const;
};

//!
//! \ingroup  Solver_Module
//! \brief    Solve A·x = b with optional preconditioner \a W (functional approach)
//!
template < typename value_t,
           typename value_pre_t >
void
kaczmarz ( const TLinearOperator< value_t > *      A,
           TVector< value_t > *                    x,
           const TVector< value_t > *              b,
           const TLinearOperator< value_pre_t > *  W         = nullptr,
           TSolverInfo *                           info      = nullptr,
           const TStopCriterion &                  stop_crit = TStopCriterion() )
{
    TKaczmarz  solver( stop_crit );

    solver.solve( A, x, b, W, info );
}

//! @endcond

}// namespace Hpro

#endif  // __HPRO_TKACZMARZ_HH
