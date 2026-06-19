#ifndef __HPRO_TAUTOSOLVER_HH
#define __HPRO_TAUTOSOLVER_HH
//
// Project     : HLIBpro
// File        : TAutoSolver.hh
// Description : class implementing a solver which automatically decides best strategy
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/solver/TSolver.hh"

namespace Hpro
{

//!
//! \ingroup  Solver_Module
//! \class    TAutoSolver
//! \brief    Implements an iterative solver automatically choosing
//!           appropriate algorithm based on matrix criteria.
//!
class TAutoSolver : public TSolver
{
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct auto solver object with corresponding stop criteria
    TAutoSolver ( const TStopCriterion &  stop_crit = TStopCriterion() );
    
    virtual ~TAutoSolver ();

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
//! \brief    Solve A·x = b with optional preconditioner \a W (functional version).
//!
template < typename value_t,
           typename value_pre_t >
void
solve ( const TLinearOperator< value_t > *      A,
        TVector< value_t > *                    x,
        const TVector< value_t > *              b,
        const TLinearOperator< value_pre_t > *  W,
        TSolverInfo *                           info      = nullptr,
        const TStopCriterion &                  stop_crit = TStopCriterion() )
{
    TAutoSolver  solver( stop_crit );

    solver.solve( A, x, b, W, info );
}

template < typename value_t,
           typename value_pre_t >
void
solve ( const TLinearOperator< value_t > &      A,
        TVector< value_t > &                    x,
        const TVector< value_t > &              b,
        const TLinearOperator< value_pre_t > &  W,
        TSolverInfo *                           info      = nullptr,
        const TStopCriterion &                  stop_crit = TStopCriterion() )
{
    solve( &A, &x, &b, &W, info, stop_crit );
}

//!
//! \ingroup  Solver_Module
//! \brief    Solve A·x = b (functional version).
//!
template < typename value_t >
void
solve ( const TLinearOperator< value_t > *  A,
        TVector< value_t > *                x,
        const TVector< value_t > *          b,
        TSolverInfo *                       info      = nullptr,
        const TStopCriterion &              stop_crit = TStopCriterion() )
{
    solve< value_t, value_t >( A, x, b, nullptr, info, stop_crit );
}

template < typename value_t >
void
solve ( const TLinearOperator< value_t > &  A,
        TVector< value_t > &                x,
        const TVector< value_t > &          b,
        TSolverInfo *                       info      = nullptr,
        const TStopCriterion &              stop_crit = TStopCriterion() )
{
    solve< value_t, value_t >( &A, &x, &b, nullptr, info, stop_crit );
}

}// namespace Hpro

#endif  // __HPRO_TAUTOSOLVER_HH
