#ifndef __HPRO_TLINEARITERATION_HH
#define __HPRO_TLINEARITERATION_HH
//
// Project     : HLIBpro
// File        : TLinearIteration.hh
// Description : Linear Iteration solver
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <ostream>
#include <list>

#include "hpro/solver/TSolver.hh"

namespace Hpro
{

//!
//! \ingroup  Solver_Module
//! \class    TLinearIteration
//! \brief    Implements linear iteration \f$x_{i+1} = x_k + N (A x_i - b)\f$
//!
//! \detail   The class TLinearIteration implements a linear iteration solver
//!           based on the second normal form \f$x_{i+1} = x_k + N (A x_i - b)\f$
//!           with the iteration matrix \f$N\f$. Please note, that the iteration
//!           matrix is defined as the preconditioner \f$W\f$, e.g., for other
//!           linear iterations like Jacobi, GS or SOR, please use the correct
//!           preconditioner.
//!
class TLinearIteration : public TSolver
{
protected:
    //! @cond

    // damping factor
    double  _damping;

    //! @endcond
    
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct linear iteration solver object 
    TLinearIteration ( const TStopCriterion &  stop_crit = TStopCriterion() );
    
    //! construct linear iteration solver object with damping
    TLinearIteration ( const double            damping,
                       const TStopCriterion &  stop_crit = TStopCriterion() );
    
    //! dtor
    virtual ~TLinearIteration ();

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
    
    //! solve A·X = B with optional preconditioner \a W
    template < typename value_t >
    void solve  ( const TLinearOperator< value_t > *  A,
                  TMatrix< value_t > *                X,
                  const TMatrix< value_t > *          B,
                  const TLinearOperator< value_t > *  W    = nullptr,
                  TSolverInfo *                       info = nullptr ) const;

    //! solve A·x = b with optional preconditioner \a W
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
linear_iteration ( const TLinearOperator< value_t > *      A,
                   TVector< value_t > *                    x,
                   const TVector< value_t > *              b,
                   const TLinearOperator< value_pre_t > *  W         = nullptr,
                   TSolverInfo *                           info      = nullptr,
                   const TStopCriterion &                  stop_crit = TStopCriterion() )
{
    TLinearIteration  solver( 1.0, stop_crit );

    solver.solve( A, x, b, W, info );
}

/////////////////////////////////////////////////////
//
// for backward compatibility
//

using TRichardson = TLinearIteration;

//!
//! \ingroup  Solver_Module
//! \brief    Solve A·x = b with optional preconditioner \a W (functional approach)
//!
template < typename value_t >
void richardson ( const TLinearOperator< value_t > *  A,
                  TVector< value_t > *                x,
                  const TVector< value_t > *          b,
                  const TLinearOperator< value_t > *  W         = nullptr,
                  TSolverInfo *                       info      = nullptr,
                  const TStopCriterion &              stop_crit = TStopCriterion() )
{
    TLinearIteration  solver( 1.0, stop_crit );

    solver.solve( A, x, b, W, info );
}

}// namespace Hpro

#endif  // __HPRO_TRICHARDSON_HH
