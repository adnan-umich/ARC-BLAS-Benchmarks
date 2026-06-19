#ifndef __HPRO_TGMRES_HH
#define __HPRO_TGMRES_HH
//
// Project     : HLIBpro
// File        : TGMRES.hh
// Description : class implementing GMRES
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/solver/TSolver.hh"

namespace Hpro
{

//!
//! \ingroup  Solver_Module
//! \class    TGMRES
//! \brief    Implements GMRES iteration with restart.
//!
class TGMRES : public TSolver
{
private:
    //! maximal dimension of Krylov space before restart
    const uint  _restart;
    
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct GMRES solver object with default restart
    TGMRES ( const TStopCriterion &  stop_crit = TStopCriterion() );
    
    //! construct GMRES solver object with given restart
    TGMRES ( const uint              restart,
             const TStopCriterion &  stop_crit = TStopCriterion() );
    
    virtual ~TGMRES ();

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
//! \fn       gmres
//! \brief    Solve A·x = b with optional preconditioner \a W (functional approach)
//!
template < typename value_t,
           typename value_pre_t >
void
gmres ( const uint                              restart,
        const TLinearOperator< value_t > *      A,
        TVector< value_t > *                    x,
        const TVector< value_t > *              b,
        const TLinearOperator< value_pre_t > *  W         = nullptr,
        TSolverInfo *                           info      = nullptr,
        const TStopCriterion &                  stop_crit = TStopCriterion() )
{
    TGMRES  solver( restart, stop_crit );

    solver.solve( A, x, b, W, info );
}

}// namespace Hpro

#endif  // __HPRO_TGMRES_HH
