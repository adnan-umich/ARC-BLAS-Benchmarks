#ifndef __HPRO_TSOLVER_HH
#define __HPRO_TSOLVER_HH
//
// Project     : HLIBpro
// File        : TSolver.hh
// Description : base-class for all iterative solvers
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <ostream>
#include <list>

#include "hpro/matrix/TLinearOperator.hh"
#include "hpro/vector/TVector.hh"

namespace Hpro
{

enum solver_status_t
{
    iterating = 0,  // solver is still iterating
    converged,      // solver has converged
    diverged,       // solver has diverged
    failed          // solver has failed
};

//!
//! \ingroup  Solver_Module
//! \struct   TSolverInfo
//! \brief    datatype to share iteration history
//!
struct TSolverInfo
{
    //! @cond
        
public:
    //!
    //! \struct  hist_t
    //! \brief   used to store history information
    //!
    struct hist_t {
        double   res_norm;
    };
        
protected:
    //! name of solver
    std::string          _solver;

    //! status of solver
    solver_status_t      _status;
        
    //! number of iterations
    size_t               _n_iter;
        
    //! average convergance rate
    double               _conv_rate;

    //! residual norm after iteration
    double               _res_norm;

    //! if true, iteration history information is stored
    bool                 _store_hist;

    //! iteration history data
    std::list< hist_t >  _history;

    //! if true, print convergence history
    bool                 _print;

    //! @endcond
        
public:
    //
    // ctor
    //
        
    //! constructor and store history if \a store_hist is true
    //! and print history if \a print is true respectively
    TSolverInfo ( const bool store_hist = false,
                  const bool print      = false );

    //
    // access data
    //
        
    //! return true if convergence information is present
    bool has_data () const { return (n_iter() != 0) || (history().size() > 0); }

    //! return solver name
    const std::string &  solver () const { return _solver; }

    //! return true, if iteration has converged
    bool           has_converged () const { return _status == converged; }

    //! return true, if iteration has diverged
    bool           has_diverged  () const { return _status == diverged; }

    //! return true, if iteration has failed
    bool           has_failed    () const { return _status == failed; }

    //! return status of solver
    solver_status_t status        () const { return _status; }

    //! return number of iteration steps
    size_t         n_iter    () const { return _n_iter; }

    //! return convergence rate
    double         conv_rate () const { return _conv_rate; }

    //! return current residual norm
    double         res_norm  () const { return _res_norm; }

    //! return iteration history data
    const std::list< hist_t > &  history () const { return _history; }

    //! return true, if history is stored
    bool           store_hist () const { return _store_hist; }
        
    //! return true, if convergence history is printed is stored
    bool           print_hist () const { return _print; }
        
    //
    // write internal data and options
    //
        
    //! append data for single iteration step to history (if set so)
    //! and update internal iteration data
    void append         ( const uint it, const double res_norm );

    //! set solver name
    void set_solver     ( const std::string & name ) { _solver = name; }

    //! set convergence status
    void set_status     ( const solver_status_t  s  ) { _status = s; }

    //! set number of iteration steps
    void set_n_iter     ( const uint     n    ) { _n_iter    = n; }

    //! set convergence rate
    void set_conv_rate  ( const double   conv ) { _conv_rate = conv; }

    //! set current residual norm
    void set_res_norm   ( const double   norm ) { _res_norm  = norm; }

    //! turn on/off printing of history
    void set_print_hist ( const bool     b    ) { _print = b; }

    //! reset all data
    void reset ();

    //
    // misc.
    //
        
    //! convert data in information object to string
    std::string  to_string () const;

    //! stream output
    friend std::ostream & operator << ( std::ostream & os, const TSolverInfo & info )
    {
        return os << info.to_string();
    }

    //! write history in Gnuplot format to stream
    void  print_gnuplot ( std::ostream &  os );
};

//!
//! \ingroup  Solver_Module
//! \struct   TStopCriterion
//! \brief     stopping criterion for iterative solvers
//!
class TStopCriterion
{
public:
    // maximal number of iterations
    uint    max_iter;
        
    // absolute reduction of residual norm
    double  abs_res_reduct;
        
    // relative reduction of residual norm compared to start residual
    double  rel_res_reduct;
        
    // maximal relative growth of residual norm compared to start residual
    double  rel_res_growth;

    // ctors
    TStopCriterion ( const uint    max_iter       = CFG::Solver::max_iter,
                     const double  abs_res_red    = CFG::Solver::abs_res_red,
                     const double  rel_res_red    = CFG::Solver::rel_res_red,
                     const double  rel_res_growth = CFG::Solver::rel_res_growth );

    TStopCriterion ( const TStopCriterion &  stop_crit );
    TStopCriterion ( TStopCriterion &&       stop_crit );

    TStopCriterion & operator =  ( const TStopCriterion &  stop_crit );
    TStopCriterion & operator =  ( TStopCriterion &&       stop_crit );
        
    //! return true if stop condition is met
    virtual bool stopped ( const uint     it, 
                           const double     norm, 
                           const double     norm0, 
                           TSolverInfo *  info ) const;

    //! convert to string
    std::string  to_string () const;
};

//!
//! \ingroup  Solver_Module
//! \brief    joins stop criteria while overloading/preferring undefined values
//!           (interprete as logical and)
//!
TStopCriterion
operator + ( const TStopCriterion &  crit1,
             const TStopCriterion &  crit2 );

//!
//! \ingroup  Solver_Module
//! \brief    sets maximal number of iteration steps
//!
TStopCriterion 
max_steps           ( const uint  steps );

//!
//! \ingroup  Solver_Module
//! \brief    sets relative reduction of residual
//!
TStopCriterion 
relative_reduction  ( const double  red );

//!
//! \ingroup  Solver_Module
//! \brief    sets absolute reduction of residual
//!
TStopCriterion 
absolute_reduction  ( const double  red );

//!
//! unions of possible matrix/vector types to enable virtual
//! solve function in TSolver and derived classes
//!
    
//!
//! \ingroup  Solver_Module
//! \class    TSolver
//! \brief    Solver base class defining interface and implementing
//!           simple solver (Richardson iteration)
//!
class TSolver
{
protected:

    //! @cond

    // stop criterion
    TStopCriterion  _stop_criterion;

    // if true, start value of iteration will be initialised by iteration method
    bool            _initialise_start_value;

    // if true, exact residual is computed (default: depends on precond./solver)
    bool            _use_exact_residual;
    
    //! @admcond
    
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct solver object with corresponding stop criteria
    TSolver ( const TStopCriterion &  stop_crit );

    //! dtor
    virtual ~TSolver ();

    //! turn initialisation of start value on/off
    void  initialise_start_value ( const bool  b );

    //! choose type of residual
    void  use_exact_residual     ( const bool  b );

    //! return choice of residual
    bool  use_exact_residual     () const { return _use_exact_residual; }
    
    ////////////////////////////////////////////////
    //
    // solving the system
    //

    //! solve A·x = b with optional preconditioner \a W
    virtual
    void solve ( any_const_operator_t  A,
                 any_vector_t          x,
                 any_const_vector_t    b,
                 any_const_operator_t  W,
                 TSolverInfo *         info = nullptr ) const = 0;
    virtual
    void solve ( any_const_operator_t  A,
                 any_vector_t          x,
                 any_const_vector_t    b,
                 TSolverInfo *         info = nullptr ) const = 0;

    //! set stop criterion
    virtual void set_stop_crit ( const TStopCriterion &  stop_crit );

    //! return true if stop condition is met
    virtual bool stopped       ( const uint     it, 
                                 const double   norm, 
                                 const double   norm0, 
                                 TSolverInfo *  info ) const;
    
    //! initialises start value of iteration
    void set_start_value ( any_vector_t          x,
                           any_const_vector_t    b ) const;
    void set_start_value ( any_vector_t          x,
                           any_const_vector_t    b,
                           any_const_operator_t  W ) const;
};

//
// to instantiate mixed precision solve method
//
#define HPRO_INST_SINGLE_SOLVE_METHOD( solverclass, type1, type2 )      \
    template void solverclass :: solve< type1, type2 > ( const TLinearOperator< type1 > *, \
                                                         TVector< type1 > *              , \
                                                         const TVector< type1 > *        , \
                                                         const TLinearOperator< type2 > *, \
                                                         TSolverInfo *                    ) const;

#define HPRO_INST_SOLVE_METHOD( solverclass )                           \
    HPRO_INST_SINGLE_SOLVE_METHOD( solverclass, float, float )          \
    HPRO_INST_SINGLE_SOLVE_METHOD( solverclass, float, double )         \
    HPRO_INST_SINGLE_SOLVE_METHOD( solverclass, double, float )         \
    HPRO_INST_SINGLE_SOLVE_METHOD( solverclass, double, double )        \
    HPRO_INST_SINGLE_SOLVE_METHOD( solverclass, std::complex< float >, std::complex< float > ) \
    HPRO_INST_SINGLE_SOLVE_METHOD( solverclass, std::complex< float >, std::complex< double > ) \
    HPRO_INST_SINGLE_SOLVE_METHOD( solverclass, std::complex< double >, std::complex< float > ) \
    HPRO_INST_SINGLE_SOLVE_METHOD( solverclass, std::complex< double >, std::complex< double > )

//
// mixed precision solve method calling template method of particular solver
//
template < typename solver_t >
void
solve_mp ( const solver_t &      solver,
           const std::string &   solver_name,
           any_const_operator_t  A,
           any_vector_t          x,
           any_const_vector_t    b,
           any_const_operator_t  W,
           TSolverInfo *         info )
{
    using std::get;
    
    if (( A.index() != x.index() ) || ( A.index() != b.index() ))
        HERROR( ERR_ARG, "(" + solver_name + ") solve", "A, x and b have different value type" );

    switch ( W.index() )
    {
        case REAL_FP32 :
            switch ( A.index() )
            {
                case REAL_FP32 : solver.solve( get< REAL_FP32 >( A ), get< REAL_FP32 >( x ), get< REAL_FP32 >( b ), get< REAL_FP32 >( W ), info ); break;
                case REAL_FP64 : solver.solve( get< REAL_FP64 >( A ), get< REAL_FP64 >( x ), get< REAL_FP64 >( b ), get< REAL_FP32 >( W ), info ); break;
                default:
                    HERROR( ERR_ARG, "(" + solver_name + ") solve", "A, x, b and W have unsupported value type" );
            }// switch
            break;
            
        case REAL_FP64 :
            switch ( A.index() )
            {
                case REAL_FP32 : solver.solve( get< REAL_FP32 >( A ), get< REAL_FP32 >( x ), get< REAL_FP32 >( b ), get< REAL_FP64 >( W ), info ); break;
                case REAL_FP64 : solver.solve( get< REAL_FP64 >( A ), get< REAL_FP64 >( x ), get< REAL_FP64 >( b ), get< REAL_FP64 >( W ), info ); break;
                default:
                    HERROR( ERR_ARG, "(" + solver_name + ") solve", "A, x, b and W have unsupported value type" );
            }// switch
            break;
            
        case COMPLEX_FP32 :
            switch ( A.index() )
            {
                case COMPLEX_FP32 : solver.solve( get< COMPLEX_FP32 >( A ), get< COMPLEX_FP32 >( x ), get< COMPLEX_FP32 >( b ), get< COMPLEX_FP32 >( W ), info ); break;
                case COMPLEX_FP64 : solver.solve( get< COMPLEX_FP64 >( A ), get< COMPLEX_FP64 >( x ), get< COMPLEX_FP64 >( b ), get< COMPLEX_FP32 >( W ), info ); break;
                default:
                    HERROR( ERR_ARG, "(" + solver_name + ") solve", "A, x, b and W have unsupported value type" );
            }// switch
            break;
            
        case COMPLEX_FP64 :
            switch ( A.index() )
            {
                case COMPLEX_FP32 : solver.solve( get< COMPLEX_FP32 >( A ), get< COMPLEX_FP32 >( x ), get< COMPLEX_FP32 >( b ), get< COMPLEX_FP64 >( W ), info ); break;
                case COMPLEX_FP64 : solver.solve( get< COMPLEX_FP64 >( A ), get< COMPLEX_FP64 >( x ), get< COMPLEX_FP64 >( b ), get< COMPLEX_FP64 >( W ), info ); break;
                default:
                    HERROR( ERR_ARG, "(" + solver_name + ") solve", "A, x, b and W have unsupported value type" );
            }// switch
            break;

        default :
            HERROR( ERR_ARG, "(" + solver_name + ") solve", "A, x, b and W have unsupported value type" );
    }// switch
}

}// namespace Hpro

#endif  // __HPRO_TSOLVER_HH
