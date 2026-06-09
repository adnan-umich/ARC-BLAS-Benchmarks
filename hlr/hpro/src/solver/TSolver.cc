//
// Project     : HLIBpro
// File        : TSolver.cc
// Description : base-class for all iterative solvers
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <iostream>
#include <cmath>
#include <sstream>

#include <boost/format.hpp>

#include "hpro/solver/TSolver.hh"

namespace Hpro
{

////////////////////////////////////////////////
//
// constructor and destructor
//

TSolver::TSolver ( const TStopCriterion &  stop_crit )
        : _stop_criterion( stop_crit )
        , _initialise_start_value( CFG::Solver::init_start_value )
        , _use_exact_residual( CFG::Solver::use_exact_residual )
{}

TSolver::~TSolver ()
{}

//
// turn initialisation of start value on/off
//
void
TSolver::initialise_start_value ( const bool  b )
{
    _initialise_start_value = b;
}

//
// choose exact/default residual
//
void
TSolver::use_exact_residual ( const bool  b )
{
    _use_exact_residual = b;
}

//
// set stop criterion
//
void
TSolver::set_stop_crit ( const TStopCriterion &  stop_crit )
{
    _stop_criterion = stop_crit;
}

//
// return true if stop condition is met
//
bool
TSolver::stopped ( const uint     it, 
                   const double     norm, 
                   const double     norm0, 
                   TSolverInfo *  info ) const
{
    return _stop_criterion.stopped( it, norm, norm0, info );
}

//
// initialises start value of iteration
//
void
TSolver::set_start_value ( any_vector_t          x,
                           any_const_vector_t    b ) const
{
    using std::get;

    if ( ! _initialise_start_value )
        return;
    
    if ( x.index() != b.index() )
        HERROR( ERR_ARG, "(TSolver) set_start_value", "x and b are of different type" );
                    
    switch ( x.index() )
    {
        case 0 : get< 0 >( x )->assign( float(1), get< 0 >( b ) ); break;
        case 1 : get< 1 >( x )->assign( double(1), get< 1 >( b ) ); break;
        case 2 : get< 2 >( x )->assign( std::complex< float >(1), get< 2 >( b ) ); break;
        case 3 : get< 3 >( x )->assign( std::complex< double >(1), get< 3 >( b ) ); break;
        default:
            HERROR( ERR_ARG, "(TSolver) set_start_value", "unknown vector value type" );
    }// switch
}
void
TSolver::set_start_value ( any_vector_t          x,
                           any_const_vector_t    b,
                           any_const_operator_t  W ) const
{
    using std::get;
    
    if ( ! _initialise_start_value )
        return;
    
    if ( x.index() != b.index() )
        HERROR( ERR_ARG, "(TSolver) set_start_value", "x and b are of different type" );

    if ( std::visit( [] ( auto &&  W_ptr ) { return ( W_ptr == nullptr ? true : false ); }, W ) )
    {
        set_start_value( x, b );
        return;
    }// if
        
    switch ( W.index() )
    {
        case REAL_FP32 :
            switch ( x.index() )
            {
                case REAL_FP32 : apply( get< REAL_FP32 >( W ), get< REAL_FP32 >( b ), get< REAL_FP32 >( x ) ); break;
                case REAL_FP64 : apply( get< REAL_FP32 >( W ), get< REAL_FP64 >( b ), get< REAL_FP64 >( x ) ); break;
                default:
                    HERROR( ERR_ARG, "(TSolver) set_start_value", "unknown vector/matrix type" );
            }// switch
            break;
            
        case REAL_FP64 :
            switch ( x.index() )
            {
                case REAL_FP32 : apply( get< REAL_FP64 >( W ), get< REAL_FP32 >( b ), get< REAL_FP32 >( x ) ); break;
                case REAL_FP64 : apply( get< REAL_FP64 >( W ), get< REAL_FP64 >( b ), get< REAL_FP64 >( x ) ); break;
                default:
                    HERROR( ERR_ARG, "(TSolver) set_start_value", "unknown vector/matrix type" );
            }// switch
            break;
            
        case COMPLEX_FP32 :
            switch ( x.index() )
            {
                case COMPLEX_FP32 : apply( get< COMPLEX_FP32 >( W ), get< COMPLEX_FP32 >( b ), get< COMPLEX_FP32 >( x ) ); break;
                case COMPLEX_FP64 : apply( get< COMPLEX_FP32 >( W ), get< COMPLEX_FP64 >( b ), get< COMPLEX_FP64 >( x ) ); break;
                default:
                    HERROR( ERR_ARG, "(TSolver) set_start_value", "unknown vector/matrix type" );
            }// switch
            break;
            
        case COMPLEX_FP64 :
            switch ( x.index() )
            {
                case COMPLEX_FP32 : apply( get< COMPLEX_FP64 >( W ), get< COMPLEX_FP32 >( b ), get< COMPLEX_FP32 >( x ) ); break;
                case COMPLEX_FP64 : apply( get< COMPLEX_FP64 >( W ), get< COMPLEX_FP64 >( b ), get< COMPLEX_FP64 >( x ) ); break;
                default:
                    HERROR( ERR_ARG, "(TSolver) set_start_value", "unknown vector/matrix type" );
            }// switch
            break;

        default :
            HERROR( ERR_ARG, "(TSolver) set_start_value", "unknown vector/matrix type" );
    };
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//
// TSolverInfo
//
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//
// constructor
//
TSolverInfo::TSolverInfo ( const bool  astore_hist,
                           const bool  aprint )
{
    reset();

    _store_hist = astore_hist;
    _print      = aprint;
}

void
TSolverInfo::append ( const uint  it,
                      const double  residual_norm )
{
    //
    // reset info if first entry
    //

    if ( it == 0 )
        reset();
    
    //
    // at least store first entry for future computations
    //

    double  old_resnorm = -1;

    if ( ! history().empty() && _store_hist )
    {
        old_resnorm = history().back().res_norm;
    }// if

    const double  reduction = residual_norm / old_resnorm;
    
    if ( history().empty() || _store_hist )
    {
        hist_t  h = { residual_norm };
        
        _history.push_back( h );
    }// if
    
    // update entries
    double  res_norm0 = history().front().res_norm;

    _n_iter    = it;
    _res_norm  = residual_norm;

    if ( it > 0 )
        _conv_rate = std::pow( residual_norm / res_norm0, double(1) / double(it) );

    
    if ( _print )
    {
        const bool  use_unicode = false;
        const bool  use_color   = false;
                              
        if ( it == 0 )
        {
            if ( use_unicode )
            {
                if ( _store_hist )
                    std::cout << boost::format( "%|=6| │ %|=13| │ %|=13| │ %|=13|" ) % "step" % "|b-Ax|" % "red." % "rate"
                              << std::endl
                              << "───────┼───────────────┼───────────────┼───────────────"
                              << std::endl;
                else
                    std::cout << boost::format( "%|=6| │ %|=13| │ %|=13|" ) % "step" % "|b-Ax|" % "rate"
                              << std::endl
                              << "───────┼───────────────┼───────────────"
                              << std::endl;
            }// if
            else
            {
                if ( _store_hist )
                    std::cout << boost::format( "%|=6| | %|=13| | %|=13| | %|=13|" ) % "step" % "|b-Ax|" % "red." % "rate"
                              << std::endl
                              << "-------+---------------+---------------+---------------"
                              << std::endl;
                else
                    std::cout << boost::format( "%|=6| | %|=13| | %|=13|" ) % "step" % "|b-Ax|" % "rate"
                              << std::endl
                              << "-------+---------------+---------------"
                              << std::endl;
            }// else
        }// if
        
        if ( use_unicode )
            std::cout << boost::format( "%|6| │ %13.4e │ " ) % it % residual_norm;
        else
            std::cout << boost::format( "%|6| | %13.4e | " ) % it % residual_norm;

        if ( it > 0 )
        {
            if ( use_color )
            {
                if ( _store_hist )
                {
                    if ( reduction <= 1 )
                        std::cout << boost::format( "%13.4e " ) % reduction;
                    else
                        std::cout << boost::format( "%13.4e " ) % reduction;

                    std::cout << "│ ";
                    
                    if ( _conv_rate <= 1 )
                        std::cout << boost::format( "%13.4e " ) % _conv_rate;
                    else
                        std::cout << boost::format( "%13.4e " ) % _conv_rate;
                }// if
                else
                {
                    if ( _conv_rate <= 1 )
                        std::cout << boost::format( "%13.4e " ) % _conv_rate;
                    else
                        std::cout << boost::format( "%13.4e " ) % _conv_rate;
                }// else
            }// if
            else
            {
                if ( use_unicode )
                {
                    if ( _store_hist )
                        std::cout << boost::format( "%13.4e │ %13.4e " ) % reduction % _conv_rate;
                    else
                        std::cout << boost::format( "%13.4e" ) % _conv_rate;
                }// if
                else
                {
                    if ( _store_hist )
                        std::cout << boost::format( "%13.4e | %13.4e " ) % reduction % _conv_rate;
                    else
                        std::cout << boost::format( "%13.4e" ) % _conv_rate;
                }// else
            }// else
        }// if
        else
        {
            if ( _store_hist )
            {
                if ( use_unicode ) std::cout << "              │               ";
                else               std::cout << "              |               ";
            }// if
        }// else

        std::cout << std::endl;
    }// if
}

//
// reset data
//
void
TSolverInfo::reset ()
{
    _status     = iterating;
    _n_iter     = 0;
    _conv_rate  = double(0);
    _res_norm   = double(0);

    _history.clear();
}
        
//
// output
//
std::string
TSolverInfo::to_string () const
{
    if ( ! has_data() )
        return "no convergence information defined";
    else
    {
        std::ostringstream  str;

        const bool  use_color   = false;

        if ( use_color )
        {
            if ( has_converged() )
                str << solver()
                    << " : "
                    << "converged"
                    << " after "
                    << n_iter() << " steps"
                    << " with "
                    << "rate = " << boost::format( "%.4e" ) % conv_rate();
            else if ( has_diverged() )
                str << solver()
                    << " : "
                    << "diverged"
                    << " after "
                    << n_iter() << " steps";
            else if ( has_failed() )
                str << solver()
                    << " : "
                    << "failed"
                    << " after "
                    << n_iter() << " steps";
            else
                str << solver()
                    << " : "
                    << "not converged"
                    << " after "
                    << n_iter() << " steps"
                    << " with "
                    << "rate = " << boost::format( "%.4e" ) % conv_rate();
        }// if
        else
        {
            if ( has_converged() )
                str << solver() << " : converged after " << n_iter() << " steps"
                    << " with rate = " << boost::format( "%.4e" ) % conv_rate();
            else if ( has_diverged() )
                str << solver() << " : diverged after " << n_iter() << " steps";
            else if ( has_failed() )
                str << solver() << " : failed after " << n_iter() << " steps";
            else
                str << solver() << " : not converged after " << n_iter() << " steps"
                    << " with rate = " << boost::format( "%.4e" ) % conv_rate();
        }// else
        
        str << " ( |Ax-b| = " << boost::format( "%.4e" ) % res_norm() << " )";

        return str.str();
    }// else
}

//
// write history in Gnuplot format to stream
//
void
TSolverInfo::print_gnuplot ( std::ostream &  os )
{
    // only print, if at least two entries are available
    if ( _history.size() <= 1 )
        return;

    //
    // just write data to file
    //

    uint  i = 0;

    os << "# See GNUPLOT manual for a description of the format" << std::endl
       << std::endl
       << "# <step>  <residual norm>" << std::endl;
    
    for ( auto  iter = _history.cbegin(); iter != _history.cend(); ++iter, ++i )
    {
        os << i << ' ' << (*iter).res_norm << std::endl;
    }// for
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//
// TStopCriterion
//
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TStopCriterion::TStopCriterion ( const uint  amax_iter,
                                 const double  aabs_res_red,
                                 const double  arel_res_red,
                                 const double  arel_res_growth )
        : max_iter( amax_iter )
        , abs_res_reduct( aabs_res_red )
        , rel_res_reduct( arel_res_red )
        , rel_res_growth( arel_res_growth )
{
}

TStopCriterion::TStopCriterion ( const TStopCriterion &  stop_crit )
{
    *this = stop_crit;
}

TStopCriterion::TStopCriterion ( TStopCriterion &&  stop_crit )
{
    *this = stop_crit;
}

TStopCriterion &
TStopCriterion::operator =  ( const TStopCriterion &  stop_crit )
{
    max_iter       = ( stop_crit.max_iter       >  0 ? stop_crit.max_iter       : CFG::Solver::max_iter );
    abs_res_reduct = ( stop_crit.abs_res_reduct >= 0 ? stop_crit.abs_res_reduct : CFG::Solver::abs_res_red );
    rel_res_reduct = ( stop_crit.rel_res_reduct >= 0 ? stop_crit.rel_res_reduct : CFG::Solver::rel_res_red );
    rel_res_growth = ( stop_crit.rel_res_growth >= 0 ? stop_crit.rel_res_growth : CFG::Solver::rel_res_growth );

    return *this;
}

TStopCriterion &
TStopCriterion::operator =  ( TStopCriterion &&  stop_crit )
{
    max_iter       = ( stop_crit.max_iter       >  0 ? stop_crit.max_iter       : CFG::Solver::max_iter );
    abs_res_reduct = ( stop_crit.abs_res_reduct >= 0 ? stop_crit.abs_res_reduct : CFG::Solver::abs_res_red );
    rel_res_reduct = ( stop_crit.rel_res_reduct >= 0 ? stop_crit.rel_res_reduct : CFG::Solver::rel_res_red );
    rel_res_growth = ( stop_crit.rel_res_growth >= 0 ? stop_crit.rel_res_growth : CFG::Solver::rel_res_growth );

    return *this;
}

//
// stop criterion
//
bool
TStopCriterion::stopped ( const uint     it,
                          const double     norm,
                          const double     norm0,
                          TSolverInfo *  info ) const
{
    if (( max_iter > 0 ) && ( it >= max_iter ))
    {
        if ( info != nullptr )
        {
            info->set_status( iterating );
        }// if

        return true;
    }// if

    if ((( abs_res_reduct > 0 ) && ( norm < abs_res_reduct )) ||
        (( rel_res_reduct > 0 ) && ( norm < (rel_res_reduct * norm0) )))
    {
        if ( info != nullptr )
        {
            info->set_status( converged );
        }// if

        return true;
    }// if

    if (( rel_res_growth > 0 ) && ( norm > (rel_res_growth * norm0) ))
    {
        if ( info != nullptr )
        {
            info->set_status( diverged );
        }// if
        
        return true;
    }// if

    return false;
}

//
// convert to string
//
std::string
TStopCriterion::to_string () const
{
    return Hpro::to_string( "stop( max steps = %d, absolute reduction = %.4e, relative reduction = %.4e )",
                            max_iter, abs_res_reduct, rel_res_reduct );
}

//
// joins stop criteria while overloading/preferring undefined values
//
TStopCriterion
operator + ( const TStopCriterion &  crit1,
             const TStopCriterion &  crit2 )
{
    TStopCriterion  res;

    if ( crit1.max_iter > 0 )
    {
        if ( crit2.max_iter > 0 )
            res.max_iter = std::max( crit1.max_iter, crit2.max_iter );
        else
            res.max_iter = crit1.max_iter;
    }// if
    else
    {
        res.max_iter = crit2.max_iter;
    }// else
    
    if ( crit1.abs_res_reduct >= 0 )
    {
        if ( crit2.abs_res_reduct >= 0 )
            res.abs_res_reduct = std::min( crit1.abs_res_reduct, crit2.abs_res_reduct );
        else
            res.abs_res_reduct = crit1.abs_res_reduct;
    }// if
    else
    {
        res.abs_res_reduct = crit2.abs_res_reduct;
    }// else
    
    if ( crit1.rel_res_reduct >= 0 )
    {
        if ( crit2.rel_res_reduct >= 0 )
            res.rel_res_reduct = std::min( crit1.rel_res_reduct, crit2.rel_res_reduct );
        else
            res.rel_res_reduct = crit1.rel_res_reduct;
    }// if
    else
    {
        res.rel_res_reduct = crit2.rel_res_reduct;
    }// else
    
    if ( crit1.rel_res_growth >= 0 )
    {
        if ( crit2.rel_res_growth >= 0 )
            res.rel_res_growth = std::max( crit1.rel_res_growth, crit2.rel_res_growth );
        else
            res.rel_res_growth = crit1.rel_res_growth;
    }// if
    else
    {
        res.rel_res_growth = crit2.rel_res_growth;
    }// else

    return res;
}

TStopCriterion 
max_steps ( const uint  steps )
{
    TStopCriterion  res;

    res.max_iter       = steps;
    res.abs_res_reduct = -1;
    res.rel_res_reduct = -1;
    res.rel_res_growth = -1;

    return res;
}

TStopCriterion 
relative_reduction  ( const double  red )
{
    TStopCriterion  res;

    res.max_iter       = 0;
    res.abs_res_reduct = -1;
    res.rel_res_reduct = red;
    res.rel_res_growth = -1;

    return res;
}

TStopCriterion 
absolute_reduction  ( const double  red )
{
    TStopCriterion  res;

    res.max_iter       = 0;
    res.abs_res_reduct = red;
    res.rel_res_reduct = -1;
    res.rel_res_growth = -1;

    return res;
}

}// namespace Hpro
