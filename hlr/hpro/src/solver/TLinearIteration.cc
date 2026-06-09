//
// Project     : HLIBpro
// File        : TLinearIteration.cc
// Description : base-class for all iterative solvers
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <cmath>
#include <sstream>

#include "hpro/solver/TLinearIteration.hh"

namespace Hpro
{

// import from TSolver
void
solve_mp ( const TSolver &       solver,
           const std::string &   solver_name,
           any_const_operator_t  A,
           any_vector_t          x,
           any_const_vector_t    b,
           any_const_operator_t  W,
           TSolverInfo *         info );

////////////////////////////////////////////////
//
// constructor and destructor
//

TLinearIteration::TLinearIteration ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
        , _damping( 1.0 )
{}

TLinearIteration::TLinearIteration ( const double            damping,
                                     const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
        , _damping( damping )
{}

TLinearIteration::~TLinearIteration ()
{}

////////////////////////////////////////////////
//
// solving the system
//

//
// the real solving
//
template < typename value_t,
           typename value_pre_t >
void
TLinearIteration::solve ( const TLinearOperator< value_t > *      A,
                          TVector< value_t > *                    x,
                          const TVector< value_t > *              b,
                          const TLinearOperator< value_pre_t > *  W,
                          TSolverInfo *                           info ) const
{
    using  real_t = real_type_t< value_t >;
        
    if ( x == nullptr ) HERROR( ERR_ARG, "(TLinearIteration) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TLinearIteration) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TLinearIteration) solve", "A = nullptr" );
    
    real_t  norm, norm0;
    auto    r = b->copy();

    if ( info != nullptr )
        info->set_solver( "Lin. Iter." );
    
    set_start_value( x, b, W );

    // r = Ax - b
    apply( A, x, r.get() );
    r->axpy( value_t(-1), b );
    
    norm = norm0 = r->norm2();
    
    if ( info != nullptr )
        info->append( 0, norm );

    //
    // standard linear iteration
    //

    uint  it = 0;
    
    while ( true )
    {
        // iteration step : x_new = x_old - W r
        if ( W != nullptr )
        {
            apply_add( - _damping, *W, *r, *x );
        }// if
        else
        {
            x->axpy( -_damping, r.get() );
        }// else

        it++;

        // update residual
        apply( A, x, r.get() );
        r->axpy( value_t(-1), b );

        norm = r->norm2();

        if ( info != nullptr )
            info->append( it, norm );

        if ( stopped( it, norm, norm0, info ) )
            break;
    }// while
}

//
// the real solving
//
template < typename value_t >
void
TLinearIteration::solve ( const TLinearOperator< value_t > *  A,
                          TMatrix< value_t > *                X,
                          const TMatrix< value_t > *          B,
                          const TLinearOperator< value_t > *  W,
                          TSolverInfo *                       info ) const
{
    using  real_t = real_type_t< value_t >;

    if ( X == nullptr ) HERROR( ERR_ARG, "(TLinearIteration) solve", "x = nullptr" );
    if ( B == nullptr ) HERROR( ERR_ARG, "(TLinearIteration) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TLinearIteration) solve", "A = nullptr" );
    
    real_t  norm, norm0;
    auto    R = B->copy();

    if ( info != nullptr )
        info->set_solver( "Lin. Iter." );

    X->scale( value_t(0) );
    // set_start_value( X, B, W );

    // R = AX - B
    R->scale( value_t(-1) );
    apply_add( value_t(1), A, X, R.get() );

    norm = norm0 = norm_F( R.get() );
    
    if ( info != nullptr )
        info->append( 0, norm );

    //
    // standard linear iteration
    //

    uint  it = 0;
    
    while ( true )
    {
        // iteration step : x_new = x_old - W r
        if ( W != nullptr )
        {
            apply_add< value_t >( - _damping, W, R.get(), X );
        }// if
        else
        {
            add( value_t(-_damping), R.get(), value_t(1), X, acc_exact );
        }// else

        it++;

        // update residual
        B->copy_to( R.get() );
        R->scale( value_t(-1) );
        apply_add( value_t(1), A, X, R.get() );

        norm = norm_F( R.get() );

        if ( info != nullptr )
            info->append( it, norm );

        if ( stopped( it, norm, norm0, info ) )
            break;
    }// while
}

//
// generic implementation for "virtual" solve method
//
void
TLinearIteration::solve ( any_const_operator_t  A,
                          any_vector_t          x,
                          any_const_vector_t    b,
                          any_const_operator_t  W,
                          TSolverInfo *         info ) const
{
    solve_mp( *this, "TLinearIteration", A, x, b, W, info );
}

template < typename value_t > using opptr_t = TLinearOperator< value_t > *;

void
TLinearIteration::solve ( any_const_operator_t  A,
                          any_vector_t          x,
                          any_const_vector_t    b,
                          TSolverInfo *         info ) const
{
    using std::get;
    
    if (( A.index() != x.index() ) ||
        ( A.index() != b.index() ))
        HERROR( ERR_ARG, "(TLinearIteration) solve", "A, x and b have different value type" );

    switch ( A.index() )
    {
        case 0: this->solve( get< 0 >( A ), get< 0 >( x ), get< 0 >( b ), opptr_t< float >( nullptr ), info ); break;
        case 1: this->solve( get< 1 >( A ), get< 1 >( x ), get< 1 >( b ), opptr_t< double >( nullptr ), info ); break;
        case 2: this->solve( get< 2 >( A ), get< 2 >( x ), get< 2 >( b ), opptr_t< std::complex< float > >( nullptr ), info ); break;
        case 3: this->solve( get< 3 >( A ), get< 3 >( x ), get< 3 >( b ), opptr_t< std::complex< double > >( nullptr ), info ); break;

        default:
            HERROR( ERR_ARG, "(TLinearIteration) solve", "A, x and b have unsupported value type" );
    }// switch
}

// instantiate solve method
HPRO_INST_SOLVE_METHOD( TLinearIteration )

}// namespace Hpro
