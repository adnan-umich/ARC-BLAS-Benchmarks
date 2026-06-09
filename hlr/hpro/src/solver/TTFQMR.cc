//
// Project     : HLIBpro
// File        : TTFQMR.cc
// Description : class implementing TFQMR
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/solver/TTFQMR.hh"

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

namespace
{

//
// actual TFQMR iteration
//
template < typename value_t,
           typename value_pre_t >
void
tfqmr ( const TLinearOperator< value_t > *      A,
        TVector< value_t > *                    x,
        const TVector< value_t > *              b,
        const TLinearOperator< value_pre_t > *  W,
        TSolverInfo *                           info,
        const TSolver *                         solver )
{
    using  real_t = real_type_t< value_t >;

    uint    it;
    real_t  norm   = 0;
    real_t  norm_0 = 0;

    if ( info != nullptr )
        info->set_solver( "TFQMR" );
    
    solver->set_start_value( x, b, W );

    auto  t  = b->copy();
    auto  t2 = b->copy();
    auto  rr = b->copy(); // for exact residual
    
    // r0 = W( b - A x )
    auto  r0 = b->copy();

    if ( W != nullptr )
    {
        apply_add( value_t(-1), A, x, r0.get() );

        if ( solver->use_exact_residual() )
            norm = norm_0 = r0->norm2();
        
        apply( W, r0.get(), t.get() );
        r0->assign( value_t(1), t.get() );

        if ( ! solver->use_exact_residual() )
            norm = norm_0 = r0->norm2();
    }// if
    else
    {
        apply_add( value_t(-1), A, x, r0.get() );
        norm = norm_0 = r0->norm2();
    }// else

    auto  w  = r0->copy();
    auto  y  = r0->copy();
    auto  d  = x->copy();

    // v = Ay
    auto  v  = x->copy();

    if ( W != nullptr )
    {
        apply( A, y.get(), d.get() );
        apply( W, d.get(), v.get() );
    }// if
    else
    {
        apply( A, y.get(), v.get() );
    }// else

    t->assign( value_t(1), v.get() );
    
    d->scale( 0 );
    
    real_t   tau   = r0->norm2();
    real_t   theta = 0;
    value_t  eta   = 0;
    value_t  rho   = dot( r0.get(), r0.get() );

    if ( info != nullptr )
        info->append( 0, norm );

    bool  converged = solver->stopped( 0, norm, norm_0, info );
    
    it = 0;
    while ( ! converged )
    {
        auto  sigma = dot( r0.get(), v.get() );

        if ( sigma == value_t(0) )
        {
            if ( info != nullptr )
                info->set_status( failed );
            else
                HERROR( ERR_DIV_ZERO, "(TTFQMR) solve", "breakdown" );
            
            return;
        }// if
        
        auto  alpha = rho / sigma;

        if ( alpha == value_t(0) )
        {
            if ( info != nullptr )
                info->set_status( failed );
            else
                HERROR( ERR_DIV_ZERO, "(TTFQMR) solve", "breakdown" );
            
            return;
        }// if
            
        for ( uint  m = 1; m <= 2; ++m )
        {
            // w = w - α A y
            if ( m == 1 )
            {
                // reuse computed product from last step
            }// if
            else
            {
                // compute Ay (also reused below)
                if ( W != nullptr )
                {
                    apply( A, y.get(), t2.get() );
                    apply( W, t2.get(), t.get() );
                }// if
                else
                {
                    apply( A, y.get(), t.get() );
                }// else
            }// else

            axpy< value_t >( -alpha, t.get(), w.get() );
            
            auto  theta_old = theta;
            auto  eta_old   = eta;

            if ( tau == real_t(0) )
            {
                if ( info != nullptr )
                    info->set_status( failed );
                else
                    HERROR( ERR_DIV_ZERO, "(TTFQMR) solve", "breakdown" );
                
                return;
            }// if
            
            theta   = w->norm2() / tau;
            auto  c = real_t(1) / std::sqrt( 1 + theta * theta );
            tau     = tau * theta * c;
            eta     = c*c * alpha;

            // d = y + (θ²η/α)d
            scale< value_t >( Math::square( theta_old ) * eta_old / alpha, d.get() );
            d->axpy( value_t(1), y.get() );

            // x = x + η d
            axpy< value_t >( eta, d.get(), x );

            // norm (approximation)
            if ( solver->use_exact_residual() )
            {
                apply( A, x, rr.get() );
                rr->axpy( value_t(-1), b );
                norm = rr->norm2();
            }// if
            else
            {
                norm = std::sqrt( 2*it + m ) * tau;
            }// else

            if ( info != nullptr )
                info->append( 2*it + m, norm );
            
            if ( solver->stopped( 2*it + m, norm, norm_0, info ) )
            {
                converged = true;
                break;
            }// if

            // update y after first step
            if ( m == 1 )
            {
                // y = y - αv
                axpy< value_t >( -alpha, v.get(), y.get() );
            }// if
        }// for

        if ( rho == value_t(0) )
        {
            if ( info != nullptr )
                info->set_status( failed );
            else
                HERROR( ERR_DIV_ZERO, "(TTFQMR) solve", "breakdown" );
            
            return;
        }// if
        
        auto  rho_new = dot( r0.get(), w.get() );
        auto  beta    = rho_new / rho;

        rho = rho_new;
        
        // y' = w + βy
        scale< value_t >( beta, y.get() );
        axpy< value_t >( value_t(1), w.get(), y.get() );
        
        // v = Ay' + β( Ay + β v ) = Ay' + βAy + β²v
        scale< value_t >( beta*beta, v.get() );
        axpy< value_t >( beta, t.get(), v.get() );

        if ( W != nullptr )
        {
            apply( A, y.get(), t2.get() );
            apply( W, t2.get(), t.get() );
        }// if
        else
        {
            apply( A, y.get(), t.get() );
        }// else

        axpy< value_t >( value_t(1), t.get(), v.get() );
        
        ++it;
    }// while
}

}// namespace anonymous

////////////////////////////////////////////////
//
// constructor and destructor
//

TTFQMR::TTFQMR ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
{}

TTFQMR::~TTFQMR ()
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
TTFQMR::solve ( const TLinearOperator< value_t > *     A,
                TVector< value_t > *                   x,
                const TVector< value_t > *             b,
                const TLinearOperator< value_pre_t > * W,
                TSolverInfo *                          info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TTFQMR) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TTFQMR) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TTFQMR) solve", "A = nullptr" );
    
    tfqmr( A, x, b, W, info, this );
}

//
// generic implementation for "virtual" solve method
//
void
TTFQMR::solve ( any_const_operator_t  A,
                any_vector_t          x,
                any_const_vector_t    b,
                any_const_operator_t  W,
                TSolverInfo *         info ) const
{
    solve_mp( *this, "TTFQMR", A, x, b, W, info );
}

template < typename value_t > using opptr_t = TLinearOperator< value_t > *;

void
TTFQMR::solve ( any_const_operator_t  A,
             any_vector_t          x,
             any_const_vector_t    b,
             TSolverInfo *         info ) const
{
    using std::get;
    
    if (( A.index() != x.index() ) ||
        ( A.index() != b.index() ))
        HERROR( ERR_ARG, "(TTFQMR) solve", "A, x and b have different value type" );

    switch ( A.index() )
    {
        case 0: this->solve( get< 0 >( A ), get< 0 >( x ), get< 0 >( b ), opptr_t< float >( nullptr ), info ); break;
        case 1: this->solve( get< 1 >( A ), get< 1 >( x ), get< 1 >( b ), opptr_t< double >( nullptr ), info ); break;
        case 2: this->solve( get< 2 >( A ), get< 2 >( x ), get< 2 >( b ), opptr_t< std::complex< float > >( nullptr ), info ); break;
        case 3: this->solve( get< 3 >( A ), get< 3 >( x ), get< 3 >( b ), opptr_t< std::complex< double > >( nullptr ), info ); break;

        default:
            HERROR( ERR_ARG, "(TTFQMR) solve", "A, x and b have unsupported value type" );
    }// switch
}

// instantiate solve method
HPRO_INST_SOLVE_METHOD( TTFQMR )

}// namespace Hpro
