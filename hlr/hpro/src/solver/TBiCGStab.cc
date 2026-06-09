//
// Project     : HLIBpro
// File        : TBiCGStab.cc
// Description : class implementing BiCG-Stab
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/System.hh"

#include "hpro/solver/TBiCGStab.hh"

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
// actual BiCG-Stab algorithm
//
template < typename value_t,
           typename value_pre_t >
void
bicgstab ( const TLinearOperator< value_t > *      A,
           TVector< value_t > *                    x,
           const TVector< value_t > *              b,
           const TLinearOperator< value_pre_t > *  W,
           TSolverInfo *                           info,
           const TSolver *                         solver )
{
    using  real_t = real_type_t< value_t >;
    
    const real_t  zero_limit = Math::square( Limits::epsilon< real_t >() );
    
    if ( info != nullptr )
        info->set_solver( "BiCGStab" );
    
    solver->set_start_value( x, b, W );

    // r = b - Ax
    auto  r  = b->copy();
    auto  t  = b->copy(); // for exact residual

    auto  norm0 = real_t(0);

    if ( W != nullptr )
    {
        apply_add( value_t(-1), A, x, t.get() );

        if ( solver->use_exact_residual() )
            norm0 = t->norm2();
        
        apply( W, t.get(), r.get() );

        if ( ! solver->use_exact_residual() )
            norm0 = r->norm2();
    }// if
    else
    {
        apply_add( value_t(-1), A, x, r.get() );
        norm0 = r->norm2();
    }// else

    auto  norm = norm0;

    // choose r0, with <r0,r> != 0
    auto  r0   = r->copy();

    // create/set up remaining aux. vectors
    auto  p    = r0->copy();
    auto  Ap   = p->copy();
    auto  Ar   = r->copy();
    
    //
    // iterate until stop-criterion reached
    //

    if ( info != nullptr )
        info->append( 0, norm0 );

    if ( solver->stopped( 0, norm, norm0, info ) )
        return;
    
    uint  it = 1;
    
    while ( true )
    {
        // α = <r_j, r_0> / <Ap_j, r_0>
        if ( W != nullptr )
        {
            apply( A, p.get(), t.get() );
            apply( W, t.get(), Ap.get() );
        }// if
        else
        {
            apply( A, p.get(), Ap.get() );
        }// else

        auto  r_r0  = dot( r.get(),  r0.get() );
        auto  Ap_r0 = dot( Ap.get(), r0.get() );
        
        if (( Math::abs( Ap_r0 ) < zero_limit ) || ( Math::abs( r_r0 ) < zero_limit ))
        {
            if ( info != nullptr ) info->set_status( failed );
            else                   HERROR( ERR_DIV_ZERO, "(TBiCGStab) solve", "breakdown" );
            
            return;
        }// if
        
        auto  alpha = r_r0 / Ap_r0;

        // x = x + αp
        axpy< value_t >( alpha, p.get(), x );

        // r = r - alpha A p_j
        axpy< value_t >( -alpha, Ap.get(), r.get() );

        if ( solver->use_exact_residual() && ( W != nullptr ))
        {
            apply( A, x, t.get() );
            t->axpy( value_t(-1), b );
            norm = t->norm2();
        }// if
        else
        {
            norm = r->norm2();
        }// else
        
        if ( info != nullptr )
            info->append( it, norm );
        
        if ( solver->stopped( it, norm, norm0, info ) )
            break;

        it++;

        // ω = <Ar_j, s_j> / <Ar_j, Ar_j>
        if ( W != nullptr )
        {
            apply( A, r.get(), t.get() );
            apply( W, t.get(), Ar.get() );
        }// if
        else
        {
            apply( A, r.get(), Ar.get() );
        }// else

        auto  Ar_r  = dot( Ar.get(), r.get() );
        auto  Ar_Ar = dot( Ar.get(), Ar.get() );
        auto  omega = Ar_r / Ar_Ar;

        if ( Math::abs( omega ) < zero_limit )
        {
            if ( info != nullptr ) info->set_status( failed );
            else                   HERROR( ERR_DIV_ZERO, "(TBiCGStab) solve", "breakdown" );
            
            return;
        }// if
        
        // x = x + ωr
        axpy< value_t >( omega, r.get(), x );

        // r = r - ωAr
        axpy< value_t >( - omega, Ar.get(), r.get() );

        if ( solver->use_exact_residual() && ( W != nullptr ))
        {
            apply( A, x, t.get() );
            t->axpy( value_t(-1), b );
            norm = t->norm2();
        }// if
        else
        {
            norm = r->norm2();
        }// else
        
        if ( info != nullptr )

            info->append( it, norm );
        
        if ( solver->stopped( it, norm, norm0, info ) )
            break;

        it++;

        // β = <r_j+1,r0> / <r_j,r0> · alpha_j / omega_j
        auto  beta = ( dot( r.get(), r0.get() ) / r_r0 ) * ( alpha / omega );

        // p = r_j+1 + β( p_j - ωAp_j )
        scale< value_t >( beta, p.get() );
        axpy< value_t >( -omega * beta, Ap.get(), p.get() );
        axpy< value_t >( value_t(1), r.get(), p.get() );
    }// while
}

}// namespace anonymous

////////////////////////////////////////////////
//
// constructor and destructor
//

TBiCGStab::TBiCGStab ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
{}

TBiCGStab::~TBiCGStab ()
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
TBiCGStab::solve ( const TLinearOperator< value_t > *     A,
                   TVector< value_t > *                   x,
                   const TVector< value_t > *             b,
                   const TLinearOperator< value_pre_t > * W,
                   TSolverInfo *                          info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TBiCGStab) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TBiCGStab) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TBiCGStab) solve", "A = nullptr" );

    bicgstab( A, x, b, W, info, this );
}

//
// generic implementation for "virtual" solve method
//
void
TBiCGStab::solve ( any_const_operator_t  A,
                   any_vector_t          x,
                   any_const_vector_t    b,
                   any_const_operator_t  W,
                   TSolverInfo *         info ) const
{
    solve_mp( *this, "TBiCGStab", A, x, b, W, info );
}

template < typename value_t > using opptr_t = TLinearOperator< value_t > *;

void
TBiCGStab::solve ( any_const_operator_t  A,
             any_vector_t          x,
             any_const_vector_t    b,
             TSolverInfo *         info ) const
{
    using std::get;
    
    if (( A.index() != x.index() ) ||
        ( A.index() != b.index() ))
        HERROR( ERR_ARG, "(TBiCGStab) solve", "A, x and b have different value type" );

    switch ( A.index() )
    {
        case 0: this->solve( get< 0 >( A ), get< 0 >( x ), get< 0 >( b ), opptr_t< float >( nullptr ), info ); break;
        case 1: this->solve( get< 1 >( A ), get< 1 >( x ), get< 1 >( b ), opptr_t< double >( nullptr ), info ); break;
        case 2: this->solve( get< 2 >( A ), get< 2 >( x ), get< 2 >( b ), opptr_t< std::complex< float > >( nullptr ), info ); break;
        case 3: this->solve( get< 3 >( A ), get< 3 >( x ), get< 3 >( b ), opptr_t< std::complex< double > >( nullptr ), info ); break;

        default:
            HERROR( ERR_ARG, "(TBiCGStab) solve", "A, x and b have unsupported value type" );
    }// switch
}

// instantiate solve method
HPRO_INST_SOLVE_METHOD( TBiCGStab )

}// namespace Hpro
