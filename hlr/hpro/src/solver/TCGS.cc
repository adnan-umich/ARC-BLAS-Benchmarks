//
// Project     : HLIBpro
// File        : TCGS.cc
// Description : class implementing CGS
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/solver/TCGS.hh"

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
// actual CGS iteration
//
template < typename value_t,
           typename value_pre_t >
void
cgs ( const TLinearOperator< value_t > *      A,
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
        info->set_solver( "CGS" );
    
    solver->set_start_value( x, b, W );

    auto  r  = b->copy();
    auto  r0 = b->copy();
    auto  rt = b->copy(); // for exact residual

    // r0 = W( b - A x )
    apply_add( value_t(-1), A, x, r0.get() );

    if ( solver->use_exact_residual() )
        norm_0 = r0->norm2();

    if ( W != nullptr )
    {
        apply( W, r0.get(), r.get() );
        r0->assign( value_t(1), r.get() );
    }// if
    else
        r->assign( value_t(1), r0.get() );

    if ( ! solver->use_exact_residual() )
        norm_0 = r0->norm2();

    auto  p  = r->copy();
    auto  u  = r->copy();
    auto  q  = r->copy();
    auto  t  = r->copy();
    
    norm = norm_0;

    if ( info != nullptr )
        info->append( 0, norm );
    
    it = 0;
    while ( true )
    {
        // t = Ap
        if ( W != nullptr )
        {
            apply( A, p.get(), q.get() );
            apply( W, q.get(), t.get() );
        }// if
        else
        {
            apply( A, p.get(), t.get() );
        }// else

        // α = <rᵢ, r₀> / <Apᵢ, r₀>
        auto  sigma = dot( t.get(), r0.get() );
        auto  rho   = dot( r.get(), r0.get() );
        auto  alpha = rho / sigma;

        // q = u - αAp
        q->assign( 1.0, u.get() );
        axpy< value_t >( -alpha, t.get(), q.get() );

        // u = u + q
        axpy< value_t >( 1.0, q.get(), u.get() );
        
        // x = x + α(u + q)
        axpy< value_t >( alpha, u.get(), x );

        // r = r - α A (u + q)
        if ( W != nullptr )
        {
            apply( A, u.get(), t.get() );
            apply_add( -alpha, W, t.get(), r.get() );
        }// if
        else
        {
            apply_add( -alpha, A, u.get(), r.get() );
        }// else

        // one step per MVM, so two steps here
        it += 2;

        if ( solver->use_exact_residual() && ( W != nullptr ))
        {
            apply( A, x, rt.get() );
            rt->axpy( value_t(-1), b );
            norm = rt->norm2();
        }// if
        else
        {
            norm = r->norm2();
        }// else
        
        if ( info != nullptr )
            info->append( it, norm );
        
        if ( solver->stopped( it, norm, norm_0, info ) )
            break;

        
        auto  beta = dot( r.get(), r0.get() ) / rho;

        // u = r + βq
        u->assign( 1.0, r.get() );
        axpy< value_t >( beta, q.get(), u.get() );
        
        // p = u + βq + β²p
        scale< value_t >( beta*beta, p.get() );
        axpy< value_t >( beta, q.get(), p.get() );
        axpy< value_t >( value_t(1), u.get(), p.get() );
    }// while
}

}// namespace anonymous

////////////////////////////////////////////////
//
// constructor and destructor
//

TCGS::TCGS ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
{}

TCGS::~TCGS ()
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
TCGS::solve ( const TLinearOperator< value_t > *      A,
              TVector< value_t > *                    x,
              const TVector< value_t > *              b,
              const TLinearOperator< value_pre_t > *  W,
              TSolverInfo *                           info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TCGS) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TCGS) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TCGS) solve", "A = nullptr" );
    
    cgs( A, x, b, W, info, this );
}

//
// generic implementation for "virtual" solve method
//
void
TCGS::solve ( any_const_operator_t  A,
              any_vector_t          x,
              any_const_vector_t    b,
              any_const_operator_t  W,
              TSolverInfo *         info ) const
{
    solve_mp( *this, "TCGS", A, x, b, W, info );
}

template < typename value_t > using opptr_t = TLinearOperator< value_t > *;

void
TCGS::solve ( any_const_operator_t  A,
             any_vector_t          x,
             any_const_vector_t    b,
             TSolverInfo *         info ) const
{
    using std::get;
    
    if (( A.index() != x.index() ) ||
        ( A.index() != b.index() ))
        HERROR( ERR_ARG, "(TCGS) solve", "A, x and b have different value type" );

    switch ( A.index() )
    {
        case 0: this->solve( get< 0 >( A ), get< 0 >( x ), get< 0 >( b ), opptr_t< float >( nullptr ), info ); break;
        case 1: this->solve( get< 1 >( A ), get< 1 >( x ), get< 1 >( b ), opptr_t< double >( nullptr ), info ); break;
        case 2: this->solve( get< 2 >( A ), get< 2 >( x ), get< 2 >( b ), opptr_t< std::complex< float > >( nullptr ), info ); break;
        case 3: this->solve( get< 3 >( A ), get< 3 >( x ), get< 3 >( b ), opptr_t< std::complex< double > >( nullptr ), info ); break;

        default:
            HERROR( ERR_ARG, "(TCGS) solve", "A, x and b have unsupported value type" );
    }// switch
}

// instantiate solve method
HPRO_INST_SOLVE_METHOD( TCGS )

}// namespace Hpro
