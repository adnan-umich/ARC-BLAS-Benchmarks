//
// Project     : HLIBpro
// File        : TCG.cc
// Description : class implementing CG
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/solver/TCG.hh"

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
// actual CG iteration
//
template < typename value_t,
           typename value_pre_t >
void
cg ( const TLinearOperator< value_t > *      A,
     TVector< value_t > *                    x,
     const TVector< value_t > *              b,
     const TLinearOperator< value_pre_t > *  W,
     TSolverInfo *                           info,
     const TSolver *                         solver )
{
    using  real_t = real_type_t< value_t >;

    uint     it;
    value_t  lambda, rho, rho_new;
    real_t   norm, norm_0;

    if ( info != nullptr )
        info->set_solver( "CG" );
    
    solver->set_start_value( x, b, W );

    auto  r = b->copy();
    auto  p = x->copy();
    auto  a = x->copy();
    auto  q = x->copy();

    // r = b - A x
    apply_add( -1, *A, *x, *r );

    // p = W r
    if ( W != nullptr ) apply( *W, *r, *p );
    else                p->assign( value_t(1), r.get() );

    // rho = < r, p >
    rho  = dot( r.get(), p.get() );
    norm = norm_0 = r->norm2(); 

    if ( info != nullptr )
        info->append( 0, norm );
    
    it = 0;
    while ( ! solver->stopped( it, norm, norm_0, info ) )
    {
        // a = A p
        apply( *A, *p, *a );
        
        // lambda = rho / < a, p >
        lambda = rho / dot( a.get(), p.get() );

        // x = x + lambda p; r = r - lambda a
        axpy(  lambda, p.get(), x );
        axpy( -lambda, a.get(), r.get() );

        it++;

        // q = W r
        if ( W != nullptr ) apply( *W, *r, *q );
        else                q->assign( value_t(1), r.get() );
        
        // rho' = < q, r >
        rho_new = dot( q.get(), r.get() );

        // check stop criteria
        if ( solver->use_exact_residual() && ( W != nullptr ))
        {
            norm = r->norm2();
        }// if
        else
        {
            norm = Math::abs( Math::sqrt( rho_new ) );
        }// else

        // p = r + (rho' / rho) p
        scale( rho_new / rho, p.get() );
        axpy( value_t(1), q.get(), p.get() );

        rho = rho_new;

        if ( info != nullptr )
            info->append( it, norm );
    }// while
}

}// namespace anonymous

////////////////////////////////////////////////
//
// constructor and destructor
//

TCG::TCG ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
{}

TCG::~TCG ()
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
TCG::solve ( const TLinearOperator< value_t > *      A,
             TVector< value_t > *                    x,
             const TVector< value_t > *              b,
             const TLinearOperator< value_pre_t > *  W,
             TSolverInfo *                           info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TCG) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TCG) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TCG) solve", "A = nullptr" );
    
    if ( ! A->is_self_adjoint() )
        HWARNING( "(TCG) solve : operator not self adjoint" );

    if (( W != nullptr ) && ! W->is_self_adjoint() )
        HWARNING( "(TCG) solve : precondition operator not self adjoint" );

    cg( A, x, b, W, info, this );
}

//
// generic implementation for "virtual" solve method
//
void
TCG::solve ( any_const_operator_t  A,
             any_vector_t          x,
             any_const_vector_t    b,
             any_const_operator_t  W,
             TSolverInfo *         info ) const
{
    solve_mp( *this, "TCG", A, x, b, W, info );
}

template < typename value_t > using opptr_t = TLinearOperator< value_t > *;

void
TCG::solve ( any_const_operator_t  A,
             any_vector_t          x,
             any_const_vector_t    b,
             TSolverInfo *         info ) const
{
    using std::get;
    
    if (( A.index() != x.index() ) ||
        ( A.index() != b.index() ))
        HERROR( ERR_ARG, "(TCG) solve", "A, x and b have different value type" );

    switch ( A.index() )
    {
        case 0: this->solve( get< 0 >( A ), get< 0 >( x ), get< 0 >( b ), opptr_t< float >( nullptr ), info ); break;
        case 1: this->solve( get< 1 >( A ), get< 1 >( x ), get< 1 >( b ), opptr_t< double >( nullptr ), info ); break;
        case 2: this->solve( get< 2 >( A ), get< 2 >( x ), get< 2 >( b ), opptr_t< std::complex< float > >( nullptr ), info ); break;
        case 3: this->solve( get< 3 >( A ), get< 3 >( x ), get< 3 >( b ), opptr_t< std::complex< double > >( nullptr ), info ); break;

        default:
            HERROR( ERR_ARG, "(TCG) solve", "A, x and b have unsupported value type" );
    }// switch
}

// instantiate solve method
HPRO_INST_SOLVE_METHOD( TCG )

}// namespace Hpro
