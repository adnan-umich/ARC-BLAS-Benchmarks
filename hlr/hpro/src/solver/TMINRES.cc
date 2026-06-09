//
// Project     : HLIBpro
// File        : TMINRES.cc
// Description : class implementing MINRES
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/solver/TMINRES.hh"

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
// epsilon for check of (almost) linear dependency in Krylov subspace
//
template < typename value_t >
constexpr real_type_t< value_t > minres_eps ()
{
    return Limits::epsilon<real_type_t< value_t >>() * real_type_t< value_t >(5e1);
}

//
// apply plane rotation (  cs   sn ) to (a)
//                      ( -sn'  cs )    (b)
//
template < typename value_t >
void
apply_rot ( const value_t  cs,
            const value_t  sn,
            value_t &      a,
            value_t &      b )
{
    const value_t  ta = a;
    const value_t  tb = b;

    a =   cs             * ta + sn * tb;
    b = - Math::conj(sn) * ta + cs * tb;
}

//
// MINRES algorithm
//
template < typename value_t,
           typename value_pre_t >
void
minres ( const TLinearOperator< value_t > *       A,
         TVector< value_t > *                     x,
         const TVector< value_t > *               b,
         const TLinearOperator< value_pre_t > *   W,
         TSolverInfo *                            info,
         const TSolver *                          solver )
{
    using  real_t = real_type_t< value_t >;
    
    solver->set_start_value( x, b, W );

    auto  p    = b->copy();
    auto  p1   = b->copy();
    auto  p2   = b->copy();
    auto  v    = b->copy();
    auto  v1   = b->copy();
    auto  vbar = b->copy();
    auto  rr   = b->copy(); // for exact residual

    //
    // v_0 = r / | r | = W ( Ax - b ) / |...|
    //
        
    real_t   norm  = 0;
    real_t   norm0 = 0;
    value_t  e;
    
    if ( W != nullptr )
    {
        apply( *A, *x, *p );
        p->axpy( value_t(-1), b );

        if ( solver->use_exact_residual() )
            norm0 = p->norm2();
        
        apply( *W, *p, *v );

        if ( ! solver->use_exact_residual() )
            norm0 = v->norm2();
    }// if
    else
    {
        apply( *A, *x, *v );
        v->axpy( value_t(-1), b );
        norm0 = v->norm2();
    }// else

    norm = norm0;

    if ( info != nullptr )
        info->append( 0, norm );
        
    // test stop condition
    if ( solver->stopped( 0, norm, norm0, info ) )
        return;

    auto  norm_v = norm;

    // adjust norm_v if neccessary
    if ( solver->use_exact_residual() && ( W != nullptr ))
        norm_v = v->norm2();
    
    v->scale( value_t(1) / norm_v );

    e = norm_v;
        
    //
    // MINRES iteration
    //
    
    uint  it  = 0;
    value_t     cs1 = value_t(0), sn1 = value_t(0);
    value_t     cs2 = value_t(0), sn2 = value_t(0);

    while ( ! solver->stopped( it, norm, norm0, info ) )
    {
        //
        // v_j+1 = W A v_j
        //
            
        if ( W != nullptr )
        {
            apply( *A, *v, *p  );
            apply( *W, *p, *vbar );
        }// if
        else
            apply( *A, *v, *vbar );

        norm = vbar->norm2();
            
        //
        // compute scalar products with previous directions
        //

        value_t  H;
        value_t  H1 = value_t(0);
        value_t  H2 = value_t(0);

        if ( it > 0 )
            H1 = dot( v1.get(), vbar.get() );
        H = dot( v.get(), vbar.get() );

        //
        // orthogonalize wrt. previously computed basis-vectors
        //
        
        if ( it > 0 ) axpy( -H1, v1.get(), vbar.get() );
        axpy( -H, v.get(), vbar.get() );
        
        value_t  Hbar = vbar->norm2();

        //
        // check if new direction is in span of previous directions
        // and stop inner iteration if so
        //
        
        if ( Math::abs( Hbar ) < norm * minres_eps< value_t >() )
            break;
        
        scale( value_t(1) / Hbar, vbar.get() );
            
        //
        // apply old Givens rotations to new column
        //
        
        if ( it > 1 ) apply_rot( cs2, sn2, H2, H1 );
        if ( it > 0 ) apply_rot( cs1, sn1, H1, H  );
            
        //
        // compute and apply Givens rotation for new element
        //

        value_t  cs, sn;
        
        Math::givens( H, Hbar, cs, sn );
        apply_rot( cs, sn, H, Hbar );

        //
        // apply Givens rotation to RHS (new pos is 0 in e1)
        //

        value_t  ebar;

        ebar = - sn * e;
        e    =   cs * e;

        //
        // adjust x
        //
        
        p->assign( value_t(1), v.get() );
        if ( it > 0 ) axpy( - H1, p1.get(), p.get() );
        if ( it > 1 ) axpy( - H2, p2.get(), p.get() );
        scale( value_t(1) / H, p.get() );

        axpy( -e, p.get(), x );

        if ( solver->use_exact_residual() && ( W != nullptr ))
        {
            apply( *A, *x, *rr );
            rr->axpy( value_t(-1), b );
            norm = rr->norm2();
        }// if
        else
        {
            // last coefficient is actual residual norm
            norm = Math::abs( ebar );
        }// else
                
        //
        // copy values for next step
        //
        
        v1->assign( value_t(1), v.get()  );
        v->assign(  value_t(1), vbar.get() );
        p2->assign( value_t(1), p1.get() );
        p1->assign( value_t(1), p.get()  );

        cs2 = cs1;
        sn2 = sn1;
        cs1 = cs;
        sn1 = sn;

        e  = ebar;

        it++;
        
        if ( info != nullptr )
            info->append( it, norm );
    }// for
}

}// namespace anonymous

////////////////////////////////////////////////
//
// constructor and destructor
//

TMINRES::TMINRES ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
{}

TMINRES::~TMINRES ()
{}

////////////////////////////////////////////////
//
// solving the system
//

//
// MINRES iteration
// - symmetric and optimized version if GMRES
//
template < typename value_t,
           typename value_pre_t >
void
TMINRES::solve ( const TLinearOperator< value_t > *      A,
                 TVector< value_t > *                    x,
                 const TVector< value_t > *              b,
                 const TLinearOperator< value_pre_t > *  W,
                 TSolverInfo *                           info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TMINRES) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TMINRES) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TMINRES) solve", "A = nullptr" );
    
    if ( info != nullptr )
        info->set_solver( "MINRES" );

    if ( ! A->is_self_adjoint() )
        HWARNING( "(TMINRES) solve : operator not self adjoint" );

    if (( W != nullptr ) && ! W->is_self_adjoint() )
        HWARNING( "(TMINRES) solve : precondition operator not self adjoint" );

    minres( A, x, b, W, info, this );
}

//
// generic implementation for "virtual" solve method
//
void
TMINRES::solve ( any_const_operator_t  A,
                 any_vector_t          x,
                 any_const_vector_t    b,
                 any_const_operator_t  W,
                 TSolverInfo *         info ) const
{
    solve_mp( *this, "TMINRES", A, x, b, W, info );
}

template < typename value_t > using opptr_t = TLinearOperator< value_t > *;

void
TMINRES::solve ( any_const_operator_t  A,
                 any_vector_t          x,
                 any_const_vector_t    b,
                 TSolverInfo *         info ) const
{
    using std::get;
    
    if (( A.index() != x.index() ) ||
        ( A.index() != b.index() ))
        HERROR( ERR_ARG, "(TMINRES) solve", "A, x and b have different value type" );

    switch ( A.index() )
    {
        case 0: this->solve( get< 0 >( A ), get< 0 >( x ), get< 0 >( b ), opptr_t< float >( nullptr ), info ); break;
        case 1: this->solve( get< 1 >( A ), get< 1 >( x ), get< 1 >( b ), opptr_t< double >( nullptr ), info ); break;
        case 2: this->solve( get< 2 >( A ), get< 2 >( x ), get< 2 >( b ), opptr_t< std::complex< float > >( nullptr ), info ); break;
        case 3: this->solve( get< 3 >( A ), get< 3 >( x ), get< 3 >( b ), opptr_t< std::complex< double > >( nullptr ), info ); break;

        default:
            HERROR( ERR_ARG, "(TMINRES) solve", "A, x and b have unsupported value type" );
    }// switch
}

// instantiate solve method
HPRO_INST_SOLVE_METHOD( TMINRES )

}// namespace Hpro
