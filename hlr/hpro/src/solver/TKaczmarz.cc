//
// Project     : HLIBpro
// File        : TKaczmarz.cc
// Description : iterative solver based on Kaczmarz algorithm
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/solver/TKaczmarz.hh"

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

TKaczmarz::TKaczmarz ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
{}

TKaczmarz::~TKaczmarz ()
{}

////////////////////////////////////////////////
//
// solving the system
//

//
// the actual Kaczmarz iteration
// (see Hackbusch: "Large iterative Systems", Sec 8.2)
//
template < typename value_t,
           typename value_pre_t >
void
TKaczmarz::solve ( const TLinearOperator< value_t > *      A,
                   TVector< value_t > *                    x,
                   const TVector< value_t > *              b,
                   const TLinearOperator< value_pre_t > *  W,
                   TSolverInfo *                           info ) const
{
    using  real_t = real_type_t< value_t >;

    if ( x == nullptr ) HERROR( ERR_ARG, "(TKaczmarz) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TKaczmarz) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TKaczmarz) solve", "A = nullptr" );
    
    real_t                                 norm, norm0;
    std::unique_ptr< TVector< value_t > >  r( b->copy() );
    std::unique_ptr< TVector< value_t > >  a( b->copy() );
    std::unique_ptr< TVector< value_t > >  e( b->copy() );
    std::unique_ptr< TVector< value_t > >  t( b->copy() );

    if ( info != nullptr )
        info->set_solver( "Kaczmarz" );

    set_start_value( x, b, W );    

    // r = A x - b
    apply( A, x, t.get() );
    t->axpy( value_t(-1), b );
    if ( W != nullptr ) apply( W, t.get(), r.get() );
    else                r->assign( value_t(1), t.get() );

    norm = norm0 = r->norm2();

    if ( info != nullptr )
        info->append( 0, norm );

    e->scale( value_t(0) );

    //
    // Kaczmarz iteration
    //

    uint        it = 0;
    const uint  n  = uint(e->size());
    
    while ( ! stopped( it, norm, norm0, info ) )
    {
        // determine <it>'th column of A
        e->set_entry( it % n, value_t(1) );

#if 1
        //
        // left transformation
        //
        
        // a = A e
        apply( A, e.get(), t.get() );
        if ( W != nullptr ) apply( W, t.get(), a.get() );
        else                a->assign( value_t(1), t.get() );
        
        // <Ax-b, a>
        const value_t  r_dot_a = dot( r.get(), a.get() );

        // <a, a>
        const value_t  a_dot_a = dot( a.get(), a.get() );

        // if <a,a> == 0, then A is not regular
        if ( abs( a_dot_a ) < 1e-32 )
            HERROR( ERR_MAT_SINGULAR, "(TKaczmarz) solve",
                    to_string( "e_%d in kern(A)", it % n ) );
        
        // x = x - e <r,a>/<a,a>
        x->axpy( - r_dot_a / a_dot_a, e.get() );
#else
        //
        // right transformation
        //
        
        // a = A^H e
        apply( A, e.get(), a.get(), MATOP_ADJ );
        
        // <Ax-b, e>
        const value_t  r_dot_e = dot( r.get(), e.get() );

        // <a, a>
        const value_t  a_dot_a = dot( a.get(), a.get() );

        // if <a,a> == 0, then A is not regular
        if ( abs( a_dot_a ) < 1e-32 )
            HERROR( ERR_MAT_SINGULAR, "(TKaczmarz) solve",
                    to_string( "e_%d in kern(A)", it+1 ) );
        
        // x = x - a <r,e>/<a,a>
        x->axpy( - r_dot_e / a_dot_a, a.get() );
#endif
        
        e->set_entry( it % n, value_t(0) );

        // update residual
        apply( A, x, t.get() );
        t->axpy( value_t(-1), b );
        if ( W != nullptr ) apply( W, t.get(), r.get() );
        else                r->assign( value_t(1), t.get() );

        norm = r->norm2();
        
        // Sys::Out::printfln( "|r_%d| = %.8e", it+1, norm );

        it++;
        if ( info != nullptr )
            info->append( it+1, norm );
    }// while
}

//
// generic implementation for "virtual" solve method
//
void
TKaczmarz::solve ( any_const_operator_t  A,
                   any_vector_t          x,
                   any_const_vector_t    b,
                   any_const_operator_t  W,
                   TSolverInfo *         info ) const
{
    solve_mp( *this, "TKaczmarz", A, x, b, W, info );
}

template < typename value_t > using opptr_t = TLinearOperator< value_t > *;

void
TKaczmarz::solve ( any_const_operator_t  A,
                   any_vector_t          x,
                   any_const_vector_t    b,
                   TSolverInfo *         info ) const
{
    using std::get;
    
    if (( A.index() != x.index() ) ||
        ( A.index() != b.index() ))
        HERROR( ERR_ARG, "(TKaczmarz) solve", "A, x and b have different value type" );

    switch ( A.index() )
    {
        case 0: this->solve( get< 0 >( A ), get< 0 >( x ), get< 0 >( b ), opptr_t< float >( nullptr ), info ); break;
        case 1: this->solve( get< 1 >( A ), get< 1 >( x ), get< 1 >( b ), opptr_t< double >( nullptr ), info ); break;
        case 2: this->solve( get< 2 >( A ), get< 2 >( x ), get< 2 >( b ), opptr_t< std::complex< float > >( nullptr ), info ); break;
        case 3: this->solve( get< 3 >( A ), get< 3 >( x ), get< 3 >( b ), opptr_t< std::complex< double > >( nullptr ), info ); break;

        default:
            HERROR( ERR_ARG, "(TKaczmarz) solve", "A, x and b have unsupported value type" );
    }// switch
}

// instantiate solve method
HPRO_INST_SOLVE_METHOD( TKaczmarz )

}// namespace Hpro
