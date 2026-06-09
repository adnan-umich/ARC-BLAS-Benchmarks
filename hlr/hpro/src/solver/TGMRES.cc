//
// Project     : HLIBpro
// File        : TGMRES.cc
// Description : class implementing GMRES
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/base/System.hh"

#include "hpro/solver/TGMRES.hh"

namespace Hpro
{

using std::unique_ptr;

// namespace abbr.
namespace B = BLAS;

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
// abbreviation for applying preconditioner or just copy vector
//
template < typename value_mat_t,
           typename value_vec_t >
void
gmres_apply ( const TLinearOperator< value_mat_t > *  M,
              const TVector< value_vec_t > &          v_src,
              TVector< value_vec_t > &                v_dest )
{
    if ( M != nullptr ) apply( M, & v_src, & v_dest );
    else                v_dest.assign( value_vec_t(1), & v_src );
}

//
// epsilon for check of (almost) linear dependency in Krylov subspace
//
template < typename value_t >
constexpr real_type_t< value_t > gmres_eps ()
{
    return Limits::epsilon<real_type_t< value_t >>() * real_type_t< value_t >(5e1);
}

////////////////////////////////////////////////
//
// apply plane rotation (  cs  sn ) to (a)
//                      ( -sn  cs )    (b)
//

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
// the real solving
//
template < typename value_t,
           typename value_pre_t >
void
gmres ( const uint                              restart,
        const TLinearOperator< value_t > *      A,
        TVector< value_t > *                    x,
        const TVector< value_t > *              b,
        const TLinearOperator< value_pre_t > *  W,
        TSolverInfo *                           info,
        const TSolver *                         solver )
{
    using  real_t = real_type_t< value_t >;

    uint                  j, it;
    real_t                norm = 0, norm_0 = 0, norm_old = 0;
    auto                  d = b->copy();
    B::Vector< value_t >  y( restart ), e1( restart+1 );
    B::Vector< value_t >  givc( restart+1 ), givs( restart+1 );
    B::Matrix< value_t >  H( restart+1, restart );
    auto                  V = std::vector< unique_ptr< TVector< value_t > > >( restart+1 );

    if ( info != nullptr )
        info->set_solver( "GMRES(" + to_string( "%d", restart ) + ")" );
    
    solver->set_start_value( x, b, W );

    //
    // GMRES iteration
    //
    
    int  ndiv = 0;     // counter for divergent inner iterations
        
    it = 0;
    
    while ( true )
    {
        //
        // v_0 = r / | r | = W ( b - Ax ) / |...|
        //
        
        if ( V[0] == nullptr )
            V[0] = b->copy();
        
        d->assign( value_t(1), b );
        apply_add( value_t(-1), A, x, d.get() );

        if ( solver->use_exact_residual() )
        {
            norm = d->norm2();
        }// if
        
        gmres_apply( W, *d, *V[0] );

        if ( ! solver->use_exact_residual() )
        {
            norm = V[0]->norm2();
        }// if

        if ( it == 0 )
        {
            norm_0 = norm_old = norm;

            if ( info != nullptr )
                info->append( 0, norm );
        }// if
        
        // test stop condition
        if ( solver->stopped( it, norm, norm_0, info ) )
            break;

        // check for divergent inner iterations
        if ( norm > norm_old ) ndiv++;
        else                   ndiv = 0;

        // stop after 10 divergent inner iterations
        if ( ndiv > 10 )
        {
            HWARNING( "in (TGMRES) solve : residual continuously growing (10 inner it.)" );
            break;
        }// if

        auto  norm_V0 = norm;
        
        // adjust norm_V0 if neccessary
        if ( solver->use_exact_residual() && ( W != nullptr ))
            norm_V0 = V[0]->norm2();
        
        if ( norm_V0 != real_t(0) ) scale( value_t(1) / norm_V0, V[0].get() );
        else                        break;

        //
        // e1 = | r | (1,0,0,...)^T
        //
        
        B::fill( value_t(0), e1 );
        e1(0) = norm_V0;
        
        //
        // iterate in Krylov-subspace
        //

        uint  dim = restart;
        
        for ( j = 0; j < restart; j++ )
        {
            //
            // v_j+1 = W A v_j
            //
            
            if ( V[j+1] == nullptr )
                V[j+1] = b->copy();
        
            gmres_apply( A, *V[j], *d );
            gmres_apply( W, *d,    *V[j+1] );

            norm = V[j+1]->norm2();
            
            //
            // orthogonalize wrt. previously computed basis-vectors
            //
            
            for ( uint i = 0; i <= j; i++ )
            {
                H(i,j) = dot( V[i].get(), V[j+1].get() );
                axpy( -H(i,j), V[i].get(), V[j+1].get() );
            }// for
            
            H(j+1,j) = V[j+1]->norm2();

            //
            // check if new direction is in span of previous directions
            // and stop inner iteration if it is the case
            //
            
            if ( Math::abs( H(j+1,j) ) < norm * gmres_eps< value_t >() )
            {
                dim = j+1;
                it++;
                break;
            }// if
            
            scale( value_t(1) / H(j+1,j), V[j+1].get() );
            
            //
            // apply old givens rotations to new column
            //
            
            for ( uint i = 0; i < j; i++ )
                apply_rot( givc(i), givs(i), H(i,j), H(i+1,j) );

            //
            // compute and apply givens rotation for new element
            //

            Math::givens( H(j,j), H(j+1,j), givc(j), givs(j) );
            apply_rot( givc(j), givs(j), H(j,j), H(j+1,j) );
            H(j+1,j) = value_t(0);

            //
            // apply givens rotation to RHS (new pos is 0 in e1)
            //
            
            e1(j+1) = - Math::conj( givs(j) ) * e1(j);
            e1(j)   =   givc(j)               * e1(j);


            if ( solver->use_exact_residual() )
            {
                //
                // solve min_y | beta e^1 - H y |
                // (H is upper triangular)
                //

                const uint      ldim = j+1;
                B::Matrix< value_t >  H_dim(  H,  B::Range( 0, ldim-1 ), B::Range( 0, ldim-1 ), copy_value );
                B::Vector< value_t >  e1_dim( e1, B::Range( 0, ldim-1 ), copy_value );
                
                B::solve_tri( B::upper_triangular, B::general_diag, H_dim, e1_dim );

                auto  x_copy = x->copy();
                
                for ( uint i = 0; i < ldim; i++ )
                    axpy( e1_dim(i), V[i].get(), x_copy.get() );

                gmres_apply( A, *x_copy, *d );
                d->axpy( value_t(-1), b );
                
                norm = d->norm2();
            }// if
            else
            {
                // last coefficient is actual residual norm
                norm = Math::abs( e1(j+1) );
            }// else

            it++;
            
            if ( info != nullptr )
                info->append( it, norm );
            
            //
            // test stop condition
            //

            if ( solver->stopped( it, norm, norm_0, info ) )
            {
                dim = j+1;
                break;
            }// if

            norm_old = norm;
        }// for

        //
        // solve min_y | beta e^1 - H y |
        // (H is upper triangular)
        //

        B::Matrix< value_t >  H_dim(  H,  B::Range( 0, dim-1 ), B::Range( 0, dim-1 ) );
        B::Vector< value_t >  e1_dim( e1, B::Range( 0, dim-1 ) );
        
        B::solve_tri( B::upper_triangular, B::general_diag, H_dim, e1_dim );
        
        //
        // adjust x
        //

        for ( uint i = 0; i < dim; i++ )
            axpy( e1(i), V[i].get(), x );

        // test stop condition
        if ( solver->stopped( it, norm, norm_0, info ) )
            break;
    }// while
}

}// namespace anonymous

////////////////////////////////////////////////
//
// constructor and destructor
//

TGMRES::TGMRES ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
        , _restart( CFG::Solver::gmres_restart )
{}

TGMRES::TGMRES ( const uint              restart,
                 const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
        , _restart( restart )
{}

TGMRES::~TGMRES ()
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
TGMRES::solve ( const TLinearOperator< value_t > *      A,
                TVector< value_t > *                    x,
                const TVector< value_t > *              b,
                const TLinearOperator< value_pre_t > *  W,
                TSolverInfo *                           info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TGMRES) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TGMRES) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TGMRES) solve", "A = nullptr" );

    gmres( _restart, A, x, b, W, info, this );
}

//
// generic implementation for "virtual" solve method
//
void
TGMRES::solve ( any_const_operator_t  A,
                any_vector_t          x,
                any_const_vector_t    b,
                any_const_operator_t  W,
                TSolverInfo *         info ) const
{
    solve_mp( *this, "TGMRES", A, x, b, W, info );
}

template < typename value_t > using opptr_t = TLinearOperator< value_t > *;

void
TGMRES::solve ( any_const_operator_t  A,
                any_vector_t          x,
                any_const_vector_t    b,
                TSolverInfo *         info ) const
{
    using std::get;
    
    if (( A.index() != x.index() ) ||
        ( A.index() != b.index() ))
        HERROR( ERR_ARG, "(TGMRES) solve", "A, x and b have different value type" );

    switch ( A.index() )
    {
        case 0: this->solve( get< 0 >( A ), get< 0 >( x ), get< 0 >( b ), opptr_t< float >( nullptr ), info ); break;
        case 1: this->solve( get< 1 >( A ), get< 1 >( x ), get< 1 >( b ), opptr_t< double >( nullptr ), info ); break;
        case 2: this->solve( get< 2 >( A ), get< 2 >( x ), get< 2 >( b ), opptr_t< std::complex< float > >( nullptr ), info ); break;
        case 3: this->solve( get< 3 >( A ), get< 3 >( x ), get< 3 >( b ), opptr_t< std::complex< double > >( nullptr ), info ); break;

        default:
            HERROR( ERR_ARG, "(TGMRES) solve", "A, x and b have unsupported value type" );
    }// switch
}

// instantiate solve method
HPRO_INST_SOLVE_METHOD( TGMRES )

}// namespace Hpro
