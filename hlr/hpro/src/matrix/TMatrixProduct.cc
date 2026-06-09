//
// Project     : HLIBpro
// File        : TMatrixProduct.cc
// Description : Represents product of two matrices (linear ops)
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/matrix/TMatrixProduct.hh"

namespace Hpro
{

namespace
{

//
// return operator corresponding to op0( op1 )
//
matop_t
apply_op ( const matop_t  op0,
           const matop_t  op1 )
{
    switch ( op0 )
    {
        case apply_normal :
        {
            switch ( op1 ) 
            {
                case apply_normal     : return apply_normal;
                case apply_conjugate  : return apply_conjugate;
                case apply_transposed : return apply_transposed;
                case apply_adjoint    : return apply_adjoint;
            }// switch
        }
        break;

        case apply_conjugate :
        {
            switch ( op1 ) 
            {
                case apply_normal     : return apply_conjugate;
                case apply_conjugate  : return apply_normal;
                case apply_transposed : return apply_adjoint;
                case apply_adjoint    : return apply_transposed;
            }// switch
        }
        break;

        case apply_transposed :
        {
            switch ( op1 ) 
            {
                case apply_normal     : return apply_transposed;
                case apply_conjugate  : return apply_adjoint;
                case apply_transposed : return apply_normal;
                case apply_adjoint    : HERROR( ERR_ARG, "apply_op", "no operator for conjugate" ); return apply_conjugate;
            }// switch
        }
        break;

        case apply_adjoint :
        {
            switch ( op1 ) 
            {
                case apply_normal     : return apply_adjoint;
                case apply_conjugate  : return apply_transposed;
                case apply_transposed : HERROR( ERR_ARG, "apply_op", "no operator for conjugate" ); return apply_conjugate; 
                case apply_adjoint    : return apply_normal;
            }// switch
        }
        break;
        
    }// switch

    HERROR( ERR_CONSISTENCY, "apply_op", "not reachable" ); 
}

}// namespace anonymous

///////////////////////////////////////////////////////////
//
// ctor
//

template < typename value_t >
TMatrixProduct< value_t >::TMatrixProduct ( const value_t                       alpha0,
                                            const TLinearOperator< value_t > *  A0,
                                            const bool                          is_owner )
        : TLinearOperator< value_t >()
        , _is_owner( is_owner )
{
    if ( A0 == nullptr )
        HERROR( ERR_ARG, "(TMatrixProduct) ctor", "matrix is NULL" );

    _factors.push_back( { A0, apply_normal, alpha0 } );
}

template < typename value_t >
TMatrixProduct< value_t >::TMatrixProduct ( const value_t                       alpha0,
                                            const matop_t                       op0,
                                            const TLinearOperator< value_t > *  A0,
                                            const bool                          is_owner )
        : TLinearOperator< value_t >()
        , _is_owner( is_owner )
{
    if ( A0 == nullptr )
        HERROR( ERR_ARG, "(TMatrixProduct) ctor", "matrix is NULL" );

    _factors.push_back( { A0, op0, alpha0 } );
}

template < typename value_t >
TMatrixProduct< value_t >::TMatrixProduct ( const value_t                       alpha0,
                                            const TLinearOperator< value_t > *  A0,
                                            const value_t                       alpha1,
                                            const TLinearOperator< value_t > *  A1,
                                            const bool                          is_owner )
        : TLinearOperator< value_t >()
        , _is_owner( is_owner )
{
    if (( A0 == nullptr ) || ( A1 == nullptr ))
        HERROR( ERR_ARG, "(TMatrixProduct) ctor", "matrix is NULL" );

    _factors.push_back( { A0, apply_normal, alpha0 } );
    _factors.push_back( { A1, apply_normal, alpha1 } );
}

template < typename value_t >
TMatrixProduct< value_t >::TMatrixProduct ( const value_t                       alpha0,
                                            const matop_t                       op0,
                                            const TLinearOperator< value_t > *  A0,
                                            const value_t                       alpha1,
                                            const matop_t                       op1,
                                            const TLinearOperator< value_t > *  A1,
                                            const bool                          is_owner )
        : TLinearOperator< value_t >()
        , _is_owner( is_owner )
{
    if (( A0 == nullptr ) || ( A1 == nullptr ))
        HERROR( ERR_ARG, "(TMatrixProduct) ctor", "matrix is NULL" );

    _factors.push_back( { A0, op0, alpha0 } );
    _factors.push_back( { A1, op1, alpha1 } );
}


template < typename value_t >
TMatrixProduct< value_t >::TMatrixProduct ( const value_t                       alpha0,
                                            const TLinearOperator< value_t > *  A0,
                                            const value_t                       alpha1,
                                            const TLinearOperator< value_t > *  A1,
                                            const value_t                       alpha2,
                                            const TLinearOperator< value_t > *  A2,
                                            const bool                          is_owner )
        : TLinearOperator< value_t >()
        , _is_owner( is_owner )
{
    if (( A0 == nullptr ) || ( A1 == nullptr ) || ( A2 == nullptr ))
        HERROR( ERR_ARG, "(TMatrixProduct) ctor", "matrix is NULL" );

    _factors.push_back( { A0, apply_normal, alpha0 } );
    _factors.push_back( { A1, apply_normal, alpha1 } );
    _factors.push_back( { A2, apply_normal, alpha2 } );
}

template < typename value_t >
TMatrixProduct< value_t >::TMatrixProduct ( const value_t                       alpha0,
                                            const matop_t                       op0,
                                            const TLinearOperator< value_t > *  A0,
                                            const value_t                       alpha1,
                                            const matop_t                       op1,
                                            const TLinearOperator< value_t > *  A1,
                                            const value_t                       alpha2,
                                            const matop_t                       op2,
                                            const TLinearOperator< value_t > *  A2,
                                            const bool                          is_owner )
        : TLinearOperator< value_t >()
        , _is_owner( is_owner )
{
    if (( A0 == nullptr ) || ( A1 == nullptr ) || ( A2 == nullptr ))
        HERROR( ERR_ARG, "(TMatrixProduct) ctor", "matrix is NULL" );

    _factors.push_back( { A0, op0, alpha0 } );
    _factors.push_back( { A1, op1, alpha1 } );
    _factors.push_back( { A2, op2, alpha2 } );
}

template < typename value_t >
TMatrixProduct< value_t >::~TMatrixProduct ()
{
    if ( _is_owner )
    {
        const auto  N = _factors.size();
        
        for ( size_t  i = 0; i < N; ++i )
        {
            auto  ptr = _factors[i].linop;
            
            delete ptr;

            // make sure, that we have no double pointers!
            for ( size_t  j = i+1; j < N; ++j )
                if ( _factors[j].linop == ptr )
                    _factors[j].linop = nullptr;
        }// for
    }// if
}

///////////////////////////////////////////////////////////
//
// linear operator properties
//

//
// return true, if field type is complex
//
template < typename value_t >
bool
TMatrixProduct< value_t >::is_complex () const
{
    return _factors[0].linop->is_complex();
}
    
//
// return true, of operator is self adjoint
//
template < typename value_t >
bool
TMatrixProduct< value_t >::is_self_adjoint () const
{
    // only test trivial case of A₁ == A₂ == A₃ ...
    auto  T = _factors[0].linop;

    for ( auto &  f : _factors )
        if (( T != f.linop ) || ( f.op != apply_normal ))
            return false;

    return T->is_self_adjoint();
}
    
///////////////////////////////////////////////////////////
//
// linear operator mapping
//
    
//
// mapping function of linear operator \f$A\f$, e.g. \f$ y := A(x)\f$.
// Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
//
template < typename value_t >
void
TMatrixProduct< value_t >::apply ( const TVector< value_t > *  x,
                                   TVector< value_t > *        y,
                                   const matop_t               op ) const
{
    const auto  N = _factors.size();
    
    if ( N == 1 )
    {
        std::unique_ptr< TVector< value_t > >  tx, ty;

        ty = ( _factors[N-1].op == apply_normal ? _factors[N-1].linop->range_vector() : _factors[N-1].linop->domain_vector() );
        _factors[N-1].linop->apply( x, ty.get(), _factors[N-1].op );
        
        for ( int  i = int(N-2); i >= 1; --i )
        {
            tx = std::move( ty );
            ty = ( op == apply_normal ? _factors[i].linop->range_vector() : _factors[i].linop->domain_vector() );
            
            _factors[i].linop->apply( tx.get(), ty.get(), _factors[i].op );
        }// for
        
        _factors[0].linop->apply( ty.get(), y, _factors[0].op );
    }// if
    else
    {
        std::unique_ptr< TVector< value_t > >  tx, ty;
        
        ty = ( apply_op( op, _factors[0].op ) == apply_normal ? _factors[0].linop->range_vector() : _factors[0].linop->domain_vector() );
        _factors[0].linop->apply( x, ty.get(), apply_op( op, _factors[0].op ) );
        
        for ( size_t  i = 1; i < N-1; ++i )
        {
            tx = std::move( ty );
            ty = ( apply_op( op, _factors[i].op ) == apply_normal ? _factors[i].linop->range_vector() : _factors[i].linop->domain_vector() );
            
            _factors[i].linop->apply( tx.get(), ty.get(), apply_op( op, _factors[i].op ) );
        }// for
        
        _factors[N-1].linop->apply( ty.get(), y, apply_op( op, _factors[N-1].op ) );
    }// else
}

//
// mapping function with update: \f$ y := y + \alpha A(x)\f$.
// Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
//
template < typename value_t >
void
TMatrixProduct< value_t >::apply_add ( const value_t               alpha,
                                       const TVector< value_t > *  x,
                                       TVector< value_t > *        y,
                                       const matop_t               op ) const
{
    const auto  N = _factors.size();
    
    if ( N == 1 )
    {
        std::unique_ptr< TVector< value_t > >  tx, ty;

        ty = ( _factors[N-1].op == apply_normal ? _factors[N-1].linop->range_vector() : _factors[N-1].linop->domain_vector() );
        _factors[N-1].linop->apply( x, ty.get(), _factors[N-1].op );
        
        for ( int  i = int(N-2); i >= 1; --i )
        {
            tx = std::move( ty );
            ty = ( op == apply_normal ? _factors[i].linop->range_vector() : _factors[i].linop->domain_vector() );
            
            _factors[i].linop->apply( tx.get(), ty.get(), _factors[i].op );
        }// for

        _factors[0].linop->apply_add( alpha, ty.get(), y, _factors[0].op );
    }// if
    else
    {
        std::unique_ptr< TVector< value_t > >  tx, ty;

        ty = ( apply_op( op, _factors[0].op ) == apply_normal ? _factors[0].linop->range_vector() : _factors[0].linop->domain_vector() );
        _factors[0].linop->apply( x, ty.get(), apply_op( op, _factors[0].op ) );
        
        for ( size_t  i = 1; i < N-1; ++i )
        {
            tx = std::move( ty );
            ty = ( apply_op( op, _factors[i].op ) == apply_normal ? _factors[i].linop->range_vector() : _factors[i].linop->domain_vector() );
            
            _factors[i].linop->apply( tx.get(), ty.get(), apply_op( op, _factors[i].op ) );
        }// for

        _factors[N-1].linop->apply_add( alpha, ty.get(), y, apply_op( op, _factors[N-1].op ) );
    }// else
}

template < typename value_t >
void
TMatrixProduct< value_t >::apply_add   ( const value_t               , // alpha,
                                         const TMatrix< value_t > *  , // X,
                                         TMatrix< value_t > *        , // Y,
                                         const matop_t                 // op
                                         ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// same as above but only the dimension of the vector spaces is tested,
// not the corresponding index sets
//
template < typename value_t >
void
TMatrixProduct< value_t >::apply_add   ( const value_t                    alpha,
                                         const BLAS::Vector< value_t > &  x,
                                         BLAS::Vector< value_t > &        y,
                                         const matop_t                    op ) const
{
    const auto  N = _factors.size();

    if ( N == 1 )
    {
        _factors[0].linop->apply_add( _factors[0].scale * alpha, x, y, apply_op( op, _factors[0].op ) );
    }// if
    else
    {
        BLAS::Vector< value_t >  tx, ty;

        auto  init_ty =
            [this,&ty,op] ( const size_t  i )
            {
                const auto &  f = _factors[i];
                
                if ( apply_op( op, f.op ) == apply_normal )
                    ty = std::move( BLAS::Vector< value_t >( f.linop->range_dim() ) );
                else
                    ty = std::move( BLAS::Vector< value_t >( f.linop->domain_dim() ) );
            };
        
        auto  linop_apply =
            [this,op] ( const size_t                     i,
                        const BLAS::Vector< value_t > &  vx,
                        BLAS::Vector< value_t > &        vy)
            {
                const auto &  f = _factors[i];
                
                f.linop->apply_add( f.scale, vx, vy, apply_op( op, f.op ) );
            };
        
        if ( op == apply_normal )
        {
            init_ty( N-1 );
            linop_apply( N-1, x, ty );
        
            for ( int  i = int(N-2); i >= 1; --i )
            {
                tx = std::move( ty );
                init_ty( i );
                linop_apply( i, tx, ty );
            }// for

            _factors[0].linop->apply_add( _factors[0].scale * alpha, ty, y, _factors[0].op );
        }// if
        else
        {
            init_ty( 0 );
            linop_apply( 0, x, ty );
        
            for ( size_t  i = 1; i < N-1; ++i )
            {
                tx = std::move( ty );
                init_ty( i );
                linop_apply( i, tx, ty );
            }// for

            _factors[N-1].linop->apply_add( alpha*_factors[N-1].scale, ty, y, apply_op( op, _factors[N-1].op ) );
        }// else
    }// else
}

template < typename value_t >
void
TMatrixProduct< value_t >::apply_add   ( const value_t                    alpha,
                                         const BLAS::Matrix< value_t > &  X,
                                         BLAS::Matrix< value_t > &        Y,
                                         const matop_t                    op ) const
{
    const auto  N = _factors.size();

    if ( N == 1 )
    {
        _factors[0].linop->apply_add( _factors[0].scale * alpha, X, Y, apply_op( op, _factors[0].op ) );
    }// if
    else
    {
        auto  TX = BLAS::Matrix< value_t >();
        auto  TY = BLAS::Matrix< value_t >();
        
        auto  init_ty = [this,&X,&TY,op] ( const size_t  i )
        {
            const auto &  f = _factors[i];
            
            if ( apply_op( op, f.op ) == apply_normal ) TY = std::move( BLAS::Matrix< value_t >( f.linop->range_dim(), X.ncols() ) );
            else                                        TY = std::move( BLAS::Matrix< value_t >( f.linop->domain_dim(), X.ncols() ) );
        };
        
        auto  linop_apply = [this,op] ( const size_t                     i,
                                        const BLAS::Matrix< value_t > &  VX,
                                        BLAS::Matrix< value_t > &        VY )
        {
            const auto &  f = _factors[i];
            
            f.linop->apply_add( f.scale, VX, VY, apply_op( op, f.op ) );
        };
    
        if ( op == apply_normal )
        {
            init_ty( N-1 );
            linop_apply( N-1, X, TY );
        
            for ( int  i = int(N-2); i >= 1; --i )
            {
                TX = std::move( TY );
                init_ty( i );
                linop_apply( i, TX, TY );
            }// for

            _factors[0].linop->apply_add( _factors[0].scale * alpha, TY, Y, _factors[0].op );
        }// if
        else
        {
            init_ty( 0 );
            linop_apply( 0, X, TY );
        
            for ( size_t  i = 1; i < N-1; ++i )
            {
                TX = std::move( TY );
                init_ty( i );
                linop_apply( i, TX, TY );
            }// for

            _factors[N-1].linop->apply_add( alpha*_factors[N-1].scale, TY, Y, apply_op( op, _factors[N-1].op ) );
        }// else
    }// else
}

///////////////////////////////////////////////////////////
//
// access to vector space elements
//

//
// return dimension of domain
//
template < typename value_t >
size_t
TMatrixProduct< value_t >::domain_dim () const
{
    auto &  f = _factors.back();
    
    if ( f.op == apply_normal ) return f.linop->domain_dim();
    else                        return f.linop->range_dim();
}

//
// return dimension of range
//
template < typename value_t >
size_t
TMatrixProduct< value_t >::range_dim  () const
{
    auto &  f = _factors.front();
    
    if ( f.op == apply_normal ) return f.linop->range_dim();
    else                        return f.linop->domain_dim();
}

//
// return vector in domain space
//
template < typename value_t >
auto
TMatrixProduct< value_t >::domain_vector  () const -> std::unique_ptr< TVector< value_t > >
{
    auto &  f = _factors.back();
    
    if ( f.op == apply_normal ) return f.linop->domain_vector();
    else                        return f.linop->range_vector();
}

//
// return vector in range space
//
template < typename value_t >
auto
TMatrixProduct< value_t >::range_vector   () const -> std::unique_ptr< TVector< value_t > >
{
    auto &  f = _factors.front();
    
    if ( f.op == apply_normal ) return f.linop->range_vector();
    else                        return f.linop->domain_vector();
}

///////////////////////////////////////////////////////////
//
// explicit instantiation
//

template class TMatrixProduct< float >;
template class TMatrixProduct< double >;
template class TMatrixProduct< std::complex< float > >;
template class TMatrixProduct< std::complex< double > >;

}// namespace Hpro
