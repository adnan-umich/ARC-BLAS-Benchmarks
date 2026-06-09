//
// Project     : HLIBpro
// File        : TMatrixSum.cc
// Description : Represents sum of two matrices (linear ops)
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/blas/Algebra.hh"
#include "hpro/matrix/TMatrixSum.hh"

namespace Hpro
{

// namespace
// {

// //
// // return operator corresponding to op0( op1 )
// //
// matop_t
// apply_op ( const matop_t  op0,
//            const matop_t  op1 )
// {
//     switch ( op0 )
//     {
//         case apply_normal :
//         {
//             switch ( op1 ) 
//             {
//                 case apply_normal     : return apply_normal;
//                 case apply_conjugate  : return apply_conjugate;
//                 case apply_transposed : return apply_transposed;
//                 case apply_adjoint    : return apply_adjoint;
//             }// switch
//         }
//         break;

//         case apply_conjugate :
//         {
//             switch ( op1 ) 
//             {
//                 case apply_normal     : return apply_conjugate;
//                 case apply_conjugate  : return apply_normal;
//                 case apply_transposed : return apply_adjoint;
//                 case apply_adjoint    : return apply_transposed;
//             }// switch
//         }
//         break;

//         case apply_transposed :
//         {
//             switch ( op1 ) 
//             {
//                 case apply_normal     : return apply_transposed;
//                 case apply_conjugate  : return apply_adjoint;
//                 case apply_transposed : return apply_normal;
//                 case apply_adjoint    : HERROR( ERR_ARG, "apply_op", "no operator for conjugate" ); return apply_conjugate;
//             }// switch
//         }
//         break;

//         case apply_adjoint :
//         {
//             switch ( op1 ) 
//             {
//                 case apply_normal     : return apply_adjoint;
//                 case apply_conjugate  : return apply_transposed;
//                 case apply_transposed : HERROR( ERR_ARG, "apply_op", "no operator for conjugate" ); return apply_conjugate; 
//                 case apply_adjoint    : return apply_normal;
//             }// switch
//         }
//         break;
        
//     }// switch

//     HERROR( ERR_CONSISTENCY, "apply_op", "not reachable" ); 
// }

// }// namespace anonymous

///////////////////////////////////////////////////////////
//
// ctor
//

template < typename value_t >
TMatrixSum< value_t >::TMatrixSum ( const value_t                       alpha0,
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
        HERROR( ERR_ARG, "(TMatrixSum) ctor", "matrix is NULL" );

    _summands.push_back( { A0, op0, alpha0 } );
    _summands.push_back( { A1, op1, alpha1 } );
}


template < typename value_t >
TMatrixSum< value_t >::TMatrixSum ( const value_t                       alpha0,
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
        HERROR( ERR_ARG, "(TMatrixSum) ctor", "matrix is NULL" );

    _summands.push_back( { A0, op0, alpha0 } );
    _summands.push_back( { A1, op1, alpha1 } );
    _summands.push_back( { A2, op2, alpha2 } );
}

template < typename value_t >
TMatrixSum< value_t >::~TMatrixSum ()
{
    if ( _is_owner )
    {
        for ( size_t  i = 0; i < _summands.size(); ++i )
        {
            auto  ptr = _summands[i].linop;
            
            delete ptr;

            // make sure, that we have no double pointers!
            for ( size_t  j = i+1; j < _summands.size(); ++j )
                if ( _summands[j].linop == ptr )
                    _summands[j].linop = nullptr;
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
TMatrixSum< value_t >::is_complex () const
{
    return _summands[0].linop->is_complex();
}
    
//
// return true, of operator is self adjoint
//
template < typename value_t >
bool
TMatrixSum< value_t >::is_self_adjoint () const
{
    // TODO: test scalars in complex case

    // only test trivial case of A₁ == A₂ == A₃ ...
    auto  T = _summands[0].linop;

    for ( auto &  s : _summands )
        if (( T != s.linop ) || ( s.op != apply_normal ))
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
TMatrixSum< value_t >::apply ( const TVector< value_t > *  x,
                               TVector< value_t > *        y,
                               const matop_t               op ) const
{
    y->fill( value_t(0) );
    
    for ( auto &  s : _summands )
        s.linop->apply_add( s.scale, x, y, BLAS::apply_op( op, s.op ) );
}

//
// mapping function with update: \f$ y := y + \alpha A(x)\f$.
// Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
//
template < typename value_t >
void
TMatrixSum< value_t >::apply_add ( const value_t               alpha,
                                   const TVector< value_t > *  x,
                                   TVector< value_t > *        y,
                                   const matop_t               op ) const
{
    for ( auto &  s : _summands )
        s.linop->apply_add( alpha*s.scale, x, y, BLAS::apply_op( op, s.op ) );
}

template < typename value_t >
void
TMatrixSum< value_t >::apply_add   ( const value_t               , // alpha,
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
TMatrixSum< value_t >::apply_add   ( const value_t                    alpha,
                                     const BLAS::Vector< value_t > &  x,
                                     BLAS::Vector< value_t > &        y,
                                     const matop_t                    op ) const
{
    BLAS::fill( value_t(0), y );
    
    for ( auto &  s : _summands )
        s.linop->apply_add( alpha*s.scale, x, y, BLAS::apply_op( op, s.op ) );
}

template < typename value_t >
void
TMatrixSum< value_t >::apply_add   ( const value_t                    alpha,
                                     const BLAS::Matrix< value_t > &  X,
                                     BLAS::Matrix< value_t > &        Y,
                                     const matop_t                    op ) const
{
    BLAS::fill( value_t(0), Y );
    
    for ( auto &  s : _summands )
        s.linop->apply_add( alpha*s.scale, X, Y, BLAS::apply_op( op, s.op ) );
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
TMatrixSum< value_t >::domain_dim () const
{
    const auto  s = _summands.front();
    
    if ( s.op == apply_normal ) return s.linop->domain_dim();
    else                        return s.linop->range_dim();
}

//
// return dimension of range
//
template < typename value_t >
size_t
TMatrixSum< value_t >::range_dim () const
{
    const auto  s = _summands.front();
    
    if ( s.op == apply_normal ) return s.linop->range_dim();
    else                        return s.linop->domain_dim();
}

//
// return vector in domain space
//
template < typename value_t >
auto
TMatrixSum< value_t >::domain_vector  () const -> std::unique_ptr< TVector< value_t > >
{
    const auto  s = _summands.front();

    if ( s.op == apply_normal ) return s.linop->domain_vector();
    else                        return s.linop->range_vector();
}

//
// return vector in range space
//
template < typename value_t >
auto
TMatrixSum< value_t >::range_vector   () const -> std::unique_ptr< TVector< value_t > >
{
    const auto  s = _summands.front();

    if ( s.op == apply_normal ) return s.linop->range_vector();
    else                        return s.linop->domain_vector();
}

///////////////////////////////////////////////////////////
//
// explicit instantiation
//

template class TMatrixSum< float >;
template class TMatrixSum< double >;
template class TMatrixSum< std::complex< float > >;
template class TMatrixSum< std::complex< double > >;

}// namespace Hpro
