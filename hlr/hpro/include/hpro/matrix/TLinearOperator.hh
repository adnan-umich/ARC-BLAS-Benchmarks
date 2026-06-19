#ifndef __HPRO_TLINEAROPERATOR_HH
#define __HPRO_TLINEAROPERATOR_HH
//
// Project     : HLIBpro
// File        : TLinearOperator.hh
// Description : baseclass for all linear operators
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <variant>

#include "hpro/base/error.hh"
#include "hpro/base/TTypeInfo.hh"
#include "hpro/blas/types.hh"
#include "hpro/blas/Vector.hh"
#include "hpro/blas/Matrix.hh"
#include "hpro/vector/TVector.hh"
#include "hpro/vector/convert.hh"

namespace Hpro
{

// local RTTI types
DECLARE_TYPE( TLinearOperator );

template < typename value_t >
class TMatrix;

//!
//! \ingroup Matrix_Module
//! \class   TLinearOperator
//! \brief   Base class for all linear operators mapping vectors to vectors.
//!
//!          Many standard arithmetic operations only depend upon a linear operator
//!          providing the mapping between vectors, e.g. iterativ solvers. Instead of
//!          requiring a full matrix and hence the need for an implementation of the
//!          full matrix algebra, an object of type TLinearOperator is fully sufficient
//!          in such cases.
//!
template < typename T_value >
class TLinearOperator : public TTypeInfo
{
public:
    using  value_t = T_value;
    using  real_t  = typename real_type< value_t >::type_t;
    
public:
    ///////////////////////////////////////////////////////////
    //
    // dtor
    //

    virtual ~TLinearOperator () {}

    ///////////////////////////////////////////////////////////
    //
    // linear operator properties
    //

    //! return true, if field type is complex valued
    virtual bool  is_complex      () const { return is_complex_type< value_t >::value; }
    
    //! return true, if field type is real valued
    virtual bool  is_real         () const { return ! is_complex(); }
    
    //! return true, of operator is self adjoint
    virtual bool  is_self_adjoint () const = 0;
    
    ///////////////////////////////////////////////////////////
    //
    // linear operator mapping
    //
    
    //!
    //! mapping function of linear operator \f$A\f$, e.g. \f$ y := A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply       ( const TVector< value_t > *  x,
                                TVector< value_t > *        y,
                                const matop_t               op = apply_normal ) const = 0;

    //!
    //! mapping function with update: \f$ y := y + \alpha A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply_add   ( const value_t               alpha,
                                const TVector< value_t > *  x,
                                TVector< value_t > *        y,
                                const matop_t               op = apply_normal ) const = 0;

    virtual void  apply_add   ( const value_t               alpha,
                                const TMatrix< value_t > *  X,
                                TMatrix< value_t > *        Y,
                                const matop_t               op = apply_normal ) const = 0;
    
    //! same as above but only the dimension of the vector spaces is tested,
    //! not the corresponding index sets
    virtual void  apply_add   ( const value_t                    alpha,
                                const BLAS::Vector< value_t > &  x,
                                BLAS::Vector< value_t > &        y,
                                const matop_t                    op = apply_normal ) const = 0;

    virtual void  apply_add   ( const value_t                    alpha,
                                const BLAS::Matrix< value_t > &  X,
                                BLAS::Matrix< value_t > &        Y,
                                const matop_t                    op = apply_normal ) const = 0;

    ///////////////////////////////////////////////////////////
    //
    // access to vector space data
    //

    //! return dimension of domain
    virtual size_t  domain_dim () const = 0;
    
    //! return dimension of range
    virtual size_t  range_dim  () const = 0;
    
    //! return vector in domain space
    virtual auto    domain_vector  () const -> std::unique_ptr< TVector< value_t > > = 0;

    //! return vector in range space
    virtual auto    range_vector   () const -> std::unique_ptr< TVector< value_t > > = 0;

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //

    HPRO_RTTI_BASE( TLinearOperator );
};

//////////////////////////////////////////////////////////
//
// variant version of linear operator
//

using any_operator_t = std::variant<
    TLinearOperator< float > *,
    TLinearOperator< double > *,
    TLinearOperator< std::complex< float > > *,
    TLinearOperator< std::complex< double > > * >;
    
using any_const_operator_t = std::variant<
    const TLinearOperator< float > *,
    const TLinearOperator< double > *,
    const TLinearOperator< std::complex< float > > *,
    const TLinearOperator< std::complex< double > > * >;
    
//////////////////////////////////////////////////////////
//
// functional form of TLinearOperator methods
//

//!
//! \ingroup Matrix_Module
//! \fn      apply
//! \brief   Functional form of TLinearOperator::apply.
//!
template < typename value_mat_t,
           typename value_vec_t >
void
apply ( const TLinearOperator< value_mat_t > *  A,
        const TVector< value_vec_t > *          x,
        TVector< value_vec_t > *                y,
        const matop_t                           op = apply_normal )
{
    if (( A == nullptr ) || ( x == nullptr ) || ( y == nullptr ))
        HERROR( ERR_ARG, "apply", "argument is NULL" );

    if constexpr ( std::is_same< value_mat_t, value_vec_t >::value )
    {
        A->apply( x, y, op );
    }// if
    else
    {
        auto  tx = convert< value_mat_t, value_vec_t >( *x );
        auto  ty = convert< value_mat_t, value_vec_t >( *y );

        A->apply( tx.get(), ty.get(), op );

        convert_to( *ty, *y );
    }// else
}

template < typename value_mat_t,
           typename value_vec_t >
void
apply ( const TLinearOperator< value_mat_t > &  A,
        const TVector< value_vec_t > &          x,
        TVector< value_vec_t > &                y,
        const matop_t                           op = apply_normal )
{
    apply( &A, &x, &y, op );
}

template < typename value_mat_t,
           typename value_vec_t >
void
apply ( const TLinearOperator< value_mat_t > &  A,
        const BLAS::Vector< value_vec_t > &     x,
        BLAS::Vector< value_vec_t > &           y,
        const matop_t                           op = apply_normal )
{
    BLAS::fill( value_vec_t(0), y );
    A.apply_add( x, y, op );
}

//!
//! \ingroup Matrix_Module
//! \fn      apply
//! \brief   Functional form of TLinearOperator::apply.
//!
template < typename value_mat_t,
           typename value_vec_t,
           typename alpha_t >
void
apply_add ( const alpha_t                           alpha,
            const TLinearOperator< value_mat_t > *  A,
            const TVector< value_vec_t > *          x,
            TVector< value_vec_t > *                y,
            const matop_t                           op = apply_normal )
{
    if (( A == nullptr ) || ( x == nullptr ) || ( y == nullptr ))
        HERROR( ERR_ARG, "apply_add", "argument is NULL" );

    if constexpr ( std::is_same< value_mat_t, value_vec_t >::value )
    {
        A->apply_add( value_mat_t(alpha), x, y, op );
    }// if
    else
    {
        auto  tx = convert< value_mat_t, value_vec_t >( *x );
        auto  ty = convert< value_mat_t, value_vec_t >( *y );

        A->apply_add( value_mat_t(alpha), tx.get(), ty.get(), op );

        convert_to( *ty, *y );
    }// else
}

template < typename value_t,
           typename alpha_t >
void
apply_add ( const alpha_t                       alpha,
            const TLinearOperator< value_t > *  A,
            const TMatrix< value_t > *          X,
            TMatrix< value_t > *                Y,
            const matop_t                       op = apply_normal )
{
    if (( A == nullptr ) || ( X == nullptr ) || ( Y == nullptr ))
        HERROR( ERR_ARG, "apply_add", "argument is NULL" );

    A->apply_add( value_t(alpha), X, Y, op );
}

template < typename value_mat_t,
           typename value_vec_t,
           typename alpha_t >
void
apply_add ( const alpha_t                           alpha,
            const TLinearOperator< value_mat_t > &  A,
            const TVector< value_vec_t > &          x,
            TVector< value_vec_t > &                y,
            const matop_t                           op = apply_normal )
{
    apply_add( alpha, &A, &x, &y, op );
}

template < typename value_t,
           typename alpha_t >
void
apply_add ( const alpha_t                       alpha,
            const TLinearOperator< value_t > &  A,
            const TMatrix< value_t > &          X,
            TMatrix< value_t > &                Y,
            const matop_t                       op = apply_normal )
{
    A.apply_add( value_t(alpha), & X, & Y, op );
}

//! return dimension of domain
template < typename value_t > size_t  domain_dim ( const TLinearOperator< value_t > *  op ) { return op->domain_dim(); }
template < typename value_t > size_t  domain_dim ( const TLinearOperator< value_t > &  op ) { return op.domain_dim(); }
    
//! return dimension of range
template < typename value_t > size_t  range_dim  ( const TLinearOperator< value_t > *  op ) { return op->range_dim(); }
template < typename value_t > size_t  range_dim  ( const TLinearOperator< value_t > &  op ) { return op.range_dim(); }
    
//! return vector in domain space
template < typename value_t > std::unique_ptr< TVector< value_t > >  domain_vector  ( const TLinearOperator< value_t > *  op ) { return op->domain_vector(); }
template < typename value_t > std::unique_ptr< TVector< value_t > >  domain_vector  ( const TLinearOperator< value_t > &  op ) { return op.domain_vector(); }

//! return vector in range space
template < typename value_t > std::unique_ptr< TVector< value_t > >  range_vector   ( const TLinearOperator< value_t > *  op ) { return op->range_vector(); }
template < typename value_t > std::unique_ptr< TVector< value_t > >  range_vector   ( const TLinearOperator< value_t > &  op ) { return op.range_vector(); }

//////////////////////////////////////////////////////////
//
// debug helpers
//

namespace DBG
{

//!
//! write linear operator to file
//!
template < typename value_t >
void
write ( const TLinearOperator< value_t > *  M,
        const std::string &                 filename,
        const std::string &                 matname );

template < typename value_t >
void
write ( const TLinearOperator< value_t > &  M,
        const std::string &                 filename,
        const std::string &                 matname )
{
    write( &M, filename, matname );
}

}// namespace DBG

}// namespace Hpro

#endif  // __HPRO_TLINEAROPERATOR_HH
