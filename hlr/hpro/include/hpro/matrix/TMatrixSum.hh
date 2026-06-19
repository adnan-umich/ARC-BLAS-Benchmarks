#ifndef __HPRO_TMATRIXSUM_HH
#define __HPRO_TMATRIXSUM_HH
//
// Project     : HLIBpro
// File        : TMatrixSum.hh
// Description : Represents sum of n matrices (linear ops)
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <deque>

#include "hpro/matrix/TLinearOperator.hh"

namespace Hpro
{

// local RTTI types
DECLARE_TYPE( TMatrixSum );

//!
//! \ingroup Matrix_Module
//! \class   TMatrixSum
//! \brief   Represents sum α₁A₁ + α₂A₂ + α₃A₃... of matrices (linear ops)
//!
template < typename T_value >
class TMatrixSum : public TLinearOperator< T_value >
{
    //! @cond
public:
 
    using  value_t = T_value;
   
private:

    struct summand_t {
        const TLinearOperator< value_t > *  linop;
        const matop_t                       op;
        const value_t                       scale;
    };
    
    // summands in matrix sum
    std::deque< summand_t >  _summands;

    // if true, local object has ownership
    const bool               _is_owner;
    
    //! @endcond
public:
    ///////////////////////////////////////////////////////////
    //
    // dtor
    //

    //! two matrices
    TMatrixSum ( const value_t                       alpha1,
                 const TLinearOperator< value_t > *  A1,
                 const value_t                       alpha2,
                 const TLinearOperator< value_t > *  A2,
                 const bool                          is_owner = false );
    
    TMatrixSum ( const value_t                       alpha1,
                 const matop_t                       op1,
                 const TLinearOperator< value_t > *  A1,
                 const value_t                       alpha2,
                 const matop_t                       op2,
                 const TLinearOperator< value_t > *  A2,
                 const bool                          is_owner = false );
    

    //! three matrices
    TMatrixSum ( const value_t                       alpha1,
                 const TLinearOperator< value_t > *  A1,
                 const value_t                       alpha2,
                 const TLinearOperator< value_t > *  A2,
                 const value_t                       alpha3,
                 const TLinearOperator< value_t > *  A3,
                 const bool                          is_owner = false );
    
    TMatrixSum ( const value_t                       alpha1,
                 const matop_t                       op1,
                 const TLinearOperator< value_t > *  A1,
                 const value_t                       alpha2,
                 const matop_t                       op2,
                 const TLinearOperator< value_t > *  A2,
                 const value_t                       alpha3,
                 const matop_t                       op3,
                 const TLinearOperator< value_t > *  A3,
                 const bool                          is_owner = false );

    virtual ~TMatrixSum ();

    ///////////////////////////////////////////////////////////
    //
    // linear operator properties
    //

    //! return true, if field type is complex
    virtual bool  is_complex      () const;
    
    //! return true, of operator is self adjoint
    virtual bool  is_self_adjoint () const;
    
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
                                const matop_t               op = apply_normal ) const;

    //!
    //! mapping function with update: \f$ y := y + \alpha A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply_add   ( const value_t               alpha,
                                const TVector< value_t > *  x,
                                TVector< value_t > *        y,
                                const matop_t               op = apply_normal ) const;

    virtual void  apply_add   ( const value_t               alpha,
                                const TMatrix< value_t > *  X,
                                TMatrix< value_t > *        Y,
                                const matop_t               op = apply_normal ) const;
    
    //! same as above but only the dimension of the vector spaces is tested,
    //! not the corresponding index sets
    virtual void  apply_add   ( const value_t                    alpha,
                                const BLAS::Vector< value_t > &  x,
                                BLAS::Vector< value_t > &        y,
                                const matop_t                    op = apply_normal ) const;

    virtual void  apply_add   ( const value_t                    alpha,
                                const BLAS::Matrix< value_t > &  X,
                                BLAS::Matrix< value_t > &        Y,
                                const matop_t                    op = apply_normal ) const;

    ///////////////////////////////////////////////////////////
    //
    // access to vector space elements
    //

    //! return dimension of domain
    virtual size_t  domain_dim () const;
    
    //! return dimension of range
    virtual size_t  range_dim  () const;
    
    //! return vector in domain space
    virtual auto    domain_vector  () const -> std::unique_ptr< TVector< value_t > >;

    //! return vector in range space
    virtual auto    range_vector   () const -> std::unique_ptr< TVector< value_t > >;

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //

    HPRO_RTTI_DERIVED( TMatrixSum, TLinearOperator< value_t > );
};

//
// functions to return matrix sum objects
//

template < typename value_t >
std::unique_ptr< TMatrixSum< value_t > >
matrix_sum ( const TLinearOperator< value_t > *  A0,
             const TLinearOperator< value_t > *  A1,
             const bool                          is_owner = false )
{
    return std::make_unique< TMatrixSum< value_t > >( value_t(1), apply_normal, A0,
                                                      value_t(1), apply_normal, A1,
                                                      is_owner );
}
template < typename value_t >
std::unique_ptr< TMatrixSum< value_t > >
matrix_sum ( const TLinearOperator< value_t > &  A0,
             const TLinearOperator< value_t > &  A1,
             const bool                          is_owner = false )
{
    return std::make_unique< TMatrixSum< value_t > >( value_t(1), apply_normal, &A0,
                                                      value_t(1), apply_normal, &A1,
                                                      is_owner );
}
                 
template < typename value_t,
           typename alpha0_t,
           typename alpha1_t >
std::unique_ptr< TMatrixSum< value_t > >
matrix_sum ( const alpha0_t                      alpha0,
             const TLinearOperator< value_t > *  A0,
             const alpha1_t                      alpha1,
             const TLinearOperator< value_t > *  A1,
             const bool                          is_owner = false )
{
    return std::make_unique< TMatrixSum< value_t > >( value_t(alpha0), apply_normal, A0,
                                                      value_t(alpha1), apply_normal, A1,
                                                      is_owner );
}
template < typename value_t,
           typename alpha0_t,
           typename alpha1_t >
std::unique_ptr< TMatrixSum< value_t > >
matrix_sum ( const alpha0_t                      alpha0,
             const TLinearOperator< value_t > &  A0,
             const alpha1_t                      alpha1,
             const TLinearOperator< value_t > &  A1,
             const bool                          is_owner = false )
{
    return std::make_unique< TMatrixSum< value_t > >( value_t(alpha0), apply_normal, &A0,
                                                      value_t(alpha1), apply_normal, &A1,
                                                      is_owner );
}
                 
template < typename value_t,
           typename alpha0_t,
           typename alpha1_t >
std::unique_ptr< TMatrixSum< value_t > >
matrix_sum ( const alpha0_t                      alpha0,
             const matop_t                       op0,
             const TLinearOperator< value_t > *  A0,
             const alpha1_t                      alpha1,
             const matop_t                       op1,
             const TLinearOperator< value_t > *  A1,
             const bool                          is_owner = false )
{
    return std::make_unique< TMatrixSum< value_t > >( value_t(alpha0), op0, A0,
                                                      value_t(alpha1), op1, A1,
                                                      is_owner );
}
template < typename value_t,
           typename alpha0_t,
           typename alpha1_t >
std::unique_ptr< TMatrixSum< value_t > >
matrix_sum ( const alpha0_t                      alpha0,
             const matop_t                       op0,
             const TLinearOperator< value_t > &  A0,
             const alpha1_t                      alpha1,
             const matop_t                       op1,
             const TLinearOperator< value_t > &  A1,
             const bool                          is_owner = false )
{
    return std::make_unique< TMatrixSum< value_t > >( value_t(alpha0), op0, &A0,
                                                      value_t(alpha1), op1, &A1,
                                                      is_owner );
}
    

template < typename value_t >
std::unique_ptr< TMatrixSum< value_t > >
matrix_sum ( const TLinearOperator< value_t > *  A0,
             const TLinearOperator< value_t > *  A1,
             const TLinearOperator< value_t > *  A2,
             const bool                          is_owner = false )
{
    return std::make_unique< TMatrixSum< value_t > >( value_t(1), apply_normal, A0,
                                                      value_t(1), apply_normal, A1,
                                                      value_t(1), apply_normal, A2,
                                                      is_owner );
}
template < typename value_t >
std::unique_ptr< TMatrixSum< value_t > >
matrix_sum ( const TLinearOperator< value_t > &  A0,
             const TLinearOperator< value_t > &  A1,
             const TLinearOperator< value_t > &  A2,
             const bool                          is_owner = false )
{
    return std::make_unique< TMatrixSum< value_t > >( value_t(1), apply_normal, &A0,
                                                      value_t(1), apply_normal, &A1,
                                                      value_t(1), apply_normal, &A2,
                                                      is_owner );
}
    
template < typename value_t,
           typename alpha0_t,
           typename alpha1_t,
           typename alpha2_t >
std::unique_ptr< TMatrixSum< value_t > >
matrix_sum ( const alpha0_t                      alpha0,
             const TLinearOperator< value_t > *  A0,
             const alpha1_t                      alpha1,
             const TLinearOperator< value_t > *  A1,
             const alpha2_t                      alpha2,
             const TLinearOperator< value_t > *  A2,
             const bool                          is_owner = false )
{
    return std::make_unique< TMatrixSum< value_t > >( value_t(alpha0), apply_normal, A0,
                                                      value_t(alpha1), apply_normal, A1,
                                                      value_t(alpha2), apply_normal, A2,
                                                      is_owner );
}
template < typename value_t,
           typename alpha0_t,
           typename alpha1_t,
           typename alpha2_t >
std::unique_ptr< TMatrixSum< value_t > >
matrix_sum ( const alpha0_t                      alpha0,
             const TLinearOperator< value_t > &  A0,
             const alpha1_t                      alpha1,
             const TLinearOperator< value_t > &  A1,
             const alpha2_t                      alpha2,
             const TLinearOperator< value_t > &  A2,
             const bool                          is_owner = false )
{
    return std::make_unique< TMatrixSum< value_t > >( value_t(alpha0), apply_normal, &A0,
                                                      value_t(alpha1), apply_normal, &A1,
                                                      value_t(alpha2), apply_normal, &A2,
                                                      is_owner );
}
                 
template < typename value_t,
           typename alpha0_t,
           typename alpha1_t,
           typename alpha2_t >
std::unique_ptr< TMatrixSum< value_t > >
matrix_sum ( const alpha0_t                      alpha0,
             const matop_t                       op0,
             const TLinearOperator< value_t > *  A0,
             const alpha1_t                      alpha1,
             const matop_t                       op1,
             const TLinearOperator< value_t > *  A1,
             const alpha2_t                      alpha2,
             const matop_t                       op2,
             const TLinearOperator< value_t > *  A2,
             const bool                          is_owner = false )
{
    return std::make_unique< TMatrixSum< value_t > >( value_t(alpha0), op0, A0,
                                                      value_t(alpha1), op1, A1,
                                                      value_t(alpha2), op2, A2,
                                                      is_owner );
}
template < typename value_t,
           typename alpha0_t,
           typename alpha1_t,
           typename alpha2_t >
std::unique_ptr< TMatrixSum< value_t > >
matrix_sum ( const alpha0_t                      alpha0,
             const matop_t                       op0,
             const TLinearOperator< value_t > &  A0,
             const alpha1_t                      alpha1,
             const matop_t                       op1,
             const TLinearOperator< value_t > &  A1,
             const alpha2_t                      alpha2,
             const matop_t                       op2,
             const TLinearOperator< value_t > &  A2,
             const bool                          is_owner = false )
{
    return std::make_unique< TMatrixSum< value_t > >( value_t(alpha0), op0, &A0,
                                                      value_t(alpha1), op1, &A1,
                                                      value_t(alpha2), op2, &A2,
                                                      is_owner );
}
                 
}// namespace Hpro

#endif  // __HPRO_TMATRIXSUM_HH
