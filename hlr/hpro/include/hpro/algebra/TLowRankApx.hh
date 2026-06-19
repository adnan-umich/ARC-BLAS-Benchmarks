#ifndef __HPRO_TLOWRANKAPX_HH
#define __HPRO_TLOWRANKAPX_HH
//
// Project     : HLib
// File        : TLowRankApx.hh
// Description : classes for computing a low rank approximation of a dense matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <deque>
#include <atomic>

#include "hpro/base/error.hh"
#include "hpro/base/TPoint.hh"
#include "hpro/blas/Matrix.hh"
#include "hpro/cluster/TCluster.hh"
#include "hpro/cluster/TGeomCluster.hh"
#include "hpro/matrix/TCoeffFn.hh"
#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"
#include "hpro/parallel/TMutex.hh"

namespace Hpro
{
    
//!
//! \{
//! \name Low Rank Approximation
//!       Classes and functions related to low rank approximation of dense matrices.
//!

//!
//! \ingroup Algebra_Module
//! \class   TLowRankApx
//! \brief   base class for all low rank approximation techniques
//!
template < typename T_value >
class TLowRankApx
{
public:

    using value_t = T_value;
    using real_t  = real_type_t< value_t >;
    
    //
    // statistical data of computations
    //
    struct stat_t
    {
        // number of requested matrix coefficients
        std::atomic< size_t >  ncoeff;

        // ctor nullifying data
        stat_t ()
                : ncoeff(0)
        {}
    };
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TLowRankApx () {}

    virtual ~TLowRankApx () {}

    //////////////////////////////////////
    //
    // build low-rank matrix
    //

    //! build low rank matrix for block cluster \a bct with
    //! rank defined by accuracy \a acc
    virtual
    std::unique_ptr< TMatrix< value_t > >
    build ( const TBlockCluster *   bct,
            const TTruncAcc &       acc ) const
    {
        if ( bct == nullptr )
            HERROR( ERR_ARG, "(TLowRankApx) build", "block cluster is null" );
        
        const TIndexSet  rowis( * ( bct->rowcl() ) );
        const TIndexSet  colis( * ( bct->colcl() ) );

        return build( bis( rowis, colis ), acc );
    }

    //! build low rank matrix for block index set \a bis with
    //! rank defined by accuracy \a acc
    virtual
    std::unique_ptr< TMatrix< value_t > >
    build ( const TBlockIndexSet &  bis,
            const TTruncAcc &       acc ) const = 0;

    //! indicate if algorithm provides statistics
    virtual bool  has_statistics () const { return false; }
};

//!
//! \ingroup Algebra_Module
//! \class   TZeroLRApx
//! \brief   Approximate all low-rank blocks by zero, e.g. for nearfield only.
//!
//!          If only the near field blocks of a matrix should be computed,
//!          TZeroLRApx will do that by approximating all far field blocks
//!          by zero, e.g. a 0-rank low rank matrix.
//!
template < typename T_value >
class TZeroLRApx : public TLowRankApx< T_value >
{
public:
    using value_t = T_value;
    using real_t  = real_type_t< value_t >;

public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TZeroLRApx () {}

    virtual ~TZeroLRApx () {}
    
    //////////////////////////////////////
    //
    // build low-rank matrix
    //

    //! return low rank matrix of rank 0
    virtual
    std::unique_ptr< TMatrix< value_t > >
    build ( const TBlockIndexSet &  bis,
            const TTruncAcc &       acc ) const;
    using TLowRankApx< value_t >::build;
};
    
//!
//! \ingroup Algebra_Module
//! \class   TDenseLRApx
//! \brief   Computes dense matrix block without approximation.
//!
//!          Instead of performing approximation for a matrix block, the whole block
//!          is computed and returned as a dense matrix.
//!
//!          This is usually used for debugging or accuracy tests.
//!
template < typename T_value >
class TDenseLRApx : public TLowRankApx< T_value >
{
public:
    //
    // template type as public member type
    //

    using  value_t    = T_value;
    using  real_t     = real_type_t< value_t >;
    using  coeff_fn_t = TCoeffFn< value_t >;
    
protected:
    // function return coefficient for index-pair
    const coeff_fn_t *  _coeff_fn;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TDenseLRApx ( const coeff_fn_t *  coeff_fn );

    virtual ~TDenseLRApx () {}
    
    //////////////////////////////////////
    //
    // build low-rank matrix
    //

    //! build low rank matrix for block index set \a bis with
    //! rank defined by accuracy \a acc
    virtual
    std::unique_ptr< TMatrix< value_t > >
    build ( const TBlockIndexSet &  bis,
            const TTruncAcc &       acc ) const;
    using TLowRankApx< value_t >::build;
};

//!
//! \ingroup Algebra_Module
//! \class   THCA
//! \brief   uses hybrid cross approximation (HCA) for computing low rank approximation
//!
//!          THCA provides a low rank approximation algorithm with a guaranteed
//!          approximation quality. It is based on the <em>generator function</em> \f$\gamma(x,y)\f$
//!          of a BEM kernel function \f$k(x,y)\f$ and it's derivatives \f$ D_x \gamma(x, y_{l}) \f$
//!          and \f$ D_y \gamma(x, y_{l}) \f$.
//!
//!          The class THCA needs a user implemented generator function of type TGeneratorFn, in which
//!          the function itself and the integrals of the corresponding derivates are defined.
//!
//!          Furthermore, HCA is based on interpolation, of which the order defines the accuracy of
//!          the final result. This interpolation order together with an accuracy for the approximation
//!          of the generator function is specific to the given problem, and hence, user defined.
//!
template < typename T_value >
class THCA : public TLowRankApx< T_value >
{
public:
    //
    // template arguments as internal types
    //

    using  value_t    = T_value;
    using  real_t     = real_type_t< value_t >;
    using  coeff_fn_t = TCoeffFn< value_t >;

    //!
    //! \class  TGeneratorFn
    //! \brief  class defining kernel generator function used by HCA
    //!
    class TGeneratorFn
    {
    public:
        //
        // constructor and destructor
        //

        TGeneratorFn () {}
        
        virtual ~TGeneratorFn () {}

        //! indicate complex nature of function
        virtual bool     is_complex    () const { return is_complex_type< value_t >::value; }
        
        //!
        //! Evaluate generator function γ at (\a x, \a y).
        //!
        virtual value_t  eval          ( const T3Point &  x,
                                         const T3Point &  y ) const = 0;

        //!
        //! Evaluate \f$ \int \phi_i(x) D_x \gamma(x, y_{l}) dx\f$ for i ∈ \a is
        //! and points \f$ y_{l} \f$ defined by \a pts.
        //! Store results in \a matrix at index (i,l).
        //!
        virtual void     integrate_dx  ( const TIndexSet &               is,
                                         const std::vector< T3Point > &  pts,
                                         BLAS::Matrix< value_t > &       matrix ) const = 0;
        //!
        //! Evaluate \f$ \int \phi_j(y) D_y \gamma(x_{l}, y) dy\f$ for j ∈ \a is 
        //! and points \f$ x_{l} \f$ defined by \a pts. 
        //! Store results in \a matrix at index (j,l).
        //!
        virtual void     integrate_dy  ( const TIndexSet &               is,
                                         const std::vector< T3Point > &  pts,
                                         BLAS::Matrix< value_t > &       matrix ) const = 0;
    };

    DISABLE_COPY_OP( THCA );
};

//!
//! \ingroup Algebra_Module
//! \class   TPermHCAGeneratorFn
//! \brief   base class for HCA generator functions using row/column permutations
//!
//!          Provides basic permutation management for evaluating the integrals
//!          over the derivatives of the kernel generator function.
//!
template < typename T_value >
class TPermHCAGeneratorFn : public THCA< T_value >::TGeneratorFn
{
public:
    //
    // template types as internal types
    //
    using  value_t = T_value;

protected:
    
    // mapping from int. to ext. numbering
    const TPermutation *  _row_perm_i2e;
    const TPermutation *  _col_perm_i2e;
    
public:
    //!
    //! constructor
    //! - \a order defines (maximal) quadrature order for
    //!   evaluating the integrals 
    //!
    TPermHCAGeneratorFn ( const TPermutation *  row_perm_i2e,
                          const TPermutation *  col_perm_i2e )
            : _row_perm_i2e( row_perm_i2e ),
              _col_perm_i2e( col_perm_i2e )
    {}
    
    //!
    //! Evaluate \f$ \int \phi_i(x) D_x \gamma(x, y_{l}) dx\f$ for i ∈ \a is
    //! and points \f$ y_{l} \f$ defined by \a pts.
    //! Store results in \a matrix at index (i,l).
    //!
    virtual void  integrate_dx  ( const TIndexSet &               is,
                                  const std::vector< T3Point > &  pts,
                                  BLAS::Matrix< value_t > &       matrix ) const
    {
        //
        // permute indices
        //

        const bool            has_perm  = (this->_row_perm_i2e != NULL);
        std::vector< idx_t >  idxs( is.size() );
        size_t                pos = 0;

        for ( auto idx : is )
        {
            const idx_t  ex_idx = (has_perm ? this->_row_perm_i2e->permute( idx ) : idx );

            idxs[ pos++ ] = ex_idx;
        }// for

        integrate_dx_perm( idxs, pts, matrix );
    }
    
    //!
    //! Evaluate \f$ \int \phi_j(y) D_y \gamma(x_{l}, y) dy\f$ for j ∈ \a is 
    //! and points \f$ x_{l} \f$ defined by \a pts. 
    //! Store results in \a matrix at index (j,l).
    //!
    virtual void  integrate_dy  ( const TIndexSet &               is,
                                  const std::vector< T3Point > &  pts,
                                  BLAS::Matrix< value_t > &       matrix ) const
    {
        //
        // permute indices
        //

        const bool            has_perm  = (this->_col_perm_i2e != NULL);
        std::vector< idx_t >  idxs( is.size() );
        size_t                pos = 0;

        for ( auto idx : is )
        {
            const idx_t  ex_idx = (has_perm ? this->_col_perm_i2e->permute( idx ) : idx );

            idxs[ pos++ ] = ex_idx;
        }// for

        integrate_dy_perm( idxs, pts, matrix );
    }

    //!
    //! Evaluate \f$ \int \phi_i(x) D_x \gamma(x, y_{l}) dx\f$ for i ∈ \a is
    //! and points \f$ y_{l} \f$ defined by \a pts.
    //! Indices in \a idxs are given in external order.
    //!
    virtual void  integrate_dx_perm  ( const std::vector< idx_t > &    idxs,
                                       const std::vector< T3Point > &  pts,
                                       BLAS::Matrix< value_t > &       matrix ) const = 0;
    //!
    //! Evaluate \f$ \int \phi_j(y) D_y \gamma(x_{l}, y) dy\f$ for j ∈ \a is 
    //! and points \f$ x_{l} \f$ defined by \a pts. 
    //! Indices in \a idxs are given in external order.
    //!
    virtual void  integrate_dy_perm  ( const std::vector< idx_t > &    idxs,
                                       const std::vector< T3Point > &  pts,
                                       BLAS::Matrix< value_t > &       matrix ) const = 0;
};

//! \}

}// namespace Hpro

#endif  // __HPRO_TLOWRANKAPX_HH
