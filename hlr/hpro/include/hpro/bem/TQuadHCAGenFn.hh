#ifndef __HPRO_TQUADHCAGENFN_HH
#define __HPRO_TQUADHCAGENFN_HH
//
// Project     : HLIBpro
// File        : TQuadHCAGenFn.hh
// Description : class providing HCA functionality for BEM bilinear forms
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/algebra/TLowRankApx.hh"
#include "hpro/bem/TFnSpace.hh"

namespace Hpro
{

//!
//! quadrature rule for single triangle in "vectorised" form
//! - size of x, y, w ≥ npts due to possible padding
//!
template < typename value_t >
struct tri_quad_rule_t
{
    // number of quadrature points
    size_t                   npts;
    
    // decoupled x and y coordinates of quadrature points
    std::vector< value_t >   x, y;
    
    // quadrature weights
    std::vector< value_t >   w;
};

//!
//! \ingroup BEM_Module
//! \class   TQuadHCAGenFn
//! \brief   base class for HCA generator functions using quadrature
//!
//!          Evaluates the integrals over the derivatives of the kernel generator
//!          function by means of quadrature.
//!
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
class TQuadHCAGenFn : public TPermHCAGeneratorFn< T_value >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = T_value;
    using  real_t     = real_type_t< value_t >;

    using  quad_rule_vec_t = std::vector< tri_quad_rule_t< real_t > >;
    
    //!
    //! statistics
    //!
    struct stat_t
    {
        //! number of calls to "eval_dx"
        size_t  n_eval_dx;

        //! number of calls to "eval_dy"
        size_t  n_eval_dy;

        //! default ctor
        stat_t ()
                : n_eval_dx(0), n_eval_dy(0)
        {}
    };

protected:
    // ansatz space
    const ansatzsp_t *  _ansatz_sp;
    
    // test space
    const testsp_t *    _test_sp;
    
    // quadrature order to use
    const uint          _quad_order;
    
    // cache for quadratur points and weights
    quad_rule_vec_t     _quad_rules_cache;

    // statistics record
    stat_t *            _stat;
    
public:
    //!
    //! constructor for HCA generator function over ansatz space \a ansatzsp
    //! test space \a testsp with permutations for row and column indices
    //! and (maximal) quadrature order \a quad_order
    //!
    TQuadHCAGenFn ( const ansatzsp_t *    ansatzsp,
                    const testsp_t *      testsp,
                    const uint            quad_order,
                    const TPermutation *  row_perm_i2e,
                    const TPermutation *  col_perm_i2e,
                    stat_t *              stat = NULL );
    
    //
    // access function spaces
    //

    const ansatzsp_t *  ansatz_space () const { return _ansatz_sp; }
    const testsp_t *    test_space   () const { return _test_sp; }


    //!
    //! Evaluate \f$ \int \phi_i(x) D_x \gamma(x, y_{l}) dx\f$ for i ∈ \a is
    //! and points \f$ y_{l} \f$ defined by \a pts.
    //!
    virtual void  integrate_dx_perm  ( const std::vector< idx_t > &    idxs,
                                       const std::vector< T3Point > &  pts,
                                       BLAS::Matrix< value_t > &       matrix ) const;
    //!
    //! Evaluate \f$ \int \phi_j(y) D_y \gamma(x_{l}, y) dy\f$ for j ∈ \a is 
    //! and points \f$ x_{l} \f$ defined by \a pts. 
    //!
    virtual void  integrate_dy_perm  ( const std::vector< idx_t > &    idxs,
                                       const std::vector< T3Point > &  pts,
                                       BLAS::Matrix< value_t > &       matrix ) const;

    //
    // handle statistics
    //

    void  set_stat ( stat_t *  stat ) { _stat = stat; }

protected:
    //
    // quadrature management
    //
    
    // return quadrature points for preset order
    const tri_quad_rule_t< real_t > *  get_quad_rule  ()  const
    {
        return  & _quad_rules_cache[_quad_order];
    }

    // return quadrature points for order \a order
    // - \a order must not be larger than order in ctor
    //
    const tri_quad_rule_t< real_t > *  get_quad_rule  ( const uint  order ) const
    {
        if ( order >= _quad_rules_cache.size() )
            HERROR( ERR_ARG, "(TQuadHCAGenFn) get_quad_rule", "order exceeds maximal order in ctor" );
    
        return  & _quad_rules_cache[ order ];
    }

    //
    // generator derivative evaluation for quadrature
    //

    //! Evaluate \f$ D_x \gamma(x, y) \f$ on with \f$x\f$ defined by quadrature
    //! points on triangle \a tri_idx. The computed values for each quadrature
    //! point i are stored on \a values[i].
    virtual void  eval_dx  ( const idx_t                        tri_idx,
                             const T3Point &                    y,
                             const tri_quad_rule_t< real_t > &  quad_rule,
                             std::vector< value_t > &           values ) const = 0;
    
    //! Evaluate \f$ D_y \gamma(x, y) \f$ on with \f$y\f$ defined by quadrature
    //! points on triangle \a tri_idx. The computed values for each quadrature
    //! point i are stored on \a values[i].
    virtual void  eval_dy  ( const T3Point &                    x,
                             const idx_t                        tri_idx,
                             const tri_quad_rule_t< real_t > &  quad_rule,
                             std::vector< value_t > &           values ) const = 0;
    
};

//!
//! \ingroup BEM_Module
//! \class   TInvarBasisQuadHCAGenFn
//! \brief   class for BEM HCA generator functions with invariant basis functions
//!
//!          Optimises evaluation of basis function for quadrature points under the
//!          assumption, that they are invariant with respect to translation and
//!          scaling, e.g. can be precomputed on reference triangles.
//!
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
class TInvarBasisQuadHCAGenFn : public TQuadHCAGenFn< T_ansatzsp, T_testsp, T_value >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = T_value;
    using  real_t     = real_type_t< value_t >;

    //
    // inherit from base class
    //
    using  stat_t     = typename TQuadHCAGenFn< ansatzsp_t, testsp_t, value_t >::stat_t;
    
    //
    // value types of basis functions
    //
    using  ansatz_value_t = value_type_t< ansatzsp_t >;
    using  test_value_t   = value_type_t< testsp_t >;
    
protected:
    
    // type for storing basis function values
    using  ansatz_store_t = std::vector< std::vector< std::vector< ansatz_value_t > > >;
    using  test_store_t   = std::vector< std::vector< std::vector< test_value_t > > >;

protected:
    // cache for precomputed basis functions
    ansatz_store_t  _ansatz_val;
    test_store_t    _test_val;
    
public:
    //!
    //! constructor for HCA generator function over ansatz space \a ansatzsp
    //! test space \a testsp with permutations for row and column indices
    //! and (maximal) quadrature order \a quad_order
    //!
    TInvarBasisQuadHCAGenFn ( const ansatzsp_t *    ansatzsp,
                              const testsp_t *      testsp,
                              const uint            quad_order,
                              const TPermutation *  row_perm_i2e,
                              const TPermutation *  col_perm_i2e,
                              stat_t *              stat = NULL );
    

    //!
    //! Evaluate \f$ \int \phi_i(x) D_x \gamma(x, y_{l}) dx\f$ for i ∈ \a is
    //! and points \f$ y_{l} \f$ defined by \a pts.
    //! Store results in \a matrix at index (i,l).
    //!
    virtual void  integrate_dx_perm  ( const std::vector< idx_t > &    idxs,
                                       const std::vector< T3Point > &  pts,
                                       BLAS::Matrix< value_t > &       matrix ) const;
    //!
    //! Evaluate \f$ \int \phi_j(y) D_y \gamma(x_{l}, y) dy\f$ for j ∈ \a is 
    //! and points \f$ x_{l} \f$ defined by \a pts. 
    //! Store results in \a matrix at index (j,l).
    //!
    virtual void  integrate_dy_perm  ( const std::vector< idx_t > &    idxs,
                                       const std::vector< T3Point > &  pts,
                                       BLAS::Matrix< value_t > &       matrix ) const;

protected:
    //
    // management of basis function values
    //

    //! return ansatz basis function values for index \a idx in triangle \a tri
    //! for quadrature points of order \a order
    const std::vector< ansatz_value_t > *
    ansatz_val ( const idx_t                idx,
                 const TGrid::triangle_t &  tri,
                 const uint                 order ) const
    {
        return & _ansatz_val[ this->ansatz_space()->triangle_index( idx, tri ) ][ order ];
    }
    
    //! same as \see ansatz_val but for test space
    const std::vector< test_value_t > *
    test_val   ( const idx_t                idx,
                 const TGrid::triangle_t &  tri,
                 const uint                 order ) const
    {
        return & _test_val[ this->test_space()->triangle_index( idx, tri ) ][ order ];
    }
};

}// namespace Hpro

#endif  // __HPRO_TQUADHCAGENFN_HH
