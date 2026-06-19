#ifndef __HPRO_THELMHOLTZBF_HH
#define __HPRO_THELMHOLTZBF_HH
//
// Project     : HLIBpro
// File        : THelmholtzBF.hh
// Description : bilinear forms for Helmholtz operator
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/algebra/TLowRankApx.hh"
#include "hpro/bem/TQuadBEMBF.hh"
#include "hpro/bem/TQuadHCAGenFn.hh"

namespace Hpro
{
    
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// standard bilinear forms for Helmholtz
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
//
// single layer potential of Helmholtz operator
//
////////////////////////////////////////////////////

//!
//! \ingroup BEM_Module
//! \class   THelmholtzSLPBF
//! \brief   Bilinear form for Helmholtz single layer potential
//!
//!          THelmholtzSLPBF implements the bilinear form for the Helmholtz
//!          single layer potential with the kernel function
//!          \f[ \frac{e^{i\kappa \|x-y\|_2}}{\|x-y\|_2} \f]
//!          i.e. for the integral equation
//!          \f[ 4 \pi \int_{\Gamma} \frac{u(y) \cdot e^{i\kappa \|x-y\|_2}}{\|x-y\|_2} dy = f(x) \f]
//!
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
class THelmholtzSLPBF : public TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, T_value >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = T_value;
    using  real_t     = real_type_t< value_t >;

private:
    //! @cond
    
    // i * wave number
    const value_t  _ikappa;

    // kernel function implementation
    void ( * _kernel_fn ) ( const TGrid::triangle_t &    tri0,
                            const TGrid::triangle_t &    tri1,
                            const tripair_quad_rule_t *  rule,
                            const value_t                ikappa,
                            const ansatzsp_t *           ansatz_sp,
                            const testsp_t *             test_sp,
                            std::vector< value_t > &     values );
    
    //! @endcond
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //
    
    THelmholtzSLPBF ( const value_t       kappa,
                      const ansatzsp_t *  aansatzsp,
                      const testsp_t *    atestsp,
                      const uint          quad_order = CFG::BEM::quad_order );

    THelmholtzSLPBF ( const value_t       kappa,
                      const ansatzsp_t *  aansatzsp,
                      const testsp_t *    atestsp,
                      const real_t        quad_error );

    virtual ~THelmholtzSLPBF () {}
    
    // return format of bilinear form, e.g. symmetric
    matform_t  format () const;

    // eval kernel function at quadrature points
    virtual void  eval_kernel  ( const idx_t                  tri0idx,
                                 const idx_t                  tri1idx,
                                 const TGrid::triangle_t &    tri0,
                                 const TGrid::triangle_t &    tri1,
                                 const tripair_quad_rule_t *  quad_rule,
                                 std::vector< value_t > &     values,
                                 const real_t                 eps = 0 ) const;
};

////////////////////////////////////////////////////
//
// double layer potential of Helmholtz operator
//
////////////////////////////////////////////////////

//!
//! \ingroup BEM_Module
//! \class   THelmholtzDLPBF
//! \brief   Bilinear form for Helmholtz double layer potential
//! 
//!          THelmholtzDLPBF implements the bilinear form for the Helmholtz
//!          double layer potential with the kernel function
//!          \f[ \frac{e^{i \cdot \kappa \|x-y\|_2} (i\cdot \kappa \|x-y\|_2 - 1) \langle n(y), y-x \rangle}{\|x-y\|_2^3} \f].
//!          or adjoint form
//!          \f[ \frac{e^{i \cdot \kappa \|x-y\|_2} (i\cdot \kappa \|x-y\|_2 - 1) \langle n(x), x-y \rangle}{\|x-y\|_2^3} \f].
//!
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
class THelmholtzDLPBF : public TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, T_value >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = T_value;
    using  real_t     = real_type_t< value_t >;

private:
    //! @cond

    // enables/disables adjoint form
    const bool     _adjoint;
    
    // i times wave number
    const value_t  _ikappa;
    
    // kernel function implementation
    void ( * _kernel_fn ) ( const idx_t                  tri_id,
                            const bool                   adjoint,
                            const TGrid::triangle_t &    tri0,
                            const TGrid::triangle_t &    tri1,
                            const tripair_quad_rule_t *  rule,
                            const value_t                ikappa,
                            const ansatzsp_t *           ansatz_sp,
                            const testsp_t *             test_sp,
                            std::vector< value_t > &     values );
    
    //! @endcond
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //
    
    THelmholtzDLPBF ( const value_t       kappa,
                      const ansatzsp_t *  aansatzsp,
                      const testsp_t *    atestsp,
                      const bool          adjoint    = false,
                      const uint          quad_order = CFG::BEM::quad_order );

    THelmholtzDLPBF ( const value_t       kappa,
                      const ansatzsp_t *  aansatzsp,
                      const testsp_t *    atestsp,
                      const bool          adjoint,
                      const real_t        quad_error );

    virtual ~THelmholtzDLPBF () {}

    // return format of bilinear form, e.g. symmetric
    matform_t  format () const;

    // eval kernel function at quadrature points
    virtual void  eval_kernel  ( const idx_t                  tri0idx,
                                 const idx_t                  tri1idx,
                                 const TGrid::triangle_t &    tri0,
                                 const TGrid::triangle_t &    tri1,
                                 const tripair_quad_rule_t *  quad_rule,
                                 std::vector< value_t > &     values,
                                 const real_t                 eps = 0 ) const;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// HCA functions for Helmholtz
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

//!
//! \ingroup BEM_Module
//! \class   THelmholtzSLPGenFn
//! \brief   kernel generator function for Helmholtz SLP
//!
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
class THelmholtzSLPGenFn : public TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, T_value >
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
    
private:

    //
    // HCA function impl.
    //

    void ( * _eval_dxy_impl ) ( const value_t &                   ikappa,
                                const tri_quad_rule_t< real_t > & quad_rule,
                                const T3Point                     x[3],
                                const T3Point &                   y,
                                std::vector< value_t > &          values );
    
protected:

    // i · wave number
    const value_t  _ikappa;
    
public:
    //
    // constructor
    //
    THelmholtzSLPGenFn ( const value_t         kappa,
                         const ansatzsp_t *    ansatzsp,
                         const testsp_t *      testsp,
                         const TPermutation *  row_perm_i2e,
                         const TPermutation *  col_perm_i2e,
                         const uint            quad_order = CFG::BEM::quad_order );
    
    //!
    //! evaluate generator function at (\a x, \a y)
    //!
    value_t  eval  ( const T3Point &  x,
                     const T3Point &  y ) const
    {
        const auto  one_over_4pi = real_t(1) / ( real_t(4) * Math::pi< real_t >() );
        const auto  dist         = real_t( ( x - y ).norm2() );
        
        return one_over_4pi * Math::exp( _ikappa * dist ) / ( dist );
    }

protected:
    //! Evaluate \f$ D_x \gamma(x, y) \f$ on with \f$x\f$ defined by quadrature
    //! points on triangle \a tri_idx. The computed values for each quadrature
    //! point i are stored on \a values[i].
    virtual void  eval_dx  ( const idx_t                       tri_idx,
                             const T3Point &                   y,
                             const tri_quad_rule_t< real_t > & quad_rule,
                             std::vector< value_t > &          values ) const
    {
        const TGrid *            grid = this->ansatz_space()->grid();
        const TGrid::triangle_t  tri  = grid->triangle( tri_idx );
        const T3Point            x[3] = { grid->vertex( tri.vtx[0] ),
                                          grid->vertex( tri.vtx[1] ),
                                          grid->vertex( tri.vtx[2] ) };
    
        _eval_dxy_impl( _ikappa, quad_rule, x, y, values );
    }
    
    //! Evaluate \f$ D_y \gamma(x, y) \f$ on with \f$y\f$ defined by quadrature
    //! points on triangle \a tri_idx. The computed values for each quadrature
    //! point i are stored on \a values[i].
    virtual void  eval_dy  ( const T3Point &                   x,
                             const idx_t                       tri_idx,
                             const tri_quad_rule_t< real_t > & quad_rule,
                             std::vector< value_t > &          values ) const
    {
        const TGrid *            grid = this->test_space()->grid();
        const TGrid::triangle_t  tri  = grid->triangle( tri_idx );
        const T3Point            y[3] = { grid->vertex( tri.vtx[0] ),
                                          grid->vertex( tri.vtx[1] ),
                                          grid->vertex( tri.vtx[2] ) };
    
        _eval_dxy_impl( _ikappa, quad_rule, y, x, values );
    }

    DISABLE_COPY_OP( THelmholtzSLPGenFn );
};

//!
//! \ingroup BEM_Module
//! \class   THelmholtzDLPGenFn
//! \brief   kernel generator function for Helmholtz DLP
//!
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
class THelmholtzDLPGenFn : public TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, T_value >
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
    
private:

    //
    // HCA function impl.
    //

    void ( * _eval_dx_impl ) ( const value_t &                   ikappa,
                               const tri_quad_rule_t< real_t > & quad_rule,
                               const T3Point                     x[3],
                               const T3Point &                   y,
                               std::vector< value_t > &          values );
    
    void ( * _eval_dy_impl ) ( const value_t &                   ikappa,
                               const tri_quad_rule_t< real_t > & quad_rule,
                               const T3Point &                   x,
                               const T3Point                     y[3],
                               const T3Point &                   normal,
                               std::vector< value_t > &          values );
    
protected:

    // i · wave number
    const value_t  _ikappa;
    
public:
    //
    // constructor
    //
    THelmholtzDLPGenFn ( const value_t         kappa,
                         const ansatzsp_t *    ansatzsp,
                         const testsp_t *      testsp,
                         const TPermutation *  row_perm_i2e,
                         const TPermutation *  col_perm_i2e,
                         const uint            quad_order = CFG::BEM::quad_order );
    
    //!
    //! evaluate generator function at (\a x, \a y)
    //!
    value_t  eval  ( const T3Point &  x,
                     const T3Point &  y ) const
    {
        const auto  one_over_4pi = real_t(1) / ( real_t(4) * Math::pi< real_t >() );
        const auto  dist         = real_t( ( x - y ).norm2() );
        
        return one_over_4pi * Math::exp( _ikappa * dist ) / ( dist );
    }

protected:
    //! Evaluate \f$ D_x \gamma(x, y) \f$ on with \f$x\f$ defined by quadrature
    //! points on triangle \a tri_idx. The computed values for each quadrature
    //! point i are stored on \a values[i].
    virtual void  eval_dx  ( const idx_t                       tri_idx,
                             const T3Point &                   y,
                             const tri_quad_rule_t< real_t > & quad_rule,
                             std::vector< value_t > &          values ) const
    {
        const TGrid *            grid = this->ansatz_space()->grid();
        const TGrid::triangle_t  tri  = grid->triangle( tri_idx );
        const T3Point            x[3] = { grid->vertex( tri.vtx[0] ),
                                          grid->vertex( tri.vtx[1] ),
                                          grid->vertex( tri.vtx[2] ) };
    
        _eval_dx_impl( _ikappa, quad_rule, x, y, values );
    }
    
    //! Evaluate \f$ D_y \gamma(x, y) \f$ on with \f$y\f$ defined by quadrature
    //! points on triangle \a tri_idx. The computed values for each quadrature
    //! point i are stored on \a values[i].
    virtual void  eval_dy  ( const T3Point &                   x,
                             const idx_t                       tri_idx,
                             const tri_quad_rule_t< real_t > & quad_rule,
                             std::vector< value_t > &          values ) const
    {
        const TGrid *            grid = this->ansatz_space()->grid();
        const TGrid::triangle_t  tri  = grid->triangle( tri_idx );
        const T3Point            y[3] = { grid->vertex( tri.vtx[0] ),
                                          grid->vertex( tri.vtx[1] ),
                                          grid->vertex( tri.vtx[2] ) };
        const T3Point            normal( grid->tri_normal( tri_idx ) );
    
        _eval_dy_impl( _ikappa, quad_rule, x, y, normal, values );
    }

    DISABLE_COPY_OP( THelmholtzDLPGenFn );
};

}// namespace

#endif  // __THELMHOLTZBF_HH
