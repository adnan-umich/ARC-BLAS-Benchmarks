#ifndef __HPRO_TMATERNCOVCOEFFFN_HH
#define __HPRO_TMATERNCOVCOEFFFN_HH
//
// Project     : HLIBpro
// File        : TMaternCovCoeffFn.hh
// Description : matrix coefficients for matern covariance kernel
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/TPoint.hh"
#include "hpro/matrix/TCoeffFn.hh"

namespace Hpro
{

//!
//! \class   TMaternCovCoeffFn
//! \brief   Matern covariance coefficient function
//!
//!          TMaternCovCoeffFn implements the Matern covarience kernel
//!          \f[ C(h, \theta ) := \frac{\sigma^2}{2^{\nu-1} \Gamma(\nu)} \left(\frac{h}{\ell}\right)^{\nu} K_{\nu}\left(\frac{h}{\ell}\right) \f]
//!          with \f$\theta=(\sigma,\ell,\nu)\f$, where \f$\sigma\f$ is the variance; \f$\nu > 0\f$ controls
//!          the smoothness of the random field, with larger values of corresponding to smoother fields;
//!          and \f$\ell\f$ is a spatial range parameter that measures how quickly the correlation of the
//!          random field decays with distance.
//!
template < typename T_point = TPoint >
class TMaternCovCoeffFn : public TCoeffFn< double >
{
public:
    using  value_t = double;
    using  point_t = T_point;
    
private:
    //! @cond
    
    // special type of Matern kernel for optimization
    enum  matern_type_t  { one_half, three_half, five_half, general };
    
    const value_t                   _length;
    const value_t                   _nu;
    const value_t                   _sigmasq;
    const value_t                   _scale_fac;
    const std::vector< point_t > &  _row_vertices;
    const std::vector< point_t > &  _col_vertices;
    const matern_type_t             _matern_type;

    //! @endcond
    
public:
    //!
    //! constructor
    //!
    TMaternCovCoeffFn ( const value_t                   sigma,
                        const value_t                   length,
                        const value_t                   nu,
                        const std::vector< point_t > &  vertices );

    TMaternCovCoeffFn ( const value_t                   sigma,
                        const value_t                   length,
                        const value_t                   nu,
                        const std::vector< point_t > &  row_vertices,
                        const std::vector< point_t > &  col_vertices );

    //!
    //! coefficient evaluation
    //!
    virtual void eval  ( const std::vector< idx_t > &  rowidxs,
                         const std::vector< idx_t > &  colidxs,
                         value_t *                     matrix ) const;
    using TCoeffFn< value_t >::eval;

    //!
    //! return format of matrix, e.g. symmetric or hermitian
    //!
    virtual matform_t  matrix_format  () const { return symmetric; }
};

}// namespace Hpro

#endif //  __HPRO_TMATERNCOVCOEFFFN_HH
