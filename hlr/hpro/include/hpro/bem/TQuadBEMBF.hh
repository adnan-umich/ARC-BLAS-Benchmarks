#ifndef __HPRO_TQUADBEMBF_HH
#define __HPRO_TQUADBEMBF_HH
//
// Project     : HLIBpro
// File        : TQuadBEMBF.hh
// Description : classes for bilinearforms in BEM-applications using quadrature
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/bem/TBEMBF.hh"

namespace Hpro
{
    
//!
//! \class  tripair_quad_rule_t
//! \brief  holds quadrature rule with points and weights for two triangles
//!
struct tripair_quad_rule_t
{
    // number of quadrature points
    size_t                  npts;
    
    // decoupled quadrature rule
    std::vector< T2Point >  pts1;    // coord. of points in triangle 1
    std::vector< T2Point >  pts2;    // coord. of points in triangle 1
    
    // "vectorised" coordinates (in single and double precision)
    std::vector< float >    x1_sgl, y1_sgl;  // coord. of points in triangle 1
    std::vector< float >    x2_sgl, y2_sgl;  // coord. of points in triangle 2

    std::vector< double >   x1_dbl, y1_dbl;
    std::vector< double >   x2_dbl, y2_dbl;

    // weight
    std::vector< float >    w_sgl;
    std::vector< double >   w_dbl;

    // resize to given number of points
    void
    resize ( const size_t  anpts )
    {
        npts = anpts;
        
        pts1.resize( npts );
        pts2.resize( npts );

        x1_sgl.resize( npts );
        y1_sgl.resize( npts );
        x2_sgl.resize( npts );
        y2_sgl.resize( npts );

        x1_dbl.resize( npts );
        y1_dbl.resize( npts );
        x2_dbl.resize( npts );
        y2_dbl.resize( npts );

        w_sgl.resize( npts );
        w_dbl.resize( npts );
    }

    // assign given data
    void
    assign ( const size_t     i,
             const T2Point &  pt1,
             const T2Point &  pt2,
             const double     wght )
    {
        pts1[i]   = pt1;
        pts2[i]   = pt2;

        x1_sgl[i] = pt1[0];
        y1_sgl[i] = pt1[1];
        x2_sgl[i] = pt2[0];
        y2_sgl[i] = pt2[1];

        x1_dbl[i] = pt1[0];
        y1_dbl[i] = pt1[1];
        x2_dbl[i] = pt2[0];
        y2_dbl[i] = pt2[1];

        w_sgl[i]  = wght;
        w_dbl[i]  = wght;
    }
    
    // return data depending on input type
    template < typename  real_t >
    const real_t *  weights () const
    {
        if constexpr ( std::is_same_v< float,  real_t > ) return w_sgl.data();
        if constexpr ( std::is_same_v< double, real_t > ) return w_dbl.data();
    }

    template < typename  real_t >
    const real_t *  x1 () const
    {
        if constexpr ( std::is_same_v< float,  real_t > ) return x1_sgl.data();
        if constexpr ( std::is_same_v< double, real_t > ) return x1_dbl.data();
    }

    template < typename  real_t >
    const real_t *  x2 () const
    {
        if constexpr ( std::is_same_v< float,  real_t > ) return x2_sgl.data();
        if constexpr ( std::is_same_v< double, real_t > ) return x2_dbl.data();
    }

    template < typename  real_t >
    const real_t *  y1 () const
    {
        if constexpr ( std::is_same_v< float,  real_t > ) return y1_sgl.data();
        if constexpr ( std::is_same_v< double, real_t > ) return y1_dbl.data();
    }

    template < typename  real_t >
    const real_t *  y2 () const
    {
        if constexpr ( std::is_same_v< float,  real_t > ) return y2_sgl.data();
        if constexpr ( std::is_same_v< double, real_t > ) return y2_dbl.data();
    }
};

//!
//! \ingroup BEM_Module
//! \class   TQuadBEMBF
//! \brief   Base class for all quadrature based bilinear forms.
//!
//!          TQuadBEMBF extends TBEMBF by providing quadrature rules for
//!          triangle pairs and defining a kernel function evaluation interface.
//!
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_value >
class TQuadBEMBF : public TBEMBF< T_ansatzsp, T_testsp, T_value >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t     = T_ansatzsp;
    using  testsp_t       = T_testsp;
    using  value_t        = T_value;
    using  real_t         = real_type_t< value_t >;
    using  ansatz_value_t = value_type_t< ansatzsp_t >;
    using  test_value_t   = value_type_t< testsp_t >;

protected:
    
    // type for storing quadrature points for different order and type
    using  quad_rules_t = std::vector< std::vector< tripair_quad_rule_t > >;

protected:
    //! @cond
    
    // quadrature order
    const uint      _quad_order;
    
    // quadrature rules for all orders and types
    quad_rules_t    _quad_rules;

    // if true, the quad. order is reduced depending on triangle distance
    bool            _quad_dist_adaptive;
    
    // relative error for adaptive order
    real_t          _ada_error;
    
    //! @endcond
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct bilinear form over function spaces \a ansatzsp and \a testsp
    //! using quadratur with order \a order; if \a dist_ada is true, the
    //! quadrature order is adaptively adjusted according to distance
    TQuadBEMBF ( const ansatzsp_t *  ansatzsp,
                 const testsp_t *    testsp,
                 const uint          order    = CFG::BEM::quad_order,
                 const bool          dist_ada = CFG::BEM::adaptive_quad_order,
                 const real_t        ada_err  = 0 );

    //! destructor
    virtual ~TQuadBEMBF () {}

    //! set accuracy for quadrature evaluation
    virtual void set_accuracy ( const double  acc )
    {
        if ( acc > 0 )
        {
            _quad_dist_adaptive = true;
            _ada_error          = acc;
        }// if
    }
    
    //////////////////////////////////////
    //
    // evaluate bilinearform
    //

    //! evaluate subblock defined by \a row_ind × \a col_ind; the indices
    //! in \a row_ind and \a col_ind can be arbitrary, e.g. must not be
    //! contiguous
    virtual void  eval  ( const std::vector< idx_t > &  row_ind,
                          const std::vector< idx_t > &  col_ind,
                          BLAS::Matrix< value_t > &     values ) const;

    //! same as above but accuracy during coefficient evaluation is set to <eps>
    virtual void  eval  ( const std::vector< idx_t > &  row_ind,
                          const std::vector< idx_t > &  col_ind,
                          BLAS::Matrix< value_t > &     values,
                          const double                  eps ) const;
    
    //
    // functions for determining quadrature rule
    //
    
    //! reorder triangle vertices such that the k common vertices are ordered
    //! from 0,…,k-1; the number of common vertices is returned
    uint  reorder_common  ( idx_t *        vtx0idxs,
                            idx_t *        vtx1idxs ) const;
    
    //! adjust quadrature order \a order depending on diameter and distance of triangles
    uint  adjust_order    ( const idx_t *  vtx0idxs,
                            const idx_t *  vtx1idxs,
                            const uint     order ) const;
    
    //! return quadrature rule for \a ncommon vertices and order \a order
    const tripair_quad_rule_t *  quad_rule ( const uint  ncommon,
                                             const uint  order ) const
    {
        return & _quad_rules[ ncommon ][ order ];
    }
    
    //
    // kernel function evaluation
    //

    //! compute kernel at quadrature points in triangles \a tri0idx and \a tri1idx 
    //! with coordinate indices \a tri0 and \a tri1 using quadrature rule \a quad_rule;
    //! the results for all points are returned in \a values
    virtual void  eval_kernel  ( const idx_t                  tri0idx,
                                 const idx_t                  tri1idx,
                                 const TGrid::triangle_t &    tri0,
                                 const TGrid::triangle_t &    tri1,
                                 const tripair_quad_rule_t *  quad_rule,
                                 std::vector< value_t > &     values,
                                 const real_t                 eps = 0 ) const = 0;

    DISABLE_COPY_OP( TQuadBEMBF );
};

//!
//! \ingroup BEM_Module
//! \class   TInvarBasisQuadBEMBF
//! \brief   Class for quadrature based bilinear forms with invariant basis functions.
//!
//!          TConstBasisQuadBEMBF extends TQuadBEMBF for function spaces with invariant
//!          basis functions, e.g. with respect to translation, scaling and rotation.
//!          This allows their precomputation for the unit triangle, improving efficiency.
//!
template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_val >
class TInvarBasisQuadBEMBF : public TQuadBEMBF< T_ansatzsp, T_testsp, T_val >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t     = T_ansatzsp;
    using  testsp_t       = T_testsp;
    using  value_t        = T_val;
    using  real_t         = real_type_t< value_t >;
    using  ansatz_value_t = value_type_t< ansatzsp_t >;
    using  test_value_t   = value_type_t< testsp_t >;

protected:
    
    // type for storing basis function values
    using  ansatz_store_t = std::vector< std::vector< std::vector< std::vector< ansatz_value_t > > > >;
    using  test_store_t   = std::vector< std::vector< std::vector< std::vector< test_value_t > > > >;

protected:
    //! @cond
    
    // ansatz function values for different quadrature rules
    ansatz_store_t  _ansatz_val;

    // test function values for different quadrature rules
    test_store_t    _test_val;

    //! @endcond
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct bilinear form over function spaces \a ansatzsp and \a testsp
    //! using quadratur with order \a order; if \a dist_ada is true, the
    //! quadrature order is adaptively adjusted according to distance
    TInvarBasisQuadBEMBF ( const ansatzsp_t *  ansatzsp,
                           const testsp_t *    testsp,
                           const uint          order    = CFG::BEM::quad_order,
                           const bool          dist_ada = CFG::BEM::adaptive_quad_order,
                           const real_t        ada_err  = 0 );

    //! destructor
    virtual ~TInvarBasisQuadBEMBF () {}

    //////////////////////////////////////
    //
    // evaluate bilinearform
    //

    //! evaluate subblock defined by \a row_ind × \a col_ind; the indices
    //! in \a row_ind and \a col_ind can be arbitrary, e.g. must not be
    //! contiguous
    virtual void  eval  ( const std::vector< idx_t > &  row_ind,
                          const std::vector< idx_t > &  col_ind,
                          BLAS::Matrix< value_t > &     values ) const;

    //! same as above but accuracy during coefficient evaluation is set to <eps>
    virtual void  eval  ( const std::vector< idx_t > &  row_ind,
                          const std::vector< idx_t > &  col_ind,
                          BLAS::Matrix< value_t > &     values,
                          const double                  eps ) const;
    
    //
    // basis function evaluation
    //

    //! compute ansatz and test basis functions for all quadrature points
    void  compute_basis_func ();
    
    //! return ansatz basis function values for index \a idx in triangle \a tri
    //! for quadrature points with \a ncommon common vertices and order \a order
    const std::vector< ansatz_value_t > *
    ansatz_val ( const idx_t                idx,
                 const TGrid::triangle_t &  tri,
                 const uint                 ncommon,
                 const uint                 order ) const
    {
        return & _ansatz_val[ this->ansatz_space()->triangle_index( idx, tri ) ][ ncommon ][ order ];
    }
    
    //! same as \see ansatz_val but for test space
    const std::vector< test_value_t > *
    test_val   ( const idx_t                idx,
                 const TGrid::triangle_t &  tri,
                 const uint                 ncommon,
                 const uint                 order ) const
    {
        return & _test_val[ this->test_space()->triangle_index( idx, tri ) ][ ncommon ][ order ];
    }

    DISABLE_COPY_OP( TInvarBasisQuadBEMBF );
};

}// namespace Hpro

#endif  // __HPRO_TBEMBF_HH
