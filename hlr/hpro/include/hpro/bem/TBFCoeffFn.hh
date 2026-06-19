#ifndef __HPRO_BEM_TBFCOEFFFN_H
#define __HPRO_BEM_TBFCOEFFFN_H
//
// Project     : HLIBpro
// File        : TBFCoeffFn.hh
// Description : matrix coefficient functions for bilinear forms
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/matrix/TCoeffFn.hh"
#include "hpro/bem/TBEMBF.hh"

namespace Hpro
{

////////////////////////////////////////////////////////////////////
//!
//! \class  TBFCoeffFn
//! \brief  Provide matrix coefficients defined by bilinear forms
//!
template < typename  T_bf > 
class TBFCoeffFn : public TCoeffFn< typename T_bf::value_t >
{
public:
    //
    // template types as internal types
    //
    using  bf_t       = T_bf;
    using  ansatzsp_t = typename bf_t::ansatzsp_t;
    using  testsp_t   = typename bf_t::testsp_t;
    using  value_t    = typename bf_t::value_t;

protected:
    //! bilinear form to evaluate
    bf_t *  _bf;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct coefficient function for bilinear form \a bf
    TBFCoeffFn ( bf_t *  bf )
            : _bf( bf )
    {
        if ( _bf == NULL )
            HERROR( ERR_ARG, "(TBFCoeffFn)", "invalid bilinear form supplied (NULL)" );
    }

    //! destructor
    virtual ~TBFCoeffFn () {}

    //! return true if function is complex valued
    bool       is_complex     () const { return _bf->is_complex(); }
        
    //! return format of matrix, e.g. symmetric or hermitian
    matform_t  matrix_format  () const { return _bf->format(); }
    
    //! set accuracy for coefficient evaluation
    virtual void set_accuracy ( const double  acc )
    {
        _bf->set_accuracy( acc );
    }
    
    //////////////////////////////////////
    //
    // evaluation of coefficient
    //

    //! evaluate matrix coefficients in \a rowis × \a colis and store values
    //! in \a matrix (real valued)
    virtual void eval ( const std::vector< idx_t > &  rowidxs,
                        const std::vector< idx_t > &  colidxs,
                        value_t *                     matrix,
                        const double                  eps = 0 ) const
    {
        //
        // call bilinear form
        //
        
        const size_t  n = rowidxs.size();
        const size_t  m = colidxs.size();
        auto          M = BLAS::Matrix< value_t >( matrix, n, 1, m, n );
        
        _bf->eval( rowidxs, colidxs, M, eps );
            
        for ( idx_t  i = 0; i < idx_t(n*m); ++i )
            matrix[i] = M.data()[i];
    }

    using TCoeffFn< value_t >::eval;

    DISABLE_COPY_OP( TBFCoeffFn );
};

}// namespace Hpro

#endif  // __HPRO_BEM_TBFCOEFFFN_H
