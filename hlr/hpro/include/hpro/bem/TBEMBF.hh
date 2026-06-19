#ifndef __HPRO_TBEMBF_HH
#define __HPRO_TBEMBF_HH
//
// Project     : HLIBpro
// File        : TBEMBF.hh
// Description : classes for bilinearforms in BEM-applications
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/blas/Matrix.hh"
#include "hpro/bem/TFnSpace.hh"
#include "hpro/matrix/TMatrix.hh"

namespace Hpro
{
    
////////////////////////////////////////////////////////////////////////////
//
// TBilinearForm
//
////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup BEM_Module
//! \class   TBilinearForm
//! \brief   Base class for all bilinear forms.
//!
//!          TBilinearForm is mainly introduced to use BEM bilinear forms
//!          without ansatz/test space arguments in evaluation, e.g. if
//!          only a generic pointer to a bilinear form is needed without
//!          any further knowledge.
//!
template < typename  T_val >
class TBilinearForm
{
public:
    //
    // template types as internal types
    //

    using  value_t = T_val;
        
public:

    //
    // dtor
    //
    virtual ~TBilinearForm () {}
    
    //////////////////////////////////////
    //
    // evaluate bilinearform
    //

    //! evaluate subblock defined by \a row_ind × \a col_ind; the indices
    //! in \a row_ind and \a col_ind can be arbitrary, e.g. must not be
    //! contiguous
    virtual void       eval       ( const std::vector< idx_t > &  row_ind,
                                    const std::vector< idx_t > &  col_ind,
                                    BLAS::Matrix< value_t > &     values ) const = 0;

    //! same as above but accuracy during coefficient evaluation is set to <eps>
    virtual void       eval       ( const std::vector< idx_t > &  row_ind,
                                    const std::vector< idx_t > &  col_ind,
                                    BLAS::Matrix< value_t > &     values,
                                    const double                  eps ) const = 0;

    //! return true if bilinear form is complex valued
    bool               is_complex () const { return is_complex_type< value_t >::value; }

    //! return format of bilinear form, e.g. symmetric
    virtual matform_t  format     () const { return MATFORM_NONSYM; }
    
    //! set accuracy for evaluation
    virtual void set_accuracy ( const double  /* acc */ ) {} // ignored by default
    
protected:
    DISABLE_COPY_OP( TBilinearForm );
};



////////////////////////////////////////////////////////////////////////////
//
// TBEMBF
//
////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup BEM_Module
//! \class   TBEMBF
//! \brief   Base class for BEM bilinear forms with ansatz and test space.
//!
//!          TBEMBF is the \em real base class for all BEM bilinear forms.
//!          It defines ansatz and test spaces in the form of template
//!          parameters.
//!
template < typename  T_ansatzsp,   // ansatz space
           typename  T_testsp,     // test space
           typename  T_val >       // value type
class TBEMBF : public TBilinearForm< T_val >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = T_val;
        
protected:
    //! function space for ansatz functions
    const ansatzsp_t *  _ansatz_sp;
    
    //! function space for test functions
    const testsp_t *    _test_sp;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct bilinear form over function spaces \a ansatzsp and \a testsp
    TBEMBF ( const ansatzsp_t *  aansatzsp,
             const testsp_t *    atestsp )
            : _ansatz_sp(aansatzsp), _test_sp(atestsp)
    {
        if (( _ansatz_sp == NULL ) || ( _test_sp == NULL ))
            HERROR( ERR_ARG, "(TBEMBF)", "invalid function spaces (NULL)" );
    }

    //! destructor
    virtual ~TBEMBF () {}

    //////////////////////////////////////
    //
    // access internal variables
    //

    //! return ansatz space
    const ansatzsp_t *  ansatz_space () const { return _ansatz_sp; }

    //! return test space
    const testsp_t *    test_space   () const { return _test_sp; }

protected:
    DISABLE_COPY_OP( TBEMBF );
};

}// namespace Hpro

#endif  // __HPRO_TBEMBF_HH
