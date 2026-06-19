#ifndef __HPRO_TMASSBF_HH
#define __HPRO_TMASSBF_HH
//
// Project     : HLIBpro
// File        : TMassBF.hh
// Description : bilinear form for mass matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/bem/TBEMBF.hh"

namespace Hpro
{
    
//!
//! \class  TMassBF
//! \brief  bilinear form for mass matrix
//!
template < typename T_ansatzsp,
           typename T_testsp,
           typename T_value >
class TMassBF : public TBEMBF< T_ansatzsp,
                               T_testsp,
                               T_value >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = T_value;

protected:
    
    // quadrature rule
    std::vector< T2Point >  _quad_pts;
    std::vector< double >   _quad_wghts;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TMassBF ( const ansatzsp_t *  aansatzsp,
              const testsp_t *    atestsp,
              const uint          quad_order = CFG::BEM::quad_order );

    virtual ~TMassBF () {}

    // return format of bilinear form, e.g. symmetric
    matform_t  format () const;

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
};

}// namespace Hpro

#endif  // __HPRO_TMASSBF_HH
