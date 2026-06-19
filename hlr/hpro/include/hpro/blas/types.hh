#ifndef __HPRO_BLAS_TYPES_HH
#define __HPRO_BLAS_TYPES_HH
//
// Project     : HLIBpro
// File        : types.hh
// Description : provide basic types for BLAS operations
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

namespace Hpro
{

//! \enum  matop_t
//! \brief matrix transformations, e.g. transposed or conjugate transpose
enum matop_t
{
    //! do not change matrix
    apply_normal     = 'N',
    MATOP_NORM       = 'N',

    //! use transposed matrix
    apply_transposed = 'T',
    apply_trans      = 'T',
    MATOP_TRANS      = 'T',

    //! use adjoint, e.g. conjugate transposed matrix
    apply_adjoint    = 'C',
    apply_adj        = 'C',
    apply_conjtrans  = 'C',       
    MATOP_ADJ        = 'C',        
    MATOP_CONJTRANS  = 'C',        

    //! use conjugate matrix (NOT WELL SUPPORTED!!!)
    apply_conjugate  = 'R',
    apply_conj       = 'R',
    MATOP_CONJ       = 'R'
};

//! \enum   approx_t
//! \brief  different types of low-rank approximation methods
enum approx_t : int
{
    use_svd      = 0,
    use_rrqr     = 1,
    use_rand     = 2,
    use_randsvd  = 2,
    use_randlr   = 3,
    use_aca      = 4,
    use_svd_pair = 5
};

namespace BLAS
{

//! \enum   diag_type_t
//! \brief  defines whether diagonal is unit or non-unit
enum diag_type_t
{
    unit_diag        = 'U',
    general_diag     = 'N'
};

//! \enum   tri_type_t
//! \brief  defines whether matrix is upper or lower triangular
enum tri_type_t
{
    lower_triangular = 'L',
    upper_triangular = 'U'
};

//! \enum   eval_side_t
//! \brief  defines whether triangular system is multiplied from left or right
enum eval_side_t
{
    from_left        = 'L',
    from_right       = 'R'
};

}// namespace BLAS

}// namespace Hpro

#endif  // __HPRO_BLAS_TYPES_HH
