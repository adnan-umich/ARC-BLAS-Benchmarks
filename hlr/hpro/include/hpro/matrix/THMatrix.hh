#ifndef __HPRO_THMATRIX_HH
#define __HPRO_THMATRIX_HH
//
// Project     : HLIBpro
// File        : THMatrix.hh
// Description : class for H-matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <list>

#include "hpro/cluster/TPermutation.hh"
#include "hpro/matrix/TBlockMatrix.hh"
#include "hpro/parallel/NET.hh"

namespace Hpro
{

// local matrix type
DECLARE_TYPE( THMatrix );

//!
//! \ingroup  Matrix_Module
//! \class    THMatrix
//! \brief    Class for an H-matrix, which extends block matrices
//!           with additional functionality, e.g. permutations.
//!
template < typename T_value >
class THMatrix : public TBlockMatrix< T_value >
{
public:
    using  value_t = T_value;
    using  real_t  = typename real_type< value_t >::type_t;
    
private:

    //! @cond

    // maximal ID of all subblocks
    int            _max_id;
    
    // mappings of indices from external to
    // internal numbering (for mul_vec)
    TPermutation   _row_perm_e2i, _col_perm_e2i;
    TPermutation   _row_perm_i2e, _col_perm_i2e;

    //! @endcond
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    THMatrix ( const TBlockCluster *   cluster = nullptr );

    THMatrix ( const TBlockIndexSet &  bis );

    THMatrix ( const TIndexSet &  arowis,
               const TIndexSet &  acolis );

    virtual ~THMatrix () {}

    /////////////////////////////////////////////////
    //
    // access internal data
    //

    // access single matrix coefficient
    virtual value_t  entry  ( const idx_t i, const idx_t j ) const;

    // compute min/max tau/sigma indices
    void comp_min_max_idx ();

    // access permutations
    void set_row_perm ( const TPermutation & perm_e2i,
                        const TPermutation & perm_i2e )
    {
        _row_perm_e2i = perm_e2i;
        _row_perm_i2e = perm_i2e;
    }

    void set_col_perm ( const TPermutation & perm_e2i,
                        const TPermutation & perm_i2e )
    {
        _col_perm_e2i = perm_e2i;
        _col_perm_i2e = perm_i2e;
    }

    const TPermutation & row_perm_e2i () const { return _row_perm_e2i; }
    const TPermutation & row_perm_i2e () const { return _row_perm_i2e; }
    const TPermutation & col_perm_e2i () const { return _col_perm_e2i; }
    const TPermutation & col_perm_i2e () const { return _col_perm_i2e; }

    // return true if _all_ mappings are present and of correct size
    bool has_perm () const
    {
        return (( _row_perm_e2i.size() == this->rows() ) &&
                ( _row_perm_i2e.size() == this->rows() ) &&
                ( _col_perm_e2i.size() == this->cols() ) &&
                ( _col_perm_i2e.size() == this->cols() ));
    }

    // return maximal ID of subblocks
    int  max_id () const { return _max_id; }
    
    /////////////////////////////////////////////////
    //
    // BLAS-routines
    //

    //! compute y ≔ β·y + α·op(M)·x, with M = this
    virtual void mul_vec ( const value_t               alpha,
                           const TVector< value_t > *  x,
                           const value_t               beta,
                           TVector< value_t > *        y,
                           const matop_t               op = MATOP_NORM ) const;
    using TMatrix< value_t >::mul_vec;

    /////////////////////////////////////////////////
    //
    // matrix operations
    //

    //! transpose matrix
    virtual void transpose ();
    
    /////////////////////////////////////////////////
    //
    // misc.
    //

    //
    // virtual constructor
    //

    // return matrix of same class (but no content)
    virtual auto  create       () const -> std::unique_ptr< TMatrix< value_t > >
    {
        return std::make_unique< THMatrix >();
    }
    
    // return copy of matrix
    virtual auto  copy         () const -> std::unique_ptr< TMatrix< value_t > >;
    
    // copy matrix wrt. given accuracy and coarsening
    virtual auto  copy         ( const TTruncAcc &     acc,
                                 const bool            coarsen = false ) const -> std::unique_ptr< TMatrix< value_t > >;

    // return structural copy of matrix
    virtual auto  copy_struct  () const -> std::unique_ptr< TMatrix< value_t > >;
    
    // copy matrix into A
    virtual void  copy_to      ( TMatrix< value_t > *  A ) const;
    virtual void  copy_to      ( TMatrix< value_t > *  A,
                                 const TTruncAcc &     acc,
                                 const bool            coarsen = false ) const;

    // copy complete structural information from given matrix
    virtual void  copy_struct_from ( const TMatrix< value_t > * M );
    
    //! copy complete structural information from given matrix with different type
    template < typename T_value_M >
    void copy_struct_from_all ( const THMatrix< T_value_M > * M )
    {
        TBlockMatrix< value_t >::copy_struct_from_all( M );
        comp_min_max_idx();
        set_row_perm( M->row_perm_e2i(), M->row_perm_i2e() );
        set_col_perm( M->col_perm_e2i(), M->col_perm_i2e() );
    }

    // return appropriate vector-types for matrix
    virtual auto  row_vector   () const -> std::unique_ptr< TVector< value_t > >;
    virtual auto  col_vector   () const -> std::unique_ptr< TVector< value_t > >;
    
    //
    // size of object
    //
    
    // return size in bytes used by this object
    virtual size_t byte_size () const;

    //
    // RTTI
    //

    HPRO_RTTI_DERIVED( THMatrix, TBlockMatrix< value_t > )
    
    //
    // serialisation
    //

    virtual void read  ( TByteStream & s );
    virtual void build ( TByteStream & s );
    virtual void write ( TByteStream & s ) const;

    // returns size of object in bytestream
    virtual size_t bs_size () const;
};

}// namespace Hpro

#endif  // __HPRO_THMATRIX_HH
