#ifndef __HPRO_TZEROMATRIX_HH
#define __HPRO_TZEROMATRIX_HH
//
// Project     : HLIBpro
// File        : TZeroMatrix.hh
// Description : class for a zero matrix, i.e. with only zero coefficients
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/cluster/TPermutation.hh"

#include "hpro/matrix/TMatrix.hh"

namespace Hpro
{

// local matrix type
DECLARE_TYPE( TZeroMatrix );

//!
//! \ingroup  Matrix_Module
//! \class    TZeroMatrix
//! \brief    Class for a null matrix with only zero coefficients
//
template < typename T_value >
class TZeroMatrix : public TMatrix< T_value >
{
public:
    using  value_t = T_value;
    using  real_t  = typename real_type< value_t >::type_t;

private:
    //! number of rows in matrix
    size_t  _rows;

    //! number of columns in matrix
    size_t  _cols;

public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct null matrix of size \a anrows × \a ancols
    TZeroMatrix ( const size_t  anrows,
                  const size_t  ancols )
            : _rows(0)
            , _cols(0)
    {
        set_size( anrows, ancols );
    }
    
    //! construct null matrix with size defined by block cluster \a bct
    TZeroMatrix ( const TBlockCluster *  bct = nullptr )
            : TMatrix< value_t >( bct )
            , _rows(0)
            , _cols(0)
    {
        if ( bct != nullptr )
            set_cluster( bct );
    }

    virtual ~TZeroMatrix () {}

    /////////////////////////////////////////////////
    //
    // access data
    //

    //! set block cluster of matrix
    virtual void  set_cluster  ( const TBlockCluster * bct );

    //! directly set dimension of matrix
    virtual void  set_size     ( const size_t  nrows,
                                 const size_t  ncols );
    
    //! return number of rows in matrix
    size_t        rows         () const { return _rows; }

    //! return number of columns in matrix
    size_t        cols         () const { return _cols; }

    //! return true, if matrix is zero
    virtual bool  is_zero      () const { return true; }
    
    /////////////////////////////////////////////////
    //
    // manage stored entries
    //

    //! return matrix coefficient a_ij (real valued)
    virtual value_t  entry  ( const idx_t,
                              const idx_t ) const { return value_t(0); }

    /////////////////////////////////////////////////
    //
    // BLAS-routines (real valued)
    //

    //! compute this ≔ α·this
    virtual void scale ( const value_t ) {}
    
    //! compute y ≔ β·y + α·op(M)·x, with M = this
    virtual void mul_vec ( const value_t               alpha,
                           const TVector< value_t > *  x,
                           const value_t               beta,
                           TVector< value_t > *        y,
                           const matop_t               op = apply_normal ) const;
    using TMatrix< value_t >::mul_vec;

    //! compute this ≔ this + α · matrix
    virtual void add ( const value_t  alpha, const TMatrix< value_t > * matrix );
        
    //! transpose matrix
    virtual void transpose ();
    
    //! conjugate matrix coefficients
    virtual void conjugate ();
    
    //! truncate matrix to given accuracy (NOT YET IMPLEMENTED)
    virtual void truncate ( const TTruncAcc & ) {}
        
    /////////////////////////////////////////////////
    //
    // serialisation
    //

    //! read data from stream \a s and copy to matrix
    virtual void read  ( TByteStream & s );
    
    //! use data from stream \a s to build matrix
    virtual void build ( TByteStream & s );

    //! write data to stream \a s
    virtual void write ( TByteStream & s ) const;

    //! returns size of object in bytestream
    virtual size_t bs_size () const;

    /////////////////////////////////////////////////
    //
    // misc.
    //

    //
    // virtual ctors
    //
    
    //! return matrix of same class (but no content)
    virtual auto  create       () const -> std::unique_ptr< TMatrix< value_t > > { return std::make_unique< TZeroMatrix< value_t > >(); }
    
    //! return copy of matrix
    virtual auto  copy         () const -> std::unique_ptr< TMatrix< value_t > >;
    using TMatrix< value_t >::copy;

    //! return structural copy of matrix
    virtual auto  copy_struct  () const -> std::unique_ptr< TMatrix< value_t > >;

    //! copy matrix into matrix \a A
    virtual void copy_to       ( TMatrix< value_t > * A ) const;
    using TMatrix< value_t >::copy_to;

    //
    // type checking
    //

    HPRO_RTTI_DERIVED( TZeroMatrix, TMatrix< value_t > );

    //! return size in bytes used by this object
    virtual size_t byte_size () const;
};

}// namespace Hpro

#endif  // __HPRO_TZEROMATRIX_HH
