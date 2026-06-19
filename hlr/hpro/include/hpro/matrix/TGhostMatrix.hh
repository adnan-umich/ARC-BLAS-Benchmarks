#ifndef __HPRO_TGHOSTMATRIX_HH
#define __HPRO_TGHOSTMATRIX_HH
//
// \file         TGhostMatrix.hh
//
// Project     : HLIBpro
// Description : class for representing remote matrix blocks
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"

namespace Hpro
{

// local matrix type
DECLARE_TYPE( TGhostMatrix );

//!
//! \ingroup Matrix_Module
//! \class   TGhostMatrix
//! \brief   The class acts as a place holder for non-local matrix blocks to
//!          access logical information, e.g. size, processor number, but
//!          can not perform any computations
//!
template < typename T_value >
class TGhostMatrix : public TMatrix< T_value >
{
public:
    using  value_t = T_value;
    using  real_t  = typename real_type< value_t >::type_t;

private:
    //! number of rows in matrix
    size_t  _n_rows;
    
    //! number of columns in matrix
    size_t  _n_cols;

public:
    ///////////////////////////////////////////////////////
    //
    // ctors and dtor
    //

    //! construct zero-sized matrix
    TGhostMatrix ()
            : TMatrix< value_t >()
            , _n_rows( 0 )
            , _n_cols(0)
    {}

    //! construct matrix defined over given block index set and
    //! on given processor set
    TGhostMatrix ( const TBlockIndexSet &  is,
                   const TProcSet &        ps )
            : TMatrix< value_t >()
            , _n_rows( is.row_is().size() )
            , _n_cols( is.col_is().size() )
    {
        this->set_ofs( is.row_is().first(), is.col_is().first() );
        this->set_procs( ps );
    }

    ////////////////////////////////////////////////////////
    //
    // structure of matrix
    //

    //! return number of rows
    virtual size_t  rows      () const { return _n_rows; }

    //! return number of columns
    virtual size_t  cols      () const { return _n_cols; }

    //! set dimension of matrix
    virtual void    set_size  ( const size_t n, const size_t m )
    {
        _n_rows = n;
        _n_cols = m;
    }

    /////////////////////////////////////////////////
    //
    // misc.
    //

    // transpose matrix
    virtual void   transpose  ()
    {
        TMatrix< value_t >::transpose();

        std::swap( _n_rows, _n_cols );
    }  
    
    // conjugate matrix coefficients
    virtual void   conjugate  ()
    {}
    
    //! truncate matrix to given accuracy
    virtual void   truncate   ( const TTruncAcc & )
    {}

    //! return matrix of same class (but no content)
    virtual auto   create     () const -> std::unique_ptr< TMatrix< value_t > >
    {
        return std::make_unique< TGhostMatrix< value_t > >();
    }

    //! return size in bytes used by this object
    virtual size_t byte_size  () const
    {
        return TMatrix< value_t >::byte_size() + sizeof(size_t) * 2;
    }

    //
    // RTTI
    //

    HPRO_RTTI_DERIVED( TGhostMatrix, TMatrix< value_t > )
};

}// namespace Hpro

#endif  // __HPRO_TGHOSTMATRIX_HH
