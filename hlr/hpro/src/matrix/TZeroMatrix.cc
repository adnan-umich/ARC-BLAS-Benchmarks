//
// Project     : HLIBpro
// File        : TZeroMatrix.cc
// Description : class for a zero matrix, i.e. with only zero coefficients
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/matrix/TZeroMatrix.hh"

namespace Hpro
{

/////////////////////////////////////////////////
//
// access data
//

//! set block cluster of matrix
template < typename value_t >
void
TZeroMatrix< value_t >::set_cluster  ( const TBlockCluster * bct )
{
    TMatrix< value_t >::set_cluster( bct );

    if ( bct != nullptr )
        set_size( bct->rowcl()->size(),
                  bct->colcl()->size() );
}

//! directly set dimension of matrix
template < typename value_t >
void
TZeroMatrix< value_t >::set_size ( const size_t  nrows,
                        const size_t  ncols )
{
    _rows = nrows;
    _cols = ncols;
}
    
/////////////////////////////////////////////////
//
// BLAS-routines
//

//! compute y ≔ β·y + α·op(M)·x, with M = this
template < typename value_t >
void
TZeroMatrix< value_t >::mul_vec ( const value_t              ,
                                  const TVector< value_t > * ,
                                  const value_t              beta,
                                  TVector< value_t >       * y,
                                  const matop_t               ) const
{
    y->scale( beta );
}

//! compute this ≔ this + α · matrix
template < typename value_t >
void
TZeroMatrix< value_t >::add ( const value_t, const TMatrix< value_t > * M )
{
    if ( ! IS_TYPE( M, TZeroMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TZeroMatrix) cadd", M->typestr() );
}

//
// transpose matrix
//
template < typename value_t >
void
TZeroMatrix< value_t >::transpose ()
{
    TMatrix< value_t >::transpose();
    
    std::swap( _rows, _cols );
}

//
// conjugate matrix coefficients
//
template < typename value_t >
void
TZeroMatrix< value_t >::conjugate ()
{
}

//
// serialisation
//

//! read data from stream \a s and copy to matrix
template < typename value_t >
void
TZeroMatrix< value_t >::read  ( TByteStream & s )
{
    TMatrix< value_t >::read( s );

    s.get( _rows );
    s.get( _cols );
}
    
//! use data from stream \a s to build matrix
template < typename value_t >
void
TZeroMatrix< value_t >::build ( TByteStream & s )
{
    TMatrix< value_t >::build( s );

    s.get( _rows );
    s.get( _cols );
}

//! write data to stream \a s
template < typename value_t >
void
TZeroMatrix< value_t >::write ( TByteStream & s ) const
{
    TMatrix< value_t >::write( s );

    s.put( _rows );
    s.put( _cols );
}

//! returns size of object in bytestream
template < typename value_t >
size_t
TZeroMatrix< value_t >::bs_size () const
{
    return TMatrix< value_t >::bs_size() + + sizeof(_rows) + sizeof(_cols);
}

/////////////////////////////////////////////////
//
// misc.
//

//
// return copy of matrix
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TZeroMatrix< value_t >::copy  () const
{
    auto           M = TMatrix< value_t >::copy();
    TZeroMatrix< value_t > *  Z = ptrcast( M.get(), TZeroMatrix< value_t > );

    copy_to( Z );

    return M;
}

//
// return structural copy of matrix
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TZeroMatrix< value_t >::copy_struct  () const
{
    auto  M = TMatrix< value_t >::copy_struct();
    auto  Z = ptrcast( M.get(), TZeroMatrix< value_t > );

    copy_to( Z );

    return M;
}

//
// copy matrix into matrix \a A
//
template < typename value_t >
void
TZeroMatrix< value_t >::copy_to ( TMatrix< value_t > * A ) const
{
    if ( ! IS_TYPE( A, TZeroMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TZeroMatrix) copy_to", A->typestr() );

    auto  M = ptrcast( A, TZeroMatrix< value_t > );

    M->set_size( rows(), cols() );
    M->set_ofs( this->row_ofs(), this->col_ofs() );
}

//
// return size in bytes used by this object
//
template < typename value_t >
size_t
TZeroMatrix< value_t >::byte_size () const
{
    return TMatrix< value_t >::byte_size() + sizeof(_rows) + sizeof(_cols);
}

template class TZeroMatrix< float >;
template class TZeroMatrix< double >;
template class TZeroMatrix< std::complex< float > >;
template class TZeroMatrix< std::complex< double > >;

}// namespace Hpro
