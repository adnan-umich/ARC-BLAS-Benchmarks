#ifndef __HPRO_BLAS_MATRIX_HH
#define __HPRO_BLAS_MATRIX_HH
//
// Project     : HLIBpro
// File        : Matrix.hh
// Description : implements dense matrix class for BLAS operations
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <type_traits>

#include "hpro/base/traits.hh"
#include "hpro/base/error.hh"

#include "hpro/blas/MemBlock.hh"
#include "hpro/blas/MatrixBase.hh"
#include "hpro/blas/Vector.hh"

#include "hpro/cluster/TPermutation.hh"

namespace Hpro
{

namespace BLAS
{

//!
//! \ingroup  BLAS_Module
//! \class    Matrix
//! \brief    Standard dense matrix in basic linear algebra, i.e. BLAS/LAPACK.
//!
//! \detail   Standard dmatrixvector class in linear algebra routines. It holds a
//!           reference to a memory block, which is controlled by this matrix,
//!           e.g. if the matrix is gone or resized, the memory block is also gone.
//!
template < typename T_value >
class Matrix : public MatrixBase< Matrix< T_value > >, public MemBlock< T_value >
{
    // ensure only floating point types (or complex version)
    static_assert( std::is_floating_point< T_value >::value || is_complex_type< T_value >::value,
                   "only floating point types supported" );
    
public:
    //! internal value type
    using  value_t = T_value;

    //! super class type
    using  super_t = MemBlock< value_t >;
    
private:
    //! dimensions of matrix
    size_t   _length[2];

    //! strides of data in memory block (rows and columns)
    size_t   _stride[2];

public:
    //
    // constructor and destructor
    //

    //! creates zero sized matrix
    Matrix () noexcept
            : MemBlock< value_t >()
            , _length{ 0, 0 }
            , _stride{ 0, 0 }
    {}

    //! creates matrix of size \a anrows × \a ancols
    Matrix ( const size_t  anrows,
             const size_t  ancols )
            : MemBlock< value_t >( anrows * ancols )
            , _length{ anrows, ancols }
            , _stride{ 1, anrows }
    {}

    //! copy constructor
    Matrix ( const Matrix &       M,
             const copy_policy_t  p = copy_reference )
            : MemBlock< value_t >()
            , _length{ 0, 0 }
            , _stride{ 0, 0 }
    {
        switch ( p )
        {
            case copy_reference :
                (*this) = M;
                break;

            case copy_value :
                _length[0] = M._length[0];
                _length[1] = M._length[1];
                _stride[0] = 1;
                _stride[1] = _length[0];
                super_t::alloc_wo_value( _length[0] * _length[1] );

                for ( idx_t j = 0; j < idx_t( _length[1] ); j++ )
                    for ( idx_t i = 0; i < idx_t( _length[0] ); i++ )
                        (*this)(i,j) = M( i, j );
                
                break;
        }// switch
    }

    //! copy constructor for other matrix types
    // template < typename T_matrix >
    // Matrix ( const T_matrix & M )
    //         : MemBlock<value_t>()
    //         , _length{ 0, 0 }
    //         , _stride{ 0, 0 }
    // {
    //     (*this) = M;
    // }
    
    //! move constructor
    Matrix ( Matrix &&  M ) noexcept
            : MemBlock< value_t >( std::move( M ) )
            , _length{ 0, 0 }
            , _stride{ 0, 0 }
    {
        M._data = nullptr;
        
        std::swap( _length, M._length );
        std::swap( _stride, M._stride );
    }

    //! creates matrix using part of \a M defined by \a r1 × \a r2
    //! \a p defines whether data is copied or referenced
    Matrix ( const Matrix &       M,
             const Range &        ar1,
             const Range &        ar2,
             const copy_policy_t  p = copy_reference )
            : MemBlock<value_t>()
            , _length{ 0, 0 }
            , _stride{ 0, 0 }
    {
        const Range  r1( ar1 == Range::all ? M.row_range() : ar1 );
        const Range  r2( ar2 == Range::all ? M.col_range() : ar2 );
        
        _length[0] = r1.size() / r1.stride();
        _length[1] = r2.size() / r2.stride();

        switch ( p )
        {
            case copy_reference :
                _stride[0] = r1.stride() * M.row_stride();
                _stride[1] = r2.stride() * M.col_stride();
            
                super_t::init( M.data() + r1.first() * M.row_stride() + r2.first() * M.col_stride() );
                break;

            case copy_value :
                super_t::alloc_wo_value( _length[0] * _length[1] );
                _stride[0] = 1;
                _stride[1] = _length[0];

                for ( idx_t j = 0; j < idx_t( _length[1] ); j++ )
                    for ( idx_t i = 0; i < idx_t( _length[0] ); i++ )
                        (*this)(i,j) = M( r1.first() + i * idx_t( r1.stride() ),
                                          r2.first() + j * idx_t( r2.stride() ) );
                break;
        }// switch
    }

    //! direct setting of matrix from raw data pointer \a ptr, \a nrows, \a rowstride, \a ncols and \a colstride
    Matrix ( value_t *            aptr,
             const size_t         anrows,
             const size_t         arowstride,
             const size_t         ancols,
             const size_t         acolstride,
             const copy_policy_t  p = copy_reference )
            : MemBlock< value_t >()
            , _length{ anrows, ancols }
            , _stride{ arowstride, acolstride }
    {
        switch ( p )
        {
            case copy_reference :
                super_t::init( aptr, false );
                break;
                
            case copy_value :
                super_t::alloc_wo_value( anrows * ancols );
                _stride[0] = 1;
                _stride[1] = anrows;
        
                for ( idx_t j = 0; j < idx_t( _length[1] ); j++ )
                    for ( idx_t i = 0; i < idx_t( _length[0] ); i++ )
                        (*this)(i,j) = aptr[ j * acolstride + i * arowstride ];
                break;
        }// switch
    }

    //! copy operator for matrices (always copy reference! for real copy, use BLAS::copy)
    Matrix & operator = ( const Matrix &  M )
    {
        super_t::init( M.data(), false );
        _length[0] = M.nrows();
        _length[1] = M.ncols();
        _stride[0] = M.row_stride();
        _stride[1] = M.col_stride();

        return *this;
    }

    //! move operator
    Matrix & operator = ( Matrix &&  M ) noexcept
    {
        if ( this != & M ) // prohibit self-moving
        {
            super_t::init( M, M._is_owner );
            
            _length[0] = M.nrows();
            _length[1] = M.ncols();
            _stride[0] = M.row_stride();
            _stride[1] = M.col_stride();

            M._data   = nullptr;
            M._length[0] = 0;
            M._length[1] = 0;
            M._stride[0] = 0;
            M._stride[1] = 0;
        }// if

        return *this;
    }
    
    //
    // data access
    //

    //! return number of rows of matrix
    size_t       nrows        () const noexcept { return _length[0]; }
    
    //! return number of columns of matrix
    size_t       ncols        () const noexcept { return _length[1]; }

    //! return BLAS matrix view of matrix object
    blasview_t   blas_view    () const noexcept { return BLAS_NORMAL; }

    //! return number of rows of actual BLAS matrix
    size_t       blas_nrows   () const noexcept { return nrows(); }

    //! return number of columns of actual BLAS matrix
    size_t       blas_ncols   () const noexcept { return ncols(); }

    //! return coefficient (i,j)
    value_t      operator ()  ( const idx_t i, const idx_t j ) const noexcept
    {
        HASSERT( i < idx_t(_length[0]) && j < idx_t(_length[1]), ERR_ARR_BOUND, "(Matrix) operator ()", "" );
        return super_t::_data[ j * _stride[1] + i * _stride[0] ];
    }

    //! return reference to coefficient (i,j)
    value_t &    operator ()  ( const idx_t i, const idx_t j ) noexcept
    {
        HASSERT( i < idx_t(_length[0]) && j < idx_t(_length[1]), ERR_ARR_BOUND, "(Matrix) operator ()", "" );
        return super_t::_data[ j * _stride[1] + i * _stride[0] ];
    }

    //! return pointer to coefficient (i,j)
    value_t *    ptr          ( const idx_t i, const idx_t j ) noexcept
    {
        HASSERT( i < idx_t(_length[0]) && j < idx_t(_length[1]), ERR_ARR_BOUND, "(Matrix) operator ()", "" );
        return super_t::_data + ( j * _stride[1] + i * _stride[0] );
    }
    const value_t *  ptr      ( const idx_t i, const idx_t j ) const noexcept
    {
        HASSERT( i < idx_t(_length[0]) && j < idx_t(_length[1]), ERR_ARR_BOUND, "(Matrix) operator ()", "" );
        return super_t::_data + ( j * _stride[1] + i * _stride[0] );
    }

    //! return pointer to internal data
    value_t *    data         () const noexcept { return super_t::_data; }

    //! return stride w.r.t. row index set
    size_t       row_stride   () const noexcept { return _stride[0]; }

    //! return stride w.r.t. column index set
    size_t       col_stride   () const noexcept { return _stride[1]; }

    //! optimised resize: only change if (n,m) != (nrows,ncols)
    void         resize       ( const size_t  n,
                                const size_t  m )
    {
        if (( _length[0] != n ) || ( _length[1] != m ))
        {
            *this = std::move( Matrix( n, m ) );
        }// if
    }
    
    //
    // construction operators
    //

    //! create real copy of matrix
    Matrix< value_t >  copy () const
    {
        Matrix< value_t >  M( *this, copy_value );

        return M;
    }
    
    //! create reference to this matrix
    Matrix< value_t >  reference () const
    {
        Matrix< value_t >  M( *this, copy_reference );

        return M;
    }
    
    //! return matrix referencing sub matrix defined by \a r1 × \a r2
    Matrix< value_t >  operator () ( const Range & r1, const Range & r2 ) const
    {
        return Matrix< value_t >( *this, r1, r2 );
    }
                          
    //! return vector referencing part of column \a j defined by \a r
    Vector< value_t >  operator () ( const Range & r, const idx_t  j ) const
    {
        return Vector< value_t >( *this, r, j );
    }
                          
    //! return vector referencing part of row \a i defined by \a r
    Vector< value_t >  operator () ( const idx_t  i, const Range & r ) const
    {
        return Vector< value_t >( *this, i, r );
    }
                          
    //! return vector referencing column \a j
    Vector< value_t >  column   ( const idx_t  j ) const
    {
        return (*this)( Range( 0, idx_t(nrows())-1 ), j );
    }
    
    //! return vector referencing row \a i
    Vector< value_t >  row      ( const idx_t  i ) const
    {
        return (*this)( i, Range( 0, idx_t(ncols())-1 ) );
    }
    
    //! return size in bytes used by this object
    size_t  byte_size () const
    {
        return sizeof( value_t ) * _length[0] * _length[1] + sizeof(_length) + sizeof(_stride) + sizeof(MemBlock< value_t >);
    }

    //! return size of (floating point) data in bytes handled by this object
    size_t  data_byte_size () const
    {
        return sizeof( value_t ) * _length[0] * _length[1];
    }

    //
    // tests
    //

    //! test data for invalid values, e.g. INF and NAN
    void  check_data  () const;
};

//
// trait for giving access to matrix properties
//
template < typename T >
struct matrix_trait< Matrix< T > >
{
    using  value_t = T;
};

//
// signals, that T is of matrix type
//
template < typename T >
struct is_matrix< Matrix< T > >
{
    static const bool  value = true;
};

//
// return real copy of given matrix
//
template < typename value_t >
Matrix< value_t >
copy ( const Matrix< value_t > &  M )
{
    return M.copy();
}

//!
//! convert between different dataypes 
//!
template < typename dest_value_t,
           typename src_value_t >
Matrix< dest_value_t >
convert ( const Matrix< src_value_t > &  M )
{
    // prevent complex to real conversion
    static_assert( ! is_complex_type< src_value_t >::value || is_complex_type< dest_value_t >::value,
                   "can not convert complex to real type" );
    
    if constexpr ( std::is_same< dest_value_t, src_value_t >::value )
        return M.copy();
    else
    {
        auto  T = Matrix< dest_value_t >( M.nrows(), M.ncols() );

        const size_t  n = M.nrows() * M.ncols();

        for ( size_t  i = 0; i < n; ++i )
            T.data()[i] = dest_value_t( M.data()[i] );

        return T;
    }// else
}

//
// return reference to matrix
//
template < typename T >
Matrix< T >
reference ( const Matrix< T > &  M )
{
    return M.reference();
}

//
// create identity matrix
//
template < typename T >
Matrix< T >
identity ( const size_t  n );

//
// create random matrix
//
template < typename T >
Matrix< T >
random ( const size_t  nrows,
         const size_t  ncols );

//
// permute matrix
//
template < typename value_t >
void
permute ( Matrix< value_t > &   M,
          const TPermutation &  row_perm,
          const TPermutation &  col_perm );

}// namespace BLAS

//////////////////////////////////////////////////////////
//
// debug helpers
//

namespace DBG
{

//!
//! write matrix to file
//!
template <typename T>
void
write ( const BLAS::Matrix< T > &  M,
        const std::string &        filename,
        const std::string &        matname );

}// namespace DBG

}// namespace Hpro

//
// include matrix views
//
#include "hpro/blas/matrix_view.hh"

#endif  // __HPRO_BLAS_MATRIX_HH
