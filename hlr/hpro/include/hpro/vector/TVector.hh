#ifndef __HPRO_TVECTOR_HH
#define __HPRO_TVECTOR_HH
//
// Project     : HLIBpro
// File        : TVector.hh
// Description : baseclass for all vector-classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <variant>

#include "hpro/base/traits.hh"
#include "hpro/base/System.hh"
#include "hpro/base/TStreamable.hh"
#include "hpro/base/TTypeInfo.hh"
#include "hpro/cluster/TIndexSet.hh"

namespace Hpro
{

//////////////////////////////////////////////////////////
//
// base vector class
//

// local matrix type
DECLARE_TYPE( TVector );

//!
//! \ingroup Vector_Module
//! \class   TVector
//! \brief   Base class for all vectors defining basic interface.
//!
template < typename T_value >
class TVector : public TStreamable, public TTypeInfo
{
public:
    using  value_t = T_value;
    using  real_t  = typename real_type< value_t >::type_t;
    
protected:
    //@cond
    
    //! first index vector represents
    idx_t  _ofs;

    //@endcond
    
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct real or complex valued vector with first index \a offset
    TVector ( const idx_t  offset      = 0 )
            : _ofs(offset)
    {}

    //! copy constructor
    TVector ( const TVector< value_t > &  v )
            : _ofs(0)
    {
        assign( 1.0, & v );
    }

    //! dtor
    virtual ~TVector () {}

    ////////////////////////////////////////////////
    //! \{
    //!
    //! \name value type information
    //!
    
    //! return true, if field type is complex valued
    virtual bool  is_complex  () const { return is_complex_type< value_t >::value; }
    
    //! return true, if field type is real valued
    virtual bool  is_real     () const { return ! is_complex(); }
    
    //!\}
    
    ////////////////////////////////////////////////
    //! \{
    //!
    //! \name index set functionality
    //!

    //! return first index (offset)
    idx_t           ofs      () const { return _ofs; }

    //! set first index (offset)
    virtual void    set_ofs  ( const idx_t  n ) { _ofs = n; }
    
    //! return size of vector
    virtual size_t  size     () const = 0;

    //! return index set
    TIndexSet       is       () const { return TIndexSet( ofs(), idx_t(ofs() + size()) - 1 ); }

    //!\}

public:
    ////////////////////////////////////////////////
    //! \{
    //!
    //! \name access entries
    //!

    //! return i'th entry
    virtual value_t  entry  ( const idx_t  i ) const;

    //! set \a i'th entry
    virtual void set_entry  ( const idx_t  i, const value_t  f );

    //! \}

    ////////////////////////////////////////////////
    //! \{
    //!
    //! \name BLAS functionality
    //!

    //! conjugate entries
    virtual void conjugate ()
    { HERROR( ERR_NOT_IMPL, "(TVector) conjugate", "" ); }
        
    //! fill with constant
    virtual void fill ( const value_t )
    { HERROR( ERR_NOT_IMPL, "(TVector) fill", "" ); }

    //! fill with random numbers
    virtual void fill_rand ( const uint )
    { HERROR( ERR_NOT_IMPL, "(TVector) fill_rand", "" ); }

    //! scale vector by constant factor
    virtual void scale ( const value_t )
    { HERROR( ERR_NOT_IMPL, "(TVector) scale", "" ); }

    //! this ≔ f · vector
    virtual void assign ( const value_t, const TVector< value_t > * )
    { HERROR( ERR_NOT_IMPL, "(TVector) assign", "" ); }

    //! copy operator for all vectors
    TVector< value_t > &  operator = ( const TVector< value_t > & v )
    {
        assign( value_t(1), & v );
        return *this;
    }
    
    //! return dot-product, \f$<x,y> = x^H · y\f$, where \f$x\f$ = this
    virtual value_t dot  ( const TVector< value_t > * ) const
    { HERROR( ERR_NOT_IMPL, "(TVector) dot", "" ); }

    //! return dot-product, \f$<x,y> = x^T · y\f$, where \f$x\f$ = this
    virtual value_t dotu ( const TVector< value_t > * ) const
    { HERROR( ERR_NOT_IMPL, "(TVector) dotu", "" ); }

    //! return euclidean norm
    virtual real_t  norm2 () const
    {
        if constexpr ( is_complex_type< value_t >::value )
            return Math::sqrt( std::real( dot( this ) ) );
        else
            return Math::sqrt( dot( this ) );
    }

    //! return infimum norm
    virtual real_t  norm_inf () const { HERROR( ERR_NOT_IMPL, "(TVector) norm_inf", "" ); }
    
    //! this ≔ this + α·x
    virtual void axpy ( const value_t, const TVector< value_t > * )
    { HERROR( ERR_NOT_IMPL, "(TVector) axpy", "" ); }
    
    //! \}
    
    ////////////////////////////////////////////////
    //! \{
    //!
    //! \name misc.
    //!

    //
    // memory consumption
    //
    
    //! return size in bytes used by this object
    virtual size_t  byte_size        () const { return sizeof(_ofs); }

    //! return size of (floating point) data in bytes handled by this object
    virtual size_t data_byte_size    () const { return 0; }

    //! return size in bytes used by this distributed object,
    //! i.e. of all distributed sub vectors
    virtual size_t  global_byte_size () const;

    //
    // virtual constructor
    //

    //! return copy of vector
    virtual auto  copy    () const -> std::unique_ptr< TVector< value_t > > = 0;
    
    //! assign local values to vector \a x
    virtual void  copy_to ( TVector< value_t > * x ) const
    {
        if ( x == nullptr )
            HERROR( ERR_ARG, "(TVector) copy_to", "vector is NULL" );
        x->assign( 1.0, this );
    }
    
    //! return object of same class
    virtual auto  create  () const -> std::unique_ptr< TVector< value_t > > = 0;

    //
    // restriction
    //

    //! create vector restricted to real part of coefficients
    virtual auto restrict_re  () const -> std::unique_ptr< TVector< real_t > >;

    //! create vector restricted to imaginary part of coefficients
    virtual auto restrict_im  () const -> std::unique_ptr< TVector< real_t > >;
    
    //! create vector restricted to absolute value of coefficients
    virtual auto restrict_abs () const -> std::unique_ptr< TVector< real_t > >;
    
    //
    // serialisation
    //

    //! read vector data from byte stream
    virtual void read  ( TByteStream & s );

    //! write vector data to byte stream
    virtual void write ( TByteStream & s ) const;

    //! returns size of object in bytestream
    virtual size_t bs_size () const;

    //
    // parallel methods
    //

    //! sum up nparts parallel copies
    //! (if bs != NULL it will be used)
    virtual void sum ( const TProcSet & p,
                       const uint       pid,
                       const uint       nparts,
                       TByteStream *    bs = NULL );

    //! same as \see sum but sum up between all processors in \a p
    virtual void sum ( const TProcSet & p );
    
    //
    // output
    //

    //! print vector to stdout
    virtual void print ( const uint ofs = 0 ) const;

    HPRO_RTTI_BASE( TVector )
    
    //! \}
};

//////////////////////////////////////////////////////////
//
// variant version of vectors
//

using any_vector_t = std::variant<
    TVector< float > *,
    TVector< double > *,
    TVector< std::complex< float > > *,
    TVector< std::complex< double > > * >;
    
using any_const_vector_t = std::variant<
    const TVector< float > *,
    const TVector< double > *,
    const TVector< std::complex< float > > *,
    const TVector< std::complex< double > > * >;

//////////////////////////////////////////////////////////
//
// wrappers for vector functions
//

//! return dot product <x,y> = x^H · y
template < typename value_t >
value_t  dot  ( const TVector< value_t > *  x,
                const TVector< value_t > *  y )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "dot", "x is null" );
    
    return x->dot( y );
}

template < typename value_t >
value_t  dot  ( const TVector< value_t > &  x,
                const TVector< value_t > &  y )
{
    return x.dot( &y );
}

//! return dot product <x,y> = x^T · y
template < typename value_t >
value_t  dotu  ( const TVector< value_t > *  x,
                 const TVector< value_t > *  y )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "dotu", "x is null" );
    
    return x->dotu( y );
}

template < typename value_t >
value_t  dotu  ( const TVector< value_t > &  x,
                 const TVector< value_t > &  y )
{
    return x.dotu( &y );
}

//! return euclidean norm of \a x
template < typename value_t >
typename real_type< value_t >::type_t
norm_2  ( const TVector< value_t > *  x )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "norm_2", "x is null" );
    
    return x->norm2();
}

template < typename value_t >
typename real_type< value_t >::type_t
norm_2  ( const TVector< value_t > &  x )
{
    return x.norm2();
}

//! return infimum norm of \a x
template < typename value_t >
typename real_type< value_t >::type_t
norm_inf  ( const TVector< value_t > *  x )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "norm_inf", "x is null" );
    
    return x->norm_inf();
}

template < typename value_t >
typename real_type< value_t >::type_t
norm_inf  ( const TVector< value_t > &  x )
{
    return x.norm_inf();
}

//////////////////////////////////////////////////////////
//
// template wrappers
//

//
// scale
//
template < typename value_t >
void
scale ( const value_t         alpha,
        TVector< value_t > *  x )
{
    x->scale( alpha );
}

template < typename value_t >
void
scale ( const value_t         alpha,
        TVector< value_t > &  x )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "scale", "x is null" );
    
    x->scale( alpha );
}

//
// axpy
//
template < typename value_t >
void
axpy ( const value_t               alpha,
       const TVector< value_t > *  x,
       TVector< value_t > *        y )
{
    if ( y == nullptr )
        HERROR( ERR_ARG, "axpy", "y is null" );
    
    y->axpy( alpha, x );
}

template < typename value_t >
void
axpy ( const value_t               alpha,
       const TVector< value_t > &  x,
       TVector< value_t > &        y )
{
    y.axpy( alpha, & x );
}

template < typename value_t >
void
add ( const value_t               alpha,
      const TVector< value_t > &  x,
      TVector< value_t > &        y )
{
    y.axpy( alpha, & x );
}

//
// copy vectors
//
template < typename value_t >
std::unique_ptr< TVector< value_t > >
copy ( const TVector< value_t > *  v )
{
    if ( v == nullptr )
        HERROR( ERR_ARG, "copy", "vector is null" );
    
    return v->copy();
}

template < typename value_t >
std::unique_ptr< TVector< value_t > >
copy ( const TVector< value_t > &  v )
{
    return v.copy();
}

//////////////////////////////////////////////////////////
//
// debug helpers
//

namespace DBG
{

//
// write vector to file
//
template < typename value_t >
void write ( const TVector< value_t > *  v,
             const std::string &         filename,
             const std::string &         vecname );

template < typename value_t >
void write ( const TVector< value_t > &  v,
             const std::string &         filename,
             const std::string &         vecname )
{
    write( &v, filename, vecname );
}

}// namespace DBG

}// namespace Hpro

#endif  // __HPRO_TVECTOR_HH
