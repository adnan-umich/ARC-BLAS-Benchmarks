#ifndef __HPRO_TSCALARVECTOR_HH
#define __HPRO_TSCALARVECTOR_HH
//
// Project     : HLIBpro
// File        : TScalarVector.hh
// Description : class for a vector of scalar type
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/types.hh"
#include "hpro/blas/Vector.hh"
#include "hpro/blas/Algebra.hh"

#include "hpro/cluster/TCluster.hh"
#include "hpro/cluster/TPermutation.hh"

#include "hpro/vector/TVector.hh"

namespace Hpro
{

// local matrix type
DECLARE_TYPE( TScalarVector );

// chunk size for updates to scalar vectors
const size_t  SCALAR_CHUNK_SIZE = 64;
    
//!
//! \ingroup Vector_Module
//! \class   TScalarVector 
//! \brief   Class for a scalar vector.
//!
template < typename T_value >
class TScalarVector : public TVector< T_value >
{
public:
    using  value_t = T_value;
    using  real_t  = typename real_type< value_t >::type_t;
    
protected:
    //! @cond
    
    //! real valued vector data
    BLAS::Vector< value_t >  _vec;

    //! size of vector
    size_t                   _size;

    //! @endcond
    
public:
    //////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //
    // constructors without copy/move
    //
    
    //! construct zero sized vector
    TScalarVector ()
            : TVector< value_t >( 0 )
            , _size( 0 ) 
    {}
    
    //! construct vector of size \a n with offset \a offset
    TScalarVector ( const size_t  n,
                    const idx_t   offset      = 0 )
            : TVector< value_t >( offset )
            , _size( 0 )
    {
        set_size( n );
    }
    
    //! construct vector with size defined by indexset \a ais
    TScalarVector ( const TIndexSet &  ais )
            : TVector< value_t >( ais.first() )
            , _size( 0 )
    {
        set_size( ais.size() );
    }
    
    //
    // copy constructors
    //
    
    //! construct vector with size defined by indexset \a ais
    //! and data defined by \a bvec
    TScalarVector ( const TIndexSet &                ais,
                    const BLAS::Vector< value_t > &  bvec )
            : TVector< value_t >( ais.first() )
            , _vec( bvec )
            , _size( bvec.length() )
    {
        init_chunk_mutices();
    }
    
    //! standard copy constructor
    TScalarVector ( const TScalarVector< value_t > &  v )
            : TVector< value_t >()
            , _size( 0 )
    {
        assign( value_t(1), & v );
    }

    //! standard move constructor
    TScalarVector ( TScalarVector< value_t > &&       v )
            : TVector< value_t >( v.ofs() )
            , _size( 0 )
    {
        _size = v._vec.length();
        _vec  = std::move( v._vec );
    }

    //
    // move constructors
    //

    //! construct vector with size defined by indexset \a ais
    //! and data defined by \a bvec
    TScalarVector ( const TIndexSet &           ais,
                    BLAS::Vector< value_t > &&  bvec )
            : TVector< value_t >( ais.first() )
            , _vec( std::move( bvec ) )
            , _size( _vec.length() )  // use _rvec, because bvec is now "empty"
    {
        init_chunk_mutices();
    }

    //! destructor
    virtual ~TScalarVector ()
    {}

    //////////////////////////////////////////////////
    //
    // access vector coefficients
    //

    //! return size of vector
    virtual size_t          size        () const
    {
        return _size;
    }

    //! access coefficent \a i (real valued)
    virtual value_t         entry       ( const idx_t  i ) const
    {
        HASSERT( i < idx_t(_size), ERR_ARR_BOUND, "(TScalarVector) entry", "" );
        return _vec(i);
    }
    
    //! set coefficient \a i to \a f (real valued)
    virtual void            set_entry   ( const idx_t  i, const value_t  f )
    {
        HASSERT( i < idx_t(_size), ERR_ARR_BOUND, "(TScalarVector) set_entry", "" );
        _vec(i) = f;
    }
    
    //! add \a f to \a i'th entry
    virtual void            add_entry   ( const idx_t  i, const value_t  f )
    {
        HASSERT( i < idx_t(_size), ERR_ARR_BOUND, "(TScalarVector) add_entry", "" );
        _vec(i) += f;
    }

    //
    // access internal vectors as BLAS vectors
    //

    //! return real valued data
    BLAS::Vector< value_t > &        blas_vec  ()       { return _vec; }

    //! return constant real valued data
    const BLAS::Vector< value_t > &  blas_vec  () const { return _vec; }

    //
    // handle index set
    //
    
    //! set size of vector
    virtual void set_size    ( const size_t  n );
    
    //! define vector by cluster
    virtual void set_cluster ( const TCluster *  c )
    {
        if ( c != NULL )
        {
            set_size( c->size() );
            this->set_ofs( c->first() );
        }// if
    }

    //! define vector by indexset
    virtual void set_is ( const TIndexSet &  ais )
    {
        set_size( ais.size() );
        this->set_ofs( ais.first() );
    }

    //
    // copy/assign methods
    //
    
    //! set internal data directly (real valued)
    virtual void set_vector ( const BLAS::Vector< value_t > &  vec,
                              const idx_t                      offset )
    {
        set_size( vec.length() );
        this->set_ofs( offset );
        BLAS::copy( vec, _vec );
    }

    //! copy from vector \a v
    virtual void copy_from ( const TScalarVector< value_t > *  v )
    {
        if ( v == NULL )
            HERROR( ERR_ARG, "(TScalarVector) copy_from", "v is NULL" );

        set_vector( v->blas_vec(), v->ofs() );
    }

    //! copy to vector \a v
    virtual void copy_to ( TScalarVector< value_t > *  v ) const
    {
        if ( v == NULL )
            HERROR( ERR_ARG, "(TScalarVector) copy_to", "v is NULL" );

        v->copy_from( this );
    }

    //
    // return reference to sub vector
    //

    TScalarVector< value_t >
    sub_vector ( const TIndexSet &  ais )
    {
        if ( ! this->is().is_sub( ais ) )
            HERROR( ERR_INDEXSET, "(TScalarVector) sub_vector",
                    "given index set is NOT a sub set of local index set" );

        return TScalarVector( ais, BLAS::Vector< value_t >( _vec, ais - this->ofs() ) );
    }
    
    const TScalarVector< value_t >
    sub_vector ( const TIndexSet &  ais ) const
    {
        if ( ! this->is().is_sub( ais ) )
            HERROR( ERR_INDEXSET, "(TScalarVector) sub_vector",
                    "given index set is NOT a sub set of local index set" );

        return TScalarVector( ais, BLAS::Vector< value_t >( _vec, ais - this->ofs() ) );
    }
    
    
    //! copy from C array \a v
    virtual void copy_from ( const value_t *  v );

    //! copy to C array \a v
    virtual void copy_to   ( value_t *        v );
    using TVector< value_t >::copy_to;

    //! standard copy operator
    TScalarVector< value_t > &  operator = ( const TScalarVector< value_t > &  v )
    {
        assign( value_t(1), & v );
        return * this;
    }

    //! permute entries according to \a perm
    void permute ( const TPermutation &  perm );

protected:

    //
    // handle chunk mutex array
    //

    //! initialise mutices for vector chunks
    void init_chunk_mutices ()
    {
        // _mutices = std::move( std::vector< mutex_t >( size() / SCALAR_CHUNK_SIZE + 1 ) );
        // // _mutices.resize( size() / SCALAR_CHUNK_SIZE + 1 ); // +1 for non-multiples
    }

public:
    //////////////////////////////////////////////////
    //
    // BLAS-routines
    //

    //! conjugate coefficients
    virtual void     conjugate  ()
    {
        if ( is_complex_type< value_t >::value )
            BLAS::conj( blas_vec() );
    }
    
    //! fill vector with constant \a α
    virtual void     fill        ( const value_t  alpha );

    //! fill vector with random numbers
    virtual void     fill_rand   ( const uint     seed );

    //! set this ≔ α · this
    virtual void     scale       ( const value_t  alpha );

    //! set this ≔ α · x
    virtual void     assign      ( const value_t               alpha,
                                   const TVector< value_t > *  x );

    //! return inner product <this, x> = this^H · x
    virtual value_t  dot         ( const TVector< value_t > *  x ) const;

    //! return inner product <this, x> = this^T · x
    virtual value_t  dotu        ( const TVector< value_t > *  x ) const;

    //! compute ‖·‖₂
    virtual real_t   norm2       () const;

    //! compute ‖·‖∞
    virtual real_t   norm_inf    () const;
    
    //! set this ≔ this + α · x
    virtual void     axpy        ( const value_t alpha, const TVector< value_t > * x );

    //! set this ≔ this + x (thread safe, is(x) ⊆ is(this))
    virtual void     add_sub_mt  ( const TScalarVector< value_t > &  x );

    //////////////////////////////////////////////////
    //
    // misc. methods
    //

    //! return size in bytes used by this object
    virtual size_t byte_size      () const;

    //! return size of (floating point) data in bytes handled by this object
    virtual size_t data_byte_size () const { return sizeof(value_t) * size(); }

    //
    // virtual constructor
    //

    //! return copy of vector
    virtual auto  copy   () const -> std::unique_ptr< TVector< value_t > > { return std::make_unique< TScalarVector< value_t > >( *this ); }
    
    //! return object of same class
    virtual auto  create () const -> std::unique_ptr< TVector< value_t > > { return std::make_unique< TScalarVector< value_t > >(); }
    
    //
    // restriction
    //

    //! return vector restricted to real part of coefficients
    virtual auto  restrict_re  () const -> std::unique_ptr< TVector< real_t > >;

    //! return vector restricted to imaginary part of coefficients
    virtual auto  restrict_im  () const -> std::unique_ptr< TVector< real_t > >;
    
    //! return vector restricted to absolute value of coefficients
    virtual auto  restrict_abs () const -> std::unique_ptr< TVector< real_t > >;
    
    //
    // stream output
    //

    //! print vector information
    virtual void print ( const uint ofs = 0 ) const;

    //
    // serialisation
    //

    //! read vector from stream
    virtual void read  ( TByteStream & s );

    //! write vector to stream
    virtual void write ( TByteStream & s ) const;

    //! returns size of object in bytestream
    virtual size_t bs_size () const;

    //
    // parallel methods
    //

    //! pointwise summation between all vectors in \a ps
    virtual void sum ( const TProcSet & ps );
    using TVector< value_t >::sum;

    HPRO_RTTI_DERIVED( TScalarVector, TVector< value_t > )
};

//////////////////////////////////////////////////////////
//
// template wrappers
//

//
// blas_rvec/_cvec const
//
template <typename value_t>       BLAS::Vector< value_t > &  blas_vec  (       TScalarVector< value_t > *  v ) { return v->blas_vec(); }
template <typename value_t> const BLAS::Vector< value_t > &  blas_vec  ( const TScalarVector< value_t > *  v ) { return v->blas_vec(); }

//!
//! functional version of TScalarVector::sub_vector
//!
template <typename value_t>
TScalarVector< value_t >
sub_vector ( TScalarVector< value_t > *  v,
             const TIndexSet &           is )
{
    return v->sub_vector( is );
}
             
template <typename value_t>
const TScalarVector< value_t >
sub_vector ( const TScalarVector< value_t > *  v,
             const TIndexSet &                 is )
{
    return v->sub_vector( is );
}

//
// type checks
//

template <typename value_t> bool is_scalar ( const TVector< value_t > &  v ) noexcept { return IS_TYPE( & v, TScalarVector ); }
template <typename value_t> bool is_scalar ( const TVector< value_t > *  v ) noexcept { return ( v != nullptr ) && IS_TYPE( v, TScalarVector ); }

}// namespace Hpro

#endif  // __HPRO_TSCALARVECTOR_HH
