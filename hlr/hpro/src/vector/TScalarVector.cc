//
// Project     : HLIBpro
// File        : TScalarVector.cc
// Description : class for a vector of scalar type
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <iostream>

#include "hpro/blas/Algebra.hh"
#include "hpro/parallel/NET.hh"

#include "TRNG.hh"

#include "hpro/vector/TScalarVector.hh"

namespace Hpro
{

// namespace abbr.
namespace B = BLAS;

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// TScalarVector
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////
//
// misc. methods
//

//
// set size of vector
//
template < typename value_t >
void
TScalarVector< value_t >::set_size ( const size_t  n )
{
    if ( n != _size )
    {
        if ( n == 0 ) _vec = B::Vector< value_t >();
        else          _vec = B::Vector< value_t >( n );

        _size = n;

        init_chunk_mutices();        
    }// if
}

//
// return size in bytes used by this object
//
template < typename value_t >
size_t
TScalarVector< value_t >::byte_size () const
{
    return TVector< value_t >::byte_size() + sizeof(size_t) + sizeof(B::Vector<value_t>) + sizeof(value_t) * size();
}

//////////////////////////////////////////////////
//
// BLAS-routines (real valued)
//

//
// fill with constant
//
template < typename value_t >
void
TScalarVector< value_t >::fill ( const value_t  f )
{
    B::fill( f, _vec );
}

//
// fill with random numbers (but prevent zero entries)
//
template <>
void
TScalarVector< float >::fill_rand ( const uint seed )
{
    TRNG  rng( seed );

    for ( idx_t  i = 0; i < idx_t(_size); i++ )
        _vec(i) = float( rng( 1000.0 ) - 500.0 );
}
    
template <>
void
TScalarVector< double >::fill_rand ( const uint seed )
{
    TRNG  rng( seed );

    for ( idx_t  i = 0; i < idx_t(_size); i++ )
        _vec(i) = double( rng( 1000.0 ) - 500.0 );
}

template <>
void
TScalarVector< std::complex< float > >::fill_rand ( const uint seed )
{
    TRNG  rng( seed );

    for ( idx_t  i = 0; i < idx_t(_size); i++ )
        _vec(i) = std::complex< float >( rng( 1000.0 ) - 500.0, rng( 1000.0 ) - 500.0 );
}
    
template <>
void
TScalarVector< std::complex< double > >::fill_rand ( const uint seed )
{
    TRNG  rng( seed );

    for ( idx_t  i = 0; i < idx_t(_size); i++ )
        _vec(i) = std::complex< double >( rng( 1000.0 ) - 500.0, rng( 1000.0 ) - 500.0 );
}

//
// this = a * this
//
template < typename value_t >
void
TScalarVector< value_t >::scale ( const value_t  f )
{
    if ( f == value_t(0) ) fill( value_t(0) );
    else                   B::scale( f, _vec );
}

//
// this = x
//
template < typename value_t >
void
TScalarVector< value_t >::assign ( const value_t  f,
                                   const TVector< value_t >* x )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) assign", "vector is nullptr" );
    
    if ( IS_TYPE( x, TScalarVector ) )
    {
        const auto  sx = cptrcast( x, TScalarVector< value_t >);
        
        set_is( sx->is() );

        B::copy( sx->_vec, _vec );

        if ( f != value_t(1) )
            scale( f );
    }// if
    else
        HERROR( ERR_VEC_TYPE, "(TScalarVector) assign", x->typestr() );
}

//
// return euclidean norm
//
template < typename value_t >
real_type_t< value_t >
TScalarVector< value_t >::norm2 () const
{
    return B::norm2( _vec );
}

//
// return infimum norm
//
template < typename value_t >
real_type_t< value_t >
TScalarVector< value_t >::norm_inf () const
{
    return Math::abs( _vec( B::max_idx( _vec ) ) );
}   

//
// this = this + a * x
//
template < typename value_t >
void
TScalarVector< value_t >::axpy ( const value_t               f,
                                 const TVector< value_t > *  x )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) axpy", "vector is nullptr" );
    
    if ( ! IS_TYPE( x, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TScalarVector) axpy", x->typestr() );

    const auto  sx = cptrcast( x, TScalarVector< value_t >);

    if ( size() != sx->size() )
        HERROR( ERR_VEC_SIZE, "(TScalarVector) axpy", "x has wrong size" );
    
    B::add( f, sx->_vec, _vec );
}

//
// this = this + x (multi-thread safe version, x may be sub vector)
//
template < typename value_t >
void
TScalarVector< value_t >::add_sub_mt ( const TScalarVector< value_t > &  x )
{
    HERROR( ERR_NOT_IMPL, "", "" );
    
    if ( ! x.is().is_subset_of( this->is() ) )
        HERROR( ERR_INDEXSET, "(TScalarVector) add_sub_mt", "given vector is not a sub-vector" );
    
    const TIndexSet  is_x( x.is() );
    const TIndexSet  is_y( this->is() );
    const idx_t      ofs_loc_glo = is_x.first() - is_y.first();
    idx_t            start_idx   = ofs_loc_glo;
    idx_t            chunk       = start_idx / SCALAR_CHUNK_SIZE;
    const idx_t      last_idx    = is_x.last() - is_y.first();
    idx_t            end_idx     = std::min< idx_t >( (chunk+1) * SCALAR_CHUNK_SIZE - 1, last_idx );
    
    while ( start_idx <= end_idx )
    {
        const B::Range        is_chunk( start_idx, end_idx );
        B::Vector< value_t >  y_chunk( blas_vec(), is_chunk );
        B::Vector< value_t >  x_chunk( x.blas_vec(), is_chunk - ofs_loc_glo );

        {
            // lock_t  lock( _mutices[ chunk ] );
                
            B::add( value_t(1), x_chunk, y_chunk );
        }

        ++chunk;
        start_idx = end_idx + 1;
        end_idx   = std::min< idx_t >( end_idx + SCALAR_CHUNK_SIZE, last_idx );
    }// while
}

//
// return inner product ( this^H * x )
//
template < typename value_t >
value_t
TScalarVector< value_t >::dot ( const TVector< value_t >* x ) const
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) dot", "vector is nullptr" );
    
    if ( IS_TYPE( x, TScalarVector ) )
    {
        return B::dot( _vec, cptrcast( x, TScalarVector< value_t >)->_vec );
    }// if
    else
        HERROR( ERR_VEC_TYPE, "(TScalarVector) dot", x->typestr() );
    
    return value_t(0);
}

//
// return inner product ( this^T * x )
//
template < typename value_t >
value_t
TScalarVector< value_t >::dotu ( const TVector< value_t >* x ) const
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) dotu", "vector is nullptr" );
    
    if ( IS_TYPE( x, TScalarVector ) )
    {
        return B::dotu( _vec, cptrcast( x, TScalarVector< value_t >)->_vec );
    }// if
    else
        HERROR( ERR_VEC_TYPE, "(TScalarVector) dotu", x->typestr() );
    
    return value_t(0);
}

//////////////////////////////////////////////////
//
// misc. methods
//

//
// create vector restricted to real/imaginary part of coefficients
//
template < typename value_t >
auto
TScalarVector< value_t >::restrict_re () const -> std::unique_ptr< TVector< real_type_t< value_t > > >
{
    auto  v = std::make_unique< TScalarVector< real_type_t< value_t > > >( this->is() );

    for ( idx_t  i = 0; i < idx_t(_size); i++ )
        v->set_entry( i, std::real( _vec(i) ) );

    return v;
}

template < typename value_t >
auto
TScalarVector< value_t >::restrict_im () const -> std::unique_ptr< TVector< real_type_t< value_t > > >
{
    auto  v = std::make_unique< TScalarVector< real_type_t< value_t > > >( this->is() );

    for ( idx_t  i = 0; i < idx_t(_size); i++ )
        v->set_entry( i, std::imag( _vec(i) ) );

    return v;
}

template < typename value_t >
auto
TScalarVector< value_t >::restrict_abs () const -> std::unique_ptr< TVector< real_type_t< value_t > > >
{
    auto  v = std::make_unique< TScalarVector< real_type_t< value_t > > >( this->is() );

    for ( idx_t  i = 0; i < idx_t(_size); i++ )
        v->set_entry( i, std::abs( _vec(i) ) );

    return v;
}

//
// copy from C array \a v
//
template < typename value_t >
void 
TScalarVector< value_t >::copy_from ( const value_t *  v )
{
    if ( v == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) copy_from", "array is nullptr" );

    for ( idx_t  i = 0; i < idx_t(_vec.length()); ++i )
        _vec(i) = v[i];
}

//
// copy to C array \a v
//
template < typename value_t >
void 
TScalarVector< value_t >::copy_to ( value_t *  v )
{
    if ( v == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) copy_to", "array is nullptr" );

    for ( idx_t  i = 0; i < idx_t(_vec.length()); ++i )
        v[i] = _vec(i);
}

//
// permute entries according to \a perm
//
template < typename value_t >
void
TScalarVector< value_t >::permute ( const TPermutation &  perm )
{
    perm.permute( this );
}

//
// stream output
//
template < typename value_t >
void
TScalarVector< value_t >::print ( const uint offset ) const
{
    for ( uint i = 0; i < offset; i++ )
        std::cout << ' ';

    std::cout << this->typestr() << " ( ";

    for ( idx_t  i = 0; i < idx_t(_size); i++ )
    {
        if ( i > 0 ) std::cout << ", ";
        std::cout << this->entry(i);
    }// for

    std::cout << " )" << std::endl;
}

//
// serialisation
//
template < typename value_t >
void
TScalarVector< value_t >::read  ( TByteStream & s )
{
    TVector< value_t >::read( s );

    size_t  vec_size;
    
    s.get( vec_size );
    set_size( vec_size );

    s.get( blas_vec().data(), sizeof(value_t) * _size );
}

template < typename value_t >
void
TScalarVector< value_t >::write ( TByteStream & s ) const
{
    TVector< value_t >::write( s );
    
    s.put( _size );
    s.put( blas_vec().data(), sizeof(value_t) * _size );
}

//
// returns size of object in bytestream
//
template < typename value_t >
size_t
TScalarVector< value_t >::bs_size () const
{
    return TVector< value_t >::bs_size() + sizeof(_size) + sizeof(value_t) * _size;
}

//
// pointwise summation between all vectors in \a ps
//
template < typename value_t >
void
TScalarVector< value_t >::sum ( const TProcSet &  ps )
{
    B::Vector< value_t >  t( _vec.length() );

    NET::reduce_all( ps, _vec.data(), t.data(), _vec.length(), NET::OP_SUM );
    B::copy( t, _vec );
}

template class TScalarVector< float >;
template class TScalarVector< double >;
template class TScalarVector< std::complex< float > >;
template class TScalarVector< std::complex< double > >;

}// namespace Hpro
