//
// Project     : HLIBpro
// File        : TBlockVector.cc
// Description : class for a vector build out of other vectors
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <iostream>

#include "hpro/vector/TScalarVector.hh"

#include "hpro/vector/TBlockVector.hh"

namespace Hpro
{

// namespace abbr.
namespace B = BLAS;

///////////////////////////////////////////////
//
// constructor and destructor
//
template < typename value_t >
TBlockVector< value_t >::~TBlockVector ()
{
    for ( uint i = 0; i < n_blocks(); i++ )
        delete block(i);
}

///////////////////////////////////////////////
//
// access internal data
//

//
// return size of vector
//
template < typename value_t >
size_t
TBlockVector< value_t >::size () const
{
    size_t  s = 0;
    
    for ( uint i = 0; i < n_blocks(); i++ )
    {
        if ( block(i) != nullptr )
            s += block(i)->size();
    }// for

    return s;
}

//
// set single vector block
//
template < typename value_t >
void
TBlockVector< value_t >::set_block ( const uint i, TVector< value_t > * v )
{
    if ( i >= n_blocks() )
        HERROR( ERR_VEC_STRUCT, "(TBlockVector) set_block", "index too high" );

    if ( block(i) != nullptr )
        delete block(i);

    _blocks[i] = v;
}
    
//
// setup blockstructure of vector
//
template < typename value_t >
void
TBlockVector< value_t >::set_block_struct ( const uint n )
{
    uint  old_n = n_blocks();

    for ( uint i = n; i < old_n; i++ )
        if ( _blocks[i] != nullptr )
            delete _blocks[i];
    
    _blocks.resize( n, nullptr );
}

////////////////////////////////////////////////
//
// access entries
//

//
// return i'th entry
//
template < typename value_t >
value_t
TBlockVector< value_t >::entry  ( const idx_t  row ) const
{
    const idx_t  vofs = this->ofs();
    
    for ( uint i = 0; i < n_blocks(); i++ )
    {
        const TVector< value_t > * v_i = block(i);
        
        if ( v_i == nullptr )
            continue;
        
        if (( v_i->ofs() - vofs <= row ) && ( row < v_i->ofs() - vofs + idx_t( v_i->size() ) ))
            return v_i->entry( row - v_i->ofs() + vofs );
    }// for

    return value_t(0);
}

//
// set i'th entry
//
template < typename value_t >
void
TBlockVector< value_t >::set_entry  ( const idx_t    row,
                           const value_t  f )
{
    const idx_t  vofs = this->ofs();
    
    for ( uint i = 0; i < n_blocks(); i++ )
    {
        TVector< value_t > * v_i = block(i);
        
        if ( v_i == nullptr )
            continue;
        
        if (( v_i->ofs() - vofs <= row ) && ( row < v_i->ofs() - vofs + idx_t( v_i->size() ) ))
        {
            v_i->set_entry( row - v_i->ofs() + vofs, f );
            return;
        }// if
    }// for
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
TBlockVector< value_t >::fill ( const value_t f )
{
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->fill( f );
}
    
//
// fill with random numbers
//
template < typename value_t >
void
TBlockVector< value_t >::fill_rand ( const uint seed )
{
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->fill_rand( seed );
}

//
// this = this + a * x
//
template < typename value_t >
void
TBlockVector< value_t >::axpy ( const value_t               f,
                                const TVector< value_t > *  x )
{
    if ( f == value_t(0) )
        return;
    
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TBlockVector) axpy", "given vector is nullptr" );

    if ( IS_TYPE( x, TBlockVector ) )
    {
        const auto  b = cptrcast( x, TBlockVector< value_t > );

        if ( b->n_blocks() != n_blocks() )
            HERROR( ERR_VEC_STRUCT, "(TBlockVector) axpy",
                    "given vector has different number of blocks" );
        
        for ( uint i = 0; i < n_blocks(); i++ )
            if ( block(i) != nullptr )
                block(i)->axpy( f, b->block(i) );
    }// if
    else if ( IS_TYPE( x, TScalarVector ) )
    {
        const auto  s = cptrcast( x, TScalarVector< value_t > );

        for ( uint i = 0; i < n_blocks(); i++ )
        {
            auto  bi = block(i);
                    
            if ( bi != nullptr )
            {
                B::Vector< value_t >      t( s->blas_vec(), bi->is() - s->ofs() );
                TScalarVector< value_t >  st( bi->is(), t );
                
                bi->axpy( f, & st );
            }// if
        }// for
    }// if
    else
    {
        HERROR( ERR_VEC_TYPE, "(TBlockVector) axpy", x->typestr() );
    }// else

}

//
// this = a * this
//
template < typename value_t >
void
TBlockVector< value_t >::scale ( const value_t  f )
{
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->scale( f );
}

//
// this = f * x
//
template < typename value_t >
void
TBlockVector< value_t >::assign ( const value_t               f,
                                  const TVector< value_t > *  x )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TBlockVector) assign", "given vector is nullptr" );
    
    if ( ! IS_TYPE( x, TBlockVector ) )
        HERROR( ERR_VEC_TYPE, "(TBlockVector) assign",
               "expected TBlockVector, got " + x->typestr() );

    const auto  b = cptrcast( x, TBlockVector< value_t > );

    this->set_ofs( b->ofs() );
    
    if ( b->n_blocks() != n_blocks() )
        HERROR( ERR_VEC_STRUCT, "(TBlockVector) assign",
               "given vector has different number of blocks" );
    
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->assign( f, b->block(i) );
}

//
// return infimum norm
//
template < typename value_t >
real_type_t< value_t >
TBlockVector< value_t >::norm_inf () const
{
    real_type_t< value_t >  f = 0.0;

    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            f = std::max( f, block(i)->norm_inf() );

    return f;
}

//
// conjugate entries
//
template < typename value_t >
void
TBlockVector< value_t >::conjugate ()
{
    for ( uint  i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->conjugate();
}

//
// return inner product ( this^H * x )
//
template < typename value_t >
value_t
TBlockVector< value_t >::dot ( const TVector< value_t > * x ) const
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TBlockVector) dot", "given vector is nullptr" );
    
    if ( ! IS_TYPE( x, TBlockVector ) )
        HERROR( ERR_VEC_TYPE, "(TBlockVector) dot",
               "expected TBlockVector, got " + x->typestr() );

    const auto  b = cptrcast( x, TBlockVector< value_t > );

    if ( b->n_blocks() != n_blocks() )
        HERROR( ERR_VEC_STRUCT, "(TBlockVector) dot", "given vector has different number of blocks" );

    value_t  f = 0.0;
    
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            f += block(i)->dot( b->block(i) );

    return f;
}

//
// return inner product ( this^T * x )
//
template < typename value_t >
value_t
TBlockVector< value_t >::dotu ( const TVector< value_t > * x ) const
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TBlockVector) dotu", "given vector is nullptr" );
    
    if ( ! IS_TYPE( x, TBlockVector ) )
        HERROR( ERR_VEC_TYPE, "(TBlockVector) dotu",
               "expected TBlockVector, got " + x->typestr() );

    const auto  b = cptrcast( x, TBlockVector< value_t > );

    if ( b->n_blocks() != n_blocks() )
        HERROR( ERR_VEC_STRUCT, "(TBlockVector) dotu", "given vector has different number of blocks" );

    value_t  f = 0.0;
    
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            f += block(i)->dotu( b->block(i) );

    return f;
}

///////////////////////////////////////////////
//
// misc.
//

//
// return size in bytes used by this object
//
template < typename value_t >
size_t
TBlockVector< value_t >::byte_size () const
{
    size_t  s = TVector< value_t >::byte_size() + sizeof(_blocks);

    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            s += block(i)->byte_size();

    return s;
}

//
// return copy of matrix
//
template < typename value_t >
std::unique_ptr< TVector< value_t > >
TBlockVector< value_t >::copy  () const
{
    auto  v = create();
    auto  x = ptrcast( v.get(), TBlockVector< value_t > );

    x->set_ofs( this->ofs() );
    x->set_block_struct( n_blocks() );

    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            x->set_block( i, block(i)->copy().release() );

    return v;
}
    
//
// stream output
//
template < typename value_t >
void
TBlockVector< value_t >::print ( const uint offset ) const
{
    for ( uint i = 0; i < offset; i++ )
        std::cout << ' ';

    std::cout << "TBlockVector [ ";
    
    for ( uint i = 0; i < n_blocks(); i++ )
    {
        if ( block(i) != nullptr )
            block(i)->print( offset+2 );
        else
        {
            for ( uint j = 0; j < offset+2; j++ )
                std::cout << ' ';
            
            std::cout << "nullptr";
        }// else

        std::cout << std::endl;
    }// if

    std::cout << " ]" << std::endl;
}

template class TBlockVector< float >;
template class TBlockVector< double >;
template class TBlockVector< std::complex< float > >;
template class TBlockVector< std::complex< double > >;

}// namespace Hpro
