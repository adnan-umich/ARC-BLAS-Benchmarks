//
// Project     : HLIBpro
// File        : vec_conv.cc
// Description : vector conversion functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "list.hh"

#include "hpro/vector/TScalarVector.hh"
#include "hpro/vector/TBlockVector.hh"

#include "hpro/vector/convert.hh"

namespace Hpro
{

// namespace abbr.
namespace B = BLAS;

namespace
{

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
// local functions
//
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
//
// conversion to scalar
//
///////////////////////////////////////////////////////////////////////

template < typename value_t >
TVector< value_t > *
block_to_scalar ( const TBlockVector< value_t > *  v )
{
    auto  s      = std::make_unique< TScalarVector< value_t > >( v->is() );
    auto  subvec = std::list< const TVector< value_t > * >();

    subvec.push_back( v );
    
    while ( ! subvec.empty() )
    {
        const auto  w = behead( subvec );

        if ( IS_TYPE( w, TBlockVector ) )
        {
            const auto  bw = cptrcast( w, TBlockVector< value_t > );

            for ( uint  i = 0; i < bw->n_blocks(); ++i )
                if ( bw->block(i) != nullptr )
                    subvec.push_back( bw->block(i) );
        }// if
        else if ( IS_TYPE( w, TScalarVector ) )
        {
            const auto  sw = cptrcast( w, TScalarVector< value_t > );

            B::Vector< value_t >  s_part( s->blas_vec(), sw->is() - s->ofs() );

            B::copy( sw->blas_vec(), s_part );
        }// if
        else
            HERROR( ERR_VEC_TYPE, "block_to_scalar", w->typestr() );
    }// while

    return s.release();
}

///////////////////////////////////////////////////////////////////////
//
// conversion to blocked
//
///////////////////////////////////////////////////////////////////////

template < typename value_t >
TVector< value_t > *
scalar_to_blocked ( const TScalarVector< value_t > *  v,
                    const TCluster *                  ct )
{
    std::unique_ptr< TVector< value_t > >  w;
        
    if ( ct->is_leaf() )
    {
        auto  v_ct = B::Vector< value_t >( v->blas_vec(), *ct - v->ofs(), copy_value );

        w = std::make_unique< TScalarVector< value_t > >( *ct, v_ct );
    }// if
    else
    {
        std::vector< TVector< value_t > * >  subblocks( ct->nsons() );
        
        for ( uint i = 0; i < ct->nsons(); ++i )
        {
            if ( ct->son(i) != nullptr )
                subblocks[i] = scalar_to_blocked( v, ct->son(i) );
        }// for

        w = std::make_unique< TBlockVector< value_t > >( *ct, subblocks );
    }// else

    return w.release();
}

}// namespace anonymous

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
// global functions
//
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template < typename dest_value_t,
           typename src_value_t >
std::unique_ptr< TVector< dest_value_t > >
convert ( const TVector< src_value_t > *  v )
{
    if ( v == nullptr )
        HERROR( ERR_ARG, "convert", "vector is null" );

    if constexpr ( std::is_same< dest_value_t, src_value_t >::value )
    {
        return v->copy();
    }// if
    else
    {
        if ( is_scalar( v ) )
        {
            auto  sv = cptrcast( v, TScalarVector< src_value_t > );
            auto  t  = BLAS::convert< dest_value_t >( sv->blas_vec() );
            
            return std::make_unique< TScalarVector< dest_value_t > >( v->is(), std::move( t ) );
        }// if
        else
            HERROR( ERR_VEC_TYPE, "convert", "vector type unsupported : " + v->typestr() );
    }// else
}

//
// convert \a v into scalar vector
//
template < typename value_t >
std::unique_ptr< TVector< value_t > >
to_scalar ( const TVector< value_t > *  v )
{
    if ( v == nullptr )
        return nullptr;

    std::unique_ptr< TVector< value_t > >  res;
    
    if ( IS_TYPE( v, TBlockVector ) )
    {
        res = std::unique_ptr< TVector< value_t > >( block_to_scalar( cptrcast( v, TBlockVector< value_t > ) ) );
    }// if
    else
        HERROR( ERR_VEC_TYPE, "to_scalar", v->typestr() );
        
    return res;
}

//
// convert \a v into blocked vector as defined by
// to cluster tree \a ct
//
template < typename value_t >
std::unique_ptr< TVector< value_t > >
to_blocked ( const TVector< value_t > *   v,
             const TCluster *             ct )
{
    if ( v == nullptr )
        return nullptr;

    if ( ct == nullptr )
        HERROR( ERR_ARG, "to_blocked", "cluster tree is nullptr" );

    std::unique_ptr< TVector< value_t > >  res;
    
    if ( IS_TYPE( v, TScalarVector ) )
        res = std::unique_ptr< TVector< value_t > >( scalar_to_blocked( cptrcast( v, TScalarVector< value_t > ), ct ) );
    else
        HERROR( ERR_VEC_TYPE, "to_blocked", v->typestr() );
        
    return res;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// explicit instantiation
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

#define INST_CONVERT( type1, type2 )                                    \
    template std::unique_ptr< TVector< type1 > > convert< type1, type2 > ( const TVector< type2 > * );
//    template void convert_to< type1, type2 > ( const TVector< type2 > *, TVector< type1 > * );

INST_CONVERT( float, float )
INST_CONVERT( double, float )
INST_CONVERT( std::complex< float >, float )
INST_CONVERT( std::complex< double >, float )

INST_CONVERT( float, double )
INST_CONVERT( double, double )
INST_CONVERT( std::complex< float >, double )
INST_CONVERT( std::complex< double >, double )

INST_CONVERT( std::complex< float >, std::complex< float > )
INST_CONVERT( std::complex< float >, std::complex< double > )

INST_CONVERT( std::complex< double >, std::complex< float > )
INST_CONVERT( std::complex< double >, std::complex< double > )

#define  INST_ALL( type )                                     \
    template std::unique_ptr< TVector< type > > to_scalar< type > ( const TVector< type > * ); \
    template std::unique_ptr< TVector< type > > to_blocked< type > ( const TVector< type > *, const TCluster * ); \

INST_ALL( float )
INST_ALL( double )
INST_ALL( std::complex< float > )
INST_ALL( std::complex< double > )

}// namespace Hpro
