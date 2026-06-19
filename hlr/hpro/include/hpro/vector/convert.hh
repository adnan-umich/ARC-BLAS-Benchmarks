#ifndef __HPRO_VEC_CONV_HH
#define __HPRO_VEC_CONV_HH
//
// Project     : HLIBpro
// File        : vec_conv.hh
// Description : vector conversion functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/vector/TVector.hh"
#include "hpro/vector/TScalarVector.hh"

namespace Hpro
{

//!
//! convert between different value types
//!
template < typename dest_value_t,
           typename src_value_t = dest_value_t >
std::unique_ptr< TVector< dest_value_t > >
convert ( const TVector< src_value_t > *  v );
       
template < typename dest_value_t,
           typename src_value_t = dest_value_t >
std::unique_ptr< TVector< dest_value_t > >
convert ( const TVector< src_value_t > &  v )
{
    return convert< dest_value_t, src_value_t >( &v );
}

//!
//! convert between different value types
//!
template < typename dest_value_t,
           typename src_value_t = dest_value_t >
void
convert_to ( const TVector< src_value_t > *  src,
             TVector< dest_value_t > *       dest )
{
    if (( src == nullptr ) || ( dest == nullptr ))
        HERROR( ERR_ARG, "convert_to", "vector is null" );

    if constexpr ( std::is_same< dest_value_t, src_value_t >::value )
    {
        dest->assign( src_value_t(1), src );
    }// if
    else
    {
        if ( is_scalar( src ) && is_scalar( dest ) )
        {
            auto  ssrc  = cptrcast( src,  TScalarVector< src_value_t > );
            auto  vsrc  = ssrc->blas_vec();
            auto  sdest =  ptrcast( dest, TScalarVector< dest_value_t > );
            auto  vdest = sdest->blas_vec();
            
            for ( size_t  i = 0; i < vsrc.length(); ++i )
                vdest(i) = dest_value_t( vsrc(i) );
        }// if
        else
            HERROR( ERR_VEC_TYPE, "convert_to", "vector type unsupported : " + src->typestr() + "/" + dest->typestr() );
    }// else
}
       
template < typename dest_value_t,
           typename src_value_t = dest_value_t >
void
convert_to ( const TVector< src_value_t > &  src,
             TVector< dest_value_t > &       dest )
{
    convert_to< dest_value_t, src_value_t >( &src, &dest );
}

//!
//! convert \a v into scalar vector
//!
template < typename value_t >
std::unique_ptr< TVector< value_t > >
to_scalar  ( const TVector< value_t > *  v );

template < typename value_t >
std::unique_ptr< TVector< value_t > >
to_scalar  ( const TVector< value_t > &  v )
{
    return to_scalar( &v );
}

//!
//! convert \a v into blocked vector as defined by
//! to cluster tree \a ct
//!
template < typename value_t >
std::unique_ptr< TVector< value_t > >
to_blocked ( const TVector< value_t > *  v,
             const TCluster *            ct );

template < typename value_t >
std::unique_ptr< TVector< value_t > >
to_blocked ( const TVector< value_t > &  v,
             const TCluster &            ct )
{
    return to_blocked( &v, &ct );
}

}// namespace Hpro

#endif // __HPRO_VEC_CONV_HH
