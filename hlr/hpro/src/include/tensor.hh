#ifndef __HPRO_TENSOR_HH
#define __HPRO_TENSOR_HH
//
// Project     : HLIBpro
// File        : matrix.hh
// Description : tensor containers based on std::vector
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

namespace Hpro
{

template < typename T >
struct tensor2
{
public:
    using  value_type      = T;
    using  pointer         = typename std::vector< value_type >::pointer;
    using  const_pointer   = typename std::vector< value_type >::const_pointer;
    using  reference       = typename std::vector< value_type >::reference;
    using  const_reference = typename std::vector< value_type >::const_reference;
    using  size_type       = typename std::vector< value_type >::size_type;
    using  difference_type = typename std::vector< value_type >::difference_type;

private:
    std::vector< value_type >  _data;
    const size_type            _dim0;

public:
    tensor2 ( const size_type  adim0,
              const size_type  adim1 )
            : _data( adim0 * adim1 )
            , _dim0( adim0 )
    {}

    tensor2 ( const size_type  adim0,
              const size_type  adim1,
              const value_type       adefault )
            : _data( adim0 * adim1, adefault )
            , _dim0( adim0 )
    {}

    tensor2 ( const tensor2 &  t )
            : _data( t._data )
            , _dim0( t._dim0 )
    {}

    tensor2 ( tensor2 &&  t )
            : _data( std::move( t._data ) )
            , _dim0( t._dim0 )
    {}

    tensor2 &
    operator = ( tensor2 &&  t )
    {
        _data = std::move( t._data );
        _dim0 = t._dim0;

        return *this;
    }

    size_type         dim0     () const { return _dim0; }
    size_type         dim1     () const { return _data.size() / _dim0; }
    
    const value_type  operator ()  ( const size_type  i, const size_type  j ) const { return _data[ j*_dim0 + i ]; }
    value_type &      operator ()  ( const size_type  i, const size_type  j )       { return _data[ j*_dim0 + i ]; }
};
         
template < typename T >
struct tensor3
{
public:
    using  value_type      = T;
    using  pointer         = typename std::vector< value_type >::pointer;
    using  const_pointer   = typename std::vector< value_type >::const_pointer;
    using  reference       = typename std::vector< value_type >::reference;
    using  const_reference = typename std::vector< value_type >::const_reference;
    using  size_type       = typename std::vector< value_type >::size_type;
    using  difference_type = typename std::vector< value_type >::difference_type;

private:
    std::vector< value_type >  data;
    const size_type            dim0;
    const size_type            dim1;

public:
    tensor3 ( const size_type  adim0,
              const size_type  adim1,
              const size_type  adim2 )
            : data( adim0 * adim1 * adim2 )
            , dim0( adim0 )
            , dim1( adim1 )
    {}

    tensor3 ( const size_type   adim0,
              const size_type   adim1,
              const size_type   adim2,
              const value_type  adefault )
            : data( adim0 * adim1 * adim2, adefault )
            , dim0( adim0 )
            , dim1( adim1 )
    {}

    tensor3 ( const tensor3 &  t )
            : data( t.data )
            , dim0( t.dim0 )
            , dim1( t.dim1 )
    {}

    tensor3 ( tensor3 &&  t )
            : data( std::move( t.data ) )
            , dim0( t.dim0 )
            , dim1( t.dim1 )
    {}

    tensor3 &
    operator = ( tensor3 &&  t )
    {
        data = std::move( t.data );
        dim0 = t.dim0;
        dim1 = t.dim1;

        return *this;
    }
    
    const value_type  operator ()  ( const size_type  i, const size_type  j, const size_type  k ) const { return data[ (k*dim1 + j)*dim0 + i ]; }
    value_type &      operator ()  ( const size_type  i, const size_type  j, const size_type  k )       { return data[ (k*dim1 + j)*dim0 + i ]; }
};

}// namespace Hpro

#endif   // __HPRO_TENSOR_HH
