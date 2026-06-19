#ifndef __HPRO_BASE_ALLOC_HH
#define __HPRO_BASE_ALLOC_HH
//
// Project     : HLIBpro
// File        : base/alloc.hh
// Description : module containing memory allocation routines
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2025. All Rights Reserved.
//

#include <vector>

#include <tbb/tbb_allocator.h>

#include <memory>

#include <hpro/config.h>

namespace Hpro
{

//
// basic memory allocator implementation used within Hpro
// (based on malloc/free by default)
//
template < typename T >
class hpro_allocator
{
public:
    using  value_type      = T;
    using  pointer         = value_type *;
    using  const_pointer   = const value_type *;
    using  reference       = value_type &;
    using  const_reference = const value_type &;
    using  size_type       = size_t;
    using  difference_type = ptrdiff_t;

    template < class U >
    struct rebind
    {
        using  other = hpro_allocator< U >;
    };

    hpro_allocator ()                           throw() {}
    hpro_allocator ( const hpro_allocator & )   throw() {}
    template < typename U >
    hpro_allocator (const hpro_allocator<U> & ) throw() {}

    pointer       address ( reference        x ) const { return & x; }
    const_pointer address ( const_reference  x ) const { return & x; }

    //! allocate memory for <n> objects.
    pointer
    allocate ( size_type       n,
               const void * /* hint */ = 0 ) // ignored
    {
        return static_cast< pointer >( ::malloc( n * sizeof(value_type) ) );
    }

    //! release allocated block of memory
    void
    deallocate ( pointer    p,
                 size_type  /* size */ = 0 )
    {
        ::free( p );
    }

    //! Largest value for which method allocate might succeed.
    size_type
    max_size () const throw()
    {
        size_type  absolutemax = static_cast< size_type >(-1) / sizeof (value_type);
        
        return (absolutemax > 0 ? absolutemax : 1);
    }
    
    void construct ( pointer             p,
                     const value_type &  value )
    {
        ::new((void*)(p)) value_type( value );
    }

    void destroy   ( pointer  p )
    {
        p->~value_type();
    }
};

//
// same as std::allocator<void>
//
template<>
class hpro_allocator<void>
{
public:
    using  pointer = void *;
    using  const_pointer = const void *;
    using  value_type    = void;
    
    template<class U>
    struct rebind
    {
        using  other = hpro_allocator<U>;
    };
};

template<typename T, typename U> inline bool operator ==  ( const hpro_allocator<T> &, const hpro_allocator<U> & ) { return true;  }
template<typename T, typename U> inline bool operator !=  ( const hpro_allocator<T> &, const hpro_allocator<U> & ) { return false; }

////////////////////////////////////////////////////////////////////////////////
//
// define internal used memory allocator
//

template < typename T > using allocator = HPRO_ALLOCATOR ;

////////////////////////////////////////////////////////////////////////////////
//
// internal versions of std::XXX using "allocator"
//

template < typename T > using vector = std::vector< T, allocator< T > >;

}// namespace Hpro

#endif // __HPRO_BASE_ALLOC_HH
