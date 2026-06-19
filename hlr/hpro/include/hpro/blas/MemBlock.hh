#ifndef __HPRO_BLAS_MEMBLOCK_HH
#define __HPRO_BLAS_MEMBLOCK_HH
//
// Project     : HLIBpro
// File        : MemBlock.hh
// Description : class for a reference countable memory block
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2025. All Rights Reserved.
//

#include <cstdlib>

#include <atomic>

#include <hpro/base/error.hh>
#include <hpro/base/alloc.hh>
#include <hpro/base/System.hh>

namespace Hpro
{

// enable (=1) for debugging
#define HPRO_MEMBLOCK_REF_COUNT  0

// enable memory debugging (see also Mem::malloc/free)
#define HPRO_DEBUG_MALLOC        0

//!
//! indicates copy policy
//!
enum copy_policy_t
{
    copy_reference, // copy pointer to data only
    copy_value      // copy actual data
};

namespace BLAS
{

//!
//! \ingroup  BLAS_Module
//! \class    MemBlock
//! \brief    Defines a reference countable memory block.
//!
template < class T_value >
class MemBlock
{
public:
    //! internal value type
    using  value_t = T_value;

    #if HPRO_MEMBLOCK_REF_COUNT == 1
    //! for debugging: reference counter type
    using  count_t = std::atomic< size_t >;
    #endif

protected:
    //! @cond

    //! pointer of data
    value_t *             _data;

    //! indicates, this object is owner of memory block
    bool                  _is_owner;

    #if HPRO_MEMBLOCK_REF_COUNT == 1
    //! for debugging: reference counter and owner of data
    count_t               _nreferences;
    MemBlock *            _owner;
    #endif

    //! memory allocation object
    allocator< value_t >  _alloc;
    
    //! @endcond

public:
    //
    // constructors destructor
    //

    //! ctor with nullptr data and 0 references
    MemBlock  () noexcept
            : _data( nullptr )
            , _is_owner( false )
            #if HPRO_MEMBLOCK_REF_COUNT == 1
            , _nreferences( 0 )
            , _owner( nullptr )
            #endif
    {}
    
    //! ctor for \a n elements of \a value_t and 0 references
    MemBlock  ( const size_t  n )
            : _data( nullptr )
            , _is_owner( false )
            #if HPRO_MEMBLOCK_REF_COUNT == 1
            , _nreferences( 0 )
            , _owner( nullptr )
            #endif
    {
        alloc( n );
    }
    
    //! copy ctor (copy reference!)
    MemBlock  ( const MemBlock &  b )
            : _data( b._data )
            , _is_owner( false )
            #if HPRO_MEMBLOCK_REF_COUNT == 1
            , _nreferences( 0 )
            , _owner( b._owner )
            #endif
    {
        #if HPRO_MEMBLOCK_REF_COUNT == 1
        add_ref_owner();
        #endif
    }
    
    //! move ctor (move ownership)
    MemBlock  ( MemBlock &&  b ) noexcept
            : _data( b._data )
            , _is_owner( b._is_owner )
            #if HPRO_MEMBLOCK_REF_COUNT == 1
            , _owner( b._owner )
            #endif
    {
        #if HPRO_MEMBLOCK_REF_COUNT == 1
        _nreferences = b._nreferences;
        
        if ( b._is_owner )
            b._owner = this;

        add_ref_owner();
        #endif

        b._is_owner = false;
    }

    //! ctor with general data pointer \a ptr and size \a n
    MemBlock  ( value_t *   aptr,
                const bool  ais_owner = false )
            : _data( nullptr )
            , _is_owner( false )
            #if HPRO_MEMBLOCK_REF_COUNT == 1
            , _nreferences( 0 )
            , _owner( nullptr )
            #endif
    {
        init( aptr, ais_owner );
    }
    
    //! dtor removing all data <b> even if references exist </b>!
    ~MemBlock  ()
    {
        if ( _is_owner )
        {
            #if HPRO_MEMBLOCK_REF_COUNT == 1
            if ( _nreferences != 0 )
                HERROR( ERR_CONSISTENCY, "(MemBlock) dtor", "freeing data with references" );
            #endif
            
            #if HPRO_DEBUG_MALLOC == 1
            Mem::free( _data );
            #else
            _alloc.deallocate( _data, 0 );
            #endif
        }// if
    }

    //! copy operator (copy reference!)
    MemBlock &  operator = ( const MemBlock &  b )
    {
        _data     = b._data;
        _is_owner = false;
        
        #if HPRO_MEMBLOCK_REF_COUNT == 1
        del_ref_owner();
        
        _owner = b._owner;
        add_ref_owner();
        #endif

        return *this;
    }

    //! move ctor (move ownership)
    MemBlock &  operator = ( MemBlock &&  b ) noexcept
    {
        _data     = b._data;
        _is_owner = b._is_owner;  // only owner if b was owner
        
        #if HPRO_MEMBLOCK_REF_COUNT == 1
        del_ref_owner();

        if ( b._is_owner )
        {
            b._owner     = this;
            _nreferences = b._nreferences;
        }// if
        
        _owner = b._owner;
        add_ref_owner();    // either for this or for b
        #endif

        b._is_owner = false;

        return *this;
    }

    //
    // initialise memory block
    //

    // initialise with raw memory pointer
    void init ( value_t *   ptr,
                const bool  ais_owner = false ) noexcept
    {
        if ( _is_owner )
        {
            #if HPRO_MEMBLOCK_REF_COUNT == 1
            if ( _nreferences != 0 )
                HERROR( ERR_CONSISTENCY, "(MemBlock) init", "freeing data with references" );

            _owner = nullptr;
            #endif
            
            #if HPRO_DEBUG_MALLOC == 1
            Mem::free( _data );
            #else
            _alloc.deallocate( _data, 0 );
            #endif
        }// if

        #if HPRO_MEMBLOCK_REF_COUNT == 1
        del_ref_owner();
        #endif

        _data     = ptr;
        _is_owner = ais_owner;

        #if HPRO_MEMBLOCK_REF_COUNT == 1
        if ( ais_owner )
        {
            _owner       = this;
            _nreferences = 0;
        }// if
        else
        {
            // no reference counting possible, since no MemBlock available
            _owner = nullptr;
        }// else
        #endif
    }

    // initialise with memory block (copy OR move)
    void init ( MemBlock &  b,
                const bool  ais_owner = false )
    {
        if ( ais_owner && ! b._is_owner )
            HERROR( ERR_CONSISTENCY,  "(MemBlock) init", "can not be owner of data NOT owned by given block" );
        
        if ( _is_owner )
        {
            #if HPRO_MEMBLOCK_REF_COUNT == 1
            if ( _nreferences != 0 )
                HERROR( ERR_CONSISTENCY, "(MemBlock) init", "freeing data with references" );

            _owner = nullptr;
            #endif
            
            #if HPRO_DEBUG_MALLOC == 1
            Mem::free( _data );
            #else
            _alloc.deallocate( _data, 0 );
            #endif
        }// if

        _data     = b._data;
        _is_owner = ais_owner;

        #if HPRO_MEMBLOCK_REF_COUNT == 1
        if ( ais_owner && b._is_owner )
        {
            b._owner     = this;
            _nreferences = b._nreferences;
        }// if

        _owner = b._owner;
        add_ref_owner();    // either for this or for b
        #endif
        
        if ( _is_owner )
            b._is_owner = false;
    }

    // allocate and initialise memory block
    void alloc ( const size_t   n,
                 const value_t  init_val = value_t(0) )
    {
        // if ( n > 0 )
        //     std::cout << "alloc" << std::endl;

        alloc_wo_value( n );

        for ( size_t  i = 0; i < n; i++ )
            _data[i] = init_val;
    }

    // only allocate memory block without intialisation
    void alloc_wo_value ( const size_t  n )
    {
        // if ( n > 0 )
        //     std::cout << "alloc" << std::endl;
        
        if ( _is_owner )
        {
            #if HPRO_MEMBLOCK_REF_COUNT == 1
            if ( _nreferences != 0 )
                HERROR( ERR_CONSISTENCY, "(MemBlock) alloc_wo_value", "freeing data with references" );

            _owner = nullptr;
            #endif
            
            #if HPRO_DEBUG_MALLOC == 1
            Mem::free( _data );
            #else
            _alloc.deallocate( _data, 0 );
            #endif
        }// if
        
        #if HPRO_DEBUG_MALLOC == 1
        _data     = static_cast< value_t * >( Mem::alloc( sizeof(value_t) * n ) );
        #else
        _data     = _alloc.allocate( n );
        #endif
        _is_owner = true;

        #if HPRO_MEMBLOCK_REF_COUNT == 1
        _owner       = this;
        _nreferences = 0;
        #endif
    }
    
    //
    // access data
    //

    //! return pointer to internal array
    value_t *       data     ()       noexcept { return _data; }

    //! return const pointer to internal array
    const value_t * data     () const noexcept { return _data; }

    //! return is_owner status
    bool            is_owner () const noexcept { return _is_owner; }
    
    #if HPRO_MEMBLOCK_REF_COUNT == 1

    //
    // reference counting
    //
    
    // add reference
    void add_ref ()
    {
        ++_nreferences;
    }

    // add reference to owner of memblock
    void add_ref_owner ()
    {
        if (( _owner != nullptr ) && ( _owner != this ))
            _owner->add_ref();
    }

    // delete reference
    void del_ref ()
    {
        --_nreferences;
    }

    // delete reference of owner of memblock
    void del_ref_owner ()
    {
        if (( _owner != nullptr ) && ( _owner != this ))
            _owner->del_ref();
    }

    #endif
};

}// namespace BLAS

}// namespace Hpro

#endif  // __HPRO_MEMBLOCK_HH
