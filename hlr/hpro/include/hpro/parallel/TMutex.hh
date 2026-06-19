#ifndef __HPRO_TMUTEX_HH
#define __HPRO_TMUTEX_HH
//
// Project     : HLIBpro
// File        : TMutex.hh
// Description : class for a mutex
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <mutex>

#include "hpro/base/types.hh"

namespace Hpro
{

//!
//! \class  TMutex
//! \brief  Wraps default mutices.
//!
class TMutex : public std::mutex
{
public:
    /////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct non-locked mutex
    TMutex ()
            : std::mutex()
    {}

    /////////////////////////////////////////////////
    //
    // misc.
    //

    //! return size in bytes used by this object
    size_t  byte_size () const { return sizeof( std::mutex ); }
};

//!
//! \class  TLockable
//! \brief  Base class for all mutex equipped classes
//!
class TLockable
{
private:
    //! @cond
    
    //! given mutex for locking
    TMutex  _mutex;

    //! @endcond
    
public:
    
    //! give access to internal mutex
    TMutex &  mutex   ()  { return _mutex; }

    //! lock local mutex
    void      lock    ()  { _mutex.lock(); }

    //! unlock local mutex
    void      unlock  ()  { _mutex.unlock(); }

    /////////////////////////////////////////////////
    //
    // misc.
    //

    //! return size in bytes used by this object
    size_t  byte_size () const
    {
        return sizeof( _mutex );
    }
};
    
//!
//! \class  TScopedLock
//! \brief  Provides automatic lock and unlock for mutices.
//!
class TScopedLock : public std::lock_guard< TMutex >
{
private:
    //! @cond
    
    // prevent copy operations
    TScopedLock ( const TScopedLock & );
    void operator = ( const TScopedLock & );

    //! @endcond
    
public:
    //! ctor: lock mutex
    explicit TScopedLock ( TMutex &  m )
            : std::lock_guard< TMutex >( m )
    {}

    //! ctor: lock lockable object
    explicit TScopedLock ( TLockable &  m )
            : std::lock_guard< TMutex >( m.mutex() )
    {}
};

#define LOCK_EXPR( mtx, expr ) { TScopedLock  __lock_##mtx( mtx ); expr; }

}// namespace Hpro

#endif  // __HPRO_TMUTEX_HH
