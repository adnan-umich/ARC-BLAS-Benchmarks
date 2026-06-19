#ifndef __HPRO_TUPDATEACCUMULATOR_HH
#define __HPRO_TUPDATEACCUMULATOR_HH
//
// Project     : HLIBpro
// File        : TUpdateAccumulator.hh
// Description : class for handling updates to matrix blocks
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <list>
#include <deque>
#include <memory>

#include "hpro/parallel/TMutex.hh"

namespace Hpro
{

// forward decl.
template < typename value_t >
class TMatrix;

//
// base class for direct (e.g. non-recursive) updates
//
template < typename value_t >
class TDirectMatrixUpdate
{
public:
    // dtor
    virtual ~TDirectMatrixUpdate () {}
    
    // apply internal update to matrix \a M with accuracy \a acc
    virtual auto         compute          ( const TTruncAcc &  acc ) -> std::unique_ptr< TMatrix< value_t > > = 0;

    // return linear operator representing update
    virtual auto         linear_op        () const -> std::unique_ptr< TLinearOperator< value_t > > = 0;
    
    // return true if result of update is dense
    virtual bool         is_dense_result  () const = 0;
    
    // return descriptive info
    virtual std::string  to_string        () const { return "TDirectMatrixUpdate"; };
};

//
// base class for recursive updates
//
template < typename value_t >
class TRecursiveMatrixUpdate
{
public:
    // dtor
    virtual ~TRecursiveMatrixUpdate () {}
    
    // apply internal update to matrix \a M with accuracy \a acc
    virtual void         apply      ( TMatrix< value_t > *  M,
                                      const TTruncAcc &     acc ) = 0;

    // return descriptive info
    virtual std::string  to_string  () const { return "TRecursiveMatrixUpdate"; };
};

//!
//! \ingroup Matrix_Module
//! \class   TUpdateAccumulator
//! \brief   Handles updates for a single matrix block by accumulating
//!          direct updates and recursive (pending) updates
//!
template < typename value_t >
class TUpdateAccumulator : public TLockable
{
    //! @cond
    
public:

    // set of direct updates
    using  direct_updates_t    = std::deque< TDirectMatrixUpdate< value_t > * >;
    
    // set of recursive updates
    using  recursive_updates_t = std::list< TRecursiveMatrixUpdate< value_t > * >;
    
private:
    // accumulated updates
    std::unique_ptr< TMatrix< value_t > >  _accumulated;

    // list of pending direct updates
    direct_updates_t                       _pending_direct;

    // list of pending recursive updates
    recursive_updates_t                    _pending_recursive;

    //! @endcond
    
public:
    //! dtor
    ~TUpdateAccumulator ()
    {
        clear_updates();
    }
    
    //! initialise matrix for accumulated updates
    void       init                   ( const TMatrix< value_t > *    M );

    //! compute and apply local direct updates;
    //! - use \a dest as hint for destination type, e.g. to choose format of accumulator
    //! - directly update \a dest if \a update_dest == true (accumulator is zero afterwards)
    void       apply_direct           ( const TTruncAcc &     acc,
                                        TMatrix< value_t > *  dest        = nullptr,
                                        const bool            update_dest = false );

    //! return true if accumulator holds any updates
    bool       has_updates            () const;
    
    //
    // access local updates
    //
    
    //! access accumulated updates
    auto       accumulated_updates    ()       ->       TMatrix< value_t > * { return _accumulated.get(); }
    auto       accumulated_updates    () const -> const TMatrix< value_t > * { return _accumulated.get(); }

    //! access set of direct pending updates
    auto       pending_direct         ()       ->       direct_updates_t &    { return _pending_direct; }
    auto       pending_direct         () const -> const direct_updates_t &    { return _pending_direct; }
    
    //! access set of recursive pending updates
    auto       pending_recursive      ()       ->       recursive_updates_t & { return _pending_recursive; }
    auto       pending_recursive      () const -> const recursive_updates_t & { return _pending_recursive; }
    
    //
    // add updates
    //
    
    //! add update matrix
    void       add_update             ( const TMatrix< value_t > *  M,
                                        const TTruncAcc &           acc,
                                        const TMatrix< value_t > *  dest = nullptr );
    
    //! add update from parent matrix
    void       add_parent_update      ( const TMatrix< value_t > *  M,
                                        const TTruncAcc &           acc );
    
    //! add update U to set of recursive pending updates
    void       add_pending_direct     ( TDirectMatrixUpdate< value_t > *  U )
    {
        if ( U == nullptr )
            return;
        
        TScopedLock  mlock( *this );

        _pending_direct.push_back( U );
    }

    //! add update U to set of recursive pending updates
    void       add_pending_recursive  ( TRecursiveMatrixUpdate< value_t > *  U )
    {
        if ( U == nullptr )
            return;
        
        TScopedLock  mlock( *this );

        _pending_recursive.push_back( U );
    }
    
    //
    // clear stored updates
    //
    
    //! remove matrix with accumulated updates
    void       clear_accumulated      ()
    {
        TScopedLock  mlock( *this );

        _accumulated.reset( nullptr );
    }
    
    //! remove list of pending updates
    void       clear_pending          ();
    
    //! clear all updates
    void       clear_updates          ()
    {
        clear_accumulated();
        clear_pending();
    }

};

}// namespace Hpro

#endif  // __HPRO_TUPDATEACCUMULATOR_HH
