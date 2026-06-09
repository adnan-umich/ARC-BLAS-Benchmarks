//
// Project     : HLIBpro
// File        : TUpdateAccumulator.cc
// Description : class for handling updates to matrix blocks
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TBlockMatrix.hh"

#include "TRNG.hh"
#include "list.hh"

#include "hpro/matrix/structure.hh"

namespace Hpro
{

namespace B = BLAS;

using std::deque;
using std::pair;
using std::list;
using std::make_pair;
using std::unique_ptr;

//
// initialise matrix for accumulated updates
//
template < typename value_t >
void
TUpdateAccumulator< value_t >::init ( const TMatrix< value_t > *  M )
{
    TScopedLock  mlock( *this );

    if ( _accumulated.get() == nullptr )
    {
        if ( is_dense( M ) )
        {
            _accumulated = std::make_unique< TDenseMatrix< value_t > >( M->row_is(), M->col_is() );
        }// if
        else
        {
            _accumulated = std::make_unique< TRkMatrix< value_t > >( M->row_is(), M->col_is() );
        }// else
    }// if
}

//
// return true if accumulator holds any updates
//
template < typename value_t >
bool
TUpdateAccumulator< value_t >::has_updates () const
{
    if ( pending_direct().empty() && pending_recursive().empty() && ( accumulated_updates() == nullptr ))
        return false;

    return true;
}
    
//
// compute sum of local direct updates
//
template < typename value_t >
void
TUpdateAccumulator< value_t >::apply_direct ( const TTruncAcc &     /* acc */,
                                              TMatrix< value_t > *  /* M */,
                                              const bool            /* update_M */ )
{
    HERROR( ERR_NOT_IMPL, "", "" );
}
    

//
// remove list of pending updates
//
template < typename value_t >
void
TUpdateAccumulator< value_t >::clear_pending ()
{
    TScopedLock  mlock( *this );
    
    for ( auto  upd : _pending_direct )
        delete upd;
    
    _pending_direct.clear();

    for ( auto  upd : _pending_recursive )
        delete upd;
    
    _pending_recursive.clear();
}

//
// add update matrix
//
template < typename value_t >
void
TUpdateAccumulator< value_t >::add_update ( const TMatrix< value_t > *  /* M */,
                                            const TTruncAcc &           /* acc */,
                                            const TMatrix< value_t > *  /* dest */ )
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// add update from parent matrix
//
template < typename value_t >
void
TUpdateAccumulator< value_t >::add_parent_update ( const TMatrix< value_t > * /* M */,
                                                   const TTruncAcc &          /* acc */ )
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// explicit instantiation
//

template class TUpdateAccumulator< float >;
template class TUpdateAccumulator< double >;
template class TUpdateAccumulator< std::complex< float > >;
template class TUpdateAccumulator< std::complex< double > >;

}// namespace Hpro
