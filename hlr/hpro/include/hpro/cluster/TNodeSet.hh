#ifndef __HPRO_TNODESET_HH
#define __HPRO_TNODESET_HH
//
// Project     : HLIBpro
// File        : TNodeSet.hh
// Description : class for representing a set of nodes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <cstring>
#include <vector>

#include "hpro/base/types.hh"

namespace Hpro
{

//
// iterate through nodeset
//
//#define FOREACH_NS( iter, set )  for ( TNodeSet::TIterator  iter = (set).begin(); ! iter.eol(); ++iter )

//!
//! \ingroup  Cluster_Module
//! \typedef  node_t
//! \brief    Basic type for storing node names.
//!
using  node_t = idx_t;

//!
//! \ingroup  Cluster_Module
//! \class    TNodeSet
//! \brief    Represents a set of nodes by an array
//!           - maximal size of nodeset, e.g. the capacity, is equal to array-size
//!           - current number of members is equal to _n_nodes
//!
class TNodeSet
{
public:
    //!
    //! \class    iterator
    //! \brief    STL iterator for TNodeSet
    //!
    struct iterator
    {
    private:
        //! referenced node set
        const TNodeSet *  _ns;
        
        //! current position
        size_t            _pos;
        
    public:
        iterator ( const TNodeSet *  ns,
                   const size_t      pos )
                : _ns(ns)
                , _pos(pos)
        {}

        iterator  operator ++ ()    noexcept { return iterator( _ns, ++_pos ); }
        iterator  operator ++ (int) noexcept { return iterator( _ns, _pos++ ); }

        iterator  operator -- ()    noexcept { return iterator( _ns, --_pos ); }
        iterator  operator -- (int) noexcept { return iterator( _ns, _pos-- ); }
        
        node_t    operator *  () const noexcept { return _ns->operator[](_pos); }

        bool      operator == ( const iterator &  it ) const noexcept { return (_pos == it._pos) && (_ns == it._ns); }
        bool      operator != ( const iterator &  it ) const noexcept { return (_pos != it._pos) || (_ns != it._ns); }
    };

    using  const_iterator = iterator;

protected:
    //! array containing nodeset
    std::vector< node_t >  _data;
    
    //! actual number of nodes in set
    size_t                 _nnodes;

public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct node set with maximal capacity \a set_size
    TNodeSet ( const size_t set_size = 0 )
            : _data( set_size )
            , _nnodes(0)
    {}

    //! copy ctor
    TNodeSet ( const TNodeSet & set )
            : _nnodes(0)
    {
        *this = set;
    }
    
    //////////////////////////////////////
    //
    // map array functions
    //

    //! return capacity of node set
    size_t    size        () const noexcept { return _data.size(); }
    
    //! return \a i'th node in set
    node_t &  operator [] ( const size_t  i )       noexcept { return _data[i]; }
    node_t    operator [] ( const size_t  i ) const noexcept { return _data[i]; }

    //! return STL iterators
    iterator        begin  () const noexcept { return iterator( this, 0 ); } 
    iterator        end    () const noexcept { return iterator( this, nnodes() ); } 

    //! return STL const iterators
    const_iterator  cbegin () const noexcept { return iterator( this, 0 ); } 
    const_iterator  cend   () const noexcept { return iterator( this, nnodes() ); } 
    
    //////////////////////////////////////
    //
    // manage node set
    //

    //! return number of stored nodes
    size_t  nnodes () const noexcept { return _nnodes; }

    //! directly set number of stored nodes
    void    set_nnodes ( const size_t n ) noexcept { _nnodes = n; }

    //! append node \a node to the set (without incrementing the array-size)
    void    append ( const node_t  node )
    {
        _data[ _nnodes++ ] = node;
    }

    //! remove all nodes from set (array-size is not changed)
    void    remove_all () noexcept
    {
        _nnodes = 0;
    }
    
    //! adjust capacity to \a n <b>and</b> node number
    void resize ( const size_t n, const bool = true )
    {
        _data.resize( n );
        _nnodes = std::min< size_t >( _nnodes, size() );
    }

    //! adjust capacity to \a n <b>and</b> node number
    void resize ( const size_t n, const uint &, const bool )
    {
        _data.resize( n );
        _nnodes = std::min< size_t >( _nnodes, size() );
    }

    //! copy operator
    TNodeSet & operator = ( const TNodeSet & set )
    {
        resize( set.size() );
        
        _nnodes = set._nnodes;

        if ( _nnodes > 0 )
            memcpy( _data.data(), set._data.data(), sizeof(node_t) * _nnodes );
        
        return *this;
    }
};

}// namespace Hpro

#endif  // __HPRO_TNODESET_HH
