//
// Project     : HLIBpro
// File        : TCluster.cc
// Description : baseclass for a cluster tree
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <iostream>

#include "treealg.hh"
#include "hpro/base/error.hh"
#include "hpro/base/System.hh"

#include "hpro/cluster/TCluster.hh"

namespace Hpro
{

//
// destructor
//

TCluster::~TCluster ()
{
    for ( size_t  i = 0; i < nsons(); i++ )
        delete _sons[i];
}

//
// set number of sons
//
void
TCluster::set_nsons ( const size_t  n )
{
    const size_t  old = nsons();
    
    _sons.resize( n );

    for ( size_t  i = old; i < n; i++ )
        _sons[i] = nullptr;
}
    
//
// set sons
//

void
TCluster::set_son ( const size_t  i,
                    TCluster *    son_ct,
                    const bool    del )
{
    if ( del && (_sons[i] != nullptr) && (_sons[i] != son_ct))
        delete _sons[i];

    _sons[i] = son_ct;

    if ( son_ct != nullptr )
        son_ct->set_parent( this );
}

void
TCluster::add_son ( TCluster *  son_ct,
                    const bool  inc_nsons )
{
    if ( son_ct == nullptr )
        return;

    //
    // look for free sons
    //

    for ( size_t  i = 0; i < nsons(); i++ )
    {
        // stop, if cluster is already a son
        if ( _sons[i] == son_ct )
            return;
        
        // stop, if free slot was found
        if ( _sons[i] == nullptr )
        {
            _sons[i] = son_ct;
            return;
        }// if
    }// for

    if ( inc_nsons )
    {
        // create free slot
        _sons.resize( nsons() + 1 );
        _sons[ nsons() - 1 ] = son_ct;
        son_ct->set_parent( this );
    }// if
    else
        HERROR( ERR_CONSISTENCY, "(TCluster) add_son", "no free slot" );
}

//
// return size of tree
//
size_t
TCluster::nnodes () const
{
    return Hpro::Tree::nnodes( this );
}

//
// return depth of tree
//
size_t
TCluster::depth () const
{
    return Hpro::Tree::depth( this );
}

//
// collect leaves in given list
//
void
TCluster::collect_leaves ( std::list< TCluster * > &  leaves,
                           const int                  tdepth,
                           const int                  level ) const
{
    if ( is_leaf() || ((tdepth != -1) && (level == tdepth)))
        leaves.push_back( const_cast< TCluster * >( this ) );
    else
    {
        for ( size_t i = 0; i < nsons(); i++ )
            if ( _sons[i] != nullptr )
                _sons[i]->collect_leaves( leaves, tdepth, level+1 );
    }// for
}

//
// virtual constructor
//
TCluster *
TCluster::copy  () const
{
    auto  c = create();

    c->set_id( id() );
    c->set_parent( _parent );
    c->set_first_last( first(), last() );
    c->set_nsons( nsons() );
    c->set_domain( is_domain() );

    for ( size_t i = 0; i < nsons(); i++ )
    {
        if ( son(i) != nullptr )
            c->set_son( i, son(i)->copy() );
    }// for

    return c;
}

//
// stream output
//
void
TCluster::print ( const uint ofs ) const
{
    for ( uint i = 0; i < ofs; i++ )
        std::cout << ' ';

    std::cout << *this << std::endl;

    for ( size_t i = 0; i < nsons(); i++ )
    {
        if ( son(i) != nullptr )
            son(i)->print( ofs + 4 );
    }// for
}

//
// return size in bytes used by this object
//
size_t
TCluster::byte_size () const
{
    size_t count = ( TIndexSet::byte_size() + sizeof(_id) + sizeof(_parent) + sizeof(_sons) + sizeof(TCluster*) * _sons.size() + sizeof(bool) );
    
    for ( size_t i = 0; i < nsons(); i++ )
    {
        if ( son(i) != nullptr )
            count += son(i)->byte_size();
    }// for
    
    return count;
}

//
// flatten hierarchy of cluster tree
//
void
flatten ( TCluster * cl )
{
    if (( cl == nullptr ) || ( cl->is_leaf() ))
        return;

    //
    // recurse and flatten in sons and collect new son clusters
    //

    std::list< TCluster * >  new_sons;
    size_t                   n_new_sons = 0;

    for ( size_t  i = 0; i < cl->nsons(); ++i )
    {
        auto  son_i = cl->son(i);
        
        if ( son_i != nullptr )
        {
            if ( son_i->is_leaf() )
            {
                new_sons.push_back( son_i );
                n_new_sons++;
                cl->set_son( i, nullptr, false );
            }// if
            else
            {
                flatten( son_i );

                for ( size_t  j = 0; j < son_i->nsons(); ++j )
                {
                    new_sons.push_back( son_i->son(j) );
                    n_new_sons++;
                }// for
                                    

                // delete son
                son_i->set_nsons( 0 );
                cl->set_son( i, nullptr );
            }// else
        }// if
    }// for

    //
    // replace sons by new sons
    //

    size_t  pos = 0;

    // std::cout << "nnew_sons = " << new_sons.size() << std::endl;
    
    cl->set_nsons( n_new_sons );

    for ( auto  new_son : new_sons )
    {
        cl->set_son( pos++, new_son );
    }// for
}

}// namespace Hpro
