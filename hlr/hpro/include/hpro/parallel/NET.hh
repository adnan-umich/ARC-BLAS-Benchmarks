#ifndef __HPRO_NET_HH
#define __HPRO_NET_HH
//
// Project     : HLIBpro
// File        : NET.hh
// Description : encapsulation of network functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <stdlib.h>

#include "hpro/config.h"

#if defined(__GNUC__) && !defined(__ICC)
#  if __GNUC__ >= 4  &&  __GNUC_MINOR__ >= 6
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wold-style-cast"
#  endif
#endif

#if HPRO_NET_TYPE == 2
#include <mpi.h>
#endif

#if defined(__GNUC__) && !defined(__ICC)
#  if __GNUC__ >= 4  &&  __GNUC_MINOR__ >= 6
#    pragma GCC diagnostic pop
#  endif
#endif

#include "hpro/base/types.hh"
#include "hpro/base/error.hh"
#include "hpro/parallel/TProcSet.hh"

namespace Hpro
{

namespace NET
{

/////////////////////////////////////////////////////
//
// reduce operations
//
enum reduce_op_t
{
    OP_MAX,     // pointwise maximum of input arrays
    OP_SUM      // pointwise summation of input arrays
};

//
// non-blocking request
//
#if HPRO_NET_TYPE == 2
using  request_t = MPI_Request;
#else
using  request_t = void *;
#endif

// standard definition for no processor
const uint NO_PROC = Limits::max<uint>();

////////////////////////////////////////
//
// initialization and ending
//

void init ();
void done ();

////////////////////////////////////////
//
// properties of parallel machine
//

//! return number of processors in parallel machins
uint nprocs ();
void set_nprocs ( uint p ); // DEBUG

//! return local processor id
uint pid    ();
void set_pid ( uint p ); // DEBUG
    
//! return true, if local processor is equal to given processor
inline
bool is_local  ( const uint  p )
{
    return ( ( p == NO_PROC ) || ( p == pid() ) );
}

//! return true, if local processor is equal to given processor set
inline
bool is_local  ( const TProcSet &  ps )
{
    return ( ( ps != PROCSET_INVALID ) && ( ps.size() == 1 ) && ( ps.first() == pid() ));
}

//! return true, if local processor is element of given processor set
inline
bool contains_local  ( const TProcSet &  ps )
{
    return ( ( ps == PROCSET_INVALID ) || ps.is_in( pid() ) );
}

////////////////////////////////////////
//
// direct (non-BSP) communication
// (blocking)
//

//
// blocking point-to-point communication
//

// direct send/recv
void    dsend   ( const uint  dest,   const void * buf, size_t bsize );
void    drecv   ( const uint  source,       void * buf, size_t bsize, int tag = 0 );

//! return size of message to receive from node \a source
size_t  dprobe  ( const uint  source, const int tag = 0 );
               
//
// nonblocking point-to-point communication
//

// internal send function
request_t
isend_intern ( uint          dest,
               const void *  buf,
               size_t        bsize,
               int           tag );

// send function
template <typename T>
request_t
isend ( const uint  dest,
        const T *   buf,
        size_t      count,
        int         tag = 0 )
{
    return isend_intern( dest, buf, sizeof(T) * count, tag );
}

// single data send
template <typename T>
request_t
isend1 ( const uint  dest,
         const T &   buf,
         int         tag = 0 )
{
    return isend_intern( dest, & buf, sizeof(T), tag );
}

// internal send function
request_t
irecv_intern ( uint          source,
               const void *  buf,
               size_t        bsize,
               int           tag );

// recieve function
template <typename T>
request_t
irecv ( const uint  source,
        const T *   buf,
        size_t      count,
        int         tag = 0 )
{
    return irecv_intern( source, buf, sizeof(T) * count, tag );
}

// single data recieve
template <typename T>
request_t
irecv1 ( const uint  source,
         const T &   buf,
         int         tag = 0 )
{
    return irecv_intern( source, & buf, sizeof(T), tag );
}

//
// wait for communication request to finish
//
void
wait ( request_t &  req );

//
// wait for one or more communication requests to finish
// - return number of finished requests
//
int
wait_some ( const std::vector< request_t > &  reqs,
            std::vector< int > &              inds );

// test single requests
bool
test ( request_t &  req );

//
// collective communication
//

//! reduce data in \a inbuf of all nodes in \a ps as
//! defined by \a op into \a outbuf of master node of \a ps
template <typename T>
void
reduce       ( const TProcSet &   ps,
               const T *          inbuf,
               T *                outbuf,
               const size_t       count,
               const reduce_op_t  op );
    
//! reduce data in \a inbuf of all nodes in \a ps as
//! defined by \a op into \a outbuf of all nodes in \a ps
template <typename T>
void
reduce_all   ( const TProcSet &   ps,
               const T *          inbuf,
               T *                outbuf,
               const size_t       count,
               const reduce_op_t  op );

//! reduce single data \a inbuf of all nodes in \a ps as
//! defined by \a op and return result 
template <typename T>
T
reduce_all   ( const TProcSet &   ps,
               const T            inbuf,
               const reduce_op_t  op )
{
    T  outbuf;

    reduce_all( ps, & inbuf, & outbuf, 1, op );

    return  outbuf;
}
    
//! internal broadcast function
void
broadcast_intern ( const uint        root,
                   const TProcSet &  ps,
                   void *            buffer,
                   const size_t      count );
    
//! send \a count elements of data in \a buffer of \a root
//! to \a buffer of all nodes in \a ps
template <typename T>
void
broadcast    ( const uint        root,
               const TProcSet &  ps,
               T *               buffer,
               const size_t      count )
{
    broadcast_intern( root, ps, buffer, sizeof(T) * count );
}
    
//! send \a count elements of data in \a buffer of master of \a ps
//! to \a buffer of all nodes in \a ps
template <typename T>
void
broadcast    ( const TProcSet &  ps,
               T *               buffer,
               const size_t      count )
{
    broadcast< T >( ps.master(), ps, buffer, count );
}
    
}// namespace NET

}// namespace Hpro

#endif  // __HPRO_NET_HH
