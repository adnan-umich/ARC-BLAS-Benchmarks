//
// Project     : HLIBpro
// File        : NET.cc
// Description : encapsulation of network functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <cstring>
#include <list>

#include "list.hh"

#include "hpro/base/System.hh"

#include "hpro/parallel/NET.hh"

namespace Hpro
{

//
// different types of network modes
//

#define NET_TYPE_SEQ  1
#define NET_TYPE_MPI  2
#define NET_TYPE_SHM  3

namespace
{

////////////////////////////////////////
//
// local data
//

// flag which indicates that network was init.
bool  net_initialised = false;

// structure of the parallel machine
int   net_pid    = 0;
int   net_nprocs = 1;

}// namespace anonymous

namespace NET
{

////////////////////////////////////////
//
// initialization and ending
//

void
init ()
{
    if ( ! net_initialised )
        net_initialised = true;
}

void
done ()
{
    if ( ! net_initialised )
        HERROR( ERR_INIT, "(NET) done", "network not initialised" );

    net_initialised = false;
}

////////////////////////////////////////
//
// properties of parallel machine
//

uint
nprocs ()
{
    return net_nprocs;
}

void
set_nprocs ( uint p )
{
    net_nprocs = p;
}

uint
pid ()
{
    return net_pid;
}

void
set_pid ( uint p )
{
    net_pid = p;
}
    
////////////////////////////////////////
//
// collective communication
//

//
// reduce data in \a inbuf of all nodes in \a ps as
// defined by \a op into \a outbuf of master node of \a ps
//
template <typename T>
void
reduce ( const TProcSet &   ps,
         const T *,
         T *,
         const size_t,
         const reduce_op_t )
{
    if ( ! net_initialised )
        HERROR( ERR_INIT, "(NET) reduce", "network not initialised" );

    if (( ps.size() > 1 ) || (( ps.size() == 1 ) && ( ps.first() != 0 )))
        HERROR( ERR_ARG, "(NET) reduce", "processor set invalid (sequential computation!)" );
}
    
#define  INST_REDUCE( type )                         \
    template                                         \
    void                                             \
    reduce< type > ( const TProcSet &   ps,          \
                     const type *       inbuf,       \
                     type *             outbuf,      \
                     const size_t       count,       \
                     const reduce_op_t  op )

INST_REDUCE( size_t );
INST_REDUCE( float );
INST_REDUCE( double );
INST_REDUCE( std::complex<float> );
INST_REDUCE( std::complex<double> );

//
// reduce data in \a inbuf of all nodes in \a ps as
// defined by \a op into \a outbuf of all nodes in \a ps
//
template <typename T>
void
reduce_all ( const TProcSet &   ps,
             const T *,
             T *,
             const size_t,
             const reduce_op_t )
{
    if ( ! net_initialised )
        HERROR( ERR_INIT, "(NET) reduce_all", "network not initialised" );

    if (( ps.size() > 1 ) || (( ps.size() == 1 ) && ( ps.first() != 0 )))
        HERROR( ERR_ARG, "(NET) reduce_all", "processor set invalid (sequential computation!)" );
}

#define  INST_REDUCE_ALL( type )                       \
    template                                           \
    void                                               \
    reduce_all< type >  ( const TProcSet &   ps,       \
                          const type *       inbuf,    \
                          type *             outbuf,   \
                          const size_t       count,    \
                          const reduce_op_t  op )

INST_REDUCE_ALL( size_t );
INST_REDUCE_ALL( float );
INST_REDUCE_ALL( double );
INST_REDUCE_ALL( std::complex<float> );
INST_REDUCE_ALL( std::complex<double> );

//
// send \a count elements of data in \a buffer of master of \a ps
// to \a buffer of all nodes in \a ps
//
void
broadcast_intern ( const uint,
                   const TProcSet &  ps,
                   void *,
                   const size_t )
{
    if ( ! net_initialised )
        HERROR( ERR_INIT, "(NET) broadcast", "network not initialised" );

    if (( ps.size() > 1 ) || (( ps.size() == 1 ) && ( ps.first() != 0 )))
        HERROR( ERR_ARG, "(NET) broadcast", "processor set invalid (sequential computation!)" );
}

////////////////////////////////////////
//
// direct (non-BSP) communication
//

void
dsend ( const unsigned int dest, const void *, const size_t )
{
    if ( ! net_initialised )
        HERROR( ERR_INIT, "(NET) dsend", "network not initialised" );

    if ( dest != 0 )
        HERROR( ERR_ARG, "(NET) dsend", "invalid destination" );
    
    if ( dest == 0 )
        HERROR( ERR_ARG, "(NET) dsend", "source = destination" );
}

void
drecv ( const unsigned int source, void *, const size_t, const int )
{
    if ( ! net_initialised )
        HERROR( ERR_INIT, "(NET) drecv", "network not initialised" );

    if ( source != 0 )
        HERROR( ERR_ARG, "(NET) drecv", "invalid source" );
    
    if ( source == 0 )
        HERROR( ERR_ARG, "(NET) drecv", "source = destination" );
}

//
// return size of message to receive
//
size_t
dprobe ( const uint source,
         const int )
{
    if ( ! net_initialised )
        HERROR( ERR_INIT, "(NET) dprobe", "network not initialised" );

    if ( source != 0 )
        HERROR( ERR_INIT, "(NET) dprobe", "invalid source" );
    
    if ( source == 0 )
        HERROR( ERR_INIT, "(NET) dprobe", "source = destination" );

    return 0;
}

//
// direct send
//
request_t
isend_intern ( uint          dest,
               const void *,
               size_t,
               int )
{
    if ( ! net_initialised )
        HERROR( ERR_INIT, "(NET) isend", "network not initialised" );

    if ( dest != 0 )
        HERROR( ERR_ARG, "(NET) isend", "invalid destination" );
    
    if ( dest == 0 )
        HERROR( ERR_ARG, "(NET) isend", "source = destination" );

    return nullptr;
}

//
// direct send
//
request_t
irecv_intern ( uint          dest,
               const void *,
               size_t,
               int )
{
    if ( ! net_initialised )
        HERROR( ERR_INIT, "(NET) irecv", "network not initialised" );

    if ( dest != 0 )
        HERROR( ERR_ARG, "(NET) irecv", "invalid destination" );
    
    if ( dest == 0 )
        HERROR( ERR_ARG, "(NET) irecv", "source = destination" );

    return nullptr;
}

//
// wait for corresponding receive
//
void
wait ( request_t & )
{
    if ( ! net_initialised )
        HERROR( ERR_INIT, "(NET) isend", "network not initialised" );

}

int
wait_some ( const std::vector< request_t > &,
            std::vector< int > & )
{
    if ( ! net_initialised )
        HERROR( ERR_INIT, "(NET) isend", "network not initialised" );

    return 0;
}

//
// test single requests
//
bool
test ( request_t & )
{
    if ( ! net_initialised )
        HERROR( ERR_INIT, "(NET) isend", "network not initialised" );

    return true;
}

}// namespace NET

}// namespace Hpro
