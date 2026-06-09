//
// Project     : HLIBpro
// File        : TStreamable.cc
// Description : baseclass for all streamable classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/error.hh"
#include "hpro/parallel/NET.hh"

#include "hpro/base/TStreamable.hh"

namespace Hpro
{

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//
// local routines (decl.)
//
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

namespace
{

//
// internal scatter routine
//
void scatter_bs ( const TProcSet & procs,
                  const uint       pid,
                  TStreamable    * S,
                  TByteStream    * bs );

}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//
// TStreamable
//
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////
//
// data distribution
//

//
// distribute local data to all processors in group
// (if bs != nullptr it will be used)
//
void
TStreamable::scatter ( const TProcSet & p,
                       const uint       pid,
                       TByteStream    * bs )
{
    //
    // write S to bytestream at local master and start
    // real distribution algorithm
    //
    
    TByteStream  * ibs;
    bool           del_bs = false;
    
    if ( bs != nullptr ) 
        ibs = bs;
    else
    {
        ibs    = new TByteStream;
        del_bs = true;
    }// else

    write( ibs->set_size( bs_size() ) );
    
    scatter_bs( p, pid, this, ibs );

    if ( del_bs )
        delete ibs;
}

void
TStreamable::scatter ( const TProcSet & procs )
{
    //
    // write S to bytestream at local master and start
    // real distribution algorithm
    //
    
    TByteStream  bs;
    
    write( bs.set_size( bs_size() ) );
    
    scatter_bs( procs, NET::pid(), this, & bs );
}

////////////////////////////////////////////////
//
// reduce operations
//

//
// sum up nparts parallel copies
// (if bs != nullptr it will be used)
//
// void
// TStreamable::sum ( const TProcSet & p,
//                    const uint       pid,
//                    const uint       nparts,
//                    TByteStream    * bs )
// {
//     HERROR( ERR_NOT_IMPL, "(TStreamable) sum", "" );
// }
    
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//
// local routines (impl.)
//
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

namespace
{

//
// scatter data of master in <procs> to all other processors
// in <procs>
//
void
scatter_bs ( const TProcSet & ps,
             const uint       pid,
             TStreamable *    S,
             TByteStream *    bs )
{
    TProcSet  procs( ps );
    
    while ( procs.size() > 1 )
    {
        //
        // send result to adjacent processor set
        //

        if ( pid == procs.master() )
            NET::dsend( procs.subset(1).master(), bs->data(), bs->size() );
        
        if ( pid == procs.subset(1).master() )
        {
            size_t  n = NET::dprobe( procs.master() );

            bs->set_size( n );
            
            NET::drecv( procs.master(), bs->data(), bs->size() );

            S->read( * bs );
        }// if

        //
        // next level of distribution
        //

        if ( procs.subset(0).is_in( pid ) ) procs = procs.subset(0);
        else                                procs = procs.subset(1);
    }// while
}

}// namespace anonymous

}// namespace Hpro
