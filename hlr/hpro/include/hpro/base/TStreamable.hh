#ifndef __HPRO_TSTREAMABLE_HH
#define __HPRO_TSTREAMABLE_HH
//
// Project     : HLIBpro
// File        : TStreamable.hh
// Description : baseclass for all streamable classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/types.hh"
#include "hpro/base/TByteStream.hh"
#include "hpro/parallel/TProcSet.hh"

namespace Hpro
{

//
// defines basic routines for reading/writing from/to bytestreams
// and for distribution of data
//
class TStreamable
{
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TStreamable () {}

    virtual ~TStreamable () {}

    ////////////////////////////////////////////////
    //
    // bytestream IO
    //
    
    // read object from stream
    virtual void read  ( TByteStream & ) {}

    // write object to stream 
    virtual void write ( TByteStream & ) const {}

    // returns size of object in bytestream
    virtual size_t bs_size () const  { return 0; }

    ////////////////////////////////////////////////
    //
    // data distribution
    //

    // distribute local data to all processors in group
    // (if bs != NULL it will be used)
    virtual void scatter ( const TProcSet & p,
                           const uint       pid,
                           TByteStream    * bs = NULL );
    virtual void scatter ( const TProcSet & procs );
    
    ////////////////////////////////////////////////
    //
    // reduce operations
    //

    // sum up nparts parallel copies
    // (if bs != NULL it will be used)
//     virtual void sum ( const TProcSet & p,
//                        const uint       pid,
//                        const uint       nparts,
//                        TByteStream    * bs = NULL );
};

}// namespace

#endif  // __HPRO_TSTREAMABLE_HH
