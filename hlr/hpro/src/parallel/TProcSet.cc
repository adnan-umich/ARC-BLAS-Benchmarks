//
// Project     : HLIBpro
// File        : TProcSet.cc
// Description : class for representing processor sets
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/System.hh"

#include "hpro/parallel/TProcSet.hh"

namespace Hpro
{

//
// represents invalid processor set
//
const TProcSet  PROCSET_INVALID( Limits::max< uint >(), Limits::max< uint >() );

///////////////////////////////////////////////////
//
// serialisation
//
///////////////////////////////////////////////////

//
// read data from stream
//
void
TProcSet::read ( TByteStream &  bs )
{
    bs.get( _first );
    bs.get( _last );
}

//
// write data to stream
//
void
TProcSet::write ( TByteStream &  bs ) const
{
    bs.put( _first );
    bs.put( _last );
}

//
// returns size of object in bytestream
//
size_t
TProcSet::bs_size () const
{
    return sizeof(_first) + sizeof(_last);
}

}// namespace
