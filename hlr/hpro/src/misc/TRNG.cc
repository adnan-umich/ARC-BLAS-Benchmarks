//
// Project     : HLIBpro
// File        : TRNG.cc
// Description : class for a random-number-generator
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <atomic>

#include <time.h>

#include "hpro/base/types.hh"

#include "TRNG.hh"

namespace Hpro
{

namespace
{

//
// return seed for RNGs
//
ulong
rand_seed ()
{
    static  std::atomic< ulong >  _seed( ulong(::time( nullptr )) );

    return _seed++;
}

}// namespace anonymous

//
// ctor
//
TRNG::TRNG ( const unsigned long seed )
        : _generator( seed == 0 ? rand_seed() : seed )
        , _distribution( 0, 1 )
        , _rng( _generator, _distribution )
{}

//
// initialise rng
//
void
TRNG::init ( const unsigned long aseed )
{
    _rng.engine().seed( aseed );
    _rng.distribution().reset();
}

}// namespace Hpro
