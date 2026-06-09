#ifndef __HPRO_TRNG_HH
#define __HPRO_TRNG_HH
//
// Project     : HLIBpro
// File        : TRNG.hh
// Description : class for a random-number-generator
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

namespace Hpro
{

//!
//! \class  TRNG
//! \brief  internal random number generator based on Boost
//!
class TRNG
{
private:
    using  generator_t = boost::mt19937;
    using  distr_t     = boost::uniform_real<>;
    using  rng_t       = boost::variate_generator< generator_t, distr_t >;
    
    // RNG data
    generator_t  _generator;
    distr_t      _distribution;
    rng_t        _rng;

public:
    ////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct RNG with given seed value
    //! (0: internal initialisation)
    TRNG ( const unsigned long seed = 0 );

    ////////////////////////////////////////////
    //
    // access rng
    //

    //! return random number in the interval [0, \a max]
    template < typename value_t >
    value_t  rand          ( const value_t        max = value_t(1) ) { return max * value_t( _rng() ); }

    //! abbr. for \see rand
    template < typename value_t >
    value_t   operator ()  ( const value_t        max = value_t(1) ) { return this->rand( max ); }

    //! initialise rng by single seed (seed must not equal 0)
    void     init          ( const unsigned long  seed );
};

}// namespace

#endif  // __HPRO_TRNG_HH
