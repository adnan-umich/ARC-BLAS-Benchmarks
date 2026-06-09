//
// Project     : HLIBpro
// File        : init.cc
// Description : initialisation and finalisation
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <stdlib.h>
#include <string.h>

#include <boost/lexical_cast.hpp>

#include "hpro/base/config.hh"
#include "hpro/base/error.hh"
#include "hpro/base/types.hh"
#include <hpro/blas/cuda.hh>
#include <hpro/blas/hip.hh>
#include "hpro/parallel/NET.hh"

namespace Hpro
{

namespace
{

//
// flag for indicating init status of HLIBpro
//
bool hlib_is_init = false;

//
// inquire environment variable and return true if 
// available or <false> otherwise
//
bool
get_environ ( const std::string &  envname,
              std::string &        content )
{
    char  * val = nullptr;
    
#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)

    size_t  val_len  = 0;

    if ( ! _dupenv_s( & val, & val_len, envname.c_str() ) && ( val != nullptr ))
    {
        content = val;
        free( val );
    }// if
#else
    if (( val = ::getenv( envname.c_str() ) ) != nullptr )
        content = val;
#endif 
    else
    {
        content = "";
        return false;
    }// else

    return true;
}

}// namespace anonymous

///////////////////////////////////////////////////////////////////
//
// global initialisation routine for HLib
//
///////////////////////////////////////////////////////////////////

void
INIT ()
{
    if ( hlib_is_init )
        return;
    
    //
    // set up logging
    //

    LOG::init();

    //
    // set up configuration management
    //

    CFG::init();

    //
    // set up distributed environment
    //
    
    NET::init();
    
    //
    // read environment variable for verbosity
    //

    std::string  value;
    int          verbosity = CFG::verbosity;
    
    if ( get_environ( "HLIB_VERBOSITY", value ) )
        verbosity = std::max( 0, boost::lexical_cast< int >( value ) );
    if ( get_environ( "HPRO_VERBOSITY", value ) )
        verbosity = std::max( 0, boost::lexical_cast< int >( value ) );

    CFG::set_verbosity( verbosity );

    if ( verbose( LOG_INFO ) )
        CFG::print_parameters();
    
    if ( verbose( LOG_DEBUG ) )
        RTTI::print_registered();
    
    HINFO( "(INIT) initial memory consumption: " + Mem::to_string() );
    
    //
    // set up GPUs
    //

    CUDA::init();
    HIP::init();

    //
    // all is ready
    //
    
    hlib_is_init = true;
}

//
// finalisation routing for HLib
//
void
DONE ()
{
    if ( ! hlib_is_init )
        return;
    
    //
    // finish parallel stuff
    //

    NET::done();

    //
    // finish config system
    //

    CFG::done();

    //
    // finish logging
    //

    LOG::done();

    //
    // print maximal memory consumption
    //
    
    HINFO( "(DONE) maximal memory consumption: " + Mem::to_string( Mem::max_usage() ) );
    
    hlib_is_init = false;
}

//
// return true, if HLIBpro is initialised
//
bool
is_init ()
{
    return hlib_is_init;
}

}// namespace
