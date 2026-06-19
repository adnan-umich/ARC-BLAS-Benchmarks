#ifndef __HPRO_INIT_HH
#define __HPRO_INIT_HH
//
// Project     : HLIBpro
// File        : init.hh
// Description : initialisation and finalisation
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

namespace Hpro
{

///////////////////////////////////////////////////////////////////
//
// initialisation and finalisation of HLIBpro
//
///////////////////////////////////////////////////////////////////

//!
//! global initialisation routine for HLIBpro
//!
void INIT ();

//!
//! finalisation routing for HLIBpro
//!
void DONE ();

//!
//! return true, if HLIBpro is initialised
//!
bool is_init ();

}// namespace

#endif // __INIT_HH
