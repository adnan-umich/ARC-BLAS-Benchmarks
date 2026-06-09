#ifndef __HPRO_SCHEDULER_HH
#define __HPRO_SCHEDULER_HH
//
// Project     : HLIBpro
// File        : TListSched.hh
// Description : provides List scheduling algorithm
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/base/types.hh"

namespace Hpro
{

//
// MultiFit scheduling based upon bin-packing
//
class TMFitSched
{
protected:
    // depth of binary search
    const uint  _search_depth;
    
public:
    ///////////////////////////////////////////////
    //
    // constructor
    //

    TMFitSched ( const uint d = 10 ) : _search_depth(d) {}
    
    virtual ~TMFitSched () {}
    
    ///////////////////////////////////////////////
    //
    // scheduling algorithms
    //

    virtual bool schedule ( const uint                     p,
                            std::vector< int > &           sched,
                            const std::vector< double > &  costs ) const;

protected:
    // bin-packing "first fit decreasing" (FFD) algorithm
    uint ffd ( double                         cmax,
               uint                           p,
               const std::vector< double > &  costs,
               std::vector< int > &           sched ) const;

    DISABLE_COPY_OP( TMFitSched );
};

}// namespace Hpro

#endif  // __HPRO_SCHEDULER_HH
