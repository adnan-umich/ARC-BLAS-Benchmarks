#ifndef __HPRO_TADMCONDITION_HH
#define __HPRO_TADMCONDITION_HH
//
// Project     : HLIBpro
// File        : TAdmCondition.hh
// Description : baseclass for all admissibility conditions
//               for block-clusters
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/cluster/TBlockCluster.hh"

namespace Hpro
{

//!
//! \ingroup  Cluster_Module
//! \class    TAdmCondition
//! \brief    Defines basic interface for admissibility conditions.
//!
class TAdmCondition
{
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor
    TAdmCondition () {}

    //! dtor
    virtual ~TAdmCondition () {}

    ///////////////////////////////////////////
    //
    // check block-cluster if admissible
    //

    //! return true if cluster \a cl is admissible
    virtual bool is_adm ( const TBlockCluster * cl ) const = 0;
};

//!
//! \ingroup  Cluster_Module
//! \class    TOffDiagAdmCond
//! \brief    everything except diagonal is admissible
//!
class TOffDiagAdmCond : public TAdmCondition
{
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor
    TOffDiagAdmCond () {}

    //! dtor
    virtual ~TOffDiagAdmCond () {}

    ///////////////////////////////////////////
    //
    // check block-cluster if admissible
    //

    //! return true if cluster \a cl is admissible
    virtual bool is_adm ( const TBlockCluster * bc ) const
    {
        const auto  rowcl = bc->rowcl();
        const auto  colcl = bc->colcl();

        if (( rowcl == colcl ) || ( *rowcl == *colcl ))
            return false;

        return true;
    }
};

}// namespace Hpro

#endif  // __HPRO_TADMCONDITION_HH
