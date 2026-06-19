#ifndef __HPRO_TALGPARTSTRAT_HH
#define __HPRO_TALGPARTSTRAT_HH
//
// Project     : HLIBpro
// File        : TAlgPartStrat.hh
// Description : partitioning strategies for algebraic clustering
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/cluster/TGraph.hh"
#include "hpro/cluster/TNodeSet.hh"

namespace Hpro
{

//!
//! \ingroup  Cluster_Module
//! \class    TAlgPartStrat
//! \brief    Base class for partitioning strategies for algebraic clustering.
//!
class TAlgPartStrat
{
public:

    // dtor
    virtual ~TAlgPartStrat () {}
    
    //!
    //! compute graph bi-partitioning of \a graph and store result in \a left and \a right
    //!
    virtual void partition ( const TGraph &  graph,
                             TNodeSet &      left,
                             TNodeSet &      right ) const = 0;
};

//!
//! \ingroup  Cluster_Module
//! \class    TMETISAlgPartStrat
//! \brief    Graph partitioning using METIS
//!
class TMETISAlgPartStrat : public TAlgPartStrat
{
public:

    // ctor
    TMETISAlgPartStrat ( const bool use_random = CFG::Cluster::METIS_random );
    
    // dtor
    virtual ~TMETISAlgPartStrat () {}
    
    //! compute graph bi-partitioning of \a graph and store result in \a left and \a right
    virtual void partition ( const TGraph &  graph,
                             TNodeSet &      left,
                             TNodeSet &      right ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TScotchAlgPartStrat
//! \brief    Graph partitioning using Scotch
//!
class TScotchAlgPartStrat : public TAlgPartStrat
{
public:
    // dtor
    virtual ~TScotchAlgPartStrat () {}
    
    //! compute graph bi-partitioning of \a graph and store result in \a left and \a right
    virtual void partition ( const TGraph & graph,
                             TNodeSet &     left,
                             TNodeSet &     right ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TChacoAlgPartStrat
//! \brief    Graph partitioning using CHACO
//!
class TChacoAlgPartStrat : public TAlgPartStrat
{
public:
    // dtor
    virtual ~TChacoAlgPartStrat () {}
    
    //! compute graph bi-partitioning of \a graph and store result in \a left and \a right
    virtual void partition ( const TGraph & graph,
                             TNodeSet &     left,
                             TNodeSet &     right ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TMongooseAlgPartStrat
//! \brief    Graph partitioning using Mongoose
//!
class TMongooseAlgPartStrat : public TAlgPartStrat
{
public:

    // ctor
    TMongooseAlgPartStrat ();
    
    // dtor
    virtual ~TMongooseAlgPartStrat () {}
    
    //! compute graph bi-partitioning of \a graph and store result in \a left and \a right
    virtual void partition ( const TGraph &  graph,
                             TNodeSet &      left,
                             TNodeSet &      right ) const;
};

}// namespace Hpro

#endif  // __HPRO_TALGPARTSTRAT_HH
