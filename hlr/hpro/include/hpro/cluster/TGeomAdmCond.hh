#ifndef __HPRO_TSTDADMCOND_HH
#define __HPRO_TSTDADMCOND_HH
//
// Project     : HLIBpro
// File        : TGeomAdmCond.hh
// Description : classes for geometry based admissibility-conditions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2024. All Rights Reserved.
//

#include <hpro/base/TPoint.hh>
#include <hpro/cluster/types.hh>

#include <hpro/cluster/TAdmCondition.hh>

namespace Hpro
{

//!
//! \ingroup  Cluster_Module
//! \class    TStrongGeomAdmCond
//! \brief    strong admissibility for FEM/BEM applications
//!           normal :  min( diam(τ), diam(σ) ) ≤ η·dist(τ,σ)
//!           use_max:  max( diam(τ), diam(σ) ) ≤ η·dist(τ,σ)
//!
class TStrongGeomAdmCond : public TAdmCondition
{
protected:
    //! parameter for ratio between diameter and distance
    const double       _eta;

    //! choose between min/max cluster diameters
    const diam_mode_t  _diam_mode;

    //! defines periodicity of coordinates
    const TPoint       _period;

public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct standard admissiblity 
    TStrongGeomAdmCond ( const double       eta       = 2.0,
                         const diam_mode_t  diam_mode = use_min_diam )
            : _eta(eta)
            , _diam_mode(diam_mode)
    {}

    //! construct standard admissiblity with periodic geometry
    TStrongGeomAdmCond ( const TPoint &     period,
                         const double       eta       = 2.0,
                         const diam_mode_t  diam_mode = use_min_diam )
            : _eta(eta)
            , _diam_mode(diam_mode)
            , _period(period)
    {}

    //! dtor
    virtual ~TStrongGeomAdmCond () {}

    ///////////////////////////////////////////
    //
    // check block-cluster if admissible
    //

    //! return true if cluster \a cl is admissible
    virtual bool is_adm ( const TBlockCluster * cl ) const;

    DISABLE_COPY_OP( TStrongGeomAdmCond );
};

// strong admissibility = standard admissibility
using TStdGeomAdmCond = TStrongGeomAdmCond;

//!
//! \ingroup  Cluster_Module
//! \class    TWeakGeomAdmCond
//! \brief    blocks are admissible if the number of dimensions with overlap of clusters is limited
//!
class TWeakGeomAdmCond : public TAdmCondition
{
private:
    // permitted number of dimensions with cluster overlap (default: 0)
    //   0 = shared vertex
    //   1 = shared edge
    //   2 = shared face  (full overlap in 2D)
    //   3 = shared cube  (full overlap in 3D)
    const uint    _noverlap;

    // permitted absolute cluster overlap per dimension (default: 0)
    const double  _hoverlap;
    
    // parameter for standard admissibility (default: 2)
    const double  _eta;

public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    TWeakGeomAdmCond ( const uint    noverlap = 0,
                       const double  hoverlap = 0.0,
                       const double  eta      = 2.0 )
            : TAdmCondition()
            , _noverlap( noverlap )
            , _hoverlap( hoverlap )
            , _eta( eta )
    {}

    virtual ~TWeakGeomAdmCond () {}

    ///////////////////////////////////////////
    //
    // check block-cluster if admissible
    //

    //! return true if \a cl is either weakly or strongly admissible
    virtual bool is_adm ( const TBlockCluster * cl ) const;
};

// strong admissibility = standard admissibility
using TWeakStdGeomAdmCond = TWeakGeomAdmCond;

//!
//! \ingroup  Cluster_Module
//! \class    THiLoFreqGeomAdmCond
//! \brief    Admissibility for high and low frequency regimes
//!
class THiLoFreqGeomAdmCond : public TAdmCondition
{
protected:
    // (modulus of) wave number 
    const double  _kappa;
    
    // parameter η of standard admissibility
    const double  _eta;

    // maximal number of wave lengths per cluster
    const double  _max_waves;

public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    THiLoFreqGeomAdmCond ( const std::complex< double >  kappa,
                           const double                  max_waves,
                           const double                  eta = 2.0 );

    virtual ~THiLoFreqGeomAdmCond () {}

    ///////////////////////////////////////////
    //
    // check block-cluster if admissible
    //

    //! return true if \a cl is admissible
    virtual bool is_adm ( const TBlockCluster * cl ) const;
};

}// namespace Hpro

#endif  // __HPRO_TGEOMADMCOND_HH
