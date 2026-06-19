#ifndef __HPRO_TCLUSTERVIS_HH
#define __HPRO_TCLUSTERVIS_HH
//
// Project     : HLIBpro
// File        : TClusterVis.hh
// Description : class for cluster tree and block cluster tree I/O
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <string>

#include "hpro/cluster/TCluster.hh"
#include "hpro/cluster/TBlockCluster.hh"
#include "hpro/cluster/TCoordinate.hh"
#include "hpro/cluster/TPermutation.hh"

namespace Hpro
{

//!
//! \ingroup  IO_Module
//! \class    TClusterVis
//! \brief    base class for cluster tree visualisation
//!
class TClusterVis
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TClusterVis () {}

    virtual ~TClusterVis () {}
    
    ///////////////////////////////////////////////
    //
    // visualisation function
    //

    //! visualise cluster tree
    //! \param  cl         cluster tree to visualise
    //! \param  filename   name of file to write visual. to
    virtual void print ( const TCluster *     cl,
                         const std::string &  filename ) const = 0;
};

//!
//! \ingroup  IO_Module
//! \class    TBlockClusterVis
//! \brief    base class for block cluster tree visualisation
//!
class TBlockClusterVis
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TBlockClusterVis () {}

    virtual ~TBlockClusterVis () {}
    
    ///////////////////////////////////////////////
    //
    // visualisation function
    //

    //! visualise  block cluster tree
    //! \param  cl         block cluster tree to visualise
    //! \param  filename   name of file to write to
    virtual void print ( const TBlockCluster *  cl,
                         const std::string &    filename ) const = 0;
};

/////////////////////////////////////////////////////////////////
//
// 2D visualisation (PostScript, SVG)
//
/////////////////////////////////////////////////////////////////

class T2DPrinter;

//
//! \ingroup  IO_Module
//! \class    T2DClusterVis
//! \brief    base class for cluster tree visualisation in 2D
//
class T2DClusterVis : public TClusterVis
{
public:
    //! @cond
    
    //
    // visualisation options (see below)
    //
    struct option_t
    {
        bool  tree;
        bool  node_procs;
        bool  id;

        option_t ();
    };

private:
    // options for visualisation
    option_t  _opt;
    
    //! @endcond
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    T2DClusterVis () {}

    virtual ~T2DClusterVis () {}
    
    ///////////////////////////////////////////////
    //
    // option management
    //
    
    //! turn on/off drawing of tree
    T2DClusterVis &  tree        ( const bool  b );
    
    //! turn on/off drawing of per-node processors
    T2DClusterVis &  node_procs  ( const bool  b );
    
    //! turn on/off drawing of node IDs
    T2DClusterVis &  id          ( const bool  b );
    
    ///////////////////////////////////////////////
    //
    // visualisation functions
    //

    virtual void print ( const TCluster *     c,
                         const std::string &  filename ) const;

protected:
    // return printer object for given format
    virtual T2DPrinter *  get_printer ( const double         width,
                                        const double         height,
                                        const std::string &  filename ) const = 0;
};

//
//! \ingroup  IO_Module
//! \class    TPSClusterVis
//! \brief    class for cluster tree visualisation in PostScript format
//
class TPSClusterVis : public T2DClusterVis
{
protected:
    // return printer object for given format
    virtual T2DPrinter *  get_printer ( const double         width,
                                        const double         height,
                                        const std::string &  filename ) const;
};

//
//! \ingroup  IO_Module
//! \fn       print_ps
//! \brief    functional version of TPSClusterVis
//
void
print_ps ( const TCluster *     cl,
           const std::string &  filename );

//
//! \ingroup  IO_Module
//! \class    T2DBlockClusterVis
//! \brief    base class for block cluster tree visualisation in 2D
//
class T2DBlockClusterVis : public TBlockClusterVis
{
public:
    //! @cond
    
    //
    // visualisation options (see below)
    //
    struct option_t
    {
        bool    border;
        bool    background;
        bool    all_nodes;
        bool    loc_is;
        bool    id;
        bool    inner_id;
        bool    procs;
        bool    single_proc;
        uint    pid;
        bool    legend;
        double  max_size_ratio;
        bool    adaptive_lw;
        double  base_line_width;
        
        option_t ();
    };

private:
    // options for visualisation
    option_t   _opt;
    
    //! @endcond
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! default ctor
    T2DBlockClusterVis () {}

    virtual ~T2DBlockClusterVis () {}
    
    ///////////////////////////////////////////////
    //
    // option managment
    //

    //! turn on/off printing of border around nodes (default: on)
    T2DBlockClusterVis &  border      ( const bool  b );
        
    //! turn on/off printing of background of nodes (default: on)
    T2DBlockClusterVis &  background  ( const bool  b );
        
    //! turn on/off showing of all nodes in tree, even inner nodes (default: on)
    T2DBlockClusterVis &  all_nodes   ( const bool  b );
        
    //! turn on/off printing local index sets (default: off)
    T2DBlockClusterVis &  loc_is      ( const bool  b );

    //! turn on/off printing local cluster ID (default: off)
    T2DBlockClusterVis &  id          ( const bool  b );

    //! turn on/off printing cluster ID for inner nodes (default: off)
    T2DBlockClusterVis &  inner_id    ( const bool  b );

    //! turn on/off printing tree for one processor (default: off)
    T2DBlockClusterVis &  procs       ( const bool  b );
        
    //! turn on/off printing tree for one processor (default: off)
    T2DBlockClusterVis &  single_proc ( const bool  b );
        
    //! set processor to restrict visualisation to (default: none);
    //! also activates "single_proc"
    T2DBlockClusterVis &  pid         ( const uint  p );

    //! turn on/off printing of legend (default: off)
    T2DBlockClusterVis &  legend      ( const bool  b );

    //! set maximal allowed ratio of block size compared to largest block size
    //! (default: 1000); if ratio is exceeded, block will not be printed to 
    //! limit file size for very large block cluster trees
    T2DBlockClusterVis &  max_size_ratio   ( const double  r );

    //! switch on/off adaptive line width, e.g. multiply lw by distance from leaf;
    //! if off, constant lw = base lw is used (default: off)
    T2DBlockClusterVis &  adaptive_lw      ( const bool    b );
    
    //! set base line width for either constant or adaptive lw (default: 1.0)
    T2DBlockClusterVis &  base_line_width  ( const double  lw );
    
    ///////////////////////////////////////////////
    //
    // visualisation functions
    //

    virtual void print ( const TBlockCluster *  c,
                         const std::string &    filename ) const;

protected:
    // return printer object for given format
    virtual T2DPrinter *  get_printer ( const double         width,
                                        const double         height,
                                        const std::string &  filename ) const = 0;
};

//
//! \ingroup  IO_Module
//! \class    TPSBlockClusterVis
//! \brief    class for block cluster tree visualisation in PostScript format
//
class TPSBlockClusterVis : public T2DBlockClusterVis
{
protected:
    // return printer object for given format
    virtual T2DPrinter *  get_printer ( const double         width,
                                        const double         height,
                                        const std::string &  filename ) const;
};

//
//! \ingroup  IO_Module
//! \fn       print_ps
//! \brief    functional version of TPSBlockClusterVis
//
void
print_ps ( const TBlockCluster *  cl,
           const std::string &    filename );

/////////////////////////////////////////////////////////////////
//
// VTK visualisation
//
/////////////////////////////////////////////////////////////////

//
//! \ingroup  IO_Module
//! \class    TVTKBlockClusterVis
//! \brief    class for block cluster tree visualisation in VTK
//
class TVTKClusterVis : public TClusterVis
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! default ctor
    TVTKClusterVis () {}

    virtual ~TVTKClusterVis () {}
    
    ///////////////////////////////////////////////
    //
    // visualisation functions
    //

    virtual void print ( const TCluster *     cl,
                         const std::string &  filename ) const;

    virtual void print ( const TCluster *      cl,
                         const TCoordinate &   coord,
                         const TPermutation &  perm_i2e,
                         const std::string &   filename ) const;
};

//
//! \ingroup  IO_Module
//! \fn       print_vtk
//! \brief    functional version of TVTKClusterVis
//
void
print_vtk ( const TCluster *      cl,
            const TCoordinate &   coord,
            const TPermutation &  perm_i2e,
            const std::string &   filename );

//
//! \ingroup  IO_Module
//! \class    TVTKBlockClusterVis
//! \brief    class for block cluster tree visualisation in VTK
//
class TVTKBlockClusterVis : public TBlockClusterVis
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! default ctor
    TVTKBlockClusterVis () {}

    virtual ~TVTKBlockClusterVis () {}
    
    ///////////////////////////////////////////////
    //
    // visualisation functions
    //

    virtual void print ( const TBlockCluster *  c,
                         const std::string &    filename ) const;
};

//
//! \ingroup  IO_Module
//! \fn       print_vtk
//! \brief    functional version of TVTKBlockClusterVis
//
void
print_vtk ( const TBlockCluster *  cl,
            const std::string &    filename );

/////////////////////////////////////////////////////////////////
//
// GraphViz visualisation
//
/////////////////////////////////////////////////////////////////

//
//! \ingroup  IO_Module
//! \class    TGVClusterVis
//! \brief    cluster tree visualisation in GraphViz format
//
class TGVClusterVis : public TClusterVis
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TGVClusterVis () {}

    virtual ~TGVClusterVis () {}
    
    ///////////////////////////////////////////////
    //
    // visualisation functions
    //

    virtual void print ( const TCluster *     c,
                         const std::string &  filename ) const;
};

//
//! \ingroup  IO_Module
//! \fn       print_gv
//! \brief    functional version of TGVClusterVis
//
void
print_gv ( const TCluster *      cl,
           const std::string &   filename );

//
//! \ingroup  IO_Module
//! \class    TGVBlockClusterVis
//! \brief    block cluster tree visualisation in GraphViz format
//
class TGVBlockClusterVis : public TBlockClusterVis
{
private:
    // processor to restrict visualisation to
    const int   _pid;
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TGVBlockClusterVis ( const int  pid = -1 )
            : _pid( pid )
    {}

    virtual ~TGVBlockClusterVis () {}
    
    ///////////////////////////////////////////////
    //
    // visualisation functions
    //

    virtual void print ( const TBlockCluster *  c,
                         const std::string &    filename ) const;
};

//
//! \ingroup  IO_Module
//! \fn       print_gv
//! \brief    functional version of TGVClusterVis
//
void
print_gv ( const TBlockCluster *  cl,
           const std::string &    filename );

}// namespace Hpro

#endif  // __TCLUSTERVIS_HH
