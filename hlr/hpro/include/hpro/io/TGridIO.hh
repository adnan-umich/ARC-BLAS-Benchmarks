#ifndef __HPRO_TGRIDIO_HH
#define __HPRO_TGRIDIO_HH
//
// Project     : HLIBpro
// File        : TGridIO.hh
// Description : contains grid input/output classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/bem/TGrid.hh"

namespace Hpro
{

//!
//! \{
//! \name Grid I/O
//! Functions for input/output of grids.
//!

//////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TGridIO
//! \brief    Base class for reading grids.
//!
class TGridIO
{
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TGridIO () {}

    virtual ~TGridIO () {}

    //////////////////////////////////////
    //
    // read in a grid from a file
    //

    //! read and return grid from file \a filename
    virtual std::unique_ptr< TGrid >
    read   ( const std::string &  filename ) const = 0;

    //! write \a grid to file \a filename
    virtual void
    write  ( const TGrid *        grid,
             const std::string &  filename ) const = 0;
};

//////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TAutoGridIO
//! \brief    Class for grid I/O with file format detection
//!
class TAutoGridIO : public TGridIO
{
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TAutoGridIO () {}

    virtual ~TAutoGridIO () {}

    //////////////////////////////////////
    //
    // read in a grid from a file
    //

    //! read and return grid from file \a filename
    virtual std::unique_ptr< TGrid >
    read   ( const std::string &  filename ) const;

    //! write \a grid to file \a filename
    virtual void
    write  ( const TGrid *        grid,
             const std::string &  filename ) const;
};

//!
//! \ingroup  IO_Module
//! \brief    Read grid from file with automatic file format detection.
//!
//!           Read grid from file \a filename while the format of the file
//!           is detected automatically based on the file content and the
//!           file name.
//!
std::unique_ptr< TGrid >
read_grid   ( const std::string &  filename );

//!
//! \ingroup  IO_Module
//! \brief    Write \a grid to file with automatic file format detection.
//!
//!           Write \a grid to file \a filename while the format of the file
//!           is detected automatically based on the file name.
//!
void
write_grid  ( const TGrid *        grid,
              const std::string &  filename );

inline
void
write_grid  ( const TGrid &        grid,
              const std::string &  filename )
{
    write_grid( & grid, filename );
}

//////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    THLibGridIO
//! \brief    Class for grid I/O in HLIB format.
//!
class THproGridIO : public TGridIO
{
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    THproGridIO () {}

    virtual ~THproGridIO () {}

    //////////////////////////////////////
    //
    // read in a grid from a file
    //

    //! read and return grid from file \a filename
    virtual std::unique_ptr< TGrid >
    read   ( const std::string &  filename ) const;

    //! write \a grid to file \a filename
    virtual void
    write  ( const TGrid *        grid,
             const std::string &  filename ) const;
};

using THLibGridIO = THproGridIO;

//////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TPlyGridIO
//! \brief    Class for grid I/O in Ply format.
//!
class TPlyGridIO : public TGridIO
{
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TPlyGridIO () {}

    virtual ~TPlyGridIO () {}

    //////////////////////////////////////
    //
    // read in a grid from a file
    //

    //! read and return grid from file \a filename
    virtual std::unique_ptr< TGrid >
    read   ( const std::string &  filename ) const;

    //! write \a grid to file \a filename
    virtual void
    write  ( const TGrid *        grid,
             const std::string &  filename ) const;
};

//////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TSurfMeshGridIO
//! \brief    Class for grid I/O in Surface Mesh format.
//!
class TSurfMeshGridIO : public TGridIO
{
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TSurfMeshGridIO () {}

    virtual ~TSurfMeshGridIO () {}

    //////////////////////////////////////
    //
    // read in a grid from a file
    //

    //! read and return grid from file \a filename
    virtual std::unique_ptr< TGrid >
    read   ( const std::string &  filename ) const;

    //! write \a grid to file \a filename
    virtual void
    write  ( const TGrid *        grid,
             const std::string &  filename ) const;
};

//////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TGMSHGridIO
//! \brief    Class for grid I/O in GMSH format.
//!
class TGMSHGridIO : public TGridIO
{
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TGMSHGridIO () {}

    virtual ~TGMSHGridIO () {}

    //////////////////////////////////////
    //
    // read in a grid from a file
    //

    //! read and return grid from file \a filename
    virtual std::unique_ptr< TGrid >
    read   ( const std::string &  filename ) const;

    //! write \a grid to file \a filename
    virtual void
    write  ( const TGrid *        grid,
             const std::string &  filename ) const;
};

//! \}

}// namespace

#endif  // __HPRO_TGRIDIO_HH
