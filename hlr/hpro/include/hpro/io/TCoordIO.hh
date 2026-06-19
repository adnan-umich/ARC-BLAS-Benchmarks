#ifndef __HPRO_TCOORDIO_HH
#define __HPRO_TCOORDIO_HH
//
// Project     : HLIBpro
// File        : TCoordIO.hh
// Description : classes for coordinate input/output
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <string>

#include "hpro/cluster/TCoordinate.hh"

namespace Hpro
{

//!
//! \{
//! \name Coordinate I/O
//! Functions for input/output of coordinates.
//!

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TCoordIO
//! \brief    Base class for coordinate I/O defining interface
//!
class TCoordIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TCoordIO () {}

    virtual ~TCoordIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write coordinates \a coo to file \a filename
    virtual void
    write ( const TCoordinate *  coo,
            const std::string &  filename ) const;
    
    //! return and return coordinates from file \a filename
    virtual std::unique_ptr< TCoordinate >
    read  ( const std::string &  filename ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TAutoCoordIO
//! \brief    Class for coordinate I/O with file format detection
//!
class TAutoCoordIO : public TCoordIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TAutoCoordIO () {}

    virtual ~TAutoCoordIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write coordinates \a coo to file \a filename
    virtual void
    write ( const TCoordinate *  coo,
            const std::string &  filename ) const;

    //! return and return coordinates from file \a filename
    virtual std::unique_ptr< TCoordinate >
    read  ( const std::string &  filename ) const;
};

//!
//! \ingroup  IO_Module
//! \brief    Read coord from file with automatic file format detection.
//!
//!           Read coord from file \a filename while the format of the file
//!           is detected automatically.
//!
std::unique_ptr< TCoordinate >
read_coord  ( const char *  filename );

//!
//! \ingroup  IO_Module
//! \brief    Write coord to file with automatic choice of file format.
//!
//!           Write coord \a coord to file \a filename with automatic choice
//!           of file format based on filename suffix:
//!           - hm     : HLIBpro format
//!           - mat, m : Matlab format
//!           - coo    : SAMG format
//!
void
write_coord  ( const TCoordinate *  coord,
               const char *         filename );

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TSAMGCoordIO
//! \brief    Class for coordinate I/O in SAMG format
//!
class TSAMGCoordIO : public TCoordIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TSAMGCoordIO () {}

    virtual ~TSAMGCoordIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write coordinates \a coo to file \a filename
    virtual void
    write ( const TCoordinate *  coo,
                                  const std::string &  filename ) const;

    //! return and return coordinates from file \a filename
    virtual std::unique_ptr< TCoordinate >
    read  ( const std::string &  filename ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TMatlabCoordIO
//! \brief    Class for coordinate I/O in Matlab format
//!
class TMatlabCoordIO : public TCoordIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TMatlabCoordIO () {}

    virtual ~TMatlabCoordIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write coordinates \a coo with name "coo" to file \a filename
    virtual void
    write ( const TCoordinate *  coo,
            const std::string &  filename ) const;

    //! write coordinates \a coo with name \a cooname to file \a filename
    virtual void
    write ( const TCoordinate *  coo,
            const std::string &  filename,
            const std::string &  cooname ) const;

    //! read and return first coordinates from file \a filename
    virtual std::unique_ptr< TCoordinate >
    read  ( const std::string &  filename ) const;

    //! read and return coordinates with name \a cooname from file \a filename
    virtual std::unique_ptr< TCoordinate >
    read  ( const std::string &  filename,
            const std::string &  cooname ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    THproCoordIO
//! \brief    Class for coordinate I/O in HLIB format
//!
class THproCoordIO : public TCoordIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    THproCoordIO () {}

    virtual ~THproCoordIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write coordinates \a coo to file \a filename
    virtual void
    write ( const TCoordinate *  coo,
            const std::string &  filename ) const;

    //! return and return coordinates from file \a filename
    virtual std::unique_ptr< TCoordinate >
    read  ( const std::string &  filename ) const;
};

using THLibCoordIO = THproCoordIO;

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TMMCoordIO
//! \brief    Class for coordinate I/O in MatrixMarket format
//!
class TMMCoordIO : public TCoordIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TMMCoordIO () {}

    virtual ~TMMCoordIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write coordinates \a coo to file \a filename
    virtual void
    write ( const TCoordinate *  coo,
            const std::string &  filename ) const;

    //! return and return coordinates from file \a filename
    virtual std::unique_ptr< TCoordinate >
    read  ( const std::string &  filename ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TPLTMGCoordIO
//! \brief    Class for coordinate I/O in PLTMG format
//!
class TPLTMGCoordIO : public TCoordIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TPLTMGCoordIO () {}

    virtual ~TPLTMGCoordIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write coordinates \a coo to file \a filename
    virtual void
    write ( const TCoordinate *  coo,
            const std::string &  filename ) const;

    //! return and return coordinates from file \a filename
    virtual std::unique_ptr< TCoordinate >
    read  ( const std::string &  filename ) const;
};

//////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TGMSHCoordIO
//! \brief    Class for coordinate I/O in GMSH format.
//!
class TGMSHCoordIO : public TCoordIO
{
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TGMSHCoordIO () {}

    virtual ~TGMSHCoordIO () {}

    //////////////////////////////////////
    //
    // read in a grid from a file
    //

    //! write coordinates \a coo to file \a filename
    virtual void
    write ( const TCoordinate *  coo,
            const std::string &  filename ) const;

    //! return and return coordinates from file \a filename
    virtual std::unique_ptr< TCoordinate >
    read  ( const std::string &  filename ) const;
};

//! \}

}// namespace Hpro

#endif  // __HPRO_TCOORDIO_HH
