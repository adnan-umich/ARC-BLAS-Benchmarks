#ifndef __HPRO_TVECTORIO_HH
#define __HPRO_TVECTORIO_HH
//
// Project     : HLIBpro
// File        : TVectorIO.hh
// Description : class for vector input/output
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/blas/Vector.hh"
#include "hpro/vector/TVector.hh"

namespace Hpro
{

//!
//! \{
//! \name Vector I/O
//! Functions for input/output of vectors.
//!

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TAutoVectorIO
//! \brief    Class for vector I/O with automatic
//!           file format detection
//!
class TAutoVectorIO 
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TAutoVectorIO () {}

    ~TAutoVectorIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write vector \a A to file \a filename
    template < typename value_t >
    void
    write ( const TVector< value_t > *  A,
            const std::string &         filename ) const;

    //! write vector \a A to file \a filename with optional
    //! vector name \a vecname (if file format has corresponding support)
    template < typename value_t >
    void
    write ( const TVector< value_t > *  A,
            const std::string &         filename,
            const std::string &         vecname ) const;

    //! read and return vector from file \a filename
    template < typename value_t >
    std::unique_ptr< TVector< value_t > >
    read  ( const std::string &  filename ) const;

    //! read and return vector from file \a filename with optional
    //! vector name \a vecname (if file format has corresponding support)
    template < typename value_t >
    std::unique_ptr< TVector< value_t > >
    read  ( const std::string &  filename,
            const std::string &  vecname ) const;
};

//!
//! \ingroup  IO_Module
//! \brief    Read vector from file with automatic file format detection.
//!
//!           Read vector from file \a filename while the format of the file
//!           is detected automatically. If the file format supports storage
//!           of several matrices (or other objects), only the first vector
//!           will be read.
//!
template < typename value_t >
std::unique_ptr< TVector< value_t > >
read_vector  ( const std::string &  filename );

//!
//! \ingroup  IO_Module
//! \brief    Read vector from file with automatic file format detection.
//!
//!           Read vector from file \a filename while the format of the file
//!           is detected automatically. If the file format supports named
//!           matrices, the vector with the name \a vecname will be read.
//!           Furthermore, if the file format supports storage of several matrices
//!           (or other objects), only the first vector will be read.
//!
template < typename value_t >
std::unique_ptr< TVector< value_t > >
read_vector  ( const std::string &  filename,
               const std::string &  vecname );

//!
//! \ingroup  IO_Module
//! \brief    Write vector to file with automatic choice of file format.
//!
//!           Write vector \a A to file \a filename while the format of the file
//!           is chosen automatically based upon the filename suffix.
//!           - hm                    : HLIBpro format
//!           - mat, m                : Matlab format
//!           - hb, rb, rua, rsa, psa : Harwell/Boeing format
//!           - rhs, sol              : SAMG format
//!           - mtx                   : VectorMarket format
//!
template < typename value_t >
void
write_vector  ( const TVector< value_t > *  A,
                const std::string &         filename );

//!
//! \ingroup  IO_Module
//! \brief    Write vector to file with automatic choice of file format.
//!
//!           Write vector \a A to file \a filename while the format of the file
//!           is chosen automatically based upon the filename suffix. If the file
//!           format supports named matrices, the vector will be stored under the
//!           name \a vecname.
//!
template < typename value_t >
void
write_vector  ( const TVector< value_t > *  A,
                const std::string &         filename,
                const std::string &         vecname );

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TSAMGVectorIO
//! \brief    Class for vector I/O in SAMG format
//!
class TSAMGVectorIO 
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TSAMGVectorIO () {}

    ~TSAMGVectorIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write vector \a v to file \a fname
    template < typename value_t >
    void
    write ( const TVector< value_t > *  v,
            const std::string &         fname ) const;

    //! read and return vector from file \a fname
    template < typename value_t >
    std::unique_ptr< TVector< value_t > >
    read  ( const std::string &  fname ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TMatlabVectorIO
//! \brief    Class for vector I/O in Matlab format
//!
class TMatlabVectorIO 
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TMatlabVectorIO () {}

    ~TMatlabVectorIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write vector \a v with name "v" to file \a fname
    template < typename value_t >
    void
    write ( const TVector< value_t > *  v,
            const std::string  &        fname ) const;
    
    //! write vector \a v with name \a vname to file \a fname
    template < typename value_t >
    void
    write ( const TVector< value_t > *  v,
            const std::string  &        fname,
            const std::string  &        vname ) const;

    //! read first vector in file \a fname
    template < typename value_t >
    std::unique_ptr< TVector< value_t > >
    read  ( const std::string & fname ) const;

    //! read and return vector named \a vname from file \a fname
    //! (if vname = "", return first vector)
    template < typename value_t >
    std::unique_ptr< TVector< value_t > >
    read  ( const std::string & fname,
            const std::string & vname ) const;

    // write BLAS vector \a v with name \a mname to file \a fname
    template < typename value_t >
    void
    write ( const BLAS::Vector< value_t > &  v,
            const std::string &              fname,
            const std::string &              mname ) const;

};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    THproVectorIO
//! \brief    Class for vector I/O in HLIB format
//!
class THproVectorIO 
{
protected:
    // use compression
    const bool  _compressed;
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    THproVectorIO ( const bool compressed = false )
            : _compressed(compressed)
    {}

    ~THproVectorIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write vector \a v to file \a fname
    template < typename value_t >
    void
    write ( const TVector< value_t > *  v,
            const std::string &         fname ) const;

    //! read and return vector from file \a fname
    template < typename value_t >
    std::unique_ptr< TVector< value_t > >
    read  ( const std::string &  fname ) const;
};

using THLibVectorIO = THproVectorIO;

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    THBVectorIO
//! \brief    Class for vector I/O in Harwell-Boeing and
//!           Rutherford-Boeing format
//!
class THBVectorIO 
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    THBVectorIO () {}

    ~THBVectorIO () {}

    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write vector \a v to file \a fname
    template < typename value_t >
    void
    write ( const TVector< value_t > *  v,
            const std::string &         fname ) const;

    //! read and return vector from file \a fname
    template < typename value_t >
    std::unique_ptr< TVector< value_t > >
    read  ( const std::string &  fname ) const;
};


///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TMMVectorIO
//! \brief    Class for vector I/O in MatrixMarket format
//!
class TMMVectorIO 
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TMMVectorIO () {}

    ~TMMVectorIO () {}

    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write vector \a v to file \a fname
    template < typename value_t >
    void
    write ( const TVector< value_t > *  v,
            const std::string &         fname ) const;

    //! read and return vector from file \a fname
    template < typename value_t >
    std::unique_ptr< TVector< value_t > >
    read  ( const std::string &  fname ) const;
};

//! \}

}// namespace

#endif  // __HPRO_TVECTORIO_HH
