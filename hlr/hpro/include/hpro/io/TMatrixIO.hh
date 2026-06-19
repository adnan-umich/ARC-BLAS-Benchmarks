#ifndef __HPRO_TMATRIXIO_HH
#define __HPRO_TMATRIXIO_HH
//
// Project     : HLIBpro
// File        : TMatrixIO.hh
// Description : class for matrix input/output
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/blas/Matrix.hh"
#include "hpro/matrix/TMatrix.hh"

namespace Hpro
{

//!
//! \{
//! \name Matrix I/O
//! Functions for input/output of matrices.
//!

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TAutoMatrixIO
//! \brief    Class for matrix I/O with automatic
//!           file format detection
//!
class TAutoMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TAutoMatrixIO () {}

    virtual ~TAutoMatrixIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a filename
    template < typename value_t >
    void
    write ( const TMatrix< value_t > *  A,
            const std::string &         filename ) const;

    //! write matrix \a A to file \a filename with optional
    //! matrix name \a matname (if file format has corresponding support)
    template < typename value_t >
    void
    write ( const TMatrix< value_t > *  A,
            const std::string &         filename,
            const std::string &         matname ) const;

    //! read and return matrix from file \a filename
    template < typename value_t >
    std::unique_ptr< TMatrix< value_t > >
    read  ( const std::string &  name ) const;

    //! read and return matrix from file \a filename with optional
    //! matrix name \a matname (if file format has corresponding support)
    template < typename value_t >
    std::unique_ptr< TMatrix< value_t > >
    read  ( const std::string &  filename,
            const std::string &  matname ) const;
};

//!
//! \ingroup  IO_Module
//! \brief    Read matrix from file with automatic file format detection.
//!
//!           Read matrix from file \a filename while the format of the file
//!           is detected automatically. If the file format supports storage
//!           of several matrices (or other objects), only the first matrix
//!           will be read.
//!
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
read_matrix  ( const std::string &  filename );

//!
//! \ingroup  IO_Module
//! \brief    Read matrix from file with automatic file format detection.
//!
//!           Read matrix from file \a filename while the format of the file
//!           is detected automatically. If the file format supports named
//!           matrices, the matrix with the name \a matname will be read.
//!           Furthermore, if the file format supports storage of several matrices
//!           (or other objects), only the first matrix will be read.
//!
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
read_matrix  ( const std::string &  filename,
               const std::string &  matname );

//!
//! \ingroup  IO_Module
//! \brief    Write matrix to file with automatic choice of file format.
//!
//!           Write matrix \a A to file \a filename while the format of the file
//!           is chosen automatically based upon the filename suffix.
//!           - hm                    : HLIBpro format
//!           - mat, m                : Matlab format
//!           - hb, rb, rua, rsa, psa : Harwell/Boeing format
//!           - amg                   : SAMG format
//!           - mtx                   : MatrixMarket format
//!           - hdf, h5               : HDF5 format
//!
template < typename value_t >
void
write_matrix  ( const TMatrix< value_t > *  A,
                const std::string &         filename );

template < typename value_t >
void
write_matrix  ( const TMatrix< value_t > &  A,
                const std::string &         filename )
{
    write_matrix( &A, filename );
}

//!
//! \ingroup  IO_Module
//! \brief    Write matrix to file with automatic choice of file format.
//!
//!           Write matrix \a A to file \a filename while the format of the file
//!           is chosen automatically based upon the filename suffix. If the file
//!           format supports named matrices, the matrix will be stored under the
//!           name \a matname.
//!
template < typename value_t >
void
write_matrix  ( const TMatrix< value_t > *  A,
                const std::string &         filename,
                const std::string &         matname );

template < typename value_t >
void
write_matrix  ( const TMatrix< value_t > &  A,
                const std::string &         filename,
                const std::string &         matname )
{
    write_matrix( &A, filename, matname );
}

//!
//! \ingroup  IO_Module
//! \brief    Read linear operator from file
//!
//!           Read linear operator from file \a filename. Currently only the
//!           HLIBpro file format supports linear operators.
//!
template < typename value_t >
std::unique_ptr< TLinearOperator< value_t > >
read_linop  ( const std::string &  filename );

//!
//! \ingroup  IO_Module
//! \brief    Write linear operator to file.
//!
//!           Write linear operator \a A to file \a filename in HLIBpro format.
//!
template < typename value_t >
void
write_linop ( const TLinearOperator< value_t > *  A,
              const std::string &                 filename );

template < typename value_t >
void
write_linop ( const TLinearOperator< value_t > &  A,
              const std::string &                 filename )
{
    write_linop( &A, filename );
}

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TOctaveMatrixIO
//! \brief    Class for matrix I/O in octave format
//!
class TOctaveMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TOctaveMatrixIO ();

    virtual ~TOctaveMatrixIO ();
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a name
    template < typename value_t >
    void
    write ( const TMatrix< value_t > *  A,
            const std::string &         name ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TSAMGMatrixIO
//! \brief    Class for matrix I/O in SAMG format
//!
class TSAMGMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TSAMGMatrixIO () {}

    virtual ~TSAMGMatrixIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a name
    template < typename value_t >
    void
    write ( const TMatrix< value_t > *  A,
            const std::string &         name ) const;

    //! read and return matrix from file \a name
    template < typename value_t >
    std::unique_ptr< TMatrix< value_t > >
    read  ( const std::string &  name ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TMatlabMatrixIO
//! \brief    Class for matrix I/O in Matlab format
//!
class TMatlabMatrixIO
{
private:
    //! if true (default), any permutation of the matrices will be
    //! applied before writing them to the file
    const bool  _permute;
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor with internal permutation set to \a perm
    TMatlabMatrixIO ( const bool perm = CFG::IO::permute_save )
            : _permute( perm )
    {}

    virtual ~TMatlabMatrixIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A with name "M" to file \a fname
    template < typename value_t >
    void
    write ( const TMatrix< value_t > *  A,
            const std::string &         fname ) const;

    // write matrix \a A with name \a mname to file \a fname
    template < typename value_t >
    void
    write ( const TMatrix< value_t > *  A,
            const std::string &         fname,
            const std::string &         mname ) const;

    // write linear operator \a A with name \a mname to file \a fname
    template < typename value_t >
    void
    write ( const TLinearOperator< value_t > *  A,
            const std::string &                 fname,
            const std::string &                 mname ) const;

    // write BLAS matrix \a A with name \a mname to file \a fname
    template < typename value_t >
    void
    write ( const BLAS::Matrix< value_t > &  A,
            const std::string &              fname,
            const std::string &              mname ) const;

    //! read and return matrix from file \a fname (assuming only one entry available)
    template < typename value_t >
    std::unique_ptr< TMatrix< value_t > >
    read  ( const std::string &  fname ) const;

    //! read matrix with name \a mname from file \a fname
    template < typename value_t >
    std::unique_ptr< TMatrix< value_t > >
    read  ( const std::string &  fname,
            const std::string &  mname ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    THproMatrixIO
//! \brief    Class for matrix I/O in HLIB format
//!
class THproMatrixIO
{
protected:
    //! flag for using compression
    bool  _compressed;
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor with compression set to \a compressed
    THproMatrixIO ( const bool compressed = false )
            : _compressed(compressed)
    {}

    virtual ~THproMatrixIO () {}

    //! set internal compression usage to \a compressed
    void set_compression ( const bool  compressed ) { _compressed = compressed; }
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a fname
    template < typename value_t >
    void
    write ( const TMatrix< value_t > *  A,
            const std::string &         fname ) const;

    //! read and return matrix from file \a fname
    template < typename value_t >
    std::unique_ptr< TMatrix< value_t > >
    read  ( const std::string &  fname ) const;

    ///////////////////////////////////////////////
    //
    // interface for linear operators
    //

    //! write linear operator \a A to file \a fname
    template < typename value_t >
    void
    write_linop ( const TLinearOperator< value_t > *  A,
                  const std::string &                 fname ) const;

    //! read and return linear operator from file \a fname
    template < typename value_t >
    std::unique_ptr< TLinearOperator< value_t > >
    read_linop  ( const std::string &  fname ) const;
};

using THLibMatrixIO = THproMatrixIO;

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TPLTMGMatrixIO
//! \brief    Class for matrix I/O in PLTMG format
//!
class TPLTMGMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TPLTMGMatrixIO () {}

    virtual ~TPLTMGMatrixIO () {}

    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a fname
    template < typename value_t >
    void
    write ( const TMatrix< value_t > *  A,
            const std::string &         fname ) const;

    //! read and return matrix from file \a fname
    template < typename value_t >
    std::unique_ptr< TMatrix< value_t > >
    read  ( const std::string &  fname ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    THBMatrixIO
//! \brief    Class for matrix I/O in Harwell-Boeing and
//!           Rutherford-Boeing format
//!
class THBMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    THBMatrixIO () {}

    virtual ~THBMatrixIO () {}

    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a fname
    template < typename value_t >
    void
    write ( const TMatrix< value_t > *  A,
            const std::string &         fname ) const;

    //! read and return matrix from file \a fname
    template < typename value_t >
    std::unique_ptr< TMatrix< value_t > >
    read  ( const std::string &  fname ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TMMMatrixIO
//! \brief    Class for matrix I/O in MatrixMarket format.
//!
class TMMMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TMMMatrixIO () {}

    virtual ~TMMMatrixIO () {}

    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a fname
    template < typename value_t >
    void
    write ( const TMatrix< value_t > *  A,
            const std::string &         fname ) const;

    //! read and return matrix from file \a fname
    template < typename value_t >
    std::unique_ptr< TMatrix< value_t > >
    read  ( const std::string &  fname ) const;
};


///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    THDF5MatrixIO
//! \brief    Class for matrix I/O in HDF5 format.
//!
class THDF5MatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    THDF5MatrixIO () {}

    virtual ~THDF5MatrixIO () {}

    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a fname
    template < typename value_t >
    void
    write ( const TMatrix< value_t > *  A,
            const std::string &         fname ) const;

    // write matrix \a A with name \a mname to file \a fname
    template < typename value_t >
    void
    write ( const TMatrix< value_t > *  A,
            const std::string &         fname,
            const std::string &         mname ) const;

    // write matrix \a A with name \a mname to file \a fname
    template < typename value_t >
    void
    write ( const BLAS::Matrix< value_t > &  A,
            const std::string &              fname,
            const std::string &              mname ) const;
    
    //! read and return matrix from file \a fname
    template < typename value_t >
    std::unique_ptr< TMatrix< value_t > >
    read  ( const std::string &  fname ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TH2LibMatrixIO
//! \brief    Class for matrix I/O in format used by H2Lib.
//!
class TH2LibMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TH2LibMatrixIO () {}

    virtual ~TH2LibMatrixIO () {}

    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a fname
    template < typename value_t >
    void
    write ( const TMatrix< value_t > *  A,
            const std::string &         fname ) const;

    //! read and return matrix from file \a fname
    template < typename value_t >
    std::unique_ptr< TMatrix< value_t > >
    read  ( const std::string &  fname ) const;
};

//! \}

}// namespace Hpro

#endif  // __HPRO_TMATRIXIO_HH
