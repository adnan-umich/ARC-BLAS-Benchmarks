#ifndef __HPRO_TCOORDINATE_HH
#define __HPRO_TCOORDINATE_HH
//
// Project     : HLIBpro
// File        : TCoordinate.hh
// Description : class for encapsulating coordinate data
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/base/TPoint.hh"

#include "hpro/cluster/TBBox.hh"

namespace Hpro
{

//! options for handling coordinate data
enum coord_data_t
{
    copy_coord_data,       // copy coordinate data
    reference_coord_data   // reference coordinate data
};
    
//!
//! \class  TCoordinate
//! \brief  stores coordinate information for indices
//!         - geometrical vectors are stored as double pointers
//!         - data can be externally defined, e.g. memory managed
//!           by user
//!
class TCoordinate
{
private:
    //! @cond
    
    // stores coordinates as array of vectors (arrays)
    std::vector< double * >  _coord;

    // additional bounding box data
    std::vector< double * >  _bbmin;
    std::vector< double * >  _bbmax;

    // dimension of the coordinates
    const uint               _dim;

    // periodicity of coordinates
    TPoint                   _period;
    
    // defines handling of given coordinate data
    const coord_data_t       _coord_data;

    //! @endcond
    
public:
    ////////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct coordinate object with \a acord holding the
    //! coordinate vectors, each of dimension \a adim; if \a ext
    //! is true, memory is handled externally
    TCoordinate ( const std::vector< double * > &  acoord,
                  const uint                       adim,
                  const coord_data_t               coord_data = copy_coord_data );
                  
    //! construct coordinate object with \a acord holding the
    //! coordinate vectors, each of dimension \a adim with additional
    //! bounding box information in \a abbmin and \a abbmax; if \a ext
    //! is true, memory is handled externally
    TCoordinate ( const std::vector< double * > &  acoord,
                  const uint                       adim,
                  const std::vector< double * > &  abbmin,
                  const std::vector< double * > &  abbmax,
                  const coord_data_t               coord_data = copy_coord_data );

    //! special versions for TPoint (always copy data)
    TCoordinate ( const std::vector< TPoint > &    acoord );
    
    TCoordinate ( const std::vector< TPoint > &    acoord,
                  const std::vector< TBBox > &     abbox );

    //! special versions for T2Point (always copy data)
    TCoordinate ( const std::vector< T2Point > &   acoord );
    
    TCoordinate ( const std::vector< T2Point > &   acoord,
                  const std::vector< TBBox > &     abbox );

    //! special versions for T3Point (always copy data)
    TCoordinate ( const std::vector< T3Point > &   acoord );
    
    TCoordinate ( const std::vector< T3Point > &   acoord,
                  const std::vector< TBBox > &     abbox );

    ~TCoordinate ();

    ////////////////////////////////////////////////////////
    //
    // give access to data
    //

    //! return number of stored coordinates
    size_t          ncoord   ()                 const { return _coord.size(); }

    //! return dimension of stored coordinates
    uint            dim      ()                 const { return _dim; }

    //! return coordinate \a i
    const double *  coord    ( const idx_t  i ) const { return _coord[i]; }

    //! return minimal coordinate of \a i'th bounding box
    const double *  bbmin    ( const idx_t  i ) const { return _bbmin[i]; }

    //! return maximal coordinate of \a i'th bounding box
    const double *  bbmax    ( const idx_t  i ) const { return _bbmax[i]; }

    //! return true if bounding box data is present
    bool            has_bbox () const
    {
        if (( _bbmin.size() == ncoord() ) && ( _bbmax.size() == ncoord() ))
            return true;
        else
            return false;
    }

    // directly access arrays
    const std::vector< double * > &  coord () const { return _coord; }
    const std::vector< double * > &  bbmin () const { return _bbmin; }
    const std::vector< double * > &  bbmax () const { return _bbmax; }

    //! return bounding box of coordinate set
    TBBox  bounding_box () const;
    
    //! set periodicity to \a p
    void set_periodicity ( const TPoint & p )
    {
        if ( ! ( p.dim() == _dim ) )
            HERROR( ERR_DIM, "(TCoordinate) set_periodicity",
                    "periodicity vector has wrong dimension" );
        
        _period = p;
    }

    //! return periodicity vector
    const TPoint & periodicity () const { return _period; }

    ////////////////////////////////////////////////////////
    //
    // misc.
    //

    //! return memory consumption
    size_t byte_size () const;

    DISABLE_COPY_OP( TCoordinate );
};

}// namespace

#endif  // __HPRO_TCOORDINATE_HH
