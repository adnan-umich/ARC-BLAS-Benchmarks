#ifndef __HPRO_TGRIDVIS_HH
#define __HPRO_TGRIDVIS_HH
//
// Project     : HLIBpro
// File        : TGridVis.hh
// Description : grid visualisation classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/bem/TGrid.hh"
#include "hpro/bem/TFnSpace.hh"
#include "hpro/vector/TScalarVector.hh"

namespace Hpro
{
    
///////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TGridVis
//! \brief    Base class for grid visualisation.
//!
class TGridVisBase
{
private:
    //! @cond

    // minimal and maximal values for grid function output
    // - if equal (default), both values are chosen based on given function
    double       _min_fn_val, _max_fn_val;

    // name of colourmap (jet, hot, copper, hsv, bone, rainbow, coolwarm)
    std::string  _cmap;
    
    //! @endcond
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TGridVisBase ();

    virtual ~TGridVisBase () {}

    //! set value interval for grid function visualisation
    //! - colour is chosen based function value and position in interval
    //! - if actual function value is out of interval, value is clipped
    void         set_func_value_interval ( const double  minval,
                                           const double  maxval );

    //! access function value interval
    double       min_func_value          () const { return _min_fn_val; }
    double       max_func_value          () const { return _max_fn_val; }

    //! set colour map
    void         set_colourmap           ( const std::string &  cmap );

    //! return name of colourmap
    std::string  colourmap               () const { return _cmap; }
};

///////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    T2DGridVis
//! \brief    Base class for 2D grid visualisation (by projection).
//!

class T2DGridVis : public TGridVisBase
{
protected:
    //! @cond
    
    // viewing direction
    T3Point  _view;

    // apply lighting
    bool     _lighting;

    // draw bounding box
    bool     _draw_bbox;

    // draw coord. axis
    bool     _draw_axis;
    
    // draw also contour of triangle (instead of just area)
    bool     _draw_contour;
    
    //! @endcond
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    T2DGridVis ();
    T2DGridVis ( const T3Point view_dir );

    virtual ~T2DGridVis () {}

    //////////////////////////////////////
    //
    // options
    //

    //! set viewing direction (default: (1, 0, 0))
    T2DGridVis &  view_dir      ( const T3Point &  view );

    //! turn on/off lighting (default: off)
    T2DGridVis &  lighting      ( const bool       b );
    
    //! turn on/off printing of bounding box (default: off)
    T2DGridVis &  draw_bbox     ( const bool       b );

    //! turn on/off printing of coord. axis (default: off)
    T2DGridVis &  draw_axis     ( const bool       b );
    
    //! turn on/off drawing of triangle contours (default: on)
    T2DGridVis &  draw_contour  ( const bool       b );
};

///////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TPSGridVis
//! \brief    Class for grid visualisation in PostScript format.
//!

class TPSGridVis : public T2DGridVis
{
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TPSGridVis () {}
    TPSGridVis ( const T3Point aview_dir ) : T2DGridVis( aview_dir ) {}

    virtual ~TPSGridVis () {}

    //////////////////////////////////////
    //
    // print grid
    //
    
    //! print \a grid to file \a filename
    void print ( const TGrid *        grid,
                 const std::string &  filename ) const;

    //! print \a grid to file \a filename with colours according to
    //! values in vector \a vec which holds coefficients for function
    //! in \a fnspace
    template < typename fnspace_t,
               typename value_t >
    void print ( const TGrid *               grid,
                 const fnspace_t *           fnspace,
                 const TVector< value_t > *  vec,
                 const std::string &         filename ) const;

    template < typename fnspace_t,
               typename value_t >
    void print ( const TGrid &               grid,
                 const fnspace_t &           fnspace,
                 const TVector< value_t > &  vec,
                 const std::string &         filename ) const
    {
        static_assert( std::is_floating_point< value_t >::value, "only real valued types supported" );
        print( &grid, &fnspace, &vec, filename );
    }
};

//
//! \ingroup  IO_Module
//! \fn       print_ps
//! \brief    functional version of TPSGridVis
//
void
print_ps ( const TGrid *        grid,
           const std::string &  filename );

inline
void
print_ps ( const TGrid &        grid,
           const std::string &  filename )
{
    print_ps( &grid, filename );
}

///////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TVTKGridVis
//! \brief    Class for grid visualisation in VTK format.
//!
class TVTKGridVis : public TGridVisBase
{
private:
    //! @cond
    
    // draw normal direction
    bool  _draw_normal;
    
    //! @endcond
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TVTKGridVis ();

    virtual ~TVTKGridVis () {}

    //////////////////////////////////////
    //
    // options
    //
    
    //! turn on/off printing of normal direction (default: off)
    TVTKGridVis &  draw_normal  ( const bool  b );
    
    //////////////////////////////////////
    //
    // print grid
    //
    
    //! print \a grid to file \a filename
    void print ( const TGrid *        grid,
                 const std::string &  filename ) const;

    //! print \a grid to file \a filename with colours according to
    //! values in vector \a vec which holds coefficients for function
    //! in \a fnspace 
    template < typename fnspace_t,
               typename value_t >
    void print ( const TGrid *               grid,
                 const fnspace_t *           fnspace,
                 const TVector< value_t > *  vec,
                 const std::string &         filename ) const;

    template < typename fnspace_t,
               typename value_t >
    void print ( const TGrid &               grid,
                 const fnspace_t &           fnspace,
                 const TVector< value_t > &  vec,
                 const std::string &         filename ) const
    {
        static_assert( std::is_floating_point< value_t >::value, "only real valued types supported" );
        print( &grid, &fnspace, &vec, filename );
    }
};

//
//! \ingroup  IO_Module
//! \fn       print_vtk
//! \brief    functional version of TVTKGridVis
//
void
print_vtk ( const TGrid *        grid,
            const std::string &  filename );

inline
void
print_vtk ( const TGrid &        grid,
            const std::string &  filename )
{
    print_vtk( &grid, filename );
}


}// namespace

#endif  // __HPRO_TGRID_HH
