#ifndef __HPRO_COLOURMAP_HH
#define __HPRO_COLOURMAP_HH
//
// Project     : HLIBpro
// File        : TColourMap.hh
// Description : provides various colourmaps
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include <hpro/base/types.hh>

#include "colour.hh"

namespace Hpro
{

//!
//! \class  TColourMap
//! \brief  Baseclass for a colourmap, defining interface.
//!
class TColourMap
{
protected:
    // vector containing RGB colours of map
    std::vector< colour_t >  _map;

public:
    //
    // constructor and destructor
    //

    TColourMap () {}

    virtual ~TColourMap () {}

    //
    // access map entries
    //

    //! return i'th entry in map
    const colour_t   ientry ( const size_t  i ) const
    {
        return _map[i];
    }

    //! map f ∈ [0,1] to i ∈ [0,_map.size()-1] and return entry
    const colour_t   fentry ( const double  f ) const
    {
        return _map[ size_t( std::min( 1.0, std::max( 0.0, f ) ) * (_map.size()-1) ) ];
    }

    //
    // initialisation
    //

    //! initialise values for a map of n entries
    virtual void init ( const size_t  n ) = 0;
};

//!
//! \fn  colourmap
//! \brief  return colourmap by name, e.g. jet, hot, coolwarm, etc.
//!
TColourMap *
colourmap ( const std::string &  name,
            const size_t         n = 1000 );


#define  DEFINE_CMAP( name )                    \
    class name : public TColourMap              \
    {                                           \
    public:                                     \
        name ( const size_t  n )                \
        {                                       \
            init( n );                          \
        }                                       \
        virtual void init ( const size_t n );   \
    }

//!
//! \class  TTangoColourMap
//! \brief  colour map from tango desktop project
//!
DEFINE_CMAP( TTangoColourMap );

//!
//! \class  TJetColourMap
//! \brief  colourmap from dark blue through blue, cyan, green, yellow, red to dark red
//!
DEFINE_CMAP( TJetColourMap );


//!
//! \class  THotColourMap
//! \brief  colourmap from black through dark red, red, orange, yellow and white
//!
DEFINE_CMAP( THotColourMap );

//!
//! \class  TCopperColourMap
//! \Brief  Copper colourmap from black to light copper tone
//!
DEFINE_CMAP( TCopperColourMap );

//!
//! \class  THSVColourMap
//! \brief  colourmap using HSV colour space
//!
DEFINE_CMAP( THSVColourMap );

//!
//! \class  TBoneColourMap
//! \brief  grey colourmap with light blue tone
//!
DEFINE_CMAP( TBoneColourMap );

//!
//! \class  TRainbowColourMap
//! \brief  colourmap from red through orange, yellow, green, blue to violet
//!
DEFINE_CMAP( TRainbowColourMap );

//!
//! \class  TCoolWarmColourMap
//! \brief  diverging colourmap from cool over white to warm colours
//!         (see http://www.sandia.gov/~kmorel/documents/ColorMaps/)
//!
DEFINE_CMAP( TCoolWarmColourMap );

//!
//! \class  TMagmaColourMap
//! \brief  perceptually uniform colourmap
//!         (definition from https://github.com/BIDS/colormap)
//!
DEFINE_CMAP( TMagmaColourMap );

//!
//! \class  TInfernoColourMap
//! \brief  perceptually uniform colourmap
//!         (definition from https://github.com/BIDS/colormap)
//!
DEFINE_CMAP( TInfernoColourMap );

//!
//! \class  TPlasmaColourMap
//! \brief  perceptually uniform colourmap
//!         (definition from https://github.com/BIDS/colormap)
//!
DEFINE_CMAP( TPlasmaColourMap );

//!
//! \class  TViridisColourMap
//! \brief  perceptually uniform colourmap
//!         (definition from https://github.com/BIDS/colormap)
//!
DEFINE_CMAP( TViridisColourMap );

//!
//! \class  TParulaColourMap
//! \brief  colour map from Matlab
//!
DEFINE_CMAP( TParulaColourMap );

//!
//! \class  TRandomColourMap
//! \brief  randomise given base colourmap
//!
class TRandomColourMap : public TColourMap
{
private:
    // base colour map
    TColourMap &  _base_cmap;
    
public:
    //
    // constructor and destructor
    //

    TRandomColourMap ( TColourMap &  base_cmap,
                       const size_t  n )
            : _base_cmap( base_cmap )
    {
        init( n );
    }

    virtual ~TRandomColourMap () {}

    //
    // initialise values for a map of n entries
    //

    virtual void init ( const size_t n );
};

//!
//! \class  TShuffleColourMap
//! \brief  shuffle given base colourmap in cyclic way
//!
class TShuffleColourMap : public TColourMap
{
private:
    // base colour map
    TColourMap &  _base_cmap;
    
public:
    //
    // constructor and destructor
    //

    TShuffleColourMap ( TColourMap &  base_cmap,
                        const size_t  n )
            : _base_cmap( base_cmap )
    {
        init( n );
    }

    virtual ~TShuffleColourMap () {}

    //
    // initialise values for a map of n entries
    //

    virtual void init ( const size_t n );
};

}// namespace Hpro

#endif  // __HPRO_COLOURMAP_HH
