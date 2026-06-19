#ifndef __HPRO_TBSHMBUILDER_HH
#define __HPRO_TBSHMBUILDER_HH
//
// Project     : HLIBpro
// File        : TBSHMBuilder.hh
// Description : class for building h-matrices out of bytestreams
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/TByteStream.hh"

#include "hpro/matrix/TMatrix.hh"

namespace Hpro
{

//
// takes a bytestream and reads in matrix
//
class TBSHMBuilder
{
public:
    //
    // constructor and destructor
    //

    TBSHMBuilder () {}

    virtual ~TBSHMBuilder () {}

    //
    // construct matrix out of given bytestream
    //

    template < typename value_t >
    std::unique_ptr< TMatrix< value_t > >  build ( TByteStream & bs ) const;

protected:
    //
    // special functions
    //

    // return matrix corresponding to given type
    template < typename value_t >
    std::unique_ptr< TMatrix< value_t > >  build_matrix ( uint t ) const;
};

}// namespace

#endif // __HPRO_TBSHMBUILDER_HH
