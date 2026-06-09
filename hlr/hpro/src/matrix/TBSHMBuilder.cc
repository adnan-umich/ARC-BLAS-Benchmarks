//
// Project     : HLIBpro
// File        : TBSHMBuilder.cc
// Description : class for building h-matrices out of bytestreams
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <iostream>

#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"
#include "hpro/matrix/TBlockMatrix.hh"
#include "hpro/matrix/THMatrix.hh"

#include "hpro/matrix/TBSHMBuilder.hh"

namespace Hpro
{

using std::make_unique;

/////////////////////////////////////////////////////////////////////////////
//
// local defines
//
/////////////////////////////////////////////////////////////////////////////

// define to do print detailed informations about computations
#define LOG( msg )   // std::cout << msg << std::endl << std::flush;

/////////////////////////////////////////////////////////////////////////////
//
// TBSHMBuilder
//
/////////////////////////////////////////////////////////////////////////////

//
// construct matrix out of given bytestream
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TBSHMBuilder::build ( TByteStream & bs ) const
{
    //
    // read type and build corresponding sub-matrix
    //

    uint  t;

    bs.get( & t, sizeof(uint) );

    std::unique_ptr< TMatrix< value_t > >  A( build_matrix< value_t >( t ) );

    // check if type of matrix is known
    if ( A.get() == nullptr )
        HERROR( ERR_NULL, "(TBSHMBuilder) build", "" );

    // use member function to build matrix
    A->build( bs );

    return A;
}

//
// return matrix corresponding to given type
//
template < typename value_t >
std::unique_ptr< TMatrix< value_t > >
TBSHMBuilder::build_matrix ( uint t ) const
{
    std::unique_ptr< TMatrix< value_t > >  M;

    if      ( t == TYPE_ID( TDenseMatrix ) ) M = make_unique< TDenseMatrix< value_t > >();
    else if ( t == TYPE_ID( TRkMatrix )    ) M = make_unique< TRkMatrix< value_t > >();
    else if ( t == TYPE_ID( TBlockMatrix ) ) M = make_unique< TBlockMatrix< value_t > >();
    else if ( t == TYPE_ID( THMatrix )     ) M = make_unique< THMatrix< value_t > >();
    else
        HERROR( ERR_MAT_TYPE, "(TBSHMBuilder) build_matrix", RTTI::id_to_type( t ) );

    LOG( "created " << M->typestr() );

    return M;
}

//
// explicit template instantiations
//

#define INST_BUILD( type )                              \
    template                                            \
    std::unique_ptr< TMatrix< type > >                  \
    TBSHMBuilder::build< type > ( TByteStream & ) const

INST_BUILD( float );
INST_BUILD( double );
INST_BUILD( std::complex< float > );
INST_BUILD( std::complex< double > );

#define INST_BUILD_MATRIX( type )                       \
    template                                            \
    std::unique_ptr< TMatrix< type > >                  \
    TBSHMBuilder::build_matrix< type > ( uint ) const

INST_BUILD_MATRIX( float );
INST_BUILD_MATRIX( double );
INST_BUILD_MATRIX( std::complex< float > );
INST_BUILD_MATRIX( std::complex< double > );

}// namespace Hpro
