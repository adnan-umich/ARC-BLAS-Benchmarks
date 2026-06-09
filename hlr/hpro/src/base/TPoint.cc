//
// Project     : HLIBpro
// File        : TPoint.cc
// Description : class for a n-dimensional vector
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <sstream>

#include "hpro/base/TPoint.hh"

namespace Hpro
{

std::string
TPoint::to_string () const
{
    std::ostringstream  tmp;
    
    tmp << "( ";
        
    for ( uint i = 0; i < dim(); i++ )
    {
        if ( i > 0 )
            tmp << ", ";
            
        tmp << (*this)[i];
    }// for

    tmp << " )";
    
    return tmp.str();
}

std::string
T3Point::to_string () const
{
    std::ostringstream  tmp;

    tmp << "( " << _data[0] << ", " << _data[1] << ", " << _data[2] << " )";

    return tmp.str();
}

std::string
T2Point::to_string () const
{
    std::ostringstream  tmp;

    tmp << "( " << _data[0] << ", " << _data[1] << " )";

    return tmp.str();
}

}// namespace Hpro
