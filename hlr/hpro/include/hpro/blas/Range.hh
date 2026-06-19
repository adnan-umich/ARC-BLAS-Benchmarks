#ifndef __HPRO_BLAS_RANGE_HH
#define __HPRO_BLAS_RANGE_HH
//
// Project     : HLIBpro
// File        : Range.hh
// Description : implements class for index set with stride
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <ostream>

#include "hpro/base/error.hh"
#include "hpro/cluster/TIndexSet.hh"

namespace Hpro
{

namespace BLAS
{

//!
//! \ingroup  BLAS_Module
//! \class    Range
//! \brief    indexset with modified ctors
//!
class Range : public TIndexSet
{
public:
    //
    // constructors
    // - last < first => indexset = ∅
    //

    //! create index set { \a afirst ... \a last } with stride \a astride
    Range ( const idx_t   afirst,
            const idx_t   alast ) noexcept
            : TIndexSet( afirst, alast )
    {
    }

    //! create index set { 0, \a size-1 }
    Range ( const idx_t  size ) noexcept
            : TIndexSet( 0, size-1 )
    {}

    //! copy constructor for TIndexSet objects
    Range ( const TIndexSet &  is ) noexcept
            : TIndexSet( is )
    {}
            
    //
    // access data
    //

    //! return stride of index set
    size_t  stride () const noexcept { return 1; }

    //
    // output
    //

    //! stream output
    friend std::ostream & operator << ( std::ostream & os, const Range & r )
    {
        os << "[ " << r.first() << " , " << r.last() << " ]";
        return os;
    }

public:
    static Range all;
};

}// namespace BLAS

}// namespace Hpro

#endif  // __HPRO_BLAS_RANGE_HH
