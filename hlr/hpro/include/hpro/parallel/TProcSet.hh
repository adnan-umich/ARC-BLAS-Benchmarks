#ifndef __HPRO_TPROCSET_HH
#define __HPRO_TPROCSET_HH
//
// Project     : HLIBpro
// File        : TProcSet.hh
// Description : class for representing processor sets
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>
#include <string>

#include "hpro/base/types.hh"
#include "hpro/base/String.hh"
#include "hpro/base/System.hh"
#include "hpro/base/TByteStream.hh"

namespace Hpro
{

//////////////////////////////////////////////////////////////////////
//!
//! \class TProcSet
//! \brief describes a processor set of continuously numbered
//!        processors
//!
//////////////////////////////////////////////////////////////////////

class TProcSet
{
public:
    // STL iterator
    struct iterator
    {
    private:
        idx_t  _pos;

    public:
        iterator () = delete;
        iterator ( const idx_t  pos ) noexcept : _pos(pos) {}

        iterator  operator ++ ()    noexcept { return iterator( ++_pos ); }
        iterator  operator ++ (int) noexcept { return iterator( _pos++ ); }

        iterator  operator -- ()    noexcept { return iterator( --_pos ); }
        iterator  operator -- (int) noexcept { return iterator( _pos-- ); }

        idx_t     operator *  () const noexcept { return _pos; }

        bool      operator == ( const iterator &  it ) const noexcept { return _pos == it._pos; }
        bool      operator != ( const iterator &  it ) const noexcept { return _pos != it._pos; }
    };
    
    using  const_iterator = iterator;
    
private:
    //!< first and last processors in set
    uint  _first, _last;

public:
    ///////////////////////////////////
    //
    // constructor and destructor
    //

    //! constructs set of size \a ps_size
    explicit TProcSet ( const uint  ps_size = 0 ) noexcept
            : _first(0)
            , _last(ps_size-1)
    {
        _last = std::max( _first, _last );
    }

    //! constructs set from \a ps_first to \a ps_last
    TProcSet ( const uint  ps_first,
               const uint  ps_last ) noexcept
            : _first(ps_first)
            , _last(ps_last)
    {
        _last = std::max( _first, _last );
    }

    //! copy constructor
    TProcSet ( const TProcSet & ps ) noexcept
    {
        *this = ps;
    }

    //
    // operations over processor sets
    //

    //! return ID of first processor in set
    uint first () const noexcept { return _first; }

    //! return ID of last processor in set
    uint last  () const noexcept { return _last;  }

    //! return i'th subset in an n-partition of the local set
    /*!
      Split local processor set into n subsets and return i'th subset.
      \param i  subset to return
      \param n  number of subsets to split into
     */
    TProcSet subset ( const uint  i,
                      const uint  n = 2 ) const noexcept
    {
        const uint  l = size() / n;

        if ( i == n-1 )
            return TProcSet( _first + i*l, _last );
        else
            return TProcSet( _first + i*l, _first + (i+1)*l - 1 );
    }

    //! split processor set in n subsets (of approx. equal size)
    /*!
      Split local processor set into n subsets of approximately the same size
      and return the partition in psets
      \param n      number of subsets to generate
      \param psets  array to store the subsets
     */
    void split ( const uint                 n,
                 std::vector< TProcSet > &  psets ) const
    {
        psets.resize( n );
        
        for ( uint i = 0; i < n; i++ )
            psets[i] = subset( i, n );
    }

    //! return union with given procsessor set
    /*!
      Construct a processor set defined by minimum of first processor
      and maximum of last processor in both sets
    */
    TProcSet  join ( const TProcSet & ps ) const noexcept
    {
        return TProcSet( std::min( first(), ps.first() ),
                         std::max( last(),  ps.last()  ) );
    }
    
    //! return true if given proc is in set
    bool  is_in    ( const uint p ) const noexcept { return ((p >= _first) && (p <= _last)); }
    bool  contains ( const uint p ) const noexcept { return is_in( p ); }

    //! return size of set
    uint  size    () const noexcept { return _last - _first + 1; }

    //! return true if processor set is empty
    bool  empty   () const noexcept { return size() == 0; }

    //! return id of master processor
    uint  master  () const noexcept { return _first; }

    //
    // STL interface
    //

    iterator        begin  () const noexcept { return iterator( _first ); }
    iterator        end    () const noexcept { return iterator( _last+1 ); }
        
    const_iterator  cbegin () const noexcept { return const_iterator( _first ); }
    const_iterator  cend   () const noexcept { return const_iterator( _last+1 ); }
        
    //
    // misc.
    //

    //! copy operator
    TProcSet &  operator =  ( const TProcSet & ps ) noexcept
    {
        _first = ps._first;
        _last  = ps._last;
        return *this;
    }
    
    //! equality operators
    bool  operator ==  ( const TProcSet & ps ) const noexcept { return (( _first == ps._first ) && ( _last == ps._last )); }

    //! inequality operators
    bool  operator !=  ( const TProcSet & ps ) const noexcept { return ! ( *this == ps ); }
    
    //! return size in bytes used by this object
    size_t  byte_size  () const noexcept { return 2*sizeof(uint); }
    
    //! return string representation
    std::string  to_string  () const
    {
        if ( _first != _last ) return Hpro::to_string( "[%d:%d]", _first, _last );
        else                   return Hpro::to_string( "[%d]", _first );
    }

    //
    // serialisation
    //

    //! read data from stream
    void    read     ( TByteStream & s );

    //! write data to stream
    void    write    ( TByteStream & s ) const;

    //! returns size of object in bytestream
    size_t  bs_size  () const;
};

//////////////////////////////////////////////////////////
//
// functional ctors
//

inline
TProcSet
ps ( const uint  size ) noexcept
{
    return TProcSet( size );
}

inline
TProcSet
ps_single ( const uint  pid ) noexcept
{
    return TProcSet( pid, pid );
}

inline
TProcSet
ps ( const uint  first,
     const uint  last ) noexcept
{
    return TProcSet( first, last );
}

//
// stream output
//
inline std::ostream &
operator << ( std::ostream &    os,
              const TProcSet &  ps )
{
    return os << ps.to_string();
}

//
// processor set representing invalid state
//
extern const TProcSet  PROCSET_INVALID;

//
// external wrappers for member functions
//
inline
std::vector< TProcSet >
split ( const TProcSet &  ps,
        const uint        n )
{
    std::vector< TProcSet >  ps_split;

    ps.split( n, ps_split );

    return ps_split;
}

inline
TProcSet
join ( const TProcSet &  ps1,
       const TProcSet &  ps2 ) noexcept
{
    return ps1.join( ps2 );
}

//!
//! Compute intersection of processor sets
//!
inline
TProcSet
intersect ( const TProcSet &  ps1,
            const TProcSet &  ps2 ) noexcept
{
    if (( ps1.size() == 0 ) || ( ps2.size() == 0 ))
        return TProcSet();
    else
        return TProcSet( std::max( ps1.first(), ps2.first() ),
                         std::min( ps1.last(), ps2.last() ) );
}

}// namespace Hpro

#endif // __HPRO_TPROCSET_HH
