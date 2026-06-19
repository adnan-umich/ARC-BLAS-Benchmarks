#ifndef __HPRO_TINDEXSET_HH
#define __HPRO_TINDEXSET_HH
//
// Project     : HLIBpro
// File        : TIndexSet.hh
// Description : class for an indexset
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/types.hh"
#include "hpro/base/String.hh"

namespace Hpro
{

//
// foreach construct to go through indexsets
//
//#define foreach_is( i, set ) for ( idx_t i = (set).first(); i <= (set).last(); ++i )

/////////////////////////////////////////////////////////////////////////////
//!
//! \class TIndexSet
//! \brief Represents an indexset with contigously numbered indices, defined
//!        by the first and last index in the set.
//!
/////////////////////////////////////////////////////////////////////////////

class TIndexSet
{
public:
    //////////////////////////////////////////////////////////
    //
    // std::map comparator
    //

    struct map_cmp_t
    {
        bool  operator ()  ( const TIndexSet &  is1,
                             const TIndexSet &  is2 ) const noexcept
        {
            return is1.is_strictly_left_of( is2 );
        }
    };

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
    //! first index in the set
    idx_t  _first;

    //! last index in the set
    idx_t  _last;

public:
    ////////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct empty index set
    TIndexSet () noexcept
    {
        _first = 0;
        _last  = -1;
    }

    //! construct indexset if size \a n
    explicit
    TIndexSet ( const size_t  n ) noexcept
    {
        _first = 0;
        _last  = idx_t(n)-1;
    }

    //! construct indexset by first and last index
    TIndexSet ( const idx_t  afirst,
                const idx_t  alast ) noexcept
    {
        set_first_last( afirst, alast );
    }

    //! copy constructor
    TIndexSet ( const TIndexSet & is ) noexcept
    {
        *this = is;
    }
    
    ////////////////////////////////////////////////////////
    //
    // access internal data
    //

    //! return first index in set
    idx_t  first () const noexcept { return _first; }
    
    //! return last index in set
    idx_t  last  () const noexcept { return _last;  }

    //! return last index in set
    size_t size  () const noexcept
    {
        auto  n = _last - _first + 1;

        return n >= 0 ? n : 0;
    }
    
    //! set indexset by first and last index
    void set_first_last ( const idx_t  afirst,
                          const idx_t  alast ) noexcept
    {
        _first = afirst;
        _last  = alast; // max( afirst, alast );
    }

    //! set indexset by first and size
    void set_first_size ( const idx_t   afirstf,
                          const size_t  asize ) noexcept
    {
        _first = afirstf;
        _last  = idx_t(afirstf + asize) - 1;
    }

    //! return true if given index is part of indexset and false otherwise
    bool is_in ( const idx_t  idx ) const noexcept
    {
        return (( idx >= _first ) && ( idx <= _last ));
    }
    
    //! return true if given indexset is subset
    bool is_sub ( const TIndexSet & is ) const noexcept
    {
        return (( is.first() >= first() ) && ( is.last() <= last() ));
    }
    bool is_subset ( const TIndexSet & is ) const noexcept
    {
        return (( is.first() >= first() ) && ( is.last() <= last() ));
    }

    //! return true if local indexset is subset of given indexset
    bool is_subset_of ( const TIndexSet & is ) const noexcept
    {
        return is.is_subset( *this );
    }
        
    //! return true if indexset is empty
    bool is_empty () const noexcept
    {
        return size() == 0;
    }
        
    ////////////////////////////////////////////////////////
    //
    // STL interface
    //

    iterator        begin  () const noexcept { return iterator( _first );  }
    iterator        end    () const noexcept { return iterator( _last+1 ); }
        
    const_iterator  cbegin () const noexcept { return iterator( _first );  }
    const_iterator  cend   () const noexcept { return iterator( _last+1 ); }
        
    ////////////////////////////////////////////////////////
    //
    // misc. methods
    //

    //! copy operator
    TIndexSet & operator = ( const TIndexSet & is ) noexcept
    {
        _first = is._first;
        _last  = is._last;

        return *this;
    }
    
    //! equality operator
    bool  operator == ( const TIndexSet & is ) const noexcept
    {
        return ((_first == is._first) && (_last == is._last));
    }
    
    //! inequality operator
    bool  operator != ( const TIndexSet & is ) const noexcept
    {
        return ((_first != is._first) || (_last != is._last));
    }
    
    //! this is strictly left of \a is iff ∀ i ∈ this, j ∈ is : i < j
    bool  is_strictly_left_of   ( const TIndexSet & is ) const noexcept
    { return _last <  is._first; }

    //! this is left or equal to \a is iff this ∖ is < is
    bool  is_left_or_equal_to   ( const TIndexSet & is ) const noexcept
    { return ( _first <= is._first ) && ( _last <= is._last ); }

    //! this is strictly right of \a is iff ∀ i ∈ this, j ∈ is : i > j
    bool  is_strictly_right_of  ( const TIndexSet & is ) const noexcept
    { return _first > is._last; }

    //! this is right or equal to \a is iff this ∖ is > is
    bool  is_right_or_equal_to  ( const TIndexSet & is ) const noexcept
    { return ( _first >= is._first ) && ( _last >= is._last ); }
    
    //! string output
    std::string to_string () const
    {
        if ( size() == 0 ) return "{}";
        if ( size() == 1 ) return Hpro::to_string( "{%d}", _first );
        else               return Hpro::to_string( "{%d:%d}", _first, _last );
    }

    //! return size in bytes used by this object
    size_t byte_size () const { return 2*sizeof(idx_t); }
};

//////////////////////////////////////////////////////////
//
// functional ctors
//

inline
TIndexSet
is ( const idx_t  first,
     const idx_t  last ) noexcept
{
    return TIndexSet( first, last );
}

//////////////////////////////////////////////////////////
//
// set operations
//

//!
//! Compute union of indexsets.
//! ATTENTION: a gap between is1 and is2 is part of the result!
//!
inline
TIndexSet
join ( const TIndexSet &  is1,
       const TIndexSet &  is2 ) noexcept
{
    return TIndexSet( std::min( is1.first(), is2.first() ), std::max( is1.last(), is2.last() ) );
}

//!
//! Compute intersection of indexsets.
//!
inline
TIndexSet
intersect ( const TIndexSet &  is1,
            const TIndexSet &  is2 ) noexcept
{
    return TIndexSet( std::max( is1.first(), is2.first() ), std::min( is1.last(), is2.last() ) );
}

//
// return true if is0 intersects with is1
//
inline
bool
is_intersecting ( const TIndexSet &  is0,
                  const TIndexSet &  is1 )
{
    return  ( intersect( is0, is1 ).size() > 0 );
}

//////////////////////////////////////////////////////////
//
// arithmetics for indexsets
//

//!
//! add offset to indexset
//!
inline
TIndexSet
operator + ( const TIndexSet &  is,
             const idx_t        ofs ) noexcept
{
    return TIndexSet( is.first() + ofs, is.last() + ofs );
}

//!
//! subtract offset from indexset
//!
inline
TIndexSet
operator - ( const TIndexSet &  is,
             const idx_t        ofs ) noexcept
{
    return TIndexSet( is.first() - ofs, is.last() - ofs );
}

//////////////////////////////////////////////////////////
//
// misc.
//

//
// stream output
//
inline std::ostream &
operator << ( std::ostream &     os,
              const TIndexSet &  is )
{
    return os << is.to_string();
}

/////////////////////////////////////////////////////////////////////////////
//!
//! \class TBlockIndexSet
//! \brief Represents a product of two indexsets.
//!
/////////////////////////////////////////////////////////////////////////////

class TBlockIndexSet
{
private:
    //! row indexset
    TIndexSet  _row_is;

    //! column indexset
    TIndexSet  _col_is;

public:
    ////////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct empty index set
    TBlockIndexSet () noexcept
    {
    }

    //! construct indexset defined by given sets
    TBlockIndexSet ( const TIndexSet &  arow_is,
                     const TIndexSet &  acol_is ) noexcept
            : _row_is( arow_is )
            , _col_is( acol_is )
    {
    }

    //! copy constructor
    TBlockIndexSet ( const TBlockIndexSet & is ) noexcept
    {
        *this = is;
    }
    
    ////////////////////////////////////////////////////////
    //
    // indexset management
    //

    //! return row indexset
    const TIndexSet &  row_is () const noexcept { return _row_is; }
    
    //! return column indexset
    const TIndexSet &  col_is () const noexcept { return _col_is; }

    //! return true if given index (i,j) is part of indexset and false otherwise
    bool is_in ( const idx_t  row_idx,
                 const idx_t  col_idx ) const noexcept
    {
        return _row_is.is_in( row_idx ) && _col_is.is_in( col_idx );
    }

    //! return true if given indexset is subset
    bool is_sub ( const TBlockIndexSet & is ) const noexcept
    {
        return _row_is.is_sub( is.row_is() ) && _col_is.is_sub( is.col_is() );
    }
    bool is_subset ( const TBlockIndexSet & is ) const noexcept
    {
        return _row_is.is_sub( is.row_is() ) && _col_is.is_sub( is.col_is() );
    }
    
    //! return true if local indexset is subset of given indexset
    bool is_subset_of ( const TBlockIndexSet & is ) const noexcept
    {
        return is.is_subset( *this );
    }
    
    //! return true if indexset is empty
    bool is_empty () const noexcept
    {
        return ( _row_is.is_empty() || _col_is.is_empty() );
    }
        
    ////////////////////////////////////////////////////////
    //
    // misc. methods
    //

    //! copy operator
    TBlockIndexSet & operator = ( const TBlockIndexSet & is ) noexcept
    {
        _row_is = is._row_is;
        _col_is = is._col_is;

        return *this;
    }
    
    //! equality operator
    bool  operator == ( const TBlockIndexSet & is ) const noexcept
    {
        return ((_row_is == is._row_is) && (_col_is == is._col_is));
    }
    
    //! inequality operator
    bool  operator != ( const TBlockIndexSet & is ) const noexcept
    {
        return ! ( *this == is );
    }
    
    //! string output
    std::string  to_string () const
    {
        return _row_is.to_string() + "×" + _col_is.to_string();
    }

    //! return size in bytes used by this object
    size_t  byte_size () const noexcept { return _row_is.byte_size() + _col_is.byte_size(); }
};

//
// functional ctor
//

inline
TBlockIndexSet
bis ( const TIndexSet &  rowis,
      const TIndexSet &  colis ) noexcept
{
    return TBlockIndexSet( rowis, colis );
}

//
// stream output
//
inline std::ostream &
operator << ( std::ostream &          os,
              const TBlockIndexSet &  is )
{
    return os << is.to_string();
}

//
// return transposed block index set
//
inline
TBlockIndexSet
transpose ( const TBlockIndexSet &  bs ) noexcept
{
    return TBlockIndexSet( bs.col_is(), bs.row_is() );
}

//
// return true if is0 intersects with is1
//
inline
bool
is_intersecting ( const TBlockIndexSet &  is0,
                  const TBlockIndexSet &  is1 )
{
    return  ( is_intersecting( is0.row_is(), is1.row_is() ) &&
              is_intersecting( is0.col_is(), is1.col_is() ) );
}

}// namespace Hpro

#endif  // __HPRO_TINDEXSET_HH
