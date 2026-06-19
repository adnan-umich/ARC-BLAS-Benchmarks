#ifndef __HPRO_TFNSPACE_HH
#define __HPRO_TFNSPACE_HH
//
// Project     : HLIBpro
// File        : TFnSpace.hh
// Description : implements function spaces over grids
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/base/TTypeInfo.hh"

#include "hpro/cluster/TCoordinate.hh"

#include "hpro/bem/TGrid.hh"

namespace Hpro
{

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// base class for a function space over grids
// - implements grid and index management
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

// local class type
DECLARE_TYPE( TFnSpace );

template < typename T_value = double >
class TFnSpace : public TTypeInfo
{
public:
    //
    // value type of basis function
    //

    using  value_t = T_value;
    
protected:
    // grid the function space lives on
    const TGrid *           _grid;

    // array holding positions of indices
    std::vector< T3Point >  _indices;

    // mapping of indices to support
    // - entries from _supp_list_ptr[i] to _supp_list_ptr[i+1]-1 in the array
    //   _supp_list defines the triangles in the support of index i
    std::vector< idx_t >    _supp_list_ptr;  
    std::vector< idx_t >    _supp_list;

    // mapping of triangles to indices
    // - entries from _tri_idx_ptr[i] to _tri_idx_ptr[i+1]-1 in the array
    //   _tri_idx defines the indices local to triangle i
    std::vector< idx_t >    _tri_idx_ptr;  
    std::vector< idx_t >    _tri_idx;

public:
    
    //
    // list of triangles in support
    //
    class TSupportSet
    {
    public:
        //
        // STL iterator
        //
        class iterator
        {
        private:
            // local function space
            const TFnSpace *  _fn_space;

            // current position in list
            idx_t             _pos;

        public:
            // ctor for specific triangle and function space
            iterator ( const TFnSpace *  fn_space,
                       const idx_t       pos )
                    : _fn_space( fn_space )
                    , _pos( pos )
          
            {
                if ( fn_space == nullptr )
                    HERROR( ERR_ARG, "(TSupportSet::iterator) ctor", "function space is NULL" );
            }

            iterator  operator ++ ()    { return iterator( _fn_space, ++_pos ); }
            iterator  operator ++ (int) { return iterator( _fn_space, _pos++ ); }

            iterator  operator -- ()    { return iterator( _fn_space, --_pos ); }
            iterator  operator -- (int) { return iterator( _fn_space, _pos-- ); }

            idx_t  operator * () const { return _fn_space->supp_list( _pos ); }

            bool  operator == ( const iterator &  it ) const
            {
                return (_pos == it._pos) && (_fn_space == it._fn_space);
            }

            bool  operator != ( const iterator &  it ) const
            {
                return ! operator == ( it );
            }
        };

        using  const_iterator = iterator;
        
    private:
        // local function space
        const TFnSpace *  _fn_space;

        // lower and upper bound for support list
        const idx_t       _lb, _ub;

    public:
        // ctor for specific triangle and function space
        TSupportSet ( const TFnSpace *  fn_space,
                      const idx_t       lb,
                      const idx_t       ub )
                : _fn_space( fn_space )
                , _lb( lb )
                , _ub( ub )
          
        {
            if ( fn_space == nullptr )
                HERROR( ERR_ARG, "(TSupportSet) ctor", "function space is NULL" );
        }

        iterator        begin  () const { return iterator( _fn_space, _lb ); }
        iterator        end    () const { return iterator( _fn_space, _ub ); }

        const_iterator  cend   () const { return const_iterator( _fn_space, _ub ); }
        const_iterator  cbegin () const { return const_iterator( _fn_space, _lb ); }

        //! number of triangles
        size_t          size   () const { return _ub - _lb; }
    };

    //
    // represents set of indices per triangle
    //
    class TTriangleIndexSet
    {
    public:
        //
        // iterator for TTriangleIndexSet
        //
        class iterator
        {
        private:
            // local function space
            const TFnSpace *  _fn_space;

            // current position in list
            idx_t             _pos;

        public:
            // ctor for specific triangle and function space
            iterator ( const TFnSpace *  fn_space,
                       const idx_t       pos )
                    : _fn_space( fn_space )
                    , _pos( pos )
          
            {
                if ( fn_space == NULL )
                    HERROR( ERR_ARG, "(TTriangleIndexIterator) ctor", "function space is NULL" );
            }

            iterator  operator ++ ()    { return iterator( _fn_space, ++_pos ); }
            iterator  operator ++ (int) { return iterator( _fn_space, _pos++ ); }

            iterator  operator -- ()    { return iterator( _fn_space, --_pos ); }
            iterator  operator -- (int) { return iterator( _fn_space, _pos-- ); }

            idx_t  operator * () const { return _fn_space->tri_idx( _pos ); }

            bool  operator == ( const iterator &  it ) const
            {
                return (_pos == it._pos) && (_fn_space == it._fn_space);
            }

            bool  operator != ( const iterator &  it ) const
            {
                return ! operator == ( it );
            }
        };
    
        using  const_iterator = iterator;
        
    private:
        // local function space
        const TFnSpace *  _fn_space;

        // lower and upper bound for index list
        const idx_t       _lb, _ub;

    public:
        // ctor for specific triangle and function space
        TTriangleIndexSet ( const TFnSpace *  fn_space,
                            const idx_t       lb,
                            const idx_t       ub )
                : _fn_space( fn_space )
                , _lb( lb )
                , _ub( ub )
        {
            if ( fn_space == NULL )
                HERROR( ERR_ARG, "(TTriangleIndexSet) ctor", "function space is NULL" );
        }

        iterator        begin  () const { return iterator( _fn_space, _lb ); }
        iterator        end    () const { return iterator( _fn_space, _ub ); }
        
        const_iterator  cend   () const { return const_iterator( _fn_space, _ub ); }
        const_iterator  cbegin () const { return const_iterator( _fn_space, _lb ); }

        //! number of indices
        size_t          size   () const { return _ub - _lb; }
    };
    
public:
    ////////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TFnSpace ( const TGrid * agrid );

    virtual ~TFnSpace ();
    
    ////////////////////////////////////////////////////////
    //
    // give access to local data
    //

    //! return underlying grid
    const TGrid *    grid           ()                 const { return _grid; }

    //! return number of indices, e.g. basis functions
    size_t           n_indices      ()                 const { return _indices.size(); }

    //! return coordinate vector of \a i'th index
    const T3Point &  index_coord    ( const idx_t  i ) const { return _indices[i]; }

    //! access index support information
    idx_t            supp_list      ( const idx_t  i ) const { return _supp_list[i]; }

    //! return set of triangles forming support for basis function \a i
    TSupportSet      support        ( const idx_t  i ) const
    {
        return TSupportSet( this, _supp_list_ptr[i], _supp_list_ptr[i+1] );
    }
        
    //! access triangle index information
    idx_t            tri_idx        ( const idx_t  i ) const { return _tri_idx[i]; }

    //! return iterator for local indices of triangle \a i
    TTriangleIndexSet  triangle_indices ( const idx_t  i ) const
    {
        return TTriangleIndexSet( this, _tri_idx_ptr[i], _tri_idx_ptr[i+1] );
    }
        
    ////////////////////////////////////////////////////////
    //
    // function evaluation
    //

    //! evaluate function \a fn at all index positions on grid
    //! and build corresponding vector
    template < typename bemfn_t >
    std::unique_ptr< TScalarVector< value_type_t< bemfn_t > > >
    eval ( const bemfn_t *  fn ) const
    {
        using fn_value_t = value_type_t< bemfn_t >;
        
        if ( fn == nullptr )
            HERROR( ERR_ARG, "(TFnSpace) eval", "given function is nullptr" );

        //
        // build vector
        //

        auto  v = std::make_unique< TScalarVector< fn_value_t > >( this->n_indices(), 0 );
    
        for ( uint i = 0; i < this->n_indices(); i++ )
        {
            // get index position
            const auto  t = index_coord( i );

            // get normal at index (only valid for constant ansatz)
            const auto  n = _grid->tri_normal( _supp_list_ptr[i] );

            // evaluate
            v->set_entry( i, fn->eval( t, n ) );
        }// for

        return v;
    }
        
    
    ////////////////////////////////////////////////////////
    //
    // misc.
    //

    //! return coordinate info for degrees of freedom (with or without bbox info)
    std::unique_ptr< TCoordinate >  build_coord ( const bool  with_bbox = true ) const;

    //! return size of function space in bytes
    virtual size_t byte_size () const;

    //
    // RTTI
    //

    HPRO_RTTI_BASE( TFnSpace );
};

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// function space for piecewise constant grid functions 
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

// local class type
DECLARE_TYPE( TConstFnSpace );

template < typename T_value = double >
class TConstFnSpace : public TFnSpace< T_value >
{
public:
    //
    // value type of basis function
    //

    using  value_t = T_value;
    
    static_assert( std::is_floating_point< value_t >::value, "only supporting floating point type" );

public:
    ////////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TConstFnSpace ( const TGrid * agrid );

    virtual ~TConstFnSpace ();

    //
    // basis evaluation
    //

    //! return number of basis functions in unit triangle
    size_t   n_unit_bases      () const
    {
        return 1;
    }
    
    //! return local index in triangle \a tri for basis function \a phi
    idx_t    triangle_index    ( const idx_t,
                                 const TGrid::triangle_t & ) const
    {
        return 0;
    }
    
    //! evaluate basis function \a i at unit coordinate ( \a x, \a y )
    value_t  eval_basis_unit   ( const idx_t,
                                 const double,
                                 const double,
                                 const idx_t * ) const
    {
        return value_t(1);
    }
    
    value_t  eval_basis_unit   ( const idx_t,
                                 const double,
                                 const double ) const
    {
        return value_t(1);
    }
    
    //
    // RTTI
    //

    HPRO_RTTI_DERIVED( TConstFnSpace, TFnSpace< value_t > );

protected:
    ////////////////////////////////////////////////////////
    //
    // private functions
    //

    // construct function space by building indices
    // and their support
    void construct ();
};

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// function space for piecewise linear grid functions 
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

DECLARE_TYPE( TLinearFnSpace );

template < typename T_value = double >
class TLinearFnSpace : public TFnSpace< T_value >
{
public:
    //
    // value type of basis function
    //

    using  value_t = T_value;

    static_assert( std::is_floating_point< value_t >::value, "only supporting floating point type" );
        
public:
    ////////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TLinearFnSpace ( const TGrid * agrid );

    virtual ~TLinearFnSpace ();

    //
    // basis evaluation
    //

    //! return number of basis functions in unit triangle
    size_t   n_unit_bases      () const
    {
        return 3;
    }
    
    //! return local index in triangle \a tri for basis function \a phi
    idx_t    triangle_index    ( const idx_t                phi,
                                 const TGrid::triangle_t &  tri ) const
    {
        if      ( tri.vtx[0] == phi ) return 0;
        else if ( tri.vtx[1] == phi ) return 1;
        else if ( tri.vtx[2] == phi ) return 2;

        HERROR( ERR_CONSISTENCY, "(TLinearFnSpace) triangle_index", "index not in triangle" );
    }
    
    //! evaluate basis function \a i at unit coordinate ( \a s, \a t )
    //! on triangle \a triidx; \a trivtx contain the triangle coordinates
    //! in changed order such that (\a s,\a t) are defined w.r.t. \a trivtx[0]
    value_t  eval_basis_unit   ( const idx_t    i,
                                 const double   s,
                                 const double   t,
                                 const idx_t *  vtxidxs ) const
    {
        if      ( vtxidxs[0] == i ) return value_t(1-s-t);
        else if ( vtxidxs[1] == i ) return value_t(s);
        else if ( vtxidxs[2] == i ) return value_t(t);

        return value_t(0);
    }
    
    //! evaluate basis function with triangle local index \a i at
    //! unit coordinate ( \a s, \a t )
    value_t  eval_basis_unit   ( const idx_t    i,
                                 const double   s,
                                 const double   t ) const
    {
        switch ( i )
        {
            case 0 : return value_t(1-s-t);
            case 1 : return value_t(s);
            case 2 : return value_t(t);
        }// switch

        return value_t(0);
    }
    
    //
    // RTTI
    //
    
    HPRO_RTTI_DERIVED( TLinearFnSpace, TFnSpace< value_t > );

protected:
    ////////////////////////////////////////////////////////
    //
    // private functions
    //

    // construct function space by building indices
    // and their support
    void construct ();
};

}// namespace

#endif // _TFNSPACE_HH
