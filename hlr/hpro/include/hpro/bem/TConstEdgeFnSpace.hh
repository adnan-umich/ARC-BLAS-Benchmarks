#ifndef __HPRO_TCONSTEDGEFNSPACE_HH
#define __HPRO_TCONSTEDGEFNSPACE_HH
//
// Project     : HLIBpro
// File        : TConstEdgeFnSpace.hh
// Description : function space for constant normal linear tangential (CN/LT) edge elements
// Author      : Ronald Kriemann, Jonas Ballani
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/bem/TFnSpace.hh"

namespace Hpro
{
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// function space for constant normal linear tangential (CN/LT) 
// edge elements
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

// local class type
DECLARE_TYPE( TConstEdgeFnSpace );

class TConstEdgeFnSpace : public TFnSpace< T2Point >
{
public:
    //
    // value type of basis function
    //
    using  value_t = T2Point;
    
    // type for storing edges
    struct edge_t
    {
        idx_t vtx[2];
    };
    
protected:
    
    //! list of all edges
    std::vector< edge_t >    _edge_list;
    
public:
    ////////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TConstEdgeFnSpace ( const TGrid * agrid );

    virtual ~TConstEdgeFnSpace ();
    
    // give access to internal data
    edge_t  edge( const idx_t  i ) const { return _edge_list[i]; } 
    
    //
    // basis evaluation
    //

    //! return number of basis functions in unit triangle
    size_t  n_unit_bases    () const
    {
        return 3;
    }
    
    //! return local index in triangle \a tri for basis function \a phi
    idx_t  triangle_index  ( const idx_t  phi,
                             const idx_t  tri ) const
    {
        if      ( _tri_idx[ _tri_idx_ptr[ tri ] ]     == phi ) return 0;
        else if ( _tri_idx[ _tri_idx_ptr[ tri ] + 1 ] == phi ) return 1;
        else if ( _tri_idx[ _tri_idx_ptr[ tri ] + 2 ] == phi ) return 2;

        HERROR( ERR_CONSISTENCY, "(TConstEdgeFnSpace) triangle_index", "index not in triangle" );
    }
    
    //! return local index in triangle \a tri for basis function \a phi
    idx_t  triangle_index  ( const idx_t                phi,
                             const TGrid::triangle_t &  tri ) const
    {
        const idx_t tri_idx0 = _supp_list[ _supp_list_ptr[phi] ];
        const idx_t tri_idx1 = _supp_list[ _supp_list_ptr[phi] + 1 ];
    
        const TGrid::triangle_t  tri0 = _grid->triangle( tri_idx0 );
    
        idx_t vtx[3], vtx0[3];
        
        for ( uint j = 0; j < 3; j++ )
        {
            vtx[j]  = tri.vtx[j];
            vtx0[j] = tri0.vtx[j];
        }// for    
        
        uint  ncommon = 0;

        for ( uint j = 0; j < 3; j++ )
        for ( uint l = 0; l < 3; l++ )
        {
            if ( vtx[j] == vtx0[l] )
            {
                std::swap( vtx[ncommon], vtx[j] );
                std::swap( vtx0[ncommon], vtx0[l] );
                ncommon++;
                break;
            }// if
        }// for
        
        if ( ncommon == 3 ) return triangle_index( phi, tri_idx0 );
        else                return triangle_index( phi, tri_idx1 );
    }
    
   
    //! evaluate basis function with triangle local index \a i at
    //! unit coordinate ( \a s, \a t )
    value_t  eval_basis_unit  ( const idx_t    i,
                                const double   s,
                                const double   t ) const
    {
        switch ( i )
        {
            case 0 : return value_t( s, t );
            case 1 : return value_t( s - 1.0, t );
            case 2 : return value_t( s, t - 1.0 );
        }// switch

        return value_t(0,0);
    }
              
    // evaluate basis function i on triangle defined by vtxidxs 
    // for all quadrature points given on the reference triangle,
    // note that the value still has to be scaled by +/-0.5 * |e|/|tri|
    // - return false if triangle not in support
    bool eval_basis ( const idx_t                     i,
                      const TGrid::triangle_t &       vtxidxs,
                      const std::vector< T2Point > &  ref_points,
                      std::vector< T3Point > &        values ) const;
                    
    // evaluate basis function i on triangle defined by vtxidxs 
    // for all quadrature points which have been already transformed
    // to the given triangle,
    // note that the value still has to be scaled by +/-0.5 * |e|/|tri|
    // - return false if triangle not in support
    bool eval_basis ( const idx_t                     i,
                      const TGrid::triangle_t &       tri,
                      const std::vector< T3Point > &  quad_points,
                      std::vector< T3Point > &        values ) const;
                     
    // scaling factor for basis function i on triangle tri = +/-0.5 * |e|/|tri|
    double  scaling_factor_basis ( const idx_t  i,
                                   const idx_t  tri ) const;
                                
    // transform quadrature points on reference triangle to general triangle
    // defined by its vertices
    void  ref_points2triangle ( const TGrid::triangle_t &       tri,
                                const std::vector< T2Point > &  ref_points, 
                                std::vector< T3Point > &        quad_points ) const;
       
    //
    // RTTI
    //

    HPRO_RTTI_DERIVED( TConstEdgeFnSpace, TFnSpace );

protected:
    
    // type for storing a vertex belonging to an edge
    using  vtx_edge_t      = struct { idx_t  vtx; idx_t  edge; };

    // type for storing a list of all edges emanating from a vertex
    using  vtx_edge_list_t = std::list< vtx_edge_t >;

    ////////////////////////////////////////////////////////
    //
    // private functions
    //

    // construct function space by building indices
    // and their support
    void construct ();
};

}// namespace

#endif // _TCONSTEDGEFNSPACE_HH
