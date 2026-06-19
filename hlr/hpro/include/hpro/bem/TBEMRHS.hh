#ifndef __HPRO_TBEMRHS_HH
#define __HPRO_TBEMRHS_HH
//
// Project     : HLIBpro
// File        : TBEMRHS.hh
// Description : classes for building RHS-vectors in BEM-applications
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/cluster/TIndexSet.hh"
#include "hpro/cluster/TPermutation.hh"

#include "hpro/bem/TFnSpace.hh"
#include "hpro/bem/TGaussQuad.hh"

#include "hpro/vector/TVector.hh"

namespace Hpro
{

////////////////////////////////////////////////////////
//
// baseclass defining interface
//
template < typename  T_fnspace,
           typename  T_rhsfn,
           typename  T_val >
class TBEMRHS
{
public:
    //! template arguments as internal types
    using  fnspace_t = T_fnspace;
    using  rhsfn_t   = T_rhsfn;
    using  value_t   = T_val;
    
public:
    ///////////////////////////////////////////////////////////////
    //
    // constructur and destructor
    //

    TBEMRHS () {}

    virtual ~TBEMRHS () {}

    ///////////////////////////////////////////////////////////////
    //!
    //! build RHS-vector for \a fnspace and RHS
    //! function \a rhs; by \a perm one can set a permutation for
    //! to final vector, e.g., from internal to external numbering
    //!
    virtual std::unique_ptr< TVector< value_t > >  build ( const fnspace_t *     fnspace,
                                                           const rhsfn_t *       rhs,
                                                           const TPermutation *  perm = NULL ) const = 0;
};

////////////////////////////////////////////////////////
//
// building RHS using quadrature
//
template < typename  T_fnspace,
           typename  T_rhsfn,
           typename  T_val >
class TQuadBEMRHS : public TBEMRHS< T_fnspace, T_rhsfn, T_val >
{
public:
    //! template arguments as internal types
    using  fnspace_t = T_fnspace;
    using  rhsfn_t   = T_rhsfn;
    using  value_t   = T_val;
    
private:
    // quadrature order
    const uint  _quad_order;
    
public:
    ///////////////////////////////////////////////////////////////
    //
    // constructur and destructor
    //

    TQuadBEMRHS ( const uint  aquad_order )
            : _quad_order( std::max( 1u, aquad_order ) )
    {}

    virtual ~TQuadBEMRHS () {}

    ///////////////////////////////////////////////////////////////
    //!
    //! build RHS-vector for \a fnspace and RHS
    //! function \a rhs; by \a perm one can set a permutation for
    //! to final vector, e.g., from internal to external numbering
    //!
    virtual std::unique_ptr< TVector< value_t > >  build ( const fnspace_t *     fnspace,
                                                           const rhsfn_t *       rhs,
                                                           const TPermutation *  perm = NULL ) const
    {
        HASSERT( fnspace != nullptr, ERR_ARG, "(TQuadBEMRHS) build", "function space is nullptr" );
        HASSERT( rhs     != nullptr, ERR_ARG, "(TQuadBEMRHS) build", "RHS is nullptr" );

        using  real_t = real_type_t< value_t >;
    
        //
        // build vector
        //

        std::vector< T2Point >  qpts;   // quadrature points for each triangle
        std::vector< double >   qwghts; // quadrature weights for each triangle
        TTriGaussQuad           quad;   // Gauss quadrature rules for triangles
        const TGrid *           grid = fnspace->grid();
        auto                    v    = std::make_unique< TScalarVector< value_t > >( fnspace->n_indices(), 0 );
    
        // get quadrature points
        quad.build( _quad_order, qpts, qwghts );

        const size_t  npts = qpts.size();
    
        for ( idx_t  idx = 0; idx < idx_t(fnspace->n_indices()); idx++ )
        {
            const idx_t  pidx = ( perm != nullptr ? perm->permute( idx ) : idx );
            value_t      val  = value_t(0);

            //
            // integrate rhs function over support of index
            //

            for ( auto  tri_idx : fnspace->support( idx ) )
            {
                const auto  tri = grid->triangle( tri_idx );
                const auto  v0  = grid->vertex( tri.vtx[0] );
                const auto  v1  = grid->vertex( tri.vtx[1] );
                const auto  v2  = grid->vertex( tri.vtx[2] );
                const auto  J   = real_t( grid->tri_size( tri_idx ) );
                const auto  n   = grid->tri_normal( tri_idx );

                //
                // evaluate at quadrature points
                //
            
                value_t  f = value_t(0);
            
                for ( size_t  j = 0; j < npts; j++ )
                {
                    // transform to triangle coordinates
                    const auto  a1 = qpts[j][0];
                    const auto  a2 = qpts[j][1];
                    const auto  a0 = 1.0 - a1 - a2;
                    const auto  x  = a0 * v0 + a1 * v1 + a2 * v2;
                
                    f += real_t( qwghts[j] ) * rhs->eval( x, n ) * value_t( fnspace->eval_basis_unit( idx, a1, a2, tri.vtx ) );
                }// for

                val += J * f;
            }// for
            
            v->set_entry( pidx, val );
        }// for

        return std::unique_ptr< TVector< value_t > >( v.release() );
    }

    DISABLE_COPY_OP( TQuadBEMRHS );
};

}// namespace

#endif  // __HPRO_TBEMRHS_HH
