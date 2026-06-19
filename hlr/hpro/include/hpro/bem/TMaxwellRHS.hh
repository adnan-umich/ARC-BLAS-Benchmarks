#ifndef __HPRO_TMAXWELLRHS_HH
#define __HPRO_TMAXWELLRHS_HH
//
// Project     : HLIBpro
// File        : TMaxwellRHS.hh
// Description : classes for building RHS-vectors in Maxwell applications
// Author      : Ronald Kriemann, Jonas Ballani
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/bem/TBEMRHS.hh"

namespace Hpro
{
//
// value type of BEM function with traits
//
struct complex3_t
{
    std::complex< double >  val[3];
};

template <>
struct real_type< complex3_t >
{
    using  type_t = double;
};

template <>
struct is_complex_type< complex3_t >
{
    static const bool value = true;
};

////////////////////////////////////////////////////////
//
// building RHS for Maxwell EFIE using quadrature
//

template < typename T_fnspace, typename T_rhsfn >
class TMaxwellEFIERHS : public TBEMRHS< T_fnspace, T_rhsfn, std::complex< double > >
{
public:
    //! template argument as internal type
    using  value_t   = std::complex< double >;
    using  rhsfn_t   = T_rhsfn;
    using  fnspace_t = T_fnspace;

private:
    // quadrature order
    const uint  _quad_order;
    
public:
    ///////////////////////////////////////////////////////////////
    //
    // constructur and destructor
    //

    TMaxwellEFIERHS ( const uint  aquad_order )
            : _quad_order( aquad_order )
    {}

    virtual ~TMaxwellEFIERHS () {}

    ///////////////////////////////////////////////////////////////
    //!
    //! build RHS-vector for \a fnspace and RHS
    //! function \a rhs; by \a perm one can set a permutation for
    //! to final vector, e.g., from internal to external numbering
    //!
    virtual
    std::unique_ptr< TVector< value_t > >
    build ( const fnspace_t *     fnspace,
            const rhsfn_t *       rhs,
            const TPermutation *  perm = NULL ) const
    {
        HASSERT( fnspace != nullptr, ERR_ARG, "(TMaxwellEFIERHS) build", "function space is nullptr" );
        HASSERT( rhs     != nullptr, ERR_ARG, "(TMaxwellEFIERHS) build", "RHS is nullptr" );

        using  real_t = real_type_t< value_t >;

        //
        // build vector
        //

        std::vector< T3Point >  t( 3 );      // holds triangle coordinates
        std::vector< T3Point >  qpts;        // quadrature points for each triangle
        std::vector< double >   qwghts;      // quadrature weights for each triangle
        TTriGaussQuad           quad;        // Gauss quadrature rules for triangles
        std::vector< T3Point >  phi_values;  // basis function values
        TGrid::triangle_t       vid;         // vertex indices of triangle
        auto                    v = std::make_unique< TScalarVector< value_t > >( fnspace->n_indices(), 0 );

        for ( idx_t  idx = 0; idx < idx_t(fnspace->n_indices()); idx++ )
        {
            const idx_t  pidx = ( perm != nullptr ? perm->permute( idx ) : idx );

            //
            // integrate rhs function over support of index
            //

            for ( auto  tri : fnspace->support( idx ) )
            {
                const auto  J = real_t( fnspace->grid()->tri_size( tri ) );

                vid = fnspace->grid()->triangle( tri );
            
                for ( uint vtx = 0; vtx < 3; vtx++ )
                    t[vtx] = fnspace->grid()->vertex( vid.vtx[vtx] );

                const auto  n = fnspace->grid()->tri_normal( tri );

                // get quadrature points
                quad.build( _quad_order, t, qpts, qwghts );

                const auto  npts = qpts.size();

                if ( phi_values.size() < npts )
                    phi_values.resize( npts );

                // evaluate basis function
                fnspace->eval_basis( idx, vid, qpts, phi_values );

                const auto  scale_phi = fnspace->scaling_factor_basis( idx, tri );

                // integrate
                value_t  f = value_t(0);

                for ( size_t  j = 0; j < npts; j++ )
                {
                    auto  fn_val  = rhs->eval( qpts[j], n );
                    auto  dot_val = ( fn_val.val[0] * real_t( phi_values[j][0] ) +
                                      fn_val.val[1] * real_t( phi_values[j][1] ) +
                                      fn_val.val[2] * real_t( phi_values[j][2] ) );
                                     
                    f += real_t(qwghts[j]) * dot_val;
                }// for
            
                v->add_entry( pidx, J * f * scale_phi );
            }// for
        }// for

        return std::unique_ptr< TVector< value_t > >( v.release() );
    }


    DISABLE_COPY_OP( TMaxwellEFIERHS );
};

////////////////////////////////////////////////////////
//
// building RHS for Maxwell MFIE using quadrature
//

template < typename T_fnspace,
           typename T_rhsfn >
class TMaxwellMFIERHS : public TBEMRHS< T_fnspace, T_rhsfn, std::complex< double > >
{
public:
    //! template argument as internal type
    using  value_t   = std::complex< double >;
    using  rhsfn_t   = T_rhsfn;
    using  fnspace_t = T_fnspace;
          
private:
    // quadrature order
    const uint  _quad_order;
    
public:
    ///////////////////////////////////////////////////////////////
    //
    // constructur and destructor
    //

    TMaxwellMFIERHS ( const uint  aquad_order )
            : _quad_order( aquad_order )
    {}

    virtual ~TMaxwellMFIERHS () {}

    ///////////////////////////////////////////////////////////////
    //!
    //! build RHS-vector for \a fnspace and RHS
    //! function \a rhs; by \a perm one can set a permutation for
    //! to final vector, e.g., from internal to external numbering
    //!
    virtual
    std::unique_ptr< TVector< value_t > >
    build ( const fnspace_t *     fnspace,
            const rhsfn_t *       rhs,
            const TPermutation *  perm = nullptr ) const
    {
        HASSERT( fnspace != nullptr, ERR_ARG, "(TMaxwellMFIERHS) build", "function space is nullptr" );
        HASSERT( rhs     != nullptr, ERR_ARG, "(TMaxwellMFIERHS) build", "RHS is nullptr" );

        using  real_t = real_type_t< value_t >;

        //
        // build vector
        //

        std::vector< T3Point >  t( 3 );      // holds triangle coordinates
        std::vector< T3Point >  qpts;        // quadrature points for each triangle
        std::vector< double >   qwghts;      // quadrature weights for each triangle
        TTriGaussQuad           quad;        // Gauss quadrature rules for triangles
        std::vector< T3Point >  phi_values;  // basis function values
        TGrid::triangle_t       vid;         // vertex indices of triangle
        auto                    v = std::make_unique< TScalarVector< value_t > >( fnspace->n_indices(), 0 );

        for ( idx_t  idx = 0; idx < idx_t(fnspace->n_indices()); idx++ )
        {
            const idx_t  pidx = ( perm != nullptr ? perm->permute( idx ) : idx );

            //
            // integrate rhs function over support of index
            //

            for ( auto  tri : fnspace->support( idx ) )
            {
                const auto  J = real_t( fnspace->grid()->tri_size( tri ) );

                vid = fnspace->grid()->triangle( tri );
            
                for ( uint vtx = 0; vtx < 3; vtx++ )
                    t[vtx] = fnspace->grid()->vertex( vid.vtx[vtx] );

                const auto  n = fnspace->grid()->tri_normal( tri );

                // get quadrature points
                quad.build( _quad_order, t, qpts, qwghts );

                const auto  npts = qpts.size();

                if ( phi_values.size() < npts )
                    phi_values.resize( npts );

                // evaluate basis function
                fnspace->eval_basis( idx, vid, qpts, phi_values );

                const auto  scale_phi = fnspace->scaling_factor_basis( idx, tri );

                // integrate
                auto  f = value_t(0);

                for ( size_t  j = 0; j < npts; j++ )
                {
                    auto        fn_val  = rhs->eval( qpts[j], n );
                    const auto  n_x_phi = cross( n, phi_values[j] );
                    auto        dot_val = ( fn_val.val[0] * real_t( n_x_phi[0] ) +
                                               fn_val.val[1] * real_t( n_x_phi[1] ) +
                                               fn_val.val[2] * real_t( n_x_phi[2] ) );

                    f += real_t(qwghts[j]) * dot_val;
                }// for

                v->add_entry( pidx, J * f * (-scale_phi) );
            }// for
        }// for

        return std::unique_ptr< TVector< value_t > >( v.release() );
    }


    DISABLE_COPY_OP( TMaxwellMFIERHS );
};

}// namespace Hpro

#endif  // __HPRO_TMAXWELLRHS_HH
