//
// Project     : HLIBpro
// File        : TMassBF.cc
// Description : bilinear form for mass matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/bem/TGaussQuad.hh"
#include "hpro/bem/TConstEdgeFnSpace.hh"

#include "hpro/bem/TMassBF.hh"

namespace Hpro
{

//
// ctor
//
template < typename T_ansatzsp,
           typename T_testsp,
           typename T_value >
TMassBF< T_ansatzsp, T_testsp, T_value >::TMassBF ( const T_ansatzsp *  aansatzsp,
                                                    const T_testsp *    atestsp,
                                                    const uint          aorder )
        : TBEMBF< T_ansatzsp, T_testsp, value_t >( aansatzsp, atestsp )
{
    //
    // build quad-points and weights
    //

    TTriGaussQuad  gaussquad;

    gaussquad.build( std::max( 1u, aorder ), _quad_pts, _quad_wghts );
}

//
// evaluate subblock defined by \a row_ind × \a col_ind; the indices
// in \a row_ind and \a col_ind can be arbitrary, e.g. must not be
// contiguous
//
template < typename T_ansatzsp,
           typename T_testsp,
           typename T_value >
void
TMassBF< T_ansatzsp, T_testsp, T_value >::eval  ( const std::vector< idx_t > &  row_ind,
                                                  const std::vector< idx_t > &  col_ind,
                                                  BLAS::Matrix< value_t > &     values ) const
{
    const size_t  nrows = row_ind.size();
    const size_t  ncols = col_ind.size();

    for ( size_t  col = 0; col < ncols; col++ )
    {
        const idx_t  j = col_ind[col];
        
        for ( size_t  row = 0; row < nrows; row++ )
        {
            const idx_t  i = row_ind[row];
        
            //
            // integrate over test space
            //

            auto        value     = value_t(0);
            const auto  ansatz_sp = this->ansatz_space();
            const auto  test_sp   = this->test_space();
            auto        support_i = ansatz_sp->support(i);
            auto        support_j = test_sp->support(j);

            if ( ansatz_sp->grid() != test_sp->grid() )
                HERROR( ERR_NOT_IMPL, "(TMassBF) eval", "unequal grids not supported" );
    
            for ( auto  tri1 : support_j )
            {
                const auto  J1 = value_t( test_sp->grid()->tri_size( tri1 ) );
        
                //
                // test if support of ansatz function is disjoint
                // and immediately go one if so
                //

                bool  intersects = false;
        
                for ( auto  tri0 : support_i )
                {
                    if ( tri0 == tri1 ) // only valid for equal grids!
                    {
                        intersects = true;
                        break;
                    }// if
                }// for

                if ( ! intersects )
                    continue;
        
                //
                // evaluate basis functions at quadrature points
                //

                idx_t         t1[3];
                const size_t  npts   = _quad_pts.size();
                auto          tvalue = value_t(0);
        
                for ( uint  vtx = 0; vtx < 3; vtx++ )
                    t1[vtx] = test_sp->grid()->triangle( tri1 ).vtx[vtx];
        
                for ( size_t  k = 0; k < npts; ++k )
                {
                    const double  a1 = _quad_pts[k][0];
                    const double  a2 = _quad_pts[k][1];
            
                    tvalue += ( value_t( _quad_wghts[k] ) *
                                ansatz_sp->eval_basis_unit( i, a1, a2, t1 ) * // TODO (see above)
                                test_sp->eval_basis_unit( j, a1, a2, t1 ) );
                }// for

                value += J1 * tvalue;
            }// for
    
            values( idx_t( row ), idx_t( col ) ) = value;
        }// for
    }// for
}

template <>
void
TMassBF< TConstEdgeFnSpace, TConstEdgeFnSpace, float >::eval  ( const std::vector< idx_t > &,
                                                                const std::vector< idx_t > &,
                                                                BLAS::Matrix< value_t > & ) const
{}

template <>
void
TMassBF< TConstEdgeFnSpace, TConstEdgeFnSpace, double >::eval  ( const std::vector< idx_t > &,
                                                                 const std::vector< idx_t > &,
                                                                 BLAS::Matrix< value_t > & ) const
{}

template <>
void
TMassBF< TConstEdgeFnSpace, TConstEdgeFnSpace, std::complex< float > >::eval  ( const std::vector< idx_t > &,
                                                                                const std::vector< idx_t > &,
                                                                                BLAS::Matrix< value_t > & ) const
{}

template <>
void
TMassBF< TConstEdgeFnSpace, TConstEdgeFnSpace, std::complex< double > >::eval  ( const std::vector< idx_t > &,
                                                                                 const std::vector< idx_t > &,
                                                                                 BLAS::Matrix< value_t > & ) const
{}

// return format of bilinear form, e.g. symmetric
template < typename T_ansatzsp,
           typename T_testsp,
           typename T_value >
matform_t
TMassBF< T_ansatzsp, T_testsp, T_value >::format () const
{
    if ( reinterpret_cast< const void * >( this->ansatz_space() ) == reinterpret_cast< const void * >( this->test_space() ) )
        return symmetric;
    else
        return unsymmetric;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

#define INST_ALL( type1, type2 )                                    \
    template class TMassBF< TConstFnSpace< type1 >,  TConstFnSpace< type1 >, type2 >; \
    template class TMassBF< TConstFnSpace< type1 >,  TLinearFnSpace< type1 >, type2 >; \
    template class TMassBF< TLinearFnSpace< type1 >, TConstFnSpace< type1 >, type2 >; \
    template class TMassBF< TLinearFnSpace< type1 >, TLinearFnSpace< type1 >, type2 >;

INST_ALL( float,  float )
INST_ALL( double, double )
INST_ALL( float,  std::complex< float > )
INST_ALL( double, std::complex< double > )

template class TMassBF< TConstEdgeFnSpace, TConstEdgeFnSpace, float >;
template class TMassBF< TConstEdgeFnSpace, TConstEdgeFnSpace, double >;
template class TMassBF< TConstEdgeFnSpace, TConstEdgeFnSpace, std::complex< float > >;
template class TMassBF< TConstEdgeFnSpace, TConstEdgeFnSpace, std::complex< double > >;

}// namespace Hpro
