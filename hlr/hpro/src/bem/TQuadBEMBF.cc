//
// Project     : HLIBpro
// File        : TBEMBF.cc
// Description : classes for bilinearforms in BEM-applications
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <set>

#include "unordered_map.hh"

#include "hpro/base/error.hh"
#include "hpro/base/config.hh"
#include "hpro/bem/TConstEdgeFnSpace.hh"
#include "hpro/bem/TSauterTriQuad.hh"

#include "hpro/bem/TQuadBEMBF.hh"

namespace Hpro
{

namespace B = BLAS;

//
// ctor
//
template < typename ansatzsp_t,
           typename testsp_t,
           typename value_t >
TQuadBEMBF< ansatzsp_t, testsp_t, value_t >::TQuadBEMBF ( const ansatzsp_t *  aansatzsp,
                                                          const testsp_t *    atestsp,
                                                          const uint          aorder,
                                                          const bool          dist_ada,
                                                          const real_t        ada_err )
        : TBEMBF< ansatzsp_t, testsp_t, value_t >( aansatzsp, atestsp )
        , _quad_order( aorder )
        , _quad_dist_adaptive( dist_ada )
        , _ada_error( ada_err )
{
    //
    // build quad-points and with up to max-order
    //

    TSauterTriQuad  quad;

    _quad_rules.resize( 4 );
    
    for ( uint  ncommon = 0; ncommon < 4; ++ncommon )
    {
        _quad_rules[ ncommon ].resize(  _quad_order + 1 );
        
        for ( uint  order = 1; order < _quad_order + 1; ++order )
        {
            auto  rule = std::vector< TSauterTriQuad::rule_t >();
            
            quad.build( ncommon, order, rule  );

            //
            // copy data into decoupled form
            //

            const size_t  npts        = rule.size();
            const size_t  padded_npts = CFG::Mach::simd_padded_size< real_t >( npts );

            _quad_rules[ ncommon ][ order ].resize( padded_npts );

            size_t  i = 0;

            // copy "real" quadrature data
            for ( ; i < npts; ++i )
            {
                _quad_rules[ ncommon ][ order ].assign( i, rule[i].pttri1, rule[i].pttri2, rule[i].wght );
                
                // _quad_rules[ ncommon ][ order ].pts1[i]  = rule[i].pttri1;
                // _quad_rules[ ncommon ][ order ].pts2[i]  = rule[i].pttri2;
                // _quad_rules[ ncommon ][ order ].w_sgl[i] = rule[i].wght;
                // _quad_rules[ ncommon ][ order ].w_dbl[i] = rule[i].wght;

                // _quad_rules[ ncommon ][ order ].x1[i]    = rule[i].pttri1[0];
                // _quad_rules[ ncommon ][ order ].y1[i]    = rule[i].pttri1[1];
                // _quad_rules[ ncommon ][ order ].x2[i]    = rule[i].pttri2[0];
                // _quad_rules[ ncommon ][ order ].y2[i]    = rule[i].pttri2[1];
            }// for

            // fill rest with zero (only weights to avoid division by zero in coordinate computations)
            for ( ; i < padded_npts; ++i )
            {
                _quad_rules[ ncommon ][ order ].assign( i, rule[0].pttri1, rule[0].pttri2, 0.0 );

                // _quad_rules[ ncommon ][ order ].pts1[i]  = rule[0].pttri1;
                // _quad_rules[ ncommon ][ order ].pts2[i]  = rule[0].pttri2;
                // _quad_rules[ ncommon ][ order ].w_sgl[i] = 0.0;
                // _quad_rules[ ncommon ][ order ].w_dbl[i] = 0.0;

                // _quad_rules[ ncommon ][ order ].x1[i]    = rule[0].pttri1[0];
                // _quad_rules[ ncommon ][ order ].y1[i]    = rule[0].pttri1[1];
                // _quad_rules[ ncommon ][ order ].x2[i]    = rule[0].pttri2[0];
                // _quad_rules[ ncommon ][ order ].y2[i]    = rule[0].pttri2[1];
            }// for
        }// for
    }// for
}

//////////////////////////////////////
//
// evaluate bilinearform
//

//
// evaluate subblock defined by \a row_ind × \a col_ind; the indices
// in \a row_ind and \a col_ind can be arbitrary, e.g. must not be
// contiguous
//
template < typename ansatzsp_t,
           typename testsp_t,
           typename value_t >
void
TQuadBEMBF< ansatzsp_t, testsp_t, value_t >::eval  ( const std::vector< idx_t > &,
                                                     const std::vector< idx_t > &,
                                                     BLAS::Matrix< value_t > & ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

template < typename ansatzsp_t,
           typename testsp_t,
           typename value_t >
void
TQuadBEMBF< ansatzsp_t, testsp_t, value_t >::eval  ( const std::vector< idx_t > &,
                                                     const std::vector< idx_t > &,
                                                     BLAS::Matrix< value_t > &,
                                                     const double ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// reorder triangle vertices such that common vertices are ordered
// from 0,…,#common-1; number of common vertices is returned
//
template < typename ansatzsp_t,
           typename testsp_t,
           typename value_t >
uint
TQuadBEMBF< ansatzsp_t, testsp_t, value_t >::reorder_common ( idx_t *  vtx0idxs,
                                                              idx_t *  vtx1idxs ) const
{
    //
    // count number of common vertices
    //
    
    uint  ncommon = 0;

    for ( uint i = 0; i < 3; i++ )
    {
        for ( uint j = 0; j < 3; j++ )
        {
            if ( vtx0idxs[i] == vtx1idxs[j] )
            {
                std::swap( vtx0idxs[ncommon], vtx0idxs[i] );
                std::swap( vtx1idxs[ncommon], vtx1idxs[j] );
                ncommon++;
                break;
            }// if
        }// for
    }// for
    
    if ( ncommon > 3 )
        HERROR( ERR_CONSISTENCY, "(TQuadBEMBF) reorder_common",
               "triangles have more than 3 common vertices" );

    return ncommon;
}

//
// adjust quadrature order depending on diameter and distance of triangles
//
template < typename ansatzsp_t,
           typename testsp_t,
           typename value_t >
uint
TQuadBEMBF< ansatzsp_t, testsp_t, value_t >::adjust_order ( const idx_t *  vtx0idxs,
                                                            const idx_t *  vtx1idxs,
                                                            const uint     order ) const
{
    const TGrid *  ansatz_grid = this->ansatz_space()->grid();
    const TGrid *  test_grid   = this->test_space()->grid();
    
    //
    // compute bounding box of vertex set
    //
    
    T3Point  bbmin0( ansatz_grid->vertex( vtx0idxs[0] ) );
    T3Point  bbmax0( bbmin0 );
    T3Point  bbmin1( test_grid->vertex( vtx1idxs[0] ) );
    T3Point  bbmax1( bbmin1 );

    for ( uint  j = 1; j < 3; j++ )
    {
        for ( uint  i = 0; i < 3; i++ )
        {
            bbmin0[i] = std::min( bbmin0[i], ansatz_grid->vertex( vtx0idxs[j] )[i] );
            bbmax0[i] = std::max( bbmax0[i], ansatz_grid->vertex( vtx0idxs[j] )[i] );
            bbmin1[i] = std::min( bbmin1[i], test_grid->vertex( vtx1idxs[j] )[i] );
            bbmax1[i] = std::max( bbmax1[i], test_grid->vertex( vtx1idxs[j] )[i] );
        }// for
    }// for
    
    const double  diam_t0 = norm2( bbmax0 - bbmin0 );
    const double  diam_t1 = norm2( bbmax1 - bbmin1 );
    double        dist    = 0.0;
    
#if 1

    //
    // determine distance between bounding boxes
    //
    
    for ( uint  i = 0; i < 3; i++ )
    {
        const double min0 = bbmin0[i];
        const double max0 = bbmax0[i];
        const double min1 = bbmin1[i];
        const double max1 = bbmax1[i];
        
        if      (min0 > max1) dist += Math::square( min0 - max1 );
        else if (max0 < min1) dist += Math::square( min1 - max0 );
    }// for
    
    dist = Math::sqrt( dist );
    
#else

    //
    // distance = ∥ dist( center0 - center1 ) ∥₂
    //
    
    dist  = norm2( 0.5 * ((bbmax0 + bbmin0) - (bbmax1 + bbmin1)) );
    dist -= 0.5 * ( diam_t0 + diam_t1 );
    
#endif

    //
    // inflate diameter and reduce order as long as distance is still
    // large compared to diameter
    //
    
    double  diam   = 5.0 * std::max( diam_t0, diam_t1 );
    uint    torder = order;
    
    while (( dist > diam ) && ( torder > 1 ))
    {
        diam *= 5.0;
        torder--;
    }// while

    return torder;
}    

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instantiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

#define INST_ALL( type1, type2 )                                        \
    template class TQuadBEMBF< TConstFnSpace< type1 >,  TConstFnSpace< type1 >,  type2 >; \
    template class TQuadBEMBF< TConstFnSpace< type1 >,  TLinearFnSpace< type1 >, type2 >; \
    template class TQuadBEMBF< TLinearFnSpace< type1 >, TConstFnSpace< type1 >,  type2 >; \
    template class TQuadBEMBF< TLinearFnSpace< type1 >, TLinearFnSpace< type1 >, type2 >;

INST_ALL( float,  float )
INST_ALL( double, double )
INST_ALL( float,  std::complex< float > )
INST_ALL( double, std::complex< double > )

template class TQuadBEMBF< TConstEdgeFnSpace, TConstEdgeFnSpace, std::complex< float > >;
template class TQuadBEMBF< TConstEdgeFnSpace, TConstEdgeFnSpace, std::complex< double > >;

}// namespace
