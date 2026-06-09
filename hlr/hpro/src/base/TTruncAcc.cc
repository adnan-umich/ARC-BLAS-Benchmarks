//
// Project     : HLIBpro
// File        : TTruncAcc.cc
// Description : defines truncation accuracy for low-rank matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/cluster/TBlockCluster.hh"
#include "hpro/matrix/TMatrix.hh"

#include "hpro/base/TTruncAcc.hh"

namespace Hpro
{

// global objects
const TTruncAcc  acc_exact( 0.0, 0.0 );

namespace
{

//////////////////////////////////////////////////////////////////
//
// local functions
//
//////////////////////////////////////////////////////////////////

//
// ratio of singular-values compared to highest considered as non-zero
//
template < typename value_t >
value_t
svd_eps ()
{
    return Limits::epsilon< value_t >();
}

}// namespace anonymous

//////////////////////////////////////////////////////////////////
//
// TTruncAcc
//
//////////////////////////////////////////////////////////////////

namespace
{

template < typename value_t >
size_t
truncation_rank ( const BLAS::Vector< value_t > &  sv,
                  const TTruncAcc &                acc )
{
    // in case of empty matrix
    if ( sv.length() > 0 && ( sv(0) == value_t(0) ))
        return 0;
        
    if ( acc.norm_mode() == spectral_norm )
    {
        //
        // with fixed precision: find k such that σ_(k+1) ≤ β
        // with β = σ₀ · ε for relative error and
        //      β = ε      for absolute error
        //
        
        auto   eps = value_t(0);
        idx_t  k   = idx_t( sv.length() );

        // initialise with either fixed rank or fixed accuracy
        if ( acc.is_fixed_rank() )
        {
            eps = svd_eps< value_t >() * Math::abs( sv(0) );
            k   = idx_t( std::min( sv.length(), acc.rank() ) );
        }// if
        else
        {
            eps = value_t( acc.rel_eps() ) * Math::abs( sv(0) );
            k   = idx_t( sv.length() );
        }// else

        // apply absolute lower limit for singular values
        eps = std::max( eps, value_t( acc.abs_eps() ) );

        // apply maximal rank
        if ( acc.has_max_rank() )
            k = std::min( k, idx_t( acc.max_rank() ) );

        // compare singular values and stop, if truncation rank was reached
        for ( idx_t  i = 0; i < k; ++i )
        {
            if ( Math::abs( sv(i) ) < eps )
            {
                k = i;
                break;
            }// if
        }// for
    
        return k;
    }// if
    else // _norm_mode == frobenius_norm
    {
        //
        // with fixed precision: find smallest k such that √(Σ_i=k^n σ_i²) ≤ β
        // with β = ε |A|_F = ε √(Σ_i=1^n σ_i²) for relative error and
        //      β = ε                           for absolute error
        //
        
        auto   eps  = value_t(0);
        idx_t  k    = idx_t( sv.length() );
        auto   norm = value_t(0);

        // compute Frobenius norm of matrix as √(Σ_i σ_i²)
        for ( idx_t  i = 0; i < idx_t(sv.length()); ++i )
            norm += Math::square( sv(i) );
        norm = Math::sqrt( norm );
            
        // initialise with either fixed rank or fixed accuracy
        if ( acc.is_fixed_rank() )
        {
            eps = svd_eps< value_t >() * norm;
            k   = idx_t( std::min( sv.length(), acc.rank() ) );
        }// if
        else
        {
            eps = value_t( acc.rel_eps() ) * norm;
            k   = idx_t( sv.length() );
        }// else

        // apply absolute lower limit for singular values
        eps = std::max( eps, value_t( acc.abs_eps() ) );

        // apply maximal rank
        if ( acc.has_max_rank() )
            k = std::min( k, idx_t( acc.max_rank() ) );

        // find smallest k such that √(Σ_k^n σ_i²) ≤ ε
        auto  rest = value_t(0);

        for ( idx_t  i = k-1; i >= 0; --i )
        {
            rest += Math::square( sv(i) );
        
            if ( Math::sqrt( rest ) > eps )
                break;

            k--;
        }// for

        return k;
    }// else
}

}// namespace anonymous

//
// return truncation rank based on internal strategy
// and given singular values \a sv
//
size_t
TTruncAcc::trunc_rank  ( const BLAS::Vector< float > &   sv ) const
{
    return truncation_rank< float >(  sv, *this );
}

size_t
TTruncAcc::trunc_rank  ( const BLAS::Vector< double > &  sv ) const
{
    return truncation_rank< double >( sv, *this );
}

// instantiate the above method
// template size_t TTruncAcc::trunc_rank< float >  ( const BLAS::Vector< float > &   sv ) const;
// template size_t TTruncAcc::trunc_rank< double > ( const BLAS::Vector< double > &  sv ) const;

//
// convert to string
//
std::string
TTruncAcc::to_string () const
{
    if ( norm_mode() == spectral_norm )
    {
        if ( is_fixed_rank() )
            return Hpro::to_string( "spectral( k = %d )", rank() );
        else
            return Hpro::to_string( "spectral( ε = %.4e )", rel_eps() );
    }// if
    else
    {
        if ( is_fixed_rank() )
            return Hpro::to_string( "frobenius( k = %d )", rank() );
        else
            return Hpro::to_string( "frobenius( ε = %.4e )", rel_eps() );
    }// else
}

//////////////////////////////////////////////////////////////////
//
// TBlockTruncAcc
//
//////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////
//
// accuracy management
//

// return accuracy description for individual subblock
const TTruncAcc &
TBlockTruncAcc::acc ( const TBlockCluster *  bc ) const
{
    if ( bc == nullptr )
        HERROR( ERR_ARG, "(TBlockTruncAcc) acc", "block cluster is nullptr" );
    
    return get_acc( * bc->rowcl(), * bc->colcl() );
}

template < typename value_t >
const TTruncAcc &
TBlockTruncAcc::acc ( const TMatrix< value_t > *  M ) const
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "(TBlockTruncAcc) acc", "matrix is nullptr" );
    
    const idx_t  row_ofs = M->row_ofs();
    const idx_t  col_ofs = M->col_ofs();
    
    return get_acc( TIndexSet( row_ofs, row_ofs + idx_t(M->rows()) - 1 ),
                    TIndexSet( col_ofs, col_ofs + idx_t(M->cols()) - 1 ) );
}

const TTruncAcc
TBlockTruncAcc::acc ( const TIndexSet &  rowis,
                      const TIndexSet &  colis ) const
{
    return get_acc( rowis, colis );
}

//
// local block accuracy method
//

const TTruncAcc &
TBlockTruncAcc::get_acc ( const TIndexSet & rowis,
                          const TIndexSet & colis ) const
{
    //
    // return accuracy for given block index set
    //

    if ( _block_acc.size() == 0 )
    {
        // no block information available, just return local
        // accuracy
        return *this;
    }// if
    else
    {
        //
        // first check, if given indexset equals global indexset
        // or if it is a real subblock
        //
        
        if (( rowis.first() == _row_idx_sets.front().first() ) &&
            ( rowis.last()  == _row_idx_sets.back().last()   ) &&
            ( colis.first() == _col_idx_sets.front().first() ) &&
            ( colis.last()  == _col_idx_sets.back().last()   ) )
            return *this;
            
        //
        // look, which subblock contains cluster
        //

        size_t  col_idx = 0;

        for ( auto  col_it = _col_idx_sets.cbegin(); col_it < _col_idx_sets.cend(); ++col_it )
        {
            if ( col_it->is_sub( colis ) )
            {
                break;
            }// if

            col_idx++;
        }// for

        if ( col_idx >= _col_idx_sets.size() )
            HERROR( ERR_CONSISTENCY, "(TBlockTruncAcc) get_acc", "not a valid subblock" );
                 
        size_t  row_idx = 0;

        for ( auto  row_it = _row_idx_sets.cbegin(); row_it < _row_idx_sets.end(); ++row_it )
        {
            if ( row_it->is_sub( rowis ) )
            {
                break;
            }// if

            row_idx++;
        }// for

        if ( row_idx >= _row_idx_sets.size() )
            HERROR( ERR_CONSISTENCY, "(TBlockTruncAcc) get_acc", "not a valid subblock" );
                 

        return _block_acc[ col_idx * _row_idx_sets.size() + row_idx ];
    }// else
}

std::string
TBlockTruncAcc::to_string () const
{
    std::string  out;
    
    for ( size_t  i = 0; i < _row_idx_sets.size(); ++i )
    {
        std::string  row_is = _row_idx_sets[i].to_string();

        for ( size_t  j = 0; j < _col_idx_sets.size(); ++j )
        {
            out += row_is + " × " + _col_idx_sets[j].to_string() + " : ";
            
            out += _block_acc[ j * _row_idx_sets.size() + i ].to_string();
            out += '\n';
        }// for
    }// for

    return out;
}

}// namespace
