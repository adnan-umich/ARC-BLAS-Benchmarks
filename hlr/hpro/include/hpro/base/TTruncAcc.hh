#ifndef __HPRO_TTRUNCACC_HH
#define __HPRO_TTRUNCACC_HH
//
// Project     : HLIBpro
// File        : TTruncAcc.hh
// Description : defines truncation accuracy for low-rank matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/base/types.hh"
#include "hpro/base/traits.hh"
#include "hpro/base/config.hh"
#include "hpro/base/System.hh"
#include "hpro/cluster/TIndexSet.hh"
#include "hpro/blas/Vector.hh"

namespace Hpro
{

//
// forward decl.
//
class TBlockCluster;

template < typename value_t >
class TMatrix;

//!
//! different norms within truncation
//!
enum trunc_norm_t
{
    spectral_norm,
    frobenius_norm
};

//!
//! \class  TTruncAcc
//! \brief  Defines accuracy for truncation of low rank blocks.
//!
class TTruncAcc
{
private:
    //! @cond
    
    //! relative truncation accuracy
    double        _rel_eps;

    //! absolute truncation accuracy
    double        _abs_eps;

    //! fixed truncation rank (negative means no limit)
    int           _rank;
    
    //! upper limit for truncation rank (negative means no limit)
    int           _max_rank;

    //! norm to be used for truncation
    trunc_norm_t  _norm_mode;
    
    //! @endcond
    
public:
    /////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //!
    //! construct accuracy object for exact truncation
    //!
    TTruncAcc ()
            : _rel_eps(0.0)
            , _abs_eps(CFG::Arith::abs_eps)
            , _rank(-1)
            , _max_rank(-1)
            , _norm_mode( spectral_norm )
    {}

    //!
    //! construct accuracy object for fixed rank truncation
    //!
    TTruncAcc ( const int     k,
                const double  absolute_eps = CFG::Arith::abs_eps )
            : _rel_eps(0.0)
            , _abs_eps(absolute_eps)
            , _rank(std::max(0,k))
            , _max_rank(-1)
            , _norm_mode( spectral_norm )
    {}

    //!
    //! construct accuracy object for fixed accuracy truncation
    //! in spectral norm
    //!
    TTruncAcc ( const double  relative_eps,
                const double  absolute_eps = CFG::Arith::abs_eps )
            : _rel_eps(relative_eps)
            , _abs_eps(absolute_eps)
            , _rank(-1)
            , _max_rank(-1)
            , _norm_mode( spectral_norm )
    {}

    //!
    //! construct accuracy object for fixed accuracy truncation
    //!
    TTruncAcc ( const trunc_norm_t  anorm_mode,
                const double        arelative_eps,
                const double        aabsolute_eps = CFG::Arith::abs_eps )
            : _rel_eps(arelative_eps)
            , _abs_eps(aabsolute_eps)
            , _rank(-1)
            , _max_rank(-1)
            , _norm_mode( anorm_mode )
    {}

    //!
    //! copy constructor
    //!
    TTruncAcc ( const TTruncAcc &  ta )
            : _rel_eps(0.0)
            , _abs_eps(CFG::Arith::abs_eps)
            , _rank(-1)
            , _max_rank(-1)
            , _norm_mode( ta._norm_mode )
    {
        *this = ta;
    }

    /////////////////////////////////////////////////
    //
    // truncation rank computation
    //

    //! return truncation rank based on internal strategy
    //! and given singular values \a sv
    virtual size_t  trunc_rank  ( const BLAS::Vector< float > &   sv ) const;
    virtual size_t  trunc_rank  ( const BLAS::Vector< double > &  sv ) const;
    
    /////////////////////////////////////////////////
    //
    // accuracy management
    //

    //! return accuracy description for individual subblock
    virtual const TTruncAcc &  acc ( const TBlockCluster * ) const
    {
        return *this;
    }

    //! return accuracy description for individual submatrix
    template < typename value_t >
    const TTruncAcc &  acc ( const TMatrix< value_t > * ) const
    {
        return *this;
    }

    //! return accuracy description for individual subblock
    virtual const TTruncAcc    acc ( const TIndexSet &  /* rowis */,
                                     const TIndexSet &  /* colis */ ) const
    {
        return *this;
    }

    //! abbreviation via () operator
    const TTruncAcc &  operator () ( const TBlockCluster *  bc ) const
    {
        return acc( bc );
    }
    
    //! abbreviation via () operator
    template < typename value_t >
    const TTruncAcc &  operator () ( const TMatrix< value_t > *  M ) const
    {
        return acc( M );
    }
    
    //! return accuracy description for individual submatrix
    const TTruncAcc    operator () ( const TIndexSet &  rowis,
                                     const TIndexSet &  colis ) const
    {
        return acc( rowis, colis );
    }

    // explicit virtual destructor
    virtual ~TTruncAcc () {}
    
    /////////////////////////////////////////////////
    //
    // access accuracy data
    //

    //! return fixed rank (nonnegative!)
    size_t  rank           () const { return std::max( 0, _rank ); }

    //! return maximal rank (nonnegative!)
    size_t  max_rank       () const { return std::max( 0, _max_rank ); }

    //! return true if maximal truncation rank was defined
    bool    has_max_rank   () const { return _max_rank >= 0; }
    
    //! return relative accuracy
    double  rel_eps        () const { return _rel_eps; }

    //! return absolute accuracy
    double  abs_eps        () const { return _abs_eps; }

    //! return true if accuracy is fixed rank
    bool    is_fixed_rank  () const { return _rank >= 0; }

    //! return true if accuracy is fixed precision
    bool    is_fixed_prec  () const { return ! is_fixed_rank(); }

    //! return true if accuracy is "exact"
    bool    is_exact       () const { return ! is_fixed_rank() && (_rel_eps == 0.0) && (_abs_eps == 0.0); }

    //! set maximal rank in truncation
    void    set_max_rank   ( const int  k )
    {
        _max_rank = std::max( 0, k );
    }

    //! return norm mode of truncation
    trunc_norm_t  norm_mode () const { return _norm_mode; }
    
    //! copy operator
    TTruncAcc & operator = ( const TTruncAcc & ta )
    {
        _rank      = ta._rank;
        _max_rank  = ta._max_rank;
        _rel_eps   = ta._rel_eps;
        _abs_eps   = ta._abs_eps;
        _norm_mode = ta._norm_mode;

        return *this;
    }

    /////////////////////////////////////////////////
    //
    // misc
    //

    //! return string representation
    virtual std::string  to_string () const;
};

/////////////////////////////////////////////////
//
// functional ctors
//

//!
//! create accuracy object with fixed (relative) precision \a relative_eps
//! using spectral norm in truncation
//!
inline
TTruncAcc
fixed_prec ( const double  relative_eps,
             const double  absolute_eps = CFG::Arith::abs_eps )
{
    return TTruncAcc( relative_eps, absolute_eps );
}

//!
//! create accuracy object with fixed (relative) precision \a relative_eps
//!
inline
TTruncAcc
fixed_prec ( const trunc_norm_t  norm_mode,
             const double        relative_eps,
             const double        absolute_eps = CFG::Arith::abs_eps )
{
    return TTruncAcc( norm_mode, relative_eps, absolute_eps );
}

//!
//! create accuracy object with relative precision \a relative_eps
//! using spectral norm in truncation
//!
inline
TTruncAcc
relative_prec ( const double  relative_eps,
                const double  absolute_eps = CFG::Arith::abs_eps )
{
    return TTruncAcc( relative_eps, absolute_eps );
}

//!
//! create accuracy object with relative precision \a relative_eps
//!
inline
TTruncAcc
relative_prec ( const trunc_norm_t  norm_mode,
                const double        relative_eps,
                const double        absolute_eps = CFG::Arith::abs_eps )
{
    return TTruncAcc( norm_mode, relative_eps, absolute_eps );
}

//!
//! create accuracy object with fixed absolute precision \a absolute_eps
//! using spectral norm in truncation
//!
inline
TTruncAcc
absolute_prec ( const double  absolute_eps )
{
    return TTruncAcc( double(0), absolute_eps );
}

//!
//! create accuracy object with fixed absolute precision \a absolute_eps
//!
inline
TTruncAcc
absolute_prec ( const trunc_norm_t  norm_mode,
                const double        absolute_eps )
{
    return TTruncAcc( norm_mode, double(0), absolute_eps );
}

//!
//! create accuracy object with fixed rank \a k
//!
inline
TTruncAcc
fixed_rank ( const int     k,
             const double  absolute_eps = CFG::Arith::abs_eps )
{
    return TTruncAcc( k, absolute_eps );
}

//!
//! \class TBlockTruncAcc
//! \brief Truncation accuracy defined blockwise for block index sets.
//!
class TBlockTruncAcc : public TTruncAcc
{
private:
    //! row indexsets defining accuracy matrix
    std::vector< TIndexSet >  _row_idx_sets;

    //! column indexsets defining accuracy matrix
    std::vector< TIndexSet >  _col_idx_sets;

    //! accuracy matrix (stored column wise)
    std::vector< TTruncAcc >  _block_acc;

public:
    /////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //!
    //! construct exact accuracy object
    //!
    TBlockTruncAcc ()
    {
    }

    //!
    //! construct accuracy object for fixed rank truncation
    //!
    TBlockTruncAcc ( const int     k,
                     const double  absolute_eps = CFG::Arith::abs_eps )
            : TTruncAcc( k, absolute_eps )
    {
    }

    //!
    //! construct accuracy object for fixed accuracy truncation
    //!
    TBlockTruncAcc ( const double  relative_eps,
                     const double  absolute_eps = CFG::Arith::abs_eps )
            : TTruncAcc( relative_eps, absolute_eps )
    {
    }

    //!
    //! construct accuracy object for block wise accuracy
    //!
    TBlockTruncAcc ( const std::vector< TIndexSet > &  row_idx_sets,
                     const std::vector< TIndexSet > &  col_idx_sets,
                     const std::vector< TTruncAcc > &  block_acc )
    {
        _row_idx_sets = row_idx_sets;
        _col_idx_sets = col_idx_sets;
        _block_acc    = block_acc;
    }

    //!
    //! copy constructor
    //!
    TBlockTruncAcc ( const TBlockTruncAcc & ta )
            : TTruncAcc()
    {
        *this = ta;
    }

    /////////////////////////////////////////////////
    //
    // accuracy management
    //

    //! return accuracy description for individual subblock defined by cluster
    virtual const TTruncAcc &  acc ( const TBlockCluster *  bc ) const;

    //! return accuracy description for individual subblock defined by matrix
    template < typename value_t >
    const TTruncAcc &  acc ( const TMatrix< value_t > *  M  ) const;
    
    //! return accuracy description for individual subblock
    virtual const TTruncAcc    acc ( const TIndexSet &      rowis,
                                     const TIndexSet &      colis ) const;

    /////////////////////////////////////////////////
    //
    // access accuracy data
    //

    //! copy operator
    TBlockTruncAcc & operator = ( const TBlockTruncAcc & ta )
    {
        // call copy operator of base class
        TTruncAcc::operator = ( ta );
        
        _row_idx_sets = ta._row_idx_sets;
        _col_idx_sets = ta._col_idx_sets;
        _block_acc    = ta._block_acc;

        return *this;
    }

    /////////////////////////////////////////////////
    //
    // misc
    //

    //! return string representation
    virtual std::string  to_string () const;
    
private:
    //
    // local block accuracy method
    //

    //! return accuracy for sub block defined by \a rowis × \a colis
    const TTruncAcc &  get_acc ( const TIndexSet & rowis,
                                 const TIndexSet & colis ) const;
};

//////////////////////////////////////////////////////////
//
// predefined accuracies
//

//! global accuracy object for exact truncation:
//!   never throw away a singular value
extern const TTruncAcc  acc_exact;

//! return accuracy object for exact truncation up to machine precision
template < typename value_t >
const TTruncAcc         acc_machine ()
{
    return TTruncAcc( Limits::epsilon< real_type_t< value_t > >(), 0.0 );
}

}// namespace Hpro

#endif  // __HPRO_TTRUNCACC_HH
