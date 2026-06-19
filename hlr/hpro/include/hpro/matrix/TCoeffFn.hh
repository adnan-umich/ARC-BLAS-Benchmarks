#ifndef __HPRO_TCOEFFFN_HH
#define __HPRO_TCOEFFFN_HH
//
// Project     : HLIBpro
// File        : TCoeffFn.hh
// Description : represents a function for matrix coefficients
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>
#include <memory>

#include "hpro/base/types.hh"
#include "hpro/base/traits.hh"

#include "hpro/cluster/TIndexSet.hh"
#include "hpro/cluster/TPermutation.hh"

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TDenseMatrix.hh"

namespace Hpro
{
    
////////////////////////////////////////////////////////////////////
//!
//! \class  TCoeffFn
//! \brief  Base class for returning coefficient for a given indexpair
//!         (i,j) in _internal_ ordering
//!
template < typename T >
class TCoeffFn
{
public:
    //
    // template type as public member type
    //

    using  value_t = T;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    //! default constructor
    TCoeffFn () {}

    //! destructor
    virtual ~TCoeffFn () {}

    //! return format of matrix, e.g. symmetric or hermitian
    virtual matform_t  matrix_format  () const { return MATFORM_NONSYM; }

    //! set accuracy for coefficient evaluation
    virtual void set_accuracy ( const double  /* acc */ ) {} // ignored by default
    
    //////////////////////////////////////
    //
    // evaluation of coefficient
    //

    //! return τ x σ subblock of real valued matrix as dense matrix block
    //! in column major format
    virtual void eval  ( const TIndexSet &  rowis,
                         const TIndexSet &  colis,
                         value_t *          matrix,
                         const double       eps = 0 ) const
    {
        const size_t          n = rowis.size();
        const size_t          m = colis.size();
        std::vector< idx_t >  tau( n );
        std::vector< idx_t >  sig( m );

        for ( auto i : rowis )
            tau[i-rowis.first()] = i;

        for ( auto i : colis )
            sig[i-colis.first()] = i;

        if ( CFG::Build::check_cb_ret )
        {
            for ( size_t  i = 0; i < m*n; i++ )
                matrix[i] = Limits::nan< value_t >();
        }// if
            
        eval( tau, sig, matrix, eps );

        if ( CFG::Build::check_cb_ret )
        {
            for ( size_t  i = 0; i < m*n; i++ )
                if ( Math::is_nan( matrix[i] ) )
                {
                    HERROR( ERR_NAN, "(TCoeffFn) eval", "matrix coefficient not set" );
                }// if
        }// if
    }

    //! return matrix coefficient for index positions defined by \a tau
    //! and \a sigma, which may be not be consecutively numbered
    virtual void  eval   ( const std::vector< idx_t > &  tau,
                           const std::vector< idx_t > &  sigma,
                           value_t *                     matrix,
                           const double                  eps = 0 ) const = 0;

    //////////////////////////////////////
    //
    // matrix block construction
    //

    //
    // construct dense matrix for \a rowis × \a colis
    //
    virtual
    std::unique_ptr< TMatrix< value_t > >
    build ( const TIndexSet &  rowis,
            const TIndexSet &  colis ) const
    {
        auto  M = std::make_unique< TDenseMatrix< value_t > >( rowis, colis );

        eval( rowis, colis, blas_mat< value_t >( M.get() ).data() );

        return M;
    }
};

////////////////////////////////////////////////////////////////////
//!
//! \class TPermCoeffFn
//! \brief Eval coefficient function with reordered indices, e.g.
//!        change internal to external ordering
//!
template < typename T >
class TPermCoeffFn : public TCoeffFn< T >
{
public:
    //
    // template type as public member type
    //

    using  value_t = T;
    using  coeff_t = TCoeffFn< value_t >;
    
protected:
    //! @cond

    // base coefficient function to call with permuted indices
    coeff_t *                         _base_coeff_fn;

    // manages coefficient function object in case of shared ownership
    const std::shared_ptr< coeff_t >  _base_coeff_ptr;
    
    // permutations to convert from internal to external numbering
    const TPermutation *              _row_perm_i2e;
    const TPermutation *              _col_perm_i2e;

    //! @endcond

public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TPermCoeffFn ( coeff_t *                           base_coeff,
                   const TPermutation *                row_perm,
                   const TPermutation *                col_perm )
            : _base_coeff_fn( base_coeff )
            , _row_perm_i2e( row_perm )
            , _col_perm_i2e( col_perm )
    {
        if ( _base_coeff_fn == nullptr )
            HERROR( ERR_ARG, "(TPermCoeffFn) ctor", "base coeff. function is NULL" );
    }

    TPermCoeffFn ( const std::shared_ptr< coeff_t > &  base_coeff,
                   const TPermutation *                row_perm,
                   const TPermutation *                col_perm )
            : _base_coeff_fn( base_coeff.get() )
            , _base_coeff_ptr( base_coeff )
            , _row_perm_i2e( row_perm )
            , _col_perm_i2e( col_perm )
    {
        if ( _base_coeff_fn == nullptr )
            HERROR( ERR_ARG, "(TPermCoeffFn) ctor", "base coeff. function is NULL" );
    }

    virtual ~TPermCoeffFn () {}

    //! return format of matrix, e.g. symmetric or hermitian
    virtual matform_t  matrix_format  () const { return _base_coeff_fn->matrix_format(); }
        
    //! set accuracy for coefficient evaluation
    virtual void set_accuracy ( const double  acc )
    {
        _base_coeff_fn->set_accuracy( acc );
    }
    
    //////////////////////////////////////
    //
    // evaluation of coefficient
    //

    // return tau x sigma subblock of matrix as dense matrix block
    // in column major format
    virtual void eval ( const TIndexSet &  rowis,
                        const TIndexSet &  colis,
                        value_t *          matrix,
                        const double       eps = 0 ) const
    {
        const auto  ptau = permute_row( rowis );
        const auto  psig = permute_col( colis );

        _base_coeff_fn->eval( ptau, psig, matrix, eps );

        if ( CFG::Build::check_cb_ret )
        {
            for ( size_t  i = 0; i < rowis.size()*colis.size(); i++ )
                if ( Math::is_nan( matrix[i] ) )
                {
                    HERROR( ERR_NAN, "(TPermCoeffFn) eval", "matrix coefficient not set" );
                }// if
        }// if
    }

    //! return matrix coefficient for index positions defined by \a tau
    //! and \a sigma, which may be not be consecutively numbered
    virtual void  eval   ( const std::vector< idx_t > &  rowidxs,
                           const std::vector< idx_t > &  colidxs,
                           value_t *                     matrix,
                           const double                  eps = 0 ) const
    {
        const auto  ptau = permute_row( rowidxs );
        const auto  psig = permute_col( colidxs );

        _base_coeff_fn->eval( ptau, psig, matrix, eps );

        if ( CFG::Build::check_cb_ret )
        {
            for ( size_t  i = 0; i < rowidxs.size()*colidxs.size(); i++ )
                if ( Math::is_nan( matrix[i] ) )
                {
                    HERROR( ERR_NAN, "(TPermCoeffFn) eval", "matrix coefficient not set" );
                }// if
        }// if
    }

    //
    // permute index sets
    //

    std::vector< idx_t >
    permute_row ( const TIndexSet & is ) const
    {
        std::vector< idx_t >  pidxs( is.size() );

        if ( _row_perm_i2e != nullptr )
        {
            for ( auto  idx : is )
                pidxs[ idx - is.first() ] = idx_t( _row_perm_i2e->permute( idx ) );
        }// if
        else
        {
            for ( auto  idx : is )
                pidxs[ idx - is.first() ] = idx;
        }// else

        return pidxs;
    }
    
    std::vector< idx_t >
    permute_row ( const std::vector< idx_t > &  idxs ) const
    {
        const size_t          n = idxs.size();
        std::vector< idx_t >  pidxs( n );

        if ( _row_perm_i2e != nullptr )
        {
            for ( size_t  i = 0; i < n; ++i )
                pidxs[ i ] = idx_t( _row_perm_i2e->permute( idxs[ i ] ) );
        }// if
        else
        {
            pidxs = idxs;
        }// else

        return pidxs;
    }
    
    std::vector< idx_t >
    permute_col ( const TIndexSet & is ) const
    {
        std::vector< idx_t >  pidxs( is.size() );

        if ( _col_perm_i2e != nullptr )
        {
            for ( auto  idx : is )
                pidxs[ idx - is.first() ] = idx_t( _col_perm_i2e->permute( idx ) );
        }// if
        else
        {
            for ( auto  idx : is )
                pidxs[ idx - is.first() ] = idx;
        }// else

        return pidxs;
    }
    
    std::vector< idx_t >
    permute_col ( const std::vector< idx_t > &  idxs ) const
    {
        const size_t          n = idxs.size();
        std::vector< idx_t >  pidxs( n );

        if ( _col_perm_i2e != nullptr )
        {
            for ( size_t  i = 0; i < n; ++i )
                pidxs[ i ] = idx_t( _col_perm_i2e->permute( idxs[ i ] ) );
        }// if
        else
        {
            pidxs = idxs;
        }// else

        return pidxs;
    }
    
    DISABLE_COPY_OP( TPermCoeffFn );
};
    
////////////////////////////////////////////////////////////////////
//!
//! \class TCBCoeffFn
//! \brief Eval real valued matrix coefficients with call-back function.
//!
template < typename T >
class TCBCoeffFn : public TCoeffFn< T >
{
public:
    //
    // template type as public member type
    //

    using  value_t = T;
    
    //
    // local function type
    //

    using  func_t = void (*) ( const size_t,
                               const int *,
                               const size_t,
                               const int *,
                               value_t *,
                               void * );

protected:
    // function to call for a coefficient
    const func_t         _fn;
    void *               _fn_arg;

    // matrix format of coefficients
    const matform_t      _matrix_format;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TCBCoeffFn ( const func_t     fn,
                 void *           arg,
                 const matform_t  mat_format )
            : _fn( fn )
            , _fn_arg( NULL )
            , _matrix_format( mat_format )
    {
        _fn_arg = arg;

        if ( _fn == NULL )
            HERROR( ERR_ARG, "(TCBCoeffFn)", "NULL function supplied" );
    }

    virtual ~TCBCoeffFn () {}

    //! return format of matrix, e.g. symmetric or hermitian
    virtual matform_t  matrix_format  () const { return _matrix_format; }
        
    //////////////////////////////////////
    //
    // evaluation of coefficient
    //

    // return tau x sigma subblock of matrix as dense matrix block
    // in column major format
    virtual void eval ( const std::vector< idx_t > &  rowidxs,
                        const std::vector< idx_t > &  colidxs,
                        value_t *                     matrix,
                        const double                  eps = 0 ) const
    {
        //
        // convert to int and call callback function
        //
            
        const size_t        n = rowidxs.size();
        const size_t        m = colidxs.size();
        std::vector< int >  itau( n );
        std::vector< int >  isig( m );

        for ( size_t  i = 0; i < n; ++i )
            itau[i] = int( rowidxs[i] );

        for ( size_t  i = 0; i < m; ++i )
            isig[i] = int( colidxs[i] );

        if ( CFG::Build::check_cb_ret )
        {
            for ( size_t  i = 0; i < m*n; i++ )
                matrix[i] = Limits::nan< value_t >();
        }// if
            
        _fn( itau.size(), itau.data(), isig.size(), isig.data(), matrix, _fn_arg );

        if ( CFG::Build::check_cb_ret )
        {
            for ( size_t  i = 0; i < m*n; i++ )
                if ( Math::is_nan( matrix[i] ) )
                {
                    HERROR( ERR_NAN, "(TCBCoeffFn) eval", "matrix coefficient not set" );
                }// if
        }// if
    }

    using TCoeffFn< value_t >::eval;
    
    DISABLE_COPY_OP( TCBCoeffFn );
};
    
////////////////////////////////////////////////////////////////////
//
// return matrix coefficients from given dense matrix in
// column major format
//
template < typename T >
class TDenseCoeffFn : public TCoeffFn< T >
{
public:
    //
    // template type as public member type
    //

    using  value_t = T;
    
protected:

    //! @cond

    // index offsets for rows/columns
    idx_t                    _row_ofs, _col_ofs;
    
    // dense matrix
    BLAS::Matrix< value_t >  _matrix;

    // matrix format
    matform_t                _matrix_format;
    
    //! @endcond
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct coefficient function with \a M as dense matrix
    //! holding all coefficients
    TDenseCoeffFn ( const BLAS::Matrix< value_t > &  M,
                    const matform_t                  mat_format )
            : _row_ofs( 0 )
            , _col_ofs( 0 )
            , _matrix( M )
            , _matrix_format( mat_format )
    {}

    //! construct coefficient function with \a M as dense matrix
    //! holding all coefficients
    TDenseCoeffFn ( const TDenseMatrix< value_t > *  M )
            : TCoeffFn< value_t >()
    {
        if ( M == nullptr )
            HERROR( ERR_ARG, "TDenseCoeffFn", "M is null" );

        _row_ofs       = M->row_ofs();
        _col_ofs       = M->col_ofs();
        _matrix        = blas_mat( M );
        _matrix_format = M->form();
    }

    virtual ~TDenseCoeffFn () {}

    //! return format of matrix, e.g. symmetric or hermitian
    virtual matform_t  matrix_format  () const
    {
        return _matrix_format;
    }
    
    //////////////////////////////////////
    //
    // evaluation of coefficient
    //

    // return \a rowidxs x \a colidxs subblock of matrix as dense matrix block
    // in column major format
    virtual void eval ( const std::vector< idx_t > &  rowidxs,
                        const std::vector< idx_t > &  colidxs,
                        value_t *                     matrix,
                        const double                  eps = 0 ) const
    {
        const size_t  n = rowidxs.size();
        const size_t  m = colidxs.size();
            
        for ( size_t  j = 0; j < m; ++j )
        {
            const idx_t  ex_j = colidxs[ j ];
                
            for ( size_t  i = 0; i < n; ++i )
            {
                const idx_t  ex_i = rowidxs[ i ];
                    
                matrix[ ( j * n ) + i ] = _matrix( ex_i - _row_ofs, ex_j - _col_ofs );
            }// for
        }// for
    }

    using TCoeffFn< value_t >::eval;
    
    DISABLE_COPY_OP( TDenseCoeffFn );
};
    
}// namespace

#endif // __HPRO_TCOEFFFN_HH
