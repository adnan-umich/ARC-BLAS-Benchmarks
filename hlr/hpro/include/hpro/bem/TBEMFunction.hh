#ifndef __HPRO_TBEMFUNCTION_HH
#define __HPRO_TBEMFUNCTION_HH
//
// Project     : HLIBpro
// File        : TBEMFunction.hh
// Description : classes for evaluating functions on grids
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/types.hh"
#include "hpro/base/TPoint.hh"

namespace Hpro
{
    
///////////////////////////////////////////////////////////////
//
// function which can be evaluated at grid points
//
///////////////////////////////////////////////////////////////

template < typename T_val >
class TBEMFunction
{
public:
    //! template argument as internal type
    using  value_t = T_val;

public:
    ///////////////////////////////////////////////////////////////
    //
    // destructor
    //
    
    virtual ~TBEMFunction () {}

    //! return true if function is complex valued and false otherwise
    bool is_complex () const { return is_complex_type< value_t >::value; }
        
    ///////////////////////////////////////////////////////////////
    //
    // evaluate function at <x> and return real or complex
    // valued result; <n> is the (unit) outward normal direction
    // at point <x>
    //
        
    virtual value_t  eval  ( const T3Point &  x,
                             const T3Point &  n ) const = 0;
};

///////////////////////////////////////////////////////////////
//
// wrapper for callback functions
//
///////////////////////////////////////////////////////////////

template < typename T_val >
class TBEMCBFunction : public TBEMFunction< T_val >
{
public:
    //! template argument as internal type
    using  value_t = T_val;
    
    // callback function types
    using  func_t  = void (*) ( const double *  x,
                                const double *  n,
                                void *          arg,
                                value_t *       res );
        
private:
    //! callback function
    const func_t  _fn;

    //! user argument to callback function
    void *        _arg;

private:
    // prevent default constructor
    TBEMCBFunction ();
        
public:
    ///////////////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TBEMCBFunction ( const func_t  afn,
                     void *        aarg )
            : _fn(afn), _arg(aarg)
    {
        if ( _fn == NULL )
            HERROR( ERR_ARG, "(TBEMCBFunction)", "invalid function supplied (NULL)" );
    }

    virtual ~TBEMCBFunction () {}

    ///////////////////////////////////////////////////////////////
    //!
    //! evaluate function at \a x and return real or complex
    //! valued result; \a n is the (unit) outward normal direction
    //! at point \a x
    //!
    virtual value_t eval ( const T3Point &  x,
                           const T3Point &  n ) const
    {
        value_t  res = value_t(0);

        if ( CFG::Build::check_cb_ret )
            res = Limits::nan<value_t>();
            
        _fn( x.vector(), n.vector(), _arg, & res );

        if ( CFG::Build::check_cb_ret && Math::is_nan( res ) )
            HERROR( ERR_NAN, "(TBEMCBFunction) eval", "result not set" );
        
        return res;
    }

    DISABLE_COPY_OP( TBEMCBFunction );
};
    
}// namespace

#endif  // __HPRO_TBEMFUNCTION_HH
