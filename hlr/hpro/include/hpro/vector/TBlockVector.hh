#ifndef __HPRO_TBLOCKVECTOR_HH
#define __HPRO_TBLOCKVECTOR_HH
//
// Project     : HLIBpro
// File        : TBlockVector.hh
// Description : class for a vector build out of other vectors
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/base/System.hh"
#include "hpro/vector/TVector.hh"

namespace Hpro
{

// local matrix type
DECLARE_TYPE( TBlockVector );

//!
//! \ingroup Vector_Module
//! \class   TBlockVector
//! \brief   Class for a blocked, scalar vector
//!
template < typename T_value >
class TBlockVector : public TVector< T_value >
{
public:
    using  value_t = T_value;
    using  real_t  = typename real_type< value_t >::type_t;
    
private:
    //! individual vector blocks
    std::vector< TVector< value_t > * >  _blocks;

public:
    
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct block vector with \a nb sub blocks
    TBlockVector ( const uint  nb = 0 )
    {
        _blocks.resize( nb, NULL );
    }

    //! construct block vector over index set \a ais and
    //! with subblocks \a ablocks
    TBlockVector ( const TIndexSet &                            ais,
                   const std::vector< TVector< value_t > * > &  ablocks )
            : TVector< value_t >( ais.first() )
            , _blocks( ablocks )
    {}

    //! dtor
    virtual ~TBlockVector ();

    ///////////////////////////////////////////////
    //
    // access internal data
    //

    //! return size of vector
    virtual size_t              size              () const;

    //! return number of blocks
    virtual uint                n_blocks          () const { return uint(_blocks.size()); }

    //! access single vector block
    TVector< value_t > *        block             ( const uint i )       { return _blocks[i]; }

    //! access single vector block
    const TVector< value_t > *  block             ( const uint i ) const { return _blocks[i]; }

    //! set single vector block
    void                        set_block         ( const uint i, TVector< value_t > * v );
    
    //! setup block structure of vector
    void                        set_block_struct  ( const uint i );

public:
    ////////////////////////////////////////////////
    //
    // access entries
    //

    //! return i'th entry
    virtual value_t  entry      ( const idx_t  i ) const;

    //! set i'th entry
    virtual void     set_entry  ( const idx_t  i, const value_t   f );

    //////////////////////////////////////////////////
    //
    // BLAS-routines
    //

    //! conjugate entries
    virtual void     conjugate ();
        
    //! fill with constant
    virtual void     fill   ( const value_t f );
    
    //! fill with random numbers
    virtual void     fill_rand ( const uint seed );

    //! set this ≔ this + α · x
    virtual void     axpy   ( const value_t alpha, const TVector< value_t > *  x );

    //! set this ≔ α · this
    virtual void     scale  ( const value_t alpha );

    //! set this ≔ α · x
    virtual void     assign ( const value_t alpha, const TVector< value_t > *  x );

    // return dot product <this,x> = this^H x
    virtual value_t  dot    ( const TVector< value_t > *  x ) const;

    // return dot product <this,x> = this^T x
    virtual value_t  dotu   ( const TVector< value_t > *  x ) const;

    //! compute ‖·‖₂
    virtual real_t   norm2  () const { return Math::sqrt( std::real( dot( this ) ) ); }

    //! compute ‖·‖∞
    virtual real_t   norm_inf () const;

    ///////////////////////////////////////////////
    //
    // misc.
    //

    //
    // memory consumption
    //
    
    //! return size in bytes used by this object
    virtual size_t byte_size () const;

    //
    // virtual constructor
    //

    //! return copy of vector
    virtual auto  copy  () const -> std::unique_ptr< TVector< value_t > >;
    
    //! return object of same class
    virtual auto  create () const -> std::unique_ptr< TVector< value_t > > { return std::make_unique< TBlockVector< value_t > >(); }

    //
    // output
    //

    //! print vector to stdout
    virtual void print ( const uint ofs = 0 ) const;

    
    HPRO_RTTI_DERIVED( TBlockVector, TVector< value_t > )
};

}// namespace

#endif  // __HPRO_TBLOCKVECTOR_HH
