#ifndef __HPRO_TPOINT_HH
#define __HPRO_TPOINT_HH
//
// Project     : HLIBpro
// File        : TPoint.hh
// Description : class for a n-dimensional vector
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <vector>

#include "hpro/base/types.hh"
#include "hpro/base/error.hh"
#include "hpro/base/System.hh"

namespace Hpro
{

//////////////////////////////////////////////////////////
//
// vector in euclidean space named as "point" since
// vectors already have different meaning in HLIBpro
//
class TPoint : public std::vector< double >
{
public:
    using  value_t = double;
    
public:
    /////////////////////////////////////////
    //
    // constructor and destructor
    //

    TPoint ( const uint  adim = 0 )
            : std::vector< double >( adim )
    {}

    TPoint ( const uint  adim, const double * coord )
            : std::vector< double >( adim )
    {
        if (( adim > 0 ) && ( coord == nullptr ))
            HERROR( ERR_ARG, "(TPoint)", "coordinate array is NULL" );
        
        for ( uint  i = 0; i < adim; i++ )
            (*this)[i] = coord[i];
    }
    
    TPoint ( const double x0 )
            : std::vector< double >{ x0 }
    {}
    
    TPoint ( const double x0,
             const double x1 )
            : std::vector< double >{ x0, x1 }
    {}
    
    TPoint ( const double x0,
             const double x1,
             const double x2 )
            : std::vector< double >{ x0, x1, x2 }
    {}
    
    TPoint ( const double x0,
             const double x1,
             const double x2,
             const double x3 )
            : std::vector< double >{ x0, x1, x2, x3 }
    {}
    
    TPoint ( const double x0,
             const double x1,
             const double x2,
             const double x3,
             const double x4 )
            : std::vector< double >{ x0, x1, x2, x3, x4 }
    {}
    
    TPoint ( const TPoint & v )
            : std::vector< double >( v )
    {}

    ~TPoint () {}
    
    /////////////////////////////////////////
    //
    // access to coordinates
    //

    uint   dim () const  { return uint(size()); }
    
    TPoint & set_dim ( const uint  n )
    {
        resize( size_t(n) );

        return *this;
    }
    
    double *        vector ()       { return this->data(); }
    const double *  vector () const { return this->data(); }

    /////////////////////////////////////////
    //
    // linear algebra
    //

    //
    // this = f * (1,1,...,1)
    //
    TPoint & fill ( const double f )
    {
        for ( auto &  x : *this )
            x = f;
        
        return *this;
    }

    //
    // assign 2- and 3-dimensional vectors
    //
    TPoint & assign ( const double x, const double y )
    {
        (*this) = { x, y };
        
        return *this;
    }
    TPoint & assign ( const double x, const double y, const double z )
    {
        (*this) = { x, y, z };
        
        return *this;
    }

    //
    // assign linear combinations
    //
    TPoint & assign ( const double f, const TPoint & v )
    {
        if ( dim() != v.dim() ) set_dim( v.dim() );

        for ( uint  i = 0; i < dim(); i++ )
            (*this)[i] = f * v[i];
    
        return *this;
    }
    TPoint & assign ( const double f1, const TPoint & v1,
                      const double f2, const TPoint & v2  )
    {
        if ( v1.dim() != v2.dim() )
            HERROR( ERR_DIM, "(TPoint) assign", "" );

        if ( dim() != v1.dim() ) set_dim( v1.dim() );

        for ( uint  i = 0; i < dim(); i++ )
            (*this)[i] = f1 * v1[i] + f2 * v2[i];
    
        return *this;
    }
    
    //
    // add linear combinations
    //
    TPoint & add ( const double f, const TPoint & v )
    {
        if ( dim() != v.dim() )
            HERROR( ERR_DIM, "(TPoint) assign", "" );

        for ( uint i = 0; i < dim(); i++ )
            (*this)[i] += f * v[i];
        
        return *this;
    }
    TPoint & add ( const double f1, const TPoint & v1,
                   const double f2, const TPoint & v2 )
    {
        if (( dim() != v1.dim() ) || ( dim() != v2.dim() ))
            HERROR( ERR_DIM, "(TPoint) assign", "" );

        for ( uint i = 0; i < dim(); i++ )
            (*this)[i] += f1 * v1[i] + f2 * v2[i];
        
        return *this;
    }

    //
    // this = a * this
    //
    TPoint & scale ( const double f )
    {
        for ( uint i = 0; i < dim(); i++ )
            (*this)[i] *= f;
        
        return *this;
    }

    //
    // return inner product ( this * v )
    //
    double dot ( const TPoint & v ) const
    {
        double  f = 0.0;
        uint  s = std::min( dim(), v.dim() );
        
        for ( uint i = 0; i < s; i++ )
            f += (*this)[i] * v[i];

        return f;
    }

    //
    // return euclidean norm
    //
    double norm2 () const
    {
        return Math::sqrt( dot( *this ) );
    }

    //
    // normalise vector w.r.t. euclidean norm
    //
    TPoint & normalise2 ()
    {
        double  f = norm2();

        if ( f != 0.0 )
            scale( 1.0 / f );
        else
            HERROR( ERR_DIV_ZERO, "(TPoint) normalise2", "" );

        return *this;
    }
    
    //
    // cross product for 2D and 3D
    //
    TPoint & assign_cross ( const TPoint & v0, const TPoint & v1 )
    {
        if ( v0.dim() != v1.dim() )
            HERROR( ERR_DIM, "(TPoint) assign_cross", "" );

        set_dim( 3 );
        
        if ( v0.dim() == 2 )
        {
            (*this)[0] = 0.0;
            (*this)[1] = 0.0;
            (*this)[2] = v0[0]*v1[1] - v0[1]*v1[0];
        }// if
        else if ( v0.dim() == 3 )
        {
            (*this)[0] = v0[1]*v1[2] - v0[2]*v1[1];
            (*this)[1] = v0[2]*v1[0] - v0[0]*v1[2];
            (*this)[2] = v0[0]*v1[1] - v0[1]*v1[0];
        }// if
        else
            HERROR( ERR_DIM, "(TPoint) assign_cross", "only 2 and 3 dimensional vectors supported" );

        return *this;
    }

    //////////////////////////////////////////////////
    //
    // misc. methods
    //

    //
    // assignment
    //
    TPoint & operator = ( const TPoint & v ) { return assign( 1.0, v ); }
                                                   
    //
    // compare operators
    //
    bool operator == ( const TPoint & v ) const
    {
        if ( dim() != v.dim() )
            return false;
        
        for ( uint i = 0; i < dim(); i++ )
            if ( (*this)[i] != v[i] )
                return false;
        
        return true;
    }
    
    bool operator != ( const TPoint & v ) const { return ! ((*this) == v); }

    //
    // return memory footprint of object
    //
    size_t byte_size () const
    {
        return sizeof(double*) + sizeof(size_t) + sizeof(double*) * dim();
    }
    
    //
    // stream output
    //
    std::string to_string () const;
};

//
// some usual functions
//
inline double norm2 ( const TPoint & v )                     { return v.norm2(); }
inline double dot   ( const TPoint & v1, const TPoint & v2 ) { return v1.dot( v2 ); }

inline TPoint
operator - ( const TPoint & v )
{
    TPoint  t( v );

    t.scale( -1 );
    
    return t;
}

inline TPoint
operator + ( const TPoint & v1, const TPoint & v2 )
{
    TPoint  t( v1 );

    t.add( 1.0, v2 );
    
    return t;
}

inline TPoint
operator - ( const TPoint & v1, const TPoint & v2 )
{
    TPoint  t( v1 );

    t.add( -1.0, v2 );
    
    return t;
}

inline TPoint
operator * ( const TPoint & v, const double f )
{
    TPoint  t( v );

    t.scale( f );
    
    return t;
}

inline TPoint
operator * ( const double f, const TPoint & v )
{
    TPoint  t( v );

    t.scale( f );
    
    return t;
}

//////////////////////////////////////////////////////////
//
// same as above but 3d version
//
class T3Point
{
public:
    using  value_t = double;
    
private:
    // the coordinates
    double  _data[3];
    
public:
    /////////////////////////////////////////
    //
    // constructor and destructor
    //

    T3Point ()
    {
        _data[0] = _data[1] = _data[2] = 0.0;
    }

    T3Point ( const double  ax,
              const double  ay,
              const double  az )
    {
        _data[0] = ax;
        _data[1] = ay;
        _data[2] = az;
    }
    
    T3Point ( const double *  v )
    {
        _data[0] = v[0];
        _data[1] = v[1];
        _data[2] = v[2];
    }
    
    T3Point ( const T3Point &  v )
    {
        *this = v;
    }

    explicit
    T3Point ( const TPoint &  v )
    {
        *this = v;
    }

    /////////////////////////////////////////
    //
    // access to coordinates
    //

    uint dim () const  { return 3; }
    
    double *        vector ()       { return _data; }
    const double *  vector () const { return _data; }

    double &  operator [] ( const uint  i )       { return _data[i]; }
    double    operator [] ( const uint  i ) const { return _data[i]; }

    double  x () const { return _data[0]; }
    double  y () const { return _data[1]; }
    double  z () const { return _data[2]; }
    
    /////////////////////////////////////////
    //
    // linear algebra
    //

    //
    // addition
    //
    void  operator += ( const T3Point & v )
    {
        _data[0] += v._data[0];
        _data[1] += v._data[1];
        _data[2] += v._data[2];
    }
    
    void  operator += ( const TPoint & v )
    {
        if ( v.dim() != dim() )
            HERROR( ERR_DIM, "(T3Point) operator +=", "given point is not 3D" );
        
        _data[0] += v[0];
        _data[1] += v[1];
        _data[2] += v[2];
    }
    
    //
    // subtraction
    //
    void  operator -= ( const T3Point & v )
    {
        _data[0] -= v._data[0];
        _data[1] -= v._data[1];
        _data[2] -= v._data[2];
    }
    
    void  operator -= ( const TPoint & v )
    {
        if ( v.dim() != dim() )
            HERROR( ERR_DIM, "(T3Point) operator -=", "given point is not 3D" );
        
        _data[0] -= v[0];
        _data[1] -= v[1];
        _data[2] -= v[2];
    }
    
    //
    // multiplication
    //
    void  operator *= ( const double f )
    {
        _data[0] *= f;
        _data[1] *= f;
        _data[2] *= f;
    }

    // elementwise multiplication; _not_ the dot product!
    void  operator *= ( const T3Point &  v )
    {
        _data[0] *= v[0];
        _data[1] *= v[1];
        _data[2] *= v[2];
    }
    
    //
    // return inner product ( this * v )
    //
    double  dot ( const T3Point & v ) const
    {
        return ( ( _data[0] * v._data[0] ) +
                 ( _data[1] * v._data[1] ) +
                 ( _data[2] * v._data[2] ) );
    }

    //
    // return euclidean norm
    //
    double  norm2 () const
    {
        return Math::sqrt( dot( *this ) );
    }

    //
    // normalise vector w.r.t. euclidean norm
    //
    T3Point &  normalise2 ()
    {
        const double  f = dot( *this );

        if ( f != 0.0 )
            *this *= Math::rsqrt( f );
        else
            HERROR( ERR_DIV_ZERO, "(T3Point) normalise2", "" );

        return *this;
    }
    
    //////////////////////////////////////////////////
    //
    // misc. methods
    //

    //
    // assignment
    //
    T3Point &  operator = ( const double      f )
    {
        _data[0] = f; _data[1] = f; _data[2] = f;
        return *this;
    }
    T3Point &  operator = ( const T3Point & v )
    {
        _data[0] = v._data[0]; _data[1] = v._data[1]; _data[2] = v._data[2];
        return *this;
    }
    T3Point &  operator = ( const TPoint  & v )
    {
        if ( v.dim() != dim() )
            HERROR( ERR_DIM, "(T3Point) operator =", "given point is not 3D" );
        
        _data[0] = v[0]; _data[1] = v[1]; _data[2] = v[2];
        return *this;
    }
                                                   
    //
    // compare operators
    //
    bool  operator == ( const T3Point & v ) const
    {
        if (( _data[0] != v._data[0] ) ||
            ( _data[1] != v._data[1] ) ||
            ( _data[2] != v._data[2] ))
            return false;
        else 
            return true;
    }
    
    bool  operator != ( const T3Point & v ) const { return ! ((*this) == v); }

    //
    // return memory footprint of object
    //
    size_t  byte_size () const { return 3*sizeof(double); }
    
    //
    // stream output
    //
    std::string  to_string () const;
};

//
// some usual functions and operators
//
inline double  norm2 ( const T3Point & v )                      { return v.norm2(); }
inline double  dot   ( const T3Point & v1, const T3Point & v2 ) { return v1.dot( v2 ); }

inline T3Point
operator - ( const T3Point & v )
{
    return T3Point( - v.x(), - v.y(), - v.z() );
}

inline T3Point
operator + ( const T3Point & v1, const T3Point & v2 )
{
    return T3Point( v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z() );
}

inline T3Point
operator - ( const T3Point & v1, const T3Point & v2 )
{
    return T3Point( v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z() );
}

inline T3Point
operator * ( const T3Point & v1, const double f )
{
    return T3Point( v1.x() * f, v1.y() * f, v1.z() * f );
}

inline T3Point
operator * ( const double f, const T3Point & v1 )
{
    return T3Point( v1.x() * f, v1.y() * f, v1.z() * f );
}

inline T3Point
cross ( const T3Point & v, const T3Point & w )
{
    return T3Point( (v.y() * w.z()) - (v.z() * w.y()),
                    (v.z() * w.x()) - (v.x() * w.z()),
                    (v.x() * w.y()) - (v.y() * w.x()) );
}

inline T3Point
spherical ( const double  theta,  // polar angle
            const double  phi,    // azimuthal angle
            const double  r )     // radial distance
{
    const auto  r_sin_theta = r * std::sin( theta );

    return T3Point( r_sin_theta * std::cos( phi ),
                    r_sin_theta * std::sin( phi ),
                    r * std::cos( theta ) );
}
            
//////////////////////////////////////////////////////////
//
// same as above but 2d version
//
class T2Point
{
public:
    using  value_t = double;
    
private:
    // the coordinates
    double  _data[2];
    
public:
    /////////////////////////////////////////
    //
    // constructor and destructor
    //

    T2Point ()
    {
        _data[0] = _data[1] = 0.0;
    }

    T2Point ( const double  ax,
              const double  ay )
    {
        _data[0] = ax;
        _data[1] = ay;
    }
    
    T2Point ( const T2Point & v )
    {
        *this = v;
    }

    explicit
    T2Point ( const TPoint & v )
    {
        *this = v;
    }

    /////////////////////////////////////////
    //
    // access to coordinates
    //

    uint dim () const  { return 2; }
    
    double       * vector ()       { return _data; }
    const double * vector () const { return _data; }

    double & operator [] ( const uint  i )       { return _data[i]; }
    double   operator [] ( const uint  i ) const { return _data[i]; }
    
    double  x () const { return _data[0]; }
    double  y () const { return _data[1]; }

    /////////////////////////////////////////
    //
    // linear algebra
    //

    //
    // addition
    //
    void  operator += ( const T2Point & v )
    {
        _data[0] += v._data[0];
        _data[1] += v._data[1];
    }
    
    void  operator += ( const TPoint & v )
    {
        if ( v.dim() != dim() )
            HERROR( ERR_DIM, "(T2Point) operator +=", "given point is not 3D" );
        
        _data[0] += v[0];
        _data[1] += v[1];
    }
    
    //
    // subtraction
    //
    void  operator -= ( const T2Point & v )
    {
        _data[0] -= v._data[0];
        _data[1] -= v._data[1];
    }
    
    void  operator -= ( const TPoint & v )
    {
        if ( v.dim() != dim() )
            HERROR( ERR_DIM, "(T2Point) operator -=", "given point is not 3D" );
        
        _data[0] -= v[0];
        _data[1] -= v[1];
    }
    
    //
    // multiplication
    //
    void  operator *= ( const double f )
    {
        _data[0] *= f;
        _data[1] *= f;
    }

    // elementwise multiplication; _not_ the dot product!
    void  operator *= ( const T3Point &  v )
    {
        _data[0] *= v[0];
        _data[1] *= v[1];
    }

    //
    // return inner product ( this * v )
    //
    double dot ( const T2Point & v ) const
    {
        return ( ( _data[0] * v._data[0] ) + ( _data[1] * v._data[1] ) );
    }

    //
    // return euclidean norm
    //
    double norm2 () const
    {
        return Math::sqrt( dot( *this ) );
    }

    //
    // normalise vector w.r.t. euclidean norm
    //
    T2Point & normalise2 ()
    {
        const double  f = dot( *this );

        if ( f != 0.0 )
            *this *= Math::rsqrt( f );
        else
            HERROR( ERR_DIV_ZERO, "(T2Point) normalise2", "" );

        return *this;
    }
    
    //////////////////////////////////////////////////
    //
    // misc. methods
    //

    //
    // assignment
    //
    T2Point &  operator = ( const double     f ) { _data[0] = f;     _data[1] = f;     return *this; }
    T2Point &  operator = ( const T2Point &  v ) { _data[0] = v.x(); _data[1] = v.y(); return *this; }
    T2Point &  operator = ( const TPoint  &  v ) { _data[0] = v[0];  _data[1] = v[1];  return *this; }
                                                   
    //
    // compare operators
    //
    bool operator == ( const T2Point & v ) const
    {
        if (( _data[0] != v._data[0] ) || ( _data[1] != v._data[1] ))
            return false;
        else 
            return true;
    }
    
    bool operator != ( const T2Point & v ) const { return ! ((*this) == v); }

    //
    // return memory footprint of object
    //
    size_t byte_size () const { return 2*sizeof(double); }
    
    //
    // stream output
    //
    std::string  to_string () const;
};

//
// some usual functions
//
inline double norm2 ( const T2Point & v )                      { return v.norm2(); }
inline double dot   ( const T2Point & v1, const T2Point & v2 ) { return v1.dot( v2 ); }

inline T2Point
operator - ( const T2Point & v )
{
    return T2Point( - v.x(), - v.y() );
}

inline T2Point
operator + ( const T2Point & v1, const T2Point & v2 )
{
    return T2Point( v1.x() + v2.x(), v1.y() + v2.y() );
}

inline T2Point
operator - ( const T2Point & v1, const T2Point & v2 )
{
    return T2Point( v1.x() - v2.x(), v1.y() - v2.y() );
}

inline T2Point
operator * ( const T2Point & v1, const double f )
{
    return T2Point( v1.x() * f, v1.y() * f );
}

inline T2Point
operator * ( const double f, const T2Point & v1 )
{
    return T2Point( v1.x() * f, v1.y() * f );
}

inline T2Point
spherical ( const double  phi,    // polar angle
            const double  r )     // radial distance
{
    return T2Point( r * std::cos( phi ),
                    r * std::sin( phi ) );
}
            
}// namespace Hpro

#endif  // __HPRO_TPOINT_HH
