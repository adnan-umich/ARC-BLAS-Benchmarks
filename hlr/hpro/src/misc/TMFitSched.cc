//
// Project     : HLIBpro
// File        : TMFitSched.cc
// Description : provides multifit scheduling algorithm
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <algorithm>

#include "hpro/base/error.hh"
#include "hpro/base/System.hh"

#include "scheduler.hh"

namespace Hpro
{

namespace
{

//
// compare object for computing sort permutations
//
template < typename T >
class TRevPermComp
{
private:
    const std::vector< T > &  _data;
    
public:

    TRevPermComp ( const std::vector< T > &  data )
            : _data( data )
    {}
    
    bool operator () ( const size_t  t1,
                       const size_t  t2 ) const
    {
        return _data[t1] > _data[t2]; // reverse!
    }
};

}// namespace anonymous

///////////////////////////////////////////////
//
// scheduling algorithm
//

//
// multifit algorithm
//
bool
TMFitSched::schedule ( const uint                     p,
                       std::vector< int > &           sched,
                       const std::vector< double > &  costs ) const
{
    //
    // handle special case
    //

    size_t  n = costs.size();
    
    sched.resize( n );

    if ( n <= p )
    {
        for ( size_t  i = 0; i < n; i++ )
            sched[i] = int(i);

        return true;
    }// if
    
    //
    // sort costs (higher first)
    //
    
    std::vector< size_t >   perm( n );
    TRevPermComp< double >  cmp( costs );

    for ( size_t  i = 0; i < n; ++i )
        perm[i] = i;

    // compute permutation for sorting
    std::sort( perm.begin(), perm.end(), cmp );

    // apply permutation to costs array
    std::vector< double >  tmp_costs( costs );

    for ( size_t  i = 0; i < n; ++i )
        tmp_costs[i] = costs[ perm[i] ];

    //
    // do binary search
    //

    std::vector< int >  best_sched( n );
    bool                found_valid = false;
    uint                proc;
    double              c, cmin, cmax;
    uint                depth = 0;

    // determine search interval
    cmax = c = 0;
    for ( size_t  i = 0; i < n; i++ )
    {
        c   += tmp_costs[i];
        cmax = std::max( cmax, tmp_costs[i] );
    }// for

    c /= double(p);

    cmin = std::max( c, cmax );
    cmax = std::max( 2*c, cmax );

    while ( depth < _search_depth )
    {
        c = (cmin + cmax) / 2.0;

        proc = ffd( c, p, tmp_costs, sched );
        
        HDEBUG( to_string( "(TMFitSched) schedule : found scheduling for %d processors and c_max = %f",
                           proc, c ) );
        
        if ( proc <= p )
        {
            // remember scheduling for this valid distribution
            if ( proc == p )
            {
                found_valid = true;

                for ( uint i = 0; i < n; i++ )
                    best_sched[ perm[i] ] = sched[i];
            }// if
            
            cmax = c;
        }// if
        else
            cmin = c;

        depth++;
    }// while

    // reset last valid scheduling
    if ( found_valid )
    {
        for ( size_t  i = 0; i < n; i++ )
            sched[i] = best_sched[i];
    }// if

    return found_valid;
}

//
// bin-packing "first fit decreasing" (FFD) algorithm
//
uint
TMFitSched::ffd ( double                         cmax,
                  uint                           p,
                  const std::vector< double > &  costs,
                  std::vector< int > &           sched ) const
{
    std::vector< double >  bins( p );
    size_t                 n = costs.size();
    uint                   last_bin = 0;

    //
    // we asume that the costs are ordered
    //
    // choose next free item and put it in first
    // bin with enough space to hold it
    // capacity of bin is defined by cmax
    //

    for ( uint i = 0; i < p; i++ )
        bins[i] = 0.0;
    
    for ( size_t  i = 0; i < n; i++ )
    {
        double  c     = costs[i];
        uint    idx   = 0;
        bool    found = false;

        //
        // look for first bin able to hold c
        //

        for ( uint j = 0; j < last_bin; j++ )
        {
            if ( bins[j] + c <= cmax )
            {
                idx   = j;
                found = true;
                break;
            }// if
        }// for

        //
        // if a bin was found, put c into it
        // else create new bin
        //

        if ( found )
        {
            bins[idx] += c;
            sched[i]   = idx;
        }// if
        else
        {
            // more than p bins are not allowed
            if ( last_bin >= p )
                return p+1;
            
            bins[ last_bin ] = c;
            sched[i]         = last_bin;
            
            last_bin++;
        }// else
    }// for

    // return number of bins used for packing
    return last_bin;
}

}// namespace Hpro
