//
// Project     : HLIBpro
// File        : TCoordVis.cc
// Description : coordinate visualisation classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <fstream>
#include <string>
#include <memory>
#include <random>

#include "baseio.hh"
#include "TPSPrinter.hh"
#include "TRNG.hh"
#include "TColourMap.hh"

#include "hpro/io/TCoordVis.hh"

namespace Hpro
{

using std::unique_ptr;
using std::string;
using std::ostream;

namespace
{

//
// colours for labeled colouring
//
const colour_t  LBL_COLOURS[] = { Tango::Aluminium4,
                                  Tango::ScarletRed1,
                                  Tango::Chameleon1,
                                  Tango::Butter1,
                                  Tango::SkyBlue1,
                                  Tango::Chocolate1,
                                  Tango::Orange1,
                                  Tango::Plum1 };

const size_t    LBL_NCOLOURS = sizeof(LBL_COLOURS) / sizeof(LBL_COLOURS[0]);

}// namespace anonymous

///////////////////////////////////////////////////////////////
//
// PostScript based coord visualisation
//
///////////////////////////////////////////////////////////////

//////////////////////////////////////
//
// constructor and destructor
//

TPSCoordVis::TPSCoordVis ()
{
    _view = T3Point( 1, 0, 0 );
}

TPSCoordVis::TPSCoordVis ( const T3Point view_dir )
        : _view(view_dir)
{
    if ( _view.norm2() < 1e-16 )
        _view = T3Point( 1, 0, 0 );
    else
        _view.normalise2();
}

//////////////////////////////////////
//
// print coord
//

namespace
{

double
min_coord_dist ( const TCoordinate * coord )
{
    //
    // try to approximate minimal distance between coordinates
    //

    std::random_device               rd{};
    std::mt19937                     generator{ rd() };
    std::uniform_int_distribution<>  distr{ 0, int(coord->ncoord()-1) };

    size_t        n_samples = 0;
    const size_t  max_samples = std::min< size_t >( 1000, coord->ncoord() );
    double        min_dist  = -1;

    do
    {
        const idx_t   i0 = distr( generator );
        const TPoint  x0( 2, coord->coord(i0) );
        bool          success = false;

        for ( uint  i = 0; i < max_samples; ++i )
        {
            const idx_t  i1 = distr( generator );

            if ( i0 != i1 )
            {
                const TPoint  x1( 2, coord->coord(i1) );
                const TPoint  d( x0 - x1 );
                const double  dist = d.norm2();
                    
                if ( std::abs( dist ) > 1e-10 )
                {
                    if (( min_dist < 0 ) || ( dist < min_dist ))
                        min_dist = dist;

                    success = true;
                }// if
            }// if
        }// for

        if ( success )
            n_samples++;
            
    } while ( n_samples < max_samples );

    return min_dist;
}

void
coord_print_ps ( const TCoordinate *          coord,
                 const std::vector< uint > *  label,
                 const std::vector< bool > *  filter,
                 const string &               filename )
{
    if ( coord->dim() == 2 )
    {
        //
        // determine bounding box of coordinates
        //

        T2Point  bbmin( 0, 0 );
        T2Point  bbmax( 1, 1 );

        if ( coord->ncoord() > 0 )
        {
            bbmin = T2Point( coord->coord(0)[0], coord->coord(0)[1] );
            bbmax = bbmin;
            
            for ( uint i = 1; i < coord->ncoord(); i++ )
            {
                bbmin[0] = std::min( bbmin[0], coord->coord(i)[0] );
                bbmax[0] = std::max( bbmax[0], coord->coord(i)[0] );

                bbmin[1] = std::min( bbmin[1], coord->coord(i)[1] );
                bbmax[1] = std::max( bbmax[1], coord->coord(i)[1] );
            }// for
        }// if
        
        //
        // determine magnification factor such that output
        // roughly fits in 500x500 window
        //
    
        const double  len_x = bbmax[0] - bbmin[0];
        const double  len_y = bbmax[1] - bbmin[1];
        double        mag   = 1.0;

        if ( std::max( len_x, len_y ) * mag > 500 )
        {
            while ( std::max( len_x, len_y ) * mag > 500 )
                mag /= 2.0;
        }// if
        else
        {
            while ( std::max( len_x, len_y ) * mag < 500 )
                mag *= 2.0;
            mag /= 2.0;
        }// else
    
        //
        // set up printer
        //

        TPSPrinter  prn( uint(mag*len_x), uint(mag*len_y), add_extension( filename, "eps" ) );
        
        prn.begin();
        prn.scale( 1.0, -1.0 );
        prn.scale( mag, mag );
        prn.translate( -bbmin[0], -bbmax[1] );

        //
        // print coordinates
        //

        const auto  min_dist = min_coord_dist( coord );
        
        prn.set_gray( 0 );
        
        for ( idx_t  i = 0; i < idx_t(coord->ncoord()); i++ )
        {
            // filter coordinates
            if (( filter != nullptr ) && ( idx_t(filter->size()) > i ))
            {
                if ( ! (*filter)[i] )
                    continue;
            }// if

            // colourise based on label
            if (( label != nullptr ) && ( idx_t(label->size()) > i ))
            {
                if ( (*label)[i] >= LBL_NCOLOURS )
                {
                    HERROR( ERR_ARG, "(TPSCoordVis) print",
                            "label larger than LBL_NCOLOURS = " + to_string( "%d", LBL_NCOLOURS ) );
                }// if
                
                const colour_t  col = LBL_COLOURS[ (*label)[i] ];
                
                prn.set_colour( col );
            }// if
            
            prn.fill_circle( coord->coord(i)[0], coord->coord(i)[1], min_dist / 5.0 );
        }// for

        prn.end();
    }// if
    else
        HERROR( ERR_NOT_IMPL, "(TPSCoordVis) print", "only 2d coordinates supported" );
}

}// namespace anonymous

//
// print <coord> to file <filename>
//
void
TPSCoordVis::print ( const TCoordinate * coord,
                     const string &      filename ) const
{
    coord_print_ps( coord, nullptr, nullptr, filename );
}

//
// print coordinates in a given cluster
//
void
TPSCoordVis::print ( const TCoordinate *          coord,
                     const std::vector< uint > &  label,
                     const string &               filename ) const
{
    coord_print_ps( coord, & label, nullptr, filename );
}

//
//! print only coordinates \a coord with filter(i)=true
//
void
TPSCoordVis::print ( const TCoordinate *          coord,
                     const std::vector< uint > &  label,
                     const std::vector< bool > &  filter,
                     const std::string &          filename ) const
{
    coord_print_ps( coord, & label, & filter, filename );
}

///////////////////////////////////////////////////////////////
//
// VTK based coord visualisation
//
///////////////////////////////////////////////////////////////

namespace
{

//
// write point set in VTK format
//
void
vtk_write_pointset ( ostream &           out,
                     const TCoordinate *  coord )
{
    out << "# vtk DataFile Version 2.0" << std::endl
        << "HLIBpro coordinates" << std::endl
        << "ASCII" << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;

    if ( coord != nullptr )
    {
        out << "POINTS " << coord->ncoord() << " FLOAT" << std::endl;
        
        for ( uint  i = 0; i < coord->ncoord(); ++i )
        {
            uint  j = 0;

            for ( ; j < coord->dim(); ++j ) out << coord->coord(i)[j] << " ";
            for ( ; j < 3;            ++j ) out << "0 ";
            out << std::endl;
        }// for

        out << "CELLS " << coord->ncoord() << " " << 2 * coord->ncoord() << std::endl;

        for ( uint  i = 0; i < coord->ncoord(); ++i )
        {
            out << "1 " << i << " ";
        }// for
        out << std::endl;
        
        out << "CELL_TYPES " << coord->ncoord() << std::endl;
        
        for ( uint  i = 0; i < coord->ncoord(); ++i )
        {
            out << "1 ";
        }// for
        out << std::endl;
    }// if
}

}// namespace anonymous

//
// print <coord> to file <filename>
//
void
TVTKCoordVis::print ( const TCoordinate * coord,
                      const string &      filename ) const
{
    unique_ptr< std::ostream >  out_ptr( open_write( add_extension( filename, "vtk" ) ) );
    std::ostream &              out = * out_ptr.get();

    if ( coord != nullptr )
    {
        vtk_write_pointset( out, coord );
    }// if
}

//
// print coordinates \a coord and colour coordinates with same
// label with same colour; the labels are defined by \label and
// must begin with 0
//
void
TVTKCoordVis::print ( const TCoordinate *          coord,
                      const std::vector< uint > &  label,
                      const string &               filename ) const
{
    unique_ptr< std::ostream >  out_ptr( open_write( add_extension( filename, "vtk" ) ) );
    std::ostream &              out = * out_ptr.get();

    if ( coord != nullptr )
    {
        vtk_write_pointset( out, coord );

        if ( label.size() < coord->ncoord() )
            return;
        
        out << "CELL_DATA " << coord->ncoord() << std::endl
            << "COLOR_SCALARS vtxcolour 3" << std::endl;
        
        for ( uint  i = 0; i < coord->ncoord(); ++i )
        {
            if ( label[i] >= LBL_NCOLOURS )
            {
                HERROR( ERR_ARG, "(TVTKCoordVis) print",
                        "label larger than LBL_NCOLOURS = " + to_string( "%d", LBL_NCOLOURS ) );
            }// if
            
            const colour_t  col = LBL_COLOURS[ label[i] ];

            out << double( col.red   ) / 255.0 << " "
                << double( col.green ) / 255.0 << " "
                << double( col.blue  ) / 255.0 << " "
                << std::endl;
        }// for
    }// if
}

// //
// // print coordinates \a coord and mark those in cluster \a bc
// // - conversion of indices from coordinate numbering to
// //   cluster numbering is performed by \a perm_i2e 
// void
// TVTKCoordVis::print ( const TCoordinate *   coord,
//                       const TBlockCluster * bc,
//                       const TPermutation *  perm_e2i,
//                       const string &        filename ) const
// {
//     ostream  out( filename.c_str() );

//     if ( coord != nullptr )
//     {
//         vtk_write_pointset( out, coord );

//         if ( bc == nullptr )
//             return;
        
//         const TCluster *  rowcl = bc->rowcl();
//         const TCluster *  colcl = bc->colcl();
        
//         if (( rowcl == nullptr ) || ( colcl == nullptr ))
//             return;
        
//         out << "CELL_DATA " << coord->ncoord() << std::endl
//             << "COLOR_SCALARS " << coord->ncoord() << " 3" << std::endl;
        
//         for ( uint  i = 0; i < coord->ncoord(); ++i )
//         {
//             const idx_t  cidx = perm_e2i->permute( i );

//             if ( rowcl->is_in( cidx ) )
//             {
//                 out << "1 0 0 ";
//             }// if
//             else if ( colcl->is_in( cidx ) )
//             {
//                 out << "0 1 0 ";
//             }// if
//             else
//             {
//                 out << "0.5 0.5 0.5 ";
//             }// if
//         }// for
//     }// if
// }

//
// print coordinates \a coord and the connections between
// them as defined by the coefficients of the sparse matrix \a S
// - the ordering of \a S is assumed to be equal to \a coord
//
template < typename value_t >
void
TVTKCoordVis::print ( const TCoordinate *               coord,
                      const TSparseMatrix< value_t > *  S,
                      const string &                    filename,
                      const coltype_t                   col_type ) const
{
    using  real_t = real_type_t< value_t >;
    
    unique_ptr< std::ostream >  out_ptr( open_write( add_extension( filename, "vtk" ) ) );
    std::ostream &              out = * out_ptr.get();

    if ( coord == nullptr )
        return;

    out << "# vtk DataFile Version 2.0" << std::endl
        << "HLIBpro coordinates" << std::endl
        << "ASCII" << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;

    out << "POINTS " << coord->ncoord() << " FLOAT" << std::endl;
        
    for ( uint  i = 0; i < coord->ncoord(); ++i )
    {
        uint  j = 0;
            
        for ( ; j < coord->dim(); ++j )
            out << coord->coord(i)[j] << " ";
        for ( ; j < 3; ++j )
            out << "0 ";
        out << std::endl;
    }// for

    //
    // determine number of "real" edges to print
    //

    const bool  is_nonsym = S->is_nonsym();
    size_t      nedges = 0;

    for ( idx_t  row = 0; row < idx_t( S->rows() ); ++row )
    {
        const auto  lb = S->rowptr( row );
        const auto  ub = S->rowptr( row+1 );

        for ( auto  j = lb; j < ub; ++j )
        {
            const auto  neigh = S->colind( j );

            //
            // only draw one edge between two indices, hence,
            // only count lower half (but check upper in case
            // of unsymmetry); always exclude diagonal
            //

            if (( neigh < row ) || (( neigh > row ) && is_nonsym && ! S->has_entry( neigh, row )))
            {
                nedges++;
            }// if
        }// for
    }// for

#if 0

    //
    // cells: only the edges
    //
    
    out << "CELLS " << nedges << " " << 3 * nedges << std::endl;

    for ( size_t  row = 0; row < S->rows(); ++row )
    {
        const idx_t  lb = S->rowptr( row );
        const idx_t  ub = S->rowptr( row+1 );

        for ( idx_t  j = lb; j < ub; ++j )
        {
            const idx_t  neigh = S->colind( j );

            if (( neigh < idx_t(row) ) || (( neigh > idx_t(row) ) && is_nonsym && ! S->has_entry( neigh, idx_t(row) )))
            {
                out << "2 " << row << " " << neigh << " ";
            }// if
        }// for
    }// for
    out << std::endl;

    //
    // cell types
    //
    
    out << "CELL_TYPES " << nedges << std::endl;
        
    for ( uint  i = 0; i < nedges; ++i )
    {
        out << "3 ";
    }// for
    out << std::endl;

    // stop, if no colour is wanted
    if ( col_type == TVTKCoordVis::NO_COLOUR )
        return;
        
    //
    // finally, define color of edges
    //

    out << "CELL_DATA " << nedges << std::endl
        << "COLOR_SCALARS cellcolour" << " 3" << std::endl;
        
    // for edge color, determine minimal and maximal value and 
    // use colourmap

    real_t  min_val = Limits::max< real_t >();
    real_t  max_val = 0.0;
        
    for ( size_t  row = 0; row < S->rows(); ++row )
    {
        const idx_t  lb = S->rowptr( row );
        const idx_t  ub = S->rowptr( row+1 );
        
        for ( idx_t  j = lb; j < ub; ++j )
        {
            const idx_t  neigh = S->colind( j );
            
            if ( idx_t(row) == neigh )
                continue;
            
            min_val = min( min_val, Math::abs( S->coeff( j ) ) );
            max_val = max( max_val, Math::abs( S->coeff( j ) ) );
        }// for
    }// for

    if ( col_type == TVTKCoordVis::LOG_COLOUR )
    {
        min_val = Math::log( min_val );
        max_val = Math::log( max_val );
    }// if
    
    THotColourMap  cmap( 1000 );
    const real_t   val_diff = max_val - min_val;

    if ( val_diff == 0.0 )
    {
        //
        // only draw gray edges if all coefficients have same value
        //
        
        for ( uint  i = 0; i < nedges; ++i )
        {
            out << "0.5 0.5 0.5" << " ";
        }// for
        out << std::endl;
    }// if
    else
    {
        //
        // colour edges using colourmap
        //
        
        for ( size_t  row = 0; row < S->rows(); ++row )
        {
            const idx_t  lb = S->rowptr( row );
            const idx_t  ub = S->rowptr( row+1 );

            for ( idx_t  j = lb; j < ub; ++j )
            {
                const idx_t  neigh = S->colind( j );

                // remark: the if below is to ensure same egde order as before
                if (( neigh < idx_t(row) ) || (( neigh > idx_t(row) ) && is_nonsym && ! S->has_entry( neigh, idx_t(row) )))
                {
                    // use maximal value of edges if two are present
                    auto  val = Math::abs( S->coeff( j ) );
                        
                    if ( S->has_entry( neigh, row ) )
                        val = max( val, Math::abs( S->entry( neigh, idx_t(row) ) ) );

                    if ( col_type == TVTKCoordVis::LOG_COLOUR )
                        val = Math::log( val );
                    
                    const colour_t  col( cmap.fentry( (val - min_val) / val_diff ) );

                    out << double( col.red   ) / 255.0 << " "
                        << double( col.green ) / 255.0 << " "
                        << double( col.blue  ) / 255.0 << " " << std::endl;
                }// if
            }// for
        }// for
        out << std::endl;
    }// else

#else
    
    //
    // cells: first the coordinates, then the edges
    //
    
    out << "CELLS " << coord->ncoord() + nedges << " " << 2 * coord->ncoord() + 3 * nedges << std::endl;

    for ( uint  i = 0; i < coord->ncoord(); ++i )
    {
        out << "1 " << i << " ";
    }// for
    out << std::endl;
        
    for ( idx_t  row = 0; row < idx_t(S->rows()); ++row )
    {
        const auto  lb = S->rowptr( row );
        const auto  ub = S->rowptr( row+1 );

        for ( auto  j = lb; j < ub; ++j )
        {
            const auto  neigh = S->colind( j );

            if (( neigh < row ) || (( neigh > row ) && is_nonsym && ! S->has_entry( neigh, row )))
            {
                out << "2 " << row << " " << neigh << " ";
            }// if
        }// for
    }// for
    out << std::endl;

    //
    // same for cell types
    //
    
    out << "CELL_TYPES " << coord->ncoord() + nedges << std::endl;
        
    for ( uint  i = 0; i < coord->ncoord(); ++i )
    {
        out << "1 " << " ";
    }// for
    out << std::endl;

    for ( uint  i = 0; i < nedges; ++i )
    {
        out << "3 " << " ";
    }// for
    out << std::endl;

    //
    // finally, define color of nodes and edges
    //

    out << "CELL_DATA " << coord->ncoord() + nedges << std::endl
        << "COLOR_SCALARS cellcolour" << " 3" << std::endl;
        
    for ( uint  i = 0; i < coord->ncoord(); ++i )
    {
        out << "0.5 0.5 0.5 ";
    }// for


    // for edge color, determine minimal and maximal value and 
    // use colourmap

    real_t  min_val = Limits::max< real_t >();
    real_t  max_val = 0.0;
        
    for ( idx_t  i = 0; i < idx_t(S->n_non_zero()); ++i )
    {
        min_val = std::min( min_val, Math::abs( S->coeff( i ) ) );
        max_val = std::max( max_val, Math::abs( S->coeff( i ) ) );
    }// if

    if ( col_type == TVTKCoordVis::LOG_COLOUR )
    {
        min_val = Math::log( min_val );
        max_val = Math::log( max_val );
    }// if
    
    THotColourMap  cmap( 1000 );
    const real_t   val_diff = max_val - min_val;

    if ( val_diff == 0.0 )
    {
        //
        // only draw gray edges if all coefficients have same value
        //
        
        for ( uint  i = 0; i < nedges; ++i )
        {
            out << "0.5 0.5 0.5" << " ";
        }// for
        out << std::endl;
    }// if
    else
    {
        //
        // colour edges using colourmap
        //
        
        for ( idx_t  row = 0; row < idx_t(S->rows()); ++row )
        {
            const auto  lb = S->rowptr( row );
            const auto  ub = S->rowptr( row+1 );

            for ( auto  j = lb; j < ub; ++j )
            {
                const auto  neigh = S->colind( j );
                auto        val   = Math::abs( S->coeff( j ) );

                // remark: the if below is to ensure same egde order as before
                if (( neigh < row ) || (( neigh > row ) && is_nonsym && ! S->has_entry( neigh, row )))
                {
                    // use maximal value of edges if two are present
                    if ( S->has_entry( neigh, row ) )
                        val = std::max( val, Math::abs( S->entry( neigh, row ) ) );

                    if ( col_type == TVTKCoordVis::LOG_COLOUR )
                        val = Math::log( val );

                    const colour_t  col( cmap.fentry( (val - min_val) / val_diff ) );

                    out << double( col.red   ) / 255.0 << " "
                        << double( col.green ) / 255.0 << " "
                        << double( col.blue  ) / 255.0 << " ";
                }// if
            }// for
        }// for
        out << std::endl;
    }// else

#endif
}

#define INST_VTKPRINT( type ) \
    template \
    void \
    TVTKCoordVis::print< type > ( const TCoordinate *           , \
                                  const TSparseMatrix< type > * , \
                                  const string &                , \
                                  const coltype_t                ) const;

INST_VTKPRINT( float )
INST_VTKPRINT( double )
INST_VTKPRINT( std::complex< float > )
INST_VTKPRINT( std::complex< double > )

//
// functional versions
//
void
print_ps ( const TCoordinate *  coord,
           const std::string &  filename )
{
    TPSCoordVis  vis;

    vis.print( coord, filename );
}

void
print_vtk ( const TCoordinate *  coord,
            const std::string &  filename )
{
    TVTKCoordVis  vis;

    vis.print( coord, filename );
}

void
print_vtk ( const TCoordinate *          coord,
            const std::vector< uint > &  label,
            const std::string &          filename )
{
    TVTKCoordVis  vis;

    vis.print( coord, label, filename );
}

}// namespace
