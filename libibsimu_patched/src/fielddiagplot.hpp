/*! \file fielddiagplot.hpp
 *  \brief %Field diagnostic plotter.
 */

/* Copyright (c) 2005-2011,2014 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * taneli.kalvas@jyu.fi.
 * 
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
 */

#ifndef FIELDDIAGPLOT_HPP
#define FIELDDIAGPLOT_HPP 1


#include "types.hpp"
#include "frame.hpp"
#include "xygraph.hpp"
#include "vec3d.hpp"
#include "geometry.hpp"
#include "scalarfield.hpp"
#include "vectorfield.hpp"


#define FIELDD_DIAG_NONE FIELD_NONE
#define FIELDD_DIAG_EPOT FIELD_EPOT
#define FIELDD_DIAG_SCHARGE FIELD_SCHARGE
#define FIELDD_DIAG_TRAJDENS FIELD_TRAJDENS
#define FIELDD_DIAG_EFIELD FIELD_EFIELD
#define FIELDD_DIAG_EFIELD_X FIELD_EFIELD_X 
#define FIELDD_DIAG_EFIELD_Y FIELD_EFIELD_Y
#define FIELDD_DIAG_EFIELD_Z FIELD_EFIELD_Z
#define FIELDD_DIAG_BFIELD FIELD_BFIELD
#define FIELDD_DIAG_BFIELD_X FIELD_BFIELD_X
#define FIELDD_DIAG_BFIELD_Y FIELD_BFIELD_Y
#define FIELDD_DIAG_BFIELD_Z FIELD_BFIELD_Z


/*! \brief Location type for field diagnostics.
 *
 *  Indicator for coordinate
 */
enum field_loc_type_e {
    FIELDD_LOC_NONE = 0, /*!< \brief Dummy location */
    FIELDD_LOC_X,        /*!< \brief Location x-coordinate */
    FIELDD_LOC_Y,        /*!< \brief Location y-coordinate */
    FIELDD_LOC_Z,        /*!< \brief Location z-coordinate */
    FIELDD_LOC_DIST      /*!< \brief Location as distance */
};


/*! \brief Field diagnostics plot.
 *
 *  A class for building xy-plots of various fields in simulations.
 */
class FieldDiagPlot {

    Frame              *_frame;

    const Geometry     *_geom;
    const ScalarField  *_epot;
    const ScalarField  *_scharge;
    const ScalarField  *_trajdens;
    const VectorField  *_efield;
    const VectorField  *_bfield;

    size_t              _N;
    Vec3D               _x1;
    Vec3D               _x2;

    field_diag_type_e   _diag[2];
    field_loc_type_e    _loc[2];

    XYGraph            *_graph[2];
    LegendEntry        *_legend[2];

    void build_data( std::vector<double> coord[4], 
		     std::vector<double> fielddata[2] ) const;
    std::string diagnostic_label( field_diag_type_e diag ) const;

public:

    /*! \brief Constructor for field diagnostics plot.
     */
    FieldDiagPlot( Frame &frame, const Geometry &geom );

    /*! \brief Destructor for field diagnostics plot.
     */
    ~FieldDiagPlot();

    /*! \brief Add pointer to electric potential.
     */
    void set_epot( const ScalarField *epot ) {
	_epot = epot;
    }

    /*! \brief Add pointer to electric field.
     */
    void set_efield( const VectorField *efield ) {
	_efield = efield;
    }

    /*! \brief Add pointer to space charge density map.
     */
    void set_scharge( const ScalarField *scharge ) {
	_scharge = scharge;
    }

    /*! \brief Add pointer to trajectory density map.
     */
    void set_trajdens( const ScalarField *trajdens ) {
	_trajdens = trajdens;
    }

    /*! \brief Add pointer to magnetic field.
     */
    void set_bfield( const VectorField *bfield ) {
	_bfield = bfield;
    }

    /*! \brief Set coordinates for field diagnostics.
     *
     *  The fields to be plotted are evaluated at \a N steps from \a
     *  x1 to \a x2. The first point of the plot is exactly at \a x1
     *  and the last at \a x2.
     */
    void set_coordinates( size_t N, const Vec3D &x1, const Vec3D &x2 ) {
	_N = N;
	_x1 = x1;
	_x2 = x2;
    }

    /*! \brief Get start coordinates of diagnostic line.
     */
    const Vec3D &start( void ) {
	return( _x1 );
    }

    /*! \brief Get end coordinates of diagnostic line.
     */
    const Vec3D &end( void ) {
	return( _x2 );
    }

    /*! \brief Get number of steps on diagnostic line.
     */
    const size_t &N( void ) {
	return( _N );
    }

    /*! \brief Set field and location plot types
     *
     *  The plot can have two x-axes and two y-axes. The diagnostic
     *  for y1-axis is set by \a diag[0] and the diagnostic for
     *  y2-axis by \a diag[1]. The x1-axis is defined by \a loc[0] and
     *  x2-axis by \a loc[1].
     */
    void set_diagnostic( const field_diag_type_e diag[2], const field_loc_type_e loc[2] ) {
	_diag[0] = diag[0];
	_diag[1] = diag[1];
	_loc[0] = loc[0];
	_loc[1] = loc[1];
    }

    /*! \brief Get diagnostic type for y-axis \a i.
     */
    const field_diag_type_e &get_diagnostic_type( int i ) {
	return( _diag[i] );
    }

    /*! \brief Get location type for x-axis \a i.
     */
    const field_loc_type_e &get_location_type( int i ) {
	return( _loc[i] );
    }

    /*! \brief Export plotted data as ASCII.
     */
    void export_data( const std::string &filename ) const;

    /*! \brief Rebuild plot.
     */
    void build_plot( void );
};



#endif






