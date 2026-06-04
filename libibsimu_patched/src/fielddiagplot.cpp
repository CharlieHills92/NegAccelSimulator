/*! \file fielddiagplot.cpp
 *  \brief %Field diagnostic plotter.
 */

/* Copyright (c) 2005-2014 Taneli Kalvas. All rights reserved.
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

#include <fstream>
#include "ibsimu.hpp"
#include "fielddiagplot.hpp"


FieldDiagPlot::FieldDiagPlot( Frame &frame, const Geometry &geom )
    : _frame(&frame), _geom(&geom), _epot(NULL), _scharge(NULL), 
      _trajdens(NULL), _efield(NULL), _bfield(NULL), _N(100)
{
    _diag[0] = FIELD_EPOT;
    _diag[1] = FIELD_NONE;
    _loc[0] = FIELDD_LOC_DIST;
    _loc[1] = FIELDD_LOC_NONE;

    _graph[0] = NULL;
    _graph[1] = NULL;

    _legend[0] = NULL;
    _legend[1] = NULL;
}


FieldDiagPlot::~FieldDiagPlot()
{
    if( _graph[0] )
	delete _graph[0];
    if( _graph[1] )
	delete _graph[1];

    if( _legend[0] )
	delete _legend[0];
    if( _legend[1] )
	delete _legend[1];
}


void FieldDiagPlot::build_data( std::vector<double> coord[4],
				std::vector<double> fielddata[2] ) const
{
    // Build coordinate points
    coord[0].reserve( _N );
    coord[1].reserve( _N );
    coord[2].reserve( _N );
    coord[3].reserve( _N );
    for( size_t a = 0; a < 3; a++ ) {
	for( size_t b = 0; b < _N; b++ ) {
	    coord[a].push_back( _x1[a] + (b/(_N-1.0))*(_x2[a]-_x1[a]) );
	}
    }

    // Build distance
    for( size_t b = 0; b < _N; b++ ) {
	Vec3D d = (b/(_N-1.0))*(_x2-_x1);
	coord[3].push_back( d.norm2() );
    }

    // Build field data 
    for( size_t a = 0; a < 2; a++ ) {
	fielddata[a].reserve( _N );
	switch( _diag[a] ) {

	case FIELD_EPOT:
	    if( _epot == NULL )
		throw( Error( ERROR_LOCATION, "epot undefined" ) );		
	    for( size_t b = 0; b < _N; b++ ) {
		fielddata[a].push_back( (*_epot)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) ) );
	    }
	    break;

	case FIELD_EFIELD:
	    if( _efield == NULL )
		throw( Error( ERROR_LOCATION, "efield undefined" ) );		
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D E = (*_efield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		fielddata[a].push_back( E.norm2() );
	    }
	    break;

	case FIELD_EFIELD_X:
	    if( _efield == NULL )
		throw( Error( ERROR_LOCATION, "efield undefined" ) );		
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D E = (*_efield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		fielddata[a].push_back( E[0] );
	    }
	    break;

	case FIELD_EFIELD_Y:
	    if( _efield == NULL )
		throw( Error( ERROR_LOCATION, "efield undefined" ) );		
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D E = (*_efield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		fielddata[a].push_back( E[1] );
	    }
	    break;

	case FIELD_EFIELD_Z:
	    if( _efield == NULL )
		throw( Error( ERROR_LOCATION, "efield undefined" ) );		
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D E = (*_efield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		fielddata[a].push_back( E[2] );
	    }
	    break;

	case FIELD_SCHARGE:
	    if( _scharge == NULL )
		throw( Error( ERROR_LOCATION, "scharge undefined" ) );
	    for( size_t b = 0; b < _N; b++ ) {
		fielddata[a].push_back( (*_scharge)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) ) );
	    }
	    break;

	case FIELD_TRAJDENS:
	    if( _trajdens == NULL )
		throw( Error( ERROR_LOCATION, "trajectory density undefined" ) );
	    for( size_t b = 0; b < _N; b++ ) {
		fielddata[a].push_back( (*_trajdens)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) ) );
	    }
	    break;

	case FIELD_BFIELD:
	    if( _bfield == NULL )
		throw( Error( ERROR_LOCATION, "bfield undefined" ) );
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D B = (*_bfield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		fielddata[a].push_back( B.norm2() );
	    }
	    break;

	case FIELD_BFIELD_X:
	    if( _bfield == NULL )
		throw( Error( ERROR_LOCATION, "bfield undefined" ) );
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D B = (*_bfield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		fielddata[a].push_back( B[0] );
	    }
	    break;

	case FIELD_BFIELD_Y:
	    if( _bfield == NULL )
		throw( Error( ERROR_LOCATION, "bfield undefined" ) );
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D B = (*_bfield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		fielddata[a].push_back( B[1] );
	    }
	    break;

	case FIELD_BFIELD_Z:
	    if( _bfield == NULL )
		throw( Error( ERROR_LOCATION, "bfield undefined" ) );
	    for( size_t b = 0; b < _N; b++ ) {
		Vec3D B = (*_bfield)( Vec3D(coord[0][b],coord[1][b],coord[2][b]) );
		fielddata[a].push_back( B[2] );
	    }
	    break;

	case FIELD_NONE:
	    // Nothing to do
	    break;

	default:
	    throw( Error( ERROR_LOCATION, "undefined request" ) );
	    break;


	}
    }
}


std::string FieldDiagPlot::diagnostic_label( field_diag_type_e diag ) const
{
    switch( diag ) {
	
    case FIELD_EPOT:
	return( "\\phi  (V)" );
	break;
	
    case FIELD_EFIELD:
	return( "|E| (V/m)" );
	break;
	
    case FIELD_EFIELD_X:
	return( "E_x (V/m)" );
	break;
	
    case FIELD_EFIELD_Y:
	if( _geom->geom_mode() == MODE_CYL )
	    return( "E_r (V/m)" );
	else
	    return( "E_y (V/m)" );
	break;
	
    case FIELD_EFIELD_Z:
	return( "E_z (V/m)" );
	break;
	
    case FIELD_SCHARGE:
	return( "\\rho  (C/m^3)" );
	break;

    case FIELD_TRAJDENS:
	return( "J  (A/m^2)" );
	break;
	
    case FIELD_BFIELD:
	return( "B (T)" );
	break;
	
    case FIELD_BFIELD_X:
	return( "B_x (T)" );
	break;
	
    case FIELD_BFIELD_Y:
	if( _geom->geom_mode() == MODE_CYL )
	    return( "B_r (T)" );
	else
	    return( "B_y (T)" );
	break;
	
    case FIELD_BFIELD_Z:
	if( _geom->geom_mode() == MODE_CYL )
	    return( "B_\\theta  (T)" );
	else
	    return( "B_z (T)" );
	break;

    case FIELD_NONE:
    default:
	return( "" );
	break;
    }

    return( "" );
}


void FieldDiagPlot::build_plot( void )
{
    // Remove old graphs
    _frame->clear_graphs();

    if( _graph[0] )
	delete _graph[0];
    if( _graph[1] )
	delete _graph[1];

    if( _legend[0] )
	delete _legend[0];
    if( _legend[1] )
	delete _legend[1];

    // Build data
    std::vector<double> coord[4];
    std::vector<double> fielddata[2];
    build_data( coord, fielddata );

    // Set x-ranges manually
    _frame->ruler_autorange_enable( PLOT_AXIS_X1, false, false );
    _frame->ruler_autorange_enable( PLOT_AXIS_X2, false, false );

    //_frame->ruler_autorange_enable( PLOT_AXIS_Y1, false, false );
    //_frame->ruler_autorange_enable( PLOT_AXIS_Y2, false, false );

    // Enable y-axes according to use
    if( _diag[0] != FIELD_NONE )
	_frame->force_enable_ruler( PLOT_AXIS_Y1, true );
    if( _diag[1] != FIELD_NONE )
	_frame->force_enable_ruler( PLOT_AXIS_Y2, true );

    // X data
    PlotAxis xaxis_use = PLOT_AXIS_X1;
    std::vector<double> *datax = NULL;
    for( size_t b = 0; b < 2; b++ ) {
	
	// Check if xaxis disabled
	if( _loc[b] == FIELDD_LOC_NONE ) {
	    //std::cout << "Disabled\n";
	    continue;
	}
	
	// Select xaxis accordingly
	PlotAxis xaxis;
	if( b == 0 )
	    xaxis = PLOT_AXIS_X1;
	else
	    xaxis = PLOT_AXIS_X2;

	// Force enable ruler
	_frame->force_enable_ruler( xaxis, true );

	// Get datax pointer and set xaxis properties
	std::vector<double> *dataxt = NULL;
	if( _loc[b] == FIELDD_LOC_DIST ) { // Distance
	    _frame->set_axis_label( xaxis, "Distance (m)" );
	    _frame->set_ranges( xaxis, 0.0, coord[3][_N-1] );
	    dataxt = &coord[3];
	    //std::cout << "Distance: " << 0.0 << " to " << coord[3][_N-1] << "\n";
	} else if( _loc[b] == FIELDD_LOC_X ) { // X
	    _frame->set_axis_label( xaxis, "x (m)" );
	    _frame->set_ranges( xaxis, _x1[0], _x2[0] );
	    dataxt = &coord[0];
	    //std::cout << "X: " << _x1[0] << " to " << _x2[0] << "\n";
	} else if( _loc[b] == FIELDD_LOC_Y ) { // Y
	    if( _geom->geom_mode() == MODE_CYL )
		_frame->set_axis_label( xaxis, "r (m)" );
	    else
		_frame->set_axis_label( xaxis, "y (m)" );
	    _frame->set_ranges( xaxis, _x1[1], _x2[1] );
	    dataxt = &coord[1];
	    //std::cout << "Y: " << _x1[1] << " to " << _x2[1] << "\n";
	} else if( _loc[b] == FIELDD_LOC_Z ) { // Z
	    _frame->set_axis_label( xaxis, "z (m)" );
	    _frame->set_ranges( xaxis, _x1[2], _x2[2] );
	    dataxt = &coord[2];
	    //std::cout << "Z: " << _x1[2] << " to " << _x2[2] << "\n";
	}

	// Use dataxt if values not constant
	if( dataxt != NULL && (*dataxt)[0] != (*dataxt)[_N-1] ) {
	    datax = dataxt;
	    xaxis_use = xaxis;
	}
    }
    if( datax == NULL )
	throw( Error( ERROR_LOCATION, "both x-axes undefined or constant" ) );


    // Y data
    for( size_t a = 0; a < 2; a++ ) {

	// Check if there is data for yaxis
	if( _diag[a] == FIELD_NONE ) {
	    _graph[a] = NULL;
	    continue;
	}

	// Select yaxis accordingly
	PlotAxis yaxis;
	if( a == 0 )
	    yaxis = PLOT_AXIS_Y1;
	else
	    yaxis = PLOT_AXIS_Y2;

	// Set yaxis label
	_frame->set_axis_label( yaxis, diagnostic_label( _diag[a] ) );

	// Make graph
	_graph[a] = new XYGraph( *datax, fielddata[a] );
	_graph[a]->set_line_style( XYGRAPH_LINE_SOLID );
	_graph[a]->set_point_style( XYGRAPH_POINT_DISABLE, true, 0.0 );
	if( a == 0 )
	    _graph[a]->set_color( Vec3D(1,0,0) );
	else
	    _graph[a]->set_color( Vec3D(0,0,1) );

	_legend[a] = new LegendEntry( *_graph[a], diagnostic_label( _diag[a] ) );
	
	// Add graph
	_frame->add_graph( xaxis_use, yaxis, _graph[a], _legend[a] );
    }

    /*
    // Save default "zoom fit" ranges
    _frame->get_ranges( PLOT_AXIS_X1, _range[0], _range[2] );
    _frame->get_ranges( PLOT_AXIS_Y1, _range[1], _range[3] );
    _frame->get_ranges( PLOT_AXIS_X2, _range[4], _range[6] );
    _frame->get_ranges( PLOT_AXIS_Y2, _range[5], _range[7] );
    */

}


void FieldDiagPlot::export_data( const std::string &filename ) const
{
    std::ofstream fstr( filename.c_str() );
    ibsimu.message( 1 ) << "Exporting field diagnostic data to \'" << filename << "\'\n";

    // Build data
    std::vector<double> coord[4];
    std::vector<double> fielddata[2];
    build_data( coord, fielddata );

    // Write header
    fstr << "# ";
    fstr << std::setw(12) << "X (m) ";
    fstr << std::setw(14) << "Y (m) ";
    fstr << std::setw(14) << "Z (m) ";
    fstr << std::setw(14) << "Distance (m) ";
    if( _diag[0] != FIELD_NONE )
	fstr << std::setw(13) << diagnostic_label( _diag[0] ) << " ";
    if( _diag[1] != FIELD_NONE )
	fstr << std::setw(13) << diagnostic_label( _diag[1] ) << " ";
    fstr << "\n";

    // Write data
    for( size_t a = 0; a < _N; a++ ) {

	// Write coordinates
	for( size_t b = 0; b < 4; b++ )
	    fstr << std::setw(13) << coord[b][a] << " ";

	// Write data
	for( size_t b = 0; b < 2; b++ ) {
	    if( _diag[b] != FIELD_NONE )
		fstr << std::setw(13) << fielddata[b][a] << " ";
	}

	fstr << "\n";
    }

    fstr.close();
}






