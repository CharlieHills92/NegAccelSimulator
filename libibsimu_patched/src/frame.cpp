/*! \file frame.cpp
 *  \brief %Frame for plots
 */

/* Copyright (c) 2005-2013 Taneli Kalvas. All rights reserved.
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

#include <iostream>
#include <cmath>
#include <limits>
#include "compmath.hpp"
#include "frame.hpp"
#include "ibsimu.hpp"


//#define DEBUG_FRAME 1


#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MAX3(X,Y,Z) ((X) > (Y) ? ((X) > (Z) ? (X) : (Z)) : ((Y) > (Z) ? (Y) : (Z)))

#define MIN(X,Y) ((X) > (Y) ? (Y) : (X))
#define MIN3(X,Y,Z) ((X) > (Y) ? ((Y) > (Z) ? (Z) : (Y)) : ((X) > (Z) ? (Z) : (X)))


Frame::Frame()
    : _offx(0), _offy(0), _width(640), _height(480), 
      _fontsize(12.0), _titlespace(10.0), _cmlspace(10.0), _bg(Vec3D(1,1,1)), _fg(Vec3D(0,0,0)),
      _legend_enable(true), _legend_pos(LEGEND_POS_TOP_RIGHT), _cm_legend(NULL),
      _cml_enable(true), _fixedaspect(PLOT_FIXED_ASPECT_DISABLED), _automargin(true)
{
    _legend.set_font_size( _fontsize );

    // Ruler X1
    _ruler[0].set_coord_index( 0 );
    _ruler[0].set_indir( true );
    _fenable[0] = false;
    _autorange[0] = true;
    _autorange[1] = true;

    // Ruler Y1
    _ruler[1].set_coord_index( 1 );
    _ruler[1].set_indir( false );
    _fenable[1] = false;
    _autorange[2] = true;
    _autorange[3] = true;

    // Ruler X2
    _ruler[2].set_coord_index( 0 );
    _ruler[2].set_indir( false );
    _fenable[2] = false;
    _autorange[4] = true;
    _autorange[5] = true;

    // Ruler Y2
    _ruler[3].set_coord_index( 1 );
    _ruler[3].set_indir( true );
    _fenable[3] = false;
    _autorange[6] = true;
    _autorange[7] = true;
}


Frame::Frame( const Frame &frame )
    : _offx(frame._offx), _offy(frame._offy), _width(frame._width), _height(frame._height),
      _fontsize(frame._fontsize), _titlespace(frame._titlespace), _cmlspace(frame._cmlspace),
      _bg(frame._bg), _fg(frame._fg), _dobj(frame._dobj), _legend_enable(frame._legend_enable),
      _legend_pos(frame._legend_pos), _legend(frame._legend), _cml_enable(frame._cml_enable),
      _title(frame._title), _fixedaspect(frame._fixedaspect), _automargin(frame._automargin)      
{
    for( size_t a = 0; a < 4; a++ ) {
	_ruler[a] = frame._ruler[a];
	_cm[a] = frame._cm[a];
	_enable[a] = frame._enable[a];
	_fenable[a] = frame._fenable[a];
	_autorange[2*a+0] = frame._autorange[2*a+0];
	_autorange[2*a+1] = frame._autorange[2*a+1];
	_range_min[a] = frame._range_min[a];
	_range_max[a] = frame._range_max[a];
	_tmargin[a] = frame._tmargin[a];
    }

    if( frame._cm_legend ) {
	_cm_legend = new ColormapLegend( *frame._cm_legend );
    } else {
	_cm_legend = NULL;
    }
}


Frame::~Frame()
{
    if( _cm_legend )
	delete _cm_legend;    
}


void Frame::set_frame_clipping( cairo_t *cairo )
{
    cairo_rectangle( cairo, _offx+_tmargin[0], 
		     _offy+_tmargin[3], 
		     _width-_tmargin[0]-_tmargin[2], 
		     _height-_tmargin[1]-_tmargin[3] );
    cairo_clip( cairo );
}


void Frame::unset_frame_clipping( cairo_t *cairo )
{
    cairo_reset_clip( cairo );
}


void Frame::set_axis_label( PlotAxis axis, const std::string &label )
{
    switch( axis ) {
    case PLOT_AXIS_X1:
	_ruler[0].set_axis_label( label );
	break;
    case PLOT_AXIS_Y1:
	_ruler[1].set_axis_label( label );
	break;
    case PLOT_AXIS_X2:
	_ruler[2].set_axis_label( label );
	break;
    case PLOT_AXIS_Y2:
	_ruler[3].set_axis_label( label );
	break;
    case PLOT_AXIS_Z:
	break;
    };
}

void Frame::set_font_size( double size ) 
{
    _ruler[0].set_font_size( size );
    _ruler[1].set_font_size( size );
    _ruler[2].set_font_size( size );
    _ruler[3].set_font_size( size );
    _title.set_font_size( size );
    _legend.set_font_size( size );
    _fontsize = size;
    _titlespace = 10.0/12.0*size;
    _cmlspace = 10.0/12.0*size;
}

void Frame::set_automargin( bool enable )
{
    _automargin = enable;
}


void Frame::set_fixed_aspect( PlotFixedMode mode )
{
    _fixedaspect = mode;
}


void Frame::add_graph( PlotAxis xaxis, PlotAxis yaxis, 
		       Graph *graph, LegendEntry *legend )
{
    _dobj.push_back( DObj( xaxis, yaxis, graph ) );
    _legend.add_entry( legend );
}


void Frame::clear_graphs( void )
{
    _dobj.clear();
    _legend.clear_entries();
    _cm_legend = NULL;
}


void Frame::enable_legend( bool enable )
{
    _legend_enable = enable;
}


void Frame::set_legend_position( legend_position_e pos )
{
    _legend_pos = pos;
}


void Frame::ruler_autorange_enable( PlotAxis axis, bool min, bool max )
{
    switch( axis ) {
    case PLOT_AXIS_X1:
	_ruler[0].set_autorange( min, max );
	break;
    case PLOT_AXIS_Y1:
	_ruler[1].set_autorange( min, max );
	break;
    case PLOT_AXIS_X2:
	_ruler[2].set_autorange( min, max );
	break;
    case PLOT_AXIS_Y2:
	_ruler[3].set_autorange( min, max );
	break;
    case PLOT_AXIS_Z:
	// Do nothing
	break;
    };
}


void Frame::set_title( const std::string &title )
{
    _title.set_text( title );
}


void Frame::get_margins( double margin[4] ) const
{
    margin[0] = _tmargin[0];
    margin[1] = _tmargin[1];
    margin[2] = _tmargin[2];
    margin[3] = _tmargin[3];
}


void Frame::get_frame_edges( double edge[4] ) const
{
    edge[0] = _offx + _tmargin[0];
    edge[1] = _offy + _tmargin[3];
    edge[2] = _offx + _width - _tmargin[2];
    edge[3] = _offy + _height - _tmargin[1];
}


Coordmapper Frame::get_coordmapper( PlotAxis xaxis, PlotAxis yaxis ) const
{
    return( Coordmapper( _cm[xaxis], _cm[yaxis] ) );
}


void Frame::set_ranges( PlotAxis axis, double min, double max )
{
    bool autorange_min, autorange_max;
    if( comp_isinf(min) )
	autorange_min = true;
    else
	autorange_min = false;
    if( comp_isinf(max) )
	autorange_max = true;
    else
	autorange_max = false;

    switch( axis ) {
    case PLOT_AXIS_X1:
	_autorange[0] = autorange_min;
	_autorange[1] = autorange_max;
	_ruler[0].set_autorange( autorange_min, autorange_max );
	_ruler[0].set_ranges( min, max );
	break;
    case PLOT_AXIS_Y1:
	_autorange[2] = autorange_min;
	_autorange[3] = autorange_max;
	_ruler[1].set_autorange( autorange_min, autorange_max );
	_ruler[1].set_ranges( min, max );
	break;
    case PLOT_AXIS_X2:
	_autorange[4] = autorange_min;
	_autorange[5] = autorange_max;
	_ruler[2].set_autorange( autorange_min, autorange_max );
	_ruler[2].set_ranges( min, max );
	break;
    case PLOT_AXIS_Y2:
	_autorange[6] = autorange_min;
	_autorange[7] = autorange_max;
	_ruler[3].set_autorange( autorange_min, autorange_max );
	_ruler[3].set_ranges( min, max );
	break;
    case PLOT_AXIS_Z:
	// Do nothing
	break;
    };
}


void Frame::get_ranges( PlotAxis axis, double &min, double &max ) const
{
    switch( axis ) {
    case PLOT_AXIS_X1:
	_ruler[0].get_ranges( min, max );
	break;
    case PLOT_AXIS_Y1:
	_ruler[1].get_ranges( min, max );
	break;
    case PLOT_AXIS_X2:
	_ruler[2].get_ranges( min, max );
	break;
    case PLOT_AXIS_Y2:
	_ruler[3].get_ranges( min, max );
	break;
    case PLOT_AXIS_Z:
	// Do nothing
	break;
    };    
}


void Frame::calculate_autoranging( void )
{
    double range[8]; // x1min, x1max, y1min, y1max, x2min, x2max, y2min, y2max

#ifdef DEBUG_FRAME
    std::cout << "\nFRAME::CALCULATE_AUTORANGING()\n\n";
#endif

    for( size_t a = 0; a < 4; a++ ) {
	range[2*a+0] = std::numeric_limits<double>::infinity();
	range[2*a+1] = -std::numeric_limits<double>::infinity();
    }

    for( size_t a = 0; a < _dobj.size(); a++ ) {

	// Get bounding box for drawable
        double bbox[4] = { std::numeric_limits<double>::infinity(),
			   std::numeric_limits<double>::infinity(),
			   -std::numeric_limits<double>::infinity(),
			   -std::numeric_limits<double>::infinity() };
        _dobj[a]._graph->get_bbox( bbox );

	// Sort to correct axes and get largest combination of bounding boxes
	if( _dobj[a]._xaxis == PLOT_AXIS_X1 ) {
	    if( bbox[0] < range[0] )
		range[0] = bbox[0];
	    if( bbox[2] > range[1] )
		range[1] = bbox[2];
	} else /* PLOT_AXIS_X2 */ {
	    if( bbox[0] < range[4] )
		range[4] = bbox[0];
	    if( bbox[2] > range[5] )
		range[5] = bbox[2];
	}
	if( _dobj[a]._yaxis == PLOT_AXIS_Y1 ) {
	    if( bbox[1] < range[2] )
		range[2] = bbox[1];
	    if( bbox[3] > range[3] )
		range[3] = bbox[3];
	} else /* PLOT_AXIS_Y2 */ {
	    if( bbox[1] < range[6] )
		range[6] = bbox[1];
	    if( bbox[3] > range[7] )
		range[7] = bbox[3];
	}
    }

    // Set ruler ranges for autoranged axes and/or adjust ranges to
    // prevent zero range span
    for( size_t a = 0; a < 4; a++ ) {

	//bool autorange_min, autorange_max;
	double range_min, range_max;
	//_ruler[a].get_autorange( autorange_min, autorange_max );

	// Get manually set/previous ranges
	_ruler[a].get_ranges( range_min, range_max );

	// Set automatic ranges
	if( _autorange[2*a+0] )
	    range_min = range[2*a+0];
	if( _autorange[2*a+1] )
	    range_max = range[2*a+1];

	// Correct if range span zero
	if( range_max - range_min == 0.0 ) {
	    range_max = 1.1*range_max + 1.0e-6;
	    range_min = 0.9*range_min - 1.0e-6;
	}

	// Set new ranges
	_range_min[a] = range_min;
	_range_max[a] = range_max;
	_ruler[a].set_ranges( range_min, range_max );

#ifdef DEBUG_FRAME
	std::cout << "ruler[" << a << "] = {" << range_min << ", " << range_max << "}\n";
#endif
    }
}


void Frame::calculate_rulers( cairo_t *cairo, bool ruler_tic_bbox_test )
{
#ifdef DEBUG_FRAME
    std::cout << "\nFRAME::CALCULATE_RULERS()\n\n";
#endif

    // Recalculate enabled rulers
    for( size_t a = 0; a < 4; a++ ) {
	if( _enable[a] )
	    _ruler[a].calculate( cairo, _cm[a], ruler_tic_bbox_test );
    }
}


void Frame::calculate_ruler_autoenable( void )
{
    // Preset to disabled if not force enabled
    for( size_t a = 0; a < 4; a++ ) {
	if( _fenable[a] )
	    _enable[a] = true;
	else
	    _enable[a] = false;
    }

    // Get used axes
    for( size_t a = 0; a < _dobj.size(); a++ ) {
	if( _dobj[a]._xaxis == PLOT_AXIS_X1 )
	    _enable[0] = true;
	if( _dobj[a]._yaxis == PLOT_AXIS_Y1 )
	    _enable[1] = true;
	if( _dobj[a]._xaxis == PLOT_AXIS_X2 )
	    _enable[2] = true;
	if( _dobj[a]._yaxis == PLOT_AXIS_Y2 )
	    _enable[3] = true;
    }
}


void Frame::force_enable_ruler( PlotAxis axis, bool force )
{
    switch( axis ) {
    case PLOT_AXIS_X1:
	_fenable[0] = force;
	break;
    case PLOT_AXIS_Y1:
	_fenable[1] = force;
	break;
    case PLOT_AXIS_X2:
	_fenable[2] = force;
	break;
    case PLOT_AXIS_Y2:
	_fenable[3] = force;
	break;
    case PLOT_AXIS_Z:
	break;
    };
}


void Frame::calculate_frame( cairo_t *cairo )
{
#ifdef DEBUG_FRAME
    std::cout << "\nFRAME::CALCULATE_FRAME()\n\n";
#endif

    // Set initial margins
    _tmargin[0] = _tmargin[1] = _tmargin[2] = _tmargin[3] = 0.0;

    // Calculate first guess plot rectangle
    double prec[4]; // xmin, ymin, xmax, ymax
    prec[0] = _offx + _tmargin[0];
    prec[1] = _offy + _tmargin[3];
    prec[2] = _offx + _width - _tmargin[2];
    prec[3] = _offy + _height - _tmargin[1];

    // Set ruler end points
    _ruler[0].set_endpoints( prec[0], prec[3], prec[2], prec[3] );
    _ruler[1].set_endpoints( prec[0], prec[3], prec[0], prec[1] );
    _ruler[2].set_endpoints( prec[0], prec[1], prec[2], prec[1] );
    _ruler[3].set_endpoints( prec[2], prec[3], prec[2], prec[1] );

    // Set ruler properties
    calculate_ruler_autoenable();
    calculate_rulers( cairo, false );
    mirror_rulers();

    // Get ruler bounding boxes
    double bbox[4][4];
    for( size_t a = 0; a < 4; a++ )
	_ruler[a].get_bbox( cairo, bbox[a], _cm[a], false );

    /*
    // Plot ruler bounding boxes
    for( size_t a = 0; a < 4; a++ ) {
	cairo_set_source_rgb( cairo, 1, 0, 0 );
	cairo_set_line_width( cairo, 1 );
	cairo_rectangle( cairo, bbox[a][0], bbox[a][1], 
			 bbox[a][2]-bbox[a][0], bbox[a][3]-bbox[a][1] );
	cairo_stroke( cairo );
    }
    */

    if( _automargin ) {
	// Increase margins to fit rulers inside allocated frame size
	_tmargin[0] += prec[0]-MIN3(bbox[0][0],bbox[1][0],bbox[2][0]);
	_tmargin[1] += MAX3(bbox[0][3],bbox[1][3],bbox[3][3])-prec[3];
	_tmargin[2] += MAX3(bbox[0][2],bbox[2][2],bbox[3][2])-prec[2];
	_tmargin[3] += prec[1]-MIN3(bbox[1][1],bbox[2][1],bbox[3][1]);

	// Increase margins to fit in colormap legend
	if( _cm_legend && _cml_enable ) {
	    double lw, lh;
	    _cm_legend->set_height( prec[3]-prec[1] );
	    _cm_legend->get_size( cairo, lw, lh );
	    _tmargin[2] += lw + _cmlspace;
	}

	// Increase margins to fit in plot title
	if( _title.get_text() != "" ) {
	    double bb[4];
	    _title.get_bbox( cairo, bb );
	    _tmargin[3] += bb[3]-bb[1] + _titlespace;
	}
    }

    bool autorbck[4];
    if( _fixedaspect != PLOT_FIXED_ASPECT_DISABLED ) {
	// Fixed aspect ratio routines assume x1 and y1 to be used.
	double range[4];
	_ruler[0].get_ranges( range[0], range[2] );
	_ruler[1].get_ranges( range[1], range[3] );

	// Disable ruler autorange temporarily
	_ruler[0].get_autorange( autorbck[0], autorbck[1] );
	_ruler[1].get_autorange( autorbck[2], autorbck[3] );
	_ruler[0].set_autorange( false, false );
	_ruler[1].set_autorange( false, false );

	if( _fixedaspect == PLOT_FIXED_ASPECT_EXTEND_RANGE ) {
	    // Modify ranges to have 1:1 aspect ratio. Algorithm keeps
	    // center point fixed and selects the axis to be modified
	    // so that the shown area does not decrease.
	    double xx = (_width - _tmargin[0] - _tmargin[2]) / (range[2]-range[0]);
	    double yy = (_height - _tmargin[1] - _tmargin[3]) / (range[3]-range[1]);
	    if( yy > xx  ) {
		// Increase y _range
		double ran = range[3]-range[1];
		double inc = (_height - _tmargin[1] - _tmargin[3])/xx - ran;
		range[3] += 0.5*inc;
		range[1] -= 0.5*inc;
		_ruler[1].set_ranges( range[1], range[3] );
	    } else {
		// Increase x _range
		double ran = range[2]-range[0];
		double inc = (_width - _tmargin[0] - _tmargin[2])/yy - ran;
		range[2] += 0.5*inc;
		range[0] -= 0.5*inc;
		_ruler[0].set_ranges( range[0], range[2] );
	    }
	} else {
	    // Increase margin either on right side or bottom to have
	    // 1:1 aspect ratio
	    double xx = (_width - _tmargin[0] - _tmargin[2]) / (range[2]-range[0]);
	    double yy = (_height - _tmargin[1] - _tmargin[3]) / (range[3]-range[1]);
	    if( yy > xx  ) {
		// Increase margin on bottom
		_tmargin[1] = _height - _tmargin[3] - xx*(range[3]-range[1]);
	    } else {
		// Increase margin on right side
		_tmargin[2] = _width - _tmargin[0] - yy*(range[2]-range[0]);
	    }
	}
    }

    if( _automargin ) {
	// Recalculate corner coordinates
	prec[0] = _offx + _tmargin[0];
	prec[1] = _offy + _tmargin[3];
	prec[2] = _offx + _width - _tmargin[2];
	prec[3] = _offy + _height - _tmargin[1];

	// Set ruler end points
	_ruler[0].set_endpoints( prec[0], prec[3], prec[2], prec[3] );
	_ruler[1].set_endpoints( prec[0], prec[3], prec[0], prec[1] );
	_ruler[2].set_endpoints( prec[0], prec[1], prec[2], prec[1] );
	_ruler[3].set_endpoints( prec[2], prec[3], prec[2], prec[1] );
    }

    // Reset autorange values
    for( int a = 0; a < 4; a++ )
	_ruler[a].set_ranges( _range_min[a], _range_max[a] );

    calculate_rulers( cairo, true );
    mirror_rulers();

    // Set title position
    _ruler[2].get_bbox( cairo, bbox[2], _cm[2], false );
    _title.set_location( _offx+0.5*_width, bbox[2][1]-_titlespace );
    _title.set_alignment( 0.5, 0.0 );

    // Return autorange settings
    if( _fixedaspect != PLOT_FIXED_ASPECT_DISABLED ) {
	_ruler[0].set_autorange( autorbck[0], autorbck[1] );
	_ruler[1].set_autorange( autorbck[2], autorbck[3] );
    }
}


void Frame::mirror_rulers( void )
{
    // Set enable/disable and do mirroring of rulers if necessary

#ifdef DEBUG_FRAME
    std::cout << "\nFRAME::MIRROR_RULERS()\n\n";

    for( size_t a = 0; a < 4; a++ ) 
	std::cout << "_enable[" << a << "] = " << _enable[a] << "\n";
#endif

    // X-axes
    if( !_enable[0] && !_enable[2] ) {
	// Both disabled
	_ruler[0].enable_labels( false );
	_ruler[2].enable_labels( false );
    } else {
	if( _enable[0] )
	    _ruler[0].enable_labels( true );
	else {
	    // Mirror 0 = 2
	    _ruler[0].copy_tics( _ruler[2] );
	    _cm[0] = _cm[2];
	    _ruler[0].enable_labels( false );
	}
	if( _enable[2] )
	    _ruler[2].enable_labels( true );
	else {
	    // Mirror 2 = 0
	    _ruler[2].copy_tics( _ruler[0] );
	    _cm[2] = _cm[0];
	    _ruler[2].enable_labels( false );
	}
    }

    // Y-axes
    if( !_enable[1] && !_enable[3] ) {
	_ruler[1].enable_labels( false );
	_ruler[3].enable_labels( false );
    } else {
	if( _enable[1] )
	    _ruler[1].enable_labels( true );
	else {
	    // Mirror 1 = 3
	    _ruler[1].copy_tics( _ruler[3] );
	    _cm[1] = _cm[3];
	    _ruler[1].enable_labels( false );
	}
	if( _enable[3] )
	    _ruler[3].enable_labels( true );
	else {
	    // Mirror 3 = 1
	    _ruler[3].copy_tics( _ruler[1] );
	    _cm[3] = _cm[1];
	    _ruler[3].enable_labels( false );
	}
    }

#ifdef DEBUG_FRAME
    for( size_t a = 0; a < 4; a++ ) {
	std::cout << "_ruler[" << a << "].debug_print():\n";
	_ruler[a].debug_print( std::cout );
	_cm[a].debug_print( std::cout );
    }
#endif
}


void Frame::draw_frame( cairo_t *cairo )
{
#ifdef DEBUG_FRAME
    std::cout << "\nFRAME::DRAW_FRAME()\n\n";
#endif

    // Draw rulers
    for( size_t a = 0; a < 4; a++ ) {
	//std::cout << "_ruler[" << a << "].debug_print() BEFORE DRAW:\n";
	//_ruler[a].debug_print( std::cout );

	_ruler[a].draw( cairo, _cm[a], false );

	//std::cout << "_ruler[" << a << "].debug_print() AFTER DRAW:\n";
	//_ruler[a].debug_print( std::cout );
    }

    // Draw title
    _title.draw( cairo );
}


void Frame::enable_colormap_legend( bool enable )
{
    _cml_enable = enable;
}


void Frame::build_colormap_legend( void )
{
    if( _cm_legend )
	delete _cm_legend;
    _cm_legend = NULL;

    size_t a;
    Colormap *cm;
    for( a = 0; a < _dobj.size(); a++ ) {
	cm = dynamic_cast<Colormap *>( _dobj[a]._graph );
	if( cm && _cml_enable )
	    break;
    }
    if( a == _dobj.size() )
	// Nothing to do
	return;

    _cm_legend = new ColormapLegend( *cm );
    _cm_legend->set_font_size( _fontsize );
}


void Frame::draw_colormap_legend( cairo_t *cairo )
{
    if( !_cm_legend || !_cml_enable )
	return;

    double bbox[4];
    _ruler[3].get_bbox( cairo, bbox, _cm[3], false );
    _cm_legend->set_height( _height-_tmargin[1]-_tmargin[3] );
    _cm_legend->plot( cairo, bbox[2] + _cmlspace,
		      _offy+_height-_tmargin[1] );
}


void Frame::draw_legend( cairo_t *cairo )
{
    //std::cout << "draw_legend()\n";

    if( !_legend_enable )
	return;

    double lw, lh;
    _legend.get_size( cairo, lw, lh );

    double x = 0.0;
    double y = 0.0;

    if( (_legend_pos & LEGEND_POS_VERTICAL_MASK) == LEGEND_POS_BOTTOM )
	y = _offy+_height-_tmargin[1];
    else if( (_legend_pos & LEGEND_POS_VERTICAL_MASK) == LEGEND_POS_MIDDLE )
	y = _offy+_tmargin[3]+0.5*(_height-_tmargin[1]-_tmargin[3]+lh);
    else if( (_legend_pos & LEGEND_POS_VERTICAL_MASK) == LEGEND_POS_TOP )
	y = _offy+_tmargin[3]+lh;

    if( (_legend_pos & LEGEND_POS_HORIZONTAL_MASK) == LEGEND_POS_LEFT )
	x = _offx+_tmargin[0];
    else if( (_legend_pos & LEGEND_POS_HORIZONTAL_MASK) == LEGEND_POS_CENTER )
	x = _offx+_tmargin[0]+0.5*(_width-_tmargin[0]-_tmargin[2]-lw);
    else if( (_legend_pos & LEGEND_POS_HORIZONTAL_MASK) == LEGEND_POS_RIGHT )
	x = _offx+_width-_tmargin[2]-lw;

    /*
    std::cout << "  lw = " << lw << "\n";
    std::cout << "  lh = " << lh << "\n";

    std::cout << "  offx = " << _offx << "\n";
    std::cout << "  offy = " << _offy << "\n";

    std::cout << "  width = " << _width << "\n";
    std::cout << "  height = " << _height << "\n";

    std::cout << "  x = " << x << "\n";
    std::cout << "  y = " << y << "\n";
    */

    _legend.plot( cairo, x, y );
}


void Frame::draw( cairo_t *cairo )
{
#ifdef DEBUG_FRAME
    std::cout << "\nFRAME::DRAW()\n\n";
    std::cout << "width = " << _width << "\n";
    std::cout << "height = " << _height << "\n";
#endif

    // Build colormap legend
    build_colormap_legend();

    // Draw background
    cairo_rectangle( cairo, _offx, _offy, _width, _height );
    cairo_set_source_rgb( cairo, _bg[0], _bg[1], _bg[2] );
    cairo_fill( cairo );

    // Get drawable bounding boxes and set ruler ranges
    calculate_autoranging();

    // Calculate frame location and size
    calculate_frame( cairo );

    // Draw contents with clipping on
    set_frame_clipping( cairo );
    for( size_t a = 0; a < _dobj.size(); a++ ) {

        cairo_save( cairo );
	PlotAxis xaxis = _dobj[a]._xaxis;
	PlotAxis yaxis = _dobj[a]._yaxis;
	double range[4];
	_ruler[xaxis].get_ranges( range[0], range[2] );
	_ruler[yaxis].get_ranges( range[1], range[3] );
	Coordmapper cm( _cm[xaxis], _cm[yaxis] );
        _dobj[a]._graph->plot( cairo, &cm, range );
        cairo_restore( cairo );
    }
    unset_frame_clipping( cairo );

    draw_legend( cairo );
    draw_colormap_legend( cairo );

    // Draw frame (on top of user drawn image)
    draw_frame( cairo );
}

