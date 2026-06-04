/*! \file xygraph.cpp
 *  \brief XY-graph
 */

/* Copyright (c) 2005-2011,2013 Taneli Kalvas. All rights reserved.
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

#include <limits>
#include <iostream>
#include <cmath>
#include "compmath.hpp"
#include "xygraph.hpp"


XYGraph::XYGraph()
    : _linewidth(1.0), _color(Vec3D(1,0,0)), _linestyle(XYGRAPH_LINE_DISABLE), 
      _pointstyle(XYGRAPH_POINT_CIRCLE),
      _point_filled(true), _point_scale(3.0),
      _histogram(false), _extend_histogram(false)
{
}


XYGraph::XYGraph( const std::vector<double> &xdata, 
		  const std::vector<double> &ydata )
    : _linewidth(1.0), _color(Vec3D(1,0,0)), _linestyle(XYGRAPH_LINE_DISABLE), 
      _pointstyle(XYGRAPH_POINT_CIRCLE),
      _point_filled(true), _point_scale(3.0), 
      _histogram(false), _extend_histogram(false),
      _xdata(xdata), _ydata(ydata)
{
}


void XYGraph::plot_point( cairo_t *cairo, double x, double y )
{
    cairo_save( cairo );
    cairo_translate( cairo, x, y );
    cairo_scale( cairo, _point_scale, _point_scale );
    
    if( _pointstyle == XYGRAPH_POINT_CIRCLE ) {
	cairo_move_to( cairo, 1.0, 0.0 );
	cairo_arc( cairo, 0.0, 0.0, 1.0, 0.0, 2.0*M_PI );
    } else if( _pointstyle == XYGRAPH_POINT_BOX ) {
	cairo_rectangle( cairo, -1.0, -1.0, 2.0, 2.0 );
    }

    cairo_restore( cairo );
    if( _point_filled )
	cairo_fill( cairo );
    else
	cairo_stroke( cairo );
}


void XYGraph::plot_standard_lines( cairo_t *cairo, const Coordmapper *cm, const double range[4], size_t N )
{
    if( _linestyle == XYGRAPH_LINE_SOLID ) {
	bool cont = false;
	for( size_t a = 0; a < N; a++ ) {
	    double x = _xdata[a];
	    double y = _ydata[a];
	    if( comp_isnan(x) || comp_isnan(y) ) {
		cont = false;
		continue;
	    }
	    cm->transform( x, y );
	    if( cont )
		cairo_line_to( cairo, x, y );
	    else
		cairo_move_to( cairo, x, y );
	    cont = true;
	}
	cairo_stroke( cairo );
    }
}


void XYGraph::plot_histogram_lines( cairo_t *cairo, const Coordmapper *cm, const double range[4], size_t N )
{
    if( N <= 1 || _linestyle != XYGRAPH_LINE_SOLID )
	return;

    double zx = range[0];
    double zy = 0.0;
    cm->transform( zx, zy );
    double ox = _xdata[0];
    double oy = _ydata[0];
    cm->transform( ox, oy );
    double x = _xdata[1];
    double y = _ydata[1];
    cm->transform( x, y );
    double w = x-ox;
    if( _extend_histogram ) {
	cairo_move_to( cairo, zx, zy );
	cairo_line_to( cairo, ox-0.5*w, zy );
    } else 
	cairo_move_to( cairo, ox-0.5*w, zy );
    cairo_line_to( cairo, ox-0.5*w, oy );

    for( size_t a = 1; a < N; a++ ) {
	x = _xdata[a];
	y = _ydata[a];
	cm->transform( x, y );
	w = x-ox;
	cairo_line_to( cairo, ox+0.5*w, oy );
	cairo_line_to( cairo, ox+0.5*w, y );
	ox = x;
	oy = y;
    }

    zx = range[2];
    zy = 0.0;
    cm->transform( zx, zy );

    cairo_line_to( cairo, x+0.5*w, y );
    cairo_line_to( cairo, x+0.5*w, zy );
    if( _extend_histogram )
	cairo_line_to( cairo, zx, zy );

    cairo_stroke( cairo );
}


void XYGraph::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4])
{
    size_t N = _xdata.size() < _ydata.size() ? _xdata.size() : _ydata.size();

    cairo_set_line_width( cairo, _linewidth );
    cairo_set_source_rgb( cairo, _color[0], _color[1], _color[2] );

    if( _histogram )
	plot_histogram_lines( cairo, cm, range, N );
    else
	plot_standard_lines( cairo, cm, range, N );

    if( _pointstyle != XYGRAPH_POINT_DISABLE ) {
	for( size_t a = 0; a < N; a++ ) {
	    double x = _xdata[a];
	    double y = _ydata[a];
	    if( comp_isnan(x) || comp_isnan(y) )
		continue;

	    cm->transform( x, y );
	    plot_point( cairo, x, y );
	}
    }
}


void XYGraph::plot_sample( cairo_t *cairo, double x, double y, double width, double height )
{
    cairo_set_line_width( cairo, _linewidth );
    cairo_set_source_rgb( cairo, _color[0], _color[1], _color[2] );

    if( _linestyle == XYGRAPH_LINE_SOLID ) {

	cairo_move_to( cairo, x, y-0.5*height );
	cairo_line_to( cairo, x+width, y-0.5*height );
	cairo_stroke( cairo );

    }
    if( _pointstyle != XYGRAPH_POINT_DISABLE ) {

	plot_point( cairo, x+0.5*width, y-0.5*height );
    }
}


void XYGraph::get_bbox( double bbox[4] )
{
    size_t N = _xdata.size() < _ydata.size() ? _xdata.size() : _ydata.size();

    bbox[0] = std::numeric_limits<double>::infinity();
    bbox[1] = std::numeric_limits<double>::infinity();
    bbox[2] = -std::numeric_limits<double>::infinity();
    bbox[3] = -std::numeric_limits<double>::infinity();
    for( size_t a = 0; a < N; a++ ) {
	if( _xdata[a] < bbox[0] )
	    bbox[0] = _xdata[a];
	if( _xdata[a] > bbox[2] )
	    bbox[2] = _xdata[a];
	if( _ydata[a] < bbox[1] )
	    bbox[1] = _ydata[a];
	if( _ydata[a] > bbox[3] )
	    bbox[3] = _ydata[a];
    }
}


void XYGraph::set_data( const std::vector<double> &xdata, 
			const std::vector<double> &ydata )
{
    _xdata = xdata;
    _ydata = ydata;
}


void XYGraph::set_color( const Vec3D &color )
{
    _color = color;
}


void XYGraph::set_line_width( double linewidth )
{
    _linewidth = fabs( linewidth );
}


void XYGraph::set_line_style( line_style_e linestyle, double linewidth )
{
    _linewidth = fabs( linewidth );
    _linestyle  = linestyle;
}


void XYGraph::set_point_style( point_style_e pointstyle, bool filled, double scale )
{
    _pointstyle = pointstyle;
    _point_filled = filled;
    _point_scale = fabs( scale );
}


void XYGraph::set_histogram( bool histo )
{
    _histogram = histo;
}


void XYGraph::extend_histogram( bool extend )
{
    _extend_histogram = extend;
}
