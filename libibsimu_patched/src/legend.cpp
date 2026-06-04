/*! \file legend.cpp
 *  \brief Plot legends
 */

/* Copyright (c) 2005-2009,2011-2013,2015 Taneli Kalvas. All rights reserved.
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

#include <math.h>
#include "legend.hpp"


// Legend sizes are relative to font size
#define LEGEND_SAMPLE_WIDTH 3.0
#define LEGEND_SAMPLE_HEIGHT 1.0
#define LEGEND_SAMPLE_LINESKIP 1.5
#define LEGEND_SAMPLE_SEPARATION 0.3
#define LEGEND_MARGIN 0.5


#define COLORMAP_LEGEND_WIDTH 1.2


void LegendEntry::plot( cairo_t *cairo, double x, double y )
{
    double scale = _label.get_font_size();
    _graph.plot_sample( cairo, x, y, scale*LEGEND_SAMPLE_WIDTH, scale*LEGEND_SAMPLE_HEIGHT );
    _label.set_location( x + scale*(LEGEND_SAMPLE_WIDTH + LEGEND_SAMPLE_SEPARATION), 
			 y );
    _label.draw( cairo );
}


void LegendEntry::get_size( cairo_t *cairo, double &width, double &height ) 
{
    double bbox[4];
    double scale = _label.get_font_size();
    _label.get_bbox( cairo, bbox );
    width = scale*(LEGEND_SAMPLE_WIDTH + LEGEND_SAMPLE_SEPARATION)
	+ (bbox[2]-bbox[0]);
    height = scale*LEGEND_SAMPLE_HEIGHT;
}


void LegendEntry::set_font_size( double fontsize )
{
    _label.set_font_size( fontsize );
}


/* ********************************************************************* */



MultiEntryLegend::MultiEntryLegend()
    : _fontsize(12.0)
{

}


void MultiEntryLegend::add_entry( LegendEntry *entry )
{
    if( entry )
	entry->set_font_size( _fontsize );
    _entry.push_back( entry );
}


void MultiEntryLegend::clear_entries( void )
{
    _entry.clear();
}


void MultiEntryLegend::plot( cairo_t *cairo, double x, double y )
{
    uint32_t k = 0;
    for( int32_t i = _entry.size()-1; i >= 0; i-- ) {
	if( _entry[i] ) {
	    _entry[i]->plot( cairo, 
			     x + _fontsize*LEGEND_MARGIN, 
			     y - _fontsize*(LEGEND_SAMPLE_LINESKIP*k+LEGEND_MARGIN) );
	    k++;
	}
    }
}


void MultiEntryLegend::get_size( cairo_t *cairo, double &width, double &height )
{
    width = height = 0.0;
    for( uint32_t i = 0; i < _entry.size(); i++ ) {
	double w, h;
	if( _entry[i] ) {
	    _entry[i]->get_size( cairo, w, h );
	    if( w > width )
		width = w;
	    if( i == _entry.size()-1 )
		height += _fontsize*LEGEND_SAMPLE_HEIGHT;
	    else
		height += _fontsize*LEGEND_SAMPLE_LINESKIP;
	}
    }

    height += _fontsize*2.0*LEGEND_MARGIN;
    width += _fontsize*2.0*LEGEND_MARGIN;
}


double MultiEntryLegend::get_font_size( void )
{
    return( _fontsize );
}


void MultiEntryLegend::set_font_size( double fontsize )
{
    _fontsize = fontsize;
    for( uint32_t i = 0; i < _entry.size(); i++ ) {
	if( _entry[i] )
	    _entry[i]->set_font_size( fontsize );
    }
}


/* ********************************************************************* */



ColormapLegend::ColormapLegend( Colormap &colormap ) 
    : _height(0.0), _fontsize(12.0), _color((Vec3D(0,0,0))),
      _ticlen_in(5.0), _ticlen_out(5.0),
      _ticspace(5.0), _colormap(colormap) 
{

}


void ColormapLegend::build_legend( double x, double y )
{
    // Make 6 tics
    _tic.clear();
    size_t N = 6;
    int sign;
    double zmin = _colormap.zscale_inv( 0.0 );
    double zmax = _colormap.zscale_inv( 1.0 );
    double zspan = zmax - zmin;
    if( zmin >= -1.0e-6*zspan && zmax >= 0.0 ) {
	sign = +1;
    } else if( zmax <= 1.0e-6*zspan && zmin <= 0.0 ) {
	sign = -1;
    } else {
	sign = 0;
    }
    for( size_t a = 0; a < N; a++ ) {

	// zval should depend on zscale settings
	double xx = x + _width + 2*_ticlen_out + _ticspace;
	double yy = y - _height*a/(N-1.0);

	char str[128];
	double zval = _colormap.zscale_inv( a/(N-1.0) );
	if( sign == +1 && fabs(zval) < 1e-6*fabs(zmax) )
	    zval = 0.0;
	else if( sign == -1 && fabs(zval) < 1e-6*fabs(zmin) )
	    zval = 0.0;
	snprintf( str, 128, "%.4g", zval );

	_tic.push_back( Tic(yy,str) );
	_tic[a]._label.set_font_size( _fontsize );
	_tic[a]._label.set_color( _color );
	_tic[a]._label.set_location( xx, yy );
	_tic[a]._label.set_alignment( 0.0, 0.5, true );
    }
}


void ColormapLegend::plot_colormap_palette_to_image_surface( cairo_surface_t *surface, int plim[4] )
{
    unsigned char *buf = cairo_image_surface_get_data( surface );
	
    int width  = cairo_image_surface_get_width( surface );
    int height = cairo_image_surface_get_height( surface );
    int stride = cairo_image_surface_get_stride( surface );

    // Check that drawing is only done to valid buffer
    if( plim[0] < 0 )
	plim[0] = 0;
    if( plim[1] < 0 )
	plim[1] = 0;
    if( plim[2] >= width )
	plim[2] = width-1;
    if( plim[3] >= height )
	plim[3] = height-1;
    if( plim[0] >= plim[2] || plim[1] >= plim[3] )
	return;
	/*
	throw( Error( ERROR_LOCATION, (std::string)"incorrect pixel limits: " + 
		      "plim[0] = " + to_string(plim[0]) + 
		      ", plim[1] = " + to_string(plim[1]) +
		      ", plim[2] = " + to_string(plim[2]) + 
		      ", plim[3] = " + to_string(plim[3]) ) );
	*/

    // Flush to ensure all writing to the image was done
    cairo_surface_flush( surface );

    const Palette &palette = _colormap.palette();
    for( int j = plim[1]; j <= plim[3]; j++ ) {

	// Transform j to color
	Vec3D c = palette( (double)(j-plim[3]) / (double)(plim[1]-plim[3]) );
	
	for( int i = plim[0]; i <= plim[2]; i++ ) {

	    buf[j*stride+4*i+0] = (unsigned char)(255*c[2]);  // Blue
	    buf[j*stride+4*i+1] = (unsigned char)(255*c[1]);  // Green
	    buf[j*stride+4*i+2] = (unsigned char)(255*c[0]);  // Red
	    buf[j*stride+4*i+3] = (unsigned char)255;         // Alpha
	}
    }

    // Mark the image dirty so cairo clears its caches.
    cairo_surface_mark_dirty( surface );
}


void ColormapLegend::plot_colomap_palette( cairo_t *cairo, int plim[4] )
{
    // Prepare by fetching surface and its parameters
    cairo_surface_t *surface = cairo_get_target( cairo );
    if( cairo_image_surface_get_data( surface ) == NULL ) {

	// Make a new image surface for plot
	int width = plim[2]-plim[0]+1;
	int height = plim[3]-plim[1]+1;
	cairo_surface_t *nsurface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32,
								width,
								height );
	int nlim[4] = { 0, 0, width-1, height-1 };
	plot_colormap_palette_to_image_surface( nsurface, nlim );

	cairo_save( cairo );
	cairo_set_source_surface( cairo, nsurface, plim[0], plim[1] );
	cairo_pattern_t *pattern = cairo_get_source( cairo );
	cairo_pattern_set_filter( pattern, CAIRO_FILTER_GOOD );
	cairo_rectangle( cairo, plim[0], plim[1], width, height );
	cairo_clip( cairo );
	cairo_paint( cairo );
	cairo_restore( cairo );
	cairo_surface_destroy( nsurface );

    } else {

	// Plot directly to image surface
	cairo_format_t format = cairo_image_surface_get_format( surface );
	if( format != CAIRO_FORMAT_ARGB32 )
	    throw( Error( ERROR_LOCATION, "cairo image surface format not supported" ) );

	plot_colormap_palette_to_image_surface( surface, plim );
    }
}


// The point (x,y) is the lower left corner of the whole legend, not
// the location of the corner of the rectangle.
//
void ColormapLegend::plot( cairo_t *cairo, double x, double y )
{
    build_legend( x, y );

    // Draw colormap palette sample
    int plim[4] = { (int)floor(x+_ticlen_out+0.5),
		    (int)floor(y-_height+0.5),
		    (int)floor(x+_ticlen_out+_width),
		    (int)floor(y) };
    plot_colomap_palette( cairo, plim );

    // Draw box and tics
    cairo_save( cairo );
    cairo_set_source_rgba( cairo, _color[0], _color[1], _color[2], 1.0 );
    cairo_set_line_width( cairo, 1.0 );

    cairo_rectangle( cairo, x + _ticlen_out, y, 
		     _width, -_height );
    for( size_t a = 0; a < _tic.size(); a++ ) {
	cairo_move_to( cairo, x, _tic[a]._loc );
	cairo_line_to( cairo, x + _ticlen_out + _ticlen_in, _tic[a]._loc );
	cairo_move_to( cairo, x + _ticlen_out + _width - _ticlen_in, _tic[a]._loc );
	cairo_line_to( cairo, x + 2*_ticlen_out + _width, _tic[a]._loc );
    }
    cairo_stroke( cairo );

    for( size_t a = 0; a < _tic.size(); a++ )
	_tic[a]._label.draw( cairo );

    cairo_restore( cairo );
}


void ColormapLegend::get_size( cairo_t *cairo, double &width, double &height )
{
    build_legend( 0.0, 0.0 );

    double maxwidth = _width + _ticlen_out*2;
    for( size_t a = 0; a < _tic.size(); a++ ) {
	double bbox[4];
	_tic[a]._label.get_bbox( cairo, bbox );
	if( bbox[2] > maxwidth )
	    maxwidth = bbox[2];
    }

    width = maxwidth;
    height = _height;
}


void ColormapLegend::set_font_size( double fontsize )
{
    _fontsize   = fontsize;
    _ticlen_in  = 5.0*fontsize/12.0;
    _ticlen_out = 5.0*fontsize/12.0;
    _ticspace   = 5.0*fontsize/12.0;
    _width      = COLORMAP_LEGEND_WIDTH*fontsize;
}


void ColormapLegend::set_height( double height )
{
    _height = height;
}

