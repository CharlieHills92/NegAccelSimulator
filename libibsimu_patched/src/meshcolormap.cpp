/*! \file meshcolormap.cpp
 *  \brief Mesh based colormap graph for plotting
 */

/* Copyright (c) 2005-2012,2014,2015 Taneli Kalvas. All rights reserved.
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

#include <cmath>
#include <limits>
#include <iomanip>
#include <iostream>
#include "compmath.hpp"
#include "meshcolormap.hpp"


//#define DEBUG_COLORMAP 1


MeshColormap::MeshColormap()
    : _interpolation(INTERPOLATION_BILINEAR), _zscale(ZSCALE_LINEAR), 
      _zmin(0.0), _zmax(0.0), _n(0), _m(0), _intrp(NULL), _zscale_prepared(false)
{

}


MeshColormap::MeshColormap( const double datarange[4], size_t n, size_t m, 
		    const std::vector<double> &data )
    : _interpolation(INTERPOLATION_BILINEAR), _zscale(ZSCALE_LINEAR), 
      _zmin(0.0), _zmax(0.0), _n(0), _m(0), _intrp(NULL), _zscale_prepared(false)
{
    set_data( datarange, n, m, data );
}


MeshColormap::MeshColormap( const MeshColormap &colormap )
    : _palette(colormap._palette), _interpolation(colormap._interpolation),
      _zscale(colormap._zscale), _zmin(colormap._zmin), _zmax(colormap._zmax),
      _n(colormap._n), _m(colormap._m), _f(colormap._f), _intrp(NULL), 
      _zscale_prepared(false)
{
    _datarange[0] = colormap._datarange[0];
    _datarange[1] = colormap._datarange[1];
    _datarange[2] = colormap._datarange[2];
    _datarange[3] = colormap._datarange[3];

    prepare_data_interpolation();
}


MeshColormap::~MeshColormap()
{
    if( _intrp )
	delete _intrp;
}


void MeshColormap::clear_data( void ) 
{
    _n = 0;
    _m = 0;
    if( _intrp )
	delete _intrp;
    _intrp = NULL;
    _f.clear();
}


void MeshColormap::set_data( const double datarange[4], size_t n, size_t m, 
			 const std::vector<double> &data )
{
#ifdef DEBUG_COLORMAP
    std::cout << "MeshColormap::set_data()\n";
#endif

    _n = n;
    _m = m;
    if( _intrp )
	delete _intrp;
    _intrp = NULL;
    if( _n == 0 || _m == 0 ) {
	_f.clear();
	return;
    }

    _datarange[0] = datarange[0];
    _datarange[1] = datarange[1];
    _datarange[2] = datarange[2];
    _datarange[3] = datarange[3];

    if( n*m != data.size() )
	throw( Error( ERROR_LOCATION, "data size not equal to n*m" ) );
    _f = data;

    // Default range is from minimum to maximum
    _zmin = std::numeric_limits<double>::infinity();
    _zmax = -std::numeric_limits<double>::infinity();
    for( size_t j = 0; j < _m; j++ ) {
	for( size_t i = 0; i < _n; i++ ) {
	    int p = i+j*_n;
	    if( _f[p] < _zmin )
		_zmin = _f[p];
	    if( _f[p] > _zmax )
		_zmax = _f[p];
	}
    }
    if( _zmin == _zmax ) {
	_zmax += 1.0;
	//_zmin -= 1.0;
    }

    _zscale_prepared = false;
    prepare_data_interpolation();
}


interpolation_e MeshColormap::get_interpolation( void ) const
{
    return( _interpolation );
}


void MeshColormap::set_interpolation( interpolation_e interpolation )
{
    _interpolation = interpolation;
    prepare_data_interpolation();
}


zscale_e MeshColormap::get_zscale( void ) const
{
    return( _zscale );
}


void MeshColormap::set_zscale( zscale_e zscale )
{
    _zscale = zscale;
    _zscale_prepared = false;
}


void MeshColormap::prepare_zscaling( void )
{
    // Error if either end is at zero with LOG scaling
    if( _zscale == ZSCALE_LOG && _zmin <= 0.0 && _zmax >= 0.0 )
	throw( Error( ERROR_LOCATION, "zmin and zmax on different sides of zero and using logscale" ) );

    double zspan = _zmax - _zmin;
    if( _zmin >= -1.0e-6*zspan && _zmax >= 0.0 ) {
	// Completely on positive side
	_sign = +1;
    } else if( _zmax <= 1.0e-6*zspan && _zmin <= 0.0 ) {
	// Completely on negative side
	_sign = -1;
    } else {
	// Both negative and positive
	_sign = 0;
    }

    // Calculate scaling coefficients
    if( _zscale == ZSCALE_LINEAR ) {
	_scale_A = _zmin;
	_scale_B = _zmax;
    } else if( _zscale == ZSCALE_LOG ) {
	if( _sign > 0 ) {
	    _scale_A = log(_zmin);
	    _scale_B = log(_zmax);
	} else {
	    _scale_A = log(-_zmin);
	    _scale_B = log(-_zmax);
	}
    } else if( _zscale == ZSCALE_RELLOG ) {
	_scale_A = log(0.001);
	_scale_B = log(1.001);
    } else {
	_scale_A = 0.0;
	_scale_B = 1.0;
    }

    _scale_C = 1.0/(_scale_B-_scale_A);
}



double MeshColormap::zscale_inv( double val )
{
    if( !_zscale_prepared )
	prepare_zscaling();
    _zscale_prepared = true;

    if( _zscale == ZSCALE_LINEAR ) {
	return( val/_scale_C + _scale_A );
    } else if( _zscale == ZSCALE_LOG ) {
	if( _sign > 0 )
	    return( exp(val/_scale_C + _scale_A) );
	else /* if( sign < 0 ) */
	    return( -exp(val/_scale_C + _scale_A) );
    } else if( _zscale == ZSCALE_RELLOG ) {
	if( _sign > 0 )
	    return( _zmin+(_zmax-_zmin)*(exp(val/_scale_C+_scale_A)-0.001) );
	else if( _sign < 0 )
	    return( _zmax+(_zmin-_zmax)*(exp((1.0-val)/_scale_C+_scale_A)-0.001) );
	else {
	    if( val > 0.5 )
		return( 0.0+(_zmax-0.0)*(exp(2.0*(val-0.5)/_scale_C+_scale_A)-0.001) );
	    else
		return( 0.0+(_zmin-0.0)*(exp((1.0-2.0*val)/_scale_C+_scale_A)-0.001) );
	}
    }

    return( val );
}


double MeshColormap::zscale( double val )
{
    if( !_zscale_prepared )
	prepare_zscaling();
    _zscale_prepared = true;

    if( _zscale == ZSCALE_LINEAR )
	return( _scale_C*(val-_scale_A) );
    else if( _zscale == ZSCALE_LOG ) {
	if( _sign > 0 )
	    return( _scale_C*(log(val)-_scale_A) );
	else /* if( sign < 0 ) */
	    return( _scale_C*(log(-val)-_scale_A) );
    } else if( _zscale == ZSCALE_RELLOG ) {
	if( _sign > 0 ) {
	    double x = (val-_zmin)/(_zmax-_zmin);
	    return( _scale_C*(log(0.001+x)-_scale_A) );
	} else if( _sign < 0 ) {
	    double x = (val-_zmax)/(_zmin-_zmax);
	    return( 1.0+_scale_C*(_scale_A-log(0.001+x)) );
	} else {
	    if( val > 0.0 ) {
		double x = (val-0.0)/(_zmax-0.0);
		return( 0.5+0.5*_scale_C*(log(0.001+x)-_scale_A) );
	    } else {
		double x = (val-0.0)/(_zmin-0.0);
		return( 0.5+0.5*_scale_C*(_scale_A-log(0.001+x)) );
	    }
	}
    }

    return( val );
}


void MeshColormap::plot_to_image_surface( cairo_surface_t *surface, const Coordmapper *cm, int plim[4] )
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

    if( !_intrp )
	throw( Error( ERROR_LOCATION, "no data available" ) );

    // Flush to ensure all writing to the image was done
    cairo_surface_flush( surface );

    for( int j = plim[1]; j <= plim[3]; j++ ) {
	for( int i = plim[0]; i <= plim[2]; i++ ) {

	    // Transform to logical coordinates
	    double x[2] = { (double)i, (double)j };
	    cm->inv_transform( x[0], x[1] );
	    double val = zscale( get_value( x[0], x[1] ) );
	    if( comp_isinf( val ) || comp_isnan( val ) )
		continue;

	    Vec3D c = _palette( val );
	    buf[j*stride+4*i+0] = (unsigned char)(255*c[2]);  // Blue
	    buf[j*stride+4*i+1] = (unsigned char)(255*c[1]);  // Green
	    buf[j*stride+4*i+2] = (unsigned char)(255*c[0]);  // Red
	    buf[j*stride+4*i+3] = (unsigned char)255;         // Alpha
	}
    }

    // Mark the image dirty so cairo clears its caches.
    cairo_surface_mark_dirty( surface );
}


void MeshColormap::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] )
{
#ifdef DEBUG_COLORMAP
    std::cout << "MeshColormap::plot()\n";
#endif

    // If colormap empty, do nothing
    if( _n == 0 || _m == 0 )
	return;

#ifdef DEBUG_COLORMAP
    std::cout << "datarange[0] = " << _datarange[0] << "\n"
	      << "datarange[1] = " << _datarange[1] << "\n"
	      << "datarange[2] = " << _datarange[2] << "\n"
	      << "datarange[3] = " << _datarange[3] << "\n\n";
    std::cout << "range[0] = " << range[0] << "\n"
	      << "range[1] = " << range[1] << "\n"
	      << "range[2] = " << range[2] << "\n"
	      << "range[3] = " << range[3] << "\n\n";
#endif

    // Range of plot limited by current plot range and data range -> final range
    double frange[4];
    if( _datarange[0] < range[0] )
	frange[0] = range[0];
    else
	frange[0] = _datarange[0];
    if( _datarange[1] < range[1] )
	frange[1] = range[1];
    else
	frange[1] = _datarange[1];

    if( _datarange[2] > range[2] )
	frange[2] = range[2];
    else
	frange[2] = _datarange[2];
    if( _datarange[3] > range[3] )
	frange[3] = range[3];
    else
	frange[3] = _datarange[3];

    // Calculate pixel ranges
    double prange[4];
    cm->transform( &prange[0], &frange[0] );
    cm->transform( &prange[2], &frange[2] );

    // Calculate pixel integer limits of drawn area, 
    // y flipped to have smaller numbers as 0 and 1, 
    // bigger as 2 and 3.
    int plim[4] = { (int)floor(prange[0]+0.5),
		    (int)floor(prange[3]+0.5),
		    (int)floor(prange[2]+0.5),
		    (int)floor(prange[1]+0.5) };

    if( plim[2] <= plim[0] || plim[3] <= plim[1] ) {
	// If drawing area inverted, do nothing
	return;
    }

    // Prepare by fetching surface and its parameters
    cairo_surface_t *surface = cairo_get_target( cairo );
    if( cairo_image_surface_get_data( surface ) == NULL ) {

	// Make a new image surface for plot
	int width = plim[2]-plim[0]+1;
	int height = plim[3]-plim[1]+1;
	cairo_surface_t *nsurface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32,
								width,
								height );
	// frange is clipped range in simulation coordinates
	// transform from simu coord to pixels in nsurface.
	double xx = (width-1)/(frange[2]-frange[0]);
	double x0 = -xx*frange[0];
	double yy = (height-1)/(frange[1]-frange[3]);
	double y0 = -yy*frange[3];
	Coordmapper *ncm = new Coordmapper( xx, x0, yy, y0 );
	int nlim[4] = { 0, 0, width-1, height-1 };
	plot_to_image_surface( nsurface, ncm, nlim );
	delete ncm;

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

	plot_to_image_surface( surface, cm, plim );
    }
}


void MeshColormap::plot_sample( cairo_t *cairo, double x, double y, double width, double height )
{

}


void MeshColormap::get_bbox( double bbox[4] )
{
    bbox[0] = _datarange[0];
    bbox[1] = _datarange[1];
    bbox[2] = _datarange[2];
    bbox[3] = _datarange[3];
}


void MeshColormap::set_palette( const Palette &palette )
{
    _palette = palette;
}


void MeshColormap::get_zrange( double &min, double &max ) const
{
    min = _zmin;
    max = _zmax;
}


void MeshColormap::set_zrange( double min, double max )
{
    if( min < max ) {
	_zmin = min;
	_zmax = max;
    } else {
	_zmin = max;
	_zmax = min;
    }
    _zscale_prepared = false;
}


void MeshColormap::prepare_data_interpolation( void )
{
    // Free old interpolation
    if( _intrp )
	delete _intrp;

    // Make a new interpolation
    switch( _interpolation ) {
    case INTERPOLATION_CLOSEST:
	_intrp = new ClosestInterpolation2D( _n, _m, _f );
	break;
    case INTERPOLATION_BILINEAR:
	_intrp = new BiLinearInterpolation2D( _n, _m, _f );
	break;
    case INTERPOLATION_BICUBIC:
	_intrp = new BiCubicInterpolation2D( _n, _m, _f );
	break;
    default:
	throw( Error( ERROR_LOCATION, "unknown interpolation type" ) );	
    }
}


double MeshColormap::get_value( double x, double y ) const
{	    
    // Calculate relative point in data
    double t = (x-_datarange[0])/(_datarange[2]-_datarange[0]);
    double u = (y-_datarange[1])/(_datarange[3]-_datarange[1]);
    
    return( (*_intrp)( t, u ) );
}

