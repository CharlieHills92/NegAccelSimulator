/*! \file softwarerenderer.hpp
 *  \brief Software 3D renderer
 */

/* Copyright (c) 2012,2013,2022 Taneli Kalvas. All rights reserved.
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
#include "softwarerenderer.hpp"
#include "ibsimu.hpp"


SoftwareRenderer::SoftwareRenderer( GtkWidget *darea )
    : _darea(darea), _surface(NULL), 
      _buf(NULL), _zbuf(NULL),
      _material_diffuse_color(0.8,0.0,0.0), 
      _material_ambient_color(0.2,0.0,0.0),
      _color(0.0,0.0,0.0)
{
    ibsimu.message( 1 ) << "Using SoftwareRenderer\n";
    ibsimu.flush();
}


SoftwareRenderer::SoftwareRenderer( cairo_surface_t *surface )
    : _darea(NULL), _surface(surface), 
      _buf(NULL), _zbuf(NULL),
      _material_diffuse_color(0.8,0.0,0.0), 
      _material_ambient_color(0.2,0.0,0.0),
      _color(0.0,0.0,0.0)
{

}

SoftwareRenderer::~SoftwareRenderer()
{
    if( _zbuf )
	delete [] _zbuf;
}


void SoftwareRenderer::set_pixel( int i, int j, const Vec3D &color )
{
    int p = j*_stride+4*i;
    _buf[p+0] = (unsigned char)(255*color[2]);  // Blue
    _buf[p+1] = (unsigned char)(255*color[1]);  // Green
    _buf[p+2] = (unsigned char)(255*color[0]);  // Red
    _buf[p+3] = (unsigned char)255;             // Alpha
}


void SoftwareRenderer::clear( const Vec3D &color )
{
    for( int j = 0; j < _height; j++ ) {
	for( int i = 0; i < _width; i++ ) {
	    set_pixel( i, j, color );
	    _zbuf[i+j*_width] = std::numeric_limits<double>::infinity();
	}
    }
}


void SoftwareRenderer::swap( int &x0, int &y0, double &z0, 
			     int &x1, int &y1, double &z1 )
{
    int xt = x0;
    int yt = y0;
    double zt = z0;
    x0 = x1;
    y0 = y1;
    z0 = z1;
    x1 = xt;
    y1 = yt;
    z1 = zt;
}


void SoftwareRenderer::line_2d( int x0, int y0, double z0,
				int x1, int y1, double z1,
				const Vec3D &color )
{
    if( abs(y1-y0) > abs(x1-x0) ) {

	if( y0 > y1 )
	    swap( x0, y0, z0, x1, y1, z1 );

	int dx = x1-x0;
	int dy = y1-y0;    
	double dz = z1-z0;
	double dyi = 1.0/dy;
	
	int y = y0;
	if( y < 0 )
	    y = 0;
	int ymax = y1;
	if( ymax >= _height )
	    ymax = _height-1;
	while( y <= ymax ) {
	    int x = (int)(dx*(y-y0)*dyi) + x0;
	    double z = dz*(y-y0)*dyi + z0;
	    if( x >= 0 && x < _width && z <= _zbuf[x+y*_width] ) {
		set_pixel( x, y, color );
		_zbuf[x+y*_width] = z;
	    }	    
	    y++;
	}

    } else {

	if( x0 > x1 )
	    swap( x0, y0, z0, x1, y1, z1 );

	int dx = x1-x0;
	int dy = y1-y0;    
	double dz = z1-z0;
	double dxi = 1.0/dx;
	
	int x = x0;
	if( x < 0 )
	    x = 0;
	int xmax = x1;
	if( xmax >= _width )
	    xmax = _width-1;
	while( x <= xmax ) {
	    int y = (int)(dy*(x-x0)*dxi) + y0;
	    double z = dz*(x-x0)*dxi + z0;
	    if( y >= 0 && y < _height && z <= _zbuf[x+y*_width] ) {
		set_pixel( x, y, color );
		_zbuf[x+y*_width] = z;
	    }	    
	    x++;
	}

    }
}


void SoftwareRenderer::flat_2d_triangle( int x0, int y0, double z0,
					 int x1, int y1, double z1,
					 int x2, int y2, double z2,
					 const Vec3D &color )
{
    // Clip by x
    if( (x0 < 0 && x1 < 0 && x2 < 0) ||
	(x0 >= _width && x1 >= _width && x2 >= _width) )
	return;

    // Sort by y
    if( y0 > y2 )
	swap( x0, y0, z0, x2, y2, z2 );
    if( y0 > y1 )
	swap( x0, y0, z0, x1, y1, z1 );
    if( y1 > y2 )
	swap( x1, y1, z1, x2, y2, z2 );

    // Clip by y
    if( y2 < 0 || y0 > _height )
	return;

    // Process from y0 to y1 (or clip)
    int y = y0;
    if( y < 0 )
	y = 0;
    int ymax = y1;
    if( ymax >= _height )
	ymax = _height-1;
    int dx1 = x1-x0;
    int dy1 = y1-y0;    
    double dy1i = 1.0/dy1;
    double dz1 = z1-z0;
    int dx2 = x2-x0;
    int dy2 = y2-y0;
    double dy2i = 1.0/dy2;
    double dz2 = z2-z0;
    int xx1, xx2;
    double zz1, zz2;
    while( y <= ymax ) {
	if( fabs(dy1) <= 0.1 ) {
	    xx1 = x0;
	    xx2 = x1;
	    zz1 = z0;
	    zz2 = z1;
	} else {
	    xx1 = (int)(dx1*(y-y0)*dy1i) + x0;
	    xx2 = (int)(dx2*(y-y0)*dy2i) + x0;
	    zz1 = dz1*(y-y0)*dy1i + z0;
	    zz2 = dz2*(y-y0)*dy2i + z0;
	}
	if( xx1 < xx2 ) {
	    double z = zz1;
	    double dz = (zz2-zz1)/(xx2-xx1);
	    if( xx1 < 0 )
		xx1 = 0;
	    if( xx2 >= _width )
		xx2 = _width-1;
	    for( int x = xx1; x < xx2; x++ ) {
		if( z <= _zbuf[x+y*_width] ) {
		    set_pixel( x, y, color );
		    _zbuf[x+y*_width] = z;
		}
		z += dz;
	    }
	} else {
	    double z = zz2;
	    double dz = (zz1-zz2)/(xx1-xx2);
	    if( xx1 >= _width )
		xx1 = _width-1;
	    if( xx2 < 0 )
		xx2 = 0;
	    for( int x = xx2; x < xx1; x++ ) {
		if( z <= _zbuf[x+y*_width] ) {
		    set_pixel( x, y, color );
		    _zbuf[x+y*_width] = z;
		}
		z += dz;
	    }
	}
	y++;
    }

    // Process from y1 to y2 (or clip)
    ymax = y2;
    if( ymax >= _height )
	ymax = _height-1;
    dx1 = x2-x1;
    dy1 = y2-y1;
    dy1i = 1.0/dy1;
    dz1 = z2-z1;
    while( y <= ymax ) {
	if( fabs(dy1) <= 0.1 ) {
	    xx1 = x1;
	    xx2 = x2;
	    zz1 = z1;
	    zz2 = z2;
	} else {
	    xx1 = (int)(dx1*(y-y1)*dy1i) + x1;
	    xx2 = (int)(dx2*(y-y0)*dy2i) + x0;
	    zz1 = dz1*(y-y1)*dy1i + z1;
	    zz2 = dz2*(y-y0)*dy2i + z0;
	}
	if( xx1 < xx2 ) {
	    double z = zz1;
	    double dz = (zz2-zz1)/(xx2-xx1);
	    if( xx1 < 0 )
		xx1 = 0;
	    if( xx2 >= _width )
		xx2 = _width-1;
	    for( int x = xx1; x < xx2; x++ ) {
		if( z <= _zbuf[x+y*_width] ) {
		    set_pixel( x, y, color );
		    _zbuf[x+y*_width] = z;
		}
		z += dz;
	    }
	} else {
	    double z = zz2;
	    double dz = (zz1-zz2)/(xx1-xx2);
	    if( xx1 >= _width )
		xx1 = _width-1;
	    if( xx2 < 0 )
		xx2 = 0;
	    for( int x = xx2; x < xx1; x++ ) {
		if( z <= _zbuf[x+y*_width] ) {
		    set_pixel( x, y, color );
		    _zbuf[x+y*_width] = z;
		}
		z += dz;
	    }
	}
	y++;
    }
}


void SoftwareRenderer::enable_view_settings( void )
{
    // Prepare matrices
    _modelview = _view*_model;
    _normalmatrix = _modelview.inverse().transpose();
    _totalmatrix = _projection*_modelview;

    // Normalize light direction (infinitely far away model)
    //_light_location = _modelview.transform_vector( _light_location );
    _light_location.normalize();
}


void SoftwareRenderer::start_rendering( void )
{
    if( _darea ) {
	// Make image surface for rendering
	_width = gtk_widget_get_allocated_width( _darea );
	_height = gtk_widget_get_allocated_height( _darea );
	_surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, _width, _height );
    } else {
	_width = cairo_image_surface_get_width( _surface );
	_height = cairo_image_surface_get_height( _surface );
    }
    _buf = cairo_image_surface_get_data( _surface );
    _stride = cairo_image_surface_get_stride( _surface );

    cairo_surface_flush( _surface );

    if( _zbuf )
	delete [] _zbuf;
    _zbuf = new double[_width*_height];

    // Clear window
    clear( Vec3D( 1.0, 1.0, 1.0 ) );
}


void SoftwareRenderer::end_rendering( cairo_t *cairo )
{
    cairo_surface_mark_dirty( _surface );

    if( _darea ) {
	// Copy image surface to window
	cairo_set_source_surface( cairo, _surface, 0, 0 );
	cairo_pattern_set_filter( cairo_get_source(cairo), CAIRO_FILTER_FAST );
	cairo_rectangle( cairo, 0, 0, _width, _height );
	cairo_clip( cairo );
	cairo_scale( cairo, 50, 50 );
	cairo_translate( cairo, 0, 0 );
	cairo_paint( cairo );

	// Free image surface
	if( _surface )
	    cairo_surface_destroy( _surface );
	_surface = NULL;
    }

    if( _zbuf )
	delete [] _zbuf;
    _zbuf = NULL;
}


void SoftwareRenderer::set_material_diffuse_color( Vec3D color )
{
    _material_diffuse_color = color;
}


void SoftwareRenderer::set_material_ambient_color( Vec3D color )
{
    _material_ambient_color = color;
}


void SoftwareRenderer::set_color( Vec3D color )
{
    _color = color;
}


void SoftwareRenderer::disable_lighting( void )
{

}


void SoftwareRenderer::enable_lighting( void )
{

}


void SoftwareRenderer::flat_triangle( const Vec3D &x0, 
				      const Vec3D &x1, 
				      const Vec3D &x2,
				      const Vec3D &n )
{
    Vec4D y0 = _totalmatrix.transform_homogenous_point( x0 );
    Vec4D y1 = _totalmatrix.transform_homogenous_point( x1 );
    Vec4D y2 = _totalmatrix.transform_homogenous_point( x2 );

    y0.homogenize();
    y1.homogenize();
    y2.homogenize();

    // Backface culling, sign of z-component of cross( y1-y0, y2-y0 )
    if( (y1[0]-y0[0])*(y2[1]-y0[1]) >= (y1[1]-y0[1])*(y2[0]-y0[0]) )
	return;

    Vec3D nn = _normalmatrix.transform_vector( n );
    double dcoef = nn*_light_location;
    Vec3D color;
    for( int a = 0; a < 3; a++ ) {
	color[a] = _light_ambient_color[a]*_material_ambient_color[a] + 
	    dcoef*_light_diffuse_color[a]*_material_diffuse_color[a];
	if( color[a] > 1.0 )
	    color[a] = 1.0;
	else if( color[a] < 0.0 )
	    color[a] = 0.0;
    }

    flat_2d_triangle( 0.5*_width*(y0[0]+1.0), 0.5*_height*(1.0-y0[1]), y0[2],
		      0.5*_width*(y1[0]+1.0), 0.5*_height*(1.0-y1[1]), y1[2],
		      0.5*_width*(y2[0]+1.0), 0.5*_height*(1.0-y2[1]), y2[2],
		      color );
}


void SoftwareRenderer::shaded_triangle( const Vec3D &x0, const Vec3D &c0,
					const Vec3D &x1, const Vec3D &c1,
					const Vec3D &x2, const Vec3D &c2,
					const Vec3D &n )
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


void SoftwareRenderer::line( const Vec3D &x0,
			     const Vec3D &x1 )
{
    Vec4D y0 = _totalmatrix.transform_homogenous_point( x0 );
    Vec4D y1 = _totalmatrix.transform_homogenous_point( x1 );

    y0.homogenize();
    y1.homogenize();

    line_2d( 0.5*_width*(y0[0]+1.0), 0.5*_height*(1.0-y0[1]), y0[2],
	     0.5*_width*(y1[0]+1.0), 0.5*_height*(1.0-y1[1]), y1[2], 
	     _color );
}
