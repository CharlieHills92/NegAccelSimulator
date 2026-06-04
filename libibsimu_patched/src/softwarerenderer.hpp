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


#ifndef SOFTWARERENDERER_HPP
#define SOFTWARERENDERER_HPP 1


#include <cairo.h>
#include <gtk/gtk.h>
#include "renderer.hpp"


/*! \brief Software 3D z-buffer renderer.
 *
 *  Capable of rendering flat shaded triangles and lines in 3D space.
 *
 *  Intended for replacing OpenGL when it is not available and for
 *  making hardcopies.
 */
class SoftwareRenderer : public Renderer {

    GtkWidget       *_darea;
    cairo_surface_t *_surface;
    unsigned char   *_buf;
    int              _width;
    int              _height;
    int              _stride;

    double          *_zbuf;

    Vec3D            _material_diffuse_color;
    Vec3D            _material_ambient_color;
    Vec3D            _color;

    Transformation   _modelview;
    Transformation   _totalmatrix;
    Transformation   _normalmatrix;

    void swap( int &x0, int &y0, double &z0, 
	       int &x1, int &y1, double &z1 );

    void set_pixel( int i, int j, const Vec3D &color );
    void clear( const Vec3D &color );
    void flat_2d_triangle( int x0, int y0, double z0,
			   int x1, int y1, double z1,
			   int x2, int y2, double z2,
			   const Vec3D &color );
    void line_2d( int x0, int y0, double z0,
		  int x1, int y1, double z1,
		  const Vec3D &color );

public:

    /*! \brief Constructor for rendering to drawing area.
     */
    SoftwareRenderer( GtkWidget *darea );

    /*! \brief Constructor for rendering to cairo image surface.
     */
    SoftwareRenderer( cairo_surface_t *surface );

    /*! \brief Destructor.
     */
    virtual ~SoftwareRenderer();

    virtual void start_rendering( void );
    virtual void end_rendering( cairo_t *cairo );

    virtual void set_material_diffuse_color( Vec3D color );
    virtual void set_material_ambient_color( Vec3D color );
    virtual void set_color( Vec3D color );

    virtual void disable_lighting( void );
    virtual void enable_lighting( void );

    virtual void enable_view_settings( void );

    virtual void flat_triangle( const Vec3D &x0, 
				const Vec3D &x1, 
				const Vec3D &x2,
				const Vec3D &n );
    virtual void shaded_triangle( const Vec3D &x0, const Vec3D &c0,
				  const Vec3D &x1, const Vec3D &c1,
				  const Vec3D &x2, const Vec3D &c2,
				  const Vec3D &n );
    virtual void line( const Vec3D &x0,
		       const Vec3D &x1 );
};


#endif
