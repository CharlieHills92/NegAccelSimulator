/*! \file glrenderer.hpp
 *  \brief OpenGL 3D renderer
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


#include <GL/gl.h>
#include "glrenderer.hpp"
#include "ibsimu.hpp"


GLRenderer::GLRenderer( GtkWidget *darea )
    : _darea(darea),
      _material_diffuse_color(0.8,0.0,0.0), 
      _material_ambient_color(0.2,0.0,0.0),
      _color(0.0,0.0,0.0)
{
    GdkGLConfigMode mode = (GdkGLConfigMode)( GDK_GL_MODE_RGBA |
					      GDK_GL_MODE_DEPTH |
					      GDK_GL_MODE_DOUBLE );
    GdkGLConfig *gl_config = gdk_gl_config_new_by_mode( mode );
    if( !gl_config )
        throw( ErrorGLInit() );

    if( !gtk_widget_set_gl_capability( darea, gl_config, NULL, TRUE,
				       GDK_GL_RGBA_TYPE ) )
        throw( ErrorGLInit() );

    ibsimu.message( 1 ) << "Using GLRenderer\n";
    ibsimu.flush();
}


GLRenderer::~GLRenderer()
{

}


void GLRenderer::start_rendering( void )
{
    // Initialize OpenGL context
    _glcontext = gtk_widget_get_gl_context( _darea );
    _gldrawable = gtk_widget_get_gl_drawable( _darea );
    if( !gdk_gl_drawable_gl_begin( _gldrawable, _glcontext ) )
	g_assert_not_reached();

    // Set OpenGL viewport
    GtkAllocation alloc;
    gtk_widget_get_allocation( _darea, &alloc );
    glViewport( 0, 0, alloc.width, alloc.height );

    // Clear
    glClearColor( 1.0, 1.0, 1.0, 0.0 );
    glClearDepth( 1.0 );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glEnable( GL_LIGHTING );
    glEnable( GL_LIGHT0 );
    glEnable( GL_DEPTH_TEST );
}


void GLRenderer::end_rendering( cairo_t *cairo )
{
    // Finish draw and close OpenGL context
    if( gdk_gl_drawable_is_double_buffered( _gldrawable) )
        gdk_gl_drawable_swap_buffers( _gldrawable ); 
    else
        glFlush();
    gdk_gl_drawable_gl_end( _gldrawable );
}


void GLRenderer::enable_view_settings( void )
{
    glMatrixMode( GL_PROJECTION );
    Transformation T = _projection.transpose();
    glLoadMatrixd( &T[0] );

    glMatrixMode( GL_MODELVIEW );
    T = _view.transpose();
    glLoadMatrixd( &T[0] );

    float l[4] = { (float)_light_location[0], 
		   (float)_light_location[1], 
		   (float)_light_location[2], 
		   0.0 };
    glLightfv( GL_LIGHT0, GL_POSITION, l );
    float c1[4] = { (float)_light_ambient_color[0], 
		    (float)_light_ambient_color[1], 
		    (float)_light_ambient_color[2], 
		    1.0 };
    glLightfv( GL_LIGHT0, GL_AMBIENT, c1 );
    float c2[4] = { (float)_light_diffuse_color[0], 
		    (float)_light_diffuse_color[1], 
		    (float)_light_diffuse_color[2], 
		    1.0 };
    glLightfv( GL_LIGHT0, GL_DIFFUSE, c2 );

    float c3[4] = { (float)_material_ambient_color[0], 
		    (float)_material_ambient_color[1], 
		    (float)_material_ambient_color[2], 
		    1.0 };
    glMaterialfv( GL_FRONT, GL_AMBIENT, c3 );
    float c4[4] = { (float)_material_diffuse_color[0], 
		    (float)_material_diffuse_color[1], 
		    (float)_material_diffuse_color[2], 
		    1.0 };
    glMaterialfv( GL_FRONT, GL_DIFFUSE, c4 );

    T = _model.transpose();
    glMultMatrixd( &T[0] );
}


void GLRenderer::set_material_diffuse_color( Vec3D color )
{
    _material_diffuse_color = color;
    float c4[4] = { (float)_material_diffuse_color[0],
		    (float)_material_diffuse_color[1],
		    (float)_material_diffuse_color[2],
		    1.0 };
    glMaterialfv( GL_FRONT, GL_DIFFUSE, c4 );
}


void GLRenderer::set_material_ambient_color( Vec3D color )
{
    _material_ambient_color = color;
    float c3[4] = { (float)_material_ambient_color[0],
		    (float)_material_ambient_color[1],
		    (float)_material_ambient_color[2],
		    1.0 };
    glMaterialfv( GL_FRONT, GL_AMBIENT, c3 );
}


void GLRenderer::set_color( Vec3D color )
{
    _color = color;
    glColor3dv( &_color[0] );
}


void GLRenderer::disable_lighting( void )
{
    glDisable( GL_LIGHTING );
    glColor3dv( &_color[0] );
}


void GLRenderer::enable_lighting( void )
{
    glEnable( GL_LIGHTING );
}


void GLRenderer::flat_triangle( const Vec3D &x0, 
				const Vec3D &x1, 
				const Vec3D &x2,
				const Vec3D &n )
{
    glBegin( GL_TRIANGLES );
    glNormal3dv( &n[0] );
    glVertex3dv( &x0[0] );
    glVertex3dv( &x1[0] );
    glVertex3dv( &x2[0] );
    glEnd();
}


void GLRenderer::shaded_triangle( const Vec3D &x0, const Vec3D &c0,
				  const Vec3D &x1, const Vec3D &c1,
				  const Vec3D &x2, const Vec3D &c2,
				  const Vec3D &n )
{
    glBegin( GL_TRIANGLES );
    glNormal3dv( &n[0] );

    float a0[4] = { (float)(0.2*c0[0]),
		    (float)(0.2*c0[1]),
		    (float)(0.2*c0[2]),
		    1.0 };
    glMaterialfv( GL_FRONT, GL_AMBIENT, a0 );
    float d0[4] = { (float)(0.8*c0[0]),
		    (float)(0.8*c0[1]),
		    (float)(0.8*c0[2]),
		    1.0 };
    glMaterialfv( GL_FRONT, GL_DIFFUSE, d0 );
    glVertex3dv( &x0[0] );

    float a1[4] = { (float)(0.2*c1[0]),
		    (float)(0.2*c1[1]),
		    (float)(0.2*c1[2]),
		    1.0 };
    glMaterialfv( GL_FRONT, GL_AMBIENT, a1 );
    float d1[4] = { (float)(0.8*c1[0]),
		    (float)(0.8*c1[1]),
		    (float)(0.8*c1[2]),
		    1.0 };
    glMaterialfv( GL_FRONT, GL_DIFFUSE, d1 );
    glVertex3dv( &x1[0] );

    float a2[4] = { (float)(0.2*c2[0]),
		    (float)(0.2*c2[1]),
		    (float)(0.2*c2[2]),
		    1.0 };
    glMaterialfv( GL_FRONT, GL_AMBIENT, a2 );
    float d2[4] = { (float)(0.8*c2[0]),
		    (float)(0.8*c2[1]),
		    (float)(0.8*c2[2]),
		    1.0 };
    glMaterialfv( GL_FRONT, GL_DIFFUSE, d2 );
    glVertex3dv( &x2[0] );

    glEnd();
}


void GLRenderer::line( const Vec3D &x0,
		       const Vec3D &x1 )
{
    glBegin( GL_LINES );
    glVertex3dv( &x0[0] );
    glVertex3dv( &x1[0] );
    glEnd();
}
