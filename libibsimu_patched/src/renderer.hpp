/*! \file renderer.hpp
 *  \brief 3D renderer base class
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


#ifndef RENDERER_HPP
#define RENDERER_HPP 1


#include <cairo.h>
#include "transformation.hpp"


/*! \brief 3D renderer base class.
 *
 *  Front-end for 3D renderer classes used in 3D geometry plotters.
 */
class Renderer {

protected:

    Vec3D            _light_diffuse_color;
    Vec3D            _light_location;
    Vec3D            _light_ambient_color;

    Transformation   _model;
    Transformation   _view;
    Transformation   _projection;

public:

    Renderer();
    virtual ~Renderer();


    void set_light_diffuse_color( Vec3D color );
    void set_light_location( Vec3D location );
    void set_light_ambient_color( Vec3D color );

    virtual void set_material_diffuse_color( Vec3D color ) = 0;
    virtual void set_material_ambient_color( Vec3D color ) = 0;
    virtual void set_color( Vec3D color ) = 0;

    void set_model_transformation( const Transformation &model );
    void set_projection_frustum( double left, double right,
				 double bottom, double top,
				 double near, double far );
    void set_view_look_at( const Vec3D &camera, 
			   const Vec3D &target,
			   const Vec3D &up );


    virtual void start_rendering( void ) = 0;
    virtual void end_rendering( cairo_t *cairo ) = 0;


    virtual void disable_lighting( void ) = 0;
    virtual void enable_lighting( void ) = 0;


    virtual void enable_view_settings( void ) = 0;


    virtual void flat_triangle( const Vec3D &x0,
				const Vec3D &x1,
				const Vec3D &x2,
				const Vec3D &n ) = 0;
    virtual void shaded_triangle( const Vec3D &x0, const Vec3D &c0,
				  const Vec3D &x1, const Vec3D &c1,
				  const Vec3D &x2, const Vec3D &c2,
				  const Vec3D &n ) = 0;
    virtual void line( const Vec3D &x0,
		       const Vec3D &x1 ) = 0;
};


#endif
