/*! \file renderer.cpp
 *  \brief 3D renderer base class
 */

/* Copyright (c) 2012 Taneli Kalvas. All rights reserved.
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


#include "renderer.hpp"


Renderer::Renderer() 
 : _light_diffuse_color(1.0,1.0,1.0), 
   _light_location(1.0,1.0,1.0),
   _light_ambient_color(1.0,1.0,1.0)
{

}


Renderer::~Renderer()
{

}


void Renderer::set_light_diffuse_color( Vec3D color )
{
    _light_diffuse_color = color;
}


void Renderer::set_light_location( Vec3D location )
{
    _light_location = location;
}


void Renderer::set_light_ambient_color( Vec3D color )
{
    _light_ambient_color = color;
}


void Renderer::set_model_transformation( const Transformation &model )
{
    _model = model;
}


void Renderer::set_projection_frustum( double left, double right,
				       double bottom, double top,
				       double near, double far )
{
    _projection = Transformation( 2.0*near/(right-left), 0.0, (right+left)/(right-left), 0.0,
				  0.0, 2*near/(top-bottom), (top+bottom)/(top-bottom), 0.0,
				  0.0, 0.0, (far+near)/(near-far), 2.0*far*near/(near-far),
				  0.0, 0.0, -1.0, 0.0 );
}


void Renderer::set_view_look_at( const Vec3D &camera, 
				 const Vec3D &target,
				 const Vec3D &up )
{
    Vec3D forward = target-camera;
    forward.normalize();
    Vec3D side = cross( forward, up );
    side.normalize();
    Vec3D realup = cross( side, forward );

    Transformation rot = Transformation( side[0], side[1], side[2], 0.0,
					 realup[0], realup[1], realup[2], 0.0,
					 -forward[0], -forward[1], -forward[2], 0.0,
					 0.0, 0.0, 0.0, 1.0 );
    Transformation trans = Transformation(  1,  0,  0, -camera[0],
					    0,  1,  0, -camera[1],
					    0,  0,  1, -camera[2],
					    0,  0,  0,    1 );
    _view = rot*trans;
}


