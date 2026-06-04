/*! \file mesh.cpp
 *  \brief Rectangular mesh definition.
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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


#include "mesh.hpp"
#include "error.hpp"
#include "ibsimu.hpp"


Mesh::Mesh()
    : _geom_mode(MODE_3D), _h(1.0)
{
    _div_h = 1.0/_h;
}


Mesh::Mesh( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h )
    : _geom_mode(geom_mode), _size(size), _origo(origo)
{
    _h = fabs(h);
    _div_h = 1.0/_h;
    _max = Vec3D( _origo(0)+_h*(_size[0]-1),
		  _origo(1)+_h*(_size[1]-1),
		  _origo(2)+_h*(_size[2]-1) );

    // Checks
    if( _h == 0.0 )
	throw( Error( ERROR_LOCATION, "zero mesh step size" ) );
    else if( _geom_mode == MODE_CYL && _origo[1] < 0.0 )
	throw( Error( ERROR_LOCATION, "negative origo in r-direction" ) );
    else if( _size[0] <= 0 || _size[1] <= 0 || _size[2] <= 0 )
    	throw( Error( ERROR_LOCATION, "zero mesh size" ) );
}


void Mesh::reset( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h )
{
    _geom_mode = geom_mode;
    _size = size;
    _origo = origo;
    _h = fabs(h);
    _div_h = 1.0/_h;
    _max = Vec3D( _origo(0)+_h*(_size(0)-1),
		  _origo(1)+_h*(_size(1)-1),
		  _origo(2)+_h*(_size(2)-1) );

    // Checks
    if( _h == 0.0 )
	throw( Error( ERROR_LOCATION, "zero mesh step size" ) );
    else if( _geom_mode == MODE_CYL && _origo[1] < 0.0 )
	throw( Error( ERROR_LOCATION, "negative origo in r-direction" ) );
    else if( _size[0] <= 0 || _size[1] <= 0 || _size[2] <= 0 )
    	throw( Error( ERROR_LOCATION, "zero mesh size" ) );
}


Mesh::Mesh( std::istream &is )
{
    _geom_mode = (geom_mode_e)read_int32( is );
    _size      = Int3D( is );
    _origo     = Vec3D( is );
    _h         = read_double( is );

    // Calculate vector max and reciprocal of h
    _div_h = 1.0/_h;
    _max = Vec3D( _origo(0)+_h*(_size[0]-1),
		  _origo(1)+_h*(_size[1]-1),
		  _origo(2)+_h*(_size[2]-1) );
}


uint32_t Mesh::dim( void ) const
{
    switch( _geom_mode ) {
    case MODE_1D:
	return( 1 );
	break;
    case MODE_2D:
	return( 2 );
	break;
    case MODE_CYL:
	return( 2 );
	break;
    default:
	return( 3 );
	break;
    }
}


Int3D Mesh::closest_node( Vec3D x ) const
{
    return( Int3D( floor((x[0]-_origo(0))/_h + 0.5),
		   floor((x[1]-_origo(1))/_h + 0.5),
		   floor((x[2]-_origo(2))/_h + 0.5) ) );
}


Int3D Mesh::mesh_number( Vec3D x ) const
{
    return( Int3D( floor((x[0]-_origo(0))/_h),
		   floor((x[1]-_origo(1))/_h),
		   floor((x[2]-_origo(2))/_h) ) );
}


Vec3D Mesh::coord_of_node( Int3D n ) const
{
    return( Vec3D( _origo(0)+n[0]*_h,
		   _origo(1)+n[1]*_h,
		   _origo(2)+n[2]*_h ) );
}


void Mesh::save( std::ostream &os ) const
{
    write_int32( os, _geom_mode );
    _size.save( os );
    _origo.save( os );
    write_double( os, _h );
}


bool Mesh::operator==( const Mesh &m ) const
{
    if( _geom_mode == m._geom_mode && 
	_size == m._size && 
	_origo.approx( m._origo, 1.0e-6 ) && 
	fabs( (_h - m._h) / _h ) < 1.0e-6 )
	return( true );
    return( false );
}


bool Mesh::operator!=( const Mesh &m ) const
{
    if( _geom_mode != m._geom_mode ||
	_size != m._size || 
	!_origo.approx( m._origo, 1.0e-6 ) || 
	_h != m._h )
	return( true );
    return( false );
}


void Mesh::debug_print( std::ostream &os ) const
{
    os << "**Mesh\n";

    switch( _geom_mode ) {
    case MODE_1D:
	os << "geom_mode = MODE_1D\n";
	break;
    case MODE_2D:
	os << "geom_mode = MODE_2D\n";
	break;
    case MODE_CYL:
	os << "geom_mode = MODE_CYL\n";
	break;
    case MODE_3D:
	os << "geom_mode = MODE_3D\n";
	break;
    }
    os << "size = (" 
       << _size[0] << ", "
       << _size[1] << ", "
       << _size[2] << ")\n";
    os << "origo = (" 
       << _origo[0] << ", "
       << _origo[1] << ", "
       << _origo[2] << ")\n";
    os << "max = (" 
       << _max[0] << ", "
       << _max[1] << ", "
       << _max[2] << ")\n";
    os << "h = " << _h << "\n";
}



