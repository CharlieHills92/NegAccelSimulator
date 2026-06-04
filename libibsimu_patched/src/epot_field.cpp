/*! \file epot_field.hpp
 *  \brief Electric potential field.
 */

/* Copyright (c) 2011 Taneli Kalvas. All rights reserved.
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


#include "epot_field.hpp"


EpotField::EpotField( const Geometry &geom )
    : MeshScalarField(geom), _geom(&geom)
{
}


EpotField::EpotField( const EpotField &f )
    : MeshScalarField(f), _geom(f._geom)
{
}


EpotField::EpotField( std::istream &s, const Geometry &geom )
    : MeshScalarField(s), _geom(&geom)
{
    // Check geometry compatibility
    if( _geom_mode != geom.geom_mode() ||
	_size != geom.size() || 
	_origo != geom.origo() || 
	_h != geom.h() ) {
	throw( Error( ERROR_LOCATION, "incompatible Geometry object" ) );
    }
}


EpotField::~EpotField()
{
}


const Geometry *EpotField::geom( void ) const
{
    return( _geom );
}


/* Evaluator for electric potential field
 *
 * Try to find virtual potentials at solid nodes to get better
 * interpolation close to solids. Is not possible like it was done
 * before. Epot must be continuous! This is difficult to achieve.
 *
 * TODO. Current implementation uses simple interpolation. No
 * knowledge of solids.
 */
double EpotField::operator()( const Vec3D &x ) const
{
    if( !_F )
	return( 0.0 );

    switch( _geom_mode ) {
    case MODE_1D:
    {
	int32_t i = (int32_t)floor( (x[0]-_origo[0])*_div_h );
	if( i < 0 )
	    i = 0;
	else if( i >= _size[0]-1 )
	    i = _size[0]-2;

	double t = _div_h*( x[0]-(i*_h+_origo[0]) );

	return( (1.0-t)*_F[i] + t*_F[i+1] );
	break;
    }
    case MODE_2D:
    case MODE_CYL:
    {
	int32_t i = (int32_t)floor( (x[0]-_origo[0])*_div_h );
	int32_t j = (int32_t)floor( (x[1]-_origo[1])*_div_h );
	if( i < 0 )
	    i = 0;
	else if( i >= _size[0]-1 )
	    i = _size[0]-2;
	if( j < 0 )
	    j = 0;
	else if( j >= _size[1]-1 )
	    j = _size[1]-2;

	double t = _div_h*( x[0]-(i*_h+_origo[0]) );
	double u = _div_h*( x[1]-(j*_h+_origo[1]) );

	int32_t ptr = _size[0]*j + i;
	return( (1.0-t)*(1.0-u)*_F[ptr] +
		(    t)*(1.0-u)*_F[ptr+1] +
		(1.0-t)*(    u)*_F[ptr+_size[0]] +
		(    t)*(    u)*_F[ptr+1+_size[0]] );
	break;
    }
    default:
    {
	int32_t i = (int32_t)floor( (x[0]-_origo[0])*_div_h );
	int32_t j = (int32_t)floor( (x[1]-_origo[1])*_div_h );
	int32_t k = (int32_t)floor( (x[2]-_origo[2])*_div_h );
	if( i < 0 )
	    i = 0;
	else if( i >= _size[0]-1 )
	    i = _size[0]-2;
	if( j < 0 )
	    j = 0;
	else if( j >= _size[1]-1 )
	    j = _size[1]-2;
	if( k < 0 )
	    k = 0;
	else if( k >= _size[2]-1 )
	    k = _size[2]-2;

	double t = _div_h*( x[0]-(i*_h+_origo[0]) );
	double u = _div_h*( x[1]-(j*_h+_origo[1]) );
	double v = _div_h*( x[2]-(k*_h+_origo[2]) );

	int32_t b  = _size[0]*_size[1];
	int32_t ptr = b*k + _size[0]*j + i;
	return( (1.0-t)*(1.0-u)*(1.0-v)*_F[ptr] +
		(    t)*(1.0-u)*(1.0-v)*_F[ptr+1] +
		(1.0-t)*(    u)*(1.0-v)*_F[ptr+_size[0]] +
		(    t)*(    u)*(1.0-v)*_F[ptr+1+_size[0]] +
		(1.0-t)*(1.0-u)*(    v)*_F[ptr+b] +
		(    t)*(1.0-u)*(    v)*_F[ptr+1+b] +
		(1.0-t)*(    u)*(    v)*_F[ptr+_size[0]+b] +
		(    t)*(    u)*(    v)*_F[ptr+1+_size[0]+b] );
	break;
    }
    }
}
