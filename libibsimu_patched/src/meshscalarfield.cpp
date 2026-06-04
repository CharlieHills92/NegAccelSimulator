/*! \file meshscalarfield.cpp
 *  \brief Mesh based scalar field.
 */

/* Copyright (c) 2005-2012 Taneli Kalvas. All rights reserved.
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

#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>
#include <limits>
#include "meshscalarfield.hpp"
#include "ibsimu.hpp"


MeshScalarField::MeshScalarField()
    : _F(NULL)
{
}


MeshScalarField::MeshScalarField( const Mesh &m )
    : Mesh(m)
{
    check_definition();

    _F = new double[_size[0]*_size[1]*_size[2]];
    memset( _F, 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
}


MeshScalarField::MeshScalarField( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h )
    : Mesh(geom_mode,size,origo,h)
{
    check_definition();

    _F = new double[_size[0]*_size[1]*_size[2]];
    memset( _F, 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
}


MeshScalarField::MeshScalarField( std::istream &s )
    : Mesh(s)
{
    check_definition();

    ibsimu.message( 1 ) << "Constructing ScalarField from stream\n";

    _F = new double[_size[0]*_size[1]*_size[2]];
    read_compressed_block( s, _size[0]*_size[1]*_size[2]*sizeof(double), (int8_t *)_F );
}


MeshScalarField::MeshScalarField( const MeshScalarField &f )
    : Mesh(f)
{
    _F = new double[_size[0]*_size[1]*_size[2]];
    memcpy( _F, f._F, _size[0]*_size[1]*_size[2]*sizeof(double) );
}


MeshScalarField::~MeshScalarField()
{
    if( _F )
	delete [] _F;
}


void MeshScalarField::check_definition()
{
    // Check that mesh size need for MeshScalarField are fullfilled
    if( _geom_mode == MODE_3D ) {
	if( _size[0] < 2 || _size[1] < 2 || _size[2] < 2 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    } else if( _geom_mode == MODE_2D || _geom_mode == MODE_CYL ) {
	if( _size[0] < 2 || _size[1] < 2 || _size[2] != 1 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    } else {
	if( _size[0] < 2 || _size[1] != 1 || _size[2] != 1 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    }
}


void MeshScalarField::clear()
{
    if( _F )
	memset( _F, 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
}


void MeshScalarField::reset( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h )
{
    Mesh::reset( geom_mode, size, origo, h );
    check_definition();

    if( _F )
	delete [] _F;
    _F = new double[_size[0]*_size[1]*_size[2]];
    memset( _F, 0, _size[0]*_size[1]*_size[2]*sizeof(double) );
}


void MeshScalarField::get_minmax( double &min, double &max ) const
{
    if( !_F ) {
	min = max = 0.0;
	return;
    }

    min = max = _F[0];

    /*
    size_t imin = 0, jmin = 0, kmin = 0;
    size_t imax = 0, jmax = 0, kmax = 0;
    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t k = 0; k < _size[2]; k++ ) {
	for( size_t j = 0; j < _size[1]; j++ ) {
	    for( size_t i = 0; i < _size[0]; i++ ) {
		double F = _F[i+(j+k*_size[1])*_size[0]];
		if( F < min ) {
		    imin = i;
		    jmin = j;
		    kmin = k;
		    min = F;
		} if( F > max ) {
		    imax = i;
		    jmax = j;
		    kmax = k;
		    max = F;
		}
	    }
	}
    }

    std::cout << "min  = " << min << "\n";
    std::cout << "imin = " << imin << "\n";
    std::cout << "jmin = " << jmin << "\n";
    std::cout << "kmin = " << kmin << "\n";

    std::cout << "max  = " << max << "\n";
    std::cout << "imax = " << imax << "\n";
    std::cout << "jmax = " << jmax << "\n";
    std::cout << "kmax = " << kmax << "\n";
    */

    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t a = 1; a < ncount; a++ ) {
	if( _F[a] < min )
	    min = _F[a];
	if( _F[a] > max )
	    max = _F[a];
    }
}


MeshScalarField &MeshScalarField::operator=( const MeshScalarField &f )
{
    (Mesh)(*this) = (Mesh)f;

    if( _F )
	delete [] _F;
    _F = new double[_size[0]*_size[1]*_size[2]];
    memcpy( _F, f._F, _size[0]*_size[1]*_size[2]*sizeof(double) );
    return( *this );
}


MeshScalarField &MeshScalarField::operator+=( const MeshScalarField &f )
{
    if( (Mesh)(*this) != (Mesh)f )
	throw( Error( ERROR_LOCATION, "non-matching fields" ) );

    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t a = 0; a < ncount; a++ )
	_F[a] += f._F[a];
    return( *this );
}


MeshScalarField &MeshScalarField::operator-=( const MeshScalarField &f )
{
    if( (Mesh)(*this) != (Mesh)f )
	throw( Error( ERROR_LOCATION, "non-matching fields" ) );

    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t a = 0; a < ncount; a++ )
	_F[a] -= f._F[a];
    return( *this );
}


MeshScalarField &MeshScalarField::operator*=( double x )
{
    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t a = 0; a < ncount; a++ )
	_F[a] *= x;
    return( *this );
}


MeshScalarField &MeshScalarField::operator/=( double x )
{
    double xi = 1.0/x;
    size_t ncount = _size[0]*_size[1]*_size[2];
    for( size_t a = 0; a < ncount; a++ )
	_F[a] *= xi;
    return( *this );
}


double MeshScalarField::operator()( const Vec3D &x ) const
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


void MeshScalarField::save( const std::string &filename ) const
{
    ibsimu.message( 1 ) << "Saving MeshScalarField to file \'" << filename << "\'.\n";

    std::ofstream os( filename.c_str(), std::ios_base::binary );
    if( !os.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\' for writing" ) );
    save( os );
    os.close();
}


void MeshScalarField::save( std::ostream &os ) const
{
    Mesh::save( os );
    write_compressed_block( os, _size[0]*_size[1]*_size[2]*sizeof(double), (int8_t *)_F );
}


void MeshScalarField::debug_print( std::ostream &os ) const
{
    Mesh::debug_print( os );

    os << "**MeshScalarField\n";
    os << "F = (";
    if( _size[0]*_size[1]*_size[2] < 10 ) {
	int a;
	for( a = 0; a < _size[0]*_size[1]*_size[2]-1; a++ )
	    os << _F[a] << ", ";
	if( a < _size[0]*_size[1]*_size[2] )
	    os << _F[a] << ")\n";
    } else {
	// Print only 10 first nodes
	for( int a = 0; a < 10; a++ )
	    os << _F[a] << ", ";
	os << "... )\n";
    }

}



