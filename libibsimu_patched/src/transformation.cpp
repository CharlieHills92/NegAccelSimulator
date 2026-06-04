/*! \file transformation.cpp
 *  \brief Affine transformation
 */

/* Copyright (c) 2010-2012,2015 Taneli Kalvas. All rights reserved.
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
#include <iomanip>
#include "transformation.hpp"


Transformation::Transformation() 
{
    x[0] = x[5] = x[10] = x[15] = 1.0; 

    x[1] = x[2] = x[3] = x[4] = x[6] = x[7] = x[8] = x[9] 
	= x[11] = x[12] = x[13] = x[14] = 0.0;
}


Transformation::Transformation( double x11, double x12, double x13, double x14,
				double x21, double x22, double x23, double x24,
				double x31, double x32, double x33, double x34,
				double x41, double x42, double x43, double x44 ) 
{ 
    x[0]  = x11;
    x[1]  = x12;
    x[2]  = x13;
    x[3]  = x14;
    x[4]  = x21;
    x[5]  = x22;
    x[6]  = x23;
    x[7]  = x24;
    x[8]  = x31;
    x[9]  = x32;
    x[10] = x33;
    x[11] = x34;
    x[12] = x41;
    x[13] = x42;
    x[14] = x43;
    x[15] = x44;
}


Transformation::Transformation( const Transformation &m )
{ 
    memcpy( x, m.x, 16*sizeof(double) );
}


Transformation::Transformation( std::istream &is )
{
    for( int i = 0; i < 16; i++ )
	x[i] = read_double( is );
}


Transformation::~Transformation()
{

}


double Transformation::determinant( void ) const
{
    return( + x[3] * x[6] * x[9] * x[12]
            - x[2] * x[7] * x[9] * x[12]
            - x[3] * x[5] * x[10] * x[12]
            + x[1] * x[7] * x[10] * x[12]
            + x[2] * x[5] * x[11] * x[12]
            - x[1] * x[6] * x[11] * x[12]
            - x[3] * x[6] * x[8] * x[13]
            + x[2] * x[7] * x[8] * x[13]
            + x[3] * x[4] * x[10] * x[13]
            - x[0] * x[7] * x[10] * x[13]
            - x[2] * x[4] * x[11] * x[13]
            + x[0] * x[6] * x[11] * x[13]
            + x[3] * x[5] * x[8] * x[14]
            - x[1] * x[7] * x[8] * x[14]
            - x[3] * x[4] * x[9] * x[14]
            + x[0] * x[7] * x[9] * x[14]
            + x[1] * x[4] * x[11] * x[14]
            - x[0] * x[5] * x[11] * x[14]
            - x[2] * x[5] * x[8] * x[15]
            + x[1] * x[6] * x[8] * x[15]
            + x[2] * x[4] * x[9] * x[15]
            - x[0] * x[6] * x[9] * x[15]
            - x[1] * x[4] * x[10] * x[15]
            + x[0] * x[5] * x[10] * x[15] );
}


Transformation Transformation::inverse( void ) const
{
    double idet = determinant();
    if( idet == 0.0 )
	throw( Error( ERROR_LOCATION, "can't invert matrix: zero determinant" ) );
    idet = 1.0/idet;

    Transformation ret( x[6]*x[11]*x[13] - x[7]*x[10]*x[13] + x[7]*x[9]*x[14] - 
			x[5]*x[11]*x[14] - x[6]*x[9]*x[15] + x[5]*x[10]*x[15],
			x[3]*x[10]*x[13] - x[2]*x[11]*x[13] - x[3]*x[9]*x[14] + 
			x[1]*x[11]*x[14] + x[2]*x[9]*x[15] - x[1]*x[10]*x[15],
			x[2]*x[7]*x[13] - x[3]*x[6]*x[13] + x[3]*x[5]*x[14] - 
			x[1]*x[7]*x[14] - x[2]*x[5]*x[15] + x[1]*x[6]*x[15],
			x[3]*x[6]*x[9] - x[2]*x[7]*x[9] - x[3]*x[5]*x[10] + 
			x[1]*x[7]*x[10] + x[2]*x[5]*x[11] - x[1]*x[6]*x[11],
			x[7]*x[10]*x[12] - x[6]*x[11]*x[12] - x[7]*x[8]*x[14] + 
			x[4]*x[11]*x[14] + x[6]*x[8]*x[15] - x[4]*x[10]*x[15],
			x[2]*x[11]*x[12] - x[3]*x[10]*x[12] + x[3]*x[8]*x[14] - 
			x[0]*x[11]*x[14] - x[2]*x[8]*x[15] + x[0]*x[10]*x[15],
			x[3]*x[6]*x[12] - x[2]*x[7]*x[12] - x[3]*x[4]*x[14] + 
			x[0]*x[7]*x[14] + x[2]*x[4]*x[15] - x[0]*x[6]*x[15],
			x[2]*x[7]*x[8] - x[3]*x[6]*x[8] + x[3]*x[4]*x[10] - 
			x[0]*x[7]*x[10] - x[2]*x[4]*x[11] + x[0]*x[6]*x[11],
			x[5]*x[11]*x[12] - x[7]*x[9]*x[12] + x[7]*x[8]*x[13] - 
			x[4]*x[11]*x[13] - x[5]*x[8]*x[15] + x[4]*x[9]*x[15],
			x[3]*x[9]*x[12] - x[1]*x[11]*x[12] - x[3]*x[8]*x[13] + 
			x[0]*x[11]*x[13] + x[1]*x[8]*x[15] - x[0]*x[9]*x[15],
			x[1]*x[7]*x[12] - x[3]*x[5]*x[12] + x[3]*x[4]*x[13] - 
			x[0]*x[7]*x[13] - x[1]*x[4]*x[15] + x[0]*x[5]*x[15],
			x[3]*x[5]*x[8] - x[1]*x[7]*x[8] - x[3]*x[4]*x[9] + 
			x[0]*x[7]*x[9] + x[1]*x[4]*x[11] - x[0]*x[5]*x[11],
			x[6]*x[9]*x[12] - x[5]*x[10]*x[12] - x[6]*x[8]*x[13] + 
			x[4]*x[10]*x[13] + x[5]*x[8]*x[14] - x[4]*x[9]*x[14],
			x[1]*x[10]*x[12] - x[2]*x[9]*x[12] + x[2]*x[8]*x[13] - 
			x[0]*x[10]*x[13] - x[1]*x[8]*x[14] + x[0]*x[9]*x[14],
			x[2]*x[5]*x[12] - x[1]*x[6]*x[12] - x[2]*x[4]*x[13] + 
			x[0]*x[6]*x[13] + x[1]*x[4]*x[14] - x[0]*x[5]*x[14],
			x[1]*x[6]*x[8] - x[2]*x[5]*x[8] + x[2]*x[4]*x[9] - 
			x[0]*x[6]*x[9] - x[1]*x[4]*x[10] + x[0]*x[5]*x[10] );
    ret *= idet;
    return( ret );
}


Transformation Transformation::transpose( void ) const
{
    Transformation ret( x[0], x[4], x[8], x[12],
			x[1], x[5], x[9], x[13],
			x[2], x[6], x[10], x[14],
			x[3], x[7], x[11], x[15] );
    return( ret );
}


bool Transformation::operator!=( const Transformation &m ) const
{
    for( size_t i = 0; i < 16; i++ ) {
        if( x[i] != m[i] )
	    return( true );
    }
    return( false );
}


const Transformation &Transformation::operator*=( double s )
{
    for( size_t i = 0; i < 16; i++ )
        x[i] *= s;
    return( *this );
}


Transformation Transformation::operator*( const Transformation &m ) const
{
    Transformation r;

    r[ 0] = x[ 0]*m[ 0] + x[ 1]*m[ 4] + x[ 2]*m[ 8] + x[ 3]*m[12];
    r[ 1] = x[ 0]*m[ 1] + x[ 1]*m[ 5] + x[ 2]*m[ 9] + x[ 3]*m[13];
    r[ 2] = x[ 0]*m[ 2] + x[ 1]*m[ 6] + x[ 2]*m[10] + x[ 3]*m[14];
    r[ 3] = x[ 0]*m[ 3] + x[ 1]*m[ 7] + x[ 2]*m[11] + x[ 3]*m[15];

    r[ 4] = x[ 4]*m[ 0] + x[ 5]*m[ 4] + x[ 6]*m[ 8] + x[ 7]*m[12];
    r[ 5] = x[ 4]*m[ 1] + x[ 5]*m[ 5] + x[ 6]*m[ 9] + x[ 7]*m[13];
    r[ 6] = x[ 4]*m[ 2] + x[ 5]*m[ 6] + x[ 6]*m[10] + x[ 7]*m[14];
    r[ 7] = x[ 4]*m[ 3] + x[ 5]*m[ 7] + x[ 6]*m[11] + x[ 7]*m[15];

    r[ 8] = x[ 8]*m[ 0] + x[ 9]*m[ 4] + x[10]*m[ 8] + x[11]*m[12];
    r[ 9] = x[ 8]*m[ 1] + x[ 9]*m[ 5] + x[10]*m[ 9] + x[11]*m[13];
    r[10] = x[ 8]*m[ 2] + x[ 9]*m[ 6] + x[10]*m[10] + x[11]*m[14];
    r[11] = x[ 8]*m[ 3] + x[ 9]*m[ 7] + x[10]*m[11] + x[11]*m[15];

    r[12] = x[12]*m[ 0] + x[13]*m[ 4] + x[14]*m[ 8] + x[15]*m[12];
    r[13] = x[12]*m[ 1] + x[13]*m[ 5] + x[14]*m[ 9] + x[15]*m[13];
    r[14] = x[12]*m[ 2] + x[13]*m[ 6] + x[14]*m[10] + x[15]*m[14];
    r[15] = x[12]*m[ 3] + x[13]*m[ 7] + x[14]*m[11] + x[15]*m[15];

    return( r );
}


Vec4D Transformation::operator*( const Vec4D &v ) const
{
    Vec4D r;

    r[0] =  x[0]*v[0] +  x[1]*v[1] +  x[2]*v[2] +  x[3]*v[3];
    r[1] =  x[4]*v[0] +  x[5]*v[1] +  x[6]*v[2] +  x[7]*v[3];
    r[2] =  x[8]*v[0] +  x[9]*v[1] + x[10]*v[2] + x[11]*v[3];
    r[3] = x[12]*v[0] + x[13]*v[1] + x[14]*v[2] + x[15]*v[3];

    return( r );
}


Vec4D Transformation::operator%( const Vec4D &v ) const
{
    Vec4D r;

    r[0] =  x[0]*v[0] +  x[4]*v[1] +  x[8]*v[2] + x[12]*v[3];
    r[1] =  x[1]*v[0] +  x[5]*v[1] +  x[9]*v[2] + x[13]*v[3];
    r[2] =  x[2]*v[0] +  x[6]*v[1] + x[10]*v[2] + x[14]*v[3];
    r[3] =  x[3]*v[0] +  x[7]*v[1] + x[11]*v[2] + x[15]*v[3];

    return( r );
}


Vec3D Transformation::inv_transform_point( const Vec3D &xin ) const
{
    Transformation t = this->inverse();
    Vec4D r = t * Vec4D( xin[0], xin[1], xin[2], 1.0 );

    return( Vec3D( r[0], r[1], r[2] ) );
}


Vec3D Transformation::inv_transform_vector( const Vec3D &xin ) const
{
    Transformation t = this->inverse();
    Vec4D r = t * Vec4D( xin[0], xin[1], xin[2], 0.0 );

    return( Vec3D( r[0], r[1], r[2] ) );
}


void Transformation::reset( void ) 
{
    Transformation t1 = unity();
    *this = t1;
}


void Transformation::translate( const Vec3D &d ) 
{
    *this = translation( d ) * (*this);
}


void Transformation::translate_before( const Vec3D &d ) 
{
    *this = (*this) * translation( d );
}


void Transformation::scale( const Vec3D &s ) 
{
    *this = scaling( s ) * (*this);
}


void Transformation::scale_before( const Vec3D &s ) 
{
    *this = (*this) * scaling( s );
}


void Transformation::rotate_x( double a ) 
{
    *this = rotation_x( a ) * (*this);
}


void Transformation::rotate_x_before( double a ) 
{
    *this = (*this) * rotation_x( a );
}


void Transformation::rotate_y( double a ) 
{
    *this = rotation_y( a ) * (*this);
}


void Transformation::rotate_y_before( double a ) 
{
    *this = (*this) * rotation_y( a );
}


void Transformation::rotate_z( double a ) 
{
    *this = rotation_z( a ) * (*this);
}


void Transformation::rotate_z_before( double a ) 
{
    *this = (*this) * rotation_z( a );
}


Transformation Transformation::unity( void ) 
{
    return( Transformation(  1,  0,  0,  0,
			     0,  1,  0,  0,
			     0,  0,  1,  0,
			     0,  0,  0,  1 ) );
}


Transformation Transformation::translation( const Vec3D &d ) 
{
    return( Transformation(  1,  0,  0, d[0],
			     0,  1,  0, d[1],
			     0,  0,  1, d[2],
			     0,  0,  0,    1 ) );
}


Transformation Transformation::scaling( const Vec3D &s ) 
{
    return( Transformation( s[0],    0,    0,  0,
			    0, s[1],    0,  0,
			    0,    0, s[2],  0,
			    0,    0,    0,  1 ) );
}
 

Transformation Transformation::rotation_x( double a ) 
{
    return( Transformation(  1,      0,       0,  0,
			     0, cos(a), -sin(a),  0,
			     0, sin(a),  cos(a),  0,
			     0,      0,       0,  1 ) );
}


Transformation Transformation::rotation_y( double a ) 
{
    return( Transformation(  cos(a),  0, sin(a),  0,
			     0,  1,      0,  0,
			     -sin(a),  0, cos(a),  0,
			     0,  0,      0,  1 ) );
}


Transformation Transformation::rotation_z( double a ) 
{
    return( Transformation( cos(a), -sin(a),  0,  0,
			    sin(a),  cos(a),  0,  0,
			    0,       0,  1,  0,
			    0,       0,  0,  1 ) );
}


void Transformation::save( const std::string &filename ) const
{
    std::ofstream os( filename.c_str(), std::ios_base::binary );
    if( !os.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\' for writing" ) );
    save( os );
    os.close();
}


void Transformation::save( std::ostream &os ) const
{
    for( int i = 0; i < 16; i++ )
	write_double( os, x[i] );
}


std::ostream &operator<<( std::ostream &os, const Transformation &t ) 
{
    for( size_t i = 0; i < 4; i++ ) {
	for( size_t j = 0; j < 4; j++ ) {
	    os << std::setw(12) << to_string(t.x[4*i+j]).substr(0,12) << " ";
	}
	os << "\n";
    }

    return( os );
}


void Transformation::debug_print( std::ostream &os ) const
{
    std::cout << "**Transformation\n";

    os << "x = \n";
    os << *this << "\n";
}

