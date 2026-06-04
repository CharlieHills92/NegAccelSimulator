/*! \file vec3d.cpp
 *  \brief Three dimensional vectors.
 */

/* Copyright (c) 2010 Taneli Kalvas. All rights reserved.
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


#include "vec3d.hpp"


Vec3D::Vec3D( const class Vec4D &vec ) 
{ 
    p[0] = vec[0];
    p[1] = vec[1];
    p[2] = vec[2];
}


bool Vec3D::operator!=( const Vec3D &x ) const 
{ 
    for( int a = 0; a < 3; a++ ) {
	if( p[a] != x.p[a] )
	    return( true );
    }
    return( false );
}


bool Vec3D::operator==( const Vec3D &x ) const
{
   for( int a = 0; a < 3; a++ ) {
       if( p[a] != x.p[a] )
	   return( false );
   }

   return( true );
}


double Vec3D::max( void ) const
{
    if( p[0] > p[1] ) {
	if( p[0] > p[2] )
	    return( p[0] );
	else
	    return( p[2] );
    } else {
	if( p[1] > p[2] )
	    return( p[1] );
	else
	    return( p[2] );
    }
}


int32_t Int3D::max( void ) const
{
    if( l[0] > l[1] ) {
	if( l[0] > l[2] )
	    return( l[0] );
	else
	    return( l[2] );
    } else {
	if( l[1] > l[2] )
	    return( l[1] );
	else
	    return( l[2] );
    }
}


bool Vec3D::approx( const Vec3D &x, double eps ) const
{
    for( int a = 0; a < 3; a++ ) {
        if( fabs( p[a] - x.p[a] ) > 1.0e-6 &&
            fabs( (p[a] - x.p[a]) / p[a] ) > 1.0e-6 ) {
            return( false );
        }
    }
    return( true );
}


int Vec3D::min_element( void ) const
{
    int imin = 0;
    double min = fabs( p[0] );
    for( int b = 1; b < 3; b++ ) {
	if( fabs( p[b] ) < min ) {
	    min = fabs( p[b] );
	    imin = b;
	}
    }
    return( imin );
}


Vec3D Vec3D::standard_basis( int i )
{
    Vec3D x;
    x[i] = 1.0;
    return( x );
}


Vec3D Vec3D::arb_perpendicular( void ) const
{
    int i = min_element();
    Vec3D x = Vec3D::standard_basis(i);
    return( cross( *this, x ) );
}


