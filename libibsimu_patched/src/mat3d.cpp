/*! \file mat3d.cpp
 *  \brief Three-by-three matrices.
 */

/* Copyright (c) 2010-2012 Taneli Kalvas. All rights reserved.
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


#include <iomanip>
#include "mat3d.hpp"
#include "error.hpp"


Mat3D::Mat3D()
{
    a[0] = a[1] = a[2] = a[3] = a[4] = a[5] = 
	a[6] = a[7] = a[8] = 0.0;
}


Mat3D::Mat3D( double a11, double a12, double a13,
	      double a21, double a22, double a23,
	      double a31, double a32, double a33 )
{
    a[0] = a11;
    a[1] = a12;
    a[2] = a13;

    a[3] = a21;
    a[4] = a22;
    a[5] = a23;

    a[6] = a31;
    a[7] = a32;
    a[8] = a33;
}


double Mat3D::determinant( void ) const
{
    double x = a[0]*(a[4]*a[8]-a[5]*a[7]);
    double y = a[1]*(a[3]*a[8]-a[5]*a[6]);
    double z = a[2]*(a[3]*a[7]-a[4]*a[6]);
    return( x-y+z );
}


Mat3D Mat3D::inverse( void ) const
{
    double idet = determinant();
    if( idet == 0.0 )
	throw( Error( ERROR_LOCATION, "can't invert matrix: zero determinant" ) );
    idet = 1.0/idet;

    Mat3D result( (a[4]*a[8]-a[5]*a[7])*idet,
		  (a[2]*a[7]-a[1]*a[8])*idet,
		  (a[1]*a[5]-a[2]*a[4])*idet,
		  (a[5]*a[6]-a[3]*a[8])*idet,
		  (a[0]*a[8]-a[2]*a[6])*idet,
		  (a[2]*a[3]-a[0]*a[5])*idet,
		  (a[3]*a[7]-a[4]*a[6])*idet,
		  (a[1]*a[6]-a[0]*a[7])*idet,
		  (a[0]*a[4]-a[1]*a[3])*idet );
    return( result );
}


Mat3D Mat3D::inverse( double det ) const
{
    double idet = 1.0/det;
    Mat3D result( (a[4]*a[8]-a[5]*a[7])*idet,
		  (a[2]*a[7]-a[1]*a[8])*idet,
		  (a[1]*a[5]-a[2]*a[4])*idet,
		  (a[5]*a[6]-a[3]*a[8])*idet,
		  (a[0]*a[8]-a[2]*a[6])*idet,
		  (a[2]*a[3]-a[0]*a[5])*idet,
		  (a[3]*a[7]-a[4]*a[6])*idet,
		  (a[1]*a[6]-a[0]*a[7])*idet,
		  (a[0]*a[4]-a[1]*a[3])*idet );
    return( result );
}


Vec3D Mat3D::operator*( const Vec3D &x ) const
{
    Vec3D result( a[0]*x[0] + a[1]*x[1] + a[2]*x[2],
		  a[3]*x[0] + a[4]*x[1] + a[5]*x[2],
		  a[6]*x[0] + a[7]*x[1] + a[8]*x[2] );
    return( result );
}


std::ostream &operator<<( std::ostream &os, const Mat3D &m ) 
{
    os << std::setw(12) << to_string(m(0)).substr(0,12) << " ";
    os << std::setw(12) << to_string(m(1)).substr(0,12) << " ";
    os << std::setw(12) << to_string(m(2)).substr(0,12) << "\n";

    os << std::setw(12) << to_string(m(3)).substr(0,12) << " ";
    os << std::setw(12) << to_string(m(4)).substr(0,12) << " ";
    os << std::setw(12) << to_string(m(5)).substr(0,12) << "\n";

    os << std::setw(12) << to_string(m(6)).substr(0,12) << " ";
    os << std::setw(12) << to_string(m(7)).substr(0,12) << " ";
    os << std::setw(12) << to_string(m(8)).substr(0,12) << "\n";
    return( os );
}

