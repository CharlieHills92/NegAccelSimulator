/*! \file vectorfield.cpp
 *  \brief %Vector field base
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


#include <limits>
#include "vectorfield.hpp"


void VectorField::get_minmax( const Mesh &mesh, double &min, double &max ) const
{
    min = std::numeric_limits<double>::infinity();
    max = -std::numeric_limits<double>::infinity();

    for( uint32_t k = 0; k < mesh.size(2); k++ ) {
	for( uint32_t j = 0; j < mesh.size(1); j++ ) {
	    for( uint32_t i = 0; i < mesh.size(0); i++ ) {

		Vec3D x( mesh.origo(0)+i*mesh.h(),
			 mesh.origo(1)+i*mesh.h(),
			 mesh.origo(2)+i*mesh.h() );
		double val = (*this)( x ).ssqr();
		if( val < min )
		    min = val;
		if( val > max )
		    max = val;
	    }
	}
    }
    min = sqrt( min );
    max = sqrt( max );
}


void VectorField::get_minmax( const Mesh &mesh, Vec3D &min, Vec3D &max ) const
{
    min = Vec3D( std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    max = Vec3D( -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity() );

    for( uint32_t k = 0; k < mesh.size(2); k++ ) {
	for( uint32_t j = 0; j < mesh.size(1); j++ ) {
	    for( uint32_t i = 0; i < mesh.size(0); i++ ) {

		Vec3D x( mesh.origo(0)+i*mesh.h(),
			 mesh.origo(1)+i*mesh.h(),
			 mesh.origo(2)+i*mesh.h() );
		Vec3D val = (*this)( x ).ssqr();
		for( size_t b = 0; b < 3; b++ ) {
		    if( val[b] < min[b] )
			min[b] = val[b];
		    if( val[b] > max[b] )
			max[b] = val[b];
		}
	    }
	}
    }
}
