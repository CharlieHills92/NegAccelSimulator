/*! \file vec4d.hpp
 *  \brief Homogenous vectors for three dimensional space.
 */

/* Copyright (c) 2005-2010 Taneli Kalvas. All rights reserved.
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

#ifndef VEC4D_HPP
#define VEC4D_HPP 1


#include <math.h>
#include <stdint.h>
#include <iostream>
#include <iostream>
#include <iomanip>
#include "vec3d.hpp"
#include "file.hpp"
#include "error.hpp"


/*! \brief Homogenous vector for three dimensional space
 *
 *  Homogenous space has 4-vectors (x,y,z,w).
 *
 *  Most operations assume the arguments are either vectors (w=0) or
 *  points (w=1). These are only checked where the algorithm is
 *  dependent on this information. Otherwise the fourth coordinate (w)
 *  is assumed what the algorithm is designed for. User is responsible
 *  for calling legal operations for corrent type of vectors.
 */
class Vec4D {

    double p[4];

public:

    Vec4D() { p[0] = 0.0; p[1] = 0.0; p[2] = 0.0; p[3] = 0.0; }
    Vec4D( double x ) { p[0] = x; p[1] = 0.0; p[2] = 0.0; p[3] = 0.0; }
    Vec4D( double x, double y ) { p[0] = x; p[1] = y; p[2] = 0.0; p[3] = 0.0; }
    Vec4D( double x, double y, double z ) { p[0] = x; p[1] = y; p[2] = z; p[3] = 0.0; }
    Vec4D( double x, double y, double z, double w ) { p[0] = x; p[1] = y; p[2] = z; p[3] = w; }

    /*! \brief Convert 3D vector to 4D vector 
     *
     *  Makes 3D vector a point.
     */
    Vec4D( const class Vec3D &vec );

    Vec4D( std::istream &s ) {
	p[0] = read_double( s );
	p[1] = read_double( s );
	p[2] = read_double( s );
	p[3] = read_double( s );
    }
    ~Vec4D() {}

    double &operator[]( int i ) { return( p[i] ); }
    const double &operator[]( int i ) const { return( p[i] ); }
    double &operator()( int i ) { return( p[i] ); }
    const double &operator()( int i ) const { return( p[i] ); }

    /*! \brief Addition.
     *
     *  Only valid for point+vector=point or
     *  vector+vector=vector. Output is of the correct type.
     */
    Vec4D operator+( const Vec4D &vec ) const { 
	return( Vec4D( p[0] + vec[0], 
		       p[1] + vec[1],
		       p[2] + vec[2],
		       (p[2] == vec[2] ? 0.0 : 1.0) ) );
    }

    /*! \brief Difference
     *
     *  Only valid for vector-vector=vector, point-vector=point or
     *  point-point=vector. Output is of the correct type.
     */
    Vec4D operator-( const Vec4D &vec ) const {
	return( Vec4D( p[0] - vec[0],
		       p[1] - vec[1],
		       p[2] - vec[2],
		       (p[2] == vec[2] ? 0.0 : 1.0) ) );
    }

    /*! \brief Accumulation
     *
     *  Only valid for point += vector or vector += vector. Output
     *  type does not change.
     */
    Vec4D &operator+=( const Vec4D &vec ) { 
	p[0] += vec[0];
	p[1] += vec[1];
	p[2] += vec[2];
	return( *this );
    }

    /*! \brief Dot product
     *
     *  Valid for vectors only.
     */
    double operator*( const Vec4D &vec ) const { 
	return( p[0] * vec[0] +
		p[1] * vec[1] +
		p[2] * vec[2] );
    }

    /*! \brief %Vector scaling.
     *
     *  Valid for points and vectors. Scaling does not affect w.
     */
    Vec4D operator*( double x ) const { 
	return( Vec4D( x*p[0], x*p[1], x*p[2], p[3] ) );
    }

    /*! \brief %Vector scaling.
     *
     *  Valid for points and vectors. Scaling does not affect w.
     */
    Vec4D &operator*=( double x ) { 
	p[0] *= x;
	p[1] *= x;
	p[2] *= x;
	return( *this );
    }

    /*! \brief %Vector scaling with divisor.
     *
     *  Valid for points and vectors. Scaling does not affect w.
     */
    Vec4D &operator/=( double x ) { 
	double div = 1.0/x;
	p[0] *= div;
	p[1] *= div;
	p[2] *= div;
	return( *this );
    }

    /*! \brief Inequality test.
     *
     *  Also tests w.
     */
    bool operator!=( const Vec4D &x ) { 
	if( p[0] != x.p[0] || p[1] != x.p[1] || p[2] != x.p[2] || p[3] != x.p[3] )
	    return( true );
	return( false ); 
    }

    /*! \brief Equality test.
     *
     *  Also tests w.
     */
    bool operator==( const Vec4D &x ) { 
	if( p[0] == x.p[0] && p[1] == x.p[1] && p[2] == x.p[2] && p[3] == x.p[3] )
	    return( true );
	return( false ); 
    }

    /*! \brief Assignment.
     */
    Vec4D &operator=( const Vec4D &x ) { 
	p[0] = x[0];
	p[1] = x[1];
	p[2] = x[2];
	p[3] = x[3];
	return( *this );
    }

    /*! \brief Homogenize vector
     *
     *  Homogenize vector by dividing all components of vector with
     *  w. Not valid for w=0.
     */
    void homogenize() {
	double inv_w = 1.0/p[3];
	p[0] *= inv_w;
	p[1] *= inv_w;
	p[2] *= inv_w;
	p[3] = 1.0;
    }

    /*! \brief Normalize vector
     *
     *  Only valid for vector. Output is guaranteed to be a vector.
     */
    void normalize() {
	double inv_norm = 1.0/sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
	p[0] *= inv_norm;
	p[1] *= inv_norm;
	p[2] *= inv_norm;
	p[3] = 0.0;
    }

    /*! \brief Returns 2-norm of vector
     *
     *  \f$ ||x||_2 = \sqrt{ \Sigma_{i=1}^3 x_i^2 } \f$
     */
    double norm2() const {
	return( sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] ) );
    }

    /*! \brief Returns square of 2-norm of vector
     *
     *  \f$ (||x||_2)^2 = \Sigma_{i=1}^3 x_i^2 \f$
     */
    double ssqr() const {
	return( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
    }

    void save( std::ostream &s ) const { 
	write_double( s, p[0] );
	write_double( s, p[1] );
	write_double( s, p[2] ); 
	write_double( s, p[3] ); 
    }

    /*! \brief Cross product
     *
     *  Only valid for vectors. Output is guaranteed to be a vector.
     */
    friend Vec4D cross( const Vec4D &vec1, const Vec4D &vec2 );

    /*! \brief Second norm of vector.
     */
    friend double norm2( const Vec4D &vec );

    /*! \brief %Vector scaling.
     *
     *  Does not affect w.
     */
    friend Vec4D operator*( double x, const Vec4D &vec );

    /*! \brief Outputting to stream.
     */
    friend std::ostream &operator<<( std::ostream &os, const Vec4D &vec );
};


inline double norm2( const Vec4D &vec ) {
    return( vec.norm2() );
}

inline Vec4D cross( const Vec4D &vec1, const Vec4D &vec2 ) { 
    return( Vec4D( vec1[1] * vec2[2] - vec1[2] * vec2[1], 
		   vec1[2] * vec2[0] - vec1[0] * vec2[2],
		   vec1[0] * vec2[1] - vec1[1] * vec2[0],
		   0.0 ) );
}


inline Vec4D operator*( double x, const Vec4D &vec )
{
    return( Vec4D( x*vec[0], x*vec[1], x*vec[2], vec[3] ) );
}


inline std::ostream &operator<<( std::ostream &os, const Vec4D &vec ) 
{
    os << std::setw(12) << to_string(vec[0]).substr(0,12) << " ";
    os << std::setw(12) << to_string(vec[1]).substr(0,12) << " ";
    os << std::setw(12) << to_string(vec[2]).substr(0,12) << " ";
    os << std::setw(12) << to_string(vec[3]).substr(0,12);
    return( os );
}


#endif



















