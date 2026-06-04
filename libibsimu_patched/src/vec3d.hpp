/*! \file vec3d.hpp
 *  \brief Three dimensional vectors.
 */

/* Copyright (c) 2005-2012,2016 Taneli Kalvas. All rights reserved.
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

#ifndef VEC3D_HPP
#define VEC3D_HPP 1


#include <math.h>
#include <stdint.h>
#include <iostream>
#include <iostream>
#include <iomanip>
#include "vec4d.hpp"
#include "file.hpp"


/*! \brief Three dimensional vector.
 */
class Vec3D {

    double p[3];

public:

    Vec3D() { p[0] = 0.0; p[1] = 0.0; p[2] = 0.0; }
    Vec3D( double x ) { p[0] = x; p[1] = 0.0; p[2] = 0.0; }
    Vec3D( double x, double y ) { p[0] = x; p[1] = y; p[2] = 0.0; }
    Vec3D( double x, double y, double z ) { p[0] = x; p[1] = y; p[2] = z; }

    Vec3D( const class Vec4D &vec );

    Vec3D( std::istream &s ) {
	p[0] = read_double( s );
	p[1] = read_double( s );
	p[2] = read_double( s );
    }
    ~Vec3D() {}

    double &operator[]( int i ) { return( p[i] ); }
    const double &operator[]( int i ) const { return( p[i] ); }
    double &operator()( int i ) { return( p[i] ); }
    const double &operator()( int i ) const { return( p[i] ); }

    /*! \brief %Vector addition
     */
    Vec3D operator+( const Vec3D &vec ) const { 
	return( Vec3D( p[0] + vec[0], 
		       p[1] + vec[1],
		       p[2] + vec[2] ) );
    }

    /*! \brief %Vector difference
     */
    Vec3D operator-( const Vec3D &vec ) const {
	return( Vec3D( p[0] - vec[0],
		       p[1] - vec[1],
		       p[2] - vec[2] ) );
    } 

    /*! \brief %Vector accumulation
     */
    Vec3D &operator+=( const Vec3D &vec ) { 
	p[0] += vec[0];
	p[1] += vec[1];
	p[2] += vec[2];
	return( *this );
    }

    /*! \brief Dot product
     */
    double operator*( const Vec3D &vec ) const { 
	return( p[0] * vec[0] +
		p[1] * vec[1] +
		p[2] * vec[2] );
    }

    /*! \brief %Vector scaling.
     */
    Vec3D operator*( double x ) const { 
	return( Vec3D( x*p[0], x*p[1], x*p[2] ) );
    }

    /*! \brief Unary minus.
     */
    Vec3D operator-( void ) const { 
	return( Vec3D( -p[0], -p[1], -p[2] ) );
    }

    /*! \brief %Vector scaling.
     */
    Vec3D &operator*=( double x ) { 
	p[0] *= x;
	p[1] *= x;
	p[2] *= x;
	return( *this );
    }

    /*! \brief %Vector scaling with divisor.
     */
    Vec3D &operator/=( double x ) { 
	double div = 1.0/x;
	p[0] *= div;
	p[1] *= div;
	p[2] *= div;
	return( *this );
    }

    /*! \brief Inequality test.
     *
     *  Require exact equality.
     */
    bool operator!=( const Vec3D &x ) const;

    /*! \brief Equality test.
     *
     *  Requires exact equality.
     */
    bool operator==( const Vec3D &x ) const;

    /*! \brief Approximate equality test.
     *
     *  Does not require exact equality, but absolute or relative
     *  error less than eps (which ever is less strict). Be careful
     *  using this function!
     */
    bool approx( const Vec3D &x, double eps = 1.0e-6 ) const;

    /*! \brief Assignment.
     */
    Vec3D &operator=( const Vec3D &x ) { 
	p[0] = x[0];
	p[1] = x[1];
	p[2] = x[2];
	return( *this );
    }

    /*! \brief Assignment of every coordinate.
     */
    Vec3D &operator=( const double &x ) { 
	p[0] = x;
	p[1] = x;
	p[2] = x;
	return( *this );
    }

    /*! \brief Calculate absolute value of each component.
     */
    void abs( void ) {
	p[0] = fabs( p[0] );
	p[1] = fabs( p[1] );
	p[2] = fabs( p[2] );
    }

    /*! \brief Normalize vector
     */
    void normalize( void ) {
	double inv_norm = 1.0/sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
	p[0] *= inv_norm;
	p[1] *= inv_norm;
	p[2] *= inv_norm;
    }

    /*! \brief Returns 2-norm of vector
     *
     *  \f$ ||x||_2 = \sqrt{ \Sigma_{i=1}^n x_i^2 } \f$
     */
    double norm2( void ) const {
	return( sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] ) );
    }

    /*! \brief Returns inf-norm of vector
     *
     *  Returns maximum component of vector.
     */
    double max( void ) const;

    /*! \brief Returns square of 2-norm of vector
     *
     *  \f$ (||x||_2)^2 = \Sigma_{i=1}^n x_i^2 \f$
     */
    double ssqr( void ) const {
	return( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
    }

    /*! \brief Returns the index of element with minimum magnitude (abs).
     */
    int min_element( void ) const;

    /*! \brief Returns arbitrary vector perpendicular to input vector.
     */
    Vec3D arb_perpendicular( void ) const;

    /*! \brief Saves data to stream \a os.
     */
    void save( std::ostream &os ) const { 
	write_double( os, p[0] );
	write_double( os, p[1] );
	write_double( os, p[2] ); 
    }

    /*! \brief Returns standard basis vector \a i.
     */
    static Vec3D standard_basis( int i );

    /*! \brief Cross product
     */
    friend Vec3D cross( const Vec3D &vec1, const Vec3D &vec2 );

    /*! \brief Second norm of vector.
     */
    friend double norm2( const Vec3D &vec );

    /*! \brief Sum of squares or square of 2-norm of vector.
     */
    friend double ssqr( const Vec3D &vec );

    /*! \brief %Vector scaling.
     */
    friend Vec3D operator*( double x, const Vec3D &vec );

    /*! \brief %Vector scaling for integer vector.
     */
    friend Vec3D operator*( double x, const class Int3D &i );

    /*! \brief Outputting to stream.
     */
    friend std::ostream &operator<<( std::ostream &os, const Vec3D &vec );
    friend std::basic_ostream<wchar_t> &operator<<( std::basic_ostream<wchar_t> &os, const Vec3D &vec );
};


inline double norm2( const Vec3D &vec ) {
    return( vec.norm2() );
}


inline double ssqr( const Vec3D &vec ) {
    return( vec.ssqr() );
}


inline Vec3D cross( const Vec3D &vec1, const Vec3D &vec2 ) { 
    return( Vec3D( vec1[1] * vec2[2] - vec1[2] * vec2[1], 
		   vec1[2] * vec2[0] - vec1[0] * vec2[2],
		   vec1[0] * vec2[1] - vec1[1] * vec2[0] ) );
}


inline Vec3D operator*( double x, const Vec3D &vec )
{
    return( Vec3D( x*vec.p[0], x*vec.p[1], x*vec.p[2] ) );
}


inline std::ostream &operator<<( std::ostream &os, const Vec3D &vec ) 
{
    os << std::setw(12) << vec[0] << " " 
       << std::setw(12) << vec[1] << " "
       << std::setw(12) << vec[2];
    //os << std::setw(12) << to_string(vec[0]).substr(0,12) << " ";
    //os << std::setw(12) << to_string(vec[1]).substr(0,12) << " ";
    //os << std::setw(12) << to_string(vec[2]).substr(0,12);
    return( os );
}


inline std::basic_ostream<wchar_t> &operator<<( std::basic_ostream<wchar_t> &os, const Vec3D &vec )
{
    os << std::setw(12) << vec[0] << " " 
       << std::setw(12) << vec[1] << " "
       << std::setw(12) << vec[2];
    //os << std::setw(12) << to_wstring(vec[0]).substr(0,12) << " ";
    //os << std::setw(12) << to_wstring(vec[1]).substr(0,12) << " ";
    //os << std::setw(12) << to_wstring(vec[2]).substr(0,12);
    return( os );
}


/*! \brief 3D Integer vector class.
 */
class Int3D {
    int32_t l[3];

public:

    Int3D() { l[0] = 0; l[1] = 0; l[2] = 0; }
    Int3D( int32_t i ) { l[0] = i; l[1] = 0; l[2] = 0; }
    Int3D( int32_t i, int32_t j ) { l[0] = i; l[1] = j; l[2] = 0; }
    Int3D( int32_t i, int32_t j, int32_t k ) { l[0] = i; l[1] = j; l[2] = k; }
    Int3D( std::istream &s ) {
	l[0] = read_int32( s );
	l[1] = read_int32( s );
	l[2] = read_int32( s );
    }
    ~Int3D() {}

    int32_t &operator[]( int i ) { return( l[i] ); }
    const int32_t &operator[]( int i ) const { return( l[i] ); }
    int32_t &operator()( int i ) { return( l[i] ); }
    const int32_t &operator()( int i ) const { return( l[i] ); }

    /*! \brief Integer vector addition
     */
    Int3D operator+( const Int3D &i ) const { 
	return( Int3D( l[0] + i[0], 
		       l[1] + i[1],
		       l[2] + i[2] ) );
    }

    /*! \brief Integer vector difference
     */
    Int3D operator-( const Int3D &i ) {
	return( Int3D( l[0] - i[0],
		       l[1] - i[1],
		       l[2] - i[2] ) );
    } 

    Int3D operator*( int i ) { 
	return( Int3D( i*l[0], i*l[1], i*l[2] ) );
    }

    Vec3D operator*( double x ) { 
	return( Vec3D( x*l[0], x*l[1], x*l[2] ) );
    }

    /*! \brief Inequality test.
     */
    bool operator!=( const Int3D &i ) const { 
	if( l[0] != i.l[0] || l[1] != i.l[1] || l[2] != i.l[2] )
	    return( true );
	return( false ); 
    }

    /*! \brief Equality test.
     */
    bool operator==( const Int3D &i ) const { 
	if( l[0] == i.l[0] && l[1] == i.l[1] && l[2] == i.l[2] )
	    return( true );
	return( false ); 
    }


    /*!  \brief Returns maximum component of vector.
     */
    int32_t max( void ) const;

    void save( std::ostream &s ) const { 
	write_int32( s, l[0] );
	write_int32( s, l[1] );
	write_int32( s, l[2] ); 
    }

    friend Vec3D operator*( double x, const Int3D &i );
    friend Int3D operator*( int x, const Int3D &i );
    friend std::ostream &operator<<( std::ostream &os, const Int3D &vec );
    friend std::basic_ostream<wchar_t> &operator<<( std::basic_ostream<wchar_t> &os, const Int3D &vec );
};


inline Vec3D operator*( double x, const Int3D &i )
{
    Vec3D res;
    res[0] = x*i.l[0];
    res[1] = x*i.l[1];
    res[2] = x*i.l[2];
    return( res );
}


inline Int3D operator*( int x, const Int3D &i )
{
    Int3D res;
    res[0] = x*i.l[0];
    res[1] = x*i.l[1];
    res[2] = x*i.l[2];
    return( res );
}


inline std::ostream &operator<<( std::ostream &os, const Int3D &vec ) 
{
    os << std::setw(12) << vec[0] << " " 
       << std::setw(12) << vec[1] << " "
       << std::setw(12) << vec[2];
    //os << std::setw(12) << to_string(vec[0]).substr(0,12) << " ";
    //os << std::setw(12) << to_string(vec[1]).substr(0,12) << " ";
    //os << std::setw(12) << to_string(vec[2]).substr(0,12);
    return( os );
}


inline std::basic_ostream<wchar_t> &operator<<( std::basic_ostream<wchar_t> &os, const Int3D &vec ) 
{
    os << std::setw(12) << vec[0] << " " 
       << std::setw(12) << vec[1] << " "
       << std::setw(12) << vec[2];
    //os << std::setw(12) << to_wstring(vec[0]).substr(0,12) << " ";
    //os << std::setw(12) << to_wstring(vec[1]).substr(0,12) << " ";
    //os << std::setw(12) << to_wstring(vec[2]).substr(0,12);
    return( os );
}


#endif

