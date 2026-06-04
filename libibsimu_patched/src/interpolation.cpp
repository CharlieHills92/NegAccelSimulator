/*! \file interpolation.cpp
 *  \brief Two dimensional interpolation
 */

/* Copyright (c) 2005-2009,2011,2012 Taneli Kalvas. All rights reserved.
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

#include <cmath>
#include <limits>


#include "interpolation.hpp"
#include "error.hpp"



Interpolation2D::Interpolation2D( size_t n, size_t m, const std::vector<double> &f )
    : _n(n), _m(m)
{
    if( _n*_m != f.size() )
	throw( Error( ERROR_LOCATION, "data size not equal to n*m" ) );

    _f = f;
}


const double &Interpolation2D::__f( int i, int j ) const 
{
    return( _f[i+j*_n] );
}


double &Interpolation2D::__f( int i, int j )
{
    return( _f[i+j*_n] );
}


/* ********************************************************************************************
 * CLOSEST
 * ******************************************************************************************** */


ClosestInterpolation2D::ClosestInterpolation2D( size_t n, size_t m, const std::vector<double> &f )
    : Interpolation2D(n,m,f)
{
}


double ClosestInterpolation2D::operator()( double x, double y ) const
{
    int i = (int)floor( x*(_n-1) + 0.5 );
    int j = (int)floor( y*(_m-1) + 0.5 );

    if( i < 0 || i >= (int)_n || j < 0 || j >= (int)_m )
	return( std::numeric_limits<double>::quiet_NaN() );

    return( _f[i+j*_n] );
}





/* ********************************************************************************************
 * BILINEAR
 * ******************************************************************************************** */


BiLinearInterpolation2D::BiLinearInterpolation2D( size_t n, size_t m, const std::vector<double> &f )
    : Interpolation2D(n,m,f)
{
}


double BiLinearInterpolation2D::operator()( double x, double y ) const
{
    double o = x*(_n-1);
    double p = y*(_m-1);
	
    int i = (int)floor( o );
    int j = (int)floor( p );
	
    o -= i;
    p -= j;
	
    if( i < 0 ) {
	i = 0;
	o = 0.0;
    } else if( i >= (int)(_n-1) ) {
	i = _n-2;
	o = 1.0;
    }
    
    if( j < 0 ) {
	j = 0;
	p = 0.0;
    } else if( j >= (int)(_m-1) ) {
	j = _m-2;
	p = 1.0;
    }
    
    return( (1-o)*(1-p)*_f[i+  j*_n] +
	       o *(1-p)*_f[i+1+j*_n] +
	    (1-o)*   p *_f[i+  (j+1)*_n] +
	       o *   p *_f[i+1+(j+1)*_n] );
}





/* ********************************************************************************************
 * BICUBIC
 * ******************************************************************************************** */


const double BiCubicInterpolation2D::wt[16][16] = {
    { 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    { 0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    {-3,  3,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    { 2, -2,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    { 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0},
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0},
    { 0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0},
    { 0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0},
    {-3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0},
    { 0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0},
    { 9, -9, -9,  9,  6,  3, -6, -3,  6, -6,  3, -3,  4,  2,  2,  1},
    {-6,  6,  6, -6, -3, -3,  3,  3, -4,  4, -2,  2, -2, -2, -1, -1},
    { 2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0},
    { 0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0},
    {-6,  6,  6, -6, -4, -2,  4,  2, -3,  3, -3,  3, -2, -1, -2, -1},
    { 4, -4, -4,  4,  2,  2, -2, -2,  2, -2,  2, -2,  1,  1,  1,  1} 
};


const double &BiCubicInterpolation2D::__fx( int i, int j ) const
{
    return( _fx[i+j*_n] );
}


double &BiCubicInterpolation2D::__fx( int i, int j )
{
    return( _fx[i+j*_n] );
}


const double &BiCubicInterpolation2D::__fy( int i, int j ) const
{
    return( _fy[i+j*_n] );
}


double &BiCubicInterpolation2D::__fy( int i, int j )
{
    return( _fy[i+j*_n] );
}


const double &BiCubicInterpolation2D::__fxy( int i, int j ) const
{
    return( _fxy[i+j*_n] );
}


double &BiCubicInterpolation2D::__fxy( int i, int j )
{
    return( _fxy[i+j*_n] );
}


void BiCubicInterpolation2D::calc_coefs( double *c, double *x )
{
    double t;

    // Multiply x with wt from the left and save the result vector to c
    for( size_t i = 0; i < 16; i++ ) {
	
	t = 0.0;
	for( size_t j = 0; j < 16; j++ )
	    t += wt[i][j] * x[j];
	c[i] = t;
    }
}


BiCubicInterpolation2D::BiCubicInterpolation2D( size_t n, size_t m, const std::vector<double> &f )
    : Interpolation2D(n,m,f)
{
    if( n < 3 || m < 3 )
	throw( Error( ERROR_LOCATION, "too small mesh for bicubic interpolation" ) );

    _fx.resize( _f.size() );
    _fy.resize( _f.size() );
    _fxy.resize( _f.size() );
    _c.resize( 16*(_n-1)*(_m-1) );
    
    // Calculate derivatives
    for( size_t i = 1; i < _n-1; i++ ) {
	for( size_t j = 1; j < _m-1; j++ ) {
	    __fx (i,j) = 0.5*(__f(i+1,j)-__f(i-1,j));
	    __fy (i,j) = 0.5*(__f(i,j+1)-__f(i,j-1));
	    __fxy(i,j) = 0.25*(__f(i+1,j+1)-__f(i-1,j+1)-__f(i+1,j-1)+__f(i-1,j-1));
	}
    }

    // Boundaries
    for( size_t i = 0; i < _n; i++ ) {
	__fx (i,0)    = 0.0;
	__fx (i,_m-1) = 0.0;
	__fy (i,0)    = 0.0;
	__fy (i,_m-1) = 0.0;
	__fxy(i,0)    = 0.0;
	__fxy(i,_m-1) = 0.0;
    }
    for( size_t j = 0; j < _m; j++ ) {
	__fx (0,   j) = 0.0;
	__fx (_n-1,j) = 0.0;
	__fy (0,   j) = 0.0;
	__fy (_n-1,j) = 0.0;
	__fxy(0,   j) = 0.0;
	__fxy(_n-1,j) = 0.0;
    }

    // Calculate interpolation 
    double x[16];
    for( size_t i = 0; i < _n-1; i++ ) {
	for( size_t j = 0; j < _m-1; j++ ) {

	    x[0]  = __f  (i,  j  );
	    x[1]  = __f  (i+1,j  );
	    x[2]  = __f  (i,  j+1);
	    x[3]  = __f  (i+1,j+1);

	    x[4]  = __fx (i,  j  );
	    x[5]  = __fx (i+1,j  );
	    x[6]  = __fx (i,  j+1);
	    x[7]  = __fx (i+1,j+1);
		
	    x[8]  = __fy (i,  j  );
	    x[9]  = __fy (i+1,j  );
	    x[10] = __fy (i,  j+1);
	    x[11] = __fy (i+1,j+1);

	    x[12] = __fxy(i,  j  );
	    x[13] = __fxy(i+1,j  );
	    x[14] = __fxy(i,  j+1);
	    x[15] = __fxy(i+1,j+1);

	    calc_coefs( &_c[16*(i+j*(_n-1))], x );
	}
    }
}


double BiCubicInterpolation2D::operator()( double x, double y ) const
{
    double o = x*(_n-1);
    double p = y*(_m-1);
	
    int i = (int)floor( o );
    int j = (int)floor( p );
	
    o -= i;
    p -= j;
	
    if( i < 0 ) {
	i = 0;
	o = 0.0;
    } else if( i >= (int)(_n-1) ) {
	i = _n-2;
	o = 1.0;
    }
    
    if( j < 0 ) {
	j = 0;
	p = 0.0;
    } else if( j >= (int)(_m-1) ) {
	j = _m-2;
	p = 1.0;
    }

    int cc = 16*(i+j*(_n-1));
    double o2 = o*o;
    double o3 = o2*o;
    double p2 = p*p;
    double p3 = p2*p;
    double val = 0.0;
    val +=   _c[cc+0] +  _c[cc+1]*o +  _c[cc+2]*o2 +  _c[cc+3]*o3;
    val += ( _c[cc+4] +  _c[cc+5]*o +  _c[cc+6]*o2 +  _c[cc+7]*o3)*p;
    val += ( _c[cc+8] +  _c[cc+9]*o + _c[cc+10]*o2 + _c[cc+11]*o3)*p2;
    val += (_c[cc+12] + _c[cc+13]*o + _c[cc+14]*o2 + _c[cc+15]*o3)*p3;

    return( val );
}

