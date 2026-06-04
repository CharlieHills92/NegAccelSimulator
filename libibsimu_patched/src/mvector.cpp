/*! \file mvector.cpp
 *  \brief N-dimensional vector
 */

/* Copyright (c) 2005-2009,2012 Taneli Kalvas. All rights reserved.
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

#include "config.hpp"
#include "mvector.hpp"
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iomanip>
#ifdef USE_BLAS
extern "C" {
#include BLAS_HEADER
}
#endif



/* ******************* *
 * VectorLA            *
 * ******************* */

Vector::VectorLA Vector::VectorLA::operator+( const Vector::VectorLA &vecla ) const 
{
    if( _refs[0]._vec->_n != vecla._refs[0]._vec->_n )
	throw( ErrorDim( ERROR_LOCATION ) );
    Vector::VectorLA res( *this );
    res._refs.insert( res._refs.end(), vecla._refs.begin(), vecla._refs.end() );
    return( res );
}


Vector::VectorLA Vector::VectorLA::operator-( const Vector::VectorLA &vecla ) const 
{
    if( _refs[0]._vec->_n != vecla._refs[0]._vec->_n )
	throw( ErrorDim( ERROR_LOCATION ) );
    Vector::VectorLA res( *this );
    int origsize = res._refs.size();
    res._refs.insert( res._refs.end(), vecla._refs.begin(), vecla._refs.end() );
    std::vector<VectorRef>::iterator itend = res._refs.end();
    for( std::vector<VectorRef>::iterator it = res._refs.begin()+origsize; it != itend; it++ )
	it->_coef *= -1.0;
    return( res );
}


Vector::VectorLA Vector::VectorLA::operator-() const 
{
    Vector::VectorLA res( *this );
    std::vector<VectorRef>::iterator itend = res._refs.end();
    for( std::vector<VectorRef>::iterator it = res._refs.begin(); it != itend; it++ )
	it->_coef *= -1.0;
    return( res );
}


Vector::VectorLA Vector::VectorLA::operator*( double x ) const 
{
    Vector::VectorLA res( *this );
    std::vector<VectorRef>::iterator itend = res._refs.end();
    for( std::vector<VectorRef>::iterator it = res._refs.begin(); it != itend; it++ )
	it->_coef *= x;
    return( res );
}


Vector::VectorLA operator*( double x, const Vector::VectorLA &vecla )
{
    Vector::VectorLA res( vecla );
    std::vector<Vector::VectorRef>::iterator itend = res._refs.end();
    for( std::vector<Vector::VectorRef>::iterator it = res._refs.begin(); it != itend; it++ )
	it->_coef *= x;
    return( res );
}




/* ******************* *
 * Vector              *
 * ******************* */

inline void Vector::allocate( void )
{
    if( !(_val = (double *)malloc( _n*sizeof(double) )) ) {
	_n = 0;
	throw( ErrorNoMem( ERROR_LOCATION ) );
    }
}


inline void Vector::callocate( void )
{
    if( !(_val = (double *)calloc( _n, sizeof(double) )) ) {
	_n = 0;
	throw( ErrorNoMem( ERROR_LOCATION ) );
    }
}


inline void Vector::reallocate( void )
{
    if( _n == 0 ) {
	_val = NULL;
	return;
    }

    double *tmp;
    if( !(tmp = (double *)realloc( _val, _n*sizeof(double) )) ) {
	free( _val );
	_n = 0;
	throw( ErrorNoMem( ERROR_LOCATION ) );
    }
    _val = tmp;
}


Vector::Vector( int n ) 
{
    _n = n;
    callocate();
}


Vector::Vector( int n, const double *val ) 
{
    _n = n;
    allocate();
    memcpy( _val, val, _n*sizeof(double) );
}


Vector::Vector( int n, double val )
{
    _n = n;
    allocate();
    for( int i = 0; i < _n; i++ )
	_val[i] = val;
}


Vector::Vector( const Vector &vec ) 
{
    _n = vec._n;
    allocate();
    memcpy( _val, vec._val, _n*sizeof(double) );
}


Vector::Vector( const Vector::VectorLA &vecla ) 
{
    int refsize = vecla._refs.size();
    _n = vecla._refs[0]._vec->_n;
    callocate();
#ifdef USE_BLAS
    for( int a = 0; a < refsize; a++ )
	BLAS(daxpy)( n, vecla._refs[a]._coef, vecla._refs[a]._vec->_val, 1, _val, 1 );
#else
    for( int i = 0; i < _n; i++ ) {
	_val[i] = vecla._refs[0]._coef * vecla._refs[0]._vec->_val[i];
	for( int a = 1; a < refsize; a++ )
	    _val[i] += vecla._refs[a]._coef * vecla._refs[a]._vec->_val[i];
    }
#endif
}


Vector::Vector( const struct MatrixMulVec &matvec )
{
    _n = 0;
    _val = NULL;
    matvec._mat->multiply_by_vector( *this, *matvec._vec );
}


Vector::~Vector() 
{
    free( _val );
}


void Vector::resize( int n )
{
    if( _n != n ) {
	_n = n;
	reallocate();
    }
}


void Vector::clear( void )
{
    memset( _val, 0, _n*sizeof(double) );    
}


void Vector::merge( Vector &vec )
{
    _n   = vec._n;
    _val = vec._val;

    vec._n   = 0;
    vec._val = NULL;
}


Vector::VectorLA Vector::operator+( const Vector::VectorLA &vecla ) const 
{
    if( _n != vecla._refs[0]._vec->_n )
	throw( ErrorDim( ERROR_LOCATION ) );
    Vector::VectorLA res( *this );
    res._refs.insert( res._refs.end(), vecla._refs.begin(), vecla._refs.end() );
    return( res );
}


Vector::VectorLA Vector::operator-( const Vector::VectorLA &vecla ) const 
{
    if( _n != vecla._refs[0]._vec->_n )
	throw( ErrorDim( ERROR_LOCATION ) );
    Vector::VectorLA res( *this );
    int origsize = res._refs.size();
    res._refs.insert( res._refs.end(), vecla._refs.begin(), vecla._refs.end() );
    std::vector<VectorRef>::iterator itend = res._refs.end();
    for( std::vector<VectorRef>::iterator it = res._refs.begin()+origsize; it != itend; it++ )
	it->_coef *= -1.0;
    return( res );
}


Vector::VectorLA Vector::operator-() const
{
    Vector::VectorLA res( *this, -1.0 );
    return( res );
}


Vector::VectorLA Vector::operator*( double x ) const
{
    Vector::VectorLA res( *this, x );
    return( res );
}


Vector &Vector::operator+=( const VectorLA &vecla )
{
    if( _n != vecla._refs[0]._vec->_n )
	throw( ErrorDim( ERROR_LOCATION ) );
#ifdef USE_BLAS
    for( int a = 0; a < vecla._refs.size(); a++ )
	BLAS(daxpy)( _n, vecla._refs[a]._coef, vecla._refs[a]._vec->_val, 1, val, 1 );
#else
    int refsize = vecla._refs.size();
    for( int i = 0; i < _n; i++ ) {
	for( int a = 0; a < refsize; a++ )
	    _val[i] += vecla._refs[a]._coef * vecla._refs[a]._vec->_val[i];
    }
#endif
    return( *this );
}


Vector &Vector::operator-=( const VectorLA &vecla )
{
    if( _n != vecla._refs[0]._vec->_n )
	throw( ErrorDim( ERROR_LOCATION ) );
#ifdef USE_BLAS
    for( int a = 0; a < vecla._refs.size(); a++ )
	BLAS(daxpy)( _n, -vecla._refs[a]._coef, vecla._refs[a]._vec->_val, 1, _val, 1 );
#else
    int refsize = vecla._refs.size();
    for( int i = 0; i < _n; i++ ) {
	for( int a = 0; a < refsize; a++ )
	    _val[i] -= vecla._refs[a]._coef * vecla._refs[a]._vec->_val[i];
    }
#endif
    return( *this );
}


Vector &Vector::operator*=( double x )
{
#ifdef USE_BLAS
    BLAS(dscal)( _n, x, _val, 1 );
#else
    for( int i = 0; i < _n; i++ )
	_val[i] *= x;
#endif
    return( *this );
}


Vector &Vector::operator=( double x )
{
    for( int i = 0; i < _n; i++ )
	_val[i] = x;
    return( *this );
}


Vector &Vector::operator=( const Vector &vec ) 
{
    if( _n != vec._n ) {
	_n = vec._n;
	reallocate();
    }
    memcpy( _val, vec._val, _n*sizeof(double) );
    return( *this );
}


Vector &Vector::operator=( const Vector::VectorLA &vecla ) 
{
    int refsize = vecla._refs.size();
    int a;
    bool tmp = 0;
    double *res;
    
    // Check if result is same as some argument
    for( a = 0; a < refsize; a++ )
	if( this == vecla._refs[a]._vec )
	    break;
    if( a != refsize ) {
	// Temp needed for result
	tmp = 1;
	if( !(res = (double *)malloc( _n*sizeof(double) )) )
	    throw( ErrorNoMem( ERROR_LOCATION ) );
    } else {
	// No temp needed
	if( _n != vecla._refs[0]._vec->_n ) {
	    _n = vecla._refs[0]._vec->_n;
	    reallocate();
	}
	res = _val;
    }

#ifdef USE_BLAS
    memset( res, 0, _n*sizeof(double) );
    for( a = 0; a < vecla._refs.size(); a++ )
	BLAS(daxpy)( _n, vecla._refs[a]._coef, vecla._refs[a]._vec->_val, 1, _res, 1 );
#else
    for( int i = 0; i < _n; i++ ) {
	res[i] = vecla._refs[0]._coef * vecla._refs[0]._vec->_val[i];
	for( a = 1; a < refsize; a++ )
	    res[i] += vecla._refs[a]._coef * vecla._refs[a]._vec->_val[i];
    }
#endif

    if( tmp ) {
	free( _val );
	_val = res;
    }

    return( *this );
}


Vector &Vector::operator=( const struct MatrixMulVec &matvec )
{
    matvec._mat->multiply_by_vector( *this, *matvec._vec );
    return( *this );
}


bool Vector::operator==( const Vector &vec ) const
{
    if( _n != vec._n )
	return( false );
    for( int i = 0; i < _n; i++ ) {
	if( _val[i] != vec._val[i] )
	    return( false );
    }
    return( true );
}


bool Vector::operator!=( const Vector &vec ) const
{
    if( _n != vec._n )
	return( true );
    for( int i = 0; i < _n; i++ ) {
	if( _val[i] != vec._val[i] )
	    return( true );
    }
    return( false );
}


Vector::VectorLA operator*( double x, Vector &vec ) 
{
    Vector::VectorLA res( vec, x );    
    return( res );
}


std::ostream &operator<<( std::ostream &os, const Vector &vec ) 
{
    for( int i = 0; i < vec._n; i++ )
	os << std::setw(6) << to_string(vec(i)).substr(0,6) << " ";
    return( os );
}


double dot_prod( const Vector &vec1, const Vector &vec2 )
{
    if( vec1._n != vec2._n )
	throw( ErrorDim( ERROR_LOCATION ) );
#ifdef USE_BLAS
    return( BLAS(ddot)( vec1._n, vec1._val, 1, vec2._val, 1 ) );
#else

    double res = 0.0;
    double *ptr1 = vec1._val;
    double *ptr1end = &vec1._val[vec1._n];
    double *ptr2 = vec2._val;
    while( ptr1 != ptr1end )
	res += (*(ptr1++)) * (*(ptr2++));
    return( res );
#endif
}


double norm1( const Vector &vec )
{
#ifdef USE_BLAS
    return( BLAS(dasum)( vec._n, vec._val, 1 ) );
#else
    double res = 0.0;
    for( int i = 0; i < vec._n; i++ )
	res += fabs( vec._val[i] );
    return( res );
#endif
}


double norm2( const Vector &vec )
{
#ifdef USE_BLAS
    return( BLAS(dnrm2)( vec._n, vec._val, 1 ) );
#else
    double res = 0.0;
    for( int i = 0; i < vec._n; i++ )
	res += vec._val[i]*vec._val[i];
    return( sqrt( res ) );
#endif
}


double ssqr( const Vector &vec )
{
    double res = 0.0;
    for( int i = 0; i < vec._n; i++ )
	res += vec._val[i]*vec._val[i];
    return( res );
}


double min( const Vector &vec )
{
    double res = std::numeric_limits<double>::infinity();
    for( int i = 0; i < vec._n; i++ ) {
	if( vec._val[i] < res )
	    res = vec._val[i];
    }
    return( res );
}


double min_abs( const Vector &vec )
{
    double res = std::numeric_limits<double>::infinity();
    double x;
    for( int i = 0; i < vec._n; i++ ) {
	x = fabs( vec._val[i] );
	if( x < res )
	    res = x;
    }
    return( res );
}


double max( const Vector &vec )
{
    double res = -std::numeric_limits<double>::infinity();
    for( int i = 0; i < vec._n; i++ ) {
	if( vec._val[i] > res )
	    res = vec._val[i];
    }
    return( res );
}


double max_abs( const Vector &vec )
{
    double res = -std::numeric_limits<double>::infinity();
    double x;
    for( int i = 0; i < vec._n; i++ ) {
	x = fabs( vec._val[i] );
	if( x > res )
	    res = x;
    }
    return( res );
}


void swap( Vector &vec1, Vector &vec2 )
{
    int tn       = vec1._n;
    double *tval = vec1._val;

    vec1._n   = vec2._n;
    vec1._val = vec2._val;

    vec2._n   = tn;
    vec2._val = tval;
}

