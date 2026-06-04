/*! \file crowmatrix.cpp
 *  \brief Compressed row sparse matrix algebra
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

#include "config.hpp"
#include "crowmatrix.hpp"
#include "ccolmatrix.hpp"
#include "coordmatrix.hpp"
#include "mvector.hpp"
#include "sort.hpp"
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iomanip>


inline void CRowMatrix::allocate( void )
{
    if( !(_ptr = (int *)malloc( (_n+1)*sizeof(int) )) ) {
	_ptr = _col = NULL;
	_val = NULL;
	_n = _m = _nz = _asize = 0;	
	throw( ErrorNoMem( ERROR_LOCATION ) );
    }

    if( _asize == 0 ) {
	_col = NULL;
	_val = NULL;
    } else {
	if( !(_col = (int *)malloc( _asize*sizeof(int) )) ) {
	    free( _ptr );
	    _ptr = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	
	if( !(_val = (double *)malloc( _asize*sizeof(double) )) ) {
	    free( _ptr );
	    free( _col );
	    _ptr = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
    }
}


inline void CRowMatrix::reallocate( void )
{
    int *tmp;
    double *tmp2;

    if( !(tmp = (int *)realloc( _ptr, (_n+1)*sizeof(int) )) ) {
	free( _ptr );
	free( _col );
	free( _val );
	_ptr = _col = NULL;
	_val = NULL;
	_n = _m = _nz = _asize = 0;
	throw( ErrorNoMem( ERROR_LOCATION ) );
    }
    _ptr = tmp;

    if( _asize == 0 ) {
	_col = NULL;
	_val = NULL;
    } else {
	if( !(tmp = (int *)realloc( _col, _asize*sizeof(int) )) ) {
	    free( _ptr );
	    free( _col );
	    free( _val );
	    _ptr = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	_col = tmp;
	if( !(tmp2 = (double *)realloc( _val, _asize*sizeof(double) )) ) {
	    free( _ptr );
	    free( _col );
	    free( _val );
	    _ptr = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	_val = tmp2;
    }
}


CRowMatrix::CRowMatrix()
    : _n(0), _m(0), _nz(0), _asize(0), _col(NULL), _val(NULL)
{
    allocate();
}


CRowMatrix::CRowMatrix( int n, int m )
{
    _n     = n;
    _m     = m;
    _nz    = 0;
    _asize = 0;
    allocate();
    memset( _ptr, 0, (_n+1)*sizeof(int) );
}


CRowMatrix::CRowMatrix( int n, int m, int nz, int asize,
			int *ptr, int *col, double *val )
{
    _n     = n;
    _m     = m;
    _nz    = nz;
    _asize = asize;
    _ptr   = ptr;
    _col   = col;
    _val   = val;
}


void CRowMatrix::build( const class CRowMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    memcpy( _ptr, mat._ptr, (_n+1)*sizeof(int) );
    memcpy( _col, mat._col, _nz*sizeof(int) );
    memcpy( _val, mat._val, _nz*sizeof(double) );
}


CRowMatrix::CRowMatrix( const CRowMatrix &mat )
    : _ptr(NULL), _col(NULL), _val(NULL)
{
    build( mat );
}


CRowMatrix &CRowMatrix::operator=( const CRowMatrix &mat )
{
    build( mat );
    return( *this );
}


void CRowMatrix::build( const class CColMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();
    
    int *c = new int[_n];
    int i, j, ind, rr;

    /* Count number of entries in each row to c. */
    memset( c, 0, _n*sizeof(int) );
    for( i = 0; i < _nz; i++ )
	c[mat._row[i]]++;

    /* Set up ptr. */
    _ptr[0] = 0;
    for( i = 1; i <= _n; i++ )
	_ptr[i] = _ptr[i-1] + c[i-1];

    /* Copy input to output using c as a column pointer. */
    for( i = 0, j = 0; i < _m; i++ ) {
	for( ; j < mat._ptr[i+1]; j++ ) {
	    rr = mat._row[j];
	    ind = _ptr[rr+1] - c[rr];
	    c[rr]--;
	    _col[ind] = i;
	    _val[ind] = mat._val[j];
	}
    }

    delete [] c;
}


CRowMatrix::CRowMatrix( const class CColMatrix &mat )
    : _ptr(NULL), _col(NULL), _val(NULL)
{
    build( mat );
}


CRowMatrix &CRowMatrix::operator=( const class CColMatrix &mat )
{
    build( mat );
    return( *this );
}


void CRowMatrix::build( const class CoordMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    int *c = new int[_n];
    int i, start, r;

    /* Count number of entries in each row to c. */
    memset( c, 0, _n*sizeof(int) );
    for( i = 0; i < _nz; i++ )
	c[mat._row[i]]++;

    /* Set up ptr. */
    _ptr[0] = 0;
    for( i = 1; i <= _n; i++ )
	_ptr[i] = _ptr[i-1] + c[i-1];

    /* Copy input to output using c as a column pointer. */
    for( i = 0; i < _nz; i++ ) {
	r = mat._row[i];
	start = _ptr[r];
	c[r]--;
	_col[start+c[r]] = mat._col[i];
	_val[start+c[r]] = mat._val[i];
    }

    delete [] c;
}


CRowMatrix::CRowMatrix( const class CoordMatrix &mat )
    : _ptr(NULL), _col(NULL), _val(NULL)
{
    build( mat );
}


CRowMatrix &CRowMatrix::operator=( const class CoordMatrix &mat )
{
    build( mat );
    return( *this );
}


CRowMatrix::CRowMatrix( const class Matrix &mat )
    : _ptr(NULL), _col(NULL), _val(NULL)
{
    const CRowMatrix  *crmat;
    const CColMatrix  *ccmat;
    const CoordMatrix *comat;

    if( (crmat = dynamic_cast<const CRowMatrix *>(&mat)) != 0 )
	build( *crmat );
    else if( (ccmat = dynamic_cast<const CColMatrix *>(&mat)) != 0 )
	build( *ccmat );
    else if( (comat = dynamic_cast<const CoordMatrix *>(&mat)) != 0 )
	build( *ccmat );
    else
        throw( ErrorUnimplemented( ERROR_LOCATION, "Couldn't convert unknown matrix type" ) );
}


CRowMatrix &CRowMatrix::operator=( const class Matrix &mat )
{
    const CRowMatrix  *crmat;
    const CColMatrix  *ccmat;
    const CoordMatrix *comat;

    if( (crmat = dynamic_cast<const CRowMatrix *>(&mat)) != 0 )
	build( *crmat );
    else if( (ccmat = dynamic_cast<const CColMatrix *>(&mat)) != 0 )
	build( *ccmat );
    else if( (comat = dynamic_cast<const CoordMatrix *>(&mat)) != 0 )
	build( *ccmat );
    else
        throw( ErrorUnimplemented( ERROR_LOCATION, "Couldn't convert unknown matrix type" ) );

    return( *this );
}


CRowMatrix::~CRowMatrix()
{
    free( _ptr );
    free( _col );
    free( _val );
}


void CRowMatrix::resize( int n, int m )
{
    free( _col );
    free( _val );
    _col = NULL;
    _val = NULL;

    if( _n != n ) {
	int *tmp;
	if( !(tmp = (int *)realloc( _ptr, (n+1)*sizeof(int) )) ) {
	    free( _ptr );
	    _ptr = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	_ptr = tmp;
    }

    _n     = n;
    _m     = m;
    _nz    = 0;
    _asize = 0;
    memset( _ptr, 0, (_n+1)*sizeof(int) );
}


void CRowMatrix::clear( void )
{
    _nz    = 0;
    _asize = 0;
    memset( _ptr, 0, (_n+1)*sizeof(int) );

    free( _col );
    free( _val );
    _col = NULL;
    _val = NULL;
}


inline void CRowMatrix::clear_no_check( int i, int j )
{
    int a;

    /* Search for the element */
    for( a = _ptr[i]; a < _ptr[i+1]; a++ )
	if( _col[a] == j )
	    break;
    /* Do nothing if the element does not exist */
    if( a == _ptr[i+1] )
	return;

    /* Move data */
    int movesize = _nz-a-1;
    memmove( &_val[a], &_val[a+1], movesize*sizeof(double) );
    memmove( &_col[a], &_col[a+1], movesize*sizeof(int) );

    /* Update pointers */
    _nz--;
    for( a = i+1; a < _n+1; a++ )
	_ptr[a]--;
}


void CRowMatrix::clear_check( int i, int j )
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    clear_no_check( i, j );
}


void CRowMatrix::reserve( int size )
{
    if( size > _asize ) {
	_asize = size;
	reallocate();
    }
}


void CRowMatrix::set_nz( int nz )
{
    if( nz > _asize ) {
	_asize = nz;
	reallocate();
    }
    _nz = nz;
}


void CRowMatrix::merge( CRowMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._asize;
    _ptr   = mat._ptr;
    _col   = mat._col;
    _val   = mat._val;

    mat._n     = 0;
    mat._m     = 0;
    mat._nz    = 0;
    mat._asize = 0;
    mat.allocate();
}


void CRowMatrix::order_ascending( void )
{
    /* Sort each row. */
    for( int i = 0; i < _n; i++ )
	insertion_sort_iv( _col, _val, _ptr[i], _ptr[i+1] );
}


bool CRowMatrix::check_ascending( void ) const
{
    /* Check each row. */
    for( int i = 0; i < _n; i++ ) {
	int prev;
	for( int j = _ptr[i]; j < _ptr[i+1]; j++ ) {
	    if( j > _ptr[i] && prev > _col[j] )
		return( false );
	    prev = _col[j];
	}
    }

    return( true );
}


inline double CRowMatrix::get_no_check( int i, int j ) const
{
    for( int a = _ptr[i]; a < _ptr[i+1]; a++ )
	if( _col[a] == j )
	    return( _val[a] );
    return( 0.0 );
}


inline double &CRowMatrix::set_no_check( int i, int j )
{
    /* Use existing element if it exists */
    for( int a = _ptr[i]; a < _ptr[i+1]; a++ )
    if( _col[a] == j )
        return( _val[a] );

    /* Reserve new space if necessary */
    if( _nz+1 > _asize )
	reserve( _asize+_m );

    /* Move existing data */
    int movesize = _nz-_ptr[i+1];
    memmove( &_val[_ptr[i+1]+1], &_val[_ptr[i+1]], movesize*sizeof(double) );
    memmove( &_col[_ptr[i+1]+1], &_col[_ptr[i+1]], movesize*sizeof(int) );

    /* Set new data */
    _val[_ptr[i+1]] = 0.0;
    _col[_ptr[i+1]] = j;

    /* Update pointers */
    _nz++;
    for( int a = i+1; a < _n+1; a++ )
	_ptr[a]++;

    return( _val[_ptr[i+1]-1] );
}


double CRowMatrix::get_check( int i, int j ) const
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    return( get_no_check( i, j ) );
}


double &CRowMatrix::set_check( int i, int j )
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    return( set_no_check( i, j ) );
}


void CRowMatrix::set_row( int i, int N, const int *col, const double *val )
{
    if( i >= _n )
	throw( ErrorRange( ERROR_LOCATION, i, _n ) );
    else if( N > _m )
	throw( Error( ERROR_LOCATION, "too many elements" ) );

    /* Reserve new space if necessary */
    int oldsize = _ptr[i+1] - _ptr[i];
    if( _nz+N-oldsize > _asize )
	reserve( _asize+_m );

    /* Move existing data */
    int offset = N-oldsize;
    int movesize = _nz-_ptr[i+1];
    memmove( &_val[_ptr[i+1]+offset], &_val[_ptr[i+1]], movesize*sizeof(double) );
    memmove( &_col[_ptr[i+1]+offset], &_col[_ptr[i+1]], movesize*sizeof(int) );

    /* Set new data */
    int err = 0;
    for( int a = 0; a < N; a++ ) {
	for( int b = a+1; b < N; b++ ) {
	    if( col[a] == col[b] )
		err = 2;
	}
	_val[_ptr[i]+a] = val[a];
	if( (_col[_ptr[i]+a] = col[a]) >= _m )
	    err = 1;
    }

    /* Update pointers */
    _nz += offset;
    for( int a = i+1; a < _n+1; a++ )
	_ptr[a] += offset;

    if( err == 1 )
	throw( Error( ERROR_LOCATION, "column index out of range" ) );
    else if( err == 2 )
	throw( Error( ERROR_LOCATION, "repeated column index" ) );
}


void CRowMatrix::construct_add( int i, int j, double val )
{
    /* Reserve new space if necessary */
    if( _nz+1 > _asize )
	reserve( _asize+_m );

    /* Set new data */
    _val[_ptr[i+1]] = val;
    _col[_ptr[i+1]] = j;

    /* Update the number of nonzeroes and only the next pointer */
    _nz++;
    _ptr[i+1]++;
    if( i+1 < _n )
	_ptr[i+2] = _ptr[i+1];
}


void CRowMatrix::debug_print( std::ostream &os ) const
{
    os << "n     = " << _n << "\n";
    os << "m     = " << _m << "\n";
    os << "nz    = " << _nz << "\n";
    os << "asize = " << _asize << "\n";

    os << "ptr[] = {";
    for( int i = 0; i < _n; i++ )
	os << _ptr[i] << ", ";
    os << _ptr[_n] << "}\n";

    os << "col[] = {";
    if( _nz <= 0 ) {
	os << "}\n";
    } else {
	for( int i = 0; i < _nz-1; i++ )
	    os << _col[i] << ", ";
	os << _col[_nz-1] << "}\n";
    }

    os << "val[] = {";
    if( _nz <= 0 ) {
	os << "}\n";
    } else {
	for( int i = 0; i < _nz-1; i++ )
	    os << _val[i] << ", ";
	os << _val[_nz-1] << "}\n";
    }
}



/* ************************************** *
 * Matrix-Vector operations               *
 * ************************************** */


void CRowMatrix::multiply_by_vector( Vector &x, const Vector &b ) const
{
    if( b.size() != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    x.resize( _n );
    x.clear();

    double sum;
    for( int i = 0; i < _n; i++ ) {
	sum = 0;
	int end = _ptr[i+1];
	for( int a = _ptr[i]; a < end; a++ )
	    sum += _val[a] * b._val[_col[a]];
	x._val[i] = sum;
    }
}


void CRowMatrix::lower_unit_solve( Vector &x, const Vector &b ) const
{
    // Make checks
    if( _n != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix not squrare" ) );
    if( b.size() != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    x = b;

    for( int i = 0; i < _n; i++ ) {
	for( int j = _ptr[i]; j < _ptr[i+1]; j++ )
	    x[i] -= _val[j]*x[_col[j]];
    }
}


void CRowMatrix::upper_diag_solve( Vector &x, const Vector &b ) const
{
    // Make checks
    if( _n != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix not square" ) );
    if( b.size() != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    x = b;

    int32_t i, j;
    for( i = _n-1; i >= 0; i-- ) {
	j = _ptr[i+1]-1;
	for( ; j > (int32_t)_ptr[i]; j-- )
	    x[i] -= _val[j] * x[_col[j]];
	x[i] /= _val[j];
    }
}


void CRowMatrix::LU_solve( Vector &x, const Vector &b ) const
{
    // Make checks
    if( _n != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix not squrare" ) );
    if( b.size() != _m )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    x = b;

    // Solve L*x=b
    for( int i = 0; i < _n; i++ ) {
	for( int j = _ptr[i]; j < _ptr[i+1]; j++ ) {
	    if( _col[j] >= i ) 
		break;
	    x[i] -= _val[j]*x[_col[j]];
	}
    }

    // Solve U*x=x
    for( int i = _n-1; i >= 0; i-- ) {
	int j = _ptr[i+1]-1;
	for( ; j > _ptr[i]; j-- ) {
	    if( _col[j] == i )
		break;
	    x[i] -= _val[j] * x[_col[j]];
	}
	x[i] /= _val[j];
    }
}

