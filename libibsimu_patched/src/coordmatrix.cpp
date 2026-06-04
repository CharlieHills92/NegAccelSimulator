/*! \file coordmatrix.cpp
 *  \brief Sparse coordinate-based sparse matrices.
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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
#include "coordmatrix.hpp"
#include "ccolmatrix.hpp"
#include "crowmatrix.hpp"
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iomanip>



inline void CoordMatrix::allocate( void )
{
    if( _asize == 0 ) {
	_row = _col = NULL;
	_val = NULL;
    } else {
	if( !(_row = (int *)malloc( _asize*sizeof(int) )) ) {
	    _row = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	if( !(_col = (int *)malloc( _asize*sizeof(int) )) ) {
	    free( _row );
	    _row = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	if( !(_val = (double *)malloc( _asize*sizeof(double) )) ) {
	    free( _row );
	    free( _col );
	    _row = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
    }
}


inline void CoordMatrix::reallocate( void )
{
    if( _asize == 0 ) {
	free( _row );	
	free( _col );
	free( _val );
	_row = _col = NULL;
	_val = NULL;
    } else {
	int *tmp;
	if( !(tmp = (int *)realloc( _row, _asize*sizeof(int) )) ) {
	    free( _row );
	    free( _col );
	    free( _val );
	    _row = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	_row = tmp;
	if( !(tmp = (int *)realloc( _col, _asize*sizeof(int) )) ) {
	    free( _row );
	    free( _col );
	    free( _val );
	    _row = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	_col = tmp;
	double *tmp2;
	if( !(tmp2 = (double *)realloc( _val, _asize*sizeof(double) )) ) {
	    free( _row );
	    free( _col );
	    free( _val );
	    _row = _col = NULL;
	    _val = NULL;
	    _n = _m = _nz = _asize = 0;
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	_val = tmp2;
    }
}


CoordMatrix::CoordMatrix( int n, int m )
{
    _n     = n;
    _m     = m;
    _nz    = 0;
    _asize = 0;
    _row   = NULL;
    _col   = NULL;
    _val   = NULL;    
}


CoordMatrix::CoordMatrix( int n, int m, int nz, 
			  const int *row, const int *col, const int *val )
{
    _n     = n;
    _m     = m;
    _nz    = nz;
    _asize = nz;
    allocate();

    memcpy( _row, row, _nz*sizeof(int) );
    memcpy( _col, col, _nz*sizeof(int) );
    memcpy( _val, val, _nz*sizeof(double) );
}


void CoordMatrix::build( const class CoordMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    memcpy( _row, mat._row, _nz*sizeof(int) );
    memcpy( _col, mat._col, _nz*sizeof(int) );
    memcpy( _val, mat._val, _nz*sizeof(double) );
}


CoordMatrix::CoordMatrix( const CoordMatrix &mat )
    : _row(NULL), _col(NULL), _val(NULL)
{
    build( mat );
}

CoordMatrix &CoordMatrix::operator=( const CoordMatrix &mat )
{
    build( mat );
    return( *this );
}


void CoordMatrix::build( const class CRowMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    int i, j;
    for( i = j = 0; i < _n; i++ ) {
	int e = mat._ptr[i+1];
	for( ; j < e; j++ ) {
	    _row[j] = i;
	    _col[j] = mat._col[j];
	    _val[j] = mat._val[j];
	}
    }
}


CoordMatrix::CoordMatrix( const class CRowMatrix &mat )
    : _row(NULL), _col(NULL), _val(NULL)
{
    build( mat );
}

CoordMatrix &CoordMatrix::operator=( const CRowMatrix &mat )
{
    build( mat );
    return( *this );
}


void CoordMatrix::build( const class CColMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._nz;
    reallocate();

    int i, j;
    for( i = j = 0; j < _m; j++ ) {
	int e = mat._ptr[j+1];
	for( ; i < e; i++ ) {
	    _row[i] = mat._row[i];
	    _col[i] = j;
	    _val[i] = mat._val[i];
	}
    }
}


CoordMatrix::CoordMatrix( const class CColMatrix &mat )
    : _row(NULL), _col(NULL), _val(NULL)
{
    build( mat );
}


CoordMatrix &CoordMatrix::operator=( const CColMatrix &mat )
{
    build( mat );
    return( *this );
}


CoordMatrix::CoordMatrix( const class Matrix &mat )
    : _row(NULL), _col(NULL), _val(NULL)
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


CoordMatrix &CoordMatrix::operator=( const Matrix &mat )
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


CoordMatrix::~CoordMatrix()
{
    free( _row );
    free( _col );
    free( _val );
}


void CoordMatrix::resize( int n, int m )
{
    _n     = n;
    _m     = m;
    _nz    = 0;
    _asize = 0;
    free( _row );
    free( _col );
    free( _val );
    _row = NULL;
    _col = NULL;
    _val = NULL;
}


void CoordMatrix::merge( CoordMatrix &mat )
{
    _n     = mat._n;
    _m     = mat._m;
    _nz    = mat._nz;
    _asize = mat._asize;
    free( _row );
    free( _col );
    free( _val );
    _row = mat._row;
    _col = mat._col;
    _val = mat._val;

    mat._n     = 0;
    mat._m     = 0;
    mat._nz    = 0;
    mat._asize = 0;
    mat._row   = NULL;
    mat._col   = NULL;
    mat._val   = NULL;
}


void CoordMatrix::clear( void )
{
    _nz    = 0;
    _asize = 0;
    free( _row );
    free( _col );
    free( _val );
    _row = NULL;
    _col = NULL;
    _val = NULL;
}


inline void CoordMatrix::clear_no_check( int i, int j )
{

}


void CoordMatrix::clear_check( int i, int j )
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    clear_no_check( i, j );
}


void CoordMatrix::reserve( int size )
{
    if( size > _asize ) {
	_asize = size;
	reallocate();
    }
}


void CoordMatrix::set_nz( int nz )
{
    if( nz > _asize ) {
	_asize = nz;
	reallocate();
    }
    _nz = nz;
}


void CoordMatrix::order_ascending_row_column( void )
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


void CoordMatrix::order_ascending_column_row( void )
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


void CoordMatrix::debug_print( std::ostream &os ) const
{
    os << "n     = " << _n << "\n";
    os << "m     = " << _m << "\n";
    os << "nz    = " << _nz << "\n";
    os << "asize = " << _asize << "\n";

    os << "row[] = {";
    if( _nz <= 0 ) {
	os << "}\n";
    } else {
	for( int i = 0; i < _nz-1; i++ )
	    os << _row[i] << ", ";
	os << _row[_nz-1] << "}\n";
    }

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


inline double CoordMatrix::get_no_check( int i, int j ) const
{
    for( int a = 0; a < _nz; a++ )
	if( _row[a] == i && _col[a] == j  )
	    return( _val[a] );
    return( 0.0 );
}


inline double &CoordMatrix::set_no_check( int i, int j )
{
    /* Use existing element if it exists */
    for( int a = 0; a < _nz; a++ )
	if( _row[a] == i && _col[a] == j  )
	    return( _val[a] );

    /* Reserve new space if necessary */
    if( _nz >= _asize )
	reserve( _asize+_n );

    /* Set new data */
    _row[_nz] = i;
    _col[_nz] = j;
    _val[_nz] = 0.0;
    _nz += 1;

    return( _val[_nz-1] );
}


double CoordMatrix::get_check( int i, int j ) const
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    return( get_no_check( i, j ) );
}


double &CoordMatrix::set_check( int i, int j )
{
    if( i >= _n || j >= _m )
	throw( ErrorRange( ERROR_LOCATION, i, _n, j, _m ) );

    return( set_no_check( i, j ) );
}


void CoordMatrix::set_no_duplicate_check( int i, int j, double val )
{
    /* Reserve new space if necessary */
    if( _nz >= _asize )
	reserve( _asize+_n );

    /* Set new data */
    _row[_nz] = i;
    _col[_nz] = j;
    _val[_nz] = val;
    _nz += 1;
}



/* ************************************** *
 * Matrix-Vector operations               *
 * ************************************** */


void CoordMatrix::multiply_by_vector( Vector &x, const Vector &b ) const
{
    if( b.size() != _m )
	throw( ErrorDim( ERROR_LOCATION, "Matrix dimension does not match vector" ) );

    x.resize( _n );
    x.clear();

    for( int i = 0; i < _nz; i++ )
	x._val[_row[i]] += _val[i] * b._val[_col[i]];
}


void CoordMatrix::lower_unit_solve( Vector &x, const Vector &b ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


void CoordMatrix::upper_diag_solve( Vector &x, const Vector &b ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}








