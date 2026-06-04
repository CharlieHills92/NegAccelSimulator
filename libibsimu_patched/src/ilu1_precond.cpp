/*! \file ilu1_precond.cpp
 *  \brief ILU1 preconditioner for sparse matrices
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


#include "ilu1_precond.hpp"
#include "sort.hpp"
#include "error.hpp"
#include <stdlib.h>


void ILU1_Precond::add_element( int *ptr, int **col, int n, int &nz, int &asize, int c )
{
    if( nz == asize ) {
	// Reallocate memory
	asize += n;
	int *ncol = (int *)realloc( (void *)(*col), asize*sizeof(int) );
	if( ncol == NULL ) {
	    free( ptr );
	    free( *col );
	    throw( ErrorNoMem( ERROR_LOCATION ) );
	}
	(*col) = ncol;
    }
    (*col)[nz] = c;
    nz++;
}


std::string ILU1_Precond::typestring( void ) const
{
    return( "ILU1" );
}
    

ILU1_Precond::ILU1_Precond()
{
    _LU = NULL;
}


ILU1_Precond::~ILU1_Precond()
{
    delete _LU;
}


void ILU1_Precond::prepare( const CRowMatrix &A )
{
    if( _LU != NULL )
	delete _LU;

    // Make checks
    if( A.columns() != A.rows() )
	throw( ErrorDim( ERROR_LOCATION, "matrix not square" ) );
    if( !A.check_ascending() )
	throw( ErrorDim( ERROR_LOCATION, "matrix not in ascending order" ) );

    // Initialize LU matrix needed for ILU1
    int n = A.rows();
    int nz = 0;
    int asize = A.nz_elements();
    int *ptr = (int *)malloc( (n+1)*sizeof(int) );
    int *col = (int *)malloc( asize*sizeof(int) );

    // Go through IKJ Gaussian elimination algorithm, 
    // only using original A elemets as input (ILU1)
    for( int i = 0; i < n; i++ ) {

	// Copy elements of row i
	ptr[i] = nz;
	for( int cp = A.ptr(i); cp < A.ptr(i+1); cp++ )
	    add_element( ptr, &col, n, nz, asize, A.col(cp) );

	// Go through columns 0 <= k < i
	for( int kp = A.ptr(i); kp < A.ptr(i+1); kp++ ) {
	    int k = A.col(kp);
	    if( k >= i )
		break;

	    // Find elements (k,j), for j > k
	    for( int jp = A.ptr(k); jp < A.ptr(k+1); jp++ ) {
		int j = A.col(jp);
		if( j <= k )
		    continue;
		
		// Find if (i,j) exists in UL
		int mp;
		for( mp = ptr[i]+kp-A.ptr(i)+1; mp < nz; mp++ )
		    if( col[mp] == j )
			break; // Exists
		if( mp == nz ) {
		    // Does not exist, add new element to LU
		    add_element( ptr, &col, n, nz, asize, j );
		}
	    }
	}
    }

    // End processing, mark number of elements
    ptr[n] = nz;

    // Sort rows to ascending column order
    for( int i = 0; i < n; i++ )
	insertion_sort_i( col, ptr[i], ptr[i+1] );

    // Create LU matrix with zero value elements
    double *val = (double *)calloc( nz, sizeof(double) );
    _LU = new CRowMatrix( n, n, nz, nz, ptr, col, val );
}


void ILU1_Precond::construct( const CRowMatrix &A )
{
    if( !A.check_ascending() )
	throw( ErrorDim( ERROR_LOCATION, "matrix not in ascending order" ) );

    // Set LU matrix element values from A
    // only using original A elements as input (ILU1)
    for( int i = 0; i < _LU->rows(); i++ ) {
	int k = A.ptr(i);
	int kmax = A.ptr(i+1);
	for( int j = _LU->ptr(i); j < _LU->ptr(i+1); j++ ) {
	    if( k == kmax ) {
		_LU->val(j) = 0.0;
		continue;
	    }
	    if( _LU->col(j) == A.col(k) ) {
		_LU->val(j) = A.val(k);
		k++;
	    } else
		_LU->val(j) = 0.0;
	}
    }

    // Go through IKJ Gaussian elimination algorithm 
    for( int i = 1; i < _LU->rows(); i++ ) {

	// Go through columns 0 <= k < i
	for( int kp = _LU->ptr(i); kp < _LU->ptr(i+1); kp++ ) {
	    int k = _LU->col(kp);
	    if( k >= i )
		break;

	    // Find diagonal (k,k)
	    int jp;
	    for( jp = _LU->ptr(k); jp < _LU->ptr(k+1); jp++ )
		if( _LU->col(jp) == k )  break;

	    double mult = _LU->val(kp)/_LU->val(jp);
	    _LU->val(kp) = mult;

	    // Process elements (k,j), for j > k
	    for( jp++; jp < _LU->ptr(k+1); jp++ ) {
		int j = _LU->col(jp);
		
		// Find (i,j)
		for( int mp = kp+1; mp < _LU->ptr(i+1); mp++ ) {
		    if( _LU->col(mp) == j ) {
			_LU->val(mp) -= mult*_LU->val(jp);
			break;
		    }
		}
	    }
	}
    }

    // Done factorizing
}


void ILU1_Precond::clear( void )
{
    if( _LU != NULL )
	delete _LU;
    _LU = NULL;
}


bool ILU1_Precond::is_prepared( void ) const
{
    return( _LU );
}


const CRowMatrix *ILU1_Precond::get_matrix( void ) const
{
    return( _LU );
}


void ILU1_Precond::debug_print( std::ostream &os ) const
{
    if( _LU == NULL )
	os << "LU = NULL\n";
    else
	os << "LU = \n" << *_LU << "\n\n";
}


void ILU1_Precond::solve( Vector &x, const Vector &b ) const
{
    _LU->LU_solve( x, b );
}

