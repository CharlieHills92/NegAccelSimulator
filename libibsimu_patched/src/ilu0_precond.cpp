/*! \file ilu0_precond.cpp
 *  \brief ILU0 preconditioner for sparse matrices
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

#include "ilu0_precond.hpp"
#include "error.hpp"
#include <stdlib.h>


std::string ILU0_Precond::typestring( void ) const
{
    return( "ILU0" );
}
    

ILU0_Precond::ILU0_Precond()
{
    _LU = NULL;
}


ILU0_Precond::~ILU0_Precond()
{
    delete _LU;
}


void ILU0_Precond::prepare( const CRowMatrix &A )
{
    if( _LU != NULL )
	delete _LU;

    _LU = new CRowMatrix( A );
}


void ILU0_Precond::construct( const CRowMatrix &A )
{
    // Copy A contents
    memcpy( &_LU->val(0), &A.val(0), A.nz_elements()*sizeof(double) );

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
		
		// Find (i,j), do nothing if not found (ILU0)
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


void ILU0_Precond::clear( void )
{
    if( _LU != NULL )
	delete _LU;
    _LU = NULL;
}


bool ILU0_Precond::is_prepared( void ) const
{
    return( _LU );
}


const CRowMatrix *ILU0_Precond::get_matrix( void ) const
{
    return( _LU );
}


void ILU0_Precond::debug_print( std::ostream &os ) const
{
    if( _LU == NULL )
	os << "LU = NULL\n";
    else
	os << "LU = \n" << *_LU << "\n\n";
}


void ILU0_Precond::solve( Vector &x, const Vector &b ) const
{
    _LU->LU_solve( x, b );
}

