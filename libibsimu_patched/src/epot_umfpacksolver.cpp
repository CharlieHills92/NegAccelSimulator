/*! \file epot_umfpacksolver.cpp
 *  \brief UMFPACK matrix solver for electric potential problem
 */

/* Copyright (c) 2005-2013 Taneli Kalvas. All rights reserved.
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


#include <umfpack.h>
#include <limits>
#include "epot_umfpacksolver.hpp"
#include "ibsimu.hpp"


EpotUMFPACKSolver::EpotUMFPACKSolver( Geometry &geom, 
				      double newton_r_eps, 
				      double newton_step_eps, 
				      uint32_t newton_imax,
				      bool gnewton )
    : EpotMatrixSolver(geom), _numeric(0), _gnewton(gnewton), _newton_r_eps(newton_r_eps), 
      _newton_step_eps(newton_step_eps), _newton_imax(newton_imax)
{

}


EpotUMFPACKSolver::EpotUMFPACKSolver( Geometry &geom, std::istream &s )
    : EpotMatrixSolver(geom,s)
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );

}


void EpotUMFPACKSolver::save( std::ostream &s ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


EpotUMFPACKSolver::~EpotUMFPACKSolver()
{
    if( _numeric )
	umfpack_di_free_numeric( &_numeric );
}


void EpotUMFPACKSolver::set_newton_imax( uint32_t newton_imax ) 
{
    _newton_imax = newton_imax;
}


void EpotUMFPACKSolver::set_newton_residual_eps( double newton_r_eps ) 
{
    _newton_r_eps = newton_r_eps;
}


void EpotUMFPACKSolver::set_newton_step_eps( double newton_step_eps ) 
{
    _newton_step_eps = newton_step_eps;
}


void EpotUMFPACKSolver::set_gnewton( bool enable )
{
    _gnewton = enable;
}


double EpotUMFPACKSolver::get_newton_residual( void ) const
{
    return( _newton_res );
}


double EpotUMFPACKSolver::get_newton_step( void ) const
{
    return( _newton_step );
}


void EpotUMFPACKSolver::reset_problem( void )
{
    if( _numeric )
	umfpack_di_free_numeric( &_numeric );
    _numeric = 0;
    reset_matrix();
}


void EpotUMFPACKSolver::umfpack_error( const std::string func, int status )
{
    if( status == UMFPACK_ERROR_n_nonpositive )
	throw( Error( ERROR_LOCATION, func + ": n less than or equal to zero" ) );
    else if( status == UMFPACK_ERROR_invalid_matrix )
	throw( Error( ERROR_LOCATION, func + ": invalid matrix" ) );
    else if( status == UMFPACK_ERROR_out_of_memory )
	throw( Error( ERROR_LOCATION, func + ": memory allocation error" ) );
    else if( status == UMFPACK_ERROR_argument_missing )
	throw( Error( ERROR_LOCATION, func + ": argument missing" ) );
    else if( status == UMFPACK_ERROR_internal_error )
	throw( Error( ERROR_LOCATION, func + ": internal error" ) );
    else if( status == UMFPACK_WARNING_singular_matrix )
	throw( Error( ERROR_LOCATION, func + ": singular matrix" ) );
    else if( status == UMFPACK_ERROR_invalid_Symbolic_object )
	throw( Error( ERROR_LOCATION, func + ": invalid symbolic object" ) );
    else if( status == UMFPACK_ERROR_invalid_system )
	throw( Error( ERROR_LOCATION, func + ": invalid system" ) );
    else if( status == UMFPACK_ERROR_different_pattern )
	throw( Error( ERROR_LOCATION, func + ": different pattern" ) );
    else if( status == UMFPACK_ERROR_invalid_Numeric_object )
	throw( Error( ERROR_LOCATION, func + ": invalid numeric object" ) );
    else if( status != UMFPACK_OK )
	throw( Error( ERROR_LOCATION, "unknown error in " + func ) );
}


void EpotUMFPACKSolver::umfpack_decompose( const CColMatrix &mat )
{
    int   status;
    void *symbolic;

    // Free old decomposition
    if( _numeric )
	umfpack_di_free_numeric( &_numeric );
    _numeric = 0;

    status = umfpack_di_symbolic( mat.columns(), mat.rows(), 
				  &mat.ptr(0), &mat.row(0), &mat.val(0), 
				  &symbolic, (double *)NULL, (double *)NULL );
    if( status != UMFPACK_OK )
	umfpack_error( "umfpack_di_symbolic", status );

    status = umfpack_di_numeric( &mat.ptr(0), &mat.row(0), &mat.val(0), 
				 symbolic, &_numeric, (double *)NULL, 
				 (double *)NULL );
    if( status != UMFPACK_OK )
	umfpack_error( "umfpack_di_numeric", status );

    // Free symbolic data
    umfpack_di_free_symbolic( &symbolic );
}


void EpotUMFPACKSolver::umfpack_solve( const CColMatrix &mat, const Vector &rhs, Vector &sol, 
				       bool force_decomposition )
{
    int status;

    // Resize solution vector
    sol.resize( rhs.size() );

    // Do decomposition if needed
    if( !_numeric || force_decomposition )
	umfpack_decompose( mat );

    status = umfpack_di_solve( UMFPACK_A, &mat.ptr(0), &mat.row(0), &mat.val(0), 
			       sol.get_data(), rhs.get_data(), _numeric, 
			       (double *)NULL, (double *)NULL );
    if( status != UMFPACK_OK )
	umfpack_error( "umfpack_di_solve", status );
}


void EpotUMFPACKSolver::subsolve( MeshScalarField &epot, const MeshScalarField &scharge )
{
    // Preprocess and set starting guess
    preprocess( epot, scharge );
    Vector X;
    set_initial_guess( epot, X );

    if( linear() ) {

	ibsimu.message(1) << "Using UMFPACK solver\n";
	ibsimu.flush();

	// Fetch matrix form of problem
	const CRowMatrix *A;
	const Vector *B;
	get_vecmat( &A, &B );

	CColMatrix Acol( *A );
	umfpack_solve( Acol, *B, X );

    } else {

	int32_t a;
        const CRowMatrix *J;
        const Vector *R;
        Vector dX;

	_newton_res = 0.0;
	_newton_step = 0.0;
	if( _gnewton ) {

	    ibsimu.message(1) << "Using Globally convergent Newton-Raphson UMFPACK solver ("
			      << "newton_imax = " << _newton_imax
			      << ", newton_r_eps = " << _newton_r_eps
			      << ", newton_step_eps = " << _newton_step_eps
			      << ")\n";
	    ibsimu.message(1) << std::setw(5)  << "Round" << " " 
			      << std::setw(14) << "Step size" << " " 
			      << std::setw(14) << "Step fac" << " " 
			      << std::setw(14) << "Residual" << "\n";
	    ibsimu.flush();
	    
	    // Globally convergent Newton-Raphson
            Vector Xold( X.size() );

            // First jacobian and residual
            get_resjac( &J, &R, X );
            double f = ssqr( *R );

	    for( a = 0; a < (int)_newton_imax; a++ ) {

		// Calculate dX = J^{-1}*R
                dX.clear();
		CColMatrix Jcol( *J );
		umfpack_solve( Jcol, *R, dX, true );

		// Search for acceptable step for which residual decreases
                double t = 2.0;
                double fold = f;
                Xold = X;
                while( f >= fold ) {

                    t *= 0.5;
                    X = Xold - t*dX;
                    get_resjac( &J, &R, X );
                    f = ssqr( *R );
                    if( t <= std::numeric_limits<double>::epsilon() )
                        break;
                }
		
		// Check for convergence
		_newton_res = max_abs( *R );
		_newton_step = max_abs( dX );
		
		ibsimu.message(1) << std::setw(5) << a << " " 
				  << std::setw(14) << _newton_step << " " 
				  << std::setw(14) << t << " " 
				  << std::setw(14) << _newton_res << "\n";
		ibsimu.flush();

		if( _newton_res < _newton_r_eps || _newton_step < _newton_step_eps )
		    break;
	    }

	} else {

	    ibsimu.message(1) << "Using Newton-Raphson UMFPACK solver("
			      << "newton_imax = " << _newton_imax
			      << ", newton_r_eps = " << _newton_r_eps
			      << ", newton_step_eps = " << _newton_step_eps
			      << ")\n";
	    ibsimu.message(1) << std::setw(5)  << "Round" << " " 
			      << std::setw(14) << "Step size" << " " 
			      << std::setw(14) << "Residual" << "\n";
	    ibsimu.flush();
	    
	    for( a = 0; a < (int)_newton_imax; a++ ) {
		// Calculate dX = J^{-1}*R
		get_resjac( &J, &R, X );
		CColMatrix Jcol( *J );
		dX.clear();
		umfpack_solve( Jcol, *R, dX, true );
		
		// Take step
		X -= dX;
		
		// Check for convergence
		_newton_res = max_abs( *R );
		_newton_step = max_abs( dX );
		
		ibsimu.message(1) << std::setw(5) << a << " " 
				  << std::setw(14) << _newton_step << " " 
				  << std::setw(14) << _newton_res << "\n";
		ibsimu.flush();

		if( _newton_res < _newton_r_eps || _newton_step < _newton_step_eps )
		    break;
	    }
	}

	if( _newton_res < _newton_r_eps || _newton_step < _newton_step_eps )
	    ibsimu.message(1) << "  Newton-Raphson converged\n";
	else
	    ibsimu.message(1) << "  Maximum number of Newton-Raphson iterations\n";
    }

    // Postprocess and set solution
    set_solution( epot, X );
    postprocess();
}

    
void EpotUMFPACKSolver::debug_print( std::ostream &os ) const
{
    EpotMatrixSolver::debug_print( os );
    os << "**EpotUMFPACKSolver\n";

    os << "newton_Reps = " << _newton_r_eps << "\n";
    os << "newton_dXeps = " << _newton_step_eps << "\n";
    os << "newton_imax = " << _newton_imax << "\n";
}
