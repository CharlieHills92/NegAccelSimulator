/*! \file epot_gssolver.cpp
 *  \brief Gauss-Seidel solver for electric potential problem
 */

/* Copyright (c) 2011-2013 Taneli Kalvas. All rights reserved.
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


#include "epot_gssolver.hpp"
#include "ibsimu.hpp"
#include "constants.hpp"
#include "compmath.hpp"
#include "statusprint.hpp"


// Minimum step factor for globally convergent Newton-Raphson
#define FMIN 1.0e-6


EpotGSSolver::EpotGSSolver( Geometry &geom )
    : EpotSolver( geom ), _epot(NULL), _rhs(NULL), _iter(0), _imax(100000), 
      _eps(1.0e-4), _step(0.0), _err(0.0), _w(1.66),
      _local_Ulim_fac(10.0), _local_Ulim(0.0), _local_imax(1), _local_eps(1.0e-6)
{
    
}


EpotGSSolver::EpotGSSolver( Geometry &geom, std::istream &s )
    : EpotSolver(geom,s)
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


void EpotGSSolver::reset_problem( void )
{
    // Do nothing
}


void EpotGSSolver::set_eps( double eps )
{
    _eps = eps;
}


double EpotGSSolver::get_potential_change_norm( void ) const
{
    return( _step );
}


double EpotGSSolver::get_error_estimate( void ) const
{
    return( _err );
}


uint32_t EpotGSSolver::get_iter( void ) const
{
    return( _iter );
}


void EpotGSSolver::set_imax( uint32_t imax )
{
    _imax = imax;
}


void EpotGSSolver::set_w( double w )
{
    _w = w;
}


void EpotGSSolver::set_plasma_solver_parameters( double Ulim_fac, uint32_t imax, double eps )
{
    _local_Ulim_fac = Ulim_fac;
    _local_imax = imax;
    _local_eps = eps;
}


double EpotGSSolver::solve_shield_potential( double epf, double cof, double rhs, double p ) const
{
    // If potential large enough, no plasma calculation needed
    if( fabs(p) > _local_Ulim )
	return( p + (epf - cof*p - rhs) / cof );

    // Globally convergent Newton-Raphson solution for node potential
    double rhst, drhst; 
    shield_newton( rhst, drhst, p );
    double F = epf - cof*p - rhs*rhst;
    double Fnew, pnew;
    for( uint32_t q = 0; q < _local_imax; q++ ) {
	double deltap = F / ( cof + rhs*drhst );
	double fac = 1.0;
	while( fac >= FMIN ) {
	    pnew = p + fac*deltap;
	    shield_newton( rhst, drhst, pnew );
	    Fnew = epf - cof*pnew - rhs*rhst;
	    if( fabs(Fnew) < fabs(F) )
		break;
	    fac *= 0.5;
	}
	if( fac < FMIN || p == pnew )
	    break;
	F = Fnew;
	p = pnew;
	if( fabs(deltap) < _local_eps )
	    break;
    }
    return( p );
}


double EpotGSSolver::solve_nsimp_potential( double epf, double cof, double rhs, double p ) const
{
    // If potential large enough, no plasma calculation needed
    if( p > _local_Ulim )
	return( p + (epf - cof*p - rhs) / cof );

    // Globally convergent Newton-Raphson solution for node potential
    double rhst, drhst; 
    nsimp_newton( rhst, drhst, p );
    double F = epf - cof*p - rhs - rhst;
    double Fnew, pnew;
    for( uint32_t q = 0; q < _local_imax; q++ ) {
	double deltap = F / ( cof + drhst );
	double fac = 1.0;
	while( fac >= FMIN ) {
	    pnew = p + fac*deltap;
	    nsimp_newton( rhst, drhst, pnew );
	    Fnew = epf - cof*pnew - rhs - rhst;
	    if( fabs(Fnew) < fabs(F) )
		break;
	    fac *= 0.5;
	}
	if( fac < FMIN || p == pnew )
	    break;
	F = Fnew;
	p = pnew;
	if( fabs(deltap) < _local_eps )
	    break;
    }
    return( p );
}


double EpotGSSolver::solve_pexp_potential( double epf, double cof, double rhs, double p ) const
{
    // If potential small enough, no plasma calculation needed
    if( p < _local_Ulim )
	return( p + (epf - cof*p - rhs) / cof );

    // Globally convergent Newton-Raphson solution for node potential
    double rhst, drhst;
    pexp_newton( rhst, drhst, p );
    double F = epf - cof*p - rhs - rhst;
    double Fnew, pnew;
    for( uint32_t q = 0; q < _local_imax; q++ ) {
	double deltap = F / ( cof + drhst );
	double fac = 1.0;
	while( fac >= FMIN ) {
	    pnew = p + fac*deltap;
	    pexp_newton( rhst, drhst, pnew );
	    Fnew = epf - cof*pnew - rhs - rhst;
	    if( fabs(Fnew) < fabs(F) )
		break;
	    fac *= 0.5;
	}
	if( fac < FMIN )
	    break;
	F = Fnew;
	p = pnew;
	if( fabs(deltap) < _local_eps )
	    break;
    }
    return( p );
}


/* *****************************************************************************
 * 3D
 */

double EpotGSSolver::gs_process_near_solid_3d( const uint8_t *nearsolid_ptr, uint32_t a,
					       uint32_t dj, uint32_t dk, uint8_t bindex ) const
{
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    double cof = 0.0;
    double epf = 0.0;

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis
    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0/(beta*beta);
	epf += 2.0/(beta*beta)*(*_epot)(a+1);
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0/(alpha*alpha);
	epf += 2.0/(alpha*alpha)*(*_epot)(a-1);
    } else {
	cof += 2.0/(alpha*beta);
	epf += 2.0/(alpha+beta)*( (*_epot)(a-1)/alpha + (*_epot)(a+1)/beta );
    }

    // Ymin direction
    alpha = 1.0;
    if( sflag & 0x04 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Ymax direction
    beta = 1.0;
    if( sflag & 0x08 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Y axis
    if( bindex & EPOT_SOLVER_BYMIN ) {
	cof += 2.0/(beta*beta);
	epf += 2.0/(beta*beta)*(*_epot)(a+dj);
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	cof += 2.0/(alpha*alpha);
	epf += 2.0/(alpha*alpha)*(*_epot)(a-dj);
    } else {
	cof += 2.0/(alpha*beta);
	epf += 2.0/(alpha+beta)*( (*_epot)(a-dj)/alpha + (*_epot)(a+dj)/beta );
    }

    // Zmin direction
    alpha = 1.0;
    if( sflag & 0x10 ) {
	alpha = *ptr/255.0;
	ptr++;
    }
    
    // Zmax direction
    beta = 1.0;
    if( sflag & 0x20 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Z axis
    if( bindex & EPOT_SOLVER_BZMIN ) {
	cof += 2.0/(beta*beta);
	epf += 2.0/(beta*beta)*(*_epot)(a+dk);
    } else if( bindex & EPOT_SOLVER_BZMAX ) {
	cof += 2.0/(alpha*alpha);
	epf += 2.0/(alpha*alpha)*(*_epot)(a-dk);
    } else {
	cof += 2.0/(alpha*beta);
	epf += 2.0/(alpha+beta)*( (*_epot)(a-dk)/alpha + (*_epot)(a+dk)/beta );
    }

    if( _plasma == PLASMA_PEXP ) {
	return( solve_pexp_potential( epf, cof, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	return( solve_nsimp_potential( epf, cof, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	return( solve_shield_potential( epf, cof, (*_rhs)(a), (*_epot)(a) ) );
    }

    return( (1.0/cof) * ( epf - (*_rhs)(a) ) );
}


double EpotGSSolver::gs_process_pure_vacuum_3d( uint32_t a, uint32_t dj, uint32_t dk ) const
{
    if( _plasma == PLASMA_PEXP ) {
	double epf = (*_epot)(a+1) + (*_epot)(a-1) 
	    + (*_epot)(a+dj) + (*_epot)(a-dj) 
	    + (*_epot)(a+dk) + (*_epot)(a-dk);
	return( solve_pexp_potential( epf, 6.0, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double epf = (*_epot)(a+1) + (*_epot)(a-1) 
	    + (*_epot)(a+dj) + (*_epot)(a-dj) 
	    + (*_epot)(a+dk) + (*_epot)(a-dk);
	return( solve_nsimp_potential( epf, 6.0, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	double epf = (*_epot)(a+1) + (*_epot)(a-1) 
	    + (*_epot)(a+dj) + (*_epot)(a-dj) 
	    + (*_epot)(a+dk) + (*_epot)(a-dk);
	return( solve_shield_potential( epf, 6.0, (*_rhs)(a), (*_epot)(a) ) );
    }

    return( (1.0/6.0) * ( (*_epot)(a+1)  + (*_epot)(a-1) +
			  (*_epot)(a+dj) + (*_epot)(a-dj) + 
			  (*_epot)(a+dk) + (*_epot)(a-dk) - (*_rhs)(a) ) );
}


double EpotGSSolver::gs_process_neumann_3d( uint32_t a, uint32_t dj, uint32_t dk, 
					    uint8_t bindex ) const
{
    double cof = 0.0;
    double epf = 0.0;

    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0;
	epf += 2.0*(*_epot)(a+1);
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0;
	epf += 2.0*(*_epot)(a-1);
    } else {
	cof += 2.0;
	epf += (*_epot)(a+1) + (*_epot)(a-1);
    }

    if( bindex & EPOT_SOLVER_BYMIN ) {
	cof += 2.0;
	epf += 2.0*(*_epot)(a+dj);
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	cof += 2.0;
	epf += 2.0*(*_epot)(a-dj);
    } else {
	cof += 2.0;
	epf += (*_epot)(a+dj) + (*_epot)(a-dj);
    }

    if( bindex & EPOT_SOLVER_BZMIN ) {
	cof += 2.0;
	epf += 2.0*(*_epot)(a+dk);
    } else if( bindex & EPOT_SOLVER_BZMAX ) {
	cof += 2.0;
	epf += 2.0*(*_epot)(a-dk);
    } else {
	cof += 2.0;
	epf += (*_epot)(a+dk) + (*_epot)(a-dk);
    }

    if( _plasma == PLASMA_PEXP ) {
	return( solve_pexp_potential( epf, cof, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	return( solve_nsimp_potential( epf, cof, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	return( solve_shield_potential( epf, cof, (*_rhs)(a), (*_epot)(a) ) );
    }

    return( (1.0/cof) * ( epf - (*_rhs)(a) ) );
}


double EpotGSSolver::gs_loop_3d( void ) const
{
    // Go through all nodes
    const double w2 = 1.0-_w;
    const uint32_t dj = _geom.size(0);
    const uint32_t dk = _geom.size(0)*_geom.size(1);

    double res = 0.0;
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    uint32_t a = k*dk+j*dj;
	    for( uint32_t i = 0; i < _geom.size(0); i++, a++ ) {

		double Vold = (*_epot)(a);
		double Vnew;
		uint32_t mesh = _geom.mesh(a);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		uint8_t bindex = boundary_index(i,j,k);

		if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
		    Vnew = gs_process_near_solid_3d( nearsolid_ptr, a, dj, dk, bindex );
		} else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		    Vnew = gs_process_pure_vacuum_3d( a, dj, dk );
		} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		    Vnew = gs_process_neumann_3d( a, dj, dk, bindex );
		} else {
		    // Dirichlet or otherwise fixed
		    continue;
		}
		Vnew = _w*Vnew + w2*Vold;
		(*_epot)(a) = Vnew;
		double dx = Vnew - Vold;
		res += dx*dx;
		if( comp_isinf(dx) ) {
		    throw( Error( ERROR_LOCATION, "Potential inf at location = " + to_string(i) + 
				  ", " + to_string(j) + ", " + to_string(k) ) );
		} else if( comp_isnan(dx) ) {
		    throw( Error( ERROR_LOCATION, "Potential NaN at location = " + to_string(i) + 
				  ", " + to_string(j) + ", " + to_string(k) ) );
		}
	    }
	}
    }

    return( sqrt(res) );
}


/* *****************************************************************************
 * CYL
 */

double EpotGSSolver::gs_process_near_solid_cyl( const uint8_t *nearsolid_ptr, 
						uint32_t i, uint32_t j, uint8_t bindex ) const
{
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    double cof = 0.0;
    double epf = 0.0;

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis
    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0/(beta*beta);
	epf += 2.0/(beta*beta)*(*_epot)(i+1,j);
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0/(alpha*alpha);
	epf += 2.0/(alpha*alpha)*(*_epot)(i-1,j);
    } else {
	cof += 2.0/(alpha*beta);
	epf += 2.0/(alpha+beta)*( (*_epot)(i-1,j)/alpha + (*_epot)(i+1,j)/beta );
    }

    // Ymin direction
    alpha = 1.0;
    if( sflag & 0x04 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Ymax direction
    beta = 1.0;
    if( sflag & 0x08 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Y axis
    if( bindex & EPOT_SOLVER_BYMIN ) {
	// On-axis
	cof += 4.0;
	epf += 4.0*(*_epot)(i,j+1);
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	cof += 2.0/(alpha*alpha);
	epf += 2.0/(alpha*alpha)*(*_epot)(i-1,j);	
    } else {
	cof += 2.0/(alpha*beta);
	epf += 1.0/(alpha+beta)*( (2.0/alpha-1.0/j)*(*_epot)(i,j-1) + (2.0/beta+1.0/j)*(*_epot)(i,j+1) );
    }

    if( _plasma == PLASMA_PEXP ) {
	return( solve_pexp_potential( epf, cof, (*_rhs)(i,j), (*_epot)(i,j) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	return( solve_nsimp_potential( epf, cof, (*_rhs)(i,j), (*_epot)(i,j) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	return( solve_shield_potential( epf, cof, (*_rhs)(i,j), (*_epot)(i,j) ) );
    }

    return( (1.0/cof) * ( epf - (*_rhs)(i,j) ) );
}


double EpotGSSolver::gs_process_pure_vacuum_cyl( uint32_t i, uint32_t j ) const
{
    if( _plasma == PLASMA_PEXP ) {
	double epf = (*_epot)(i+1,j) + (*_epot)(i-1,j) 
	    + (1.0+0.5/j)*(*_epot)(i,j+1) 
	    + (1.0-0.5/j)*(*_epot)(i,j-1);
	return( solve_pexp_potential( epf, 4.0, (*_rhs)(i,j), (*_epot)(i,j) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double epf = (*_epot)(i+1,j) + (*_epot)(i-1,j) 
	    + (1.0+0.5/j)*(*_epot)(i,j+1) 
	    + (1.0-0.5/j)*(*_epot)(i,j-1);
	return( solve_nsimp_potential( epf, 4.0, (*_rhs)(i,j), (*_epot)(i,j) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	double epf = (*_epot)(i+1,j) + (*_epot)(i-1,j) 
	    + (1.0+0.5/j)*(*_epot)(i,j+1) 
	    + (1.0-0.5/j)*(*_epot)(i,j-1);
	return( solve_shield_potential( epf, 4.0, (*_rhs)(i,j), (*_epot)(i,j) ) );
    }
 
    return( 0.25 * ( (*_epot)(i+1,j) + (*_epot)(i-1,j)
		     + (1.0+0.5/j)*(*_epot)(i,j+1)
		     + (1.0-0.5/j)*(*_epot)(i,j-1) 
		     - (*_rhs)(i,j) ) );
}


double EpotGSSolver::gs_process_neumann_cyl( uint32_t i, uint32_t j, uint8_t bindex ) const
{
    double cof = 0.0;
    double epf = 0.0;

    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0;
	epf += 2.0*(*_epot)(i+1,j);
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0;
	epf += 2.0*(*_epot)(i-1,j);
    } else {
	cof += 2.0;
	epf += (*_epot)(i+1,j) + (*_epot)(i-1,j);
    }

    if( bindex & EPOT_SOLVER_BYMIN ) {
	// On-axis, special case
	cof += 4.0;
	epf += 4.0*(*_epot)(i,j+1);
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	cof += 2.0;
	epf += 2.0*(*_epot)(i,j-1);
    } else {
	cof += 2.0;
	epf += (1.0+0.5/j)*(*_epot)(i,j+1) + (1.0-0.5/j)*(*_epot)(i,j-1);
    }

    if( _plasma == PLASMA_PEXP ) {
	return( solve_pexp_potential( epf, cof, (*_rhs)(i,j), (*_epot)(i,j) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	return( solve_nsimp_potential( epf, cof, (*_rhs)(i,j), (*_epot)(i,j) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	return( solve_shield_potential( epf, cof, (*_rhs)(i,j), (*_epot)(i,j) ) );
    }

    return( (1.0/cof) * ( epf - (*_rhs)(i,j) ) );
}


double EpotGSSolver::gs_loop_cyl( void ) const
{
    // Go through all nodes
    const double w2 = 1.0-_w;
    const uint32_t dj = _geom.size(0);

    double res = 0.0;
    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	uint32_t a = j*dj;
	for( uint32_t i = 0; i < _geom.size(0); i++, a++ ) {
	    
	    double Vold = (*_epot)(a);
	    double Vnew;
	    uint32_t mesh = _geom.mesh(a);
	    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	    uint8_t bindex = boundary_index(i,j);

	    if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
		Vnew = gs_process_near_solid_cyl( nearsolid_ptr, i, j, bindex );
	    } else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		Vnew = gs_process_pure_vacuum_cyl( i, j );
	    } else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		Vnew = gs_process_neumann_cyl( i, j, bindex );
	    } else {
		// Dirichlet
		continue;
	    }
	    Vnew = _w*Vnew + w2*Vold;
	    (*_epot)(a) = Vnew;
	    double dx = Vnew - Vold;
	    res += dx*dx;
	    if( comp_isinf(dx) ) {
		throw( Error( ERROR_LOCATION, "Potential inf at location = " + to_string(i) + 
			      ", " + to_string(j) ) );
	    } else if( comp_isnan(dx) ) {
		throw( Error( ERROR_LOCATION, "Potential NaN at location = " + to_string(i) + 
			      ", " + to_string(j) ) );
	    }
	}
    }

    return( sqrt(res) );
}


/* *****************************************************************************
 * 2D
 */

double EpotGSSolver::gs_process_near_solid_2d( const uint8_t *nearsolid_ptr, 
					       uint32_t a, uint32_t dj,
					       uint8_t bindex ) const
{
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    double cof = 0.0;
    double epf = 0.0;

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis
    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0/(beta*beta);
	epf += 2.0/(beta*beta)*(*_epot)(a+1);
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0/(alpha*alpha);
	epf += 2.0/(alpha*alpha)*(*_epot)(a-1);
    } else {
	cof += 2.0/(alpha*beta);
	epf += 2.0/(alpha+beta)*( (*_epot)(a-1)/alpha + (*_epot)(a+1)/beta );
    }

    // Ymin direction
    alpha = 1.0;
    if( sflag & 0x04 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Ymax direction
    beta = 1.0;
    if( sflag & 0x08 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Y axis
    if( bindex & EPOT_SOLVER_BYMIN ) {
	cof += 2.0/(beta*beta);
	epf += 2.0/(beta*beta)*(*_epot)(a+dj);
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	cof += 2.0/(alpha*alpha);
	epf += 2.0/(alpha*alpha)*(*_epot)(a-dj);
    } else {
	cof += 2.0/(alpha*beta);
	epf += 2.0/(alpha+beta)*( (*_epot)(a-dj)/alpha + (*_epot)(a+dj)/beta );
    }

    if( _plasma == PLASMA_PEXP ) {
	return( solve_pexp_potential( epf, cof, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	return( solve_nsimp_potential( epf, cof, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	return( solve_shield_potential( epf, cof, (*_rhs)(a), (*_epot)(a) ) );
    }

    return( (1.0/cof) * ( epf - (*_rhs)(a) ) );
}


double EpotGSSolver::gs_process_pure_vacuum_2d( uint32_t a, uint32_t dj ) const
{
    if( _plasma == PLASMA_PEXP ) {
	double epf = (*_epot)(a+1) + (*_epot)(a-1) 
	    + (*_epot)(a+dj) + (*_epot)(a-dj);
	return( solve_pexp_potential( epf, 4.0, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double epf = (*_epot)(a+1) + (*_epot)(a-1) 
	    + (*_epot)(a+dj) + (*_epot)(a-dj);
	return( solve_nsimp_potential( epf, 4.0, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	double epf = (*_epot)(a+1) + (*_epot)(a-1) 
	    + (*_epot)(a+dj) + (*_epot)(a-dj);
	return( solve_shield_potential( epf, 4.0, (*_rhs)(a), (*_epot)(a) ) );
    }

    return( 0.25 * ( (*_epot)(a+1) + (*_epot)(a-1) +
		     (*_epot)(a+dj) + (*_epot)(a-dj) - (*_rhs)(a) ) );
}


double EpotGSSolver::gs_process_neumann_2d( uint32_t a, uint32_t dj, uint8_t bindex ) const
{
    double cof = 0.0;
    double epf = 0.0;

    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0;
	epf += 2.0*(*_epot)(a+1);
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0;
	epf += 2.0*(*_epot)(a-1);
    } else {
	cof += 2.0;
	epf += (*_epot)(a+1) + (*_epot)(a-1);
    }

    if( bindex & EPOT_SOLVER_BYMIN ) {
	cof += 2.0;
	epf += 2.0*(*_epot)(a+dj);
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	cof += 2.0;
	epf += 2.0*(*_epot)(a-dj);
    } else {
	cof += 2.0;
	epf += (*_epot)(a+dj) + (*_epot)(a-dj);
    }

    if( _plasma == PLASMA_PEXP ) {
	return( solve_pexp_potential( epf, cof, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	return( solve_nsimp_potential( epf, cof, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	return( solve_shield_potential( epf, cof, (*_rhs)(a), (*_epot)(a) ) );
    }

    return( (1.0/cof) * ( epf - (*_rhs)(a) ) );
}


double EpotGSSolver::gs_loop_2d( void ) const
{
    // Go through all nodes
    const double w2 = 1.0-_w;
    const uint32_t dj = _geom.size(0);

    double res = 0.0;
    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	uint32_t a = j*dj;
	for( uint32_t i = 0; i < _geom.size(0); i++, a++ ) {
	    
	    double Vold = (*_epot)(a);
	    double Vnew;
	    uint32_t mesh = _geom.mesh(a);
	    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	    uint8_t bindex = boundary_index(i,j);

	    if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
		Vnew = gs_process_near_solid_2d( nearsolid_ptr, a, dj, bindex );
	    } else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		Vnew = gs_process_pure_vacuum_2d( a, dj );
	    } else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		Vnew = gs_process_neumann_2d( a, dj, bindex );
	    } else {
		// Dirichlet
		continue;
	    }
	    Vnew = _w*Vnew + w2*Vold;
	    (*_epot)(a) = Vnew;
	    double dx = Vnew - Vold;
	    res += dx*dx;
	    if( comp_isinf(dx) ) {
		throw( Error( ERROR_LOCATION, "Potential inf at location = " + to_string(i) + 
			      ", " + to_string(j) ) );
	    } else if( comp_isnan(dx) ) {
		throw( Error( ERROR_LOCATION, "Potential NaN at location = " + to_string(i) + 
			      ", " + to_string(j) ) );
	    }
	}
    }

    return( sqrt(res) );
}


/* *****************************************************************************
 * 1D
 */

double EpotGSSolver::gs_process_near_solid_1d( const uint8_t *nearsolid_ptr, 
					       uint32_t i, uint8_t bindex ) const
{
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis
    double cof, epf;
    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof = 2.0/(beta*beta);
	epf = 2.0/(beta*beta)*(*_epot)(i+1);
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof = 2.0/(alpha*alpha);
	epf = 2.0/(alpha*alpha)*(*_epot)(i-1);
    } else {
	cof = 2.0/(alpha*beta);
	epf = 2.0/(alpha+beta)*( (*_epot)(i-1)/alpha + (*_epot)(i+1)/beta );
    }

    if( _plasma == PLASMA_PEXP ) {
	return( solve_pexp_potential( epf, cof, (*_rhs)(i), (*_epot)(i) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	return( solve_nsimp_potential( epf, cof, (*_rhs)(i), (*_epot)(i) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	return( solve_shield_potential( epf, cof, (*_rhs)(i), (*_epot)(i) ) );
    }

    return( (1.0/cof) * ( epf - (*_rhs)(i) ) );
}


double EpotGSSolver::gs_process_pure_vacuum_1d( uint32_t i ) const
{
    if( _plasma == PLASMA_PEXP ) {
	double epf = (*_epot)(i+1) + (*_epot)(i-1);
	return( solve_pexp_potential( epf, 2.0, (*_rhs)(i), (*_epot)(i) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double epf = (*_epot)(i+1) + (*_epot)(i-1);
	return( solve_nsimp_potential( epf, 2.0, (*_rhs)(i), (*_epot)(i) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	double epf = (*_epot)(i+1) + (*_epot)(i-1);
	return( solve_shield_potential( epf, 2.0, (*_rhs)(i), (*_epot)(i) ) );
    }

    return( 0.5 * ( (*_epot)(i+1) + (*_epot)(i-1) - (*_rhs)(i) ) );
}


double EpotGSSolver::gs_process_neumann_1d( uint32_t i, uint8_t bindex ) const
{
    double epf;
    if( bindex & EPOT_SOLVER_BXMIN )
	epf = 2.0*(*_epot)(i+1);
    else
	epf = 2.0*(*_epot)(i-1);

    if( _plasma == PLASMA_PEXP ) {
	return( solve_pexp_potential( epf, 2.0, (*_rhs)(i), (*_epot)(i) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	return( solve_nsimp_potential( epf, 2.0, (*_rhs)(i), (*_epot)(i) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	return( solve_shield_potential( epf, 2.0, (*_rhs)(i), (*_epot)(i) ) );
    }

    // (2*phi_i+1 - 2*phi_i) = -h^2*rho/eps + 2*h*f_N
    return( 0.5*(epf - (*_rhs)(i)) );
}


double EpotGSSolver::gs_loop_1d( void ) const
{
    // Go through all nodes
    double w2 = 1.0-_w;

    double res = 0.0;
    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
	    
	double Vold = (*_epot)(i);
	double Vnew;
	uint32_t mesh = _geom.mesh(i);
	uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	uint8_t bindex = boundary_index(i);

	if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
	    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
	    Vnew = gs_process_near_solid_1d( nearsolid_ptr, i, bindex );
	} else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
	    Vnew = gs_process_pure_vacuum_1d( i );
	} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
	    Vnew = gs_process_neumann_1d( i, bindex );
	} else {
	    // Dirichlet
	    continue;
	}
	Vnew = _w*Vnew + w2*Vold;
	(*_epot)(i) = Vnew;
	double dx = Vnew - Vold;
	res += dx*dx;
	if( comp_isinf(dx) ) {
	    throw( Error( ERROR_LOCATION, "Potential inf at location = " + to_string(i) ) );
	} else if( comp_isnan(dx) ) {
	    throw( Error( ERROR_LOCATION, "Potential NaN at location = " + to_string(i) ) );
	}
    }

    return( sqrt(res) );
}


/* *****************************************************************************
 * Common
 */


double EpotGSSolver::near_solid_neumann_rhs_contribution( uint32_t i, uint32_t j, uint32_t k, 
							  uint8_t bindex, const Vec3D &x ) const
{
    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( _geom.mesh(i,j,k) & SMESH_NEAR_SOLID_INDEX_MASK );
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    double rhs = 0.0;

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis
    if( bindex & EPOT_SOLVER_BXMIN )
	rhs += 2.0*_geom.h()*_geom.get_boundary(1).value(x) / beta;
    else if( bindex & EPOT_SOLVER_BXMAX )
	rhs += 2.0*_geom.h()*_geom.get_boundary(2).value(x) / alpha;

    // Ymin direction
    alpha = 1.0;
    if( sflag & 0x04 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Ymax direction
    beta = 1.0;
    if( sflag & 0x08 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Y axis
    if( _geom.geom_mode() == MODE_CYL ) {
	// Cylindrical geometry
	if( bindex & EPOT_SOLVER_BYMAX )
	    rhs += (1.0)/(2.0*alpha)*(2.0/alpha+1.0/j)*
		2.0*_geom.h()*_geom.get_boundary(4).value(x);
    } else {
	if( bindex & EPOT_SOLVER_BYMIN )
	    rhs += 2.0*_geom.h()*_geom.get_boundary(3).value(x) / beta;
	else if( bindex & EPOT_SOLVER_BYMAX )
	    rhs += 2.0*_geom.h()*_geom.get_boundary(4).value(x) / alpha;
    }

    // Zmin direction
    alpha = 1.0;
    if( sflag & 0x10 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Zmax direction
    beta = 1.0;
    if( sflag & 0x20 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Z axis
    if( bindex & EPOT_SOLVER_BZMIN )
	rhs += 2.0*_geom.h()*_geom.get_boundary(5).value(x) / beta;
    else if( bindex & EPOT_SOLVER_BZMAX )
	rhs += 2.0*_geom.h()*_geom.get_boundary(6).value(x) / alpha;

    return( rhs );
}


void EpotGSSolver::preprocess( const MeshScalarField &scharge )
{
    EpotSolver::preprocess( *_epot );

    // Build right-hand-side
    if( _rhs )
	delete _rhs;
    _rhs = new MeshScalarField( (const Mesh)scharge );

    // Build rhs
    Vec3D x;
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	x[2] = _geom.origo(2) + _geom.h()*k;
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    x[1] = _geom.origo(1) + _geom.h()*j;
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		x[0] = _geom.origo(0) + _geom.h()*i;

		uint32_t mesh = _geom.mesh(i,j,k);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		uint8_t bindex = boundary_index_general(i,j,k);

		if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {

		    (*_rhs)(i,j,k) = -scharge(i,j,k)*_geom.h()*_geom.h()/EPSILON0;

		    if( bindex )
			(*_rhs)(i,j,k) += near_solid_neumann_rhs_contribution( i, j, k, bindex, x );

		} else if( node_id == SMESH_NODE_ID_NEUMANN ) {

		    (*_rhs)(i,j,k) = -scharge(i,j,k)*_geom.h()*_geom.h()/EPSILON0;

		    if( bindex & EPOT_SOLVER_BXMIN )
			(*_rhs)(i,j,k) += 2.0*_geom.h()*_geom.get_boundary(1).value(x);
		    else if( bindex & EPOT_SOLVER_BXMAX )
			(*_rhs)(i,j,k) += 2.0*_geom.h()*_geom.get_boundary(2).value(x);
		    if( _geom.geom_mode() == MODE_CYL ) {
			if( bindex & EPOT_SOLVER_BYMAX )
			    (*_rhs)(i,j,k) += (1.0+0.5/j)*2.0*_geom.h()*_geom.get_boundary(4).value(x);
		    } else {
			if( bindex & EPOT_SOLVER_BYMIN )
			    (*_rhs)(i,j,k) += 2.0*_geom.h()*_geom.get_boundary(3).value(x);
			else if( bindex & EPOT_SOLVER_BYMAX )
			    (*_rhs)(i,j,k) += 2.0*_geom.h()*_geom.get_boundary(4).value(x);
		    }
		    if( bindex & EPOT_SOLVER_BZMIN )
			(*_rhs)(i,j,k) += 2.0*_geom.h()*_geom.get_boundary(5).value(x);
		    else if( bindex & EPOT_SOLVER_BZMAX )
			(*_rhs)(i,j,k) += 2.0*_geom.h()*_geom.get_boundary(6).value(x);

		} else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {

		    (*_rhs)(i,j,k) = -scharge(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
		}
	    }
	}
    }
}


void EpotGSSolver::postprocess( void )
{
    delete _rhs;
    _rhs = NULL;

    EpotSolver::postprocess();
}


double EpotGSSolver::error_scale( void )
{
    double maxsize = _geom.size().max();
    if( _geom.geom_mode() == MODE_3D ) {
	// Coefficients from a fit to spherical condenser 
        // test data with 200x200x200 resolution
	const double a =  0.0642162;
	const double b = -0.0821098;
	const double c =  0.0377858;
	const double d = -0.00640262;
	return( (a + _w*(b + _w*(c + _w*d) ) ) * sqrt(maxsize) );
    } else {
	// Coefficients from a fit to cylindrical condenser 
        // test data with 400x400 resolution
	const double a =  0.0473260;
	const double b = -0.0611217;
	const double c =  0.0284808;
	const double d = -0.00488292;
	return( (a + _w*(b + _w*(c + _w*d) ) ) * maxsize );
    }
}


void EpotGSSolver::prepare_local_gnewton_settings( void )
{
    if( _plasma == PLASMA_PEXP ) {
	// Limit calculation to X times electron temp from plasma potential
	_local_Ulim = _Up - _local_Ulim_fac*_Te;
	ibsimu.message( 1 ) << "Limiting plasma calculation to U > " << _local_Ulim << " V\n";
    } else if( _plasma == PLASMA_NSIMP ) {
	_local_Ulim = 0.0;
	for( uint32_t a = 0; a < _Ei.size(); a++ ) {
	    if( _Ei[a] > _local_Ulim )
		_local_Ulim = _Ei[a];
	}
	// Limit calculation to 10 times maximum energy
	_local_Ulim = _local_Ulim_fac*_local_Ulim;
	ibsimu.message( 1 ) << "Limiting plasma calculation to U < " << _local_Ulim << " V\n";
    } else if( _plasma == PLASMA_SHIELD ) {
	_local_Ulim = fabs(_Up)+fabs(_local_Ulim_fac*_Te);
	ibsimu.message( 1 ) << "Limiting plasma calculation to |U| < " << fabs(_local_Ulim) << " V\n";
    }
}


void EpotGSSolver::subsolve( MeshScalarField &epot, const MeshScalarField &scharge )
{
    //StatusPrint sp( ibsimu.message( 1 ) );
    StatusPrint sp;
    ibsimu.message( 1 ) << "Using Gauss-Seidel solver (" 
			<< "w = " << _w
			<< ", imax = " << _imax
			<< ", eps = " << _eps
			<< ")\n";

    if( !linear() )
	prepare_local_gnewton_settings();

    std::stringstream ss;
    ss << std::setw(5) << 0 << " " << std::scientific << std::setw(20) << 0;
    sp.print( ss.str() );
    
    double errscale = error_scale();

    // Set epot pointer and preprocess
    _epot = &epot;
    preprocess( scharge );

    // Loop until converged
    _iter = 0;
    while( _iter < _imax ) {
	if( _geom.geom_mode() == MODE_3D )
	    _step = gs_loop_3d();
	else if( _geom.geom_mode() == MODE_2D )
	    _step = gs_loop_2d();
	else if( _geom.geom_mode() == MODE_CYL )
	    _step = gs_loop_cyl();
	else if( _geom.geom_mode() == MODE_1D )
	    _step = gs_loop_1d();
	else
	    throw( ErrorUnimplemented( ERROR_LOCATION ) );

	_err = errscale*_step;

	std::stringstream ss;
	ss << "  " << std::setw(5) << _iter << " " << std::scientific << std::setw(20) << _err;
	sp.print( ss.str() );

	_iter++;
	if( _err < _eps )
	    break;
	if( comp_isinf(_err) || comp_isnan(_err) )
	    break;
    }

    // Postprocess
    postprocess();

    std::stringstream ss2;
    ss2 << "  " << std::setw(5) << _iter << " " << std::scientific << std::setw(20) << _err;
    sp.print( ss2.str(), true );
    ibsimu.message( 1 ) << "\n";
    if( _iter == _imax )
	ibsimu.message( 1 ) << "Maximum number of iteration rounds done.\n";
    ibsimu.message( 1 ) << "error estimate = " << _err << "\n";
    ibsimu.message( 1 ) << "iterations = " << _iter << "\n";
}


void EpotGSSolver::save( std::ostream &s ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


void EpotGSSolver::debug_print( std::ostream &os ) const 
{
    EpotSolver::debug_print( os );
    os << "**EpotGSSolver\n";
    os << "imax = " << _imax << "\n";
    os << "eps = " << _eps << "\n";
    os << "w = " << _w << "\n";

}

