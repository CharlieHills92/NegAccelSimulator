/*! \file epot_mgsolver.cpp
 *  \brief Multigrid solver for electric potential problem
 */

/* Copyright (c) 2011-2013,2021 Taneli Kalvas. All rights reserved.
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


#include "epot_mgsubsolver.hpp"
#include "ibsimu.hpp"
#include "constants.hpp"
#include "compmath.hpp"


// Minimum step factor for globally convergent Newton-Raphson
#define FMIN 1.0e-6


/* *****************************************************************************
 * EpotMGSubSolver
 */

EpotMGSubSolver::EpotMGSubSolver( const EpotSolver &epsolver, Geometry &geom,
				  double Ulim, uint32_t imax, double eps )
    : EpotSolver( epsolver, geom ), _Ulim(Ulim), _imax(imax), _eps(eps)
{

}


double EpotMGSubSolver::error_scale_mg( void ) const
{
    double maxsize = _geom.size().max();
    if( _geom.geom_mode() == MODE_3D ) {
	// Coefficients from a fit to spherical condenser
        // test data 11-321 resolutions
	const double a = 0.00347906;
	const double b = 0.234629;
	return( a*pow(maxsize,b) );
    } else {
	// Coefficients from a fit to cylindrical condenser
        // test data with 21-2561 resolutions
	const double a = 0.0034861;
	const double b = 0.767609;
	return( a*pow(maxsize,b) );
    }
}


double EpotMGSubSolver::error_scale( double w ) const
{
    double maxsize = _geom.size().max();
    if( _geom.geom_mode() == MODE_3D ) {
	// Coefficients from a fit to spherical condenser 
        // test data with 200x200x200 resolution
	const double a =  0.0642162;
	const double b = -0.0821098;
	const double c =  0.0377858;
	const double d = -0.00640262;
	return( (a + w*(b + w*(c + w*d) ) ) * sqrt(maxsize) );
    } else {
	// Coefficients from a fit to cylindrical condenser
        // test data with 400x400 resolution
	const double a =  0.0473260;
	const double b = -0.0611217;
	const double c =  0.0284808;
	const double d = -0.00488292;
	return( (a + w*(b + w*(c + w*d) ) ) * maxsize );
    }
}


double EpotMGSubSolver::mg_smooth( MeshScalarField *epot, const MeshScalarField *rhs )
{
    _epot = epot;
    _rhs  = rhs;

    switch( _geom.geom_mode() ) {
    case MODE_1D:
	return( rbgs_loop_1d() );
	break;
    case MODE_2D:
	return( rbgs_loop_2d() );
	break;
    case MODE_CYL:
	return( rbgs_loop_cyl() );
	break;
    case MODE_3D:
	return( rbgs_loop_3d() );
	break;
    default:
	break;
    }

    throw( ErrorAssert( ERROR_LOCATION ) );
}


double EpotMGSubSolver::mg_solve( MeshScalarField *epot, const MeshScalarField *rhs, double w )
{
    _epot = epot;
    _rhs  = rhs;

    switch( _geom.geom_mode() ) {
    case MODE_1D:
	return( sor_loop_1d( w ) );
	break;
    case MODE_2D:
	return( sor_loop_2d( w ) );
	break;
    case MODE_CYL:
	return( sor_loop_cyl( w ) );
	break;
    case MODE_3D:
	return( sor_loop_3d( w ) );
	break;
    default:
	break;
    }

    throw( ErrorAssert( ERROR_LOCATION ) );
}


void EpotMGSubSolver::defect( MeshScalarField *defect, MeshScalarField *epot, const MeshScalarField *rhs,
			      bool after_smooth )
{
    _defect = defect;
    _epot   = epot;
    _rhs    = rhs;

    switch( _geom.geom_mode() ) {
    case MODE_1D:
	defect_1d();
	break;
    case MODE_2D:
	defect_2d();
	break;
    case MODE_CYL:
	defect_cyl();
	break;
    case MODE_3D:
	defect_3d( after_smooth );
	break;
    default:
	break;
    }
}


void EpotMGSubSolver::preprocess( MeshScalarField &epot )
{
    EpotSolver::preprocess( epot );
}


void EpotMGSubSolver::postprocess( void )
{
    EpotSolver::postprocess();
}


double EpotMGSubSolver::solve_shield_potential( double epf, double cof, double rhs, double p ) const
{
    // If potential large enough, no plasma calculation needed
    if( p > _Ulim )
	return( p + (epf - cof*p - rhs) / cof );

    // Globally convergent Newton-Raphson solution for node potential
    double rhst, drhst; 
    shield_newton( rhst, drhst, p );
    double F = epf - cof*p - rhs - rhst;
    double Fnew, pnew;
    for( uint32_t q = 0; q < _imax; q++ ) {
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
	if( fabs(deltap) < _eps )
	    break;
    }
    return( p );
}


double EpotMGSubSolver::solve_nsimp_potential( double epf, double cof, double rhs, double p ) const
{
    // If potential large enough, no plasma calculation needed
    if( p > _Ulim )
	return( p + (epf - cof*p - rhs) / cof );

    // Globally convergent Newton-Raphson solution for node potential
    double rhst, drhst; 
    nsimp_newton( rhst, drhst, p );
    double F = epf - cof*p - rhs - rhst;
    double Fnew, pnew;
    for( uint32_t q = 0; q < _imax; q++ ) {
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
	if( fabs(deltap) < _eps )
	    break;
    }
    return( p );
}


double EpotMGSubSolver::solve_pexp_potential( double epf, double cof, double rhs, double p ) const
{
    // If potential small enough, no plasma calculation needed
    if( p < _Ulim )
	return( p + (epf - cof*p - rhs) / cof );

    // Globally convergent Newton-Raphson solution for node potential
    double rhst, drhst;
    pexp_newton( rhst, drhst, p );
    double F = epf - cof*p - rhs - rhst;
    double Fnew, pnew;
    for( uint32_t q = 0; q < _imax; q++ ) {
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
	if( fabs(deltap) < _eps )
	    break;
    }
    return( p );
}




/* *****************************************************************************
 * 1D Gauss-Seidel
 */


double EpotMGSubSolver::gs_process_near_solid_1d( const uint8_t *nearsolid_ptr, 
						  uint32_t a, uint8_t bindex ) const
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
	epf = 2.0/(beta*beta)*(*_epot)(a+1);
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof = 2.0/(alpha*alpha);
	epf = 2.0/(alpha*alpha)*(*_epot)(a-1);
    } else {
	cof = 2.0/(alpha*beta);
	epf = 2.0/(alpha+beta)*( (*_epot)(a-1)/alpha + (*_epot)(a+1)/beta );
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


double EpotMGSubSolver::gs_process_pure_vacuum_1d( uint32_t a ) const
{
    // (phi_i-1 - 2*phi_i + phi_i+1) / h^2 = rho_i/eps_0
    if( _plasma == PLASMA_PEXP ) {
	double epf = (*_epot)(a+1) + (*_epot)(a-1);
	return( solve_pexp_potential( epf, 2.0, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	double epf = (*_epot)(a+1) + (*_epot)(a-1);
	return( solve_nsimp_potential( epf, 2.0, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	double epf = (*_epot)(a+1) + (*_epot)(a-1);
	return( solve_shield_potential( epf, 2.0, (*_rhs)(a), (*_epot)(a) ) );
    }

    return( (1.0/2.0) * ( (*_epot)(a+1) + (*_epot)(a-1) - (*_rhs)(a) ) );
}


double EpotMGSubSolver::gs_process_neumann_1d( uint32_t a, uint8_t bindex ) const
{
    double epf;
    if( bindex & EPOT_SOLVER_BXMIN )
	epf = 2.0*(*_epot)(a+1);
    else
	epf = 2.0*(*_epot)(a-1);

    if( _plasma == PLASMA_PEXP ) {
	return( solve_pexp_potential( epf, 2.0, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_NSIMP ) {
	return( solve_nsimp_potential( epf, 2.0, (*_rhs)(a), (*_epot)(a) ) );
    } else if( _plasma == PLASMA_SHIELD ) {
	return( solve_shield_potential( epf, 2.0, (*_rhs)(a), (*_epot)(a) ) );
    }

    // (2*phi_i+1 - 2*phi_i) = -h^2*rho/eps + 2*h*f_N
    return( 0.5*(epf - (*_rhs)(a)) );
}


double EpotMGSubSolver::rbgs_loop_1d( void ) const
{
    // Go through internal nodes once using Red-Black ordering
    double res = 0.0;
    for( uint32_t rb = 0; rb < 2; rb++ ) {
	    
	// Odd first, even second
	uint32_t i = 0;
	if( rb != 0 )
	    i = 0;
	else
	    i = 1;

	for( ; i < _geom.size(0); i+=2 ) {
	    
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
	    (*_epot)(i) = Vnew;
	    double dx = Vnew - Vold;
	    res += dx*dx;
	    if( comp_isinf(dx) ) {
		throw( Error( ERROR_LOCATION, "Potential inf at location = " + to_string(i) ) );
	    } else if( comp_isnan(dx) ) {
		throw( Error( ERROR_LOCATION, "Potential NaN at location = " + to_string(i) ) );
	    }
	}
    }

    return( 2.0*sqrt(res) );
}


double EpotMGSubSolver::sor_loop_1d( double w ) const
{
    const double w2 = 1.0-w;
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
	Vnew = w*Vnew + w2*Vold;
	(*_epot)(i) = Vnew;
	double dx = Vnew - Vold;
	res += dx*dx;
	if( comp_isinf(dx) ) {
	    throw( Error( ERROR_LOCATION, "Potential inf at location = " + to_string(i) ) );
	} else if( comp_isnan(dx) ) {
	    throw( Error( ERROR_LOCATION, "Potential NaN at location = " + to_string(i) ) );
	}
    }

    return( 2.0*sqrt(res) );
}


/* *****************************************************************************
 * 1D Defect
 */


double EpotMGSubSolver::defect_near_solid_1d( const uint8_t *nearsolid_ptr, uint32_t a, uint8_t bindex ) const
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
	epf = 2.0/(beta*beta)*(*_epot)(a+1);
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof = 2.0/(alpha*alpha);
	epf = 2.0/(alpha*alpha)*(*_epot)(a-1);
    } else {
	cof = 2.0/(alpha*beta);
	epf = 2.0/(alpha+beta)*( (*_epot)(a-1)/alpha + (*_epot)(a+1)/beta );
    }

    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_SHIELD ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	shield_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a)*rhst );
    }

    return( epf - cof*(*_epot)(a) - (*_rhs)(a) );
}


double EpotMGSubSolver::defect_pure_vacuum_1d( uint32_t a ) const
{
    // (phi_i-1 - 2*phi_i + phi_i+1) / h^2 = rho_i/eps_0
    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( (*_epot)(a+1) + (*_epot)(a-1) - 2.0*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( (*_epot)(a+1) + (*_epot)(a-1) - 2.0*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_SHIELD ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	shield_newton( rhst, drhst, p );
	return( (*_epot)(a+1) + (*_epot)(a-1) - 2.0*p - (*_rhs)(a)*rhst );
    }

    return( (*_epot)(a+1) + (*_epot)(a-1) - 2.0*(*_epot)(a) - (*_rhs)(a) );
}


double EpotMGSubSolver::defect_neumann_1d( uint32_t a, uint8_t bindex ) const
{
    double epf;
    if( bindex & EPOT_SOLVER_BXMIN )
	epf = 2.0*(*_epot)(a+1);
    else
	epf = 2.0*(*_epot)(a-1);

    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( epf - 2.0*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( epf - 2.0*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_SHIELD ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	shield_newton( rhst, drhst, p );
	return( epf - 2.0*p - (*_rhs)(a)*rhst );
    }

    // (2*phi_i+1 - 2*phi_i) = -h^2*rho/eps + 2*h*f_N
    return( epf - 2.0*(*_epot)(a) - (*_rhs)(a) );
}


void EpotMGSubSolver::defect_1d( void ) const
{
    // Go through all nodes
    for( uint32_t i = 0; i < _geom.size(0); i++ ) {

	uint32_t mesh = _geom.mesh(i);
	uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	uint8_t bindex = boundary_index(i);

	double D;
	if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
	    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
	    D = defect_near_solid_1d( nearsolid_ptr, i, bindex );
	} else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
	    D = defect_pure_vacuum_1d( i );
	} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
	    D = defect_neumann_1d( i, bindex );
	} else {
	    D = 0.0;
	}
	(*_defect)(i) = D;
    }
}


/* *****************************************************************************
 * 2D Gauss-Seidel
 */


double EpotMGSubSolver::gs_process_near_solid_2d( const uint8_t *nearsolid_ptr, 
						  uint32_t a, uint32_t dj, uint8_t bindex ) const
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


double EpotMGSubSolver::gs_process_pure_vacuum_2d( uint32_t a, uint32_t dj ) const
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

    return( (1.0/4.0) * ( (*_epot)(a+1) + (*_epot)(a-1) +
			  (*_epot)(a+dj) + (*_epot)(a-dj) - (*_rhs)(a) ) );
}


double EpotMGSubSolver::gs_process_neumann_2d( uint32_t a, uint32_t dj, uint8_t bindex ) const
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


double EpotMGSubSolver::rbgs_loop_2d( void ) const
{
    // Go through all nodes once using Red-Black ordering
    double res = 0.0;
    for( uint32_t rb = 0; rb < 2; rb++ ) {
	const uint32_t dj = _geom.size(0);
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    
	    // Odd first, even second
	    uint32_t i = j % 2;
	    if( i != rb )
		i = 0;
	    else
		i = 1;

	    for( ; i < _geom.size(0); i+=2 ) {
	    
		uint32_t a = j*dj+i;
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
		(*_epot)(a) = Vnew;
		double dx = Vnew - Vold;
		res += dx*dx;
		if( comp_isinf(dx) ) {
		    throw( Error( ERROR_LOCATION, "Potential inf at location = (" + to_string(i) + 
				  ", " + to_string(j) + ")" ) );
		} else if( comp_isnan(dx) ) {
		    throw( Error( ERROR_LOCATION, "Potential NaN at location = (" + to_string(i) + 
				  ", " + to_string(j) + ")" ) );
		}
	    }
	}
    }

    return( 4.0*sqrt(res) );
}


double EpotMGSubSolver::sor_loop_2d( double w ) const
{
    const double w2 = 1.0-w;
    const uint32_t dj = _geom.size(0);

    double res = 0.0;
    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	for( uint32_t i = 0; i < _geom.size(0); i++ ) {

	    uint32_t a = j*dj+i;
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
	    Vnew = w*Vnew + w2*Vold;
	    (*_epot)(a) = Vnew;
	    double dx = Vnew - Vold;
	    res += dx*dx;
	    if( comp_isinf(dx) ) {
		throw( Error( ERROR_LOCATION, "Potential inf at location = (" + to_string(i) + 
			      ", " + to_string(j) + ")" ) );
	    } else if( comp_isnan(dx) ) {
		throw( Error( ERROR_LOCATION, "Potential NaN at location = (" + to_string(i) + 
			      ", " + to_string(j) + ")" ) );
	    }
	}
    }

    return( 4.0*sqrt(res) );
}


/* *****************************************************************************
 * 2D Defect
 */


double EpotMGSubSolver::defect_near_solid_2d( const uint8_t *nearsolid_ptr, 
					      uint32_t a, uint32_t dj, uint8_t bindex ) const
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
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_SHIELD ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	shield_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a)*rhst );
    }

    return( epf - cof*(*_epot)(a) - (*_rhs)(a) );
}



double EpotMGSubSolver::defect_pure_vacuum_2d( uint32_t a, uint32_t dj ) const
{
    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( (*_epot)(a+1) + (*_epot)(a-1) + (*_epot)(a+dj) 
		+ (*_epot)(a-dj) - 4.0*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( (*_epot)(a+1) + (*_epot)(a-1) + (*_epot)(a+dj) 
		+ (*_epot)(a-dj) - 4.0*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_SHIELD ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	shield_newton( rhst, drhst, p );
	return( (*_epot)(a+1) + (*_epot)(a-1) + (*_epot)(a+dj) 
		+ (*_epot)(a-dj) - 4.0*p - (*_rhs)(a)*rhst );
    }

    return( (*_epot)(a+1) + (*_epot)(a-1) +
	    (*_epot)(a+dj) + (*_epot)(a-dj) - 4.0*(*_epot)(a) - (*_rhs)(a) );
}


double EpotMGSubSolver::defect_neumann_2d( uint32_t a, uint32_t dj, uint8_t bindex ) const
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
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_SHIELD ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	shield_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a)*rhst );
    }

    return( epf - cof*(*_epot)(a) - (*_rhs)(a) );
}


void EpotMGSubSolver::defect_2d( void ) const
{
    // Go through all nodes
    const uint32_t dj = _geom.size(0);
    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	for( uint32_t i = 0; i < _geom.size(0); i++ ) {

	    uint32_t a = j*dj+i;
	    uint32_t mesh = _geom.mesh(a);
	    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	    uint8_t bindex = boundary_index(i,j);

	    double D;
	    if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
		D = defect_near_solid_2d( nearsolid_ptr, a, dj, bindex );
	    } else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		D = defect_pure_vacuum_2d( a, dj );
	    } else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		D = defect_neumann_2d( a, dj, bindex );
	    } else {
		D = 0.0;
	    }
	    (*_defect)(a) = D;
	}
    }
}



/* *****************************************************************************
 * CYL Gauss-Seidel
 */


double EpotMGSubSolver::gs_process_near_solid_cyl( const uint8_t *nearsolid_ptr, 
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
	epf += 2.0/(alpha*alpha)*(*_epot)(i,j-1);
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


double EpotMGSubSolver::gs_process_pure_vacuum_cyl( uint32_t i, uint32_t j ) const
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
 
    return( (1.0/4.0) * ( (*_epot)(i+1,j) + (*_epot)(i-1,j)
			  + (1.0+0.5/j)*(*_epot)(i,j+1)
			  + (1.0-0.5/j)*(*_epot)(i,j-1) 
			  - (*_rhs)(i,j) ) );
}


double EpotMGSubSolver::gs_process_neumann_cyl( uint32_t i, uint32_t j, uint8_t bindex ) const
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


double EpotMGSubSolver::rbgs_loop_cyl( void ) const
{
    const uint32_t dj = _geom.size(0);

    // Go through all nodes once using Red-Black ordering
    double res = 0.0;
    for( uint32_t rb = 0; rb < 2; rb++ ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    
	    // Odd first, even second
	    uint32_t i = j % 2;
	    if( i != rb )
		i = 0;
	    else
		i = 1;

	    for( ; i < _geom.size(0); i+=2 ) {
	    
		uint32_t a = j*dj+i;
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
		(*_epot)(a) = Vnew;
		double dx = Vnew - Vold;
		res += dx*dx;
		if( comp_isinf(dx) ) {
		    throw( Error( ERROR_LOCATION, "Potential inf at location = (" + to_string(i) + 
				  ", " + to_string(j) + ")" ) );
		} else if( comp_isnan(dx) ) {
		    throw( Error( ERROR_LOCATION, "Potential NaN at location = (" + to_string(i) + 
				  ", " + to_string(j) + ")" ) );
		}
	    }
	}
    }

    return( 4.0*sqrt(res) );
}


double EpotMGSubSolver::sor_loop_cyl( double w ) const
{
    const double w2 = 1.0-w;
    const uint32_t dj = _geom.size(0);

    double res = 0.0;
    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	for( uint32_t i = 0; i < _geom.size(0); i++ ) {

	    uint32_t a = j*dj+i;
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
	    Vnew = w*Vnew + w2*Vold;
	    (*_epot)(a) = Vnew;
	    double dx = Vnew - Vold;
	    res += dx*dx;
	    if( comp_isinf(dx) ) {
		throw( Error( ERROR_LOCATION, "Potential inf at location = (" + to_string(i) + 
			      ", " + to_string(j) + ")" ) );
	    } else if( comp_isnan(dx) ) {
		throw( Error( ERROR_LOCATION, "Potential NaN at location = (" + to_string(i) + 
			      ", " + to_string(j) + ")" ) );
	    }
	}
    }

    return( 4.0*sqrt(res) );
}


/* *****************************************************************************
 * CYL Defect
 */


double EpotMGSubSolver::defect_near_solid_cyl( const uint8_t *nearsolid_ptr, 
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
	epf += 2.0/(alpha*alpha)*(*_epot)(i,j-1);	
    } else {
	cof += 2.0/(alpha*beta);
	epf += 1.0/(alpha+beta)*( (2.0/alpha-1.0/j)*(*_epot)(i,j-1) + (2.0/beta+1.0/j)*(*_epot)(i,j+1) );
    }

    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(i,j);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(i,j) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(i,j);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(i,j) - rhst );
    } else if( _plasma == PLASMA_SHIELD ) {
	double p = (*_epot)(i,j);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(i,j)*rhst );
    }

    return( epf - cof*(*_epot)(i,j) - (*_rhs)(i,j) );
}



double EpotMGSubSolver::defect_pure_vacuum_cyl( uint32_t i, uint32_t j ) const
{
    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(i,j);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( (*_epot)(i+1,j) + (*_epot)(i-1,j)
		+ (1.0+0.5/j)*(*_epot)(i,j+1) + (1.0-0.5/j)*(*_epot)(i,j-1) 
		- 4.0*p - (*_rhs)(i,j) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(i,j);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( (*_epot)(i+1,j) + (*_epot)(i-1,j) 
		+ (1.0+0.5/j)*(*_epot)(i,j+1) + (1.0-0.5/j)*(*_epot)(i,j-1) 
		- 4.0*p - (*_rhs)(i,j) - rhst );
    } else if( _plasma == PLASMA_SHIELD ) {
	double p = (*_epot)(i,j);
	double rhst, drhst;
	shield_newton( rhst, drhst, p );
	return( (*_epot)(i+1,j) + (*_epot)(i-1,j) 
		+ (1.0+0.5/j)*(*_epot)(i,j+1) + (1.0-0.5/j)*(*_epot)(i,j-1) 
		- 4.0*p - (*_rhs)(i,j)*rhst );
    }
 
    return( (*_epot)(i+1,j) + (*_epot)(i-1,j)
	    + (1.0+0.5/j)*(*_epot)(i,j+1) + (1.0-0.5/j)*(*_epot)(i,j-1) 
	    - 4.0*(*_epot)(i,j) - (*_rhs)(i,j) );
}


double EpotMGSubSolver::defect_neumann_cyl( uint32_t i, uint32_t j, uint8_t bindex ) const
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
	double p = (*_epot)(i,j);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(i,j) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(i,j);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(i,j) - rhst );
    } else if( _plasma == PLASMA_SHIELD ) {
	double p = (*_epot)(i,j);
	double rhst, drhst;
	shield_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(i,j)*rhst );
    }

    return( epf - cof*(*_epot)(i,j) - (*_rhs)(i,j) );
}


void EpotMGSubSolver::defect_cyl( void ) const
{
    // Go through all nodes
    const uint32_t dj = _geom.size(0);
    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	for( uint32_t i = 0; i < _geom.size(0); i++ ) {

	    uint32_t a = j*dj+i;
	    uint32_t mesh = _geom.mesh(a);
	    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
	    uint8_t bindex = boundary_index(i,j);

	    double D;
	    if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
		D = defect_near_solid_cyl( nearsolid_ptr, i, j, bindex );
	    } else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		D = defect_pure_vacuum_cyl( i, j );
	    } else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		D = defect_neumann_cyl( i, j, bindex );
	    } else {
		D = 0.0;
	    }
	    (*_defect)(a) = D;
	}
    }
}


/* *****************************************************************************
 * 3D Gauss-Seidel
 */


double EpotMGSubSolver::gs_process_near_solid_3d( const uint8_t *nearsolid_ptr, 
						  uint32_t a, uint32_t dj, uint32_t dk, uint8_t bindex ) const
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


double EpotMGSubSolver::gs_process_pure_vacuum_3d( uint32_t a, uint32_t dj, uint32_t dk ) const
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

    return( (1.0/6.0) * ( (*_epot)(a+1) + (*_epot)(a-1) + 
			  (*_epot)(a+dj) + (*_epot)(a-dj) +
			  (*_epot)(a+dk) + (*_epot)(a-dk) - (*_rhs)(a) ) );
}


double EpotMGSubSolver::gs_process_neumann_3d( uint32_t a, uint32_t dj, uint32_t dk, uint8_t bindex ) const
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


double EpotMGSubSolver::rbgs_loop_3d( void ) const
{
    const uint32_t dj = _geom.size(0);
    const uint32_t dk = _geom.size(0)*_geom.size(1);

    // Go through all nodes once using Red-Black ordering
    double res = 0.0;
    for( uint32_t rb = 0; rb < 2; rb++ ) {
	for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    
		// Odd first, even second
		uint32_t i = (k+j) % 2;
		if( i != rb )
		    i = 0;
		else
		    i = 1;

		for( ; i < _geom.size(0); i+=2 ) {
	    
		    uint32_t a = k*dk+j*dj+i;
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
			// Dirichlet
			continue;
		    }
		    (*_epot)(a) = Vnew;
		    double dx = Vnew - Vold;
		    res += dx*dx;
		    if( comp_isinf(dx) ) {
			throw( Error( ERROR_LOCATION, "Potential inf at location = (" + to_string(i) + 
				      ", " + to_string(j) + ", " + to_string(k) + ")" ) );
		    } else if( comp_isnan(dx) ) {
			throw( Error( ERROR_LOCATION, "Potential NaN at location = (" + to_string(i) + 
				      ", " + to_string(j) + ", " + to_string(k) + ")" ) );
		    }
		}
	    }
	}
    }

    return( 6.0*sqrt(res) );
}


double EpotMGSubSolver::sor_loop_3d( double w ) const
{
    const double w2 = 1.0-w;
    const uint32_t dj = _geom.size(0);
    const uint32_t dk = _geom.size(0)*_geom.size(1);

    double res = 0.0;
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {

		uint32_t a = k*dk+j*dj+i;
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
		    // Dirichlet
		    continue;
		}
		Vnew = w*Vnew + w2*Vold;
		(*_epot)(a) = Vnew;
		double dx = Vnew - Vold;
		res += dx*dx;
		if( comp_isinf(dx) ) {
		    throw( Error( ERROR_LOCATION, "Potential inf at location = (" + to_string(i) + 
				  ", " + to_string(j) + ", " + to_string(k) + ")" ) );
		} else if( comp_isnan(dx) ) {
		    throw( Error( ERROR_LOCATION, "Potential NaN at location = (" + to_string(i) + 
				  ", " + to_string(j) + ", " + to_string(k) + ")" ) );
		}
	    }
	}
    }

    return( 6.0*sqrt(res) );
}


/* *****************************************************************************
 * 3D Defect
 */


double EpotMGSubSolver::defect_near_solid_3d( const uint8_t *nearsolid_ptr, 
					      uint32_t a, uint32_t dj, uint32_t dk, uint8_t bindex ) const
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
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_SHIELD ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	shield_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a)*rhst );
    }

    return( epf - cof*(*_epot)(a) - (*_rhs)(a) );
}



double EpotMGSubSolver::defect_pure_vacuum_3d( uint32_t a, uint32_t dj, uint32_t dk ) const
{
    if( _plasma == PLASMA_PEXP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( (*_epot)(a+1) + (*_epot)(a-1) 
		+ (*_epot)(a+dj) + (*_epot)(a-dj) 
		+ (*_epot)(a+dk) + (*_epot)(a-dk) 
		- 6.0*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( (*_epot)(a+1) + (*_epot)(a-1) 
		+ (*_epot)(a+dj) + (*_epot)(a-dj) 
		+ (*_epot)(a+dk) + (*_epot)(a-dk) 
		- 6.0*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_SHIELD ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	shield_newton( rhst, drhst, p );
	return( (*_epot)(a+1) + (*_epot)(a-1) 
		+ (*_epot)(a+dj) + (*_epot)(a-dj) 
		+ (*_epot)(a+dk) + (*_epot)(a-dk) 
		- 6.0*p - (*_rhs)(a)*rhst );
    }

    return( (*_epot)(a+1) + (*_epot)(a-1) 
	    + (*_epot)(a+dj) + (*_epot)(a-dj) 
	    + (*_epot)(a+dk) + (*_epot)(a-dk) 
	    - 6.0*(*_epot)(a) - (*_rhs)(a) );
}


double EpotMGSubSolver::defect_neumann_3d( uint32_t a, uint32_t dj, 
					   uint32_t dk, uint8_t bindex ) const
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
	double p = (*_epot)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a) - rhst );
    } else if( _plasma == PLASMA_SHIELD ) {
	double p = (*_epot)(a);
	double rhst, drhst;
	shield_newton( rhst, drhst, p );
	return( epf - cof*p - (*_rhs)(a)*rhst );
    }

    return( epf - cof*(*_epot)(a) - (*_rhs)(a) );
}


void EpotMGSubSolver::defect_3d( bool after_smooth ) const
{
    // Go through all nodes
    const uint32_t dj = _geom.size(0);
    const uint32_t dk = _geom.size(0)*_geom.size(1);
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		
		uint32_t a = k*dk+j*dj+i;
		// RBGS smoother does even pass last. Even nodes have
		// zero defect.
		if( after_smooth && (i+k+j) % 2 == 0 ) {
		    (*_defect)(a) = 0.0;
		    continue;
		}

		uint32_t mesh = _geom.mesh(a);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		uint8_t bindex = boundary_index(i,j,k);

		double D;
		if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( mesh & SMESH_NEAR_SOLID_INDEX_MASK );
		    D = defect_near_solid_3d( nearsolid_ptr, a, dj, dk, bindex );
		} else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		    D = defect_pure_vacuum_3d( a, dj, dk );
		} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		    D = defect_neumann_3d( a, dj, dk, bindex );
		} else {
		    D = 0.0;
		}
		(*_defect)(a) = D;
	    }
	}
    }
}



