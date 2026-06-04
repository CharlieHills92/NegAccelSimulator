/*! \file epot_solver.cpp
 *  \brief Poisson equation problem for solving electric potential.
 */

/* Copyright (c) 2005-2013,2017,2021 Taneli Kalvas. All rights reserved.
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

#include "epot_solver.hpp"
#include "timer.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "compmath.hpp"
#include "constants.hpp"



/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */


EpotSolver::EpotSolver( Geometry &geom ) 
    : _geom(geom), _plasma(PLASMA_NONE), 
      _rhoe(0.0), _Te(0.0), _Up(0.0), 
      _force_pot(0.0), _force_pot_func(0), _force_pot_func2(0), 
      _init_plasma_func(0), _plasma_calc_func(0)
{

}


EpotSolver::EpotSolver( const EpotSolver &epsolver, Geometry &geom )
    : _geom(geom), _plasma(epsolver._plasma),
      _rhoe(epsolver._rhoe), _Te(epsolver._Te), _Up(epsolver._Up),
      _rhoi(epsolver._rhoi), _Ei(epsolver._Ei), 
      _force_pot(epsolver._force_pot), _force_pot_func(epsolver._force_pot_func),
      _force_pot_func2(epsolver._force_pot_func2),
      _init_plasma_func(epsolver._init_plasma_func),
      _plasma_calc_func(epsolver._plasma_calc_func)
{

}


EpotSolver::EpotSolver( Geometry &geom, std::istream &s )
    : _geom(geom)
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


/* ************************************** *
 * EpotSolver constructing               *
 * ************************************** */


void EpotSolver::set_parameters( const EpotSolver &epsolver )
{
    _plasma = epsolver._plasma;
    _rhoe = epsolver._rhoe;
    _Te = epsolver._Te;
    _Up = epsolver._Up;
    _rhoi = epsolver._rhoi;
    _Ei = epsolver._Ei;
    _force_pot = epsolver._force_pot;
    _force_pot_func = epsolver._force_pot_func;
    _force_pot_func2 = epsolver._force_pot_func2;
    _init_plasma_func = epsolver._init_plasma_func;
    _plasma_calc_func = epsolver._plasma_calc_func;
}


void EpotSolver::set_forced_potential_volume( double force_pot, 
					      CallbackFunctorB_V *force_pot_func )
{
    _force_pot = force_pot;
    _force_pot_func = force_pot_func;
}


void EpotSolver::set_forced_potential_volume( CallbackFunctorD_V *force_pot_func )
{
    _force_pot_func2 = force_pot_func;
}


void EpotSolver::set_initial_plasma( double Up, 
				     CallbackFunctorB_V *init_plasma_func )
{
    _plasma     = PLASMA_INITIAL;
    _Up         = Up;
    if( !init_plasma_func )
	throw( Error( ERROR_LOCATION, "NULL initial plasma function" ) );
    _init_plasma_func = init_plasma_func;
    reset_problem();
}


void EpotSolver::set_plasma_calc_region( CallbackFunctorB_V *plasma_calc_func )
{
    _plasma_calc_func = plasma_calc_func;
    reset_problem();
}


void EpotSolver::set_pexp_plasma( double rhoe, double Te, double Up )
{
    _init_plasma_func = NULL;
    _plasma     = PLASMA_PEXP;
    _rhoe       = -fabs(rhoe); // Ensure correct sign of charge density
    _Te         = Te;
    _Up         = Up;
    reset_problem();
}



void EpotSolver::set_nsimp_initial_plasma( CallbackFunctorB_V *init_plasma_func )
{
    _plasma     = PLASMA_INITIAL;
    _Up         = 0.0;
    if( !init_plasma_func )
	throw( Error( ERROR_LOCATION, "NULL initial plasma function" ) );
    _init_plasma_func = init_plasma_func;
    reset_problem();
}


void EpotSolver::set_nsimp_plasma( double rhop, double Ep, 
				   std::vector<double> rhoi, std::vector<double> Ei )
{
    _init_plasma_func = NULL;
    _plasma = PLASMA_NSIMP;
    _rhoi.clear();
    _Ei.clear();

    if( rhop < 0 )
	ibsimu.message( MSG_WARNING, 1 ) << "Warning: negative fast species charge density\n";
    _rhoi.push_back( rhop );
    _Ei.push_back( Ep );

    if(  rhoi.size() != Ei.size() )
	throw( Error( ERROR_LOCATION, "different size rhoi and Ep vectors" ) );
    for( size_t a = 0; a < rhoi.size(); a++ ) {
	if( rhoi[a] < 0 )
	    ibsimu.message( MSG_WARNING, 1 ) << "Warning: negative thermal species " << a << " charge density\n";
	_rhoi.push_back( rhoi[a] );
	_Ei.push_back( Ei[a] );
    }

    reset_problem();
}


void EpotSolver::set_shield_plasma( double Tm, double Um )
{
    _init_plasma_func = NULL;
    _plasma     = PLASMA_SHIELD;
    _Te        = Tm;
    _Up        = Um;
    reset_problem();
}


void EpotSolver::shield_newton( double &rhs, double &drhs, double epot ) const
{
    double K = tanh(_plA*epot - _plB);
    rhs  = 0.5*(1.0+K);
    drhs = 0.5*_plA*(1.0-K*K);
}


void EpotSolver::pexp_newton( double &rhs, double &drhs, double epot ) const
{
    rhs  = _plA*exp( _plB*epot - _plC );
    drhs = _plB*rhs;
}


void EpotSolver::nsimp_newton( double &rhs, double &drhs, double epot ) const
{
    // fast
    double r = _plB*epot;
    if( r > 0.0 ) {
	rhs  = _plA*( 1.0 + 2.0/sqrt(M_PI)*r );
	drhs = _plC;
    } else {
	rhs  = _plA*( 1.0 + erf( r ) );
	drhs = _plC*exp( -r*r );
    }

    // thermal
    for( size_t i = 0; i < _plD.size(); i++ ) {
	double k = _plE[i]*epot;
	if( k > 10.0 ) {
	    // Prevent overruns with exponential function during iteration
	    // by changing to linear slope
	    double exp10 = exp( 10.0 );
	    double d = _plD[i]*exp10;
	    double w = _plD[i]*exp10 + d*(k - 10.0);
	    rhs  += w;
	    drhs += _plE[i]*d;
	} else {
	    double w = _plD[i]*exp( k );
	    rhs  += w;
	    drhs += _plE[i]*w;
	}
    }
}


void EpotSolver::preprocess( MeshScalarField &epot )
{
    if( _plasma == PLASMA_PEXP ) {
	// Calculate plasma parameters for positive ion extraction. 
	// Calculation done as rhs + A*exp(B*x-C)
	_plA = -_rhoe*_geom.h()*_geom.h()/EPSILON0;
	_plB = 1.0/_Te;
	_plC = _Up/_Te;
    } else if( _plasma == PLASMA_NSIMP ) {
	// Calculate plasma parameters for negative ion extraction. 
	// The total right hand side is rhs + r_fast + r_thermal, where
	// r_fast = A(1+erf(B*x)) and
	// r_thermal = D*exp(E*x), for i >= 0. The derivatives are
	// r_fast' = A*B*2/sqrt(pi)*exp(-B^2*x^2)
	// r_thermal' = D*E*exp(E*x)
	_plA = -_rhoi[0]*_geom.h()*_geom.h()/EPSILON0;
	_plB = -1.0/_Ei[0];
	_plC = _plA*_plB*2.0/sqrt(M_PI);
	_plD.clear();
	_plE.clear();
	for( uint32_t i = 1; i < _rhoi.size(); i++ ) {
	    _plD.push_back( -_rhoi[i]*_geom.h()*_geom.h()/EPSILON0 );
	    _plE.push_back( -1.0/_Ei[i] );
	}
    } else if( _plasma == PLASMA_SHIELD ) {
	_plA = 1.0/_Te;
	_plB = _Up/_Te;
    }


    // Set forced vacuum nodes and dirichlet nodes to correct
    // potential. Mark fixed vacuum nodes with a tag.
    Vec3D x;
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	x[2] = _geom.origo(2) + _geom.h()*k;
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    x[1] = _geom.origo(1) + _geom.h()*j;
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		x[0] = _geom.origo(0) + _geom.h()*i;

		uint32_t mesh = _geom.mesh(i,j,k);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		if( node_id == SMESH_NODE_ID_NEAR_SOLID ||
		    node_id == SMESH_NODE_ID_PURE_VACUUM ) {

		    // Vacuum
		    double val;
		    if( _force_pot_func2 && 
			comp_isfinite( (val = (*_force_pot_func2)( x ))) ) {

			// Mark as fixed vacuum
			_geom.mesh(i,j,k) |= SMESH_NODE_FIXED;
			epot(i,j,k) = val;

		    } else if( _force_pot_func && (*_force_pot_func)( x ) ) {

			// Mark as fixed vacuum
			_geom.mesh(i,j,k) |= SMESH_NODE_FIXED;
			epot(i,j,k) = _force_pot;

		    } else if( (_plasma == PLASMA_PEXP_INITIAL ||
				_plasma == PLASMA_NSIMP_INITIAL) && 
			       _init_plasma_func && (*_init_plasma_func)( x ) ) {

			// Mark as fixed vacuum
			_geom.mesh(i,j,k) |= SMESH_NODE_FIXED;
			epot(i,j,k) = _Up;

		    }

		} else if( node_id == SMESH_NODE_ID_NEUMANN ) {

		    if( _force_pot_func && (*_force_pot_func)( x ) ) {

			// Mark as fixed vacuum
			_geom.mesh(i,j,k) = SMESH_NODE_ID_PURE_VACUUM_FIX;
			epot(i,j,k) = _force_pot;

		    } else if( (_plasma == PLASMA_PEXP_INITIAL ||
				_plasma == PLASMA_NSIMP_INITIAL) && 
			       _init_plasma_func && (*_init_plasma_func)( x ) ) {

			// Mark as fixed vacuum
			_geom.mesh(i,j,k) = SMESH_NODE_ID_PURE_VACUUM_FIX;
			epot(i,j,k) = _Up;

		    }

		} else if( node_id == SMESH_NODE_ID_DIRICHLET ) {
		    
		    // Dirichlet
		    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		    epot(i,j,k) = _geom.get_boundary( boundary ).value( x );

		} else {

		    throw( ErrorUnimplemented( ERROR_LOCATION ) );		    

		}
	    }
	}
    }
}


void EpotSolver::postprocess( void )
{
    // Remove fixed vacuum tags
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {

		uint32_t mesh = _geom.mesh(i,j,k);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		if( node_id == SMESH_NODE_ID_NEAR_SOLID_FIX ) {
		    // Change to near solid node, keeping index pointer
		    uint32_t index = SMESH_NEAR_SOLID_INDEX_MASK & mesh;
		    _geom.mesh(i,j,k) = SMESH_NODE_ID_NEAR_SOLID | index;
		} else if( node_id == SMESH_NODE_ID_PURE_VACUUM_FIX ) {
		    // Change to vacuum node if not on boundary
		    if( i == 0 || i == _geom.size(0)-1 )
			continue;
		    if( (_geom.geom_mode() == MODE_2D || _geom.geom_mode() == MODE_CYL ||
			 _geom.geom_mode() == MODE_3D) && (j == 0 || j == _geom.size(1)-1) )
			continue;
		    if( _geom.geom_mode() == MODE_3D && (k == 0 || k == _geom.size(2)-1) )
			continue;
		    _geom.mesh(i,j,k) = SMESH_NODE_ID_PURE_VACUUM;
		}
	    }
	}
    }

    // Change PURE_VACUUM_FIX nodes on Neumann boundaries back to
    // Neumann nodes.

    // Xmin and Xmax
    for( uint32_t bound = 1; bound <= 2; bound++ ) {
	uint32_t i = 0;
	if( bound == 2 ) i = _geom.size(0)-1;
	for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
		uint32_t mesh = _geom.mesh(i,j,k);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		if( node_id == SMESH_NODE_ID_PURE_VACUUM_FIX &&
		    _geom.get_boundary(bound).type() == BOUND_NEUMANN ) {
		    _geom.mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | bound;
		}
	    }
	}
    }
    if( _geom.geom_mode() == MODE_2D || _geom.geom_mode() == MODE_CYL ||
	_geom.geom_mode() == MODE_3D ) {
	// Ymin and Ymax
	for( uint32_t bound = 3; bound <= 4; bound++ ) {
	    uint32_t j = 0;
	    if( bound == 4 ) j = _geom.size(1)-1;
	    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
		for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		    uint32_t mesh = _geom.mesh(i,j,k);
		    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		    if( node_id == SMESH_NODE_ID_PURE_VACUUM_FIX &&
			_geom.get_boundary(bound).type() == BOUND_NEUMANN ) {
			_geom.mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | bound;
		    }
		}
	    }
	}
    }
    if( _geom.geom_mode() == MODE_3D ) {
	// Zmin and Zmax
	for( uint32_t bound = 5; bound <= 6; bound++ ) {
	    uint32_t k = 0;
	    if( bound == 6 ) k = _geom.size(2)-1;
	    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
		for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		    uint32_t mesh = _geom.mesh(i,j,k);
		    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		    if( node_id == SMESH_NODE_ID_PURE_VACUUM_FIX &&
			_geom.get_boundary(bound).type() == BOUND_NEUMANN ) {
			_geom.mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | bound;
		    }
		}
	    }
	}
    }
}


bool EpotSolver::linear( void ) const
{
    switch( _plasma ) {
    case PLASMA_NONE:
    case PLASMA_PEXP_INITIAL:
    case PLASMA_NSIMP_INITIAL:
	return( true );
	break;
    case PLASMA_PEXP:
    case PLASMA_NSIMP:
    case PLASMA_SHIELD:
	return( false );
	break;
    }

    throw( ErrorAssert( ERROR_LOCATION ) );
}


uint8_t EpotSolver::boundary_index( uint32_t i ) const
{
    uint8_t ret = 0;
    if( i == 0 )
	ret += EPOT_SOLVER_BXMIN;
    else if( i == _geom.size(0)-1 )
	ret += EPOT_SOLVER_BXMAX;
    return( ret );
}


uint8_t EpotSolver::boundary_index( uint32_t i, uint32_t j ) const
{
    uint8_t ret = 0;
    if( i == 0 )
	ret += EPOT_SOLVER_BXMIN;
    else if( i == _geom.size(0)-1 )
	ret += EPOT_SOLVER_BXMAX;
    if( j == 0 )
	ret += EPOT_SOLVER_BYMIN;
    else if( j == _geom.size(1)-1 )
	ret += EPOT_SOLVER_BYMAX;
    return( ret );
}


uint8_t EpotSolver::boundary_index( uint32_t i, uint32_t j, uint32_t k ) const
{
    uint8_t ret = 0;
    if( i == 0 )
	ret += EPOT_SOLVER_BXMIN;
    else if( i == _geom.size(0)-1 )
	ret += EPOT_SOLVER_BXMAX;
    if( j == 0 )
	ret += EPOT_SOLVER_BYMIN;
    else if( j == _geom.size(1)-1 )
	ret += EPOT_SOLVER_BYMAX;
    if( k == 0 )
	ret += EPOT_SOLVER_BZMIN;
    else if( k == _geom.size(2)-1 )
	ret += EPOT_SOLVER_BZMAX;
    return( ret );
}


uint8_t EpotSolver::boundary_index_general( uint32_t i, uint32_t j, uint32_t k ) const
{
    uint8_t ret = 0;
    if( i == 0 )
	ret += EPOT_SOLVER_BXMIN;
    else if( i == _geom.size(0)-1 )
	ret += EPOT_SOLVER_BXMAX;
    if( _geom.geom_mode() == MODE_1D )
	return( ret );

    if( j == 0 )
	ret += EPOT_SOLVER_BYMIN;
    else if( j == _geom.size(1)-1 )
	ret += EPOT_SOLVER_BYMAX;
    if( _geom.geom_mode() == MODE_2D || _geom.geom_mode() == MODE_CYL )
	return( ret );

    if( k == 0 )
	ret += EPOT_SOLVER_BZMIN;
    else if( k == _geom.size(2)-1 )
	ret += EPOT_SOLVER_BZMAX;
    return( ret );
}


/* ************************************** *
 * Solver                                 *
 * ************************************** */


MeshScalarField *EpotSolver::evaluate_scharge( const ScalarField &__scharge ) const
{
    // If scharge not defined in same points as geometry, evaluate scharge at nodes
    MeshScalarField *scharge = new MeshScalarField( _geom );
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	double z = k*_geom.h()+_geom.origo(2);
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    double y = j*_geom.h()+_geom.origo(1);
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		double x = i*_geom.h()+_geom.origo(0);
		(*scharge)(i,j,k) = __scharge( Vec3D(x,y,z) );
	    }
	}
    }

    return( scharge );
}


void EpotSolver::solve( MeshScalarField &epot, const ScalarField &__scharge )
{
    Timer t;

    ibsimu.message( 1 ) << "Solving problem\n";
    ibsimu.inc_indent();

    // Check if geometry is built
    if( !_geom.built() )
	throw( Error( ERROR_LOCATION, "geomerty not built" ) );

    // Set scharge to be used by solver
    bool scharge_internal = false;
    const MeshScalarField *scharge = dynamic_cast<const MeshScalarField *>( &__scharge );
    if( scharge == 0 || *scharge != _geom ) {
	if( *scharge != _geom )
	    ibsimu.message( 1 ) << "Converting space charge density to match geometry.\n";
	scharge_internal = true;
	scharge = evaluate_scharge( __scharge );
    }

    // Resize epot if necessary, otherwise leave old solution as a starting point
    if( epot != _geom )
	epot.reset( _geom.geom_mode(), _geom.size(), _geom.origo(), _geom.h() );

    // Solve problem
    subsolve( epot, *scharge );

    // Free scharge if allocated
    if( scharge_internal )
	delete scharge;

    // End timer
    t.stop();

    ibsimu.message( 1 ) << "time used = " << t << "\n";
    ibsimu.dec_indent();
}


/* ************************************** *
 * Misc                                   *
 * ************************************** */


const Geometry &EpotSolver::geometry( void ) const
{
    return( _geom );
}


void EpotSolver::debug_print( std::ostream &os ) const 
{
    os << "**EpotSolver\n";
    os << "plasma = " << _plasma << "\n";
    os << "rhoe = " << _rhoe << "\n";
    os << "Tc = " << _Te << "\n";
    os << "Up = " << _Up << "\n";
}

