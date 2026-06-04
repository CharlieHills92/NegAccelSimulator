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


#include "epot_mgsolver.hpp"
#include "ibsimu.hpp"
#include "constants.hpp"
#include "compmath.hpp"
#include "statusprint.hpp"


//#define DEBUG_MGSOLVER 1


/* *****************************************************************************
 * EpotMGSolver
 */


EpotMGSolver::EpotMGSolver( Geometry &geom )
    : EpotSolver( geom ), _geom_prepared(false), _levels(1), _npre(5), 
      _npost(5), _mgcycmax(100), _mgcyc(0), _mgeps(1.0e-4), _mgerr(0.0),
      _gamma(1), _step(0.0), _coarse_eps(1.0e-10), _coarse_err(0.0), 
      _w(1.7), _imax(10000), _coarse_steps(0), 
      _local_Ulim(0.0), _local_imax(1), _local_eps(1.0e-6)
{

}


EpotMGSolver::EpotMGSolver( Geometry &geom, std::istream &s )
    : EpotSolver(geom,s)
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


EpotMGSolver::~EpotMGSolver()
{
    reset_problem();
}


void EpotMGSolver::reset_problem( void )
{
    //std::cout << "Reset!\n";

    // Free geom and subsolver arrays made in prepare_mg_geom()
    for( uint32_t a = 0; a < _geomv.size(); a++ ) {
	if( a >= 1 )
	    delete _geomv[a];
	delete _epotsolverv[a];
    }

    // Free other arrays made in preprocess()
    for( uint32_t a = 0; a < _epotv.size(); a++ ) {
	delete _rhsv[a];
	delete _workv[a];
	if( _work2v[a] )
	    delete _work2v[a];
	if( a >= 1 )
	    delete _epotv[a];
    }

    _epotv.clear();
    _geomv.clear();
    _epotsolverv.clear();
    _rhsv.clear();
    _workv.clear();
    _work2v.clear();

    _step = 0.0;
    _mgerr = 0.0;
    _coarse_err = 0.0;
    _coarse_steps = 0;
    _geom_prepared = false;
}


void EpotMGSolver::set_eps( double eps )
{
    _coarse_eps = eps;
}


void EpotMGSolver::set_w( double w )
{
    _w = w;
}


void EpotMGSolver::set_imax( uint32_t imax )
{
    _imax = imax;
}


void EpotMGSolver::set_local_imax( uint32_t local_imax )
{
    _local_imax = local_imax;
}


double EpotMGSolver::get_potential_change_norm( void ) const
{
    return( _step );
}


double EpotMGSolver::get_error_estimate( void ) const
{
    return( _mgerr );
}


uint32_t EpotMGSolver::get_mgcyc( void ) const
{
    return( _mgcyc );
}


void EpotMGSolver::set_levels( uint32_t levels )
{
    if( levels < 1 )
	throw( Error( ERROR_LOCATION, "Invalid number of multigrid levels" ) ); 
    _levels = levels;
    reset_problem();
}


void EpotMGSolver::set_mgcycmax( uint32_t mgcycmax )
{
    _mgcycmax = mgcycmax;
}


void EpotMGSolver::set_mgeps( double mgeps )
{
    _mgeps = mgeps;
}


void EpotMGSolver::set_gamma( uint32_t gamma )
{
    if( gamma < 1 )
	throw( Error( ERROR_LOCATION, "Invalid multigrid cycle coefficient" ) ); 
    _gamma = gamma;
}


void EpotMGSolver::set_npre( uint32_t npre )
{
    if( npre < 1 )
	throw( Error( ERROR_LOCATION, "Invalid number of multigrid pre-smoothing cycles" ) ); 
    _npre = npre;
}


void EpotMGSolver::set_npost( uint32_t npost )
{
    if( npost < 1 )
	throw( Error( ERROR_LOCATION, "Invalid number of multigrid post-smoothing cycles" ) ); 
    _npost = npost;
}


/* *****************************************************************************
 * Common
 */


/* Create hierarchy of geometries and solvers.
 */
void EpotMGSolver::prepare_mg_geom( void )
{
    ibsimu.message( 1 ) << "Preparing multigrid geometries\n";

    geom_mode_e geom_mode = _geom.geom_mode();
    Int3D size = _geom.size();
    Vec3D origo = _geom.origo();
    double h = _geom.h();

    ibsimu.inc_indent();
    ibsimu.set_message_threshold( MSG_VERBOSE, ibsimu.get_message_threshold( MSG_VERBOSE ) - 1 );
    for( uint32_t a = 0; a < _levels; a++ ) {

	if( a == 0 ) {

	    // First (finest) level
	    _geomv.push_back( &_geom );
	    EpotMGSubSolver *mgss = new EpotMGSubSolver( *this, _geom, 
							 _local_Ulim, _local_imax, _local_eps );
	    _epotsolverv.push_back( mgss );

	} else {

	    // Build mesh density	    
	    for( uint32_t b = 0; b < _geom.dim(); b++ ) {
		if( size[b] % 2 == 0 )
		    throw( Error( ERROR_LOCATION, "Incorrect parity of mesh size " + 
				  to_string(size[b]) + " in direction " + to_string(b) + 
				  " at level " + to_string(a) ) );		    
		size[b] = (size[b]+1)/2;
	    }
	    h *= 2.0;

	    // Build geometry
	    Geometry *geom = new Geometry( geom_mode, size, origo, h );
	    for( uint32_t b = 7; b <= _geom.number_of_boundaries(); b++ )
		geom->set_solid( b, _geom.get_solid(b) );
	    for( uint32_t b = 1; b <= _geom.number_of_boundaries(); b++ ) {
		// Boundaries are of same type but value is zero.
		Bound bound = _geom.get_boundary(b);
		bound.set_value( 0.0 );
		geom->set_boundary( b, bound );
	    }
	    geom->build_mesh();
	    _geomv.push_back( geom );

	    // Make sub-solver copying parameters from this
	    EpotMGSubSolver *mgss = new EpotMGSubSolver( *this, *_geomv.back(), 
							 _local_Ulim, _local_imax, _local_eps );
	    _epotsolverv.push_back( mgss );
	}
    }
    ibsimu.set_message_threshold( MSG_VERBOSE, ibsimu.get_message_threshold( MSG_VERBOSE ) + 1 );
    ibsimu.dec_indent();
    ibsimu.message( 1 ) << "Done\n";
    ibsimu.flush();

    _geom_prepared = true;
}


double EpotMGSolver::near_solid_neumann_rhs_contribution( uint32_t i, uint32_t j, uint32_t k, 
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


void EpotMGSolver::preprocess( MeshScalarField &epot, const MeshScalarField &scharge )
{
    ibsimu.message( 1 ) << "Preprocessing\n";

    // Run preprocess geometry for all EpotSolvers, allocate rhs and epot fields 
    // and one working field
    for( uint32_t a = 0; a < _levels; a++ ) {

	if( a == 0 ) {
	    _epotv.push_back( &epot );
	} else {
	    MeshScalarField *ep = new MeshScalarField( (const Mesh)(*_geomv[a]) );
	    _epotv.push_back( ep );
	}

	MeshScalarField *rhs = new MeshScalarField( (const Mesh)(*_geomv[a]) );
	_rhsv.push_back( rhs );

	MeshScalarField *work = new MeshScalarField( (const Mesh)(*_geomv[a]) );
	_workv.push_back( work );

	if( !linear() ) {
	    if( a == 0 ) {
		_work2v.push_back( NULL );
	    } else {
		MeshScalarField *work2 = new MeshScalarField( (const Mesh)(*_geomv[a]) );
		_work2v.push_back( work2 );
	    }
	}

	// Preprocess solid meshes
	_epotsolverv[a]->preprocess( *_epotv[a] );
    }

    // Build rhs for top level
    Vec3D x;
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	x[2] = _geom.origo(2) + _geom.h()*k;
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    x[1] = _geom.origo(1) + _geom.h()*j;
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		x[0] = _geom.origo(0) + _geom.h()*i;

		uint32_t mesh = _geomv[0]->mesh(i,j,k);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		uint8_t bindex = boundary_index_general(i,j,k);

		if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {

		    (*_rhsv[0])(i,j,k) = -scharge(i,j,k)*_geom.h()*_geom.h()/EPSILON0;

		    if( bindex )
			(*_rhsv[0])(i,j,k) += near_solid_neumann_rhs_contribution( i, j, k, bindex, x );

		} else if( node_id == SMESH_NODE_ID_NEUMANN ) {

		    (*_rhsv[0])(i,j,k) = -scharge(i,j,k)*_geom.h()*_geom.h()/EPSILON0;

		    if( bindex & EPOT_SOLVER_BXMIN )
			(*_rhsv[0])(i,j,k) += 2.0*_geom.h()*_geom.get_boundary(1).value(x);
		    else if( bindex & EPOT_SOLVER_BXMAX )
			(*_rhsv[0])(i,j,k) += 2.0*_geom.h()*_geom.get_boundary(2).value(x);
		    if( _geom.geom_mode() == MODE_CYL ) {
			if( bindex & EPOT_SOLVER_BYMAX )
			    (*_rhsv[0])(i,j,k) += (1.0+0.5/j)*2.0*_geom.h()*_geom.get_boundary(4).value(x);
		    } else {
			if( bindex & EPOT_SOLVER_BYMIN )
			    (*_rhsv[0])(i,j,k) += 2.0*_geom.h()*_geom.get_boundary(3).value(x);
			else if( bindex & EPOT_SOLVER_BYMAX )
			    (*_rhsv[0])(i,j,k) += 2.0*_geom.h()*_geom.get_boundary(4).value(x);
		    }
		    if( bindex & EPOT_SOLVER_BZMIN )
			(*_rhsv[0])(i,j,k) += 2.0*_geom.h()*_geom.get_boundary(5).value(x);
		    else if( bindex & EPOT_SOLVER_BZMAX )
			(*_rhsv[0])(i,j,k) += 2.0*_geom.h()*_geom.get_boundary(6).value(x);

		} else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {

		    (*_rhsv[0])(i,j,k) = -scharge(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
		}
	    }
	}
    }
}


void EpotMGSolver::postprocess( void )
{
    ibsimu.message( 1 ) << "Postprocessing\n";

    for( uint32_t a = 0; a < _levels; a++ ) {

	if( a != 0 )
	    delete _epotv[a];
	delete _rhsv[a];
	delete _workv[a];
	if( !linear() && a != 0 )
	    delete _work2v[a];

	// Postprocess solid meshes
	_epotsolverv[a]->postprocess();
    }

    _epotv.clear();
    _rhsv.clear();
    _workv.clear();
    _work2v.clear();
}


void EpotMGSolver::restrict_3d( MeshScalarField *out, const MeshScalarField *in, bool defect )
{
    // Go through internal nodes of rougher level (out)
    int32_t s = out->size(0)-1;
    int32_t t = out->size(1)-1;
    int32_t u = out->size(2)-1;
    for( int32_t k = 1; k < u; k++ ) {
	for( int32_t j = 1; j < t; j++ ) {
	    for( int32_t i = 1; i < s; i++ ) {

		int32_t ii = 2*i;
		int32_t jj = 2*j;
		int32_t kk = 2*k;

		// RBGS smoother does even pass last. Even nodes have
		// zero defect.
		if( defect ) {
		    double D = 
			1.0/16.0*(*in)( ii-1, jj-1, kk-1 ) + 
			1.0/16.0*(*in)( ii-1, jj-1, kk+1 ) + 
			1.0/4.0 *(*in)( ii-1, jj,   kk   ) + 
			1.0/16.0*(*in)( ii-1, jj+1, kk-1 ) + 
			1.0/16.0*(*in)( ii-1, jj+1, kk+1 ) + 

			1.0/4.0 *(*in)( ii,   jj-1, kk   ) + 
			1.0/4.0 *(*in)( ii,   jj,   kk-1 ) + 
			1.0/4.0 *(*in)( ii,   jj,   kk+1 ) + 
			1.0/4.0 *(*in)( ii,   jj+1, kk   ) + 

			1.0/16.0*(*in)( ii+1, jj-1, kk-1 ) + 
			1.0/16.0*(*in)( ii+1, jj-1, kk+1 ) + 
			1.0/4.0 *(*in)( ii+1, jj,   kk   ) + 
			1.0/16.0*(*in)( ii+1, jj+1, kk-1 ) + 
			1.0/16.0*(*in)( ii+1, jj+1, kk+1 );
		    (*out)(i,j,k) = D;
		} else {
		    double D = 
			1.0/64.0*(*in)( ii-1, jj-1, kk-1 ) + 
			1.0/32.0*(*in)( ii,   jj-1, kk-1 ) + 
			1.0/64.0*(*in)( ii+1, jj-1, kk-1 ) + 
		    
			1.0/32.0*(*in)( ii-1, jj,   kk-1 ) + 
			1.0/16.0*(*in)( ii,   jj,   kk-1 ) + 
			1.0/32.0*(*in)( ii+1, jj,   kk-1 ) + 
		    
			1.0/64.0*(*in)( ii-1, jj+1, kk-1 ) + 
			1.0/32.0*(*in)( ii,   jj+1, kk-1 ) + 
			1.0/64.0*(*in)( ii+1, jj+1, kk-1 ) + 

			//
		    
			1.0/32.0*(*in)( ii-1, jj-1, kk   ) + 
			1.0/16.0*(*in)( ii,   jj-1, kk   ) + 
			1.0/32.0*(*in)( ii+1, jj-1, kk   ) + 
		    
			1.0/16.0*(*in)( ii-1, jj,   kk   ) + 
			1.0/8.0 *(*in)( ii,   jj,   kk   ) + 
			1.0/16.0*(*in)( ii+1, jj,   kk   ) + 
		    
			1.0/32.0*(*in)( ii-1, jj+1, kk   ) + 
			1.0/16.0*(*in)( ii,   jj+1, kk   ) + 
			1.0/32.0*(*in)( ii+1, jj+1, kk   ) + 

			//

			1.0/64.0*(*in)( ii-1, jj-1, kk+1 ) + 
			1.0/32.0*(*in)( ii,   jj-1, kk+1 ) + 
			1.0/64.0*(*in)( ii+1, jj-1, kk+1 ) + 
		    
			1.0/32.0*(*in)( ii-1, jj,   kk+1 ) + 
			1.0/16.0*(*in)( ii,   jj,   kk+1 ) + 
			1.0/32.0*(*in)( ii+1, jj,   kk+1 ) + 
		    
			1.0/64.0*(*in)( ii-1, jj+1, kk+1 ) + 
			1.0/32.0*(*in)( ii,   jj+1, kk+1 ) + 
			1.0/64.0*(*in)( ii+1, jj+1, kk+1 );
		    (*out)(i,j,k) = D;
		}
	    }
	}
    }

    // FACES

    // i=0 boundary nodes of level+1
    for( int32_t k = 1; k < u; k++ ) {
	for( int32_t j = 1; j < t; j++ ) {

	int32_t i = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/16.0*(*in)( ii,   jj-1, kk-1 ) + 
	    1.0/8.0 *(*in)( ii,   jj,   kk-1 ) + 
	    1.0/16.0*(*in)( ii,   jj+1, kk-1 ) + 

	    1.0/8.0 *(*in)( ii,   jj-1, kk   ) + 
	    1.0/4.0 *(*in)( ii,   jj,   kk   ) + 
	    1.0/8.0 *(*in)( ii,   jj+1, kk   ) + 

	    1.0/16.0*(*in)( ii,   jj-1, kk+1 ) + 
	    1.0/8.0 *(*in)( ii,   jj,   kk+1 ) + 
	    1.0/16.0*(*in)( ii,   jj+1, kk+1 );
	
	(*out)(i,j,k) = D;
	}
    }

    // i=size(0)-1 boundary nodes of level+1
    for( int32_t k = 1; k < u; k++ ) {
	for( int32_t j = 1; j < t; j++ ) {

	int32_t i = s;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/16.0*(*in)( ii,   jj-1, kk-1 ) + 
	    1.0/8.0 *(*in)( ii,   jj,   kk-1 ) + 
	    1.0/16.0*(*in)( ii,   jj+1, kk-1 ) + 

	    1.0/8.0 *(*in)( ii,   jj-1, kk   ) + 
	    1.0/4.0 *(*in)( ii,   jj,   kk   ) + 
	    1.0/8.0 *(*in)( ii,   jj+1, kk   ) + 

	    1.0/16.0*(*in)( ii,   jj-1, kk+1 ) + 
	    1.0/8.0 *(*in)( ii,   jj,   kk+1 ) + 
	    1.0/16.0*(*in)( ii,   jj+1, kk+1 );
	
	(*out)(i,j,k) = D;
	}
    }

    // j=0 boundary nodes of level+1
    for( int32_t k = 1; k < u; k++ ) {
	for( int32_t i = 1; i < s; i++ ) {

	int32_t j = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/16.0*(*in)( ii-1, jj,   kk-1 ) + 
	    1.0/8.0 *(*in)( ii,   jj,   kk-1 ) + 
	    1.0/16.0*(*in)( ii+1, jj,   kk-1 ) + 

	    1.0/8.0 *(*in)( ii-1, jj,   kk   ) + 
	    1.0/4.0 *(*in)( ii,   jj,   kk   ) + 
	    1.0/8.0 *(*in)( ii+1, jj,   kk   ) + 

	    1.0/16.0*(*in)( ii-1, jj,   kk+1 ) + 
	    1.0/8.0 *(*in)( ii,   jj,   kk+1 ) + 
	    1.0/16.0*(*in)( ii+1, jj,   kk+1 );
	
	(*out)(i,j,k) = D;
	}
    }

    // j=size(1)-1 boundary nodes of level+1
    for( int32_t k = 1; k < u; k++ ) {
	for( int32_t i = 1; i < s; i++ ) {

	int32_t j = t;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/16.0*(*in)( ii-1, jj,   kk-1 ) + 
	    1.0/8.0 *(*in)( ii,   jj,   kk-1 ) + 
	    1.0/16.0*(*in)( ii+1, jj,   kk-1 ) + 

	    1.0/8.0 *(*in)( ii-1, jj,   kk   ) + 
	    1.0/4.0 *(*in)( ii,   jj,   kk   ) + 
	    1.0/8.0 *(*in)( ii+1, jj,   kk   ) + 

	    1.0/16.0*(*in)( ii-1, jj,   kk+1 ) + 
	    1.0/8.0 *(*in)( ii,   jj,   kk+1 ) + 
	    1.0/16.0*(*in)( ii+1, jj,   kk+1 );
	
	(*out)(i,j,k) = D;
	}
    }

    // k=0 boundary nodes of level+1
    for( int32_t j = 1; j < t; j++ ) {
	for( int32_t i = 1; i < s; i++ ) {

	int32_t k = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/16.0*(*in)( ii-1, jj-1, kk   ) + 
	    1.0/8.0 *(*in)( ii,   jj-1, kk   ) + 
	    1.0/16.0*(*in)( ii+1, jj-1, kk   ) + 

	    1.0/8.0 *(*in)( ii-1, jj,   kk   ) + 
	    1.0/4.0 *(*in)( ii,   jj,   kk   ) + 
	    1.0/8.0 *(*in)( ii+1, jj,   kk   ) + 

	    1.0/16.0*(*in)( ii-1, jj+1, kk   ) + 
	    1.0/8.0 *(*in)( ii,   jj+1, kk   ) + 
	    1.0/16.0*(*in)( ii+1, jj+1, kk   );
	
	(*out)(i,j,k) = D;
	}
    }

    // k=size(2)-1 boundary nodes of level+1
    for( int32_t j = 1; j < t; j++ ) {
	for( int32_t i = 1; i < s; i++ ) {

	int32_t k = u;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/16.0*(*in)( ii-1, jj-1, kk   ) + 
	    1.0/8.0 *(*in)( ii,   jj-1, kk   ) + 
	    1.0/16.0*(*in)( ii+1, jj-1, kk   ) + 

	    1.0/8.0 *(*in)( ii-1, jj,   kk   ) + 
	    1.0/4.0 *(*in)( ii,   jj,   kk   ) + 
	    1.0/8.0 *(*in)( ii+1, jj,   kk   ) + 

	    1.0/16.0*(*in)( ii-1, jj+1, kk   ) + 
	    1.0/8.0 *(*in)( ii,   jj+1, kk   ) + 
	    1.0/16.0*(*in)( ii+1, jj+1, kk   );
	
	(*out)(i,j,k) = D;
	}
    }

    // EDGES

    // i = 0, j = 0
    for( int32_t k = 1; k < u; k++ ) {

	int32_t i = 0;
	int32_t j = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/4.0*(*in)( ii,   jj,   kk-1 ) + 
	    1.0/2.0*(*in)( ii,   jj,   kk   ) + 
	    1.0/4.0*(*in)( ii,   jj,   kk+1 );
	
	(*out)(i,j,k) = D;
    }

    // i = 0, j = t
    for( int32_t k = 1; k < u; k++ ) {

	int32_t i = 0;
	int32_t j = t;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/4.0*(*in)( ii,   jj,   kk-1 ) + 
	    1.0/2.0*(*in)( ii,   jj,   kk   ) + 
	    1.0/4.0*(*in)( ii,   jj,   kk+1 );
	
	(*out)(i,j,k) = D;
    }

    // i = s, j = 0
    for( int32_t k = 1; k < u; k++ ) {

	int32_t i = s;
	int32_t j = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/4.0*(*in)( ii,   jj,   kk-1 ) + 
	    1.0/2.0*(*in)( ii,   jj,   kk   ) + 
	    1.0/4.0*(*in)( ii,   jj,   kk+1 );
	
	(*out)(i,j,k) = D;
    }

    // i = s, j = t
    for( int32_t k = 1; k < u; k++ ) {

	int32_t i = s;
	int32_t j = t;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/4.0*(*in)( ii,   jj,   kk-1 ) + 
	    1.0/2.0*(*in)( ii,   jj,   kk   ) + 
	    1.0/4.0*(*in)( ii,   jj,   kk+1 );
	
	(*out)(i,j,k) = D;
    }

    // i = 0, k = 0
    for( int32_t j = 1; j < t; j++ ) {

	int32_t i = 0;
	int32_t k = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/4.0*(*in)( ii,   jj-1, kk   ) + 
	    1.0/2.0*(*in)( ii,   jj,   kk   ) + 
	    1.0/4.0*(*in)( ii,   jj+1, kk   );
	
	(*out)(i,j,k) = D;
    }

    // i = 0, k = u
    for( int32_t j = 1; j < t; j++ ) {

	int32_t i = 0;
	int32_t k = u;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/4.0*(*in)( ii,   jj-1, kk   ) + 
	    1.0/2.0*(*in)( ii,   jj,   kk   ) + 
	    1.0/4.0*(*in)( ii,   jj+1, kk   );
	
	(*out)(i,j,k) = D;
    }

    // i = s, k = 0
    for( int32_t j = 1; j < t; j++ ) {

	int32_t i = s;
	int32_t k = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/4.0*(*in)( ii,   jj-1, kk   ) + 
	    1.0/2.0*(*in)( ii,   jj,   kk   ) + 
	    1.0/4.0*(*in)( ii,   jj+1, kk   );
	
	(*out)(i,j,k) = D;
    }

    // i = s, k = u
    for( int32_t j = 1; j < t; j++ ) {

	int32_t i = s;
	int32_t k = u;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/4.0*(*in)( ii,   jj-1, kk   ) + 
	    1.0/2.0*(*in)( ii,   jj,   kk   ) + 
	    1.0/4.0*(*in)( ii,   jj+1, kk   );
	
	(*out)(i,j,k) = D;
    }

    // j = 0, k = 0
    for( int32_t i = 1; i < s; i++ ) {

	int32_t j = 0;
	int32_t k = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/4.0*(*in)( ii-1, jj,   kk   ) + 
	    1.0/2.0*(*in)( ii,   jj,   kk   ) + 
	    1.0/4.0*(*in)( ii+1, jj,   kk   );
	
	(*out)(i,j,k) = D;
    }

    // j = 0, k = u
    for( int32_t i = 1; i < s; i++ ) {

	int32_t j = 0;
	int32_t k = u;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/4.0*(*in)( ii-1, jj,   kk   ) + 
	    1.0/2.0*(*in)( ii,   jj,   kk   ) + 
	    1.0/4.0*(*in)( ii+1, jj,   kk   );
	
	(*out)(i,j,k) = D;
    }

    // j = t, k = 0
    for( int32_t i = 1; i < s; i++ ) {

	int32_t j = t;
	int32_t k = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/4.0*(*in)( ii-1, jj,   kk   ) + 
	    1.0/2.0*(*in)( ii,   jj,   kk   ) + 
	    1.0/4.0*(*in)( ii+1, jj,   kk   );
	
	(*out)(i,j,k) = D;
    }

    // j = t, k = u
    for( int32_t i = 1; i < s; i++ ) {

	int32_t j = t;
	int32_t k = u;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	int32_t kk = 2*k;
	double D = 
	    1.0/4.0*(*in)( ii-1, jj,   kk   ) + 
	    1.0/2.0*(*in)( ii,   jj,   kk   ) + 
	    1.0/4.0*(*in)( ii+1, jj,   kk   );
	
	(*out)(i,j,k) = D;
    }

    // CORNERS

    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    double D = (*in)( 2*i, 2*j, 2*k );
    (*out)(i,j,k) = D;

    i = s;
    D = (*in)( 2*i, 2*j, 2*k );
    (*out)(i,j,k) = D;

    j = t;
    D = (*in)( 2*i, 2*j, 2*k );
    (*out)(i,j,k) = D;

    i = 0;
    D = (*in)( 2*i, 2*j, 2*k );
    (*out)(i,j,k) = D;

    k = u;
    D = (*in)( 2*i, 2*j, 2*k );
    (*out)(i,j,k) = D;

    i = s;
    D = (*in)( 2*i, 2*j, 2*k );
    (*out)(i,j,k) = D;

    j = 0;
    D = (*in)( 2*i, 2*j, 2*k );
    (*out)(i,j,k) = D;

    i = 0;
    D = (*in)( 2*i, 2*j, 2*k );
    (*out)(i,j,k) = D;
}


void EpotMGSolver::restrict_cyl( MeshScalarField *out, const MeshScalarField *in, bool defect )
{
    // Go through internal nodes of rougher level (out)
    int32_t s = out->size(0)-1;
    int32_t t = out->size(1)-1;
    for( int32_t j = 1; j < t; j++ ) {
        for( int32_t i = 1; i < s; i++ ) {

	    int32_t ii = 2*i;
	    int32_t jj = 2*j;

	    // RBGS smoother does even pass last. Even nodes have
	    // zero defect.
	    if( defect ) {
		double D = 
		    1.0/2.0 *(*in)( ii-1, jj   ) + 
		    1.0/2.0 *(*in)( ii,   jj-1 ) + 
		    1.0/2.0 *(*in)( ii,   jj+1 ) + 
		    1.0/2.0 *(*in)( ii+1, jj   );
		(*out)(i,j) = D;
	    } else {
		double D = 
		    1.0/16.0*(*in)( ii-1, jj-1 ) + 
		    1.0/8.0 *(*in)( ii-1, jj   ) + 
		    1.0/16.0*(*in)( ii-1, jj+1 ) + 
		    
		    1.0/8.0 *(*in)( ii,   jj-1 ) + 
		    1.0/4.0 *(*in)( ii,   jj   ) + 
		    1.0/8.0 *(*in)( ii,   jj+1 ) + 
		    
		    1.0/16.0*(*in)( ii+1, jj-1 ) + 
		    1.0/8.0 *(*in)( ii+1, jj   ) + 
		    1.0/16.0*(*in)( ii+1, jj+1 );
		(*out)(i,j) = D;
	    }
	}
    }

    // i=0 boundary nodes of level+1
    for( int32_t j = 1; j < t; j++ ) {

	int32_t i = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*in)( ii, jj-1 ) + 
	    1.0/2.0*(*in)( ii, jj   ) + 
	    1.0/4.0*(*in)( ii, jj+1 );
	
	(*out)(i,j) = D;
    }

    // i=size(0)-1 boundary nodes of level+1
    for( int32_t j = 1; j < t; j++ ) {

	int32_t i = s;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*in)( ii, jj-1 ) + 
	    1.0/2.0*(*in)( ii, jj   ) + 
	    1.0/4.0*(*in)( ii, jj+1 );
	
	(*out)(i,j) = D;
    }

    // j=0 boundary nodes of level+1 (axis)
    for( int32_t i = 1; i < s; i++ ) {

	int32_t j = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*in)( ii-1, jj ) + 
	    1.0/2.0*(*in)( ii,   jj ) + 
	    1.0/4.0*(*in)( ii+1, jj );
	
	(*out)(i,j) = D;
    }

    // j=size(1)-1 boundary nodes of level+1
    for( int32_t i = 1; i < s; i++ ) {

	int32_t j = t;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*in)( ii-1, jj ) + 
	    1.0/2.0*(*in)( ii,   jj ) + 
	    1.0/4.0*(*in)( ii+1, jj );
	
	(*out)(i,j) = D;
    }

    // And corners
    int32_t i = 0;
    int32_t j = 0;
    double D = (*in)( 2*i, 2*j );
    (*out)(i,j) = D;

    i = s;
    D = (*in)( 2*i, 2*j );
    (*out)(i,j) = D;

    j = t;
    D = (*in)( 2*i, 2*j );
    (*out)(i,j) = D;

    i = 0;
    D = (*in)( 2*i, 2*j );
    (*out)(i,j) = D;
}


void EpotMGSolver::restrict_2d( MeshScalarField *out, const MeshScalarField *in, bool defect )
{
    // Go through internal nodes of rougher level (out)
    int32_t s = out->size(0)-1;
    int32_t t = out->size(1)-1;
    for( int32_t j = 1; j < t; j++ ) {
        for( int32_t i = 1; i < s; i++ ) {

	    int32_t ii = 2*i;
	    int32_t jj = 2*j;

	    // RBGS smoother does even pass last. Even nodes have
	    // zero defect.
	    if( defect ) {
		double D = 
		    1.0/2.0 *(*in)( ii-1, jj   ) + 
		    1.0/2.0 *(*in)( ii,   jj-1 ) + 
		    1.0/2.0 *(*in)( ii,   jj+1 ) + 
		    1.0/2.0 *(*in)( ii+1, jj   );
		(*out)(i,j) = D;
	    } else {
		double D = 
		    1.0/16.0*(*in)( ii-1, jj-1 ) + 
		    1.0/8.0 *(*in)( ii-1, jj   ) + 
		    1.0/16.0*(*in)( ii-1, jj+1 ) + 
		    
		    1.0/8.0 *(*in)( ii,   jj-1 ) + 
		    1.0/4.0 *(*in)( ii,   jj   ) + 
		    1.0/8.0 *(*in)( ii,   jj+1 ) + 
		    
		    1.0/16.0*(*in)( ii+1, jj-1 ) + 
		    1.0/8.0 *(*in)( ii+1, jj   ) + 
		    1.0/16.0*(*in)( ii+1, jj+1 );
		(*out)(i,j) = D;
	    }
	}
    }

    // i=0 boundary nodes of level+1
    for( int32_t j = 1; j < t; j++ ) {

	int32_t i = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*in)( ii, jj-1 ) + 
	    1.0/2.0*(*in)( ii, jj   ) + 
	    1.0/4.0*(*in)( ii, jj+1 );
	
	(*out)(i,j) = D;
    }

    // i=size(0)-1 boundary nodes of level+1
    for( int32_t j = 1; j < t; j++ ) {

	int32_t i = s;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*in)( ii, jj-1 ) + 
	    1.0/2.0*(*in)( ii, jj   ) + 
	    1.0/4.0*(*in)( ii, jj+1 );
	
	(*out)(i,j) = D;
    }

    // j=0 boundary nodes of level+1
    for( int32_t i = 1; i < s; i++ ) {

	int32_t j = 0;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*in)( ii-1, jj ) + 
	    1.0/2.0*(*in)( ii,   jj ) + 
	    1.0/4.0*(*in)( ii+1, jj );
	
	(*out)(i,j) = D;
    }

    // j=size(1)-1 boundary nodes of level+1
    for( int32_t i = 1; i < s; i++ ) {

	int32_t j = t;
	int32_t ii = 2*i;
	int32_t jj = 2*j;
	double D = 
	    1.0/4.0*(*in)( ii-1, jj ) + 
	    1.0/2.0*(*in)( ii,   jj ) + 
	    1.0/4.0*(*in)( ii+1, jj );
	
	(*out)(i,j) = D;
    }

    // And corners
    int32_t i = 0;
    int32_t j = 0;
    double D = (*in)( 2*i, 2*j );
    (*out)(i,j) = D;

    i = s;
    D = (*in)( 2*i, 2*j );
    (*out)(i,j) = D;

    j = t;
    D = (*in)( 2*i, 2*j );
    (*out)(i,j) = D;

    i = 0;
    D = (*in)( 2*i, 2*j );
    (*out)(i,j) = D;
}


void EpotMGSolver::restrict_1d( MeshScalarField *out, const MeshScalarField *in, bool defect )
{
    // Go through internal nodes of rougher level (out)
    int32_t s = out->size(0)-1;
    for( int32_t i = 1; i < s; i++ ) {

	int32_t ii = 2*i;

	// RBGS smoother does even pass last. Even nodes have
	// zero defect.
	if( defect ) {
	    double D = 
		(*in)( ii-1 ) + 
		(*in)( ii+1 );
	    (*out)(i) = D;
	} else {
	    double D = 
		1.0/4.0*(*in)( ii-1 ) + 
		1.0/2.0*(*in)( ii   ) + 
		1.0/4.0*(*in)( ii+1 );
	    (*out)(i) = D;
	}
    }

    // Boundary nodes completely separated
    (*out)(0) = (*in)(0);
    (*out)(s) = (*in)(2*s);
}


void EpotMGSolver::restrict( MeshScalarField *out, const MeshScalarField *in, bool defect )
{
    switch( _geom.geom_mode() ) {
    case MODE_1D:
	restrict_1d( out, in, defect );
	break;
    case MODE_2D:
	restrict_2d( out, in, defect );
	break;
    case MODE_CYL:
	restrict_cyl( out, in, defect );
	break;
    case MODE_3D:
	restrict_3d( out, in, defect );
	break;
    default:
	break;
    }
}


void EpotMGSolver::prolong_add_3d( MeshScalarField *out, int32_t i, int32_t j, int32_t k, double C )
{
    if( i >= 0 && i < (int32_t)out->size(0) &&
	j >= 0 && j < (int32_t)out->size(1) &&
	k >= 0 && k < (int32_t)out->size(2) )
	(*out)(i,j,k) += C;
}


void EpotMGSolver::prolong_3d( MeshScalarField *out, const MeshScalarField *in )
{
    // Clear output field
    out->clear();

    // Loop through all input nodes
    int32_t s = in->size(0);
    int32_t t = in->size(1);
    int32_t u = in->size(2);
    for( int32_t k = 0; k < u; k++ ) {
	for( int32_t j = 0; j < t; j++ ) {
	    for( int32_t i = 0; i < s; i++ ) {

		int32_t ii = 2*i;
		int32_t jj = 2*j;
		int32_t kk = 2*k;
		double C = (*in)(i,j,k);

		// No need to prolong correction to odd points. Odd
		// field will be overwritten by first RBGS loop.
		prolong_add_3d( out, ii-1, jj-1, kk,   1.0/4.0*C );
		prolong_add_3d( out, ii-1, jj,   kk-1, 1.0/4.0*C );
		prolong_add_3d( out, ii-1, jj,   kk+1, 1.0/4.0*C );
		prolong_add_3d( out, ii-1, jj+1, kk,   1.0/4.0*C );

		prolong_add_3d( out, ii,   jj-1, kk-1, 1.0/4.0*C );
		prolong_add_3d( out, ii,   jj-1, kk+1, 1.0/4.0*C );
		prolong_add_3d( out, ii,   jj,   kk,   1.0*C );
		prolong_add_3d( out, ii,   jj+1, kk-1, 1.0/4.0*C );
		prolong_add_3d( out, ii,   jj+1, kk+1, 1.0/4.0*C );

		prolong_add_3d( out, ii+1, jj-1, kk,   1.0/4.0*C );
		prolong_add_3d( out, ii+1, jj,   kk-1, 1.0/4.0*C );
		prolong_add_3d( out, ii+1, jj,   kk+1, 1.0/4.0*C );
		prolong_add_3d( out, ii+1, jj+1, kk,   1.0/4.0*C );
	    }
	}
    }    
}


void EpotMGSolver::prolong_add_cyl( MeshScalarField *out, int32_t i, int32_t j, double C )
{
    if( i >= 0 && i < (int32_t)out->size(0) &&
	j >= 0 && j < (int32_t)out->size(1) )
	(*out)(i,j) += C;
}


void EpotMGSolver::prolong_cyl( MeshScalarField *out, const MeshScalarField *in )
{
    // Clear output field
    out->clear();

    // Loop through all input nodes
    int32_t s = in->size(0);
    int32_t t = in->size(1);
    for( int32_t j = 0; j < t; j++ ) {
	for( int32_t i = 0; i < s; i++ ) {

	    int32_t ii = 2*i;
	    int32_t jj = 2*j;
	    double C = (*in)(i,j);

	    // No need to prolong correction to odd points. Odd
	    // field will be overwritten by first RBGS loop.
	    prolong_add_cyl( out, ii-1, jj-1, 1.0/4.0*C );
	    //prolong_add_cyl( out, ii-1, jj,   1.0/2.0*C );
	    prolong_add_cyl( out, ii-1, jj+1, 1.0/4.0*C );

	    //prolong_add_cyl( out, ii,   jj-1, 1.0/2.0*C );
	    prolong_add_cyl( out, ii,   jj,   1.0*C );
	    //prolong_add_cyl( out, ii,   jj+1, 1.0/2.0*C );

	    prolong_add_cyl( out, ii+1, jj-1, 1.0/4.0*C );
	    //prolong_add_cyl( out, ii+1, jj,   1.0/2.0*C );
	    prolong_add_cyl( out, ii+1, jj+1, 1.0/4.0*C );
	}
    }    
}


void EpotMGSolver::prolong_add_2d( MeshScalarField *out, int32_t i, int32_t j, double C )
{
    if( i >= 0 && i < (int32_t)out->size(0) &&
	j >= 0 && j < (int32_t)out->size(1) )
	(*out)(i,j) += C;
}


void EpotMGSolver::prolong_2d( MeshScalarField *out, const MeshScalarField *in )
{
    // Clear output field
    out->clear();

    // Loop through all input nodes
    int32_t s = in->size(0);
    int32_t t = in->size(1);
    for( int32_t j = 0; j < t; j++ ) {
	for( int32_t i = 0; i < s; i++ ) {

	    int32_t ii = 2*i;
	    int32_t jj = 2*j;
	    double C = (*in)(i,j);

	    // No need to prolong correction to odd points. Odd
	    // field will be overwritten by first RBGS loop.
	    prolong_add_2d( out, ii-1, jj-1, 1.0/4.0*C );
	    //prolong_add_2d( out, ii-1, jj,   1.0/2.0*C );
	    prolong_add_2d( out, ii-1, jj+1, 1.0/4.0*C );

	    //prolong_add_2d( out, ii,   jj-1, 1.0/2.0*C );
	    prolong_add_2d( out, ii,   jj,   1.0*C );
	    //prolong_add_2d( out, ii,   jj+1, 1.0/2.0*C );

	    prolong_add_2d( out, ii+1, jj-1, 1.0/4.0*C );
	    //prolong_add_2d( out, ii+1, jj,   1.0/2.0*C );
	    prolong_add_2d( out, ii+1, jj+1, 1.0/4.0*C );
	}
    }    
}


void EpotMGSolver::prolong_add_1d( MeshScalarField *out, int32_t i, double C )
{
    if( i >= 0 && i < (int32_t)out->size(0) )
	(*out)(i) += C;
}


void EpotMGSolver::prolong_1d( MeshScalarField *out, const MeshScalarField *in )
{
    // Clear output field
    out->clear();

    // Loop through all input nodes
    int32_t s = in->size(0);
    for( int32_t i = 0; i < s; i++ ) {
	
	int32_t ii = 2*i;
	double C = (*in)(i);

	// No need to prolong correction to odd points. Odd
	// field will be overwritten by first RBGS loop.
	//prolong_add_1d( out, ii-1, 1.0/2.0*C );
	prolong_add_1d( out, ii,       1.0*C );
	//prolong_add_1d( out, ii+1, 1.0/2.0*C );
    }    
}


void EpotMGSolver::prolong( MeshScalarField *out, const MeshScalarField *in )
{
    switch( _geom.geom_mode() ) {
    case MODE_1D:
	prolong_1d( out, in );
	break;
    case MODE_2D:
	prolong_2d( out, in );
	break;
    case MODE_CYL:
	prolong_cyl( out, in );
	break;
    case MODE_3D:
	prolong_3d( out, in );
	break;
    default:
	break;
    }
}


// Make correction Xnew=X+V
void EpotMGSolver::correct( const Geometry *geom, MeshScalarField *sol, const MeshScalarField *corr )
{
    // Loop through all nodes, only correct non-fixed even points.
    // Odd field will be overwritten by first RBGS loop.
    int32_t s = geom->size(0);
    int32_t t = geom->size(1);
    int32_t u = geom->size(2);
    for( int32_t k = 0; k < u; k++ ) {
	for( int32_t j = 0; j < t; j++ ) {
	    for( int32_t i = (k+j) % 2; i < s; i+=2 ) {

		if( !(geom->mesh(i,j,k) & SMESH_NODE_FIXED) )
		    (*sol)(i,j,k) += (*corr)(i,j,k);
	    }
	}
    }
}


void EpotMGSolver::mg_recurse( uint32_t level )
{
    ibsimu.inc_indent();

    // Linear case: clear last field, expected result closer to zero than result 
    // from last round. Never clear top level field.
    if( linear() && level != 0 )
	_epotv[level]->clear();

    if( level == _levels-1 ) {
 
	ibsimu.message( 2 ) << "Roughest level solution\n";
	
	// Last level, solve the roughest problem until convergence
	uint32_t a = 0;
	while( a < _imax ) {
	    _step = _epotsolverv[level]->mg_solve( _epotv[level], _rhsv[level], _w );
	    _coarse_err = _epotsolverv[level]->error_scale( _w ) * _step;
	    a++;
	    if( _coarse_err < _coarse_eps )
		break;
	}
	_coarse_steps += a;

	if( ibsimu.get_message_threshold( MSG_VERBOSE ) > 1 ) {
	    ibsimu.message( 2 ) << a << " iterations done\n";
	    ibsimu.message( 2 ) << _coarse_err << " accuracy reached\n";
	    if( a >= _imax ) {
		ibsimu.message( 2 ) << "maximum number of iterations done\n";
	    }
	}

	ibsimu.dec_indent();
	return;
    } 

    // Do gamma cycles of next level
    for( uint32_t gcyc = 0; gcyc < _gamma; gcyc++ ) {

	if( _gamma > 1 )
	    ibsimu.message( 2 ) << "Doing gamma cycle " << gcyc+1 << "/" 
				<< _gamma << " at level " << level << "\n";
	else
	    ibsimu.message( 2 ) << "Doing cycle at level " << level << "\n";

	// Pre smoothing
	for( uint32_t a = 0; a < _npre; a++ )
	    _epotsolverv[level]->mg_smooth( _epotv[level], _rhsv[level] );

	// Defect
	_epotsolverv[level]->defect( _workv[level], _epotv[level], _rhsv[level], true );

	if( linear() ) {

	    restrict( _rhsv[level+1], _workv[level], true );
	    (*_rhsv[level+1]) *= -1.0;

	} else {

	    restrict( _workv[level+1], _workv[level], true );
	    restrict( _epotv[level+1], _epotv[level], false );

	    // Store solution
	    (*_work2v[level+1]) = (*_epotv[level+1]);

	    // Calculate defect
	    _epotsolverv[level+1]->defect( _rhsv[level+1], _epotv[level+1], _workv[level+1], false );
	}

	// Recurse to next level to solve for correction V (epot)
	// A*V=-D
	mg_recurse( level+1 );
	
	if( linear() ) {

	    // Prolong correction (epot) to work field
	    prolong( _workv[level], _epotv[level+1] );
	    
	} else {

	    // Compute correction
	    (*_workv[level+1]) = (*_epotv[level+1]);
	    (*_workv[level+1]) -= (*_work2v[level+1]);
	    prolong( _workv[level], _workv[level+1] );
	}

	// Make correction Xnew=X+V
	correct( _geomv[level], _epotv[level], _workv[level] );

	// Post smoothing
	for( uint32_t a = 0; a < _npost; a++ )
	    _step = _epotsolverv[level]->mg_smooth( _epotv[level], _rhsv[level] );
	_mgerr = _epotsolverv[level]->error_scale_mg() * _step;

	if( ibsimu.get_message_threshold( MSG_VERBOSE ) > 1 )
	    ibsimu.message( 2 ) << _mgerr << " accuracy reached\n";
     }

    ibsimu.dec_indent();
}


void EpotMGSolver::print_field( const MeshScalarField *F )
{
    for( size_t i = 0; i < F->size(0); i++ ) {
	ibsimu.message( 1 ) << std::setw(8) << (*F)(i) << " ";
    }
    ibsimu.message( 1 ) << "\n";
}


void EpotMGSolver::prepare_local_gnewton_settings( void )
{
    if( _plasma == PLASMA_PEXP ) {
	// Limit calculation to 10 times electron temp from plasma potential
	_local_Ulim = _Up - 10.0*_Te;
	ibsimu.message( 1 ) << "Limiting plasma calculation to U > " << _local_Ulim << " V\n";
    } else if( _plasma == PLASMA_NSIMP ) {
	_local_Ulim = 0.0;
	for( uint32_t a = 0; a < _Ei.size(); a++ ) {
	    if( _Ei[a] > _local_Ulim )
		_local_Ulim = _Ei[a];
	}
	// Limit calculation to 10 times maximum energy
	_local_Ulim = 10.0*_local_Ulim; 
	ibsimu.message( 1 ) << "Limiting plasma calculation to U < " << _local_Ulim << " V\n";
    } else if( _plasma == PLASMA_SHIELD ) {
	_local_Ulim = fabs(_Up)+fabs(10.0*_Te);
	ibsimu.message( 1 ) << "Limiting plasma calculation to |U| < " << fabs(_local_Ulim) << " V\n";
    }
}


void EpotMGSolver::subsolve( MeshScalarField &epot, const MeshScalarField &scharge )
{
    //ibsimu.message( 1 ) << "\n";

    ibsimu.message( 1 ) << "Using Multigrid solver (" 
			<< "levels = " << _levels
			<< ", npre = " << _npre
			<< ", npost = " << _npost
			<< ", mgcycmax = " << _mgcycmax
			<< ", mgeps = " << _mgeps
			<< ", gamma = " << _gamma
			<< ", w = " << _w
			<< ", coarse_eps = " << _coarse_eps
			<< ", imax = " << _imax
			<< ")\n";
    if( linear() )
	ibsimu.message( 1 ) << "Linear problem\n";
    else {
	ibsimu.message( 1 ) << "Nonlinear problem\n";
	prepare_local_gnewton_settings();
    }

    if( !_geom_prepared )
	prepare_mg_geom();
    preprocess( epot, scharge );

    // Status print
    StatusPrint sp;

    // Call iterator mgcyc times
    uint32_t a;
    ibsimu.message( 1 ) << "Taking multigrid cycles\n";
    for( a = 0; a < _mgcycmax; a++ ) {
	if( ibsimu.get_message_threshold( MSG_VERBOSE ) > 1 )
	    ibsimu.message( 2 ) << "Cycle " << a+1 << "\n";

	_coarse_steps = 0;
	mg_recurse( 0 );

	// Print statistics
	if( ibsimu.get_message_threshold( MSG_VERBOSE ) == 1 ) {
	    std::stringstream ss;
	    ss << std::setw(5) << a+1 << ", "
	       << "rough " << _coarse_steps << " " << _coarse_err << ", "
	       << "fine " << _mgerr << "      ";
	    sp.print( ss.str(), true );
	}

	if( _mgerr < _mgeps )
	    break;
    }
    ibsimu.message( 1 ) << "\nDone\n";
    _mgcyc = a;

    if( a == _mgcycmax )
	ibsimu.message( 1 ) << "Maximum number of cycles done\n";
    else
	ibsimu.message( 1 ) << "Convergence criterion reached\n";

    postprocess();
}


void EpotMGSolver::save( std::ostream &s ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}




void EpotMGSolver::debug_print( std::ostream &os ) const 
{
    EpotSolver::debug_print( os );
    os << "**EpotMGSolver\n";
    os << "geom_prepared = " << _geom_prepared << "\n";
    os << "levels = " << _levels << "\n";
    os << "npre = " << _npre << "\n";
    os << "npost = " << _npost << "\n";
    os << "mgcycmax = " << _mgcycmax << "\n";
    os << "mgcyc = " << _mgcyc << "\n";
    os << "mgeps = " << _mgeps << "\n";
    os << "mgerr = " << _mgerr << "\n";
    os << "gamma = " << _gamma << "\n";
    os << "step = " << _step << "\n";
    os << "coarse_eps = " << _coarse_eps << "\n";
    os << "coarse_err = " << _coarse_err << "\n";
    os << "w = " << _w << "\n";
    os << "imax = " << _imax << "\n";
    os << "coarse_steps = " << _coarse_steps << "\n";
}

