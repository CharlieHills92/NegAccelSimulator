/*! \file particleiterator.hpp
 *  \brief %Particle iterator
 */

/* Copyright (c) 2005-2013,2018 Taneli Kalvas, Tobin Jones. All rights reserved.
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

#ifndef PARTICLEITERATOR_HPP
#define PARTICLEITERATOR_HPP 1


#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_poly.h>
#include "geometry.hpp"
#include "mat3d.hpp"
#include "compmath.hpp"
#include "trajectory.hpp"
#include "particles.hpp"
#include "vectorfield.hpp"
#include "meshscalarfield.hpp"
#include "scharge.hpp"
#include "scheduler.hpp"
#include "polysolver.hpp"
#include "particledatabase.hpp"
#include "cfifo.hpp"
#include "ibsimu.hpp"


//#define DEBUG_PARTICLE_ITERATOR 1


#ifdef DEBUG_PARTICLE_ITERATOR 
#define DEBUG_MESSAGE(x) ibsimu.message(MSG_DEBUG_GENERAL,1) << x
#define DEBUG_INC_INDENT() ibsimu.inc_indent()
#define DEBUG_DEC_INDENT() ibsimu.dec_indent()
#else
#define DEBUG_MESSAGE(x) do {} while(0)
#define DEBUG_INC_INDENT() do {} while(0)
#define DEBUG_DEC_INDENT() do {} while(0)
#endif


#define COLLISION_EPS 1.0e-6


/*! \brief %Particle iterator type.
 */
enum particle_iterator_type_e {
    PARTICLE_ITERATOR_ADAPTIVE = 0,
    PARTICLE_ITERATOR_FIXED_STEP_LEN
};



/*! \brief %Mesh intersection (collision) coordinate data
 *
 * Contains one coordinate data and the direction of particle for one
 * intersection between particle trajectory and mesh plane. Templated
 * for particle point type (see ParticlePBase).
 *
 * Templated for particle point type (see ParticlePBase).
 */
template <class PP> class ColData {
public:
    PP                _x;         /*!< \brief %Mesh intersection coordinates. */
    int               _dir;       /*!< \brief Direction of particle at intersection.
				   *  i: -1/+1, j: -2/+2, k: -3:/+3. */
    
    /*! \brief Default constructor.
     */
    ColData() : _dir(0) {}

    /*! \brief Constructor for collision at \a x into direction \a dir.
     */
    ColData( PP x, int dir ) : _x(x), _dir(dir) {}
    
    /*! \brief Compare coldata entry times.
     *
     *  Used for sorting coldata entries.
     */
    bool operator<( const ColData &cd ) const {
	return( _x[0] < cd._x[0] );
    }

    /*! \brief Find mesh intersections of linearly interpolated
     *  particle trajectory segment.
     *
     *  Makes a linear interpolation between points \a x1 and \a x2
     *  and searches intersection points of this line and \a
     *  mesh. Intersection points are saved to vector \a coldata in
     *  increasing time order.
     */
    static void build_coldata_linear( std::vector<ColData> &coldata, const Mesh &mesh,
				      const PP &x1, const PP &x2 ) {
	
	DEBUG_MESSAGE( "Building coldata using linear interpolation\n" );
	DEBUG_INC_INDENT();

	coldata.clear();

	for( size_t a = 0; a < PP::dim(); a++ ) {
	    
            int a1 = (int)floor( (x1[2*a+1]-mesh.origo(a))/mesh.h() );
            int a2 = (int)floor( (x2[2*a+1]-mesh.origo(a))/mesh.h() );
            if( a1 > a2 ) {
                int a = a2;
                a2 = a1;
                a1 = a;
            }

            for( int b = a1+1; b <= a2; b++ ) {

                // Save intersection coordinates
                double K = (b*mesh.h() + mesh.origo(a) - x1[2*a+1]) / 
                    (x2[2*a+1] - x1[2*a+1]);
                if( K < 0.0 ) K = 0.0;
                else if( K > 1.0 ) K = 1.0;

		DEBUG_MESSAGE( "Adding point " << x1 + (x2-x1)*K << "\n" );

                if( x2[2*a+1] > x1[2*a+1] )
                    coldata.push_back( ColData( x1 + (x2-x1)*K, a+1 ) );
                else
                    coldata.push_back( ColData( x1 + (x2-x1)*K, -a-1 ) );
            }
        }

	// Sort intersections in increasing time order
	sort( coldata.begin(), coldata.end() );

	DEBUG_DEC_INDENT();
	DEBUG_MESSAGE( "Coldata built\n" );
    }

    /*! \brief Find mesh intersections of polynomially interpolated
        particle trajectory segment.
     *
     *  Makes a polynomial interpolation between points \ə x1 and \a
     *  x2 and searches intersection points of this line and \a
     *  mesh. Intersection points are saved to vector \a coldata in
     *  increasing time order.
     */
    static void build_coldata_poly( std::vector<ColData> &coldata, const Mesh &mesh,
				    const PP &x1, const PP &x2 ) {
	
	DEBUG_MESSAGE( "Building coldata using polynomial interpolation\n" );
	DEBUG_INC_INDENT();

	coldata.clear();

	// Construct trajectory representation
	TrajectoryRep1D *traj = new TrajectoryRep1D[PP::dim()];
	for( size_t a = 0; a < PP::dim(); a++ ) {
	    traj[a].construct( x2[0]-x1[0], 
			       x1[2*a+1], x1[2*a+2], 
			       x2[2*a+1], x2[2*a+2] );
	    DEBUG_MESSAGE( "Trajectory polynomial " << a << " order = " 
			   << traj[a].get_representation_order() << "\n" );
	}

	// Solve trajectory intersections
	for( size_t a = 0; a < PP::dim(); a++ ) {

	    // Mesh number of x1 (start point)
	    int i = (int)floor( (x1[2*a+1]-mesh.origo(a))/mesh.h() );
	    
	    // Search to negative (dj = -1) and positive (dj = +1) mesh directions
	    for( int dj = -1; dj <= 1; dj += 2 ) {
		int j = i;
		if( dj == +1 )
		    j = i+1;
		int Kcount;  // Solution counter
		double K[3]; // Solution array
		while( 1 ) {

		    // Intersection point
		    double val = mesh.origo(a) + mesh.h() * j;
		    if( val < mesh.origo(a)-mesh.h() )
			break;
		    else if( val > mesh.max(a)+mesh.h() )
			break;

		    DEBUG_MESSAGE( "Searching intersections at coord(" << a << ") = " << val << "\n" );

		    Kcount = traj[a].solve( K, val );
		    if( Kcount == 0 )
			break; // No valid roots

		    // Save roots to coldata
		    DEBUG_MESSAGE( "Found " << Kcount << " valid roots\n" );
		    DEBUG_INC_INDENT();
		    for( int b = 0; b < Kcount; b++ ) {
			PP xcol;
			double x, v;
			xcol(0) = x1[0] + K[b]*(x2[0]-x1[0]);
			for( size_t c = 0; c < PP::dim(); c++ ) {
			    traj[c].coord( x, v, K[b] );
			    if( a == c )
				xcol[2*c+1] = val; // limit numerical inaccuracy
			    else
				xcol[2*c+1] = x;
			    xcol[2*c+2] = v;
			}
			if( mesh.geom_mode() == MODE_CYL )
			    xcol[5] = x1[5] + K[b]*(x2[5]-x1[5]);

			DEBUG_MESSAGE( "K = " << K[b] << "\n" );
			DEBUG_MESSAGE( "Adding point " << xcol << "\n" );

			if( xcol[2*a+2] >= 0.0 )
			    coldata.push_back( ColData( xcol, a+1 ) );
			else
			    coldata.push_back( ColData( xcol, -a-1 ) );
		    }
		    DEBUG_DEC_INDENT();

		    j += dj;
		}
	    }
	}

	// Sort intersections in increasing time order
	sort( coldata.begin(), coldata.end() );

	delete[] traj;

	DEBUG_DEC_INDENT();
	DEBUG_MESSAGE( "Coldata built\n" );
    }

};


/*! \brief %Particle iterator class for continuous Vlasov-type iteration.
 *
 * Templated for particle point type (see ParticlePBase).
 *
 *  \todo Detailed documentation needed.
 *  \todo PIC style iterator needed.
 */
template <class PP> class ParticleIterator {

    gsl_odeiv_system           _system;        /**< \brief GSL ODE integrator system. */
    gsl_odeiv_step            *_step;          /**< \brief GSL ODE integrator stepper. */
    gsl_odeiv_control         *_control;       /**< \brief GSL ODE integrator constrol. */
    gsl_odeiv_evolve          *_evolve;        /**< \brief GSL ODE integrator integrator. */

    particle_iterator_type_e   _type;          /**< \brief Iteratory type. */

    trajectory_interpolation_e _intrp;         /*!< \brief Interpolation type. */
    scharge_deposition_e       _scharge_dep;   /*!< \brief Space charge deposition type. */
    double                     _epsabs;        /*!< \brief Absolute error limit. */
    double                     _epsrel;        /*!< \brief Relative error limit. */
    uint32_t                   _maxsteps;      /*!< \brief Maximum number of simulation steps for particle. */
    double                     _maxt;          /*!< \brief Maximum particle lifetime. */
    bool                       _save_points;   /*!< \brief Save all points? */
    uint32_t                   _trajdiv;       /*!< \brief Divisor for saved trajectories,
					        * if 3, every third trajectory is saved. */
    bool                       _mirror[6];     /*!< \brief Is particle mirrored on boundary? */
    bool                       _surface_collision;

    ParticleIteratorData       _pidata;        /*!< \brief User data provided to PP::get_derivatives(). */
    TrajectoryHandlerCallback *_thand_cb;      /*!< \brief Trajectory handler callback. */
    TrajectoryEndCallback     *_tend_cb;       /*!< \brief Trajectory end callback. */
    TrajectorySurfaceCollisionCallback *_tsur_cb;   /*!< \brief Trajectory surface collision callback. */
    const TrajectoryEndCallback *_bsup_cb;     /*!< \brief B-field plasma suppression callback. */
    const TrajectoryEndCallback *_esup_cb;     /*!< \brief E-field plasma suppression callback. */
    ParticleDataBase          *_pdb;           /*!< \brief Particle database pointer for adding secondary particles. */
    pthread_mutex_t           *_scharge_mutex; /*!< \brief Space charge mutex. */

    PP                         _xi;            /*!< \brief Previous mesh intersection coordinates 
					        *   or starting point. */
    std::vector<PP>            _traj;          /*!< \brief %Particle trajectory data for current trajectory. */
    std::vector<ColData<PP> >  _coldata;       /*!< \brief %Mesh intersection coordinate data. */
    CFiFo<PP,4>                _cdpast;        /*!< \brief Past three intersection coords. */

    ParticleStatistics         _stat;          /*!< \brief Particle statistics. */


    /*! \brief Save trajectory point \a x.
     *
     *  Throw error if run out of memory.
     */
    void save_trajectory_point( PP x ) {

	try {
	    _traj.push_back( x );
	} catch( std::bad_alloc ) {
	    throw( ErrorNoMem( ERROR_LOCATION, "Out of memory saving trajectory" ) );
	}
    }

    /*! \brief Return solid number from nodes around cube (i,j,k).
     */
    uint32_t get_solid( int i, int j, int k ) {
	uint32_t node;
	if( PP::dim() == 2 ) {
	    node = _pidata._geom->mesh( i, j );
	    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
		(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
		return( node & SMESH_BOUNDARY_NUMBER_MASK );
	    node = _pidata._geom->mesh( i+1, j );
	    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
		(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
		return( node & SMESH_BOUNDARY_NUMBER_MASK );
	    node = _pidata._geom->mesh( i, j+1 );
	    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
		(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
		return( node & SMESH_BOUNDARY_NUMBER_MASK );
	    node = _pidata._geom->mesh( i+1, j+1 );
	    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
		(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
		return( node & SMESH_BOUNDARY_NUMBER_MASK );
	} else {
	    node = _pidata._geom->mesh( i, j, k );
	    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
		(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
		return( node & SMESH_BOUNDARY_NUMBER_MASK );
	    node = _pidata._geom->mesh( i+1, j, k );
	    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
		(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
		return( node & SMESH_BOUNDARY_NUMBER_MASK );
	    node = _pidata._geom->mesh( i, j+1, k );
	    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
		(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
		return( node & SMESH_BOUNDARY_NUMBER_MASK );
	    node = _pidata._geom->mesh( i+1, j+1, k );
	    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
		(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
		return( node & SMESH_BOUNDARY_NUMBER_MASK );
	    node = _pidata._geom->mesh( i, j, k+1 );
	    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
		(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
		return( node & SMESH_BOUNDARY_NUMBER_MASK );
	    node = _pidata._geom->mesh( i+1, j, k+1 );
	    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
		(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
		return( node & SMESH_BOUNDARY_NUMBER_MASK );
	    node = _pidata._geom->mesh( i, j+1, k+1 );
	    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
		(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
		return( node & SMESH_BOUNDARY_NUMBER_MASK );
	    node = _pidata._geom->mesh( i+1, j+1, k+1 );
	    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
		(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
		return( node & SMESH_BOUNDARY_NUMBER_MASK );
	}
	return( 0 );
    }

    /*! \brief Check for particle collision with triangulated solid surface.
     *
     *  Checks for intersection between trajectory segment and surface
     *  triangles in grid cell i.
     */
    bool check_collision_surface( Particle<PP> &particle, const PP &x1, const PP &x2, 
				  PP &status_x, const int32_t i[3] ) {

	DEBUG_MESSAGE( "Checking collisions with solids - surface check\n" );
	DEBUG_INC_INDENT();

	if( i[0] < 0 )
	    return( true );
	else if( i[0] >= (int32_t)_pidata._geom->size(0)-1 )
	    return( true );
	else if( i[1] < 0 )
	    return( true );
	else if( i[1] >= (int32_t)_pidata._geom->size(1)-1 )
	    return( true );
	else if( i[2] < 0 )
	    return( true );
	else if( i[2] >= (int32_t)_pidata._geom->size(2)-1 )
	    return( true );

	Vec3D v1 = x1.location();
	Vec3D v2 = x2.location();

	DEBUG_MESSAGE( "v1 = " << v1 << "\n" );
	DEBUG_MESSAGE( "v2 = " << v2 << "\n" );

	int32_t tric = _pidata._geom->surface_trianglec( i[0], i[1], i[2] );
	int32_t ptr = _pidata._geom->surface_triangle_ptr( i[0], i[1], i[2] );

	DEBUG_MESSAGE( "tric = " << tric << "\n" );

	// Go through surface triangles at mesh cube i
	for( int32_t a = 0; a < tric; a++ ) {
	    const VTriangle &tri = _pidata._geom->surface_triangle( ptr+a );
	    const Vec3D &va = _pidata._geom->surface_vertex( tri[0] );
	    const Vec3D &vb = _pidata._geom->surface_vertex( tri[1] );
	    const Vec3D &vc = _pidata._geom->surface_vertex( tri[2] );

	    DEBUG_MESSAGE( "a = " << a << "\n" );
	    DEBUG_MESSAGE( "tri[" << a << "][0] = " << va << "\n" );
	    DEBUG_MESSAGE( "tri[" << a << "][1] = " << vb << "\n" );
	    DEBUG_MESSAGE( "tri[" << a << "][2] = " << vc << "\n" );

	    // Solve for intersection between trajectory segment and surface triangle
	    // r1 + (r2-r1)*K[0] = ra + (rb-ra)*K[1] + (rc-ra)*K[2]
	    Mat3D m( v2[0]-v1[0], -vb[0]+va[0], -vc[0]+va[0],
		     v2[1]-v1[1], -vb[1]+va[1], -vc[1]+va[1],
		     v2[2]-v1[2], -vb[2]+va[2], -vc[2]+va[2] );
	    double mdet = m.determinant();
	    if( mdet == 0.0 )
		continue;
	    Mat3D minv = m.inverse( mdet );
	    Vec3D off( -v1[0]+va[0], -v1[1]+va[1], -v1[2]+va[2] );
	    Vec3D K = minv*off;
	    double K3 = K[1]+K[2];
	    DEBUG_MESSAGE( "K = " << K << "\n" );

	    // Check for intersection at valid ranges
	    // Allow COLLISION_EPS amount of overlap, double inclusion is 
	    // not an issue here, missing an intersection is a problem.
	    if( K[0] > -COLLISION_EPS && K[0] < 1.0+COLLISION_EPS && 
		K[1] > -COLLISION_EPS && K[1] < 1.0+COLLISION_EPS && 
		K[2] > -COLLISION_EPS && K[2] < 1.0+COLLISION_EPS && 
		K3 > -COLLISION_EPS && K3 < 1.0+COLLISION_EPS ) {

		DEBUG_MESSAGE( "Intersection found\n" );

		// Found intersection, set collision coordinates
		for( uint32_t b = 0; b < PP::size(); b++ )
		    status_x[b] = x1[b] + K[0]*(x2[b]-x1[b]);

		// Remove all points from trajectory after time status_x[0].
		// Does this ever happen???
		for( int32_t b = _traj.size()-1; b > 0; b-- ) {
		    if( _traj[b][0] > status_x[0] )
			_traj.pop_back();
		    else
			break;
		}

		// Update status and statistics
		particle.set_status( PARTICLE_COLL );
		uint32_t solid = get_solid( i[0], i[1], i[2] );
		_stat.add_bound_collision( solid, particle.IQ() );

		if( _tsur_cb )
		    (*_tsur_cb)( &particle, &status_x, ptr+a, K[1], K[2] );

		DEBUG_MESSAGE( "Solid collision detected\n" );
		DEBUG_DEC_INDENT();

		return( false );
	    }
	}

	DEBUG_MESSAGE( "No collisions\n" );
	DEBUG_DEC_INDENT();

	return( true );
    }

    /*! \brief Check for particle collision with solid surface with inside().
     *
     *  Particle propagates from x1 to x2, where x1 is in
     *  vacuum. Checks if x2 is inside solid and if it is, the
     *  collision point is bracketed between x1 and x2. Particle
     *  coordinates at status_x are set to collision coordinates.
     *  Returns false if particle collided.
     */
    bool check_collision_solid( Particle<PP> &particle, const PP &x1, const PP &x2, 
				PP &status_x ) {

	DEBUG_MESSAGE( "Checking collisions with solids - inside check\n" );
	DEBUG_INC_INDENT();

	// If inside solid, bracket for collision point
	Vec3D v2 = x2.location();
	int32_t bound = _pidata._geom->inside( v2 );
	if( bound < 7 ) {
	    DEBUG_MESSAGE( "No collisions\n" );
	    DEBUG_DEC_INDENT();
	    return( true );
	}
	Vec3D vc;
	Vec3D v1 = x1.location();
	double K = _pidata._geom->bracket_surface( bound, v2, v1, vc );

	// Calculate new PP
	for( size_t a = 0; a < PP::size(); a++ )
	    status_x[a] = x2[a] + K*(x1[a]-x2[a]);

	// Remove all points from trajectory after time status_x[0].
	for( size_t a = _traj.size()-1; a > 0; a-- ) {
	    if( _traj[a][0] > status_x[0] )
		_traj.pop_back();
	    else
		break;
	}

	// Save last trajectory point and update status
	//save_trajectory_point( status_x );
	particle.set_status( PARTICLE_COLL );

	// Update collision statistics for boundary
	_stat.add_bound_collision( bound, particle.IQ() );

	DEBUG_MESSAGE( "Solid collision detected\n" );
	DEBUG_DEC_INDENT();

	return( false );
    }

    /*! \brief Check for particle collision with solid
     *
     *  Two different algorithms.
     */
    bool check_collision( Particle<PP> &particle, const PP &x1, const PP &x2, 
			  PP &status_x, const int i[3] ) {

	if( _surface_collision )
	    return( check_collision_surface( particle, x1, x2, status_x, i ) );
	return( check_collision_solid( particle, x1, x2, status_x ) );
    }


    /*! \brief Mirror trajectory.
     *
     *  Trajectory is mirrored at \a _coldata[c] on axis \a at \a
     *  border, where -1 is the negative side and +1 is the positive
     *  side.. Already saved trajectory points are checked back to
     *  _xi.
     */
    void handle_mirror( size_t c, int i[3], size_t a, int border, PP &x2 ) {

	DEBUG_MESSAGE( "Mirror trajectory\n" );
	DEBUG_INC_INDENT();

	double xmirror;
	if( border < 0 ) {
	    xmirror = _pidata._geom->origo(a);
	    i[a] = -i[a]-1;
	} else {
	    xmirror = _pidata._geom->max(a);
	    i[a] = 2*_pidata._geom->size(a)-i[a]-3;
	}

	DEBUG_MESSAGE( "xmirror = " << xmirror << "\n" );
	DEBUG_MESSAGE( "i = (" << i[0] << ", " << i[1] << ", " << i[2] << ")\n" );
	DEBUG_MESSAGE( "xi = " << _xi << "\n" );
	
	// Check if found edge at first encounter
	bool caught_at_boundary = false;
	if( _coldata[c]._dir == border*((int)a+1) && 
	    ( i[a] == 0 || i[a] == (int)_pidata._geom->size(a)-2 ) ) {
	    caught_at_boundary = true;
	    DEBUG_MESSAGE( "caught_at_boundary\n" );
	}

	// Mirror traj back to _xi
	if( caught_at_boundary ) {
	    save_trajectory_point( _coldata[c]._x );
	} else {
	    for( int b = _traj.size()-1; b > 0; b-- ) {
		if( _traj[b][0] >= _xi[0] ) {
		    
		    DEBUG_MESSAGE( "mirroring traj[" << b << "] = " << _traj[b] << "\n" );
		    _traj[b][2*a+1] = 2.0*xmirror - _traj[b][2*a+1];
		    _traj[b][2*a+2] *= -1.0;
		} else
		    break;
	    }
	}

	// Mirror rest of the coldata
	for( size_t b = c; b < _coldata.size(); b++ ) {
	    if( (size_t)abs(_coldata[b]._dir) == a+1 )
		_coldata[b]._dir *= -1;
	    _coldata[b]._x[2*a+1] = 2.0*xmirror - _coldata[b]._x[2*a+1];
	    _coldata[b]._x[2*a+2] *= -1.0;
	}

	if( caught_at_boundary )
	    save_trajectory_point( _coldata[c]._x );

	// Mirror calculation point
	x2[2*a+1] = 2.0*xmirror - x2[2*a+1];
	x2[2*a+2] *= -1.0;
	
	// Coordinates changed, reset integrator
	gsl_odeiv_step_reset( _step );
	gsl_odeiv_evolve_reset( _evolve );

	DEBUG_DEC_INDENT();
    }


    void handle_collision( Particle<PP> &particle, uint32_t bound, size_t c, PP &status_x ) {

	DEBUG_MESSAGE( "Handle collision\n" );

	//save_trajectory_point( _coldata[c]._x );
	status_x = _coldata[c]._x;
	particle.set_status( PARTICLE_OUT );
	_stat.add_bound_collision( bound, particle.IQ() );
    }


    /*! \brief Return if node (i,j) is a solid node
     */
    bool is_solid( int i, int j ) {
	uint32_t node = _pidata._geom->mesh_check( i, j );
	if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	    (node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	    return( true );
	return( false );
    }

    /*! \brief Return if node (i,j,k) is a solid node
     */
    bool is_solid( int i, int j, int k ) {
	uint32_t node = _pidata._geom->mesh_check( i, j, k );
	if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	    (node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	    return( true );
	return( false );
    }

    /*! \brief Handle particle mesh intersection.
     *
     *  Particle mesh coordinates \a i are advanced through
     *  intersection described in \a _coldata[c]. Collisions with
     *  solids and particle mirroring on boundaries are handled. Final
     *  particle coordinates \a x2 are updated accordingly. Returns
     *  false if particle collided.
     */
    bool handle_trajectory_advance( Particle<PP> &particle, size_t c, int i[3], PP &x2 ) {

	bool surface_collision = false;
	DEBUG_MESSAGE( "Handle trajectory advance\n" );
	DEBUG_INC_INDENT();

	// Check for collisions with solids and advance coordinates i.
	if( PP::dim() == 2 ) {
	    if( _coldata[c]._dir == -1 ) {
		if( (is_solid(i[0],i[1]) || is_solid(i[0], i[1]+1)) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2, i ) )
		    surface_collision = true;
		i[0]--;
	    } else if( _coldata[c]._dir == +1 ) {
		if( (is_solid(i[0]+1,i[1]) || is_solid(i[0]+1,i[1]+1)) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2, i ) )
		    surface_collision = true;
		i[0]++;
	    } else if( _coldata[c]._dir == -2 ) {
		if( (is_solid(i[0],i[1]) || is_solid(i[0]+1,i[1])) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2, i ) )
		    surface_collision = true;
		i[1]--;
	    } else {
		if( (is_solid(i[0],  i[1]+1) || is_solid(i[0]+1,i[1]+1)) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2, i ) )
		    surface_collision = true;
		i[1]++;
	    }
	} else if( PP::dim() == 3 ) {
	    if( _coldata[c]._dir == -1 ) {
		if( (is_solid(i[0],  i[1],  i[2]  ) || 
		     is_solid(i[0],  i[1]+1,i[2]  ) ||
		     is_solid(i[0],  i[1],  i[2]+1) ||
		     is_solid(i[0],  i[1]+1,i[2]+1)) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2, i ) )
		    surface_collision = true;
		i[0]--;
	    } else if( _coldata[c]._dir == +1 ) {
		if( (is_solid(i[0]+1,i[1],  i[2]  ) || 
		     is_solid(i[0]+1,i[1]+1,i[2]  ) ||
		     is_solid(i[0]+1,i[1],  i[2]+1) ||
		     is_solid(i[0]+1,i[1]+1,i[2]+1)) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2, i ) )
		    surface_collision = true;
		i[0]++;
	    } else if( _coldata[c]._dir == -2 ) {
		if( (is_solid(i[0],  i[1],i[2]  ) || 
		     is_solid(i[0]+1,i[1],i[2]  ) ||
		     is_solid(i[0],  i[1],i[2]+1) ||
		     is_solid(i[0]+1,i[1],i[2]+1)) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2, i ) )
		    surface_collision = true;
		i[1]--;
	    } else if( _coldata[c]._dir == +2 ) {
		if( (is_solid(i[0],  i[1]+1,i[2]  ) || 
		     is_solid(i[0]+1,i[1]+1,i[2]  ) ||
		     is_solid(i[0],  i[1]+1,i[2]+1) ||
		     is_solid(i[0]+1,i[1]+1,i[2]+1)) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2, i ) )
		    surface_collision = true;
		i[1]++;
	    } else if( _coldata[c]._dir == -3 ) {
		if( (is_solid(i[0],  i[1],  i[2]) || 
		     is_solid(i[0]+1,i[1],  i[2]) ||
		     is_solid(i[0],  i[1]+1,i[2]) ||
		     is_solid(i[0]+1,i[1]+1,i[2])) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2, i ) )
		    surface_collision = true;
		i[2]--;
	    } else {
		if( (is_solid(i[0],  i[1],  i[2]+1) || 
		     is_solid(i[0]+1,i[1],  i[2]+1) ||
		     is_solid(i[0],  i[1]+1,i[2]+1) ||
		     is_solid(i[0]+1,i[1]+1,i[2]+1)) &&
		    !check_collision( particle, _xi, _coldata[c]._x, x2, i ) )
		    surface_collision = true;
		i[2]++;
	    }
	} else {
	    throw( Error( ERROR_LOCATION, "unsupported dimension number" ) );
	}
	
	if( surface_collision ) {
	    DEBUG_MESSAGE( "Surface collision!\n" );
	    DEBUG_DEC_INDENT();
	    return( false );
	}

	// Check for collisions/mirroring with simulation boundary. Here
	// coordinates i are already advanced to next mesh.
	for( size_t a = 0; a < PP::dim(); a++ ) {

	    if( i[a] < 0 ) {
		DEBUG_MESSAGE( "Boundary collision at boundary " << 2*a+0 << "\n" );
		if( _mirror[2*a] )
		    handle_mirror( c, i, a, -1, x2 );
		else {
		    handle_collision( particle, 1+2*a, c, x2 );
		    DEBUG_DEC_INDENT();
		    return( false );
		}
	    } else if( i[a] >= (int32_t)(_pidata._geom->size(a)-1) ) {
		DEBUG_MESSAGE( "Boundary collision at boundary " << 2*a+1 << "\n" );
		if( _mirror[2*a+1] )
		    handle_mirror( c, i, a, +1, x2 );
		else {
		    handle_collision( particle, 2+2*a, c, x2 );
		    DEBUG_DEC_INDENT();
		    return( false );
		}
	    }
	}

	DEBUG_MESSAGE( "No collision.\n" );
	DEBUG_DEC_INDENT();
	return( true );
    }

    /*! \brief Limit trajectory advance from \a x1 to \a x2 to double
     *  the simulation box.
     *
     *  Return true if limitation is done, false if data is left untouched.
     */
    bool limit_trajectory_advance( const PP &x1, PP &x2 ) {

	bool touched = false;
	for( size_t a = 0; a < PP::dim(); a++ ) {

	    double lim1 = _pidata._geom->origo(a) - 
		(_pidata._geom->size(a)-1)*_pidata._geom->h();
	    double lim2 = _pidata._geom->origo(a) + 
		2*(_pidata._geom->size(a)-1)*_pidata._geom->h();

	    if( x2[2*a+1] < lim1 ) {
		double K = (lim1 - x1[2*a+1]) / (x2[2*a+1] - x1[2*a+1]);
		x2 = x1 + K*(x2-x1);
		touched = true;
		DEBUG_MESSAGE( "Limiting step to " << x2 << "\n" );
	    } else if(x2[2*a+1] > lim2 ) {
		double K = (lim2 - x1[2*a+1]) / (x2[2*a+1] - x1[2*a+1]);
		x2 = x1 + K*(x2-x1);
		touched = true;
		DEBUG_MESSAGE( "Limiting step to " << x2 << "\n" );
	    }
	}

	return( touched );
    }

    /*! \brief Build coldata.
     *
     *  Throw error if run out of memory.
     */
    void build_coldata( bool force_linear, const PP &x1, const PP &x2 ) {

	try {
	    if( _intrp == TRAJECTORY_INTERPOLATION_POLYNOMIAL && !force_linear )
		ColData<PP>::build_coldata_poly( _coldata, *_pidata._geom, x1, x2 );
	    else
		ColData<PP>::build_coldata_linear( _coldata, *_pidata._geom, x1, x2 );
	} catch( std::bad_alloc ) {
	    throw( ErrorNoMem( ERROR_LOCATION, "out of memory building collision data" ) );
	}
    }


    /*! \brief Handle particle iteration step from coordinates \a x1 to \a x2.
     *
     *  Searches mesh intersections between points \a x1 and \a x2 and
     *  builds ColData. Checks for collisions with solids and
     *  boundaries and sets space charge. If particle collides with
     *  mirroring boundary, the \a x2 coordinates are changed and GSL
     *  ODE integrator is resetted. Space charge is deposited at each
     *  mesh line crossing.
     *
     *  If \a force_linear is true, linear interpolation of trajectory
     *  is used regardless of interpolation settings.
     *
     *  If \a first_step is true, this is the first step and the
     *  particle is allowed to get into the simulation box if the
     *  definition point (x1) was outside the geometry.
     *
     *  Return true if particle status is PARTICLE_OK after trajectory
     *  step, false otherwise.
     */
    bool handle_trajectory( Particle<PP> &particle, const PP &x1, PP &x2, 
			    bool force_linear, bool first_step ) {

	DEBUG_MESSAGE( "Handle trajectory from x1 to x2:\n" <<
		       "  x1 = " << x1 << "\n" << 
		       "  x2 = " << x2 << "\n" );
	DEBUG_INC_INDENT();

	// Limit trajectory advance to double the simulation box
	// If limitation done, force to linear interpolation
	if( limit_trajectory_advance( x1, x2 ) )
	    force_linear = true;

	build_coldata( force_linear, x1, x2 );

	// TODO
	// Remove entrance to geometry if coming from outside or make
	// code to skip the collision detection for these particles

	// No intersections, nothing to do
  //NdH: I added also here a TrajectoryHandlerCallback, and a check to exit this function with false if the particle is neutralized. 
	if( _coldata.size() == 0 ) 
  	{
	    DEBUG_MESSAGE( "No coldata\n" );
	    DEBUG_DEC_INDENT();
		if( _thand_cb )
		{
			if( _pidata._geom->inside( x2.location() ) == 0 ) //NdH: Only call TrajectoryHandlerCallback if the end point is a vacuum node. Not sure why this is not checked by IBSimu itself.
			{
				PP x1temp = PP(); 
				x1temp.clear(); 
				x1temp = x1temp + x1; 
				(*_thand_cb)( &particle, &x1temp, &x2 ); // CP 20240521 back to previous definition of the trajectory handler operator
          	}
      
			if( particle.get_status() != PARTICLE_OK ) // Particle can be charge-exchanged: status = (2, PARTICLE_COLL). So build in detection here.
			{
				// ibsimu.message( 1 ) << "Coldata length zero, and CX event: special case! Size of pdb: " << to_string( _pdb->size() ) << "\n";
				save_trajectory_point( x2 );
				return( false );
          	}
      	}
	    return( true );
	}

	// Starting mesh index
	int i[3] = {0, 0, 0};
	for( size_t cdir = 0; cdir < PP::dim(); cdir++ )
	    i[cdir] = (int)floor( (x1[2*cdir+1]-_pidata._geom->origo(cdir))*_pidata._geom->div_h() );

	// Process intersection points
	DEBUG_MESSAGE( "Process coldata points:\n" );
	DEBUG_INC_INDENT();
	for( size_t a = 0; a < _coldata.size(); a++ ) {

	    if( _save_points )
		save_trajectory_point( _coldata[a]._x );

	    DEBUG_MESSAGE( "Coldata " << a << "\n" <<
			   "  x = " << _coldata[a]._x << "\n" <<
			   "  i = (" << std::setw(3) << i[0] << ", "
			   << std::setw(3) << i[1] << ", "
			   << std::setw(3) << i[2] << ")\n" << 
			   "  dir = " << _coldata[a]._dir << "\n" );

	    // Update space charge for mesh volume i
	    if( _scharge_dep == SCHARGE_DEPOSITION_LINEAR && _pidata._scharge ) {
		_cdpast.push( _coldata[a]._x );
		scharge_add_from_trajectory_linear( *_pidata._scharge, _scharge_mutex, particle.IQ(),
						    _coldata[a]._dir, _cdpast, i );
	    }

	    // Advance particle in mesh, update i, check for possible
	    // collisions and do mirroring for trajectory and coldata.
	    handle_trajectory_advance( particle, a, i, x2 );

	    if( _scharge_dep == SCHARGE_DEPOSITION_PIC && _pidata._scharge ) {
		scharge_add_from_trajectory_pic( *_pidata._scharge, _scharge_mutex, particle.IQ(), 
						 _xi, _coldata[a]._x );
	    }

      	if( (_thand_cb) && (particle.get_status() == PARTICLE_OK) )  //As a result of this being 'in the loop', the TrajectoryHandlerCallback is called every intersection with the grid. Maybe to speed up, this can be moved outwards later.
      	{                //But then the step size need to be limited, and/or the trajectory integrated. It also might cause problems with the space charge deposition.
          if( a == 0 ) //Since for _coldata.size() of n, there are n+1 trajectory pieces, the first piece gets treated differently.
          {
              PP x1temp = PP(); 
              x1temp.clear();   //Not sure if this is needed.
              x1temp = x1temp + x1; 

              (*_thand_cb)( &particle, &x1temp, &_coldata[a]._x );  // CP 20240521 back to previous definition of the trajectory handler operator
          }
          else 
          {
              (*_thand_cb)( &particle, &_coldata[a - 1]._x, &_coldata[a]._x ); // CP 20240521 back to previous definition of the trajectory handler operator
          }

          if( (a == (_coldata.size() - 1 )) && ( _pidata._geom->inside( x2.location() ) == 0 ) && (particle.get_status() == PARTICLE_OK) ) //Note the if instead of the else if -> there could be only 1 coldata point in between the first and last point of the trajectory!
          {   //NdH: Only call TrajectoryHandlerCallback if the end point is a vacuum node. Not sure why this is not checked by IBSimu itself.
              (*_thand_cb)( &particle, &_coldata[a]._x, &x2 );  // CP 20240521 back to previous definition of the trajectory handler operator
          }

      	} // (*_thand_cb)( &particle, &_coldata[a]._x, &x2 ); //This was the original function call.

	    // Clear coldata and exit if particle collided.
	    if( particle.get_status() != PARTICLE_OK ) {
			// save_trajectory_point( x2 );
			save_trajectory_point( _coldata[a]._x ); // CP 20240522 save the location at which the particle has changed status
			_coldata.clear();
			DEBUG_DEC_INDENT();
			DEBUG_MESSAGE( "Interrupting coldata processing\n" );
			DEBUG_DEC_INDENT();
			return( false );
	    }

	    // Update last accepted intersection point xi.
	    _xi = _coldata[a]._x;
	}

	_coldata.clear();
	DEBUG_DEC_INDENT();
	DEBUG_MESSAGE( "Coldata done\n" );
	DEBUG_DEC_INDENT();
	return( true );
    }


    /*! \brief Is particle mirroring required at axis in cylindrical symmetry?
     */
     bool axis_mirror_required( const PP &x2 ) {
	 return( _pidata._geom->geom_mode() == MODE_CYL && 
		 x2[4] < 0.0 && 
		 x2[3] <= 0.01*_pidata._geom->h() &&
		 x2[3]*fabs(x2[5]) <= 1.0e-9*fabs(x2[4]) );
		 
     }


    /*! \brief Handle particle mirrored at axis in cylindrical symmetry.
     *
     *  Return true if particle status is PARTICLE_OK after trajectory
     *  step, false otherwise.
     */
    bool handle_axis_mirror_step( Particle<PP> &particle, const PP &x1, PP &x2 ) {

	// Get acceleration at x2
	double dxdt[5];
	PP::get_derivatives( x2[0], &x2[1], dxdt, (void *)&_pidata );

	// Calculate crossover point assuming zero acceleration in
	// r-direction and constant acceleration in x-direction
	double dt = -x2[3]/x2[4];
	PP xc;
	xc.clear();
	xc[0] = x2[0]+dt;
	xc[1] = x2[1]+(x2[2]+0.5*dxdt[1]*dt)*dt;
	xc[2] = x2[2];
	xc[3] = x2[3]+x2[4]*dt;
	xc[4] = x2[4];
	xc[5] = x2[5];

	// Mirror x2 to x3
	PP x3 = 2*xc - x2;
	x3[3] *= -1.0;
	x3[4] *= -1.0;
	x3[5] *= -1.0;

	DEBUG_MESSAGE( "Particle mirror:\n" <<
		       "  x1: " << x1 << "\n" <<
		       "  x2: " << x2 << "\n" << 
		       "  xc: " << xc << "\n" << 
		       "  x3: " << x3 << "\n" );
	
	// Handle step with linear interpolation to avoid going to r<=0
	if( !handle_trajectory( particle, x2, x3, true, false ) )
	    return( false ); // Particle done

	// Save trajectory calculation points
	save_trajectory_point( x2 );
	save_trajectory_point( xc );
	xc[4] *= -1.0;
	xc[5] *= -1.0;
	save_trajectory_point( xc );
	
	// Next step not a continuation of previous one, reset
	// integrator
	gsl_odeiv_step_reset( _step );
	gsl_odeiv_evolve_reset( _evolve );

	// Continue iteration at mirrored point
	x2 = x3;
	return( true );
    }
    
    /*! \brief Check particle definition.
     *
     *  False is returned if particle is defined out of simulation
     *  area or inside a solid. False is also returned if illegal
     *  input (NaN) is received.
     *
     *  Removed feature: False is also returned if a particle
     *  is defined on the simulation border going out of simulation
     *  with the corresponding border not having mirroring enabled. If
     *  the mirroring is enabled, the particle is mirrored and true is
     *  returned. Otherwise the particle definition is ok and true is
     *  returned.
     */
    bool check_particle_definition( PP &x ) {

	DEBUG_MESSAGE( "Particle defined at x = " << x << "\n" );
	DEBUG_INC_INDENT();

	// Check for NaN
	for( size_t a = 0; a < PP::size(); a++ ) {
	    if( comp_isnan( x[a] ) ) {
		DEBUG_MESSAGE( "Particle coordinate NaN\n" );
		DEBUG_DEC_INDENT();
    // ibsimu.message( 1 ) << "NdH: BADDEF detected, Particle coordinate NaN\n";
		return( false );
	    }
	}

	// Check if inside solids or outside simulation geometry box.
	if( _surface_collision ) {
	    if( _pidata._geom->surface_inside( x.location() ) ) {
		DEBUG_MESSAGE( "Particle inside solid or outside simulation\n" );
		DEBUG_DEC_INDENT();
    // ibsimu.message( 1 ) << "NdH: BADDEF detected, Particle inside solid or outside simulation " << to_string(_pidata._geom->surface_inside( x.location() )) << "\n";
		return( false );
	    }
	} else { 
	    if( _pidata._geom->inside( x.location() ) ) {
		DEBUG_MESSAGE( "Particle inside solid or outside simulation\n" );
		DEBUG_DEC_INDENT();
    // ibsimu.message( 1 ) << "NdH: BADDEF detected, Particle inside solid or outside simulation " << to_string( _pidata._geom->inside( x.location() ) )<< "\n";
		return( false );
	    }
	}

	// Check if particle on simulation geometry border and directed outwards
	/*
	for( size_t a = 0; a < PP::dim(); a++ ) {
	    if( x[2*a+1] == _pidata._geom->origo(a) && x[2*a+2] < 0.0 ) {
		if( _mirror[2*a] ) {
		    x[2*a+2] *= -1.0;
#ifdef DEBUG_PARTICLE_ITERATOR
	ibsimu.message(MSG_DEBUG_GENERAL,1) << "Mirroring to:\n";
	ibsimu.message(MSG_DEBUG_GENERAL,1) << "  x = " << x << "\n";
#endif
		} else {
		    return( false );
		}

	    } else if( x[2*a+1] == _pidata._geom->max(a) & x[2*a+2] > 0.0 ) {
		if( _mirror[2*a+1] ) {
		    x[2*a+2] *= -1.0;
#ifdef DEBUG_PARTICLE_ITERATOR
	ibsimu.message(MSG_DEBUG_GENERAL,1) << "Mirroring to:\n";
	ibsimu.message(MSG_DEBUG_GENERAL,1) << "  x = " << x << "\n";
#endif
		} else {
		    return( false );
		}
	    }
	}

#ifdef DEBUG_PARTICLE_ITERATOR
	ibsimu.message(MSG_DEBUG_GENERAL,1) << "Definition ok\n";
#endif

	*/
	DEBUG_MESSAGE( "Ok\n" );
	DEBUG_DEC_INDENT();
	return( true );
    }
    
    double calculate_dt( const PP &x, const double *dxdt ) {

	double spd = 0.0, acc = 0.0;

	for( size_t a = 0; a < (PP::size()-1)/2; a++ ) {
	    //ibsimu.message(MSG_DEBUG_GENERAL,1) << "spd += " << dxdt[2*a]*dxdt[2*a] << "\n";
	    spd += dxdt[2*a]*dxdt[2*a];
	    //ibsimu.message(MSG_DEBUG_GENERAL,1) << "acc += " << dxdt[2*a+1]*dxdt[2*a+1] << "\n";
	    acc += dxdt[2*a+1]*dxdt[2*a+1];
	}
	if( _pidata._geom->geom_mode() == MODE_CYL ) {
	    //ibsimu.message(MSG_DEBUG_GENERAL,1) << "MODE_CYL\n";
	    //ibsimu.message(MSG_DEBUG_GENERAL,1) << "spd += " << x[3]*x[3]*x[5]*x[5] << "\n";
	    spd += x[3]*x[3]*x[5]*x[5];
	    //ibsimu.message(MSG_DEBUG_GENERAL,1) << "acc += " << x[3]*x[3]*dxdt[4]*dxdt[4] << "\n";
	    acc += x[3]*x[3]*dxdt[4]*dxdt[4];
	}
	double v2 = spd;
	//ibsimu.message(MSG_DEBUG_GENERAL,1) << "spd = " << sqrt(spd) << "\n";
	//ibsimu.message(MSG_DEBUG_GENERAL,1) << "acc = " << sqrt(acc) << "\n";
	spd = _pidata._geom->h() / sqrt(spd);
	acc = sqrt( 2.0*_pidata._geom->h() / sqrt(acc) );
	double dt = ( spd < acc ? spd : acc );

	if( _pidata._relativistic && _pidata._bfield ) {
	    Vec3D xc = x.location();
	    Vec3D B = (*_pidata._bfield)( xc );
	    if( _pidata._bsup_cb ) {
		double factor = (*_pidata._bsup_cb)( xc );
		B *= factor;
	    }
	    double b2 = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
	    double qm = fabs( _pidata._qm );
	    if( comp_isfinite(v2) && v2 >= 0.0 && v2 < SPEED_C2 &&
		comp_isfinite(b2) && b2 > 0.0 && qm > 0.0 ) {
		// CPI patch on 20260604: resolve the first relativistic trial step below
		// one twentieth of the local gyroperiod in strong B fields.
		double gamma = 1.0/sqrt( 1.0 - v2/SPEED_C2 );
		double gyro_dt = 2.0*M_PI*gamma / ( 20.0*qm*sqrt(b2) );
		if( comp_isfinite(gyro_dt) && gyro_dt > 0.0 )
		    dt = ( dt < gyro_dt ? dt : gyro_dt );
	    }
	}

	return( dt );
    }

public:

    /*! \brief Constructor for new particle iterator.
     *
     *  New particle iterator is initialized with given settings.
     *
     *  \param type Particle iterator type used
     *  \param epsabs Absolute error limit in iteration
     *  \param epsrel Relative error limit in iteration
     *  \param intrp Interpolation type.
     *  \param scharge_dep Space charge deposition type.
     *  \param maxsteps Maximum number of steps to take before particle is killed
     *  \param maxt Maximum flight time for a particle
     *  \param save_points Flag for saving all intersection points of trajectories
     *  \param trajdiv Trajectory saving divisor. Only every trajdiv:th particle 
     *  trajectory saved.
     *  \param mirror %Particle mirroring on surfaces
     *  \param scharge Space charge field to save to
     *  \param scharge_mutex Space charge write mutex
     *  \param efield Electric field in the geometry
     *  \param bfield Magnetic field in the geometry
     *  \param geom %Geometry definition
     *
     *  The particle iterator is given the settings for calculation
     *  and geometry, electric field and space charge map to
     *  build. Pointer to first particle in the particle database
     *  vector is used to calculate the particle number from the
     *  particle memory location.
     */
    ParticleIterator( particle_iterator_type_e type, double epsabs, double epsrel, 
		      trajectory_interpolation_e intrp, scharge_deposition_e scharge_dep, 
		      uint32_t maxsteps, double maxt, bool save_points,
		      uint32_t trajdiv, bool mirror[6], MeshScalarField *scharge, 
		      pthread_mutex_t *scharge_mutex,
		      const VectorField *efield, const VectorField *bfield, 
		      const Geometry *geom ) 
	: _type(type), _intrp(intrp), _scharge_dep(scharge_dep), _epsabs(epsabs), _epsrel(epsrel), 
	  _maxsteps(maxsteps), _maxt(maxt), _save_points(save_points), _trajdiv(trajdiv), 
	  _surface_collision(false), _pidata(scharge,efield,bfield,geom), 
	  _thand_cb(0), _tend_cb(0), _tsur_cb(0), _bsup_cb(0), _esup_cb(0), _pdb(0), _scharge_mutex(scharge_mutex), 
	  _stat(geom->number_of_boundaries()) {
	
	// Initialize mirroring
	_mirror[0] = mirror[0];
	_mirror[1] = mirror[1];
	_mirror[2] = mirror[2];
	_mirror[3] = mirror[3];
	_mirror[4] = mirror[4];
	_mirror[5] = mirror[5];

	// Initialize system of ordinary differential equations (ODE)
	_system.jacobian  = NULL;
	_system.params    = (void *)&_pidata;
	_system.function  = PP::get_derivatives;
	_system.dimension = PP::size()-1; // Time is not part of differential equation dimensions

	// Make scale
	// 2D:  x vx y vy
	// Cyl: x vx r vr omega
	// 3D:  x vx y vy z vz
	double scale_abs[PP::size()-1];
	for( uint32_t a = 0; a < (uint32_t)PP::size()-2; a+=2 ) {
	    scale_abs[a+0] = 1.0;
	    scale_abs[a+1] = 1.0e6;
	}
	if( _pidata._geom->geom_mode() == MODE_CYL )
	    scale_abs[4] = 1.0;

	// Initialize ODE solver
	_step    = gsl_odeiv_step_alloc( gsl_odeiv_step_rkck, _system.dimension );
	//_control = gsl_odeiv_control_standard_new( _epsabs, _epsrel, 1.0, 1.0 );
	_control = gsl_odeiv_control_scaled_new( _epsabs, _epsrel, 1.0, 1.0, scale_abs, PP::size()-1 );
	_evolve  = gsl_odeiv_evolve_alloc( _system.dimension );
    }


    /*! \brief Destructor.
     */
    ~ParticleIterator() {
	gsl_odeiv_evolve_free( _evolve );
	gsl_odeiv_control_free( _control );
	gsl_odeiv_step_free( _step );
    }


    /*! \brief Enable/disable surface collision model.
     */
    void set_surface_collision( bool surface_collision ) {
	if( surface_collision && _pidata._geom->geom_mode() == MODE_2D )
	    throw( Error( ERROR_LOCATION, "2D surface collision not supported" ) );
	if( surface_collision && !_pidata._geom->surface_built() )
	    throw( Error( ERROR_LOCATION, "surface model not built" ) );
	_surface_collision = surface_collision;
    }


    /*! \brief Set trajectory handler callback. 
     */
    void set_trajectory_handler_callback( TrajectoryHandlerCallback *thand_cb, ParticleDataBase *pdb ) {
	_thand_cb = thand_cb;
  _pdb = pdb;
    }


    /*! \brief Set trajectory end callback. 
     */
    void set_trajectory_end_callback( TrajectoryEndCallback *tend_cb, ParticleDataBase *pdb ) {
	_tend_cb = tend_cb;
	_pdb = pdb;
    }


    /*! \brief Set trajectory surface collision callback. 
     */
    void set_trajectory_surface_collision_callback( TrajectorySurfaceCollisionCallback *tsur_cb ) {
	_tsur_cb = tsur_cb;
    }


    /*! \brief Set B-field potential dependent suppression callback.
     */
    void set_bfield_suppression_callback( const CallbackFunctorD_V *bsup_cb ) {
	_pidata.set_bfield_suppression_callback( bsup_cb );
    }

    /*! \brief Set E-field potential dependent suppression callback.
     */
    void set_efield_suppression_callback( const CallbackFunctorD_V *esup_cb ) {
    _pidata.set_efield_suppression_callback( esup_cb );
    }

    /*! \brief Set relativistic particle iteration.
     */
    void set_relativistic( bool enable ) {
	_pidata.set_relativistic( enable );
    }

    /*! \brief Get particle iterator statistics.
     */
    const ParticleStatistics &get_statistics( void ) const {
	return( _stat );
    }

    
    /*! \brief Iterate a particle from start to end.
     *
     *  Iterate particle \a particle from start to end. This function
     *  is called by the Scheduler \a scheduler, which provides
     *  particles to be solved. Reference to \a scheduler is provided
     *  for the possibility to add secondary particles to particle
     *  database.
     */
    void operator()( Particle<PP> *particle, uint32_t pi ) {

	DEBUG_MESSAGE( "Calculating particle " << pi << "\n" );
	DEBUG_INC_INDENT();

	// Check particle status
	if( particle->get_status() != PARTICLE_OK ) {
	    DEBUG_DEC_INDENT();
	    return;
	}

	// Copy starting point to x and 
	PP x = particle->x();

	// Check particle definition
	if( !check_particle_definition( x ) ) {
      // ibsimu.message( 1 ) << "NdH: BADDEF encountered in particle " << to_string(pi) << "\n";
	    particle->set_status( PARTICLE_BADDEF );
      // particle->set_status( PARTICLE_TIME ); //NdH: THIS IS ONLY DONE FOR DEBUGGING PURPOSES, SO SEE THE PARTICLE TRACK! IT SHOULD BE REVERTED IMMEDIATELY AFTERWARDS!!1!
	    _stat.inc_end_baddef();
	    DEBUG_DEC_INDENT();
	    return;
	}
	particle->x() = x;

	// Reset trajectory and save first trajectory point.
	_traj.clear();
	save_trajectory_point( x );
	_pidata._qm = particle->qm();
	_xi = x;
	if( _scharge_dep == SCHARGE_DEPOSITION_LINEAR ) {
	    _cdpast.clear();
	    _cdpast.push( x );
	}

	// Reset integrator
	gsl_odeiv_step_reset( _step );
	gsl_odeiv_evolve_reset( _evolve );
	
	// Make initial guess for step size
	double dxdt[PP::size()-1];
	PP::get_derivatives( 0.0, &x[1], dxdt, (void *)&_pidata );
	double dt = calculate_dt( x, dxdt );

	DEBUG_MESSAGE( "Starting values for iteration:\n" << 
		       "  dxdt = " );
	for( size_t a = 0; a < PP::size()-1; a++ )
	    DEBUG_MESSAGE( dxdt[a] << " " );
	DEBUG_MESSAGE( "\n" <<
		       "  dt = " << dt << "\n" );
	DEBUG_MESSAGE( "Starting iteration\n" );
	
	// Iterate ODEs until maximum steps are done, time is used 
	// or particle collides.
	PP x2;
	size_t nstp = 0; // Steps taken
	while( nstp < _maxsteps && x[0] < _maxt ) {

	    // Take a step.
	    x2 = x;

	    DEBUG_MESSAGE( "\nStep *** *** *** ****\n" <<
			   "  x  = " << x2 << "\n" <<
			   "  dt = " << dt << " (proposed)\n" );

	    while( true ) {
		int retval = gsl_odeiv_evolve_apply( _evolve, _control, _step, &_system, 
						     &x2[0], _maxt, &dt, &x2[1] );
		if( retval == IBSIMU_DERIV_ERROR ) {
		    DEBUG_MESSAGE( "Step rejected\n" <<
				   "  x2 = " << x2 << "\n" <<
				   "  dt = " << dt << "\n" );
		    x2[0] = x[0]; // Reset time (this shouldn't be necessary - there 
		                  // is a bug in GSL-1.12, report has been sent)
		    dt *= 0.5;
		    if( dt == 0.0 )
			throw( Error( ERROR_LOCATION, "too small step size" ) );
		    //nstp++;
		    continue;
		} else if( retval == GSL_SUCCESS ) {
		    break;
		} else {
		    throw( Error( ERROR_LOCATION, "gsl_odeiv_evolve_apply failed" ) );
		}
	    }
	    
	    // Check step count number and step size validity
	    if( nstp >= _maxsteps )
	        break;
	    if( x2[0] == x[0] ) {
	        // Print failed trajectory
		ibsimu.message( 1 ) << "Particle calculation failed. Coordinates:\n";
		for( size_t a = 0; a < _traj.size(); a++ )
		    ibsimu.message( 1 ) << a << " " << _traj[a] << "\n";
	        throw( Error( ERROR_LOCATION, "too small step size calculating particle " + to_string(pi) ) );
	    }
	    
	    // Increase step count.
	    nstp++;

	    DEBUG_MESSAGE( "Step accepted from x1 to x2:\n" <<
			   "  x1 = " << x << "\n" <<
			   "  x2 = " << x2 << "\n" );

	    // Handle collisions and space charge of step.
	    if( !handle_trajectory( *particle, x, x2, false, nstp == 0 ) ) {
		x = x2;
		break; // Particle done
	    }

	    // Check if particle mirroring is required to avoid 
	    // singularity at symmetry axis.
	    if( axis_mirror_required( x2 ) ) {
		if( !handle_axis_mirror_step( *particle, x, x2 ) )
		    break; // Particle done
	    }

	    // Propagate coordinates
	    x = x2;

	    // Save trajectory point
	    save_trajectory_point( x2 );
	}

	DEBUG_MESSAGE( "Done iterating\n" );

	// Check if step count or time limited 
	if( nstp == _maxsteps ) {
	    particle->set_status( PARTICLE_NSTP );
	    _stat.inc_end_step();
	} else if( x[0] >= _maxt ) {
	    particle->set_status( PARTICLE_TIME );
	    _stat.inc_end_time();
	}

	// Save step count
	_stat.inc_sum_steps( nstp );

	// Save trajectory of current particle
	if( _trajdiv != 0 && pi % _trajdiv == 0 )
	    particle->copy_trajectory( _traj );

	// Save last particle location
	particle->x() = x;

	// Call trajectory end callback
	if( _tend_cb )
	    (*_tend_cb)( particle, _pdb );

	DEBUG_DEC_INDENT();
    }

};


#endif

