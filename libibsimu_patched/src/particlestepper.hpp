/*! \file particlestepper.hpp
 *  \brief %Particle stepper using Boris leap-frog
 */

/* Copyright (c) 2015,2022 Taneli Kalvas. All rights reserved.
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


#ifndef PARTICLESTEPPER_HPP
#define PARTICLESTEPPER_HPP 1


#include "geometry.hpp"
#include "particles.hpp"
#include "vectorfield.hpp"
#include "meshscalarfield.hpp"
#include "scharge.hpp"
#include "particledatabase.hpp"


template<class PP> class ParticleStepper {

    double                     _dt;
    bool                       _surface_collision;
    uint32_t                   _trajdiv;       /*!< \brief Divisor for saved trajectories,
					        * if 3, every third particle is saved. */
    bool                       _mirror[6];     /*!< \brief Is particle mirrored on boundary? */
    MeshScalarField           *_scharge;
    const VectorField         *_efield;
    const VectorField         *_bfield;
    const Geometry            *_geom;

public:


    ParticleStepper( double dt, uint32_t trajdiv, bool mirror[6], 
		     MeshScalarField *scharge, const VectorField *efield, 
		     const VectorField *bfield, const Geometry *geom )
	: _dt(dt), _surface_collision(false), _trajdiv(trajdiv), _scharge(scharge), 
	  _efield(efield), _bfield(bfield), _geom(geom) {
	// Initialize mirroring
	_mirror[0] = mirror[0];
	_mirror[1] = mirror[1];
	_mirror[2] = mirror[2];
	_mirror[3] = mirror[3];
	_mirror[4] = mirror[4];
	_mirror[5] = mirror[5];
    }


    ~ParticleStepper() {

    }


    /*! \brief Initialize particle stepping velocity backwards by 0.5*dt
     */
    void initialize( Particle<PP> *particle, uint32_t pi ) {
	if( PP::geom_mode() == MODE_3D ) {
	    Vec3D E, B;
	    Vec3D x = particle->x().location();
	    Vec3D v = particle->x().velocity();
	    if( _efield )
		E = (*_efield)( x );
	    if( _bfield )
		B = (*_bfield)( x );
	    Vec3D a = particle->qm()*(E+cross(v,B));
	    v = v - 0.5*a*_dt;
	    (*particle)[2] = v[0];
	    (*particle)[4] = v[1];
	    (*particle)[6] = v[2];

	} else {
	    throw( ErrorUnimplemented( ERROR_LOCATION, "Particle stepping for geometry modes other than MODE_3D unimplemented" ) );
	}
    }


    /*! \brief Take one dt step forward for particle
     */
    void step( Particle<PP> *particle, uint32_t pi ) {

	// Check particle status
	if( particle->get_status() != PARTICLE_OK )
	    return;

	Vec3D x;
	if( PP::geom_mode() == MODE_3D ) {
	    Vec3D E, B;
	    x = particle->x().location();
	    if( _efield )
		E = (*_efield)( x );
	    if( _bfield )
		B = (*_bfield)( x );
	    Vec3D vminus = particle->x().velocity() + 0.5*particle->qm()*E*_dt;
	    Vec3D t = 0.5*particle->qm()*B*_dt;
	    Vec3D vprime = vminus + cross(vminus,t);
	    Vec3D s = 2.0/(1+t.ssqr())*t;
	    Vec3D vplus = vminus + cross(vprime,s);
	    Vec3D v = vplus + 0.5*particle->qm()*E*_dt;
	    x += v*_dt;

	    for( int a = 0; a < 3; a++ ) {
		if( _mirror[2*a] && x[a] < _geom->origo(a) ) {
		    x[a] = 2.0*_geom->origo(a) - x[a];
		    v[a] *= -1;
		}
		if( _mirror[2*a+1] && x[a] > _geom->max(a) ) {
		    x[a] = 2.0*_geom->max(a) - x[a];
		    v[a] *= -1;
		}
	    }

	    (*particle)[0] += _dt; // Advance time, time follows position, not velocity
	    (*particle)[1] = x[0];
	    (*particle)[2] = v[0];
	    (*particle)[3] = x[1];
	    (*particle)[4] = v[1];
	    (*particle)[5] = x[2];
	    (*particle)[6] = v[2];

	} else {
	    throw( ErrorUnimplemented( ERROR_LOCATION, "Particle stepping for geometry modes other than MODE_3D unimplemented" ) );
	}

	// Collision detection
	if( _surface_collision ) {
	    if( _geom->surface_inside( x ) ) {
		particle->set_status( PARTICLE_COLL );
		return;
	    }
	} else { 
	    if( _geom->inside( x ) ) {
		particle->set_status( PARTICLE_COLL );
		return;
	    }
	}

	// Space charge deposition
	scharge_add_step_pic( *_scharge, particle->IQ(), x );

	// Save trajectory data
	if( _trajdiv != 0 && pi % _trajdiv == 0 )
	    particle->add_trajectory_point( particle->x() );
    }



    /*! \brief Enable/disable surface collision model.
     */
    void set_surface_collision( bool surface_collision ) {
	if( surface_collision && _geom->geom_mode() == MODE_2D )
	    throw( Error( ERROR_LOCATION, "2D surface collision not supported" ) );
	if( surface_collision && !_geom->surface_built() )
	    throw( Error( ERROR_LOCATION, "surface model not built" ) );
	_surface_collision = surface_collision;
    }
};


#endif

