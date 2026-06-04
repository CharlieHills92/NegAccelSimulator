/*! \file particles.cpp
 *  \brief %Particle and particle point objects
 */

/* Copyright (c) 2005-2014 Taneli Kalvas. All rights reserved.
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

#include "particles.hpp"
#include "constants.hpp"
#include "trajectory.hpp"
#include "compmath.hpp"
#include "mat3d.hpp"
#include "ibsimu.hpp"
#include <iostream>
#include <iomanip>


//#define DEBUG_PARTICLE_DERIVATIVES 1


#ifdef DEBUG_PARTICLE_DERIVATIVES
#define DEBUG_MESSAGE(x) ibsimu.message(MSG_DEBUG_GENERAL,1) << x
#define DEBUG_INC_INDENT() ibsimu.inc_indent()
#define DEBUG_DEC_INDENT() ibsimu.dec_indent()
#else
#define DEBUG_MESSAGE(x) do {} while(0)
#define DEBUG_INC_INDENT() do {} while(0)
#define DEBUG_DEC_INDENT() do {} while(0)
#endif


int ParticleP2D::get_derivatives( double t, const double *x, double *dxdt, void *data )
{
    DEBUG_MESSAGE( "get_derivatives: " + to_string(x[0]) + " " + to_string(x[2]) + "\n" );

    Vec3D E, B, xc( x[0], x[2], 0.0 );
    ParticleIteratorData *pidata = (ParticleIteratorData *)data;

    // Prevent sampling field data outside simulation box, allow h tolerance
    for( uint32_t a = 0; a < 2; a++ ) {
	if( xc[a] < pidata->_geom->origo(a)-pidata->_geom->h() )
	    return( IBSIMU_DERIV_ERROR );
	else if( xc[a] > pidata->_geom->max(a)+pidata->_geom->h() )
	    return( IBSIMU_DERIV_ERROR );
    }
    
    if( pidata->_efield )
	E = (*pidata->_efield)( xc );
    if( pidata->_bfield )
	B = (*pidata->_bfield)( xc );

    if( comp_isnan(E[0]) || comp_isnan(E[1]) ) {
	DEBUG_MESSAGE( "return NaN\n" );
	return( IBSIMU_DERIV_ERROR );
    }

    /* Positions: dx/dt = vx, dy/dt = vy */
    dxdt[0] = x[1];
    dxdt[2] = x[3];

    /* Velocities dvx/dt = ax, dvy/dt = ay */
    if( pidata->_relativistic ) {
	double pxm = pidata->_qm * (E[0] + x[3]*B[2]);
	double pym = pidata->_qm * (E[1] - x[1]*B[2]);
	double v2 = x[1]*x[1] + x[3]*x[3];
	double gamma = 1.0/sqrt( 1.0 - v2/SPEED_C2 );
	double gamma3 = gamma*gamma*gamma;
	double a = gamma3*x[1]*x[1]/SPEED_C2 + gamma;
	double b = gamma3*x[1]*x[3]/SPEED_C2; // = c
	double d = gamma3*x[3]*x[3]/SPEED_C2 + gamma;
	double idet = 1.0/( a*d - b*b );
	dxdt[1] = idet*( d*pxm - b*pym );
	dxdt[3] = idet*( -b*pxm + a*pym );
    } else {
	dxdt[1] = pidata->_qm * (E[0] + x[3]*B[2]);
	dxdt[3] = pidata->_qm * (E[1] - x[1]*B[2]);
    }
	
    DEBUG_MESSAGE( "return dxdt = " + to_string(dxdt[0]) + " " + to_string(dxdt[1]) + " " + 
		   to_string(dxdt[2]) + " " + to_string(dxdt[3]) + " " + "\n" );
    return( GSL_SUCCESS );
}


int ParticleP2D::trajectory_intersections_at_plane( std::vector<ParticleP2D> &intsc, 
						    int crd, double val,
						    const ParticleP2D &x1, const ParticleP2D &x2,
						    int extrapolate )
{
    // Construct trajectory interpolation
    double dt = x2[0] - x1[0];
    TrajectoryRep1D trep[2];
    trep[0].construct( dt, x1[1], x1[2], x2[1], x2[2] );
    trep[1].construct( dt, x1[3], x1[4], x2[3], x2[4] );

    // Solve for intersections
    double K[3];
    int nroots = trep[crd].solve( K, val, extrapolate );

    // Save intersection points
    double x[2], v[2];
    for( int b = 0; b < nroots; b++ ) {
	trep[0].coord( x[0], v[0], K[b] );
	trep[1].coord( x[1], v[1], K[b] );
	intsc.push_back( ParticleP2D( x1[0]+dt*K[b], x[0], v[0], x[1], v[1] ) );
    }

    return( nroots );
}


int ParticlePCyl::get_derivatives( double t, const double *x, double *dxdt, void *data )
{
    Vec3D E, B, xc( x[0], x[2], 0.0 );
    ParticleIteratorData *pidata = (ParticleIteratorData *)data;

    // Prevent sampling field data outside simulation box, allow h
    // tolerance except at r=0 boundary
    if( xc[0] < pidata->_geom->origo(0)-pidata->_geom->h() )
	return( IBSIMU_DERIV_ERROR );
    else if( xc[0] > pidata->_geom->max(0)+pidata->_geom->h() )
	return( IBSIMU_DERIV_ERROR );
    else if( xc[1] <= 0.0 )
	return( IBSIMU_DERIV_ERROR );
    else if( xc[1] > pidata->_geom->max(1)+pidata->_geom->h() )
	return( IBSIMU_DERIV_ERROR );
    
    if( pidata->_efield )
	E = (*pidata->_efield)( xc );
    if( pidata->_bfield ) {
	B = (*pidata->_bfield)( xc );
	if( pidata->_bsup_cb ) {
	    double factor = (*pidata->_bsup_cb)( xc );
	    B *= factor;
	}
    }

    if( comp_isnan(E[0]) || comp_isnan(E[1]) )
	return( IBSIMU_DERIV_ERROR );

    /* Positions: dx/dt = vx, dr/dt = vr */
    dxdt[0] = x[1];
    dxdt[2] = x[3];
    
    /* Velocities:
     * dvx/dt = ax
     * dvr/dt = ar+r*(dtheta/dt)^2 
     * d^2theta/dt^2 = (a_theta-dr/dt*dtheta/dt)/r 
     */
    if( pidata->_relativistic ) {
	throw( ErrorUnimplemented( ERROR_LOCATION, "Relativistic particle iteration unimplemented" ) );
    } else {
	dxdt[1] = pidata->_qm * (E[0] + x[3]*B[2] - x[2]*x[4]*B[1]);
	dxdt[3] = pidata->_qm * (E[1] + x[2]*x[4]*B[0] - x[1]*B[2]) + x[2]*x[4]*x[4];
	if( x[2] == 0.0 )
	    dxdt[4] = 0.0;
	else
	    dxdt[4] = (pidata->_qm * (x[1]*B[1] - x[3]*B[0]) - 2.0*x[3]*x[4]) / x[2];
    }
    
    return( GSL_SUCCESS );
}


int ParticlePCyl::trajectory_intersections_at_plane( std::vector<ParticlePCyl> &intsc, 
						     int crd, double val,
						     const ParticlePCyl &x1, const ParticlePCyl &x2,
						     int extrapolate )
{
    // Construct trajectory interpolation
    double dt = x2[0] - x1[0];
    TrajectoryRep1D trep[2];
    trep[0].construct( dt, x1[1], x1[2], x2[1], x2[2] );
    trep[1].construct( dt, x1[3], x1[4], x2[3], x2[4] );

    // Solve for intersections
    double K[3];
    int nroots = trep[crd].solve( K, val, extrapolate );
    
    // Save intersection points
    double x[2], v[2];
    for( int b = 0; b < nroots; b++ ) {
	trep[0].coord( x[0], v[0], K[b] );
	trep[1].coord( x[1], v[1], K[b] );
	intsc.push_back( ParticlePCyl( x1[0]+K[b]*dt, x[0], v[0], x[1], v[1], x1[5]+K[b]*(x2[5]-x1[5]) ) );
    }

    return( nroots );
}


int ParticleP3D::get_derivatives( double t, const double *x, double *dxdt, void *data )
{
    Vec3D E, B, xc( x[0], x[2], x[4] );
    ParticleIteratorData *pidata = (ParticleIteratorData *)data;

    // Prevent sampling field data outside simulation box, allow h tolerance
    for( uint32_t a = 0; a < 3; a++ ) {
	if( xc[a] < pidata->_geom->origo(a)-pidata->_geom->h() )
	    return( IBSIMU_DERIV_ERROR );
	else if( xc[a] > pidata->_geom->max(a)+pidata->_geom->h() )
	    return( IBSIMU_DERIV_ERROR );
    }
    
    if( pidata->_efield )
	E = (*pidata->_efield)( xc );
    if( pidata->_bfield ) {
	B = (*pidata->_bfield)( xc );
	if( pidata->_bsup_cb ) {
	    double factor = (*pidata->_bsup_cb)( xc );
	    B *= factor;
	}
    }

    if( comp_isnan(E[0]) || comp_isnan(E[1]) || comp_isnan(E[2]) )
	return( IBSIMU_DERIV_ERROR );
    if( !comp_isfinite(B[0]) || !comp_isfinite(B[1]) || !comp_isfinite(B[2]) )
	return( IBSIMU_DERIV_ERROR );

    /* Positions: dx/dt = vx, dy/dt = vy, dz/dt = vz */
    dxdt[0] = x[1];
    dxdt[2] = x[3];
    dxdt[4] = x[5];
    
    /* Velocities: dvx/dt = ax, dvy/dt = ay, dvz/dt = az */
    if( pidata->_relativistic ) {
	// CPI patch on 20260604: avoid NaN timestep in relativistic tracing.
	double v2 = x[1]*x[1] + x[3]*x[3] + x[5]*x[5];
    if( !comp_isfinite(v2) || v2 >= SPEED_C2 ) {
        return( IBSIMU_DERIV_ERROR );
    }
    double inv_gamma2 = 1.0 - v2/SPEED_C2;
    if( !comp_isfinite(inv_gamma2) || inv_gamma2 <= 0.0 ) {
        return( IBSIMU_DERIV_ERROR );
    }
    double gamma = 1.0/sqrt( inv_gamma2 );
    if( !comp_isfinite(gamma) ) {
        return( IBSIMU_DERIV_ERROR );
    }
	double gamma3c = gamma*gamma*gamma/SPEED_C2;
    if( !comp_isfinite(gamma3c) ) {
        return( IBSIMU_DERIV_ERROR );
    }
	double a12 = gamma3c*x[1]*x[3];
	double a13 = gamma3c*x[1]*x[5];
	double a23 = gamma3c*x[3]*x[5];
	Mat3D m( gamma3c*x[1]*x[1] + gamma, a12, a13,
		 a12, gamma3c*x[3]*x[3] + gamma, a23,
		 a13, a23, gamma3c*x[5]*x[5] + gamma );
    double mdet = m.determinant();
    if( !comp_isfinite(mdet) || mdet == 0.0 ) {
        return( IBSIMU_DERIV_ERROR );
    }
    Mat3D minv = m.inverse( mdet );
	Vec3D pm( pidata->_qm * (E[0] + x[3]*B[2] - x[5]*B[1]),
		  pidata->_qm * (E[1] + x[5]*B[0] - x[1]*B[2]),
		  pidata->_qm * (E[2] + x[1]*B[1] - x[3]*B[0]) );
    if( !comp_isfinite(pm[0]) || !comp_isfinite(pm[1]) || !comp_isfinite(pm[2]) ) {
        return( IBSIMU_DERIV_ERROR );
    }
	Vec3D dvdt = minv*pm;
    if( !comp_isfinite(dvdt[0]) || !comp_isfinite(dvdt[1]) || !comp_isfinite(dvdt[2]) ) {
        return( IBSIMU_DERIV_ERROR );
    }
	dxdt[1] = dvdt(0);
	dxdt[3] = dvdt(1);
	dxdt[5] = dvdt(2);
    } else {
	dxdt[1] = pidata->_qm * (E[0] + x[3]*B[2] - x[5]*B[1]);
	dxdt[3] = pidata->_qm * (E[1] + x[5]*B[0] - x[1]*B[2]);
	dxdt[5] = pidata->_qm * (E[2] + x[1]*B[1] - x[3]*B[0]);
    }

    return( GSL_SUCCESS );
}


int ParticleP3D::trajectory_intersections_at_plane( std::vector<ParticleP3D> &intsc, 
						    int crd, double val,
						    const ParticleP3D &x1, const ParticleP3D &x2,
						    int extrapolate )
{
    // Construct trajectory interpolation
    double dt = x2[0] - x1[0];
    TrajectoryRep1D trep[3];
    trep[0].construct( dt, x1[1], x1[2], x2[1], x2[2] );
    trep[1].construct( dt, x1[3], x1[4], x2[3], x2[4] );
    trep[2].construct( dt, x1[5], x1[6], x2[5], x2[6] );

    // Solve for intersections
    double K[3];
    int nroots = trep[crd].solve( K, val, extrapolate );
    
    // Save intersection points
    double x[3], v[3];
    for( int b = 0; b < nroots; b++ ) {
	trep[0].coord( x[0], v[0], K[b] );
	trep[1].coord( x[1], v[1], K[b] );
	trep[2].coord( x[2], v[2], K[b] );
	intsc.push_back( ParticleP3D( x1[0]+dt*K[b], x[0], v[0], x[1], v[1], x[2], v[2] ) );
    }

    return( nroots );
}

