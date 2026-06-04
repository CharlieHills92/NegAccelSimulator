/*! \file scharge.cpp
 *  \brief Space charge deposition functions
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

#include <limits>

#include "trajectory.hpp"
#include "scharge.hpp"
#include "ibsimu.hpp"


//#define DEBUG_SCHARGE 1


#ifdef DEBUG_SCHARGE
#define DEBUG_MESSAGE(x) ibsimu.message(MSG_DEBUG_GENERAL,1) << x
#define DEBUG_INC_INDENT() ibsimu.inc_indent()
#define DEBUG_DEC_INDENT() ibsimu.dec_indent()
#else
#define DEBUG_MESSAGE(x) do {} while(0)
#define DEBUG_INC_INDENT() do {} while(0)
#define DEBUG_DEC_INDENT() do {} while(0)
#endif


void scharge_finalize_pic( MeshScalarField &scharge )
{
    ibsimu.message( 1 ) << "Finalizing space charge density map (PIC method)\n";
    ibsimu.inc_indent();

    switch( scharge.geom_mode() ) {
    case MODE_2D:
    {
	// Convert charge map to space charge density map
	scharge /= (scharge.h()*scharge.h());

	// Correct boundaries
	for( uint32_t i = 0; i < scharge.size(0); i++ ) {
	    scharge( i, 0 ) *= 2.0;
	    scharge( i, scharge.size(1)-1 ) *= 2.0;
	}
	for( uint32_t j = 0; j < scharge.size(1); j++ ) {
	    scharge( 0, j ) *= 2.0;
	    scharge( scharge.size(0)-1, j ) *= 2.0;
	}
	break;
    }
    case MODE_CYL:
    {
	// Convert charge map to space charge density map
 	for( uint32_t i = 0; i < scharge.size(0); i++ ) {
	    for( uint32_t j = 0; j < scharge.size(1); j++ ) {
		if( j == 0 ) {
                    double rj2 = scharge.h()+scharge.origo(1);
                    scharge( i, j ) /= (M_PI*scharge.h()*(rj2*rj2));
                } else {
                    double rj1 = (j-0.5)*scharge.h()+scharge.origo(1);
                    double rj2 = (j+0.5)*scharge.h()+scharge.origo(1);
                    scharge( i, j ) /= (M_PI*scharge.h()*(rj2*rj2-rj1*rj1));
                }
	    }
	}

	// Correct boundaries
	for( uint32_t i = 0; i < scharge.size(0); i++ ) {
	    scharge( i, 0 ) *= 2.0;
	    scharge( i, scharge.size(1)-1 ) *= 2.0;
	}
	for( uint32_t j = 0; j < scharge.size(1); j++ ) {
	    scharge( 0, j ) *= 2.0;
	    scharge( scharge.size(0)-1, j ) *= 2.0;
	}
	break;
    }
    case MODE_3D:
    {
	// Convert charge map to space charge density map
	scharge /= (scharge.h()*scharge.h()*scharge.h());

	// Correct boundaries
 	for( uint32_t i = 0; i < scharge.size(0); i++ ) {
	    for( uint32_t j = 0; j < scharge.size(1); j++ ) {
		scharge( i, j, 0 ) *= 2.0;
		scharge( i, j, scharge.size(2)-1 ) *= 2.0;
	    }
	}
	for( uint32_t i = 0; i < scharge.size(0); i++ ) {
	    for( uint32_t k = 0; k < scharge.size(2); k++ ) {
		scharge( i, 0, k ) *= 2.0;
		scharge( i, scharge.size(1)-1, k ) *= 2.0;
	    }
	}	
	for( uint32_t j = 0; j < scharge.size(1); j++ ) {
	    for( uint32_t k = 0; k < scharge.size(2); k++ ) {
		scharge( 0, j, k ) *= 2.0;
		scharge( scharge.size(0)-1, j, k ) *= 2.0;
	    }
	}	
	break;
    }
    default:
    {
	throw( Error( ERROR_LOCATION, "unsupported dimension number" ) );
    }
    }

    ibsimu.dec_indent();
}


void scharge_finalize_linear( MeshScalarField &scharge )
{
    ibsimu.message( 1 ) << "Finalizing space charge density map (LINEAR method)\n";
    ibsimu.inc_indent();

    switch( scharge.geom_mode() ) {
    case MODE_2D:
    {
	// Convert to space charge density map
	scharge /= scharge.h();
	break;
    }
    default:
    {
	throw( Error( ERROR_LOCATION, "unsupported geometry mode" ) );
    }
    }

    ibsimu.dec_indent();
}

	
void scharge_add_from_trajectory_pic( MeshScalarField &scharge, pthread_mutex_t *mutex, 
				      double I, const ParticleP2D &x1, const ParticleP2D &x2 )
{
    TrajectoryRep1D traj[2];
    double x[2];
    double v[2];
    double t[2];
    int32_t i[2];

    DEBUG_MESSAGE( "Calculating space charge\n" );
    DEBUG_INC_INDENT();
    DEBUG_MESSAGE( "x1 = " << x1 << "\n" );
    DEBUG_MESSAGE( "x2 = " << x2 << "\n" );

    double dt = x2[0]-x1[0];
    for( size_t a = 0; a < 2; a++ ) {
	traj[a].construct( dt, x1[2*a+1], x1[2*a+2], x2[2*a+1], x2[2*a+2] );
	traj[a].coord( x[a], v[a], 0.5 );
	i[a] = (int32_t)floor( ( x[a]-scharge.origo(a) ) * scharge.div_h() );
	t[a] = ( x[a]-(i[a]*scharge.h()+scharge.origo(a)) ) * scharge.div_h();

	DEBUG_MESSAGE( "a = " << a << "\n" );
	DEBUG_MESSAGE( "x = " << x[a] << "\n" );
	DEBUG_MESSAGE( "i = " << i[a] << "\n" );
	DEBUG_MESSAGE( "t = " << t[a] << "\n" );

	// Add charge to boundaries when over simulation area
	if( i[a] < 0 ) {
	    i[a] = 0;
	    t[a] = 0.0;
	} else if( i[a] >= (int32_t)scharge.size(a)-1 ) {
	    i[a] = scharge.size(a)-2;
	    t[a] = 1.0;
	}
    }

    double Q = I*dt;
    int p = scharge.size(0)*i[1] + i[0];
    pthread_mutex_lock( mutex );
    scharge( p )                   += (1.0-t[0])*(1.0-t[1])*Q;
    scharge( p+scharge.size(0) )   += (1.0-t[0])*t[1]*Q;
    scharge( p+1 )                 += t[0]*(1.0-t[1])*Q;
    scharge( p+1+scharge.size(0) ) += t[0]*t[1]*Q;
    pthread_mutex_unlock( mutex );

    DEBUG_DEC_INDENT();
}


void scharge_add_from_trajectory_pic( MeshScalarField &scharge, pthread_mutex_t *mutex, 
				      double I, const ParticlePCyl &x1, const ParticlePCyl &x2 )
{
    TrajectoryRep1D traj[2];
    double x[2];
    double v[2];
    double t[2];
    int32_t i[2];

    double dt = x2[0]-x1[0];

    // x-direction
    traj[0].construct( dt, x1[1], x1[2], x2[1], x2[2] );
    traj[0].coord( x[0], v[0], 0.5 );
    i[0] = (int32_t)floor( ( x[0]-scharge.origo(0) ) * scharge.div_h() );
    t[0] = ( x[0]-(i[0]*scharge.h()+scharge.origo(0)) ) * scharge.div_h();

    // r-direction
    traj[1].construct( dt, x1[3], x1[4], x2[3], x2[4] );
    traj[1].coord( x[1], v[1], 0.5 );
    i[1] = (int32_t)floor( ( x[1]-scharge.origo(1) ) * scharge.div_h() );
    double rj1 = i[1]*scharge.h()+scharge.origo(1);
    double rj2 = rj1+scharge.h();
    rj1 = rj1*rj1;
    rj2 = rj2*rj2;
    t[1] = (x[1]*x[1]-rj1) / (rj2-rj1);

    for( size_t a = 0; a < 2; a++ ) {
	// Add charge to boundaries when over simulation area
	if( i[a] < 0 ) {
	    i[a] = 0;
	    t[a] = 0.0;
	} else if( i[a] >= (int32_t)scharge.size(a)-1 ) {
	    i[a] = scharge.size(a)-2;
	    t[a] = 1.0;
	}
    }

    double Q = I*dt;
    int p = scharge.size(0)*i[1] + i[0];
    pthread_mutex_lock( mutex );
    scharge( p )                   += (1.0-t[0])*(1.0-t[1])*Q;
    scharge( p+scharge.size(0) )   += (1.0-t[0])*t[1]*Q;
    scharge( p+1 )                 += t[0]*(1.0-t[1])*Q;
    scharge( p+1+scharge.size(0) ) += t[0]*t[1]*Q;    
    pthread_mutex_unlock( mutex );
}


void scharge_add_from_trajectory_pic( MeshScalarField &scharge, pthread_mutex_t *mutex, 
				      double I, const ParticleP3D &x1, const ParticleP3D &x2 )
{
    TrajectoryRep1D traj[3];
    double x[3];
    double v[3];
    double t[3];
    int32_t i[3];

    double dt = x2[0]-x1[0];
    for( size_t a = 0; a < 3; a++ ) {
	traj[a].construct( dt, x1[2*a+1], x1[2*a+2], x2[2*a+1], x2[2*a+2] );
	traj[a].coord( x[a], v[a], 0.5 );
	//x[a] = 0.5*(x1[2*a+1]+x2[2*a+1]);
	i[a] = (int32_t)floor( ( x[a]-scharge.origo(a) ) * scharge.div_h() );
	t[a] = ( x[a]-(i[a]*scharge.h()+scharge.origo(a)) ) * scharge.div_h();

	// Add charge to boundaries when over simulation area
	if( i[a] < 0 ) {
	    i[a] = 0;
	    t[a] = 0.0;
	} else if( i[a] >= (int32_t)scharge.size(a)-1 ) {
	    i[a] = scharge.size(a)-2;
	    t[a] = 1.0;
	}
    }

    double Q = I*dt;
    int p = scharge.size(0)*scharge.size(1)*i[2] + scharge.size(0)*i[1] + i[0];

    pthread_mutex_lock( mutex );
    scharge( p )                   += (1.0-t[0])*(1.0-t[1])*(1.0-t[2])*Q;
    scharge( p+scharge.size(0) )   += (1.0-t[0])*t[1]*(1.0-t[2])*Q;
    scharge( p+1 )                 += t[0]*(1.0-t[1])*(1.0-t[2])*Q;
    scharge( p+1+scharge.size(0) ) += t[0]*t[1]*(1.0-t[2])*Q;

    p += scharge.size(0)*scharge.size(1);
    scharge( p )                   += (1.0-t[0])*(1.0-t[1])*t[2]*Q;
    scharge( p+scharge.size(0) )   += (1.0-t[0])*t[1]*t[2]*Q;
    scharge( p+1 )                 += t[0]*(1.0-t[1])*t[2]*Q;
    scharge( p+1+scharge.size(0) ) += t[0]*t[1]*t[2]*Q;
    pthread_mutex_unlock( mutex );
}


double distance( const double x[2], const ParticleP2D &p )
{
    double t1 = x[0]-p[1];
    double t2 = x[1]-p[3];
    return( sqrt( t1*t1 + t2*t2 ) );
}


/* Find closest distance from x to line segment between particle
 * points p1 and p2. Return the distance and set interpolated particle
 * coordinates at this location in intrp. Uses linear interpolation.
 */
double closest_point( ParticleP2D &intrp, const double x[2], 
		      const ParticleP2D &p1, const ParticleP2D &p2 )
{
    double t1 = p1[1]-p2[1];
    double t2 = p1[3]-p2[3];
    double l2 = t1*t1 + t2*t2;  // l2=|p1-p2|^2
    if( l2 == 0.0 ) {
	intrp = p1;
	return( distance( x, p1 ) );
    }
		
    // Use parametric presentation of line: p1 + t*(p2-p1)
    // Projection of point x to line is at t = [(x-p1) . (p2-p1)] / |p2-p1|^2
    double d1[2] = { x[0]-p1[1], x[1]-p1[3] };
    double d2[2] = { p2[1]-p1[1], p2[3]-p1[3] };
    double t = ( d1[0]*d2[0] + d1[1]*d2[1] ) / l2;
    if( t < 0.0 ) {
	intrp = p1;
	return( distance( x, p1 ) );
    } else if( t > 1.0 ) {
	intrp = p2;
	return( distance( x, p2 ) );
    }
    intrp = p1 + t*(p2-p1);
    return( distance( x, intrp ) );
}


double distance( const double x[3], const ParticleP3D &p )
{
    double t1 = x[0]-p[1];
    double t2 = x[1]-p[3];
    double t3 = x[2]-p[5];
    return( sqrt( t1*t1 + t2*t2 + t3*t3 ) );
}


/* Find closest distance from x to line segment between particle
 * points p1 and p2. Return the distance and set interpolated particle
 * coordinates at this location in intrp. Uses linear interpolation.
 */
double closest_point( ParticleP3D &intrp, const double x[3], 
		      const ParticleP3D &p1, const ParticleP3D &p2 )
{
    double t1 = p1[1]-p2[1];
    double t2 = p1[3]-p2[3];
    double t3 = p1[5]-p2[5];
    double l2 = t1*t1 + t2*t2 + t3*t3;  // l2=|p1-p2|^2
    if( l2 == 0.0 ) {
	intrp = p1;
	return( distance( x, p1 ) );
    }
    
    // Use parametric presentation of line: p1 + t*(p2-p1)
    // Projection of point x to line is at t = [(x-p1) . (p2-p1)] / |p2-p1|^2
    double d1[3] = { x[0]-p1[1], x[1]-p1[3], x[2]-p1[5] };
    double d2[3] = { p2[1]-p1[1], p2[3]-p1[3], p2[5]-p1[5] };
    double t = ( d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2] ) / l2;
    if( t < 0.0 ) {
	intrp = p1;
	return( distance( x, p1 ) );
    } else if( t > 1.0 ) {
	intrp = p2;
	return( distance( x, p2 ) );
    }
    intrp = p1 + t*(p2-p1);
    return( distance( x, intrp ) );
}


void scharge_add_from_trajectory_linear( MeshScalarField &scharge, pthread_mutex_t *mutex, 
					 double I, int dir, const CFiFo<ParticleP2D,4> &cdpast, const int i[3] )
{
#ifdef DEBUG_SCHARGE
    std::cout << "scharge_add_from_trajectory\n";
    std::cout << "I = " << I << "\n";
    std::cout << "dir = " << dir << "\n";
    std::cout << "i = { " << i[0] << ", " << i[1] << ", " << i[2] << " }\n";
#endif

    // Node indices to check
    int ni[2][2];
    if( dir == -1 ) {
	ni[0][0] = i[0]+1;
	ni[0][1] = i[1];
	ni[1][0] = i[0]+1;
	ni[1][1] = i[1]+1;
    } else if( dir == +1 ) {
	ni[0][0] = i[0];
	ni[0][1] = i[1];
	ni[1][0] = i[0];
	ni[1][1] = i[1]+1;
    } else if( dir == -2 ) {
	ni[0][0] = i[0];
	ni[0][1] = i[1]+1;
	ni[1][0] = i[0]+1;
	ni[1][1] = i[1]+1;
    } else {
	ni[0][0] = i[0];
	ni[0][1] = i[1];
	ni[1][0] = i[0]+1;
	ni[1][1] = i[1];
    }
    
    // Process nodes
    ParticleP2D p_closest( 0, 0, 0, 0, 0 );
    for( int b = 0; b < 2; b++ ) {
	double x[2] = { scharge.origo(0)+ni[b][0]*scharge.h(),
		        scharge.origo(1)+ni[b][1]*scharge.h() };
#ifdef DEBUG_SCHARGE
	std::cout << "Finding closest point to node " << ni[b][0] 
		  << ", " << ni[b][1] << "\n";
	std::cout << "x = " << x[0] << ", " << x[1] << "\n";
#endif
	// Find closest point to past trajectory
	double d_closest = std::numeric_limits<double>::infinity();
	for( int c = 1; c < cdpast.size(); c++ ) {
	    ParticleP2D p;
	    double d = closest_point( p, x, cdpast[c-1], cdpast[c] );
	    if( d < d_closest ) {
		p_closest = p;
		d_closest = d;
	    }
#ifdef DEBUG_SCHARGE
	    std::cout << "Segment " << c << "\n";
	    std::cout << "p1 = " << cdpast[c-1] << "\n";
	    std::cout << "p2 = " << cdpast[c] << "\n";
	    std::cout << "d = " << d << "\n";
	    std::cout << "p = " << p << "\n";
#endif
	}

#ifdef DEBUG_SCHARGE
	std::cout << "Done\n";
	std::cout << "d_closest = " << d_closest << "\n";
	std::cout << "p_closest = " << p_closest << "\n";
#endif
	// Use closest point to calculate rho*h=I/v
	if( d_closest > scharge.h() )
	    continue;
	double v = sqrt( p_closest[2]*p_closest[2] + p_closest[4]*p_closest[4] );
	double Q = I*(scharge.h()-d_closest)/(scharge.h()*v);
	pthread_mutex_lock( mutex );
	scharge( ni[b][0], ni[b][1] ) += Q;
	pthread_mutex_unlock( mutex );	
#ifdef DEBUG_SCHARGE
	std::cout << "v = " << v << "\n";
	std::cout << "Q = " << Q << "\n";
#endif
    }
}


void scharge_add_from_trajectory_linear( MeshScalarField &scharge, pthread_mutex_t *mutex, 
					 double I, int dir, const CFiFo<ParticleP3D,4> &cdpast, const int i[3] )
{
#ifdef DEBUG_SCHARGE
    std::cout << "scharge_add_from_trajectory\n";
    std::cout << "I = " << I << "\n";
    std::cout << "dir = " << dir << "\n";
    std::cout << "i = { " << i[0] << ", " << i[1] << ", " << i[2] << " }\n";
#endif

    // Node indices to check
    int ni[4][3];
    if( dir == -1 ) {
	ni[0][0] = i[0]+1;
	ni[0][1] = i[1];
	ni[0][2] = i[2];

	ni[1][0] = i[0]+1;
	ni[1][1] = i[1]+1;
	ni[1][2] = i[2];

	ni[2][0] = i[0]+1;
	ni[2][1] = i[1];
	ni[2][2] = i[2]+1;

	ni[3][0] = i[0]+1;
	ni[3][1] = i[1]+1;
	ni[3][2] = i[2]+1;
    } else if( dir == +1 ) {
	ni[0][0] = i[0];
	ni[0][1] = i[1];
	ni[0][2] = i[2];

	ni[1][0] = i[0];
	ni[1][1] = i[1]+1;
	ni[1][2] = i[2];

	ni[2][0] = i[0];
	ni[2][1] = i[1];
	ni[2][2] = i[2]+1;

	ni[3][0] = i[0];
	ni[3][1] = i[1]+1;
	ni[3][2] = i[2]+1;
    } else if( dir == -2 ) {
	ni[0][0] = i[0];
	ni[0][1] = i[1]+1;
	ni[0][2] = i[2];

	ni[1][0] = i[0]+1;
	ni[1][1] = i[1]+1;
	ni[1][2] = i[2];

	ni[2][0] = i[0];
	ni[2][1] = i[1]+1;
	ni[2][2] = i[2]+1;

	ni[3][0] = i[0]+1;
	ni[3][1] = i[1]+1;
	ni[3][2] = i[2]+1;
    } else if( dir == +2 ) {
	ni[0][0] = i[0];
	ni[0][1] = i[1];
	ni[0][2] = i[2];

	ni[1][0] = i[0]+1;
	ni[1][1] = i[1];
	ni[1][2] = i[2];

	ni[2][0] = i[0];
	ni[2][1] = i[1];
	ni[2][2] = i[2]+1;

	ni[3][0] = i[0]+1;
	ni[3][1] = i[1];
	ni[3][2] = i[2]+1;
    } else if( dir == -3 ) {
	ni[0][0] = i[0];
	ni[0][1] = i[1];
	ni[0][2] = i[2]+1;

	ni[1][0] = i[0]+1;
	ni[1][1] = i[1];
	ni[1][2] = i[2]+1;

	ni[2][0] = i[0];
	ni[2][1] = i[1]+1;
	ni[2][2] = i[2]+1;

	ni[3][0] = i[0]+1;
	ni[3][1] = i[1]+1;
	ni[3][2] = i[2]+1;
    } else {
	ni[0][0] = i[0];
	ni[0][1] = i[1];
	ni[0][2] = i[2];

	ni[1][0] = i[0]+1;
	ni[1][1] = i[1];
	ni[1][2] = i[2];

	ni[2][0] = i[0];
	ni[2][1] = i[1]+1;
	ni[2][2] = i[2];

	ni[3][0] = i[0]+1;
	ni[3][1] = i[1]+1;
	ni[3][2] = i[2];
    }
    
    // Process nodes
    ParticleP3D p_closest( 0, 0, 0, 0, 0, 0, 0 );
    for( int b = 0; b < 4; b++ ) {
	double x[3] = { scharge.origo(0)+ni[b][0]*scharge.h(),
		        scharge.origo(1)+ni[b][1]*scharge.h(),
			scharge.origo(2)+ni[b][2]*scharge.h() };
#ifdef DEBUG_SCHARGE
	std::cout << "Finding closest point to node " << ni[b][0] 
		  << ", " << ni[b][1]
		  << ", " << ni[b][2] << "\n";
	std::cout << "x = " << x[0] << ", " << x[1] << ", " << x[2] << "\n";
#endif
	// Find closest point to past trajectory
	double d_closest = std::numeric_limits<double>::infinity();
	for( int c = 1; c < cdpast.size(); c++ ) {
	    ParticleP3D p;
	    double d = closest_point( p, x, cdpast[c-1], cdpast[c] );
	    if( d < d_closest ) {
		p_closest = p;
		d_closest = d;
	    }
#ifdef DEBUG_SCHARGE
	    std::cout << "Segment " << c << "\n";
	    std::cout << "p1 = " << cdpast[c-1] << "\n";
	    std::cout << "p2 = " << cdpast[c] << "\n";
	    std::cout << "d = " << d << "\n";
	    std::cout << "p = " << p << "\n";
#endif
	}

#ifdef DEBUG_SCHARGE
	std::cout << "Done\n";
	std::cout << "d_closest = " << d_closest << "\n";
	std::cout << "p_closest = " << p_closest << "\n";
#endif
	// Use closest point to calculate rho*h^2=I/v
	// [I/v] = As/m = C/m
	// Is circular ok or should we do bilinear here?
	if( d_closest > scharge.h() )
	    continue;
	double v = sqrt( p_closest[2]*p_closest[2] + 
			 p_closest[4]*p_closest[4] + 
			 p_closest[6]*p_closest[6] );
	double Q = I*(scharge.h()-d_closest)/(scharge.h()*v);
	pthread_mutex_lock( mutex );
	scharge( ni[b][0], ni[b][1], ni[b][2] ) += Q;
	pthread_mutex_unlock( mutex );	
#ifdef DEBUG_SCHARGE
	std::cout << "v = " << v << "\n";
	std::cout << "Q = " << Q << "\n";
#endif
    }    
}


void scharge_add_from_trajectory_linear( MeshScalarField &scharge, pthread_mutex_t *mutex, 
					 double I, int dir, const CFiFo<ParticlePCyl,4> &cdpast, const int i[3] )
{
    throw( Error( ERROR_LOCATION, "unsupported geometry mode" ) );
}



void scharge_finalize_step_pic( MeshScalarField &scharge )
{
    // Convert charge map to space charge density map
    scharge /= (scharge.h()*scharge.h()*scharge.h());

    // Correct boundaries
    for( uint32_t i = 0; i < scharge.size(0); i++ ) {
	for( uint32_t j = 0; j < scharge.size(1); j++ ) {
	    scharge( i, j, 0 ) *= 2.0;
	    scharge( i, j, scharge.size(2)-1 ) *= 2.0;
	}
    }
    for( uint32_t i = 0; i < scharge.size(0); i++ ) {
	for( uint32_t k = 0; k < scharge.size(2); k++ ) {
	    scharge( i, 0, k ) *= 2.0;
	    scharge( i, scharge.size(1)-1, k ) *= 2.0;
	}
    }	
    for( uint32_t j = 0; j < scharge.size(1); j++ ) {
	for( uint32_t k = 0; k < scharge.size(2); k++ ) {
	    scharge( 0, j, k ) *= 2.0;
	    scharge( scharge.size(0)-1, j, k ) *= 2.0;
	}
    }
}


void scharge_add_step_pic( MeshScalarField &scharge, double Q, const Vec3D &x )
{
    int i[3];
    double t[3];

    for( size_t a = 0; a < 3; a++ ) {
	i[a] = (int32_t)floor( ( x[a]-scharge.origo(a) ) * scharge.div_h() );
	t[a] = ( x[a]-(i[a]*scharge.h()+scharge.origo(a)) ) * scharge.div_h();

	// Add charge to boundaries when over simulation area
	if( i[a] < 0 ) {
	    i[a] = 0;
	    t[a] = 0.0;
	} else if( i[a] >= (int32_t)scharge.size(a)-1 ) {
	    i[a] = scharge.size(a)-2;
	    t[a] = 1.0;
	}
    }

    int p = scharge.size(0)*scharge.size(1)*i[2] + scharge.size(0)*i[1] + i[0];

    scharge( p )                   += (1.0-t[0])*(1.0-t[1])*(1.0-t[2])*Q;
    scharge( p+scharge.size(0) )   += (1.0-t[0])*t[1]*(1.0-t[2])*Q;
    scharge( p+1 )                 += t[0]*(1.0-t[1])*(1.0-t[2])*Q;
    scharge( p+1+scharge.size(0) ) += t[0]*t[1]*(1.0-t[2])*Q;

    p += scharge.size(0)*scharge.size(1);
    scharge( p )                   += (1.0-t[0])*(1.0-t[1])*t[2]*Q;
    scharge( p+scharge.size(0) )   += (1.0-t[0])*t[1]*t[2]*Q;
    scharge( p+1 )                 += t[0]*(1.0-t[1])*t[2]*Q;
    scharge( p+1+scharge.size(0) ) += t[0]*t[1]*t[2]*Q;
}
