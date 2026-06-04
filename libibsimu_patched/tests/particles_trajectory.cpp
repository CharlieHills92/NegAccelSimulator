/*! \file particles_trajectory.cpp 
 *  \brief Test particle iterator in constant electric field.
 *
 *  \test Test particle iterator in constant electric field.
 */

/*
 *  The test case has a 1000 V potential difference between
 *  x=-0.05 m and x=0.05 m planes. The electric field is therefore
 *  -1000 V / 0.1 m = -10000 V/m. The particle trajectories should
 *  follow trajectories defined by constant acceleration motion
 *  \f[ x = x_0 + v_0 t + \frac{1}{2} a_x t^2, \f]
 *  where the acceleration is
 *  \f[ a_x = \frac{F_x}{m} = \frac{Eq}{m}. \f]
 *  In this test case the particles have q = 1e, m = 1u. Therefore 
 *  \f$a_x = -9.64853082148e11\mathrm{~m/s}^2 \f$
 *  The particles start at origo with different starting velocity 
 *  vectors. Key points and energy conservation are checked.
 *
 */


#include <iostream>
#include <cstdlib>
#include <iomanip>
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "meshvectorfield.hpp"
#include "epot_bicgstabsolver.hpp"
#include "epot_efield.hpp"
#include "particles.hpp"
#include "particledatabase.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "config.h"

#ifdef GTK3
#include "gtkplotter.hpp"
#endif

using namespace std;


const double ax = -9.64853082148e11;


void check_particle( Particle2D &p, double vx, double vy )
{
    for( uint32_t b = 0; b < p.traj_size(); b++ ) {
	double t = p.traj(b)(0);
	double x = vx * t + 0.5*ax*t*t;
	double y = vy * t;
	if( fabs( p.traj(b)(1) - x ) > 1e-6 ||
	    fabs( p.traj(b)(3) - y ) > 1e-6 ) {
	    std::cout << "Error: trajectory differs from theory\n";
	    exit( 1 );
	}
	    
    }

}


void test( int argc, char **argv )
{
    Geometry geom( MODE_2D, Int3D(11,11,1), Vec3D(-0.05,-0.05,0.0), 0.01 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 1000.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,      0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,      0.0) );
    geom.build_mesh();

    EpotField epot( geom );
    MeshVectorField bfield;
    MeshScalarField scharge( geom );

    EpotBiCGSTABSolver solver( geom );
    solver.solve( epot, scharge );

    EpotEfield efield( epot );

    ParticleDataBase2D pdb( geom );
    pdb.set_thread_count( 1 );
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0, 0.0,  0.0, 0.0,  1e5 ) );
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0, 0.0,  1e5, 0.0,  0.0 ) );
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0, 0.0,  0.0, 0.0, -1e5 ) );
    pdb.iterate_trajectories( scharge, efield, bfield );

    // Check particle trajectory points
    check_particle( pdb.particle(0),  0.0,  1e5 );
    check_particle( pdb.particle(1),  1e5,  0.0 );
    check_particle( pdb.particle(2),  0.0, -1e5 );

    GeomPlotter geomplotter( geom );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.plot_png( "particles_trajectory.png" );

#ifdef GTK3
    if( false ) {
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.set_scharge( &scharge );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
    }
#endif
}


