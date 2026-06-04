/*! \file particles_mirror.cpp
 *  \brief Test particle mirroring
 *
 *  \test Test particle mirroring
 *
 *  The test case has a 1000 V potential difference between
 *  x=-0.05 m and x=0.05 m planes. The electric field is therefore
 *  -1000 V / 0.1 m = -10000 V/m. The particle trajectories should
 *  follow trajectories defined by constant acceleration motion
 *  \f[ x = x_0 + v_0 t + \frac{1}{2} a_x t^2, \f]
 *  where the acceleration is
 *  \f[ a_x = \frac{F_x}{m} = \frac{Eq}{m}. \f]
 *  In this test case the particles have q = 1e, m = 1u. Therefore 
 *  \f$a_x = -9.64853082148e11\mathrm{~m/s}^2 \f$
 */


#include <iostream>
#include <iomanip>
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "meshvectorfield.hpp"
#include "epot_bicgstabsolver.hpp"
#include "epot_efield.hpp"
#include "particles.hpp"
#include "error.hpp"
#include "ibsimu.hpp"


using namespace std;


const double ax = -9.64853082148e11;
double maxerr = 0.0;
bool err = false;
bool err_energy = false;


void check_particle( Particle2D &p, double x0, double vx0, double y0, double vy0, EpotField &epot )
{
    // Check energy conservation
    double Eref = 0.0;
    //std::cout << "Checking energy conservation\n";
    for( uint32_t b = 0; b < p.traj_size(); b++ ) {
	double Ek = 0.5*MASS_U*( p.traj(b)(2)*p.traj(b)(2) + 
				 p.traj(b)(4)*p.traj(b)(4) );
	double Ep = CHARGE_E*epot( Vec3D(p.traj(b)(1), p.traj(b)(3), 0.0) );
	if( b == 0 )
	    Eref = Ek + Ep;
	else if( fabs( (Ek + Ep - Eref)/Eref ) > 4.0e-4 ) {
	    err_energy = true;
	}
	//std::cout << "x = " << p.traj(b)(1) << "\tfabs((Ek+Ep-Eref)/Eref) = " << fabs( (Ek + Ep - Eref)/Eref ) << "\n";
    }

    // Check trajectory
    //std::cout << "\nChecking trajectory\n";
    for( uint32_t b = 0; b < p.traj_size(); b++ ) {
	double t = p.traj(b)(0);
	if( t < 3.21935897726e-7 - 0.1e-7 ) {
	    double x = x0 + vx0 * t + 0.5*ax*t*t;
	    //double vx = vx0 + ax*t;
	    double y = y0 + vy0 * t;
	    //std::cout << x << "\t" << p.traj(b)(1) << "\t"
	    //<< y << "\t" << p.traj(b)(3) << "\t"
	    //<< vx << "\t" << p.traj(b)(2) << "\n";
	    double xerr = fabs( p.traj(b)(1) - x );
	    double yerr = fabs( p.traj(b)(3) - y );
	    //double vxerr = fabs( p.traj(b)(2) - vx );
	    if( xerr > 5e-5 || yerr > 5e-5 )
		err = true;
	    //std::cout << "xerr = " << xerr << "\n";
	    //std::cout << "yerr = " << yerr << "\n";
	    //std::cout << "vxerr = " << vxerr << "\n";
	} else if( t > 3.21935897726e-7 + 0.1e-7 ) {
	    double t2 = 3.21935897726e-7;
	    double new_x0 = x0 + vx0 * t2 + 0.5*ax*t2*t2;
	    double new_y0 = y0 + vy0 * t2;
	    double new_vx0 = -vx0 - ax*t2;
	    t = t - t2;
	    double x = new_x0 + new_vx0 * t + 0.5*ax*t*t;
	    //double vx = new_vx0 + ax*t;
	    double y = new_y0 + vy0 * t;
	    double xerr = fabs( p.traj(b)(1) - x );
	    double yerr = fabs( p.traj(b)(3) - y );
	    //double vxerr = fabs( p.traj(b)(2) - vx );
	    //std::cout << x << "\t" << p.traj(b)(1) << "\t"
	    //<< y << "\t" << p.traj(b)(3) << "\t"
	    //<< vx << "\t" << p.traj(b)(2) << "\n";
	    if( xerr > 5e-5 || yerr > 5e-5 ) 
		err = true;
	    //std::cout << "xerr = " << xerr << "\n";
	    //std::cout << "yerr = " << yerr << "\n";
	    //std::cout << "vxerr = " << vxerr << "\n";
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
    //g.debug_print();

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;

    EpotBiCGSTABSolver solver( geom );
    solver.solve( epot, scharge );

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_ANTIMIRROR, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase2D pdb( geom );
    pdb.set_thread_count( 1 );
    bool pmirror[6] = { true, false, false, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_accuracy( 1.0e-6, 1.0e-6 ); // 1.0e-6: 4e-4 error in energy
                                        // 1.0e-7: 1e-5 error in energy
                                        // 1.0e-8: 8e-7 error in energy
    pdb.set_save_all_points( true );
    pdb.set_trajectory_interpolation( TRAJECTORY_INTERPOLATION_POLYNOMIAL );
    pdb.add_particle( 0.0, 1.0, 1.0, ParticleP2D( 0.0, 0.0, 0.0, -0.04, 1e5 ) );
    pdb.iterate_trajectories( scharge, efield, bfield );
    //pdb.debug_print( cout );

    // Check particle trajectory points
    check_particle( pdb.particle(0), 0.0, 0.0, -0.04, 1e5, epot );

    GeomPlotter geomplotter( geom );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.plot_png( "particles_mirror.png" );

    if( err )
	std::cout << "Error: trajectory differs from theory\n";
    if( err_energy )
	std::cout << "Error: energy not conserved\n";
    if( err || err_energy )
	exit( 1 );
}


