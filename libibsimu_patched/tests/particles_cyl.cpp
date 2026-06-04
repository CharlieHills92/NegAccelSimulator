/*! \file particles_cyl.cpp 
 *  \brief Test particle iterator in zero field in cylindrical coordinates.
 *
 *  \test Test particle iterator in zero field in cylindrical coordinates.
 */

/*
 *  Test particle iterator in zero field in cylindrical coordinates.
 *  The trajectory is linear in (x,y) coordinates and is compared to cylindrical
 *  coordinates, which are
 *  \f[ r = \sqrt{ x^2 + y^2 }, \f]
 *  and
 *  \f[ \frac{\partial r}{\partial t} = \frac{2x \frac{\partial x}{\partial t} + 
    2y \frac{\partial y}{\partial t}}{2 \sqrt{ x^2 + y^2 }}. \f]
 *  Key points and energy conservation are checked.
 *
 */


#include <cstdlib>
#include <fstream>
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


using namespace std;


double vrfunc( double y, double vy, double z, double vz )
{
    return( (y*vy+z*vz)/sqrt(y*y+z*z) );
}


double wfunc( double y, double vy, double z, double vz )
{
    return( (z*vy-y*vz)/(y*y+z*z) );
}


void test( int argc, char **argv )
{
    Geometry geom( MODE_CYL, Int3D(11,11,1), Vec3D(-0.05,0.0,0.0), 0.01 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,    0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,      0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,    0.0) );
    geom.build_mesh();

    EpotField epot( geom );
    MeshScalarField scharge( geom );

    EpotBiCGSTABSolver solver( geom );
    solver.solve( epot, scharge );

    EpotEfield efield( epot );
    MeshVectorField bfield;

    ParticleDataBaseCyl pdb( geom );
    pdb.set_accuracy( 1.0e-6, 1.0e-6 ); // Extremely sensitive on error requirement
    pdb.set_thread_count( 1 );
    pdb.set_max_steps( 100 );
    double y0 = 0.003;
    double vy = -1.5e5;
    double z0 = 0.002;
    double vz = -2.0e5;
    double r0 = sqrt(y0*y0 + z0*z0);
    double vr = vrfunc( y0, vy, z0, vz );
    double w = wfunc( y0, vy, z0, vz );
    pdb.add_particle( 0.0, 1.0, 1.0, ParticlePCyl( 0, 0.0, 0.0, r0, vr, w ) );
    pdb.iterate_trajectories( scharge, efield, bfield );
    //pdb.debug_print();

    ParticleCyl &part = pdb.particle(0);
    ofstream ostr( "particles_cyl.dat" );
    //ostream &ostr = cout;
    ostr << "# "
	 << setw(62) << "Analytic" << " " 
	 << setw(38) << "Simulated" << "\n"; 
    ostr << "# "
	 << setw(10) << "y (m)" << " " 
	 << setw(12) << "z (m)" << " " 
	 << setw(12) << "r (m)" << " "
	 << setw(12) << "vr (m/s)" << " "
	 << setw(12) << "w (rad/s)" << " "

	 << setw(12) << "r (m)" << " "
	 << setw(12) << "vr (m/s)" << " "
	 << setw(12) << "w (rad/s)" << "\n";
    bool err = false;
    for( uint32_t b = 0; b < part.traj_size(); b++ ) {
	double t = part.traj(b)(0);
	double y = y0 + vy * t;
	double z = z0 + vz * t;
	ostr << setw(12) << y << " " 
	     << setw(12) << z << " " 
	     << setw(12) << sqrt(y*y+z*z) << " " 
	     << setw(12) << vrfunc(y,vy,z,vz) << " "
	     << setw(12) << wfunc(y,vy,z,vz) << " ";

	double r = part.traj(b)(3);
	double vr2 = part.traj(b)(4);
	double w2 = part.traj(b)(5);
	ostr << setw(12) << r << " "
	     << setw(12) << vr2 << " "
	     << setw(12) << w2 << "\n";
	
	if( fabs( sqrt(y*y+z*z) - r ) > 1e-6 ) {
	    err = true;
	}
    }

    ostr.close();
    if( err ) {
	std::cout << "Error: trajectory differs from theory\n";
	exit( 1 );
    }
}


