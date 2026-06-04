/*! \file particles_bfieldcyl.cpp 
 *  \brief Test particle iterator in constant magnetic field in 
 *  cylindrical geometry.
 *
 *  \test Test particle iterator in constant magnetic field in 
 *  cylindrical geometry.
 */


#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "particledatabase.hpp"
#include "epot_bicgstabsolver.hpp"
#include "meshvectorfield.hpp"
#include "epot_efield.hpp"
#include "particles.hpp"
#include "error.hpp"
#include "ibsimu.hpp"


using namespace std;


void test( int argc, char **argv )
{
    Geometry geom( MODE_CYL, Int3D(21,11,1), Vec3D(-0.05,0.0,0.0), 0.01 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,    0.0) );
    geom.set_boundary( 3, Bound(BOUND_DIRICHLET,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,    0.0) );
    geom.build_mesh();
    //g.debug_print();

    EpotField epot( geom );
    MeshScalarField scharge( geom );

    EpotBiCGSTABSolver solver( geom );
    solver.solve( epot, scharge );

    EpotEfield efield( epot );

    // Make magnetic field
    bool fout[3] = {true,false,false};
    MeshVectorField bfield( MODE_CYL, fout, Int3D(2,2,1), Vec3D(-0.05,0.0,0.0), 0.2 );
    bfield.set( 0, 0, Vec3D(1,0,0) );
    bfield.set( 1, 0, Vec3D(1,0,0) );
    bfield.set( 0, 1, Vec3D(1,0,0) );
    bfield.set( 1, 1, Vec3D(1,0,0) );

    ParticleDataBaseCyl pdb( geom );
    pdb.set_thread_count( 1 );
    pdb.set_max_steps( 1000 );
    //pdb.set_max_time( 260e-9 );
    pdb.add_particle( 0.0, 1.0, 1.0, ParticlePCyl( 0, 
						   -0.045, 2e6, 
						   0.020728544449, 0.0, 
						   -2e6/0.020728544449 ) );
    pdb.add_particle( 0.0, 1.0, 1.0, ParticlePCyl( 0, 
						   -0.045, 2e6, 
						   0.04, 0.0, 
						   -2e6/0.04 ) );
    pdb.add_particle( 0.0, 1.0, 1.0, ParticlePCyl( 0, 
						   -0.045, 2e6, 
						   0.01, 0.0, 
						   -2e6/0.01 ) );
    pdb.iterate_trajectories( scharge, efield, bfield );
    //pdb.debug_print();

    // Check particle trajectory points
    bool err = false;
    ofstream ostr( "particles_bfieldcyl.dat" );
    ostr << "# "
	 << setw(12) << "Time (s)" << " "
	 << setw(14) << "x (m)" << " "
	 << setw(14) << "r (m)" << "\n";
    ParticleCyl &prt = pdb.particle(0);
    for( uint32_t b = 0; b < prt.traj_size(); b++ ) {
	double t = prt.traj(b)(0);
	double x = prt.traj(b)(1);
	double r = prt.traj(b)(3);
	ostr << setw(14) << t << " "
	     << setw(14) << x << " "
	     << setw(14) << r << "\n";
	if( fabs( r - 0.020728544449 ) > 1e-5 ) {
	    std::cout << "Error: trajectory differs from theory\n";
	    exit( 1 );
	}
    }

    GeomPlotter geomplotter( geom );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.set_particle_div( 1 );
    geomplotter.plot_png( "particles_bfieldcyl.png" );

    ostr.close();
    if( err ) {
	std::cout << "Error: trajectory differs from theory\n";
	exit( 1 );
    }
}


