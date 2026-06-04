/*! \file particles_edip.cpp 
 *  \brief Test particle iterator in circular electrostatic dipole.
 *
 *  \test Test particle iterator in circular electrostatic dipole.
 */

/*
 *  Test particle iterator in circular electrostatic dipole.
 *  The trajectory is circular while
 *  \f[ F = ma = \frac{mv^2}{r}, \f]
 *  where \f$r\f$ is the radius of turn in the dipole. The force of 
 *  electrostatic field
 *  \f[ F = qE = q \frac{A}{r}, \f]
 *  where 
 *  \f[ A = \frac{\phi_a - \phi_b}{\log\frac{r_a}{r_b}}, \f]
 *  in a coaxial system where the potential distribution is
 *  \f[ \phi = A \log r + B. \f]
 *  It can be solved that
 *  \f[ A = \frac{2 U}{\log(\frac{r_a}{r_b})}, \f]
 *  where +U is the potential of the other electrode and -U is
 *  the potential of the other electrode on the dipole. Therefore
 *  \f[ U = \frac{m}{q} \frac{v^2}{2} \log(\frac{r_a}{r_b}). \f]
 *  For a system with \f$r_a = 6\mathrm{~cm}\f$ and \f$r_b = 4\mathrm{~cm}\f$,
 *  and particle of Neon-20, q=+6 accelerated with 10 kV extraction 
 *  voltage, the potential on the dipole is about 4055 V.
 *
 */


#include <fstream>
#include <iomanip>
#include "geometry.hpp"
#include "epot_bicgstabsolver.hpp"
#include "func_solid.hpp"
#include "meshvectorfield.hpp"
#include "epot_efield.hpp"
#include "particles.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "geomplotter.hpp"
#include "config.h"

#ifdef GTK3
#include "gtkplotter.hpp"
#endif


using namespace std;


bool inside_electrode( double x, double y, double z )
{
    return( x*x + y*y < 0.04*0.04 );
}

bool outside_electrode( double x, double y, double z )
{
    return( x*x + y*y > 0.06*0.06 );
}


void test( int argc, char **argv )
{
    Geometry geom( MODE_2D, Int3D(71,71,1), Vec3D(0.0,0.0,0.0), 0.001 );

    Solid *solid_in = new FuncSolid( inside_electrode );
    geom.set_solid( 7, solid_in );

    Solid *solid_out = new FuncSolid( outside_electrode );
    geom.set_solid( 8, solid_out );

    // Speed of particle
    double vy = sqrt(2*10e3*6.0*CHARGE_E/(19.9924*MASS_U));
    //std::cout << "Velocity = " << vy << "\n";

    // Deflection pot
    double U = 19.9924*MASS_U*vy*vy*log(6.0/4.0)/(6.0*CHARGE_E*2.0);
    //std::cout << "U = " << U << "\n";

    geom.set_boundary( 1, Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,   -U) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET ,  +U) );
    geom.build_mesh();
    //geom.debug_print();

    EpotField epot( geom );
    MeshScalarField scharge( geom );

    EpotBiCGSTABSolver solver( geom );
    solver.solve( epot, scharge );

    EpotEfield efield( epot );
    MeshVectorField bfield;

    ParticleDataBase2D pdb( geom );
    //pdb.set_accuracy( 1.0e-8, 1.0e-8 );
    pdb.set_thread_count( 1 );
    pdb.set_max_steps( 100 );

    // Add 20-Ne, q=6+
    pdb.add_particle( 0.0, 6.0, 19.9924, ParticleP2D( 0, 0.05, 0.0, 0.0, vy ) );
    pdb.iterate_trajectories( scharge, efield, bfield );
    //pdb.debug_print();

    /*
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
    */

    GeomPlotter geomplotter( geom );
    geomplotter.set_size( 1024, 768 );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.plot_png( "particles_edip.png" );

#ifdef GTK3
    if( false ) {
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_scharge( &scharge );
	plotter.set_epot( &epot );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
    }
#endif
}



