/*! \file particles_relativistic.cpp 
 *  \brief Test relativistic particle iterator in 2D.
 *
 *  \test  Test relativistic particle iterator in 2D.
 */


#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "epot_bicgstabsolver.hpp"
#include "meshvectorfield.hpp"
#include "epot_efield.hpp"
#include "func_solid.hpp"
#include "geomplotter.hpp"
#include "particles.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "ibsimutest.hpp"
#include "config.h"

#ifdef GTK3
#include "gtkplotter.hpp"
#endif


using namespace std;


// 2D: relativistic electron in magnetic field, test for radius
// B*r=p/q=gamma*m0*v/q
void test1( int argc, char **argv )
{
    ibsimu.message(1) << "------------ Test1 ------------\n";

    Geometry geom( MODE_2D, Int3D(11,11,1), Vec3D(-0.05,-0.05,0.0), 0.01 );
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

    double v = 1.5e8;
    double gamma = 1.0/sqrt(1.0-v*v/SPEED_C2);
    //std::cout << "gamma = " << gamma << "\n";
    double Q = -1.0;
    double q = Q*CHARGE_E;
    double M = 1.0/1822.88;
    double m = M*MASS_U;
    double B = -0.05;
    double r = gamma*m*v/(B*q);
    //std::cout << "r(relativistic) = " << r << "\n";
    //std::cout << "r(classic) = " << m*v/(B*q) << "\n";

    bool fout[3] = {false,false,true};
    MeshVectorField bfield( MODE_2D, fout, Int3D(2,2,1), Vec3D(-0.05,-0.05,0.0), 0.1 );
    bfield.set( 0, 0, -Vec3D(0,0,B) );
    bfield.set( 1, 0, -Vec3D(0,0,B) );
    bfield.set( 0, 1, -Vec3D(0,0,B) );
    bfield.set( 1, 1, -Vec3D(0,0,B) );

    ParticleDataBase2D pdb( geom );
    pdb.set_thread_count( 1 );
    pdb.set_relativistic( true );
    pdb.set_max_steps( 100 );
    pdb.add_particle( 0.0, Q, M, ParticleP2D( 0, r, 0.0, 0.0, v ) );
    pdb.iterate_trajectories( scharge, efield, bfield );
    //pdb.debug_print( std::cout );

    TrajectoryDiagnosticData tdata;
    std::vector<trajectory_diagnostic_e> diagnostics;
    diagnostics.push_back( DIAG_X );
    pdb.trajectories_at_plane( tdata, AXIS_Y, 0.0, diagnostics );
    int rcount = 0;
    double rave = 0.0;
    for( uint32_t a = 0; a < tdata.traj_size()-1; a+=2 ) {
	double rr = 0.5*fabs(tdata(a,0)-tdata(a+1,0));
	//std::cout << "y1 = " << tdata(a,0) << "\n";
	//std::cout << "y2 = " << tdata(a+1,0) << "\n";
	//std::cout << "rr = " << rr << "\n\n";
	rave += rr;
	rcount++;
    }
    rave = rave/rcount;
    //std::cout << "rave = " << rave << "\n";

    // Check particle trajectory points
    bool err = false;
    ofstream ostr( "particles_relativistic.dat" );
    ostr << "# "
	 << setw(12) << "Time (s)" << " "
	 << setw(14) << "x (m)" << " "
	 << setw(14) << "y (m)" << " "
	 << setw(14) << "r (m)" << "\n";
    Particle2D &prt = pdb.particle(0);
    for( uint32_t b = 0; b < prt.traj_size(); b++ ) {
	double t = prt.traj(b)(0);
	double x = prt.traj(b)(1);
	double y = prt.traj(b)(3);
	double rmeas = sqrt(x*x+y*y);
	ostr << setw(14) << t << " "
	     << setw(14) << x << " "
	     << setw(14) << y << " "
	     << setw(14) << rmeas << "\n";
	if( fabs( rmeas - r ) > 1e-5 )
	    err = true;
    }

    GeomPlotter geomplotter( geom );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.plot_png( "particles_relativistic.png" );

    ostr.close();
    if( err )
	throw( ErrorTest( ERROR_LOCATION, "trajectory differs from theory" ) );
}


bool electrode1_func( double x, double y, double z )
{
    double r = sqrt(y*y + z*z);
    return( x < 0.01 && r > 0.01 && 
            x < r-0.004 );
}


bool electrode2_func( double x, double y, double z )
{
    double r = sqrt(y*y + z*z);
    return( x > 0.09 && r > 0.01 &&
            0.1-x < r-0.004 );
}


// 3D: accelerate electrons to 500 keV, check final energy
void test2( int argc, char **argv )
{
    ibsimu.message(1) << "------------ Test2 ------------\n";

    double Q = -1.0;
    //double q = Q*CHARGE_E;
    double M = 1.0/1822.88;
    //double m = M*MASS_U;

    //Geometry geom( MODE_3D, Int3D(26,26,26), Vec3D(0,-0.05,-0.05), 0.004 );
    Geometry geom( MODE_3D, Int3D(51,51,51), Vec3D(0,-0.05,-0.05), 0.002 );

    Solid *solid1 = new FuncSolid( electrode1_func );
    geom.set_solid( 7, solid1 );
    Solid *solid2 = new FuncSolid( electrode2_func );
    geom.set_solid( 8, solid2 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,    3.0e3) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,  500.0e3) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,      0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,      0.0) );
    geom.set_boundary( 5, Bound(BOUND_NEUMANN,      0.0) );
    geom.set_boundary( 6, Bound(BOUND_NEUMANN,      0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,    3.0e3) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET,  500.0e3) );
    geom.build_mesh();

    EpotField epot( geom );
    MeshScalarField scharge( geom );

    EpotBiCGSTABSolver solver( geom );

    EpotEfield efield( epot );
    MeshVectorField bfield;
    ParticleDataBase3D pdb( geom );
    pdb.set_relativistic( true );
    pdb.set_thread_count( 4 );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );
    for( size_t a = 0; a < 4; a++ ) {

        solver.solve( epot, scharge );
	efield.recalculate();

        pdb.clear();
        pdb.add_cylindrical_beam_with_energy( 4000, 50.0e3, Q, M, 
                                              3.0e3, 0.0, 0.5,
                                              Vec3D(0,0,0), // center
                                              Vec3D(0,1,0), // dir1
                                              Vec3D(0,0,1), // dir2
                                              0.005 );      // radius
        pdb.iterate_trajectories( scharge, efield, bfield );
    }

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

    GeomPlotter gplotter( geom );
    gplotter.set_scharge( &scharge );
    gplotter.set_epot( &epot );
    gplotter.set_particle_database( &pdb );
    gplotter.set_particle_div( 11 );
    gplotter.set_font_size( 16 );
    gplotter.set_view( VIEW_XY );
    gplotter.plot_png( "particles_relativistic2_xy.png" );
    //gplotter.set_view( VIEW_YZ, 0 );
    //gplotter.plot_png( "particles_relativistic2_yz.png" );

    TrajectoryDiagnosticData tdata;
    std::vector<trajectory_diagnostic_e> diag;
    diag.push_back( DIAG_EK );
    pdb.trajectories_at_plane( tdata, AXIS_X, 0.1, diag );
    double Ek = 0.0;
    for( uint32_t a = 0; a < tdata.traj_size(); a++ )
	Ek += tdata(a,0);
    Ek = Ek/tdata.traj_size();
    if( fabs(Ek-500e3) > 50.0 ) {
	ofstream of( "particles_relativistic_out.txt" );
	for( uint32_t a = 0; a < tdata.traj_size(); a++ )
	    of << tdata(a,0) << "\n";
	of.close();
	throw( ErrorTest( ERROR_LOCATION, "incorrect particle energy" ) );
    }
}


void test( int argc, char **argv )
{
    test1( argc, argv );
    test2( argc, argv );
}
