/*! \file nsimp_plasma3d.cpp
 *  \brief Test with a simple negative ion plasma in 3d
 *  electrode configuration.
 *
 *  \test Test with a simple negative ion plasma in 3d
 *  electrode configuration.
 */


#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "epot_bicgstabsolver.hpp"
#include "epot_umfpacksolver.hpp"
#include "epot_gssolver.hpp"
#include "epot_mgsolver.hpp"
#include "particledatabase.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "particlediagplotter.hpp"
#include "geomplotter.hpp"
#include "config.h"

#ifdef GTK3
#include "gtkplotter.hpp"
#endif


using namespace std;


bool solid1( double x, double y, double z )
{
    double r = sqrt(y*y+z*z);

    if( x <= 0.0 && r <= 0.005 )
        return( false );

    return( x <= 0.004 && r >= 0.0030 && 
	    r >= 4.8*x - 0.003 );
}


bool solid2( double x, double y, double z )
{
    double r = sqrt(y*y+z*z);

    return( (x >= 0.025 || r <= 3.25*x - 0.06925) &&
	    x >= 0.0231 && r >= 0.0055 );
}


bool initial_plasma( double x, double y, double z )
{
    return( x <= 0.00055 );
}


void test( int argc, char **argv )
{
    Geometry geom( MODE_3D, Int3D(61,81,81), Vec3D(-0.001,-2.0e-2,-2.0e-2), 5.0e-4 );

    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0 ) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, +6.0e3) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 5, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 6, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0)  );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, +6.0e3) );
    geom.build_mesh();

    //EpotMGSolver solver( geom );
    //solver.set_levels( 3 );
    //EpotGSSolver solver( geom );
    //solver.set_w( 1.93 );
    //EpotUMFPACKSolver solver( geom );
    EpotBiCGSTABSolver solver( geom );
    //solver.set_newton_imax( 1 );
    InitialPlasma initp( AXIS_X, 0.0006 );
    solver.set_nsimp_initial_plasma( &initp );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;
    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase3D pdb( geom );
    bool pmirror[6] = { false, false, true, false, true, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );

    for( size_t i = 0; i < 8; i++ ) {

	if( i == 1 ) {
	    // 50%/50%
	    std::vector<double> Ei;
	    std::vector<double> rhoi;
	    Ei.push_back( 1.0 );
            double rhoneg = fabs(pdb.get_rhosum());
            rhoi.push_back( 0.5*rhoneg );
	    solver.set_nsimp_plasma( 0.5*rhoneg, 10.0, rhoi, Ei );

	    // 100% fast
	    /*
	    std::vector<double> Ei;
	    std::vector<double> rhoi;
            double rhoneg = fabs(pdb.get_rhosum());
	    solver.set_nsimp_plasma( rhoneg, 10.0, rhoi, Ei );
	    */

	    // 100% thermal
	    /*
	    std::vector<double> Ei;
	    std::vector<double> rhoi;
	    Ei.push_back( 1.0 );
            double rhoneg = fabs(pdb.get_rhosum());
            rhoi.push_back( rhoneg );
	    solver.set_nsimp_plasma( 0.0, 10.0, rhoi, Ei );
	    */
	}

	solver.solve( epot, scharge );
	efield.recalculate();

	pdb.clear();
	// H-
	// 1 mA total -> 35.37 A/m2
	pdb.add_cylindrical_beam_with_energy( 50000, -35.37, -1.0, 1.0, 
                                              5.0, 0.0, 0.5, 
                                              Vec3D(-1e-3,0,0),
                                              Vec3D(0,1,0),
                                              Vec3D(0,0,1), 0.005 );

	// e-
	// 20 mA total -> 707.4 A/m2
	pdb.add_cylindrical_beam_with_energy( 50000, -707.4, -1.0, 1.0/1836.15, 
                                              5.0, 0.0, 0.5, 
                                              Vec3D(-1e-3,0,0),
                                              Vec3D(0,1,0),
                                              Vec3D(0,0,1), 0.005 );
	pdb.iterate_trajectories( scharge, efield, bfield );
    }

#ifdef GTK3
    if( false ) {
	MeshScalarField tdens( geom );
	pdb.build_trajectory_density_field( tdens );
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.set_efield( &efield );
	plotter.set_trajdens( &tdens );
	plotter.set_scharge( &scharge );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
    }
#endif

    GeomPlotter gplotter( geom );
    gplotter.set_size( 1024, 768 );
    gplotter.set_epot( &epot );
    std::vector<double> eqlines;
    eqlines.push_back( +1.0 );
    eqlines.push_back( +2.0 );
    eqlines.push_back( +4.0 );
    eqlines.push_back( +8.0 );
    gplotter.set_eqlines_manual( eqlines );
    gplotter.set_particle_database( &pdb );
    gplotter.plot_png( "nsimp_plasma3d.png" );
}

