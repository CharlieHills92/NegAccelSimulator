/*! \file plasma3d.cpp
 *  \brief Test with a plasma in 3d electrode configuration.
 *
 *  \test Test with a plasma in 3d electrode configuration.
 */


#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "epot_gssolver.hpp"
#include "epot_mgsolver.hpp"
#include "epot_bicgstabsolver.hpp"
#include "particledatabase.hpp"
#include "geometry.hpp"
#include "convergence.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
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
    if( x <= 0.0003 && r <= 0.001 )
	return( false );
    return( x <= 0.00187 && r >= 0.00054 && r >= 2.28*x - 0.0010 &&
	    (r >= 0.00054 || r >= 0.0015) );
}


bool solid2( double x, double y, double z )
{
    double r = sqrt(y*y+z*z);
    return( x >= 0.0095 && r >= 0.0023333 && r >= 0.01283 - x );
}


void test( int argc, char **argv )
{
    Geometry geom( MODE_3D, Int3D(121,73,73), Vec3D(0,0,0), 0.0001 );
    //Geometry geom( MODE_3D, Int3D(241,145,145), Vec3D(0,0,0), 0.00005 );
    //Geometry geom( MODE_3D, Int3D(481,289,289), Vec3D(0,0,0), 0.000025 );

    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0 ) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, -8.0e3) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 5, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 6, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0)  );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, -8.0e3) );
    geom.build_mesh();
    geom.build_surface();

    //EpotMGSolver solver( geom );
    //solver.set_levels( 4 );
    EpotBiCGSTABSolver solver( geom );
    InitialPlasma initp( AXIS_X, 0.0006 );
    solver.set_initial_plasma( 5.0, &initp );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;
    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
				     FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase3D pdb( geom );
    bool pmirror[6] = { false, false, true, false, true, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( false );
    pdb.set_surface_collision( true );

    Emittance emit;

    Convergence conv;
    conv.add_epot( epot );
    conv.add_scharge( scharge );
    conv.add_emittance( 0, emit );
	
    for( size_t i = 0; i < 3; i++ ) {

	if( i == 1 ) {
	    double rhoe = pdb.get_rhosum();
	    solver.set_pexp_plasma( -rhoe, 5.0, 5.0 );
	}

	solver.solve( epot, scharge );
	efield.recalculate();

	pdb.clear();
	pdb.add_cylindrical_beam_with_energy( 5000, 600.0, 1.0, 1.0, 
					      5.0, 0.0, 0.5, 
					      Vec3D(0,0,0),
					      Vec3D(0,1,0),
					      Vec3D(0,0,1), 0.001 );
	pdb.iterate_trajectories( scharge, efield, bfield );

	ParticleDiagPlotter pplotter( geom, pdb, AXIS_X, 0.0119, PARTICLE_DIAG_PLOT_SCATTER,
				      DIAG_Y, DIAG_YP );
	emit = pplotter.calculate_emittance();
	conv.evaluate_iteration();

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

    ofstream ofconv( "plasma3d_conv.dat" );
    conv.print_history( ofconv );
    ofconv.close();

    MeshScalarField tdens( geom );
    pdb.build_trajectory_density_field( tdens );

    GeomPlotter gplotter( geom );
    gplotter.set_size( 1024, 768 );
    gplotter.set_view( VIEW_XY, 0 );
    gplotter.set_epot( &epot );
    std::vector<double> eqlines;
    eqlines.push_back( -4.0 );
    eqlines.push_back( -2.0 );
    eqlines.push_back( 0.01 );
    eqlines.push_back( +2.0 );
    eqlines.push_back( +4.0 );
    gplotter.set_eqlines_manual( eqlines );
    gplotter.set_particle_database( &pdb );
    gplotter.set_particle_div( 0 );
    gplotter.set_trajdens( &tdens );
    gplotter.set_fieldgraph_plot( FIELD_TRAJDENS );
    gplotter.fieldgraph()->set_zscale( ZSCALE_RELLOG );
    gplotter.plot_png( "plasma3d.png" );
}

