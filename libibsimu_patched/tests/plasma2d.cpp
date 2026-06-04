/*! \file plasma2d.cpp
 *  \brief Test with a plasma in 2d electrode configuration.
 *
 *  \test Test with a plasma in 2d electrode configuration.
 */


#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "epot_mgsolver.hpp"
#include "epot_gssolver.hpp"
#include "epot_umfpacksolver.hpp"
#include "epot_bicgstabsolver.hpp"
#include "particledatabase.hpp"
#include "geometry.hpp"
#include "convergence.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "particlediagplotter.hpp"
#include "fielddiagplotter.hpp"
#include "geomplotter.hpp"
#include "config.h"

#ifdef GTK3
#include "gtkplotter.hpp"
#endif

using namespace std;


bool solid1( double x, double y, double z )
{
    return( x <= 2.0e-3 && y >= 0.5e-3 && y >= 2.0*x - 1.0e-3 &&
	    (x >= 0.5e-3 || y >= 1.5e-3) );
}


bool solid2( double x, double y, double z )
{
    return( x >= 10.0e-3 && y >= 1.5e-3 && y >= 12.0e-3 - x );
}


void test( int argc, char **argv )
{
    //Geometry geom( MODE_2D, Int3D(121,71,1), Vec3D(0,0,0), 1e-4 );
    //Geometry geom( MODE_2D, Int3D(241,141,1), Vec3D(0,0,0), 5e-5 );
    //Geometry geom( MODE_2D, Int3D(1201,705,1), Vec3D(0,0,0), 1e-5 );
    Geometry geom( MODE_2D, Int3D(76,45,1), Vec3D(0,0,0), 1.6e-4 );

    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0 ) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, -8.0e3) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0)  );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, -8.0e3) );
    geom.build_mesh();
    
    //EpotMGSolver solver( geom );
    //solver.set_mgcycmax( 10 );
    //solver.set_levels( 4 );
    //EpotUMFPACKSolver solver( geom );
    EpotGSSolver solver( geom );
    //EpotBiCGSTABSolver solver( geom );
    InitialPlasma initp( AXIS_X, 0.0006 );
    solver.set_initial_plasma( 5.0, &initp );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;
    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase2D pdb( geom );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );

    Convergence conv;
    conv.add_epot( epot );
    conv.add_scharge( scharge );
    Emittance emit;
    conv.add_emittance( 0, emit );

    for( size_t i = 0; i < 5; i++ ) {

	if( i == 1 ) {
	    double rhoe = pdb.get_rhosum();
	    //solver.set_w( 1.0 );
	    //solver.set_imax( 1000 );
	    solver.set_pexp_plasma( -rhoe, 5.0, 5.0 );
	}

	solver.solve( epot, scharge );
	efield.recalculate();

	pdb.clear();
	pdb.add_2d_beam_with_energy( 50000, 600.0, 1.0, 1.0, 
				     5.0, 0.0, 0.5, 
				     0.0, 0.0, 
				     0.0, 1.5e-3 );
	pdb.iterate_trajectories( scharge, efield, bfield );

	ParticleDiagPlotter pp( geom, pdb, AXIS_X, 11.90e-3, 
				PARTICLE_DIAG_PLOT_SCATTER, 
				DIAG_Y, DIAG_YP );
	emit = pp.calculate_emittance();

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

    ofstream ofconv( "plasma2d_conv.dat" );
    conv.print_history( ofconv );
    ofconv.close();

    ParticleDiagPlotter pplotter1( geom, pdb, AXIS_X, 1e-6, 
				   PARTICLE_DIAG_PLOT_HISTO2D, 
				   DIAG_Y, DIAG_YP );
    //pplotter1.set_ranges( -0.008, -0.15001, 0.004, 0.05 );
    pplotter1.set_font_size( 20 );
    pplotter1.set_size( 800, 600 );
    pplotter1.plot_png( "plasma2d_emit1.png" );

    ParticleDiagPlotter pplotter2( geom, pdb, AXIS_X, 2.0e-3, 
				   PARTICLE_DIAG_PLOT_HISTO2D, 
				   DIAG_Y, DIAG_YP );
    //pplotter2.set_ranges( -0.008, -0.15001, 0.004, 0.05 );
    pplotter2.set_font_size( 20 );
    pplotter2.set_size( 800, 600 );
    pplotter2.plot_png( "plasma2d_emit2.png" );

    ParticleDiagPlotter pplotter3( geom, pdb, AXIS_X, 6.0e-3, 
				   PARTICLE_DIAG_PLOT_HISTO2D, 
				   DIAG_Y, DIAG_YP );
    //pplotter3.set_ranges( -0.008, -0.15001, 0.004, 0.05 );
    pplotter3.set_font_size( 20 );
    pplotter3.set_size( 800, 600 );
    pplotter3.plot_png( "plasma2d_emit3.png" );

    ParticleDiagPlotter pplotter4( geom, pdb, AXIS_X, 11.90e-3, 
				   PARTICLE_DIAG_PLOT_HISTO2D, 
				   DIAG_Y, DIAG_YP );
    //pplotter4.set_ranges( -0.008, -0.15001, 0.004, 0.05 );
    pplotter4.set_font_size( 20 );
    pplotter4.set_size( 800, 600 );
    pplotter4.plot_png( "plasma2d_emit4.png" );

    MeshScalarField tdens( geom );
    pdb.build_trajectory_density_field( tdens );
    GeomPlotter gplotter( geom );
    gplotter.set_size( 800, 600 );
    gplotter.set_font_size( 20 );

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
    gplotter.plot_png( "plasma2d.png" );

    //gplotter.set_ranges( 0, 0, 0.002, 0.002 );
    //gplotter.plot_png( "plasma2d_zoom.png" );
    /*
    ParticleDiagPlotter pplotter( geom, pdb, AXIS_X, 0.01, 
				  PARTICLE_DIAG_PLOT_HISTO2D, DIAG_Y, DIAG_YP );
    pplotter.plot_png( "plasma2d_emit.png" );

    FieldDiagPlotter fplotter( geom );
    fplotter.set_scharge( &scharge );
    fplotter.set_epot( &epot );
    fplotter.set_coordinates( 100, Vec3D(0.006,0,0), Vec3D(0.006,0.007,0) );
    field_diag_type_e diag[2] = {FIELDD_DIAG_EPOT, FIELDD_DIAG_SCHARGE};
    field_loc_type_e loc[2] = {FIELDD_LOC_Y, FIELDD_LOC_NONE};
    fplotter.set_diagnostic( diag, loc );
    fplotter.plot_png( "plasma2d_field.png" );
    */
}

