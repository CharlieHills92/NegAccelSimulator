/*! \file plasmacyl.cpp
 *  \brief Test with a plasma in cylindrical electrode configuration.
 *
 *  \test Test with a plasma in cylindrical electrode configuration.
 */


#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "epot_gssolver.hpp"
#include "epot_bicgstabsolver.hpp"
#include "epot_mgsolver.hpp"
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
    Geometry geom( MODE_CYL, Int3D(241,141,1), Vec3D(0,0,0), 5e-5 );

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
    
    EpotBiCGSTABSolver solver( geom );
    //EpotUMFPACKSolver solver( geom );
    //EpotMGSolver solver( geom );
    //solver.set_levels( 3 );
    //solver.set_mgcycmax( 1 );
    //EpotGSSolver solver( geom );
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

    ParticleDataBaseCyl pdb( geom );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );

    Emittance emit;

    Convergence conv;
    conv.add_epot( epot );
    conv.add_scharge( scharge );
    conv.add_emittance( 0, emit );

    for( size_t i = 0; i < 3; i++ ) {

	if( i == 1 ) {
	    double rhoe = pdb.get_rhosum();
	    //solver.set_w( 1.0 );
	    //solver.set_imax( 1000 );
	    solver.set_pexp_plasma( -rhoe, 5.0, 5.0 );
	}

	solver.solve( epot, scharge );
	efield.recalculate();

	pdb.clear();
	pdb.add_2d_beam_with_energy( 5000, 600.0, 1.0, 1.0, 
				     5.0, 0.0, 2.0, 
				     0.0, 0.0, 
				     0.0, 1.5e-3 );
	pdb.iterate_trajectories( scharge, efield, bfield );

	ParticleDiagPlotter pplotter( geom, pdb, AXIS_X, 0.0119, 
				      PARTICLE_DIAG_PLOT_SCATTER,
				      DIAG_R, DIAG_RP );
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

    ofstream ofconv( "plasmacyl_conv.dat" );
    conv.print_history( ofconv );
    ofconv.close();

    /*
    std::vector<trajectory_diagnostic_e> diagnostics;
    diagnostics.push_back( DIAG_R );
    diagnostics.push_back( DIAG_RP );
    diagnostics.push_back( DIAG_AP );
    diagnostics.push_back( DIAG_CURR );
    diagnostics.push_back( DIAG_NO );
    TrajectoryDiagnosticData tdata;
    pdb.trajectories_at_plane( tdata, AXIS_X, 1.0e-6, diagnostics );
    tdata.export_data( "plasmacyl_diag.txt" );
    */

    ParticleDiagPlotter pplotter1( geom, pdb, AXIS_X, 2e-6, 
				   PARTICLE_DIAG_PLOT_HISTO2D, 
				   DIAG_Y, DIAG_YP );
    pplotter1.set_font_size( 20 );
    pplotter1.set_size( 800, 600 );
    pplotter1.plot_png( "plasmacyl_emit1.png" );
    pplotter1.export_data( "plasmacyl_emit1.txt" );

    /*
    ParticleDiagPlotter pplotter2( geom, pdb, AXIS_X, 2.0e-3, 
				   PARTICLE_DIAG_PLOT_HISTO2D, 
				   DIAG_Y, DIAG_YP );
    //pplotter2.set_ranges( -0.008, -0.15001, 0.004, 0.05 );
    pplotter2.set_font_size( 20 );
    pplotter2.set_size( 800, 600 );
    pplotter2.plot_png( "plasmacyl_emit2.png" );

    ParticleDiagPlotter pplotter3( geom, pdb, AXIS_X, 6.0e-3, 
				   PARTICLE_DIAG_PLOT_HISTO2D, 
				   DIAG_Y, DIAG_YP );
    //pplotter3.set_ranges( -0.008, -0.15001, 0.004, 0.05 );
    pplotter3.set_font_size( 20 );
    pplotter3.set_size( 800, 600 );
    pplotter3.plot_png( "plasmacyl_emit3.png" );

    ParticleDiagPlotter pplotter4( geom, pdb, AXIS_X, 11.90e-3, 
				   PARTICLE_DIAG_PLOT_HISTO2D, 
				   DIAG_Y, DIAG_YP );
    //pplotter4.set_ranges( -0.008, -0.15001, 0.004, 0.05 );
    pplotter4.set_font_size( 20 );
    pplotter4.set_size( 800, 600 );
    pplotter4.plot_png( "plasmacyl_emit4.png" );

    */

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
    gplotter.plot_png( "plasmacyl.png" );
}

