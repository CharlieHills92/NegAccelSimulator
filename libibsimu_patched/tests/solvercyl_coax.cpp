/*! \file solvercyl_coax.cpp 
 *  \brief Test solver with a cylindrical problem with coaxial
 *  electrodes.
 *
 *  \test Test solver with a cylindrical problem with coaxial
 *  electrodes.
 */


#include <fstream>
#include <iomanip>
#include "epot_bicgstabsolver.hpp"
#include "epot_mgsolver.hpp"
#include "epot_gssolver.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "config.h"

#ifdef GTK3
#include "gtkplotter.hpp"
#endif


using namespace std;


const double r1 = 0.021;
const double r2 = 0.07;
const double V1 = 0.0;
const double V2 = 10.0;
const double A = (V1-V2)/log(r1/r2);
const double B = V1-A*log(r1);


bool solid1( double x, double y, double z )
{
    return( y <= r1 );
}


bool solid2( double x, double y, double z )
{
    return( y >= r2 );
}


double phi( double r )
{
    return( A*log(r)+B );
}


void test( int argc, char **argv )
{
    //Geometry geom( MODE_CYL, Int3D(5,81,1), Vec3D(0,0,0), 0.001 );
    Geometry geom( MODE_CYL, Int3D(33,81,1), Vec3D(0,0,0), 0.001 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET, 10.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, 10.0) );
    geom.build_mesh();

    //EpotGSSolver solver( geom );
    //solver.set_imax( 100000 );
    //EpotBiCGSTABSolver solver( 
    EpotMGSolver solver( geom );
    solver.set_levels( 4 );
    //solver.set_mgcycmax( 1 );
    EpotField epot( geom );
    MeshScalarField scharge( geom );
    solver.solve( epot, scharge );

    double maxerr = 0.0;
    uint32_t loci = 0;
    uint32_t locj = 0;
    ofstream ostr( "solvercyl_coax.dat" );
    ostr << "# "
         << setw(12) << "x (m)" << " " 
         << setw(14) << "r (m)" << " " 
         << setw(14) << "potential (V)" << " "
         << setw(14) << "theory (V)" << "\n";
    for( uint32_t i = 0; i < geom.size(0); i++ ) {
        for( uint32_t j = 0; j < geom.size(1); j++ ) {
            double x = i*geom.h();
            double r = j*geom.h();
	    double err = fabs( epot(i,j) - phi(r) );
            if( r > r1 && r < r2 && err > maxerr ) {
		maxerr = err;
		loci = i;
		locj = j;
	    }
            ostr << setw(14) << x << " " 
                 << setw(14) << r << " " 
                 << setw(14) << epot(i,j) << " "
                 << setw(14) << phi(r) << "\n";
        }
	ostr << "\n\n";
    }
    ostr.close();

#ifdef GTK3
    if( false ) {
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_scharge( &scharge );
	plotter.set_epot( &epot );
	plotter.new_geometry_plot_window();
	plotter.run();
    }
#endif

    if( maxerr > 0.0022 ) {
	std::cout << "Maximum error = " << maxerr << "\n";
	std::cout << " at (" 
		  << loci << ", " 
		  << locj << ")\n";
	std::cout << "Error: solved potential differs from theory\n";
	exit( 1 );
    }
}
