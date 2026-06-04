/*! \file solver3d_sphere.cpp 
 *  \brief Test solver with a 3d problem made of two concentric spheres.
 *
 *  \test Test solver with a 3d problem made of two concentric spheres.
 */


#include <fstream>
#include <iomanip>
#include "epot_gssolver.hpp"
#include "epot_mgsolver.hpp"
#include "epot_bicgstabsolver.hpp"
#include "epot_umfpacksolver.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "config.h"

#ifdef GTK3
#include "gtkplotter.hpp"
#endif


using namespace std;


double phi_a = 0.0;
double phi_b = 10.0;
double r_a = 0.02;
double r_b = 0.07;
double A = (phi_a - phi_b)/(1.0/r_a - 1.0/r_b);
double B = phi_a + A/r_a;


bool solid1( double x, double y, double z )
{
    return( x*x+y*y+z*z <= 0.02*0.02 );
}


bool solid2( double x, double y, double z )
{
    return( x*x+y*y+z*z >= 0.07*0.07 );
}


double phi( double r )
{
    if( r <= 0.02 )
	return( phi_a );
    else if( r >= 0.07 )
	return( phi_b );
    return( A/r - B );
}

/* Red-Black Gauss-Seidel
 *
 * w = 1.00, iter = 7347
 *
 *
 */

/* Optimize w at h=0.001
 *
 * w = 1.00, iter = 7380
 *
 * w = 1.85, iter = 755
 * w = 1.90, iter = 484
 * w = 1.91, iter = 424
 * w = 1.92, iter = 359
 * w = 1.93, iter = 279
 * w = 1.94, iter = 315
 * w = 1.95, iter = 375
 * w = 1.96, iter = 465
 * w = 1.99, iter = 1794
 *
 */

/* Speed test
 *
 * w = 1.93, -O2 compiler optimization
 *
 * fys80:
 * Intel(R) Core(TM)2 Quad CPU Q6700
 * 2666.67 MHz, 4096 KB
 *
 * golgata:
 * AMD Athlon(tm) 64 X2 Dual Core Processor 4800+
 * 2500.000 MHz, 512 KB cache
 *
 * h        iter    t/fys80       t/golgata
 * ---------------------------------------------------------
 * 0.002     260       0.26            0.48  
 * 0.001     279       2.01            3.44
 * 0.0005   1272      68.58          120.33
 * 0.00025  4542    1993.02
 *
 * Hand optimization: (i,j,k) -> (a,dj,dk)
 *
 * h        iter    t/fys80       t/golgata
 * ---------------------------------------------------------
 * 0.002     260          -            0.31
 * 0.001     279          -            2.41
 * 0.0005   1272          -           85.81
 * 0.00025                -
 *
 * Inline keyword: no effect
 * Compiler: -O3
 *
 * h        iter    t/fys80       t/golgata
 * ---------------------------------------------------------
 * 0.002     260          -            
 * 0.001     279          -            2.62
 * 0.0005   1272          -           85.41
 * 0.00025                -
 * 
 * Typed into one function (force inline)
 *
 * h        iter    t/fys80       t/golgata
 * ---------------------------------------------------------
 * 0.002     260          -            0.32
 * 0.001     279          -            2.51
 * 0.0005   1272          -           87.36
 * 0.00025                -
 * 
 * Back to separate functions and red-black sor algorithm, w = 1.93
 *
 * h        iter    t/fys80       t/golgata    max err
 * ---------------------------------------------------------
 * 0.002     258          -            0.48   0.293490
 * 0.001     258          -            3.12   0.138084
 * 0.0005   1200          -          112.24  0.0670977
 * 0.00025                -
 * 
 * Comparison same as previous, disabling red-black
 *
 * h        iter    t/fys80       t/golgata    max err
 * ---------------------------------------------------------
 * 0.002     260          -            0.30    0.29349
 * 0.001     279          -            2.39   0.138079
 * 0.0005   1272          -           83.66  0.0670906
 * 0.00025                -
 * 
 * Analysis. about 35 % overhead from algorithm change. Iteration
 * count does not increase -> suggests need for adjusting w.
 *
 */

void test_simu( int argc, char **argv )
{
    //double h = 0.0008;
    //int32_t size = (int32_t)ceil(0.08/h) + 1;
    int32_t size = 96+1;
    double h = 0.08/(size-1);
    Geometry g( MODE_3D, Int3D(size,size,size), Vec3D(0,0,0), h );
    Solid *s1 = new FuncSolid( solid1 );
    g.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    g.set_solid( 8, s2 );
    g.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 2, Bound(BOUND_DIRICHLET, 10.0) );
    g.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 4, Bound(BOUND_DIRICHLET, 10.0) );
    g.set_boundary( 5, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 6, Bound(BOUND_DIRICHLET, 10.0) );
    g.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 8, Bound(BOUND_DIRICHLET, 10.0) );
    g.build_mesh();
    
    EpotMGSolver solver( g );
    solver.set_levels( 4 );
    //EpotGSSolver solver( g );
    //EpotBiCGSTABSolver solver( g );
    //EpotUMFPACKSolver solver( g );
    //EpotRBGSSolver solver( g );
    EpotField epot( g );
    MeshScalarField scharge( g );

    //solver.set_imax( 1000000 );
    //solver.set_eps( 1e-6 );
    //solver.set_w( 1.93 );
    solver.solve( epot, scharge );

    /*
    ofstream osepot( "solver3d_sphere_h=0.0005.dat" );
    epot.save( osepot );
    osepot.close();
    return;



    ifstream isepot( "solver3d_sphere_h=0.0005.dat" );
    MeshScalarField epot_ref( isepot );
    isepot.close();

    ofstream dout( "solver3d_sphere_convergence_h=0.002.dat" );

    for( uint32_t iter = 0; iter < 100000; iter++ ) {
	// Take one iteration step
	solver.set_imax( 1 );
	solver.set_eps( 1e-20 );
	solver.set_w( 1.00 );
	solver.solve( epot, scharge );
	double res = solver.get_residual();

	// Calculate absolute error
	double abserr = 0.0;
	for( uint32_t a = 0; a < epot.nodecount(); a++ ) {
	    double err = fabs(epot(a)-epot_ref(a));
	    if( err > abserr )
		abserr = err;
	}

	// Store statistics
	dout << setw(12) << iter << " "
	     << setw(12) << res << " "
	     << setw(12) << res << "\n";

	if( res <= 1.0e-8 )
	    break;
    }

    dout.close();
    */

    bool err = false;
    double maxerr = 0.0;
    int maxerrl[3] = {0,0,0};
    ofstream ostr( "solver3d_sphere.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "y (m)" << " " 
	 << setw(14) << "z (m)" << " " 
	 << setw(14) << "r (m)" << " " 
	 << setw(14) << "potential (V)" << " "
	 << setw(14) << "theory (V)" << "\n";
    for( uint32_t a = 0; a < g.size(0); a++ ) {
	for( uint32_t b = 0; b < g.size(1); b++ ) {
	    for( uint32_t c = 0; c < g.size(2); c++ ) {
		double x = a*g.h();
		double y = b*g.h();
		double z = c*g.h();
		double r = sqrt(x*x + y*y + z*z);
		double e = epot(a,b,c);
		double error = fabs( epot(a,b,c) - phi(r) );
		if( r < 0.02 )
		    e = 0.0;
		else if( r > 0.07 )
		    e = 10.0;
		else if( r > 0.025 && error > 0.32 )
		    err = true;
		if( r > 0.02 && r < 0.07 && error > maxerr ) {
		    maxerr = error;
		    maxerrl[0] = a;
		    maxerrl[1] = b;
		    maxerrl[2] = c;
		}
		ostr << setw(14) << x << " " 
		     << setw(14) << y << " " 
		     << setw(14) << z << " " 
		     << setw(14) << r << " " 
		     << setw(14) << e << " "
		     << setw(14) << phi(r) << "\n";
	    }
	}
    }

    ostr.close();

    /*
    std::cout << "Maximum error = " << maxerr << "\n";
    std::cout << " at (" 
	      << maxerrl[0] << ", " 
	      << maxerrl[1] << ", " 
	      << maxerrl[2] << ")\n";
    */

#ifdef GTK3
    if( false ) {
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &g );
	plotter.set_scharge( &scharge );
	plotter.set_epot( &epot );
	plotter.new_geometry_plot_window();
	plotter.run();
    }
#endif

    if( err ) {
	std::cout << "Error: solved potential differs from theory\n";
	std::cout << "Maximum error = " << maxerr << "\n";
	std::cout << " at (" 
		  << maxerrl[0] << ", " 
		  << maxerrl[1] << ", " 
		  << maxerrl[2] << ")\n";
	exit( 1 );
    }
}


void test( int argc, char **argv )
{
    test_simu( argc, argv );

}
