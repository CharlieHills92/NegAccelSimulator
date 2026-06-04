/*! \file solver2d_coax.cpp 
 *  \brief Test solver with a 2d problem made of two concentric cylinders.
 *
 *  \test Test solver with a 2d problem made of two concentric cylinders.
 */

/*
 *  Analytically this problem is solved by using cylindrical coordinates.
 *  The Poisson equation 
 *  \f[ \nabla^2 \phi = -\frac{\rho}{\epsilon} \f]
 *  becomes
 *  \f[ \frac{\partial^2 \phi}{\partial r^2} + 
    \frac{1}{r} \frac{\partial \phi}{\partial r} +
    \frac{1}{r^2} \frac{\partial^2 \phi}{\partial \theta^2} +
    \frac{\partial^2 \phi}{\partial z^2} = -\frac{\rho}{\epsilon}. \f]
 *  The third and fourth terms are zero in this symmetric case and in case
 *  of no space charge, the differential equation becomes
 *  \f[ \frac{\partial^2 \phi}{\partial r^2} + 
    \frac{1}{r} \frac{\partial \phi}{\partial r} = 0. \f]
 *  By setting \f$ \frac{\partial \phi}{\partial r} = y \f$ and integrating
 *  the resulting separable first order DE we get  
 *  \f[ \frac{\partial \phi}{\partial r} = A \frac{1}{r}, \f]
 *  which can again be integrated to get the final solution
 *  \f[ \phi = A \ln r + B, \f]
 *  where A and B are integration constants set by boundary conditions.
 *  If we have Dirichlet boundaries \f$ \phi=0\mathrm{~V} \f$ at 
 *  \f$ x=0.02\mathrm{~m} \f$ and \f$ \phi=10\mathrm{~V} \f$ at 
 *  \f$ x=0.07\mathrm{~m} \f$, the constants become
 *  \f[ A=7.98235600148\mathrm{~V} \mathrm{~and~} B=31.2271603153 \f]
 *
 */


#include <fstream>
#include <iomanip>
#include "epot_bicgstabsolver.hpp"
#include "epot_gssolver.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_field.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "ibsimutest.hpp"
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
    return( x*x+y*y <= r1*r1 );
}


bool solid2( double x, double y, double z )
{
    return( x*x+y*y >= r2*r2 );
}


double phi( double r )
{
    if( r < r1 )
	return( V1 );
    else if( r > r2 )
	return( V2 );
    return( A*log(r)+B );
}


void test( int argc, char **argv )
{
    Geometry geom( MODE_2D, Int3D(41,41,1), Vec3D(0,0,0), 0.002 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 10.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET, 10.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, 10.0) );
    geom.build_mesh();

    EpotBiCGSTABSolver solver( geom );
    EpotField epot( geom );
    MeshScalarField scharge( geom );
    solver.solve( epot, scharge );

    bool err = false;
    ofstream ostr( "solver2d_coax.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "y (m)" << " " 
	 << setw(14) << "r (m)" << " " 
	 << setw(14) << "potential (V)" << " "
	 << setw(14) << "theory (V)" << "\n";
    for( uint32_t a = 0; a < geom.size(0); a++ ) {
	for( uint32_t b = 0; b < geom.size(1); b++ ) {
	    double x = a*geom.h();
	    double y = b*geom.h();
	    double r = sqrt(x*x + y*y);
	    if( r > r1 && r < r2 && fabs( epot(a,b) - phi(r) ) > 0.15  )
		err = true;
	    ostr << setw(14) << x << " " 
		 << setw(14) << y << " " 
		 << setw(14) << r << " " 
		 << setw(14) << epot(a,b) << " "
		 << setw(14) << phi(r) << "\n";
	}
    }

    ostr.close();

#ifdef GTK3
    if( false ) {
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.new_geometry_plot_window();
	plotter.run();
    }
#endif

    if( err )
	throw( ErrorTest( ERROR_LOCATION, "Error: solved potential differs from theory" ) );
}


