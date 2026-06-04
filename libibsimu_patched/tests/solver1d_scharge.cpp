/*! \file solver1d_scharge.cpp
 *  \brief Test solver with a 1d problem with space charge.
 *
 *  \test The simple 1d problems with constant space charge are easily
 *  solved analytically. The Poisson equation 
 *  \f[ \nabla^2 \phi = -\frac{\rho}{\epsilon} \f]
 *  can be integrated twice to get
 *  \f[ \phi = -\frac{\rho}{2 \epsilon} x^2 + Ax + B, \f]
 *  where A and B are integration constants set by boundary conditions.
 *  If we have set space charge \f$ \rho = 1\cdot10^{-4}\mathrm{~C/m}^3 \f$
 *  and Dirichlet boundaries \f$ \phi=0 \f$ at x=0 and x=0.1 m, the 
 *  constants become
 *  \f[ A=\frac{\rho}{2 \epsilon}\cdot 0.1\mathrm{~m}=564.705\mathrm{~V/m} \mathrm{~and~} B=0 \f]
 *
 */


#include <fstream>
#include <iomanip>
#include "epot_bicgstabsolver.hpp"
#include "meshscalarfield.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"


using namespace std;


const double eps0 = 8.85418781762e-12;
double phi( double x )
{
    return( -1.0e-4*x*x/(2.0*eps0) + 1.0e-4/(2.0*eps0)*0.1*x + 0 );
}


void test( int argc, char **argv )
{
    bool err = false;

    Geometry geom( MODE_1D, Int3D(11,1,1), Vec3D(0,0,0), 0.01 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 0.0) );
    geom.build_mesh();

    MeshScalarField epot( geom );
    MeshScalarField scharge( geom );
    for( uint32_t a = 0; a < geom.size(0); a++ )
	scharge(a) = 1.0e-4;

    EpotBiCGSTABSolver solver( geom );
    solver.solve( epot, scharge );

    ofstream ostr( "solver1d_scharge.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "potential (V)" << "\n";
    for( uint32_t a = 0; a < geom.size(0); a++ ) {
	if( fabs( epot(a) - phi(a*geom.h()) ) > 1.0e-4  )
	    err = true;
	ostr << setw(14) << a*geom.h() << " " 
	     << setw(14) << epot(a) << "\n";
    }
    ostr << "\n\n";
    ostr.close();

    if( err ) {
	std::cout << "Error: solved potential differs from theory\n";
	exit( 1 );
    }
}


