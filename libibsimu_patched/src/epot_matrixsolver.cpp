/*! \file epot_matrixsolver.cpp
 *  \brief Matrix solver for electric potential problem
 */

/* Copyright (c) 2011,2012,2017,2021 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * taneli.kalvas@jyu.fi.
 * 
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
 */


#include "epot_matrixsolver.hpp"
#include "constants.hpp"
#include "ibsimu.hpp"


EpotMatrixSolver::Node2DoF::Node2DoF() 
  : _size(0), _n2d(0) 
{

}


EpotMatrixSolver::Node2DoF::Node2DoF( Int3D size ) 
  : _size(size) 
{
    _n2d = new uint32_t[_size[0]*_size[1]*_size[2]];
}


EpotMatrixSolver::Node2DoF::~Node2DoF() 
{
    if( _n2d )
	delete _n2d; 
}


void EpotMatrixSolver::Node2DoF::clear( void ) 
{
    if( _n2d )
	delete _n2d;
    _n2d = 0;
    _size = Int3D(0,0,0);
}


void EpotMatrixSolver::Node2DoF::resize( Int3D size ) 
{
    _size = size;
    if( _n2d )
	delete _n2d;
    _n2d = new uint32_t[_size[0]*_size[1]*_size[2]];
}


void EpotMatrixSolver::Node2DoF::debug_print( std::ostream &os ) const
{
    uint32_t a;
    uint32_t nc = _size[0]*_size[1]*_size[2];
    os << "**Node2DoF\n";
    os << "size = ("
        << _size[0] << ", "
        << _size[1] << ", "
        << _size[2] << ")\n";
    os << "n2d = (";
    if( _n2d ) {
	if( nc <= 10 ) {
	    for( a = 0; a < nc-1; a++ )
		os << _n2d[a] << ", ";
	    os << _n2d[a] << ")\n";
	} else {
	    for( a = 0; a < 10; a++ )
		os << _n2d[a] << ", ";
	    os << " ...)\n";
	}
    } else {
	os << ")\n";
    }
}


EpotMatrixSolver::EpotMatrixSolver( Geometry &geom )
    : EpotSolver(geom), _dof(0), _fd_mat(0), _fd_vec(0), _d_vec(0)
{

}


EpotMatrixSolver::EpotMatrixSolver( Geometry &geom, std::istream &s )
    : EpotSolver(geom,s)
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


EpotMatrixSolver::~EpotMatrixSolver() 
{
    if( _fd_mat )
        delete _fd_mat;
    if( _fd_vec )
        delete _fd_vec;
}


/* Set a link in matrix system.
 *
 * Makes the node a dependent on node b. If b is positive this means
 * that element (a,b) of matrix A is set to value val. If b is
 * negative, the node has a fixed potential and the dependence should
 * be added to the vector side of the system of equations on row
 * a. The value of potential at node b is stored in epot(-b).
 */
void EpotMatrixSolver::set_link( uint32_t a, uint32_t b, double val )
{
    //std::cout << "set_link( a = " << a << ", b = " << b << ", val = " << val << ")\n";

    if( (b & N2D_TYPE_MASK) == N2D_TYPE_FIXED ) {
	//std::cout << "epot = " << (*_epot)(b & N2D_INDEX_MASK) << "\n";
        (*_fd_vec)(a) += -val * (*_epot)(b & N2D_INDEX_MASK);
    } else {
	//std::cout << "matrix construct\n";
        (*_fd_mat).construct_add( a & N2D_INDEX_MASK, 
				  b & N2D_INDEX_MASK, val );
    }
}


void EpotMatrixSolver::add_vacuum_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;

    switch( _geom.geom_mode() ) {
    case MODE_1D:
        set_link( a, _n2d(i-1,j,k), 1.0 );
        set_link( a, _n2d(i,j,k), -2.0 );
        set_link( a, _n2d(i+1,j,k), 1.0 );
        break;
    case MODE_2D:
        set_link( a, _n2d(i,j-1,k), 1.0 );
        set_link( a, _n2d(i-1,j,k), 1.0 );
        set_link( a, _n2d(i,j,k), -4.0 );
        set_link( a, _n2d(i+1,j,k), 1.0 );
        set_link( a, _n2d(i,j+1,k), 1.0 );
        break;
    case MODE_CYL:
        set_link( a, _n2d(i,j-1,k), 1.0-0.5/j );
        set_link( a, _n2d(i-1,j,k), 1.0 );
        set_link( a, _n2d(i,j,k), -4.0 );
        set_link( a, _n2d(i+1,j,k), 1.0 );
        set_link( a, _n2d(i,j+1,k), 1.0+0.5/j );
        break;
    case MODE_3D:
        set_link( a, _n2d(i,j,k-1), 1.0 );
        set_link( a, _n2d(i,j-1,k), 1.0 );
        set_link( a, _n2d(i-1,j,k), 1.0 );
        set_link( a, _n2d(i,j,k), -6.0 );
        set_link( a, _n2d(i+1,j,k), 1.0 );
        set_link( a, _n2d(i,j+1,k), 1.0 );
        set_link( a, _n2d(i,j,k+1), 1.0 );
        break;
    }

    if( _plasma == PLASMA_PEXP || _plasma == PLASMA_NSIMP || _plasma == PLASMA_SHIELD ) {

	bool inplasma = true;
	if( _plasma_calc_func )
	    // Test if within plasma calculation region
	    inplasma = (*_plasma_calc_func)(x);

	if( !inplasma ) {
	    (*_d_vec)(a) = 0; // No plasma calculation
	} else if( _plasma == PLASMA_PEXP ) {
	    double p = (*_sol)(a);
	    double rhst, drhst;
	    pexp_newton( rhst, drhst, p );
	    (*_fd_vec)(a) += rhst;
	    (*_d_vec)(a) = drhst;
	} else if( _plasma == PLASMA_NSIMP ) {
	    double p = (*_sol)(a);
	    double rhst, drhst;
	    nsimp_newton( rhst, drhst, p );
	    (*_fd_vec)(a) += rhst;
	    (*_d_vec)(a) = drhst;
	} else if( _plasma == PLASMA_SHIELD ) {
	    double p = (*_sol)(a);
	    double rhst, drhst;
	    shield_newton( rhst, drhst, p );
	    double R = -(*_scharge)(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
	    (*_fd_vec)(a) += R*rhst;
	    (*_d_vec)(a) = R*drhst;
	}
    }

    // Right hand side
    if( _plasma != PLASMA_SHIELD )
	(*_fd_vec)(a) += -(*_scharge)(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
}


void EpotMatrixSolver::add_near_solid_node_1d( uint32_t i, const Vec3D &x )
{
    uint32_t a = _n2d(i) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i);

    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( _geom.mesh(i) & SMESH_NEAR_SOLID_INDEX_MASK );
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis
    if( bindex & EPOT_SOLVER_BXMIN ) {
	set_link( a, _n2d(i), -2.0/(beta*beta) );
	set_link( a, _n2d(i+1), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	set_link( a, _n2d(i), -2.0/(alpha*alpha) );
	set_link( a, _n2d(i-1), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value(x) / alpha;
    } else {
	set_link( a, _n2d(i), -2.0/(alpha*beta) );
	set_link( a, _n2d(i-1), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i+1), 2.0/((alpha+beta)*beta) );
    }
}


void EpotMatrixSolver::add_near_solid_node_2d( uint32_t i, uint32_t j, const Vec3D &x )
{
    uint32_t a = _n2d(i,j) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j);

    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( _geom.mesh(i,j) & SMESH_NEAR_SOLID_INDEX_MASK );
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    double cof = 0.0;

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis
    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i+1,j), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i-1,j), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value(x) / alpha;
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i-1,j), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i+1,j), 2.0/((alpha+beta)*beta) );
    }

    // Ymin direction
    alpha = 1.0;
    if( sflag & 0x04 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Ymax direction
    beta = 1.0;
    if( sflag & 0x08 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Y axis
    if( bindex & EPOT_SOLVER_BYMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i,j+1), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(3).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i,j-1), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(4).value(x) / alpha;
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i,j-1), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i,j+1), 2.0/((alpha+beta)*beta) );
    }

    // Middle node
    set_link( a, _n2d(i,j), -cof );
}


void EpotMatrixSolver::add_near_solid_node_cyl( uint32_t i, uint32_t j, const Vec3D &x )
{
    uint32_t a = _n2d(i,j) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j);

    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( _geom.mesh(i,j) & SMESH_NEAR_SOLID_INDEX_MASK );
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    double cof = 0.0;

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis
    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i+1,j), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i-1,j), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value(x) / alpha;
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i-1,j), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i+1,j), 2.0/((alpha+beta)*beta) );
    }

    // Ymin direction
    alpha = 1.0;
    if( sflag & 0x04 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Ymax direction
    beta = 1.0;
    if( sflag & 0x08 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Y axis
    if( bindex & EPOT_SOLVER_BYMIN ) {
	// On-axis
	cof += 4.0;
	set_link( a, _n2d(i,j+1), 4.0 );
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i,j-1), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += (1.0)/(2.0*alpha)*(2.0/alpha+1.0/j)*
	    2.0*_geom.h()*_geom.get_boundary(4).value(x) / alpha;
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i,j-1), 1.0/(alpha+beta)*(2.0/alpha-1.0/j) );
	set_link( a, _n2d(i,j+1), 1.0/(alpha+beta)*(2.0/beta+1.0/j) );
    }

    // Middle node
    set_link( a, _n2d(i,j), -cof );
}


void EpotMatrixSolver::add_near_solid_node_3d( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j,k);

    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( _geom.mesh(i,j,k) & SMESH_NEAR_SOLID_INDEX_MASK );
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    double cof = 0.0;

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis
    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i+1,j,k), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i-1,j,k), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value(x) / alpha;
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i-1,j,k), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i+1,j,k), 2.0/((alpha+beta)*beta) );
    }

    // Ymin direction
    alpha = 1.0;
    if( sflag & 0x04 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Ymax direction
    beta = 1.0;
    if( sflag & 0x08 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Y axis
    if( bindex & EPOT_SOLVER_BYMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i,j+1,k), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(3).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i,j-1,k), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(4).value(x) / alpha;
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i,j-1,k), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i,j+1,k), 2.0/((alpha+beta)*beta) );
    }

    // Zmin direction
    alpha = 1.0;
    if( sflag & 0x10 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Zmax direction
    beta = 1.0;
    if( sflag & 0x20 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Z axis
    if( bindex & EPOT_SOLVER_BZMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i,j,k+1), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(5).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BZMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i,j,k-1), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(6).value(x) / alpha;
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i,j,k-1), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i,j,k+1), 2.0/((alpha+beta)*beta) );
    }

    // Middle node
    set_link( a, _n2d(i,j,k), -cof );
}


void EpotMatrixSolver::add_near_solid_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;

    switch( _geom.geom_mode() ) {
    case MODE_1D:
	add_near_solid_node_1d(i,x);
	break;
    case MODE_2D:
	add_near_solid_node_2d(i,j,x);
	break;
    case MODE_CYL:
	add_near_solid_node_cyl(i,j,x);
	break;
    case MODE_3D:
	add_near_solid_node_3d(i,j,k,x);
	break;
    }

    if( _plasma == PLASMA_PEXP || _plasma == PLASMA_NSIMP || _plasma == PLASMA_SHIELD ) {

	bool inplasma = true;
	if( _plasma_calc_func )
	    // Test if within plasma calculation region
	    inplasma = (*_plasma_calc_func)(x);

	if( !inplasma ) {
	    (*_d_vec)(a) = 0; // No plasma calculation
	} else if( _plasma == PLASMA_PEXP ) {
	    double p = (*_sol)(a);
	    double rhst, drhst;
	    pexp_newton( rhst, drhst, p );
	    (*_fd_vec)(a) += rhst;
	    (*_d_vec)(a) = drhst;
	} else if( _plasma == PLASMA_NSIMP ) {
	    double p = (*_sol)(a);
	    double rhst, drhst;
	    nsimp_newton( rhst, drhst, p );
	    (*_fd_vec)(a) += rhst;
	    (*_d_vec)(a) = drhst;
	} else if( _plasma == PLASMA_SHIELD ) {
	    double p = (*_sol)(a);
	    double rhst, drhst;
	    shield_newton( rhst, drhst, p );
	    double R = -(*_scharge)(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
	    (*_fd_vec)(a) += R*rhst;
	    (*_d_vec)(a) = R*drhst;
	}
    }

    // Right hand side
    if( _plasma != PLASMA_SHIELD )
	(*_fd_vec)(a) += -(*_scharge)(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
}


void EpotMatrixSolver::add_neumann_node_1d( uint32_t i, const Vec3D &x )
{
    uint32_t a = _n2d(i) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i);

    if( bindex & EPOT_SOLVER_BXMIN ) {
	set_link( a, _n2d(i), -2.0 );
	set_link( a, _n2d(i+1), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value( x );
    } else {
	set_link( a, _n2d(i-1), 2.0 );
	set_link( a, _n2d(i), -2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value( x );
    }
}


void EpotMatrixSolver::add_neumann_node_2d( uint32_t i, uint32_t j, const Vec3D &x )
{
    uint32_t a = _n2d(i,j) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j);

    if( bindex & EPOT_SOLVER_BXMIN ) {
	set_link( a, _n2d(i+1,j), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value( x );
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	set_link( a, _n2d(i-1,j), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value( x );
    } else {
	set_link( a, _n2d(i-1,j), 1.0 );
	set_link( a, _n2d(i+1,j), 1.0 );
    }
    
    if( bindex & EPOT_SOLVER_BYMIN ) {
	set_link( a, _n2d(i,j+1), 2.0 );
	(*_fd_vec)(a) += -2.0*_geom.h()*_geom.get_boundary(3).value( x );
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	set_link( a, _n2d(i,j-1), 2.0 );
	(*_fd_vec)(a) += -2.0*_geom.h()*_geom.get_boundary(4).value( x );
    } else {
	set_link( a, _n2d(i,j-1), 1.0 );
	set_link( a, _n2d(i,j+1), 1.0 );
    }

    set_link( a, _n2d(i,j), -4.0 );
}


void EpotMatrixSolver::add_neumann_node_cyl( uint32_t i, uint32_t j, const Vec3D &x )
{
    uint32_t a = _n2d(i,j) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j);

    if( bindex & EPOT_SOLVER_BYMIN ) {
	// On-axis
	if( bindex & EPOT_SOLVER_BXMIN ) {
	    set_link( a, _n2d(i+1,j), 2.0 );
	    (*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value( x );
	} else if( bindex & EPOT_SOLVER_BXMAX ) {
	    set_link( a, _n2d(i-1,j), 2.0 );
	    (*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value( x );
	} else {
	    set_link( a, _n2d(i-1,j), 1.0 );
	    set_link( a, _n2d(i+1,j), 1.0 );
	}
	set_link( a, _n2d(i,j+1), 4.0 );
	set_link( a, _n2d(i,j), -6.0 );
    } else {
	// Off-axis
	if( bindex & EPOT_SOLVER_BXMIN ) {
	    set_link( a, _n2d(i+1,j), 2.0 );
	    (*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value( x );
	} else if( bindex & EPOT_SOLVER_BXMAX ) {
	    set_link( a, _n2d(i-1,j), 2.0 );
	    (*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value( x );
	} else {
	    set_link( a, _n2d(i-1,j), 1.0 );
	    set_link( a, _n2d(i+1,j), 1.0 );
	}
	
	if( bindex & EPOT_SOLVER_BYMAX ) {
	    set_link( a, _n2d(i,j-1), 2.0 );
	    (*_fd_vec)(a) += (1.0+0.5/j)*2.0*_geom.h()*_geom.get_boundary(2).value( x );
	} else {
	    set_link( a, _n2d(i,j-1), 1.0-0.5/j );
	    set_link( a, _n2d(i,j+1), 1.0+0.5/j );
	}
	
	set_link( a, _n2d(i,j), -4.0 );
    }
}


void EpotMatrixSolver::add_neumann_node_3d( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j,k);

    if( bindex & EPOT_SOLVER_BXMIN ) {
	set_link( a, _n2d(i+1,j,k), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value( x );
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	set_link( a,_n2d(i-1,j,k), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value( x );
    } else {
	set_link( a, _n2d(i-1,j,k), 1.0 );
	set_link( a, _n2d(i+1,j,k), 1.0 );
    }

    if( bindex & EPOT_SOLVER_BYMIN ) {
	set_link( a, _n2d(i,j+1,k), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(3).value( x );
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	set_link( a, _n2d(i,j-1,k), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(4).value( x );
    } else {
	set_link( a, _n2d(i,j-1,k), 1.0 );
	set_link( a, _n2d(i,j+1,k), 1.0 );
    }
    
    if( bindex & EPOT_SOLVER_BZMIN ) {
	set_link( a, _n2d(i,j,k+1), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(5).value( x );
    } else if( bindex & EPOT_SOLVER_BZMAX ) {
	set_link( a, _n2d(i,j,k-1), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(6).value( x );
    } else {
	set_link( a, _n2d(i,j,k-1), 1.0 );
	set_link( a, _n2d(i,j,k+1), 1.0 );
    }
    
    set_link( a, _n2d(i,j,k), -6.0 );
}


void EpotMatrixSolver::add_neumann_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;

    switch( _geom.geom_mode() ) {
    case MODE_1D:
	add_neumann_node_1d( i, x );
        break;
    case MODE_2D:
	add_neumann_node_2d( i, j, x );
        break;
    case MODE_CYL:
	add_neumann_node_cyl( i, j, x );
	break;
    case MODE_3D:
	add_neumann_node_3d( i, j, k, x );
        break;
    }

    if( _plasma == PLASMA_PEXP || _plasma == PLASMA_NSIMP || _plasma == PLASMA_SHIELD ) {

	bool inplasma = true;
	if( _plasma_calc_func )
	    // Test if within plasma calculation region
	    inplasma = (*_plasma_calc_func)(x);

	if( !inplasma ) {
	    (*_d_vec)(a) = 0; // No plasma calculation
	} else if( _plasma == PLASMA_PEXP ) {
	    double p = (*_sol)(a);
	    double rhst, drhst;
	    pexp_newton( rhst, drhst, p );
	    (*_fd_vec)(a) += rhst;
	    (*_d_vec)(a) = drhst;
	} else if( _plasma == PLASMA_NSIMP ) {
	    double p = (*_sol)(a);
	    double rhst, drhst;
	    nsimp_newton( rhst, drhst, p );
	    (*_fd_vec)(a) += rhst;
	    (*_d_vec)(a) = drhst;
	} else if( _plasma == PLASMA_SHIELD ) {
	    double p = (*_sol)(a);
	    double rhst, drhst;
	    shield_newton( rhst, drhst, p );
	    double R = -(*_scharge)(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
	    (*_fd_vec)(a) += R*rhst;
	    (*_d_vec)(a) = R*drhst;
	}
    }

    // Right hand side
    if( _plasma != PLASMA_SHIELD )
	(*_fd_vec)(a) += -(*_scharge)(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
}


void EpotMatrixSolver::reset_matrix( void )
{
    if( _fd_mat )
        delete _fd_mat;
    if( _fd_vec )
        delete _fd_vec;
    _fd_mat = 0;
    _fd_vec = 0;
    _dof = 0;

    _n2d.clear();
}


void EpotMatrixSolver::set_initial_guess( const MeshScalarField &epot, Vector &X ) const
{
    if( _dof == 0  )
	throw( Error( ERROR_LOCATION, "preprocess not done" ) );

    X.resize( _dof );
    for( int32_t a = 0; a < (int32_t)_geom.nodecount(); a++ ) {
	if( (_n2d(a) & N2D_TYPE_MASK) == N2D_TYPE_FREE )
	    X( _n2d(a) & N2D_INDEX_MASK ) = epot(a);
    }
}


void EpotMatrixSolver::set_solution( MeshScalarField &epot, const Vector &X ) const
{ 
    if( _dof == 0  )
	throw( Error( ERROR_LOCATION, "preprocess not done" ) );

    for( uint32_t a = 0; a < _geom.nodecount(); a++ ) {
	uint32_t b = _n2d(a);
        if( (b & N2D_TYPE_MASK) == N2D_TYPE_FREE )
            epot(a) = X(b);
    }
}


void EpotMatrixSolver::preprocess( MeshScalarField &epot, const MeshScalarField &scharge )
{
    _epot = &epot;
    _scharge = &scharge;
    EpotSolver::preprocess( epot );

    reset_matrix();

    // Build n2d array and calculate degrees of freedom.
    _n2d.resize( _geom.size() );
    _dof = 0;
    for( int32_t a = 0; a < (int32_t)_geom.nodecount(); a++ ) {

	uint32_t mesh = _geom.mesh(a);
	bool fixed = mesh & SMESH_NODE_FIXED;

	if( fixed )
	    // Fixed node
	    _n2d(a) = N2D_TYPE_FIXED | a;
	else { 
	    // Free node
	    _n2d(a) = N2D_TYPE_FREE | _dof;
	    _dof++;
	}
    }

    ibsimu.message( 1 ) << "dof = " << _dof << "\n";
    if( _dof == 0 )
        throw( Error( ERROR_LOCATION, "zero degrees of freedom" ) );
    
    // Allocate problem matrix and vector
    _fd_mat = new CRowMatrix( _dof, _dof );
    _fd_vec = new Vector( _dof );
}


void EpotMatrixSolver::build_mat_vec( void )
{
    if( _dof == 0  )
	throw( Error( ERROR_LOCATION, "preprocess not done" ) );
    _fd_mat->clear();
    _fd_vec->clear();

    // Build matrix and rhs vector
    Vec3D x;
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	x[2] = _geom.origo(2) + _geom.h()*k;
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    x[1] = _geom.origo(1) + _geom.h()*j;
            for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		x[0] = _geom.origo(0) + _geom.h()*i;

		uint32_t a = (k*_geom.size(1)+j)*_geom.size(0)+i;
		uint32_t mesh = _geom.mesh(a);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		bool fixed = mesh & SMESH_NODE_FIXED;

		if( fixed ) {
		    //ibsimu.message( 1 ) << "Fixed\n";
		    continue;
		} else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
		    //ibsimu.message( 1 ) << "Vacuum\n";
		    add_vacuum_node( i, j, k, x );
		} else if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		    //ibsimu.message( 1 ) << "Near solid\n";
		    add_near_solid_node( i, j, k, x );
		} else if( node_id == SMESH_NODE_ID_NEUMANN ) {
		    //ibsimu.message( 1 ) << "Neumann\n";
		    add_neumann_node( i, j, k, x );
		}
	    }
        }
    }

    // Order matrix
    _fd_mat->order_ascending();
}


void EpotMatrixSolver::postprocess( void )
{ 
    reset_matrix();
    EpotSolver::postprocess();
}


void EpotMatrixSolver::get_vecmat( const CRowMatrix **A, const Vector **B )
{
    build_mat_vec();
    *A = _fd_mat;
    *B = _fd_vec;
}


void EpotMatrixSolver::get_resjac( const CRowMatrix **J, const Vector **R, const Vector &X )
{
    // Build residual and jacobian from linear matric and vector.
    // Calculate R = J0*X - B(X) and J = J0 + I*D(X)
    _d_vec = new Vector( _dof );

    // Construct whole right hand side to _fd_vec, nonlinear component of diagonal 
    // to _d_vec and linear part of jacobian to _fd_mat.
    _sol = &X;
    build_mat_vec();

    // Linear part
    Vector w = (*_fd_mat) * X;

    for( uint32_t a = 0; a < _dof; a++ ) {

	(*_fd_vec)(a) = w(a) - (*_fd_vec)(a);
	_fd_mat->set(a,a) -= (*_d_vec)(a);
    }

    *J = _fd_mat;
    *R = _fd_vec;

    delete _d_vec;
}


/*
bool EpotMatrixSolver::linear( void ) const
{
    if( _plasma == PLASMA_NONE || _plasma == PLASMA_PEXP_INITIAL || _plasma == PLASMA_NSIMP_INITIAL )
        return( true );
    else
        return( false );
}
*/


void EpotMatrixSolver::debug_print( std::ostream &os ) const
{
    EpotSolver::debug_print( os );
    os << "**EpotMatrixSolver\n";
    os << "dof = " << _dof << "\n";
    _n2d.debug_print( os );
    if( _fd_mat )
	os << "fd_mat = \n" << *_fd_mat << "\n";
    else
	os << "fd_mat = NULL\n";
    if( _fd_vec )
	os << "fd_vec = \n" << *_fd_vec << "\n";
    else
	os << "fd_vec = NULL\n";
}


void EpotMatrixSolver::save( std::ostream &s ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}
