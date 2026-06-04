/*! \file epot_efield.cpp
 *  \brief Electric potential base electric field.
 */

/* Copyright (c) 2005-2011,2013,2014 Taneli Kalvas. All rights reserved.
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

#include <limits>
#include <cmath>
#include "epot_efield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"


EpotEfield::EpotEfield( const EpotField &epot )
    : _epot(epot), _geom(epot.geom())
{
    _F[0] = NULL;
    _F[1] = NULL;
    _F[2] = NULL;

    _extrpl[0] = FIELD_EXTRAPOLATE;
    _extrpl[1] = FIELD_EXTRAPOLATE;
    _extrpl[2] = FIELD_EXTRAPOLATE;
    _extrpl[3] = FIELD_EXTRAPOLATE;
    _extrpl[4] = FIELD_EXTRAPOLATE;
    _extrpl[5] = FIELD_EXTRAPOLATE;

    precalc();
}    


void EpotEfield::copy_1d( const EpotEfield &efield )
{
    uint32_t n = _epot.size(0)+1;

    _F[0] = new double[n];
    _F[1] = NULL;
    _F[2] = NULL;

    memcpy( _F[0], efield._F[0], n*sizeof(double) );
}


void EpotEfield::copy_2d( const EpotEfield &efield )
{
    uint32_t n = _epot.size(0)+1;
    uint32_t m = _epot.size(1)+1;

    _F[0] = new double[n*(m-1)];
    _F[1] = new double[(n-1)*m];
    _F[2] = NULL;

    memcpy( _F[0], efield._F[0], n*(m-1)*sizeof(double) );
    memcpy( _F[1], efield._F[1], (n-1)*m*sizeof(double) );
}


void EpotEfield::copy_3d( const EpotEfield &efield )
{
    uint32_t n = _epot.size(0)+1;
    uint32_t m = _epot.size(1)+1;
    uint32_t o = _epot.size(2)+1;

    _F[0] = new double[n*(m-1)*(o-1)];
    _F[1] = new double[(n-1)*m*(o-1)];
    _F[2] = new double[(n-1)*(m-1)*o];

    memcpy( _F[0], efield._F[0], n*(m-1)*(o-1)*sizeof(double) );
    memcpy( _F[1], efield._F[1], (n-1)*m*(o-1)*sizeof(double) );
    memcpy( _F[2], efield._F[2], (n-1)*(m-1)*o*sizeof(double) );
}


EpotEfield::EpotEfield( const EpotEfield &efield )
    : _epot(efield._epot), _geom(efield._geom)
{
    memcpy( _extrpl, efield._extrpl, 6*sizeof(field_extrpl_e) );

    switch( _epot.geom_mode() ) {
    case MODE_1D:
	copy_1d( efield );
	break;
    case MODE_2D:
    case MODE_CYL:
	copy_2d( efield );
	break;
    case MODE_3D:
	copy_3d( efield );
	break;
    }
}


EpotEfield::~EpotEfield()
{
    for( size_t i = 0; i < 3; i++ ) {
	if( _F[i] != NULL )
	    delete [] _F[i];
    }
}


void EpotEfield::set_extrapolation( field_extrpl_e extrpl[6] ) 
{
    memcpy( _extrpl, extrpl, 6*sizeof(field_extrpl_e) );
    precalc();
}


void EpotEfield::recalculate( void )
{
    precalc();
}


uint8_t EpotEfield::solid_dist( uint32_t node, uint32_t dir ) const
{
    const uint8_t *nptr = _geom->nearsolid_ptr(node & SMESH_NEAR_SOLID_INDEX_MASK);
    uint8_t neighbours = nptr[0];
    nptr++;
    uint32_t a = 0;
    while( a < dir ) {
	if( neighbours & 0x01 )
	    nptr++;
	neighbours = neighbours >> 1;
	a++;
    }
    if( (neighbours & 0x01) == 0x00 )
	throw( Error( ERROR_LOCATION, (const std::string)"no near neighbour in selected direction"
		      ", dir = " + to_string(dir) 
		      + ", neighbours = " + 
		      to_string((int)*_geom->nearsolid_ptr(node & SMESH_NEAR_SOLID_INDEX_MASK)) ) );

    return( *nptr );
}


void EpotEfield::precalc_1d( void )
{
    double h = _epot.h();
    uint32_t n = _epot.size(0)+1;

    _F[0] = new double[n];

    for( uint32_t i = 1; i < n-1; i++ ) {
	
	uint32_t node1 = _geom->mesh(i-1);
	uint32_t node1id = node1 & SMESH_NODE_ID_MASK;
	uint32_t node2 = _geom->mesh(i);
	uint32_t node2id = node2 & SMESH_NODE_ID_MASK;

	if( node1id == SMESH_NODE_ID_NEAR_SOLID &&
	    node2id == SMESH_NODE_ID_DIRICHLET && 
	    (node2 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {

	    // Between near solid node and solid node
	    uint8_t dist = solid_dist( node1, 1 );
	    _F[0][i] = 255.0*( _epot(i-1) - _epot(i) ) / (dist*h);

	} else if( node2id == SMESH_NODE_ID_NEAR_SOLID &&
		   node1id == SMESH_NODE_ID_DIRICHLET && 
		   (node1 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {

	    // Between near solid node and solid node
	    uint8_t dist = solid_dist( node2, 0 );
	    _F[0][i] = 255.0*( _epot(i-1) - _epot(i) ) / (dist*h);

	} else if( node1id == SMESH_NODE_ID_DIRICHLET && 
		   (node1 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 &&
		   node2id == SMESH_NODE_ID_DIRICHLET && 
		   (node2 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
	    
	    // Inside solid, search for near solid in x-direction
	    uint32_t node0 = _geom->mesh_check(i-2);
	    uint32_t node0id = node0 & SMESH_NODE_ID_MASK;
	    uint32_t node3 = _geom->mesh_check(i+1);
	    uint32_t node3id = node3 & SMESH_NODE_ID_MASK;
	    if( node0id == SMESH_NODE_ID_NEAR_SOLID &&
		node3id == SMESH_NODE_ID_NEAR_SOLID ) {
		
		// No surface nearby
		_F[0][i] = 0.0;
		
	    } else if( node0id == SMESH_NODE_ID_NEAR_SOLID ) {
		
		// Between near solid node and solid node
		uint8_t dist = solid_dist( node0, 1 );
		_F[0][i] = 255.0*( _epot(i-2) - _epot(i-1) ) / (dist*h);
		
	    } else if( node3id == SMESH_NODE_ID_NEAR_SOLID ) {
		
		// Between near solid node and solid node
		uint8_t dist = solid_dist( node3, 0 );
		_F[0][i] = 255.0*( _epot(i) - _epot(i+1) ) / (dist*h);
		
	    } else {
		
		// Ambiguous case, surface in both directions
		_F[0][i] = 0.0;
	    }
	    
	} else {
	    
	    // Free space
	    _F[0][i] = (_epot(i-1) - _epot(i)) / h;

	}
    }

    // Extrapolate boundary from last two nodes.
    // If SYMMETRIC_POTENTIAL extrapolation, force E-field interpolation to
    // zero when approaching boundary, otherwise extrapolate from last
    // two nodes.
    if( _extrpl[0] == FIELD_SYMMETRIC_POTENTIAL )
	_F[0][0] = -_F[0][1];
    else
	_F[0][0] = 2.0*_F[0][1]-_F[0][2];

    if( _extrpl[1] == FIELD_SYMMETRIC_POTENTIAL )
	_F[0][n-1] = -_F[0][n-2];
    else
	_F[0][n-1] = 2.0*_F[0][n-2]-_F[0][n-3];
}


void EpotEfield::precalc_2d( void )
{
    double h = _epot.h();
    uint32_t n = _epot.size(0)+1;
    uint32_t m = _epot.size(1)+1;

    _F[0] = new double[n*(m-1)];
    _F[1] = new double[(n-1)*m];

    // Do Ex-field
    for( uint32_t j = 0; j < m-1; j++ ) {
	for( uint32_t i = 1; i < n-1; i++ ) {
	
	    uint32_t node1 = _geom->mesh(i-1,j);
	    uint32_t node1id = node1 & SMESH_NODE_ID_MASK;
	    uint32_t node2 = _geom->mesh(i,j);
	    uint32_t node2id = node2 & SMESH_NODE_ID_MASK;
		
	    if( node1id == SMESH_NODE_ID_NEAR_SOLID &&
		node2id == SMESH_NODE_ID_DIRICHLET && 
		(node2 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		
		// Between near solid node and solid node
		uint8_t dist = solid_dist( node1, 1 );
		_F[0][i+j*n] = 255.0*( _epot(i-1,j) - _epot(i,j) ) / (dist*h);
		
	    } else if( node2id == SMESH_NODE_ID_NEAR_SOLID &&
		       node1id == SMESH_NODE_ID_DIRICHLET && 
		       (node1 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		
		// Between near solid node and solid node
		uint8_t dist = solid_dist( node2, 0 );
		_F[0][i+j*n] = 255.0*( _epot(i-1,j) - _epot(i,j) ) / (dist*h);
		
	    } else if( node1id == SMESH_NODE_ID_DIRICHLET && 
		       (node1 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 &&
		       node2id == SMESH_NODE_ID_DIRICHLET && 
		       (node2 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		
		// Inside solid, search for near solid in x-direction
		uint32_t node0 = _geom->mesh_check(i-2,j);
		uint32_t node0id = node0 & SMESH_NODE_ID_MASK;
		uint32_t node3 = _geom->mesh_check(i+1,j);
		uint32_t node3id = node3 & SMESH_NODE_ID_MASK;
		if( node0id == SMESH_NODE_ID_NEAR_SOLID &&
		    node3id == SMESH_NODE_ID_NEAR_SOLID ) {

		    // No surface nearby
		    _F[0][i+j*n] = 0.0;

		} else if( node0id == SMESH_NODE_ID_NEAR_SOLID ) {
		    
		    // Between near solid node and solid node
		    uint8_t dist = solid_dist( node0, 1 );
		    _F[0][i+j*n] = 255.0*( _epot(i-2,j) - _epot(i-1,j) ) / (dist*h);

		} else if( node3id == SMESH_NODE_ID_NEAR_SOLID ) {
		    
		    // Between near solid node and solid node
		    uint8_t dist = solid_dist( node3, 0 );
		    _F[0][i+j*n] = 255.0*( _epot(i,j) - _epot(i+1,j) ) / (dist*h);

		} else {
		    
		    // Ambiguous case, surface in both directions
		    _F[0][i+j*n] = 0.0;
		}

	    } else {
		
		// Free space
		_F[0][i+j*n] = (_epot(i-1,j) - _epot(i,j)) / h;
	    }
	}
	
	// If SYMMETRIC_POTENTIAL extrapolation, force E-field interpolation to
	// zero when approaching boundary, otherwise extrapolate from last
	// two nodes.
	if( _extrpl[0] == FIELD_SYMMETRIC_POTENTIAL )
	    _F[0][0+j*n] = -_F[0][1+j*n];
	else
	    _F[0][0+j*n] = 2.0*_F[0][1+j*n]-_F[0][2+j*n];
	
	if( _extrpl[1] == FIELD_SYMMETRIC_POTENTIAL )
	    _F[0][n-1+j*n] = -_F[0][n-2+j*n];
	else
	    _F[0][n-1+j*n] = 2.0*_F[0][n-2+j*n]-_F[0][n-3+j*n];
    }

    // Do Ey-field
    for( uint32_t i = 0; i < n-1; i++ ) {
	for( uint32_t j = 1; j < m-1; j++ ) {
	
	    uint32_t node1 = _geom->mesh(i,j-1);
	    uint32_t node1id = node1 & SMESH_NODE_ID_MASK;
	    uint32_t node2 = _geom->mesh(i,j);
	    uint32_t node2id = node2 & SMESH_NODE_ID_MASK;

	    if( node1id == SMESH_NODE_ID_NEAR_SOLID &&
		node2id == SMESH_NODE_ID_DIRICHLET && 
		(node2 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		
		// Between near solid node and solid node
		uint8_t dist = solid_dist( node1, 3 );
		_F[1][i+j*_epot.size(0)] = 255.0*( _epot(i,j-1) - _epot(i,j) ) / (dist*h);

	    } else if( node2id == SMESH_NODE_ID_NEAR_SOLID &&
		       node1id == SMESH_NODE_ID_DIRICHLET && 
		       (node1 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		
		// Between near solid node and solid node
		uint8_t dist = solid_dist( node2, 2 );
		_F[1][i+j*_epot.size(0)] = 255.0*( _epot(i,j-1) - _epot(i,j) ) / (dist*h);
		
	    } else if( node1id == SMESH_NODE_ID_DIRICHLET && 
		       (node1 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 &&
		       node2id == SMESH_NODE_ID_DIRICHLET && 
		       (node2 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		
		// Inside solid, search for near solid in y-direction
		uint32_t node0 = _geom->mesh_check(i,j-2);
		uint32_t node0id = node0 & SMESH_NODE_ID_MASK;
		uint32_t node3 = _geom->mesh_check(i,j+1);
		uint32_t node3id = node3 & SMESH_NODE_ID_MASK;
		if( node0id == SMESH_NODE_ID_NEAR_SOLID &&
		    node3id == SMESH_NODE_ID_NEAR_SOLID ) {

		    // No surface nearby
		    _F[1][i+j*_epot.size(0)] = 0.0;

		} else if( node0id == SMESH_NODE_ID_NEAR_SOLID ) {
		    
		    // Between near solid node and solid node
		    uint8_t dist = solid_dist( node0, 3 );
		    _F[1][i+j*_epot.size(0)] = 255.0*( _epot(i,j-2) - _epot(i,j-1) ) / (dist*h);

		} else if( node3id == SMESH_NODE_ID_NEAR_SOLID ) {
		    
		    // Between near solid node and solid node
		    uint8_t dist = solid_dist( node3, 2 );
		    _F[1][i+j*_epot.size(0)] = 255.0*( _epot(i,j) - _epot(i,j+1) ) / (dist*h);

		} else {
		    
		    // Ambiguous case, surface in both directions
		    _F[1][i+j*_epot.size(0)] = 0.0;
		}

	    } else {
		
		// Free space
		_F[1][i+j*_epot.size(0)] = (_epot(i,j-1) - _epot(i,j)) / h;
	    }
	}

	// If SYMMETRIC_POTENTIAL extrapolation, force E-field interpolation to
	// zero when approaching boundary, otherwise extrapolate from last
	// two nodes.
	if( _extrpl[2] == FIELD_SYMMETRIC_POTENTIAL )
	    _F[1][i] = -_F[1][i+_epot.size(0)];
	else
	    _F[1][i] = 2.0*_F[1][i+_epot.size(0)]-_F[1][i+2*_epot.size(0)];
	
	if( _extrpl[3] == FIELD_SYMMETRIC_POTENTIAL )
	    _F[1][i+(m-1)*_epot.size(0)] = -_F[1][i+(m-2)*_epot.size(0)];
	else
	    _F[1][i+(m-1)*_epot.size(0)] = 2.0*_F[1][i+(m-2)*_epot.size(0)]-_F[1][i+(m-3)*_epot.size(0)];
    }
}


void EpotEfield::precalc_3d( void )
{
    double h = _epot.h();
    uint32_t n = _epot.size(0)+1;
    uint32_t m = _epot.size(1)+1;
    uint32_t o = _epot.size(2)+1;

    _F[0] = new double[n*(m-1)*(o-1)];
    _F[1] = new double[(n-1)*m*(o-1)];
    _F[2] = new double[(n-1)*(m-1)*o];

    // Do Ex-field
    for( uint32_t k = 0; k < o-1; k++ ) {
	for( uint32_t j = 0; j < m-1; j++ ) {
	    for( uint32_t i = 1; i < n-1; i++ ) {
	
		uint32_t node1 = _geom->mesh(i-1,j,k);
		uint32_t node1id = node1 & SMESH_NODE_ID_MASK;
		uint32_t node2 = _geom->mesh(i,j,k);
		uint32_t node2id = node2 & SMESH_NODE_ID_MASK;

		if( node1id == SMESH_NODE_ID_NEAR_SOLID &&
		    node2id == SMESH_NODE_ID_DIRICHLET && 
		    (node2 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		    
		    // Between near solid node and solid node
		    uint8_t dist = solid_dist( node1, 1 );
		    _F[0][i+(j+k*_epot.size(1))*n] = 
			255.0*( _epot(i-1,j,k) - _epot(i,j,k) ) / (dist*h);
		    
		} else if( node2id == SMESH_NODE_ID_NEAR_SOLID &&
			   node1id == SMESH_NODE_ID_DIRICHLET && 
			   (node1 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		    
		    // Between near solid node and solid node
		    uint8_t dist = solid_dist( node2, 0 );
		    _F[0][i+(j+k*_epot.size(1))*n] = 
			255.0*( _epot(i-1,j,k) - _epot(i,j,k) ) / (dist*h);
		    
		}  else if( node1id == SMESH_NODE_ID_DIRICHLET && 
			    (node1 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 &&
			    node2id == SMESH_NODE_ID_DIRICHLET && 
			    (node2 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		
		    // Inside solid, search for near solid in x-direction
		    uint32_t node0 = _geom->mesh_check(i-2,j,k);
		    uint32_t node0id = node0 & SMESH_NODE_ID_MASK;
		    uint32_t node3 = _geom->mesh_check(i+1,j,k);
		    uint32_t node3id = node3 & SMESH_NODE_ID_MASK;
		    if( node0id == SMESH_NODE_ID_NEAR_SOLID &&
			node3id == SMESH_NODE_ID_NEAR_SOLID ) {
			
			// No surface nearby
			_F[0][i+(j+k*_epot.size(1))*n] = 0.0;
			
		    } else if( node0id == SMESH_NODE_ID_NEAR_SOLID ) {
			
			// Between near solid node and solid node
			uint8_t dist = solid_dist( node0, 1 );
			_F[0][i+(j+k*_epot.size(1))*n] = 255.0*( _epot(i-2,j,k) - _epot(i-1,j,k) ) / (dist*h);
			
		    } else if( node3id == SMESH_NODE_ID_NEAR_SOLID ) {
			
			// Between near solid node and solid node
			uint8_t dist = solid_dist( node3, 0 );
			_F[0][i+(j+k*_epot.size(1))*n] = 255.0*( _epot(i,j,k) - _epot(i+1,j,k) ) / (dist*h);
			
		    } else {
			
			// Ambiguous case, surface in both directions
			_F[0][i+(j+k*_epot.size(1))*n] = 0.0;
		    }
		    
		} else {
		    
		    // Free space
		    _F[0][i+(j+k*_epot.size(1))*n] = (_epot(i-1,j,k) - _epot(i,j,k)) / h;
		}
	    }
	    
	    // If SYMMETRIC_POTENTIAL extrapolation, force E-field interpolation to
	    // zero when approaching boundary, otherwise extrapolate from last
	    // two nodes.
	    if( _extrpl[0] == FIELD_SYMMETRIC_POTENTIAL )
		_F[0][0+(j+k*_epot.size(1))*n] = -_F[0][1+(j+k*_epot.size(1))*n];
	    else
		_F[0][0+(j+k*_epot.size(1))*n] = 2.0*_F[0][1+(j+k*_epot.size(1))*n]-_F[0][2+(j+k*_epot.size(1))*n];
	    
	    if( _extrpl[1] == FIELD_SYMMETRIC_POTENTIAL )
		_F[0][n-1+(j+k*_epot.size(1))*n] = -_F[0][n-2+(j+k*_epot.size(1))*n];
	    else
		_F[0][n-1+(j+k*_epot.size(1))*n] = 2.0*_F[0][n-2+(j+k*_epot.size(1))*n]-_F[0][n-3+(j+k*_epot.size(1))*n];
	}
    }

    // Do Ey-field
    for( uint32_t k = 0; k < o-1; k++ ) {
	for( uint32_t i = 0; i < n-1; i++ ) {
	    for( uint32_t j = 1; j < m-1; j++ ) {
	
		uint32_t node1 = _geom->mesh(i,j-1,k);
		uint32_t node1id = node1 & SMESH_NODE_ID_MASK;
		uint32_t node2 = _geom->mesh(i,j,k);
		uint32_t node2id = node2 & SMESH_NODE_ID_MASK;

		if( node1id == SMESH_NODE_ID_NEAR_SOLID &&
		    node2id == SMESH_NODE_ID_DIRICHLET && 
		    (node2 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		
		    // Between near solid node and solid node
		    uint8_t dist = solid_dist( node1, 3 );
		    _F[1][i+(j+k*m)*_epot.size(0)] = 255.0*( _epot(i,j-1,k) - _epot(i,j,k) ) / (dist*h);

		} else if( node2id == SMESH_NODE_ID_NEAR_SOLID &&
			   node1id == SMESH_NODE_ID_DIRICHLET && 
			   (node1 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		
		    // Between near solid node and solid node
		    uint8_t dist = solid_dist( node2, 2 );
		    _F[1][i+(j+k*m)*_epot.size(0)] = 
			255.0*( _epot(i,j-1,k) - _epot(i,j,k) ) / (dist*h);
		
		} else if( node1id == SMESH_NODE_ID_DIRICHLET && 
			   (node1 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 &&
			   node2id == SMESH_NODE_ID_DIRICHLET && 
			   (node2 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		
		    // Inside solid, search for near solid in y-direction
		    uint32_t node0 = _geom->mesh_check(i,j-2,k);
		    uint32_t node0id = node0 & SMESH_NODE_ID_MASK;
		    uint32_t node3 = _geom->mesh_check(i,j+1,k);
		    uint32_t node3id = node3 & SMESH_NODE_ID_MASK;
		    if( node0id == SMESH_NODE_ID_NEAR_SOLID &&
			node3id == SMESH_NODE_ID_NEAR_SOLID ) {
			
			// No surface nearby
			_F[1][i+(j+k*m)*_epot.size(0)] = 0.0;
			
		    } else if( node0id == SMESH_NODE_ID_NEAR_SOLID ) {
			
			// Between near solid node and solid node
			uint8_t dist = solid_dist( node0, 3 );
			_F[1][i+(j+k*m)*_epot.size(0)] = 255.0*( _epot(i,j-2,k) - _epot(i,j-1,k) ) / (dist*h);
			
		    } else if( node3id == SMESH_NODE_ID_NEAR_SOLID ) {
			
			// Between near solid node and solid node
			uint8_t dist = solid_dist( node3, 2 );
			_F[1][i+(j+k*m)*_epot.size(0)] = 255.0*( _epot(i,j,k) - _epot(i,j+1,k) ) / (dist*h);
			
		    } else {
			
			// Ambiguous case, surface in both directions
			_F[1][i+(j+k*m)*_epot.size(0)] = 0.0;
		    }
		    
		} else {
		    
		    // Free space
		    _F[1][i+(j+k*m)*_epot.size(0)] = (_epot(i,j-1,k) - _epot(i,j,k)) / h;
		}
	    }

	    // If SYMMETRIC_POTENTIAL extrapolation, force E-field interpolation to
	    // zero when approaching boundary, otherwise extrapolate from last
	    // two nodes.
	    if( _extrpl[2] == FIELD_SYMMETRIC_POTENTIAL )
		_F[1][i+(k*m)*_epot.size(0)] = -_F[1][i+(1+k*m)*_epot.size(0)];
	    else
		_F[1][i+(k*m)*_epot.size(0)] = 2.0*_F[1][i+(1+k*m)*_epot.size(0)]-_F[1][i+(2+k*m)*_epot.size(0)];
	    
	    if( _extrpl[3] == FIELD_SYMMETRIC_POTENTIAL )
		_F[1][i+(m-1+k*m)*_epot.size(0)] = -_F[1][i+(m-2+k*m)*_epot.size(0)];
	    else
		_F[1][i+(m-1+k*m)*_epot.size(0)] = 2.0*_F[1][i+(m-2+k*m)*_epot.size(0)]-_F[1][i+(m-3+k*m)*_epot.size(0)];

	}
    }

    // Do Ez-field
    for( uint32_t i = 0; i < n-1; i++ ) {
	for( uint32_t j = 0; j < m-1; j++ ) {
	    for( uint32_t k = 1; k < o-1; k++ ) {
	
		uint32_t node1 = _geom->mesh(i,j,k-1);
		uint32_t node1id = node1 & SMESH_NODE_ID_MASK;
		uint32_t node2 = _geom->mesh(i,j,k);
		uint32_t node2id = node2 & SMESH_NODE_ID_MASK;

		if( node1id == SMESH_NODE_ID_NEAR_SOLID &&
		    node2id == SMESH_NODE_ID_DIRICHLET && 
		    (node2 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		
		    // Between near solid node and solid node
		    uint8_t dist = solid_dist( node1, 5 );
		    _F[2][i+(j+k*_epot.size(1))*_epot.size(0)] = 255.0*( _epot(i,j,k-1) - _epot(i,j,k) ) / (dist*h);

		} else if( node2id == SMESH_NODE_ID_NEAR_SOLID &&
			   node1id == SMESH_NODE_ID_DIRICHLET && 
			   (node1 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		
		    // Between near solid node and solid node
		    uint8_t dist = solid_dist( node2, 4 );
		    _F[2][i+(j+k*_epot.size(1))*_epot.size(0)] = 255.0*( _epot(i,j,k-1) - _epot(i,j,k) ) / (dist*h);
		
		} else if( node1id == SMESH_NODE_ID_DIRICHLET && 
			   (node1 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 &&
			   node2id == SMESH_NODE_ID_DIRICHLET && 
			   (node2 & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) {
		
		    // Inside solid, search for near solid in z-direction
		    uint32_t node0 = _geom->mesh_check(i,j,k-2);
		    uint32_t node0id = node0 & SMESH_NODE_ID_MASK;
		    uint32_t node3 = _geom->mesh_check(i,j,k+1);
		    uint32_t node3id = node3 & SMESH_NODE_ID_MASK;
		    if( node0id == SMESH_NODE_ID_NEAR_SOLID &&
			node3id == SMESH_NODE_ID_NEAR_SOLID ) {
			
			// No surface nearby
			_F[2][i+(j+k*_epot.size(1))*_epot.size(0)] = 0.0;
			
		    } else if( node0id == SMESH_NODE_ID_NEAR_SOLID ) {
			
			// Between near solid node and solid node
			uint8_t dist = solid_dist( node0, 5 );
			_F[2][i+(j+k*_epot.size(1))*_epot.size(0)] = 255.0*( _epot(i,j,k-2) - _epot(i,j,k-1) ) / (dist*h);
			
		    } else if( node3id == SMESH_NODE_ID_NEAR_SOLID ) {
			
			// Between near solid node and solid node
			uint8_t dist = solid_dist( node3, 4 );
			_F[2][i+(j+k*_epot.size(1))*_epot.size(0)] = 255.0*( _epot(i,j,k) - _epot(i,j,k+1) ) / (dist*h);
			
		    } else {
			
			// Ambiguous case, surface in both directions
			_F[2][i+(j+k*_epot.size(1))*_epot.size(0)] = 0.0;
		    }
		    
		} else {
		    
		    // Free space
		    _F[2][i+(j+k*_epot.size(1))*_epot.size(0)] = (_epot(i,j,k-1) - _epot(i,j,k)) / h;
		}
	    }

	    // If SYMMETRIC_POTENTIAL extrapolation, force E-field interpolation to
	    // zero when approaching boundary, otherwise extrapolate from last
	    // two nodes.
	    if( _extrpl[4] == FIELD_SYMMETRIC_POTENTIAL )
		_F[2][i+j*_epot.size(0)] = -_F[2][i+(j+_epot.size(1))*_epot.size(0)];
	    else
		_F[2][i+j*_epot.size(0)] = 2.0*_F[2][i+(j+_epot.size(1))*_epot.size(0)]-
		    _F[2][i+(j+2*_epot.size(1))*_epot.size(0)];
	    
	    if( _extrpl[5] == FIELD_SYMMETRIC_POTENTIAL )
		_F[2][i+(j+(o-1)*_epot.size(1))*_epot.size(0)] = -_F[2][i+(j+(o-2)*_epot.size(1))*_epot.size(0)];
	    else
		_F[2][i+(j+(o-1)*_epot.size(1))*_epot.size(0)] = 2.0*_F[2][i+(j+(o-2)*_epot.size(1))*_epot.size(0)]-
		    _F[2][i+(j+(o-3)*_epot.size(1))*_epot.size(0)];
	}
    }
}


/* Calculate electric field at points between potential nodes
 *
 *
 */
void EpotEfield::precalc( void )
{
    ibsimu.message( 1 ) << "Calculating efield\n";
    ibsimu.inc_indent();

    // Delete old fields
    for( size_t i = 0; i < 3; i++ ) {
	if( _F[i] != NULL )
	    delete [] _F[i];
    }

    switch( _epot.geom_mode() ) {
    case MODE_1D:
	precalc_1d();
	break;
    case MODE_2D:
    case MODE_CYL:
	precalc_2d();
	break;
    case MODE_3D:
	precalc_3d();
	break;
    }

    ibsimu.dec_indent();
}


void EpotEfield::debug_print( std::ostream &os ) const
{
    os << "**EpotEfield\n";
    
    os << "extrpl = (";
    for( uint32_t a = 0; a < 6; a++ ) {
	switch( _extrpl[a] ) {
	case FIELD_EXTRAPOLATE:
	    os << "FIELD_EXTRAPOLATE";
	    break;
	case FIELD_MIRROR:
	    os << "FIELD_MIRROR";
	    break;
	case FIELD_ANTIMIRROR:
	    os << "FIELD_ANTIMIRROR";
	    break;
	case FIELD_SYMMETRIC_POTENTIAL:
	    os << "FIELD_SYMMETRIC_POTENTIAL";
	    break;
	case FIELD_ZERO:
	    os << "FIELD_ZERO";
	    break;
	case FIELD_NAN:
	    os << "FIELD_NAN";
	    break;
	}
	if( a != 5 ) os << ", ";
    }
    os << ")\n";

    uint32_t nodecount = _geom->nodecount();
    for( uint32_t b = 0; b < 3; b++ ) {
	os << "F[" << b << "] = ";
	if( _F[b] == NULL ) {
	    os << "NULL\n";
	    continue;
	}
	os << "(";
	if( nodecount < 10 ) {
	    uint32_t a;
	    for( a = 0; a < nodecount-1; a++ )
		os << _F[b][a] << ", ";
	    if( a < nodecount )
		os << _F[b][a] << ")\n";
	} else {
	    // Print only 10 first nodes
	    for( uint32_t a = 0; a < 10; a++ )
		os << _F[b][a] << ", ";
	    os << "... )\n";
	}
    }
}


const Vec3D EpotEfield::operator()( const Vec3D &x ) const
{
    Vec3D R, X(x);
    Vec3D sign( 1.0, 1.0, 1.0 );

    switch( _geom->geom_mode() ) {
    case MODE_1D:
    {
	if( !_F[0] )
	    break;

	if( X[0] < _geom->origo(0) ) {
	    if( _extrpl[0] == FIELD_ZERO ) {
		R[0] = 0.0;
		break;
	    } else if( _extrpl[0] == FIELD_NAN ) {
		R[0] = std::numeric_limits<double>::quiet_NaN();
		break;
	    } else if( X[0] < _geom->origo(0)-_geom->size(0)*_geom->h() ) {
		// Outside double the simulation box: return NaN
		R[0] = std::numeric_limits<double>::quiet_NaN();
		break;
	    } else if( _extrpl[0] == FIELD_MIRROR ) {
		X[0] = 2.0*_geom->origo(0) - X[0];
	    } else if( _extrpl[0] == FIELD_ANTIMIRROR || 
		       _extrpl[0] == FIELD_SYMMETRIC_POTENTIAL ) {
		sign[0] *= -1.0;
		X[0] = 2.0*_geom->origo(0) - X[0];
	    }
	} else if( X[0] > _geom->max(0) ) {
	    if( _extrpl[1] == FIELD_ZERO ) {
		R[0] = 0.0;
		break;
	    } else if( _extrpl[1] == FIELD_NAN ) {
		R[0] = std::numeric_limits<double>::quiet_NaN();
		break;
	    } else if( X[0] > _geom->origo(0)+2.0*_geom->size(0)*_geom->h() ) {
		// Outside double the simulation box: return NaN
		R[0] = std::numeric_limits<double>::quiet_NaN();
		break;
	    } else if( _extrpl[1] == FIELD_MIRROR ) {
		X[0] = 2.0*_geom->max(0) - X[0];
	    } else if( _extrpl[1] == FIELD_ANTIMIRROR || 
		       _extrpl[1] == FIELD_SYMMETRIC_POTENTIAL ) {
		sign[0] *= -1.0;
		X[0] = 2.0*_geom->max(0) - X[0];
	    }
	}

	int32_t i = (int32_t)floor( (X[0]-_geom->origo(0))*_geom->div_h() + 0.5 );
	if( i < 0 )
	    i = 0;
	else if( i >= (int32_t)_geom->size(0) )
	    i = _geom->size(0)-1;

	double t = _geom->div_h()*( X[0]-((i-0.5)*_geom->h()+_geom->origo(0)) );

	R[0] = sign[0]*( (1.0-t)*_F[0][i] + t*_F[0][i+1] );
	break;
    }
    case MODE_2D:
    case MODE_CYL:
    {
	if( !_F[0] || !_F[1] )
	    break;

	for( int a = 0; a < 2; a++ ) {
	    if( X[a] < _geom->origo(a) ) {
		if( _extrpl[2*a] == FIELD_ZERO ) {
		    break;
		} else if( _extrpl[2*a] == FIELD_NAN ) {
		    R[0] = std::numeric_limits<double>::quiet_NaN();
		    R[1] = std::numeric_limits<double>::quiet_NaN();
		    break;
		} else if( X[a] < _geom->origo(a)-_geom->size(a)*_geom->h() ) {
		    // Outside double the simulation box: return NaN
		    R[0] = std::numeric_limits<double>::quiet_NaN();
		    R[1] = std::numeric_limits<double>::quiet_NaN();
		    break;
		} else if( _extrpl[2*a] == FIELD_MIRROR ) {
		    X[a] = 2.0*_geom->origo(a) - X[a];
		} else if( _extrpl[2*a] == FIELD_ANTIMIRROR || 
			   _extrpl[2*a] == FIELD_SYMMETRIC_POTENTIAL ) {
		    sign[a] *= -1.0;
		    X[a] = 2.0*_geom->origo(a) - X[a];
		}
	    } else if( X[a] > _geom->max(a) ) {
		if( _extrpl[2*a+1] == FIELD_ZERO ) {
		    break;
		} else if( _extrpl[2*a+1] == FIELD_NAN ) {
		    R[0] = std::numeric_limits<double>::quiet_NaN();
		    R[1] = std::numeric_limits<double>::quiet_NaN();
		    break;
		} else if( X[a] > _geom->origo(a)+2.0*_geom->size(a)*_geom->h() ) {
		    // Outside double the simulation box: return NaN
		    R[0] = std::numeric_limits<double>::quiet_NaN();
		    R[1] = std::numeric_limits<double>::quiet_NaN();
		    break;
		} else if( _extrpl[2*a+1] == FIELD_MIRROR ) {
		    X[a] = 2.0*_geom->max(a) - X[a];
		} else if( _extrpl[2*a+1] == FIELD_ANTIMIRROR || 
			   _extrpl[2*a+1] == FIELD_SYMMETRIC_POTENTIAL ) {
		    sign[a] *= -1.0;
		    X[a] = 2.0*_geom->max(a) - X[a];
		}
	    }
	}

	// Ex
	if( true ) {
	    int32_t i = (int32_t)floor( (X[0]-_geom->origo(0))*_geom->div_h() + 0.5 );
	    int32_t j = (int32_t)floor( (X[1]-_geom->origo(1))*_geom->div_h() );
	    if( i < 0 )
		i = 0;
	    else if( i >= (int32_t)_geom->size(0) )
		i = _geom->size(0)-1;
	    if( j < 0 )
		j = 0;
	    else if( j >= (int32_t)_geom->size(1)-1 )
		j = _geom->size(1)-2;
	    
	    double t = _geom->div_h()*( X[0]-((i-0.5)*_geom->h()+_geom->origo(0)) );
	    double u = _geom->div_h()*( X[1]-(j*_geom->h()+_geom->origo(1)) );
	    
	    size_t b = _geom->size(0)+1;
	    size_t c = i+j*b;
	    R[0] = sign[0]*( (1.0-u)*(1.0-t)*_F[0][c]   + 
		             (1.0-u)*     t *_F[0][c+1] +
		                  u *(1.0-t)*_F[0][c+b] +
			          u *     t *_F[0][c+b+1] );
	}

	// Ey
	if( true ) {
	    int32_t i = (int32_t)floor( (X[0]-_geom->origo(0))*_geom->div_h() );
	    int32_t j = (int32_t)floor( (X[1]-_geom->origo(1))*_geom->div_h() + 0.5 );
	    if( i < 0 )
		i = 0;
	    else if( i >= (int32_t)_geom->size(0)-1 )
		i = _geom->size(0)-2;
	    if( j < 0 )
		j = 0;
	    else if( j >= (int32_t)_geom->size(1) )
		j = _geom->size(1)-1;
	    
	    double t = _geom->div_h()*( X[0]-(i*_geom->h()+_geom->origo(0)) );
	    double u = _geom->div_h()*( X[1]-((j-0.5)*_geom->h()+_geom->origo(1)) );

	    size_t b = _geom->size(0);
	    size_t c = i+j*b;
	    R[1] = sign[1]*( (1.0-u)*(1.0-t)*_F[1][c]   + 
			     (1.0-u)*     t *_F[1][c+1] +
			          u *(1.0-t)*_F[1][c+b] +
			          u *     t *_F[1][c+b+1] );
	}
	break;
    }
    case MODE_3D:
    {
	if( !_F[0] || !_F[1] || !_F[2] )
	    break;

	for( int a = 0; a < 3; a++ ) {
	    if( X[a] < _geom->origo(a) ) {
		if( _extrpl[2*a] == FIELD_ZERO ) {
		    break;
		} else if( _extrpl[2*a] == FIELD_NAN ) {
		    R[0] = std::numeric_limits<double>::quiet_NaN();
		    R[1] = std::numeric_limits<double>::quiet_NaN();
		    R[2] = std::numeric_limits<double>::quiet_NaN();
		    break;
		} else if( X[a] < _geom->origo(a)-_geom->size(a)*_geom->h() ) {
		    // Outside double the simulation box: return NaN
		    R[0] = std::numeric_limits<double>::quiet_NaN();
		    R[1] = std::numeric_limits<double>::quiet_NaN();
		    R[2] = std::numeric_limits<double>::quiet_NaN();
		    break;
		} else if( _extrpl[2*a] == FIELD_MIRROR ) {
		    X[a] = 2.0*_geom->origo(a) - X[a];
		} else if( _extrpl[2*a] == FIELD_ANTIMIRROR || 
			   _extrpl[2*a] == FIELD_SYMMETRIC_POTENTIAL ) {
		    sign[a] *= -1.0;
		    X[a] = 2.0*_geom->origo(a) - X[a];
		}
	    } else if( X[a] > _geom->max(a) ) {
		if( _extrpl[2*a+1] == FIELD_ZERO ) {
		    break;
		} else if( _extrpl[2*a+1] == FIELD_NAN ) {
		    R[0] = std::numeric_limits<double>::quiet_NaN();
		    R[1] = std::numeric_limits<double>::quiet_NaN();
		    R[2] = std::numeric_limits<double>::quiet_NaN();
		    break;
		} else if( X[a] > _geom->origo(a)+2.0*_geom->size(a)*_geom->h() ) {
		    // Outside double the simulation box: return NaN
		    R[0] = std::numeric_limits<double>::quiet_NaN();
		    R[1] = std::numeric_limits<double>::quiet_NaN();
		    R[2] = std::numeric_limits<double>::quiet_NaN();
		    break;
		} else if( _extrpl[2*a+1] == FIELD_MIRROR ) {
		    X[a] = 2.0*_geom->max(a) - X[a];
		} else if( _extrpl[2*a+1] == FIELD_ANTIMIRROR || 
			   _extrpl[2*a+1] == FIELD_SYMMETRIC_POTENTIAL ) {
		    sign[a] *= -1.0;
		    X[a] = 2.0*_geom->max(a) - X[a];
		}
	    }
	}

	// Ex
	if( true ) {
	    int32_t i = (int32_t)floor( (X[0]-_geom->origo(0))*_geom->div_h() + 0.5 );
	    int32_t j = (int32_t)floor( (X[1]-_geom->origo(1))*_geom->div_h() );
	    int32_t k = (int32_t)floor( (X[2]-_geom->origo(2))*_geom->div_h() );
	    if( i < 0 )
		i = 0;
	    else if( i >= (int32_t)_geom->size(0) )
		i = _geom->size(0)-1;
	    if( j < 0 )
		j = 0;
	    else if( j >= (int32_t)_geom->size(1)-1 )
		j = _geom->size(1)-2;
	    if( k < 0 )
		k = 0;
	    else if( k >= (int32_t)_geom->size(2)-1 )
		k = _geom->size(2)-2;
	    
	    double t = _geom->div_h()*( X[0]-((i-0.5)*_geom->h()+_geom->origo(0)) );
	    double u = _geom->div_h()*( X[1]-(j*_geom->h()+_geom->origo(1)) );
	    double v = _geom->div_h()*( X[2]-(k*_geom->h()+_geom->origo(2)) );
	    
	    size_t a = _geom->size(0)+1;
	    size_t b = a*_geom->size(1);
	    size_t c = i+j*a+k*b;
	    R[0] = sign[0]*( (1.0-v)*(1.0-u)*(1.0-t)*_F[0][c]   + 
		             (1.0-v)*(1.0-u)*     t *_F[0][c+1] +
		             (1.0-v)*     u *(1.0-t)*_F[0][c+a] +
			     (1.0-v)*     u *     t *_F[0][c+a+1] +
			          v *(1.0-u)*(1.0-t)*_F[0][c+b]   + 
		                  v *(1.0-u)*     t *_F[0][c+b+1] +
		                  v *     u *(1.0-t)*_F[0][c+b+a] +
			          v *     u *     t *_F[0][c+b+a+1] );
	}

	// Ey
	if( true ) {
	    int32_t i = (int32_t)floor( (X[0]-_geom->origo(0))*_geom->div_h() );
	    int32_t j = (int32_t)floor( (X[1]-_geom->origo(1))*_geom->div_h() + 0.5 );
	    int32_t k = (int32_t)floor( (X[2]-_geom->origo(2))*_geom->div_h() );
	    if( i < 0 )
		i = 0;
	    else if( i >= (int32_t)_geom->size(0)-1 )
		i = _geom->size(0)-2;
	    if( j < 0 )
		j = 0;
	    else if( j >= (int32_t)_geom->size(1) )
		j = _geom->size(1)-1;
	    if( k < 0 )
		k = 0;
	    else if( k >= (int32_t)_geom->size(2)-1 )
		k = _geom->size(2)-2;
	    
	    double t = _geom->div_h()*( X[0]-(i*_geom->h()+_geom->origo(0)) );
	    double u = _geom->div_h()*( X[1]-((j-0.5)*_geom->h()+_geom->origo(1)) );
	    double v = _geom->div_h()*( X[2]-(k*_geom->h()+_geom->origo(2)) );

	    size_t a = _geom->size(0);
	    size_t b = a*(_geom->size(1)+1);
	    size_t c = i+j*a+k*b;
	    R[1] = sign[1]*( (1.0-v)*(1.0-u)*(1.0-t)*_F[1][c]   + 
		             (1.0-v)*(1.0-u)*     t *_F[1][c+1] +
		             (1.0-v)*     u *(1.0-t)*_F[1][c+a] +
			     (1.0-v)*     u *     t *_F[1][c+a+1] +
			          v *(1.0-u)*(1.0-t)*_F[1][c+b]   + 
		                  v *(1.0-u)*     t *_F[1][c+b+1] +
		                  v *     u *(1.0-t)*_F[1][c+b+a] +
			          v *     u *     t *_F[1][c+b+a+1] );
	}

	// Ez
	if( true ) {
	    int32_t i = (int32_t)floor( (X[0]-_geom->origo(0))*_geom->div_h() );
	    int32_t j = (int32_t)floor( (X[1]-_geom->origo(1))*_geom->div_h() );
	    int32_t k = (int32_t)floor( (X[2]-_geom->origo(2))*_geom->div_h() + 0.5 );
	    if( i < 0 )
		i = 0;
	    else if( i >= (int32_t)_geom->size(0)-1 )
		i = _geom->size(0)-2;
	    if( j < 0 )
		j = 0;
	    else if( j >= (int32_t)_geom->size(1)-1 )
		j = _geom->size(1)-2;
	    if( k < 0 )
		k = 0;
	    else if( k >= (int32_t)_geom->size(2) )
		k = _geom->size(2)-1;
	    
	    double t = _geom->div_h()*( X[0]-(i*_geom->h()+_geom->origo(0)) );
	    double u = _geom->div_h()*( X[1]-(j*_geom->h()+_geom->origo(1)) );
	    double v = _geom->div_h()*( X[2]-((k-0.5)*_geom->h()+_geom->origo(2)) );

	    size_t a = _geom->size(0);
	    size_t b = a*_geom->size(1);
	    size_t c = i+j*a+k*b;
	    R[2] = sign[2]*( (1.0-v)*(1.0-u)*(1.0-t)*_F[2][c]   + 
		             (1.0-v)*(1.0-u)*     t *_F[2][c+1] +
		             (1.0-v)*     u *(1.0-t)*_F[2][c+a] +
			     (1.0-v)*     u *     t *_F[2][c+a+1] +
			          v *(1.0-u)*(1.0-t)*_F[2][c+b]   + 
		                  v *(1.0-u)*     t *_F[2][c+b+1] +
		                  v *     u *(1.0-t)*_F[2][c+b+a] +
			          v *     u *     t *_F[2][c+b+a+1] );
	}
	break;
    }
    }

    return( R );
}



