/*! \file solidgraph.cpp
 *  \brief %Graph for plotting solids.
 */

/* Copyright (c) 2005-2012 Taneli Kalvas. All rights reserved.
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

#include <vector>
#include <limits>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include "compmath.hpp"
#include "solidgraph.hpp"
#include "vec3d.hpp"
#include "lineclip.hpp"
#include "ibsimu.hpp"


//#define DEBUG_SOLIDGRAPH 1


void SolidGraph::SolidPoints::add_point( Point x )
{
#ifdef DEBUG_SOLIDGRAPH
    std::cout << "  Add point (" << x[0] << ", " << x[1] << ")\n";
#endif

    _p.push_back( x );
}


bool SolidGraph::is_edge( uint32_t N, uint32_t node1, uint32_t node2 )
{
    return( (N == node1 && N != node2) ||
	    (N != node1 && N == node2) );
}


void SolidGraph::step( int32_t &nextx, int32_t &nexty, int32_t dir )
{
    if( dir == 0 ) {
	nexty = nexty+2;
    } else if( dir == 1 ) {
	if( nexty%2 == 0 ) {
	    nexty = nexty+1;
	} else {
	    nextx = nextx-1;
	    nexty = nexty+1;
	}
    } else if( dir == 2 ) {
	nextx = nextx-1;
    } else if( dir == 3 ) {
	if( nexty%2 == 0 ) {
	    nexty = nexty-1;
	} else {
	    nextx = nextx-1;
	    nexty = nexty-1;
	}
    } else if( dir == 4 ) {
	nexty = nexty-2;
    } else if( dir == 5 ) {
	if( nexty%2 == 0 ) {
	    nextx = nextx+1;
	    nexty = nexty-1;
	} else {
	    nexty = nexty-1;
	}
    } else if( dir == 6 ) {
	nextx = nextx+1;
    } else if( dir == 7 ) {
	if( nexty%2 == 0 ) {
	    nextx = nextx+1;
	    nexty = nexty+1;
	} else {
	    nexty = nexty+1;
	}
    }
}

/* Loop around solid and save points to \a solid. Uses near solid data
 * from Geometry for solid surface location when available, otherwise
 * draws on the solid edge node. Only does drawing in coordinate axes
 * directions from the surface. Direction of travel includes 45 degree
 * directions:
 *
 *   1 0 7
 *   2 X 6
 *   3 4 5
 */
void SolidGraph::loop( SolidPoints *solid, int32_t x, int32_t y,
		       char *done, uint32_t sizex, uint32_t sizey )
{
#ifdef DEBUG_SOLIDGRAPH    
    std::cout << "\nSolidGraph::loop()\n";
#endif

    // Construct node id
    uint32_t nodeid = SMESH_NODE_ID_DIRICHLET | solid->N();

    // Save starting point
    int32_t startx = x;
    int32_t starty = y;

    // Initialize mesh coordinates
    int32_t i1[3];
    int32_t i2[3];
    i1[_vb[2]] = _level;
    i2[_vb[2]] = _level;

    // Get mesh coordinates
    bool even = (y%2 == 0);
    i1[_vb[0]] = x-1;
    i1[_vb[1]] = y/2-1;
    if( even ) {
	i2[_vb[0]] = i1[_vb[0]]+1;
	i2[_vb[1]] = i1[_vb[1]];
    } else {
	i2[_vb[0]] = i1[_vb[0]];
	i2[_vb[1]] = i1[_vb[1]]+1;
    }
    // Get node IDs
    uint32_t node1 = _geom.mesh_check( i1[0], i1[1], i1[2] );
    uint32_t node2 = _geom.mesh_check( i2[0], i2[1], i2[2] );

    // Loop around solid
    do {

#ifdef DEBUG_SOLIDGRAPH
	std::cout << "  At (x,y) = (" << x << ", " << y << ")\n";
	std::cout << "        i1 = (" << i1[0] << ", " << i1[1] << ", " << i1[2] << ")\n";
	std::cout << "        i2 = (" << i2[0] << ", " << i2[1] << ", " << i2[2] << ")\n";
	if( node1 == nodeid )
	    std::cout << "     node1 = " << node1 << " (inside)\n";
	else
	    std::cout << "     node1 = " << node1 << " (outside)\n";
	if( node2 == nodeid )
	    std::cout << "     node2 = " << node2 << " (inside)\n";
	else
	    std::cout << "     node2 = " << node2 << " (outside)\n";
#endif

	// Mark done
	done[x+y*sizex] = 1;

	// Save the point and get next dir
	int32_t dir;
	if( even ) {
	    if( node1 == nodeid ) {
		// Solid on left, go up-right
#ifdef DEBUG_SOLIDGRAPH
		std::cout << "  Solid on left\n";
#endif
		dir = 7; 
		if( (node2 & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEAR_SOLID ) {
		    int32_t soliddir = 2*_vb[0];
#ifdef DEBUG_SOLIDGRAPH
		    std::cout << "  soliddir = " << soliddir << "\n";
#endif
		    double dist = _geom.solid_dist( i2[0], i2[1], i2[2], soliddir )/255.0;
		    solid->add_point( Point( (i2[_vb[0]]-dist)*_geom.h()+_geom.origo(_vb[0]),
					     i2[_vb[1]]*_geom.h()+_geom.origo(_vb[1]) ) );
		} else {
		    solid->add_point( Point( (i2[_vb[0]]-0.5)*_geom.h()+_geom.origo(_vb[0]),
					     i2[_vb[1]]*_geom.h()+_geom.origo(_vb[1]) ) );
		}
	    } else {
		// Solid on right, go down-left
#ifdef DEBUG_SOLIDGRAPH
		std::cout << "  Solid on right\n";
#endif
		dir = 3;
		if( (node1 & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEAR_SOLID ) {
		    int32_t soliddir = 2*_vb[0]+1;
#ifdef DEBUG_SOLIDGRAPH
		    std::cout << "  soliddir = " << soliddir << "\n";
#endif
		    double dist = _geom.solid_dist( i1[0], i1[1], i1[2], soliddir )/255.0;
		    solid->add_point( Point( (i1[_vb[0]]+dist)*_geom.h()+_geom.origo(_vb[0]),
					     i1[_vb[1]]*_geom.h()+_geom.origo(_vb[1]) ) );
		} else {
		    solid->add_point( Point( (i1[_vb[0]]+0.5)*_geom.h()+_geom.origo(_vb[0]),
					     i1[_vb[1]]*_geom.h()+_geom.origo(_vb[1]) ) );
		}
	    }
	} else {
	    if( node1 == nodeid ) {
		// Solid on bottom, go up-left
#ifdef DEBUG_SOLIDGRAPH
		std::cout << "  Solid on bottom\n";
#endif
		dir = 1;
		if( (node2 & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEAR_SOLID ) {
		    int32_t soliddir = 2*_vb[1];
#ifdef DEBUG_SOLIDGRAPH
		    std::cout << "  soliddir = " << soliddir << "\n";
#endif
		    double dist = _geom.solid_dist( i2[0], i2[1], i2[2], soliddir )/255.0;
		    solid->add_point( Point( i2[_vb[0]]*_geom.h()+_geom.origo(_vb[0]),
					     (i2[_vb[1]]-dist)*_geom.h()+_geom.origo(_vb[1]) ) );
		} else {
		    solid->add_point( Point( i2[_vb[0]]*_geom.h()+_geom.origo(_vb[0]),
					     (i2[_vb[1]]-0.5)*_geom.h()+_geom.origo(_vb[1]) ) );
		}
	    } else {
		// Solid on top, go down-right
#ifdef DEBUG_SOLIDGRAPH
		std::cout << "  Solid on top\n";
#endif
		dir = 5;
		if( (node1 & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEAR_SOLID ) {
		    int32_t soliddir = 2*_vb[1]+1;
#ifdef DEBUG_SOLIDGRAPH
		    std::cout << "  soliddir = " << soliddir << "\n";
#endif
		    double dist = _geom.solid_dist( i1[0], i1[1], i1[2], soliddir )/255.0;
		    solid->add_point( Point( i1[_vb[0]]*_geom.h()+_geom.origo(_vb[0]),
					     (i1[_vb[1]]+dist)*_geom.h()+_geom.origo(_vb[1]) ) );
		} else {
		    solid->add_point( Point( i1[_vb[0]]*_geom.h()+_geom.origo(_vb[0]),
					     (i1[_vb[1]]+0.5)*_geom.h()+_geom.origo(_vb[1]) ) );
		}
	    }
	}

#ifdef DEBUG_SOLIDGRAPH
	    std::cout << "  Find next step\n";
#endif

	// Find next step
	for( int32_t a = 0; a < 8; a++ ) {

	    int32_t nextx = x;
	    int32_t nexty = y;
	    step( nextx, nexty, dir );

	    // Get mesh coordinates
	    even = (nexty%2 == 0);
	    i1[_vb[0]] = nextx-1;
	    i1[_vb[1]] = nexty/2-1;
	    if( even ) {
		i2[_vb[0]] = i1[_vb[0]]+1;
		i2[_vb[1]] = i1[_vb[1]];
	    } else {
		i2[_vb[0]] = i1[_vb[0]];
		i2[_vb[1]] = i1[_vb[1]]+1;
	    }
	    // Get node id's    
	    node1 = _geom.mesh_check( i1[0], i1[1], i1[2] );
	    node2 = _geom.mesh_check( i2[0], i2[1], i2[2] );

#ifdef DEBUG_SOLIDGRAPH
	    std::cout << "    Probing dir = " << dir << "\n";
	    std::cout << "      (x,y) = (" << nextx << ", " << nexty << ")\n";	
	    std::cout << "         i1 = (" << i1[0] << ", " << i1[1] << ", " << i1[2] << ")\n";
	    std::cout << "         i2 = (" << i2[0] << ", " << i2[1] << ", " << i2[2] << ")\n";
	    if( node1 == nodeid )
		std::cout << "        node1 = " << node1 << " (inside)\n";
	    else
		std::cout << "        node1 = " << node1 << " (outside)\n";
	    if( node2 == nodeid )
		std::cout << "        node2 = " << node2 << " (inside)\n";
	    else
		std::cout << "        node2 = " << node2 << " (outside)\n";
#endif

	    // Accept if edge
	    if( is_edge( nodeid, node1, node2 ) ) {
#ifdef DEBUG_SOLIDGRAPH
		std::cout << "    Accepting step\n";
#endif
		x = nextx;
		y = nexty;
		break;
	    }
	    // Otherwise proceed to next direction
	    dir++;
	    if( dir > 7 )
		dir -= 8;
	}

    } while( x != startx || y != starty );

}


void SolidGraph::clear_data( void )
{
    for( size_t a = 0; a < _solid.size(); a++  ) {
	delete _solid[a];
    }
    _solid.clear();
}


void SolidGraph::build_solid( uint32_t N, char *done, uint32_t sizex, uint32_t sizey )
{
#ifdef DEBUG_SOLIDGRAPH
    std::cout << "build_solid()\n";
#endif
    
    // Create points container
    SolidPoints *solid = new SolidPoints( N );
    uint32_t nodeid = SMESH_NODE_ID_DIRICHLET | N;

    // Go through done array
    int32_t i1[3];  // Mesh coordinates
    int32_t i2[3];  // Mesh coordinates
    i1[_vb[2]] = _level;
    i2[_vb[2]] = _level;
    for( int32_t y = 0; y < (int32_t)sizey; y++ ) {
	bool even = (y%2 == 0);
	uint32_t xcount = (even ? sizex-1 : sizex);
	for( int32_t x = 0; x < (int32_t)xcount; x++ ) {

#ifdef DEBUG_SOLIDGRAPH
	    std::cout << "Processing x = " << x 
		      << ", y = " << y << "\n";
#endif

	    // Skip if already processed
	    if( done[x+y*sizex] )
		continue;

	    i1[_vb[0]] = x-1;
	    i1[_vb[1]] = y/2-1;
	    if( even ) {
		i2[_vb[0]] = i1[_vb[0]]+1;
		i2[_vb[1]] = i1[_vb[1]];
	    } else {
		i2[_vb[0]] = i1[_vb[0]];
		i2[_vb[1]] = i1[_vb[1]]+1;
	    }

#ifdef DEBUG_SOLIDGRAPH
	    std::cout << "Processing i1 = (" 
		      << i1[0] << "," 
		      << i1[1] << "," 
		      << i1[2] << "), i2 = ("
		      << i2[0] << "," 
		      << i2[1] << "," 
		      << i2[2] << ")\n";
#endif

	    uint32_t node1 = _geom.mesh_check( i1[0], i1[1], i1[2] );
	    uint32_t node2 = _geom.mesh_check( i2[0], i2[1], i2[2] );
	    if( !is_edge( nodeid, node1, node2 ) )
		continue;

#ifdef DEBUG_SOLIDGRAPH
	    std::cout << "Found edge\n";
#endif

	    if( solid->size() != 0 ) {
		// Add break before new loop
		solid->add_point( Point(std::numeric_limits<double>::quiet_NaN(),
					std::numeric_limits<double>::quiet_NaN()) );
	    }
	    loop( solid, x, y, done, sizex, sizey );
	}
    }

    // Save points container
    _solid.push_back( solid );
}


SolidGraph::SolidGraph( const Geometry &geom ) 
    : Graph3D(geom), _geom(geom), _color(Vec3D(0.2,0.2,1.0)), _cache(true)
{
}


SolidGraph::~SolidGraph()
{
    for( size_t a = 0; a < _solid.size(); a++ )
	delete _solid[a];
}


void SolidGraph::disable_cache( void )
{
    _cache = false;
}


void SolidGraph::build_data( void )
{
#ifdef DEBUG_SOLIDGRAPH
    std::cout << "SolidGraph::build_data()\n";
#endif

    // Clear old data
    clear_data();

    // Array for marking processed nodes
    // size + 2 in both directions for marking over the boundary
    // 2*size-1 in y direction for half steps
    uint32_t sizex = _geom.size(_vb[0])+2;
    uint32_t sizey = 2*(_geom.size(_vb[1])+2) - 1;

    // Reserve done array
    char *done = new char[sizex*sizey];

    // Build data for all solids
    for( size_t a = 0; a < _geom.number_of_solids(); a++ ) {

	// Clear done array
	memset( done, 0, sizex*sizey*sizeof(char) );

	build_solid( 7+a, done, sizex, sizey );
    }

    // Free done array
    delete [] done;
}


void SolidGraph::plot_data( cairo_t *cairo, LineClip &lc, const Coordmapper *cm )
{
#ifdef DEBUG_SOLIDGRAPH
    std::cout << "SolidGraph::plot_data()\n";
#endif

    // Plot all solids
    for( size_t a = 0; a < _solid.size(); a++ ) {

	SolidPoints &sp = *_solid[a];
	if( sp.size() == 0 )
	    continue;

	double xout[2];
	cm->transform( xout, sp[0].x );
	lc.move_to( xout[0], xout[1] );

	for( size_t b = 1; b < sp.size(); b++ ) {
	    if( comp_isnan( sp[b].x[0] ) ) {
		// Break in path, separate paths
		do b++;
		while( b != sp.size() && comp_isnan( sp[b].x[0] ) );
		if( b == sp.size() )
		    break;
		lc.close_path();
		cm->transform( xout, sp[b].x );
		lc.move_to( xout[0], xout[1] );
	    } else {
		cm->transform( xout, sp[b].x );
		lc.line_to( xout[0], xout[1] );
	    }
	}

	lc.fill();
    }
}


void SolidGraph::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] )
{
#ifdef DEBUG_SOLIDGRAPH
    std::cout << "SolidGraph::plot()\n";
#endif

    if( !_cache || _solid.size() == 0 || _oview != _view || _olevel != _level ) {
	// First round or change happened
	build_data();
    }
    _oview = _view;
    _olevel = _level;

    // Set drawing properties
    cairo_set_source_rgba( cairo, _color[0], _color[1], _color[2], 1.0 );
    cairo_set_line_width( cairo, 1.0 );

    // Set clipping ranges
    double clip[4];
    cm->transform( &clip[0], &range[0] );
    cm->transform( &clip[2], &range[2] );
    LineClip lc( cairo );
    lc.set( clip[0], clip[1], clip[2], clip[3] );

    plot_data( cairo, lc, cm );
}


void SolidGraph::plot_sample( cairo_t *cairo, double x, double y, double width, double height )
{

}


void SolidGraph::get_bbox( double bbox[4] )
{
    bbox[0] = _geom.origo( _vb[0] );
    bbox[1] = _geom.origo( _vb[1] );
    bbox[2] = _geom.max( _vb[0] );
    bbox[3] = _geom.max( _vb[1] );
}
