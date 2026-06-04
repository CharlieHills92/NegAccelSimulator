/*! \file eqpotgraph.cpp
 *  \brief %Graph for plotting equipotential lines
 */

/* Copyright (c) 2005-2013 Taneli Kalvas. All rights reserved.
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

#include "eqpotgraph.hpp"
#include "lineclip.hpp"
#include "ibsimu.hpp"


EqPotGraph::EqPotGraph( const MeshScalarField &epot, const Geometry &geom )
    : Graph3D(geom), _color(Vec3D(0.2,1,0.2)), _epot(epot), _geom(geom), 
      _data_built(false), _cache(true)
{
}


EqPotGraph::~EqPotGraph()
{
    for( size_t a = 0; a < _lines.size(); a++ )
	delete _lines[a];
}


void EqPotGraph::set_eqlines_manual( const std::vector<double> &pot )
{
    _data_built = false;
    _eqlines_manual = pot;
}


void EqPotGraph::set_eqlines_auto( size_t N )
{
    _data_built = false;
    _eqlines_auto = N;
}


bool EqPotGraph::eqline_exists( double pot1, uint32_t sol1, 
				double pot2, uint32_t sol2, 
				double pot ) const 
{
    if( pot1 != pot2 && ((pot1 >= pot && pot2 <= pot) ||
			 (pot1 <= pot && pot2 >= pot)) ) {
	return( true );
    }
    return( false );
}


bool EqPotGraph::is_solid( uint32_t sol ) const
{
    return( (sol & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET && 
	    (sol & SMESH_BOUNDARY_NUMBER_MASK) >= 7 );
}


bool EqPotGraph::is_near( uint32_t sol ) const
{
    return( (sol & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEAR_SOLID );
}


void EqPotGraph::build_data( void )
{
    _data_built = true;

    // Clear old data
    for( size_t a = 0; a < _lines.size(); a++ )
	delete _lines[a];
    _lines.clear();

    // Build vector of required manual potential lines
    for( size_t a = 0; a < _eqlines_manual.size(); a++ )
	_lines.push_back( new EqPotLines( _eqlines_manual[a] ) );

    // Add automatic lines
    double min, max;
    _epot.get_minmax( min, max );
    for( size_t a = 0; a < _eqlines_auto; a++ )
	_lines.push_back( new EqPotLines( min + (a+0.5)*(max-min)/(_eqlines_auto) ) );

    // Go through mesh
    size_t i[3];
    size_t dx, dy;
    i[_vb[2]] = _level;
    if( _vb[0] == 0 ) dx = 1;
    else if( _vb[0] == 1 ) dx = _geom.size(0);
    else dx = _geom.size(1)*_geom.size(0);
    if( _vb[1] == 0 ) dy = 1;
    else if( _vb[1] == 1 ) dy = _geom.size(0);
    else dy = _geom.size(1)*_geom.size(0);
    size_t sizex = _geom.size(_vb[0])-1;
    size_t sizey = _geom.size(_vb[1])-1;
    for( i[_vb[1]] = 0; i[_vb[1]] < sizey; i[_vb[1]]++ ) {
	for( i[_vb[0]] = 0; i[_vb[0]] < sizex; i[_vb[0]]++ ) {
	    
	    //std::cout << "Scanning at " << i[_vb[0]] << " " << i[_vb[1]] << "\n";
	    size_t ptr = _geom.size(0)*_geom.size(1)*i[2] + _geom.size(0)*i[1] + i[0];

	    // Potentials and solid mesh numbers counterclockwise
	    // starting from lower left corner in vb[0], vb[1] coordinates
	    size_t ptr1       = ptr;
	    double pot1       = _epot( ptr1 );
	    uint32_t sol1     = _geom.mesh( ptr1 );

	    size_t ptr2       = ptr+dx;
	    double pot2       = _epot( ptr2 );
	    uint32_t sol2     = _geom.mesh( ptr2 );

	    size_t ptr3       = ptr+dx+dy;
	    double pot3       = _epot( ptr3 );
	    uint32_t sol3     = _geom.mesh( ptr3 );

	    size_t ptr4       = ptr+dy;
	    double pot4       = _epot( ptr4 );
	    uint32_t sol4     = _geom.mesh( ptr4 );

	    // Check for existance of each equipotential line
	    for( size_t a = 0; a < _lines.size(); a++ ) {

		size_t b = 0;
		double pot = _lines[a]->pot;
		double x[8];

		//std::cout << "pot = " << pot << "\n";

		// Bottom
		if( eqline_exists( pot1, sol1, pot2, sol2, pot ) ) {
		    if( is_near(sol1) && is_solid(sol2) ) {
			// Get solid distance for sol1 to direction of sol2 (+vb[0])
			uint8_t dist = _geom.solid_dist( ptr1, 2*_vb[0]+1 );
			double t = dist*(pot-pot1) / (255*(pot2-pot1));
			x[b+0] = ( i[_vb[0]] + t ) * _geom.h() + _geom.origo(_vb[0]);
		    } else if( is_near(sol2) && is_solid(sol1) ) {
			// Get solid distance for sol2 to direction of sol1 (-vb[0])
			uint8_t dist = _geom.solid_dist( ptr2, 2*_vb[0] );
			double t = dist*(pot-pot2) / (255*(pot1-pot2));
			x[b+0] = ( i[_vb[0]] + 1 - t ) * _geom.h() + _geom.origo(_vb[0]);
		    } else {
			double t = (pot-pot1) / (pot2-pot1);
			x[b+0] = ( i[_vb[0]] + t ) * _geom.h() + _geom.origo(_vb[0]);
		    }
		    x[b+1] = i[_vb[1]]*_geom.h() + _geom.origo(_vb[1]);
		    b += 2;
		}

		// Right
		if( eqline_exists( pot2, sol2, pot3, sol3, pot ) ) {
		    if( is_near(sol2) && is_solid(sol3) ) {
			// Get solid distance for sol2 to direction of sol3 (+vb[1])
			uint8_t dist = _geom.solid_dist( ptr2, 2*_vb[1]+1 );
			double t = dist*(pot-pot2) / (255*(pot3-pot2));
			x[b+1] = ( i[_vb[1]] + t ) * _geom.h() + _geom.origo(_vb[1]);
		    } else if( is_near(sol3) && is_solid(sol2) ) {
			// Get solid distance for sol3 to direction of sol2 (-vb[1])
			uint8_t dist = _geom.solid_dist( ptr3, 2*_vb[1] );
			double t = dist*(pot-pot3) / (255*(pot2-pot3));
			x[b+1] = ( i[_vb[1]] + 1 - t ) * _geom.h() + _geom.origo(_vb[1]);
		    } else {
			double t = (pot-pot2) / (pot3-pot2);
			x[b+1] = ( i[_vb[1]] + t ) * _geom.h() + _geom.origo(_vb[1]);
		    }
		    x[b+0] = (i[_vb[0]]+1)*_geom.h() + _geom.origo(_vb[0]);
		    b += 2;
		}
		
		// Top
		if( eqline_exists( pot4, sol4, pot3, sol3, pot ) ) {
		    if( is_near(sol4) && is_solid(sol3) ) {
			// Get solid distance for sol4 to direction of sol3 (+vb[0])
			uint8_t dist = _geom.solid_dist( ptr4, 2*_vb[0]+1 );
			double t = dist*(pot-pot4) / (255*(pot3-pot4));
			x[b+0] = ( i[_vb[0]] + t ) * _geom.h() + _geom.origo(_vb[0]);
		    } else if( is_near(sol3) && is_solid(sol4) ) {
			// Get solid distance for sol3 to direction of sol4 (-vb[0])
			uint8_t dist = _geom.solid_dist( ptr3, 2*_vb[0] );
			double t = dist*(pot-pot3) / (255*(pot4-pot3));
			x[b+0] = ( i[_vb[0]] + 1 - t ) * _geom.h() + _geom.origo(_vb[0]);
		    } else {
			double t = (pot-pot4) / (pot3-pot4);
			x[b+0] = ( i[_vb[0]] + t ) * _geom.h() + _geom.origo(_vb[0]);
		    }
		    x[b+1] = (i[_vb[1]]+1)*_geom.h() + _geom.origo(_vb[1]);
		    b += 2;
		}

		// Left
		if( eqline_exists( pot1, sol1, pot4, sol4, pot ) ) {
		    if( is_near(sol1) && is_solid(sol4) ) {
			// Get solid distance for sol1 to direction of sol4 (+vb[1])
			uint8_t dist = _geom.solid_dist( ptr1, 2*_vb[1]+1 );
			double t = dist*(pot-pot1) / (255*(pot4-pot1));
			x[b+1] = ( i[_vb[1]] + t ) * _geom.h() + _geom.origo(_vb[1]);
		    } else if( is_near(sol4) && is_solid(sol1) ) {
			// Get solid distance for sol4 to direction of sol1 (-vb[1])
			uint8_t dist = _geom.solid_dist( ptr4, 2*_vb[1] );
			double t = dist*(pot-pot4) / (255*(pot1-pot4));
			x[b+1] = ( i[_vb[1]] + 1 - t ) * _geom.h() + _geom.origo(_vb[1]);
		    } else {
			double t = (pot-pot1) / (pot4-pot1);
			x[b+1] = ( i[_vb[1]] + t ) * _geom.h() + _geom.origo(_vb[1]);
		    }
		    x[b+0] = i[_vb[0]]*_geom.h() + _geom.origo(_vb[0]);
		    b += 2;
		}

		if( b == 4 ) {
		    _lines[a]->x.push_back( Line( x[0], x[1], x[2], x[3] ) );
		} else if( b == 8 ) {
		    // Six possible connections possible for four
		    // points We don't want to draw across. 
		    // Equipotential lines don't cross -> out
		    // of four lines we draw two
		    _lines[a]->x.push_back( Line( x[0], x[1], x[2], x[3] ) );
		    _lines[a]->x.push_back( Line( x[4], x[5], x[6], x[7] ) );
		}
	    }
	}
    }
}


void EqPotGraph::disable_cache( void )
{
    _cache = false;
}


void EqPotGraph::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] )
{
    if( !_data_built || _oview != _view || _olevel != _level || !_cache ) {
	// First round or change happened
	build_data();
    }
    _oview = _view;
    _olevel = _level;

    // Set drawing properties
    cairo_set_source_rgba( cairo, _color[0], _color[1], _color[2], 1.0 );
    cairo_set_line_width( cairo, 1.0 );
    LineClip lc( cairo );

    for( size_t a = 0; a < _lines.size(); a++ ) {
	for( size_t b = 0; b < _lines[a]->x.size(); b++ ) {

	    double xout[2];
	    cm->transform( xout, &_lines[a]->x[b][0] );
	    //std::cout << "From:"  << xout[0] << " " << xout[1] << "\n";
	    lc.move_to( xout[0], xout[1] );
	    cm->transform( xout, &_lines[a]->x[b][2] );
	    //std::cout << "To:  " << xout[0] << " " << xout[1] << "\n";
	    lc.line_to( xout[0], xout[1] );
	    cairo_stroke( cairo );
	}
    }
}


void EqPotGraph::plot_sample( cairo_t *cairo, double x, double y, double width, double height )
{

}


void EqPotGraph::get_bbox( double bbox[4] )
{
    bbox[0] = _geom.origo( _vb[0] );
    bbox[1] = _geom.origo( _vb[1] );
    bbox[2] = _geom.max( _vb[0] );
    bbox[3] = _geom.max( _vb[1] );
}
