/*! \file particlegraph.cpp
 *  \brief %Graph for particle plots
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

/*! \file particlegraph.cpp
 *  \brief Particle trajectory plotter.
 */

#include <string.h>
#include <iostream>
#include "particlegraph.hpp"
#include "ibsimu.hpp"


ParticleGraph::ParticleGraph( const Geometry &geom, const ParticleDataBase &pdb, 
			      uint32_t particle_div, uint32_t particle_offset,
			      bool qm_discr )
    : Graph3D(geom), _geom(geom), _pdb(pdb), 
      _particle_div(particle_div), _particle_offset(particle_offset),
      _coordsize(3), _qm_discr(qm_discr)
{
    // Allocate work space
    _coord = new double[2*_coordsize];

    // Add default colors
    _color.push_back( Vec3D( 1.0, 0.2, 0.2 ) ); // Red

    _color.push_back( Vec3D( 0.2, 1.0, 1.0 ) ); // Cyan
    _color.push_back( Vec3D( 1.0, 0.2, 1.0 ) ); // Magenta
    _color.push_back( Vec3D( 1.0, 1.0, 0.2 ) ); // Yellow

    _color.push_back( Vec3D( 1.0, 0.5, 0.2 ) ); // Orange
    _color.push_back( Vec3D( 0.5, 0.2, 1.0 ) ); // Purple
    _color.push_back( Vec3D( 0.2, 1.0, 0.5 ) ); // Bluish green
    _color.push_back( Vec3D( 1.0, 0.2, 0.5 ) ); // Pink

}


ParticleGraph::~ParticleGraph()
{
    delete [] _coord;
}


void ParticleGraph::set_particle_div( uint32_t particle_div, uint32_t particle_offset )
{
    _particle_div = particle_div;
    _particle_offset = particle_offset;
}


void ParticleGraph::set_qm_discretation( bool qm_discr )
{
    _qm_discr = qm_discr;
}


/*  Subroutine for drawing linear interpolation line between two
 *  particle trajectory data points.  Here x = (x,y,vx,vy,t) is a
 *  point in the trajectory. Parameter first indicates first point of
 *  the trajectory.
 */
void ParticleGraph::draw_linear( const Coordmapper *cm, LineClip &lc, 
				 double x[5], bool first ) const
{
    double xout[2];
    cm->transform( xout, x );
    if( first )
	lc.move_to( xout[0], xout[1] );
    else
	lc.line_to( xout[0], xout[1] );
}


void ParticleGraph::get_point( const Coordmapper *cm, double *coord, double s, 
			       double Ax, double Bx, double Cx, double Dx, 
			       double Ay, double By, double Cy, double Dy ) const
{
    double x[2] = { ((Ax*s + Bx)*s + Cx)*s + Dx, 
		    ((Ay*s + By)*s + Cy)*s + Dy };
    cm->transform( coord, x );
}


/*  Subroutine for drawing curved line (consisting of several straight
 *  lines) between two particle trajectory data points.  Here x =
 *  (x,y,vx,vy,t) is a point in the trajectory. Parameter first
 *  indicates first point of the trajectory.
 */
void ParticleGraph::draw_curve( const Coordmapper *cm, LineClip &lc, 
				double x[5], bool first )
{
    if( first ) {
	// Save first point
	memcpy( _ox, x, 5*sizeof(double) );
	return;
    }

    // Calculate polynomial coefficients
    double dt = x[4]-_ox[4];

    double Ax = (x[2]+_ox[2])*dt + 2.0*(_ox[0]-x[0]);
    double Bx = 3.0*(x[0]-_ox[0]) - (2.0*_ox[2]+x[2])*dt;
    double Cx = _ox[2]*dt;
    double Dx = _ox[0];

    double Ay = (x[3]+_ox[3])*dt + 2.0*(_ox[1]-x[1]);
    double By = 3.0*(x[1]-_ox[1]) - (2.0*_ox[3]+x[3])*dt;
    double Cy = _ox[3]*dt;
    double Dy = _ox[1];

    // Fill database with three points
    size_t coord_N = 3;
    cm->transform( &_coord[0], &_ox[0] );
    get_point( cm, &_coord[2], 0.5, Ax, Bx, Cx, Dx, Ay, By, Cy, Dy );
    cm->transform( &_coord[4], &x[0] );
    
    // Bisect curve until desired accuracy is achieved
    size_t a;
    double maxerr, err, xt, yt;
    while( 1 ) {

	//for( a = 0; a < coord_N; a++ ) {
	//    std::cout << std::setw(8) << _coord[2*a+0]
	//	      << std::setw(8) << _coord[2*a+1] << "\n";
	//}

	// Calculate maximum error between the most recent points and
	// the linear interpolation of points from the last round.
	maxerr = 0.0;
	for( a = 1; a < coord_N; a += 2 ) {
	    xt = 0.5*(_coord[2*(a-1)] + _coord[2*(a+1)]);
	    yt = 0.5*(_coord[2*(a-1)+1] + _coord[2*(a+1)+1]);
	    xt -= _coord[2*a];
	    yt -= _coord[2*a+1];
	    err = sqrt( xt*xt + yt*yt );
	    if( err > maxerr )
		maxerr = err;
	}
	
	// If maximum error is less than 2.5 pixels, it is good enough
	// Or if number of coordinates is 33 or more
	if( maxerr < 2.5 || coord_N >= 33 ) {
	    break;
	}

	// Add more points
	coord_N += coord_N-1;

	// Allocate more space if necessary
	if( _coordsize < coord_N ) {
	    double *ct = new double[2*coord_N];
	    memcpy( ct, _coord, _coordsize*2*sizeof(double) );
	    delete [] _coord;
	    _coordsize = coord_N;
	    _coord = ct;
	}

	// Fill coordinates
	for( a = coord_N-1; a > 0; a-- ) {
	    
	    // Even
	    _coord[2*a]   = _coord[a];
	    _coord[2*a+1] = _coord[a+1];
	    a--;
      
	    // Odd
	    get_point( cm, &_coord[2*a], (double)a/(coord_N-1), 
		       Ax, Bx, Cx, Dx, Ay, By, Cy, Dy );
	}
    }
    
    // Draw lines
    lc.move_to( _coord[0], _coord[1] );
    for( a = 1; a < coord_N; a++ )
	lc.line_to( _coord[2*a+0], _coord[2*a+1] );

    memcpy( _ox, x, 5*sizeof(double) );
}


void ParticleGraph::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] )
{
    // No plotting
    if( _particle_div == 0 )
	return;

    // Q/M discriminator set
    std::vector<double> qm_set;
    if( !_qm_discr )
	cairo_set_source_rgb( cairo, _color[0][0], _color[0][1], _color[0][2] );
    cairo_set_line_width( cairo, 1.0 );

    // Set clipping ranges
    double clip[4];
    cm->transform( &clip[0], &range[0] );
    cm->transform( &clip[2], &range[2] );
    LineClip lc( cairo );
    lc.set( clip[0], clip[1], clip[2], clip[3] );

    // Loop through all particles
    for( size_t a = _particle_offset; a < _pdb.size(); a += _particle_div ) {

	// No plotting if one or less trajectory points
	if( _pdb.traj_size( a ) <= 1 )
	    continue;

	// Select color for particle
	if( _qm_discr ) {
	    size_t c;
	    const ParticleBase &p = _pdb.particle(a);
	    for( c = 0; c < qm_set.size(); c++ ) {
		volatile double pqm = p.qm();
		if( pqm == qm_set[c] )
		    break;
	    }
	    if( c == qm_set.size() )
		qm_set.push_back( p.qm() ); // New q/m

	    // Set color
	    size_t s = c%_color.size();
	    cairo_set_source_rgb( cairo, _color[s][0], _color[s][1], _color[s][2] );
	}

	// Loop through all particle trajectory points
	for( size_t b = 0; b < _pdb.traj_size( a ); b++ ) {

	    double t;
	    Vec3D loc, vel;
	    _pdb.trajectory_point( t, loc, vel, a, b );
	    double x[5] = { loc(_vb[0]), loc(_vb[1]), 
			    vel(_vb[0]), vel(_vb[1]), t };
	    if( _pdb.get_polyint() )
		draw_curve( cm, lc, x, b == 0 );
	    else
		draw_linear( cm, lc, x, b == 0 );

	}
	cairo_stroke( cairo );
    }
}


void ParticleGraph::plot_sample( cairo_t *cairo, double x, double y, double width, double height )
{

}


void ParticleGraph::get_bbox( double bbox[4] )
{
    bbox[0] = _geom.origo( _vb[0] );
    bbox[1] = _geom.origo( _vb[1] );
    bbox[2] = _geom.max( _vb[0] );
    bbox[3] = _geom.max( _vb[1] );
}


void ParticleGraph::add_color( const Vec3D &color )
{
    _color.push_back( color );
}


void ParticleGraph::clear_colors( void )
{
    _color.clear();
}
