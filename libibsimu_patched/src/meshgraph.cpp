/*! \file meshgraph.cpp
 *  \brief %Graph of rectangular mesh for geometry plots
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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

#include "meshgraph.hpp"
#include "ibsimu.hpp"


MeshGraph::MeshGraph( const Geometry &geom )
    : Graph3D(geom), _geom(geom)
{

}


MeshGraph::~MeshGraph()
{

}


void MeshGraph::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] )
{
    cairo_set_source_rgb( cairo, 0.5, 0.5, 0.5 );
    cairo_set_line_width( cairo, 1.0 );

    int irange[4];
    irange[0] = (int)floor((range[0]-_geom.origo(_vb[0])) / _geom.h());
    irange[1] = (int)floor((range[1]-_geom.origo(_vb[1])) / _geom.h());
    irange[2] = (int)ceil((range[2]-_geom.origo(_vb[0])) / _geom.h());
    irange[3] = (int)ceil((range[3]-_geom.origo(_vb[1])) / _geom.h());
    
    double xin[2];
    double xout[2];

    for( int i = irange[0]; i <= irange[2]; i++ ) {
	
	//std::cout << "mesh line at x=" << _geom.origo(_vb[0]) + _geom.h() * i << "\n";
	
	xin[0] = _geom.origo(_vb[0]) + _geom.h() * i;
	xin[1] = range[1];
	cm->transform( xout, xin );
	cairo_move_to( cairo, xout[0], xout[1] );

	xin[0] = _geom.origo(_vb[0]) + _geom.h() * i;
	xin[1] = range[3];
	cm->transform( xout, xin );
	cairo_line_to( cairo, xout[0], xout[1] );
    }

    for( int j = irange[1]; j <= irange[3]; j++ ) {
	
	//std::cout << "mesh line at y=" << _geom.origo(_vb[1]) + _geom.h() * j << "\n";
	
	xin[0] = range[0];
	xin[1] = _geom.origo(_vb[1]) + _geom.h() * j;
	cm->transform( xout, xin );
	cairo_move_to( cairo, xout[0], xout[1] );

	xin[0] = range[2];
	xin[1] = _geom.origo(_vb[1]) + _geom.h() * j;
	cm->transform( xout, xin );
	cairo_line_to( cairo, xout[0], xout[1] );
    }

    cairo_stroke( cairo );
}


void MeshGraph::plot_sample( cairo_t *cairo, double x, double y, double width, double height )
{

}


void MeshGraph::get_bbox( double bbox[4] )
{
    bbox[0] = _geom.origo( _vb[0] );
    bbox[1] = _geom.origo( _vb[1] );
    bbox[2] = _geom.max( _vb[0] );
    bbox[3] = _geom.max( _vb[1] );
}
