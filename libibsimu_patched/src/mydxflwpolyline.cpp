/*! \file mydxflwpolyline.cpp
 *  \brief DXF lwpolyline entity
 */

/* Copyright (c) 2010-2012,2014 Taneli Kalvas. All rights reserved.
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
#include "mydxflwpolyline.hpp"


MyDXFLWPolyline::MyDXFLWPolyline( class MyDXFFile *dxf )
{
#ifdef MYDXF_DEBUG
    std::cout << "  Reading entity LWPOLYLINE\n";
#endif

    int32_t N = 0;
    int32_t i = -1;
    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with entity

	else if( dxf->group_get_code() == 10 ) {
	    i++;
	    if( i == N )
		throw Error( ERROR_LOCATION, "Incorrect amount of vertices in LWPOLYLINE"
			     " on line " + dxf->linec() );
	    _p[i][0] = dxf->group_get_double();
	} else if( dxf->group_get_code() == 20 )
	    _p[i][1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 42 )
	    _p[i][2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 70 )
	    _flags = dxf->group_get_int16();
	else if( dxf->group_get_code() == 90 ) {
	    // Resize to number of vertices.
	    if( i != -1 )
		throw Error( ERROR_LOCATION, "Number of vertices specified after first "
			     "vertedin LWPOLYLINE on line " + dxf->linec() );
	    N = dxf->group_get_int32();
	    _p.resize( N );
	}

	else
	    process_group( dxf );
    }

    if( i != N-1 )
	throw Error( ERROR_LOCATION, "Incorrect amount of vertices in LWPOLYLINE"
		     " on line " + dxf->linec() );

#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


void MyDXFLWPolyline::explode( MyDXFEntities *ent, class MyDXFFile *dxf, const Transformation *t ) const
{
    MyDXFLWPolyline *line = new MyDXFLWPolyline( *this );

    // Transform points
    for( uint32_t a = 0; a < line->_p.size(); a++ )
	line->_p[a] = t->transform_point( line->_p[a] );

    // Add to entities
    ent->add_entity( line );
}


void MyDXFLWPolyline::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "LWPOLYLINE" );
    write_common( dxf, ostr );

    dxf->write_group( 90, (int32_t)_p.size() );

    for( uint32_t i = 0; i < _p.size(); i++ ) {
	dxf->write_group( 10, _p[i][0] );
	dxf->write_group( 20, _p[i][1] );
	dxf->write_group( 42, _p[i][2] );
    }

    dxf->write_group( 70, _flags );
}


void MyDXFLWPolyline::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
			    const Transformation *t, const double range[4] ) const
{
    if( _p.size() == 0 )
	return;

    Vec4D x = t->transform( _p[0] );
    cairo_move_to( cairo, x[0], x[1] );
    for( uint32_t a = 1; a < _p.size(); a++ ) {
	x = t->transform( _p[a] );
	cairo_line_to( cairo, x[0], x[1] );
    }
    if( closed() ) {
	x = t->transform( _p[0] );
	cairo_line_to( cairo, x[0], x[1] );
    }

    cairo_stroke( cairo );
}


void MyDXFLWPolyline::get_bbox( Vec3D &min, Vec3D &max, 
				const class MyDXFFile *dxf, const Transformation *t ) const
{
    min = Vec3D( std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    max = Vec3D( -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity() );

    for( uint32_t a = 0; a < _p.size(); a++ ) {
	Vec4D x = t->transform( _p[a] );
	bbox_ppoint( min, max, x );
    }

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "Polyline bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFLWPolyline::scale( class MyDXFFile *dxf, double s  )
{
    for( uint32_t a = 0; a < _p.size(); a++ )
	_p[a] *= s;
}


void MyDXFLWPolyline::translate( class MyDXFFile *dxf, const Vec3D &dx  )
{
    for( uint32_t a = 0; a < _p.size(); a++ )
	_p[a] += dx;
}


void MyDXFLWPolyline::rotate_z( class MyDXFFile *dxf, double a )
{
    Transformation t = Transformation::rotation_z( a );
    for( uint32_t a = 0; a < _p.size(); a++ )
	_p[a] = t.transform_point( _p[a] );
}


int MyDXFLWPolyline::ray_cross( double x, double y ) const
{
#ifdef MYDXF_DEBUG
    std::cout << "ray_cross( x = " << x << ", y = " << y << " )\n";
    for( uint32_t a = 0; a < _p.size(); a++ )
	std::cout << "  p[" << a << "] = " << _p[a] << "\n";
    std::cout << "\n";
#endif

    Vec3D p1 = _p[0];

    // Uncertainty region around p1
    if( fabs( x - p1[0] ) < MYDXF_PERT_EPS && y >= p1[1]-MYDXF_PERT_EPS )
	return( 2 );

    // Go through all lines 
    int c = 0;
    for( uint32_t a = 1; a < _p.size(); a++ ) {
	Vec3D p2 = _p[a];

#ifdef MYDXF_DEBUG
	std::cout << "  Test against:\n";
	std::cout << "  p1 = " << p1 << "\n";
	std::cout << "  p2 = " << p2 << "\n";
#endif

	// Uncertainty region around p2
	if( fabs( x - p2[0] ) < MYDXF_PERT_EPS && y >= p2[1]-MYDXF_PERT_EPS )
	    return( 2 );
	
	// Ray within the line ends in x-direction
	if( (x > p1[0] && x < p2[0]) || 
	    (x < p1[0] && x > p2[0]) ) {

	    // Calculate crossing y-coordinate.
	    double t = (x-p1[0])/(p2[0]-p1[0]);
	    double cy = (1.0-t)*p1[1] + t*p2[1];
	    // Boundary case y == cy is considered crossing.
	    if( y >= cy ) {
		c = !c; // Crossing
	    } else {
		// No crossing
	    }
	}

	p1 = p2;
    }

    return( c );
}


void MyDXFLWPolyline::debug_print( std::ostream &os ) const
{
    os << "  LWPOLYLINE\n";
    MyDXFFile::debug_print_format( os, "flags", _flags );
    for( uint32_t a = 0; a < _p.size(); a++ )
	MyDXFFile::debug_print_format( os, "p", _p[a] );
}


void MyDXFLWPolyline::set_start( const Vec3D &s )
{
    _p[0] = s;
}


void MyDXFLWPolyline::set_end( const Vec3D &e )
{
    _p[_p.size()-1] = e;
}


bool MyDXFLWPolyline::geom_same( const MyDXFLWPolyline &line, double eps ) const
{
    if( _p.size() != line._p.size() )
	return( 1 );

    for( uint32_t a = 0; a < _p.size(); a++ ) {
	if( norm2( _p[a] - line._p[a] ) < eps )
	    return( 1 );
    }

    return( 0 );
}


