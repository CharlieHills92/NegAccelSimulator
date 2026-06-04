/*! \file mydxfline.hpp
 *  \brief DXF line entity
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


#include "mydxfline.hpp"


//#define DEBUG_INSIDE_LOOP 1


MyDXFLine::MyDXFLine( class MyDXFFile *dxf )
{
#ifdef MYDXF_DEBUG
    std::cout << "  Reading entity LINE\n";
#endif

    // Default values
    _p1[0] = _p1[1] = _p1[2] = 0.0;
    _p2[0] = _p2[1] = _p2[2] = 0.0;

    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with entity

	else if( dxf->group_get_code() == 10 )
	    _p1[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 )
	    _p1[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 )
	    _p1[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 11 )
	    _p2[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 21 )
	    _p2[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 31 )
	    _p2[2] = dxf->group_get_double();
	else
	    process_group( dxf );
    }

#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


void MyDXFLine::explode( MyDXFEntities *ent, class MyDXFFile *dxf, const Transformation *t ) const
{
    MyDXFLine *line = new MyDXFLine( *this );

    // Transform points
    line->_p1 = t->transform_point( line->_p1 );
    line->_p2 = t->transform_point( line->_p2 );

    // Add to entities
    ent->add_entity( line );
}


void MyDXFLine::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "LINE" );
    write_common( dxf, ostr );

    dxf->write_group( 10, _p1[0] );
    dxf->write_group( 20, _p1[1] );
    dxf->write_group( 30, _p1[2] );

    dxf->write_group( 11, _p2[0] );
    dxf->write_group( 21, _p2[1] );
    dxf->write_group( 31, _p2[2] );
}


void MyDXFLine::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		      const Transformation *t, const double range[4] ) const
{
#ifdef MYDXF_DEBUG_PLOT
    std::cout << "MyDXFLine::plot()\n";
#endif

    Vec3D x1 = t->transform_point( _p1 );
    Vec3D x2 = t->transform_point( _p2 );

#ifdef MYDXF_DEBUG_PLOT
    std::cout << "MyDXFLine::plot(): plotting (" << x1[0] << ", " << x1[1] << ")\n";
    std::cout << "MyDXFLine::plot(): plotting (" << x2[0] << ", " << x2[1] << ")\n";
#endif

    cairo_move_to( cairo, x1[0], x1[1] );
    cairo_line_to( cairo, x2[0], x2[1] );

    cairo_stroke( cairo );
}


void MyDXFLine::get_bbox( Vec3D &min, Vec3D &max, 
			  const class MyDXFFile *dxf, const Transformation *t ) const
{
    Vec4D x1 = t->transform( _p1 );
    Vec4D x2 = t->transform( _p2 );

    for( int a = 0; a < 3; a++ ) {
	if( x1[a] < x2[a] ) {
	    min[a] = x1[a];
	    max[a] = x2[a];
	} else {
	    min[a] = x2[a];
	    max[a] = x1[a];
	}
    }	

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "Line bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFLine::scale( class MyDXFFile *dxf, double s  )
{
    _p1 *= s;
    _p2 *= s;
}


void MyDXFLine::translate( class MyDXFFile *dxf, const Vec3D &dx  )
{
    _p1 += dx;
    _p2 += dx;
}


void MyDXFLine::rotate_z( class MyDXFFile *dxf, double a )
{
    Transformation t = Transformation::rotation_z( a );
    _p1 = t.transform_point( _p1 );
    _p2 = t.transform_point( _p2 );
}


int MyDXFLine::ray_cross( double x, double y ) const
{
#ifdef DEBUG_INSIDE_LOOP
    std::cout << "  MyDXFLine::ray_cross( x = " << x << ", y = " << y << " )\n";
    std::cout << "  p1 = " << _p1 << "\n";
    std::cout << "  p2 = " << _p2 << "\n";
#endif

    // Uncertainty regions at line ends
    if( (fabs( x - _p1[0] ) < MYDXF_PERT_EPS && y >= _p1[1]-MYDXF_PERT_EPS) ||
	(fabs( x - _p2[0] ) < MYDXF_PERT_EPS && y >= _p2[1]-MYDXF_PERT_EPS) )
	return( 2 );

    // Ray within the line ends in x-direction
    if( (x > _p1[0] && x < _p2[0]) || 
	(x < _p1[0] && x > _p2[0]) ) {

	// Calculate crossing y-coordinate.
	double t = (x-_p1[0])/(_p2[0]-_p1[0]);
        double cy = (1.0-t)*_p1[1] + t*_p2[1];

#ifdef DEBUG_INSIDE_LOOP
	std::cout << "  cy = " << cy << "\n";
#endif

        // Boundary case y == cy is considered crossing.
        if( y >= cy )
            return( 1 );
    }

#ifdef DEBUG_INSIDE_LOOP
	std::cout << "  no crossing\n";
#endif

    // No crossing.
    return( 0 );
}


void MyDXFLine::debug_print( std::ostream &os ) const
{
    os << "  LINE\n";
    MyDXFFile::debug_print_format( os, "p1", _p1 );
    MyDXFFile::debug_print_format( os, "p2", _p2 );
}


void MyDXFLine::set_start( const Vec3D &s )
{
    _p1 = s;
}


void MyDXFLine::set_end( const Vec3D &e )
{
    _p2 = e;
}


bool MyDXFLine::geom_same( const MyDXFLine &line, double eps ) const
{
    return( (norm2( _p1 - line._p1 ) < eps &&
	     norm2( _p2 - line._p2 ) < eps) ||
	    (norm2( _p2 - line._p1 ) < eps &&
	     norm2( _p1 - line._p2 ) < eps) );
}


