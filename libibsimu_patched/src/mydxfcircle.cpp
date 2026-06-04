/*! \file mydxfcircle.hpp
 *  \brief DXF circle entity
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
#include "mydxfcircle.hpp"
#include "mydxfline.hpp"
#include "polysolver.hpp"


#define MAX_ARC_CHOP 4


MyDXFCircle::MyDXFCircle( class MyDXFFile *dxf )
{
#ifdef MYDXF_DEBUG
    std::cout << "  Reading entity CIRCLE\n";
#endif

    // Default values
    _pc[0] = _pc[1] = _pc[2] = 0.0;
    _r = 1.0;

    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with entity

	else if( dxf->group_get_code() == 10 )
	    _pc[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 )
	    _pc[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 )
	    _pc[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 40 )
	    _r = dxf->group_get_double();

	else
	    process_group( dxf );
    }

#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


void MyDXFCircle::set_center( const Vec3D &c )
{
    _pc = c;
}


void MyDXFCircle::set_radius( double r )
{
    _r = r;
}


void MyDXFCircle::explode( MyDXFEntities *ent, class MyDXFFile *dxf, const Transformation *t ) const
{
    if( (*t)[0] == (*t)[5] && (*t)[0] == (*t)[10] ) {
	// Representable as an circle
	MyDXFCircle *circle = new MyDXFCircle( *this );
	
	// Transform points
	Vec3D c = t->transform_point( _pc );
	Vec3D r = t->transform_point( _pc + Vec3D(_r,0,0) );
	circle->_pc = c;
	circle->_r = norm2( c - r );

	// Add to entities
	ent->add_entity( circle );
    } else {
	// Not representable as circle

	// Chop circle to N pieces
	int N = 16;
	Vec3D old;
	for( int i = 0; i < N; i++ ) {
	    double a = i*2.0*M_PI/(N-1.0);
	    Vec3D x = t->transform_point( Vec3D( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] ) );

	    // Add line to entities
	    if( i > 0 ) {
		MyDXFLine *line = new MyDXFLine( *this );
		line->set_start( old );
		line->set_end( x );
		ent->add_entity( line );
	    }
	    old = x;
	}
    }
}


void MyDXFCircle::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "CIRCLE" );
    write_common( dxf, ostr );

    dxf->write_group( 10, _pc[0] );
    dxf->write_group( 20, _pc[1] );
    dxf->write_group( 30, _pc[2] );

    dxf->write_group( 40, _r );
}


void MyDXFCircle::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
			const Transformation *t, const double range[4] ) const
{
#ifdef MYDXF_DEBUG_PLOT
    std::cout << "MyDXFCircle::plot()\n";
#endif

    std::vector<double> cx;
    std::vector<double> cy;

    // Divide circle to linear pieces and keep length of linear pieces
    // to less than 4 px

    // Start with three points
    double cf = 2.0*M_PI/3.0;
    for( int i = 0; i < 3; i++ ) {
	double a = i*cf;
	Vec4D p = t->transform( Vec3D( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] ) );
	cx.push_back( p[0] );
	cy.push_back( p[1] );
    }

    // Loop until error is small enough
    for( int iter = 0; iter < MAX_ARC_CHOP; iter++ ) {

	//std::cout << "LOOP: N = " << cx.size() << "\n";
	//for( int i = 0; i < (int)cx.size(); i++ )
	//std::cout << "  " << cx[i] << "\t" << cy[i] << "\n";

	// Calculate maximum step length
	double mstep = 0.0;
	for( int i = 0; i < (int)cx.size()-1; i++ ) {
	    double sx = cx[i+1] - cx[i];
	    double sy = cy[i+1] - cy[i];
	    double step = sx*sx + sy*sy;
	    if( step > mstep )
		mstep = step;
	}
	// last gap
	double sx = cx[0] - cx[cx.size()-1];
	double sy = cy[0] - cy[cx.size()-1];
	double step = sx*sx + sy*sy;
	if( step > mstep )
	    mstep = step;

	// Calculate maximum error between most recent points and the
	// linear interpolation from the last round
	/*
	double maxerr = 0.0;
	for( int i = 1; i < (int)cx.size(); i += 2 ) {
	    double tx = 0.5*(cx[i-1] + cx[i+1]);
	    double ty = 0.5*(cy[i-1] + cy[i+1]);
	    tx -= cx[i];
	    ty -= cy[i];
	    double err = tx*tx + ty*ty; // square of error
	    if( err > maxerr )
		maxerr = err;
	}
	*/

	//std::cout << "mstep = " << mstep << "\n";

	// If step < 4.0
	if( mstep < 16.0 )
	    break;

	// Add more points
	cx.resize( 2*cx.size() );
	cy.resize( 2*cy.size() );
	double cf = 2.0*M_PI/cx.size();
	for( int i = cx.size()-1; i > 0; i-- ) {
	    
	    // New point
	    double a = i*cf;
	    Vec3D p = t->transform_point( Vec3D( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] ) );
	    cx[i] = p[0];
	    cy[i] = p[1];
	    i--;

	    // Old point
	    cx[i] = cx[i/2];
	    cy[i] = cy[i/2];
	}

	//std::cout << "\n";
    }

#ifdef MYDXF_DEBUG_PLOT
    for( int i = 0; i < (int)cx.size(); i++ )
	std::cout << "MyDXFCircle::plot(): plotting (" << cx[i] << ", " << cy[i] << ")\n";
#endif

    // Error small enough -> draw lines
    cairo_move_to( cairo, cx[0], cy[0] );
    for( int i = 1; i < (int)cx.size(); i++ )
	cairo_line_to( cairo, cx[i], cy[i] );
    cairo_line_to( cairo, cx[0], cy[0] );
    cairo_stroke( cairo );
}


void MyDXFCircle::get_bbox( Vec3D &min, Vec3D &max, 
			    const class MyDXFFile *dxf, const Transformation *t ) const
{
    min = Vec3D( std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    max = Vec3D( -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity() );

    // Chop circle to 8 pieces for bbox
    double cf = 2.0*M_PI/8.0;
    for( int i = 0; i < 8; i++ ) {
	double a = i*cf;
	Vec3D p( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] );
	Vec4D x = t->transform( Vec3D( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] ) );
	bbox_ppoint( min, max, x );
    }

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "Circle bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFCircle::scale( class MyDXFFile *dxf, double s )
{
    _pc *= s;
    _r  *= s;
}


void MyDXFCircle::translate( class MyDXFFile *dxf, const Vec3D &dx )
{
    _pc += dx;
}


void MyDXFCircle::rotate_z( class MyDXFFile *dxf, double a )
{
    Transformation t = Transformation::rotation_z( a );
    _pc = t.transform_point( _pc );
}


int MyDXFCircle::ray_cross( double x, double y ) const
{
    double dx = x-_pc[0];

    // Uncertainty regions at circle bounds in x-direction and close to arc ends
    if( (fabs( dx-_r ) < MYDXF_PERT_EPS || fabs( dx+_r ) < MYDXF_PERT_EPS) && 
	y >= _pc[1]-MYDXF_PERT_EPS )
	return( 2 );

    // Outside circle bounds in x-direction
    if( (dx <= -_r || dx >= _r) )
	return( 0 );	// No crossing

    double dy = sqrt( _r*_r - dx*dx );
    double cy1 = _pc[1] + dy;
    double cy2 = _pc[1] - dy;

    if( y <= cy2 )
	return( 0 );
    else if( y >= cy1 )
	return( 0 );

    return( 1 );
}


void MyDXFCircle::debug_print( std::ostream &os ) const
{
    os << "  CIRCLE\n";
    MyDXFFile::debug_print_format( os, "pc", _pc );
    MyDXFFile::debug_print_format( os, "r", _r );
}


bool MyDXFCircle::geom_same( const MyDXFCircle &circle, double eps ) const
{
    return( norm2( _pc - circle._pc ) < eps &&
	    fabs( _r - circle._r ) < eps );
}




