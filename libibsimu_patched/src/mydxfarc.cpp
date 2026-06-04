/*! \file mydxfarc.hpp
 *  \brief DXF arc entity
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
#include "mydxfarc.hpp"
#include "mydxfline.hpp"
#include "polysolver.hpp"


#define MAX_ARC_CHOP 16
///#define DEBUG_INSIDE_LOOP 1


MyDXFArc::MyDXFArc( class MyDXFFile *dxf )
{
#ifdef MYDXF_DEBUG
    std::cout << "  Reading entity ARC\n";
#endif

    // Default values
    _pc[0] = _pc[1] = _pc[2] = 0.0;
    _r = 1.0;
    _ang1 = 0.0;
    _ang2 = 2.0*M_PI;

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

	else if( dxf->group_get_code() == 50 ) {
	    double ang1 = M_PI*dxf->group_get_double()/180.0;
	    // Enforce between 0 and 2 pi
	    _ang1 = ang1 - 2.0*M_PI*floor( ang1/(2.0*M_PI) );
	} else if( dxf->group_get_code() == 51 ) {
	    double ang2 = M_PI*dxf->group_get_double()/180.0;
	    // Enforce between 0 and 2 pi
	    _ang2 = ang2 - 2.0*M_PI*floor( ang2/(2.0*M_PI) );
	}

	else
	    process_group( dxf );
    }

#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


void MyDXFArc::explode( MyDXFEntities *ent, class MyDXFFile *dxf, const Transformation *t ) const
{
    if( (*t)[0] == (*t)[5] && (*t)[0] == (*t)[10] ) {
	// Representable as an arc
	MyDXFArc *arc = new MyDXFArc( *this );
	
	// Transform points
	Vec3D s = t->transform_point( start() );
	Vec3D e = t->transform_point( end() );
	Vec3D c = t->transform_point( center() );
	arc->set_center_and_ends( c, s, e );

	// Add to entities
	ent->add_entity( arc );
    } else {
	// Not representable as arc

	// Chop arc to N pieces
	double adiff = _ang2-_ang1;
	if( adiff < 0.0 )
	    adiff = 2.0*M_PI + adiff;
	
	Vec3D old;
	int N = 16;
	for( int i = 0; i < N; i++ ) {
	    double a = _ang1 + i*adiff/(N-1.0);
	    if( a > 2.0*M_PI )
		a -= 2.0*M_PI;
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


void MyDXFArc::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "ARC" );
    write_common( dxf, ostr );

    dxf->write_group( 10, _pc[0] );
    dxf->write_group( 20, _pc[1] );
    dxf->write_group( 30, _pc[2] );

    dxf->write_group( 40, _r );
    dxf->write_group( 50, 180*_ang1/M_PI );
    dxf->write_group( 51, 180*_ang2/M_PI );
}


void MyDXFArc::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		     const Transformation *t, const double range[4] ) const
{
    std::vector<double> cx;
    std::vector<double> cy;

    // Divide arc to linear pieces and keep length of linear pieces
    // to less than 4 px

    //std::cout << "\n\nARC plot\n";
    //std::cout << "  ang1 = " << _ang1 << "\n";
    //std::cout << "  ang2 = " << _ang2 << "\n";

    // Start with three points: ends and one in between
    double adiff = _ang2-_ang1;
    if( adiff < 0.0 )
	adiff = 2.0*M_PI + adiff;
    //std::cout << "  adiff = " << adiff << "\n";

    for( int i = 0; i < 3; i++ ) {
	double a = _ang1 + i*adiff/2.0;
	if( a > 2.0*M_PI )
	    a -= 2.0*M_PI;
	Vec3D p = t->transform_point( Vec3D( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] ) );
	//std::cout << "  a = " << a << "\n";
	//std::cout << "  p = " << p[0] << "\t" << p[1] << "\n";
	cx.push_back( p[0] );
	cy.push_back( p[1] );
    }

    // Loop until error is small enough
    for( int iter = 0; iter < MAX_ARC_CHOP; iter++ ) {

	//std::cout << "LOOP: N = " << cx.size() << "\n";
	//for( int i = 0; i < (int)cx.size(); i++ )
	//    std::cout << "  " << cx[i] << "\t" << cy[i] << "\n";

	// Calculate maximum step length
	double mstep = 0.0;
	for( int i = 0; i < (int)cx.size()-1; i++ ) {
	    double sx = cx[i+1] - cx[i];
	    double sy = cy[i+1] - cy[i];
	    double step = sx*sx + sy*sy;
	    if( step > mstep )
		mstep = step;
	}
	//std::cout << "mstep = " << mstep << "\n";
	//std::cout << "\n";

	// If step < 4.0
	if( mstep < 4.0 )
	    break;

	// Add more points
	cx.resize( 2*cx.size()-1 );
	cy.resize( 2*cy.size()-1 );
	double cf = adiff/(cx.size()-1);
	for( int i = cx.size()-1; i > 1; i-- ) {

	    // Old point
	    cx[i] = cx[i/2];
	    cy[i] = cy[i/2];
	    i--;

	    // New point
	    double a = _ang1 + i*cf;
	    if( a > 2.0*M_PI )
		a -= 2.0*M_PI;
	    Vec3D p = t->transform_point( Vec3D( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] ) );
	    //std::cout << "  a = " << a << "\n";
	    //std::cout << "  p = " << p[0] << "\t" << p[1] << "\n";
	    cx[i] = p[0];
	    cy[i] = p[1];
	}
    }

    // Error small enough -> draw lines
    cairo_move_to( cairo, cx[0], cy[0] );
    for( int i = 0; i < (int)cx.size(); i++ )
	cairo_line_to( cairo, cx[i], cy[i] );
    cairo_stroke( cairo );
}


void MyDXFArc::get_bbox( Vec3D &min, Vec3D &max, 
			 const class MyDXFFile *dxf, const Transformation *t ) const
{
    min = Vec3D( std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    max = Vec3D( -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity() );

    // Chop arc to 8 pieces for bbox
    double adiff = _ang2-_ang1;
    if( adiff < 0.0 )
	adiff = 2.0*M_PI + adiff;

    for( int i = 0; i < 8; i++ ) {
	double a = _ang1 + i*adiff/7.0;
	if( a > 2.0*M_PI )
	    a -= 2.0*M_PI;
	Vec3D x = t->transform_point( Vec3D( _pc[0]+_r*cos(a), _pc[1]+_r*sin(a), _pc[2] ) );
	bbox_ppoint( min, max, x );
    }

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "Arc bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFArc::scale( class MyDXFFile *dxf, double s )
{
    _pc *= s;
    _r  *= s;
}


void MyDXFArc::translate( class MyDXFFile *dxf, const Vec3D &dx )
{
    _pc += dx;
}


void MyDXFArc::rotate_z( class MyDXFFile *dxf, double a )
{
    Transformation t = Transformation::rotation_z( a );
    Vec3D s = t.transform_point( start() );
    Vec3D e = t.transform_point( end() );
    Vec3D c = t.transform_point( center() );
    set_center_and_ends( c, s, e );
}


int MyDXFArc::ray_cross( double x, double y ) const
{
#ifdef DEBUG_INSIDE_LOOP
    std::cout << "  MyDXFArc::ray_cross( x = " << x << ", y = " << y << " )\n";
    std::cout << "  pc = " << _pc << "\n";
    std::cout << "  r  = " << _r << "\n";
#endif

    // Quick test for points definitely under the circle -> no crossing
    if( y <= _pc[1]-_r-MYDXF_PERT_EPS )
	return( 0 );

    double dx = x-_pc[0];

    // Uncertainty regions at circle bounds in x-direction and close to arc ends
    if( ((fabs( dx-_r ) < MYDXF_PERT_EPS || fabs( dx+_r ) < MYDXF_PERT_EPS) && 
	 y >= _pc[1]-MYDXF_PERT_EPS) ||
	fabs( dx-_r*sin(_ang1) ) < MYDXF_PERT_EPS || 
	fabs( dx-_r*sin(_ang2) ) < MYDXF_PERT_EPS )
	return( 2 );

    // Outside circle bounds in x-direction
    if( (dx <= -_r || dx >= _r) )
	return( 0 );	// No crossing

    double dy = sqrt( _r*_r - dx*dx );
    double cy1 = _pc[1] + dy;
    double cy2 = _pc[1] - dy;

    if( y <= cy2 )
	// Under the circle -> no crossing
	return( 0 );

    double a1 = atan2( dy,dx );
    if( a1 < 0.0 )
	a1 += 2.0*M_PI;

    double a2 = atan2( -dy,dx );
    if( a2 < 0.0 )
	a2 += 2.0*M_PI;

#ifdef DEBUG_INSIDE_LOOP
    std::cout << "  a1 = " << a1 << "\n";
    std::cout << "  a2 = " << a2 << "\n";
#endif

    int c = 0;
    if( _ang1 < _ang2 ) {
	if( a1 > _ang1 && a1 < _ang2 && y > cy1 )
	    c = !c;
	if( a2 > _ang1 && a2 < _ang2 && y > cy2 )
	    c = !c;
    } else {
	if( (a1 < _ang2 || a1 > _ang1) && y > cy1 )
	    c = !c;
	if( (a2 < _ang2 || a2 > _ang1) && y > cy2 )
	    c = !c;
    }

    return( c );
}

void MyDXFArc::set_ang1( double ang1 ) 
{
    // Enforce between 0 and 2 pi
    _ang1 = ang1 - 2.0*M_PI*floor( ang1/(2.0*M_PI) );
}


void MyDXFArc::set_ang2( double ang2 )
{ 
    // Enforce between 0 and 2 pi
    _ang2 = ang2 - 2.0*M_PI*floor( ang2/(2.0*M_PI) );
}


void MyDXFArc::debug_print( std::ostream &os ) const
{
    os << "  ARC\n";
    MyDXFFile::debug_print_format( os, "pc", _pc );
    MyDXFFile::debug_print_format( os, "r", _r );
    MyDXFFile::debug_print_format( os, "ang1", _ang1 );
    MyDXFFile::debug_print_format( os, "ang2", _ang2 );
}


bool MyDXFArc::geom_same( const MyDXFArc &arc, double eps ) const
{
    return( norm2( _pc - arc._pc ) < eps &&
	    fabs( _r - arc._r ) < eps &&
	    norm2( start() - arc.start() ) < eps &&
	    norm2( end() - arc.end() ) < eps );
}


void MyDXFArc::set_start( const Vec3D &s )
{
    Vec3D e = end();
    set_center_point( s, e );
}


void MyDXFArc::set_end( const Vec3D &e )
{
    Vec3D s = start();
    set_center_point( s, e );
}


void MyDXFArc::set_center_and_ends( const Vec3D &c, const Vec3D &s, const Vec3D &e )
{
    _pc = c;
    _r = norm2(c-s);
    double dx = s[0]-c[0];
    double dy = s[1]-c[1];
    _ang1 = atan2( dy, dx );
    if( _ang1 < 0.0 )
	_ang1 += 2.0*M_PI;
    dx = e[0]-c[0];
    dy = e[1]-c[1];
    _ang2 = atan2( dy, dx );
    if( _ang2 < 0.0 )
	_ang2 += 2.0*M_PI;
}


void MyDXFArc::set_center_point( const Vec3D &s, const Vec3D &e )
{
    if( norm2( e-s ) < 2.0*_r ) {
	Vec3D u = 0.5*(e-s);
	Vec3D v( -u[1], u[0] );

	double A = v[0]*v[0] + v[1]*v[1];
	double B = 2.0*( v[0]*u[0] + v[1]*u[1] );
	double C = u[0]*u[0] + u[1]*u[1] - _r*_r;

	double t1, t2;
	uint32_t n = solve_quadratic( A, B, C, &t1, &t2 );
	if( n != 2 )
	    throw Error( ERROR_LOCATION, "Less than two roots found" );
	
	// solver gives solutions in ascending order and t2 = -t1
	_pc = s+u+t2*v;
    } else {
	_pc = 0.5*(s+e);
	_r = 0.5*norm2( e-s );
    }

    Vec3D dd = s-_pc;
    double dx = dd[0];
    double dy = dd[1];
    _ang1 = atan2( dy,dx );
    if( _ang1 < 0.0 )
	_ang1 += 2.0*M_PI;

    dd = e-_pc;
    dx = dd[0];
    dy = dd[1];
    _ang2 = atan2( dy,dx );
    if( _ang2 < 0.0 )
	_ang2 += 2.0*M_PI;
}





