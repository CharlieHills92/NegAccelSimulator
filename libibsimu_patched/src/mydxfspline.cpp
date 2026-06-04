/*! \file mydxfspline.hpp
 *  \brief DXF spline entity
 */

/* Copyright (c) 2011,2012,2014 Taneli Kalvas. All rights reserved.
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
#include "mydxfspline.hpp"


MyDXFSpline::MyDXFSpline( class MyDXFFile *dxf )
{
#ifdef MYDXF_DEBUG
    std::cout << "  Reading entity SPLINE\n";
#endif

    int32_t Nknot = 0;
    int32_t Ncont = 0;
    int32_t Nfit = 0;
    int32_t iknot = -1;
    int32_t icont = -1;
    int32_t ifit = -1;
    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with entity

	else if( dxf->group_get_code() == 12 )
	    _tangent0[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 22 )
	    _tangent0[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 32 )
	    _tangent0[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 13 )
	    _tangent1[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 23 )
	    _tangent1[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 33 )
	    _tangent1[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 210 )
	    _normal[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 220 )
	    _normal[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 230 )
	    _normal[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 42 )
	    _knot_tol = dxf->group_get_double();
	else if( dxf->group_get_code() == 43 )
	    _cont_tol = dxf->group_get_double();
	else if( dxf->group_get_code() == 44 )
	    _fit_tol = dxf->group_get_double();

	else if( dxf->group_get_code() == 70 )
	    _flags = dxf->group_get_int16();
	else if( dxf->group_get_code() == 71 )
	    _degree = dxf->group_get_int16();
	else if( dxf->group_get_code() == 72 ) {
	    if( iknot != -1 )
		throw Error( ERROR_LOCATION, "Number of knots specified after first "
			     "knot in SPLINE on line " + dxf->linec() );
	    Nknot = dxf->group_get_int16();
	    _knot.resize( Nknot );
	} else if( dxf->group_get_code() == 73 ) {
	    if( icont != -1 )
		throw Error( ERROR_LOCATION, "Number of control points specified after first "
			     "control point in SPLINE on line " + dxf->linec() );
	    Ncont = dxf->group_get_int16();
	    _cont.resize( Ncont );

	} else if( dxf->group_get_code() == 74 ) {
	    if( ifit != -1 )
		throw Error( ERROR_LOCATION, "Number of fit points specified after first "
			     "fir point in SPLINE on line " + dxf->linec() );
	    Nfit = dxf->group_get_int16();
	    _fit.resize( Nfit );
	    
	} else if( dxf->group_get_code() == 40 ) {
	    iknot++;
	    if( iknot == Nknot )
		throw Error( ERROR_LOCATION, "Incorrect amount of knots in SPLINE"
			     " on line " + dxf->linec() );
	    _knot[iknot] = dxf->group_get_double();

	} else if( dxf->group_get_code() == 10 ) {
	    icont++;
	    if( icont == Ncont )
		throw Error( ERROR_LOCATION, "Incorrect amount of control points in SPLINE"
			     " on line " + dxf->linec() );
	    _cont[icont][0] = dxf->group_get_double();
	} else if( dxf->group_get_code() == 20 )
	    _cont[icont][1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 )
	    _cont[icont][2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 11 ) {
	    ifit++;
	    if( ifit == Nfit )
		throw Error( ERROR_LOCATION, "Incorrect amount of fit points in SPLINE"
			     " on line " + dxf->linec() );
	    _fit[ifit][0] = dxf->group_get_double();
	} else if( dxf->group_get_code() == 21 )
	    _fit[ifit][1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 31 )
	    _fit[ifit][2] = dxf->group_get_double();

	else
	    process_group( dxf );
    }

    // Check that the knot, control and fit vectors got filled
    if( iknot != Nknot-1 )
	throw( Error( ERROR_LOCATION, "Incorrect amount of knots in SPLINE"
		      " on line " + dxf->linec() ) );
    if( icont != Ncont-1 )
	throw( Error( ERROR_LOCATION, "Incorrect amount of control points in SPLINE"
		      " on line " + dxf->linec() ) );
    if( ifit != Nfit-1 )
	throw( Error( ERROR_LOCATION, "Incorrect amount of fit points in SPLINE"
		      " on line " + dxf->linec() ) );

    // Check that sizes match
    if( _knot.size() != _cont.size() + _degree + 1 )
	throw( Error( ERROR_LOCATION, "Incorrect number knots (" + to_string(_knot.size()) + ") vs."
		      " control points (" + to_string(_cont.size()) + ") in SPLINE of degree " + 
		      to_string(_degree) + " on line " + to_string(dxf->linec()) ) );

    if( _cont.size() < (size_t)(_degree+1) )
	throw( Error( ERROR_LOCATION, "Too few control points (" + to_string(_cont.size()) + ") for"
		      " SPLINE of degree " + to_string(_degree) ) );

    check_knots();

    // QCAD makes all closed splines as periodic. Haven't had a
    // possibility test non-periodic closed splines
    if( periodic() ) {

	//std::cout << "ORIG:\n";
	//for( size_t a = 0; a < _knot.size(); a++ )
	//    std::cout << "knot[" << a << "] = " << _knot[a] << "\n";
	//for( size_t a = 0; a < _cont.size(); a++ )
	//    std::cout << "cont[" << a << "] = " << _cont[a] << "\n";

	// Add _degree amount of points
	size_t addp = _degree;

	size_t a;
	for( a = 0; a < addp/2; a++ )
	    add_to_end( _cont[a] );

	//std::cout << "ADDED END:\n";
	//for( size_t a = 0; a < _knot.size(); a++ )
	//    std::cout << "knot[" << a << "] = " << _knot[a] << "\n";
	//for( size_t a = 0; a < _cont.size(); a++ )
	//    std::cout << "cont[" << a << "] = " << _cont[a] << "\n";

	//size_t added = a;
	//std::cout << "added = " << added << "\n";
	for( ; a < _degree; a++ )
	    add_to_start( _cont[_cont.size()-1-a] );

	//std::cout << "ADDED START:\n";
	//for( size_t a = 0; a < _knot.size(); a++ )
	//    std::cout << "knot[" << a << "] = " << _knot[a] << "\n";
	//for( size_t a = 0; a < _cont.size(); a++ )
	//    std::cout << "cont[" << a << "] = " << _cont[a] << "\n";

	make_cyclic();

	//std::cout << "MADE CYCLIC:\n";
	//for( size_t a = 0; a < _knot.size(); a++ )
	//    std::cout << "knot[" << a << "] = " << _knot[a] << "\n";
	//for( size_t a = 0; a < _cont.size(); a++ )
	//    std::cout << "cont[" << a << "] = " << _cont[a] << "\n";
    }

    if( rational() )
	throw( ErrorUnimplemented( ERROR_LOCATION, "Rational spline unimplemented" ) );

    // Make polyline approximation of spline
    build_polyline();
    
#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


void MyDXFSpline::make_cyclic( void )
{
    // Start
    int32_t ref1 = _degree;
    double w1 = _knot[ref1+1] - _knot[ref1];
    for( size_t a = 1; a < _degree; a++ )
	_knot[ref1-a] -= a*w1;

    // End
    int32_t ref2 = _knot.size()-1-_degree;
    double w2 = _knot[ref2] - _knot[ref2-1];
    for( size_t a = 1; a < _degree; a++ ) 
	_knot[ref2+a] += a*w2;

}


// Add a new control point to end of spline
void MyDXFSpline::add_to_end( const Vec3D &x )
{
    //std::cout << "add_to_end( x = " << x << " )\n";

    // Modify knots
    int32_t ref = _knot.size()-_degree-1;
    double weight = _knot[ref] - _knot[ref-1];
    //std::cout << "weight = " << weight << "\n";
    for( size_t a = ref+1; a < _knot.size(); a++ )
	_knot[a] += weight;
    _knot.push_back( _knot.back() );

    // Add control point to end
    _cont.push_back( x );
}


// Add a new control point to start of spline
void MyDXFSpline::add_to_start( const Vec3D &x )
{
    //std::cout << "add_to_start( x = " << x << " )\n";

    // Modify knots
    int32_t ref = _degree;
    double weight = _knot[ref+1] - _knot[ref];
    //std::cout << "weight = " << weight << "\n";
    _knot.push_back( _knot.back() );
    int32_t a;
    for( a = _knot.size()-2; a > ref; a-- )
	_knot[a] = _knot[a-1];
    for( ; a >= 0; a-- )
	_knot[a] -= weight;

    // Add control point to start
    _cont.insert( _cont.begin(), x );
}


void MyDXFSpline::check_knots( void ) const
{
    for( size_t a = 0; a < _knot.size()-1; a++ ) {
	if( _knot[a] > _knot[a+1] ) {
	    throw( Error( ERROR_LOCATION, "Knot vector not increasing" ) );
	}
    }
}


void MyDXFSpline::explode( class MyDXFEntities *ent, MyDXFFile *dxf, const Transformation *t ) const
{
    MyDXFSpline *spline = new MyDXFSpline( *this );

    // Transform control and fit points
    for( uint32_t a = 0; a < spline->_cont.size(); a++ )
	spline->_cont[a] = t->transform_point( spline->_cont[a] );
    for( uint32_t a = 0; a < spline->_fit.size(); a++ )
	spline->_fit[a] = t->transform_point( spline->_fit[a] );

    // Add to entities
    ent->add_entity( spline );
}


void MyDXFSpline::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "SPLINE" );
    write_common( dxf, ostr );

    dxf->write_group( 210, _normal[0] );
    dxf->write_group( 220, _normal[1] );
    dxf->write_group( 230, _normal[2] );

    dxf->write_group( 70, _flags );
    dxf->write_group( 71, (int16_t)_degree );
    dxf->write_group( 72, (int16_t)_knot.size() );
    dxf->write_group( 73, (int16_t)_cont.size() );
    dxf->write_group( 74, (int16_t)_fit.size() );

    dxf->write_group( 42, _knot_tol );
    dxf->write_group( 43, _cont_tol );
    dxf->write_group( 44, _fit_tol );

    dxf->write_group( 12, _tangent0[0] );
    dxf->write_group( 22, _tangent0[1] );
    dxf->write_group( 32, _tangent0[2] );

    dxf->write_group( 13, _tangent1[0] );
    dxf->write_group( 23, _tangent1[1] );
    dxf->write_group( 33, _tangent1[2] );

    for( uint32_t i = 0; i < _knot.size(); i++ ) {
	dxf->write_group( 40, _knot[i] );
    }

    for( uint32_t i = 0; i < _cont.size(); i++ ) {
	dxf->write_group( 10, _cont[i][0] );
	dxf->write_group( 20, _cont[i][1] );
	dxf->write_group( 30, _cont[i][2] );
    }

    for( uint32_t i = 0; i < _fit.size(); i++ ) {
	dxf->write_group( 11, _fit[i][0] );
	dxf->write_group( 21, _fit[i][1] );
	dxf->write_group( 31, _fit[i][2] );
    }
}


Vec3D MyDXFSpline::start( void ) const
{
    return( _polyline[0] );
}


Vec3D MyDXFSpline::end( void ) const
{
    return( _polyline.back() );
}


void MyDXFSpline::set_start( const Vec3D &s )
{
    _polyline[0] = s;
}


void MyDXFSpline::set_end( const Vec3D &s )
{
    _polyline.back() = s;
}


int MyDXFSpline::ray_cross( double x, double y ) const
{
#ifdef MYDXF_DEBUG
    std::cout << "ray_cross( x = " << x << ", y = " << y << " )\n";
    for( uint32_t a = 0; a < _polyline.size(); a++ )
	std::cout << "  p[" << a << "] = " << _polyline[a] << "\n";
    std::cout << "\n";
#endif

    Vec3D p1 = _polyline[0];

    // Uncertainty region around p1
    if( fabs( x - p1[0] ) < MYDXF_PERT_EPS && y >= p1[1]-MYDXF_PERT_EPS )
	return( 2 );

    // Go through all lines in polyline representation
    int c = 0;
    for( uint32_t a = 1; a < _polyline.size(); a++ ) {
	Vec3D p2 = _polyline[a];

#ifdef MYDXF_DEBUG
	std::cout << "  Test against:\n";
	std::cout << "  p1 = " << p1 << "\n";
	std::cout << "  p2 = " << p2 << "\n";
#endif

	// Uncertainty region around p2
	if( fabs( x - p2[0] ) < MYDXF_PERT_EPS && y >= p2[1]-MYDXF_PERT_EPS )
	    return( 2 );

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


bool MyDXFSpline::geom_same( const MyDXFSpline &spline, double eps ) const
{
    return( false );
}


void MyDXFSpline::debug_print( std::ostream &os ) const
{
    os << "  SPLINE\n";
    MyDXFFile::debug_print_format( os, "flags", _flags );
    MyDXFFile::debug_print_format( os, "degree", _degree );
    MyDXFFile::debug_print_format( os, "knot_tol", _knot_tol );
    MyDXFFile::debug_print_format( os, "cont_tol", _cont_tol );
    MyDXFFile::debug_print_format( os, "fit_tol", _fit_tol );
    MyDXFFile::debug_print_format( os, "normal", _normal );
    MyDXFFile::debug_print_format( os, "tangent0", _tangent0 );
    MyDXFFile::debug_print_format( os, "tangent1", _tangent1 );

    for( uint32_t a = 0; a < _knot.size(); a++ )
	MyDXFFile::debug_print_format( os, "knot", _knot[a] );

    for( uint32_t a = 0; a < _cont.size(); a++ )
	MyDXFFile::debug_print_format( os, "cont", _cont[a] );

    for( uint32_t a = 0; a < _fit.size(); a++ )
	MyDXFFile::debug_print_format( os, "fit", _fit[a] );
}


// Return a point from the spline with parameter t in [0,1]
Vec3D MyDXFSpline::point( double t ) const
{
    //std::cout << "point( t = " << t << " )\n";

    // Scale parameter to range
    t = t*_knot[_knot.size()-1] + (1.0-t)*_knot[0];

    // Choose knot interval, where t is located
    uint32_t i;
    for( i = _degree; i < _knot.size()-2; i++ ) {
	if( t <= _knot[i+1] )
	    return( deboor( t, i, _degree ) );
    }
    return( deboor( t, _knot.size()-1, _degree ) );
}


// Recursive subroutine for point()
// Implements deBoor's algorithm
Vec3D MyDXFSpline::deboor( double t, int i, int k ) const
{
    /*
    for( size_t a = k; a < _degree; a++ )
	std::cout << "  ";
    std::cout << "deboor( t = " << t << ", i = " << i << ", k = " << k <<  " )\n";

    Vec3D ret;
    if( k == 0 ) {
	ret = _cont[i];
    } else {
	std::cout << "i+_degree+1-k = " <<  i+_degree+1-k << "\n";
	std::cout << "_knot[i+_degree+1-k] = " <<  _knot[i+_degree+1-k] << "\n";
	std::cout << "i = " <<  i << "\n";
	std::cout << "_knot[i] = " <<  _knot[i] << "\n";
	double diff = _knot[i+_degree+1-k] - _knot[i];
	double alpha;
	if( diff == 0.0 )
	    alpha = 0.0;
	else
	    alpha = ( t - _knot[i] ) / diff;
	for( size_t a = k; a < _degree+1; a++ )
	    std::cout << "  ";
	std::cout << "1.0-alpha = " << 1.0-alpha << "\n";
	ret = deboor( t, i-1, k-1 ) * (1.0-alpha);
	for( size_t a = k; a < _degree+1; a++ )
	    std::cout << "  ";
	std::cout << "alpha = " << alpha << "\n";
	ret += deboor( t, i,   k-1 ) * alpha;
    }

    for( size_t a = k; a < _degree; a++ )
	std::cout << "  ";
    std::cout << "ret = " << ret << "\n";

    return( ret );
    */

    if( k == 0 ) {
	return( _cont[i] );
    } else {
	double diff = _knot[i+_degree+1-k] - _knot[i];
	double alpha;
	if( diff == 0.0 )
	    alpha = 0.0;
	else
	    alpha = ( t - _knot[i] ) / diff;
	return( deboor( t, i-1, k-1 ) * (1.0-alpha) + 
		deboor( t, i,   k-1 ) * alpha );
    }
}


void MyDXFSpline::plot_knot_points( cairo_t *cairo, const Transformation *t ) const
{
    //std::cout << "plot_knot_points()\n";

    cairo_save( cairo );

    cairo_set_source_rgb( cairo, 1, 0, 0 );

    Vec4D x;
    const double r = 3.0;
    for( uint32_t a = 0; a < _knot.size(); a++ ) {

	double k = _knot[a];
	double kk = (k-_knot[0]) / (_knot.back()-_knot[0]);
	x = t->transform( point( kk ) );	
	//x.homogenize();
	//std::cout << "knot[" << a << "] = " << _knot[a] << ": x = " << x << "\n";
	cairo_move_to( cairo, x[0]+r, x[1] );
	cairo_arc( cairo, x[0], x[1], r, 0.0, 2.0*M_PI );
	cairo_fill( cairo );
    }

    cairo_restore( cairo );
}


void MyDXFSpline::plot_polyline_points( cairo_t *cairo, const Transformation *t ) const
{
    //std::cout << "plot_polyline_points()\n";

    cairo_save( cairo );

    cairo_set_source_rgb( cairo, 0, 0, 1 );

    // Plot points
    /*
    Vec3D x;
    const double r = 3.0;
    for( uint32_t a = 0; a < _polyline.size(); a++ ) {

	x = t->transform_point( _polyline[a] );
	x.homogenize();
	std::cout << "polyline[" << a << "] = " << _polyline[a] << ": x = " << x << "\n";
	cairo_move_to( cairo, x[0]+r, x[1] );
	cairo_arc( cairo, x[0], x[1], r, 0.0, 2.0*M_PI );
	cairo_fill( cairo );
    }
    */

    // Plot lines
    Vec3D x;
    for( uint32_t a = 0; a < _polyline.size(); a++ ) {

	x = t->transform_point( _polyline[a] );
	//x.homogenize();
	//std::cout << "polyline[" << a << "] = " << _polyline[a] << ": x = " << x << "\n";
	if( a == 0 )
	    cairo_move_to( cairo, x[0], x[1] );
	else
	    cairo_line_to( cairo, x[0], x[1] );
    }
    cairo_stroke( cairo );

    cairo_restore( cairo );
}


void MyDXFSpline::plot_control_points( cairo_t *cairo, const Transformation *t ) const
{
    cairo_save( cairo );

    cairo_set_source_rgb( cairo, 0, 1, 0 );

    Vec3D x;
    const double r = 3.0;
    for( uint32_t a = 0; a < _cont.size(); a++ ) {
	x = t->transform_point( _cont[a] );
	//x.homogenize();
	cairo_move_to( cairo, x[0]+r, x[1] );
	cairo_arc( cairo, x[0], x[1], r, 0.0, 2.0*M_PI );
	cairo_fill( cairo );
    }

    cairo_restore( cairo );
}


void MyDXFSpline::build_polyline( void )
{
    _polyline.clear();

    double kstart = _knot[0];
    double kend = _knot.back();

    // Build a polyline with _degree-1 midpoints and endpoints at each knot
    for( uint32_t a = 0; a < _knot.size()-1; a++ ) {

	double k1 = _knot[a];
	double k2 = _knot[a+1];

	// Avoid using connection nodes with periodic spline
	if( k2 < kstart || k1 < kstart )
	    continue;

	// Make ptc many points, +1 for last
	uint32_t ptc = _degree*2;
	uint32_t ptl = ptc;
	if( k2 == kend )
	    ptl += 1;

	for( uint32_t b = 0; b < ptl; b++ ) {
	    double k = k1+(k2-k1)*b/ptc;
	    double kk = (k-kstart) / (kend-_knot[0]);
	    _polyline.push_back( point( kk ) );
	}

	// Avoid using connection nodes with periodic spline
	if( k2 == kend )
	    break;
    }
}


void MyDXFSpline::plot( const class MyDXFFile *dxf, cairo_t *cairo,
			const Transformation *t, const double range[4] ) const
{
    Vec3D x;

    //plot_control_points( cairo, t );
    //plot_knot_points( cairo, t );
    //plot_polyline_points( cairo, t );
    //return;

    // Plot using deBoor's algorithm
    x = t->transform_point( point( 0.0 ) );
    //x.homogenize();
    cairo_move_to( cairo, x[0], x[1] );
    //std::cout << "x = ( " << x[0] << ", " << x[1] << " )\n";

    for( uint32_t a = 1; a < 101; a++ ) {
	double tt = a/100.0;
	x = t->transform_point( point( tt ) );
	//x.homogenize();
	cairo_line_to( cairo, x[0], x[1] );
	//std::cout << "x = ( " << x[0] << ", " << x[1] << " )\n";
    }

    /*
    if( _degree == 1 ) {

	// Must have at least two points to draw
	if( _cont.size() <= 2 )
	    return;

	Vec3D x = t->transform_point( _cont[0] );
	x.homogenize();
	cairo_move_to( cairo, x[0], x[1] );
	for( uint32_t a = 1; a < _cont.size(); a++ ) {
	    x = t->transform_point( _cont[a] );
	    x.homogenize();
	    cairo_line_to( cairo, x[0], x[1] );
	}
	if( closed() ) {
	    x = t->transform_point( _cont[0] );
	    x.homogenize();
	    cairo_line_to( cairo, x[0], x[1] );
	}
    } 
    */

    cairo_stroke( cairo );
}


void MyDXFSpline::get_bbox( Vec3D &min, Vec3D &max, 
			    const class MyDXFFile *dxf, const Transformation *t ) const
{
    // Bounding box estimate for the curve is the bounding box for all control points
    min = Vec3D( std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    max = Vec3D( -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity() );

    for( uint32_t a = 0; a < _cont.size(); a++ ) {
	Vec3D x = t->transform_point( _cont[a] );
	//x.homogenize();
	bbox_ppoint( min, max, x );
    }
}


void MyDXFSpline::scale( class MyDXFFile *dxf, double s )
{
    // Scale control, fit and polyline points
    for( uint32_t a = 0; a < _cont.size(); a++ )
	_cont[a] *= s;
    for( uint32_t a = 0; a < _fit.size(); a++ )
	_fit[a] *= s;
    for( uint32_t a = 0; a < _polyline.size(); a++ )
	_polyline[a] *= s;
}


void MyDXFSpline::translate( class MyDXFFile *dxf, const Vec3D &dx )
{
    // Translate control, fit and polyline points
    for( uint32_t a = 0; a < _cont.size(); a++ )
	_cont[a] += dx;
    for( uint32_t a = 0; a < _fit.size(); a++ )
	_fit[a] += dx;
    for( uint32_t a = 0; a < _polyline.size(); a++ )
	_polyline[a] += dx;
}


void MyDXFSpline::rotate_z( class MyDXFFile *dxf, double a )
{
    Transformation t = Transformation::rotation_z( a );

    // Transform control and fit points
    for( uint32_t a = 0; a < _cont.size(); a++ )
	_cont[a] = t.transform_point( _cont[a] );
    for( uint32_t a = 0; a < _fit.size(); a++ )
	_fit[a] = t.transform_point( _fit[a] );
}



