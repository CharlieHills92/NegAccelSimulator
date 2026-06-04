/*! \file ruler.cpp
 *  \brief Rulers for plot frames
 */

/* Copyright (c) 2005-2010,2012,2013 Taneli Kalvas. All rights reserved.
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

#include <cmath>
#include "ruler.hpp"
#include "compmath.hpp"


//#define DEBUG_RULER 1


Ruler::Ruler()
    : _color((Vec3D(0,0,0))), _ticlen_in(5.0), _ticlen_out(5.0), _labelspace(5.0),
      _fontsize(12.0), _label_enabled(true), _indir(true), _cind(0)
{
    _range[0] = 0.0;
    _range[1] = 1.0;

    _autorange[0] = true;
    _autorange[1] = true;

    _endpt[0] = 0.0;
    _endpt[1] = 0.0;
    _endpt[2] = 100.0;
    _endpt[3] = 0.0;
}


Ruler::Ruler( int cind )
    : _color((Vec3D(0,0,0))), _ticlen_in(5.0), _ticlen_out(5.0), _labelspace(5.0),
      _fontsize(12.0), _label_enabled(true), _indir(true), _cind(cind)
{
    _range[0] = 0.0;
    _range[1] = 1.0;

    _autorange[0] = true;
    _autorange[1] = true;

    if( cind == 0 ) {
	_endpt[0] = 0.0;
	_endpt[1] = 0.0;
	_endpt[2] = 100.0;
	_endpt[3] = 0.0;
    } else {
	_endpt[0] = 0.0;
	_endpt[1] = 0.0;
	_endpt[2] = 0.0;
	_endpt[3] = 100.0;
    }
}


Ruler::Ruler( const Ruler &ruler )
    : _color(ruler._color), _ticlen_in(ruler._ticlen_in), _ticlen_out(ruler._ticlen_out),
      _labelspace(ruler._labelspace), _fontsize(ruler._fontsize), _axislabel(ruler._axislabel), 
      _label_enabled(ruler._label_enabled), 
      _indir(ruler._indir), _cind(ruler._cind)
{
    _range[0] = ruler._range[0];
    _range[1] = ruler._range[1];

    _autorange[0] = ruler._autorange[0];
    _autorange[1] = ruler._autorange[1];

    _endpt[0] = ruler._endpt[0];
    _endpt[1] = ruler._endpt[1];
    _endpt[2] = ruler._endpt[2];
    _endpt[3] = ruler._endpt[3];

    for( size_t a = 0; a < ruler._tic.size(); a++ )
	_tic.push_back( ruler._tic[a] );
}


Ruler &Ruler::operator=( const Ruler &ruler )
{
    if( this != &ruler ) {
	_color = ruler._color;
	_ticlen_in = ruler._ticlen_in;
	_ticlen_out = ruler._ticlen_out;
	_labelspace = ruler._labelspace;
	_fontsize = ruler._fontsize;
	_axislabel = ruler._axislabel;
	_label_enabled = ruler._label_enabled;
	_indir = ruler._indir;
	_cind = ruler._cind;

	_range[0] = ruler._range[0];
	_range[1] = ruler._range[1];
	
	_autorange[0] = ruler._autorange[0];
	_autorange[1] = ruler._autorange[1];
	
	_endpt[0] = ruler._endpt[0];
	_endpt[1] = ruler._endpt[1];
	_endpt[2] = ruler._endpt[2];
	_endpt[3] = ruler._endpt[3];
	
	_tic.clear();
	for( size_t a = 0; a < ruler._tic.size(); a++ )
	    _tic.push_back( ruler._tic[a] );
    }

    return( *this );
}


void Ruler::copy_tics( const Ruler &ruler )
{
    // Copy range
    _range[0] = ruler._range[0];
    _range[1] = ruler._range[1];

    // Copy tics
    _tic.clear();
    for( size_t a = 0; a < ruler._tic.size(); a++ )
	_tic.push_back( ruler._tic[a] );
}


void Ruler::set_font_size( double size )
{
    _ticlen_in = size*5.0/12.0;
    _ticlen_out = size*5.0/12.0;
    _labelspace = size*5.0/12.0;
    _fontsize = size;
    _axislabel.set_font_size( size );
}


void Ruler::set_color( const Vec3D &color )
{
    _color = color;
}


void Ruler::set_ticlen( double inlen, double outlen )
{
    _ticlen_in  = inlen > 0.0 ? inlen : 0.0;
    _ticlen_out = outlen > 0.0 ? outlen : 0.0;
}


void Ruler::set_autorange( bool autorange_min, bool autorange_max )
{
    _autorange[0] = autorange_min;
    _autorange[1] = autorange_max;
}


void Ruler::get_autorange( bool &autorange_min, bool &autorange_max ) const
{
    autorange_min = _autorange[0];
    autorange_max = _autorange[1];
}


void Ruler::set_ranges( double min, double max )
{
#ifdef DEBUG_RULER
    std::cout << "Ruler::set_ranges( " << min << ", " << max << ")\n";
#endif

    if( comp_isinf(min) || comp_isnan(min) ) {
	if( comp_isinf(max) || comp_isnan(max) ) {
	    _range[0] = 0.0;
	    _range[1] = 1.0;
	} else {
	    _range[0] = max-1.0;
	    _range[1] = max;
	}    
    } else if( comp_isinf(max) || comp_isnan(max) ) {
	_range[0] = min;
	_range[1] = min+1.0;
    } else {
	_range[0] = min;
	_range[1] = max;
    }
}


void Ruler::get_ranges( double &min, double &max ) const
{
    min = _range[0];
    max = _range[1];
}


void Ruler::set_endpoints( double x1, double y1, double x2, double y2 )
{
    _endpt[0] = x1;
    _endpt[1] = y1;
    _endpt[2] = x2;
    _endpt[3] = y2;
}


void Ruler::set_axis_label( const std::string &label )
{
    _axislabel.set_text( label );
}


void Ruler::enable_labels( bool enable )
{
    _label_enabled = enable;
}


void Ruler::set_indir( bool ccw )
{
    _indir = ccw;
}


void Ruler::set_coord_index( int cind )
{
    _cind = cind;
}


// Test for bounding box crash. Return true if crash.
// bbox = {xmin, ymin, xmax, ymax}
bool Ruler::tic_labels_bbox_crash_x( const double bbox[4], const double obbox[4] )
{
    if( (bbox[0] <= obbox[2] && bbox[0] >= obbox[0]) ||
	(bbox[2] <= obbox[2] && bbox[2] >= obbox[0]) )
	return( true );

    return( false );
}


// Test for bounding box crash. Return true if crash.
// bbox = {xmin, ymin, xmax, ymax}
bool Ruler::tic_labels_bbox_crash_y( const double bbox[4], const double obbox[4] )
{
    if( (bbox[1] <= obbox[3] && bbox[1] >= obbox[1]) ||
	(bbox[3] <= obbox[3] && bbox[3] >= obbox[1]) )
	return( true );

    return( false );
}


// Add a tic and set its properties including location calculated
// using cm. Make bbox crash test if ruler_tic_bbox_test is
// enabled. Return true if test reports a crash. The tic is still
// built as a whole regardless of the reports result. Also maximum
// width for y-rulers and maximum height for x-rulers is built in
// maxsize. The bounding box of the tic is returned in obbox if it is
// calculated.
//
bool Ruler::add_tic( double x, cairo_t *cairo, const Coordmapper1D &cm, 
		     bool ruler_tic_bbox_test, double &maxsize, double obbox[4] ) 
{
    // Indir sign
    double sign;
    if( _indir ) sign = 1.0;
    else sign = -1.0;

    // Prevent small non-zero numbers because of rounding errors
    if( fabs(x/(_range[1]-_range[0])) < 1e-6 ) x = 0;
    char str[128];
    snprintf( str, 128, "%g", x );

    // Add tic
    _tic.push_back( Tic( x, str ) );
    size_t n = _tic.size()-1;

    // Set tic parameters
    _tic[n]._label.set_font_size( _fontsize );
    if( _cind == 0 ) {

	// xruler
	if( _indir == true ) {
	    // x1
	    _tic[n]._label.set_alignment( 0.5, 1.0, true );
	} else {
	    // x2
	    _tic[n]._label.set_alignment( 0.5, 0.0, true );
	}

	cm.transform( x );
	double y = _endpt[1]+(_endpt[3]-_endpt[1])*(x-_endpt[0])/(_endpt[2]-_endpt[0]);
	_tic[n]._label.set_location( x, y+sign*(_ticlen_out+_labelspace) );

	if( _label_enabled ) {
	    double bbox[4];
	    _tic[n]._label.get_bbox( cairo, bbox );
	    if( bbox[3]-bbox[1] > maxsize )
		maxsize = bbox[3]-bbox[1];
	    if( ruler_tic_bbox_test && n > 1 ) {
		// Make bbox crash test
		if( tic_labels_bbox_crash_x( bbox, obbox ) )
		    return( true );
	    }
	    obbox[0] = bbox[0];
	    obbox[1] = bbox[1];
	    obbox[2] = bbox[2];
	    obbox[3] = bbox[3];
	}
	
    } else /* _cind == 1 */ {

	// yruler
	if( _indir == false ) {
	    // y1
	    _tic[n]._label.set_alignment( 1.0, 0.5, true );
	} else {
	    // y2
	    _tic[n]._label.set_alignment( 0.0, 0.5, true );
	}

	cm.transform( x );
	double y = _endpt[0]+(_endpt[2]-_endpt[0])*(x-_endpt[1])/(_endpt[3]-_endpt[1]);
	_tic[n]._label.set_location( y+sign*(_ticlen_out+_labelspace), x );
	    
	if( _label_enabled ) {
	    double bbox[4];
	    _tic[n]._label.get_bbox( cairo, bbox );
	    if( bbox[2]-bbox[0] > maxsize )
		maxsize = bbox[2]-bbox[0];
	    if( ruler_tic_bbox_test && n > 1 ) {
		// Make bbox crash test
		if( tic_labels_bbox_crash_y( bbox, obbox ) )
		    return( true );
	    }
	    obbox[0] = bbox[0];
	    obbox[1] = bbox[1];
	    obbox[2] = bbox[2];
	    obbox[3] = bbox[3];
	}
    }

    return( false );
}


/* Calculate pow( 10, x ) accurately
 */
double accpow10( int x )
{
    if( x > 0 ) {

	double res = 1.0;
	for( int a = 0; a < x; a++ )
	    res *= 10.0;
	return( res );

    } else if( x < 0 ) {

	double res = 1.0;
	for( int a = 0; a < -x; a++ )
	    res *= 10.0;
	return( 1.0/res );

    }

    return( 1.0 );
}


void Ruler::calculate( cairo_t *cairo, Coordmapper1D &cm, bool ruler_tic_bbox_test ) 
{
    int num = 8;
    double step;
    double big  = fabs(_range[1]-_range[0]) / num;
    //double bigl = pow( 10, floor( log10( fabs(_range[1]-_range[0]) / num ) ) );
    double bigl = accpow10( (int)floor( log10( fabs(_range[1]-_range[0]) / num ) ) );
    int fac;
    if( bigl >= big ) 
        fac = 1;
    else if( 2*bigl >= big ) 
        fac = 2;
    else if( 5*bigl >= big )
        fac = 5;
    else 
        fac = 10;

    // Step direction
    bool stepdir = true;
    if( _range[1] < _range[0] )
	stepdir = false;

#ifdef DEBUG_RULER
    std::cout << "Ruler::calculate()\n";
    std::cout << "  min = " << _range[0] << "\n";
    std::cout << "  max = " << _range[1] << "\n";
    std::cout << "  big = " << big << "\n";
    std::cout << "  bigl = " << bigl << "\n";
    std::cout << "  fac = " << fac << "\n";
    std::cout << "  stepdir = " << stepdir << "\n";
#endif

    double maxsize = 0.0;
    double orig_range[2] = { _range[0], _range[1] };
    for( size_t a = 0; a < 2; a++ ) {

#ifdef DEBUG_RULER
	std::cout << "Ruler::calculate(): round = " << a << "\n";
#endif

	// Calculate step
	if( stepdir ) step = fac*bigl;
	else step = -fac*bigl;

	// Disable bbox tests on last round
	if( a == 1 )
	    ruler_tic_bbox_test = false;

	// Restore original range
	if( a != 0 ) {
	    _range[0] = orig_range[0];
	    _range[1] = orig_range[1];
	}

	// Do autoranges according to tic spacing
	double x = 0.0;
	if( stepdir )
	    // Allow 1% of step size error for inaccurate arithmetic
	    x = step*floor(_range[0]/step-0.01); 
	else
	    x = step*ceil(_range[0]/step+0.01);
	if( _autorange[0] )
	    _range[0] = x;
#ifdef DEBUG_RULER
	std::cout << "Ruler::calculate(): starting x = " << x << "\n";
	std::cout << "Ruler::calculate(): comp1 = " << _range[1]+0.01*step << "\n";
	std::cout << "Ruler::calculate(): comp2 = " << _range[1]-0.01*step << "\n";
#endif
	if( _autorange[1] ) {
	    int i = 0;
	    while( (stepdir && x <= _range[1]+0.01*step) || (!stepdir && x >= _range[1]-0.01*step) ) {
		i++;
		x = _range[0] + i*step;
#ifdef DEBUG_RULER
		std::cout << "Ruler::calculate(): trying x = " << x << "\n";
#endif

	    }
	    _range[1] = x;
	}

#ifdef DEBUG_RULER
	std::cout << "  step = " << step << "\n";
	std::cout << "  autoranged min = " << _range[0] << "\n";
	std::cout << "  autoranged max = " << _range[1] << "\n";
#endif

	// Calculate coordmapper with current ranges
	if( _cind == 0 ) {
	    double xx = (_endpt[2]-_endpt[0]) / (_range[1]-_range[0]);
	    double x0 = _endpt[0] - _range[0]*xx;
	    cm.set_transformation( xx, x0 );
	} else {
	    double yy = (_endpt[3]-_endpt[1]) / (_range[1]-_range[0]);
	    double y0 = _endpt[1] - _range[0]*yy;
	    cm.set_transformation( yy, y0 );
	}

	// Make tics and make bbox test if required
	maxsize = 0.0;
	double obbox[4];
	bool bbox_test_crash = false;
	_tic.clear();
	double startx = step*floor(_range[0]/step);
	x = startx;
	int i = 0;
	while( (stepdir && x < _range[0]-0.01*step) || (!stepdir && x > _range[0]+0.01*step) ) {
	    i++;
	    x = startx + i*step;
	}
	while( (stepdir && x <= _range[1]+0.01*step) || (!stepdir && x >= _range[1]-0.01*step) ) {

	    // Make a tic and set parameters and make bounding box
	    // test enabled. Returns true if bbox test reports a crash.
	    if( add_tic( x, cairo, cm, ruler_tic_bbox_test, maxsize, obbox ) ) {
		bbox_test_crash = true;
		break;
	    }

	    i++;
	    x = startx + i*step;
	}
	if( bbox_test_crash ) {
	    // Increase fac according to 1, 2, 5, 10 sequence
	    int base = 1;
	    while( fac >= 10 ) {
		fac /= 10;
		base *= 10;
	    }
	    if( fac == 1 )
		fac = 2*base;
	    else if( fac == 2 )
		fac = 5*base;
	    else
		fac = 10*base;
	    
	    continue;
	}

	// Tics accepted 
	break;
    }

    // Indir sign
    double sign;
    if( _indir ) sign = 1.0;
    else sign = -1.0;

    // Set axis label properties
    if( _label_enabled && _axislabel.get_text() != "" ) {
	if( _cind == 0 ) {
	    // xruler
	    _axislabel.set_location( 0.5*(_endpt[0]+_endpt[2]), 
				     0.5*(_endpt[1]+_endpt[3])+sign*(_ticlen_out+2.0*_labelspace+maxsize) );
	    if( _indir )
		_axislabel.set_alignment( 0.5, 1.0 );
	    else
		_axislabel.set_alignment( 0.5, 0.0 );
	    _axislabel.set_rotation( 0.0 );
	} else {
	    // yruler
	    _axislabel.set_location( 0.5*(_endpt[0]+_endpt[2])+sign*(_ticlen_out+2.0*_labelspace+maxsize), 
				     0.5*(_endpt[1]+_endpt[3]) );
	    if( _indir )
		_axislabel.set_alignment( 0.5, 1.0 );
	    else
		_axislabel.set_alignment( 0.5, 0.0 );
	    _axislabel.set_rotation( 0.5*M_PI );
	}
    }
}


void Ruler::draw( cairo_t *cairo, Coordmapper1D &cm, bool recalculate )
{
    if( recalculate )
	calculate( cairo, cm, true );

    // Set cairo parameters
    cairo_set_source_rgba( cairo, _color[0], _color[1], _color[2], 1.0 );
    cairo_set_line_width( cairo, 1.0 );

    // Draw baseline
    cairo_move_to( cairo, _endpt[0], _endpt[1] );
    cairo_line_to( cairo, _endpt[2], _endpt[3] );
    cairo_stroke( cairo );

    // Draw tics
    double sign;
    if( _indir ) sign = 1.0;
    else sign = -1.0;
    for( size_t a = 0; a < _tic.size(); a++ ) {

	if( _cind == 0 ) {
	    // xruler
	    double x = _tic[a]._x;
	    cm.transform( x );
	    double y = _endpt[1]+(_endpt[3]-_endpt[1])*(x-_endpt[0])/(_endpt[2]-_endpt[0]);

	    // Tic
	    cairo_move_to( cairo, x, y-sign*_ticlen_in );
	    cairo_line_to( cairo, x, y+sign*_ticlen_out );
	    cairo_stroke( cairo );

	} else {
	    // yruler
	    double y = _tic[a]._x;
	    cm.transform( y );
	    double x = _endpt[0]+(_endpt[2]-_endpt[0])*(y-_endpt[1])/(_endpt[3]-_endpt[1]);

	    // Tic
	    cairo_move_to( cairo, x-sign*_ticlen_in,  y );
	    cairo_line_to( cairo, x+sign*_ticlen_out, y );
	    cairo_stroke( cairo );
	}
    }

    if( _label_enabled ) {
	// Draw tics labels
	for( size_t a = 0; a < _tic.size(); a++ )
	    _tic[a]._label.draw( cairo );

	// Draw axis labels
	_axislabel.draw( cairo );
    }
}


void Ruler::get_bbox( cairo_t *cairo, double bbox[4], Coordmapper1D &cm, bool recalculate )
{
    double bb[4];

    if( recalculate )
	calculate( cairo, cm, true );

    if( _cind == 0 ) {
	// xruler
	bbox[0] = _endpt[0] < _endpt[2] ? _endpt[0] : _endpt[2];
	bbox[2] = _endpt[0] > _endpt[2] ? _endpt[0] : _endpt[2];
	if( _endpt[1] < _endpt[3] ) {
	    bbox[1] = _endpt[1];
	    bbox[3] = _endpt[3];
	} else {
	    bbox[1] = _endpt[3];
	    bbox[3] = _endpt[1];
	}

	// Tics, tic labels and axis label (y-direction)
	if( _indir ) {
	    bbox[1] -= _ticlen_in;
	    bbox[3] += _ticlen_out;
	    if( _label_enabled ) {
		for( size_t a = 0; a < _tic.size(); a++ ) {
		    _tic[a]._label.get_bbox( cairo, bb );
		    if( bbox[3] < bb[3] )
			bbox[3] = bb[3];
		}
		if( _axislabel.get_text() != "" ) {
		    _axislabel.get_bbox( cairo, bb );
		    if( bbox[3] < bb[3] )
			bbox[3] = bb[3];
		}
	    }
	} else {
	    bbox[1] -= _ticlen_out;
	    bbox[3] += _ticlen_in;
	    if( _label_enabled ) {
		for( size_t a = 0; a < _tic.size(); a++ ) {
		    _tic[a]._label.get_bbox( cairo, bb );
		    if( bbox[1] > bb[1] )
			bbox[1] = bb[1];
		}
		if( _axislabel.get_text() != "" ) {
		    _axislabel.get_bbox( cairo, bb );
		    if( bbox[1] > bb[1] )
			bbox[1] = bb[1];
		}
	    }
	}

	// Tic labels (x-direction)
	if( _label_enabled && _tic.size() > 0 ) {
	    // Minimum
	    _tic[0]._label.get_bbox( cairo, bb );
	    if( bb[0] < bbox[0] )
		bbox[0] = bb[0];

	    // Maximum
	    _tic[_tic.size()-1]._label.get_bbox( cairo, bb );
	    if( bb[2] > bbox[2] )
		bbox[2] = bb[2];
	}

    } else {
	// yruler
	bbox[1] = _endpt[1] < _endpt[3] ? _endpt[1] : _endpt[3];
	bbox[3] = _endpt[1] > _endpt[3] ? _endpt[1] : _endpt[3];
	if( _endpt[0] < _endpt[2] ) {
	    bbox[0] = _endpt[0];
	    bbox[2] = _endpt[2];
	} else {
	    bbox[0] = _endpt[2];
	    bbox[2] = _endpt[0];
	}

	// Tics, tic labels and axis label (x-direction)
	if( _indir ) {
	    bbox[0] -= _ticlen_in;
	    bbox[2] += _ticlen_out;
	    if( _label_enabled ) {
		for( size_t a = 0; a < _tic.size(); a++ ) {
		    _tic[a]._label.get_bbox( cairo, bb );
		    if( bb[2] > bbox[2] )
			bbox[2] = bb[2];
		}
		if( _axislabel.get_text() != "" ) {
		    _axislabel.get_bbox( cairo, bb );
		    if( bb[2] > bbox[2] )
			bbox[2] = bb[2];
		}
	    }
	} else {
	    bbox[0] -= _ticlen_out;
	    bbox[2] += _ticlen_in;
	    if( _label_enabled ) {
		for( size_t a = 0; a < _tic.size(); a++ ) {
		    _tic[a]._label.get_bbox( cairo, bb );
		    if( bb[0] < bbox[0] )
			bbox[0] = bb[0];
		}
		if( _axislabel.get_text() != "" ) {
		    _axislabel.get_bbox( cairo, bb );
		    if( bb[0] < bbox[0] )
			bbox[0] = bb[0];
		}
	    }
	}

	// Tic labels (y-direction)
	if( _label_enabled && _tic.size() > 0 ) {
	    // Minimum
	    _tic[0]._label.get_bbox( cairo, bb );
	    if( bb[3] > bbox[3] )
		bbox[3] = bb[3];

	    // Maximum
	    _tic[_tic.size()-1]._label.get_bbox( cairo, bb );
	    if( bb[1] < bbox[1] )
		bbox[1] = bb[1];
	}
    }
}


void Ruler::debug_print( std::ostream &os ) const 
{
    os << "**Ruler\n";
    os << "  range = (" << _range[0] << ", " << _range[1] << ")\n";
    os << "  autorange = (" << _autorange[0] << ", " << _autorange[1] << ")\n";
    os << "  endpt = (" << _endpt[0] << ", " << _endpt[1] << ", "
       << _endpt[2] << ", " << _endpt[3] << ")\n";
    os << "  axis_label = " << _axislabel << "\n";
    os << "  indir = " << _indir << "\n";
    os << "  cind = " << _cind << "\n";
    os << "  label_enabled = " << _label_enabled << "\n";
    os << "  tic_count = " << _tic.size() << "\n";
    for( size_t a = 0; a < _tic.size(); a++ ) {
	os << "    tic[" << a << "] = {" << _tic[a]._x << ", \"" << _tic[a]._label << "\"}\n";
    }
}

