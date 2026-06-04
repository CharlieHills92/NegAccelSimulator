/*! \file label.cpp
 *  \brief Plot labels
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

#include "label.hpp"
#include "fonts.hpp"
#include "compmath.hpp"
#include <cmath>


#define LINESPREAD 1.12


Label::Label()
    : _size(12.0), _family("Times"), _slant(CAIRO_FONT_SLANT_NORMAL), 
      _weight(CAIRO_FONT_WEIGHT_NORMAL), _color(Vec3D(0,0,0)), 
      _xalign(0.0), _yalign(0.0), _yzeroext(false), _rotation(0.0), 
      _xlocation(0.0), _ylocation(0.0)
{

}


Label::Label( const Label &label )
    : _text(label._text), _size(label._size), _family(label._family), _slant(label._slant), 
      _weight(label._weight), _color(label._color), 
      _xalign(label._xalign), _yalign(label._yalign), _yzeroext(label._yzeroext),
      _rotation(label._rotation), 
      _xlocation(label._xlocation), _ylocation(label._ylocation)
{
}


Label::Label( const std::string &text )
    : _text(text), _size(12.0), _family("Times"), _slant(CAIRO_FONT_SLANT_NORMAL), 
      _weight(CAIRO_FONT_WEIGHT_NORMAL), _color(Vec3D(0,0,0)), 
      _xalign(0.0), _yalign(0.0), _yzeroext(false), _rotation(0.0), 
      _xlocation(0.0), _ylocation(0.0)
{

}


Label &Label::operator=( const Label &label )
{
    _text      = label._text;
    _size      = label._size;
    _family    = label._family;
    _slant     = label._slant;
    _weight    = label._weight;
    _color     = label._color;
    _xalign    = label._xalign;
    _yalign    = label._yalign;
    _yzeroext  = label._yzeroext;
    _rotation  = label._rotation;
    _xlocation = label._xlocation;
    _ylocation = label._ylocation;

    return( *this );
}


Label::~Label()
{

}


void Label::set_font_size( double size )
{
    _size = size;
}


double Label::get_font_size( void ) const
{
    return( _size );
}


void Label::set_font_family( const std::string &family )
{
    _family = family;
}


void Label::set_font_slant(cairo_font_slant_t slant )
{
    _slant = slant;
}
    

void Label::set_font_weight( cairo_font_weight_t weight )
{
    _weight = weight;
}


void Label::set_color( const Vec3D &color )
{
    _color = color;
}


void Label::set_location( double x, double y )
{
    _xlocation = x;
    _ylocation = y;
}


void Label::set_rotation( double angle )
{
    _rotation = angle;
}


void Label::set_alignment( double x, double y, bool yzeroext )
{
    _xalign = x;
    _yalign = y;
    _yzeroext = yzeroext;
}


void Label::set_text( const std::string &text )
{
    _text = text;
}


std::string Label::get_text( void ) const
{
    return( _text );
}


void Label::process_parsed( cairo_t *cairo, const std::string &text, cairo_text_extents_t *extents0, 
			    double x0, double y0, double &x, double &y ) const
{
    if( extents0 ) {
	// Get extents
	cairo_text_extents_t extents1;
	//std::cout << "Getting extents for \'" << text << "\'\n";
	fontlib.text_extents( cairo, text, &extents1 );
	// Combine extents
	//std::cout << "0:\n";
	//std::cout << "x_advance = " << extents0->x_advance << "\n";
	//std::cout << "y_advance = " << extents0->y_advance << "\n";
	//std::cout << "width     = " << extents0->width     << "\n";
	//std::cout << "height    = " << extents0->height    << "\n";
	//std::cout << "x_bearing = " << extents0->x_bearing << "\n";
	//std::cout << "y_bearing = " << extents0->y_bearing << "\n";
	FontLib::combine_extents( extents0, x0, y0, &extents1, x, y );
	//std::cout << "1:\n";
	//std::cout << "x_advance = " << extents1.x_advance << "\n";
	//std::cout << "y_advance = " << extents1.y_advance << "\n";
	//std::cout << "width     = " << extents1.width     << "\n";
	//std::cout << "height    = " << extents1.height    << "\n";
	//std::cout << "x_bearing = " << extents1.x_bearing << "\n";
	//std::cout << "y_bearing = " << extents1.y_bearing << "\n";
	x += extents1.x_advance;
	y += extents1.y_advance;
    } else {
	fontlib.draw_text( cairo, text, x, y );
    }
}


/* Recursive use for this function!!
 * (x0,y0) is the starting point for the label cursor, provided for calculating extents
 * (x,y) is the current cursor position
 *
 */
void Label::parse_latex( cairo_t *cairo, const std::string &text, cairo_text_extents_t *extents0, 
			 double x0, double y0, double &x, double &y ) const
{
    std::string s;
    
    size_t a = 0;
    while( a < text.length() ) {
	char c = text[a];
	if( c == '\n' || c == '\\' || c == '_' || c == '^' ) {
	    // Control character next

	    // Stroke plain text so far
	    if( s != "" )
		//std::cout << "Stroke plain text\n";
		process_parsed( cairo, s, extents0, x0, y0, x, y );
	    s = "";

	    if( c == '\n' || (c == '\\' && a+1 < text.length() && text[a+1] == '\\') ) {

		// Newline
		cairo_matrix_t matrix;
		cairo_get_font_matrix( cairo, &matrix );
		x = x0;
		y += LINESPREAD*matrix.yy;
		a++;
		if( c == '\\' )
		    a++;

	    } else if( c == '\\' ) {

		// Command
		//std::cout << "Read in command\n";
		a++;
		while( a < text.length() ) {
		    c = text[a];
		    if( c == ' ' || c == '\n' || c == '\\' || c == '_' || c == '^' || 
			c == '=' || c == '-' || c == '+' || c == '/' || c == '*' ) {
			// Process command
			//if( s == "" )
			//std::cout << "s = \'" << s << "\'\n";
			int i = 0;
			while( FontLib::symbols[i].name != NULL ) {
			    if( s == FontLib::symbols[i].name )
				break;
			    i++;
			}
			if( FontLib::symbols[i].name != NULL ) {
			    // Stroke symbol
			    std::string ss = FontLib::symbols[i].ucode;
			    //std::cout << "Stroke symbol\n";
			    process_parsed( cairo, ss, extents0, x0, y0, x, y );
			    s = "";
			} else {
			    // Symbol not found, stroke symbol name
			    //std::cout << "Stroke symbol name\n";
			    process_parsed( cairo, s, extents0, x0, y0, x, y );
			    s = "";
			}
			if( c == ' ' ) {
			    // Dump if space
			    a++;
			}
			break;
		    }
		    s += text[a];
		    a++;
		}
	    } else {
		a++;
	    }
	} else {
	    // Plain text continues
	    s += text[a];
	    a++;
	}
    }

    // Stroke plain text at end
    if( s != "" ) {
	//std::cout << "Stroke plain text at end\n";
	process_parsed( cairo, s, extents0, x0, y0, x, y );
    }
}


void Label::draw( cairo_t *cairo )
{
    if( _text == "" )
	return;

    // Color
    cairo_set_source_rgba( cairo, _color[0], _color[1], _color[2], 1.0 );

    // Font size
    cairo_set_font_size( cairo, _size );

    // Process text with alignment
    fontlib.push_font( _family, _slant, _weight );
    if( comp_isnan(_xalign) || comp_isnan(_yalign) || 
	comp_isinf(_xalign) || comp_isinf(_yalign) ) {

	cairo_save( cairo );
	cairo_translate( cairo, _xlocation, _ylocation );
	cairo_rotate( cairo, _rotation );
	cairo_move_to( cairo, 0.0, 0.0 );
	double x = 0.0, y = 0.0;
	parse_latex( cairo, _text, NULL, 0, 0, x, y );
	cairo_restore( cairo );

    } else {

	cairo_text_extents_t ext;
	ext.x_bearing = 0.0;
	ext.y_bearing = 0.0;
	ext.width     = 0.0;
	ext.height    = 0.0;
	ext.x_advance = 0.0;
	ext.y_advance = 0.0;
	double x = 0.0, y = 0.0;
	parse_latex( cairo, _text, &ext, 0, 0, x, y );

	if( _yzeroext ) {
	    cairo_text_extents_t extz;
	    extz.x_bearing = 0.0;
	    extz.y_bearing = 0.0;
	    extz.width     = 0.0;
	    extz.height    = 0.0;
	    extz.x_advance = 0.0;
	    extz.y_advance = 0.0;
	    double x = 0.0, y = 0.0;
	    parse_latex( cairo, "0", &extz, 0, 0, x, y );
	    ext.y_bearing = extz.y_bearing;
	    ext.height    = extz.height;
	}

	cairo_save( cairo );
	cairo_translate( cairo, _xlocation, _ylocation );
	cairo_rotate( cairo, -_rotation );
	cairo_translate( cairo,
			 -ext.x_bearing-_xalign*ext.width,
			 -ext.y_bearing-(1.0-_yalign)*ext.height );
	cairo_move_to( cairo, 0.0, 0.0 );
	x = 0.0; y = 0.0;
	parse_latex( cairo, _text, NULL, 0, 0, x, y );
	cairo_restore( cairo );

    }

    // Done with font
    fontlib.pop_font();
}


void Label::get_extents( cairo_t *cairo, cairo_text_extents_t *extents )
{
    if( extents == NULL )
	return;

    // Font size
    cairo_set_font_size( cairo, _size );

    // Process text
    fontlib.push_font( _family, _slant, _weight );
    extents->x_bearing = 0.0;
    extents->y_bearing = 0.0;
    extents->width     = 0.0;
    extents->height    = 0.0;
    extents->x_advance = 0.0;
    extents->y_advance = 0.0;
    double x = 0.0, y = 0.0;
    parse_latex( cairo, _text, extents, 0, 0, x, y );
    fontlib.pop_font();
}


void Label::get_bbox( cairo_t *cairo, double bbox[4] ) const
{
    if( _text == "" ) {
	bbox[0] = bbox[2] = _xlocation;
	bbox[1] = bbox[3] = _ylocation;
	return;
    }

    // Font size
    cairo_set_font_size( cairo, _size );

    // Process text
    fontlib.push_font( _family, _slant, _weight );
    cairo_text_extents_t extents;
    extents.x_bearing = 0.0;
    extents.y_bearing = 0.0;
    extents.width     = 0.0;
    extents.height    = 0.0;
    extents.x_advance = 0.0;
    extents.y_advance = 0.0;
    double x = 0.0, y = 0.0;
    parse_latex( cairo, _text, &extents, 0, 0, x, y );

    double baseline = 0.0;
    if( _yzeroext ) {
	cairo_text_extents_t extz;
	extz.x_bearing = 0.0;
	extz.y_bearing = 0.0;
	extz.width     = 0.0;
	extz.height    = 0.0;
	extz.x_advance = 0.0;
	extz.y_advance = 0.0;
	x = 0.0, y = 0.0;
	parse_latex( cairo, "0", &extz, 0, 0, x, y );
	baseline = _ylocation-extz.y_bearing-(1.0-_yalign)*extz.height;
    }
    fontlib.pop_font();

    // Calculate bounding box
    double bb[4];
    bb[0] = _xlocation-_xalign*extents.width;
    bb[2] = _xlocation+(1.0-_xalign)*extents.width;
    if( _yzeroext ) {
	bb[1] = baseline+extents.y_bearing;
 	bb[3] = baseline+extents.y_bearing+extents.height;
    } else {
	bb[1] = _ylocation-(1.0-_yalign)*extents.height;
	bb[3] = _ylocation+_yalign*extents.height;
    }

    if( _rotation == 0.0 ) {
	bbox[0] = bb[0];
	bbox[1] = bb[1];
	bbox[2] = bb[2];
	bbox[3] = bb[3];
	return;
    }

    // Rotate bounding box corners
    double corner[8] = { bb[0], bb[1],
			 bb[0], bb[3],
			 bb[2], bb[3],
			 bb[2], bb[1] };
    double sinr = sin( -_rotation );
    double cosr = cos( -_rotation );
    for( size_t a = 0; a < 4; a++ ) {
	double dx = corner[2*a+0] - _xlocation;
	double dy = corner[2*a+1] - _ylocation;
	corner[2*a+0] = cosr*dx - sinr*dy + _xlocation;
	corner[2*a+1] = sinr*dx + cosr*dy + _ylocation;
    }

    /*
    for( size_t a = 0; a < 4; a++ ) {
	cairo_set_source_rgb( cairo, 1, 0, 0 );
	cairo_set_line_width( cairo, 1 );
	cairo_arc( cairo, corner[2*a+0], corner[2*a+1], 2, 0, 2.0*M_PI );
	cairo_fill( cairo );
    }
    */

    // Find minima and maxima
    bbox[0] = corner[0];
    bbox[1] = corner[1];
    bbox[2] = corner[0];
    bbox[3] = corner[1];
    for( size_t a = 1; a < 4; a++ ) {
	if( corner[2*a+0] < bbox[0] )
	    bbox[0] = corner[2*a+0];
	if( corner[2*a+0] > bbox[2] )
	    bbox[2] = corner[2*a+0];
	if( corner[2*a+1] < bbox[1] )
	    bbox[1] = corner[2*a+1];
	if( corner[2*a+1] > bbox[3] )
	    bbox[3] = corner[2*a+1];
    }

    //cairo_set_source_rgb( cairo, 1, 0, 0 );
    //cairo_set_line_width( cairo, 1 );
    //cairo_rectangle( cairo, bbox[0], bbox[1], bbox[2]-bbox[0], bbox[3]-bbox[1] );
    //cairo_stroke( cairo );
}


std::ostream &operator<<( std::ostream &os, const Label &label ) 
{
    os << label._text;
    return( os );
}
