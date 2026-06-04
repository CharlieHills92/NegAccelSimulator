/*! \file mydxffont.cpp
 *  \brief DXF text/font handling
 */

/* Copyright (c) 2012 Taneli Kalvas. All rights reserved.
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


#include "mydxffont.hpp"


// Font operation types
#define FONT_OP_MOVE 0
#define FONT_OP_LINE 1
#define FONT_OP_CURVE1 2
#define FONT_OP_CURVE2 3
#define FONT_OP_CURSOR 4


/* ************************** */


MyDXFFont::Glyph::Glyph()
    : _c(0)
{

}


MyDXFFont::Glyph::Glyph( uint32_t c )
    : _c(c)
{

}


MyDXFFont::Glyph::~Glyph()
{

}


void MyDXFFont::Glyph::reset( uint32_t c )
{
    _c = c;
    _fontops.clear();
}


void MyDXFFont::Glyph::add_op( unsigned char op, double x, double y )
{
    _fontops.push_back( FontOp( op, x, y ) );
}


void MyDXFFont::Glyph::draw( cairo_t *cairo, const Transformation *t, Vec3D &x ) const
{
    Vec3D y1;
    for( uint32_t i = 0; i < _fontops.size(); i++ ) {

	if( _fontops[i].op() == FONT_OP_MOVE ) {
	    y1 = t->transform_point(x+_fontops[i].v());
	    cairo_move_to( cairo, y1[0], y1[1] );

	} else if( _fontops[i].op() == FONT_OP_LINE ) {
	    y1 = t->transform_point(x+_fontops[i].v());
	    cairo_line_to( cairo, y1[0], y1[1] );

	} else if( _fontops[i].op() == FONT_OP_CURVE1 ) {
	    y1 = t->transform_point(x+_fontops[i].v());

	} else if( _fontops[i].op() == FONT_OP_CURVE2 ) {
	    Vec3D y2 = t->transform_point(x+_fontops[i].v());
	    cairo_curve_to( cairo, y1[0], y1[1], y1[0], y1[1], y2[0], y2[1] );

	} else if( _fontops[i].op() == FONT_OP_CURSOR ) {
	    x += _fontops[i].v();

	}
    }
    cairo_stroke( cairo );
}


void MyDXFFont::Glyph::cursor_advance( Vec3D &x ) const
{
    for( uint32_t i = 0; i < _fontops.size(); i++ ) {
	if( _fontops[i].op() == FONT_OP_CURSOR )
	    x += _fontops[i].v();
    }
}


/* ************************** */


MyDXFFont::MyDXFFont()
{
    // Build glyph database

    // Special 'zero' for unknown glyphs
    Glyph g( 0 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_LINE,   0.440, 0.660 );
    g.add_op( FONT_OP_LINE,   0.000, 0.660 );
    g.add_op( FONT_OP_LINE,   0.000, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );
    
    // Capital alphabet

    g.reset( 'A' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.333, 1.000 );
    g.add_op( FONT_OP_LINE,   0.666, 0.000 );
    g.add_op( FONT_OP_MOVE,   0.333*0.280, 0.280 );
    g.add_op( FONT_OP_LINE,   0.666-0.333*0.280, 0.280 );
    g.add_op( FONT_OP_CURSOR, 1.000, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'B' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.550 );
    g.add_op( FONT_OP_LINE,   0.290, 0.550 );
    g.add_op( FONT_OP_CURVE1, 0.290+0.275, 0.550 );
    g.add_op( FONT_OP_CURVE2, 0.290+0.275, 0.275 );
    g.add_op( FONT_OP_CURVE1, 0.290+0.275, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.290, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.280, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.280+0.225, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.280+0.225,1.000-0.225 );
    g.add_op( FONT_OP_CURVE1, 0.280+0.225, 1.000-0.450 );
    g.add_op( FONT_OP_CURVE2, 0.280, 1.000-0.450 );
    g.add_op( FONT_OP_CURSOR, 0.880, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'C' );
    g.add_op( FONT_OP_MOVE,   0.440, 1.000 );
    g.add_op( FONT_OP_LINE,   0.220, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.780 );
    g.add_op( FONT_OP_LINE,   0.000, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.000 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'D' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.500 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.330, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.550, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.550, 0.780 );
    g.add_op( FONT_OP_LINE,   0.550, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.550, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.500 );
    g.add_op( FONT_OP_CURSOR, 0.890, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'E' );
    g.add_op( FONT_OP_MOVE,   0.440, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.440, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.500 );
    g.add_op( FONT_OP_LINE,   0.440, 0.500 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'F' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.440, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.500 );
    g.add_op( FONT_OP_LINE,   0.440, 0.500 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'G' );
    g.add_op( FONT_OP_MOVE,   0.390, 0.550 );
    g.add_op( FONT_OP_LINE,   0.560, 0.550 );
    g.add_op( FONT_OP_LINE,   0.560, 0.000 );
    g.add_op( FONT_OP_LINE,   0.220, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.220 );
    g.add_op( FONT_OP_LINE,   0.000, 0.780 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 1.000 );
    g.add_op( FONT_OP_LINE,   0.560, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.890, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'H' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.500 );
    g.add_op( FONT_OP_LINE,   0.560, 0.500 );
    g.add_op( FONT_OP_MOVE,   0.560, 0.000 );
    g.add_op( FONT_OP_LINE,   0.560, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.890, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'I' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.330, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'J' );
    g.add_op( FONT_OP_MOVE,   0.330, 1.000 );
    g.add_op( FONT_OP_LINE,   0.330, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.330, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.110, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.660, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'K' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.380 );
    g.add_op( FONT_OP_LINE,   0.550, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.190, 0.380+0.190/0.550*0.620 );
    g.add_op( FONT_OP_LINE,   0.550, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'L' );
    g.add_op( FONT_OP_MOVE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'M' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.330, 0.440 );
    g.add_op( FONT_OP_LINE,   0.660, 1.000 );
    g.add_op( FONT_OP_LINE,   0.660, 0.000 );
    g.add_op( FONT_OP_CURSOR, 1.000, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'N' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.550, 0.000 );
    g.add_op( FONT_OP_LINE,   0.550, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.890, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'O' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.500 );
    g.add_op( FONT_OP_LINE,   0.000, 0.780 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 1.000 );
    g.add_op( FONT_OP_LINE,   0.330, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.550, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.550, 0.780 );
    g.add_op( FONT_OP_LINE,   0.550, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.550, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.000 );
    g.add_op( FONT_OP_LINE,   0.220, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.220 );
    g.add_op( FONT_OP_LINE,   0.000, 0.500 );
    g.add_op( FONT_OP_CURSOR, 0.890, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'P' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.330, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.550, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.550, 0.780 );
    g.add_op( FONT_OP_LINE,   0.550, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.550, 0.440 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.440 );
    g.add_op( FONT_OP_LINE,   0.000, 0.440 );
    g.add_op( FONT_OP_CURSOR, 0.890, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'Q' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.500 );
    g.add_op( FONT_OP_LINE,   0.000, 0.780 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 1.000 );
    g.add_op( FONT_OP_LINE,   0.330, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.550, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.550, 0.780 );
    g.add_op( FONT_OP_LINE,   0.550, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.550, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.000 );
    g.add_op( FONT_OP_LINE,   0.220, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.220 );
    g.add_op( FONT_OP_LINE,   0.000, 0.500 );
    g.add_op( FONT_OP_MOVE,   0.660, 0.000 );
    g.add_op( FONT_OP_LINE,   0.330, 0.220 );
    g.add_op( FONT_OP_CURSOR, 1.000, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'R' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.330, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.550, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.550, 0.780 );
    g.add_op( FONT_OP_LINE,   0.550, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.550, 0.440 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.440 );
    g.add_op( FONT_OP_LINE,   0.000, 0.440 );
    g.add_op( FONT_OP_MOVE,   0.330, 0.440 );
    g.add_op( FONT_OP_LINE,   0.550, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.890, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'S' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.330, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.550, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.550, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.550, 0.500 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.500 );
    g.add_op( FONT_OP_LINE,   0.220, 0.500 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.500 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.720 );
    g.add_op( FONT_OP_LINE,   0.000, 0.780 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 1.000 );
    g.add_op( FONT_OP_LINE,   0.550, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.890, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'T' );
    g.add_op( FONT_OP_MOVE,   0.330, 0.000 );
    g.add_op( FONT_OP_LINE,   0.330, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.660, 1.000 );
    g.add_op( FONT_OP_CURSOR, 1.000, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'U' );
    g.add_op( FONT_OP_MOVE,   0.550, 1.000 );
    g.add_op( FONT_OP_LINE,   0.550, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.550, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.000 );
    g.add_op( FONT_OP_LINE,   0.220, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.220 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.890, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'V' );
    g.add_op( FONT_OP_MOVE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.330, 0.000 );
    g.add_op( FONT_OP_LINE,   0.660, 1.000 );
    g.add_op( FONT_OP_CURSOR, 1.000, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'X' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.660, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.660, 0.000 );
    g.add_op( FONT_OP_CURSOR, 1.000, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'Y' );
    g.add_op( FONT_OP_MOVE,   0.330, 0.000 );
    g.add_op( FONT_OP_LINE,   0.330, 0.550 );
    g.add_op( FONT_OP_LINE,   0.660, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.330, 0.550 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_CURSOR, 1.000, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'Z' );
    g.add_op( FONT_OP_MOVE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.550, 1.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.550, 0.000 );
    g.add_op( FONT_OP_CURSOR, 1.000, 0.000 );
    _glyphs.push_back( g );

    g.reset( 0xC5 ); // A with ring
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.333, 1.000 );
    g.add_op( FONT_OP_LINE,   0.666, 0.000 );
    g.add_op( FONT_OP_MOVE,   0.333*0.280, 0.280 );
    g.add_op( FONT_OP_LINE,   0.666-0.333*0.280, 0.280 );
    g.add_op( FONT_OP_MOVE,   0.333, 1.280 );
    g.add_op( FONT_OP_CURVE1, 0.413, 1.280 );
    g.add_op( FONT_OP_CURVE2, 0.413, 1.200 );
    g.add_op( FONT_OP_CURVE1, 0.413, 1.120 );
    g.add_op( FONT_OP_CURVE2, 0.333, 1.120 );
    g.add_op( FONT_OP_CURVE1, 0.253, 1.120 );
    g.add_op( FONT_OP_CURVE2, 0.253, 1.200 );
    g.add_op( FONT_OP_CURVE1, 0.253, 1.280 );
    g.add_op( FONT_OP_CURVE2, 0.333, 1.280 );
    g.add_op( FONT_OP_CURSOR, 1.000, 0.000 );
    _glyphs.push_back( g );

    g.reset( 0xC4 ); // A with umlauts
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.333, 1.000 );
    g.add_op( FONT_OP_LINE,   0.666, 0.000 );
    g.add_op( FONT_OP_MOVE,   0.333*0.280, 0.280 );
    g.add_op( FONT_OP_LINE,   0.666-0.333*0.280, 0.280 );
    g.add_op( FONT_OP_MOVE,   0.100, 1.000 );
    g.add_op( FONT_OP_LINE,   0.100, 1.060 );
    g.add_op( FONT_OP_MOVE,   0.566, 1.000 );
    g.add_op( FONT_OP_LINE,   0.566, 1.060 );
    g.add_op( FONT_OP_CURSOR, 1.000, 0.000 );
    _glyphs.push_back( g );

    g.reset( 0xD6 ); // O with umlauts
    g.add_op( FONT_OP_MOVE,   0.000, 0.500 );
    g.add_op( FONT_OP_LINE,   0.000, 0.780 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 1.000 );
    g.add_op( FONT_OP_LINE,   0.330, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.550, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.550, 0.780 );
    g.add_op( FONT_OP_LINE,   0.550, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.550, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.000 );
    g.add_op( FONT_OP_LINE,   0.220, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.220 );
    g.add_op( FONT_OP_LINE,   0.000, 0.500 );
    g.add_op( FONT_OP_MOVE,   0.100, 1.100 );
    g.add_op( FONT_OP_LINE,   0.100, 1.160 );
    g.add_op( FONT_OP_MOVE,   0.450, 1.100 );
    g.add_op( FONT_OP_LINE,   0.450, 1.160 );
    g.add_op( FONT_OP_CURSOR, 0.890, 0.000 );
    _glyphs.push_back( g );

    // Small alphabet

    g.reset( 'a' );
    g.add_op( FONT_OP_MOVE,   0.060, 0.660 );
    g.add_op( FONT_OP_LINE,   0.280, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.380 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_LINE,   0.160, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.160 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.330 );
    g.add_op( FONT_OP_CURVE2, 0.160, 0.330 );
    g.add_op( FONT_OP_LINE,   0.440, 0.330 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'b' );
    g.add_op( FONT_OP_MOVE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.280, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.160 );
    g.add_op( FONT_OP_LINE,   0.440, 0.390 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.550 );
    g.add_op( FONT_OP_CURVE2, 0.280, 0.550 );
    g.add_op( FONT_OP_LINE,   0.000, 0.550 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'c' );
    g.add_op( FONT_OP_MOVE,   0.330, 0.000 );
    g.add_op( FONT_OP_LINE,   0.160, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.160 );
    g.add_op( FONT_OP_LINE,   0.000, 0.500 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.160, 0.660 );
    g.add_op( FONT_OP_LINE,   0.330, 0.660 );
    g.add_op( FONT_OP_CURSOR, 0.550, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'd' );
    g.add_op( FONT_OP_MOVE,   0.440, 1.000 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_LINE,   0.160, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.160 );
    g.add_op( FONT_OP_LINE,   0.000, 0.500 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.160, 0.660 );
    g.add_op( FONT_OP_LINE,   0.440, 0.660 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'e' );
    g.add_op( FONT_OP_MOVE,   0.440, 0.000 );
    g.add_op( FONT_OP_LINE,   0.160, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.160 );
    g.add_op( FONT_OP_LINE,   0.000, 0.440 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.440 );
    g.add_op( FONT_OP_LINE,   0.440, 0.330 );
    g.add_op( FONT_OP_LINE,   0.000, 0.330 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'f' );
    g.add_op( FONT_OP_MOVE,   0.110, 0.000 );
    g.add_op( FONT_OP_LINE,   0.110, 0.840 );
    g.add_op( FONT_OP_CURVE1, 0.110, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.270, 1.000 );
    g.add_op( FONT_OP_LINE,   0.330, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.660 );
    g.add_op( FONT_OP_LINE,   0.330, 0.660 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'g' );
    g.add_op( FONT_OP_MOVE,   0.000,-0.330 );
    g.add_op( FONT_OP_LINE,   0.280,-0.330 );
    g.add_op( FONT_OP_CURVE1, 0.440,-0.330 );
    g.add_op( FONT_OP_CURVE2, 0.440,-0.170 );
    g.add_op( FONT_OP_LINE,   0.440, 0.660 );
    g.add_op( FONT_OP_LINE,   0.160, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.500 );
    g.add_op( FONT_OP_LINE,   0.000, 0.160 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.160, 0.000 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'h' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.660 );
    g.add_op( FONT_OP_LINE,   0.280, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.500 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'i' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.660 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.940 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.330, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'j' );
    g.add_op( FONT_OP_MOVE,   0.000,-0.330 );
    g.add_op( FONT_OP_LINE,   0.060,-0.440 );
    g.add_op( FONT_OP_CURVE1, 0.220,-0.440 );
    g.add_op( FONT_OP_CURVE2, 0.220,-0.060 );
    g.add_op( FONT_OP_LINE,   0.220, 0.660 );
    g.add_op( FONT_OP_MOVE,   0.220, 0.940 );
    g.add_op( FONT_OP_LINE,   0.220, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.550, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'k' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.380 );
    g.add_op( FONT_OP_LINE,   0.440, 0.660 );
    g.add_op( FONT_OP_MOVE,   0.150, 0.380+0.150/0.440*0.280 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.550, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'l' );
    g.add_op( FONT_OP_MOVE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.110 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.110, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.440, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'm' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.660 );
    g.add_op( FONT_OP_LINE,   0.500, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.660, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.660, 0.500 );
    g.add_op( FONT_OP_LINE,   0.660, 0.000 );
    g.add_op( FONT_OP_MOVE,   0.330, 0.000 );
    g.add_op( FONT_OP_LINE,   0.330, 0.660 );
    g.add_op( FONT_OP_CURSOR, 1.000, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'n' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.660 );
    g.add_op( FONT_OP_LINE,   0.280, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.500 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'o' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.330 );
    g.add_op( FONT_OP_LINE,   0.000, 0.440 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.440 );
    g.add_op( FONT_OP_LINE,   0.440, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.220 );
    g.add_op( FONT_OP_LINE,   0.000, 0.330 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'p' );
    g.add_op( FONT_OP_MOVE,   0.000,-0.330 );
    g.add_op( FONT_OP_LINE,   0.000, 0.660 );
    g.add_op( FONT_OP_LINE,   0.280, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.500 );
    g.add_op( FONT_OP_LINE,   0.440, 0.160 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.280, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'q' );
    g.add_op( FONT_OP_MOVE,   0.440,-0.330 );
    g.add_op( FONT_OP_LINE,   0.440, 0.660 );
    g.add_op( FONT_OP_LINE,   0.160, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.500 );
    g.add_op( FONT_OP_LINE,   0.000, 0.160 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.160, 0.000 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'r' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.660 );
    g.add_op( FONT_OP_LINE,   0.220, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.330, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.550 );
    g.add_op( FONT_OP_CURSOR, 0.660, 0.000 );
    _glyphs.push_back( g );

    g.reset( 's' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.170, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.330, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.160 );
    g.add_op( FONT_OP_CURVE1, 0.330, 0.330 );
    g.add_op( FONT_OP_CURVE2, 0.170, 0.330 );
    g.add_op( FONT_OP_LINE,   0.160, 0.330 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.330 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.490 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.160, 0.660 );
    g.add_op( FONT_OP_LINE,   0.330, 0.660 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 't' );
    g.add_op( FONT_OP_MOVE,   0.160, 0.000 );
    g.add_op( FONT_OP_LINE,   0.160, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.660 );
    g.add_op( FONT_OP_LINE,   0.330, 0.660 );
    g.add_op( FONT_OP_CURSOR, 0.660, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'u' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.660 );
    g.add_op( FONT_OP_LINE,   0.000, 0.160 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.160, 0.000 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_LINE,   0.440, 0.660 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'v' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.660 );
    g.add_op( FONT_OP_LINE,   0.220, 0.000 );
    g.add_op( FONT_OP_LINE,   0.440, 0.660 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'x' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.440, 0.660 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.660 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'y' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.660 );
    g.add_op( FONT_OP_LINE,   0.220, 0.000 );
    g.add_op( FONT_OP_MOVE,   0.440, 0.660 );
    g.add_op( FONT_OP_LINE,   0.110,-0.330 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 'z' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.660 );
    g.add_op( FONT_OP_LINE,   0.440, 0.660 );
    g.add_op( FONT_OP_LINE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 0xE5 ); // a with ring
    g.add_op( FONT_OP_MOVE,   0.060, 0.660 );
    g.add_op( FONT_OP_LINE,   0.280, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.380 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_LINE,   0.160, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.160 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.330 );
    g.add_op( FONT_OP_CURVE2, 0.160, 0.330 );
    g.add_op( FONT_OP_LINE,   0.440, 0.330 );
    g.add_op( FONT_OP_MOVE,   0.250, 0.980 );
    g.add_op( FONT_OP_CURVE1, 0.330, 0.980 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.900 );
    g.add_op( FONT_OP_CURVE1, 0.330, 0.820 );
    g.add_op( FONT_OP_CURVE2, 0.250, 0.820 );
    g.add_op( FONT_OP_CURVE1, 0.170, 0.820 );
    g.add_op( FONT_OP_CURVE2, 0.170, 0.900 );
    g.add_op( FONT_OP_CURVE1, 0.170, 0.980 );
    g.add_op( FONT_OP_CURVE2, 0.250, 0.980 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 0xE4 ); // a with umlauts
    g.add_op( FONT_OP_MOVE,   0.060, 0.660 );
    g.add_op( FONT_OP_LINE,   0.280, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.380 );
    g.add_op( FONT_OP_LINE,   0.440, 0.000 );
    g.add_op( FONT_OP_LINE,   0.160, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.160 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.330 );
    g.add_op( FONT_OP_CURVE2, 0.160, 0.330 );
    g.add_op( FONT_OP_LINE,   0.440, 0.330 );
    g.add_op( FONT_OP_MOVE,   0.060, 0.800 );
    g.add_op( FONT_OP_LINE,   0.060, 0.860 );
    g.add_op( FONT_OP_MOVE,   0.380, 0.800 );
    g.add_op( FONT_OP_LINE,   0.380, 0.860 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 0xF6 ); // o with umlauts
    g.add_op( FONT_OP_MOVE,   0.000, 0.330 );
    g.add_op( FONT_OP_LINE,   0.000, 0.440 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.660 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.440 );
    g.add_op( FONT_OP_LINE,   0.440, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.220 );
    g.add_op( FONT_OP_LINE,   0.000, 0.330 );
    g.add_op( FONT_OP_MOVE,   0.060, 0.800 );
    g.add_op( FONT_OP_LINE,   0.060, 0.860 );
    g.add_op( FONT_OP_MOVE,   0.380, 0.800 );
    g.add_op( FONT_OP_LINE,   0.380, 0.860 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    // Numbers

    g.reset( '0' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.500 );
    g.add_op( FONT_OP_LINE,   0.000, 0.780 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.440, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.780 );
    g.add_op( FONT_OP_LINE,   0.440, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.220 );
    g.add_op( FONT_OP_LINE,   0.000, 0.500 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.330 );
    g.add_op( FONT_OP_LINE,   0.440, 0.670 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '1' );
    g.add_op( FONT_OP_MOVE,   0.220, 0.000 );
    g.add_op( FONT_OP_LINE,   0.220, 1.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.780 );
    g.add_op( FONT_OP_CURSOR, 0.550, 0.000 );
    _glyphs.push_back( g );

    g.reset( '2' );
    g.add_op( FONT_OP_MOVE,   0.440, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.440, 0.780 );
    g.add_op( FONT_OP_CURVE1, 0.440, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.780 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '3' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.280, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.160 );
    g.add_op( FONT_OP_LINE,   0.440, 0.390 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.550 );
    g.add_op( FONT_OP_CURVE2, 0.280, 0.550 );
    g.add_op( FONT_OP_LINE,   0.110, 0.550 );
    g.add_op( FONT_OP_MOVE,   0.280, 0.550 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.550 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.710 );
    g.add_op( FONT_OP_LINE,   0.440, 0.840 );
    g.add_op( FONT_OP_CURVE1, 0.440, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.280, 1.000 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '4' );
    g.add_op( FONT_OP_MOVE,   0.220, 1.000 );
    g.add_op( FONT_OP_LINE,   0.000, 0.220 );
    g.add_op( FONT_OP_LINE,   0.550, 0.220 );
    g.add_op( FONT_OP_MOVE,   0.390, 0.000 );
    g.add_op( FONT_OP_LINE,   0.390, 0.440 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '5' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.220, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.500 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.500 );
    g.add_op( FONT_OP_LINE,   0.000, 0.500 );
    g.add_op( FONT_OP_LINE,   0.000, 0.780 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.550, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '6' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.550 );
    g.add_op( FONT_OP_LINE,   0.220, 0.550 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.550 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.330 );
    g.add_op( FONT_OP_LINE,   0.440, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.220 );
    g.add_op( FONT_OP_LINE,   0.000, 0.780 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 1.000 );
    g.add_op( FONT_OP_LINE,   0.440, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '7' );
    g.add_op( FONT_OP_MOVE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.440, 1.000 );
    g.add_op( FONT_OP_LINE,   0.160, 0.000 );
    g.add_op( FONT_OP_MOVE,   0.220, 0.550 );
    g.add_op( FONT_OP_LINE,   0.440, 0.550 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '8' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.275 );
    g.add_op( FONT_OP_LINE,   0.000, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.220 );
    g.add_op( FONT_OP_LINE,   0.440, 0.330 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.550 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.550 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.550 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.330 );
    g.add_op( FONT_OP_LINE,   0.000, 0.275 );
    g.add_op( FONT_OP_MOVE,   0.030, 0.775 );
    g.add_op( FONT_OP_LINE,   0.030, 0.740 );
    g.add_op( FONT_OP_CURVE1, 0.030, 0.550 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.550 );
    g.add_op( FONT_OP_CURVE1, 0.410, 0.550 );
    g.add_op( FONT_OP_CURVE2, 0.410, 0.740 );
    g.add_op( FONT_OP_LINE,   0.410, 0.810 );
    g.add_op( FONT_OP_CURVE1, 0.410, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.030, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.030, 0.810 );
    g.add_op( FONT_OP_LINE,   0.030, 0.775 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '9' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.220, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.440, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.440, 0.220 );
    g.add_op( FONT_OP_LINE,   0.440, 0.780 );
    g.add_op( FONT_OP_CURVE1, 0.440, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.780 );
    g.add_op( FONT_OP_LINE,   0.000, 0.660 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.440 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.440 );
    g.add_op( FONT_OP_LINE,   0.440, 0.440 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    // Special

    g.reset( ' ' );
    g.add_op( FONT_OP_CURSOR, 0.440, 0.000 );
    _glyphs.push_back( g );

    g.reset( '+' );
    g.add_op( FONT_OP_MOVE,   0.275, 0.250 );
    g.add_op( FONT_OP_LINE,   0.275, 0.750 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.500 );
    g.add_op( FONT_OP_LINE,   0.550, 0.500 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '-' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.500 );
    g.add_op( FONT_OP_LINE,   0.550, 0.500 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( 0xB1 ); // plus/minus
    g.add_op( FONT_OP_MOVE,   0.275, 0.350 );
    g.add_op( FONT_OP_LINE,   0.275, 0.850 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.600 );
    g.add_op( FONT_OP_LINE,   0.550, 0.600 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.150 );
    g.add_op( FONT_OP_LINE,   0.550, 0.150 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '=' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.200 );
    g.add_op( FONT_OP_LINE,   0.550, 0.200 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.400 );
    g.add_op( FONT_OP_LINE,   0.550, 0.400 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '#' );
    g.add_op( FONT_OP_MOVE,   0.200, 0.800 );
    g.add_op( FONT_OP_LINE,   0.200, 0.000 );
    g.add_op( FONT_OP_MOVE,   0.400, 0.800 );
    g.add_op( FONT_OP_LINE,   0.400, 0.000 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.600 );
    g.add_op( FONT_OP_LINE,   0.600, 0.600 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.200 );
    g.add_op( FONT_OP_LINE,   0.600, 0.200 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '\\' );
    g.add_op( FONT_OP_MOVE,   0.000, 1.000 );
    g.add_op( FONT_OP_LINE,   0.330, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '/' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.330, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '(' );
    g.add_op( FONT_OP_MOVE,   0.220, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.220 );
    g.add_op( FONT_OP_LINE,   0.000, 0.780 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.440, 0.000 );
    _glyphs.push_back( g );

    g.reset( ')' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.220, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 0.220 );
    g.add_op( FONT_OP_LINE,   0.220, 0.780 );
    g.add_op( FONT_OP_CURVE1, 0.220, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.440, 0.000 );
    _glyphs.push_back( g );

    g.reset( ',' );
    g.add_op( FONT_OP_MOVE,   0.110, 0.060 );
    g.add_op( FONT_OP_LINE,   0.110, 0.000 );
    g.add_op( FONT_OP_LINE,   0.000,-0.330 );
    g.add_op( FONT_OP_CURSOR, 0.440, 0.000 );
    _glyphs.push_back( g );

    g.reset( '.' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.060 );
    g.add_op( FONT_OP_LINE,   0.000, 0.000 );
    g.add_op( FONT_OP_CURSOR, 0.330, 0.000 );
    _glyphs.push_back( g );

    g.reset( '!' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.060 );
    g.add_op( FONT_OP_LINE,   0.000, 0.000 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.330 );
    g.add_op( FONT_OP_LINE,   0.000, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.330, 0.000 );
    _glyphs.push_back( g );

    g.reset( '?' );
    g.add_op( FONT_OP_MOVE,   0.220, 0.060 );
    g.add_op( FONT_OP_LINE,   0.220, 0.000 );
    g.add_op( FONT_OP_MOVE,   0.220, 0.330 );
    g.add_op( FONT_OP_LINE,   0.220, 0.440 );
    g.add_op( FONT_OP_LINE,   0.440, 0.780 );
    g.add_op( FONT_OP_CURVE1, 0.440, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.780 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );

    g.reset( '$' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.330, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.550, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.550, 0.220 );
    g.add_op( FONT_OP_CURVE1, 0.550, 0.500 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.500 );
    g.add_op( FONT_OP_LINE,   0.220, 0.500 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.500 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.720 );
    g.add_op( FONT_OP_LINE,   0.000, 0.780 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.220, 1.000 );
    g.add_op( FONT_OP_LINE,   0.550, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.275,-0.150 );
    g.add_op( FONT_OP_LINE,   0.275, 1.150 );
    g.add_op( FONT_OP_CURSOR, 0.890, 0.000 );
    _glyphs.push_back( g );

    g.reset( '%' );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.840, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.165, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.330, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.835 );
    g.add_op( FONT_OP_CURVE1, 0.330, 0.670 );
    g.add_op( FONT_OP_CURVE2, 0.165, 0.670 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.670 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.835 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.165, 1.000 );
    g.add_op( FONT_OP_MOVE,   0.675, 0.330 );
    g.add_op( FONT_OP_CURVE1, 0.840, 0.330 );
    g.add_op( FONT_OP_CURVE2, 0.840, 0.165 );
    g.add_op( FONT_OP_CURVE1, 0.840, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.675, 0.000 );
    g.add_op( FONT_OP_CURVE1, 0.510, 0.000 );
    g.add_op( FONT_OP_CURVE2, 0.510, 0.165 );
    g.add_op( FONT_OP_CURVE1, 0.510, 0.330 );
    g.add_op( FONT_OP_CURVE2, 0.675, 0.330 );
    g.add_op( FONT_OP_CURSOR, 1.180, 0.000 );
    _glyphs.push_back( g );

    g.reset( 0xB0 ); // degrees
    g.add_op( FONT_OP_MOVE,   0.165, 1.000 );
    g.add_op( FONT_OP_CURVE1, 0.330, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.330, 0.835 );
    g.add_op( FONT_OP_CURVE1, 0.330, 0.670 );
    g.add_op( FONT_OP_CURVE2, 0.165, 0.670 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.670 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.835 );
    g.add_op( FONT_OP_CURVE1, 0.000, 1.000 );
    g.add_op( FONT_OP_CURVE2, 0.165, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.550, 0.000 );
    _glyphs.push_back( g );

    g.reset( 0xF8 ); // degrees
    g.add_op( FONT_OP_MOVE,   0.000, 0.500 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.250 );
    g.add_op( FONT_OP_CURVE2, 0.250, 0.250 );
    g.add_op( FONT_OP_CURVE1, 0.500, 0.250 );
    g.add_op( FONT_OP_CURVE2, 0.500, 0.500 );
    g.add_op( FONT_OP_CURVE1, 0.500, 0.750 );
    g.add_op( FONT_OP_CURVE2, 0.250, 0.750 );
    g.add_op( FONT_OP_CURVE1, 0.000, 0.750 );
    g.add_op( FONT_OP_CURVE2, 0.000, 0.500 );
    g.add_op( FONT_OP_MOVE,   0.000, 0.000 );
    g.add_op( FONT_OP_LINE,   0.500, 1.000 );
    g.add_op( FONT_OP_CURSOR, 0.780, 0.000 );
    _glyphs.push_back( g );


}


MyDXFFont::~MyDXFFont()
{

}


void MyDXFFont::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		      const Transformation *t, const double range[4],
		      uint32_t c, Vec3D &x ) const
{
    uint32_t index;
    for( index = 0; index < _glyphs.size(); index++ ) {
	if( _glyphs[index].ch() == c ) {
	    _glyphs[index].draw( cairo, t, x );
	    break;
	}	
    }
    if( index == _glyphs.size() ) {
	// Didn't find matching glyph, print glyph zero
	if( dxf->wlevel() >= 2 )
	    std::cout << "Missing glyph \'" << (char)(c&&0xFF) << "\' (" << std::hex << c << ")\n";
	_glyphs[0].draw( cairo, t, x );
    }
}


void MyDXFFont::size( Vec3D &min, Vec3D &max, const class MyDXFFile *dxf, uint32_t c, Vec3D &x ) const
{
    uint32_t index;
    for( index = 0; index < _glyphs.size(); index++ ) {
	if( _glyphs[index].ch() == c ) {
	    Vec3D xold = x;
	    _glyphs[index].cursor_advance( x );
	    // In x-direction
	    if( xold[0] < min[0] )
		min[0] = xold[0];
	    if( x[0] < min[0] )
		min[0] = x[0];
	    if( xold[0] > max[0] )
		max[0] = xold[0];
	    if( x[0] > max[0] )
		max[0] = x[0];
	    // In y-direction
	    if( xold[1] < min[1] )
		min[1] = xold[1];
	    if( x[1] < min[1] )
		min[1] = x[1];
	    if( xold[1]+1.0 > max[1] )
		max[1] = xold[1]+1.0;
	    if( x[1]+1.0 > max[1] )
		max[1] = x[1]+1.0;
	    break;
	}	
    }
    if( index == _glyphs.size() ) {
	// Didn't find matching glyph, print glyph zero
	_glyphs[0].cursor_advance( x );
    }
}

