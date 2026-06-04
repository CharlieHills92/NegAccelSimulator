/*! \file mydxffont.hpp
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


#ifndef MY_DXF_FONT_HPP
#define MY_DXF_FONT_HPP 1


#include "mydxffile.hpp"


/*! \brief Class for drawing text characters in MyDXFFile.
 */
class MyDXFFont 
{
    class FontOp
    {
	unsigned char _op;
	double _x;
	double _y;
    public:

	FontOp( unsigned char op, double x, double y )
	    : _op(op), _x(x), _y(y) {}
	~FontOp() {}

	unsigned char op( void ) const { return( _op ); }
	Vec3D v( void ) const { return( Vec3D(_x,_y,0.0) ); }
    };

    class Glyph
    {
	uint32_t _c;
	std::vector<FontOp> _fontops;
    public:

	Glyph();
	Glyph( uint32_t c );
	~Glyph();

	void reset( uint32_t c );
	void add_op( unsigned char op, double x, double y );

	uint32_t ch( void ) const { return( _c ); }
	void draw( cairo_t *cairo, const Transformation *t, Vec3D &x ) const;
	void cursor_advance( Vec3D &x ) const;
    };

    std::vector<Glyph> _glyphs;

public:

    MyDXFFont();

    ~MyDXFFont();

    /*! \brief Plot character \a c with cairo.
     *
     *  Character \a c is drawn using transformation \t a to \a cairo
     *  surface with current style at location \a x. Location is advanced
     *  to next character location.
     */
    void plot( const class MyDXFFile *dxf, cairo_t *cairo, 
	       const Transformation *t, const double range[4],
	       uint32_t c, Vec3D &x ) const;


    void size( Vec3D &min, Vec3D &max, const class MyDXFFile *dxf, uint32_t c, Vec3D &x ) const;
};

#endif
