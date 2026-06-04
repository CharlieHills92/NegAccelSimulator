/*! \file mydxfmtext.cpp
 *  \brief DXF mtext entity
 */

/* Copyright (c) 2010-2012,2014,2022 Taneli Kalvas. All rights reserved.
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


#include "mydxfmtext.hpp"
#include "mydxffont.hpp"


MyDXFMText::MyDXFMText() 
    : _text_height(1.0), _rect_width(0.0), _text_width(0.0), _vert_height(0.0),
      _spacing_fac(1.0), _rotation(0.0), _attachment_point(ATTACHMENT_POINT_TOP_LEFT),
      _drawing_direction(DRAWING_DIRECTION_LEFT_TO_RIGHT), _line_spacing(2),
      _style("STANDARD"), _extrusion(Vec3D(0,0,1)), _xaxis(1,0,0)
{

}


MyDXFMText::MyDXFMText( class MyDXFFile *dxf )
    : _text_height(1.0), _rect_width(0.0), _text_width(0.0), _vert_height(0.0),
      _spacing_fac(1.0), _rotation(0.0), _attachment_point(ATTACHMENT_POINT_TOP_LEFT),
      _drawing_direction(DRAWING_DIRECTION_LEFT_TO_RIGHT), _line_spacing(2),
      _style("STANDARD"), _extrusion(Vec3D(0,0,1)), _xaxis(1,0,0)
{
    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with entity

	else if( dxf->group_get_code() == 10 )
	    _p[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 )
	    _p[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 )
	    _p[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 40 )
	    _text_height = dxf->group_get_double();
	else if( dxf->group_get_code() == 41 )
	    _rect_width = dxf->group_get_double();
	else if( dxf->group_get_code() == 42 )
	    _text_width = dxf->group_get_double();
	else if( dxf->group_get_code() == 43 )
	    _vert_height = dxf->group_get_double();
	else if( dxf->group_get_code() == 44 )
	    _spacing_fac = dxf->group_get_double();

	else if( dxf->group_get_code() == 71 )
	    _attachment_point = dxf->group_get_int16();
	else if( dxf->group_get_code() == 72 )
	    _drawing_direction = dxf->group_get_int16();
	else if( dxf->group_get_code() == 73 )
	    _line_spacing = dxf->group_get_int16();

	else if( dxf->group_get_code() == 7 )
	    _style = dxf->group_get_string();

	else if( dxf->group_get_code() == 210 )
	    _extrusion[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 220 )
	    _extrusion[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 230 )
	    _extrusion[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 11 )
	    _xaxis[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 12 )
	    _xaxis[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 13 ) {
	    _xaxis[2] = dxf->group_get_double();
	    // Set rotation according to xaxis
	    _rotation = atan2(_xaxis[1],_xaxis[0]);
	} else if( dxf->group_get_code() == 50 ) {
	    _rotation = M_PI*dxf->group_get_double()/180.0;
	    // Set xaxis according to rotation
	    _xaxis[0] = sin(_rotation);
	    _xaxis[1] = cos(_rotation);
	    _xaxis[2] = 0.0;
	}

	else if( dxf->group_get_code() == 7 )
	    _style = dxf->group_get_string();


	else if( dxf->group_get_code() == 3 || dxf->group_get_code() == 1 )
	    _text += dxf->group_get_string();
	
	else
	    process_group( dxf );
    }

    // Normalize direction vectors
    _xaxis.normalize();
    _extrusion.normalize();
}


void MyDXFMText::explode( MyDXFEntities *ent, class MyDXFFile *dxf, const Transformation *t ) const
{
    MyDXFMText *text = new MyDXFMText( *this );

    // Transform points
    Vec3D p = t->transform_point( text->_p );
    Vec3D h = text->_p + Vec3D(0,text->_text_height,0);
    text->_p = p;
    text->_text_height = norm2( h-p );

    // Add to entities
    ent->add_entity( text );
}


void MyDXFMText::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "MTEXT" );
    write_common( dxf, ostr );

    // Chop into 250 character chunks. Code 1 for last (or only
    // group), others code 3.
    int done = 0;
    int remaining = _text.length();
    while( remaining ) {
	int code = 1;
	if( remaining > 250 )
	    code = 3;
	std::string s =_text.substr( done, 250 );
	done += s.length();
	remaining -= s.length();

	dxf->write_group( code, s.c_str() );
    }

    dxf->write_group( 10, _p[0] );
    dxf->write_group( 20, _p[1] );
    dxf->write_group( 30, _p[2] );

    dxf->write_group( 40, _text_height );
    dxf->write_group( 41, _rect_width );
    dxf->write_group( 42, _text_width );
    dxf->write_group( 43, _vert_height );
    dxf->write_group( 44, _spacing_fac );
    dxf->write_group( 50, _rotation );

    dxf->write_group( 71, _attachment_point );
    dxf->write_group( 72, _drawing_direction );
    dxf->write_group( 73, _line_spacing );

    dxf->write_group(  7, _style.c_str() );

    dxf->write_group( 210, _extrusion[0] );
    dxf->write_group( 220, _extrusion[1] );
    dxf->write_group( 230, _extrusion[2] );

    dxf->write_group( 11, _xaxis[0] );
    dxf->write_group( 21, _xaxis[1] );
    dxf->write_group( 31, _xaxis[2] );
}


void MyDXFMText::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		       const Transformation *t, const double range[4] ) const
{
    MyDXFFont mfont;
    if( dxf->wlevel() >= 2 )
	std::cout << "Warning: plotting for MText entity not fully implemented\n";

    // Check size
    Vec3D x, min, max;
    for( uint32_t i = 0; i < _text.length(); i++ ) {
	uint32_t c = _text[i];
	if( _text.substr(i,3) == "\\U+" ) {
	    // Unicode input
	    c = strtol( _text.substr(i+3,4).c_str(), NULL, 16 );
	    i+=6;
	} else if( _text.substr(i,2) == "\\~" ) {
	    c = ' ';
	    i++;
	} else if( _text.substr(i,3) == "%%%" ) {
	    c = '%';
	    i+=2;
	} else if( _text.substr(i,3) == "%%c" ) {
	    c = 0xF8; // diameter symbol
	    i+=2;
	} else if( _text.substr(i,3) == "%%d" ) {
	    c = 0xB0; // degrees symbol
	    i+=2;
	} else if( _text.substr(i,3) == "%%p" ) {
	    c = 0xB1; // plus/minus symbol
	    i+=2;
	}
	mfont.size( min, max, dxf, c, x );
    }    

    //std::cout << "min = " << min << "\n";
    //std::cout << "max = " << max << "\n";

    cairo_save( cairo );
    cairo_set_line_width( cairo, 0.5 );
    cairo_set_line_join( cairo, CAIRO_LINE_JOIN_BEVEL );
    Transformation t2 = *t;
    t2.translate_before( _p );
    t2.rotate_z_before( _rotation );
    t2.scale_before( Vec3D(_text_height,_text_height,_text_height) );

    // cursor starting location
    switch( _attachment_point ) {
    case ATTACHMENT_POINT_TOP_LEFT:
	x = Vec3D(0,-max[1],0);
	break;
    case ATTACHMENT_POINT_TOP_CENTER:
	x = Vec3D(-0.5*(max[0]+min[0]),-max[1],0);
	break;
    case ATTACHMENT_POINT_TOP_RIGHT:
	x = Vec3D(-max[0],-max[1],0);
	break;
    case ATTACHMENT_POINT_MIDDLE_LEFT:
	x = Vec3D(0,-0.5*(min[1]+max[1]),0);
	break;
    case ATTACHMENT_POINT_MIDDLE_CENTER:
	x = Vec3D(-0.5*(max[0]+min[0]),-0.5*(min[1]+max[1]),0);
	break;
    case ATTACHMENT_POINT_MIDDLE_RIGHT:
	x = Vec3D(-max[0],-0.5*(min[1]+max[1]),0);
	break;
    case ATTACHMENT_POINT_BOTTOM_LEFT:
	x = Vec3D(0,0,0);
	break;
    case ATTACHMENT_POINT_BOTTOM_CENTER:
	x = Vec3D(-0.5*(max[0]+min[0]),0,0);
	break;
    case ATTACHMENT_POINT_BOTTOM_RIGHT:
	x = Vec3D(-max[0],0,0);
	break;
    default:
	if( dxf->wlevel() >= 2 )
	    std::cout << "Warning: unknown attachment point in MText\n";	
	break;
    }

    for( uint32_t i = 0; i < _text.length(); i++ ) {
	uint32_t c = _text[i];
	if( _text.substr(i,3) == "\\U+" ) {
	    // Unicode input
	    c = strtol( _text.substr(i+3,4).c_str(), NULL, 16 );
	    i+=6;
	} else if( _text.substr(i,2) == "\\~" ) {
	    c = ' ';
	    i++;
	} else if( _text.substr(i,3) == "%%%" ) {
	    c = '%';
	    i+=2;
	} else if( _text.substr(i,3) == "%%c" ) {
	    c = 0xF8; // diameter symbol
	    i+=2;
	} else if( _text.substr(i,3) == "%%d" ) {
	    c = 0xB0; // degrees symbol
	    i+=2;
	} else if( _text.substr(i,3) == "%%p" ) {
	    c = 0xB1; // plus/minus symbol
	    i+=2;
	}
	mfont.plot( dxf, cairo, &t2, range, c, x );
    }
    cairo_restore( cairo );
}


void MyDXFMText::get_bbox( Vec3D &min, Vec3D &max, 
			   const class MyDXFFile *dxf, const Transformation *t ) const
{
    if( dxf->wlevel() >= 2 )
	std::cout << "Warning: bounding box for MText entity not fully implemented\n";
    min = _p;
    max = _p;

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "MText bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFMText::scale( class MyDXFFile *dxf, double s )
{
    _p *= s;
    _text_height *= s;
    _rect_width *= s;
    _text_width *= s;
    _vert_height *= s;
}


void MyDXFMText::translate( class MyDXFFile *dxf, const Vec3D &dx )
{
    _p += dx;
}


void MyDXFMText::rotate_z( class MyDXFFile *dxf, double a )
{
    Transformation t = Transformation::rotation_z( a );
    _p = t.transform_point( _p );

    _rotation += a;
    // Enforce between 0 and 2 pi
    _rotation = _rotation - 2.0*M_PI*floor( _rotation/(2.0*M_PI) );
}


void MyDXFMText::debug_print( std::ostream &os ) const
{
    os << "  MTEXT\n";
    MyDXFFile::debug_print_format( os, "text", _text );
    MyDXFFile::debug_print_format( os, "p", _p );
    MyDXFFile::debug_print_format( os, "text_height", _text_height );
    MyDXFFile::debug_print_format( os, "rect_width", _rect_width );
    MyDXFFile::debug_print_format( os, "text_width", _text_width );
    MyDXFFile::debug_print_format( os, "vert_height", _vert_height );
    MyDXFFile::debug_print_format( os, "spacing_fac", _spacing_fac );
    MyDXFFile::debug_print_format( os, "rotation", _rotation );
    MyDXFFile::debug_print_format( os, "attachment_point", _attachment_point );
    MyDXFFile::debug_print_format( os, "drawing_direction", _drawing_direction );
    MyDXFFile::debug_print_format( os, "line_spacing", _line_spacing );
    MyDXFFile::debug_print_format( os, "style", _style );
    MyDXFFile::debug_print_format( os, "extrusion", _extrusion );
    MyDXFFile::debug_print_format( os, "xaxis", _xaxis );
}

