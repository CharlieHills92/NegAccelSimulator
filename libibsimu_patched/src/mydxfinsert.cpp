/*! \file mydxfinsert.cpp
 *  \brief DXF insert entity
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
#include "mydxfinsert.hpp"


MyDXFInsert::MyDXFInsert()
    : _rotation(0.0), _col_count(1), _row_count(1), _col_spacing(0.0), _row_spacing(0.0)
{
    // Default values
    _scale[0] = _scale[1] = _scale[2] = 1.0;
}



MyDXFInsert::MyDXFInsert( class MyDXFFile *dxf )
    : _rotation(0.0), _col_count(1), _row_count(1), _col_spacing(0.0), _row_spacing(0.0)
{
#ifdef MYDXF_DEBUG
    std::cout << "  Reading entity INSERT\n";
#endif

    // Default values
    _scale[0] = _scale[1] = _scale[2] = 1.0;

    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with entity

	else if( dxf->group_get_code() == 2 )
	    _block_name = dxf->group_get_string();

	else if( dxf->group_get_code() == 10 )
	    _p[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 )
	    _p[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 )
	    _p[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 41 )
	    _scale[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 42 )
	    _scale[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 43 )
	    _scale[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 44 )
	    _col_spacing = dxf->group_get_double();
	else if( dxf->group_get_code() == 45 )
	    _row_spacing = dxf->group_get_double();

	else if( dxf->group_get_code() == 50 ) {
	    double rot = M_PI*dxf->group_get_double()/180.0;
	    // Enforce between 0 and 2 pi
	    _rotation = rot - 2.0*M_PI*floor( rot/(2.0*M_PI) );
	}

	else if( dxf->group_get_code() == 70 )
	    _col_count = dxf->group_get_int16();
	else if( dxf->group_get_code() == 71 )
	    _row_count = dxf->group_get_int16();

	else
	    process_group( dxf );
    }

#ifdef MYDXF_DEBUG
    std::cout << *this;
#endif
}


void MyDXFInsert::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "INSERT" );
    write_common( dxf, ostr );

    dxf->write_group( 2, _block_name.c_str() );

    dxf->write_group( 10, _p[0] );
    dxf->write_group( 20, _p[1] );
    dxf->write_group( 30, _p[2] );

    dxf->write_group( 41, _scale[0] );
    dxf->write_group( 42, _scale[1] );
    dxf->write_group( 43, _scale[2] );

    dxf->write_group( 44, _col_spacing );
    dxf->write_group( 45, _row_spacing );

    dxf->write_group( 50, 180*_rotation/M_PI );

    dxf->write_group( 70, _col_count );
    dxf->write_group( 71, _row_count );
}


void MyDXFInsert::plot( const class MyDXFFile *dxf, cairo_t *cairo, 
			const Transformation *t, const double range[4] ) const
{
#ifdef MYDXF_DEBUG_PLOT
    std::cout << "MyDXFInsert::plot()\n";
#endif

    // Fetch block data
    const MyDXFBlocks *blocks = dxf->get_blocks();
    const MyDXFBlock *b = blocks->get_by_name( _block_name );
    if( !b )
	return;

    Transformation t2 = *t;
    t2.translate_before( _p );
    t2.rotate_z_before( _rotation );
    t2.scale_before( _scale );

    for( int16_t col = 0; col < _col_count; col++ ) {
	for( int16_t row = 0; row < _row_count; row++ ) {
	    Transformation t3 = t2;
	    t3.translate_before( Vec3D( col*_col_spacing, row*_row_spacing, 0.0 ) );
	    b->plot( dxf, cairo, &t3, range );
	}
    }

#ifdef MYDXF_DEBUG_PLOT
    std::cout << "MyDXFInsert::plot() done\n";
#endif
}


void MyDXFInsert::get_bbox( Vec3D &min, Vec3D &max, 
			    const class MyDXFFile *dxf, const Transformation *t ) const
{
    min = Vec3D( std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    max = Vec3D( -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity() );

    // Fetch block data
    const MyDXFBlocks *blocks = dxf->get_blocks();
    const MyDXFBlock *b = blocks->get_by_name( _block_name );
    if( !b )
	return;

    Transformation t2 = *t;
    t2.translate_before( _p );
    t2.rotate_z_before( _rotation );
    t2.scale_before( _scale );

    for( int16_t col = 0; col < _col_count; col++ ) {
	for( int16_t row = 0; row < _row_count; row++ ) {
	    Transformation t3 = t2;
	    t3.translate_before( Vec3D( col*_col_spacing, row*_row_spacing, 0.0 ) );
	    Vec3D mi, ma;
	    b->get_bbox( mi, ma, dxf, &t3 );

	    for( int b = 0; b < 3; b++ ) {
		if( mi[b] < min[b] )
		    min[b] = mi[b];
		if( ma[b] > max[b] )
		    max[b] = ma[b];
	    }
	}
    }

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "Insert bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFInsert::scale( class MyDXFFile *dxf, double s )
{
    _scale *= s;

    // Fetch block data
    //MyDXFBlocks *blocks = dxf->get_blocks();
    //MyDXFBlock *b = blocks->get_by_name( _block_name );
    //if( !b )
    //return;
    //
    //b->scale( dxf, s );
}


void MyDXFInsert::translate( class MyDXFFile *dxf, const Vec3D &dx )
{
    _p += dx;

    // Fetch block data
    //MyDXFBlocks *blocks = dxf->get_blocks();
    //MyDXFBlock *b = blocks->get_by_name( _block_name );
    //if( !b )
    //return;
    //
    //b->translate( dxf, dx );
}


void MyDXFInsert::rotate_z( class MyDXFFile *dxf, double a )
{
    _rotation += a;
}


void MyDXFInsert::explode( class MyDXFEntities *ent, MyDXFFile *dxf, const Transformation *t ) const
{
    // Fetch block data
    const MyDXFBlocks *blocks = dxf->get_blocks();
    const MyDXFBlock *b = blocks->get_by_name( _block_name );
    if( !b )
	return;

    Transformation t2 = *t;
    t2.translate_before( _p );
    t2.rotate_z_before( _rotation );
    t2.scale_before( _scale );

    for( int16_t col = 0; col < _col_count; col++ ) {
	for( int16_t row = 0; row < _row_count; row++ ) {
	    Transformation t3 = t2;
	    t3.translate_before( Vec3D( col*_col_spacing, row*_row_spacing, 0.0 ) );
	    b->explode( ent, dxf, &t3 );
	}
    }    
}


void MyDXFInsert::debug_print( std::ostream &os ) const
{
    os << "  INSERT\n";
    MyDXFFile::debug_print_format( os, "block_name", _block_name );
    MyDXFFile::debug_print_format( os, "p", _p );
    MyDXFFile::debug_print_format( os, "scale", _scale );
    MyDXFFile::debug_print_format( os, "rotation", _rotation );
    MyDXFFile::debug_print_format( os, "col_count", _col_count );
    MyDXFFile::debug_print_format( os, "row_count", _row_count );
    MyDXFFile::debug_print_format( os, "col_spacing", _col_spacing );
    MyDXFFile::debug_print_format( os, "row_spacing", _row_spacing );
}



