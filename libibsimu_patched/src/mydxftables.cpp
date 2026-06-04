/*! \file mydxftables.cpp
 *  \brief DXF Tables 
 */

/* Copyright (c) 2010-2012 Taneli Kalvas. All rights reserved.
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


#include <iostream>
#include <iomanip>
#include "mydxftables.hpp"
#include "error.hpp"



/* ************************************************************************** *
 * DXF Entry                                                                  *
 * ************************************************************************** */

MyDXFTableEntry::MyDXFTableEntry()
{

}


void MyDXFTableEntry::process_group( class MyDXFFile *dxf )
{
    if( dxf->group_get_code() == 5 || dxf->group_get_code() == 105 )
	_handle = dxf->group_get_string();
    if( dxf->group_get_code() == 330 )
	_handle_to_owner = dxf->group_get_string();
}


void MyDXFTableEntry::write_common( class MyDXFFile *dxf, std::ofstream &ostr )
{
    //if( dynamic_cast<MyDXFTableEntryDimStyle *>( this ) )
    //dxf->write_group( 105, _handle.c_str() );
    //else
    dxf->write_group( 5, _handle.c_str() );
    //dxf->write_group( 330, _handle_to_owner.c_str() );
}


void MyDXFTableEntry::debug_print_common( std::ostream &os ) const
{
    MyDXFFile::debug_print_format( os, "handle", _handle );
    MyDXFFile::debug_print_format( os, "handle_to_owner", _handle_to_owner );
}


std::ostream &operator<<( std::ostream &os, const MyDXFTableEntry &e )
{
    e.debug_print( os );
    e.debug_print_common( os );
    return( os );
}


/* ************************************************************************** *
 * DXF Entry BlockRecord                                                      *
 * ************************************************************************** */

MyDXFTableEntryBlockRecord::MyDXFTableEntryBlockRecord( class MyDXFFile *dxf )
    : _units(0), _explodability(0), _scalability(0)
{
    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with record

	else if( dxf->group_get_code() == 2 )
	    _name = dxf->group_get_string();
	else if( dxf->group_get_code() == 70 )
	    _units = dxf->group_get_int16();
	else if( dxf->group_get_code() == 340 )
	    _handle_to_layout = dxf->group_get_string();
	else if( dxf->group_get_code() == 280 )
	    _explodability = dxf->group_get_int8();
	else if( dxf->group_get_code() == 281 )
	    _scalability = dxf->group_get_int8();

	else
	    process_group( dxf );
    }

}


MyDXFTableEntryBlockRecord::~MyDXFTableEntryBlockRecord()
{

}


void MyDXFTableEntryBlockRecord::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "BLOCK_RECORD" );
    write_common( dxf, ostr );

    dxf->write_group( 100, "AcDbSymbolTableRecord" );
    dxf->write_group( 100, "AcDbBlockTableRecord" );

    dxf->write_group( 2, _name.c_str() );
    dxf->write_group( 70, _units );
    dxf->write_group( 280, _explodability );
    dxf->write_group( 281, _scalability );
    dxf->write_group( 340, _handle_to_layout.c_str() );
}


void MyDXFTableEntryBlockRecord::debug_print( std::ostream &os ) const
{
    MyDXFFile::debug_print_format( os, "name", _name );
    MyDXFFile::debug_print_format( os, "units", _units );
    MyDXFFile::debug_print_format( os, "explodability", (int)_explodability );
    MyDXFFile::debug_print_format( os, "scalability", (int)_scalability );
    MyDXFFile::debug_print_format( os, "handle_to_layout", _handle_to_layout );
}


/* ************************************************************************** *
 * DXF Entry Layer                                                            *
 * ************************************************************************** */

MyDXFTableEntryLayer::MyDXFTableEntryLayer( class MyDXFFile *dxf )
    : _flags(0), _color(0), _plotting(true), _lineweight(0)
{
    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with record

	else if( dxf->group_get_code() == 2 )
	    _name = dxf->group_get_string();
	else if( dxf->group_get_code() == 6 )
	    _linetype = dxf->group_get_string();
	else if( dxf->group_get_code() == 70 )
	    _flags = dxf->group_get_int16();
	else if( dxf->group_get_code() == 62 )
	    _color = dxf->group_get_int16();
	else if( dxf->group_get_code() == 290 )
	    _plotting = dxf->group_get_bool();
	else if( dxf->group_get_code() == 370 )
	    _lineweight = dxf->group_get_int8();
	else if( dxf->group_get_code() == 390 )
	    _handle_to_plot_style_name = dxf->group_get_string();
	else if( dxf->group_get_code() == 347 )
	    _handle_to_material = dxf->group_get_string();

	else
	    process_group( dxf );
    }

}


MyDXFTableEntryLayer::~MyDXFTableEntryLayer()
{

}


void MyDXFTableEntryLayer::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "LAYER" );
    write_common( dxf, ostr );

    dxf->write_group( 100, "AcDbSymbolTableRecord" );
    dxf->write_group( 100, "AcDbLayerTableRecord" );

    dxf->write_group( 2, _name.c_str() );
    dxf->write_group( 6, _linetype.c_str() );
    dxf->write_group( 70, _flags );
    dxf->write_group( 62, _color );
    dxf->write_group( 290, _plotting );
    dxf->write_group( 370, _lineweight );
    dxf->write_group( 390, _handle_to_plot_style_name.c_str() );
    dxf->write_group( 347, _handle_to_material.c_str() );
}


void MyDXFTableEntryLayer::debug_print( std::ostream &os ) const
{
    MyDXFFile::debug_print_format( os, "name", _name );
    MyDXFFile::debug_print_format( os, "linetype", _linetype );
    MyDXFFile::debug_print_format( os, "flags", _flags );
    MyDXFFile::debug_print_format( os, "color", _color );
    MyDXFFile::debug_print_format( os, "plotting", _plotting );
    MyDXFFile::debug_print_format( os, "lineweight", (int)_lineweight );
    MyDXFFile::debug_print_format( os, "handle_to_plot_style_name", _handle_to_plot_style_name );
    MyDXFFile::debug_print_format( os, "handle_to_material", _handle_to_material );
}


/* ************************************************************************** *
 * DXF Entry Vport                                                            *
 * ************************************************************************** */

MyDXFTableEntryVport::MyDXFTableEntryVport( class MyDXFFile *dxf )
    : _flags(0), _ucs_type(0), _elevation(0.0)
{
    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 )
	    break; // Done with record

	else if( dxf->group_get_code() == 2 )
	    _name = dxf->group_get_string();
	else if( dxf->group_get_code() == 70 )
	    _flags = dxf->group_get_int16();

	else if( dxf->group_get_code() == 10 )
	    _vmin[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 )
	    _vmin[1] = dxf->group_get_double();

	else if( dxf->group_get_code() == 11 )
	    _vmax[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 21 )
	    _vmax[1] = dxf->group_get_double();

	else if( dxf->group_get_code() == 12 )
	    _vcenter[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 22 )
	    _vcenter[1] = dxf->group_get_double();

	else if( dxf->group_get_code() == 13 )
	    _snap_base[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 23 )
	    _snap_base[1] = dxf->group_get_double();

	else if( dxf->group_get_code() == 14 )
	    _snap_spacing[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 24 )
	    _snap_spacing[1] = dxf->group_get_double();

	else if( dxf->group_get_code() == 15 )
	    _grid_spacing[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 25 )
	    _grid_spacing[1] = dxf->group_get_double();

	else if( dxf->group_get_code() == 16 )
	    _view_direction[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 26 )
	    _view_direction[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 36 )
	    _view_direction[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 17 )
	    _view_target[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 27 )
	    _view_target[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 37 )
	    _view_target[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 110 )
	    _ucs_origin[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 120 )
	    _ucs_origin[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 130 )
	    _ucs_origin[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 111 )
	    _ucs_x[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 121 )
	    _ucs_x[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 131 )
	    _ucs_x[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 112 )
	    _ucs_y[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 122 )
	    _ucs_y[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 132 )
	    _ucs_y[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 79 )
	    _ucs_type = dxf->group_get_int16();

	else if( dxf->group_get_code() == 146 )
	    _elevation = dxf->group_get_double();

	else
	    process_group( dxf );
    }

}


MyDXFTableEntryVport::~MyDXFTableEntryVport()
{

}


void MyDXFTableEntryVport::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "VPORT" );
    write_common( dxf, ostr );

    dxf->write_group( 100, "AcDbSymbolTableRecord" );
    dxf->write_group( 100, "AcDbViewportTableRecord" );

    dxf->write_group( 2, _name.c_str() );
    dxf->write_group( 70, _flags );
    dxf->write_group( 10, _vmin[0] );
    dxf->write_group( 20, _vmin[1] );
    dxf->write_group( 11, _vmax[0] );
    dxf->write_group( 21, _vmax[1] );
    dxf->write_group( 12, _vcenter[0] );
    dxf->write_group( 22, _vcenter[1] );
    dxf->write_group( 13, _snap_base[0] );
    dxf->write_group( 23, _snap_base[1] );
    dxf->write_group( 14, _snap_spacing[0] );
    dxf->write_group( 24, _snap_spacing[1] );
    dxf->write_group( 15, _grid_spacing[0] );
    dxf->write_group( 25, _grid_spacing[1] );
    dxf->write_group( 16, _view_direction[0] );
    dxf->write_group( 26, _view_direction[1] );
    dxf->write_group( 36, _view_direction[2] );
    dxf->write_group( 17, _view_target[0] );
    dxf->write_group( 27, _view_target[1] );
    dxf->write_group( 37, _view_target[2] );

    dxf->write_group( 110, _ucs_origin[0] );
    dxf->write_group( 120, _ucs_origin[1] );
    dxf->write_group( 130, _ucs_origin[2] );
    dxf->write_group( 111, _ucs_x[0] );
    dxf->write_group( 121, _ucs_x[1] );
    dxf->write_group( 131, _ucs_x[2] );
    dxf->write_group( 112, _ucs_y[0] );
    dxf->write_group( 122, _ucs_y[1] );
    dxf->write_group( 132, _ucs_y[2] );
    dxf->write_group( 79, _ucs_type );
    dxf->write_group( 146, _elevation );
}


void MyDXFTableEntryVport::debug_print( std::ostream &os ) const
{
    MyDXFFile::debug_print_format( os, "name", _name );
    MyDXFFile::debug_print_format( os, "flags", _flags );
    MyDXFFile::debug_print_format( os, "vmin", _vmin );
    MyDXFFile::debug_print_format( os, "vmax", _vmax );
    MyDXFFile::debug_print_format( os, "vcenter", _vcenter );
    MyDXFFile::debug_print_format( os, "snap_base", _snap_base );
    MyDXFFile::debug_print_format( os, "snap_spacing", _snap_spacing );
    MyDXFFile::debug_print_format( os, "grid_spacing", _grid_spacing );
    MyDXFFile::debug_print_format( os, "view_direction", _view_direction );
    MyDXFFile::debug_print_format( os, "view_target", _view_target );
    MyDXFFile::debug_print_format( os, "ucs_origin", _ucs_origin );
    MyDXFFile::debug_print_format( os, "ucs_x", _ucs_x );
    MyDXFFile::debug_print_format( os, "ucs_y", _ucs_y );
    MyDXFFile::debug_print_format( os, "ucs_type", _ucs_type );
    MyDXFFile::debug_print_format( os, "elevation", _elevation );
}


/* ************************************************************************** *
 * DXF Table                                                                  *
 * ************************************************************************** */

MyDXFTable::MyDXFTable( const std::string &name, class MyDXFFile *dxf )
    : _name(name)
{
#ifdef MYDXF_DEBUG
    std::cout << "Reading Table " << _name << "\n";
#endif

    while( dxf->group_get_code() != -1 ) {
	
	if( dxf->group_get_code() == 0 ) {
	    if( dxf->group_get_string() == "ENDSEC" || dxf->group_get_string() == "ENDTAB" )
		break; // Done with table
	    if( dxf->group_get_string() == _name ) {
		if( dxf->group_get_string() == "BLOCK_RECORD" ) 
		    _entries.push_back( new MyDXFTableEntryBlockRecord( dxf ) );
		else if( dxf->group_get_string() == "LAYER" ) 
		    _entries.push_back( new MyDXFTableEntryLayer( dxf ) );
		else if( dxf->group_get_string() == "VPORT" ) 
		    _entries.push_back( new MyDXFTableEntryVport( dxf ) );
		else
		    throw Error( ERROR_LOCATION, "Error reading table entry on line " + 
				 to_string(dxf->linec()) );
	    }
	} else {

	    if( dxf->group_get_code() == 5 )
		_handle = dxf->group_get_string();
	    else if( dxf->group_get_code() == 330 )
		_handle_to_owner = dxf->group_get_string();

	    dxf->read_group();
	}
    }
}



MyDXFTable::~MyDXFTable()
{
    for( uint32_t i = 0; i < _entries.size(); i++ )
	delete _entries[i];
}


void MyDXFTable::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "TABLE" );
    dxf->write_group( 2, _name.c_str() );

    dxf->write_group( 5, _handle.c_str() );

    dxf->write_group( 100, "AcDbSymbolTable" );

    //dxf->write_group( 330, _handle_to_owner.c_str() );
    dxf->write_group( 70, (int16_t)_entries.size() );

    for( uint32_t i = 0; i < _entries.size(); i++ )
	_entries[i]->write( dxf, ostr );

    dxf->write_group( 0, "ENDTAB" );
}




void MyDXFTable::debug_print( std::ostream &os ) const
{
    os << "*** Table " << _name << " ****************************************\n";

    MyDXFFile::debug_print_format( os, "handle", _handle );
    MyDXFFile::debug_print_format( os, "handle_to_owner", _handle_to_owner );
    os << "\n";

    for( uint32_t i = 0; i < _entries.size(); i++ )
	os << (*_entries[i]) << "\n";
}






/* ************************************************************************** *
 * DXF Tables                                                                 *
 * ************************************************************************** */


MyDXFTables::MyDXFTables( class MyDXFFile *dxf )
    : _blockrecord(0), _layer(0), _vport(0)
{
#ifdef MYDXF_DEBUG
    std::cout << "Reading section TABLES\n";
#endif

    // Read TABLES section
    //
    // Ends in ENDSEC
    dxf->read_group();

    while( dxf->group_get_code() != -1 ) {
	
	if( dxf->group_get_code() != 0 ) {
	    dxf->read_group();
	    continue;
	}

	else if( dxf->group_get_string() == "ENDSEC" )
	    break; // Done with tables
	else if( dxf->group_get_string() == "TABLE" ) {
	    if( dxf->read_group() != 2 )
		throw Error( ERROR_LOCATION, "Error at start of table on line " + 
			     to_string(dxf->linec()) );

	    if( dxf->group_get_string() == "BLOCK_RECORD" )
		_blockrecord = new MyDXFTable( "BLOCK_RECORD", dxf );
	    else if( dxf->group_get_string() == "LAYER" )
		_layer = new MyDXFTable( "LAYER", dxf );
	    else if( dxf->group_get_string() == "VPORT" )
		_vport = new MyDXFTable( "VPORT", dxf );
	    else {
		// Unknown table
		if( dxf->wlevel() >= 2 )
		    std::cout << "Skipping unknown table \'" << dxf->group_get_string() << "\'\n";

		while( dxf->read_group() != -1 ) {
		    if( dxf->group_get_code() == 0 && dxf->group_get_string() == "ENDTAB" )
			break;
		}
		dxf->read_group();
	    }
	} else {
	    dxf->read_group();
	}
    }
}



MyDXFTables::~MyDXFTables()
{
    if( _layer ) 
	delete _layer;
    if( _blockrecord ) 
	delete _blockrecord;
    if( _vport ) 
	delete _vport;
}


void MyDXFTables::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "SECTION" );
    dxf->write_group( 2, "TABLES" );

    if( _blockrecord ) 
	_blockrecord->write( dxf, ostr );
    if( _layer ) 
	_layer->write( dxf, ostr );
    if( _vport ) 
	_vport->write( dxf, ostr );

    dxf->write_group( 0, "ENDSEC" );
}




void MyDXFTables::debug_print( std::ostream &os ) const
{
    os << "*** Section TABLES ****************************************\n";

    if( _blockrecord ) {
	_blockrecord->debug_print( os );
	os << "\n";
    }

    if( _layer ) {
	_layer->debug_print( os );
	os << "\n";
    }

    if( _vport ) {
	_vport->debug_print( os );
	os << "\n";
    }
}


