/*! \file mydxffile.cpp
 *  \brief DXF file.
 */

/* Copyright (c) 2010,2012 Taneli Kalvas. All rights reserved.
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


#include <string.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "mydxffile.hpp"
#include "error.hpp"


#define DXF_DEBUG_PRINT_FIELD_WIDTH 30


#define CODE_STRING(x) \
    (							 \
	((x) >= 0 && (x) <= 9) ||			 \
        (x) == 100 || (x) == 101 || (x) == 102 || (x) == 105 ||	 \
        ((x) >= 300 && (x) <= 369) ||			 \
        ((x) >= 390 && (x) <= 399) ||			 \
        ((x) >= 410 && (x) <= 419) ||			 \
        ((x) >= 430 && (x) <= 439) ||			 \
        ((x) >= 470 && (x) <= 481) ||			 \
        ((x) >= 999 && (x) <= 1009)			 \
    )

#define CODE_BOOL(x) \
    (				  \
	(x) >= 290 && (x) <= 299  \
    )

#define CODE_INT8(x) \
    (					\
	((x) >= 280 && (x) <= 289) ||	\
	((x) >= 370 && (x) <= 389)	\
    )

#define CODE_INT16(x) \
    (					\
	((x) >= 60 && (x) <= 79) ||	\
        ((x) >= 170 && (x) <= 179) ||   \
        ((x) >= 270 && (x) <= 289) ||	\
        ((x) >= 400 && (x) <= 409) ||   \
        ((x) >= 1060 && (x) <= 1070)    \
    )

#define CODE_INT32(x) \
    (						\
	((x) >= 90 && (x) <= 99) ||		\
	((x) >= 420 && (x) <= 429) ||		\
	((x) >= 440 && (x) <= 459) ||		\
	(x) == 1071				\
    )

#define CODE_INT64(x) \
    (				  \
	(x) >= 160 && (x) <= 169  \
    )

#define CODE_DOUBLE(x) \
    ( \
        ((x) >= 10 && (x) <= 59) || \
        ((x) >= 110 && (x) <= 149) || \
        ((x) >= 210 && (x) <= 239) || \
        ((x) >= 460 && (x) <= 469) || \
	((x) >= 1010 && (x) <= 1059) \
    )

	
MyDXFFile::MyDXFFile()
    : _wlevel(0), _header(0), _tables(0), _blocks(0), _entities(0)
{

}


void MyDXFFile::read( const std::string &filename )
{
    // Free old data
    if( _header )
	delete _header;
    if( _tables )
	delete _tables;
    if( _blocks )
	delete _blocks;
    if( _entities )
	delete _entities;

    // Open file
    _istr.open( filename.c_str() );
    if( !_istr.good() )
	throw Error( ERROR_LOCATION, "Couldn't open file \'" + filename + "\'" );

    // Check if binary type
    char buf[22];
    _istr.get( buf, 22 );
    if( !strncmp( buf, "AutoCAD Binary DXF\x0d\x0a\x1a\0", 22 ) )
	_ascii = false;
    else {
	_linec = 0;
	_ascii = true;
    }
    _istr.seekg( 0 );
    
    // Search for known section types
    while( read_group() != -1 ) {

#if MYDXF_DEBUG >= 2
	std::cout << "Searching for section start\n";
#endif

	if( group_get_code() == 0 && group_get_string() == "SECTION" ) {
	    if( read_group() != 2 )
		throw Error( ERROR_LOCATION, "Error at start of section on line " + 
			     to_string(_linec) );	

	    if( group_get_string() == "HEADER" ) {
		_header = new MyDXFHeader( this );
	    } else if( group_get_string() == "TABLES" ) {
		_tables = new MyDXFTables( this );
	    } else if( group_get_string() == "BLOCKS" ) {
		_blocks = new MyDXFBlocks( this );
	    } else if( group_get_string() == "ENTITIES" ) {
		_entities = new MyDXFEntities( this, false );
	    } else {
		// Unknown section
		if( wlevel() >= 2 )
		    std::cout << "Skipping unknown section \'" << group_get_string() << "\'\n";

		while( read_group() != -1 ) {
		    if( group_get_code() == 0 && group_get_string() == "ENDSEC" )
			break;
		}
	    }
	}
    }

    // Close file
    _istr.close();
}


MyDXFFile::MyDXFFile( const std::string &filename )
    : _wlevel(0), _header(0), _tables(0), _blocks(0), _entities(0)
{
    read( filename );
}


MyDXFFile::~MyDXFFile()
{
    if( _header )
	delete _header;
    if( _tables )
	delete _tables;
    if( _blocks )
	delete _blocks;
    if( _entities )
	delete _entities;
}


void MyDXFFile::write( const std::string &filename )
{
    _ostr.open( filename.c_str(), std::ios_base::trunc );

    if( _header )
	_header->write( this, _ostr );
    if( _tables )
	_tables->write( this, _ostr );
    if( _blocks )
	_blocks->write( this, _ostr );
    if( _entities )
	_entities->write( this, _ostr );

    write_group( 0, "EOF" );

    // Close file
    _ostr.close();
}


int MyDXFFile::group_get_code( void ) const
{
    return( _group_code );
}


#define GROUP_TYPE_UNKNOWN 0
#define GROUP_TYPE_STRING  1
#define GROUP_TYPE_DOUBLE  2
#define GROUP_TYPE_BOOL    3
#define GROUP_TYPE_INT8    4
#define GROUP_TYPE_INT16   5
#define GROUP_TYPE_INT32   6
#define GROUP_TYPE_INT64   7


std::string MyDXFFile::group_get_string( void ) const
{
    if( _group_type != GROUP_TYPE_STRING )
	throw Error( ERROR_LOCATION, "Wrong group type on line " + 
		     to_string(_linec) );	
    return( _group_string );
}


double MyDXFFile::group_get_double( void ) const
{
    if( _group_type != GROUP_TYPE_DOUBLE )
	throw Error( ERROR_LOCATION, "Wrong group type on line " + 
		     to_string(_linec) );	
    return( _group_double );
}


bool MyDXFFile::group_get_bool( void ) const
{
    if( _group_type != GROUP_TYPE_BOOL )
	throw Error( ERROR_LOCATION, "Wrong group type on line " + 
		     to_string(_linec) );	
    return( _group_bool );
}


int8_t MyDXFFile::group_get_int8( void ) const
{
    if( _group_type != GROUP_TYPE_INT8 )
	throw Error( ERROR_LOCATION, "Wrong group type on line " + 
		     to_string(_linec) );	
    return( _group_int8 );
}


int16_t MyDXFFile::group_get_int16( void ) const
{
    if( _group_type != GROUP_TYPE_INT16 )
	throw Error( ERROR_LOCATION, "Wrong group type on line " + 
		     to_string(_linec) );	
    return( _group_int16 );
}


int32_t MyDXFFile::group_get_int32( void ) const
{
    if( _group_type != GROUP_TYPE_INT32 )
	throw Error( ERROR_LOCATION, "Wrong group type on line " + 
		     to_string(_linec) );	
    return( _group_int32 );
}


int64_t MyDXFFile::group_get_int64( void ) const
{
    if( _group_type != GROUP_TYPE_INT64 )
	throw Error( ERROR_LOCATION, "Wrong group type on line " + 
		     to_string(_linec) );	
    return( _group_int64 );
}


int MyDXFFile::read_group( void )
{
    char buf[257];
    char *endptr;

    if( !_istr.is_open() )
	throw Error( ERROR_LOCATION, "No open file" );

    if( _ascii ) {

	// Read group code
	_linec++;
	_istr.getline( buf, 256 );
	if( _istr.eof() ) {
	    _group_code = -1;
	    return( -1 );
	}

	_group_code = strtol( buf, &endptr, 10 );
	if( endptr == buf ) {
	    throw Error( ERROR_LOCATION, "Error reading group code on line " + 
			 to_string(_linec) );
	}
	while( isspace( *endptr ) ) endptr++;
	if( *endptr != '\0' ) {
	    throw Error( ERROR_LOCATION, "Error reading group code on line " + 
			 to_string(_linec) );
	}

#if MYDXF_DEBUG >= 2
	std::cout << "Group code: " << _group_code << "\n";
#endif

	// Read group value
	_linec++;
	_istr.getline( buf, 256 );
	if( _istr.eof() ) {
	    throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			 to_string(_linec) + ", premature end of file" );
	}
	if( CODE_STRING(_group_code) ) {
	    // String value, remove possible leading end-of-line characters
	    _group_type = GROUP_TYPE_STRING;
	    int a = strlen( buf )-1;
	    while( a > 0 && (buf[a] == '\n' || buf[a] == '\r' || buf[a] == '\f') )
		a--;
	    buf[a+1] = '\0';
	    _group_string = buf;
#if MYDXF_DEBUG >= 2
	std::cout << "Group value (string): \'" << _group_string << "\'\n";
#endif

	} else if( CODE_BOOL(_group_code) ) {
	    // bool
	    int value = strtol( buf, &endptr, 10 );
	    if( endptr == buf ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    while( isspace( *endptr ) ) endptr++;
	    if( *endptr != '\0' ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    _group_type = GROUP_TYPE_BOOL;
	    _group_bool = value;
#if MYDXF_DEBUG >= 2
	std::cout << "Group value (bool): " << _group_bool << "\n";
#endif

	} else if( CODE_INT8(_group_code) ) {
	    // 8-bit integer value
	    int value = strtol( buf, &endptr, 10 );
	    if( endptr == buf ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    while( isspace( *endptr ) ) endptr++;
	    if( *endptr != '\0' ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    _group_type = GROUP_TYPE_INT8;
	    _group_int8 = value;
#if MYDXF_DEBUG >= 2
	std::cout << "Group value (int8): " << _group_int8 << "\n";
#endif

	} else if( CODE_INT16(_group_code) ) {
	    // 16-bit integer value
	    int value = strtol( buf, &endptr, 10 );
	    if( endptr == buf ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    while( isspace( *endptr ) ) endptr++;
	    if( *endptr != '\0' ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    _group_type = GROUP_TYPE_INT16;
	    _group_int16 = value;
#if MYDXF_DEBUG >= 2
	std::cout << "Group value (int16): " << _group_int16 << "\n";
#endif

	} else if( CODE_INT32(_group_code) ) { 
	    // 32-bit integer value
	    int value = strtol( buf, &endptr, 10 );
	    if( endptr == buf ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    while( isspace( *endptr ) ) endptr++;
	    if( *endptr != '\0' ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    _group_type = GROUP_TYPE_INT32;
	    _group_int32 = value;
#if MYDXF_DEBUG >= 2
	std::cout << "Group value (int32): " << _group_int32 << "\n";
#endif

	} else if( CODE_INT64(_group_code) ) { 
	    // 64-bit integer value
	    int64_t value = strtoll( buf, &endptr, 10 );
	    if( endptr == buf ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    while( isspace( *endptr ) ) endptr++;
	    if( *endptr != '\0' ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    _group_type = GROUP_TYPE_INT64;
	    _group_int64 = value;
#if MYDXF_DEBUG >= 2
	std::cout << "Group value (int64): " << _group_int64 << "\n";
#endif

	} else if( CODE_DOUBLE(_group_code) ) { 
	    // Double
	    double value = strtod( buf, &endptr );
	    if( endptr == buf ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    while( isspace( *endptr ) ) endptr++;
	    if( *endptr != '\0' ) {
		throw Error( ERROR_LOCATION, "Error reading group value on line " + 
			     to_string(_linec) );
	    }
	    _group_type = GROUP_TYPE_DOUBLE;
	    _group_double = value;
#if MYDXF_DEBUG >= 2
	std::cout << "Group value (double): " << _group_double << "\n";
#endif

	} else {
	    throw Error( ERROR_LOCATION, "Unknown group code " + to_string(_group_code)
			 + " on line " + to_string(_linec) );	
	}
	
    } else { // Read binary

	throw Error( ERROR_LOCATION, "Binary DXF file format unsupported" );

    }

    return( _group_code );
}


void MyDXFFile::write_group( int code, const char *data )
{
    if( !_ostr.is_open() )
	throw Error( ERROR_LOCATION, "No open output file" );

    if( _ascii ) {

	// Check code
	if( !CODE_STRING(code) )
	    throw Error( ERROR_LOCATION, "Incorrect code for data type" );

	_ostr << code << "\n";
	_ostr << data << "\n";

    } else {

	throw Error( ERROR_LOCATION, "Binary DXF file format unsupported" );

    }
}

void MyDXFFile::write_group( int code, double data )
{
    if( !_ostr.is_open() )
	throw Error( ERROR_LOCATION, "No open output file" );

    if( _ascii ) {

	// Check code
	if( !CODE_DOUBLE(code) )
	    throw Error( ERROR_LOCATION, "Incorrect code for data type" );

	_ostr << code << "\n";
	_ostr << std::scientific << std::setprecision(9) << data << "\n";

    } else {

	throw Error( ERROR_LOCATION, "Binary DXF file format unsupported" );

    }
}

void MyDXFFile::write_group( int code, bool data )
{
    if( !_ostr.is_open() )
	throw Error( ERROR_LOCATION, "No open output file" );

    if( _ascii ) {

	// Check code
	if( !CODE_BOOL(code) )
	    throw Error( ERROR_LOCATION, "Incorrect code for data type" );

	_ostr << code << "\n";
	_ostr << data << "\n";

    } else {

	throw Error( ERROR_LOCATION, "Binary DXF file format unsupported" );

    }
}

void MyDXFFile::write_group( int code, int8_t data )
{
    if( !_ostr.is_open() )
	throw Error( ERROR_LOCATION, "No open output file" );

    if( _ascii ) {

	// Check code
	if( !CODE_INT8(code) )
	    throw Error( ERROR_LOCATION, "Incorrect code for data type" );

	_ostr << code << "\n";
	_ostr << (int)data << "\n";

    } else {

	throw Error( ERROR_LOCATION, "Binary DXF file format unsupported" );

    }
}

void MyDXFFile::write_group( int code, int16_t data )
{
    if( !_ostr.is_open() )
	throw Error( ERROR_LOCATION, "No open output file" );

    if( _ascii ) {

	// Check code
	if( !CODE_INT16(code) )
	    throw Error( ERROR_LOCATION, "Incorrect code for data type" );

	_ostr << code << "\n";
	_ostr << data << "\n";

    } else {

	throw Error( ERROR_LOCATION, "Binary DXF file format unsupported" );

    }
}

void MyDXFFile::write_group( int code, int32_t data )
{
    if( !_ostr.is_open() )
	throw Error( ERROR_LOCATION, "No open output file" );

    if( _ascii ) {

	// Check code
	if( !CODE_INT32(code) )
	    throw Error( ERROR_LOCATION, "Incorrect code for data type" );

	_ostr << code << "\n";
	_ostr << data << "\n";

    } else {

	throw Error( ERROR_LOCATION, "Binary DXF file format unsupported" );

    }
}

void MyDXFFile::write_group( int code, int64_t data )
{
    if( !_ostr.is_open() )
	throw Error( ERROR_LOCATION, "No open output file" );

    if( _ascii ) {

	// Check code
	if( !CODE_INT64(code) )
	    throw Error( ERROR_LOCATION, "Incorrect code for data type" );

	_ostr << code << "\n";
	_ostr << data << "\n";

    } else {

	throw Error( ERROR_LOCATION, "Binary DXF file format unsupported" );

    }
}



void MyDXFFile::debug_print( std::ostream &os ) const
{
    os << "DXF File debug print:\n\n";
    
    if( _header )
	_header->debug_print( os );
    if( _tables )
	_tables->debug_print( os );
    if( _blocks )
	_blocks->debug_print( os );
    if( _entities )
	_entities->debug_print( os );
}


void MyDXFFile::debug_print_format( std::ostream &os, 
				    const std::string &fieldname, 
				    const std::string &val )
{
    os << "  " 
       << std::setw(DXF_DEBUG_PRINT_FIELD_WIDTH) << std::left << fieldname 
       << " = \'" << val << "\'\n";
}


void MyDXFFile::debug_print_format( std::ostream &os, 
				    const std::string &fieldname, 
				    double val )
{
    os << "  " 
       << std::setw(DXF_DEBUG_PRINT_FIELD_WIDTH) << std::left << fieldname 
       << " = " << val << "\n";
}


void MyDXFFile::debug_print_format( std::ostream &os, 
				    const std::string &fieldname, 
				    int val )
{
    os << "  " 
       << std::setw(DXF_DEBUG_PRINT_FIELD_WIDTH) << std::left << fieldname 
       << " = " << val << "\n";
}


void MyDXFFile::debug_print_format( std::ostream &os, 
				    const std::string &fieldname, 
				    const Vec3D &val )
{
    os << "  " 
       << std::setw(DXF_DEBUG_PRINT_FIELD_WIDTH) << std::left << fieldname 
       << " = (" << val[0] << ", " << val[1] << ", " << val[2] << ")\n";
}
