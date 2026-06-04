/*! \file readascii.cpp
 *  \brief ASCII file reading tool.
 */

/* Copyright (c) 2011-2012,2015 Taneli Kalvas. All rights reserved.
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
#include <fstream>
#include "readascii.hpp"
#include "error.hpp"


ReadAscii::ReadAscii()
    : _comment('#'), _N(0), _data(NULL)
{
}


ReadAscii::ReadAscii( const std::string &filename, int columns )
    : _comment('#'), _N(0), _data(NULL)
{
    read( filename, columns );
}


ReadAscii::~ReadAscii()
{
    clear();
}


void ReadAscii::set_comment_character( char c )
{
    _comment = c;
}


void ReadAscii::read_data_line( const std::string &str, int linec )
{
    // Parse line
    const char *ptr = str.c_str();
	    
    // Skip leading white space
    while( isspace(*ptr) ) ptr++;

    // Check if comment
    if( *ptr == _comment )
	return;

    // Check if line contained only white space
    if( *ptr == '\n' || *ptr == '\r' || *ptr == '\0' )
	return;

    // Read line data
    for( uint32_t i = 0; i < _N; i++ ) {

	char *endptr;
	double val = strtod( ptr, &endptr );
	if( endptr == ptr ) {
	    if( *ptr == '\n' || *ptr == '\r' || *ptr == '\0' )
		throw( Error( ERROR_LOCATION, "unexpected end of line reading file \'" 
			      + _filename + "\' on line " + to_string(linec) 
			      + ". Missing data column on line?" ) );
	    throw( Error( ERROR_LOCATION, "unexpected input reading file \'" 
			  + _filename + "\' on line " + to_string(linec) ) );
	}
	_data[i]->push_back( val );
	ptr = endptr;

	// Skip white space
	while( isspace(*ptr) ) ptr++;
    }

    // Check if line done
    if( *ptr != '\n' && *ptr != '\r' && *ptr != '\0' )
	throw( Error( ERROR_LOCATION, "unexpected input reading file \'" 
		      + _filename + "\' on line " + to_string(linec)
		      + ". More input than expected. Variable number of columns in file?") );
}


void ReadAscii::read( const std::string &filename, int columns )
{
    // Clear old data
    for( uint32_t i = 0; i < _N; i++ )
	delete _data[i];
    delete _data;

    _filename = filename;
    std::ifstream fin( filename.c_str() );
    if( !fin.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\'" ) );

    // Analyze first valid record of file if number of columns is not known
    if( columns == -1 ) {
	int linec = 0;
	std::string str;
	const char *ptr;
	while( !fin.eof() ) {

	    // Read line
	    std::getline( fin, str );
	    linec++;

	    // Parse line
	    ptr = str.c_str();
	    
	    // Skip leading white space
	    while( isspace(*ptr) ) ptr++;

	    // Check if comment
	    if( *ptr == _comment )
		continue;

	    // Check if line contained only white space
	    if( *ptr == '\n' || *ptr == '\r' || *ptr == '\0' )
		continue;

	    break;
	}
	if( fin.eof() ) {
	    // EOF before meaningful input, return empty file
	    _N = 0;
	    fin.close();
	    return;
	}
	
	// Count columns
	columns = 0;
	while( 1 ) {
	    char *endptr;
	    double val = strtod( ptr, &endptr );
	    if( val == 0.0 && endptr == ptr )
		throw( Error( ERROR_LOCATION, "unexpected input reading file \'" 
			      + filename + "\' on line " + to_string(linec) ) );
	    ptr = endptr;
	    columns++;

	    // Skip white space
	    while( isspace(*ptr) ) ptr++;
	    
	    // Check if line done
	    if( *ptr == '\n' || *ptr == '\r' || *ptr == '\0' )
		break;
	}
    }

    // Reserve space for data
    _N = columns;
    _data = new std::vector<double> *[_N];
    for( uint32_t i = 0; i < _N; i++ )
	_data[i] = new std::vector<double>;

    // Read data
    int linec = 0;
    fin.clear();
    fin.seekg( 0 );
    while( !fin.eof() ) {
	
	// Read line
	std::string str;
	std::getline( fin, str );
	linec++;

	read_data_line( str, linec );
    }

    fin.close();

    // Most probably these data vectors are only used for temporary
    // storage and therefore they are not reserved to exact
    // size. Excess memory use does not matter in this case.
}


void ReadAscii::clear( void )
{
    for( uint32_t i = 0; i < _N; i++ )
	delete _data[i];
    delete _data;
    _data = NULL;
    _N = 0;
}


const std::vector<double> &ReadAscii::operator[]( uint32_t i ) const
{
    if( i < _N )
	return( *_data[i] );
    throw( ErrorRange( ERROR_LOCATION, i, _N ) );
}

