/*! \file stl_solid.cpp
 *  \brief Stereolithography CAD file handling
 */

/* Copyright (c) 2011-2012,2016 Taneli Kalvas. All rights reserved.
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

#include <cstddef>
#include <limits>
#include <sstream>
#include <ctype.h>
#include <string.h>
#include "stlfile.hpp"
#include "ibsimu.hpp"


#define DEBUG_STL 1


bool ciscomp( const char *str1, const char *str2, size_t n )
{
    size_t str1len = strlen( str1 );
    size_t str2len = strlen( str2 );

    if( str1len < n ) {
	if( str1len != str2len )
	    return( false );
	n = str1len;
    }
    if( str2len < n )
	return( false );

    for( size_t i = 0; i < n ; i++ ) {
        if( tolower(str1[i]) != tolower(str2[i]) )
            return( false );
    }
    return( true );
}


/* ******************** *
 * Triangle             *
 * ******************** */


void STLFile::Triangle::read_binary_float_vector( Vec3D &x, std::ifstream &ifstr )
{
    for( size_t a = 0; a < 3; a++ ) {
	float r;
	ifstr.read( (char *)&r, 4 );
	x[a] = r;
    }
}


void STLFile::Triangle::read_ascii_float_vector( Vec3D &x, const char *buf, 
						 const std::string &filename, int linec )
{
    char *endptr;

    // Skip preceding whitespace
    while( isspace( *buf ) )
	buf++;

    for( size_t i = 0; i < 3; i++ ) {
	x[i] = strtod( buf, &endptr );
	if( endptr == buf ) {
	    throw( Error( ERROR_LOCATION, "unexpected end of line reading file \'" 
			  + filename + "\' on line " + to_string(linec) ) );
	}
	buf = endptr;
	while( isspace( *buf ) )
	    buf++;
    }

}


STLFile::Triangle::Triangle( std::ifstream &ifstr )
{
    read_binary_float_vector( _normal, ifstr );
    read_binary_float_vector( _p[0], ifstr );
    read_binary_float_vector( _p[1], ifstr );
    read_binary_float_vector( _p[2], ifstr );

    uint16_t attr;
    ifstr.read( (char *)&attr, 2 );
    _attr = attr;
}


STLFile::Triangle::Triangle( std::ifstream &ifstr, const char *buf, const std::string &filename, int &linec )
{
    if( !ciscomp( buf, "facet normal", 12 ) )
	throw( Error( ERROR_LOCATION, "Unexpected input on line " + to_string(linec) +
		      ", expecting \'facet normal\'." ) );
    buf += 12;
    read_ascii_float_vector( _normal, buf, filename, linec );

    // Tag "outer loop"
    std::string str;
    std::getline( ifstr, str );
    buf = str.c_str();
    linec++;

    // Skip whitespace
    while( isspace( *buf ) )
	buf++;

    if( !ciscomp( buf, "outer loop", 10 ) )
	throw( Error( ERROR_LOCATION, "Unexpected input on line " + to_string(linec) +
		      ", expecting \'outer loop\'." ) );

    // Read 3 vertices
    for( size_t i = 0; i < 3; i++ ) {
	std::getline( ifstr, str );
	buf = str.c_str();
	linec++;
	
	// Skip whitespace
	while( isspace( *buf ) )
	    buf++;
	
	if( !ciscomp( buf, "vertex", 6 ) )
	    throw( Error( ERROR_LOCATION, "Unexpected input on line " + to_string(linec) +
			  ", expecting \'vertex\'." ) );

	buf += 6;
	read_ascii_float_vector( _p[i], buf, filename, linec );
    }

    // Tag "endloop"
    std::getline( ifstr, str );
    buf = str.c_str();
    linec++;

    // Skip whitespace
    while( isspace( *buf ) )
	buf++;

    if( !ciscomp( buf, "endloop", 7 ) )
	throw( Error( ERROR_LOCATION, "Unexpected input on line " + to_string(linec) +
		      ", expecting \'endloop\'." ) );

    // Tag "endfacet"
    std::getline( ifstr, str );
    buf = str.c_str();
    linec++;

    // Skip whitespace
    while( isspace( *buf ) )
	buf++;

    if( !ciscomp( buf, "endfacet", 8 ) )
	throw( Error( ERROR_LOCATION, "Unexpected input on line " + to_string(linec) +
		      ", expecting \'endfacet\'." ) );


}


STLFile::Triangle::~Triangle()
{

}


uint16_t STLFile::Triangle::attr( void ) const
{
    return( _attr );
}


const Vec3D &STLFile::Triangle::normal( void ) const
{
    return( _normal );
}


void STLFile::Triangle::debug_print( std::ostream &os ) const
{
    os << "**Triangle\n";
    os << "  normal = " << _normal << "\n";
    os << "  p1     = " << _p[0] << "\n";
    os << "  p2     = " << _p[1] << "\n";
    os << "  p3     = " << _p[2] << "\n";
    os << "  attr   = " << _attr << "\n";
}



/* ******************** *
 * STLFile              *
 * ******************** */


void STLFile::read_binary( std::ifstream &ifstr )
{
    // Skip header
    char buf[80];
    ifstr.read( buf, 80 );

    // Read number of triangles
    uint32_t tcount;
    ifstr.read( (char *)(&tcount), 4 );

    // Check if sensible
    if( tcount > 10000000 )
	throw( ErrorUnimplemented( ERROR_LOCATION, "Too high number of triangles" ) );
    _tri.reserve( tcount );

    // Read triangles
    for( uint32_t a = 0; a < tcount; a++ )
	_tri.push_back( Triangle( ifstr ) );
}


void STLFile::read_ascii( std::ifstream &ifstr )
{
    // Read first line (header)
    std::string header;
    std::getline( ifstr, header );

    // Read line-by-line
    int linec = 0;
    while( !ifstr.eof() ) {
        
        std::string str;
        std::getline( ifstr, str );
	const char *buf = str.c_str();
        linec++;

	// Skip whitespace
	while( isspace( *buf ) )
	    buf++;

	if( ciscomp( buf, "endsolid", 8 ) )
	    break;

	_tri.push_back( Triangle( ifstr, buf, _filename, linec ) );
    }
}


STLFile::STLFile( const std::string &filename,
		  double vertex_matching_eps, 
		  double signed_volume_eps )
    : _filename(filename), _solid(vertex_matching_eps, signed_volume_eps)
{
    ibsimu.message( 1 ) << "Reading STL-file \'" << filename << "\'\n";
    ibsimu.inc_indent();

    std::ifstream ifstr( filename.c_str(), std::ios_base::binary );
    if( !ifstr.good() )
	throw( Error( ERROR_LOCATION, "Couldn't open file \'" + filename + "\'" ) );

    /* Check if ascii
     * Ascii files start with "solid XXX\n", on first line where XXX is free
     * form text. Second line starts with "facet normal" with possible leading whitespace
     * Binary files have free form 80 byte header
     */
    _ascii = false;
    char buf[1024];
    ifstr.read( buf, 1024 );
    std::streamsize c = ifstr.gcount();
    if( c > 6 && !strncasecmp( buf, "solid ", 6 ) ) {
	// Might be ascii, seek next line
	int a = 6;
	while( a < c && buf[a] != '\n' ) a++;
	while( a < c && isspace(buf[a]) ) a++;
	if( a+12 < c && !strncasecmp( &buf[a], "facet normal", 12 ) )
	    _ascii = true;
    }

    // Read triangle data
    ifstr.clear(); // Clear possible eofbit/failbit
    ifstr.seekg( 0 );
    if( _ascii )
	read_ascii( ifstr );
    else
	read_binary( ifstr );
    ifstr.close();

    // Convert to vertex triangle data
    build_vtriangle_data();
    _tri.clear();

    // Report
    ibsimu.message(1) << "vertexc = " << _solid.vertexc() << "\n";
    ibsimu.message(1) << "trianglec = " << _solid.trianglec() << "\n";

    // Check
    _solid.check_data();

    ibsimu.dec_indent();
}


STLFile::~STLFile()
{
    
}


void STLFile::save( const std::string &filename, bool ascii ) const
{
    if( !ascii )
	throw( ErrorUnimplemented( ERROR_LOCATION, "Binary file save not implemented" ) );	

    std::ofstream os( filename.c_str(), std::ios_base::binary );
    if( !os.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\' for writing" ) );

    os << "solid " << filename << "\n";

    for( uint32_t a = 0; a < _solid.trianglec(); a++ ) {    

	os << "facet normal 0 0 0\n";
	os << "outer loop\n";
	const VTriangle &tri = _solid.triangle(a);
	for( size_t b = 0; b < 3; b++ ) {    
	    uint32_t v = tri[b];
	    os << "vertex " <<  _solid.vertex(v)[0]
	       << " " << _solid.vertex(v)[1]
	       << " " << _solid.vertex(v)[2] << "\n";
	}
	os << "endloop\n";
	os << "endfacet\n";
    }
    
    os << "endsolid\n";

    os.close();
}


void STLFile::build_vtriangle_data( void )
{
    ibsimu.message(1) << "Making vertex connections\n";

    // Go through all triangles and add coordinates to vertex list if
    // they are not already there. Add vertex triangles in the same
    // time using made vertices.
    _solid.clear();
    for( size_t a = 0; a < _tri.size(); a++ )
	_solid.add_triangle( _tri[a][0], _tri[a][1], _tri[a][2] );

    _solid.prepare_for_inside();
}


void STLFile::debug_print( std::ostream &os ) const
{
    os << "**STLFile\n";

    os << "  vertexc = " << _solid.vertexc() << "\n";
    for( size_t a = 0; a < _solid.vertexc(); a++ )
	os << "  vertex[" << a << "] = " << _solid.vertex(a) << "\n";

    os << "  vtric = " << _solid.trianglec() << "\n";
    for( size_t a = 0; a < _solid.trianglec(); a++ ) {
	const VTriangle &tri = _solid.triangle(a);
	os << "  vtriangle[" << a << "]:\n";
	os << "    v1 = " << tri[0] << "\n";
	os << "    v2 = " << tri[1] << "\n";
	os << "    v3 = " << tri[2] << "\n";
    }
}


