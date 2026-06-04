/*! \file file.cpp
 *  \brief Bindary file writing and reading tools.
 */

/* Copyright (c) 2005-2012 Taneli Kalvas. All rights reserved.
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


#include "file.hpp"
#include <zlib.h>
#include "error.hpp"


//#define DEBUG_FILE 1


void write_int8( std::ostream &s, int8_t value )
{
#ifdef DEBUG_FILE
    std::cout << "write_int8: " << (int)value << "\n";
#endif
    s.write( (const char *)(&value), sizeof(int8_t) );
}


void write_int16( std::ostream &s, int16_t value )
{
#ifdef DEBUG_FILE
    std::cout << "write_int16: " << value << "\n";
#endif
    s.write( (const char *)(&value), sizeof(int16_t) );
}


void write_int32( std::ostream &s, int32_t value )
{
#ifdef DEBUG_FILE
    std::cout << "write_int32: " << value << "\n";
#endif
    s.write( (const char *)(&value), sizeof(int32_t) );
}


void write_uint32( std::ostream &s, uint32_t value )
{
#ifdef DEBUG_FILE
    std::cout << "write_uint32: " << value << "\n";
#endif
    s.write( (const char *)(&value), sizeof(uint32_t) );
}


void write_double( std::ostream &s, double value )
{
#ifdef DEBUG_FILE
    std::cout << "write_double: " << value << "\n";
#endif
    s.write( (const char *)(&value), sizeof(double) );
}


void write_compressed_block( std::ostream &s, uint32_t len, const int8_t *data )
{
    uLongf  buflen = uLongf(len*1.1)+12 ;
    Bytef  *buf = new Bytef[buflen];
    compress( buf, &buflen, (Bytef *)data, (uLong)len );
    write_uint32( s, len );
    write_uint32( s, buflen );
    s.write( (const char *)buf, buflen );
#ifdef DEBUG_FILE
    std::cout << "write_compressed_block: datalen = " << len
	      << ", compressedlen = " << buflen << "\n";
#endif
    delete [] buf;
}


int8_t read_int8( std::istream &s )
{
    int8_t value;
    s.read( (char *)(&value), sizeof(int8_t) );
#ifdef DEBUG_FILE
    std::cout << "read_int8: " << (int)value << "\n";
#endif
    return( value );
}


int16_t read_int16( std::istream &s )
{
    int16_t value;
    s.read( (char *)(&value), sizeof(int16_t) );
#ifdef DEBUG_FILE
    std::cout << "read_int16: " << value << "\n";
#endif
    return( value );
}


int32_t read_int32( std::istream &s )
{
    int32_t value;
    s.read( (char *)(&value), sizeof(int32_t) );
#ifdef DEBUG_FILE
    std::cout << "read_int32: " << value << "\n";
#endif
    return( value );
}


uint32_t read_uint32( std::istream &s )
{
    uint32_t value;
    s.read( (char *)(&value), sizeof(uint32_t) );
#ifdef DEBUG_FILE
    std::cout << "read_uint32: " << value << "\n";
#endif
    return( value );
}


double read_double( std::istream &s )
{
    double value;
    s.read( (char *)(&value), sizeof(double) );
#ifdef DEBUG_FILE
    std::cout << "read_double: " << value << "\n";
#endif
    return( value );
}


uint32_t read_compressed_block( std::istream &s, uint32_t len, int8_t *dest )
{
    uint32_t datalen;
    uint32_t compressedlen;
    datalen = read_uint32( s );
    compressedlen = read_uint32( s );
#ifdef DEBUG_FILE
    std::cout << "read_compressed_block: datalen = " << datalen 
	      << ", compressedlen = " << compressedlen << "\n";
#endif

    if( datalen > len )
	throw( Error( ERROR_LOCATION, "compressed data length longer than expected (" + 
		      to_string(datalen) + " > " + to_string(len) ) );

    Bytef  *buf = new Bytef [compressedlen];
    uLongf  tmp = datalen;
    s.read( (char *)buf, compressedlen );
    uncompress( (Bytef *)dest, &tmp, buf, compressedlen );
    delete [] buf;

#ifdef DEBUG_FILE
    std::cout << "done\n";
#endif
    return( datalen );
}

