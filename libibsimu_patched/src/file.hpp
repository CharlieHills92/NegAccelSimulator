/*! \file file.hpp
 *  \brief Bindary file writing and reading tools.
 */

/* Copyright (c) 2005-2010 Taneli Kalvas. All rights reserved.
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

#ifndef FILE_HPP
#define FILE_HPP 1


#include <iostream>
#include <stdint.h>

#define FILEID_NULL             0

#define FILEID_GEOMETRY      1001

#define FILEID_FUNCSOLID     2001
#define FILEID_CSGSOLID      2002
#define FILEID_DXFSOLID      2003
#define FILEID_STLSOLID      2004

#define FILEID_PARTICLEDB2D  3001
#define FILEID_PARTICLEDBCYL 3002
#define FILEID_PARTICLEDB3D  3003

#define FILEID_SCALARFIELD   4001

#define FILEID_VECTORFIELD   5001


/* **************** *
 * Write            *
 * **************** */


/*! \brief Write int8_t \a value into stream \a os.
 */
void write_int8( std::ostream &os, int8_t value );


/*! \brief Write int16_t \a value into stream \a os.
 */
void write_int16( std::ostream &os, int16_t value );


/*! \brief Write int32_t \a value into stream \a os.
 */
void write_int32( std::ostream &os, int32_t value );


/*! \brief Write uint32_t \a value into stream \a os.
 */
void write_uint32( std::ostream &os, uint32_t value );


/*! \brief Write double \a value into stream \a os.
 */
void write_double( std::ostream &os, double value );


/*! \brief Write data block \a data of length \a len bytes into stream
 *  \a os in compressed form.
 */
void write_compressed_block( std::ostream &os, uint32_t len, const int8_t *data );


/* **************** *
 * Read             *
 * **************** */


/*! \brief Read int8_t from stream \a is.
 */
int8_t read_int8( std::istream &is );


/*! \brief Read int16_t from stream \a is.
 */
int16_t read_int16( std::istream &is );


/*! \brief Read int32_t from stream \a is.
 */
int32_t read_int32( std::istream &is );


/*! \brief Read uint32_t from stream \a is.
 */
uint32_t read_uint32( std::istream &is );


/*! \brief Readd double from stream \a is.
 */
double read_double( std::istream &is );


/*! \brief Read compressed data block of length \a len bytes from stream \a is.
 */
uint32_t read_compressed_block( std::istream &is, uint32_t len, int8_t *dest );


#endif

