/*! \file readascii.hpp
 *  \brief ASCII file reading tool.
 */

/* Copyright (c) 2011,2015 Taneli Kalvas. All rights reserved.
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

#ifndef READASCII_HPP
#define READASCII_HPP 1


#include <vector>
#include <string>
#include <stdint.h>


/*! \brief Class for reading ASCII data files.
 *
 *  Reads ASCII data file formatted into \a N columns skipping empty
 *  lines and lines starting with a comment character (assignable,
 *  defaults to '#'). The number of columns may be given to read
 *  function, in which case it is checked to be correct or the number
 *  of columns may be determined from the start of the file. The
 *  number of columns must be fixed and may not change during the
 *  file. Data is stored into \a N double-type vectors, which can be
 *  read from the object after reading a datafile.
 */
class ReadAscii {

    std::string           _filename; /*!< \brief Name of read file. */
    char                  _comment;  /*!< \brief Comment character. */
    uint32_t              _N;        /*!< \brief Number of columns. */
    std::vector<double> **_data;     /*!< \brief Column data vectors. */

    void read_data_line( const std::string &str, int linec );

public:

    /*! \brief Constructor for empty class.
     */
    ReadAscii();

    /*! \brief Constructor for class from file.
     *
     *  Read ASCII data file from \a filename. If \a columns is -1 the
     *  number of data columns is determined from the file. If the
     *  number of columns is given the file is checked to have the
     *  columns.
     */
    ReadAscii( const std::string &filename, int columns = -1 );

    /*! \brief Destructor.
     */
    ~ReadAscii();

    /*! \brief Read ASCII data file.
     *
     *  Read ASCII data file from \a filename. If \a columns is -1 the
     *  number of data columns is determined from the file. If the
     *  number of columns is given the file is checked to have the
     *  columns.
     */
    void read( const std::string &filename, int columns = -1 );

    void set_comment_character( char c );

    /*! \brief Clear data.
     */
    void clear( void );

    /*! \brief Return number of columns in data.
     */
    uint32_t columns( void ) const {
	return( _N );
    }

    /*! \brief Return number of rows in data.
     */
    uint32_t rows( void ) const {
	if( _N > 0 )
	    return( _data[0]->size() );
	return( 0 );
    }

    /*! \brief Return const reference to the vector containing column \a i.
     */
    const std::vector<double> &operator[]( uint32_t i ) const;
};


#endif

