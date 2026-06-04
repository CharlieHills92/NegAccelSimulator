/*! \file cfifo.hpp
 *  \brief First-in first-out container with cyclic memory
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

#ifndef CFIFO_HPP
#define CFIFO_HPP 1


#include "error.hpp"


/*! \brief Cyclic memory first-in first-out container
 *
 *  The container has space for \a N class \a T type objects. The
 *  container is operated as first-in first-out. There is not need to
 *  empty the fifo as the oldest objects are overwritten when the fifo
 *  gets full.
 */
template <class T, int N>
class CFiFo {

    int      _size;
    int      _next;
    T        _data[N];

public:

    CFiFo() 
	: _size(0), _next(0) {
    }

    ~CFiFo() {

    }

    int size( void ) const {
	return( _size );
    }

    void clear( void ) {
	_size = 0;
	_next = 0;
    }

    void push( const T &x ) {
	_data[_next] = x;
	_next++;
	if( _next == N )
	    _next = 0;
	if( _size < N )
	    _size++;
    }

    /*! \brief Access \a i:th object from the input.
     *
     *  Object 0 is the newest one in the fifo. Throws an ErrorRange
     *  exception if trying to access non-existing object.
     */
    const T &operator[]( int i ) const {
	if( i < 0 )
	    throw( ErrorRange( ERROR_LOCATION, "Trying to access negative index" ) );
	else if( i >= _size )
	    throw( ErrorRange( ERROR_LOCATION, "Trying to access beyond fifo end" ) );
	int access = _next-i-1;
	if( access < 0 )
	    access += _size;
	return( _data[access] );
    };

};


#endif
