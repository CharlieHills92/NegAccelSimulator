/*! \file timer.hpp
 *  \brief %Timer for cputime and realtime
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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

#ifndef TIMER_HPP
#define TIMER_HPP 1


#include <iostream>
#include <time.h>
#include <sys/time.h>


/*! \brief Class for measuring code runtime in cpu time and realtime.
 */
class Timer {

    clock_t         _cpu1, _cpu2;
    struct timeval  _time1, _time2;

public:

    /*! \brief Construct and start timer.
     */
    Timer();

    /*! \brief Destruct timer.
     */
    ~Timer() {}

    /*! \brief Start timer.
     */
    void start( void );

    /*! \brief Stop timer.
     */
    void stop( void );
    
    /*! \brief Get real time elapsed in seconds.
     */
    double get_real_time( void ) const;

    /*! \brief Get cpu time elapsed in seconds.
     */
    double get_cpu_time( void ) const;

    /*! \brief Print report of elapsed time to stream.
     */    
    friend std::ostream &operator<<( std::ostream &os, const Timer &t );
};


#endif
