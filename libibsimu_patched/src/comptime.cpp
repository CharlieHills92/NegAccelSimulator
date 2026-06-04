/*! \file comptime.cpp
 *  \brief Compatible time functions
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


#include "comptime.hpp"


#if defined(WIN32) || defined(__MINGW32__)
#include <windows.h>
#endif


#if !defined(HAVE_GETTIMEOFDAY)
#if defined(WIN32) || defined(__MINGW32__)
#define EPOCHFILETIME (116444736000000000LL)
int ibs_gettimeofday( struct timeval *tv, struct timezone *tz )
{
    FILETIME        ft;
    LARGE_INTEGER   li;
    __int64         t;

    if( tv ) {
        GetSystemTimeAsFileTime(&ft);
        li.LowPart  = ft.dwLowDateTime;
        li.HighPart = ft.dwHighDateTime;
        t  = li.QuadPart;       /* In 100-nanosecond intervals */
        t -= EPOCHFILETIME;     /* Offset to the Epoch time */
        t /= 10;                /* In microseconds */
        tv->tv_sec  = (long)(t / 1000000);
        tv->tv_usec = (long)(t % 1000000);
    }

    return( 0 );
}
#else
int ibs_gettimeofday( struct timeval *tv, struct timezone *tz )
{
    tv->tv_sec = (long)time(NULL);
    tv->tv_usec = 0;
    return( 0 );
}
#endif
#else
int ibs_gettimeofday( struct timeval *tv, struct timezone *tz )
{
    return( gettimeofday( tv, tz ) );
}
#endif



#if !defined(HAVE_CLOCK_GETTIME)
#if defined(WIN32) || defined(__MINGW32__)
#define EPOCHFILETIME (116444736000000000LL)
int ibs_clock_gettime( clockid_t clk_id, struct timespec *tp )
{
    FILETIME        ft;
    LARGE_INTEGER   li;
    __int64         t;

    if( tp ) {
        GetSystemTimeAsFileTime(&ft);
        li.LowPart  = ft.dwLowDateTime;
        li.HighPart = ft.dwHighDateTime;
        t  = li.QuadPart;       /* In 100-nanosecond intervals */
        t -= EPOCHFILETIME;     /* Offset to the Epoch time */
        tp->tv_sec  = (long)(t / 10000000);
        tp->tv_nsec = 100L*((long)(t % 10000000));
    }

    return( 0 );
}
#else

#if !defined(HAVE_GETTIMEOFDAY)
int ibs_clock_gettime( clockid_t clk_id, struct timespec *tp )
{
    struct timeval tv;
    gettimeofday( &tv, NULL );
    tp->tv_sec  = tv.tv_sec;
    tp->tv_nsec = 1000L*tv.tv_usec;
    return( 0 );
}
#else
/* using time() in 1 sec intervals */
int ibs_clock_gettime( clockid_t clk_id, struct timespec *tp )
{
    tp->tv_sec = (long)time(NULL);
    tp->tv_nsec = 0;
    return( 0 );
}
#endif

#endif
#else
int ibs_clock_gettime( clockid_t clk_id, struct timespec *tp )
{
    return( clock_gettime( clk_id, tp ) );
}
#endif

