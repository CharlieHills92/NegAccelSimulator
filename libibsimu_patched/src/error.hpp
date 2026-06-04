/*! \file error.hpp
 *  \brief %Error classes and handling
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

#ifndef ERROR_HPP
#define ERROR_HPP 1


#include "config.h"
#include <string>
#include <cstring>
#include <stdlib.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <errno.h>
#include <signal.h>


/*! \brief Function for converting a type to string.
 */
template <class T>
inline std::string to_string( const T& t )
{
    std::stringstream ss;
    ss << t;
    return( ss.str() );
}


/*! \brief Function for converting a type to string.
 */
template <class T>
inline std::wstring to_wstring( const T& t )
{
    std::wstringstream ss;
    ss << t;
    return( ss.str() );
}


/*! \brief Macro for setting error location when throwing errors.
 */
#define ERROR_LOCATION ErrorLocation( __FILE__, __LINE__, __func__ )


/*! \brief %Error location class.
 *
 *  Container to store the location (source file name, line number and
 *  function name) where the error happened. Used for debugging
 *  purposes. Macro ERROR_LOCATION is defined for convenient use of
 *  class.
 */
class ErrorLocation {

    const char *_file;
    int         _line;
    const char *_func;

public:

    /*! \brief Default constructor for error location.
     *
     *  Stores a null location.
     */
    ErrorLocation();

    /*! \brief Constructor for setting error location.
     *
     *  This constructor is conveniently called with ERROR_LOCATION macro.
     */
    ErrorLocation( const char *file, int line, const char *func );

    /*! \brief Return file name of location.
     */
    std::string file( void );

    /*! \brief Return line number of location.
     */
    int line( void );

    /*! \brief Return function name of location.
     */
    std::string func( void );
};


/*! \brief Exception backtrace.
 *
 *  Saves the backtrace when constructed. Uses GNU extensions and
 *  therefore only works when compiled with GNU system.
 */
class ExceptionTracer {

    void  *_traceaddress[25];
    int    _tracecount;

public:

    /*! \brief Default constructor for exception tracer. Saves the
     *  backtrace of the program at this location for printing it when
     *  the error is caught.
     */
    ExceptionTracer();

    /*! \brief Print the backtrace to \a os.
     */
    void print_trace( std::ostream &os );
};


/*! \brief Basic error class.
 */
class Error : public ExceptionTracer {

    ErrorLocation  _loc;

protected:
    
    std::string    _error_str;

public:

    /*! \brief Default constructor for error class.
     */
    Error();

    /*! \brief Constructor for error class with error message.
     */
    Error( const std::string &str );

    /*! \brief Constructor for error class with location information.
     */
    Error( const ErrorLocation &loc );

    /*! \brief Constructor for error class with location information
     *  and error message.
     */
    Error( const ErrorLocation &loc, const std::string &str );

    /*! \brief Return error message.
     */
    std::string get_error_message( void );

    /*! \brief Print a standard error message to \a os.
     *
     *  Print error message string and the error location is source
     *  code. Also prints the exception back trace if \a print_trace
     *  is true.
     */
    void print_error_message( std::ostream &os, bool traceprint = true );
};


/*! \brief %Error class for memory allocation errors.
 */
class ErrorNoMem : public Error {

public:

    /*! \brief Constructor for memory allocation error with standard error message.
     *
     *  The error message is "memory allocation error".
     */
    ErrorNoMem( const ErrorLocation &loc );

    /*! \brief Constructor for memory allocation error with custom error message.
     */
    ErrorNoMem( const ErrorLocation &loc, const std::string &str );
};


/*! \brief %Error class for C-style errno errors.
 */
class ErrorErrno : public Error {

    int _ierrno;

public:

    /*! \brief Constructor for errno based error with standard error
     *  message from errno database.
     */
    ErrorErrno( const ErrorLocation &loc );
};


/*! \brief %Error class to use if requested feature is unimplemented.
 */
class ErrorUnimplemented : public Error {

public:

   /*! \brief Constructor for unimplemented feature error with standard error message.
     *
     *  The error message is "feature unimplemented".
     */
    ErrorUnimplemented( const ErrorLocation &loc );

    /*! \brief Constructor for unimplemented feature error with custom error message.
     */
    ErrorUnimplemented( const ErrorLocation &loc, const std::string &str );
};


/*! \brief %Error class to use if impossible things happen
 */
class ErrorAssert : public Error {

public:

   /*! \brief Constructor for assert error with standard error message.
     *
     *  The error message is "assertion failed".
     */
    ErrorAssert( const ErrorLocation &loc );

    /*! \brief Constructor for assert error with custom error message.
     */
    ErrorAssert( const ErrorLocation &loc, const std::string &str );
};


/*! \brief %Error class for dimension mismatch errors.
 */
class ErrorDim : public Error {

public:

    /*! \brief Constructor for dimension mismatch error with standard error message.
     *
     *  The error message is "dimension mismatch".
     */
    ErrorDim( const ErrorLocation &loc );

    /*! \brief Constructor for dimension mismatch error with custom error message.
     */
    ErrorDim( const ErrorLocation &loc, const std::string &str );
};


/*! \brief %Error class for index range checking errors.
 */
class ErrorRange : public Error {

public:

    /*! \brief Constructor for error message for two dimensional indexing error.
     *
     *  The index \a i is supposed to be smaller than \a n and \a j
     *  smaller than \a m.
     */
    ErrorRange( const ErrorLocation &loc, 
		uint32_t i, uint32_t n, 
		uint32_t j, uint32_t m );

    /*! \brief Constructor for error message for one dimensional indexing error.
     *
     *  The index \a i is supposed to be smaller than \a n.
     */
    ErrorRange( const ErrorLocation &loc, uint32_t i, uint32_t n );

    /*! \brief Constructor for error message for one dimensional indexing error.
     *
     *  The index \a i is supposed to be smaller than \a n. Additional
     *  message string.
     */
    ErrorRange( const ErrorLocation &loc, uint32_t i, uint32_t n, const std::string &str );

    /*! \brief Constructor for error message for indexing error.
     */
    ErrorRange( const ErrorLocation &loc, const std::string &str );
};


/*! \brief Signal handler
 */
class SignalHandler {

public:

#if defined(_GNU_SOURCE) && defined(HAVE_SIGINFO_T)
    /*! \brief Signal handler function for SIGSEGV.
     */
    static void signal_handler_SIGSEGV( int signum, siginfo_t *info, void *ptr );

    /*! \brief Signal handler function for SIGTERM.
     */
    static void signal_handler_SIGTERM( int signum, siginfo_t *info, void *ptr );
#else
    /*! \brief Signal handler function for SIGTERM.
     */
    static void signal_handler_SIGTERM( int signum );
#endif
};


#endif

