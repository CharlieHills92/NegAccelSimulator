/*! \file ibsimu.hpp
 *  \brief Ion Beam Simulator global settings
 */

/* Copyright (c) 2010-2012 Taneli Kalvas. All rights reserved.
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


#ifndef IBSIMU_HPP
#define IBSIMU_HPP 1


#include <stdint.h>
#include <stdarg.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>


class Timer;


/*! \brief Message type count.
 *
 *  Must match number of types of message_type_e.
 */
#define MSG_COUNT 5

/*! \brief Output type.
 */
enum message_type_e {
    MSG_VERBOSE = 0,
    MSG_WARNING,
    MSG_ERROR,
    MSG_DEBUG_GENERAL,
    MSG_DEBUG_DXF
};


/*! \brief Output type.
 */
enum rng_type_e {
    RNG_SOBOL = 0,
    RNG_MT
};


/*! \brief Main class for %IBSimu.
 *
 *  Used to store global settings. One instance of the class is
 *  initialized globally with the name \a ibsimu.
 */
class IBSimu 
{
    Timer        *_t;

    bool          _hello;
    uint32_t      _threadcount;
    rng_type_e    _rng;
    
    struct nullstream: std::ostream {
	struct nullbuf: std::streambuf {
	    int overflow( int c ) { return( traits_type::not_eof(c) ); }
	} m_sbuf;
	nullstream(): std::ios(&m_sbuf), std::ostream(&m_sbuf) {}
    };

    std::stringstream  _ss;
    nullstream         _ns;
    std::ofstream      _fo;
    std::ostream      *_os;             // Message output stream

    int32_t       _indent;
    int32_t       _message_threshold[MSG_COUNT];

    /*! \brief Copy constructor
     *
     *  Private. Prevent copying
     */
    IBSimu( const IBSimu &ibs );

    const IBSimu &operator=( const IBSimu &ibs ) { return( *this ); }

    size_t convert_stringstream_to_lines( const std::stringstream &ss,
					  size_t indentation,
					  std::string &target,
					  bool &lineunfinished );

public:

    /*! \brief Default constructor.
     */
    IBSimu();

    /*! \brief Default destructor.
     */
    ~IBSimu();

    /*! \brief Set message output to stream \a os.
     *
     *  Returns a reference to the old stream or nullstream if
     *  previous output file was a file opened by IBSimu.
     */
    std::ostream &set_message_output( std::ostream &os );

    /*! \brief Set message output to file \a filename.
     *
     *  Returns a reference to the old output stream or nullstream if
     *  previous output file was a file opened by IBSimu. If output
     *  stream is redefined after using this function, the file will be
     *  automatically closed.  The file will also be closed when the
     *  IBSimu object is destructed.
     */
    std::ostream &set_message_output( const std::string &filename );

    /*! \brief Print message output.
     *
     *  Returns a reference to stream for output \a type and
     *  importance \a level. If level is larger than threshold the
     *  function returns a reference to nullstream and nothing will be
     *  printed.
     */
    std::ostream &message( message_type_e type, int32_t level );

    std::ostream &message( int32_t level );

    void flush( bool finishlines = true );

    /*! \brief Increase message indentation.
     */
    void inc_indent( void );

    /*! \brief Decrease message indentation.
     */
    void dec_indent( void );

    /*! \brief Return if message output file is stdout.
     */
    bool output_is_cout();

    /*! \brief Set message threshold level.
     *
     *  Only messages with level lower than or equal to the threshold
     *  will be printed. A value of 1 is used for standard amount of
     *  output. Value of 2 for extended amount of output. Message
     *  threshold for MSG_VERBOSE defaults to 0 (no output), MSG_ERROR
     *  and MSG_WARNING default to 1. For enabling the debug messages,
     *  the library has to be compiled with debugging enabled.
     */
    void set_message_threshold( message_type_e type, int32_t level );

    /*! \brief Get message threshold level.
     *
     *  Values less than or equal to zero mean no output will be printed.
     */
    int32_t get_message_threshold( message_type_e type );




    /*! \brief Set the number of threads used for calculation.
     */
    void set_thread_count( uint32_t threadcount );

    /*! \brief Get the number of threads used for calculation.
     */
    uint32_t get_thread_count( void ) { return( _threadcount ); }



    /*! \brief Set the random number generator to use.
     *
     *  The default rng type RNG_SOBOL is a quasi random number
     *  generator based on use of Sobol sequence.  This setting
     *  affects particle distribution definitions.
     */
    void set_rng_type( rng_type_e type );

    /*! \brief Get the random number generator type.
     */
    rng_type_e get_rng_type( void ) { return( _rng ); }
    


    /*! \brief Halt execution
     *
     *  This function is called by the error handler in case of SIGTERM.
     */
    void halt( void );
};


/*! \brief Global instance of class %IBSimu.
 */
extern IBSimu ibsimu;


#endif

