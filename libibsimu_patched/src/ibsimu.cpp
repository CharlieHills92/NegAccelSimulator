/*! \file ibsimu.cpp
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


#include "config.h"
#include "id.hpp"
#include <iostream>
#include <signal.h>
#include "ibsimu.hpp"
#include "timer.hpp"
#include "error.hpp"



IBSimu::IBSimu( const IBSimu &ibs )
    : _t(0), _rng(ibs._rng), _os(&_ns)
{

}


IBSimu::IBSimu()
    : _hello(false), _threadcount(1), _rng(RNG_SOBOL), _os(&std::cout), _indent(0)
{
    // Set message level thresholds to defaults
    for( int a = 0; a < MSG_COUNT; a++ )
	_message_threshold[a] = 0;
    _message_threshold[MSG_WARNING] = 1;
    _message_threshold[MSG_ERROR] = 1;

#if defined(_GNU_SOURCE) && defined(HAVE_SIGINFO_T)
    // Set a catch for segmentation fault
    struct sigaction act_sigsegv;
    act_sigsegv.sa_sigaction = SignalHandler::signal_handler_SIGSEGV;
    sigemptyset( &act_sigsegv.sa_mask );
    act_sigsegv.sa_flags = SA_SIGINFO;
    sigaction( SIGSEGV, &act_sigsegv, NULL );

    // Set a catch for terminate/kill/int
    struct sigaction act_sigterm;
    act_sigterm.sa_sigaction = SignalHandler::signal_handler_SIGTERM;
    sigemptyset( &act_sigterm.sa_mask );
    act_sigterm.sa_flags = SA_SIGINFO;
    sigaction( SIGTERM, &act_sigterm, NULL );
    sigaction( SIGQUIT, &act_sigterm, NULL );
    sigaction( SIGINT, &act_sigterm, NULL );
#else
    signal( SIGTERM, SignalHandler::signal_handler_SIGTERM );
    signal( SIGINT, SignalHandler::signal_handler_SIGTERM );
#endif

    // Start timer for whole simulation
    _t = new Timer;
}


IBSimu::~IBSimu()
{
    _t->stop();
    message( 1 ) << "Ending simulation\n";
    message( 1 ) << "  time used = " << *_t << "\n";
    flush();

    // Close the message output file if it is open
    if( _fo.is_open() )
	_fo.close();

    delete _t;
}


std::ostream &IBSimu::message( message_type_e type, int32_t level )
{
    if( type >= MSG_COUNT || type < 0 )
	throw( Error( ERROR_LOCATION, "invalid output type" ) );

    if( level > _message_threshold[type] )
	return( _ns );

    //return( *_os );
    return( _ss );
}


std::ostream &IBSimu::message( int32_t level )
{
    if( level > _message_threshold[MSG_VERBOSE] )
	return( _ns );

    flush();

    //return( *_os );
    return( _ss );
}


void IBSimu::inc_indent( void )
{
    flush( false );
    _indent++;
}


void IBSimu::dec_indent( void )
{
    flush( false );
    _indent--;
}


std::ostream &IBSimu::set_message_output( std::ostream &os )
{
    if( _fo.is_open() ) {
	_fo.close();
	_os = &os;
	return( _ns );
    } else {
	std::ostream *osold = _os;
	_os = &os;
	return( *osold );
    }
}


std::ostream &IBSimu::set_message_output( const std::string &filename )
{
    if( _fo.is_open() ) {
	_fo.close();
	_fo.open( filename.c_str() );
	if( !_fo.is_open() )
	    throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\' for writing" ) );
	_os = &_fo;
	return( _ns );
    } else {
	std::ostream *osold = _os;
	_fo.open( filename.c_str() );
	if( !_fo.is_open() )
	    throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\' for writing" ) );
	_os = &_fo;
	return( *osold );
    }
}


bool IBSimu::output_is_cout()
{
    return( _os == &std::cout ); 
}


void IBSimu::set_message_threshold( message_type_e type, int32_t level )
{
    if( type >= MSG_COUNT || type < 0 )
	throw( Error( ERROR_LOCATION, "invalid output message type" ) );

    // Set message threshold level
    _message_threshold[type] = level;

    // If setting to verbose mode and no greeting has yet been shown
    if( type == MSG_VERBOSE && level > 0 && !_hello ) {
	_hello = true;
	message( 1 ) << "Ion Beam Simulator " VERSION " (" IBSIMU_GIT_ID ")\n";
    }
}


int32_t IBSimu::get_message_threshold( message_type_e type )
{
    if( type >= MSG_COUNT || type < 0 )
	throw( Error( ERROR_LOCATION, "invalid output type" ) );

    return( _message_threshold[type] );
}


void IBSimu::flush( bool finishlines )
{
    std::string stmp;
    bool        lineunfinished;
    convert_stringstream_to_lines( _ss,
				   _indent*2,
				   stmp,
				   lineunfinished );
    *_os << stmp;
    *_os << std::flush;

    if( lineunfinished && finishlines )
	*_os << "\n";

    _ss.str( "" );   
}


size_t IBSimu::convert_stringstream_to_lines( const std::stringstream &ss,
					      size_t indentation,
					      std::string &target,
					      bool &lineunfinished )
{
    std::string stmp( ss.str() );

    lineunfinished = false;

    if( stmp.size() ) {
	for( size_t i = 0; i < indentation; i++ )
	    target += ' ';
	lineunfinished = true;
    }
    
    size_t pos = 0;
    size_t linecount = 0;
    while( pos < stmp.size() ) {
	if( stmp[pos] == '\n' ) {
	    lineunfinished = false;
	    ++linecount;
	    
	    target += '\n'; 
	} else {
	    if( lineunfinished == false )
		for( size_t i = 0; i < indentation; i++ )
		    target += ' ';
	    
	    target += stmp[pos]; 
	    lineunfinished = true;
	}
	
	pos++;
    }

    return( linecount );
}


void IBSimu::set_thread_count( uint32_t threadcount ) 
{
    if( threadcount <= 0 )
	throw( Error( ERROR_LOCATION, "invalid parameter" ) );
    _threadcount = threadcount;
}


void IBSimu::set_rng_type( rng_type_e type ) 
{
    if( type != RNG_SOBOL && type != RNG_MT )
	throw( Error( ERROR_LOCATION, "unknown rng type" ) );
    _rng = type;
}


/* Global instance */
IBSimu ibsimu;


