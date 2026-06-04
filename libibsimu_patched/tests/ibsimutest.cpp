/*! \file ibsimutest.cpp
 *  \brief Common functions for tests.
 */


#include <iostream>
#include <string>
#include "ibsimutest.hpp"
#include "ibsimu.hpp"


ErrorTest::ErrorTest( const ErrorLocation &loc, const std::string &str )
    : Error(loc)
{
    _error_str = "Test failure: " + str;
}


ErrorTest::ErrorTest( const ErrorLocation &loc )
    : Error(loc)
{
    _error_str = "Test failure\n";
}


/* Main function for tests
 */
int main( int argc, char **argv )
{
    try {
	// Libtool binary name is "lt-blah", strip the main part
	std::string filename;
	std::string bin = argv[0];
	size_t loc = bin.rfind( "lt-" );
	if( loc != std::string::npos )
	    filename = bin.substr( loc+3 ) + "_vout.txt";
	else
	    filename = bin + "_vout.txt";

	ibsimu.set_message_output( filename );
	//ibsimu.set_message_threshold( MSG_DEBUG_GENERAL, 1 );
	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	test( argc, argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
	exit( 1 );
    }

    return( 0 );
}
