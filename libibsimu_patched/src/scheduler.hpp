/*! \file scheduler.hpp
 *  \brief Job scheduler for parallel processing
 */

/* Copyright (c) 2005-2009,2011 Taneli Kalvas, Jan Sarén. All rights reserved.
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

#ifndef SCHEDULER_HPP
#define SCHEDULER_HPP 1


#include <pthread.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <deque>
#include <time.h>
#include "comptime.hpp"


//#define SCHEDULER_DEBUG 1


//pthread_mutex_t cout_mutex = PTHREAD_MUTEX_INITIALIZER;


/*! \brief %Scheduler class for implementing consumer-producer
 *  threading.
 *
 *  %Scheduler uses a manager thread and a given number of working
 *  threads for solving problems. %Scheduler is templated with Solver,
 *  Problem and Error classes. Solver class has to provide an operator
 *  \code void operator()( Problem *p, int32_t pi ) \endcode 
 *  to solve problem p with index location pi.
 *
 *  Error class has to have a default constructor. Scheduler catches
 *  the errors of this type from the working threads and saves caught
 *  errors in a container. If an error is caught, all the working
 *  threads are interrupted and problem solving is finished. The
 *  scheduler does indicate the error state by returning false from
 *  finish(). Error state can also be queried with is_error(). Errors
 *  can be fetched from the internal containers with get_errors().
 *
 *  The %Scheduler can be used for static and dynamic problems. All
 *  threads can add problems. The problem vector is a shared resource
 *  so it must be protected with lock_mutex() and unlock_mutex().
 *
 *  Parallel processing is started with run() function and the end of
 *  processing can be waited with finish().
 */
template <class Solv, class Prob, class Err>
class Scheduler {

    class Consumer {

	/*
	enum consumer_status_e {
	    CONSUMER_CREATED = 0,
	    CONSUMER_RUNNING,
	    CONSUMER_FINISHED
	};
	*/

	//pthread_mutex_t      _mutex;             //!< \brief Mutex for active check
	pthread_t            _thread;            //!< \brief POSIX thread handle
	Solv                *_solver;            //!< \brief Solver of this consumer
	Scheduler           *_scheduler;         //!< \brief Pointer to scheduler
	//struct timeval       _t0;
	//std::vector<struct timeval> _t;
    
	void *consumer_main( void ) {
	    //struct timeval t;
	    
#ifdef SCHEDULER_DEBUG
	    std::cout << "Consumer main entrance\n";
#endif
	    //pthread_mutex_lock( &_mutex );
	    //_status = CONSUMER_RUNNING;
	    //pthread_mutex_unlock( &_mutex );

	    Prob *p;
	    uint32_t pi;
	    while( (p = _scheduler->get_next_problem( pi )) ) {
		try {
		    //gettimeofday( &t, NULL );
		    //_t.push_back( t );
		    (*_solver)( p, pi );
		    //gettimeofday( &t, NULL );
		    //_t.push_back( t );
		} catch( Err e ) {
		    //std::cout << "on_error\n";
		    // Handle error and stop solving
		    _scheduler->on_error( e, pi );
		    break;
		};
		_scheduler->inc_solved_problem();
	    }

#ifdef SCHEDULER_DEBUG
	    std::cout << "Exiting consumer\n";
#endif
	    //pthread_mutex_lock( &_mutex );
	    //_status = CONSUMER_FINISHED;
	    //pthread_mutex_unlock( &_mutex );
	    return( NULL );
	}
    
    public:

	static void *consumer_entry( void *data ) {
	    Consumer *consumer = (Consumer *)data;
	    return( consumer->consumer_main() );
	}

	Consumer( Solv *solver, Scheduler *scheduler ) : _solver(solver), _scheduler(scheduler) { 

	    //pthread_mutex_init( &_mutex, NULL );
#ifdef SCHEDULER_DEBUG
	    std::cout << "Consumer constructor\n";
#endif
	    //gettimeofday( &_t0, NULL );
	}

	~Consumer() {
	    //pthread_mutex_lock( &cout_mutex );
#ifdef SCHEDULER_DEBUG
	    std::cout << "Consumer destructor\n";
#endif
	    //for( size_t a = 0; a < _t.size(); a++ ) {
	    //std::cout << (_t[a].tv_sec-_t0.tv_sec) + 
	    //(_t[a].tv_usec-_t0.tv_usec)/1e6 << "\n";
	    //a++;
	    //std::cout << (_t[a].tv_sec-_t0.tv_sec) + 
	    //(_t[a].tv_usec-_t0.tv_usec)/1e6 << "\n\n\n";
	    //}
	    //pthread_mutex_unlock( &cout_mutex );
	}

	void run( void ) {
	    pthread_create( &_thread, NULL, consumer_entry, (void *)this );
	}

	void join( void ) {

#ifdef SCHEDULER_DEBUG
	    std::cout << "Consumer join\n";
#endif
	    //pthread_mutex_lock( &_mutex );
	    //if( _status == CONSUMER_FINISHED ) {
	    //pthread_mutex_unlock( &_mutex );
	    //return;
	    //} else if( _status == CONSUMER_CREATED ) {
	    //
	    //}
	    //pthread_mutex_unlock( &_mutex );
	    pthread_join( _thread, NULL );
	}

    };


    pthread_mutex_t         _mutex;            //!< \brief Mutex for all shared data
    pthread_cond_t          _scheduler_cond;   //!< \brief Wake up scheduler
    pthread_cond_t          _producer_cond;    //!< \brief Wake up producer
    pthread_cond_t          _consumer_cond;    //!< \brief Wake up consumer

    //size_t                  _problems_in_c;    //!< \brief Total problems in count
    //size_t                  _problems_err_c;   //!< \brief Total error problems out count
    //std::deque<Prob*>       _problems_out;     //!< \brief Problems already solved

    uint32_t                _read_c;           //!< \brief Read problems count
    uint32_t                _solved_c;         //!< \brief Solved problems count
    std::vector<Prob *>    &_problems;         //!< \brief Vector of problems

    pthread_t               _scheduler_thread; //!< \brief Scheduler main thread
    std::vector<Consumer *> _consumers;        //!< \brief Consumer objects

    bool                    _join;             //!< \brief Is join needed?
    bool                    _running;          //!< \brief Are we running
    bool                    _error;            //!< \brief Finish as soon as possible
    bool                    _done;             //!< \brief Exit after current problem
    bool                    _finish;           //!< \brief Finish all problems and exit
    std::vector<Err>        _err;              //!< \brief Errors encountered
    std::vector<int32_t>    _eprob;            //!< \brief Indices for problems causing errors


    /*! \brief %Error handler
     *
     *  Saves caught error \a e and index \a pi of problem causing the
     *  error and broadcasts other threads signalling about the
     *  error condition.
     */
    void on_error( Err &e, uint32_t pi ) {
	pthread_mutex_lock( &_mutex );
	_err.push_back( e );
	_eprob.push_back( pi );
	_error = true;
	pthread_cond_broadcast( &_scheduler_cond );
	pthread_mutex_unlock( &_mutex );
    }


    /*! \brief Increase the counter of solved problems
     */
    void inc_solved_problem( void ) {
	pthread_mutex_lock( &_mutex );
	_solved_c++;
	pthread_mutex_unlock( &_mutex );
    }

    /*! \brief Get pointer to next problem.
     *
     *  Returns pointer to the next problem from vector of
     *  problems. Increases read problems counter. Returns NULL if
     *  no problem returned. Also returns problem index in \a pi
     */
    Prob *get_next_problem( uint32_t &pi ) {
#ifdef SCHEDULER_DEBUG
	    std::cout << "get_next_problem()\n";
#endif
	pthread_mutex_lock( &_mutex );
    
	if( _done || _error ) {
	    pthread_mutex_unlock( &_mutex );
#ifdef SCHEDULER_DEBUG
	    std::cout << "get_next_problem(): Returning NULL\n";
#endif
	    pi = -1;
	    return( NULL );
	}

	if( _problems.size() == _read_c ) {
#ifdef SCHEDULER_DEBUG
	    std::cout << "get_next_problem(): No problem to return... waiting\n";
#endif
	    // Signal producer that problems are spent
	    pthread_cond_signal( &_scheduler_cond );
	    while( _problems.size() == _read_c ) {
		// Wait for new problems
		pthread_cond_wait( &_consumer_cond, &_mutex );
		if( _done || _error ) {
		    pthread_mutex_unlock( &_mutex );
#ifdef SCHEDULER_DEBUG
		    std::cout << "get_next_problem(): Returning NULL\n";
#endif
		    pi = -1;
		    return( NULL );
		}
	    }
	}

	// Return next problem
	pi = _read_c++;
	Prob *ret = _problems[pi];

#ifdef SCHEDULER_DEBUG
	std::cout << "get_next_problem(): Returning problem " << pi << "\n";
#endif

	pthread_mutex_unlock( &_mutex );
	return( ret );
    }

    
    /*! \brief Scheduler main.
     */
    void *scheduler_main( void ) {

#ifdef SCHEDULER_DEBUG
	std::cout << "Running scheduler_main()\n";
#endif

	// Start consumer threads
	for( size_t a = 0; a < _consumers.size(); a++ )
	    _consumers[a]->run();

	pthread_mutex_lock( &_mutex );

	while( 1 ) {
	    // Wait until all consumers are done with all problems or error occurs
	    while( !(_problems.size() == _solved_c || _done || _error) ) {
		//std::cout << "scheduler_main(): scheduler_cond wait 1\n";
		pthread_cond_wait( &_scheduler_cond, &_mutex );
	    }

	    if( (_finish && _problems.size() == _solved_c) || _done || _error )
		break;

	    // Problems temporarily done
	    pthread_cond_wait( &_scheduler_cond, &_mutex );
	    //std::cout << "scheduler_main(): prob_in = " << _problems_in_c
	    //<< " prob_out = " << _problems_out_c << "\n";
	    //std::cout << "scheduler_main(): scheduler_cond wait 2\n";

	    // Signal consumers to wake up
	    pthread_cond_broadcast( &_consumer_cond );
	}

	// Broadcast: done
	_done = true;
	_running = false;
	pthread_cond_broadcast( &_consumer_cond );
	pthread_mutex_unlock( &_mutex );

	// Join all consumers
	//std::cout << "scheduler_main(): Scheduler waiting in join\n";
	for( size_t a = 0; a < _consumers.size(); a++ )
	    _consumers[a]->join();

	pthread_cond_broadcast( &_producer_cond );
	//std::cout << "scheduler_main(): Exiting scheduler\n";
	return( NULL );
    }


    /*! \brief Scheduler main loop entry point.
     */
    static void *scheduler_entry( void *data ) {
	Scheduler *scheduler = (Scheduler *)data;
	return( scheduler->scheduler_main() );
    }

    /*! \brief Private copy constructor
     */
    Scheduler( const Scheduler &s ) {}

    /*! \brief Private copy operator
     */
    const Scheduler &operator=( const Scheduler &s ) {
	return( *this );
    }

public:


    /*! \brief Constructor
     *
     *  Constructor for scheduler solving problems in vector \a prob.
     */
    Scheduler( std::vector<Prob *> &prob ) 
	: _read_c(0), _solved_c(0), _problems(prob), _join(false), _running(false), 
	  _error(false), _done(false), _finish(false) {

	// Initialize pthread objects
	pthread_mutex_init( &_mutex, NULL );
	pthread_cond_init( &_scheduler_cond, NULL );
	pthread_cond_init( &_consumer_cond, NULL );
	pthread_cond_init( &_producer_cond, NULL );
    }


    /*! \brief Destructor
     */
    ~Scheduler() {

	// Force exit
	_done = true;
	finish();

	pthread_mutex_destroy( &_mutex );
	pthread_cond_destroy( &_scheduler_cond );
	pthread_cond_destroy( &_consumer_cond );
	pthread_cond_destroy( &_producer_cond );
    }


    /*! \brief Return true on errors.
     */
    bool is_error( void ) {
	// No mutex needed for one bit read
	return( _error );
    }


    /*! \brief Return true if scheduler is running.
     */
    bool is_running( void ) {
	// No mutex needed for one bit read
	return( _running );
    }


    /*! \brief Return number of solved problems.
     */
    uint32_t get_solved_count( void ) {
	pthread_mutex_lock( &_mutex );
	uint32_t ret = _solved_c;
	pthread_mutex_unlock( &_mutex );
	return( ret );
    }


    /*! \brief Return number of problems.
     */
    uint32_t get_problem_count( void ) {
	pthread_mutex_lock( &_mutex );
	uint32_t ret = _problems.size();
	pthread_mutex_unlock( &_mutex );
	return( ret );
    }


    /*! \brief Fetch errors and indices of corresponding problems.
     *
     *  \param e Container where errors are appended
     *  \param pi Container where problems are appended
     *  \return Number of error problems added to containers
     *
     *  Errors are removed from %Scheduler.
     */
    template <class Cont1, class Cont2>
    size_t get_errors( Cont1 &e, Cont2 &pi ) {
	pthread_mutex_lock( &_mutex );
	size_t r = _err.size();
	for( size_t a = 0; a < _err.size(); a++ ) {
	    e.push_back( _err[a] );
	    pi.push_back( _eprob[a] );
	}
	_err.clear();
	_eprob.clear();
	pthread_mutex_unlock( &_mutex );
	return( r );
    }    


    /*! \brief Run threads with \a N solvers.
     *
     *  Returns immediately after creating working threads. Use
     *  finish() or destructor of class to wait for work to be
     *  completed.
     */
    void run( std::vector<Solv *> solv ) {

	// Do nothing if already running
	if( _running )
	    return;

	// Create consumer threads
	for( size_t a = 0; a < solv.size(); a++ )
	    _consumers.push_back( new Consumer( solv[a], this ) );

	_read_c = 0;
	_solved_c = 0;
	_join = true;
	_running = true;
	_error = false;
	_done = false;
	_finish = false;
	_err.clear();
	_eprob.clear();

	// Create scheduler thread
	pthread_create( &_scheduler_thread, NULL, scheduler_entry, (void *)this );
    }


    /*! \brief Lock mutex for adding problems.
     */
    void lock_mutex( void ) {

	pthread_mutex_lock( &_mutex );
    }

    
    /*! \brief Unlock mutex.
     */
    void unlock_mutex( void ) {

	pthread_cond_broadcast( &_scheduler_cond );	
	pthread_mutex_unlock( &_mutex );
    }
    

    /*! \brief Force exit from scheduler.
     *
     *  Waits for the solvers to finish the problems under progress
     *  after, which all threads are stopped and resources
     *  freed. Returns false if any error occured during last
     *  run. Otherwise returns true.
     */
    bool force_exit( void ) {

	_done = true;
	return( finish() );
    }

    /*! \brief Call for solvers to finish when problems are solved.
     *
     *  Waits for solvers to finish for 1 sec. Returns true if solver
     *  finished, false if not. The finish() function needs to be
     *  called after true to free the resources.
     */
    bool wait_finish( void ) {

	pthread_mutex_lock( &_mutex );
	if( _running ) {
	    _finish = true;
	    pthread_cond_broadcast( &_scheduler_cond );

	    struct timespec ts;
	    ibs_clock_gettime( CLOCK_REALTIME, &ts );
	    ts.tv_sec += 1;
	    int rc = pthread_cond_timedwait( &_producer_cond, &_mutex, &ts );
	    if( rc == ETIMEDOUT ) {
		pthread_mutex_unlock( &_mutex );
		return( false );
	    }
	}
	pthread_mutex_unlock( &_mutex );
	return( true );
    }

    /*! \brief Wait for all problems to be solved.
     *
     *  Waits for all solvers to finish. Stops all threads and frees
     *  resources. Returns false if any error occured during last
     *  run. Otherwise returns true.
     *
     *  Solver is prematurely exited if an error occurs.
     */
    bool finish( void ) {

	pthread_mutex_lock( &_mutex );
	if( _running ) {
	    _finish = true;
	    //std::cout << "finish(): scheduler_cond broadcast\n";
	    pthread_cond_broadcast( &_scheduler_cond );
	
	    //std::cout << "finish(): producer_cond wait\n";
	    pthread_cond_wait( &_producer_cond, &_mutex );
	}
	pthread_mutex_unlock( &_mutex );

	if( _join ) {
	    // Delete consumer threads
	    for( size_t a = 0; a < _consumers.size(); a++ )
		delete _consumers[a];
	    _consumers.clear();

	    pthread_join( _scheduler_thread, NULL );
	    _join = false;
	}
	
	if( _error )
	    return( false );

	return( true );
    }

    friend class Consumer;
};



#endif

