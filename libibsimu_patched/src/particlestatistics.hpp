/*! \file particlestatistics.hpp
 *  \brief %Particle statistics
 */

/* Copyright (c) 2010-2011 Taneli Kalvas. All rights reserved.
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

#ifndef PARTICLESTATISTICS_HPP
#define PARTICLESTATISTICS_HPP 1


#include <stdint.h>
#include <ostream>
#include <vector>


/*! \brief %Particle iteration statistics.
 *
 *  Stores statistics about the particle histories.
 *
 */
class ParticleStatistics {

    uint32_t      _end_time;          /*!< \brief Number of time limited particle iterations. */
    uint32_t      _end_step;          /*!< \brief Number of step count limited particle iterations. */
    uint32_t      _end_baddef;        /*!< \brief Number of bad particle definitions. */
    uint32_t      _sum_steps;         /*!< \brief Total number of steps taken. */

    std::vector<uint32_t> _bound_collisions;  /*!< \brief Number of particles collided with electrodes. */
    std::vector<double>   _bound_current;     /*!< \brief Amount of current collided with electrodes. */

public:

    ParticleStatistics();
    ParticleStatistics( const ParticleStatistics &stat );
    ParticleStatistics( uint32_t nboundaries );

    /*! \brief Constructor for loading particle statistics from a file.
     */
    ParticleStatistics( std::istream &s );

    ~ParticleStatistics();

    const ParticleStatistics &operator=( const ParticleStatistics &stat );
    const ParticleStatistics &operator+=( const ParticleStatistics &stat );

    void clear( void );
    void reset( uint32_t nboundaries );
    
    uint32_t end_time( void ) const;
    uint32_t end_step( void ) const;
    uint32_t end_baddef( void ) const;
    uint32_t sum_steps( void ) const;

    uint32_t number_of_boundaries( void ) const;
    uint32_t bound_collisions( uint32_t bound ) const;
    uint32_t bound_collisions( void ) const;
    double bound_current( uint32_t bound ) const;
    double bound_current( void ) const;

    void inc_end_time( void ) { _end_time++; }
    void inc_end_step( void ) { _end_step++; }
    void inc_end_baddef( void ) { _end_baddef++; }
    void inc_sum_steps( void ) { _sum_steps++; }
    void inc_sum_steps( uint32_t i ) { _sum_steps += i; }

    void add_bound_collision( uint32_t bound, double IQ );

    /*! \brief Saves data to stream.
     */
    void save( std::ostream &s ) const;
};


#endif

