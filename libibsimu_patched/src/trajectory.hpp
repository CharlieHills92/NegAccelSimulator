/*! \file trajectory.hpp
 *  \brief Trajectory interpolation solver
 */

/* Copyright (c) 2005-2009,2011,2013 Taneli Kalvas. All rights reserved.
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

#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP 1


#include <iostream>
#include <stdint.h>


enum trajectory_rep_e {
    TRAJ_EMPTY = 0,
    TRAJ_LINEAR,
    TRAJ_QUADRATIC,
    TRAJ_CUBIC
};


/*! \brief Trajectory representation between two calculated points in 1d.
 *
 *  Uses either linear-linear, quadratic-linear or cubic-quadratic
 *  representation for location \a x and velocity \a v. Time is
 *  presented as a parametric, scaled variable \a K ranging from 0 to 1.
 */
class TrajectoryRep1D {

    trajectory_rep_e _rep;
    double           _A, _B, _C, _D, _E;

    static bool in( double K, int extrapolate );

public:

    /*! \brief Default constructor for empty representation.
     */
    TrajectoryRep1D() : _rep(TRAJ_EMPTY) {}

    /*! \brief Constructor for representation of trajectory from \a
     *  (x1,v1) to \a (x2,v2) in time \a dt.
     *
     *  Can be forced to use a specified representation type by
     *  setting \a force. Defaults to TRAJ_EMPTY, which means that the
     *  highest numerically stable polynomial is automatically used.
     */
    TrajectoryRep1D( double dt, double x1, double v1, double x2, double v2,
		     trajectory_rep_e force = TRAJ_EMPTY );

    /*! \brief Destructor.
     */
    ~TrajectoryRep1D() {}

    /*! \brief Construct representation of trajectory from \a (x1,v1)
     *  to \a (x2,v2) in time \a dt.
     *
     *  Can be forced to use a specified representation type by
     *  setting \a force. Defaults to TRAJ_EMPTY, which means that the
     *  highest numerically stable polynomial is automatically used.
     */
    void construct( double dt, double x1, double v1, double x2, double v2,
		    trajectory_rep_e force = TRAJ_EMPTY );

    /*! \brief Calculate location \a x and velocity \a v at parametric
     *  time \a K.
     */
    void coord( double &x, double &v, double K );

    /*! \brief Solves for trajectory intersection with location.
     *
     *  Solves the trajectory intersection with location \a x. Saves
     *  the valid solutions to array \a K in increasing order and
     *  returns the number of solutions saved. The accepted limit for
     *  parametric time \a K is 0 < K =< 1 if extrapolate is 0,
     *  -1.0e-6 < K =< 1 if extrapolate is less than 0 or 0.0 < K <=
     *  1+1.0e-6 if extrapolate is more than 0.
     */
    int solve( double K[3], double x, int extrapolate = 0 );

    /*! \brief Returns order of polynomial used.
     *
     *  Returns order 0-3.
     */
    uint32_t get_representation_order( void ) const;

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;
};


#endif



















