/*! \file scharge.hpp
 *  \brief Space charge deposition functions
 */

/* Copyright (c) 2005-2013 Taneli Kalvas. All rights reserved.
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

#ifndef SCHARGE_HPP
#define SCHARGE_HPP 1


#include <pthread.h>
#include "scalarfield.hpp"
#include "particles.hpp"
#include "cfifo.hpp"


/*! \brief Finalize space charge calculation.
 *
 *  Converts charge map built by calls to
 *  scharge_add_from_trajectory() to space charge density map suitable
 *  for potential solvers. Also corrects boundary space charge values.
 *
 */
void scharge_finalize_pic( MeshScalarField &scharge );


/*! \brief Function for adding charge to space charge density map from
 *  particle trajectory in 2d simulation.
 *
 *  Adds contribution to space charge density map \a scharge from a
 *  particle carrying current per length \a I travelling from \a x1 to
 *  \a x2. Charge density is added to the mesh assuming that it is
 *  localized at the point where particle was at time
 *  0.5*(t1+t2). This function is supposed to be used so that \a x1
 *  and \a x2 are at the mesh intersection points, through which
 *  particle trajectory has passed.
 */
void scharge_add_from_trajectory_pic( MeshScalarField &scharge, pthread_mutex_t *mutex, 
				      double I, const ParticleP2D &x1, const ParticleP2D &x2 );


/*! \brief Function for adding charge to space charge density map from
 *  particle trajectory in cylindrically symmetric simulation.
 *
 *  Same as for scharge_add_from_trajectory(). Now I is real current
 *  (A).
 */
void scharge_add_from_trajectory_pic( MeshScalarField &scharge, pthread_mutex_t *mutex, 
				      double I, const ParticlePCyl &x1, const ParticlePCyl &x2 );


/*! \brief Function for adding charge to space charge density map from
 *  particle trajectory in 3d simulation.
 *
 *  Same as for scharge_add_from_trajectory(). Now I is real current
 *  (A).
 */
void scharge_add_from_trajectory_pic( MeshScalarField &scharge, pthread_mutex_t *mutex, 
				      double I, const ParticleP3D &x1, const ParticleP3D &x2 );



void scharge_finalize_linear( MeshScalarField &scharge );

void scharge_add_from_trajectory_linear( MeshScalarField &scharge, pthread_mutex_t *mutex, 
					 double I, int dir, const CFiFo<ParticleP2D,4> &cdpast, const int i[3] );

void scharge_add_from_trajectory_linear( MeshScalarField &scharge, pthread_mutex_t *mutex, 
					 double I, int dir, const CFiFo<ParticleP3D,4> &cdpast, const int i[3] );

void scharge_add_from_trajectory_linear( MeshScalarField &scharge, pthread_mutex_t *mutex, 
					 double I, int dir, const CFiFo<ParticlePCyl,4> &cdpast, const int i[3] );










void scharge_finalize_step_pic( MeshScalarField &scharge );
void scharge_add_step_pic( MeshScalarField &scharge, double Q, const Vec3D &x );


#endif
