/*! \file epot_efield.hpp
 *  \brief Electric potential base electric field.
 */

/* Copyright (c) 2005-2014 Taneli Kalvas. All rights reserved.
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

#ifndef EPOT_EFIELD_HPP
#define EPOT_EFIELD_HPP 1


#include "vectorfield.hpp"
#include "epot_field.hpp"
#include "vec3d.hpp"
#include "types.hpp"


/*! \brief %Vector field based on interpolation of electric
 *  potential.
 *
 *  %EpotEfield contains pointers to Geometry and to MeshScalarField
 *  electric potential (epot). The electric field is defined in points
 *  between the nodes using equation Ex=(Ex0-Ex1)/h. The Ex, Ey and Ez
 *  are not defined in same positions as the nodes of the electric
 *  potential because of this. The field evaluator uses linear
 *  interpolation to return smoothly varying field values. Close to
 *  the solids, special consideration is used to return as good
 *  estimate of the field as possible.
 *
 *  The function recalculate() has to be called to remake the electric
 *  field.
 *
 *  The behaviour of the interpolation function outside mesh points
 *  can be programmed with set_extrapolation() function. Behaviour
 *  defaults to extrapolation using closest electric potential points. 
 */
class EpotEfield : public VectorField {

    field_extrpl_e    _extrpl[6];   /*!< \brief What to return outside geometry. */
    const EpotField  &_epot;        /*!< \brief Reference to electric potential. */
    const Geometry   *_geom;        /*!< \brief Pointer to geometry. */

    double           *_F[3];        /*!< \brief Vector field data in three components */

    uint8_t solid_dist( uint32_t node, uint32_t dir ) const;

    void copy_1d( const EpotEfield &efield );
    void copy_2d( const EpotEfield &efield );
    void copy_3d( const EpotEfield &efield );

    void precalc_1d( void );
    void precalc_2d( void );
    void precalc_3d( void );
    void precalc( void );

public:

    /*! \brief Constructor.
     */
    EpotEfield( const EpotField &epot );

    /*! \brief Copy constructor.
     */
    EpotEfield( const EpotEfield &efield );

    /*! \brief Destructor.
     */
    ~EpotEfield();

    /*! \brief Set the behaviour of electric field interpolation
     *  outside mesh points (extrapolation).
     *
     *  The interpolation function behaviour can be set separately for
     *  each boundary. This is done by setting the desired properties
     *  to the \a extrpl array. The interpolation function can use an
     *  extrapolation of the last three electric potential values for
     *  calculation of electric field (FIELD_EXTRAPOLATE) or it can
     *  return the mirror of the electric field across the mesh
     *  boundary like E_x(x)=E_x(-x) (FIELD_MIRROR), or it can return
     *  the mirror of the electric field across the mesh boundary like
     *  E_x(x,y,z)=-E_x(-x,y,z) (FIELD_ANTIMIRROR), or it can return a
     *  zero electric field outside the mesh (FIELD_ZERO) or NaN,
     *  which prevents the use of field evaluations in particle
     *  iterator (FIELD_NAN).
     *
     *  For simulation of symmetric cases, where electric potential is
     *  symmetric across the boundary like phi(x)=phi(-x) the E-field
     *  can be evaluated like with FIELD_ANTIMIRROR setting, but
     *  additionally the E-field is forced to be zero at the
     *  boundary. This extrapolation behaviour is selected with
     *  FIELD_SYMMETRIC_POTENTIAL. This gives second order correct
     *  electric field next to the boundary, which is essential to get
     *  physical results in cases where beam propagates next the the
     *  boundary.
     *
     *  Very far (double the size of the simulation box) the field
     *  evaluator will always return NaN.
     */
    void set_extrapolation( field_extrpl_e extrpl[6] );

    /*! \brief Recalculate electric field from potential.
     */
    void recalculate( void );

    /*! \brief Operator for getting interpolated electric field value
     *  at \a x.
     */
    virtual const Vec3D operator()( const Vec3D &x ) const;

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;
};


#endif

