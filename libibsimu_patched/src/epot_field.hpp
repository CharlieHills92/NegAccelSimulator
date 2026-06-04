/*! \file epot_field.hpp
 *  \brief Electric potential field.
 */

/* Copyright (c) 2011 Taneli Kalvas. All rights reserved.
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

#ifndef EPOT_FIELD_HPP
#define EPOT_FIELD_HPP 1


#include "geometry.hpp"
#include "meshscalarfield.hpp"


/*! \brief Electric potential field.
 *
 *  Electric potential field based on MeshScalarField. Contains an evaluator
 */
class EpotField : public MeshScalarField {

    const Geometry *_geom;      /*!< \brief Pointer to geometry. */

    // Forbidden
    EpotField( std::istream &s ) {}

public:

    /*! \brief Constructor.
     *
     *  The field is created using the mesh geometry and solid
     *  definitions of \a geom.
     */
    EpotField( const Geometry &geom );

    /*! \brief Copy constructor.
     */
    EpotField( const EpotField &f );

    /*! \brief Constructor for loading EpotField field data from a
     *  file and using solid definitions from \a geom.
     */
    EpotField( std::istream &s, const Geometry &geom );

    /*! \brief Destructor.
     */
    virtual ~EpotField();

    /*! \brief Get a pointer to geometry.
     */
    const Geometry *geom( void ) const;

    using MeshScalarField::operator();

    /*! \brief Operator for getting linearly interpolated field values.
     *
     *  The field is interpolated linearly to get the field value at
     *  \a x. If \a x is outside the mesh, the field is extrapolated
     *  linearly using the field points to \a x. This provides correct
     *  field values also close to the mesh boundaries.
     *
     *  The near solid points are evaluated taking in account the
     *  distance to the solid surface.
     */
    virtual double operator()( const Vec3D &x ) const;
};


#endif
