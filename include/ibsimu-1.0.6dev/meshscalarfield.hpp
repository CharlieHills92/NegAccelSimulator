/*! \file meshscalarfield.hpp
 *  \brief Mesh based scalar fields
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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

#ifndef MESHSCALARFIELD_HPP
#define MESHSCALARFIELD_HPP 1


#include <iostream>
#include "mesh.hpp"
#include "scalarfield.hpp"


/*! \brief Scalar field class.
 *
 *  %ScalarField is intended to be used as a container for scalar
 *  field data. The IBSimu Package uses scalar fields for storing
 *  electric potential and space charge density fields. A separate
 *  Efield exists for electric fields. The scalar field class
 *  provides a subset of vector operations to copy, sum and scale
 *  scalar fields. The class also includes operators for indexed
 *  access to mesh elements and interpolation functions for linearly
 *  interpolated smooth field data. The interpolation function uses
 *  linear extrapolation of field outside the defined mesh. This way
 *  it can be ensured that %ScalarField returns sensible values even
 *  close by to the edges of the geometry.
 *
 *  The scalar field mesh size has to be at least 2 nodes in active
 *  directions for the selected geometry mode. The mesh size has to be
 *  1 in the other (inactive) directions.
 */
class MeshScalarField : public ScalarField, public Mesh {

    /*! \brief Check mesh definition validity.
     */
    void check_definition();

protected:

    double     *_F;         /*!< \brief Scalar field data */

public:

    /*! \brief Default constructor.
     *
     *  The field made with the default constructor sets geometry mode
     *  to MODE3D, mesh cell size \a h to 1, mesh size \a size to
     *  (0,0,0) and origo \a origo to (0,0,0). The field evaluator
     *  returns always zero.
     */
    MeshScalarField();

    /*! \brief Constructor for scalarfield from \a mesh.
     *
     *  Returns a new scalar field with geometry parameters (including
     *  mesh size) set from \a m. The field is set to zero in all
     *  locations.
     */
    MeshScalarField( const Mesh &m );

    /*! \brief Constructor for set geometry.
     *
     *  Returns a new scalar field with geometry set according to
     *  parameters. The field is set to zero in all locations.
     */
    MeshScalarField( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h );

    /*! \brief Constructor for loading scalar field from a file.
     */
    MeshScalarField( std::istream &s );

    /*! \brief Copy constructor.
     */
    MeshScalarField( const MeshScalarField &f );

    /*! \brief Destructor.
     */
    virtual ~MeshScalarField();

    /*! \brief Clears the field.
     */
    void clear();

    /*! \brief Resets the field geometry.
     *
     *  Sets the field geometry according to the parameters and clears
     *  the field to zero in all locations.
     */
    void reset( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h );

    /*! \brief Search minimum and maximum values of scalar field.
     */
    void get_minmax( double &min, double &max ) const;

    /*! \brief Copy operator.
     */
    MeshScalarField &operator=( const MeshScalarField &f );

    /*! \brief Accumulation operator.
     */
    MeshScalarField &operator+=( const MeshScalarField &f );

    /*! \brief Inverse accumulation operator.
     */
    MeshScalarField &operator-=( const MeshScalarField &f );

    /*! \brief Scaling operator.
     */
    MeshScalarField &operator*=( double x );

    /*! \brief Inverse scaling operator.
     */
    MeshScalarField &operator/=( double x );

    /*! \brief Operator for pointing to elements of the field mesh.
     */
    inline double &operator()( int32_t i ) { 
	return( _F[i] ); }

    /*! \brief Operator for pointing to elements of the field mesh.
     */
    inline double &operator()( int32_t i, int32_t j ) {
	return( _F[i + j*_size(0)] ); }

    /*! \brief Operator for pointing to elements of the field mesh.
     */
    inline double &operator()( int32_t i, int32_t j, int32_t k ) {
	return( _F[i + (j + k*_size(1))*_size(0)] ); }

    /*! \brief Operator for pointing to elements of the field mesh.
     */
    inline const double &operator()( int32_t i ) const { 
	return( _F[i] ); }

    /*! \brief Operator for pointing to elements of the field mesh.
     */
    inline const double &operator()( int32_t i, int32_t j ) const {
	return( _F[i + j*_size(0)] ); }

    /*! \brief Operator for pointing to elements of the field mesh.
     */
    inline const double &operator()( int32_t i, int32_t j, int32_t k ) const {
	return( _F[i + (j + k*_size(1))*_size(0)] ); }

    /*! \brief Operator for getting linearly interpolated field values.
     *
     *  The field is interpolated linearly to get the field value at
     *  \a x. If \a x is outside the mesh, the field is extrapolated
     *  linearly using the field points to \a x. This provides correct
     *  field values also close to the mesh boundaries.
     */
    virtual double operator()( const Vec3D &x ) const;

    /*! \brief Saves scalar field data to file \a filename.
     */
    void save( const std::string &filename ) const;

    /*! \brief Saves scalar field data to stream \a os.
     */
    void save( std::ostream &os ) const;

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;
};


#endif
