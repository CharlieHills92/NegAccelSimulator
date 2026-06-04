/*! \file multimeshvectorfield.hpp
 *  \brief %Vector field using multiple meshes.
 */

/* Copyright (c) 2011,2012,2014 Taneli Kalvas. All rights reserved.
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

#ifndef MULTIMESHVECTORFIELD_HPP
#define MULTIMESHVECTORFIELD_HPP 1


#include <vector>
#include "meshvectorfield.hpp"


/*! \brief %Vector field based on multiple meshes.
 *
 *  The MultiMeshVectorField is a vector field implementation storing
 *  vector data in multiple rectangular meshes. This makes it possible
 *  to have higher precision in some areas and lower precision in
 *  other areas of the simulation. The multiple mesh vector field has
 *  several underlying mesh based vector fields in a certain
 *  order. When the field value is queried at a location, the
 *  underlying fields are queries in the defined order until a field
 *  covering the point of query is found. Is the point is outside all
 *  of the underlying fields, the last (largest and coarsest) field is
 *  used to extrapolate the field value (according to the
 *  extrapolation settings).
 *
 *  The first (largest and coarsest) field of the multiple mesh vector
 *  field is defined first either with a constructor or using reset()
 *  function. The detailed meshes are then added with add_mesh().
 *
 *  The MultiMeshVectorField can't be yet used with electric potential
 *  solvers. This functionality has been planned to be added in the
 *  future for possibility of having a higher mesh density in critical
 *  areas. Currently the MultiMeshVectorFields can only be used as
 *  magnetic fields.
 *
 */
class MultiMeshVectorField : public VectorField {

    std::vector<MeshVectorField *> _field; /*!< \brief List of underlying fields. 
					    *
					    *   First one is the coarsest and must exist.
					    */
public:

    /*! \brief Default constructor.
     *
     *  The field made with the default constructor sets geometry mode
     *  to MODE3D, mesh cell size \a h to 1, mesh size \a size to
     *  (0,0,0), origo \a origo to (0,0,0). The field evaluator
     *  returns always zero.
     */
    MultiMeshVectorField();

    /*! \brief Constructor for vector field from \a m.
     *
     *  Returns a new vector field with geometry parameters (including
     *  mesh size) set from \a m. The field is set to zero in all
     *  locations.
     */
    MultiMeshVectorField( const Mesh &m, const bool fout[3] );

    /*! \brief Constructor for set geometry.
     *
     *  Returns a new vector field with geometry set according to
     *  parameters: \a geom_mode is the geometry mode, \a size is the
     *  size of the mesh, \a origo is the location of mesh point
     *  (0,0,0) and \a h is the mesh cell size. The vector field
     *  components marked \a true in array fout are to be defined in
     *  the vector field. Components marked \a false are always
     *  zero. The field is initially set to zero in all locations.
     */
    MultiMeshVectorField( geom_mode_e geom_mode, const bool fout[3], Int3D size, 
			  Vec3D origo, double h );

    /*! \brief Constructor for vector field from ascii file.
     *
     *  The vector field for geometry mode \a geom_mode is read in
     *  from file \a filename. The lines starting with # are
     *  skipped. After that the data is read in line-by-line with one
     *  data point per line. The data columns are separated by white
     *  space. The coordinate data are \a (x, y) in 2D, \a (x, r) in
     *  Cyl and \a (x, y, z) in 3D. The field data to be read are
     *  enabled by user with \a fout. The enabled field data
     *  components are read in after the coordinate data from the data
     *  line.
     *
     *  The data points are expected to appear in coordinate sorted
     *  order because the mesh step h is determined from the spatial
     *  difference of first two data points. Spatial coordinates are
     *  multiplied with \a xscale and field components with \a fscale
     *  while read in.
     *
     *  For magnetic fields the particle iterator assumes vector 
     *  field in the following formats:
     *  In 2D:  (x, y, Bz)
     *  In Cyl: (x, r, Bx, Br, Btheta)
     *  In 3D: (x, y, z, Bx, By, Bz)
     *
     */
    MultiMeshVectorField( geom_mode_e geom_mode, const bool fout[3], double xscale, 
			  double fscale, const std::string &filename );

    /*! \brief Copy constructor.
     */
    MultiMeshVectorField( const MultiMeshVectorField &f );

    /*! \brief Constructor for loading vector field from a file.
     */
    MultiMeshVectorField( std::istream &s );

    /*! \brief Destructor.
     */
    virtual ~MultiMeshVectorField();

    /*! \brief Set the behaviour of field interpolation outside mesh
     *  points (extrapolation).
     *
     *  The interpolation function behaviour can be set separately for
     *  each boundary. This is done by setting the desired properties
     *  to the \a extrpl array. The interpolation function can use an
     *  extrapolation of the last two field values (\a
     *  FIELD_EXTRAPOLATE) or it can return the mirror of the field
     *  across the mesh boundary (\a FIELD_MIRROR), can return a zero
     *  field (\a FIELD_ZERO) or it can return a NaN (\a FIELD_NAN)
     *  outside the mesh. The \a FIELD_EXTRAPOLATE is the default behaviour.
     *
     *  The use of \a FIELD_MIRROR in case of symmetric cases, where
     *  beam is traversing next to the geometry boundary, is necessary
     *  to get physical results.
     *
     *  Very far (double the size of the simulation box) the field
     *  evaluator will always return zero.
     */
    void set_extrapolation( const field_extrpl_e extrpl[6] );

    /*! \brief Translate field in coordinate system.
     */
    void translate( Vec3D x );

    /*! \brief Scale field in coordinate system.
     */
    void scale( double s );

    /*! \brief Rotate field in coordinate system around x-axis for \a a radians.
     */
    void rotate_x( double a );

    /*! \brief Rotate field in coordinate system around y-axis for \a a radians.
     */
    void rotate_y( double a );

    /*! \brief Rotate field in coordinate system around z-axis for \a a radians.
     */
    void rotate_z( double a );

    /*! \brief Clears the field.
     */
    void clear();

    /*! \brief Resets the field geometry.
     *
     *  Sets the field geometry according to the parameters and clears
     *  the field to zero in all locations. Only the coarsest field
     *  will be defined using these parameters.
     */
    void reset( geom_mode_e geom_mode, const bool fout[3], Int3D size, 
		Vec3D origo, double h );

      
    /*! \brief Add a new mesh field to the multiple mesh vector field.
     *
     *  The new field, which is an internal copy of the field given as
     *  a parameter to this function is appended to the list of
     *  fields. The field defined last is checked first by the
     *  evaluator when searching for field value at a point. The new
     *  field has to have the same defined field components and same
     *  geometry mode as the first field defined.
     */  
    void add_mesh( const MeshVectorField &field );

    /*! \brief Add a new mesh field to the multiple mesh vector field.
     *
     *  The new field is appended to the list of fields. The field
     *  defined last is checked first by the evaluator when searching
     *  for field value at a point. The new field has same defined
     *  field components and same geometry mode as the first field
     *  defined. The \a size, \a origo and mesh step \a h can differ.
     */
    void add_mesh( Int3D size, Vec3D origo, double h );
    
    /*! \brief Add a new mesh field from ascii file.
     *
     *  The new field is appended to the list of fields. The field
     *  defined last is checked first by the evaluator when searching
     *  for field value at a point. The new field has same defined
     *  field components and same geometry mode as the first field
     *  defined. The \a size, \a origo and mesh step \a h can differ.
     */
    void add_mesh( double xscale, double fscale, const std::string &filename );

    /*! \brief Search minimum and maximum vector length values of
     *  vector field.
     */
    void get_minmax( double &min, double &max ) const;

    /*! \brief Get which field components are defined.
     */
    void get_defined_components( bool fout[3] ) const;

    /*! \brief Copy operator.
     */
    MultiMeshVectorField &operator=( const MultiMeshVectorField &f );

    /*! \brief Return const reference to subfield \a i.
     */
    const MeshVectorField &operator[]( int i ) const {
	return( *_field[i] );
    }

    /*! \brief Return reference to subfield \a i.
     */
    MeshVectorField &operator[]( int i ) {
	return( *_field[i] );
    }

    /*! \brief Operator for getting linearly interpolated field value
     *  at \a x.
     */
    virtual const Vec3D operator()( const Vec3D &x ) const;

    /*! \brief Saves data to a new file \a filename.
     */
    void save( const std::string &filename ) const;

    /*! \brief Saves vector field data to stream.
     */
    void save( std::ostream &s ) const;

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;
};


#endif

