/*! \file meshvectorfield.hpp
 *  \brief %Mesh based vector fields
 */

/* Copyright (c) 2005-2011,2014,2015 Taneli Kalvas. All rights reserved.
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

#ifndef MESHVECTORFIELD_HPP
#define MESHVECTORFIELD_HPP 1


#include "transformation.hpp"
#include "vectorfield.hpp"
#include "mesh.hpp"
#include "types.hpp"


/*! \brief %Mesh based vector field.
 *
 *  The mesh based vector field is a vector field implementation
 *  storing vector data in a even, rectangular mesh. The %IBSimu
 *  Package uses vector fields for electric and magnetic fields. The
 *  vector field class provides a subset of vector operations to copy,
 *  sum and scale vector fields. The class also includes operators for
 *  indexed access to mesh elements and interpolation functions for
 *  linearly interpolated smooth field data. The behaviour of the
 *  interpolation function can be selected outside the defined
 *  mesh. The default behaviour is to use linear extrapolation. This
 *  way it can be ensured that %VectorField returns sensible values
 *  even close by to the edges of the mesh.
 *
 *  If the size of mesh is 1 in some direction, then the field is
 *  constant in that direction. Otherwise linear interpolation is
 *  used.
 *
 *  The mesh of the magnetic field can be selected independently of
 *  the electric field in ParticleDataBase. This allows minimization
 *  of the memory use in the case where electric field needs high
 *  resolution, but magnetic field is relatively smooth.
 *
 *  The mesh based bector field can be rotated, scaled and
 *  translated. This is implemented by storing the active
 *  transformation in the %MeshVectorField. Each time the field value
 *  is queried at a point, the coordinates are transformed to the
 *  local coordinate system of the %MeshVectorField using the active
 *  transformation.
 *  
 */
class MeshVectorField : public VectorField, public Mesh {

    field_extrpl_e   _extrpl[6]; /*!< \brief What to return outside geometry. */
    double          *_F[3];      /*!< \brief Vector field data in three components
				  *
				  *   If pointer in array is NULL the component 
				  *   is not stored.
				  */

    Transformation   _T;         /*!< \brief Field transformation. */
    Transformation   _Tinv;      /*!< \brief Field transformation inversed. */

    void check_definition();

    void convert_3d_to_3d( const MeshVectorField &fin );
    void convert_cyl_to_3d( const MeshVectorField &fin );

public:


    /*! \brief Default constructor.
     *
     *  The field made with the default constructor sets geometry mode
     *  to MODE3D, mesh cell size \a h to 1, mesh size \a size to
     *  (0,0,0) and origo \a origo to (0,0,0). The field evaluator
     *  returns always zero.
     */
    MeshVectorField();

    /*! \brief Constructor for vector field from \a m.
     *
     *  Returns a new vector field with geometry parameters (including
     *  mesh size) set from \a m. The field is set to zero in all
     *  locations.
     */
    MeshVectorField( const Mesh &m, const bool fout[3] );

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
    MeshVectorField( geom_mode_e geom_mode, const bool fout[3], Int3D size, 
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
    MeshVectorField( geom_mode_e geom_mode, const bool fout[3], double xscale, 
		     double fscale, const std::string &filename );

    /*! \brief Conversion constructor.
     *
     *  Returns a new vector field with geometry set according to
     *  parameters: \a geom_mode is the geometry mode, \a size is the
     *  size of the mesh, \a origo is the location of mesh point
     *  (0,0,0) and \a h is the mesh cell size. The vector field
     *  components marked \a true in array fout are to be defined in
     *  the vector field. Components marked \a false are always
     *  zero. The field content is copied from another mesh based
     *  vector field \a fin.
     *
     *  Currently supports 
     *
     *  1. Conversion from cylindrical (x,r,Bx,Br,Btheta) to 3d 
     *  (x,y,z,Bx,By,Bz), where x -> z and r -> (x,y).
     *
     *  2. Conversion from 3d to 3d. Allows change of mesh density and 
     *  size of field, change of extrapolation effect, etc.
     *
     *  Conversion algorithm uses field evaluator of the input field
     *  and therefore the extrapolation settings and the
     *  transformation in the field affect the newly created field.
     */
    MeshVectorField( geom_mode_e geom_mode, const bool fout[3], Int3D size, 
		     Vec3D origo, double h, const MeshVectorField &fin );

    /*! \brief Copy constructor.
     */
    MeshVectorField( const MeshVectorField &f );

    /*! \brief Constructor for loading vector field from stream \a is.
     */
    MeshVectorField( std::istream &is );

    /*! \brief Destructor.
     */
    virtual ~MeshVectorField();

    /*! \brief Set the behaviour of field interpolation outside mesh
     *  points (extrapolation).
     *
     *  The interpolation function behaviour can be set separately for
     *  each boundary. This is done by setting the desired properties
     *  to the \a extrpl array. The interpolation function can use an
     *  extrapolation of the last two field values (\a
     *  FIELD_EXTRAPOLATE) or it can return the mirror of the field
     *  across the boundary as if the boundary was a perfect magnetic
     *  conductor, i.e. the field lines are perpendicular to the
     *  boundary (\a FIELD_MIRROR). It can also return the mirror of
     *  the field across the mesh boundary as if the boundary was a
     *  magnetic insulator, i.e. the field lines are parallel to the
     *  boundary (\a FIELD_ANTIMIRROR) or it can return a zero field
     *  (\a FIELD_ZERO) NaN (\a FIELD_NAN) outside the mesh. The \a
     *  FIELD_EXTRAPOLATE is the default behaviour.
     *
     *  The \a FIELD_MIRROR extrapolation mirrors the field for a boundary \a x=0 as
     *
     *      B_x(-x,y,z) =  B_x(x,y,z),
     *      B_y(-x,y,z) = -B_y(x,y,z),
     *      B_z(-x,y,z) = -B_z(x,y,z)
     *
     *  The \a FIELD_ANTIMIRROR extrapolation mirrors the field for a boundary \a x=0 as
     *
     *      B_x(-x,y,z) = -B_x(x,y,z),  
     *      B_y(-x,y,z) =  B_y(x,y,z),  
     *      B_z(-x,y,z) =  B_z(x,y,z)
     *
     *  The use of \a FIELD_MIRROR in case of symmetric cases, where
     *  beam is traversing next to the geometry boundary, is necessary
     *  to get physical results.
     *
     *  Very far (double the size of the mesh volume) the field
     *  evaluator will always return zero.
     */
    void set_extrapolation( const field_extrpl_e extrpl[6] );

    /*! \brief Set transformation to unity.
     */
    void reset_transformation( void );

    /*! \brief Set transformation as a copy of \a T.
     */
    void set_transformation( const Transformation &T );

    /*! \brief Translate field.
     */
    void translate( const Vec3D &dx );

    /*! \brief Scale field.
     */
    void scale( const Vec3D &sx );

    /*! \brief Rotate field in coordinate system around y-axis for \a a radians.
     */
    void rotate_x( double a );

    /*! \brief Rotate field in coordinate system around y-axis for \a a radians.
     */
    void rotate_y( double a );

    /*! \brief Rotate field in coordinate system around y-axis for \a a radians.
     */
    void rotate_z( double a );

    /*! \brief Clears the field.
     */
    void clear();

    /*! \brief Resets the field geometry.
     *
     *  Sets the field mesh geometry according to the parameters,
     *  clears the field to zero in all locations and resets the
     *  transformation.
     */
    void reset( geom_mode_e geom_mode, const bool fout[3], Int3D size, 
		Vec3D origo, double h );

    /*! \brief Search minimum and maximum vector length values of
     *  vector field.
     */
    void get_minmax( double &min, double &max ) const;

    /*! \brief Search minimum and maximum for each vector component.
     *
     *  Return vector min with minimum components found and vector max
     *  with maximum components found.
     */
    void get_minmax( Vec3D &min, Vec3D &max ) const;

    /*! \brief Get which field components are defined.
     */
    void get_defined_components( bool fout[3] ) const;

    /*! \brief Copy operator.
     */
    MeshVectorField &operator=( const MeshVectorField &f );

    /*! \brief Accumulation operator.
     *
     *  This operator is only capable of accumulating fields with
     *  matching meshes and tranformations.
     */
    MeshVectorField &operator+=( const MeshVectorField &f );

    /*! \brief Scaling operator for field values.
     */
    MeshVectorField &operator*=( double x );

    /*! \brief Inverse scaling operator.
     */
    MeshVectorField &operator/=( double x );

    /*! \brief Operator for getting elements of the field mesh.
     */
    const Vec3D operator()( int32_t i ) const;

    /*! \brief Operator for getting elements of the field mesh.
     */
    const Vec3D operator()( int32_t i, int32_t j ) const;

    /*! \brief Operator for getting elements of the field mesh.
     */
    const Vec3D operator()( int32_t i, int32_t j, int32_t k ) const;

    /*! \brief Operator for setting element (i) of the field mesh.
     */
    void set( int32_t i, const Vec3D &v );
    
    /*! \brief Operator for setting element (i,j) of the field mesh.
     */
    void set( int32_t i, int32_t j, const Vec3D &v );

    /*! \brief Operator for setting element (i,j,k) of the field mesh.
     */
    void set( int32_t i, int32_t j, int32_t k, const Vec3D &v );

    /*! \brief Operator for getting linearly interpolated field value
     *  at \a x.
     */
    virtual const Vec3D operator()( const Vec3D &x ) const;

    /*! \brief Saves data to a new file \a filename.
     */
    void save( const std::string &filename ) const;

    /*! \brief Saves vector field data to stream \a os.
     */
    void save( std::ostream &os ) const;

    /*! \brief Print debugging information to stream \a os.
     */
    void debug_print( std::ostream &os ) const;
};


#endif

