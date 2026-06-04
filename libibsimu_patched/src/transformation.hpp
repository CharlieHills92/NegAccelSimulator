/*! \file transformation.hpp
 *  \brief Full transformation for three dimensional homogenous space
 */

/* Copyright (c) 2010-2012,2015 Taneli Kalvas. All rights reserved.
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

#ifndef TRANSFORMATION_HPP
#define TRANSFORMATION_HPP 1


#include <string.h>
#include "vec3d.hpp"
#include "vec4d.hpp"


/*! \brief %Transformation for homogenous three dimensional space.
 *
 *  Transformation for homogenous three dimensional space operates of
 *  4-vectors of type Vec4D. The transformation contains convenience
 *  functions for making affine transformations on 3-vectors of type
 *  Vec3D.
 *
 *  The Transformation object is a 4x4 matrix with convenience
 *  functions for transformation related operations. The matrix is
 *  stored in row first order:
 *  \code
 *   0  1  2  3
 *   4  5  6  7
 *   8  9 10 11
 *  12 13 14 15
 *  \endcode
 */
class Transformation
{

    double x[16];

public:

    /*! \brief Constructor for identity transformation.
     */
    Transformation();

    /*! \brief Constructor for preset transformation matrix.
     */
    Transformation( double x11, double x12, double x13, double x14,
		    double x21, double x22, double x23, double x24,
		    double x31, double x32, double x33, double x34,
		    double x41, double x42, double x43, double x44 );

    /*! \brief Copy constructor.
     */
    Transformation( const Transformation &m );

    /*! \brief Constructor for loading transformation from stream \a is.
     */
    Transformation( std::istream &is );

    /*! \brief Destructor.
     */
    ~Transformation();

    /*! \brief Indexing for transformation matrix.
     */
    double &operator[]( int i ) {
        return( x[i] );
    }

    /*! \brief Indexing for constant transformation matrix.
     */
    const double &operator[]( int i ) const {
        return( x[i] );
    }
    
    /*! \brief Return transpose matrix.
     */
    Transformation transpose( void ) const;

    /*! \brief Return determinant of matrix.
     */
    double determinant( void ) const;

    /*! \brief Return inverse matrix.
     *
     *  Throws an error if determinant is zero.
     */
    Transformation inverse( void ) const;

    /*! \brief Non-equality test
     */
    bool operator!=( const Transformation &m ) const;

    /*! \brief Multiplication of tranformation matrix by scalar.
     */
    const Transformation &operator*=( double s );

    /*! \brief Multiplication of transformation matrices for combining
     *  transformations.
     *
     *  Tranformation is done by multiplying the matrix with a vector
     *  from the right. Therefore the multiplication of transformation
     *  matrices has the effect that the right-hand-side
     *  transformation is applied first and left-hand-side second.
     */
    Transformation operator*( const Transformation &m ) const;

    /*! \brief Multiplication of tranformation matrix by vector.
     *
     *  Makes a full transformation on the homogenous vector \a v.
     */
    Vec4D operator*( const Vec4D &v ) const;

    /*! \brief Multiplication of the transpose of the tranformation
        matrix by vector.
     *
     *  Used for transforming surface normal vectors.
     */
    Vec4D operator%( const Vec4D &v ) const;

    /*! \brief Transform homogenous vector \a xin.
     */
    Vec4D transform( const Vec4D &xin ) const {
	return( *this * xin );
    }

    /*! \brief Transform point \a xin to homogenous space.
     *
     *  Homogenization of output vector is not done.
     */
    Vec4D transform_homogenous_point( const Vec3D &xin ) const {
	return( *this * Vec4D( xin[0], xin[1], xin[2], 1.0 ) );
    }
    
    /*! \brief Transform point \a xin.
     *
     *  Assumes the transformation is affine. Homogenization of output
     *  vector is not done.
     */
    Vec3D transform_point( const Vec3D &xin ) const {
	Vec4D r = *this * Vec4D( xin[0], xin[1], xin[2], 1.0 );	
	return( Vec3D( r[0], r[1], r[2] ) );
    }

    /*! \brief Inverse transform point \a xin.
     *
     *  Assumes the transformation is affine. Homogenization of output
     *  vector is not done.
     *
     *  This is a convenience function to inverting a transformation
     *  matrix and then doing a transform. If more than one transform
     *  is done inverse() and transform() functions should be used.
     */
    Vec3D inv_transform_point( const Vec3D &xin ) const;

    /*! \brief Transform vector \a xin.
     *
     *  Assumes the transformation is affine. Homogenization of output
     *  vector is not done.
     */
    Vec3D transform_vector( const Vec3D &xin ) const {
	Vec4D r = *this * Vec4D( xin[0], xin[1], xin[2], 0.0 );	
	return( Vec3D( r[0], r[1], r[2] ) );
    }

    /*! \brief Inverse transform vector \a xin.
     *
     *  Assumes the transformation is affine. Homogenization of output
     *  vector is not done.
     *
     *  This is a convenience function to inverting a transformation
     *  matrix and then doing a transform. If more than one transform
     *  is done inverse() and transform_vector() functions should be used.
     */
    Vec3D inv_transform_vector( const Vec3D &xin ) const;

    /*! \brief Reset transformation.
     *
     *  Reset transformation to unity.
     */
    void reset( void );

    /*! \brief Translate transformation.
     *
     *  The effect of the new transformation is to first do the old
     *  transformation and then do the translation.
     */
    void translate( const Vec3D &d );

    /*! \brief Translate transformation.
     *
     *  The effect of the new transformation is to first do the
     *  translation and then do the old transformation.
     */
    void translate_before( const Vec3D &d );

    /*! \brief Scale transformation.
     *
     *  The effect of the new transformation is to first do the old
     *  transformation and then do the scaling.
     */
    void scale( const Vec3D &s );

    /*! \brief Scale transformation.
     *
     *  The effect of the new transformation is to first do the
     *  scaling and then do the old transformation.
     */
    void scale_before( const Vec3D &s );

    /*! \brief Rotate transformation around x-axis.
     *
     *  Rotate around x-axis for \a a radians.
     *
     *  The effect of the new transformation is to first do the old
     *  transformation and then do the rotation.
     */
    void rotate_x( double a );

    /*! \brief Rotate transformation around x-axis.
     *
     *  Rotate around x-axis for \a a radians.
     *
     *  The effect of the new transformation is to first do the
     *  rotation and then do the old transformation.
     */
    void rotate_x_before( double a );

    /*! \brief Rotate transformation around y-axis.
     *
     *  Rotate around y-axis for \a a radians.
     *
     *  The effect of the new transformation is to first do the old
     *  transformation and then do the rotation.
     */
    void rotate_y( double a );

    /*! \brief Rotate transformation around y-axis.
     *
     *  Rotate around y-axis for \a a radians.
     *
     *  The effect of the new transformation is to first do the
     *  rotation and then do the old transformation.
     */
    void rotate_y_before( double a );

    /*! \brief Rotate transformation around z-axis.
     *
     *  Rotate around z-axis for \a a radians.
     *
     *  The effect of the new transformation is to first do the old
     *  transformation and then do the rotation.
     */
    void rotate_z( double a );

    /*! \brief Rotate transformation around z-axis.
     *
     *  Rotate around z-axis for \a a radians.
     *
     *  The effect of the new transformation is to first do the
     *  rotation and then do the old transformation.
     */
    void rotate_z_before( double a );

    /*! \brief Return unity transformation.
     */
    static Transformation unity( void );

    /*! \brief Return translation transformation.
     */
    static Transformation translation( const Vec3D &d );

    /*! \brief Return scaling transformation.
     */
    static Transformation scaling( const Vec3D &s );
 
    /*! \brief Return rotation transformation rotating around x-axis.
     */
    static Transformation rotation_x( double a );

    /*! \brief Return rotation transformation rotating around y-axis.
     */
    static Transformation rotation_y( double a );

    /*! \brief Return rotation transformation rotating around z-axis.
     */
    static Transformation rotation_z( double a );

    /*! \brief Outputting to stream.
     */
    friend std::ostream &operator<<( std::ostream &os, const Transformation &t );

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
