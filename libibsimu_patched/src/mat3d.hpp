/*! \file mat3d.hpp
 *  \brief Three-by-three matrices.
 */

/* Copyright (c) 2010-2012 Taneli Kalvas. All rights reserved.
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

#ifndef MAT3D_HPP
#define MAT3D_HPP 1


#include "vec3d.hpp"


/*! \brief Three-by-three matrix.
 *
 *  Dense three-by-three matrix. Matrix data in row first order: \code
 *      | 0 1 2 |
 *  A = | 3 4 5 |
 *      | 6 7 8 |
 *  \endcode
 */
class Mat3D {

    double a[9];

public:

    /*! \brief Constructor for zero matrix.
     */
    Mat3D();

    /*! \brief Constructor for matrix with preset elements.
     */
    Mat3D( double a11, double a12, double a13,
	   double a21, double a22, double a23,
	   double a31, double a32, double a33 );

    /*! \brief Destructor.
     */
    ~Mat3D() {}


    /*! \brief Indexing for matrix.
     *
     *  Returns a reference to \a i element, where \a i is the index
     *  in row-first order. No checking performed - not a safe
     *  function.
     */
    double &operator()( int i ) {
        return( a[i] );
    }

    /*! \brief Indexing for constant matrix.
     *
     *  Returns a reference to \a i element, where \a i is the index
     *  in row-first order. No checking performed - not a safe
     *  function.
     */
    const double &operator()( int i ) const {
        return( a[i] );
    }
    
    /*! \brief Indexing for matrix.
     *
     *  Returns a reference to \a (i,j) element, where \a i is the row
     *  and \a j is the column number. No checking performed - not a
     *  safe function.
     */
    double &operator()( int i, int j ) {
        return( a[3*i+j] );
    }

    /*! \brief Indexing for matrix.
     *
     *  Returns a reference to \a (i,j) element, where \a i is the row
     *  and \a j is the column number. No checking performed - not a
     *  safe function.
     */
    const double &operator()( int i, int j ) const {
        return( a[3*i+j] );
    }

    /*! \brief Return determinant of matrix.
     */
    double determinant( void ) const;

    /*! \brief Return inverse matrix.
     *
     *  Throws an error if determinant is zero. For more thorough
     *  checking (determinant close to zero, for example), see the
     *  other form of inverse().
     */
    Mat3D inverse( void ) const;

    /*! \brief Return inverse matrix.
     *
     *  Calculate inverse matrix. Use pre-calculated determinant \a det.
     */
    Mat3D inverse( double det ) const;

    /*! \brief Matrix-vector multiplication.
     */
    Vec3D operator*( const Vec3D &x ) const;

    /*! \brief Outputting to stream.
     */
    friend std::ostream &operator<<( std::ostream &os, const Mat3D &m );
};


#endif

