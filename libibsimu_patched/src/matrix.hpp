/*! \file matrix.hpp
 *  \brief Basis for matrix implementations
 */

/* Copyright (c) 2005-2010,2012 Taneli Kalvas. All rights reserved.
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

#ifndef MATRIX_HPP
#define MATRIX_HPP 1


#include <iostream>
#include "mvector.hpp"


/*! \brief Container object for matrix-vector multiplication operation.
 *
 *  This container object is used to store a matrix-vector
 *  multiplication operation. For more information about the use of
 *  MatrixMulVec, see MCRowMatrix.
 */
struct MatrixMulVec {
    const class Matrix  *_mat; //!< Pointer to matrix.
    const class Vector  *_vec; //!< Pointer to vector.

    /*! \brief Constructor for MMatrixMulVec with matrix \a mat and
     *  vector \a vec.
     */
    MatrixMulVec( const Matrix &mat, const class Vector &vec ) : 
	_mat(&mat), _vec(&vec) {}

    friend class Vector;
};


/*! \brief Base matrix class.
 *
 *  The matrix class is an abstract class designed to be used as a
 *  base class for different matrix implementations.
 */
class Matrix {
    virtual double get_check( int i, int j ) const = 0;
    virtual double &set_check( int i, int j ) = 0;
    virtual double get_no_check( int i, int j ) const = 0;
    virtual double &set_no_check( int i, int j ) = 0;

public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Virtual destructor.
     */
    virtual ~Matrix() {}

/* ************************************** *
 * Access and information                 *
 * ************************************** */

    /*! \brief Returns the number of columns of the matrix.
     */
    virtual int columns( void ) const = 0;

    /*! \brief Returns the number of rows of the matrix.
     */
    virtual int rows( void ) const = 0;

    /*! \brief Returns the number of rows \a n and the number of
     *  columns \a m of the matrix.
     */
    virtual void size( int &n, int &m ) const = 0;

/* ************************************** *
 * User level control                     *
 * ************************************** */

    /*! \brief Resizes the matrix to \a nn x \a mm.
     */
    virtual void resize( int n, int m ) = 0;

    //virtual void merge( Matrix &mat ) = 0;

    /*! \brief Clears the matrix (sets all element to zero).
     */
    virtual void clear( void ) = 0;

/* ************************************** *
 * User level matrix element access       *
 * ************************************** */

    /*! \brief Function to get a matrix element value at (i,j).
     */
    inline double get( int i, int j ) const;

    /*! \brief Function to get a reference to matrix element value at
     *  (i,j).
     */
    inline double &set( int i, int j );

/* ************************************** *
 * Matrix-Vector operations               *
 * ************************************** */

    /*! \brief Operator for matrix-vector multiplication.
     */
    MatrixMulVec operator*( const class Vector &vec ) const;

    /*  \brief Calculates \a x = \a A*b.
     *
     *  Called by Vector through MatrixMulVec for efficient use of
     *  operator* in matrix-vector multiplication.
     */
    virtual void multiply_by_vector( Vector &res, const Vector &rhs ) const = 0;

    /*! \brief Solves \a A*x = \a b for lower unit diagonal matrix.
     *
     *  Matrix has to have elements only in the lower triangle. Unit
     *  diagonal is implied, it is not to be saved to matrix.
     */
    virtual void lower_unit_solve( Vector &y, const Vector &b ) const = 0;

    /*! \brief Solves \a A*x = \a b for upper diagonal matrix.
     *
     *  Matrix has to have only upper diagonal elements. The diagonal
     *  element has to be the first entry on each column (sorted
     *  ascending order).
     */
    virtual void upper_diag_solve( Vector &x, const Vector &y ) const = 0;

    friend class Vector;
};


inline double Matrix::get( int i, int j ) const
{
#ifdef SPM_RANGE_CHECK
    return( get_check( i, j ) );
#else
    return( get_no_check( i, j ) );
#endif
}    


inline double &Matrix::set( int i, int j )
{
#ifdef SPM_RANGE_CHECK
    return( set_check( i, j ) );
#else
    return( set_no_check( i, j ) );
#endif
}    


std::ostream &operator<<( std::ostream &os, const Matrix &mat );


#endif
