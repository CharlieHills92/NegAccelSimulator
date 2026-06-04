/*! \file ccolmatrix.hpp
 *  \brief Compressed column sparse matrix algebra
 */

/* Copyright (c) 2005-2012 Taneli Kalvas. All rights reserved.
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

#ifndef CCOLMATRIX_HPP
#define CCOLMATRIX_HPP 1


#include <cstdlib>
#include <iostream>
#include "matrix.hpp"
#include "error.hpp"


/*! \brief Compressed column sparse matrix class.
 *
 *  The matrix is stored in the standard compressed column sparse matrix
 *  storage mode.  In compressed column storage method all non-zero
 *  matrix elements are stored in array \a val in column-by-column
 *  order. The corresponding row indices are stored in another
 *  array \a row in the same order. The third array \a ptr contains
 *  "pointer" indices indicating start and end of each row in the
 *  first two arrays. The format itself does not require a certain
 *  ordering of elements, but some implementation might need/be faster
 *  on some ordering. Our example matrix \code
 *      | 1  2  0  0  3|
 *      | 4  5  6  0  0|
 *  A = | 0  7  8  0  9|
 *      | 0  0  0 10  0|
 *      |11  0  0  0 12|
 *  \endcode is represented in compressed column sparse matrix class as: \code
 *  ptr[] = {0, 3,  6, 8,  9, 12}
 *  val[] = {1, 4, 11, 2,  5,  7, 6, 8, 10, 3, 9, 12}
 *  row[] = {0, 1,  4, 0,  1,  2, 1, 2,  3, 0, 2,  4}
 *  \endcode
 */
class CColMatrix : public Matrix {
    int       _n;      //!< Number of rows.
    int       _m;      //!< Number of columns.
    int       _nz;     //!< Number of nonzero elements.
    int       _asize;  //!< Allocation size of col and val.
    int      *_ptr;    //!< Column pointers, size m+1.
    int      *_row;    //!< Row indices(i), nz elements in use, asize elements allocated.
    double   *_val;    //!< Element values, nz elements in use, asize elements allocated.
    
    void reallocate( void );
    void allocate( void );

    double get_check( int i, int j ) const;
    double &set_check( int i, int j );
    double get_no_check( int i, int j ) const;
    double &set_no_check( int i, int j );

    void clear_check( int i, int j );
    void clear_no_check( int i, int j );

    void build( const class CColMatrix &mat );
    void build( const class CRowMatrix &mat );
    void build( const class CoordMatrix &mat );

public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Default constructor.
     */
    CColMatrix();

    /*! \brief Constructor to make empty \a n x \a m matrix.
     */
    CColMatrix( int n, int m );

    /*! \brief Constructor to make \a n x \a m matrix from compressed
     *  column matrix data.
     *
     *  The constructed matrix uses \a ptr, \a row and \a val as its
     *  internal data. The arrays are not copied! For memory
     *  allocation compatibility reasons, arrays \a ptr, \a row and \a
     *  val should be allocated using \a malloc and/or \a realloc.
     */
    CColMatrix( int n, int m, int nz, int asize,
		int *ptr, int *row, double *val );

    /*! \brief Copy constructor.
     */
    CColMatrix( const CColMatrix &mat );

    /*! \brief Constructor for conversion from compressed row matrix.
     *
     *  The compressed column matrix will be built in ascending row
     *  index order.
     */
    CColMatrix( const class CRowMatrix &mat );

    /*! \brief Constructor for conversion from coordinate matrix.
     */
    CColMatrix( const class CoordMatrix &mat );

    /*! \brief Constructor for conversion from unknown matrix type.
     */
    CColMatrix( const class Matrix &mat );

    /*! \brief Destructor.
     */
    ~CColMatrix();

/* ************************************** *
 * Access and information                 *
 * ************************************** */

    /*! \brief Returns the number of columns in the matrix.
     */
    int columns( void ) const { return( _m ); }

    /*! \brief Returns the number of rows in the matrix.
     */
    int rows( void ) const { return( _n ); }

    /*! \brief Returns the number of columns and number of columns in \a n and \a m.
     */
    void size( int &n, int &m ) const { n = _n; m = _m; }

    /*! \brief Returns the number of non-zero elements in the matrix.
     */
    int nz_elements( void ) const { return( _nz ); }

    /*! \brief Returns the number of elements allocated for matrix.
     */
    int capacity( void ) const { return( _asize ); }

/* ************************************** *
 * User level control                     *
 * ************************************** */

    /*! \brief Resizes the matrix to size \a n x \a m.
     *
     *  All existing non-zero elements are cleared.
     */
    void resize( int n, int m );

    /*! \brief Merges matrix \a mat into the matrix leaving \a mat empty.
     *
     *  Copies contents of matrix \a mat into the matrix and sets
     *  contents of matrix \a mat to \a n = 0 and \a m = 0.
     *  \param mat Matrix to copy from.
     */
    void merge( CColMatrix &mat );

    /*! \brief Clear non-zero matrix elements, set all elements to zero.
     */
    void clear( void );

    /*! \brief Clear matrix element (i,j).
     *
     *  Removes element (i,j) from the list of non-zero matrix elements.
     */
    void clear( int i, int j );

    /*! \brief Reserve memory for \a size matrix elements.
     */
    void reserve( int size );

    /*! \brief Order (sort) matrix data in ascending row index
     *  order within each column.
     */
    void order_ascending( void );

    /*! \brief Check if matrix data is in ascending row index
     *  order within each column.
     */
    bool check_ascending( void ) const;

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;

/* ************************************** *
 * User level matrix element access       *
 * ************************************** */

    /*! \brief Function to get a matrix element value at (i,j).
     *
     *  Range checking is done for \a i and \a j if \c SPM_RANGE_CHECK
     *  is defined. Throws ErrorRange exception on range checking errors.
     */
    double get( int i, int j ) const;

    /*! \brief Function to get a reference to matrix element value at (i,j).
     *
     *  This function can be used to set or modify matrix element
     *  value. See following examples: \code
     *  A.set(0,0) = 1.2
     *  A.set(0,1) *= 2
     *  A.set(0,1) += 0.1 
     *  \endcode Note that a
     *  reference is actually a pointer to the memory location of the
     *  element and therefore unexpected things can happen if matrix
     *  is modified while using set, for example \code 
     *  A.set(0,0) = A.set(0,1) = A.set(0,2) = 5.0 
     *  \endcode does not do what you
     *  would expect it to do.
     *
     *  Range checking is done for \a i and \a j if \c SPM_RANGE_CHECK
     *  is defined. Throws ErrorRange exception on range checking errors.
     */
    double &set( int i, int j );

    /*! \brief Function to set matrix column elements.
     *
     *  The existing elements of the column \a j are replaced by \a N
     *  elements at row indices described in array \a row and with
     *  values described in array \a val. Range checking is always
     *  done for indexes \a i and \a j. Throws ErrorRange exception on
     *  range checking errors.
     */
    void set_column( int j, int N, const int *row, const double *val );

    /*! \brief Adds an element to matrix while constructing the whole
        matrix.
     *
     *  This is a special function for constructing the whole matrix
     *  column-by-column in ascending column order. The elements
     *  within a column can be defined in any order. With this
     *  function, every column of the matrix has to be defined. No
     *  checking of the definitions are made.
     *
     *  Using this function for defining a large matrix gains
     *  drasticly in speed. The function leaves all but the next
     *  column pointers unmodified. Therefore the matrix is unvalid
     *  and should not be accessed with other functions before all
     *  columns have been defined.
     */
    void construct_add( int i, int j, double val );

/* ************************************** *
 * Low level access                       *
 * ************************************** */

    /*! \brief Returns a reference to the to the internal pointer index
     *  data \a ptr of the matrix.
     */
    int &ptr( int i ) { return( _ptr[i] ); }

    /*! \brief Returns a reference to the to the internal column data
     *  of the matrix.
     */
    int &row( int i ) { return( _row[i] ); }

    /*! \brief Returns a reference to the to the internal value data
     *  of the matrix.
     */
    double &val( int i ) { return( _val[i] ); }

    /*! \brief Returns a const reference to the to the internal
     *  pointer index data \a ptr of the matrix.
     */
    const int &ptr( int i ) const { return( _ptr[i] ); }

    /*! \brief Returns a const reference to the to the internal row
     *  data of the matrix.
     */
    const int &row( int i ) const { return( _row[i] ); }

    /*! \brief Returns a const reference to the to the internal value
     *  data of the matrix.
     */
    const double &val( int i ) const { return( _val[i] ); }

    /*! \brief Set number of non-zero elements in the matrix.
     *
     *  This function is to be used with low level access
     *  functions. The number of non-zero elements should be set to
     *  the same value as \a ptr[m]. Internal arrays are resized if \a
     *  nz is larger than the allocated size.
     */
    void set_nz( int nz );

/* ************************************** *
 * Assignent operators                    *
 * ************************************** */

    /*! \brief Assignment operator.
     */
    CColMatrix &operator=( const CColMatrix &mat );

    /*! \brief Assignment operator for conversion from compressed row
     *  matrix.
     *
     *  The compressed column matrix will be built in ascending row
     *  index order.
     */
    CColMatrix &operator=( const class CRowMatrix &mat );

    /*! \brief Assignment operator for conversion from coordinate
     *  matrix.
     */
    CColMatrix &operator=( const class CoordMatrix &mat );

    /*! \brief Assignment operator for conversion from unknown matrix
     *  type.
     */
    CColMatrix &operator=( const Matrix &mat );

/* ************************************** *
 * Matrix-Vector operations               *
 * ************************************** */

    /*  \brief Calculates \a x = \a A*b.
     */
    void multiply_by_vector( Vector &x, const Vector &b ) const;

    /*! \brief Solves \a A*x = \a b for lower unit diagonal matrix.
     *
     *  Matrix has to have elements only in the lower triangle. Unit
     *  diagonal is implied, it is not to be saved to matrix.
     */
    void lower_unit_solve( Vector &x, const Vector &b ) const;

    /*! \brief Solves \a A*x = \a b for upper diagonal matrix.
     *
     *  Matrix has to have only upper diagonal elements. The diagonal
     *  element has to be the first entry on each column (sorted
     *  ascending order).
     */
    void upper_diag_solve( Vector &x, const Vector &b ) const;


    friend class CRowMatrix;
    friend class CoordMatrix;
    friend class Vector;
    friend class HBIO;
};


inline double CColMatrix::get( int i, int j ) const
{
#ifdef SPM_RANGE_CHECK
    return( get_check( i, j ) );
#else
    return( get_no_check( i, j ) );
#endif
}    


inline double &CColMatrix::set( int i, int j )
{
#ifdef SPM_RANGE_CHECK
    return( set_check( i, j ) );
#else
    return( set_no_check( i, j ) );
#endif
}    


inline void CColMatrix::clear( int i, int j )
{
#ifdef SPM_RANGE_CHECK
    clear_check( i, j );
#else
    clear_no_check( i, j );
#endif
}


#endif






















