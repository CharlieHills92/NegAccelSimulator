/*! \file mvector.hpp
 *  \brief N-dimensional vector
 */

/* Copyright (c) 2005-2010 Taneli Kalvas. All rights reserved.
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

#ifndef MVECTOR_HPP
#define MVECTOR_HPP 1


#include <vector>
#include <iostream>
#include "matrix.hpp"
#include "error.hpp"


/*! \brief Dense math vector class.
 *
 *  %Vector is intended to be used for mathematical vectors with large
 *  number of dimensions (n>20). Basic algebra operations are defined
 *  for these vectors, including addition, subtraction, scaling,
 *  different norms and dot product.
 *
 *  The %Vector class does linear algebra operations using %VectorLA
 *  container object to build and store a list of coefficients and
 *  vectors. %VectorLA is used to construct "formulas" from the vector
 *  addition (and subtraction) and multiplication operations. These
 *  formulas are stored without calculating them until a destination
 *  of the result of the operation is known. This is is done to avoid
 *  excess copying and use of temporary variables during calculation.
 *
 *  Uses BLAS if preprocessor directive USE_BLAS is defined, see
 *  config.hpp.
 */
class Vector {
    int       _n;    //!< Number of elements.
    double   *_val;  //!< Element values.

    void allocate( void );
    void callocate( void );
    void reallocate( void );

public:

    /*! \brief Container object for coefficient-vector pairs.
     *
     *  This container object is used to store coefficient and vector
     *  data for linear algebra object VectorLA.
     */
    struct VectorRef {
	const Vector *_vec;  //!< Pointer to vector.
	double        _coef; //!< Coefficient for vector \a vec.

	/*! \brief Constructor for empty VectorRef.
	 */ 
	VectorRef() : _vec(NULL), _coef(1.0) {}

	/*! \brief Constructor for VectorRef with vector \a vec and coefficient \a coef. 
	 */ 
	VectorRef( const Vector *vec, double coef ) : _vec(vec), _coef(coef) {}
    };

    /*! \brief Container object for linear algebra operations.
     *
     *  This container object is used to build and store a list of
     *  coefficients and vectors for linear algebra operations. For
     *  more information about the use of VectorLA, see Vector.
     */
    struct VectorLA {
	std::vector<VectorRef> _refs; //!< List of linear algebra operations.

	/*! \brief Default constructor.
	 */
	VectorLA() {}

	/*! \brief Copy constructor.
	 */
	VectorLA( const VectorLA &vecla ) : _refs(vecla._refs) {}

	/*! \brief Constructor for VectorLA with vector \a vec with coefficient 1.
	 */
	VectorLA( const Vector &vec ) {
	    VectorRef ref( &vec, 1.0 );
	    _refs.push_back( ref );
	}

	/*! \brief Constructor for VectorLA with vector \a vec with coefficient \a coef.
	 */
	VectorLA( const Vector &vec, double coef ) {
	    VectorRef ref( &vec, coef );
	    _refs.push_back( ref );
	}

	/*! \brief Operator for pointing to elements of linear algebra operations.
	 *
	 *  Range checking is done for \a i if \c SPM_RANGE_CHECK is
	 *  defined. Throws MErrorRange exception on range checking
	 *  errors.
	 *
	 *  This operator can be used to calculate algebra for only
	 *  one coordinate of vectors without any excess
	 *  calculation. \code double x = (A-1.3*B)[2] \endcode
	 */
	double operator[]( int i ) const;

	/*! \brief Operator for pointing to elements of linear algebra operations.
	 *
	 *  Range checking is done for \a i if \c SPM_RANGE_CHECK is
	 *  defined. Throws MErrorRange exception on range checking
	 *  errors.
	 */
	double operator()( int i ) const;

	/*! \brief Operator for adding vectors.
	 */
	VectorLA operator+( const VectorLA &vecla ) const;

	/*! \brief Operator for subtracting vectors.
	 */
	VectorLA operator-( const VectorLA &vecla ) const;

	/*! \brief Operator for unary minus.
	 */
	VectorLA operator-() const;

	/*! \brief Operator for multiplying vector with a constant.
	 */
	VectorLA operator*( double x ) const;

	/*! \brief Operator for multiplying vector with a constant.
	 */
	friend VectorLA operator*( double x, const VectorLA &vecla );
    };

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Default constructor.
     *
     *  Returns an empty vector with 0 dimensions.
     */
    Vector() : _n(0), _val(NULL) {}

    /*! \brief Constructor with set dimensionality.
     *
     *  Returns a zero vector with nn dimensions.
     *  \param n Number of dimensions.
     */
    Vector( int n );

    /*! \brief Constructor with preset dimensionality and coordinate values.
     *
     *  Returns a vector with \a n dimensions and coordinate values
     *  copied from array \a val.
     *  \param n Number of dimensions.
     *  \param val Array of coordinate values.
     */
    Vector( int n, const double *val );

    /*! \brief Constructor with preset dimensionality and a value for all coordinates.
     *
     *  Returns a vector with \a n dimensions and all coordinate values
     *  set to \a val.
     *  \param n Number of dimensions.
     *  \param val Number for coordinate values.
     */
    Vector( int n, double val );

    /*! \brief Copy constructor.
     *
     *  Returns a copy of vector \a vec.
     *  \param vec Vector to copy from.
     */
    Vector( const Vector &vec );

    /*! \brief Constructor for setting vector to the result of linear
     *  algebra.
     */
    Vector( const VectorLA &vecla );

    /*! \brief Constructor for setting vector to the result of matrix
     *  multiplication.
     */
    Vector( const struct MatrixMulVec &matvec );

    /*! \brief Destructor for vectors.
     */
    ~Vector();

    /*! \brief Returns the size of vector.
     */
    int size( void ) const { return( _n ); }

    /*! \brief Resizes a vector.
     *
     *  Resizes a vector to \a n dimensions. Contents of vector data
     *  are preserved for where ever it is possible. Newly allocated
     *  areas have uninitialized content. \param n Number of
     *  dimensions for resized vector.
     */
    void resize( int n );

    /*! \brief Clears the vector.
     *
     *  Logically the same as doing \code X = 0 \endcode but clear()
     *  is probably faster.
     */
    void clear( void );

    /*! \brief Merges vector \a vec into the vector leaving \a vec empty.
     *
     *  Copies contents of vector \a vec into the vector and sets
     *  contents of vector \a vec to \a n = 0 and \a val = NULL.
     *  \param vec Vector to copy from.
     */
    void merge( Vector &vec );

    /*! \brief Returns a pointer to the coordinate value data of the vector.
     *
     *  The data is stored in a linear double array, which starts from
     *  the first coordinate and ends to the last one. A non-const
     *  pointer to the first element of the array is returned. 
     */
    double *get_data( void ) { return( _val ); }

    /*! \brief Returns a const pointer to the coordinate value data of
     *  the vector.
     */
    const double *get_data( void ) const { return( _val ); }

    /*! \brief Operator for adding vectors.
     */
    VectorLA operator+( const VectorLA &vecla ) const;

    /*! \brief Operator for subtracting vectors.
     */
    VectorLA operator-( const VectorLA &vecla ) const;

    /*! \brief Operator for unary minus.
     */
    VectorLA operator-() const;

    /*! \brief Operator for multiplying vector with a constant.
     */
    VectorLA operator*( double x ) const;

    /*! \brief Operator for adding vectors.
     */
    Vector &operator+=( const VectorLA &vecla );

    /*! \brief Operator for subtracting vectors.
     */
    Vector &operator-=( const VectorLA &vecla );

    /*! \brief Operator for multiplying vector with a constant.
     */
    Vector &operator*=( double x );

    /*! \brief Operator for setting all vector elements to value \a x.
     */
    Vector &operator=( double x );

    /*! \brief Operator for making vector a copy of another vector.
     */
    Vector &operator=( const Vector &vec );

    /*! \brief Operator for setting vector to the result of linear algebra.
     */
    Vector &operator=( const VectorLA &vecla );

    /*! \brief Operator for setting vector to the result of matrix multiplication
     */
    Vector &operator=( const struct MatrixMulVec &matvec );

    /*! \brief Operator for comparing vectors.
     *
     *  This operator should not be used for general comparison in
     *  applications because of the nature of floating point
     *  numbers. For testing the equality of vectors \a A and \a B you
     *  should do \code max_abs(A - B) < eps \endcode where \a eps is
     *  a small number.
     */
    bool operator==( const Vector &vec ) const;

    /*! \brief Operator for comparing vectors.
     *
     *  This operator should not be used for general comparison in
     *  applications because of the nature of floating point
     *  numbers. For testing the inequality of vectors \a A and \a B
     *  you should do \code max_abs(A - B) > eps \endcode where \a
     *  eps is a small number.
     */
    bool operator!=( const Vector &vec ) const;

    /*! \brief Operator for pointing to vector elements.
     *
     *  Range checking is done for \a i if \c SPM_RANGE_CHECK is
     *  defined. Throws ErrorRange exception on range checking
     *  errors.
     */
    double &operator[]( int i );

    /*! \brief Operator for pointing to vector elements.
     *
     *  Range checking is done for \a i if \c SPM_RANGE_CHECK is
     *  defined. Throws ErrorRange exception on range checking
     *  errors.
     */
    double &operator()( int i );

    /*! \brief Operator for pointing to vector elements.
     *
     *  Range checking is done for \a i if \c SPM_RANGE_CHECK is
     *  defined. Throws ErrorRange exception on range checking
     *  errors.
     */
    double operator[]( int i ) const;

    /*! \brief Operator for pointing to vector elements.
     *
     *  Range checking is done for \a i if \c SPM_RANGE_CHECK is
     *  defined. Throws ErrorRange exception on range checking
     *  errors.
     */
    double operator()( int i ) const;

    friend class HBIO;
    friend class Matrix;
    friend class CRowMatrix;
    friend class CColMatrix;
    friend class CoordMatrix;

    /*! \brief Operator for multiplying vector with a constant.
     */
    friend VectorLA operator*( double x, Vector &vec );

    /*! \brief Operator for printing a vector.
     */
    friend std::ostream &operator<<( std::ostream &os, const Vector &vec );

    /*! \brief Returns dot product of vector \a vec1 and vector \a vec2.
     */
    friend double dot_prod( const Vector &vec1, const Vector &vec2 );

    /*! \brief Returns 1-norm of vector.
     *
     *  \f$ ||x||_1 = \Sigma_{i=1}^n |x_i| \f$
     */
    friend double norm1( const Vector &vec );

    /*! \brief Returns 2-norm of vector
     *
     *  \f$ ||x||_2 = \sqrt{ \Sigma_{i=1}^n x_i^2 } \f$
     */
    friend double norm2( const Vector &vec );

    /*! \brief Returns square of 2-norm of vector
     *
     *  \f$ (||x||_2)^2 = \Sigma_{i=1}^n x_i^2 \f$
     */
    friend double ssqr( const Vector &vec );

    /*! \brief Returns the minimum vector element value.
     */
    friend double min( const Vector &vec );

    /*! \brief Returns the minimum vector element absolute value.
     */
    friend double min_abs( const Vector &vec );

    /*! \brief Returns the maximum vector element value.
     */
    friend double max( const Vector &vec );

    /*! \brief Returns the maximum vector element absolute value.
     */
    friend double max_abs( const Vector &vec );

    /*! \brief Swaps contents of vector \a vec1 and vector \a vec2.
     */
    friend void swap( Vector &vec1, Vector &vec2 );
};


inline double Vector::VectorLA::operator[]( int i ) const {
#ifdef SPM_RANGE_CHECK
    if( i >= _refs[0]._vec->_n )
	throw( ErrorRange( ERROR_LOCATION, i, _refs[0]._vec->_n ) );
#endif
    double res = 0.0;
    std::vector<VectorRef>::const_iterator itend = _refs.end();
    for( std::vector<VectorRef>::const_iterator it = _refs.begin(); it != itend; it++ )
	res += (it->_coef) * (*(it->_vec))[i];
    return( res );
}


inline double Vector::VectorLA::operator()( int i ) const {
#ifdef SPM_RANGE_CHECK
    if( i >= _refs[0]._vec->_n )
	throw( ErrorRange( ERROR_LOCATION, i, _refs[0]._vec->_n ) );
#endif
    double res = 0.0;
    std::vector<VectorRef>::const_iterator itend = _refs.end();
    for( std::vector<VectorRef>::const_iterator it = _refs.begin(); it != itend; it++ )
	res += (it->_coef) * (*(it->_vec))[i];
    return( res );
}


inline double &Vector::operator[]( int i ) {
#ifdef SPM_RANGE_CHECK
    if( i >= _n )
	throw( ErrorRange( ERROR_LOCATION, i, _n ) );
#endif
    return( _val[i] );
}


inline double &Vector::operator()( int i ) {
#ifdef SPM_RANGE_CHECK
    if( i >= _n )
	throw( ErrorRange( ERROR_LOCATION, i, _n ) );
#endif
    return( _val[i] );
}


inline double Vector::operator[]( int i ) const {
#ifdef SPM_RANGE_CHECK
    if( i >= _n )
	throw( ErrorRange( ERROR_LOCATION, i, _n ) );
#endif
    return( _val[i] );
}


inline double Vector::operator()( int i ) const {
#ifdef SPM_RANGE_CHECK
    if( i >= _n )
	throw( ErrorRange( ERROR_LOCATION, i, _n ) );
#endif
    return( _val[i] );
}





#endif




















