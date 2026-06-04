/*! \file interpolation.hpp
 *  \brief Two dimensional interpolation
 */

/* Copyright (c) 2005-2012,2014 Taneli Kalvas. All rights reserved.
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

#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP 1


#include <vector>
#include <cstddef>


/*! \brief Base class for 2d interpolation.
 *
 *  Provides an interpolation of a function defined at mesh points \a
 *  (x,y), where 0 <= (x,y) <= 1. 
 */
class Interpolation2D
{
protected:

    size_t                    _n;    /*!< \brief Size of first coordinate of mesh. */
    size_t                    _m;    /*!< \brief Size of second coordinate of mesh. */
    std::vector<double>       _f;    /*!< \brief Function data mesh. */

    /*! \brief Constructor.
     *
     *  Makes an independent object for interpolation of data. An
     *  internal copy of \a f is created. Data mesh is assumed to be
     *  accessed with indexing \a i+j*n.
     */
    Interpolation2D( size_t n, size_t m, const std::vector<double> &f );

    const double &__f( int i, int j ) const;
    double &__f( int i, int j );

public:

    /*! \brief Virtual destructor.
     */
    virtual ~Interpolation2D() {}

    /*! \brief Operator for getting interpolation at \a (x,y).
     *
     *  Returns an interpolated value of the function at \a (x,y),
     *  where 0 <= (x,y) <= 1. Returns NaN outside this area.
     */
    virtual double operator()( double x, double y ) const = 0;
};


/*! \brief Closest point 2d interpolation.
 *
 *  Not really an interpolation. Just returns the closest point of
 *  original data.
 */
class ClosestInterpolation2D : public Interpolation2D
{

public:

    /*! \brief Constructor.
     *
     *  Makes an independent object for interpolation of data. An
     *  internal copy of \a f is created. Data mesh is assumed to be
     *  accessed with indexing \a i+j*n.
     */
    ClosestInterpolation2D( size_t n, size_t m, const std::vector<double> &f );

    /*! \brief Destructor.
     */
    virtual ~ClosestInterpolation2D() {}

    /*! \brief Operator for getting interpolation at \a (x,y).
     *
     *  Returns an interpolated value of the function at \a (x,y),
     *  where 0 <= (x,y) <= 1. Returns NaN outside this area.
     */
    virtual double operator()( double x, double y ) const;
};


/*! \brief BiLinear 2d interpolation.
 */
class BiLinearInterpolation2D : public Interpolation2D
{

public:

    /*! \brief Constructor.
     *
     *  Makes an independent object for interpolation of data. An
     *  internal copy of \a f is created. Data mesh is assumed to be
     *  accessed with indexing \a i+j*n.
     */
    BiLinearInterpolation2D( size_t n, size_t m, const std::vector<double> &f );

    /*! \brief Destructor.
     */
    virtual ~BiLinearInterpolation2D() {}

    /*! \brief Operator for getting interpolation at \a (x,y).
     *
     *  Returns an interpolated value of the function at \a (x,y),
     *  where 0 <= (x,y) <= 1. Returns NaN outside this area.
     */
    virtual double operator()( double x, double y ) const;
};


/*! \brief BiCubic 2d interpolation.
 *
 *  Calculates the derivatives of the function at mesh points using
 *  central finite differences. Zero derivatives are assumed at boundaries.
 */
class BiCubicInterpolation2D : public Interpolation2D
{

    std::vector<double>       _fx;   /*!< \brief X-derivative of function, data mesh. */
    std::vector<double>       _fy;   /*!< \brief Y-derivative of function, data mesh. */
    std::vector<double>       _fxy;  /*!< \brief XY-derivative of function, data mesh. */
    std::vector<double>       _c;    /*!< \brief Precalculated coefficients for interpolation. 
				     *
				     * 16 numbers per mesh square. Totalling 16*(n-1)*(m-1).
				     */

    const double &__fx( int i, int j ) const;
    const double &__fy( int i, int j ) const;
    const double &__fxy( int i, int j ) const;

    double &__fx( int i, int j );
    double &__fy( int i, int j );
    double &__fxy( int i, int j );

    static void calc_coefs( double *c, double *x );
    static const double wt[16][16];

public:

    /*! \brief Constructor.
     *
     *  Makes an independent object for interpolation of data. An
     *  internal copy of \a f is created. Data mesh is assumed to be
     *  accessed with indexing \a i+j*n.
     */
    BiCubicInterpolation2D( size_t n, size_t m, const std::vector<double> &f );

    /*! \brief Destructor.
     */
    virtual ~BiCubicInterpolation2D() {}

    /*! \brief Operator for getting interpolation at \a (x,y).
     *
     *  Returns an interpolated value of the function at \a (x,y),
     *  where 0 <= (x,y) <= 1. Returns NaN outside this area.
     */
    virtual double operator()( double x, double y ) const;
};


#endif
