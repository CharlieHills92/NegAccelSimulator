/*! \file random.hpp
 *  \brief %Random number generators
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

#ifndef RANDOM_HPP
#define RANDOM_HPP 1


#include <stdint.h>
#include <pthread.h>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_roots.h>



/*! \brief Base class for non-uniform random variate transformation.
 */
class Random_Variate_Transformation {

public:

    /*! \brief Virtual destructor.
     */
    virtual ~Random_Variate_Transformation() {}

    /*! \brief Return a newly allocated copy of object.
     */
    virtual Random_Variate_Transformation *copy( void ) const = 0;

    /*! \brief Returns number from distribution, transformed from
     *  uniformly distributed R, where 0 <= R <= 1.
     */
    virtual double transform( double R ) = 0;
};



/*! \brief Uniform transformation.
 *
 *  Unity transformation.
 */
class Uniform_Transformation : public Random_Variate_Transformation {

public:

    /*! \brief Constructor.
     */
    Uniform_Transformation() {}

    /*! \brief Copy constructor.
     */
    Uniform_Transformation( const Uniform_Transformation &trans ) {}

    /*! \brief Virtual destructor.
     */
    virtual ~Uniform_Transformation() {}

    /*! \brief Assignment
     */
    Uniform_Transformation &operator=( const Uniform_Transformation &trans ) { return( *this ); }

    /*! \brief Return a newly allocated copy of object.
     */
    virtual Uniform_Transformation *copy( void ) const;

    /*! \brief Returns the random variate with no transformation.
     */
    virtual double transform( double R ) { return( R ); }
};



/*! \brief Random variate transformation for Gaussian distribution.
 */
class Gaussian_Transformation : public Random_Variate_Transformation {

    gsl_function_fdf    _fdf;      /*!< \brief Function to solve for gaussian transformation. */
    gsl_root_fdfsolver *_solver;   /*!< \brief Solver for gaussian transformation. */

    static const double _rgauss_const_A;
    static const double _rgauss_const_B;

    static double rgauss_f_func( double x, void *params );
    static double rgauss_df_func( double x, void *params );
    static void rgauss_fdf_func( double x, void *params, double *f, double *df );

public:

    /*! \brief Constructor.
     */
    Gaussian_Transformation();

    /*! \brief Copy constructor.
     */
    Gaussian_Transformation( const Gaussian_Transformation &trans );

    /*! \brief Virtual destructor.
     */
    virtual ~Gaussian_Transformation();

    /*! \brief Assignment
     */
    Gaussian_Transformation &operator=( const Gaussian_Transformation &trans ) { return( *this ); }

    /*! \brief Return a newly allocated copy of object.
     */
    virtual Gaussian_Transformation *copy( void ) const;

    /*! \brief Returns number from distribution, transformed from
     *  uniformly distributed R, where 0 <= R <= 1.
     */
    virtual double transform( double R );
};



/*! \brief Random variate transformation for raised cosine distribution.
 */
class Cosine_Transformation : public Random_Variate_Transformation {

    gsl_function_fdf    _fdf;      /*!< \brief Function to solve for cosine transformation. */
    gsl_root_fdfsolver *_solver;   /*!< \brief Solver for cosine transformation. */

    static double rcosine_f_func( double x, void *params );
    static double rcosine_df_func( double x, void *params );
    static void rcosine_fdf_func( double x, void *params, double *f, double *df );

public:

    /*! \brief Constructor.
     */
    Cosine_Transformation();

    /*! \brief Copy constructor.
     */
    Cosine_Transformation( const Cosine_Transformation &trans );

    /*! \brief Virtual destructor.
     */
    virtual ~Cosine_Transformation();

    /*! \brief Assignment
     */
    Cosine_Transformation &operator=( const Cosine_Transformation &trans ) { return( *this ); }

    /*! \brief Return a newly allocated copy of object.
     */
    virtual Cosine_Transformation *copy( void ) const;

    /*! \brief Returns number from distribution, transformed from
     *  uniformly distributed R, where 0 <= R <= 1.
     */
    virtual double transform( double R );
};



/*! \brief Random variate transformation for gamma distribution.
 */
class Gamma_Transformation : public Random_Variate_Transformation {

    gsl_function_fdf    _fdf;      /*!< \brief Function to solve for cosine transformation. */
    gsl_root_fdfsolver *_solver;   /*!< \brief Solver for cosine transformation. */
    
    double _R;
    double _k;
    double _theta;
    double _rgamma_const_A;
    double _rgamma_const_B;

    static double rgamma_f_func( double x, void *params );
    static double rgamma_df_func( double x, void *params );
    static void rgamma_fdf_func( double x, void *params, double *f, double *df );

public:

    /*! \brief Constructor.
     */
    Gamma_Transformation( double k, double theta );

    /*! \brief Copy constructor.
     */
    Gamma_Transformation( const Gamma_Transformation &trans );

    /*! \brief Virtual destructor.
     */
    virtual ~Gamma_Transformation();

    /*! \brief Assignment
     */
    Gamma_Transformation &operator=( const Gamma_Transformation &trans );

    /*! \brief Return a newly allocated copy of object.
     */
    virtual Gamma_Transformation *copy( void ) const;

    /*! \brief Returns number from distribution, transformed from
     *  uniformly distributed R, where 0 <= R <= 1.
     */
    virtual double transform( double R );
};



/*! \brief %Random number generator for N dimensions.
 *
 *  This random number generator can produce random numbers in N
 *  independent dimensions. The RNG defaults to return uniformly
 *  distributed numbers between 0 and 1. Other distributions can be
 *  sampled by setting transformations. The transformations can be set
 *  independently for each dimension.
 *
 *  Random number generators should be mutex protected.
 */
class Random {

protected:

    mutable pthread_mutex_t _mutex;
    std::vector<Random_Variate_Transformation *> _transformation;

public:

    /*! \brief Constructor for \a N dimensional RNG.
     */
    Random( size_t N );

    /*! \brief Destructor.
     */
    virtual ~Random();

    /*! \brief Set random variate transformation for dimension \a i.
     *
     *  All transformations default to Uniform_Transformation.
     */
    void set_transformation( size_t i, const Random_Variate_Transformation &trans );

    /*! \brief Get a set of random numbers.
     * 
     *  Get next sampling from random number generator to array \a x.
     *  The array x must have space for \a N numbers. 
     */
    virtual void get( double *x ) const = 0;
};


/*! \brief Quasi random number generator for N dimensions.
 *
 *  This random number generator can produce quasi random numbers in N
 *  independent dimensions using the Sobol sequence. Implementation
 *  from GSL library.
 *
 *  Mutex protected.
 */
class QRandom : public Random {

    size_t    _N;        /*!< \brief Number of dimensions. */
    gsl_qrng *_qrng;     /*!< \brief Random number generator from gsl. */

    /*! \brief Prevent copying.
     */
    QRandom( const QRandom &qrng ) : Random(0) {}

public:

    /*! \brief Constructor for random number generator in \a N independent dimensions.
     */
    QRandom( size_t N );

    /*! \brief Virtual destructor.
     */
    virtual ~QRandom();

    /*! \brief Get a set of random numbers.
     * 
     *  Get next sampling from random number generator to array \a x.
     *  The array x must have space for \a N numbers. 
     */
    virtual void get( double *x ) const;
};


/*! \brief Mersenne Twister random number generator for N dimensions.
 *
 *  This RNG can produce random numbers in N independent
 *  dimensions. RNG includes functions to return uniformly distributed
 *  numbers between 0 and 1 and numbers from a gaussian distribution.
 *  Based on MT19937 generator from gsl.
 *
 *  Mutex protected.
 */
class MTRandom : public Random {

    size_t    _N;        /*!< \brief Number of dimensions. */
    gsl_rng  *_rng;      /*!< \brief Random number generator from gsl. */

    /*! \brief Prevent copying
     */
    MTRandom( const MTRandom &rng ) : Random(0) {}

public:

    /*! \brief Constructor for random number generator in \a N independent dimensions.
     */
    MTRandom( size_t N );

    /*! \brief Virtual destructor.
     */
    virtual ~MTRandom();

    /*! \brief Get a set of random numbers.
     * 
     *  Get next sampling from random number generator to array \a x.
     *  The array x must have space for \a N numbers. 
     */
    virtual void get( double *x ) const;
};


#endif

