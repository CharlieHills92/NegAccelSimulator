/*! \file random.cpp
 *  \brief %Random number generators
 */

/* Copyright (c) 2005-2009,2011 Taneli Kalvas. All rights reserved.
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


#include <iostream>
#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include "random.hpp"
#include "error.hpp"


Uniform_Transformation *Uniform_Transformation::copy( void ) const
{
    return( new Uniform_Transformation( *this) );
}



/* ********************************************************************************* */



const double Gaussian_Transformation::_rgauss_const_A = 1.0/sqrt(2.0);
const double Gaussian_Transformation::_rgauss_const_B = 1.0/sqrt(2.0*M_PI);


Gaussian_Transformation::Gaussian_Transformation()
{
    _fdf.f   = &rgauss_f_func;
    _fdf.df  = &rgauss_df_func;
    _fdf.fdf = &rgauss_fdf_func;
    _solver  = gsl_root_fdfsolver_alloc( gsl_root_fdfsolver_newton );
}


Gaussian_Transformation::Gaussian_Transformation( const Gaussian_Transformation &trans )
{
    _fdf.f   = &rgauss_f_func;
    _fdf.df  = &rgauss_df_func;
    _fdf.fdf = &rgauss_fdf_func;
    _solver  = gsl_root_fdfsolver_alloc( gsl_root_fdfsolver_newton );
}


Gaussian_Transformation::~Gaussian_Transformation()
{
    gsl_root_fdfsolver_free( _solver );
}


double Gaussian_Transformation::rgauss_f_func( double x, void *params )
{
    double *R = (double *)params;
    return( 0.5*( 1.0+erf(x*_rgauss_const_A) ) - *R );
}


double Gaussian_Transformation::rgauss_df_func( double x, void *params )
{
    return( _rgauss_const_B*exp( -0.5*x*x ) );
}


/* Subroutine providing cumulative distribution function of gaussian
 * and its derivative at x. */
void Gaussian_Transformation::rgauss_fdf_func( double x, void *params, double *f, double *df )
{
    double *R = (double *)params;

    *f  = 0.5*( 1.0+erf(x*_rgauss_const_A) ) - *R;
    *df = _rgauss_const_B*exp( -0.5*x*x );
}


Gaussian_Transformation *Gaussian_Transformation::copy( void ) const
{
    return( new Gaussian_Transformation( *this) );
}


double Gaussian_Transformation::transform( double R )
{
    int iter = 0;
    double x, x0;

    x = 0.0;
    _fdf.params = (void *)&R;
    gsl_root_fdfsolver_set( _solver, &_fdf, x );
    do {
        gsl_root_fdfsolver_iterate( _solver );
        x0 = x;
        x = gsl_root_fdfsolver_root( _solver );
        if( fabs(x - x0) < 1e-6 )
	    break;
	iter++;
    } while( iter < 100 );
    if( iter == 100 )
        throw( Error( ERROR_LOCATION, "too many iterations" ) );
    return( x );
}



/* ********************************************************************************* */



Cosine_Transformation::Cosine_Transformation()
{
    _fdf.f   = &rcosine_f_func;
    _fdf.df  = &rcosine_df_func;
    _fdf.fdf = &rcosine_fdf_func;
    _solver  = gsl_root_fdfsolver_alloc( gsl_root_fdfsolver_newton );
}


Cosine_Transformation::Cosine_Transformation( const Cosine_Transformation &trans )
{
    _fdf.f   = &rcosine_f_func;
    _fdf.df  = &rcosine_df_func;
    _fdf.fdf = &rcosine_fdf_func;
    _solver  = gsl_root_fdfsolver_alloc( gsl_root_fdfsolver_newton );
}


Cosine_Transformation::~Cosine_Transformation()
{
    gsl_root_fdfsolver_free( _solver );
}


double Cosine_Transformation::rcosine_f_func( double x, void *params )
{
    double *R = (double *)params;

    if( x < -1.0 )
	return( -(*R) );
    else if( x > 1.0 )
	return( 1.0 - (*R) );
    return( 0.5*( 1.0 + x + M_1_PI*sin(x*M_PI) ) - *R );
}


double Cosine_Transformation::rcosine_df_func( double x, void *params )
{
    if( x < -1.0 )
	return( 0.0 );
    else if( x > 1.0 )
	return( 0.0 );
    return( 0.5*( 1.0 + cos(x*M_PI) ) );
}


void Cosine_Transformation::rcosine_fdf_func( double x, void *params, double *f, double *df )
{
    double *R = (double *)params;

    if( x < -1.0 ) {
	*f = -(*R);
	*df = 0.0;
    } else if( x > 1.0 ) {
	*f = 1.0 - (*R);
	*df = 0.0;
    } else {
	*f  = 0.5*( 1.0 + x + M_1_PI*sin(x*M_PI) ) - *R;
	*df = 0.5*( 1.0 + cos(x*M_PI) );
    }
}


Cosine_Transformation *Cosine_Transformation::copy( void ) const
{
    return( new Cosine_Transformation( *this) );
}


double Cosine_Transformation::transform( double R )
{
    int iter = 0;
    double x, x0;

    x = R - 0.5; // Initial guess
    _fdf.params = (void *)&R;
    gsl_root_fdfsolver_set( _solver, &_fdf, x );
    do {
        gsl_root_fdfsolver_iterate( _solver );
        x0 = x;
        x = gsl_root_fdfsolver_root( _solver );
        if( fabs(x - x0) < 1e-6 )
	    break;
	iter++;
    } while( iter < 100 );
    if( iter == 100 )
        throw( Error( ERROR_LOCATION, "too many iterations" ) );
    return( x );
}



/* ********************************************************************************* */



Gamma_Transformation::Gamma_Transformation( double k, double theta )
{
    _k = k;
    _theta = theta;
    _rgamma_const_A = 1.0 / ( gsl_sf_gamma( _k ) * pow( _theta, _k ) );
    _rgamma_const_B = 1.0 / _theta;

    _fdf.f   = &rgamma_f_func;
    _fdf.df  = &rgamma_df_func;
    _fdf.fdf = &rgamma_fdf_func;
    _solver  = gsl_root_fdfsolver_alloc( gsl_root_fdfsolver_newton );
}


Gamma_Transformation::Gamma_Transformation( const Gamma_Transformation &trans )
{
    _k = trans._k;
    _theta = trans._theta;
    _rgamma_const_A = trans._rgamma_const_A;
    _rgamma_const_B = trans._rgamma_const_B;

    _fdf.f   = &rgamma_f_func;
    _fdf.df  = &rgamma_df_func;
    _fdf.fdf = &rgamma_fdf_func;
    _solver  = gsl_root_fdfsolver_alloc( gsl_root_fdfsolver_newton );
}


Gamma_Transformation::~Gamma_Transformation()
{
    gsl_root_fdfsolver_free( _solver );
}


Gamma_Transformation &Gamma_Transformation::operator=( const Gamma_Transformation &trans )
{
    _k = trans._k;
    _theta = trans._theta;
    _rgamma_const_A = trans._rgamma_const_A;
    _rgamma_const_B = trans._rgamma_const_B;

    return( *this );
}


double Gamma_Transformation::rgamma_f_func( double x, void *params )
{
    Gamma_Transformation *trans = (Gamma_Transformation *)params;
    double R = trans->_R;

    if( x < 0.0 )
	return( -R );
    return( gsl_sf_gamma_inc_P( trans->_k, x * trans->_rgamma_const_B ) - R );
}


double Gamma_Transformation::rgamma_df_func( double x, void *params )
{
    Gamma_Transformation *trans = (Gamma_Transformation *)params;

    if( x < 0.0 )
	return( 0.0 );
    return( pow( x, trans->_k-1.0 ) * exp( -x * trans->_rgamma_const_B ) * trans->_rgamma_const_A );
}


void Gamma_Transformation::rgamma_fdf_func( double x, void *params, double *f, double *df )
{
    Gamma_Transformation *trans = (Gamma_Transformation *)params;
    double R = trans->_R;
    
    if( x < 0.0 ) {
	*f = -R;
	*df = 0.0;
    } else {
	*f  = gsl_sf_gamma_inc_P( trans->_k, x * trans->_rgamma_const_B ) - R;
	*df = pow( x, trans->_k-1.0 ) * exp( -x * trans->_rgamma_const_B ) * trans->_rgamma_const_A;
    }
}


Gamma_Transformation *Gamma_Transformation::copy( void ) const
{
    return( new Gamma_Transformation( *this) );
}


double Gamma_Transformation::transform( double R )
{
    int iter = 0;
    double x, x0;

    //std::cout << "Transform( R = " << R << ")\n";

    x = _k*_theta; // Initial guess: mean
    _R = R;
    _fdf.params = (void *)this;
    gsl_root_fdfsolver_set( _solver, &_fdf, x );
    do {
        gsl_root_fdfsolver_iterate( _solver );
        x0 = x;
	//std::cout << "  iter " << iter << ": ";
	//std::cout << "  x = " << x << " ";
        x = gsl_root_fdfsolver_root( _solver );
	//std::cout << ",  err = " << fabs(x - x0) << "\n";
        if( fabs(x - x0) < 1e-6 )
	    break;
	iter++;
    } while( iter < 100 );
    if( iter == 100 )
        throw( Error( ERROR_LOCATION, "too many iterations" ) );
    return( x );
}



/* ********************************************************************************* */



Random::Random( size_t N )
{
    pthread_mutex_init( &_mutex, NULL );
    for( size_t a = 0; a < N; a++ )
	_transformation.push_back( new Uniform_Transformation );
}


Random::~Random()
{
    pthread_mutex_destroy( &_mutex );
    for( size_t a = 0; a < _transformation.size(); a++ )
	delete( _transformation[a] );
}


void Random::set_transformation( size_t i, const Random_Variate_Transformation &trans )
{
    if( i >= _transformation.size() )
	throw( ErrorRange( ERROR_LOCATION, i, _transformation.size() ) );

    delete( _transformation[i] );
    _transformation[i] = trans.copy();
}


/* ********************************************************************************* */



QRandom::QRandom( size_t N )
    : Random(N), _N(N)
{
    _qrng = gsl_qrng_alloc( gsl_qrng_sobol, _N );
}


QRandom::~QRandom()
{
    gsl_qrng_free( _qrng );
}


void QRandom::get( double *x ) const
{
    pthread_mutex_lock( &_mutex );
    gsl_qrng_get( _qrng, x );
    for( size_t a = 0; a < _N; a++ )
	x[a] = _transformation[a]->transform( x[a] );
    pthread_mutex_unlock( &_mutex );
}




/* ********************************************************************************* */



MTRandom::MTRandom( size_t N )
    : Random(N), _N(N)
{
    _rng = gsl_rng_alloc( gsl_rng_mt19937 );
}


MTRandom::~MTRandom()
{
    gsl_rng_free( _rng );
}


void MTRandom::get( double *x ) const
{
    pthread_mutex_lock( &_mutex );
    for( size_t a = 0; a < _N; a++ )
	x[a] = _transformation[a]->transform( gsl_rng_uniform( _rng ) );
    pthread_mutex_unlock( &_mutex );
}

