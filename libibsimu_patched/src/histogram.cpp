/*! \file histogram.cpp
 *  \brief %Histogram data handling for 1D and 2D
 */

/* Copyright (c) 2005-2011,2014 Taneli Kalvas. All rights reserved.
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

#include <cstring>
#include <cmath>
#include <limits>
#include <iostream>
#include "histogram.hpp"
#include "error.hpp"


Histogram1D::Histogram1D( uint32_t n, const double range[2] )
    : _n(n), _data(n,0.0)
{
    if( _n < 4 )
	throw( Error( ERROR_LOCATION, "too small histogram size" ) );

    _range[0] = range[0];
    _range[1] = range[1];

    _step = (_range[1]-_range[0]) / (_n-1.0);
}


Histogram1D::Histogram1D( uint32_t n, 
			  const std::vector<double> &xdata,
			  histogram_accumulation_e type )
    : _n(n), _data(n,0.0)
{
    if( _n < 4 )
	throw( Error( ERROR_LOCATION, "too small histogram size" ) );

    uint32_t N = xdata.size();
    if( N == 0 ) {
	// No input data -> return empty histogram
	_range[0] = -1.0;
	_range[1] = +1.0;
	_step = (_range[1]-_range[0]) / (_n-1.0);
	return;
    }

    // Find range limits so that the furthest points will receive no
    // contribution from data.
    double bbox[2];
    bbox[0] = std::numeric_limits<double>::infinity();
    bbox[1] = -std::numeric_limits<double>::infinity();
    for( uint32_t a = 0; a < N; a++ ) {
	if( xdata[a] < bbox[0] )
	    bbox[0] = xdata[a];
	if( xdata[a] >bbox[1] )
	    bbox[1] = xdata[a];
    }

    if( bbox[0] != bbox[1] ) {
	// Widen range by one bin width in each direction
	_range[0] = bbox[0] - (bbox[1]-bbox[0])/(_n-3.0);
	_range[1] = bbox[1] + (bbox[1]-bbox[0])/(_n-3.0);
    } else if( bbox[0] != 0.0 ) {
	// Widen range by 10 % of value in each direction
	_range[0] = bbox[0] * 0.9;
	_range[1] = bbox[1] * 1.1;
    } else {
	// Widen range by 1.0e-6
	_range[0] = bbox[0] - 1.0e-6;
	_range[1] = bbox[1] + 1.0e-6;
    }

    // Calculate step size
    _step = (_range[1]-_range[0]) / (_n-1.0);

    // Add data
    if( type == HISTOGRAM_ACCUMULATION_CLOSEST ) {
	for( uint32_t a = 0; a < N; a++ )
	    accumulate_closest( xdata[a], 1.0 );
    } else {
	for( uint32_t a = 0; a < N; a++ )
	    accumulate_linear( xdata[a], 1.0 );
    }
}


Histogram1D::Histogram1D( uint32_t n, 
			  const std::vector<double> &xdata, 
			  const std::vector<double> &wdata,
			  histogram_accumulation_e type )
    : _n(n), _data(n,0.0)
{
    if( _n < 4 )
	throw( Error( ERROR_LOCATION, "too small histogram size" ) );

    uint32_t N = xdata.size() < wdata.size() ? xdata.size() : wdata.size();
    if( N == 0 ) {
	// No input data -> return empty histogram
	_range[0] = -1.0;
	_range[1] = +1.0;
	_step = (_range[1]-_range[0]) / (_n-1.0);
	return;
    }

    // Find range limits so that the furthest points will receive no
    // contribution from data.
    double bbox[2];
    bbox[0] = std::numeric_limits<double>::infinity();
    bbox[1] = -std::numeric_limits<double>::infinity();
    for( uint32_t a = 0; a < N; a++ ) {
	if( xdata[a] < bbox[0] )
	    bbox[0] = xdata[a];
	if( xdata[a] >bbox[1] )
	    bbox[1] = xdata[a];
    }

    if( bbox[0] != bbox[1] ) {
	// Widen range by one bin width in each direction
	_range[0] = bbox[0] - (bbox[1]-bbox[0])/(_n-3.0);
	_range[1] = bbox[1] + (bbox[1]-bbox[0])/(_n-3.0);
    } else if( bbox[0] != 0.0 ) {
	// Widen range by 10 % of value in each direction
	_range[0] = bbox[0] * 0.9;
	_range[1] = bbox[1] * 1.1;
    } else {
	// Widen range by 1.0e-6
	_range[0] = bbox[0] - 1.0e-6;
	_range[1] = bbox[1] + 1.0e-6;
    }

    // Calculate step size
    _step = (_range[1]-_range[0]) / (_n-1.0);

    // Add data
    if( type == HISTOGRAM_ACCUMULATION_CLOSEST ) {
	for( uint32_t a = 0; a < N; a++ )
	    accumulate_closest( xdata[a], wdata[a] );
    } else {
	for( uint32_t a = 0; a < N; a++ )
	    accumulate_linear( xdata[a], wdata[a] );
    }
}


Histogram1D::~Histogram1D()
{

}


void Histogram1D::get_bin_range( double &min, double &max ) const
{
    min = std::numeric_limits<double>::infinity();
    max = -std::numeric_limits<double>::infinity();

    for( uint32_t a = 0; a < _n; a++ ) {
	if( _data[a] < min ) 
	    min = _data[a];
	if( _data[a] > max ) 
	    max = _data[a];
    }
}


void Histogram1D::accumulate_closest( double x, double weight )
{
    int32_t i = (int)floor( (x-_range[0]) / _step + 0.5 );

    if( i < (int32_t)_n && i >= 0 )
	_data[i] += weight;
}


void Histogram1D::accumulate_closest( const std::vector<double> &xdata )
{
    uint32_t N = xdata.size();
    for( uint32_t a = 0; a < N; a++ )
	accumulate_closest( xdata[a], 1.0 );
}


void Histogram1D::accumulate_closest( const std::vector<double> &xdata,
				      const std::vector<double> &wdata )
{
    uint32_t N = xdata.size() < wdata.size() ? xdata.size() : wdata.size();
    for( uint32_t a = 0; a < N; a++ )
	accumulate_closest( xdata[a], wdata[a] );
}


void Histogram1D::accumulate_linear( const std::vector<double> &xdata )
{
    uint32_t N = xdata.size();
    for( uint32_t a = 0; a < N; a++ )
	accumulate_linear( xdata[a], 1.0 );
}


void Histogram1D::accumulate_linear( const std::vector<double> &xdata,
				     const std::vector<double> &wdata )
{
    uint32_t N = xdata.size() < wdata.size() ? xdata.size() : wdata.size();
    for( uint32_t a = 0; a < N; a++ )
	accumulate_linear( xdata[a], wdata[a] );
}


void Histogram1D::accumulate_linear( double x, double weight )
{
    //std::cout << "(x,y) = (" << x << "," << y << ")\n";
    int32_t i = (int)floor( (x-_range[0]) / _step );
    //std::cout << "(i,j) = (" << i << "," << j << ")\n";
    double t = (x - _range[0])/_step - i;
    //std::cout << "(t,u) = (" << t << "," << u << ")\n";

    if( i < (int32_t)_n-1 && i >= 0 ) {
	// No checks necessary 
	_data[i  ] += weight*(1.0-t);
	_data[i+1] += weight*t;
	return;
    }
	
    // Thorough checks needed
    if( i >= 0 && i < (int32_t)_n )
	_data[i]   += weight*(1.0-t);
    if( i+1 >= 0 && i+1 < (int32_t)_n )
	_data[i+1] += weight*t;
}


double Histogram1D::coord( uint32_t i ) const 
{
    return( _range[0] + i*(_range[1]-_range[0]) / (_n-1.0) ); 
}


void Histogram1D::get_range( double range[2] ) const 
{ 
    range[0] =  _range[0];
    range[1] =  _range[1];
}


void Histogram1D::convert_to_density( void )
{
    // Scale with inverse of area of histogram bin
    double scale = 1.0/_step;

    *this *= scale;
}


const Histogram1D &Histogram1D::operator*=( double x )
{
    // Go through the histogram
    for( uint32_t i = 0; i < _n; i++ )
	_data[i] *= x;

    return( *this );
}




/* ****************************************************************************
 *
 */


Histogram2D::Histogram2D( uint32_t n, uint32_t m, const double range[4] )
    : _n(n), _m(m), _data(n*m,0.0)
{
    if( _n < 4 || _m < 4 )
	throw( Error( ERROR_LOCATION, "too small histogram size" ) );

    _range[0] = range[0];
    _range[1] = range[1];
    _range[2] = range[2];
    _range[3] = range[3];

    _nstep = (_range[2]-_range[0]) / (_n-1.0);
    _mstep = (_range[3]-_range[1]) / (_m-1.0);
}


Histogram2D::Histogram2D( uint32_t n, uint32_t m, 
			  const std::vector<double> &xdata,
			  const std::vector<double> &ydata,
			  histogram_accumulation_e type )
    : _n(n), _m(m), _data(n*m,0.0)
{
    if( _n < 4 || _m < 4 )
	throw( Error( ERROR_LOCATION, "too small histogram size" ) );

    uint32_t N = xdata.size() < ydata.size() ? xdata.size() : ydata.size();
    if( N == 0 ) {
	// No input data -> return empty histogram
	_range[0] = _range[1] = -1.0;
	_range[2] = _range[3] = +1.0;
	_nstep = (_range[2]-_range[0]) / (_n-1.0);
	_mstep = (_range[3]-_range[1]) / (_m-1.0);
	return;
    }

    // Find range limits so that the furthest points will receive no
    // contribution from data.
    double bbox[4];
    bbox[0] = std::numeric_limits<double>::infinity();
    bbox[1] = std::numeric_limits<double>::infinity();
    bbox[2] = -std::numeric_limits<double>::infinity();
    bbox[3] = -std::numeric_limits<double>::infinity();
    for( uint32_t a = 0; a < N; a++ ) {
	if( xdata[a] < bbox[0] )
	    bbox[0] = xdata[a];
	if( xdata[a] >bbox[2] )
	    bbox[2] = xdata[a];

	if( ydata[a] < bbox[1] )
	    bbox[1] = ydata[a];
	if( ydata[a] >bbox[3] )
	    bbox[3] = ydata[a];
    }

    if( bbox[0] != bbox[2] ) {
	// Widen range by one bin width in each direction
	_range[0] = bbox[0] - (bbox[2]-bbox[0])/(_n-3.0);
	_range[2] = bbox[2] + (bbox[2]-bbox[0])/(_n-3.0);
    } else if( bbox[0] != 0.0 ) {
	// Widen range by 10 % of value in each direction
	_range[0] = bbox[0] * 0.9;
	_range[2] = bbox[2] * 1.1;
    } else {
	// Widen range by 1.0e-6
	_range[0] = bbox[0] - 1.0e-6;
	_range[2] = bbox[2] + 1.0e-6;
    }

    if( bbox[1] != bbox[3] ) {
	// Widen range by one bin width in each direction
	_range[1] = bbox[1] - (bbox[3]-bbox[1])/(_m-3.0);
	_range[3] = bbox[3] + (bbox[3]-bbox[1])/(_m-3.0);
    } else if( bbox[1] != 0.0 ) {
	// Widen range by 10 % of value in each direction
	_range[1] = bbox[1] * 0.9;
	_range[3] = bbox[3] * 1.1;
    } else {
	// Widen range by 1.0e-6
	_range[1] = bbox[1] - 1.0e-6;
	_range[3] = bbox[3] + 1.0e-6;
    }

    // Calculate step size
    _nstep = (_range[2]-_range[0]) / (_n-1.0);
    _mstep = (_range[3]-_range[1]) / (_m-1.0);

    // Add data
    if( type == HISTOGRAM_ACCUMULATION_CLOSEST ) {
	for( uint32_t a = 0; a < N; a++ )
	    accumulate_closest( xdata[a], ydata[a], 1.0 );
    } else {
	for( uint32_t a = 0; a < N; a++ )
	    accumulate_linear( xdata[a], ydata[a], 1.0 );
    }
}


Histogram2D::Histogram2D( uint32_t n, uint32_t m, 
			  const std::vector<double> &xdata,
			  const std::vector<double> &ydata,
			  const std::vector<double> &wdata,
			  histogram_accumulation_e type )
    : _n(n), _m(m), _data(n*m,0.0)
{
    if( _n < 4 || _m < 4 )
	throw( Error( ERROR_LOCATION, "too small histogram size" ) );

    uint32_t N = xdata.size() < ydata.size() ? 
	(xdata.size() < wdata.size() ? xdata.size() : wdata.size()) :
	(ydata.size() < wdata.size() ? ydata.size() : wdata.size());
    if( N == 0 ) {
	// No input data -> return empty histogram
	_range[0] = _range[1] = -1.0;
	_range[2] = _range[3] = +1.0;
	_nstep = (_range[2]-_range[0]) / (_n-1.0);
	_mstep = (_range[3]-_range[1]) / (_m-1.0);
	return;
    }

    // Find range limits so that the furthest points will receive no
    // contribution from data.
    double bbox[4];
    bbox[0] = std::numeric_limits<double>::infinity();
    bbox[1] = std::numeric_limits<double>::infinity();
    bbox[2] = -std::numeric_limits<double>::infinity();
    bbox[3] = -std::numeric_limits<double>::infinity();
    for( uint32_t a = 0; a < N; a++ ) {
	if( xdata[a] < bbox[0] )
	    bbox[0] = xdata[a];
	if( xdata[a] >bbox[2] )
	    bbox[2] = xdata[a];

	if( ydata[a] < bbox[1] )
	    bbox[1] = ydata[a];
	if( ydata[a] >bbox[3] )
	    bbox[3] = ydata[a];
    }

    if( bbox[0] != bbox[2] ) {
	// Widen range by one bin width in each direction
	_range[0] = bbox[0] - (bbox[2]-bbox[0])/(_n-3.0);
	_range[2] = bbox[2] + (bbox[2]-bbox[0])/(_n-3.0);
    } else if( bbox[0] != 0.0 ) {
	// Widen range by 10 % of value in each direction
	_range[0] = bbox[0] * 0.9;
	_range[2] = bbox[2] * 1.1;
    } else {
	// Widen range by 1.0e-6
	_range[0] = bbox[0] - 1.0e-6;
	_range[2] = bbox[2] + 1.0e-6;
    }

    if( bbox[1] != bbox[3] ) {
	// Widen range by one bin width in each direction
	_range[1] = bbox[1] - (bbox[3]-bbox[1])/(_m-3.0);
	_range[3] = bbox[3] + (bbox[3]-bbox[1])/(_m-3.0);
    } else if( bbox[1] != 0.0 ) {
	// Widen range by 10 % of value in each direction
	_range[1] = bbox[1] * 0.9;
	_range[3] = bbox[3] * 1.1;
    } else {
	// Widen range by 1.0e-6
	_range[1] = bbox[1] - 1.0e-6;
	_range[3] = bbox[3] + 1.0e-6;
    }

    // Calculate step size
    _nstep = (_range[2]-_range[0]) / (_n-1.0);
    _mstep = (_range[3]-_range[1]) / (_m-1.0);

    // Add data
    if( type == HISTOGRAM_ACCUMULATION_CLOSEST ) {
	for( uint32_t a = 0; a < N; a++ )
	    accumulate_closest( xdata[a], ydata[a], wdata[a] );
    } else {
	for( uint32_t a = 0; a < N; a++ )
	    accumulate_linear( xdata[a], ydata[a], wdata[a] );
    }
}


Histogram2D::~Histogram2D()
{

}


double Histogram2D::icoord( uint32_t i ) const 
{
    return( _range[0] + i*(_range[2]-_range[0]) / (_n-1.0) ); 
}


double Histogram2D::jcoord( uint32_t j ) const 
{ 
    return( _range[1] + j*(_range[3]-_range[1]) / (_m-1.0) ); 
}


void Histogram2D::get_range( double range[4] ) const 
{ 
    range[0] =  _range[0];
    range[1] =  _range[1];
    range[2] =  _range[2];
    range[3] =  _range[3];
}


void Histogram2D::convert_to_density( void )
{
    // Scale with inverse of area of histogram bin
    double scale = 1.0/(_nstep*_mstep);

    *this *= scale;
}


const Histogram2D &Histogram2D::operator*=( double x )
{
    // Go through the histogram
    uint32_t size = _n*_m;
    for( uint32_t i = 0; i < size; i++ )
	_data[i] *= x;

    return( *this );
}


void Histogram2D::get_bin_range( double &min, double &max ) const
{
    min = std::numeric_limits<double>::infinity();
    max = -std::numeric_limits<double>::infinity();

    uint32_t size = _n*_m;
    for( uint32_t a = 0; a < size; a++ ) {
	if( _data[a] < min )
	    min = _data[a];
	if( _data[a] > max )
	    max = _data[a];
    }
}    


void Histogram2D::accumulate_closest( double x, double y, double weight )
{
    int32_t i = (int)floor( (x-_range[0]) / _nstep + 0.5 );
    int32_t j = (int)floor( (y-_range[1]) / _mstep + 0.5 );

    if( i < (int32_t)_n && i >= 0 &&
	j < (int32_t)_m && j >= 0 )
	_data[i+j*_n] += weight;
}


void Histogram2D::accumulate_linear( double x, double y, double weight )
{
    //std::cout << "(x,y) = (" << x << "," << y << ")\n";
    int32_t i = (int32_t)floor( (x-_range[0]) / _nstep );
    int32_t j = (int32_t)floor( (y-_range[1]) / _mstep );
    //std::cout << "(i,j) = (" << i << "," << j << ")\n";
    double t = (x - _range[0])/_nstep - i;
    double u = (y - _range[1])/_mstep - j;
    //std::cout << "(t,u) = (" << t << "," << u << ")\n";

    if( i < (int32_t)_n-1 && j < (int32_t)_m-1 && i >= 0 && j >= 0 ) {
	// No checks necessary 
	_data[i  +    j*_n] += weight*(1.0-t)*(1.0-u);
	_data[i+1+    j*_n] += weight*t*(1.0-u);
	_data[i  +(j+1)*_n] += weight*(1.0-t)*u;
	_data[i+1+(j+1)*_n] += weight*t*u;
	return;
    }
	
    // Thorough checks needed
    if( i >= 0 && j >= 0 && i < (int32_t)_n && j < (int32_t)_m )
	_data[i  +    j*_n] += weight*(1.0-t)*(1.0-u);
    if( i+1 >= 0 && j >= 0 && i+1 < (int32_t)_n && j < (int32_t)_m )
	_data[i+1+    j*_n] += weight*t*(1.0-u);
    if( i >= 0 && j+1 >= 0 && i < (int32_t)_n && j+1 < (int32_t)_m )
	_data[i  +(j+1)*_n] += weight*(1.0-t)*u;
    if( i+1 >= 0 && j+1 >= 0 && i+1 <(int32_t) _n && j+1 < (int32_t)_m )
	_data[i+1+(j+1)*_n] += weight*t*u;
}


