/*! \file trajectorydiagnostics.cpp
 *  \brief Trajectory diagnostics
 */

/* Copyright (c) 2005-2014 Taneli Kalvas. All rights reserved.
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
#include <iomanip>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "compmath.hpp"
#include "trajectorydiagnostics.hpp"
#include "ibsimu.hpp"
#include "error.hpp"


//#define DEBUG_EMITTANCECONV 1


void TrajectoryDiagnosticColumn::mirror( coordinate_axis_e axis, double level )
{
    size_t size = _data.size();
    _data.reserve( 2*size );

    // Handle X-axis
    if( axis == AXIS_X ) {
	if( _diag == DIAG_X ) {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( 2*level-_data[a] );
	} else if( _diag == DIAG_VX || _diag == DIAG_XP ) {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( -_data[a] );
	} else {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( _data[a] );
	}
    }

    // Handle Y-axis
    else if( axis == AXIS_Y || axis == AXIS_R ) {
	if( _diag == DIAG_Y || _diag == DIAG_R ) {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( 2*level-_data[a] );
	} else if( _diag == DIAG_VY || _diag == DIAG_VR || _diag == DIAG_YP || _diag == DIAG_RP || _diag == DIAG_AP ) {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( -_data[a] );
	} else {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( _data[a] );
	}
    }

    // Handle Z-axis
    else if( axis == AXIS_Z ) {
	if( _diag == DIAG_Z ) {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( 2*level-_data[a] );
	} else if( _diag == DIAG_VZ || _diag == DIAG_ZP ) {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( -_data[a] );
	} else {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( _data[a] );
	}
    }
}


void TrajectoryDiagnosticData::export_data( const std::string &filename )
{
    std::ofstream fstr( filename.c_str() );
    if( !fstr.is_open() ) {
	throw( Error( ERROR_LOCATION, 
		      "couldn't open file \'" + filename + "\' for writing" ) );
    }
    ibsimu.message( 1 ) << "Exporting diagnostic data to \'" << filename << "\'\n";

    size_t n = traj_size();
    size_t m = diag_size();

    fstr << "# ";
    for( size_t b = 0; b < m; b++ ) {
	if( b == 0 )
	    fstr << std::setw(11) << trajectory_diagnostic_string_with_unit[_column[b].diagnostic()] << " ";
	else
	    fstr << std::setw(13) << trajectory_diagnostic_string_with_unit[_column[b].diagnostic()] << " ";
    }
    fstr << "\n";

    for( size_t a = 0; a < n; a++ ) {
	for( size_t b = 0; b < m; b++ ) {
	    fstr << std::setw(13) << _column[b][a] << " ";
	}
	fstr << "\n";
    }

    fstr.close();
}


Emittance::Emittance()
    : _Isum(0.0), _xave(0.0), _xpave(0.0), _x2(0.0), _xp2(0.0), _xxp(0.0), 
      _alpha(0.0), _beta(0.0), _gamma(0.0), _epsilon(0.0)
{

}


Emittance::Emittance( const std::vector<double> &x,
		      const std::vector<double> &xp,
		      const std::vector<double> &I )
{
    size_t N = x.size() < xp.size() ? 
	(x.size() < I.size() ? x.size() : I.size()) : 
	(xp.size() < I.size() ? xp.size() : I.size());

    // Calculate averages
    _Isum  = 0.0;
    _xave  = 0.0;
    _xpave = 0.0;
    for( size_t a = 0; a < N; a++ ) {
	_Isum  += I[a];
	_xave  += x[a]*I[a];
	_xpave += xp[a]*I[a];
    }
    _xave  = _xave  / _Isum;
    _xpave = _xpave / _Isum;

    // Calculate expectation values
    _x2  = 0.0;
    _xp2 = 0.0;
    _xxp = 0.0;
    for( size_t a = 0; a < N; a++ ) {
	_x2 += (x[a]-_xave)*(x[a]-_xave)*I[a];
	_xp2 += (xp[a]-_xpave)*(xp[a]-_xpave)*I[a];
	_xxp += (x[a]-_xave)*(xp[a]-_xpave)*I[a];
    }
    _x2  = _x2  / _Isum;
    _xp2 = _xp2 / _Isum;
    _xxp = _xxp / _Isum;

    // Calculate Twiss parameters
    _epsilon = sqrt( _xp2*_x2 - _xxp*_xxp );
    _alpha   = -_xxp/_epsilon;
    _beta    = _x2/_epsilon;
    _gamma   = _xp2/_epsilon;

    // Calculate axes and angle
    _angle = 0.5*atan2( -2.0*_alpha, _beta - _gamma );
    double H = 0.5*(_beta+_gamma);
    _rmajor = sqrt( 0.5*_epsilon ) * ( sqrt(H+1.0)+sqrt(H-1.0) );
    _rminor = sqrt( 0.5*_epsilon ) * ( sqrt(H+1.0)-sqrt(H-1.0) );

#if DEBUG_EMITTANCECONV >= 1
    std::cout << "xave    = " << _xave << "\n";
    std::cout << "xpave   = " << _xpave << "\n";
    std::cout << "x2      = " << _x2 << "\n";
    std::cout << "xp2     = " << _xp2 << "\n";
    std::cout << "xxp     = " << _xxp << "\n";
    std::cout << "epsilon = " << _epsilon << "\n";
    std::cout << "alpha   = " << _alpha << "\n";
    std::cout << "beta    = " << _beta << "\n";
    std::cout << "gamma   = " << _gamma << "\n";
    std::cout << "angle   = " << _angle << "\n";
    std::cout << "H       = " << H << "\n";
    std::cout << "rmajor  = " << _rmajor << "\n";
    std::cout << "rminor  = " << _rminor << "\n";
#endif
}


Emittance::Emittance( const std::vector<double> &x,
		      const std::vector<double> &xp )
{
    size_t N = x.size() > xp.size() ? x.size() : xp.size();

    // Calculate averages
    _xave  = 0.0;
    _xpave = 0.0;
    for( size_t a = 0; a < N; a++ ) {
	_xave  += x[a];
	_xpave += xp[a];
    }
    _xave  = _xave  / (double)N;
    _xpave = _xpave / (double)N;

    // Calculate expectation values
    _x2  = 0.0;
    _xp2 = 0.0;
    _xxp = 0.0;
    for( size_t a = 0; a < N; a++ ) {
	_x2 += (x[a]-_xave)*(x[a]-_xave);
	_xp2 += (xp[a]-_xpave)*(xp[a]-_xpave);
	_xxp += (x[a]-_xave)*(xp[a]-_xpave);
    }
    _x2  = _x2  / (double)N;
    _xp2 = _xp2 / (double)N;
    _xxp = _xxp / (double)N;

    // Calculate Twiss parameters
    _epsilon = sqrt( _xp2*_x2 - _xxp*_xxp );
    _alpha   = -_xxp/_epsilon;
    _beta    = _x2/_epsilon;
    _gamma   = _xp2/_epsilon;

    // Calculate axes and angle
    _angle = 0.5*atan2( -2.0*_alpha, _beta - _gamma );
    double H = 0.5*(_beta+_gamma);
    _rmajor = sqrt( 0.5*_epsilon ) * ( sqrt(H+1.0)+sqrt(H-1.0) );
    _rminor = sqrt( 0.5*_epsilon ) * ( sqrt(H+1.0)-sqrt(H-1.0) );

#if DEBUG_EMITTANCECONV >= 1
    std::cout << "xave    = " << _xave << "\n";
    std::cout << "xpave    = " << _xpave << "\n";
    std::cout << "x2      = " << _x2 << "\n";
    std::cout << "xp2     = " << _xp2 << "\n";
    std::cout << "xxp     = " << _xxp << "\n";
    std::cout << "epsilon = " << _epsilon << "\n";
    std::cout << "alpha   = " << _alpha << "\n";
    std::cout << "beta    = " << _beta << "\n";
    std::cout << "gamma   = " << _gamma << "\n";
    std::cout << "angle   = " << _angle << "\n";
    std::cout << "H       = " << H << "\n";
    std::cout << "rmajor  = " << _rmajor << "\n";
    std::cout << "rminor  = " << _rminor << "\n";
#endif
}


Emittance::Emittance( size_t xsize, size_t xpsize, const double range[4],
		      const std::vector<double> &I )
{
    size_t N = I.size();
    if( xsize*xpsize != N )
	throw( Error( ERROR_LOCATION, 
		      "intensity data size does not match mesh size" ) );

    // Calculate averages
    _Isum  = 0.0;
    _xave  = 0.0;
    _xpave = 0.0;
    for( size_t j = 0; j < xpsize; j++ ) {
	double Qxp = range[1] + ((double)j/(xpsize-1.0))*(range[3]-range[1]);
	for( size_t i = 0; i < xsize; i++ ) {
	    double Qx = range[0] + ((double)i/(xsize-1.0))*(range[2]-range[0]);
	    double QI = I[i+j*xsize];
	    _Isum  += QI;
	    _xave  += Qx*QI;
	    _xpave += Qxp*QI;
	}
    }
    _xave  = _xave  / _Isum;
    _xpave = _xpave / _Isum;

    // Calculate expectation values
    _x2  = 0.0;
    _xp2 = 0.0;
    _xxp = 0.0;
    for( size_t j = 0; j < xpsize; j++ ) {
	double Qxp = range[1] + ((double)j/(xpsize-1.0))*(range[3]-range[1]);
	for( size_t i = 0; i < xsize; i++ ) {
	    double Qx = range[0] + ((double)i/(xsize-1.0))*(range[2]-range[0]);
	    double QI = I[i+j*xsize];
	    _x2  += (Qx-_xave)*(Qx-_xave)*QI;
	    _xp2 += (Qxp-_xpave)*(Qxp-_xpave)*QI;
	    _xxp += (Qx-_xave)*(Qxp-_xpave)*QI;
	}
    }
    _x2  = _x2  / _Isum;
    _xp2 = _xp2 / _Isum;
    _xxp = _xxp / _Isum;

    // Calculate Twiss parameters
    _epsilon = sqrt( _xp2*_x2 - _xxp*_xxp );
    _alpha   = -_xxp/_epsilon;
    _beta    = _x2/_epsilon;
    _gamma   = _xp2/_epsilon;

    // Calculate axes and angle
    _angle = 0.5*atan2( -2.0*_alpha, _beta - _gamma );
    double H = 0.5*(_beta+_gamma);
    _rmajor = sqrt( 0.5*_epsilon ) * ( sqrt(H+1.0)+sqrt(H-1.0) );
    _rminor = sqrt( 0.5*_epsilon ) * ( sqrt(H+1.0)-sqrt(H-1.0) );

#if DEBUG_EMITTANCECONV >= 1
    std::cout << "xave    = " << _xave << "\n";
    std::cout << "xpave   = " << _xpave << "\n";
    std::cout << "x2      = " << _x2 << "\n";
    std::cout << "xp2     = " << _xp2 << "\n";
    std::cout << "xxp     = " << _xxp << "\n";
    std::cout << "epsilon = " << _epsilon << "\n";
    std::cout << "alpha   = " << _alpha << "\n";
    std::cout << "beta    = " << _beta << "\n";
    std::cout << "gamma   = " << _gamma << "\n";
    std::cout << "angle   = " << _angle << "\n";
    std::cout << "H       = " << H << "\n";
    std::cout << "rmajor  = " << _rmajor << "\n";
    std::cout << "rminor  = " << _rminor << "\n";
#endif
}


void Emittance::debug_print( std::ostream &os ) const
{
    os << "xave    = " << _xave << "\n";
    os << "xpave   = " << _xpave << "\n";
    os << "x2      = " << _x2 << "\n";
    os << "xp2     = " << _xp2 << "\n";
    os << "xxp     = " << _xxp << "\n";
    os << "epsilon = " << _epsilon << "\n";
    os << "alpha   = " << _alpha << "\n";
    os << "beta    = " << _beta << "\n";
    os << "gamma   = " << _gamma << "\n";
    os << "angle   = " << _angle << "\n";
    os << "rmajor  = " << _rmajor << "\n";
    os << "rminor  = " << _rminor << "\n";
}


int min( int n1, int n2, int n3, int n4 )
{
    int n;

    if( n1 < n2 )
	n = n1;
    else
	n = n2;

    if( n3 < n )
	n = n3;

    if( n4 < n )
	n = n4;
    
    return( n );
}


EmittanceConv::EmittanceConv( uint32_t n, uint32_t m,
			      const std::vector<double> &r,
			      const std::vector<double> &rp,
			      const std::vector<double> &ap,
			      const std::vector<double> &I,
			      uint32_t rotn,
			      double xmin, double xpmin, 
			      double xmax, double xpmax )
{
    int N = min( r.size(), rp.size(), ap.size(), I.size() );

    // Find ranges
    double range[4] = { xmin, xpmin, xmax, xpmax };

    // xmax = rmax, xmin = -rmax
    double rmax = 0.0;
    for( int a = 0; a < N; a++ ) {
	if( r[a] > rmax )
	    rmax = r[a];
    }
    if( comp_isnan(range[0]) )
	range[0] = -rmax;
    if( comp_isnan(range[2]) )
	range[2] = rmax;

    // xpmax
    double _xpmax = 0.0;
    for( int a = 0; a < N; a++ ) {
	double xpm = sqrt( rp[a]*rp[a] + ap[a]*ap[a] );
	if( xpm > _xpmax )
	    _xpmax = xpm;
    }
    if( comp_isnan(range[1]) )
	range[1] = -_xpmax;
    if( comp_isnan(range[3]) )
	range[3] = _xpmax;

    // Make grid
    _grid = new Histogram2D( n, m, range );

#if DEBUG_EMITTANCECONV >= 2
    std::cout << "Building emittance conversion\n";
    std::cout << "  N     = " << N << "\n";
    std::cout << "  n     = " << n << "\n";
    std::cout << "  m     = " << m << "\n";
    std::cout << "  range = {" 
	      << range[0] << ", "
	      << range[1] << ", "
	      << range[2] << ", "
	      << range[3] << "}\n";
#endif

    // For each particle
    _Isum = 0.0;
    _x2 = _xp2 = _xxp = 0.0;
    _xave = _xpave = 0.0;
    for( int a = 0; a < N; a++ ) {

#if DEBUG_EMITTANCECONV >= 2
	std::cout << "  Particle " << a << "\n";
	std::cout << "    r  = " << r[a] << "\n";
	std::cout << "    rp = " << rp[a] << "\n";
	std::cout << "    ap = " << ap[a] << "\n";
	std::cout << "    I  = " << I[a] << "\n\n";
#endif

	// Skip zero r
	if( r[a] == 0.0 )
	    continue;

	// Rotate around xy-plane
	_Isum += I[a];
	double x2 = 0.0; // Avoid numerical noise
	double xp2 = 0.0;
	double xxp = 0.0;
	double dI = I[a]/rotn;
	for( uint32_t b = 0; b < rotn; b++ ) {

	    double rnd  = ((double)rand())/RAND_MAX;
	    double ang  = 2.0*M_PI*(b+rnd)/rotn;
	    double sint = sin( ang );
	    double cost = cos( ang );

	    double x  = r[a]*sint;
	    double xp = rp[a]*sint + ap[a]*cost;

	    x2 += x*x*dI;
	    xp2 += xp*xp*dI;
	    xxp += x*xp*dI;

	    _grid->accumulate_closest( x, xp, dI );

	    /*
	    int i = (int)floor( (x -range[0]) / _grid->nstep() + 0.5 );
	    int j = (int)floor( (xp-range[1]) / _grid->mstep() + 0.5 );

	    if( i >= 0 && i < n && j >= 0 && j < m )
		_grid->accumulate( i, j, dI );
	    */
	}

	_x2 += x2;
	_xp2 += xp2;
	_xxp += xxp;
    }

    _x2  = _x2  / _Isum;
    _xp2 = _xp2 / _Isum;
    _xxp = _xxp / _Isum;

    // Calculate Twiss parameters
    _epsilon = sqrt( _xp2*_x2 - _xxp*_xxp );
    _alpha   = -_xxp/_epsilon;
    _beta    = _x2/_epsilon;
    _gamma   = _xp2/_epsilon;

    // Calculate axes and angle
    _angle = 0.5*atan2( -2.0*_alpha, _beta - _gamma );
    double H = 0.5*(_beta+_gamma);
    _rmajor = sqrt( 0.5*_epsilon ) * ( sqrt(H+1.0)+sqrt(H-1.0) );
    _rminor = sqrt( 0.5*_epsilon ) * ( sqrt(H+1.0)-sqrt(H-1.0) );

#if DEBUG_EMITTANCECONV >= 1
    std::cout << "xave    = " << _xave << "\n";
    std::cout << "xpave    = " << _xpave << "\n";
    std::cout << "x2      = " << _x2 << "\n";
    std::cout << "xp2     = " << _xp2 << "\n";
    std::cout << "xxp     = " << _xxp << "\n";
    std::cout << "epsilon = " << _epsilon << "\n";
    std::cout << "alpha   = " << _alpha << "\n";
    std::cout << "beta    = " << _beta << "\n";
    std::cout << "gamma   = " << _gamma << "\n";
    std::cout << "angle   = " << _angle << "\n";
    std::cout << "H       = " << H << "\n";
    std::cout << "rmajor  = " << _rmajor << "\n";
    std::cout << "rminor  = " << _rminor << "\n";
#endif
}


/*	
EmittanceConv::EmittanceConv( int n, int m,
			      const std::vector<double> &r,
			      const std::vector<double> &rp,
			      const std::vector<double> &ap,
			      const std::vector<double> &I )
{
    int N = min( r.size(), rp.size(), ap.size(), I.size() );

    // Find ranges
    double range[4] = { 0.0, 0.0, 0.0, 0.0 };

    // xmax = rmax, xmin = -rmax
    double rmax = 0.0;
    for( int a = 0; a < N; a++ ) {
	if( r[a] > rmax )
	    rmax = r[a];
    }
    range[0] = -rmax;
    range[2] = rmax;

    // xpmax
    double xpmax = 0.0;
    for( int a = 0; a < N; a++ ) {
	double t = ap[a]/rp[a];
	double xpm = rp[a]*sqrt( 1.0 + t*t );
	if( xpm > xpmax )
	    xpmax = xpm;
    }
    range[1] = -xpmax;
    range[3] = xpmax;

    // Make grid
    _grid = new Histogram2D( n, m, range );
    double dx = _grid->nstep();
    double dxp = _grid->mstep();

#ifdef DEBUG_EMITTANCECONV 
    std::cout << "Building emittance conversion\n";
    std::cout << "  N     = " << N << "\n";
    std::cout << "  n     = " << n << "\n";
    std::cout << "  m     = " << m << "\n";
    std::cout << "  dx    = " << dx << "\n";
    std::cout << "  dxp   = " << dxp << "\n";
    std::cout << "  range = {" 
	      << range[0] << ", "
	      << range[1] << ", "
	      << range[2] << ", "
	      << range[3] << "}\n";
#endif

    // For each particle
    for( int a = 0; a < N; a++ ) {

#ifdef DEBUG_EMITTANCECONV 
	std::cout << "  Particle " << a << "\n";
	std::cout << "    r  = " << r[a] << "\n";
	std::cout << "    rp = " << rp[a] << "\n";
	std::cout << "    ap = " << ap[a] << "\n";
	std::cout << "    I  = " << I[a] << "\n\n";
#endif

	// Skip zero r
	if( r[a] == 0.0 )
	    continue;

	    // Draw a line on grid from -r to r
	int i = (int)floor( (-r[a]-range[0]) / dx + 0.5 );

	bool done = false;
	while( !done ) {

	    // Current density to this grid point
	    double xi = _grid->icoord(i);
	    double x1 = xi-0.5*dx;
	    double x2 = xi+0.5*dx;

#ifdef DEBUG_EMITTANCECONV 
	    std::cout << "    i  = " << i << "\n";
	    std::cout << "    xi = " << xi << "\n";
	    std::cout << "    x1 = " << x1 << "\n";
	    std::cout << "    x2 = " << x2 << "\n";
#endif

	    // Check validity of points
	    if( x1 < -r[a] ) {
		x1 = -r[a];
#ifdef DEBUG_EMITTANCECONV 
		std::cout << "    x1 = " << x1 << " (limit)\n";
#endif
	    }
	    if( x2 >= r[a] ) {
		x2 = r[a];
		done = true;
#ifdef DEBUG_EMITTANCECONV 
		std::cout << "    x2 = " << x2 << " (limit)\n";
#endif
	    }
	    if( xi < -r[a] ) {
		xi = -r[a];
#ifdef DEBUG_EMITTANCECONV 
		std::cout << "    xi = " << xi << " (limit)\n";
#endif
	    }
	    if( xi > r[a] ) {
		xi = r[a];
#ifdef DEBUG_EMITTANCECONV 
		std::cout << "    xi = " << xi << " (limit)\n";
#endif
	    }

	    double dI = I[a]*( asin(x2/r[a]) - asin(x1/r[a]) ) / ( M_PI*r[a] );

	    // Deposit to closest xp(j)
	    double t = xi/r[a];
	    double xp = rp[a]*xi/r[a] - ap[a]*sqrt( 1.0 - t*t );
	    int j = (int)floor( (xp-range[1]) / _grid->mstep() + 0.5 );

#ifdef DEBUG_EMITTANCECONV 
	    std::cout << "    dI = " << dI << "\n";
	    std::cout << "    t  = " << t << "\n";
	    std::cout << "    xp = " << xp << "\n";
	    std::cout << "    j  = " << j << "\n\n";
#endif

	    if( i >= 0 && i < n && j >= 0 && j < m )
		_grid->accumulate( i, j, dI );

	    i++;
	}
    }

    // Build Emittance statistics from grid

    // Calculate averages
    for( int j = 0; j < m; j++ ) {
	for( int i = 0; i < n; i++ ) {
	    double I = (*_grid)(i,j);
	    _Isum  += I;
	    _xave  += _grid->icoord(i)*I;
	    _xpave += _grid->jcoord(j)*I;
	}
    }
    _xave  = _xave  / _Isum;
    _xpave = _xpave / _Isum;

    // Calculate expectation values
    for( int j = 0; j < m; j++ ) {
	for( int i = 0; i < n; i++ ) {
	    double I = (*_grid)(i,j);
	    double x = _grid->icoord(i);
	    double xp = _grid->jcoord(j);
	    _x2 += (x-_xave)*(x-_xave)*I;
	    _xp2 += (xp-_xpave)*(xp-_xpave)*I;
	    _xxp += (x-_xave)*(xp-_xpave)*I;
	}
    }    
    _x2  = _x2  / _Isum;
    _xp2 = _xp2 / _Isum;
    _xxp = _xxp / _Isum;

    // Calculate Twiss parameters
    _epsilon = sqrt( _xp2*_x2 - _xxp*_xxp );
    _alpha   = -_xxp/_epsilon;
    _beta    = _x2/_epsilon;
    _gamma   = _xp2/_epsilon;

    // Calculate axes and angle
    _angle = 0.5*atan2( (-2.0*_alpha) , (_beta - _gamma) );
    double H = 0.5*(_beta+_gamma);
    _rmajor = sqrt( 0.5*_epsilon ) * ( sqrt(H+1.0)+sqrt(H-1.0) );
    _rminor = sqrt( 0.5*_epsilon ) * ( sqrt(H+1.0)-sqrt(H-1.0) );
}
*/


EmittanceConv::~EmittanceConv()
{
    if( _grid )
	delete _grid;
}

