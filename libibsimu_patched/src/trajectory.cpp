/*! \file trajectory.cpp
 *  \brief Trajectory interpolation solver
 */

/* Copyright (c) 2005-2011,2013 Taneli Kalvas. All rights reserved.
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

#include <cmath>
#include <iostream>
#include "trajectory.hpp"
#include "error.hpp"
#include "polysolver.hpp"


//#define DEBUG_TRAJECTORY 1


void TrajectoryRep1D::construct( double dt, double x1, double v1, 
				 double x2, double v2,
				 trajectory_rep_e force )
{
    double sca = 1.0/(x2-x1);

    switch( force ) {
    case TRAJ_EMPTY:
	_A = (v2+v1)*dt + 2.0*(x1-x2);
	if( fabs(_A*sca) < 1.0e-6 ) {
	    _A = x2-x1-v1*dt;
	    if( fabs(_A*sca) < 1.0e-6 ) {
		_A = x2-x1;
		_B = x1;
		_C = v2-v1;
		_D = v1;
		_rep = TRAJ_LINEAR;
		break;
	    }
	    _B = v1*dt;
	    _C = x1;
	    _D = v2-v1;
	    _E = v1;
	    _rep = TRAJ_QUADRATIC;
	    break;
	}
	_B = 3.0*(x2-x1) - (2.0*v1+v2)*dt;
	_C = v1*dt;
	_D = x1;
	_E = dt;
	_rep = TRAJ_CUBIC;
	break;
    case TRAJ_LINEAR:
	_A = x2-x1;
	_B = x1;
	_C = v2-v1;
	_D = v1;
	_rep = TRAJ_LINEAR;
	break;
    case TRAJ_QUADRATIC:
	_A = x2-x1-v1*dt;
	_B = v1*dt;
	_C = x1;
	_D = v2-v1;
	_E = v1;
	_rep = TRAJ_QUADRATIC;
	break;
    case TRAJ_CUBIC:
	_A = (v2+v1)*dt + 2.0*(x1-x2);
	_B = 3.0*(x2-x1) - (2.0*v1+v2)*dt;
	_C = v1*dt;
	_D = x1;
	_E = dt;
	_rep = TRAJ_CUBIC;
	break;
    };

#ifdef DEBUG_TRAJECTORY
    std::cout << "Constructing trajectory representation\n";
    switch( _rep ) {
    case TRAJ_EMPTY:
	std::cout << "  rep = TRAJ_EMPTY\n";
	break;
    case TRAJ_LINEAR:
	std::cout << "  rep = TRAJ_LINEAR\n";
	std::cout << "  A = " << _A << "\n";
	std::cout << "  B = " << _B << "\n";
	std::cout << "  C = " << _C << "\n";
	std::cout << "  D = " << _D << "\n";
	break;
    case TRAJ_QUADRATIC:
	std::cout << "  rep = TRAJ_QUADRATIC\n";
	std::cout << "  A = " << _A << "\n";
	std::cout << "  B = " << _B << "\n";
	std::cout << "  C = " << _C << "\n";
	std::cout << "  D = " << _D << "\n";
	std::cout << "  E = " << _E << "\n";
	break;
    case TRAJ_CUBIC:
	std::cout << "  rep = TRAJ_CUBIC\n";
	std::cout << "  A = " << _A << "\n";
	std::cout << "  B = " << _B << "\n";
	std::cout << "  C = " << _C << "\n";
	std::cout << "  D = " << _D << "\n";
	std::cout << "  E = " << _E << "\n";
	break;
    };    
#endif
}


TrajectoryRep1D::TrajectoryRep1D( double dt, double x1, double v1, 
				  double x2, double v2,
				  trajectory_rep_e force )
{
    construct( dt, x1, v1, x2, v2, force );
}


void TrajectoryRep1D::coord( double &x, double &v, double K )
{
    switch( _rep ) {
    case TRAJ_EMPTY:
	throw( Error( ERROR_LOCATION, "empty representation" ) );
	break;
    case TRAJ_LINEAR:
	x = _A*K + _B;
	v = _C*K + _D;
	break;
    case TRAJ_QUADRATIC:
	x = (_A*K + _B)*K + _C;
	v = _D*K + _E;
	break;
    case TRAJ_CUBIC:
	x = ((_A*K + _B)*K + _C)*K + _D;
	v = ((3.0*_A*K + 2.0*_B)*K + _C)/_E;
	break;
    };
}


bool TrajectoryRep1D::in( double K, int extrapolate )
{
    if( extrapolate == 0 && K > 0.0 && K <= 1.0 )
	return( true );
    else if( extrapolate < 0 && K > -1.0e-6 && K <= 1.0 )
	return( true );
    else if( K > 0.0 && K <= 1.0+1.0e-6 )
	return( true );
    return( false );
}


int TrajectoryRep1D::solve( double K[3], double x, int extrapolate )
{
#ifdef DEBUG_TRAJECTORY
    std::cout << "solve( x = " << x << " ):\n";
    switch( _rep ) {
    case TRAJ_EMPTY:
	std::cout << "  rep = TRAJ_EMPTY\n";
	break;
    case TRAJ_LINEAR:
	std::cout << "  rep = TRAJ_LINEAR\n";
	break;
    case TRAJ_QUADRATIC:
	std::cout << "  rep = TRAJ_QUADRATIC\n";
	break;
    case TRAJ_CUBIC:
	std::cout << "  rep = TRAJ_CUBIC\n";
	break;
    };
#endif

    switch( _rep ) {
    case TRAJ_EMPTY:
	throw( Error( ERROR_LOCATION, "empty representation" ) );
	break;
    case TRAJ_LINEAR:
    {
	if( _A == 0.0 )
	    return( 0 );
	// Some tests fail here of spotting K equal to one
	// Rounding to 64-bit from 80-bit help
	// This is really not the way to go... temporary
	volatile double KT = (x-_B) / _A;
	K[0] = KT;
	if( in( K[0], extrapolate ) )
	    return( 1 );
	break;
    }
    case TRAJ_QUADRATIC:
    {
	int nroots = solve_quadratic( _A, _B, _C-x, &K[0], &K[1] );
#ifdef DEBUG_TRAJECTORY
	std::cout << "  nroots = " << nroots << "\n";
	for( int a = 0; a < nroots; a++ )
	    std::cout << "  K[" << a << "] = " << K[a] <<"\n";
#endif
	if( nroots == 0 ) {
	    return( 0 );
	} else if( nroots == 1 ) {
	    if( in( K[0], extrapolate ) ) {
		return( 1 );
	    } else {
		return( 0 );
	    }
	} else /* nroots = 2 */ {
	    if( in( K[0], extrapolate ) ) {
		if( in( K[1], extrapolate ) ) {
		    return( 2 );
		} else {
		    return( 1 );
		}
	    } else {
		if( in( K[1], extrapolate ) ) {
		    K[0] = K[1];
		    return( 1 );
		} else {
		    return( 0 );
		}
	    }
	}
	break;
    }
    case TRAJ_CUBIC:
    {
	int nroots = solve_cubic( _A, _B, _C, _D-x, &K[0], &K[1], &K[2] );
	if( nroots == 0 ) {
	    return( 0 );
	} else if( nroots == 1 ) {
	    if( in( K[0], extrapolate ) ) {
		return( 1 );
	    }
	} else if( nroots == 2 ) {
	    if( in( K[0], extrapolate ) ) {
		if( in( K[1], extrapolate ) ) {
		    return( 2 );
		} else {
		    return( 1 );
		}
	    } else {
		if( in( K[1], extrapolate ) ) {
		    K[0] = K[1];
		    return( 1 );
		} else {
		    return( 0 );
		}
	    }
	} else /* nroots = 3 */ {
	    if( in( K[0], extrapolate ) ) {
		if( in( K[1], extrapolate ) ) {
		    if( in( K[2], extrapolate ) ) {
			return( 3 );
		    } else {
			return( 2 );
		    }
		} else {
		    if( in( K[2], extrapolate ) ) {
			K[1] = K[2];
			return( 2 );
		    } else {
			return( 1 );
		    }
		}
	    } else {
		if( in( K[1], extrapolate ) ) {
		    if( in( K[2], extrapolate ) ) {
			K[0] = K[1];
			K[1] = K[2];
			return( 2 );
		    } else {
			K[0] = K[1];
			return( 1 );
		    }
		} else {
		    if( in( K[2], extrapolate ) ) {
			K[0] = K[2];
			return( 1 );
		    } else {
			return( 0 );
		    }
		}
	    }
	}
	break;
    }
    };

    return( 0 );
}


uint32_t TrajectoryRep1D::get_representation_order( void ) const
{
    if( _rep == TRAJ_EMPTY )
	return( 0 );
    else if( _rep == TRAJ_LINEAR )
	return( 1 );
    else if( _rep == TRAJ_QUADRATIC )
	return( 2 );
    else if( _rep == TRAJ_CUBIC )
	return( 3 );
    
    return( (uint32_t)(-1) );
}


void TrajectoryRep1D::debug_print( std::ostream &os ) const
{
    switch( _rep ) {
    case TRAJ_EMPTY:
	os << "rep = EMPTY\n";
	break;
    case TRAJ_LINEAR:
	os << "rep = LINEAR\n";
	os << "x   = A*K + B\n";
	os << "v   = C*K + D\n";
	os << "A   = " << _A << "\n";
	os << "B   = " << _B << "\n";
	os << "C   = " << _C << "\n";
	os << "D   = " << _D << "\n";
	break;
    case TRAJ_QUADRATIC:
	os << "rep = QUADRATIC\n";
	os << "x   = A*K^2 + B*K + C\n";
	os << "v   = D*K + E\n";
	os << "A   = " << _A << "\n";
	os << "B   = " << _B << "\n";
	os << "C   = " << _C << "\n";
	os << "D   = " << _D << "\n";
	os << "E   = " << _E << "\n";
	break;
    case TRAJ_CUBIC:
	os << "rep = CUBIC\n";
	os << "x   = A*K^3 + B*K^2 + C*K + D\n";
	os << "v   = (3*A*K^2 + 2*B*K + C)/E\n";
	os << "A   = " << _A << "\n";
	os << "B   = " << _B << "\n";
	os << "C   = " << _C << "\n";
	os << "D   = " << _D << "\n";
	os << "E   = " << _E << "\n";
	break;
    }
}

