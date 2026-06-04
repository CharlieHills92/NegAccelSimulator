/*! \file axisymmetricvectorfield.cpp
 *  \brief Axisymmetric magnetic field
 */

/* Copyright (c) 2012-2014 Taneli Kalvas. All rights reserved.
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


#include "axisymmetricvectorfield.hpp"


AxisymmetricVectorField::AxisymmetricVectorField( geom_mode_e geom_mode, 
						  const std::vector<double> &z,
						  const std::vector<double> &Bz )
{
    if( geom_mode != MODE_3D && geom_mode != MODE_CYL )
	throw( Error( ERROR_LOCATION, "unsupported geometry mode" ) );
    _geom_mode = geom_mode;
    if( z.size() != Bz.size() )
	throw( Error( ERROR_LOCATION, "non-matching Bz and z arrays" ) );
    _spline = gsl_spline_alloc( gsl_interp_cspline, z.size() );
    _accel = gsl_interp_accel_alloc();
    gsl_spline_init( _spline, &z[0], &Bz[0], z.size() );
}


AxisymmetricVectorField::AxisymmetricVectorField( geom_mode_e geom_mode, 
						  double origo, double h,
						  const std::vector<double> &Bz,
						  uint32_t order )
{
    if( Bz.size() < 7 )
	throw( Error( ERROR_LOCATION, "not enough nodes" ) );
    if( geom_mode != MODE_3D && geom_mode != MODE_CYL )
	throw( Error( ERROR_LOCATION, "unsupported geometry mode" ) );
    _geom_mode = geom_mode;
    _spline = NULL;
    _accel = NULL;
    _origo = origo;
    _h = h;
    _Bz = Bz;
    _order = order;
}


AxisymmetricVectorField::AxisymmetricVectorField( const AxisymmetricVectorField &f )
{
    _geom_mode = f._geom_mode;
    if( f._spline ) {
	_spline = gsl_spline_alloc( gsl_interp_cspline, f._spline->size );
	_accel = gsl_interp_accel_alloc();
	gsl_spline_init( _spline, f._spline->x, f._spline->y, f._spline->size );
    } else {
	_origo = f._origo;
	_h = f._h;
	_Bz = f._Bz;
    }
}


AxisymmetricVectorField &AxisymmetricVectorField::operator=( const AxisymmetricVectorField &f )
{
    if( _spline ) {
	gsl_spline_free( _spline );
	gsl_interp_accel_free( _accel );
    } else {
	_Bz.clear();
    }

    _geom_mode = f._geom_mode;
    if( f._spline ) {
	_spline = gsl_spline_alloc( gsl_interp_cspline, f._spline->size );
	_accel = gsl_interp_accel_alloc();
	gsl_spline_init( _spline, f._spline->x, f._spline->y, f._spline->size );
    } else {
	_origo = f._origo;
	_h = f._h;
	_Bz = f._Bz;
    }
    return( *this );
}


AxisymmetricVectorField::~AxisymmetricVectorField()
{
    gsl_spline_free( _spline );
    gsl_interp_accel_free( _accel );
}


const Vec3D AxisymmetricVectorField::eval_spline( double z, double r ) const
{
    if( z < _spline->x[0] || z > _spline->x[_spline->size-1] )
	return( 0.0 );
    double Bz0  = gsl_spline_eval( _spline, z, _accel );
    double Bzp  = gsl_spline_eval_deriv( _spline, z, _accel );
    double Bzpp = gsl_spline_eval_deriv2( _spline, z, _accel );
    
    double Bz = Bz0 - r*r*Bzpp/4.0;
    double Br = -r*Bzp/2.0;
    return( Vec3D( Bz, Br, 0.0 ) );
}


void AxisymmetricVectorField::get_fdm_derivatives( double Bder[7], double z ) const
{
    int32_t i = floor( (z-_origo)/_h );
    if( i < 3 )
	i = 3;
    else if( i > (int32_t)_Bz.size()-5 )
	i = _Bz.size()-5;
    double s = (z-(_origo+_h*i))/_h;
    double is = 1.0-s;
    double Bz[7];

    // Calculate function values at z, z+/-h, z+/-2h and z+/-3h
    for( int32_t a = 0; a <= 6; a++ ) {
	int32_t j = i+a-3;
	Bz[a] = is*_Bz[j] + s*_Bz[j+1];
    }

    // Calculate derivatives
    Bder[0] = Bz[3];

    Bder[1] = (-Bz[0] + 9.0*Bz[1] - 45.0*Bz[2] + 45.0*Bz[4] - 9.0*Bz[5] + Bz[6]) / (60.0*_h);

    double h = _h*_h;
    Bder[2] = (2.0*Bz[0] - 27.0*Bz[1] + 270.0*Bz[2] - 490.0*Bz[3] + 270.0*Bz[4] - 27.0*Bz[5] + 2.0*Bz[6]) / (180.0*h);

    h *= _h;
    Bder[3] = (Bz[0] - 8.0*Bz[1] + 13.0*Bz[2] - 13.0*Bz[4] + 8.0*Bz[5] - Bz[6]) / (8.0*h);

    h *= _h;
    Bder[4] = (-Bz[0] + 12.0*Bz[1] - 39.0*Bz[2] + 56.0*Bz[3] - 39.0*Bz[4] + 12.0*Bz[5] - Bz[6]) / (6.0*h);

    h *= _h;
    Bder[5] = (-Bz[0] + 4.0*Bz[1] - 5.0*Bz[2] + 5.0*Bz[4] - 4.0*Bz[5] + Bz[6]) / (2.0*h);

    h *= _h;
    Bder[6] = (Bz[0] - 6.0*Bz[1] + 15.0*Bz[2] - 20.0*Bz[3] + 15.0*Bz[4] - 6.0*Bz[5] + Bz[6]) / h;
}


const Vec3D AxisymmetricVectorField::eval_fdm( double z, double r ) const
{
    double Bder[7];
    double r2 = r*r;
    double r3 = r2*r;
    double r4 = r2*r2;
    double r5 = r4*r;
    double r6 = r3*r3;
    get_fdm_derivatives( Bder, z );

    double Bz, Br;
    if( _order == 0 ) {
	Bz = Bder[0];
	Br = 0.0;
    } else if(_order == 1 ) {
	Bz = Bder[0];
	Br = -r*Bder[1]/2.0;
    } else if(_order == 2 ) {
	Bz = Bder[0] - r2*Bder[2]/4.0;
	Br = -r*Bder[1]/2.0;
    } else if(_order == 3 ) {
	Bz = Bder[0] - r2*Bder[2]/4.0;
	Br = -r*Bder[1]/2.0 + r3*Bder[3]/16.0;
    } else if(_order == 4 ) {
	Bz = Bder[0] - r2*Bder[2]/4.0 + r4*Bder[4]/64.0;
	Br = -r*Bder[1]/2.0 + r3*Bder[3]/16.0;
    } else if(_order == 5 ) {
	Bz = Bder[0] - r2*Bder[2]/4.0 + r4*Bder[4]/64.0;
	Br = -r*Bder[1]/2.0 + r3*Bder[3]/16.0 - r5*Bder[5]/384.0;
    } else if(_order >= 6 ) {
	Bz = Bder[0] - r2*Bder[2]/4.0 + r4*Bder[4]/64.0 - r6*Bder[6]/2304.0;
	Br = -r*Bder[1]/2.0 + r3*Bder[3]/16.0 - r5*Bder[5]/384.0;
    }

    return( Vec3D( Bz, Br, 0.0 ) );
}


const Vec3D AxisymmetricVectorField::operator()( const Vec3D &x ) const
{
    if( _geom_mode == MODE_3D ) {
	Vec3D B;
	if( _spline )
	    B = eval_spline( x[2], sqrt(x[0]*x[0] + x[1]*x[1]) );
	else
	    B = eval_fdm( x[2], sqrt(x[0]*x[0] + x[1]*x[1]) );
	double theta = atan2( x[1], x[0] );
	return( Vec3D( B[1]*cos(theta), B[1]*sin(theta), B[0] ) );
    } else {
	if( _spline )
	    return( eval_spline( x[0], x[1] ) );
	else
	    return( eval_fdm( x[0], x[1] ) );
    }
}


