/*! \file particledatabaseimp.cpp
 *  \brief Particle database implementation
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


#include <fstream>
#include "particledatabaseimp.hpp"
#include "compmath.hpp"
#include "random.hpp"
#include "constants.hpp"
#include "mat3d.hpp"



/* ******************************************************************************************* *
 * ParticleDataBaseImp                                                                         *
 * ******************************************************************************************* */


ParticleDataBaseImp::ParticleDataBaseImp( ParticleDataBase *pdb, const Geometry &geom )
    : _geom(geom), _epsabs(1e-6), _epsrel(1e-6), _intrp(TRAJECTORY_INTERPOLATION_POLYNOMIAL), 
      _scharge_dep(SCHARGE_DEPOSITION_PIC), _maxsteps(1000), 
      _maxt(1e-3), _save_points(false), _trajdiv(1), _rhosum(0.0), _iteration(0), 
      _relativistic(false), _surface_collision(false),
      _bsup_cb(NULL), _thand_cb(NULL), _tend_cb(NULL), _tsur_cb(NULL),
      _pdb(pdb)
{
    for( size_t a = 0; a < 6; a++ )
	_mirror[a] = false;
}


ParticleDataBaseImp::ParticleDataBaseImp( ParticleDataBase *pdb, std::istream &s, const Geometry &geom )
    : _geom(geom), _pdb(pdb)
{
    ibsimu.message( 1 ) << "Constructing ParticleDataBase from stream\n";
    ibsimu.inc_indent();

    _epsabs = read_double( s );
    _epsrel = read_double( s );
    _intrp = (trajectory_interpolation_e)read_int8( s );
    _scharge_dep = (scharge_deposition_e)read_int8( s );
    _maxsteps = read_int32( s );
    _maxt = read_double( s );
    _save_points = read_int8( s );
    _trajdiv = read_int32( s );
    for( uint32_t a = 0; a < 6; a++ )
	_mirror[a] = read_int8( s );
    _rhosum = read_double( s );
    _stat = ParticleStatistics( s );
    _iteration = read_int32( s );
    _relativistic = false; // Should be loaded
    _surface_collision = false; // Should be loaded
    _bsup_cb = NULL;
    _thand_cb = NULL;
    _tend_cb = NULL;
    _tsur_cb = NULL;
}


ParticleDataBaseImp::ParticleDataBaseImp( const ParticleDataBaseImp &pdb )
    : _geom(pdb._geom), _epsabs(pdb._epsabs), _epsrel(pdb._epsrel), _intrp(pdb._intrp), 
      _scharge_dep(pdb._scharge_dep), _maxsteps(pdb._maxsteps), _maxt(pdb._maxt), 
      _save_points(pdb._save_points), _trajdiv(pdb._trajdiv), _rhosum(pdb._rhosum), 
      _stat(pdb._stat), _iteration(0), _relativistic(pdb._relativistic), 
      _surface_collision(pdb._surface_collision), _bsup_cb(pdb._bsup_cb), 
      _thand_cb(pdb._thand_cb), _tend_cb(pdb._tend_cb), _tsur_cb(pdb._tsur_cb),
      _pdb(pdb._pdb)
{
    for( uint32_t a = 0; a < 6; a++ )
	_mirror[a] = pdb._mirror[a];
}


const ParticleDataBaseImp &ParticleDataBaseImp::operator=( const ParticleDataBaseImp &pdb )
{
    return( *this );
}


ParticleDataBaseImp::~ParticleDataBaseImp()
{

}


double ParticleDataBaseImp::energy_to_velocity( double E, double m )
{
    if( E < 5.0e-9*m*SPEED_C2 ) {
	return( sqrt(2.0*E/m) );
    } else {
	double x = E/(m*SPEED_C2)+1.0;
	double beta = sqrt( 1.0 - 1.0/(x*x) );
	return( beta*SPEED_C );
    }
}


void ParticleDataBaseImp::set_accuracy( double epsabs, double epsrel )
{
    _epsabs = epsabs;
    _epsrel = epsrel;
}


void ParticleDataBaseImp::set_bfield_suppression( const CallbackFunctorD_V *bsup_cb )
{
    _bsup_cb = bsup_cb;
}


void ParticleDataBaseImp::set_trajectory_handler_callback( TrajectoryHandlerCallback *thand_cb )
{
    _thand_cb = thand_cb;
}


void ParticleDataBaseImp::set_trajectory_end_callback( TrajectoryEndCallback *tend_cb )
{
    _tend_cb = tend_cb;
}


void ParticleDataBaseImp::set_trajectory_surface_collision_callback( TrajectorySurfaceCollisionCallback *tsur_cb )
{
    _tsur_cb = tsur_cb;
}


void ParticleDataBaseImp::set_relativistic( bool enable )
{
    _relativistic = enable;
}


void ParticleDataBaseImp::set_surface_collision( bool surface_collision )
{
    _surface_collision = surface_collision;
}


void ParticleDataBaseImp::set_polyint( bool polyint )
{
    if( polyint )
	_intrp = TRAJECTORY_INTERPOLATION_POLYNOMIAL;
    else
	_intrp = TRAJECTORY_INTERPOLATION_LINEAR;
}

    
bool ParticleDataBaseImp::get_polyint( void ) const
{
    if( _intrp == TRAJECTORY_INTERPOLATION_POLYNOMIAL )
	return( true );
    else
	return( false );
}


void ParticleDataBaseImp::set_trajectory_interpolation( trajectory_interpolation_e intrp )
{
    _intrp = intrp;
}


trajectory_interpolation_e ParticleDataBaseImp::get_trajectory_interpolation( void ) const
{
    return( _intrp );
}


void ParticleDataBaseImp::set_scharge_deposition( scharge_deposition_e type )
{
    _scharge_dep = type;
}


scharge_deposition_e ParticleDataBaseImp::get_scharge_deposition( void ) const
{
    return( _scharge_dep );
}

    
void ParticleDataBaseImp::set_max_steps( uint32_t maxsteps )
{
    if( maxsteps <= 0 )
        throw( Error( ERROR_LOCATION, "invalid parameter" ) );
    _maxsteps = maxsteps;
}


void ParticleDataBaseImp::set_max_time( double maxt )
{
    if( maxt <= 0.0 )
        throw( Error( ERROR_LOCATION, "invalid parameter" ) );
    _maxt = maxt;
}


void ParticleDataBaseImp::set_save_all_points( bool save_points )
{
    _save_points = save_points;
}


void ParticleDataBaseImp::set_save_trajectories( uint32_t div )
{
    _trajdiv = div;
}


uint32_t ParticleDataBaseImp::get_save_trajectories( void ) const
{
    return( _trajdiv );
}


void ParticleDataBaseImp::set_mirror( const bool mirror[6] )
{
    _mirror[0] = mirror[0];
    _mirror[1] = mirror[1];
    _mirror[2] = mirror[2];
    _mirror[3] = mirror[3];
    if( geom_mode() == MODE_3D ) {
        _mirror[4] = mirror[4];
        _mirror[5] = mirror[5];
    } else {
        _mirror[4] = false;
        _mirror[5] = false;
    }
}


void ParticleDataBaseImp::get_mirror( bool mirror[6] ) const
{
    mirror[0] = _mirror[0];
    mirror[1] = _mirror[1];
    mirror[2] = _mirror[2];
    mirror[3] = _mirror[3];
    mirror[4] = _mirror[4];
    mirror[5] = _mirror[5];
}


int ParticleDataBaseImp::get_iteration_number( void ) const
{
    return( _iteration );
}


double ParticleDataBaseImp::get_rhosum( void ) const
{
    return( _rhosum );
}


void ParticleDataBaseImp::set_rhosum( double rhosum )
{
    _rhosum = rhosum;
}


const ParticleStatistics &ParticleDataBaseImp::get_statistics( void ) const
{
    return( _stat );
}


void ParticleDataBaseImp::save( std::ostream &os ) const
{
    write_double( os, _epsabs );
    write_double( os, _epsrel );
    write_int8( os, _intrp );
    write_int8( os, _scharge_dep );
    write_int32( os, _maxsteps );
    write_double( os, _maxt );
    write_int8( os, _save_points );
    write_int32( os, _trajdiv );
    for( uint32_t a = 0; a < 6; a++ )
	write_int8( os, _mirror[a] );
    write_double( os, _rhosum );
    _stat.save( os );
    write_int32( os, _iteration );
    // relativistic missing
    // surface collision missing
}


void ParticleDataBaseImp::debug_print( std::ostream &os ) const 
{
    os << "epsabs = "      << _epsabs << "\n";
    os << "epsrel = "      << _epsrel << "\n";
    os << "intrp = "       << _intrp << "\n";
    os << "scharge_dep = " << _scharge_dep << "\n";
    os << "maxsteps = "    << _maxsteps << "\n";
    os << "maxt = "        << _maxt << "\n";
    os << "save_points = " << _save_points << "\n";
    os << "trajdiv = "     << _trajdiv << "\n";
    os << "mirror = (";
    for( uint32_t a = 0; a < 5; a++ )
	os << _mirror[a] << ", ";
    os << _mirror[5] << ")\n";
    os << "rhosum = "      << _rhosum << "\n";
    os << "iteration = "   << _iteration << "\n";
    os << "relativistic = "   << _relativistic << "\n";
    os << "surface_collision = "   << _surface_collision << "\n";
}


/* ******************************************************************************************* *
 * ParticleDataBase2DImp                                                                       *
 * ******************************************************************************************* */


ParticleDataBase2DImp::ParticleDataBase2DImp( ParticleDataBase *pdb, const Geometry &geom )
    : ParticleDataBasePPImp<ParticleP2D>(pdb,geom)
{
    
}


ParticleDataBase2DImp::ParticleDataBase2DImp( const ParticleDataBase2DImp &pdb )
    : ParticleDataBasePPImp<ParticleP2D>(pdb)
{

}


ParticleDataBase2DImp::ParticleDataBase2DImp( ParticleDataBase *pdb, std::istream &s, const Geometry &geom )
    : ParticleDataBasePPImp<ParticleP2D>(pdb,s,geom)
{
}


ParticleDataBase2DImp::~ParticleDataBase2DImp()
{

}


const ParticleDataBase2DImp &ParticleDataBase2DImp::operator=( const ParticleDataBase2DImp &pdb )
{
    return( *this );
}


void ParticleDataBase2DImp::add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
						       double v, double dvp, double dvt, 
						       double x1, double y1, double x2, double y2 )
{
    ibsimu.message( 1 ) << "Defining a 2d beam\n";

    _particles.reserve( _particles.size()+N );

    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 2 );
    else
	rng = new MTRandom( 2 );
    rng->set_transformation( 0, Gaussian_Transformation() );
    rng->set_transformation( 1, Gaussian_Transformation() );

    m *= MASS_U;
    q *= CHARGE_E;
    double s = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    _rhosum += J/v;
    double IQ = J*s/N; // IQ is in units A/m for 2D
    double vt[2], dv[2];
    ParticleP2D x;
    x[0] = 0.0;
    Vec3D transverse( x2-x1, y2-y1, 0.0 );
    transverse /= transverse.norm2();
    Vec3D parallel( transverse[1], -transverse[0], 0.0 );

    double Isum = 0.0;
    for( uint32_t a = 0; a < N; a++ ) {
	x[1] = x1 + (x2-x1)*(a+0.5)/((double)N);
	x[3] = y1 + (y2-y1)*(a+0.5)/((double)N);

	rng->get( vt );
	dv[0] = transverse[0]*dvt*vt[0] + parallel[0]*dvp*vt[1];
	dv[1] = transverse[1]*dvt*vt[0] + parallel[1]*dvp*vt[1];
	x[2] = parallel[0]*v + dv[0];
	x[4] = parallel[1]*v + dv[1];

	add_particle( Particle2D( IQ, q, m, x ) );
	Isum += IQ;
    }

    delete( rng );
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A/m\n";
}


void ParticleDataBase2DImp::add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
						     double E, double Tp, double Tt, 
						     double x1, double y1, double x2, double y2 )
{
    ibsimu.message( 1 ) << "Defining a 2d beam\n";

    _particles.reserve( _particles.size()+N );

    // Convert input parameters
    m *= MASS_U;
    q *= CHARGE_E;
    E *= CHARGE_E;
    Tp *= CHARGE_E;
    Tt *= CHARGE_E;
    double v = energy_to_velocity(E,m);
    double dvp = sqrt(Tp/m);
    double dvt = sqrt(Tt/m);

    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 2 );
    else
	rng = new MTRandom( 2 );
    rng->set_transformation( 0, Gaussian_Transformation() );
    rng->set_transformation( 1, Gaussian_Transformation() );

    double s = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    _rhosum += J/v;
    double IQ = J*s/N; // IQ is in units A/m for 2D
    double vt[2];
    ParticleP2D x;
    x[0] = 0.0;
    Vec3D transverse( x2-x1, y2-y1, 0.0 );
    transverse /= transverse.norm2();
    Vec3D parallel( transverse[1], -transverse[0], 0.0 );

    double Isum = 0.0;
    for( uint32_t a = 0; a < N; a++ ) {
	x[1] = x1 + (x2-x1)*(a+0.5)/((double)N);
	x[3] = y1 + (y2-y1)*(a+0.5)/((double)N);

	rng->get( vt );
	double pveld = dvp*vt[1];
	double pvel = sqrt( 2.0*E/m + pveld*pveld );
	x[2] = transverse[0]*dvt*vt[0] + parallel[0]*pvel;
	x[4] = transverse[1]*dvt*vt[0] + parallel[1]*pvel;

	add_particle( Particle2D( IQ, q, m, x ) );
	Isum += IQ;
    }

    delete( rng );
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A/m\n";
}


void ParticleDataBase2DImp::add_2d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
							   double a, double b, double e,
							   double Ex, double x0, double y0 )
{
    ibsimu.message( 1 ) << "Defining a 2d beam using Twiss parameters\n";

    _particles.reserve( _particles.size()+N );

    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 2 );
    else
	rng = new MTRandom( 2 );

    e *= 4.0; // KV-distributed emittance limited inside 4*e_rms ellipse

    m *= MASS_U;
    q *= CHARGE_E;

    double g = (1.0 + a*a)/b;
    double ymax = sqrt( b*e );
    double ypmax = sqrt( g*e );
    double rn[2];
    double IQ = I/N;

    ParticleP2D x;
    x[0] = 0.0;
    x[1] = x0;
    x[2] = sqrt(2.0*Ex*CHARGE_E/m);

    uint32_t n = 0;
    while( n < N ) {

	rng->get( rn );

	// Randomize (y,y')
	double y = -ymax + 2.0*ymax*rn[0];
	double yp = -ypmax + 2.0*ypmax*rn[1];

	// Check if inside ellipse
	double yp1, yp2;
	int nroots = solve_quadratic( b, 2.0*a*y, g*y*y-e, &yp1, &yp2 );
	if( nroots != 2 || yp < yp1 || yp > yp2 )
	    continue;
	
	// Set up particle
	x[3] = y0 + y;
	x[4] = x[2]*yp;

	add_particle( Particle2D( IQ, q, m, x ) );
	n++;
    }

    delete( rng );
}


void ParticleDataBase2DImp::add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
								 double a, double b, double e,
								 double Ex, double x0, double y0 )
{
    ibsimu.message( 1 ) << "Defining a 2d beam using Twiss parameters\n";

    _particles.reserve( _particles.size()+N );

    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 2 );
    else
	rng = new MTRandom( 2 );
    rng->set_transformation( 0, Gaussian_Transformation() );
    rng->set_transformation( 1, Gaussian_Transformation() );

    m *= MASS_U;
    q *= CHARGE_E;

    double g = (1.0 + a*a)/b;
    double h = 0.5*(b+g);
    double rmaj = sqrt(0.5*e)*(sqrt(h+1.0)+sqrt(h-1.0));
    double rmin = sqrt(0.5*e)*(sqrt(h+1.0)-sqrt(h-1.0));
    double theta = 0.5*atan2( -2.0*a, b-g );
    double w[2], rn[2];
    double IQ = I/N;

    ParticleP2D x;
    x[0] = 0.0;
    x[1] = x0;
    x[2] = sqrt(2.0*Ex*CHARGE_E/m);

    uint32_t n = 0;
    while( n < N ) {

	// Randomize point from gaussian distribution
	rng->get( rn );
	w[0] = rmaj*rn[0];
	w[1] = rmin*rn[1];

	// Rotate to correct angle
	double y  = w[0]*cos(theta) - w[1]*sin(theta);
	double yp = w[0]*sin(theta) + w[1]*cos(theta);

	// Set up particle
	x[3] = y0 + y;
	x[4] = x[2]*yp;

	add_particle( Particle2D( IQ, q, m, x ) );
	n++;
    }

    delete( rng );
}


void ParticleDataBase2DImp::add_tdens_from_segment( MeshScalarField &tdens, double IQ,
						    ParticleP2D &x1, ParticleP2D &x2 ) const
{
    double dx = 0.0;
    double x[2];
    double t[2];
    int32_t i[2];

    for( size_t a = 0; a < 2; a++ ) {
	double xx = ( x2[2*a+1] - x1[2*a+1] );
	dx  += xx*xx;
	x[a] = 0.5*( x1[2*a+1] + x2[2*a+1] );
	i[a] = (int32_t)floor( ( x[a]-tdens.origo(a) ) * tdens.div_h() );
	t[a] = ( x[a]-(i[a]*tdens.h()+tdens.origo(a)) ) * tdens.div_h();
	
	if( i[a] < 0 || i[a] >= (int32_t)(tdens.size(a)-1) )
	    return;
    }
    dx = sqrt( dx );

    double J = IQ*dx; // A
    int p = tdens.size(0)*i[1] + i[0];
    tdens( p )                 += (1.0-t[0])*(1.0-t[1])*J;
    tdens( p+tdens.size(0) )   += (1.0-t[0])*t[1]*J;
    tdens( p+1 )               += t[0]*(1.0-t[1])*J;
    tdens( p+1+tdens.size(0) ) += t[0]*t[1]*J;
}


void ParticleDataBase2DImp::build_trajectory_density_field( MeshScalarField &tdens ) const
{
    tdens.clear();

    // Go through all particle trajectories
    std::vector< ColData<ParticleP2D> > coldata;
    for( size_t a = 0; a < _particles.size(); a++ ) {
	size_t N = _particles[a]->traj_size();
	if( N < 2 )
	    continue;
	ParticleP2D x1 = _particles[a]->traj(0);
	ParticleP2D xlast = x1;
	double IQ = _particles[a]->IQ();
	for( size_t b = 1; b < N; b++ ) {
	    ParticleP2D x2 = _particles[a]->traj(b);

	    // Mesh collision points from x1 to x2
	    ColData<ParticleP2D>::build_coldata_linear( coldata, tdens, x1, x2 );
	    if( coldata.size() == 0 )
		continue;

	    // Process trajectory mesh collision points
	    for( size_t c = 0; c < coldata.size(); c++ ) {

		// Add contribution from segment from xlast to coldata[c]
		add_tdens_from_segment( tdens, IQ, xlast, coldata[c]._x );	    
		xlast = coldata[c]._x;
	    }

	    // Next trajectory line
	    x1 = x2;
	}
    }

    // Normalize to trajectory density, A/m2
    tdens /= (tdens.h()*tdens.h());

    // Fix boundaries, which only get half of the contribution they should
    // This should depend on mirroring properties!
    for( uint32_t i = 0; i < tdens.size(0); i++ ) {
	tdens( i, 0 ) *= 2.0;
	tdens( i, tdens.size(1)-1 ) *= 2.0;
    }
    for( uint32_t j = 0; j < tdens.size(1); j++ ) {
	tdens( 0, j ) *= 2.0;
	tdens( tdens.size(0)-1, j ) *= 2.0;
    }
}


void ParticleDataBase2DImp::save( std::ostream &os ) const
{
    ParticleDataBasePPImp<ParticleP2D>::save( os );
}


void ParticleDataBase2DImp::debug_print( std::ostream &os ) const 
{
    os << "**ParticleDataBase2D\n";
    ParticleDataBasePPImp<ParticleP2D>::debug_print( os );
}


/* ******************************************************************************************* *
 * ParticleDataBaseCylImp                                                                      *
 * ******************************************************************************************* */


ParticleDataBaseCylImp::ParticleDataBaseCylImp( ParticleDataBase *pdb, const Geometry &geom )
    : ParticleDataBasePPImp<ParticlePCyl>(pdb,geom)
{
    
}


ParticleDataBaseCylImp::ParticleDataBaseCylImp( const ParticleDataBaseCylImp &pdb )
    : ParticleDataBasePPImp<ParticlePCyl>(pdb)
{

}


ParticleDataBaseCylImp::ParticleDataBaseCylImp( ParticleDataBase *pdb, std::istream &s, const Geometry &geom )
    : ParticleDataBasePPImp<ParticlePCyl>(pdb,s,geom)
{
}


ParticleDataBaseCylImp::~ParticleDataBaseCylImp()
{

}


const ParticleDataBaseCylImp &ParticleDataBaseCylImp::operator=( const ParticleDataBaseCylImp &pdb )
{
    return( *this );
}


void ParticleDataBaseCylImp::add_tdens_from_segment( MeshScalarField &tdens, double IQ,
						     ParticlePCyl &x1, ParticlePCyl &x2 ) const
{
    double dx = 0.0;
    double x[2];
    double t[2];
    int32_t i[2];

    // x-direction
    double xx = ( x2[1] - x1[1] );
    xx = xx*xx;
    dx += xx;
    x[0] = 0.5*( x1[1] + x2[1] );
    i[0] = (int32_t)floor( ( x[0]-tdens.origo(0) ) * tdens.div_h() );
    t[0] = ( x[0]-(i[0]*tdens.h()+tdens.origo(0)) ) * tdens.div_h();

    // r-direction
    xx = ( x2[3] - x1[3] );
    xx = xx*xx;
    dx += xx;
    x[1] = 0.5*( x1[3] + x2[3] );
    i[1] = (int32_t)floor( ( x[1]-tdens.origo(1) ) * tdens.div_h() );
    double rj1 = i[1]*tdens.h()+tdens.origo(1);
    double rj2 = rj1+tdens.h();
    rj1 = rj1*rj1;
    rj2 = rj2*rj2;
    t[1] = (x[1]*x[1]-rj1) / (rj2-rj1);
    
    dx = sqrt( dx );

    // Check that no deposition is done over field size
    for( size_t a = 0; a < 2; a++ ) {
	if( i[a] < 0 || i[a] >= (int32_t)(tdens.size(a)-1) )
	    return;
	/*
	// Add charge to boundaries when over simulation area
	if( i[a] < 0 ) {
	    i[a] = 0;
	    t[a] = 0.0;
	} else if( i[a] >= (int32_t)tdens.size(a)-1 ) {
	    i[a] = tdens.size(a)-2;
	    t[a] = 1.0;
	}
	*/
    }

    double J = IQ*dx; // A*m
    int p = tdens.size(0)*i[1] + i[0];
    tdens( p )                 += (1.0-t[0])*(1.0-t[1])*J;
    tdens( p+tdens.size(0) )   += (1.0-t[0])*t[1]*J;
    tdens( p+1 )               += t[0]*(1.0-t[1])*J;
    tdens( p+1+tdens.size(0) ) += t[0]*t[1]*J;
}


void ParticleDataBaseCylImp::build_trajectory_density_field( MeshScalarField &tdens ) const
{
    tdens.clear();

    // Go through all particle trajectories
    std::vector< ColData<ParticlePCyl> > coldata;
    for( size_t a = 0; a < _particles.size(); a++ ) {
	size_t N = _particles[a]->traj_size();
	if( N < 2 )
	    continue;
	ParticlePCyl x1 = _particles[a]->traj(0);
	ParticlePCyl xlast = x1;
	double IQ = _particles[a]->IQ();
	for( size_t b = 1; b < N; b++ ) {
	    ParticlePCyl x2 = _particles[a]->traj(b);

	    // Mesh collision points from x1 to x2
	    ColData<ParticlePCyl>::build_coldata_linear( coldata, tdens, x1, x2 );
	    if( coldata.size() == 0 )
		continue;

	    // Process trajectory mesh collision points
	    for( size_t c = 0; c < coldata.size(); c++ ) {

		// Add contribution from segment from xlast to coldata[c]
		add_tdens_from_segment( tdens, IQ, xlast, coldata[c]._x );	    
		xlast = coldata[c]._x;
	    }

	    // Next trajectory line
	    x1 = x2;
	}
    }    

    // Normalize to trajectory density, A/m2
    for( uint32_t i = 0; i < tdens.size(0); i++ ) {
	for( uint32_t j = 0; j < tdens.size(1); j++ ) {
	    if( j == 0 ) {
		double rj2 = tdens.h()+tdens.origo(1);
		tdens( i, j ) /= (M_PI*tdens.h()*(rj2*rj2));
	    } else {
		double rj1 = (j-0.5)*tdens.h()+tdens.origo(1);
		double rj2 = (j+0.5)*tdens.h()+tdens.origo(1);
		tdens( i, j ) /= (M_PI*tdens.h()*(rj2*rj2-rj1*rj1));
	    }
	}
    }

    // Fix boundaries, which only get half of the contribution they should
    // This should depend on mirroring properties!
    for( uint32_t i = 0; i < tdens.size(0); i++ ) {
	tdens( i, 0 ) *= 2.0;
	tdens( i, tdens.size(1)-1 ) *= 2.0;
    }
    for( uint32_t j = 0; j < tdens.size(1); j++ ) {
	tdens( 0, j ) *= 2.0;
	tdens( tdens.size(0)-1, j ) *= 2.0;
    }
}

void ParticleDataBaseCylImp::add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
							double v, double dvp, double dvt, 
							double x1, double y1, double x2, double y2 )
{
    ibsimu.message( 1 ) << "Defining a cylindrical beam\n";

    _particles.reserve( _particles.size()+N );

    // 0: temperature parallel to defining line (transverse to beam)
    // 1: temperature perpendicular to line (parallel to beam)
    // 2: temperature in angular direction (skew velocity)
    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 3 );
    else
	rng = new MTRandom( 3 );
    rng->set_transformation( 0, Gaussian_Transformation() );
    rng->set_transformation( 1, Gaussian_Transformation() );
    rng->set_transformation( 2, Gaussian_Transformation() ); 

    m *= MASS_U;
    q *= CHARGE_E;
    double s = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    _rhosum += J/v;
    double IQ = 2.0*M_PI*J*s/N;
    double vt[3], dv[3];
    ParticlePCyl x;
    x[0] = 0.0;
    Vec3D transverse( x2-x1, y2-y1, 0.0 );
    transverse /= transverse.norm2();
    Vec3D parallel( transverse[1], -transverse[0], 0.0 );

    double Isum = 0.0;
    for( uint32_t a = 0; a < N; a++ ) {
	x[1] = x1 + (x2-x1)*(a+0.5)/((double)N);
	x[3] = y1 + (y2-y1)*(a+0.5)/((double)N);

	rng->get( vt );
	dv[0] = transverse[0]*dvt*vt[0] + parallel[0]*dvp*vt[1];
	dv[1] = transverse[1]*dvt*vt[0] + parallel[1]*dvp*vt[1];
	dv[2] = dvt*vt[2];
	x[2] = parallel[0]*fabs(v) + dv[0];
	x[4] = parallel[1]*fabs(v) + dv[1];
	if( x[3] == 0.0 )
	    x[5] = 0.0;
	else
	    x[5] = dv[2]/x[3];

	add_particle( ParticleCyl( IQ*x[3], q, m, x ) );
	Isum += IQ*x[3];
    }

    delete( rng );
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBaseCylImp::add_2d_beam_with_total_energy( uint32_t N, double J, double q, double m, 
							    double Etot, const ScalarField &epot, 
							    double Tp, double Tt, 
							    double x1, double y1, double x2, double y2 )
{
    ibsimu.message( 1 ) << "Defining a cylindrical beam\n";

    _particles.reserve( _particles.size()+N );

    // Convert input parameters
    m *= MASS_U;
    q *= CHARGE_E;
    Etot *= CHARGE_E;
    Tp *= CHARGE_E;
    Tt *= CHARGE_E;
    double dvp = sqrt(Tp/m);
    double dvt = sqrt(Tt/m);

    // 0: temperature parallel to defining line (transverse to beam)
    // 1: temperature perpendicular to line (parallel to beam)
    // 2: temperature in angular direction (skew velocity)
    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 3 );
    else
	rng = new MTRandom( 3 );
    rng->set_transformation( 0, Gaussian_Transformation() );
    rng->set_transformation( 1, Gaussian_Transformation() );
    rng->set_transformation( 2, Gaussian_Transformation() ); 

    double s = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    double I = 2.0*M_PI*J*s/N;
    double vt[3];
    ParticlePCyl x;
    x[0] = 0.0;
    Vec3D transverse( x2-x1, y2-y1, 0.0 );
    transverse /= transverse.norm2();
    Vec3D parallel( transverse[1], -transverse[0], 0.0 );

    double Isum = 0.0;
    for( uint32_t a = 0; a < N; a++ ) {
	x[1] = x1 + (x2-x1)*(a+0.5)/((double)N);
	x[3] = y1 + (y2-y1)*(a+0.5)/((double)N);

	rng->get( vt );
	double pveld = dvp*vt[1];
	Vec3D loc( x[1], x[3], 0.0 );
	double U = epot( loc );
        double pvel = sqrt( 2.0*(Etot-q*U)/m + pveld*pveld );
	x[2] = transverse[0]*dvt*vt[0] + parallel[0]*pvel;
	x[4] = transverse[1]*dvt*vt[0] + parallel[1]*pvel;
	if( x[3] == 0.0 )
	    x[5] = 0.0;
	else
	    x[5] = dvt*vt[2]/x[3];

	add_particle( ParticleCyl( I*x[3], q, m, x ) );
	Isum += I*x[3];
    }

    delete( rng );
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBaseCylImp::add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
						      double E, double Tp, double Tt, 
						      double x1, double y1, double x2, double y2 )
{
    ibsimu.message( 1 ) << "Defining a cylindrical beam\n";

    _particles.reserve( _particles.size()+N );

    // Convert input parameters
    m *= MASS_U;
    q *= CHARGE_E;
    E *= CHARGE_E;
    Tp *= CHARGE_E;
    Tt *= CHARGE_E;
    double v = energy_to_velocity(E,m);
    double dvp = sqrt(Tp/m);
    double dvt = sqrt(Tt/m);

    // 0: temperature parallel to defining line (transverse to beam)
    // 1: temperature perpendicular to line (parallel to beam)
    // 2: temperature in angular direction (skew velocity)
    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 3 );
    else
	rng = new MTRandom( 3 );
    rng->set_transformation( 0, Gaussian_Transformation() );
    rng->set_transformation( 1, Gaussian_Transformation() );
    rng->set_transformation( 2, Gaussian_Transformation() ); 

    double s = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    _rhosum += J/v;
    double IQ = 2.0*M_PI*J*s/N;
    double vt[3];
    ParticlePCyl x;
    x[0] = 0.0;
    Vec3D transverse( x2-x1, y2-y1, 0.0 );
    transverse /= transverse.norm2();
    Vec3D parallel( transverse[1], -transverse[0], 0.0 );

    double Isum = 0.0;
    for( uint32_t a = 0; a < N; a++ ) {
	x[1] = x1 + (x2-x1)*(a+0.5)/((double)N);
	x[3] = y1 + (y2-y1)*(a+0.5)/((double)N);

	rng->get( vt );
	double pveld = dvp*vt[1];
        double pvel = sqrt( 2.0*E/m + pveld*pveld );
	x[2] = transverse[0]*dvt*vt[0] + parallel[0]*pvel;
	x[4] = transverse[1]*dvt*vt[0] + parallel[1]*pvel;
	if( x[3] == 0.0 )
	    x[5] = 0.0;
	else
	    x[5] = dvt*vt[2]/x[3];

	add_particle( ParticleCyl( IQ*x[3], q, m, x ) );
	Isum += IQ*x[3];
    }

    delete( rng );
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBaseCylImp::add_2d_full_gaussian_beam( uint32_t N, double I, double q, double m,
							double Ex, double Tp, double Tt, 
							double x0, double dr )
{
    ibsimu.message( 1 ) << "Defining a cylindrical full gaussian beam.\n";

    _particles.reserve( _particles.size()+N );

    // Convert input parameters
    m *= MASS_U;
    q *= CHARGE_E;
    Ex *= CHARGE_E;
    Tp *= CHARGE_E;
    Tt *= CHARGE_E;
    //double v = energy_to_velocity(Ex,m);
    double dvp = sqrt(Tp/m);
    double dvt = sqrt(Tt/m);

    // Random number generator for vx, y, vy, z, vz
    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 5 );
    else
	rng = new MTRandom( 5 );
    rng->set_transformation( 0, Gaussian_Transformation() );
    rng->set_transformation( 1, Gaussian_Transformation() );
    rng->set_transformation( 2, Gaussian_Transformation() );
    rng->set_transformation( 3, Gaussian_Transformation() );
    rng->set_transformation( 4, Gaussian_Transformation() );
    double rn[5];

    double IQ = I/N;
    double Isum = 0.0;

    ParticlePCyl x;
    x[0] = 0.0;
    x[1] = x0;

    uint32_t n = 0;
    while( n < N ) {

	rng->get( rn );
	if( rn[1] < 0 || rn[3] < 0 )
	    continue;

	double pveld = dvp*rn[0];
        double vx = sqrt( 2.0*Ex/m + pveld*pveld );
	double y = dr*rn[1];
	double vy = dvt*rn[2];
	double z = dr*rn[3];
	double vz = dvt*rn[4];

	// Convert to cylindrical coordinates
	double r  = sqrt( y*y + z*z );
	if( r == 0.0 )
	    continue;
	double alpha = atan2( z, y );
	double sina = sin(alpha);
	double cosa = cos(alpha);
	double vr = vy*cosa + vz*sina;
	double va = -vy*sina + vz*cosa;

	// Set up particle
	x[2] = vx;
	x[3] = r;
	x[4] = vr;
	x[5] = va/r;

	add_particle( ParticleCyl( IQ, q, m, x ) );
	Isum += IQ;
	n++;
    }

    delete( rng );
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBaseCylImp::add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
								  double a, double b, double e,
								  double Ex, double x0 )
{
    ibsimu.message( 1 ) << "Defining a cylindrical beam using Twiss parameters\n";

    _particles.reserve( _particles.size()+N );

    m *= MASS_U;
    q *= CHARGE_E;

    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 4 );
    else
	rng = new MTRandom( 4 );
    rng->set_transformation( 0, Gaussian_Transformation() );
    rng->set_transformation( 1, Gaussian_Transformation() );
    rng->set_transformation( 2, Gaussian_Transformation() );
    rng->set_transformation( 3, Gaussian_Transformation() );

    double w[4], rn[4];

    double g = (1.0 + a*a)/b;
    double h = 0.5*(b+g);
    double rmaj = sqrt(0.5*e)*(sqrt(h+1.0)+sqrt(h-1.0));
    double rmin = sqrt(0.5*e)*(sqrt(h+1.0)-sqrt(h-1.0));
    double theta = 0.5*atan2( -2.0*a, b-g );

    double IQ = I/N;

    ParticlePCyl x;
    x[0] = 0.0;
    x[1] = x0;
    x[2] = sqrt(2.0*Ex*CHARGE_E/m);

    uint32_t n = 0;
    while( n < N ) {

	// Randomize point from gaussian distribution
	rng->get( rn );
	w[0] = rmaj*rn[0];
	w[1] = rmin*rn[1];
	w[2] = rmaj*rn[2];
	w[3] = rmin*rn[3];

	// Rotate to correct angle
	double sint = sin(theta);
	double cost = cos(theta);
	double y  = w[0]*cost - w[1]*sint;
	double yp = w[0]*sint + w[1]*cost;
	double z  = w[2]*cost - w[3]*sint;
	double zp = w[2]*sint + w[3]*cost;

	// Convert to cylindrical coordinates
	double r  = sqrt( y*y + z*z );
	if( r == 0.0 ) // reject center point
	    continue;
	double alpha = atan2( z, y );
	double sina = sin(alpha);
	double cosa = cos(alpha);
	double rp = yp*cosa + zp*sina;
	double ap = -yp*sina + zp*cosa;

	// Set up particle
	x[3] = r;
	x[4] = x[2]*rp;
	x[5] = x[2]*ap/r;

	add_particle( ParticleCyl( IQ, q, m, x ) );
	n++;
    }

    delete( rng );
}


// Find gap where x is between cum[b] and cum[b+1] by forking
uint32_t ParticleDataBaseCylImp::bisect_cumulative_array( const std::vector<double> &cum, double x )
{
    uint32_t b;
    uint32_t lowi = 0; 
    if( cum.size() <= 1 )
	return( 0 );
    uint32_t highi = cum.size()-1;

    while( 1 ) {

	// Fork
	b = (lowi+highi)/2;
	double val = cum[b];

	// Update limits
	if( x < val )
	    highi = b;
	else
	    lowi = b;

	// Test for convergence
	if( lowi+1 == highi )
	    return( lowi );
    }
}


void ParticleDataBaseCylImp::export_path_manager_data( const std::string &filename, 
						       double ref_E, double ref_q, double ref_m, 
						       double val, uint32_t Np ) const
{
    // Calculate reference particle
    ref_m *= MASS_U;
    double ref_v = sqrt(2.0*ref_E*CHARGE_E/ref_m);
    double ref_p = ref_m*ref_v;
    const double m_to_gevc2 = SPEED_C2/CHARGE_E/1.0e9;
    const double p_to_gevc = SPEED_C/CHARGE_E/1.0e9;

    ibsimu.message( 1 ) << "Making trajectory diagnostics at plane for Path Manager output\n" 
			<< "  ref_E = " << ref_E << " eV\n"
			<< "  ref_q = " << ref_q << " e\n"
			<< "  ref_m = " << ref_m/MASS_U << " u\n"
			<< "  x = " << val << "\n";

    ibsimu.message( 1 ) << "  Outputting to \'" << filename << "\'\n";

    std::ofstream ofile( filename.c_str() );
    if( !ofile.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\'" ) );

    // Scan through particle trajectory points and build data
    std::vector<ParticleCyl> tdata;
    std::vector<ParticlePCyl> intsc;
    std::vector<double> Icum;
    Icum.push_back( 0.0 );
    double Isum = 0.0;
    for( size_t a = 0; a < _particles.size(); a++ ) {
	size_t N = _particles[a]->traj_size();
	if( N < 2 )
	    continue;
	ParticlePCyl x1 = _particles[a]->traj(0);
	for( size_t b = 1; b < N; b++ ) {
	    ParticlePCyl x2 = _particles[a]->traj(b);
	    intsc.clear();
	    size_t nintsc;
	    if( b == 1 )
		nintsc = ParticlePCyl::trajectory_intersections_at_plane( intsc, 0, val, x1, x2, -1 );
	    else if( b == N-1 )
		nintsc = ParticlePCyl::trajectory_intersections_at_plane( intsc, 0, val, x1, x2, +1 );
	    else
		nintsc = ParticlePCyl::trajectory_intersections_at_plane( intsc, 0, val, x1, x2, 0 );
	    for( size_t c = 0; c < nintsc; c++ ) {
		Isum += _particles[a]->IQ();
		Icum.push_back( Icum.back()+_particles[a]->IQ() );
		tdata.push_back( ParticleCyl( _particles[a]->IQ(), _particles[a]->q(),
					      _particles[a]->m(), intsc[c] ) );
	    }
	    
	    x1 = x2;
	}
    }

    // Write header
    time_t tt = time(NULL);
    char *timebuf = ctime( &tt );
    ofile << "Beam data from IBSimu, " << timebuf; // contains newline
    ofile << "Total beam: " << Isum << " A\n";
    ofile << std::setw(17) << ref_p*p_to_gevc  << " !REFERENCE MOMENTUM [GeV/c]\n";
    ofile << std::setw(17) << 0.0              << " !REFERENCE PHASE [rad]\n";
    ofile << std::setw(17) << 0.0              << " !FREQUENCY [Hz]\n";
    ofile << std::setw(17) << ref_m*m_to_gevc2 << " !REFERENCE MASS [GeV/c2]\n";
    ofile << std::setw(17) << ref_q            << " !REFERENCE CHARGE STATE\n";
    ofile << std::setw(17) << Np               << " !NUMBER OF PARTICLES\n";

    for( uint32_t a = 0; a < Np; a++ ) {	

	// Find gap where x is between Icum[b] and Icum[b+1]
	double s = Isum*a/(Np-1); // From 0 to Isum
	uint32_t b = bisect_cumulative_array( Icum, s );

	// Randomize azimuthal angle and phase
	double ang = -M_PI + 2.0*M_PI*rand() / RAND_MAX;
	double phi = -M_PI + 2.0*M_PI*rand() / RAND_MAX;
	double r = tdata[b][3];
	double rp = tdata[b][4]/tdata[b][2];
	double ap = r*tdata[b][5]/tdata[b][2];
	double sint = sin( ang );
	double cost = cos( ang );
	double x = r*cost;
	double xp = rp*cost - ap*sint;
	double y = r*sint;
	double yp = rp*sint + ap*cost;
	double pq = tdata[b].q();
	double pm = tdata[b].m();
	double dpp = (pm*tdata[b][2]-ref_p)/ref_p;
	double q = pq/CHARGE_E;
	double m = pm*m_to_gevc2;

	ofile << std::setw(17) << a << " ";                                              // index
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << x << " ";   // x
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << xp << " ";  // x'
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << y << " ";   // y
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << yp << " ";  // y'
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << phi << " "; // phase
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << dpp << " "; // dp/p
	ofile << std::setw(17) << 0 << " ";                                              // flag
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << q << " ";   // q
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << m << " ";   // m
	ofile << "\n";
    }

    ibsimu.message( 1 ) << "  number of trajectories input = " << tdata.size() << "\n";
    ibsimu.message( 1 ) << "  number of trajectories ouput = " << Np << "\n";
    ibsimu.message( 1 ) << "  total current = " << Isum << " A\n";
}


void ParticleDataBaseCylImp::save( std::ostream &os ) const
{
    ParticleDataBasePPImp<ParticlePCyl>::save( os );
}


void ParticleDataBaseCylImp::debug_print( std::ostream &os ) const 
{
    os << "**ParticleDataBaseCyl\n";
    ParticleDataBasePPImp<ParticlePCyl>::debug_print( os );
}


/* ******************************************************************************************* *
 * ParticleDataBase3DImp                                                                       *
 * ******************************************************************************************* */


ParticleDataBase3DImp::ParticleDataBase3DImp( ParticleDataBase *pdb, const Geometry &geom )
    : ParticleDataBasePPImp<ParticleP3D>(pdb,geom)
{
    
}


ParticleDataBase3DImp::ParticleDataBase3DImp( const ParticleDataBase3DImp &pdb )
    : ParticleDataBasePPImp<ParticleP3D>(pdb)
{

}


ParticleDataBase3DImp::ParticleDataBase3DImp( ParticleDataBase *pdb, std::istream &s, const Geometry &geom )
    : ParticleDataBasePPImp<ParticleP3D>(pdb,s,geom)
{
}


ParticleDataBase3DImp::~ParticleDataBase3DImp()
{

}


const ParticleDataBase3DImp &ParticleDataBase3DImp::operator=( const ParticleDataBase3DImp &pdb )
{
    return( *this );
}


void ParticleDataBase3DImp::add_cylindrical_beam_with_velocity( uint32_t N, double J, double q, double m, 
								double v, double dvp, double dvt, Vec3D c, 
								Vec3D dir1, Vec3D dir2, double r )
{
    ibsimu.message( 1 ) << "Defining a cylindrical beam\n";

    _particles.reserve( _particles.size()+N );

    // Random number generator for two positions and three velocities (gaussian)
    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 5 );
    else
	rng = new MTRandom( 5 );
    rng->set_transformation( 2, Gaussian_Transformation() );
    rng->set_transformation( 3, Gaussian_Transformation() );
    rng->set_transformation( 4, Gaussian_Transformation() );
    double qx[5];
    double px[6];

    m *= MASS_U;
    q *= CHARGE_E;
    _rhosum += J/v;
    double IQ = M_PI*r*r*J/N;

    // Calculate and check base vectors
    Vec3D dir3 = cross( dir1, dir2 );
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( comp_isnan( dir1[0] ) || comp_isnan( dir2[0] ) || comp_isnan( dir3[0] ) )
	throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );

    // Prepare particle
    ParticleP3D x;
    x[0] = 0.0;

    double Isum = 0.0;
    uint32_t a = 0;
    while( a < N ) {

	rng->get( qx );

	// Calculate in natural (dir1,dir2,dir3) coordinates
	px[0] = -r + 2.0*r*qx[0];
	px[1] = -r + 2.0*r*qx[1];
	px[2] = 0.0;
	px[3] = dvt*qx[2];
	px[4] = dvt*qx[3];
	px[5] = dvp*qx[4] + v;

	// Reject points outside unit radius circle
	if( px[0]*px[0] + px[1]*px[1] > r*r )
	    continue;

	// Map to world coordinates
	x[1] = dir1[0]*px[0] + dir2[0]*px[1] + dir3[0]*px[2] + c[0];
	x[2] = dir1[0]*px[3] + dir2[0]*px[4] + dir3[0]*px[5];
	x[3] = dir1[1]*px[0] + dir2[1]*px[1] + dir3[1]*px[2] + c[1];
	x[4] = dir1[1]*px[3] + dir2[1]*px[4] + dir3[1]*px[5];
	x[5] = dir1[2]*px[0] + dir2[2]*px[1] + dir3[2]*px[2] + c[2];
	x[6] = dir1[2]*px[3] + dir2[2]*px[4] + dir3[2]*px[5];

	add_particle( Particle3D( IQ, q, m, x ) );
	Isum += IQ;
	a++;
    }

    delete( rng );
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBase3DImp::add_cylindrical_beam_with_energy( uint32_t N, double J, double q, double m, 
							      double E, double Tp, double Tt, Vec3D c,
							      Vec3D dir1, Vec3D dir2, double r )
{
    ibsimu.message( 1 ) << "Defining a cylindrical beam\n";

    _particles.reserve( _particles.size()+N );

    // Convert input parameters
    m *= MASS_U;
    q *= CHARGE_E;
    E *= CHARGE_E;
    Tp *= CHARGE_E;
    Tt *= CHARGE_E;
    double v = energy_to_velocity(E,m);
    double dvp = sqrt(Tp/m);
    double dvt = sqrt(Tt/m);

    // Random number generator for two positions and three velocities (gaussian)
    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 5 );
    else
	rng = new MTRandom( 5 );
    rng->set_transformation( 2, Gaussian_Transformation() );
    rng->set_transformation( 3, Gaussian_Transformation() );
    rng->set_transformation( 4, Gaussian_Transformation() );
    double qx[5];
    double px[6];

    _rhosum += J/v;
    double IQ = M_PI*r*r*J/N;

    // Calculate and check base vectors
    Vec3D dir3 = cross( dir1, dir2 );
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( comp_isnan( dir1[0] ) || comp_isnan( dir2[0] ) || comp_isnan( dir3[0] ) )
	throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );
	
    // Prepare particle
    ParticleP3D x;
    x[0] = 0.0;

    double Isum = 0.0;
    uint32_t a = 0;
    while( a < N ) {

	rng->get( qx );

	// Position in natural (dir1,dir2,dir3) coordinates
	px[0] = -r + 2.0*r*qx[0];
	px[1] = -r + 2.0*r*qx[1];
	px[2] = 0.0;

	// Reject points outside radius r circle
 	if( px[0]*px[0] + px[1]*px[1] > r*r )
	    continue;

	// Velocities in natural (dir1,dir2,dir3) coordinates
	px[3] = dvt*qx[2];
	px[4] = dvt*qx[3];
	double pveld = dvp*qx[4];
        px[5] = sqrt( 2.0*E/m + pveld*pveld );

	// Map to world coordinates
	x[1] = dir1[0]*px[0] + dir2[0]*px[1] + dir3[0]*px[2] + c[0];
	x[2] = dir1[0]*px[3] + dir2[0]*px[4] + dir3[0]*px[5];
	x[3] = dir1[1]*px[0] + dir2[1]*px[1] + dir3[1]*px[2] + c[1];
	x[4] = dir1[1]*px[3] + dir2[1]*px[4] + dir3[1]*px[5];
	x[5] = dir1[2]*px[0] + dir2[2]*px[1] + dir3[2]*px[2] + c[2];
	x[6] = dir1[2]*px[3] + dir2[2]*px[4] + dir3[2]*px[5];

	add_particle( Particle3D( IQ, q, m, x ) );
	Isum += IQ;
	a++;
    }

    delete( rng );
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBase3DImp::add_cylindrical_beam_with_total_energy( uint32_t N, double J, double q, double m, 
								    double Etot, const ScalarField &epot, 
								    double Tp, double Tt, Vec3D c,
								    Vec3D dir1, Vec3D dir2, double r )
{
    ibsimu.message( 1 ) << "Defining a cylindrical beam\n";

    _particles.reserve( _particles.size()+N );

    // Convert input parameters
    m *= MASS_U;
    q *= CHARGE_E;
    Etot *= CHARGE_E;
    Tp *= CHARGE_E;
    Tt *= CHARGE_E;
    double dvp = sqrt(Tp/m);
    double dvt = sqrt(Tt/m);

    // Random number generator for two positions and three velocities (gaussian)
    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 5 );
    else
	rng = new MTRandom( 5 );
    rng->set_transformation( 2, Gaussian_Transformation() );
    rng->set_transformation( 3, Gaussian_Transformation() );
    rng->set_transformation( 4, Gaussian_Transformation() );
    double qx[5];
    double px[6];

    double I = M_PI*r*r*J/N;

    // Calculate and check base vectors
    Vec3D dir3 = cross( dir1, dir2 );
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( comp_isnan( dir1[0] ) || comp_isnan( dir2[0] ) || comp_isnan( dir3[0] ) )
	throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );

    // Prepare particle
    ParticleP3D x;
    x[0] = 0.0;

    double Isum = 0.0;
    uint32_t a = 0;
    while( a < N ) {

	rng->get( qx );

	// Position in natural (dir1,dir2,dir3) coordinates
	px[0] = -r + 2.0*r*qx[0];
	px[1] = -r + 2.0*r*qx[1];
	px[2] = 0.0;

	// Reject points outside radius r circle
 	if( px[0]*px[0] + px[1]*px[1] > r*r )
	    continue;

	// Map location to world coordinates
	x[1] = dir1[0]*px[0] + dir2[0]*px[1] + dir3[0]*px[2] + c[0];
	x[3] = dir1[1]*px[0] + dir2[1]*px[1] + dir3[1]*px[2] + c[1];
	x[5] = dir1[2]*px[0] + dir2[2]*px[1] + dir3[2]*px[2] + c[2];

	// Velocities in natural (dir1,dir2,dir3) coordinates
	px[3] = dvt*qx[2];
	px[4] = dvt*qx[3];
	double pveld = dvp*qx[4];
	Vec3D loc( x[1], x[3], x[5] );
	double U = epot( loc );
        px[5] = sqrt( 2.0*(Etot-q*U)/m + pveld*pveld );

	// Map velocity to world coordinates
	x[4] = dir1[1]*px[3] + dir2[1]*px[4] + dir3[1]*px[5];
	x[2] = dir1[0]*px[3] + dir2[0]*px[4] + dir3[0]*px[5];
	x[6] = dir1[2]*px[3] + dir2[2]*px[4] + dir3[2]*px[5];

	add_particle( Particle3D( I, q, m, x ) );
	Isum += I;
	a++;
    }

    delete( rng );
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBase3DImp::add_rectangular_beam_with_velocity( uint32_t N, double J, double q, double m, 
								double v, double dvp, double dvt, Vec3D c, 
								Vec3D dir1, Vec3D dir2, 
								double size1, double size2 )
{
    ibsimu.message( 1 ) << "Defining a rectangular beam\n";

    _particles.reserve( _particles.size()+N );

    // Random number generator for two positions and three velocities (gaussian)
    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 5 );
    else
	rng = new MTRandom( 5 );
    rng->set_transformation( 2, Gaussian_Transformation() );
    rng->set_transformation( 3, Gaussian_Transformation() );
    rng->set_transformation( 4, Gaussian_Transformation() );
    double qx[5];
    double px[6];

    m *= MASS_U;
    q *= CHARGE_E;
    _rhosum += J/v;
    double IQ = 4.0*size1*size2*J/N;

    // Calculate and check base vectors
    Vec3D dir3 = cross( dir1, dir2 );
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( comp_isnan( dir1[0] ) || comp_isnan( dir2[0] ) || comp_isnan( dir3[0] ) )
	throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );

    // Prepare particle
    ParticleP3D x;
    x[0] = 0.0;

    double Isum = 0.0;
    uint32_t a = 0;
    while( a < N ) {

	rng->get( qx );

	// Position in natural (dir1,dir2,dir3) coordinates
	px[0] = size1*(2.0*qx[0]-1.0);
	px[1] = size2*(2.0*qx[1]-1.0);
	px[2] = 0.0;

	// Velocities in natural (dir1,dir2,dir3) coordinates
	px[3] = dvt*qx[2];
	px[4] = dvt*qx[3];
	px[5] = dvp*qx[4] + v;

	// Map to world coordinates
	x[1] = dir1[0]*px[0] + dir2[0]*px[1] + dir3[0]*px[2] + c[0];
	x[2] = dir1[0]*px[3] + dir2[0]*px[4] + dir3[0]*px[5];
	x[3] = dir1[1]*px[0] + dir2[1]*px[1] + dir3[1]*px[2] + c[1];
	x[4] = dir1[1]*px[3] + dir2[1]*px[4] + dir3[1]*px[5];
	x[5] = dir1[2]*px[0] + dir2[2]*px[1] + dir3[2]*px[2] + c[2];
	x[6] = dir1[2]*px[3] + dir2[2]*px[4] + dir3[2]*px[5];

	add_particle( Particle3D( IQ, q, m, x ) );
	Isum += IQ;
	a++;
    }

    delete( rng );
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBase3DImp::add_rectangular_beam_with_energy( uint32_t N, double J, double q, double m, 
							      double E, double Tp, double Tt, Vec3D c, 
							      Vec3D dir1, Vec3D dir2, 
							      double size1, double size2 )
{
    ibsimu.message( 1 ) << "Defining a rectangular beam\n";

    _particles.reserve( _particles.size()+N );

    // Convert input parameters
    m *= MASS_U;
    q *= CHARGE_E;
    E *= CHARGE_E;
    Tp *= CHARGE_E;
    Tt *= CHARGE_E;
    double v = energy_to_velocity(E,m);
    double dvp = sqrt(Tp/m);
    double dvt = sqrt(Tt/m);

    // Random number generator for two positions and three velocities (gaussian)
    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 5 );
    else
	rng = new MTRandom( 5 );
    rng->set_transformation( 2, Gaussian_Transformation() );
    rng->set_transformation( 3, Gaussian_Transformation() );
    rng->set_transformation( 4, Gaussian_Transformation() );
    double qx[5];
    double px[6];

    // Charge
    _rhosum += J/v;
    double IQ = 4.0*size1*size2*J/N;

    // Calculate and check base vectors
    Vec3D dir3 = cross( dir1, dir2 );
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( comp_isnan( dir1[0] ) || comp_isnan( dir2[0] ) || comp_isnan( dir3[0] ) )
	throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );

    // Prepare particle
    ParticleP3D x;
    x[0] = 0.0;

    double Isum = 0.0;
    uint32_t a = 0;
    while( a < N ) {

	rng->get( qx );

	// Calculate in natural (dir1,dir2,dir3) coordinates
	px[0] = size1*(2.0*qx[0]-1.0);
	px[1] = size2*(2.0*qx[1]-1.0);
	px[2] = 0.0;
	px[3] = dvt*qx[2];
	px[4] = dvt*qx[3];
	double pveld = dvp*qx[4];
        px[5] = sqrt( 2.0*E/m + pveld*pveld );

	// Map to world coordinates
	x[1] = dir1[0]*px[0] + dir2[0]*px[1] + dir3[0]*px[2] + c[0];
	x[2] = dir1[0]*px[3] + dir2[0]*px[4] + dir3[0]*px[5];
	x[3] = dir1[1]*px[0] + dir2[1]*px[1] + dir3[1]*px[2] + c[1];
	x[4] = dir1[1]*px[3] + dir2[1]*px[4] + dir3[1]*px[5];
	x[5] = dir1[2]*px[0] + dir2[2]*px[1] + dir3[2]*px[2] + c[2];
	x[6] = dir1[2]*px[3] + dir2[2]*px[4] + dir3[2]*px[5];

	add_particle( Particle3D( IQ, q, m, x ) );
	Isum += IQ;
	a++;
    }

    delete( rng );
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBase3DImp::add_3d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
							   double E0,
							   double a1, double b1, double e1,
							   double a2, double b2, double e2,
							   Vec3D c, Vec3D dir1, Vec3D dir2 )
{
    ibsimu.message( 1 ) << "Defining a 3d KV beam using Twiss parameters\n";

    _particles.reserve( _particles.size()+N );

    e1 *= 4.0; // KV-distributed emittance limited inside 4*e_rms ellipse
    e2 *= 4.0;

    m *= MASS_U;
    q *= CHARGE_E;

    // Calculate and check base vectors
    Vec3D dir3 = cross( dir1, dir2 );
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( comp_isnan( dir1[0] ) || comp_isnan( dir2[0] ) || comp_isnan( dir3[0] ) )
	throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );

    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 4 );
    else
	rng = new MTRandom( 4 );
    double rn[4];

    double g1 = (1.0 + a1*a1)/b1;
    double h1 = 0.5*(b1+g1);
    double rmaj1 = sqrt(0.5*e1)*(sqrt(h1+1.0)+sqrt(h1-1.0));
    double rmin1 = sqrt(0.5*e1)*(sqrt(h1+1.0)-sqrt(h1-1.0));
    double theta1 = 0.5*atan2( -2.0*a1, b1-g1 );

    double g2 = (1.0 + a2*a2)/b2;
    double h2 = 0.5*(b2+g2);
    double rmaj2 = sqrt(0.5*e2)*(sqrt(h2+1.0)+sqrt(h2-1.0));
    double rmin2 = sqrt(0.5*e2)*(sqrt(h2+1.0)-sqrt(h2-1.0));
    double theta2 = 0.5*atan2( -2.0*a2, b2-g2 );

    double IQ = I/N;

    ParticleP3D x, px;
    px[0] = x[0] = 0.0;
    px[5] = 0.0;
    px[6] = sqrt(2.0*E0*CHARGE_E/m);

    uint32_t n = 0;
    while( n < N ) {

	rng->get( rn );

	// Randomize point inside unit sphere
	double rrmaj1 = -1.0 + 2.0*rn[0];
	double rrmin1 = -1.0 + 2.0*rn[1];
	double rrmaj2 = -1.0 + 2.0*rn[2];
	double rrmin2 = -1.0 + 2.0*rn[3];

	// Check if inside sphere
	double radius2 = rrmaj1*rrmaj1 + rrmin1*rrmin1 + rrmaj2*rrmaj2 + rrmin2*rrmin2;
	if( radius2 == 0.0 || radius2 > 1.0 )
	    continue;

	// Project to surface
	double norm = 1.0/sqrt(radius2);
	rrmaj1 *= norm;
	rrmin1 *= norm;
	rrmaj2 *= norm;
	rrmin2 *= norm;

	// Scale to ellipsoid
	rrmaj1 *= rmaj1;
	rrmin1 *= rmin1;
	rrmaj2 *= rmaj2;
	rrmin2 *= rmin2;

	// Rotate to correct angle: p1, v1, p2, v2, where
	// v1 = vz*x1' and v2 = vz*x2'
	px[1] = rrmaj1*cos(theta1) - rrmin1*sin(theta1);
	px[2] = px[6]*(rrmaj1*sin(theta1) + rrmin1*cos(theta1));
	px[3] = rrmaj2*cos(theta2) - rrmin2*sin(theta2);
	px[4] = px[6]*(rrmaj2*sin(theta2) + rrmin2*cos(theta2));

	// Map to world coordinates x, vx, y, vy, z, vz 
	x[1] = dir1[0]*px[1] + dir2[0]*px[3] + dir3[0]*px[5] + c[0];
	x[2] = dir1[0]*px[2] + dir2[0]*px[4] + dir3[0]*px[6];
	x[3] = dir1[1]*px[1] + dir2[1]*px[3] + dir3[1]*px[5] + c[1];
	x[4] = dir1[1]*px[2] + dir2[1]*px[4] + dir3[1]*px[6];
	x[5] = dir1[2]*px[1] + dir2[2]*px[3] + dir3[2]*px[5] + c[2];
	x[6] = dir1[2]*px[2] + dir2[2]*px[4] + dir3[2]*px[6];

	add_particle( Particle3D( IQ, q, m, x ) );
	n++;
    }

    delete( rng );
}


void ParticleDataBase3DImp::add_3d_waterbag_beam_with_emittance( uint32_t N, double I, double q, double m,
								 double E0,
								 double a1, double b1, double e1,
								 double a2, double b2, double e2,
								 Vec3D c, Vec3D dir1, Vec3D dir2 )
{
    ibsimu.message( 1 ) << "Defining a 3d Waterbag beam using Twiss parameters\n";

    _particles.reserve( _particles.size()+N );

    e1 *= 6.0; // Waterbag distributed emittance limited inside 6*e_rms ellipse
    e2 *= 6.0;

    m *= MASS_U;
    q *= CHARGE_E;

    // Calculate and check base vectors
    Vec3D dir3 = cross( dir1, dir2 );
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( comp_isnan( dir1[0] ) || comp_isnan( dir2[0] ) || comp_isnan( dir3[0] ) )
	throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );

    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 4 );
    else
	rng = new MTRandom( 4 );
    double rn[4];

    double g1 = (1.0 + a1*a1)/b1;
    double h1 = 0.5*(b1+g1);
    double rmaj1 = sqrt(0.5*e1)*(sqrt(h1+1.0)+sqrt(h1-1.0));
    double rmin1 = sqrt(0.5*e1)*(sqrt(h1+1.0)-sqrt(h1-1.0));
    double theta1 = 0.5*atan2( -2.0*a1, b1-g1 );

    double g2 = (1.0 + a2*a2)/b2;
    double h2 = 0.5*(b2+g2);
    double rmaj2 = sqrt(0.5*e2)*(sqrt(h2+1.0)+sqrt(h2-1.0));
    double rmin2 = sqrt(0.5*e2)*(sqrt(h2+1.0)-sqrt(h2-1.0));
    double theta2 = 0.5*atan2( -2.0*a2, b2-g2 );

    double IQ = I/N;

    ParticleP3D x, px;
    px[0] = x[0] = 0.0;
    px[5] = 0.0;
    px[6] = sqrt(2.0*E0*CHARGE_E/m);

    uint32_t n = 0;
    while( n < N ) {

	rng->get( rn );

	// Randomize point inside unit sphere
	double rrmaj1 = -1.0 + 2.0*rn[0];
	double rrmin1 = -1.0 + 2.0*rn[1];
	double rrmaj2 = -1.0 + 2.0*rn[2];
	double rrmin2 = -1.0 + 2.0*rn[3];

	// Check if inside sphere
	double radius2 = rrmaj1*rrmaj1 + rrmin1*rrmin1 + rrmaj2*rrmaj2 + rrmin2*rrmin2;
	if( radius2 > 1.0 )
	    continue;

	// Scale to ellipsoid
	rrmaj1 *= rmaj1;
	rrmin1 *= rmin1;
	rrmaj2 *= rmaj2;
	rrmin2 *= rmin2;

	// Rotate to correct angle: p1, v1, p2, v2, where
	// v1 = vz*x1' and v2 = vz*x2'
	px[1] = rrmaj1*cos(theta1) - rrmin1*sin(theta1);
	px[2] = px[6]*(rrmaj1*sin(theta1) + rrmin1*cos(theta1));
	px[3] = rrmaj2*cos(theta2) - rrmin2*sin(theta2);
	px[4] = px[6]*(rrmaj2*sin(theta2) + rrmin2*cos(theta2));

	// Map to world coordinates x, vx, y, vy, z, vz 
	x[1] = dir1[0]*px[1] + dir2[0]*px[3] + dir3[0]*px[5] + c[0];
	x[2] = dir1[0]*px[2] + dir2[0]*px[4] + dir3[0]*px[6];
	x[3] = dir1[1]*px[1] + dir2[1]*px[3] + dir3[1]*px[5] + c[1];
	x[4] = dir1[1]*px[2] + dir2[1]*px[4] + dir3[1]*px[6];
	x[5] = dir1[2]*px[1] + dir2[2]*px[3] + dir3[2]*px[5] + c[2];
	x[6] = dir1[2]*px[2] + dir2[2]*px[4] + dir3[2]*px[6];

	add_particle( Particle3D( IQ, q, m, x ) );
	n++;
    }

    delete( rng );
}


void ParticleDataBase3DImp::add_3d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
								 double E0, 
								 double a1, double b1, double e1,
								 double a2, double b2, double e2,
								 Vec3D c, Vec3D dir1, Vec3D dir2 )
{
    ibsimu.message( 1 ) << "Defining a 3d gaussian beam using Twiss parameters\n";

    _particles.reserve( _particles.size()+N );

    m *= MASS_U;
    q *= CHARGE_E;

    // Calculate and check base vectors
    Vec3D dir3 = cross( dir1, dir2 );
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( comp_isnan( dir1[0] ) || comp_isnan( dir2[0] ) || comp_isnan( dir3[0] ) )
	throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );

    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 4 );
    else
	rng = new MTRandom( 4 );
    rng->set_transformation( 0, Gaussian_Transformation() );
    rng->set_transformation( 1, Gaussian_Transformation() );
    rng->set_transformation( 2, Gaussian_Transformation() );
    rng->set_transformation( 3, Gaussian_Transformation() );
    double w[4], rn[4];

    double g1 = (1.0 + a1*a1)/b1;
    double h1 = 0.5*(b1+g1);
    double rmaj1 = sqrt(0.5*e1)*(sqrt(h1+1.0)+sqrt(h1-1.0));
    double rmin1 = sqrt(0.5*e1)*(sqrt(h1+1.0)-sqrt(h1-1.0));
    double theta1 = 0.5*atan2( -2.0*a1, b1-g1 );

    double g2 = (1.0 + a2*a2)/b2;
    double h2 = 0.5*(b2+g2);
    double rmaj2 = sqrt(0.5*e2)*(sqrt(h2+1.0)+sqrt(h2-1.0));
    double rmin2 = sqrt(0.5*e2)*(sqrt(h2+1.0)-sqrt(h2-1.0));
    double theta2 = 0.5*atan2( -2.0*a2, b2-g2 );

    double IQ = I/N;

    ParticleP3D x, px;
    px[0] = x[0] = 0.0;
    px[5] = 0.0;
    px[6] = sqrt(2.0*E0*CHARGE_E/m);

    uint32_t n = 0;
    while( n < N ) {

	// Randomize point from gaussian distribution
	rng->get( rn );
	w[0] = rmaj1*rn[0];
	w[1] = rmin1*rn[1];
	w[2] = rmaj2*rn[2];
	w[3] = rmin2*rn[3];

	// Rotate to correct angle p1, v1, p2, v2, where
	// v1 = vz*x1' and v2 = vz*x2'
	px[1] = w[0]*cos(theta1) - w[1]*sin(theta1);
	px[2] = px[6]*(w[0]*sin(theta1) + w[1]*cos(theta1));
	px[3] = w[2]*cos(theta2) - w[3]*sin(theta2);
	px[4] = px[6]*(w[2]*sin(theta2) + w[3]*cos(theta2));

	// Map to world coordinates x, vx, y, vy, z, vz 
	x[1] = dir1[0]*px[1] + dir2[0]*px[3] + dir3[0]*px[5] + c[0];
	x[2] = dir1[0]*px[2] + dir2[0]*px[4] + dir3[0]*px[6];
	x[3] = dir1[1]*px[1] + dir2[1]*px[3] + dir3[1]*px[5] + c[1];
	x[4] = dir1[1]*px[2] + dir2[1]*px[4] + dir3[1]*px[6];
	x[5] = dir1[2]*px[1] + dir2[2]*px[3] + dir3[2]*px[5] + c[2];
	x[6] = dir1[2]*px[2] + dir2[2]*px[4] + dir3[2]*px[6];

	add_particle( Particle3D( IQ, q, m, x ) );
	n++;
    }

    delete( rng );
}


void ParticleDataBase3DImp::add_3d_plane_cosine_beam_with_partps( uint32_t N, double I, double m, double T, 
                          Vec3D c, Vec3D dir1, Vec3D dir2, double size1, double size2 )
{   //2016-11-21 This function should generate particles leaving a rectangular plane with a cosine distribution to facilitate neutral particle modeling.
    m *= MASS_U;           //The mass is converted from atomic mass units to kg.
    T *= KBOLTZ;           //The temperature is converted from K to J.
    double q = CHARGE_E;   //The particles carry a positive charge of 1 e.
    double IQ = I/N;       //Each particle carries a current IQ.
    double JQ = IQ/(4.0*size1*size2); // This is the current density per particle in A/m2. It is only used for the mysterious rhosum variable.
    // double dv = sqrt(T/m); //Each velocity component has a Gaussian distribution in a Maxwell-Boltzmann distribution.
    double dv = T/m;       //The square is needed so this is more convenient.

    ibsimu.message( 1 ) << "Defining a 3d rectangular beam with cosine distribution\n";

    _particles.reserve( _particles.size()+N );

    Random *rng; // Random number generator for two uniform random positions, and three Gaussian velocity components (to get a Maxwell-Boltzmann distribution for the total velocity).
    if( ibsimu.get_rng_type() == RNG_SOBOL )  //In addition, two more numbers are needed for the cosine distribution and the angles.
        rng = new QRandom( 7 );
    else
        rng = new MTRandom( 7 );
    rng->set_transformation( 0, Uniform_Transformation() ); //Ok, the first two components are default, this is just for clarity.
    rng->set_transformation( 1, Uniform_Transformation() );
    rng->set_transformation( 2, Gaussian_Transformation() );
    rng->set_transformation( 3, Gaussian_Transformation() );
    rng->set_transformation( 4, Gaussian_Transformation() );
    rng->set_transformation( 5, Uniform_Transformation() ); //Ok, the last two components are default, this is just for clarity.
    rng->set_transformation( 6, Uniform_Transformation() );
    double qx[7]; //Will hold the random numbers
    double px[6]; //Will hold the location and direction vectors in the frame of the particle
    double v;     //Will hold the total velocity

    rng->get( qx ); //The first numbers seem to be all zero from the rng, I use this useless line to initialize the rng. (NdH 2017-04-13)

    Vec3D dir3 = cross( dir1, dir2 ); // Calculate and check basis vectors
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( comp_isnan( dir1[0] ) || comp_isnan( dir2[0] ) || comp_isnan( dir3[0] ) )
        throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );

    ParticleP3D x; // Prepare particle (double t, double x, double vx, double y, double vy, double z, double vz)
    x[0] = 0.0;

    for( uint32_t a = 0; a < N; a++ ) //Loop over the particles to generate.
    {
        rng->get( qx );
        v = sqrt(qx[2]*qx[2]*dv + qx[3]*qx[3]*dv + qx[4]*qx[4]*dv); //Make Maxwell-Boltzmann distribution by sampling individual velocity components.

        px[0] = size1*(2.0*qx[0]-1.0); // Calculate in natural (dir1,dir2,dir3) coordinates
        px[1] = size2*(2.0*qx[1]-1.0);
        px[2] = 0.0;

        double sintheta = sqrt( qx[5] ); //These lines generate a cosine distribution.
        double costheta = sqrt( 1.0 - qx[5] );
        double psiangle = qx[6]*2.0*M_PI;         

        px[3] = v*sintheta*cos(psiangle); //Direction 3 is the direction normal to the plane!
        px[4] = v*sintheta*sin(psiangle);
        px[5] = v*costheta;               

        x[1] = dir1[0]*px[0] + dir2[0]*px[1] + dir3[0]*px[2] + c[0]; // Map to world coordinates
        x[2] = dir1[0]*px[3] + dir2[0]*px[4] + dir3[0]*px[5];
        x[3] = dir1[1]*px[0] + dir2[1]*px[1] + dir3[1]*px[2] + c[1];
        x[4] = dir1[1]*px[3] + dir2[1]*px[4] + dir3[1]*px[5];
        x[5] = dir1[2]*px[0] + dir2[2]*px[1] + dir3[2]*px[2] + c[2];
        x[6] = dir1[2]*px[3] + dir2[2]*px[4] + dir3[2]*px[5];

        add_particle( Particle3D( IQ, q, m, x ) );

        _rhosum += JQ/v;  //NOTE -> check later if this is correct!

    }

    delete( rng );

    ibsimu.message( 1 ) << "Defined a 3d rectangular beam with cosine distribution\n";
}

void ParticleDataBase3DImp::add_3d_cylindrical_cosine_beam_with_partps( uint32_t N, double I, double m, double T, 
                          Vec3D c, Vec3D dir1, Vec3D dir2, double r )
{   //2017-02-06 This function should generate particles leaving a circular plane with a cosine distribution to facilitate neutral particle modeling.
    m *= MASS_U;           //The mass is converted from atomic mass units to kg.
    T *= KBOLTZ;           //The temperature is converted from K to J.
    double q = CHARGE_E;   //The particles carry a positive charge of 1 e.
    double IQ = I/N;       //Each particle carries a current IQ.
    double JQ = IQ/(M_PI*r*r); // This is the current density per particle in A/m2. It is only used for the mysterious rhosum variable.
    // double dv = sqrt(T/m); //Each velocity component has a Gaussian distribution in a Maxwell-Boltzmann distribution.
    double dv = T/m;       //The square is needed so this is more convenient.

    ibsimu.message( 1 ) << "Defining a 3d cylindrical beam with cosine distribution\n";

    _particles.reserve( _particles.size()+N );

    Random *rng; // Random number generator for two uniform random positions, and three Gaussian velocity components (to get a Maxwell-Boltzmann distribution for the total velocity).
    if( ibsimu.get_rng_type() == RNG_SOBOL )  //In addition, two more numbers are needed for the cosine distribution and the angles.
        rng = new QRandom( 7 );
    else
        rng = new MTRandom( 7 );
    rng->set_transformation( 0, Uniform_Transformation() ); //Ok, the first two components are default, this is just for clarity.
    rng->set_transformation( 1, Uniform_Transformation() );
    rng->set_transformation( 2, Gaussian_Transformation() );
    rng->set_transformation( 3, Gaussian_Transformation() );
    rng->set_transformation( 4, Gaussian_Transformation() );
    rng->set_transformation( 5, Uniform_Transformation() ); //Ok, the last two components are default, this is just for clarity.
    rng->set_transformation( 6, Uniform_Transformation() );
    double qx[7]; //Will hold the random numbers
    double px[6]; //Will hold the location and direction vectors in the frame of the particle
    double v;     //Will hold the total velocity

    Vec3D dir3 = cross( dir1, dir2 ); // Calculate and check basis vectors
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( comp_isnan( dir1[0] ) || comp_isnan( dir2[0] ) || comp_isnan( dir3[0] ) )
        throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );

    ParticleP3D x; // Prepare particle (double t, double x, double vx, double y, double vy, double z, double vz)
    x[0] = 0.0;

    double Isum = 0.0;
    uint32_t a = 0;
    rng->get( qx ); //To initialize the random number generator, it is used once before, seems to be needed. 
    while( a < N ) //A while loop is used here because the particles that are outside the circle are skipped, and then a should not be increased.
    {
        rng->get( qx );

        px[0] = -r + 2.0*r*qx[0]; // Position in natural (dir1,dir2,dir3) coordinates
        px[1] = -r + 2.0*r*qx[1];
        px[2] = 0.0;

        if( px[0]*px[0] + px[1]*px[1] > r*r ) // Reject points outside radius r circle
            continue;

        v = sqrt(qx[2]*qx[2]*dv + qx[3]*qx[3]*dv + qx[4]*qx[4]*dv); //Make Maxwell-Boltzmann distribution by sampling individual velocity components.
        double sintheta = sqrt( qx[5] ); //These lines generate a cosine distribution.
        double costheta = sqrt( 1.0 - qx[5] );
        double psiangle = qx[6]*2.0*M_PI;         

        px[3] = v*sintheta*cos(psiangle); //Direction 3 is the direction normal to the plane!
        px[4] = v*sintheta*sin(psiangle);
        px[5] = v*costheta;       

        x[1] = dir1[0]*px[0] + dir2[0]*px[1] + dir3[0]*px[2] + c[0]; // Map to world coordinates
        x[2] = dir1[0]*px[3] + dir2[0]*px[4] + dir3[0]*px[5];
        x[3] = dir1[1]*px[0] + dir2[1]*px[1] + dir3[1]*px[2] + c[1];
        x[4] = dir1[1]*px[3] + dir2[1]*px[4] + dir3[1]*px[5];
        x[5] = dir1[2]*px[0] + dir2[2]*px[1] + dir3[2]*px[2] + c[2];
        x[6] = dir1[2]*px[3] + dir2[2]*px[4] + dir3[2]*px[5];

        add_particle( Particle3D( IQ, q, m, x ) );

        _rhosum += JQ/v;  //NOTE -> check later if this is correct!
        Isum += IQ;
        a++;
    }

    delete( rng );

    ibsimu.message( 1 ) << "Defined a 3d cylindrical beam with cosine distribution\n";
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";

}

void ParticleDataBase3DImp::add_3d_spherical_beam_with_energy( uint32_t N, double I, double q, double m, 
             double E, double Tp, double Tt, Vec3D cent, Vec3D pole, double angl)
{   //2018-06-21 This function should generate particles leaving a spherical surface with certain energy. Meaning of the input parameters:
    //uint32_t      Number of particles
    //double I      Total current carried by the particles (A)
    //double q      Charge of each particle (elementary charges)
    //double m      Mass of each particle (atomic mass units)
    //double E      Starting energy (eV)
    //double Tp     Temperature in parallel direction (eV)
    //double Tt     Temperature in transverse direction (eV)
    //Vec3D cent    Location of the center of the spere
    //Vec3D pole    Location of the middle of the emitting spherical cap
    //double angl   Polar angle between the rays from the center of the sphere to the apex of the cap (the pole) and the edge of the disk forming the base of the cap (degrees)

    // Convert input parameters
    m *= MASS_U;            //The mass is converted from atomic mass units to kg.
    q *= CHARGE_E;          //The particles carry charges that are a multiple of e.
    E *= CHARGE_E;
    Tp *= CHARGE_E;
    Tt *= CHARGE_E;
    angl *= (M_PI/180.0);
    double v = energy_to_velocity(E,m);
    double dvp = sqrt(Tp/m);
    double dvt = sqrt(Tt/m);
    Vec3D ndrx = pole - cent; //A vector pointing from the center to the pole
    double r = ndrx.norm2(); //Radius of sphere
    double A = 2.0*M_PI*r*r*(1.0 - cos(angl)); //Area of cap.
    double IQ = I/N;       //Each particle carries a current IQ.
    // double JQ = IQ/A;      //This is the current density per particle in A/m2. It is only used for the mysterious rhosum variable.

    ibsimu.message( 1 ) << "Defining a spherical beam\n";
    _particles.reserve( _particles.size()+N );
    _rhosum += (I/A)/v; //This is just copied from earlier routines, there it is J/V. J is the current density, or (I/A).

    ndrx.normalize(); //These vectors are the normal vectors of the sphere. The y and z direction are arbitrary, but this should be fine for a sphere. 
    Vec3D ndry = ndrx.arb_perpendicular(); //The x vector points from the center to the pole. 
    Vec3D ndrz = cross(ndry, ndrx);

    // Random number generator for two angles and three velocities (gaussian)
    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
    rng = new QRandom( 5 );
    else
    rng = new MTRandom( 5 );
    rng->set_transformation( 2, Gaussian_Transformation() );
    rng->set_transformation( 3, Gaussian_Transformation() );
    rng->set_transformation( 4, Gaussian_Transformation() );
    double qx[5];
    double vx[3];

    ParticleP3D x; // Prepare particle (double t, double x, double vx, double y, double vy, double z, double vz)
    x[0] = 0.0;

    double Isum = 0.0;
    uint32_t a = 0;
    rng->get( qx ); //To initialize the random number generator, it is used once before, seems to be needed. 
    while( a < N ) //Loop over the particles that are requested.
    {   //See https://math.stackexchange.com/questions/1585975/how-to-generate-random-points-on-a-sphere for an excellent explanation.
        rng->get( qx ); //The first two uniformly distributed numbers are used for the angle (location). First they are converted to latitude and longitude (on the full sphere)
        double lat = acos(2.0*qx[0] - 1.0) - 0.5*M_PI; //This angle is now between -pi/2 and +pi/2
        double lon = 2.0*M_PI*qx[1]; //This angle is between 0 and 2 pi.

        Vec3D plo(r*cos(lat)*cos(lon), r*cos(lat)*sin(lon), r*sin(lat)); //Position of point on sphere with radius r. 
        if( acos(plo[0]/r) > angl ) // Reject points outside spherical cap. Sounds inefficient, but compressing the whole sphere on the cap doesn't give a homogeneous distribution. 
        {
            continue; //This makes use of the fact that the points are generated in the x direction, the dot product of plo*(1,0,0) is just plo[0]. The length is always r.
        }
        else
        {
            x[1] = ndrx[0]*plo[0] + ndry[0]*plo[1] + ndrz[0]*plo[2] + cent[0]; //Transforms the basis directions. https://en.wikipedia.org/wiki/Transformation_matrix
            x[3] = ndrx[1]*plo[0] + ndry[1]*plo[1] + ndrz[1]*plo[2] + cent[1]; //Translate so that the center is at center. Now the position is made. Time to make the basis vectors of the random point.
            x[5] = ndrx[2]*plo[0] + ndry[2]*plo[1] + ndrz[2]*plo[2] + cent[2]; //The positions are now properly taken care of.

            Vec3D ntpx(cent[0] - x[1], cent[1] - x[3], cent[2] - x[5]); //Vector pointing from transformed position to the center.
            ntpx.normalize();
            Vec3D ntpy = ntpx.arb_perpendicular();
            Vec3D ntpz = cross(ntpy, ntpx);

            double pveld = dvp*qx[2]; //Velocities in natural (nx,ny,nz) coordinates. This is copied from other routines in the ParticleDataBase3DImp class. Note, nx is assumed to be the beam direction!
            vx[0] = sqrt( 2.0*E/m + pveld*pveld );
            vx[1] = dvt*qx[3];
            vx[2] = dvt*qx[4];

            x[2] = ntpx[0]*vx[0] + ntpy[0]*vx[1] + ntpz[0]*vx[2]; //Map velocity to world coordinates
            x[4] = ntpx[1]*vx[0] + ntpy[1]*vx[1] + ntpz[1]*vx[2];
            x[6] = ntpx[2]*vx[0] + ntpy[2]*vx[1] + ntpz[2]*vx[2];

            add_particle( Particle3D( IQ, q, m, x ) );
            Isum += IQ;
            a++;
        }
    }

    delete( rng );

    ibsimu.message( 1 ) << "Defined a 3d spherical beam with energy\n";
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";

}


void ParticleDataBase3DImp::trajectories_at_free_plane( TrajectoryDiagnosticData &tdata, 
							Vec3D c, Vec3D o, Vec3D p,
							const std::vector<trajectory_diagnostic_e> &diagnostics ) const
{
    // Check input
    for( size_t a = 0; a < diagnostics.size(); a++ ) {
	if( diagnostics[a] != DIAG_T && 
	    diagnostics[a] != DIAG_X && diagnostics[a] != DIAG_VX && 
	    diagnostics[a] != DIAG_Y && diagnostics[a] != DIAG_VY && 
	    diagnostics[a] != DIAG_Z && diagnostics[a] != DIAG_VZ && 
	    diagnostics[a] != DIAG_O && diagnostics[a] != DIAG_VO && 
	    diagnostics[a] != DIAG_P && diagnostics[a] != DIAG_VP && 
	    diagnostics[a] != DIAG_Q && diagnostics[a] != DIAG_VQ && 
	    diagnostics[a] != DIAG_OP && diagnostics[a] != DIAG_PP && 
	    diagnostics[a] != DIAG_CURR && diagnostics[a] != DIAG_EK && 
	    diagnostics[a] != DIAG_QM && diagnostics[a] != DIAG_CHARGE && 
	    diagnostics[a] != DIAG_MASS )
	    throw( Error( ERROR_LOCATION, "invalid diagnostic type " + to_string(diagnostics[a]) ) );
    }

    // Forward vector (normal to plane)
    Vec3D q = cross(o,p);
    if( q.norm2() == 0.0 )
	throw( Error( ERROR_LOCATION, "invalid vectors defining plane" ) );

    // Orthogonalize vector p
    p = cross(q,o);

    // Normalize vectors
    o.normalize();
    p.normalize();
    q.normalize();

    ibsimu.message( 1 ) << "Making trajectory diagnostics at plane\n" 
			<< "  c = " << c << "\n"
			<< "  o = " << o << "\n"
			<< "  p = " << p << "\n"
			<< "  q = " << q << "\n";

    // Prepare output vector
    tdata.clear();
    for( size_t a = 0; a < diagnostics.size(); a++ ) {
	tdata.add_data_column( diagnostics[a] );
    }

    // Scan through particle trajectory points
    double Isum = 0.0;
    for( size_t a = 0; a < _particles.size(); a++ ) {
	size_t N = _particles[a]->traj_size();
	if( N < 2 )
	    continue;
	ParticleP3D x1 = _particles[a]->traj(0);
	for( size_t b = 1; b < N; b++ ) {
	    ParticleP3D x2 = _particles[a]->traj(b);

	    // Solve trajectory crossing point using linear interpolation of position
	    // K[0]: parametric distance from x1 to x2
 	    // K[1]: distance in direction o
	    // K[2]: distance in direction p
	    Mat3D m( x2[1]-x1[1], -o[0], -p[0],
		     x2[3]-x1[3], -o[1], -p[1],
		     x2[5]-x1[5], -o[2], -p[2] );
	    double mdet = m.determinant();
	    if( mdet == 0.0 )
		continue;
	    Mat3D minv = m.inverse( mdet );
	    Vec3D off( c[0]-x1[1], c[1]-x1[3], c[2]-x1[5] );
	    Vec3D K = minv*off;
	    
	    if( K[0] >= 0.0 && K[0] < 1.0 ) {
		// Crossing found between x1 and x2, accumulate total current
		Isum += _particles[a]->IQ();

		// Position and velocity in xyz coordinates
		double t = K[0];
		double nt = 1.0-K[0];
		Vec3D pos( nt*x1[1]+t*x2[1], nt*x1[3]+t*x2[3], nt*x1[5]+t*x2[5] );
		Vec3D vel( nt*x1[2]+t*x2[2], nt*x1[4]+t*x2[4], nt*x1[6]+t*x2[6] );

		// Position and velocity in opq coordinates
		Vec3D pos_opq( K[1], K[2], 0.0 );
		Vec3D vel_opq( vel*o, vel*p, vel*q );

		// Fill diagnostic data
		for( size_t c = 0; c < tdata.diag_size(); c++ ) {

		    switch( tdata.diagnostic( c ) ) {
		    case DIAG_T:
			tdata.add_data( c, nt*x1[0]+t*x2[0] );
			break;

		    case DIAG_X:
			tdata.add_data( c, pos[0] );
			break;
		    case DIAG_VX:
			tdata.add_data( c, vel[0] );
			break;
		    case DIAG_Y:
			tdata.add_data( c, pos[1] );
			break;
		    case DIAG_VY:
			tdata.add_data( c, vel[1] );
			break;
		    case DIAG_Z:
			tdata.add_data( c, pos[2] );
			break;
		    case DIAG_VZ:
			tdata.add_data( c, vel[2] );
			break;

		    case DIAG_O:
			tdata.add_data( c, pos_opq[0] );
			break;
		    case DIAG_VO:
			tdata.add_data( c, vel_opq[0] );
			break;
		    case DIAG_P:
			tdata.add_data( c, pos_opq[1] );
			break;
		    case DIAG_VP:
			tdata.add_data( c, vel_opq[1] );
			break;
		    case DIAG_Q:
			tdata.add_data( c, pos_opq[2] );
			break;
		    case DIAG_VQ:
			tdata.add_data( c, vel_opq[2] );
			break;

		    case DIAG_OP:
			tdata.add_data( c, vel_opq[0]/vel_opq[2] );
			break;
		    case DIAG_PP:
			tdata.add_data( c, vel_opq[1]/vel_opq[2] );
			break;

		    case DIAG_CURR:
			tdata.add_data( c, _particles[a]->IQ() );
			break;
		    case DIAG_EK:
			if( vel.ssqr() < 1.0e-4*SPEED_C )
			    tdata.add_data( c, 0.5*_particles[a]->m()*vel.ssqr()/CHARGE_E );
			else {
			    double beta = vel.ssqr()/SPEED_C;
			    double gamma = 1.0 / sqrt( 1.0 - beta*beta );
			    double Ek = _particles[a]->m()*SPEED_C2*( gamma - 1.0 );
			    tdata.add_data( c, Ek/CHARGE_E );
			}
			break;
		    case DIAG_QM:
			tdata.add_data( c, (_particles[a]->q()/CHARGE_E) / (_particles[a]->m()/MASS_U) );
			break;
		    case DIAG_CHARGE:
			tdata.add_data( c, _particles[a]->q()/CHARGE_E );
			break;
		    case DIAG_MASS:
			tdata.add_data( c, _particles[a]->m()/MASS_U );
			break;
		    case DIAG_NO:
			tdata.add_data( c, a );
			break;
		    default:
			throw( ErrorUnimplemented( ERROR_LOCATION ) );
			break;
		    }
		}
	    }

	    // Next trajectory step
	    x1 = x2;
	}
    }

    ibsimu.message( 1 ) << "  number of trajectories = " << tdata.traj_size() << "\n";
    ibsimu.message( 1 ) << "  total current = " << Isum << " A\n";
}


/* Returns if axism type mirroring is enabled. The index axism
 * indicates mirroring across the planes in the list
 * (-,x,y,z,xy,xz,yz,xyz), starting from zero (no mirroring).
 */
bool ParticleDataBase3DImp::free_plane_mirror_enabled( uint32_t axism ) const
{
    switch( axism ) {
    case 1: // x
	return( _mirror[0] || _mirror[1] );
	break;
    case 2: // y
	return( _mirror[2] || _mirror[3] );
	break;
    case 3: // z
	return( _mirror[4] || _mirror[5] );
	break;
    case 4: // xy
	return( (_mirror[0] || _mirror[1]) && (_mirror[2] || _mirror[3]) );
	break;
    case 5: // xz
	return( (_mirror[0] || _mirror[1]) && (_mirror[4] || _mirror[5]) );
	break;
    case 6: // yz
	return( (_mirror[2] || _mirror[3]) && (_mirror[4] || _mirror[5]) );
	break;
    case 7: // xyz
	return( (_mirror[0] || _mirror[1]) && (_mirror[2] || _mirror[3]) && (_mirror[4] || _mirror[5]) );
	break;
    }

    // case 0: -
    return( true );
}


ParticleP3D ParticleDataBase3DImp::free_plane_mirror( const ParticleP3D &p, uint32_t axism ) const
{
    switch( axism ) {
    case 0: // -
	return( p );
    case 1: // x
	if( _mirror[0] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.origo(0)-p[1], -p[2], 
				 p[3], p[4], 
				 p[5], p[6] ) );
	else if( _mirror[1] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.max(0)-p[1], -p[2], 
				 p[3], p[4], 
				 p[5], p[6] ) );
	break;
    case 2: // y
	if( _mirror[2] )
	    return( ParticleP3D( p[0], 
				 p[1], p[2], 
				 2.0*_geom.origo(1)-p[3], -p[4], 
				 p[5], p[6] ) );
	else if( _mirror[3] )
	    return( ParticleP3D( p[0], 
				 p[1], p[2], 
				 2.0*_geom.max(1)-p[3], -p[4], 
				 p[5], p[6] ) );
	break;
    case 3: // z
	if( _mirror[2] )
	    return( ParticleP3D( p[0], 
				 p[1], p[2], 
				 p[3], p[4], 
				 2.0*_geom.origo(2)-p[5], -p[6] ) );
	else if( _mirror[3] )
	    return( ParticleP3D( p[0], 
				 p[1], p[2], 
				 p[3], p[4], 
				 2.0*_geom.max(2)-p[5], -p[6] ) );
	break;
    case 4: // xy
	if( _mirror[0] && _mirror[2] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.origo(0)-p[1], -p[2], 
				 2.0*_geom.origo(1)-p[3], -p[4], 
				 p[5], p[6] ) );
	else if( _mirror[0] && _mirror[3] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.origo(0)-p[1], -p[2], 
				 2.0*_geom.max(1)-p[3], -p[4], 
				 p[5], p[6] ) );
	else if( _mirror[1] && _mirror[2] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.max(0)-p[1], -p[2], 
				 2.0*_geom.origo(1)-p[3], -p[4], 
				 p[5], p[6] ) );
	else if( _mirror[1] && _mirror[3] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.max(0)-p[1], -p[2], 
				 2.0*_geom.max(1)-p[3], -p[4], 
				 p[5], p[6] ) );
	break;
    case 5: // xz
	if( _mirror[0] && _mirror[4] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.origo(0)-p[1], -p[2], 
				 p[3], p[4], 
				 2.0*_geom.origo(2)-p[5], -p[6] ) );
	else if( _mirror[0] && _mirror[5] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.origo(0)-p[1], -p[2], 
				 p[3], p[4], 
				 2.0*_geom.max(2)-p[5], -p[6] ) );
	else if( _mirror[1] && _mirror[4] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.max(0)-p[1], -p[2], 
				 p[3], p[4], 
				 2.0*_geom.origo(2)-p[5], -p[6] ) );
	else if( _mirror[1] && _mirror[5] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.max(0)-p[1], -p[2], 
				 p[3], p[4], 
				 2.0*_geom.max(2)-p[5], -p[6] ) );
	break;
    case 6: // yz
	if( _mirror[2] && _mirror[4] )
	    return( ParticleP3D( p[0], 
				 p[1], p[2], 
				 2.0*_geom.origo(1)-p[3], -p[4], 
				 2.0*_geom.origo(2)-p[5], -p[6] ) );
	else if( _mirror[2] && _mirror[5] )
	    return( ParticleP3D( p[0], 
				 p[1], p[2], 
				 2.0*_geom.origo(1)-p[3], -p[4], 
				 2.0*_geom.max(2)-p[5], -p[6] ) );
	else if( _mirror[3] && _mirror[4] )
	    return( ParticleP3D( p[0], 
				 p[1], p[2], 
				 2.0*_geom.max(1)-p[3], -p[4], 
				 2.0*_geom.origo(2)-p[5], -p[6] ) );
	else if( _mirror[3] && _mirror[5] )
	    return( ParticleP3D( p[0], 
				 p[1], p[2], 
				 2.0*_geom.max(1)-p[3], -p[4], 
				 2.0*_geom.max(2)-p[5], -p[6] ) );
	break;
    case 7: // xyz
	if( _mirror[0] && _mirror[2] && _mirror[4] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.origo(0)-p[1], -p[2], 
				 2.0*_geom.origo(1)-p[3], -p[4], 
				 2.0*_geom.origo(2)-p[5], -p[6] ) );
	else if( _mirror[0] && _mirror[2] && _mirror[5] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.origo(0)-p[1], -p[2], 
				 2.0*_geom.origo(1)-p[3], -p[4], 
				 2.0*_geom.max(2)-p[5], -p[6] ) );
	else if( _mirror[0] && _mirror[3] && _mirror[4] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.origo(0)-p[1], -p[2], 
				 2.0*_geom.max(1)-p[3], -p[4], 
				 2.0*_geom.origo(2)-p[5], -p[6] ) );
	else if( _mirror[0] && _mirror[3] && _mirror[5] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.origo(0)-p[1], -p[2], 
				 2.0*_geom.max(1)-p[3], -p[4], 
				 2.0*_geom.max(2)-p[5], -p[6] ) );
	else if( _mirror[1] && _mirror[2] && _mirror[4] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.max(0)-p[1], -p[2], 
				 2.0*_geom.origo(1)-p[3], -p[4], 
				 2.0*_geom.origo(2)-p[5], -p[6] ) );
	else if( _mirror[1] && _mirror[2] && _mirror[5] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.max(0)-p[1], -p[2], 
				 2.0*_geom.origo(1)-p[3], -p[4], 
				 2.0*_geom.max(2)-p[5], -p[6] ) );
	else if( _mirror[1] && _mirror[3] && _mirror[4] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.max(0)-p[1], -p[2], 
				 2.0*_geom.max(1)-p[3], -p[4], 
				 2.0*_geom.origo(2)-p[5], -p[6] ) );
	else if( _mirror[1] && _mirror[3] && _mirror[5] )
	    return( ParticleP3D( p[0], 
				 2.0*_geom.max(0)-p[1], -p[2], 
				 2.0*_geom.max(1)-p[3], -p[4], 
				 2.0*_geom.max(2)-p[5], -p[6] ) );
	break;
    }
    
    throw( ErrorAssert( ERROR_LOCATION ) );
}


void ParticleDataBase3DImp::export_path_manager_data( const std::string &filename, 
						      double ref_E, double ref_q, double ref_m, 
						      Vec3D c, Vec3D o, Vec3D p ) const
{
    // Forward vector (normal to plane)
    Vec3D q = cross(o,p);
    if( q.norm2() == 0.0 )
	throw( Error( ERROR_LOCATION, "invalid vectors defining plane" ) );

    // Orthogonalize vector p
    p = cross(q,o);

    // Normalize vectors
    o.normalize();
    p.normalize();
    q.normalize();

    // Calculate reference particle
    ref_m *= MASS_U;
    double ref_v = sqrt(2.0*ref_E*CHARGE_E/ref_m);
    double ref_p = ref_m*ref_v;
    const double m_to_gevc2 = SPEED_C2/CHARGE_E/1.0e9;
    const double p_to_gevc = SPEED_C/CHARGE_E/1.0e9;

    ibsimu.message( 1 ) << "Making trajectory diagnostics at plane for Path Manager output\n" 
			<< "  ref_E = " << ref_E << " eV\n"
			<< "  ref_q = " << ref_q << " e\n"
			<< "  ref_m = " << ref_m/MASS_U << " u\n"
			<< "  c = " << c << "\n"
			<< "  o = " << o << "\n"
			<< "  p = " << p << "\n"
			<< "  q = " << q << "\n";
    ibsimu.message( 1 ) << "  Outputting to \'" << filename << "\'\n";

    std::ofstream ofile( filename.c_str() );
    if( !ofile.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\'" ) );

    // Scan through particle trajectory points and build data
    std::vector<double> data[7]; // x, x', y, y', dp/p, q, m
    bool first = true;
    double IQ_first = 0.0;
    bool IQ_warning = false;
    double Isum = 0.0;
    for( size_t a = 0; a < _particles.size(); a++ ) {

	double pq = _particles[a]->q();
	double pm = _particles[a]->m();
	size_t N = _particles[a]->traj_size();
	if( N < 2 )
	    continue;

	// Mirroring of trajectory (-,x,y,z,xy,xz,yz,xyz)
	for( uint32_t axism = 0; axism < 8; axism++ ) {
	    if( !free_plane_mirror_enabled( axism ) )
		continue;

	    ParticleP3D x1 = free_plane_mirror( _particles[a]->traj(0), axism );
	    for( size_t b = 1; b < N; b++ ) {
		ParticleP3D x2 = free_plane_mirror( _particles[a]->traj(b), axism );

		// Solve trajectory crossing point using linear interpolation of position
		// K[0]: parametric distance from x1 to x2
		// K[1]: distance in direction o
		// K[2]: distance in direction p
		Vec3D K;
		try {
		    Mat3D m( x2[1]-x1[1], -o[0], -p[0],
			     x2[3]-x1[3], -o[1], -p[1],
			     x2[5]-x1[5], -o[2], -p[2] );
		    Mat3D minv = m.inverse();
		    Vec3D off( c[0]-x1[1], c[1]-x1[3], c[2]-x1[5] );
		    K = minv*off;
		} catch(...) {
		    continue;
		}
		
		if( K[0] >= 0.0 && K[0] < 1.0 ) {
		    // Crossing found between x1 and x2, accumulate total current
		    Isum += _particles[a]->IQ();
		    
		    // Position and velocity in xyz coordinates
		    double t = K[0];
		    double nt = 1.0-K[0];
		    Vec3D pos( nt*x1[1]+t*x2[1], nt*x1[3]+t*x2[3], nt*x1[5]+t*x2[5] );
		    Vec3D vel( nt*x1[2]+t*x2[2], nt*x1[4]+t*x2[4], nt*x1[6]+t*x2[6] );
		    
		    // Position and velocity in opq coordinates
		    Vec3D pos_opq( K[1], K[2], 0.0 );
		    Vec3D vel_opq( vel*o, vel*p, vel*q );
		    
		    // Check IQ fixed
		    if( first ) {
			IQ_first = _particles[a]->IQ();
		    } else if( IQ_first == 0.0 ) {
			if( !IQ_warning && _particles[a]->IQ() != 0.0 ) {
			    ibsimu.message( 1 ) << "  WARNING: Trajectory current is not constant.\n";
			    IQ_warning = true;
			}
		    } else {
			if( !IQ_warning && fabs( ( IQ_first - _particles[a]->IQ() ) / IQ_first ) > 1.0e-6 ) {
			    ibsimu.message( 1 ) << "  WARNING: Trajectory current is not constant.\n";
			    IQ_warning = true;
			}
		    }

		    // Save point
		    data[0].push_back( pos_opq[0] );                  // x (m)
		    data[1].push_back( vel_opq[0]/vel_opq[2] );       // x' (rad)
		    data[2].push_back( pos_opq[1] );                  // y (m)
		    data[3].push_back( vel_opq[1]/vel_opq[2] );       // y' (rad)
		    data[4].push_back( (pm*vel_opq[2]-ref_p)/ref_p ); // (p_z-p_ref)/p_ref
		    data[5].push_back( pq/CHARGE_E );                 // q (in e)
		    data[6].push_back( pm*m_to_gevc2 );               // m (in GeV/c2)
		    
		    first = false;
		}

		// Next trajectory step
		x1 = x2;
	    }
	}
    }

    // Write header
    time_t tt = time(NULL);
    char *timebuf = ctime( &tt );
    ofile << "Beam data from IBSimu, " << timebuf; // contains newline
    ofile << "Total beam: " << Isum << " A\n";
    ofile << std::setw(17) << ref_p*p_to_gevc  << " !REFERENCE MOMENTUM [GeV/c]\n";
    ofile << std::setw(17) << 0.0              << " !REFERENCE PHASE [rad]\n";
    ofile << std::setw(17) << 0.0              << " !FREQUENCY [Hz]\n";
    ofile << std::setw(17) << ref_m*m_to_gevc2 << " !REFERENCE MASS [GeV/c2]\n";
    ofile << std::setw(17) << ref_q            << " !REFERENCE CHARGE STATE\n";
    ofile << std::setw(17) << data[0].size()   << " !NUMBER OF PARTICLES\n";

    for( size_t a = 0; a < data[0].size(); a++ ) {	
	ofile << std::setw(17) << a << " ";       // particle number
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[0][a] << " "; // x
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[1][a] << " "; // x'
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[2][a] << " "; // y
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[3][a] << " "; // y'
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << -M_PI + 2.0*M_PI*rand() / RAND_MAX << " "; // phase
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[4][a] << " "; // dp/p
	ofile << std::setw(17) << 0 << " "; // flag (alive=0)
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[5][a] << " "; // q
	ofile << std::setw(17) << std::scientific << std::setprecision(9) << data[6][a] << " "; // m
	ofile << "\n";
    }

    ibsimu.message( 1 ) << "  number of trajectories = " << data[0].size() << "\n";
    ibsimu.message( 1 ) << "  total current = " << Isum << " A\n";
}


void ParticleDataBase3DImp::add_tdens_from_segment( MeshScalarField &tdens, double IQ,
						    ParticleP3D &x1, ParticleP3D &x2 ) const
    
{
    double dx = 0.0;
    double x[3];
    double t[3];
    int32_t i[3];

    for( size_t a = 0; a < 3; a++ ) {
	double xx = ( x2[2*a+1] - x1[2*a+1] );
	dx  += xx*xx;
	x[a] = 0.5*( x1[2*a+1] + x2[2*a+1] );
	i[a] = (int32_t)floor( ( x[a]-tdens.origo(a) ) * tdens.div_h() );
	t[a] = ( x[a]-(i[a]*tdens.h()+tdens.origo(a)) ) * tdens.div_h();
	
	if( i[a] < 0 || i[a] >= (int32_t)(tdens.size(a)-1) )
	    return;
    }
    dx = sqrt( dx );

    double J = IQ*dx; // A*m
    int p = tdens.size(0)*tdens.size(1)*i[2] + tdens.size(0)*i[1] + i[0];
    tdens( p )                 += (1.0-t[0])*(1.0-t[1])*(1.0-t[2])*J;
    tdens( p+tdens.size(0) )   += (1.0-t[0])*t[1]*(1.0-t[2])*J;
    tdens( p+1 )               += t[0]*(1.0-t[1])*(1.0-t[2])*J;
    tdens( p+1+tdens.size(0) ) += t[0]*t[1]*(1.0-t[2])*J;
    
    p += tdens.size(0)*tdens.size(1);
    tdens( p )                 += (1.0-t[0])*(1.0-t[1])*t[2]*J;
    tdens( p+tdens.size(0) )   += (1.0-t[0])*t[1]*t[2]*J;
    tdens( p+1 )               += t[0]*(1.0-t[1])*t[2]*J;
    tdens( p+1+tdens.size(0) ) += t[0]*t[1]*t[2]*J;
}


void ParticleDataBase3DImp::build_trajectory_density_field( MeshScalarField &tdens ) const
{
    tdens.clear();

    // Go through all particle trajectories
    std::vector< ColData<ParticleP3D> > coldata;
    for( size_t a = 0; a < _particles.size(); a++ ) {
	size_t N = _particles[a]->traj_size();
	if( N < 2 )
	    continue;
	ParticleP3D x1 = _particles[a]->traj(0);
	ParticleP3D xlast = x1;
	double IQ = _particles[a]->IQ();
	for( size_t b = 1; b < N; b++ ) {
	    ParticleP3D x2 = _particles[a]->traj(b);

	    // Mesh collision points from x1 to x2
	    ColData<ParticleP3D>::build_coldata_linear( coldata, tdens, x1, x2 );
	    if( coldata.size() == 0 )
		continue;

	    // Process trajectory mesh collision points
	    for( size_t c = 0; c < coldata.size(); c++ ) {

		// Add contribution from segment from xlast to coldata[c]
		add_tdens_from_segment( tdens, IQ, xlast, coldata[c]._x );	    
		xlast = coldata[c]._x;
	    }

	    // Next trajectory line
	    x1 = x2;
	}
    }    

    // Normalize to trajectory density, A/m2
    tdens /= (tdens.h()*tdens.h()*tdens.h());

    // Fix boundaries, which only get half of the contribution they should
    // This should depend on mirroring properties!
    for( uint32_t i = 0; i < tdens.size(0); i++ ) {
	for( uint32_t j = 0; j < tdens.size(1); j++ ) {
	    tdens( i, j, 0 ) *= 2.0;
	    tdens( i, j, tdens.size(2)-1 ) *= 2.0;
	}
    }
    for( uint32_t i = 0; i < tdens.size(0); i++ ) {
	for( uint32_t k = 0; k < tdens.size(2); k++ ) {
	    tdens( i, 0, k ) *= 2.0;
	    tdens( i, tdens.size(1)-1, k ) *= 2.0;
	}
    }	
    for( uint32_t j = 0; j < tdens.size(1); j++ ) {
	for( uint32_t k = 0; k < tdens.size(2); k++ ) {
	    tdens( 0, j, k ) *= 2.0;
	    tdens( tdens.size(0)-1, j, k ) *= 2.0;
	}
    }
}


void ParticleDataBase3DImp::save( std::ostream &os ) const
{
    ParticleDataBasePPImp<ParticleP3D>::save( os );
}


void ParticleDataBase3DImp::debug_print( std::ostream &os ) const 
{
    os << "**ParticleDataBase3D\n";
    ParticleDataBasePPImp<ParticleP3D>::debug_print( os );
}



