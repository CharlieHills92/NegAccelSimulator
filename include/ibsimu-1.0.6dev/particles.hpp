/*! \file particles.hpp
 *  \brief %Particle and particle point objects
 */

/* Copyright (c) 2005-2012,2015 Taneli Kalvas. All rights reserved.
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

#ifndef PARTICLES_HPP
#define PARTICLES_HPP 1


#include <vector>
#include <string>
#include <string.h>
#include <gsl/gsl_errno.h>
#include "geometry.hpp"
#include "meshscalarfield.hpp"
#include "vectorfield.hpp"
#include "vec3d.hpp"
#include "callback.hpp"
#include "constants.hpp"


/* Integer error value that is supposed to differ from internal GSL
 * error values */
#define IBSIMU_DERIV_ERROR 201


/*! \brief %Particle status enum.
 *
 *  %Particle status can be either OK (\a PARTICLE_OK), out of geometry
 *  (\a PARTICLE_OUT), collided (\a PARTICLE_COLL), badly defined (\a
 *  PARTICLE_BADDEF), simulated to time limit (\a PARTICLE_TIME) or
 *  simulated to step limit (\a PARTICLE_NSTP).
 */
enum particle_status_e {
    PARTICLE_OK = 0,
    PARTICLE_OUT,
    PARTICLE_COLL,
    PARTICLE_BADDEF,
    PARTICLE_TIME,
    PARTICLE_NSTP,
    PARTICLE_STRIP,
    PARTICLE_PERIODIC
};



/* ************************************************************************************* *
 * Particle point classes                                                                *
 * ************************************************************************************* */


/*! \brief %Particle point base class
 *
 *  A parent class for a common handle to %Particle point classes.
 */
class ParticlePBase
{

};


/*! \brief %Particle point class for 2D.
 *
 *  A collection of particle coordinates (t,x,vx,y,vy) to describe
 *  time, location and velocity in 2D.
 */
class ParticleP2D : public ParticlePBase 
{
    double         _x[5];   /*!< \brief Point coordinates */

public:

    /*! \brief Default constuctor.
     */
    ParticleP2D() {}

    /*! \brief Constructor for 2D particle point.
     */
    ParticleP2D( double t, double x, double vx, double y, double vy ) {
	_x[0] = t; _x[1] = x; _x[2] = vx; _x[3] = y; _x[4] = vy;
    }

    /*! \brief Constructor for loading particle point from a file.
     */
    ParticleP2D( std::istream &s ) {
	_x[0] = read_double( s );
	_x[1] = read_double( s );
	_x[2] = read_double( s );
	_x[3] = read_double( s );
	_x[4] = read_double( s );
    }

    /*! \brief Returns geometry mode.
     */
    static geom_mode_e geom_mode() { return(MODE_2D); }

    /*! \brief Returns number of dimensions for geometry.
     */
    static size_t dim() { return(2); }

    /*! \brief Returns number of coordinates used for particle point.
     */
    static size_t size() { return(5); }

    /*! \brief Returns time derivatives \a dxdt of coordinates at time
     *  \a t and coordinates \a x = (x,vx,y,vy) for one particle.
     *
     *  The calculation of particle trajectory is done by integrating
     *  the Lorentz equation in a form of a set of ordinary
     *  differential equations. In the case of 2D the set is:
     *  \f{eqnarray*} \frac{dx}{dt}   &=& v_x \\ 
                      \frac{dy}{dt}   &=& v_y \\ 
                      \frac{dv_x}{dt} &=& a_x = \frac{q}{m}(E_x + v_y B_z) \\ 
                      \frac{dv_y}{dt} &=& a_y = \frac{q}{m}(E_y - v_x B_z) \f}
     */
    static int get_derivatives( double t, const double *x, double *dxdt, void *data );

    /*! \brief Return the number of trajectory intersections with
     *  plane \a crd = val on the trajectory from \a x1 to \a
     *  x2. Intersection points are appended to vector \a intsc.
     */
    static int trajectory_intersections_at_plane( std::vector<ParticleP2D> &intsc, 
						  int crd, double val,
						  const ParticleP2D &x1, 
						  const ParticleP2D &x2,
						  int extrapolate = 0 );

    /*! \brief Return string representation for unit of current
     *
     *  Returns "A/m".
     */
    static const std::string IQ_unit() { return( "A/m" ); }

    /*! \brief Clear point to zero.
     */
    void clear() { memset( (void *)_x, 0, 5*sizeof(double) ); }

    /*! \brief Returns the location of particle point in Vec3D.
     */
    Vec3D location() const { return( Vec3D( _x[1], _x[3], 0.0 ) ); }

    /*! \brief Returns the velocity of particle point in Vec3D.
     */
    Vec3D velocity() const { return( Vec3D( _x[2], _x[4], 0.0 ) ); }

    /*! \brief Returns speed of particle.
     */
    double speed() { return( sqrt(_x[2]*_x[2] + _x[4]*_x[4]) ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator[]( int i ) { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator[]( int i ) const { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator()( int i ) { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator()( int i ) const { return( _x[i] ); }

    ParticleP2D operator+( const ParticleP2D &pp ) const { 
	ParticleP2D res;
	res[0] = _x[0] + pp[0];
	res[1] = _x[1] + pp[1];
	res[2] = _x[2] + pp[2];
	res[3] = _x[3] + pp[3];
	res[4] = _x[4] + pp[4];
	return( res );
    }

    ParticleP2D operator-( const ParticleP2D &pp ) const { 
	ParticleP2D res;
	res[0] = _x[0] - pp[0];
	res[1] = _x[1] - pp[1];
	res[2] = _x[2] - pp[2];
	res[3] = _x[3] - pp[3];
	res[4] = _x[4] - pp[4];
	return( res );
    }

    ParticleP2D operator*( double x ) const { 
	ParticleP2D res;
	res[0] = _x[0]*x;
	res[1] = _x[1]*x;
	res[2] = _x[2]*x;
	res[3] = _x[3]*x;
	res[4] = _x[4]*x;
	return( res );
    }

    /*! \brief Saves data to stream.
     */
    void save( std::ostream &s ) const {
	write_double( s, _x[0] );
	write_double( s, _x[1] );
	write_double( s, _x[2] );
	write_double( s, _x[3] );
	write_double( s, _x[4] );
    }
    
    friend ParticleP2D operator*( double x, const ParticleP2D &pp );
};


inline std::ostream &operator<<( std::ostream &os, const ParticleP2D &pp )
{
    os << "("
       << std::setw(12) << pp(0) << ", "
       << std::setw(12) << pp(1) << ", "
       << std::setw(12) << pp(2) << ", "
       << std::setw(12) << pp(3) << ", "
       << std::setw(12) << pp(4) << ")";
    return( os );
}


inline ParticleP2D operator*( double x, const ParticleP2D &pp )
{
    ParticleP2D res;
    res[0] = pp[0]*x;
    res[1] = pp[1]*x;
    res[2] = pp[2]*x;
    res[3] = pp[3]*x;
    res[4] = pp[4]*x;
    return( res );
}


/*! \brief %Particle point class for cylindrical coordinates.
 *
 *  A collection of particle coordinates (t,x,vx,r,vr,w) to
 *  describe time, location and velocity in cylindrical coordinates Here \f$ w = \frac{d\theta}{dt} \f$.
 */
class ParticlePCyl : public ParticlePBase 
{
    double         _x[6];   /*!< \brief Point coordinates */

public:

    /*! \brief Default constuctor.
     */
    ParticlePCyl() {}

    /*! \brief Constructor for cylindrical particle point.
     */
    ParticlePCyl( double t, double x, double vx, double r, double vr, double w ) {
	_x[0] = t; _x[1] = x; _x[2] = vx; _x[3] = r; _x[4] = vr; _x[5] = w;
    }

    /*! \brief Constructor for loading particle point from a file.
     */
    ParticlePCyl( std::istream &s ) {
	_x[0] = read_double( s );
	_x[1] = read_double( s );
	_x[2] = read_double( s );
	_x[3] = read_double( s );
	_x[4] = read_double( s );
	_x[5] = read_double( s );
    }

    /*! \brief Returns geometry mode.
     */
    static geom_mode_e geom_mode() { return(MODE_CYL); }

    /*! \brief Returns number of dimensions for geometry.
     */
    static size_t dim() { return(2); }

    /*! \brief Returns number of coordinates used for particle point.
     */
    static size_t size() { return(6); }

    /*! \brief Returns time derivatives \a dxdt of coordinates at time
     *  \a t and coordinates \a x = (x,vx,r,vr,w) for one particle.
     *
     *  The calculation of particle trajectory is done by integrating
     *  the Lorentz equation in a form of a set of ordinary
     *  differential equations. In the case of cylindrical coordinates 
     *  the set is:
     *  \f{eqnarray*} \frac{dx}{dt}   &=& v_x \\ 
                      \frac{dr}{dt}   &=& v_r \\ 
                      \frac{dv_x}{dt} &=& a_x = \frac{q}{m}(E_x + v_r B_{\theta} - v_\theta B_r) \\ 
                      \frac{dv_r}{dt} &=& a_r + r\left(\frac{d\theta}{dt}\right)^2 = 
                      \frac{q}{m}(E_y + v_{\theta} B_x - v_x B_{\theta}) + r \left(\frac{d\theta}{dt}\right)^2 \\
                      \frac{d^2\theta}{dt^2} &=& \frac{1}{r} \left(a_\theta-\frac{dr}{dt} \frac{d\theta}{dt}\right) =
		      \frac{1}{r} \left(\frac{q}{m}(v_x B_r - v_r B_x)-2 \frac{dr}{dt} \frac{d\theta}{dt}\right), \f}
     *  where \f$ v_\theta = r \frac{d\theta}{dt} \f$.
     */
    static int get_derivatives( double t, const double *x, double *dxdt, void *data );

    /*! \brief Return the number of trajectory intersections with
     *  plane \a crd = val on the trajectory from \a x1 to \a
     *  x2. Intersection points are appended to vector \a intsc.
     */
    static int trajectory_intersections_at_plane( std::vector<ParticlePCyl> &intsc, 
						  int crd, double val,
						  const ParticlePCyl &x1,
						  const ParticlePCyl &x2,
						  int extrapolate = 0 );

    /*! \brief Return string representation for unit of current
     *
     *  Returns "A".
     */
    static const std::string IQ_unit() { return( "A" ); }

    /*! \brief Clear point to zero.
     */
    void clear() { memset( (void *)_x, 0, 6*sizeof(double) ); }

    /*! \brief Returns the location of particle point in Vec3D.
     */
    Vec3D location() const { return( Vec3D( _x[1], _x[3], 0.0 ) ); }

    /*! \brief Returns the velocity of particle point in Vec3D.
     */
    Vec3D velocity() const { return( Vec3D( _x[2], _x[4], _x[5]*_x[3] ) ); }

    /*! \brief Returns speed of particle.
     */
    double speed() { return( sqrt(_x[2]*_x[2] + _x[4]*_x[4] + _x[3]*_x[3]*_x[5]*_x[5]) ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator[]( int i ) { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator[]( int i ) const { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator()( int i ) { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator()( int i ) const { return( _x[i] ); }

    ParticlePCyl operator+( const ParticlePCyl &pp ) const { 
	ParticlePCyl res;
	res[0] = _x[0] + pp[0];
	res[1] = _x[1] + pp[1];
	res[2] = _x[2] + pp[2];
	res[3] = _x[3] + pp[3];
	res[4] = _x[4] + pp[4];
	res[5] = _x[5] + pp[5];
	return( res );
    }

    ParticlePCyl operator-( const ParticlePCyl &pp ) const { 
	ParticlePCyl res;
	res[0] = _x[0] - pp[0];
	res[1] = _x[1] - pp[1];
	res[2] = _x[2] - pp[2];
	res[3] = _x[3] - pp[3];
	res[4] = _x[4] - pp[4];
	res[5] = _x[5] - pp[5];
	return( res );
    }

    ParticlePCyl operator*( double x ) const { 
	ParticlePCyl res;
	res[0] = _x[0]*x;
	res[1] = _x[1]*x;
	res[2] = _x[2]*x;
	res[3] = _x[3]*x;
	res[4] = _x[4]*x;
	res[5] = _x[5]*x;
	return( res );
    }

    /*! \brief Saves data to stream.
     */
    void save( std::ostream &s ) const {
	write_double( s, _x[0] );
	write_double( s, _x[1] );
	write_double( s, _x[2] );
	write_double( s, _x[3] );
	write_double( s, _x[4] );
	write_double( s, _x[5] );
    }
    
    friend ParticlePCyl operator*( double x, const ParticlePCyl &pp );
};


inline std::ostream &operator<<( std::ostream &os, const ParticlePCyl &pp )
{
    os << "("
       << std::setw(12) << pp(0) << ", "
       << std::setw(12) << pp(1) << ", "
       << std::setw(12) << pp(2) << ", "
       << std::setw(12) << pp(3) << ", "
       << std::setw(12) << pp(4) << ", "
       << std::setw(12) << pp(5) << ")";
    return( os );
}


inline ParticlePCyl operator*( double x, const ParticlePCyl &pp )
{
    ParticlePCyl res;
    res[0] = pp[0]*x;
    res[1] = pp[1]*x;
    res[2] = pp[2]*x;
    res[3] = pp[3]*x;
    res[4] = pp[4]*x;
    res[5] = pp[5]*x;
    return( res );
}


/*! \brief %Particle point class for 3D.
 *
 *  A collection of particle coordinates (t,x,vx,y,vy,z,vz)to describe
 *  time, location and velocity in 3D.
 */
class ParticleP3D : public ParticlePBase 
{
    double         _x[7];   /*!< \brief Point coordinates */

public:

    /*! \brief Default constuctor.
     */
    ParticleP3D() {}

    /*! \brief Constructor for 3D particle point.
     */
    ParticleP3D( double t, double x, double vx, double y, double vy, double z, double vz ) {
	_x[0] = t; _x[1] = x; _x[2] = vx; _x[3] = y; _x[4] = vy; _x[5] = z; _x[6] = vz;
    }

    /*! \brief Constructor for loading particle point from a file.
     */
    ParticleP3D( std::istream &s ) {
	_x[0] = read_double( s );
	_x[1] = read_double( s );
	_x[2] = read_double( s );
	_x[3] = read_double( s );
	_x[4] = read_double( s );
	_x[5] = read_double( s );
	_x[6] = read_double( s );
    }

    /*! \brief Returns geometry mode.
     */
    static geom_mode_e geom_mode() { return(MODE_3D); }

    /*! \brief Returns number of dimensions for geometry.
     */
    static size_t dim() { return(3); }

    /*! \brief Returns number of coordinates used for particle point.
     */
    static size_t size() { return(7); }

    /*! \brief Returns time derivatives \a dxdt of coordinates at time
     *  \a t and coordinates \a x = (x,vx,y,vy,z,vz) for one particle.
     *
     *  The calculation of particle trajectory is done by integrating
     *  the Lorentz equation in a form of a set of ordinary
     *  differential equations. In the case of 3D the set is:
     *  \f{eqnarray*} \frac{dx}{dt}   &=& v_x \\ 
                      \frac{dy}{dt}   &=& v_y \\ 
                      \frac{dz}{dt}   &=& v_z \\ 
                      \frac{dv_x}{dt} &=& a_x = \frac{q}{m}(E_x + v_y B_z - v_z B_y) \\ 
                      \frac{dv_y}{dt} &=& a_y = \frac{q}{m}(E_y + v_z B_x - v_x B_z) \\
                      \frac{dv_z}{dt} &=& a_z = \frac{q}{m}(E_z + v_x B_y - v_y B_x) \f}
     */
    static int get_derivatives( double t, const double *x, double *dxdt, void *data );

    /*! \brief Return the number of trajectory intersections with
     *  plane \a crd = val on the trajectory from \a x1 to \a
     *  x2. Intersection points are appended to vector \a intsc.
     */
    static int trajectory_intersections_at_plane( std::vector<ParticleP3D> &intsc, 
						  int crd, double val,
						  const ParticleP3D &x1, 
						  const ParticleP3D &x2, 
						  int extrapolate = 0 );

    /*! \brief Return string representation for unit of current
     *
     *  Returns "A".
     */
    static const std::string IQ_unit() { return( "A" ); }

    /*! \brief Clear point to zero.
     */
    void clear() { memset( (void *)_x, 0, 7*sizeof(double) ); }

    /*! \brief Returns the location of particle point in Vec3D.
     */
    Vec3D location() const { return( Vec3D( _x[1], _x[3], _x[5] ) ); }

    /*! \brief Returns the velocity of particle point in Vec3D.
     */
    Vec3D velocity() const { return( Vec3D( _x[2], _x[4], _x[6] ) ); }

    /*! \brief Returns speed of particle.
     */
    double speed() { return( sqrt(_x[2]*_x[2] + _x[4]*_x[4] + _x[6]*_x[6]) ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator[]( int i ) { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator[]( int i ) const { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator()( int i ) { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator()( int i ) const { return( _x[i] ); }

    ParticleP3D operator+( const ParticleP3D &pp ) const { 
	ParticleP3D res;
	res[0] = _x[0] + pp[0];
	res[1] = _x[1] + pp[1];
	res[2] = _x[2] + pp[2];
	res[3] = _x[3] + pp[3];
	res[4] = _x[4] + pp[4];
	res[5] = _x[5] + pp[5];
	res[6] = _x[6] + pp[6];
	return( res );
    }

    ParticleP3D operator-( const ParticleP3D &pp ) const { 
	ParticleP3D res;
	res[0] = _x[0] - pp[0];
	res[1] = _x[1] - pp[1];
	res[2] = _x[2] - pp[2];
	res[3] = _x[3] - pp[3];
	res[4] = _x[4] - pp[4];
	res[5] = _x[5] - pp[5];
	res[6] = _x[6] - pp[6];
	return( res );
    }

    ParticleP3D operator*( double x ) const { 
	ParticleP3D res;
	res[0] = _x[0]*x;
	res[1] = _x[1]*x;
	res[2] = _x[2]*x;
	res[3] = _x[3]*x;
	res[4] = _x[4]*x;
	res[5] = _x[5]*x;
	res[6] = _x[6]*x;
	return( res );
    }

    /*! \brief Saves data to stream.
     */
    void save( std::ostream &s ) const {
	write_double( s, _x[0] );
	write_double( s, _x[1] );
	write_double( s, _x[2] );
	write_double( s, _x[3] );
	write_double( s, _x[4] );
	write_double( s, _x[5] );
	write_double( s, _x[6] );
    }
    
    friend ParticleP3D operator*( double x, const ParticleP3D &pp );
};


inline std::ostream &operator<<( std::ostream &os, const ParticleP3D &pp )
{
    os << "("
       << std::setw(12) << pp(0) << ", "
       << std::setw(12) << pp(1) << ", "
       << std::setw(12) << pp(2) << ", "
       << std::setw(12) << pp(3) << ", "
       << std::setw(12) << pp(4) << ", "
       << std::setw(12) << pp(5) << ", "
       << std::setw(12) << pp(6) << ")";
    return( os );
}


inline ParticleP3D operator*( double x, const ParticleP3D &pp )
{
    ParticleP3D res;
    res[0] = pp[0]*x;
    res[1] = pp[1]*x;
    res[2] = pp[2]*x;
    res[3] = pp[3]*x;
    res[4] = pp[4]*x;
    res[5] = pp[5]*x;
    res[6] = pp[6]*x;
    return( res );
}




/* ************************************************************************************** *
 * Particle classes                                                                       *
 * ************************************************************************************** */


/*! \brief %Particle base class
 *
 *  A parent class for common (non-templated) parts of %Particle classes.
 */
class ParticleBase
{
protected:

    particle_status_e   _status;       /*!< \brief Status of particle */
    double              _IQ;           /*!< \brief Current or charge of particle
					*
					*   The sign is forced to be same as the sign of q/m.
					*
					*   In time-dependent simulations particles are 
					*   localized and IQ is charge, in time-independent
					*   simulation IQ is current. In case of 2D 
					*   simulations this is in units A/m or C/m.
					*   In 3D and cylindrically symmetric simulations 
					*   unit is A or C.
					*/
    double              _q;            /*!< \brief Charge q [C] */
    double              _m;            /*!< \brief Mass m [kg] */
    int                 _gen;            /*!< \brief Particle generation [#] */

    ParticleBase( double IQ, double q, double m ) 
	: _status(PARTICLE_OK), _q(q) {
	_m = fabs(m);
	if( _q < 0 )
	    _IQ = -fabs(IQ);
	else
	    _IQ = fabs(IQ);
    _gen=0;
    }

    ParticleBase( double IQ, double q, double m, int gen ) 
	: _status(PARTICLE_OK), _q(q) {
	_m = fabs(m);
	if( _q < 0 )
	    _IQ = -fabs(IQ);
	else
	    _IQ = fabs(IQ);
    _gen = gen; // Constructor for adding also the generation of the particles
    }

    /*! \brief Constructor for loading particle from a file.
     */
    ParticleBase( std::istream &s ) {
	_status = (particle_status_e)read_int32( s );
	_IQ = read_double( s );
	_q = read_double( s );
	_m = read_double( s );	
	_gen = read_int32( s );	
    }

    ~ParticleBase() {}

public:

    /*! \brief Return particle status.
     */
    particle_status_e get_status() { return( _status ); }

    /*! \brief Set particle status.
     */
    void set_status( particle_status_e status ) { _status = status; }

    /*! \brief Return current or charge carried by trajectory or
     *  particle cloud [A/C].
     *
     *  With MODE_2D, this in A/m or C/m.
     */
    double IQ() const { return( _IQ ); }

    /*! \brief Return particle charge (q) [C].
     */
    double q() const { return( _q ); }

    /*! \brief Return particle mass (m) [kg].
     */
    double m() const { return( _m ); }

    /*! \brief Return charge per mass ratio (q/m) [C/kg].
     */
    double qm() const { return( _q/_m ); }

    /*! \brief Return charge per mass ratio (q/m) [C/kg].
     */
    int gen() const { return( _gen ); }

    /*! \brief Saves data to stream.
     */
    void save( std::ostream &s ) const {
	write_int32( s, _status );
	write_double( s, _IQ );
	write_double( s, _q );
	write_double( s, _m );
	write_int32( s, _gen );
    }
};


/*! \brief %Particle class in some geometry.
 *
 *  Contains a templated description of one particle in some
 *  geometry. The geometry is selected throught templation, where PP
 *  is either ParticleP2D, ParticlePCyl or ParticleP3D. The particle
 *  classes for defined geometry modes can be used through defined
 *  types Particle2D, ParticleCyl and Particle3D.
 */
template<class PP> class Particle : public ParticleBase
{
    std::vector<PP>     _trajectory;   /*!< \brief Container of trajectory points */
    PP                  _x;            /*!< \brief Current coordinates of particle */

public:

    /*! \brief Constructor for particle.
     *
     *  \param IQ Current or charge of the macroparticle (A/C).
     *  \param q  Charge of microscopic particle (C).
     *  \param m  Mass of microscopic particle (kg).
     *  \param x  Particle coordinates (m, m/s, m, m/s, ...).
     */
    Particle( double IQ, double q, double m, const PP &x ) 
	: ParticleBase(IQ,q,m), _x(x) {}
    Particle( double IQ, double q, double m, int gen, const PP &x ) 
	: ParticleBase(IQ,q,m,gen), _x(x) {}

    /*! \brief Constructor for loading particle from a file.
     */
    Particle( std::istream &s ) 
	: ParticleBase( s ) {

	uint32_t N  = read_int32( s );
	_trajectory.reserve( N );
	for( uint32_t a = 0; a < N; a++ )
	    _trajectory.push_back( PP( s ) );
	_x = PP( s );
    }

    /*! \brief Destructor.
     */
    ~Particle() {}

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator()( int i ) { return( _x(i) ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator()( int i ) const { return( _x(i) ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator[]( int i ) { return( _x(i) ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator[]( int i ) const { return( _x(i) ); }

    /*! \brief Returns the location of particle in Vec3D.
     */
    Vec3D location() const { return( _x.location() ); }

    /*! \brief Returns the velocity of particle in Vec3D.
     */
    Vec3D velocity() const { return( _x.velocity() ); }

    /*! \brief Return reference to coordinate data.
     */
    PP &x() { return( _x ); }

    /*! \brief Return const reference to coordinate data.
     */
    const PP &x() const { return( _x ); }

    /*! \brief Return reference to trajectory data.
     */
    PP &traj( int i ) { return( _trajectory[i] ); }

    /*! \brief Return const reference to trajectory data.
     */
    const PP &traj( int i ) const { return( _trajectory[i] ); }

    /*! \brief Return number of trajectory points of particle.
     */
    size_t traj_size( void ) const { return( _trajectory.size() ); }

    /*! \brief Add trajectory point to the end of the trajectory.
     */
    void add_trajectory_point( const PP &x ) { _trajectory.push_back( x ); }

    /*! \brief Define trajectory by copying.
     */
    void copy_trajectory( const std::vector<PP> &traj ) { _trajectory = traj; }

    /*! \brief Clears the particle trajectory.
     */
    void clear_trajectory( void ) { _trajectory.clear(); }

    /*! \brief Resets the particle trajectory.
     *
     *  Resets the particle coordinates to the start of the
     *  trajectory, clears the trajectory and sets particle status as
     *  PARTICLE_OK.
     */
    void reset_trajectory( void ) { 
	_status = PARTICLE_OK;
	if( _trajectory.size() != 0 ) {
	    _x = _trajectory[0];
	    _trajectory.clear(); 
	}
    }

    /*! \brief Saves data to stream.
     */
    void save( std::ostream &s ) const {
	ParticleBase::save( s );
	write_int32( s, _trajectory.size() );
	for( uint32_t a = 0; a < _trajectory.size(); a++ )
	    _trajectory[a].save( s );
	_x.save( s );
    }
    
    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const {
	size_t a;
        os << "**Particle\n";
	switch( _status ) {
	case PARTICLE_OK:
	    os << "stat = PARTICLE_OK\n";
	    break;
	case PARTICLE_OUT:
	    os << "stat = PARTICLE_OUT\n";
	    break;
	case PARTICLE_COLL:
	    os << "stat = PARTICLE_COLL\n";
	    break;
	case PARTICLE_BADDEF:
	    os << "stat = PARTICLE_BADDEF\n";
	    break;
	case PARTICLE_TIME:
	    os << "stat = PARTICLE_TIME\n";
	    break;
	case PARTICLE_NSTP:
	    os << "stat = PARTICLE_NSTP\n";
	    break;
	case PARTICLE_STRIP:
	    os << "stat = PARTICLE_STRIP\n";
	    break;
	}
	os << "IQ   = " << _IQ << "\n";
	os << "q    = " << _q << "\n";
	os << "m    = " << _m << "\n";
	os << "x    = " << _x << "\n";
	os << "Trajectory:\n";
	for( a = 0; a < _trajectory.size(); a++ )
	    os << "x[" << a << "] = " << _trajectory[a] << "\n";

	/*
	std::cout << "Trajectory:\n";
	for( a = 0; a < _trajectory.size(); a++ ) {
	    std::cout << "x[" << a << "] = (";
	    uint32_t b;
	    const PP &tp = _trajectory[a];
	    if( tp.size() > 0 ) {
		for( b = 0; b < tp.size()-1; b++ )
		    std::cout << tp[b] << ", ";
		std::cout << tp[b] << ")\n";
	    } else {
		std::cout << ")\n";
	    }
	}
	*/
    }
};


/*! \brief Particle class in 2D.
 *
 *  A typedef for templated class Particle.
 */
typedef Particle<ParticleP2D>  Particle2D;


/*! \brief Particle class in Cylindrical symmetry.
 *
 *  A typedef for templated class Particle.
 */
typedef Particle<ParticlePCyl> ParticleCyl;


/*! \brief Particle class in 3D.
 *
 *  A typedef for templated class Particle.
 */
typedef Particle<ParticleP3D>  Particle3D;



/*! \brief Temporary data bundle for particle iterators.
 */
struct ParticleIteratorData {
    MeshScalarField          *_scharge;  /*!< \brief Space charge field or NULL. */
    const VectorField        *_efield;   /*!< \brief Electric field or NULL. */
    const VectorField        *_bfield;   /*!< \brief Magnetic field or NULL. */
    const Geometry           *_geom;     /*!< \brief Geometry. */
    double                    _qm;       /*!< \brief Precalculated q/m. */
    const CallbackFunctorD_V *_bsup_cb;  /*!< \brief B-field plasma suppression callback. */
    const CallbackFunctorD_V *_esup_cb;  /*!< \brief E-field plasma suppression callback. */ // CP 20240517
    bool                      _relativistic; /*!< \brief Do relativistic particle calc? */

    ParticleIteratorData( MeshScalarField *scharge, const VectorField *efield, 
			  const VectorField *bfield, const Geometry *geom ) 
	: _scharge(scharge), _efield(efield), _bfield(bfield), 
	  _geom(geom), _qm(0.0), _bsup_cb(0), _relativistic(false) {}

    /*! \brief Set B-field potential dependent suppression callback.
     */
    void set_bfield_suppression_callback( const CallbackFunctorD_V *bsup_cb ) {
	_bsup_cb = bsup_cb;
    }

    // CP 20240517
    /*! \brief Set E-field potential dependent suppression callback.
     */
    void set_efield_suppression_callback( const CallbackFunctorD_V *esup_cb ) {
	_esup_cb = esup_cb;
    }

    /*! \brief Set relativistic particle iteration.
     */
    void set_relativistic( bool enable ) {
	_relativistic = enable;
    }

};


#endif

