/*! \file particledatabase.hpp
 *  \brief %Particle databases
 */

/* Copyright (c) 2005-2015 Taneli Kalvas. All rights reserved.
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

#ifndef PARTICLEDATABASE_HPP
#define PARTICLEDATABASE_HPP 1


#include "scalarfield.hpp"
#include "vectorfield.hpp"
#include "particles.hpp"
#include "trajectorydiagnostics.hpp"
#include "particlestatistics.hpp"
#include "constants.hpp"
#include "types.hpp"


/*! \brief Trajectory handler callback.
 */
class TrajectoryHandlerCallback {
public:

    /*! \brief Virtual destructor.
     */
    virtual ~TrajectoryHandlerCallback() {}

    virtual void operator()( ParticleBase *particle, ParticlePBase *xcur, ParticlePBase *xend ) = 0;
};


/*! \brief Trajectory end callback.
 */
class TrajectoryEndCallback {
public:

    /*! \brief Virtual destructor.
     */
    virtual ~TrajectoryEndCallback() {}

    /*! \brief Operator called when particle calculation ends.
     */
    virtual void operator()( ParticleBase *particle, class ParticleDataBase *pdb ) = 0;
};


/*! \brief Trajectory surface collision callback.
 */
class TrajectorySurfaceCollisionCallback {
public:

    /*! \brief Virtual destructor.
     */
    virtual ~TrajectorySurfaceCollisionCallback() {}

    /*! \brief Operator called when particle collides with surface.
     *
     *  Only with triangulated surfaces. The colliding \a particle,
     *  collision location \a x, triangle index \a tri and the
     *  parametric location \a (s,t) are given as parameters in the
     *  call.
     *
     *  The coordinate \a s is the parametric distance from vertex 0
     *  to vertex 1 and coordinate \a t is the parametric distance
     *  from vertex 0 to vertex 2 inside the triangle.
     */
    virtual void operator()( ParticleBase *particle, ParticlePBase *x, uint32_t tri, 
			     double s, double t ) = 0;
};


/*! \brief Magnetic field plasma suppression for positive ion extraction.
 *
 *  Defines a magnetic field suppression, which is dependent on
 *  electric potential. If the electric potential at the location \a x
 *  is \b larger than the defined potential limit, the magnetic field
 *  will be suppressed to zero. Otherwise, no magnetic field
 *  suppression will be made. This produces a hard boundary for the
 *  magnetic field suppression.
 */
class PPlasmaBfieldSuppression : public CallbackFunctorD_V {

    double                 _phi;    /*!< \brief Limit for potential. */
    const MeshScalarField &_epot;   /*!< \brief Electric potential field. */
    
public:

    /*! \brief Constructor setting electric potential field and potential limit.
     */
    PPlasmaBfieldSuppression( const MeshScalarField &epot, double phi ) 
	: _phi(phi), _epot(epot) {}

    /*! \brief Destructor.
     */
    ~PPlasmaBfieldSuppression() {}

    /*! \brief Suppression function.
     */
    virtual double operator()( const Vec3D &x ) const {
	if( _epot( x ) > _phi )
	    return( 0.0 );
	return( 1.0 );
    }
};


/*! \brief Magnetic field plasma suppression for negative ion extraction.
 *
 *  Defines a magnetic field suppression, which is dependent on
 *  electric potential. If the electric potential at the location \a x
 *  is \b smaller than the defined potential limit, the magnetic field
 *  will be suppressed to zero. Otherwise, no magnetic field
 *  suppression will be made. This produces a hard boundary for the
 *  magnetic field suppression.
 */
class NPlasmaBfieldSuppression : public CallbackFunctorD_V {

    double             _phi;    /*!< \brief Limit for potential. */
    const MeshScalarField &_epot;   /*!< \brief Electric potential field. */
    
public:

    /*! \brief Constructor setting electric potential field and potential limit.
     */
    NPlasmaBfieldSuppression( const MeshScalarField &epot, double phi ) 
	: _phi(phi), _epot(epot) {}

    /*! \brief Destructor.
     */
    ~NPlasmaBfieldSuppression() {}
    
    /*! \brief Suppression function.
     */
    virtual double operator()( const Vec3D &x ) const {
	if( _epot( x ) < _phi )
	    return( 0.0 );
	return( 1.0 );
    }
};


/* ******************************************************************************************* *
 * ParticleDataBase classes                                                                    *
 * ******************************************************************************************* */


/*! \brief %Particle database base class.
 *
 * %Particle database base class holds the definitions of particle
 * iteration parameters. Base class also provides a possibility for
 * general pointer to particle database and virtual functions for
 * accessing particles.
 */
class ParticleDataBase {

    class ParticleDataBaseImp *_imp;
    
protected:

/* ************************************** *
 * Constructors                           *
 * ************************************** */

    /*! \brief Constructor.
     */
    ParticleDataBase();

    /*! \brief Copy constructor.
     */
    ParticleDataBase( const ParticleDataBase &pdb );

    /*! \brief Assignment.
     */
    const ParticleDataBase &operator=( const ParticleDataBase &pdb );

    /*! \brief Set particle database implementation pointer.
     *
     *  Used by child constructors.
     */
    void set_implementation_pointer( class ParticleDataBaseImp *imp );

public:

/* ************************************** *
 * Destructor                             *
 * ************************************** */

    /*! \brief Virtual destructor.
     */
    virtual ~ParticleDataBase();

/* ****************************************** *
 * Particle iteration settings and statictics *
 * ****************************************** */

    /*! \brief Set the number of threads used for calculation.
     *
     *  \deprecated This function is deprecated (it does nothing) and
     *  is replaced by global IBSimu::set_thread_count().
     */
    void set_thread_count( uint32_t threadcount ) {}

    /*! \brief Set the accuracy requirement for calculation.
     *
     *  Accuracy requirements default to \a epsabs = 1.0e-6 and \a
     *  epsrel = 1.0e-6.
     */
    void set_accuracy( double epsabs, double epsrel );

    /*! \brief Set magnetic field suppression location depedent
     *  callback functor.
     *
     *  Can be used to suppress the magnetic field effect on particles
     *  inside the plasma by making an object with a pointer to
     *  electric potential. Comparing to local potential value the
     *  magnetic field suppression can be localized to a volume with
     *  selected potential range.
     */
    void set_bfield_suppression( const CallbackFunctorD_V *functor );

    /*! \brief Set trajectory handler callback. 
     */
    void set_trajectory_handler_callback( TrajectoryHandlerCallback *thand_cb );

    /*! \brief Set trajectory end callback. 
     */
    void set_trajectory_end_callback( TrajectoryEndCallback *tend_cb );

    /*! \brief Set trajectory surface collision callback. 
     */
    void set_trajectory_surface_collision_callback( TrajectorySurfaceCollisionCallback *tsur_cb );

    /*! \brief Set relativistic particle iteration.
     */
    void set_relativistic( bool enable );

    /*! \brief Set surface collision model to be used.
     *
     *  Defaults to false, which enables the inside-test for collisions.
     */
    void set_surface_collision( bool surface_collision );

    /*! \brief Set the interpolation type to polynomial(true) or linear(false).
     *
     *  \deprecated This function is deprecated. It can still be used
     *  to set polynomial or linear interpolation. Please use set_trajectory_interpolation().
     */
    void set_polyint( bool polyint );
    
    /*! \brief Get current interpolation type.
     *
     *  \deprecated This function is deprecated. It can still be used
     *  to get the interpolation type for polynomial(true) or
     *  linear(false) interpolation. If interpolation type is
     *  something else, the function will return false. Please use
     *  get_trajectory_interpolation().
     */
    bool get_polyint( void ) const;
    
    /*! \brief Set trajectory interpolation type.
     *
     *  Defaults to polynomial interpolation.
     */
    void set_trajectory_interpolation( trajectory_interpolation_e intrp );

    /*! \brief Get trajectory interpolation type.
     */
    trajectory_interpolation_e get_trajectory_interpolation( void ) const;

    /*! \brief Set space charge deposition type.
     *
     *  Defaults to SCHARGE_DEPOSITION_PIC.
     */
    void set_scharge_deposition( scharge_deposition_e type );

    /*! \brief Get space charge deposition type.
     */
    scharge_deposition_e get_scharge_deposition( void ) const;

    /*! \brief Set maximum number of steps to iterate.
     *
     *  One thousand (1000) steps is the default
     */
    void set_max_steps( uint32_t maxsteps );

    /*! \brief Set maximum lifetime of particle in simulation.
     *
     *  One millisecond (1e-3 sec) is the default
     */
    void set_max_time( double maxt );

    /*! \brief Set trajectory saving.
     *
     *  If \a save_points is true, all trajectory intersection with
     *  mesh boundaries are saved. If \a save_points is false
     *  (default), only the calculation points are saved.
     */
    void set_save_all_points( bool save_points );

    /*! \brief Set trajectory saving.
     *
     *  If \a div is zero, no trajectories are saved.
     *  If \a div is one, every trajectory is saved.
     *  If \a div N>1, every Nth trajectory is saved.
     *
     *  By default, all trajectories are saved.
     */
    void set_save_trajectories( uint32_t div );

    /*! \brief Get trajectory saving.
     *
     *  If \a div is zero, no trajectories are saved.
     *  If \a div is one, every trajectory is saved.
     *  If \a div N>1, every Nth trajectory is saved.
     */
    uint32_t get_save_trajectories( void ) const;

    /*! \brief Set particle mirroring on boundaries.
     *
     *  Mirroring is set for (xmin,xmax,ymin,ymax,zmin,zmax)
     *  borders. Mirroring is always false for directions that do not
     *  exist.
     *
     *  By default, mirroring is disabled for all boundaries.
     */
    void set_mirror( const bool mirror[6] );

    /*! \brief Get particle mirroring on boundaries.
     *
     *  Mirroring is read for (xmin,xmax,ymin,ymax,zmin,zmax)
     *  borders. Mirroring is always false for directions that do not
     *  exist.
     */
    void get_mirror( bool mirror[6] ) const;

    /*! \brief Get the number of iteration rounds done.
     */
    int get_iteration_number( void ) const;

    /*! \brief Return sum of defined beam space charge density.
     *
     *  Returns the summed beam space charge density for all beams
     *  defined with "add_beam" functions. Does not work with
     *  individually added particle trajectories. \a rhosum is cleared
     *  with particle database clearing function clear().
     *
     *  Can be used to program plasma electron density with
     *  EpotProblem::set_pexp_plasma() for example. Please note that
     *  it gives to accumulated charge density which might be
     *  incorrect for multi-beam extraction simulation defined with
     *  several calls to "add_beam" functions.
     */
    double get_rhosum( void ) const;

    /*! \brief Set sum of defined beam space charge density.
     */
    void set_rhosum( double rhosum );

    /*! \brief Get particle iterator statistics.
     */
    const ParticleStatistics &get_statistics( void ) const;

/* ************************************** *
 * Information and queries                *
 * ************************************** */

    /*! \brief Return geometry mode.
     */
    geom_mode_e geom_mode() const;

    /*! \brief Returns particle count.
     */
    size_t size( void ) const;

    /*! \brief Returns a reference to particle \a i.
     */
    virtual ParticleBase &particle( uint32_t i ) = 0;

    /*! \brief Returns a const reference to particle \a i.
     */
    virtual const ParticleBase &particle( uint32_t i ) const = 0;

    /*! \brief Returns the length of trajectory \a i.
     */
    double traj_length( uint32_t i ) const;
    
    /*! \brief Returns number of trajectory points for particle \a i.
     */
    size_t traj_size( uint32_t i ) const;
    
    /*! \brief Gets the particle \a i trajectory point \a j as particle point.
     */
    virtual const ParticlePBase &trajectory_point( uint32_t i, uint32_t j ) const = 0;

    /*! \brief Gets the particle \a i trajectory point \a j into \a vel, \a loc and \a t.
     */
    void trajectory_point( double &t, Vec3D &loc, Vec3D &vel, uint32_t i, uint32_t j ) const;

    /*! \brief Gets trajectory diagnostic \a diagnostics at plane \a
     *  axis = \a val in trajectory diagnostic data object \a tdata.
     */
    void trajectories_at_plane( TrajectoryDiagnosticData &tdata, 
				coordinate_axis_e axis,
				double val,
				const std::vector<trajectory_diagnostic_e> &diagnostics ) const;

    /*! \brief Build trajectory density field.
     *
     *  This method builds a current density map (A/m2) from the
     *  particles. It is similar to space charge map, but the particle
     *  velocities do not affect the value, only the trajectory
     *  density and the current carried by the trajectories.
     *
     *  For historical reasons this function is called
     *  build_trajectory_density_field() even though it produces a
     *  current density map. The scalar field \a tdens mesh can differ
     *  from geometry mesh.
     */
    void build_trajectory_density_field( MeshScalarField &tdens ) const;

/* ************************************** *
 * Particle and trajectory clearing       *
 * ************************************** */

    /*! \brief Clears the particle database of all particles.
     *
     *  Clears the database of particles. Also clears beam space
     *  charge sum.
     */
    void clear( void );

    /*! \brief Clears the particle trajectory database.
     *
     *  The particle properties are otherwise conserved, but existing
     *  trajectories are cleared.
     */
    void clear_trajectories( void );

    /*! \brief Clears particle \a a in the particle trajectory database.
     *
     *  The particle properties are otherwise conserved, but existing
     *  trajectory are cleared.
     */
    void clear_trajectory( size_t a );

    /*! \brief Clears the particle trajectory database and set particles to initial coordinates
     *
     *  The particles are reset with initial coordinates and trajectories are cleared.
     *  Sets particle status for all particles as PARTICLE_OK.
     */
    void reset_trajectories( void );

    /*! \brief Clears the particle trajectory and set the particle to initial coordinates
     *
     *  The particle \a a is resetted with initial coordinates and trajectory is cleared.
     *  Sets particle status as PARTICLE_OK.
     */
    void reset_trajectory( size_t a );

/* ************************************** *
 * Particle definition                    *
 * ************************************** */

    /*! \brief Reserve memory for \a size particles.
     */
    void reserve( size_t size );

/* ************************************** *
 * Particle iterators                     *
 * ************************************** */

    /*! \brief Iterate particles through the geometry.
     *
     *  The particles defined in particle database \a pdb are iterated
     *  through electric field \a efield and magnetic field \a bfield 
     *  in geometry \a geom. Space charge density field \a scharge is set
     *  from the particle trajectories.
     *
     *  The meshes of the geometry, the space charge field and the
     *  electric field have to be equal. The magnetic field mesh can
     *  be selected independently. This allows minimization of the
     *  memory use in the case where electric field needs high
     *  resolution, but magnetic field is relatively smooth.
     */
    void iterate_trajectories( MeshScalarField &scharge, const VectorField &efield, 
			       const VectorField &bfield );

    /*! \brief Step particles forward by time step dt.
     *
     *  The particles defined in particle database \a pdb are stepped
     *  forward one time step \a dt through electric field \a efield and
     *  magnetic field \a bfield in geometry \a geom.
     */
    void step_particles( MeshScalarField &scharge, const VectorField &efield, 
			 const VectorField &bfield, double dt );

/* ************************************** *
 * Debugging, plotting and saving         *
 * ************************************** */

    /*! \brief Saves data to a new file \a filename.
     */
    virtual void save( const std::string &filename ) const = 0;

    /*! \brief Saves data to stream.
     */
    virtual void save( std::ostream &os ) const = 0;

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const = 0;
};


/*! \brief %Particle database class for two dimensions.
 *
 *  %Particle database holds the definitions of particles and possibly
 *  the trajectories of particles it the particle iterator has saved
 *  them. %ParticleDataBase2D provides a variety of convenience
 *  functions for defining particle beams.
 *
 *  Particles are always stored in the database in the order they are
 *  defined. When reading back the simulation results, the order can
 *  be used to identify the particles.
 */
class ParticleDataBase2D : public ParticleDataBase {

    class ParticleDataBase2DImp *_imp;

public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Constructor.
     */
    ParticleDataBase2D( const Geometry &geom );

    /*! \brief Copy constructor.
     */
    ParticleDataBase2D( const ParticleDataBase2D &pdb );

    /*! \brief Constructor for loading particle statistics from a file.
     */
    ParticleDataBase2D( std::istream &s, const Geometry &geom );

    /*! \brief Destructor.
     */
    ~ParticleDataBase2D();

    /*! \brief Assignment.
     */
    const ParticleDataBase2D &operator=( const ParticleDataBase2D &pdb );
    
/* ************************************** *
 * Information and queries                *
 * ************************************** */

    /*! \brief Returns a reference to particle \a i.
     */
    virtual Particle2D &particle( uint32_t i );

    /*! \brief Returns a const reference to particle \a i.
     */
    virtual const Particle2D &particle( uint32_t i ) const;
    
    /*! \brief Gets the particle \a i trajectory point \a j as particle point.
     */
    virtual const ParticleP2D &trajectory_point( uint32_t i, uint32_t j ) const;

    using ParticleDataBase::trajectory_point;

    /*! \brief Gets trajectory data at plane \a axis = \a val as
     *  particles in vector \a tdata.
     */
    void trajectories_at_plane( std::vector<Particle2D> &tdata,
				coordinate_axis_e axis,
				double val ) const;

    using ParticleDataBase::trajectories_at_plane;

/* ************************************** *
 * Particle definition                    *
 * ************************************** */

    /*! \brief Add one particle.
     *
     *  Adds one particle to particle database. Particle properties
     *  are: \a IQ is the current (A) in time-independent or charge
     *  (C) in time-dependent simulations carried by the particle
     *  cloud that the simulated particle represents, \a q is the
     *  charge state of the microscopic particle (in multiples of e),
     *  \a m is the mass of the microscopic particle (u) and \a x
     *  contains the time, position (m) and velocity (m/s) of the
     *  particle.
     */
    void add_particle( double IQ, double q, double m, const ParticleP2D &x );

    /*! \brief Add one particle.
     *
     *  Adds one particle to database.
     */
    void add_particle( const Particle2D &p );

/* ************************************** *
 * Particle beam definition               *
 * ************************************** */

    /*! \brief Add a 2d beam with energies.
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particle is \a q (in
     *  multiples of e), mass is \a m (u). The beam is defined on a
     *  line from (\a x1, \a y1) to (\a x2, \a y2). The beam
     *  propagates into a direction 90 degrees clockwise from the
     *  direction of vector pointing from (\a x1, \a y1) to (\a x2, \a
     *  y2) with a mean energy \a E (eV). The beam has parallel
     *  temperature \a Tp (eV) and transverse temperature \a Tt (eV).
     *
     *  The particle speeds of the beam in direction \a i are sampled
     *  from a gaussian distribution with standard deviation dv_i =
     *  sqrt(T_i*e/m), where \a T_i is the beam temperature in
     *  direction \a (eV), \a e is electron charge (C) and m is the
     *  mass of the ion (kg).
     *
     *  Space charge J/v is added to the \a rhosum variable.
     */
    void add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
				  double E, double Tp, double Tt, 
				  double x1, double y1, double x2, double y2 );

    /*! \brief Add a 2d beam with velocities.
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particle is \a q (in
     *  multiples of e), mass is \a m (u). The beam is defined on a
     *  line from (\a x1, \a y1) to (\a x2, \a y2). The beam
     *  propagates into a direction 90 degrees clockwise from the
     *  direction of vector pointing from (\a x1, \a y1) to (\a x2, \a
     *  y2) with a mean velocity \a v (m/s). The beam has parallel
     *  gaussian velocity distribution with standard deviation \a dvp
     *  (m/s) and transverse gaussian velocity distribution with
     *  standard deviation \a dvt (m/s).
     *
     *  Space charge J/v is added to the \a rhosum variable.
     */
    void add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
				    double v, double dvp, double dvt, 
				    double x1, double y1, double x2, double y2 );

    /*! \brief Add a 2d beam with defined KV emittance.
     *
     *  Adds a beam consisting of \a N particles and a total current
     *  of \a I (A). The particles are defined to have equal currents,
     *  charge \a q (in multiples of e) and mass \a m (u). The
     *  starting energy of the beam is \a Ex (eV) and the starting
     *  location \a (x0,y0) (center point). The beam propagates to
     *  positive x-direction. The beam is made to match Twiss
     *  parameters \f$ \alpha \f$ (a) and \f$ \beta \f$ (b) in
     *  projectional direction (y,y'). The rms-emittance of the beam
     *  is made to match \f$ \epsilon_{\mathrm{rms}} \f$ (e).
     *
     *  The beam spread in the projectional space is made according to
     *  KV/hard-edged (Kapchinsky-Vladimirsky) distribution.
     *
     */
    void add_2d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
					double a, double b, double e,
					double Ex, double x0, double y0 );

    /*! \brief Add a 2d beam with defined gaussian emittance.
     *
     *  Adds a beam consisting of \a N particles and a total current
     *  of \a I (A). The particles are defined to have equal currents,
     *  charge \a q (in multiples of e) and mass \a m (u). The
     *  starting energy of the beam is \a Ex (eV) and the starting
     *  location \a (x0,y0) (center point). The beam propagates to
     *  positive x-direction. The beam is made to match Twiss
     *  parameters \f$ \alpha \f$ (a) and \f$ \beta \f$ (b) in
     *  projectional direction (y,y'). The rms-emittance of the beam
     *  is made to match \f$ \epsilon_{\mathrm{rms}} \f$ (e).
     *
     *  The beam spread in the projectional space is made according to
     *  Gaussian distribution.
     *
     */
    void add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
					      double a, double b, double e,
					      double Ex, double x0, double y0 );

/* ************************************** *
 * Debugging, plotting and saving         *
 * ************************************** */

    /*! \brief Saves data to a new file \a filename.
     */
    virtual void save( const std::string &filename ) const;

    /*! \brief Saves data to stream.
     */
    virtual void save( std::ostream &s ) const;

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;
};




/*! \brief %Particle database class for cylindrical systems.
 *
 *  %Particle database holds the definitions of particles and possibly
 *  the trajectories of particles it the particle iterator has saved
 *  them. %ParticleDataBaseCyl provides a variety of convenience
 *  functions for defining particle beams.
 *
 *  Particles are always stored in the database in the order they are
 *  defined. When reading back the simulation results, the order can
 *  be used to identify the particles.
 */
class ParticleDataBaseCyl : public ParticleDataBase {

    class ParticleDataBaseCylImp *_imp;

public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Constructor.
     */
    ParticleDataBaseCyl( const Geometry &geom );

    /*! \brief Copy constructor.
     */
    ParticleDataBaseCyl( const ParticleDataBaseCyl &pdb );

    /*! \brief Constructor for loading particle statistics from a file.
     */
    ParticleDataBaseCyl( std::istream &s, const Geometry &geom );

    /*! \brief Destructor.
     */
    ~ParticleDataBaseCyl();

    /*! \brief Assignment.
     */
    const ParticleDataBaseCyl &operator=( const ParticleDataBaseCyl &pdb );

/* ************************************** *
 * Information and queries                *
 * ************************************** */

    /*! \brief Returns a reference to particle \a i.
     */
    virtual ParticleCyl &particle( uint32_t i );

    /*! \brief Returns a const reference to particle \a i.
     */
    virtual const ParticleCyl &particle( uint32_t i ) const;
    
    /*! \brief Gets the particle \a i trajectory point \a j as particle point.
     */
    virtual const ParticlePCyl &trajectory_point( uint32_t i, uint32_t j ) const;

    using ParticleDataBase::trajectory_point;

    /*! \brief Gets trajectory data at plane \a axis = \a val as
     *  particles in vector \a tdata.
     */
    void trajectories_at_plane( std::vector<ParticleCyl> &tdata,
				coordinate_axis_e axis,
				double val ) const;

    using ParticleDataBase::trajectories_at_plane;

    /*! \brief Export particle data as Path Manager data
     *
     *  Makes trajectory diagnostics on a plane x = \a val. Each
     *  particle in cylindrical symmetry simulations corresponds to a
     *  ring in 3D and has a different current (according to starting
     *  radius). In Path Manager particles have equal currents.
     *  Therefore the particles... are converted in a complicated way.
     *
     *  The particle properties on the diagnostics plane are written
     *  to file \a filename in the format required by the Path Manager
     *  program. Reference particle for the output has energy \a ref_E
     *  (in eV), charge state \a ref_q (in electron charges) and mass
     *  \a ref_m (in atomic mass units).
     *
     *  The file contains a header with information about the
     *  reference particle. After the header, the file has a line for
     *  each particle, which contain the following columns: \a x, \a
     *  x', \a y, \a y', \a phase, \a (p_z-p_ref)/p_ref, \a 0-flag, \a
     *  q, \a m. The particles are mirrored according the particle
     *  database mirroring settings. The phase is randomized for each
     *  particle between 0 and 2 pi corresponding to continuous beam.
     */
    void export_path_manager_data( const std::string &filename, 
				   double ref_E, double ref_q, double ref_m, 
				   double val, uint32_t Np ) const;

/* ************************************** *
 * Particle definition                    *
 * ************************************** */

    /*! \brief Add one particle.
     *
     *  Adds one particle to particle database. Particle properties
     *  are: \a IQ is the current (A) in time-independent or charge
     *  (C) in time-dependent simulations carried by the particle
     *  cloud that the simulated particle represents, \a q is the
     *  charge state of the microscopic particle (in multiples of e),
     *  \a m is the mass of the microscopic particle (u) and \a x
     *  contains the time, position (m) and velocity (m/s) of the
     *  particle.
     */
    void add_particle( double IQ, double q, double m, const ParticlePCyl &x );

    /*! \brief Add one particle.
     *
     *  Adds one particle to database.
     */
    void add_particle( const ParticleCyl &p );

 /* ************************************** *
 * Particle beam definition               *
 * ************************************** */

    /*! \brief Add a 2d beam with energies.
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particle is \a q (in
     *  multiples of e), mass is \a m (u). The beam is defined on a
     *  line from (\a x1, \a y1) to (\a x2, \a y2). The beam
     *  propagates into a direction 90 degrees clockwise from the
     *  direction of vector pointing from (\a x1, \a y1) to (\a x2, \a
     *  y2) with a mean energy \a E (eV). The beam has parallel
     *  temperature \a Tp (eV) and transverse temperature \a Tt (eV).
     *
     *  The particle speeds of the beam in direction \a i are sampled
     *  from a gaussian distribution with standard deviation dv_i =
     *  sqrt(T_i*e/m), where \a T_i is the beam temperature in
     *  direction \a (eV), \a e is electron charge (C) and m is the
     *  mass of the ion (kg).
     *
     *  Space charge J/v is added to the \a rhosum variable.
     */
    void add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
				  double E, double Tp, double Tt, 
				  double x1, double y1, double x2, double y2 );

    /*! \brief Add a 2d beam with total energy
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particle is \a q (in
     *  multiples of e), mass is \a m (u). The beam is defined on a
     *  line from (\a x1, \a y1) to (\a x2, \a y2). The beam
     *  propagates into a direction 90 degrees clockwise from the
     *  direction of vector pointing from (\a x1, \a y1) to (\a x2, \a
     *  y2) with a mean total energy \a Etot (eV). The total energy is
     *  calculated as Etot = 0.5mv^2 + qU. The electric potential U is
     *  evaluated from \a epot. The beam also has parallel temperature
     *  \a Tp (eV) and transverse temperature \a Tt (eV).
     *
     *  The particle speeds of the beam in direction \a i are sampled
     *  from a gaussian distribution with standard deviation dv_i =
     *  sqrt(T_i*e/m), where \a T_i is the beam temperature in
     *  direction \a (eV), \a e is electron charge (C) and m is the
     *  mass of the ion (kg).
     */
    void add_2d_beam_with_total_energy( uint32_t N, double J, double q, double m, 
					double Etot, const ScalarField &epot, 
					double Tp, double Tt, 
					double x1, double y1, double x2, double y2 );

    /*! \brief Add a 2d beam with velocities.
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particle is \a q (in
     *  multiples of e), mass is \a m (u). The beam is defined on a
     *  line from (\a x1, \a y1) to (\a x2, \a y2). The beam
     *  propagates into a direction 90 degrees clockwise from the
     *  direction of vector pointing from (\a x1, \a y1) to (\a x2, \a
     *  y2) with a mean velocity \a v (m/s). The beam has parallel
     *  gaussian velocity distribution with standard deviation \a dvp
     *  (m/s) and transverse gaussian velocity distribution with
     *  standard deviation \a dvt (m/s).
     *
     *  Space charge J/v is added to the \a rhosum variable.
     */
    void add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
				    double v, double dvp, double dvt, 
				    double x1, double y1, double x2, double y2 );

    /*! \brief Add a 2d beam with gaussian profile and velocity distribution
     */
    void add_2d_full_gaussian_beam( uint32_t N, double I, double q, double m,
				    double Ex, double Tp, double Tt, 
				    double x0, double dr );

    /*! \brief Add a 2d beam with defined gaussian emittance.
     *
     *  Adds a beam consisting of \a N particles and a total current
     *  of \a I (A). The particles are defined to have equal currents,
     *  charge \a q (in multiples of e) and mass \a m (u). The
     *  starting energy of the beam is \a Ex (eV) and the starting
     *  location \a x0 (center point). The beam propagates to positive
     *  x-direction. The beam is made to match Twiss parameters \f$
     *  \alpha \f$ (a), \f$ \beta \f$ (b), \f$ \epsilon \f$ (e) in
     *  projectional directions (y,y') and (z,z').
     *
     *  The even beam current distribution is different from other
     *  beam definitions for cylindrical coordinates. This might
     *  change soon.
     *
     *  \todo Redo particle distribution definition from emittance for
     *  cylindrical symmetry systems.
     */
    void add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
					      double a, double b, double e,
					      double Ex, double x0 );

/* ************************************** *
 * Debugging, plotting and saving         *
 * ************************************** */

    /*! \brief Saves data to a new file \a filename.
     */
    virtual void save( const std::string &filename ) const;

    /*! \brief Saves data to stream.
     */
    virtual void save( std::ostream &s ) const;

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;
};



/*! \brief %Particle database class for three dimensions.
 *
 *  %Particle database holds the definitions of particles and possibly
 *  the trajectories of particles it the particle iterator has saved
 *  them. %ParticleDataBase3D provides a variety of convenience
 *  functions for defining particle beams.
 *
 *  Particles are always stored in the database in the order they are
 *  defined. When reading back the simulation results, the order can
 *  be used to identify the particles.
 */
class ParticleDataBase3D : public ParticleDataBase {

    class ParticleDataBase3DImp *_imp;

public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Constructor.
     */
    ParticleDataBase3D( const Geometry &geom );

    /*! \brief Copy constructor.
     */
    ParticleDataBase3D( const ParticleDataBase3D &pdb );

    /*! \brief Constructor for loading particle statistics from a file.
     */
    ParticleDataBase3D( std::istream &s, const Geometry &geom );

    /*! \brief Destructor.
     */
    ~ParticleDataBase3D();

    /*! \brief Assignment.
     */
    const ParticleDataBase3D &operator=( const ParticleDataBase3D &pdb );

/* ************************************** *
 * Information and queries                *
 * ************************************** */

    /*! \brief Returns a reference to particle \a i.
     */
    virtual Particle3D &particle( uint32_t i );

    /*! \brief Returns a const reference to particle \a i.
     */
    virtual const Particle3D &particle( uint32_t i ) const;
    
    /*! \brief Gets the particle \a i trajectory point \a j as particle point.
     */
    virtual const ParticleP3D &trajectory_point( uint32_t i, uint32_t j ) const;

    using ParticleDataBase::trajectory_point;

    /*! \brief Gets trajectory data at plane \a axis = \a val as
     *  particles in vector \a tdata.
     */
    void trajectories_at_plane( std::vector<Particle3D> &tdata,
				coordinate_axis_e axis,
				double val ) const;

    using ParticleDataBase::trajectories_at_plane;

/* ************************************** *
 * Particle definition                    *
 * ************************************** */

    /*! \brief Add one particle.
     *
     *  Adds one particle to particle database. Particle properties
     *  are: \a IQ is the current (A) in time-independent or charge
     *  (C) in time-dependent simulations carried by the particle
     *  cloud that the simulated particle represents, \a q is the
     *  charge state of the microscopic particle (in multiples of e),
     *  \a m is the mass of the microscopic particle (u) and \a x
     *  contains the time, position (m) and velocity (m/s) of the
     *  particle.
     */
    void add_particle( double IQ, double q, double m, const ParticleP3D &x );
    void add_particle( double IQ, double q, double m, int gen, const ParticleP3D &x );

    /*! \brief Add one particle.
     *
     *  Adds one particle to database.
     */
    void add_particle( const Particle3D &p );

/* ************************************** *
 * Particle beam definition               *
 * ************************************** */


    /*! \brief Add a cylindrical beam with energies.
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particles is \a q (in
     *  multiples of e), mass is \a m (u). The beam starting surface
     *  is a disc of radius \a r centered at \a c. The normal
     *  direction of the disc is \a dir3 = \a dir1 x \a dir2. The
     *  first tangent direction is \a dir1 and the second is \a dir1 x
     *  \a dir3. If you want beam to go to positive x-direction, \a
     *  dir1 could be (0,1,0) and \a dir2 (0,0,1) for example.  The
     *  beam total energy \a Etot (eV) is defined in the normal
     *  direction. The total energy is calculated as Etot = 0.5mv^2 +
     *  qU. Temperatures are defined in normal (parallel) direction
     *  and transverse direction as \a Tp (eV) and \a Tt (eV),
     *  respectively.
     *
     *  The particle speeds of the beam in direction \a i are sampled
     *  from a gaussian distribution with standard deviation dv_i =
     *  sqrt(T_i*e/m), where \a T_i is the beam temperature in
     *  direction \a (eV), \a e is electron charge (C) and m is the
     *  mass of the ion (kg).
     */
    void add_cylindrical_beam_with_total_energy( uint32_t N, double J, double q, double m, 
						 double Etot, const ScalarField &epot, 
						 double Tp, double Tt, Vec3D c, 
						 Vec3D dir1, Vec3D dir2, double r );

    /*! \brief Add a cylindrical beam with energies.
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particles is \a q (in
     *  multiples of e), mass is \a m (u). The beam starting surface
     *  is a disc of radius \a r centered at \a c. The normal
     *  direction of the disc is \a dir3 = \a dir1 x \a dir2. The
     *  first tangent direction is \a dir1 and the second is \a dir1 x
     *  \a dir3. If you want beam to go to positive x-direction, \a
     *  dir1 could be (0,1,0) and \a dir2 (0,0,1) for example.  The
     *  beam energy \a E (eV) is defined in the normal
     *  direction. Temperatures are defined in normal (parallel)
     *  direction and transverse direction as \a Tp (eV) and \a Tt
     *  (eV), respectively.
     *
     *  The particle speeds of the beam in direction \a i are sampled
     *  from a gaussian distribution with standard deviation dv_i =
     *  sqrt(T_i*e/m), where \a T_i is the beam temperature in
     *  direction \a (eV), \a e is electron charge (C) and m is the
     *  mass of the ion (kg).
     *
     *  Space charge J/v is added to the \a rhosum variable.
     */
    void add_cylindrical_beam_with_energy( uint32_t N, double J, double q, double m, 
					   double E, double Tp, double Tt, Vec3D c, 
					   Vec3D dir1, Vec3D dir2, double r );

    /*! \brief Add a cylindrical beam with velocities.
     *
     *  Adds a beam consisting of \a N particles. The beam current
     *  density is \a J (A/m^2), charge of beam particles is \a q (in
     *  multiples of e), mass is \a m (u). The beam starting surface
     *  is a disc of radius \a r centered at \a c. The normal
     *  direction of the disc is \a dir3 = \a dir1 x \a dir2. The
     *  first tangent direction is \a dir1 and the second is \a dir1 x
     *  \a dir3. The beam velocity \a v (m/s) is defined in the normal
     *  direction. The velocity distribution standard deviation are
     *  defined in normal (parallel) direction and transverse
     *  directions as \a dvp (m/s) and \a dvt (m/s), respectively
     *
     *  The particle velocities of the beam in direction \a i are
     *  sampled from a gaussian distribution with standard deviation
     *  dv_i, which is related to beam temperature by dv_i =
     *  sqrt(T_i*e/m), where \a T_i is the beam temperature in
     *  direction \a (eV), \a e is electron charge (C) and m is the
     *  mass of the ion (kg).
     *
     *  Space charge J/v is added to the \a rhosum variable.
     */
    void add_cylindrical_beam_with_velocity( uint32_t N, double J, double q, double m, 
					     double v, double dvp, double dvt, Vec3D c, 
					     Vec3D dir1, Vec3D dir2, double r );

    /*! \brief Add a rectangular beam with energies.
     *
     *  Very much like add_cylindrical_beam_with_energy(), but the
     *  beam shape is rectangular. The beam size is defined with two
     *  variables, \a sizex and \a sizey. These sizes are half the
     *  range, i.e. the beam ranges in the first direction from center
     *  point - size1 to center point + size1.
     */
    void add_rectangular_beam_with_energy( uint32_t N, double J, double q, double m, 
					   double E, double Tp, double Tt, Vec3D c, 
					   Vec3D dir1, Vec3D dir2, double size1, double size2 );

    /*! \brief Add a rectangular beam with velocities.
     *
     *  Very much like add_cylindrical_beam_with_velocity(), but the
     *  beam shape is rectangular.  The beam size is defined with two
     *  variables, \a sizex and \a sizey. These sizes are half the
     *  range, i.e. the beam ranges in the first direction from center
     *  point - size1 to center point + size1.
     */
    void add_rectangular_beam_with_velocity( uint32_t N, double J, double q, double m, 
					     double v, double dvp, double dvt, Vec3D c, 
					     Vec3D dir1, Vec3D dir2, double size1, double size2 );


    /*! \brief Add a 3d beam with defined KV emittance.
     *
     *  Adds a beam consisting of \a N particles and a total current
     *  of \a I (A). The particles are defined to have equal currents,
     *  charge \a q (in multiples of e) and mass \a m (u). The
     *  starting energy of the beam is \a E0 (eV) and the starting
     *  location \a c (center point). The beam propagates to the
     *  direction \a dir3 = \a dir1 x \a dir2, where \a dir1 and \a
     *  dir2 are the emittance axes. The beam is made to match Twiss
     *  parameters \f$ \alpha_1 \f$ (a1), \f$ \beta_1 \f$ (b1), \f$
     *  \epsilon_{1,\mathrm{rms}} \f$ (e1) in the first projectional
     *  direction \a dir1 and \f$ \alpha_2 \f$ (a2), \f$ \beta_2 \f$
     *  (b2), \f$ \epsilon_{2,\mathrm{rms}} \f$ (e2) in the second
     *  projectional direction \a dir2.
     *
     *  The beam spread in the transverse phase space is made
     *  according to KV/hard-edged (Kapchinsky-Vladimirsky)
     *  distribution.
     */
    void add_3d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
					double E0,
					double a1, double b1, double e1,
					double a2, double b2, double e2,
					Vec3D c, Vec3D dir1, Vec3D dir2 );

    /*! \brief Add a 3d beam with defined waterbag emittance.
     *
     *  Adds a beam consisting of \a N particles and a total current
     *  of \a I (A). The particles are defined to have equal currents,
     *  charge \a q (in multiples of e) and mass \a m (u). The
     *  starting energy of the beam is \a E0 (eV) and the starting
     *  location \a c (center point). The beam propagates to the
     *  direction \a dir3 = \a dir1 x \a dir2, where \a dir1 and \a
     *  dir2 are the emittance axes. The beam is made to match Twiss
     *  parameters \f$ \alpha_1 \f$ (a1), \f$ \beta_1 \f$ (b1), \f$
     *  \epsilon_{1,\mathrm{rms}} \f$ (e1) in the first projectional
     *  direction \a dir1 and \f$ \alpha_2 \f$ (a2), \f$ \beta_2 \f$
     *  (b2), \f$ \epsilon_{2,\mathrm{rms}} \f$ (e2) in the second
     *  projectional direction \a dir2.
     *
     *  The beam spread in the transverse phase space is made
     *  according to waterbag distribution.
     */
    void add_3d_waterbag_beam_with_emittance( uint32_t N, double I, double q, double m,
					      double E0,
					      double a1, double b1, double e1,
					      double a2, double b2, double e2,
					      Vec3D c, Vec3D dir1, Vec3D dir2 );
    
    /*! \brief Add a 3d beam with defined gaussian emittance.
     *
     *  Adds a beam consisting of \a N particles and a total current
     *  of \a I (A). The particles are defined to have equal currents,
     *  charge \a q (in multiples of e) and mass \a m (u). The
     *  starting energy of the beam is \a E0 (eV) and the starting
     *  location \a c (center point). The beam propagates to the
     *  direction \a dir3 = \a dir1 x \a dir2, where \a dir1 and \a
     *  dir2 are the emittance axes. The beam is made to match Twiss
     *  parameters \f$ \alpha_1 \f$ (a1), \f$ \beta_1 \f$ (b1), \f$
     *  \epsilon_{1,\mathrm{rms}} \f$ (e1) in the first projectional
     *  direction \a dir1 and \f$ \alpha_2 \f$ (a2), \f$ \beta_2 \f$
     *  (b2), \f$ \epsilon_{2,\mathrm{rms}} \f$ (e2) in the second
     *  projectional direction \a dir2.
     *
     *  The beam spread in the transverse phase space is made
     *  according to Gaussian distribution.
     */
    void add_3d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
					      double E0, 
					      double a1, double b1, double e1,
					      double a2, double b2, double e2,
					      Vec3D c, Vec3D dir1, Vec3D dir2 );


/* ************************************** *
 * Information and queries                *
 * ************************************** */


    /*! \brief Get trajectory diagnostic on a plane.
     *
     *  Builds trajectory diagnostics on a plane defined in three
     *  dimensional space by center point vector \a c and two basis
     *  vectors \a o and \a p. The vectors \a o and \a p are
     *  normalized and p is adjusted to be orthogonal to \a o. This is
     *  done by calculating \a q = \a cross(o,p) and \a p = \a
     *  cross(q,o). 
     *
     *  The diagnostic gathered is defined by \a diagnostics. Data is
     *  returned in \a tdata.
     *
     *  Valid diagnostic types are DIAG_T, DIAG_X, DIAG_VX, DIAG_Y,
     *  DIAG_VY, DIAG_Z, DIAG_VZ, DIAG_O, DIAG_VO, DIAG_P, DIAG_VP,
     *  DIAG_Q, DIAG_VQ, DIAG_OP, DIAG_PP, DIAG_CURR, DIAG_EK and
     *  DIAG_QM, DIAG_MASS and DIAG_CHARGE.
     */
    void trajectories_at_free_plane( TrajectoryDiagnosticData &tdata, 
				     const Vec3D &c, const Vec3D &o, const Vec3D &p,
				     const std::vector<trajectory_diagnostic_e> &diagnostics ) const;


    /*! \brief Export particle data as Path Manager data
     *
     *  Makes trajectory diagnostics on a plane defined in three
     *  dimensional space by center point vector \a c and two basis
     *  vectors \a o and \a p. The vectors \a o and \a p are
     *  normalized and p is adjusted to be orthogonal to \a o. This is
     *  done by calculating \a q = \a cross(o,p) and \a p = \a
     *  cross(q,o).
     *
     *  The particle properties on the diagnostics plane are written
     *  to file \a filename in the format required by the Path Manager
     *  program. Reference particle for the output has energy \a ref_E
     *  (in eV), charge state \a ref_q (in electron charges) and mass
     *  \a ref_m (in atomic mass units).
     *
     *  The file contains a header with information about the
     *  reference particle. After the header, the file has a line for
     *  each particle, which contain the following columns: \a x, \a
     *  x', \a y, \a y', \a phase, \a (p_z-p_ref)/p_ref, \a 0-flag, \a
     *  q, \a m. The particles are mirrored according the particle
     *  database mirroring settings. The phase is randomized for each
     *  particle between 0 and 2 pi corresponding to continuous beam.
     */
    void export_path_manager_data( const std::string &filename, 
				   double ref_E, double ref_q, double ref_m, 
				   const Vec3D &c, const Vec3D &o, const Vec3D &p ) const;

/* ************************************** *
 * Debugging, plotting and saving         *
 * ************************************** */

    /*! \brief Saves data to a new file \a filename.
     */
    virtual void save( const std::string &filename ) const;

    /*! \brief Saves data to stream.
     */
    virtual void save( std::ostream &s ) const;

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;
};

#endif

