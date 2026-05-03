/*! \file particledatabaseimp.hpp
 *  \brief Particle database implementation
 */

/* Copyright (c) 2005-2013,2015 Taneli Kalvas. All rights reserved.
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

#ifndef PARTICLEDATABASEIMP_HPP
#define PARTICLEDATABASEIMP_HPP 1


#include "ibsimu.hpp"
#include "statusprint.hpp"
#include "timer.hpp"
#include "file.hpp"
#include "particles.hpp"
#include "particleiterator.hpp"
#include "trajectorydiagnostics.hpp"


class ParticleDataBaseImp {

protected:

    const Geometry             &_geom;
    double                     _epsabs;       /*!< \brief Absolute error limit for calculation. */
    double                     _epsrel;       /*!< \brief Relative error limit for calculation. */
    trajectory_interpolation_e _intrp;        /*!< \brief Trajectory interpolation type. */
    scharge_deposition_e       _scharge_dep;  /*!< \brief Space charge deposition type. */
    uint32_t                   _maxsteps;     /*!< \brief Maximum number of steps to calculate. */
    double                     _maxt;         /*!< \brief Maximum particle time in simulation. */
    bool                       _save_points;  /*!< \brief Save all points? */
    uint32_t                   _trajdiv;      /*!< \brief Divisor for saved trajectories,
					       * if 3, every third trajectory is saved. */
    bool                       _mirror[6];    /*!< \brief Boundary particle mirroring. */

    double                     _rhosum;       /*!< \brief Sum of space charge density in defined beams (C/m3). */

    ParticleStatistics         _stat;         /*!< \brief Particle statistics. */

    uint32_t                   _iteration;    /*!< \brief Iteration number. */    
    bool                       _relativistic; /*!< \brief Relativistic particle iteration. */
    bool                       _surface_collision; /*!< \brief Surface collision model. */

    const CallbackFunctorD_V  *_bsup_cb;      /*!< \brief Location dependent magnetic field suppression. */
    TrajectoryHandlerCallback *_thand_cb;     /*!< \brief Trajectory handler callback. */
    TrajectoryEndCallback     *_tend_cb;      /*!< \brief Trajectory collision callback. */
    TrajectorySurfaceCollisionCallback *_tsur_cb;    /*!< \brief Trajectory surface collision callback. */
    ParticleDataBase          *_pdb;          /*!< \brief Particle database pointer. */

    ParticleDataBaseImp( ParticleDataBase *pdb, const Geometry &geom );

    ParticleDataBaseImp( ParticleDataBase *pdb, std::istream &s, const Geometry &geom );

    ParticleDataBaseImp( const ParticleDataBaseImp &pdb );

    const ParticleDataBaseImp &operator=( const ParticleDataBaseImp &pdb );

    /*! \brief Convert energy to velocity.
     *
     *  Energy \a E given in Joules and mass \a m in kg.
     */
    static double energy_to_velocity( double E, double m );
    
    void save( std::ostream &os ) const;

public:

    virtual ~ParticleDataBaseImp();

    void set_accuracy( double epsabs, double epsrel );

    void set_bfield_suppression( const CallbackFunctorD_V *functor );

    void set_trajectory_handler_callback( TrajectoryHandlerCallback *thand_cb );

    void set_trajectory_end_callback( TrajectoryEndCallback *tend_cb );

    void set_trajectory_surface_collision_callback( TrajectorySurfaceCollisionCallback *tsur_cb );

    void set_relativistic( bool enable );

    void set_surface_collision( bool surface_collision );

    void set_polyint( bool polyint );
    
    bool get_polyint( void ) const;
    
    void set_trajectory_interpolation( trajectory_interpolation_e intrp );

    trajectory_interpolation_e get_trajectory_interpolation( void ) const;

    void set_scharge_deposition( scharge_deposition_e type );

    scharge_deposition_e get_scharge_deposition( void ) const;

    void set_max_steps( uint32_t maxsteps );

    void set_max_time( double maxt );

    void set_save_all_points( bool save_points );

    void set_save_trajectories( uint32_t div );

    uint32_t get_save_trajectories( void ) const;

    void set_mirror( const bool mirror[6] );

    void get_mirror( bool mirror[6] ) const;

    int get_iteration_number( void ) const;

    double get_rhosum( void ) const;

    void set_rhosum( double rhosum );
    
    const ParticleStatistics &get_statistics( void ) const;

    virtual geom_mode_e geom_mode() const = 0;

    virtual size_t size( void ) const = 0;

    virtual size_t traj_size( uint32_t i ) const = 0;

    virtual double traj_length( uint32_t i ) const = 0;

    virtual void trajectory_point( double &t, Vec3D &loc, Vec3D &vel, uint32_t i, uint32_t j ) const = 0;

    virtual void trajectories_at_plane( TrajectoryDiagnosticData &tdata, 
					coordinate_axis_e axis,
					double val,
					const std::vector<trajectory_diagnostic_e> &diagnostics ) const = 0;

    virtual void clear( void ) = 0;

    virtual void clear_trajectories( void ) = 0;

    virtual void clear_trajectory( size_t a ) = 0;

    virtual void reset_trajectories( void ) = 0;

    virtual void reset_trajectory( size_t a ) = 0;

    virtual void reserve( size_t size ) = 0;

    virtual void build_trajectory_density_field( MeshScalarField &tdens ) const = 0;

    virtual void iterate_trajectories( MeshScalarField &scharge, const VectorField &efield, 
				       const VectorField &bfield ) = 0;

    virtual void step_particles( MeshScalarField &scharge, const VectorField &efield, 
				 const VectorField &bfield, double dt ) = 0;

    void debug_print( std::ostream &os ) const;
};


template<class PP> class ParticleDataBasePPImp : public ParticleDataBaseImp {

    /*! \brief Add requested diagnostics to \a tdata from particle point \a x.
     */
    static void add_diagnostics( TrajectoryDiagnosticData &tdata, const PP &x, 
				 const Particle<PP> &p, int crd, int index ) {
	//std::cout << "add_diagnostics():\n";
	for( size_t a = 0; a < tdata.diag_size(); a++ ) {
	    //std::cout << "  diagnostic[" << a << "] = " << tdata.diagnostic(a) << "\n";
	    
	    double data = 0.0;
	    switch( tdata.diagnostic( a ) ) {
	    case DIAG_NONE:
		data = 0.0;
		break;
	    case DIAG_T:
		data = x[0];
		break;
	    case DIAG_X:
		data = x[1];
		break;
	    case DIAG_VX:
		data = x[2];
		break;
	    case DIAG_Y:
	    case DIAG_R:
		data = x[3];
		break;
	    case DIAG_VY:
	    case DIAG_VR:
		data = x[4];
		break;
	    case DIAG_Z:
		data = x[5];
		break;
	    case DIAG_VZ:
		data = x[6];
		break;
	    case DIAG_W:
		data = x[5];
		break;
	    case DIAG_VTHETA:
		data = x[5]*x[3];
		break;
	    case DIAG_XP:
		data = x[2]/x[2*crd+2];
		break;
	    case DIAG_YP:
	    case DIAG_RP:
		data = x[4]/x[2*crd+2];
		break;
	    case DIAG_AP:
		data = x[3]*x[5]/x[2*crd+2];
		break;
	    case DIAG_ZP:
		data = x[6]/x[2*crd+2];
		break;
	    case DIAG_CURR:
		data = p.IQ();
		break;
	    case DIAG_QM:
		data = (p.q()/CHARGE_E) / (p.m()/MASS_U);
		break;
	    case DIAG_CHARGE:
		data = p.q()/CHARGE_E;
		break;
	    case DIAG_MASS:
		data = p.m()/MASS_U;
		break;
	    case DIAG_EK:
	    {
		double beta = x.velocity().norm2()/SPEED_C;
		if( beta < 1.0e-4 )
		    data = 0.5*p.m()*x.velocity().ssqr()/CHARGE_E;
		else if( beta >= 1.0 )
		    throw( ErrorUnimplemented( ERROR_LOCATION, "particle velocity over light speed" ) );
		else {
		    double gamma = 1.0 / sqrt( 1.0 - beta*beta );
		    double Ek = p.m()*SPEED_C2*( gamma - 1.0 );
		    data = Ek/CHARGE_E;
		}
		break;
	    }
	    case DIAG_NO:
		data = index;
		break;
	    default:
		throw( ErrorUnimplemented( ERROR_LOCATION ) );
		break;
	    }
	    //std::cout << "  adding data = " << data << "\n";
	    tdata.add_data( a, data );
	}
    }

protected:

    std::vector<Particle<PP> *>                        _particles; /*!< \brief Particles. */
    Scheduler<ParticleIterator<PP>,Particle<PP>,Error> _scheduler; /*!< \brief Scheduler for solver. */

    /*! \brief Constructor, using API pdb
     */
    ParticleDataBasePPImp( ParticleDataBase *pdb, const Geometry &geom )
	: ParticleDataBaseImp(pdb,geom), _scheduler(_particles) {
	if( _geom.geom_mode() != PP::geom_mode() )
	    throw( Error( ERROR_LOCATION, "Differing geometry modes" ) );
    }

    /*! \brief Constructor from stream, using API pdb
     */
    ParticleDataBasePPImp( ParticleDataBase *pdb, std::istream &s, const Geometry &geom ) 
	: ParticleDataBaseImp(pdb,s,geom), _scheduler(_particles) {

	uint32_t N = read_int32( s );
	ibsimu.message( 1 ) << "Reading " << N << " particles.\n";
	_particles.reserve( N );
	for( uint32_t a = 0; a < N; a++ )
	    _particles.push_back( new Particle<PP>( s ) );

	ibsimu.dec_indent();
    }

    /*! \brief Copy constructor
     */
    ParticleDataBasePPImp( const ParticleDataBasePPImp &pdb ) 
	: ParticleDataBaseImp(pdb), _scheduler(_particles) {

	_particles.reserve( pdb._particles.size() );
	for( size_t a = 0; a < pdb._particles.size(); a++ )
	    _particles.push_back( new Particle<PP>( *pdb._particles[a] ) );
    }

    /*! \brief Copy assignment operator
     */
    const ParticleDataBasePPImp &operator=( const ParticleDataBasePPImp &pdb ) {

	for( size_t a; a < pdb._particles.size(); a++ )
	    delete _particles[a];
	_particles.clear();
	_particles.reserve( pdb._particles.size() );
	for( size_t a = 0; a < pdb._particles.size(); a++ )
	    _particles.push_back( new Particle<PP>( *pdb._particles[a] ) );
	return( *this );
    }

public:

    virtual ~ParticleDataBasePPImp() {
	// Delete particles
	for( size_t a = 0; a < _particles.size(); a++ )
	    delete( _particles[a] );
    }

    virtual geom_mode_e geom_mode() const { 
	return( PP::geom_mode() ); 
    }

    virtual size_t size( void ) const { 
	return( _particles.size() ); 
    }

    Particle<PP> &particle( uint32_t i ) { 
	return( *_particles[i] ); 
    }

    const Particle<PP> &particle( uint32_t i ) const { 
	return( *_particles[i] ); 
    }
    
    virtual double traj_length( uint32_t i ) const {

	size_t N = _particles[i]->traj_size();
	double len = 0.0;
	if( N < 2 )
	    return( 0.0 );
	Vec3D x1 = _particles[i]->traj(0).location();
	for( size_t b = 1; b < N; b++ ) {
	    Vec3D x2 = _particles[i]->traj(b).location();
	    len += norm2(x2-x1);
	    x1 = x2;
	}

	return( len );
    }

    virtual size_t traj_size( uint32_t i ) const { 
	return( _particles[i]->traj_size() ); 
    }

    const PP &trajectory_point( uint32_t i, uint32_t j ) const {
	return( _particles[i]->traj(j) );
    }

    virtual void trajectory_point( double &t, Vec3D &loc, Vec3D &vel, uint32_t i, uint32_t j ) const {
	PP x = _particles[i]->traj(j); 
	t = x[0];
	loc = x.location();
	vel = x.velocity();
    }

    virtual void trajectories_at_plane( std::vector< Particle<PP> > &tdata,
					coordinate_axis_e axis,
					double val ) const {

	ibsimu.message( 1 ) << "Making trajectory diagnostics at " 
			    << coordinate_axis_string[axis] << " = " << val << "\n";
	ibsimu.inc_indent();

	// Check query
	switch( PP::geom_mode() ) {
	case MODE_1D:
	    throw( Error( ERROR_LOCATION, "unsupported dimension number" ) );
	    break;
	case MODE_2D:
	    if( axis == AXIS_R || axis == AXIS_Z )
		throw( Error( ERROR_LOCATION, "nonexistent axis" ) );
	    break;
	case MODE_CYL:
	    if( axis == AXIS_Y || axis == AXIS_Z )
		throw( Error( ERROR_LOCATION, "nonexistent axis" ) );
	    break;
	case MODE_3D:
	    if( axis == AXIS_R )
		throw( Error( ERROR_LOCATION, "nonexistent axis" ) );
	    break;
	default:
	    throw( Error( ERROR_LOCATION, "unsupported dimension number" ) );
	}

	// Prepare output vector
	tdata.clear();

	// Set coordinate index
	int crd;
	switch( axis ) {
	case AXIS_X:
	    crd = 0;
	    break;
	case AXIS_Y:
	case AXIS_R:
	    crd = 1;
	    break;
	case AXIS_Z:
	    crd = 2;
	    break;
	default:
	    throw( Error( ERROR_LOCATION, "unsupported axis" ) );
	}

	// Scan through particle trajectory points
	double Isum = 0.0;
	std::vector<PP> intsc;
	for( size_t a = 0; a < _particles.size(); a++ ) {
	    size_t N = _particles[a]->traj_size();
	    if( N < 2 )
		continue;
	    PP x1 = _particles[a]->traj(0);
	    for( size_t b = 1; b < N; b++ ) {
		PP x2 = _particles[a]->traj(b);
		intsc.clear();
		size_t nintsc;
		if( b == 1 )
		    nintsc = PP::trajectory_intersections_at_plane( intsc, crd, val, x1, x2, -1 );
		else if( b == N-1 )
		    nintsc = PP::trajectory_intersections_at_plane( intsc, crd, val, x1, x2, +1 );
		else
		    nintsc = PP::trajectory_intersections_at_plane( intsc, crd, val, x1, x2, 0 );
		for( size_t c = 0; c < nintsc; c++ ) {
		    Isum += _particles[a]->IQ();
		    tdata.push_back( Particle<PP>( _particles[a]->IQ(), _particles[a]->q(),
						   _particles[a]->m(), intsc[c] ) );
		}

		x1 = x2;
	    }
	}

	ibsimu.message( 1 ) << "number of trajectories = " << tdata.size() << "\n";
	if( PP::geom_mode() == MODE_2D )
	    ibsimu.message( 1 ) << "total current = " << Isum << " A/m\n";
	else
	    ibsimu.message( 1 ) << "total current = " << Isum << " A\n";
	ibsimu.dec_indent();	
    }

    virtual void trajectories_at_plane( TrajectoryDiagnosticData &tdata, 
					coordinate_axis_e axis,
					double val,
					const std::vector<trajectory_diagnostic_e> &diagnostics ) const {

	ibsimu.message( 1 ) << "Making trajectory diagnostics at " 
			    << coordinate_axis_string[axis] << " = " << val << "\n";
	ibsimu.inc_indent();

	// Check query
	switch( PP::geom_mode() ) {
	case MODE_1D:
	    throw( Error( ERROR_LOCATION, "unsupported dimension number" ) );
	    break;
	case MODE_2D:
	    if( axis == AXIS_R || axis == AXIS_Z )
		throw( Error( ERROR_LOCATION, "nonexistent axis" ) );
	    break;
	case MODE_CYL:
	    if( axis == AXIS_Y || axis == AXIS_Z )
		throw( Error( ERROR_LOCATION, "nonexistent axis" ) );
	    break;
	case MODE_3D:
	    if( axis == AXIS_R )
		throw( Error( ERROR_LOCATION, "nonexistent axis" ) );
	    break;
	default:
	    throw( Error( ERROR_LOCATION, "unsupported dimension number" ) );
	}

	// Check diagnostics query validity
	for( size_t a = 0; a < diagnostics.size(); a++ ) {
	    if( diagnostics[a] == DIAG_NONE )
		throw( Error( ERROR_LOCATION, "invalid diagnostics query \'DIAG_NONE\'" ) );
	    else if( PP::geom_mode() != MODE_CYL && (diagnostics[a] == DIAG_R ||
						     diagnostics[a] == DIAG_VR ||
						     diagnostics[a] == DIAG_RP ||
						     diagnostics[a] == DIAG_W ||
						     diagnostics[a] == DIAG_VTHETA ||
						     diagnostics[a] == DIAG_AP) )
		throw( Error( ERROR_LOCATION, "invalid diagnostics query for geometry type" ) );
	    else if( PP::geom_mode() != MODE_3D && (diagnostics[a] == DIAG_Z ||
						    diagnostics[a] == DIAG_VZ ||
						    diagnostics[a] == DIAG_ZP) )
		throw( Error( ERROR_LOCATION, "invalid diagnostics query for geometry type" ) );
	    else if( diagnostics[a] == DIAG_O ||
		     diagnostics[a] == DIAG_VO ||
		     diagnostics[a] == DIAG_OP ||
		     diagnostics[a] == DIAG_P ||
		     diagnostics[a] == DIAG_VP ||
		     diagnostics[a] == DIAG_PP ||
		     diagnostics[a] == DIAG_Q ||
		     diagnostics[a] == DIAG_VQ )
		throw( Error( ERROR_LOCATION, "invalid diagnostics query for trajectories_at_plane()\n"
			      "use trajectories_at_free_plane()" ) );
	}

	// Prepare output vector
	tdata.clear();
	for( size_t a = 0; a < diagnostics.size(); a++ ) {
	    tdata.add_data_column( diagnostics[a] );
	}
	
	// Set coordinate index
	int crd;
	switch( axis ) {
	case AXIS_X:
	    crd = 0;
	    break;
	case AXIS_Y:
	case AXIS_R:
	    crd = 1;
	    break;
	case AXIS_Z:
	    crd = 2;
	    break;
	default:
	    throw( Error( ERROR_LOCATION, "unsupported axis" ) );
	}

	// Scan through particle trajectory points
	double Isum = 0.0;
	std::vector<PP> intsc;
	for( size_t a = 0; a < _particles.size(); a++ ) {
	    size_t N = _particles[a]->traj_size();
	    if( N < 2 )
		continue;
	    PP x1 = _particles[a]->traj(0);
	    for( size_t b = 1; b < N; b++ ) {
		PP x2 = _particles[a]->traj(b);
		intsc.clear();
		size_t nintsc;
		if( b == 1 )
		    nintsc = PP::trajectory_intersections_at_plane( intsc, crd, val, x1, x2, -1 );
		else if( b == N-1 )
		    nintsc = PP::trajectory_intersections_at_plane( intsc, crd, val, x1, x2, +1 );
		else
		    nintsc = PP::trajectory_intersections_at_plane( intsc, crd, val, x1, x2, 0 );
		for( size_t c = 0; c < nintsc; c++ ) {
		    Isum += _particles[a]->IQ();
		    add_diagnostics( tdata, intsc[c], *_particles[a], crd, a );
		}

		x1 = x2;
	    }
	}

	ibsimu.message( 1 ) << "number of trajectories = " << tdata.traj_size() << "\n";
	if( PP::geom_mode() == MODE_2D )
	    ibsimu.message( 1 ) << "total current = " << Isum << " A/m\n";
	else
	    ibsimu.message( 1 ) << "total current = " << Isum << " A\n";
	ibsimu.dec_indent();
    }

    virtual void clear( void ) { 
	for( size_t a = 0; a < _particles.size(); a++ )
	    delete( _particles[a] );
	_particles.clear();
	_rhosum = 0.0;
    }

    virtual void clear_trajectories( void ) {
	for( uint32_t a = 0; a < _particles.size(); a++ )
	    _particles[a]->clear_trajectory();
    }

    virtual void clear_trajectory( size_t a ) {
	if( a >= _particles.size() )
	    throw( ErrorRange( ERROR_LOCATION, a, _particles.size() ) );	    
	_particles[a]->clear_trajectory();
    }

    virtual void reset_trajectories( void ) {
	for( uint32_t a = 0; a < _particles.size(); a++ )
	    _particles[a]->reset_trajectory();
    }

    virtual void reset_trajectory( size_t a ) {
	if( a >= _particles.size() )
	    throw( ErrorRange( ERROR_LOCATION, a, _particles.size() ) );	    
	_particles[a]->reset_trajectory();
    }

    virtual void reserve( size_t size ) { 
	_particles.reserve( size ); 
    }

    void add_particle( const Particle<PP> &pp ) {
	_scheduler.lock_mutex();
	try {
	    _particles.push_back( new Particle<PP>( pp ) );
	} catch( std::bad_alloc ) {
	    throw( ErrorNoMem( ERROR_LOCATION, "Out of memory adding particle" ) );
	}
	_scheduler.unlock_mutex();
    }

    void add_particle( double IQ, double q, double m, const PP &x ) {	
	add_particle( Particle<PP>( IQ, CHARGE_E*q, MASS_U*m, x ) );
    }

    virtual void iterate_trajectories( MeshScalarField &scharge, const VectorField &efield, 
				       const VectorField &bfield ) {

	Timer t;
	ibsimu.message( 1 ) << "Calculating particle trajectories\n";
	ibsimu.inc_indent();
	if( _surface_collision )
	    ibsimu.message( 1 ) << "Using surface collision model\n";
	else
	    ibsimu.message( 1 ) << "Using solid collision model\n";
	if( _relativistic )
	    ibsimu.message( 1 ) << "Using relativistic iterator\n";
	else
	    ibsimu.message( 1 ) << "Using non-relativistic iterator\n";
	_iteration++;

	StatusPrint sp;
	if( ibsimu.output_is_cout() ) {
	    std::stringstream ss;
	    ss << "  " << "0 / " << _particles.size();
	    sp.print( ss.str() );
	}

	// Check geometry mode
	if( _geom.geom_mode() != PP::geom_mode() )
	    throw( Error( ERROR_LOCATION, "Differing geometry modes" ) );

	// Clear space charge
	scharge.clear();

	// Reset statistics
	_stat.reset( _geom.number_of_boundaries() );

	// Check number of particles
	if( _particles.size() == 0 ) {
	    ibsimu.message( 1 ) << "no particles to calculate\n";
	    ibsimu.dec_indent();
	    return;
	}

	// Make solvers
	pthread_mutex_t                      scharge_mutex = PTHREAD_MUTEX_INITIALIZER;
	std::vector<ParticleIterator<PP> *>  iterators;
	for( uint32_t a = 0; a < ibsimu.get_thread_count(); a++ ) {

	    iterators.push_back( new ParticleIterator<PP>( PARTICLE_ITERATOR_ADAPTIVE, _epsabs, _epsrel, 
							   _intrp, _scharge_dep, _maxsteps, _maxt, 
							   _save_points, _trajdiv, _mirror, &scharge, 
							   &scharge_mutex, &efield, &bfield, &_geom ) );
	    // iterators[a]->set_trajectory_handler_callback( _thand_cb );
      iterators[a]->set_trajectory_handler_callback( _thand_cb, _pdb );
	    iterators[a]->set_trajectory_end_callback( _tend_cb, _pdb );
	    iterators[a]->set_trajectory_surface_collision_callback( _tsur_cb );
	    iterators[a]->set_bfield_suppression_callback( _bsup_cb );
	    iterators[a]->set_relativistic( _relativistic );
	    iterators[a]->set_surface_collision( _surface_collision );
	}

	// Run scheduler
	_scheduler.run( iterators );

	// Print statistics
	if( ibsimu.output_is_cout() ) {
	    while( !_scheduler.wait_finish() ) {
		std::stringstream ss;
		ss << "  " << _scheduler.get_solved_count() << " / " 
		   << _scheduler.get_problem_count();
		sp.print( ss.str() );
	    }
	}

	// Finish scheduler
	_scheduler.finish();

	// Print final statistics
	if( ibsimu.output_is_cout() ) {
	    std::stringstream ss;
	    ss << "  " << _scheduler.get_solved_count() << " / " 
	       << _scheduler.get_problem_count() << " Done\n";
	    sp.print( ss.str(), true );
	}

	if( _scheduler.is_error() ) {
	    // Throw the error
	    std::vector<Error> err;
	    std::vector<int32_t> part;
	    _scheduler.get_errors( err, part );
	    throw( err[0] );
	}

	// Collect statistics. Free all allocated memory.
	for( uint32_t a = 0; a < ibsimu.get_thread_count(); a++ ) {
	    ParticleStatistics stat = iterators[a]->get_statistics();
	    _stat += stat;
	    delete iterators[a];
	}

	if( _scharge_dep == SCHARGE_DEPOSITION_LINEAR )
	    scharge_finalize_linear( scharge );
	else
	    scharge_finalize_pic( scharge );
	
	t.stop();
	ibsimu.message( 1 ) << "Particle histories (" << _particles.size() << " total):\n";
	ibsimu.message( 1 ) << "  flown = " << _stat.bound_collisions() << "\n";
	ibsimu.message( 1 ) << "  time limited = " << _stat.end_time() << "\n";
	ibsimu.message( 1 ) << "  step count limited = " << _stat.end_step() << "\n";
	ibsimu.message( 1 ) << "  bad definitions = " << _stat.end_baddef() << "\n";
	for( size_t a = 1; a <= _stat.number_of_boundaries(); a++ ) {
	    ibsimu.message( 1 ) << "  beam to boundary " << a << " = " << _stat.bound_current(a)
				<< " " << PP::IQ_unit() << " (" << _stat.bound_collisions(a) << " particles)" << "\n";
	}
	ibsimu.message( 1 ) << "  total steps = " << _stat.sum_steps() << "\n";
	ibsimu.message( 1 ) << "  steps per particle (ave) = " << 
	    _stat.sum_steps()/(double)_particles.size() << "\n";
	ibsimu.message( 1 ) << "time used = " << t << "\n";
	ibsimu.dec_indent();
    }

    virtual void step_particles( MeshScalarField &scharge, const VectorField &efield, 
				 const VectorField &bfield, double dt ) {

	throw( ErrorUnimplemented( ERROR_LOCATION ) );
    }


    void save( std::ostream &os ) const {

	ParticleDataBaseImp::save( os );

	write_int32( os, _particles.size() );
	for( uint32_t a = 0; a < _particles.size(); a++ )
	    _particles[a]->save( os );
    }


    void debug_print( std::ostream &os ) const {
	ParticleDataBaseImp::debug_print( os );

	for( uint32_t a = 0; a < _particles.size(); a++ ) {
	    os << "Particle " << a << ":\n";
	    _particles[a]->debug_print( os );
	}
    }

};


/*! \brief ParticleDataBase2D implementation.
 *
 *  Not part of the public interface of the library.
 */
class ParticleDataBase2DImp : public ParticleDataBasePPImp<ParticleP2D> {

    void add_tdens_from_segment( MeshScalarField &tdens, double IQ,
				 ParticleP2D &x1, ParticleP2D &x2 ) const;

public:

    ParticleDataBase2DImp( ParticleDataBase *pdb, const Geometry &geom );

    ParticleDataBase2DImp( const ParticleDataBase2DImp &pdb );

    ParticleDataBase2DImp( ParticleDataBase *pdb, std::istream &s, const Geometry &geom );

    virtual ~ParticleDataBase2DImp();

    const ParticleDataBase2DImp &operator=( const ParticleDataBase2DImp &pdb );

    virtual void build_trajectory_density_field( MeshScalarField &tdens ) const;

    void add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
				    double v, double dvp, double dvt, 
				    double x1, double y1, double x2, double y2 );

    void add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
				  double E, double Tp, double Tt, 
				  double x1, double y1, double x2, double y2 );

    void add_2d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
					double a, double b, double e,
					double Ex, double x0, double y0 );

    void add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
					      double a, double b, double e,
					      double Ex, double x0, double y0 );

    void save( std::ostream &os ) const;

    void debug_print( std::ostream &os ) const;
};



/*! \brief ParticleDataBaseCyl implementation.
 *
 *  Not part of the public interface of the library.
 */
class ParticleDataBaseCylImp : public ParticleDataBasePPImp<ParticlePCyl> {

    void add_tdens_from_segment( MeshScalarField &tdens, double IQ,
				 ParticlePCyl &x1, ParticlePCyl &x2 ) const;

    static uint32_t bisect_cumulative_array( const std::vector<double> &cum, double x );

public:

    ParticleDataBaseCylImp( ParticleDataBase *pdb, const Geometry &geom );

    ParticleDataBaseCylImp( const ParticleDataBaseCylImp &pdb );
    
    ParticleDataBaseCylImp( ParticleDataBase *pdb, std::istream &s, const Geometry &geom );

    virtual ~ParticleDataBaseCylImp();

    const ParticleDataBaseCylImp &operator=( const ParticleDataBaseCylImp &pdb );

    virtual void build_trajectory_density_field( MeshScalarField &tdens ) const;

    void add_2d_beam_with_total_energy( uint32_t N, double J, double q, double m, 
					double Etot, const ScalarField &epot, 
					double Tp, double Tt, 
					double x1, double y1, double x2, double y2 );

    void add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
				  double E, double Tp, double Tt, 
				  double x1, double y1, double x2, double y2 );

    void add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
				    double v, double dvp, double dvt, 
				    double x1, double y1, double x2, double y2 );

    void add_2d_full_gaussian_beam( uint32_t N, double I, double q, double m,
				    double Ex, double Tp, double Tt, 
				    double x0, double dr );

    void add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
					      double a, double b, double e,
					      double Ex, double x0 );

    void export_path_manager_data( const std::string &filename, 
				   double ref_E, double ref_q, double ref_m, 
				   double val, uint32_t Np ) const;

    void save( std::ostream &os ) const;

    void debug_print( std::ostream &os ) const;
};


/*! \brief ParticleDataBase3D implementation.
 *
 *  Not part of the public interface of the library.
 */
class ParticleDataBase3DImp : public ParticleDataBasePPImp<ParticleP3D> {

    void add_tdens_from_segment( MeshScalarField &tdens, double IQ,
				 ParticleP3D &x1, ParticleP3D &x2 ) const;

    bool free_plane_mirror_enabled( uint32_t axism ) const;
    ParticleP3D free_plane_mirror( const ParticleP3D &p, uint32_t axism ) const;
public:

    ParticleDataBase3DImp( ParticleDataBase *pdb, const Geometry &geom );

    ParticleDataBase3DImp( const ParticleDataBase3DImp &pdb );

    ParticleDataBase3DImp( ParticleDataBase *pdb, std::istream &s, const Geometry &geom );

    virtual ~ParticleDataBase3DImp();

    const ParticleDataBase3DImp &operator=( const ParticleDataBase3DImp &pdb );

    virtual void build_trajectory_density_field( MeshScalarField &tdens ) const;

    void add_cylindrical_beam_with_total_energy( uint32_t N, double J, double q, double m, 
						 double Etot, const ScalarField &epot, 
						 double Tp, double Tt, Vec3D c, 
						 Vec3D dir1, Vec3D dir2, double r );

    void add_cylindrical_beam_with_energy( uint32_t N, double J, double q, double m, 	
					   double E, double Tp, double Tt, Vec3D c, 
					   Vec3D dir1, Vec3D dir2, double r );

    void add_cylindrical_beam_with_velocity( uint32_t N, double J, double q, double m, 
					     double v, double dvp, double dvt, Vec3D c, 
					     Vec3D dir1, Vec3D dir2, double r );

    void add_rectangular_beam_with_energy( uint32_t N, double J, double q, double m, 
					   double E, double Tp, double Tt, Vec3D c, 
					   Vec3D dir1, Vec3D dir2, double size1, double size2 );

    void add_rectangular_beam_with_velocity( uint32_t N, double J, double q, double m, 
					     double v, double dvp, double dvt, Vec3D c, 
					     Vec3D dir1, Vec3D dir2, double size1, double size2 );

    void add_3d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
					double E0,
					double a1, double b1, double e1,
					double a2, double b2, double e2,
					Vec3D c, Vec3D dir1, Vec3D dir2 );

    void add_3d_waterbag_beam_with_emittance( uint32_t N, double I, double q, double m,
					      double E0,
					      double a1, double b1, double e1,
					      double a2, double b2, double e2,
					      Vec3D c, Vec3D dir1, Vec3D dir2 );

    void add_3d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
					      double E0, 
					      double a1, double b1, double e1,
					      double a2, double b2, double e2,
					      Vec3D c, Vec3D dir1, Vec3D dir2 );

    void add_3d_plane_cosine_beam_with_partps( uint32_t N, double I, double m, double T, 
                          Vec3D c, Vec3D dir1, Vec3D dir2, double size1, double size2 );

    void add_3d_cylindrical_cosine_beam_with_partps( uint32_t N, double I, double m, double T, 
                          Vec3D c, Vec3D dir1, Vec3D dir2, double r );

    void add_3d_spherical_beam_with_energy( uint32_t N, double I, double q, double m, 
             double E, double Tp, double Tt, Vec3D cent, Vec3D pole, double angl);

    void trajectories_at_free_plane( TrajectoryDiagnosticData &tdata, 
				     Vec3D c, Vec3D o, Vec3D p,
				     const std::vector<trajectory_diagnostic_e> &diagnostics ) const;

    void export_path_manager_data( const std::string &filename, 
				   double ref_E, double ref_q, double ref_m, 
				   Vec3D c, Vec3D o, Vec3D p ) const;    

    void save( std::ostream &os ) const;

    void debug_print( std::ostream &os ) const;
};


#endif

