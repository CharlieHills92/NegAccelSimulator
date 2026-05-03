/*! \file particledatabase.cpp
 *  \brief %Particle databases
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


#include "particledatabase.hpp"
#include "particledatabaseimp.hpp"



ParticleDataBase::ParticleDataBase()
    : _imp(0)
{
}


ParticleDataBase::ParticleDataBase( const ParticleDataBase &pdb )
    : _imp(0)
{
}


ParticleDataBase::~ParticleDataBase()
{
}


const ParticleDataBase &ParticleDataBase::operator=( const ParticleDataBase &pdb )
{
    return( *this );
}


void ParticleDataBase::set_implementation_pointer( class ParticleDataBaseImp *imp )
{
    _imp = imp;
}


void ParticleDataBase::set_accuracy( double epsabs, double epsrel )
{
    _imp->set_accuracy( epsabs, epsrel );
}


void ParticleDataBase::set_bfield_suppression( const CallbackFunctorD_V *functor )
{
    _imp->set_bfield_suppression( functor );
}


void ParticleDataBase::set_trajectory_handler_callback( TrajectoryHandlerCallback *thand_cb )
{
    _imp->set_trajectory_handler_callback( thand_cb );
}


void ParticleDataBase::set_trajectory_end_callback( TrajectoryEndCallback *tend_cb )
{
    _imp->set_trajectory_end_callback( tend_cb );
}


void ParticleDataBase::set_trajectory_surface_collision_callback( TrajectorySurfaceCollisionCallback *tsur_cb )
{
    _imp->set_trajectory_surface_collision_callback( tsur_cb );
}


void ParticleDataBase::set_relativistic( bool enable )
{
    _imp->set_relativistic( enable );
}


void ParticleDataBase::set_surface_collision( bool surface_collision )
{
    _imp->set_surface_collision( surface_collision );
}


void ParticleDataBase::set_polyint( bool polyint ) 
{
    _imp->set_polyint( polyint );
}


bool ParticleDataBase::get_polyint( void ) const 
{
    return( _imp->get_polyint() );
}


void ParticleDataBase::set_trajectory_interpolation( trajectory_interpolation_e intrp )
{
    _imp->set_trajectory_interpolation( intrp );
}


trajectory_interpolation_e ParticleDataBase::get_trajectory_interpolation( void ) const
{
    return( _imp->get_trajectory_interpolation() );
}


void ParticleDataBase::set_scharge_deposition( scharge_deposition_e type )
{
    _imp->set_scharge_deposition( type );
}


scharge_deposition_e ParticleDataBase::get_scharge_deposition( void ) const
{
    return( _imp->get_scharge_deposition() );
}


void ParticleDataBase::set_max_steps( uint32_t maxsteps ) 
{
    _imp->set_max_steps( maxsteps );
}


void ParticleDataBase::set_max_time( double maxt ) 
{
    _imp->set_max_time( maxt );
}


void ParticleDataBase::set_save_all_points( bool save_points )
{
    _imp->set_save_all_points( save_points );
}


void ParticleDataBase::set_save_trajectories( uint32_t div ) 
{
    _imp->set_save_trajectories( div );
}


uint32_t ParticleDataBase::get_save_trajectories( void ) const
{
    return( _imp->get_save_trajectories() );
}


void ParticleDataBase::set_mirror( const bool mirror[6] )
{
    _imp->set_mirror( mirror );
}


void ParticleDataBase::get_mirror( bool mirror[6] ) const 
{
    _imp->get_mirror( mirror );
}


int ParticleDataBase::get_iteration_number( void ) const 
{
    return( _imp->get_iteration_number() );
}


double ParticleDataBase::get_rhosum( void ) const
{
    return( _imp->get_rhosum() );
}


void ParticleDataBase::set_rhosum( double rhosum )
{
    _imp->set_rhosum( rhosum );
}


const ParticleStatistics &ParticleDataBase::get_statistics( void ) const
{
    return( _imp->get_statistics() );
}


geom_mode_e ParticleDataBase::geom_mode() const
{
    return( _imp->geom_mode() );
}


size_t ParticleDataBase::size( void ) const
{
    return( _imp->size() );
}


double ParticleDataBase::traj_length( uint32_t i ) const
{
    return( _imp->traj_length( i ) );
}


size_t ParticleDataBase::traj_size( uint32_t i ) const
{
    return( _imp->traj_size( i ) );
}


void ParticleDataBase::trajectory_point( double &t, Vec3D &loc, Vec3D &vel, uint32_t i, uint32_t j ) const
{
    _imp->trajectory_point( t, loc, vel, i, j );
}


void ParticleDataBase::trajectories_at_plane( TrajectoryDiagnosticData &tdata, 
					      coordinate_axis_e axis,
					      double val,
					      const std::vector<trajectory_diagnostic_e> &diagnostics ) const
{
    _imp->trajectories_at_plane( tdata, axis, val, diagnostics );
}


void ParticleDataBase::clear( void )
{
    _imp->clear();
}


void ParticleDataBase::clear_trajectories( void )
{
    _imp->clear_trajectories();
}


void ParticleDataBase::clear_trajectory( size_t a )
{
    _imp->clear_trajectory( a );
}


void ParticleDataBase::reset_trajectories( void )
{
    _imp->reset_trajectories();
}


void ParticleDataBase::reset_trajectory( size_t a )
{
    _imp->reset_trajectory( a );
}


void ParticleDataBase::reserve( size_t size )
{
    _imp->reserve( size );
}


void ParticleDataBase::build_trajectory_density_field( MeshScalarField &tdens ) const
{
    _imp->build_trajectory_density_field( tdens );
}


void ParticleDataBase::iterate_trajectories( MeshScalarField &scharge, const VectorField &efield, 
					     const VectorField &bfield )
{
    _imp->iterate_trajectories( scharge, efield, bfield );
}


void ParticleDataBase::step_particles( MeshScalarField &scharge, const VectorField &efield, 
				       const VectorField &bfield, double dt )
{
    _imp->step_particles( scharge, efield, bfield, dt );
}



/* ******************************************************************************************* *
 * ParticleDataBase2D                                                                          *
 * ******************************************************************************************* */


ParticleDataBase2D::ParticleDataBase2D( const Geometry &geom )
{
    _imp = new ParticleDataBase2DImp( this, geom );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBase2D::ParticleDataBase2D( const ParticleDataBase2D &pdb )
{
    _imp = new ParticleDataBase2DImp( *pdb._imp );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBase2D::ParticleDataBase2D( std::istream &s, const Geometry &geom )
{
    int32_t fileid = read_int32( s );
    if( fileid != FILEID_PARTICLEDB2D ) {
	if( fileid == FILEID_PARTICLEDB3D || fileid == FILEID_PARTICLEDBCYL )
	    throw( Error( ERROR_LOCATION, "unknown input reading ParticleDataBase2D from stream\n"
			  "incorrect geometry type" ) );
	else
	    throw( Error( ERROR_LOCATION, "unknown input reading ParticleDataBase2D from stream" ) );
    }
    _imp = new ParticleDataBase2DImp( this, s, geom );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBase2D::~ParticleDataBase2D()
{
    delete _imp;
}


const ParticleDataBase2D &ParticleDataBase2D::operator=( const ParticleDataBase2D &pdb )
{
    delete _imp;
    _imp = new ParticleDataBase2DImp( *pdb._imp );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );

    return( *this );
}


Particle2D &ParticleDataBase2D::particle( uint32_t i )
{
    return( _imp->particle( i ) );
}


const Particle2D &ParticleDataBase2D::particle( uint32_t i ) const
{
    return( _imp->particle( i ) ) ;
}


const ParticleP2D &ParticleDataBase2D::trajectory_point( uint32_t i, uint32_t j ) const
{
    return( _imp->trajectory_point( i, j ) );
}


void ParticleDataBase2D::trajectories_at_plane( std::vector<Particle2D> &tdata,
						coordinate_axis_e axis,
						double val ) const
{
    _imp->trajectories_at_plane( tdata, axis, val );
}


void ParticleDataBase2D::add_particle( double IQ, double q, double m, const ParticleP2D &x )
{
    _imp->add_particle( IQ, q, m, x );
}


void ParticleDataBase2D::add_particle( const Particle2D &p )
{
    _imp->add_particle( p );
}


void ParticleDataBase2D::add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
						    double v, double dvp, double dvt, 
						    double x1, double y1, double x2, double y2 )
{
    _imp->add_2d_beam_with_velocity( N, J, q, m, v, dvp, dvt, x1, y1, x2, y2 );
}


void ParticleDataBase2D::add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
						  double E, double Tp, double Tt, 
						  double x1, double y1, double x2, double y2 )
{
    _imp->add_2d_beam_with_energy( N, J, q, m, E, Tp, Tt, x1, y1, x2, y2 );
}


void ParticleDataBase2D::add_2d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
							double a, double b, double e,
							double Ex, double x0, double y0 )
{
    _imp->add_2d_KV_beam_with_emittance( N, I, q, m, a, b, e, Ex, x0, y0 );
}


void ParticleDataBase2D::add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
							      double a, double b, double e,
							      double Ex, double x0, double y0 )
{
    _imp->add_2d_gaussian_beam_with_emittance( N, I, q, m, a, b, e, Ex, x0, y0 );
}


void ParticleDataBase2D::save( const std::string &filename ) const
{
    ibsimu.message( 1 ) << "Saving ParticleDataBase2D to file \'" << filename << "\'.\n";

    std::ofstream os( filename.c_str(), std::ios_base::binary );
    if( !os.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\' for writing" ) );
    write_int32( os, FILEID_PARTICLEDB2D );
    _imp->save( os );
    os.close();
}


void ParticleDataBase2D::save( std::ostream &os ) const
{
    _imp->save( os );
}


void ParticleDataBase2D::debug_print( std::ostream &os ) const
{
    _imp->debug_print( os );
}


/* ******************************************************************************************* *
 * ParticleDataBaseCyl                                                                         *
 * ******************************************************************************************* */


ParticleDataBaseCyl::ParticleDataBaseCyl( const Geometry &geom )
{
    _imp = new ParticleDataBaseCylImp( this, geom );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBaseCyl::ParticleDataBaseCyl( const ParticleDataBaseCyl &pdb )
{
    _imp = new ParticleDataBaseCylImp( *pdb._imp );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBaseCyl::ParticleDataBaseCyl( std::istream &s, const Geometry &geom )
{
    int32_t fileid = read_int32( s );
    if( fileid != FILEID_PARTICLEDBCYL ) {
	if( fileid == FILEID_PARTICLEDB2D || fileid == FILEID_PARTICLEDB3D )
	    throw( Error( ERROR_LOCATION, "unknown input reading ParticleDataBaseCyl from stream\n"
			  "incorrect geometry type" ) );
	else
	    throw( Error( ERROR_LOCATION, "unknown input reading ParticleDataBaseCyl from stream" ) );
    }
    _imp = new ParticleDataBaseCylImp( this, s, geom );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBaseCyl::~ParticleDataBaseCyl()
{
    delete _imp;
}


const ParticleDataBaseCyl &ParticleDataBaseCyl::operator=( const ParticleDataBaseCyl &pdb )
{
    delete _imp;
    _imp = new ParticleDataBaseCylImp( *pdb._imp );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );

    return( *this );
}


ParticleCyl &ParticleDataBaseCyl::particle( uint32_t i )
{
    return( _imp->particle( i ) );
}


const ParticleCyl &ParticleDataBaseCyl::particle( uint32_t i ) const
{
    return( _imp->particle( i ) ) ;
}


const ParticlePCyl &ParticleDataBaseCyl::trajectory_point( uint32_t i, uint32_t j ) const
{
    return( _imp->trajectory_point( i, j ) );
}


void ParticleDataBaseCyl::trajectories_at_plane( std::vector<ParticleCyl> &tdata,
						 coordinate_axis_e axis,
						 double val ) const
{
    _imp->trajectories_at_plane( tdata, axis, val );
}


void ParticleDataBaseCyl::export_path_manager_data( const std::string &filename, 
						    double ref_E, double ref_q, double ref_m, 
						    double val, uint32_t Np ) const
{
    _imp->export_path_manager_data( filename, ref_E, ref_q, ref_m, val, Np );
}


void ParticleDataBaseCyl::add_particle( double IQ, double q, double m, const ParticlePCyl &x )
{
    _imp->add_particle( IQ, q, m, x );
}


void ParticleDataBaseCyl::add_particle( const ParticleCyl &p )
{
    _imp->add_particle( p );
}


void ParticleDataBaseCyl::add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
						     double v, double dvp, double dvt, 
						     double x1, double y1, double x2, double y2 )
{
    _imp->add_2d_beam_with_velocity( N, J, q, m, v, dvp, dvt, x1, y1, x2, y2 );
}


void ParticleDataBaseCyl::add_2d_beam_with_total_energy( uint32_t N, double J, double q, double m, 
							 double Etot, const ScalarField &epot, 
							 double Tp, double Tt, 
							 double x1, double y1, double x2, double y2 )
{
    _imp->add_2d_beam_with_total_energy( N, J, q, m, Etot, epot, Tp, Tt, x1, y1, x2, y2 );
}


void ParticleDataBaseCyl::add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
						   double E, double Tp, double Tt, 
						   double x1, double y1, double x2, double y2 )
{
    _imp->add_2d_beam_with_energy( N, J, q, m, E, Tp, Tt, x1, y1, x2, y2 );
}


void ParticleDataBaseCyl::add_2d_full_gaussian_beam( uint32_t N, double I, double q, double m,
						     double Ex, double Tp, double Tt, 
						     double x0, double dr )
{
    _imp->add_2d_full_gaussian_beam( N, I, q, m, Ex, Tp, Tt, x0, dr );
}


void ParticleDataBaseCyl::add_2d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
							       double a, double b, double e,
							       double Ex, double x0 )
{
    _imp->add_2d_gaussian_beam_with_emittance( N, I, q, m, a, b, e, Ex, x0 );
}


void ParticleDataBaseCyl::save( const std::string &filename ) const
{
    ibsimu.message( 1 ) << "Saving ParticleDataBaseCyl to file \'" << filename << "\'.\n";

    std::ofstream os( filename.c_str(), std::ios_base::binary );
    if( !os.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\' for writing" ) );
    write_int32( os, FILEID_PARTICLEDBCYL );
    _imp->save( os );
    os.close();
}


void ParticleDataBaseCyl::save( std::ostream &os ) const
{
    _imp->save( os );
}


void ParticleDataBaseCyl::debug_print( std::ostream &os ) const
{
    _imp->debug_print( os );
}


/* ******************************************************************************************* *
 * ParticleDataBase3D                                                                          *
 * ******************************************************************************************* */


ParticleDataBase3D::ParticleDataBase3D( const Geometry &geom )
{
    _imp = new ParticleDataBase3DImp( this, geom );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBase3D::ParticleDataBase3D( const ParticleDataBase3D &pdb )
{
    _imp = new ParticleDataBase3DImp( *pdb._imp );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBase3D::ParticleDataBase3D( std::istream &s, const Geometry &geom )
{
    int32_t fileid = read_int32( s );
    if( fileid != FILEID_PARTICLEDB3D ) {
	if( fileid == FILEID_PARTICLEDB2D || fileid == FILEID_PARTICLEDBCYL )
	    throw( Error( ERROR_LOCATION, "unknown input reading ParticleDataBase3D from stream\n"
			  "incorrect geometry type" ) );
	else
	    throw( Error( ERROR_LOCATION, "unknown input reading ParticleDataBase3D from stream" ) );
    }
    _imp = new ParticleDataBase3DImp( this, s, geom );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );
}


ParticleDataBase3D::~ParticleDataBase3D()
{
    delete _imp;
}


const ParticleDataBase3D &ParticleDataBase3D::operator=( const ParticleDataBase3D &pdb )
{
    delete _imp;
    _imp = new ParticleDataBase3DImp( *pdb._imp );
    set_implementation_pointer( (ParticleDataBaseImp *)_imp );

    return( *this );
}


Particle3D &ParticleDataBase3D::particle( uint32_t i )
{
    return( _imp->particle( i ) );
}


const Particle3D &ParticleDataBase3D::particle( uint32_t i ) const
{
    return( _imp->particle( i ) ) ;
}


const ParticleP3D &ParticleDataBase3D::trajectory_point( uint32_t i, uint32_t j ) const
{
    return( _imp->trajectory_point( i, j ) );
}


void ParticleDataBase3D::trajectories_at_plane( std::vector<Particle3D> &tdata,
						coordinate_axis_e axis,
						double val ) const
{
    _imp->trajectories_at_plane( tdata, axis, val );
}


void ParticleDataBase3D::add_particle( double IQ, double q, double m, const ParticleP3D &x )
{
    _imp->add_particle( IQ, q, m, x );
}


void ParticleDataBase3D::add_particle( const Particle3D &p )
{
    _imp->add_particle( p );
}


void ParticleDataBase3D::add_cylindrical_beam_with_velocity( uint32_t N, double J, double q, double m, 
							     double v, double dvp, double dvt, Vec3D c, 
							     Vec3D dir1, Vec3D dir2, double r )
{
    _imp->add_cylindrical_beam_with_velocity( N, J, q, m, v, dvp, dvt, c, dir1, dir2, r );
}


void ParticleDataBase3D::add_cylindrical_beam_with_total_energy( uint32_t N, double J, double q, double m, 
								 double Etot, const ScalarField &epot, 
								 double Tp, double Tt, Vec3D c, 
								 Vec3D dir1, Vec3D dir2, double r )
{
    _imp->add_cylindrical_beam_with_total_energy( N, J, q, m, Etot, epot, Tp, Tt, c, dir1, dir2, r );
}

void ParticleDataBase3D::add_cylindrical_beam_with_energy( uint32_t N, double J, double q, double m, 
							   double E, double Tp, double Tt, Vec3D c,
							   Vec3D dir1, Vec3D dir2, double r )
{
    _imp->add_cylindrical_beam_with_energy( N, J, q, m, E, Tp, Tt, c, dir1, dir2, r );
}


void ParticleDataBase3D::add_rectangular_beam_with_velocity( uint32_t N, double J, double q, double m, 
							     double v, double dvp, double dvt, Vec3D c, 
							     Vec3D dir1, Vec3D dir2, double size1, double size2 )
{
    _imp->add_rectangular_beam_with_velocity( N, J, q, m, v, dvp, dvt, c, dir1, dir2, size1, size2 );
}


void ParticleDataBase3D::add_rectangular_beam_with_energy( uint32_t N, double J, double q, double m, 
							   double E, double Tp, double Tt, Vec3D c, 
							   Vec3D dir1, Vec3D dir2, double size1, double size2 )
{
    _imp->add_rectangular_beam_with_energy( N, J, q, m, E, Tp, Tt, c, dir1, dir2, size1, size2 );
}


void ParticleDataBase3D::add_3d_KV_beam_with_emittance( uint32_t N, double I, double q, double m,
							double E0,
							double a1, double b1, double e1,
							double a2, double b2, double e2,
							Vec3D c, Vec3D dir1, Vec3D dir2 )
{
    _imp->add_3d_KV_beam_with_emittance( N, I, q, m, E0, a1, b1, e1, a2, b2, e2, c, dir1, dir2 );
}


void ParticleDataBase3D::add_3d_waterbag_beam_with_emittance( uint32_t N, double I, double q, double m,
							      double E0,
							      double a1, double b1, double e1,
							      double a2, double b2, double e2,
							      Vec3D c, Vec3D dir1, Vec3D dir2 )
{
    _imp->add_3d_waterbag_beam_with_emittance( N, I, q, m, E0, a1, b1, e1, a2, b2, e2, c, dir1, dir2 );
}


void ParticleDataBase3D::add_3d_gaussian_beam_with_emittance( uint32_t N, double I, double q, double m,
							      double E0, 
							      double a1, double b1, double e1,
							      double a2, double b2, double e2,
							      Vec3D c, Vec3D dir1, Vec3D dir2 )
{
    _imp->add_3d_gaussian_beam_with_emittance( N, I, q, m, E0, a1, b1, e1, a2, b2, e2, c, dir1, dir2 );
}


void ParticleDataBase3D::trajectories_at_free_plane( TrajectoryDiagnosticData &tdata, 
						     const Vec3D &c, const Vec3D &o, const Vec3D &p,
						     const std::vector<trajectory_diagnostic_e> &diagnostics ) const
{
    _imp->trajectories_at_free_plane( tdata, c, o, p, diagnostics );
}


void ParticleDataBase3D::export_path_manager_data( const std::string &filename, 
						   double ref_E, double ref_q, double ref_m, 
						   const Vec3D &c, const Vec3D &o, const Vec3D &p ) const
{
    _imp->export_path_manager_data( filename, ref_E, ref_q, ref_m, c, o, p );
}


void ParticleDataBase3D::save( const std::string &filename ) const
{
    ibsimu.message( 1 ) << "Saving ParticleDataBase3D to file \'" << filename << "\'.\n";

    std::ofstream os( filename.c_str(), std::ios_base::binary );
    if( !os.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\' for writing" ) );
    write_int32( os, FILEID_PARTICLEDB3D );
    _imp->save( os );
    os.close();
}


void ParticleDataBase3D::save( std::ostream &os ) const
{
    _imp->save( os );
}


void ParticleDataBase3D::debug_print( std::ostream &os ) const
{
    _imp->debug_print( os );
}


