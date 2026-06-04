/*! \file multimeshvectorfield.cpp
 *  \brief %Vector field using multiple meshes.
 */

/* Copyright (c) 2011,2014 Taneli Kalvas. All rights reserved.
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


#include <limits>
#include <fstream>
#include "multimeshvectorfield.hpp"
#include "ibsimu.hpp"
#include "compmath.hpp"
#include "error.hpp"


MultiMeshVectorField::MultiMeshVectorField()
{
    _field.push_back( new MeshVectorField );
}


MultiMeshVectorField::MultiMeshVectorField( const Mesh &m, const bool fout[3] )
{
    _field.push_back( new MeshVectorField( m, fout ) );
}


MultiMeshVectorField::MultiMeshVectorField( geom_mode_e geom_mode, const bool fout[3], Int3D size, 
					    Vec3D origo, double h )
{
    _field.push_back( new MeshVectorField( geom_mode, fout, size, origo, h ) );
}


MultiMeshVectorField::MultiMeshVectorField( geom_mode_e geom_mode, const bool fout[3], double xscale, 
					    double fscale, const std::string &filename )
{
    _field.push_back( new MeshVectorField( geom_mode, fout, xscale, fscale, filename ) );
}


MultiMeshVectorField::MultiMeshVectorField( const MultiMeshVectorField &f )
{
    for( size_t i = 0; i < f._field.size(); i++ )
	_field.push_back( new MeshVectorField( *f._field[i] ) );
}


MultiMeshVectorField::MultiMeshVectorField( std::istream &s )
{
    uint32_t N = (uint32_t)read_int32( s );
    for( size_t i = 0; i < N; i++ )
	_field.push_back( new MeshVectorField( s ) );
}


MultiMeshVectorField::~MultiMeshVectorField()
{
    for( size_t i = 0; i < _field.size(); i++ )
	delete _field[i];
}


void MultiMeshVectorField::set_extrapolation( const field_extrpl_e extrpl[6] )
{
    _field[0]->set_extrapolation( extrpl );
}


void MultiMeshVectorField::translate( Vec3D x )
{
    for( size_t i = 0; i < _field.size(); i++ )
	_field[i]->translate( x );
}


void MultiMeshVectorField::scale( double s )
{
    for( size_t i = 0; i < _field.size(); i++ )
	_field[i]->scale( s );
}


void MultiMeshVectorField::rotate_x( double a )
{
    for( size_t i = 0; i < _field.size(); i++ )
	_field[i]->rotate_x( a );
}


void MultiMeshVectorField::rotate_y( double a )
{
    for( size_t i = 0; i < _field.size(); i++ )
	_field[i]->rotate_y( a );
}


void MultiMeshVectorField::rotate_z( double a )
{
    for( size_t i = 0; i < _field.size(); i++ )
	_field[i]->rotate_z( a );
}


void MultiMeshVectorField::clear()
{
    for( size_t i = 0; i < _field.size(); i++ )
	_field[i]->clear();
}


void MultiMeshVectorField::reset( geom_mode_e geom_mode, const bool fout[3], Int3D size, 
				  Vec3D origo, double h )
{
    for( size_t i = 1; i < _field.size(); i++ )
	delete _field[i];
    _field[0]->reset( geom_mode, fout, size, origo, h );
}


void MultiMeshVectorField::add_mesh( const MeshVectorField &field )
{
    // Check geometry type
    if( field.geom_mode() != _field[0]->geom_mode() )
	throw( Error( ERROR_LOCATION, "incompatible geometries in fields" ) );

    // Check defined field components
    bool fout0[3];
    _field[0]->get_defined_components( fout0 );
    bool fout[3];
    field.get_defined_components( fout );
    for( uint32_t a = 0; a < 3; a++ ) {
	if( fout[a] != fout0[a] )
	    throw( Error( ERROR_LOCATION, "incompatible field components in fields" ) );
    }

    _field.push_back( new MeshVectorField( field ) );
    
    // Set NaN extrapolation setting
    field_extrpl_e extrpl[6] = { FIELD_NAN, FIELD_NAN, FIELD_NAN, 
				 FIELD_NAN, FIELD_NAN, FIELD_NAN };
    _field.back()->set_extrapolation( extrpl );
}
    

void MultiMeshVectorField::add_mesh( Int3D size, Vec3D origo, double h )
{
    bool fout[3];
    _field[0]->get_defined_components( fout );
    _field.push_back( new MeshVectorField( _field[0]->geom_mode(), fout, size, origo, h ) );
    
    // Set NaN extrapolation setting
    field_extrpl_e extrpl[6] = { FIELD_NAN, FIELD_NAN, FIELD_NAN, 
				 FIELD_NAN, FIELD_NAN, FIELD_NAN };
    _field.back()->set_extrapolation( extrpl );
}
    

void MultiMeshVectorField::add_mesh( double xscale, double fscale, const std::string &filename )
{
    bool fout[3];
    _field[0]->get_defined_components( fout );
    _field.push_back( new MeshVectorField( _field[0]->geom_mode(), fout, xscale, fscale, filename ) );

    // Set NaN extrapolation setting
    field_extrpl_e extrpl[6] = { FIELD_NAN, FIELD_NAN, FIELD_NAN, 
				 FIELD_NAN, FIELD_NAN, FIELD_NAN };
    _field.back()->set_extrapolation( extrpl );
}


void MultiMeshVectorField::get_minmax( double &min, double &max ) const
{
    min = std::numeric_limits<double>::infinity();
    max = -std::numeric_limits<double>::infinity();
    for( size_t i = 0; i < _field.size(); i++ ) {

	double mmin, mmax;
	_field[i]->get_minmax( mmin, mmax );

	if( mmin < min )
	    min = mmin;
	if( mmax > max )
	    max = mmax;
    }
}


void MultiMeshVectorField::get_defined_components( bool fout[3] ) const
{
    _field[0]->get_defined_components( fout );
}


MultiMeshVectorField &MultiMeshVectorField::operator=( const MultiMeshVectorField &f )
{
    // Delete old
    for( size_t i = 0; i < _field.size(); i++ )
	delete _field[i];
    _field.clear();

    // Copy new
    for( size_t i = 0; i < f._field.size(); i++ )
	_field.push_back( new MeshVectorField( *f._field[i] ) );

    return( *this );
}


const Vec3D MultiMeshVectorField::operator()( const Vec3D &x ) const
{
    // Scan fields for field evaluation
    for( int32_t i = _field.size()-1; i > 0; i-- ) {
	Vec3D R = (*_field[i])( x );
	if( !comp_isnan( R[0] ) )
	    return( R );
    }

    // Use coarsest field 
    Vec3D R = (*_field[0])( x );
    return( R );
}


void MultiMeshVectorField::save( const std::string &filename ) const
{
    ibsimu.message( 1 ) << "Saving MultiMeshVectorField to file \'" << filename << "\'.\n";

    std::ofstream os( filename.c_str() );
    if( !os.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\' for writing" ) );
    save( os );
    os.close();
}


void MultiMeshVectorField::save( std::ostream &s ) const
{
    write_int32( s, _field.size() );
    for( size_t i = 0; i < _field.size(); i++ )
	_field[i]->save( s );
}


void MultiMeshVectorField::debug_print( std::ostream &os ) const
{
    os << "**MultiMeshVectorField\n";

    os << "n = " << _field.size() << "\n";
    for( size_t i = 0; i < _field.size(); i++ )
	_field[i]->debug_print( os );
}


