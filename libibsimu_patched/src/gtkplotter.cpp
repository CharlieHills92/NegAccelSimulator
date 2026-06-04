/*! \file gtkplotter.cpp
 *  \brief GTK based plotters.
 */

/* Copyright (c) 2005-2013 Taneli Kalvas. All rights reserved.
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


#include "config.h"

#ifdef OPENGL
#include <gtk/gtkgl.h>
#endif

#include "gtkplotter.hpp"
#include "gtkgeomwindow.hpp"
#include "gtkgeom3dwindow.hpp"
#include "gtkparticlediagwindow.hpp"
#include "gtkfielddiagwindow.hpp"
#include "ibsimu.hpp"
#include "error.hpp"



bool GTKPlotter::_gtk_initialized = false;
bool GTKPlotter::_opengl = false;


GTKPlotter::GTKPlotter( int *argc, char ***argv )
    : _sdata(NULL), _geom(NULL), _epot(NULL), _efield(NULL), _scharge(NULL), _tdens(NULL), 
      _bfield(NULL), _pdb(NULL)
{
    if( !_gtk_initialized ) {

	// Initialize gtk
	if( gtk_init_check( argc, argv ) == FALSE )
	    throw( Error( ERROR_LOCATION, "Couldn't initialize GTK" ) );
	_gtk_initialized = true;

#ifdef OPENGL
	_opengl = true;
	// Initialize OpenGL
	if( gtk_gl_init_check( argc, argv ) == FALSE ||
	    gdk_gl_query_extension() == FALSE )
	    _opengl = false;
#endif
    }
}


GTKPlotter::~GTKPlotter()
{
    
}


void GTKPlotter::force_software_renderer( void ) const
{
    _opengl = false;
}


bool GTKPlotter::opengl( void ) const
{
    return( _opengl );
}


void GTKPlotter::run()
{
    ibsimu.message( 1 ) << "Running GTKPlotter\n";
    ibsimu.inc_indent();

    gtk_main();

    ibsimu.message( 1 ) << "Done\n";
    ibsimu.dec_indent();
}


GTKWindow *GTKPlotter::new_geometry_plot_window( void )
{
    GTKWindow *window = new GTKGeomWindow( *this, *_geom, _epot, _efield, _scharge, _tdens, _bfield, _pdb );
    _windows.push_back( window );

    return( window );
}


GTKWindow *GTKPlotter::new_geometry_3d_plot_window( void )
{
    GTKWindow *window = new GTKGeom3DWindow( *this, *_geom, _pdb, _sdata );
    _windows.push_back( window );

    return( window );
}


GTKWindow *GTKPlotter::new_particle_plot_window( coordinate_axis_e axis, double level, 
						 particle_diag_plot_type_e type,
						 trajectory_diagnostic_e diagx, 
						 trajectory_diagnostic_e diagy )
{
    if( !_pdb )
	throw( Error( ERROR_LOCATION, "Particle database not defined" ) );
    if( !_geom )
	throw( Error( ERROR_LOCATION, "Geometry not defined" ) );
    GTKWindow *window = new GTKParticleDiagWindow( *this, *_pdb, *_geom, axis, level, type, diagx, diagy );
    _windows.push_back( window );

    return( window );
}


GTKWindow *GTKPlotter::new_field_plot_window( size_t N, const Vec3D &x1, const Vec3D &x2,
					      const field_diag_type_e diag[2], 
					      const field_loc_type_e loc[2] )
{
    if( !_geom )
	throw( Error( ERROR_LOCATION, "Geometry not defined" ) );
    GTKWindow *window = new GTKFieldDiagWindow( *this, *_geom, N, x1, x2, diag, loc );
    _windows.push_back( window );

    return( window );
}



void GTKPlotter::delete_window( GTKWindow *window )
{
    for( size_t a = 0; a < _windows.size(); a++ ) {
	if( _windows[a] == window ) {
	    delete( window );
	    _windows.erase( _windows.begin()+a );
	    break;
	}
    }

    if( _windows.size() == 0 ) {
	gtk_main_quit();
    }
}


const std::vector<double> *GTKPlotter::get_surface_triangle_data( void ) const
{
    return( _sdata );
}


const Geometry *GTKPlotter::get_geometry( void ) const
{
    return( _geom );
}


const EpotField *GTKPlotter::get_epot( void ) const
{
    return( _epot );
}


const EpotEfield *GTKPlotter::get_efield( void ) const
{
    return( _efield );
}


const MeshScalarField *GTKPlotter::get_scharge( void ) const
{
    return( _scharge );
}


const MeshScalarField *GTKPlotter::get_trajdens( void ) const
{
    return( _tdens );
}


const VectorField *GTKPlotter::get_bfield( void ) const
{
    return( _bfield );
}


const ParticleDataBase *GTKPlotter::get_particledatabase( void ) const
{
    return( _pdb );
}


void GTKPlotter::set_surface_triangle_data( const std::vector<double> *data )
{
    _sdata = data;
}


void GTKPlotter::set_geometry( const Geometry *geom )
{
    _geom = geom;
}


void GTKPlotter::set_epot( const EpotField *epot )
{
    _epot = epot;
}


void GTKPlotter::set_efield( const EpotEfield *efield )
{
    _efield = efield;
}


void GTKPlotter::set_scharge( const MeshScalarField *scharge )
{
    _scharge = scharge;
}


void GTKPlotter::set_trajdens( const MeshScalarField *tdens )
{
    _tdens = tdens;
}


void GTKPlotter::set_bfield( const VectorField *bfield )
{
    _bfield = bfield;
}


void GTKPlotter::set_particledatabase( const ParticleDataBase *pdb )
{
    _pdb = pdb;
}

