/*! \file gtkfielddiagwindow.cpp
 *  \brief Field diagnostic window
 */

/* Copyright (c) 2005-2012,2016 Taneli Kalvas. All rights reserved.
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
#include "gtkfielddiagwindow.hpp"
#include "gtkfielddiagexportdialog.hpp"
#include "fonts.hpp"


GTKFieldDiagWindow::GTKFieldDiagWindow( GTKPlotter &plotter, const Geometry &geom, size_t N, 
					const Vec3D &x1, const Vec3D &x2,
					const field_diag_type_e diag[2], 
					const field_loc_type_e loc[2] )
    : GTKFrameWindow(plotter), _geom(geom), _plot(_frame,geom)
{
    _plot.set_epot( plotter.get_epot() );
    _plot.set_efield( plotter.get_efield() );
    _plot.set_scharge( plotter.get_scharge() );
    _plot.set_bfield( plotter.get_bfield() );

    _plot.set_coordinates( N, x1, x2 );
    _plot.set_diagnostic( diag, loc );

    // Set window title
    gtk_window_set_title( GTK_WINDOW(_window), "Field diagnostics" );

    // Add export menu item
    GtkWidget *item_export = gtk_menu_item_new_with_mnemonic( "_Export" );
    gtk_menu_shell_prepend( GTK_MENU_SHELL(_menu_file), item_export );
    g_signal_connect( G_OBJECT(item_export), "activate",
		      G_CALLBACK(menuitem_export_signal),
		      (gpointer)this );

    try {
	_plot.build_plot();
    } catch( Error e ) {
	// Destroy window
	gtk_widget_destroy( _window );
	//plotter->delete_window( this ); Didn't work, still got window and pressing 'X' caused seg fault
	throw e;
    }
    _frame.get_ranges( PLOT_AXIS_X1, _x1min, _x1max );
    _frame.get_ranges( PLOT_AXIS_X1, _x2min, _x2max );
    show();
}


GTKFieldDiagWindow::~GTKFieldDiagWindow()
{
}


void GTKFieldDiagWindow::export_data( void )
{
    GTKFieldDiagExportDialog dialog( _window, &_plot );
    dialog.run();
}


void GTKFieldDiagWindow::menuitem_export_signal( GtkToolButton *button,
						 gpointer object )
{
    GTKFieldDiagWindow *window = (GTKFieldDiagWindow *)object;
    window->export_data();
}


void GTKFieldDiagWindow::zoom_fit( void )
{
    double min = -std::numeric_limits<double>::infinity();
    double max = std::numeric_limits<double>::infinity();
    _frame.set_ranges( PLOT_AXIS_X1, _x1min, _x1max );
    _frame.set_ranges( PLOT_AXIS_Y1, min, max );
    _frame.set_ranges( PLOT_AXIS_X2, _x2min, _x2max );
    _frame.set_ranges( PLOT_AXIS_Y2, min, max );
    _frame.ruler_autorange_enable( PLOT_AXIS_X1, false, false );
    _frame.ruler_autorange_enable( PLOT_AXIS_Y1, false, false );
    _frame.ruler_autorange_enable( PLOT_AXIS_X2, false, false );
    _frame.ruler_autorange_enable( PLOT_AXIS_Y2, false, false );

    draw_and_expose();
}


void *GTKFieldDiagWindow::build_preferences( GtkWidget *notebook )
{
    return( NULL );
}


void GTKFieldDiagWindow::read_preferences( GtkWidget *notebook, void *pdata )
{

}


std::string GTKFieldDiagWindow::track_text( double x, double y )
{
    std::wstringstream ss;
    std::locale mylocale("");
    ss.imbue( mylocale );

    // Calculate relative track location on x-axis
    double xmin, xmax;
    _frame.get_ranges( PLOT_AXIS_X1, xmin, xmax );
    double t = (x - xmin)/(xmax-xmin);

    // Calculate absolute coordinates of track location
    Vec3D start = _plot.start();
    Vec3D end = _plot.end();
    Vec3D xc = start + t*(end-start);
    double dist = norm2( xc-start );

    if( _geom.geom_mode() == MODE_3D ) {
	ss << "x = " << xc[0] << " m\n"
	   << "y = " << xc[1] << " m\n"
	   << "z = " << xc[2] << " m\n"
	   << "d = " << dist << " m\n";
    } else if( _geom.geom_mode() == MODE_2D ) {
	ss << "x = " << xc[0] << " m\n"
	   << "y = " << xc[1] << " m\n"
	   << "d = " << dist << " m\n";
    } else if( _geom.geom_mode() == MODE_CYL ) {
	ss << "x = " << xc[0] << " m\n"
	   << "r = " << xc[1] << " m\n"
	   << "d = " << dist << " m\n";
    }

    // Process Y-axes
    double x2 = _track_px;
    double y2 = _track_py;
    Coordmapper cm = _frame.get_coordmapper( PLOT_AXIS_X2, PLOT_AXIS_Y2 );
    cm.inv_transform( x2, y2 );
    for( int i = 0; i < 2; i++ ) {

	double val;
	if( i == 0 ) {
	    val = y;
	} else {
	    val = y2;
	}

	switch( _plot.get_diagnostic_type( i ) ) {
	case FIELD_EPOT:
	    ss << "epot = " << val << "V\n";
	    break;
	    
	case FIELD_EFIELD:
	    ss << "|E| = " << val << "V/m\n";
	    break;
	    
	case FIELD_EFIELD_X:
	    ss << "Ex = " << val << "V/m\n";
	    break;
	    
	case FIELD_EFIELD_Y:
	    ss << "Ey = " << val << "V/m\n";
	    break;
	    
	case FIELD_EFIELD_Z:
	    ss << "Ez = " << val << "V/m\n";
	    break;
	    
	case FIELD_SCHARGE:
	    ss << "rho = " << val << "C/m^3\n";
	    break;
	    
	case FIELD_BFIELD:
	    ss << "|B| = " << val << "T\n";
	    break;
	    
	case FIELD_BFIELD_X:
	    ss << "Bx = " << val << "T\n";
	    break;
	    
	case FIELD_BFIELD_Y:
	    ss << "By = " << val << "T\n";
	    break;
	    
	case FIELD_BFIELD_Z:
	    ss << "Bz = " << val << "T\n";
	    break;
	    
	case FIELD_NONE:
	default:
	    break;
	};
    }
    
    return( wstring_to_utf8( ss.str() ) );
}

