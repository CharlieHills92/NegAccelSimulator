/*! \file gtkparticlediagwindow.cpp
 *  \brief %Particle diagnostic window.
 */

/* Copyright (c) 2005-2012,2014,2016,2022 Taneli Kalvas. All rights reserved.
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

#include <sstream>
#include <limits>
#include <locale>
#include "gtkparticlediagwindow.hpp"
#include "gtkparticlediagexportdialog.hpp"
#include "histogram.hpp"
#include "meshcolormap.hpp"
#include "fonts.hpp"


GTKParticleDiagWindow::GTKParticleDiagWindow( GTKPlotter &plotter, const ParticleDataBase &pdb, 
					      const Geometry &geom,
					      coordinate_axis_e axis, double level, 
					      particle_diag_plot_type_e type,
					      trajectory_diagnostic_e diagx, 
					      trajectory_diagnostic_e diagy )
    : GTKFrameWindow(plotter), _plot(_frame, geom, pdb, axis, level, type, diagx, diagy), _prefdata(NULL)
{
    // Set window title
    gtk_window_set_title( GTK_WINDOW(_window), "Particle diagnostics" );

    // Add export menu item
    GtkWidget *item_export = gtk_menu_item_new_with_mnemonic( "_Export" );
    gtk_menu_shell_prepend( GTK_MENU_SHELL(_menu_file), item_export );
    g_signal_connect( G_OBJECT(item_export), "activate",
                      G_CALLBACK(menuitem_export_signal),
                      (gpointer)this );

    _plot.build_plot();
    show();
}


GTKParticleDiagWindow::~GTKParticleDiagWindow()
{
    if( _prefdata )
	delete _prefdata;
}



void GTKParticleDiagWindow::plot_2d_type_toggled_process( void )
{
    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->plot_scatter_radio) ) ) {
	gtk_widget_set_sensitive( _prefdata->histo_n_spin, false );
	gtk_widget_set_sensitive( _prefdata->histo_m_spin, false );
	gtk_widget_set_sensitive( _prefdata->histo_acc_closest_radio, false );
	gtk_widget_set_sensitive( _prefdata->histo_acc_bilinear_radio, false );
	gtk_widget_set_sensitive( _prefdata->int_closest_radio, false );
	gtk_widget_set_sensitive( _prefdata->int_bilinear_radio, false );
	gtk_widget_set_sensitive( _prefdata->int_bicubic_radio, false );
	gtk_widget_set_sensitive( _prefdata->dot_size_spin, true );
    } else {
	gtk_widget_set_sensitive( _prefdata->histo_n_spin, true );
	gtk_widget_set_sensitive( _prefdata->histo_m_spin, true );
	gtk_widget_set_sensitive( _prefdata->histo_acc_closest_radio, true );
	gtk_widget_set_sensitive( _prefdata->histo_acc_bilinear_radio, true );
	gtk_widget_set_sensitive( _prefdata->int_closest_radio, true );
	gtk_widget_set_sensitive( _prefdata->int_bilinear_radio, true );
	gtk_widget_set_sensitive( _prefdata->int_bicubic_radio, true );
	gtk_widget_set_sensitive( _prefdata->dot_size_spin, false );
    }
}


void GTKParticleDiagWindow::plot_2d_type_toggled( GtkToggleButton *togglebutton,
					  gpointer         user_data )
{
   GTKParticleDiagWindow *pdwindow = (GTKParticleDiagWindow *)user_data;
    if( gtk_toggle_button_get_active( togglebutton ) )
	pdwindow->plot_2d_type_toggled_process();
}


void GTKParticleDiagWindow::build_preferences_2d( GtkWidget *notebook )
{
    // ****************************************************************************
    // Preferences notebook

    GtkWidget *grid = gtk_grid_new();
    gtk_grid_set_column_homogeneous( GTK_GRID(grid), TRUE );
    gtk_widget_set_hexpand( grid, TRUE );
    uint32_t yl = 0;

    // ****************************************************************************

    // Plot type
    GtkWidget *label = gtk_label_new( "Plot type" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    _prefdata->plot_scatter_radio = gtk_radio_button_new_with_label_from_widget( NULL,
										 "Scatter" );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->plot_scatter_radio, 1, yl, 1, 1 );
    yl++;
    _prefdata->plot_colormap_radio = gtk_radio_button_new_with_label_from_widget( 
	GTK_RADIO_BUTTON(_prefdata->plot_scatter_radio), "Colormap" );
    particle_diag_plot_type_e type = _plot.get_type();
    if( type == PARTICLE_DIAG_PLOT_SCATTER )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->plot_scatter_radio), true );
    if( type == PARTICLE_DIAG_PLOT_HISTO2D )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->plot_colormap_radio), true );
    g_signal_connect( _prefdata->plot_scatter_radio, "toggled",
		      G_CALLBACK(plot_2d_type_toggled), (gpointer)this );
    g_signal_connect( _prefdata->plot_colormap_radio, "toggled",
		      G_CALLBACK(plot_2d_type_toggled), (gpointer)this );
    if( type == PARTICLE_DIAG_PLOT_HISTO1D ) {
	gtk_widget_set_sensitive( _prefdata->plot_scatter_radio, false );
	gtk_widget_set_sensitive( _prefdata->plot_colormap_radio, false );
    }
    gtk_grid_attach( GTK_GRID(grid), _prefdata->plot_colormap_radio, 1, yl, 1, 1 );
    yl++;

    // Histogram size
    label = gtk_label_new( "Bin count x" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkAdjustment *histo_n_adj = gtk_adjustment_new( _plot.get_histogram_n(), 0, 1000, 1, 10, 0 );
    _prefdata->histo_n_spin = gtk_spin_button_new( histo_n_adj, 1, 0 );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->histo_n_spin, 1, yl, 1, 1 );
    yl++;

    // Histogram size
    label = gtk_label_new( "Bin count y" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkAdjustment *histo_m_adj = gtk_adjustment_new( _plot.get_histogram_m(), 0, 1000, 1, 10, 0 );
    _prefdata->histo_m_spin = gtk_spin_button_new( histo_m_adj, 1, 0 );
    if( type == PARTICLE_DIAG_PLOT_HISTO1D )
	gtk_widget_set_sensitive( _prefdata->histo_m_spin, false );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->histo_m_spin, 1, yl, 1, 1 );
    yl++;

    // Histogram accumulation type
    histogram_accumulation_e accu = _plot.get_histogram_accumulation();
    label = gtk_label_new( "Histogram accumulation" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    _prefdata->histo_acc_closest_radio = gtk_radio_button_new_with_label_from_widget( NULL,
										      "Closest" );
    if( accu == HISTOGRAM_ACCUMULATION_CLOSEST )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->histo_acc_closest_radio), true );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->histo_acc_closest_radio, 1, yl, 1, 1 );
    yl++;
    _prefdata->histo_acc_bilinear_radio = gtk_radio_button_new_with_label_from_widget( 
	GTK_RADIO_BUTTON(_prefdata->histo_acc_closest_radio), "Bilinear" );
    if( accu == HISTOGRAM_ACCUMULATION_LINEAR )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->histo_acc_bilinear_radio), true );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->histo_acc_bilinear_radio, 1, yl, 1, 1 );
    yl++;

    // Colormap interpolation style
    interpolation_e interpolation =  _plot.get_colormap_interpolation();
    label = gtk_label_new( "Colormap interpolation" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    _prefdata->int_closest_radio = gtk_radio_button_new_with_label_from_widget( NULL,
										"Closest" );
    if( type == PARTICLE_DIAG_PLOT_HISTO1D )
	gtk_widget_set_sensitive( _prefdata->int_closest_radio, false );
    if( interpolation == INTERPOLATION_CLOSEST )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->int_closest_radio), true );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->int_closest_radio, 1, yl, 1, 1 );
    yl++;
    _prefdata->int_bilinear_radio = gtk_radio_button_new_with_label_from_widget( 
	GTK_RADIO_BUTTON(_prefdata->int_closest_radio), "Bilinear" );
    if( type == PARTICLE_DIAG_PLOT_HISTO1D )
	gtk_widget_set_sensitive( _prefdata->int_bilinear_radio, false );
    if( interpolation == INTERPOLATION_BILINEAR )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->int_bilinear_radio), true );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->int_bilinear_radio, 1, yl, 1, 1 );
    yl++;
    _prefdata->int_bicubic_radio = gtk_radio_button_new_with_label_from_widget( 
	GTK_RADIO_BUTTON(_prefdata->int_closest_radio), "Bicubic" );
    if( type == PARTICLE_DIAG_PLOT_HISTO1D )
	gtk_widget_set_sensitive( _prefdata->int_bicubic_radio, false );
    if( interpolation == INTERPOLATION_BICUBIC )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->int_bicubic_radio), true );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->int_bicubic_radio, 1, yl, 1, 1 );
    yl++;

    // Dot size
    label = gtk_label_new( "Dot size" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkAdjustment *dot_size_adj = gtk_adjustment_new( _plot.get_dot_size(), 0.1, 10.0, 0.1, 1, 0 );
    _prefdata->dot_size_spin = gtk_spin_button_new( dot_size_adj, 0.1, 1 );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->dot_size_spin, 1, yl, 1, 1 );
    yl++;

    // Ellipse fit
    label = gtk_label_new( "Emittance ellipse fit" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    _prefdata->ellipse_check = gtk_check_button_new_with_label( "on/off" );
    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->ellipse_check), _plot.get_emittance_ellipse() );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->ellipse_check, 1, yl, 1, 1 );
    yl++;

    // Add notebook page
    label = gtk_label_new( "Particles" );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), grid, label );

    // Set active switches
    plot_2d_type_toggled_process();
}


void GTKParticleDiagWindow::build_preferences_1d( GtkWidget *notebook )
{
    // ****************************************************************************
    // Preferences notebook

    GtkWidget *grid = gtk_grid_new();
    gtk_grid_set_column_homogeneous( GTK_GRID(grid), TRUE );
    gtk_widget_set_hexpand( grid, TRUE );
    uint32_t yl = 0;

    // ****************************************************************************

    GtkWidget *label = gtk_label_new( "Bin count" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkAdjustment *histo_n_adj = gtk_adjustment_new( _plot.get_histogram_n(), 0, 1000, 1, 10, 0 );
    _prefdata->histo_n_spin = gtk_spin_button_new( histo_n_adj, 1, 0 );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->histo_n_spin, 1, yl, 1, 1 );
    yl++;

    // Histogram accumulation type
    histogram_accumulation_e accu = _plot.get_histogram_accumulation();
    label = gtk_label_new( "Histogram accumulation" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    _prefdata->histo_acc_closest_radio = gtk_radio_button_new_with_label_from_widget( 
	NULL, "Closest" );
    if( accu == HISTOGRAM_ACCUMULATION_CLOSEST )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->histo_acc_closest_radio), true );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->histo_acc_closest_radio, 1, yl, 1, 1 );
    yl++;
    _prefdata->histo_acc_bilinear_radio = gtk_radio_button_new_with_label_from_widget( 
	GTK_RADIO_BUTTON(_prefdata->histo_acc_closest_radio), "Linear" );
    if( accu == HISTOGRAM_ACCUMULATION_LINEAR )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->histo_acc_bilinear_radio), true );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->histo_acc_bilinear_radio, 1, yl, 1, 1 );
    yl++;

    // Plot presentation style
    bool histo_style =  _plot.get_histogram_style();
    label = gtk_label_new( "Plot style" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    _prefdata->style_histo_radio = gtk_radio_button_new_with_label_from_widget( 
	NULL, "Histogram" );
    if( histo_style )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->style_histo_radio), true );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->style_histo_radio, 1, yl, 1, 1 );
    yl++;
    _prefdata->style_line_radio = gtk_radio_button_new_with_label_from_widget( 
	GTK_RADIO_BUTTON(_prefdata->style_histo_radio), "Lines" );
    if( !histo_style )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->style_line_radio), true );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->style_line_radio, 1, yl, 1, 1 );
    yl++;

    // Add notebook page
    label = gtk_label_new( "Particles" );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), grid, label );
}


void *GTKParticleDiagWindow::build_preferences( GtkWidget *notebook )
{
    if( _prefdata )
	delete _prefdata;
    _prefdata = new PreferencesData;

    if( _plot.get_type() == PARTICLE_DIAG_PLOT_HISTO1D )
	build_preferences_1d( notebook );
    else
	build_preferences_2d( notebook );

    return( NULL );
}


void GTKParticleDiagWindow::read_preferences( GtkWidget *notebook, void *pdata )
{
    if( _plot.get_type() == PARTICLE_DIAG_PLOT_HISTO1D ) {

	// Histogram size
	_plot.set_histogram_n( gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(_prefdata->histo_n_spin) ) );

	// Histogram accumulation type
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->histo_acc_closest_radio) ) )
	    _plot.set_histogram_accumulation( HISTOGRAM_ACCUMULATION_CLOSEST );
	else
	    _plot.set_histogram_accumulation( HISTOGRAM_ACCUMULATION_LINEAR );

	// Plot presentation style
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->style_histo_radio) ) )
	    _plot.set_histogram_style( true );
	else
	    _plot.set_histogram_style( false );

    } else {

	particle_diag_plot_type_e type;
	particle_diag_plot_type_e otype = _plot.get_type();
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->plot_scatter_radio) ) )
	    type = PARTICLE_DIAG_PLOT_SCATTER;
	else
	    type = PARTICLE_DIAG_PLOT_HISTO2D;
	if( type != otype )
	    _plot.set_type( type );

	// Histogram size
	_plot.set_histogram_n( gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(_prefdata->histo_n_spin) ) );
	_plot.set_histogram_m( gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(_prefdata->histo_m_spin) ) );

	// Histogram accumulation type
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->histo_acc_closest_radio) ) )
	    _plot.set_histogram_accumulation( HISTOGRAM_ACCUMULATION_CLOSEST );
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->histo_acc_bilinear_radio) ) )
	    _plot.set_histogram_accumulation( HISTOGRAM_ACCUMULATION_LINEAR );

	// Interpolation style
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->int_closest_radio) ) )
	    _plot.set_colormap_interpolation( INTERPOLATION_CLOSEST );
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->int_bilinear_radio) ) )
	    _plot.set_colormap_interpolation( INTERPOLATION_BILINEAR );
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->int_bicubic_radio) ) )
	    _plot.set_colormap_interpolation( INTERPOLATION_BICUBIC );

	// Dot size
	//_plot.set_dot_size( gtk_spin_button_get_value_as_float( GTK_SPIN_BUTTON(_prefdata->dot_size_spin) ) );
	_plot.set_dot_size( gtk_spin_button_get_value( GTK_SPIN_BUTTON(_prefdata->dot_size_spin) ) );

	// Ellipse fit
	_plot.set_emittance_ellipse( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( _prefdata->ellipse_check) ) );
    }

    // Range setting, draw and expose done by calling function GTKPreferences::run()
    //zoom_fit();
    _plot.build_plot();
    //draw_and_expose();
}


void GTKParticleDiagWindow::export_data( void )
{
    GTKParticleDiagExportDialog dialog( _window, &_plot );
    dialog.run();
}


void GTKParticleDiagWindow::menuitem_export_signal( GtkToolButton *button,
						    gpointer object )
{
    GTKParticleDiagWindow *window = (GTKParticleDiagWindow *)object;
    window->export_data();
}



std::string GTKParticleDiagWindow::track_text( double x, double y )
{
    std::wstringstream ss;
    std::locale mylocale("");
    ss.imbue( mylocale );

    particle_diag_plot_type_e type;
    trajectory_diagnostic_e diagx, diagy;
    _plot.get_plot( type, diagx, diagy );

    for( int i = 0; i < 2; i++ ) {

	double val;
	trajectory_diagnostic_e diag;
	if( i == 0 ) {
	    val = x;
	    diag = diagx;
	} else {
	    val = y;
	    diag = diagy;
	}

	if( type == PARTICLE_DIAG_PLOT_HISTO1D && i == 1 ) {
	    ss << "Intensity = " << val << "\n";
	} else {
	    ss << trajectory_diagnostic_string[diag] << " = " 
	       << val << " " 
	       << trajectory_diagnostic_string_unit[diag] << "\n";
	}
    }

    if( type == PARTICLE_DIAG_PLOT_HISTO2D ) {
	if( (diagx == DIAG_X || diagx == DIAG_Y || diagx == DIAG_R || diagx == DIAG_Z) && 
	    (diagy == DIAG_X || diagy == DIAG_Y || diagy == DIAG_R || diagy == DIAG_Z) ) {
	    // Profile plot
	    const MeshColormap *cmap = _plot.get_colormap();
	    double val = cmap->get_value( x, y );
	    ss << "J = " << val << " A/m2\n";
	} else {
	    // Emittance plot
	    const MeshColormap *cmap = _plot.get_colormap();
	    double val = cmap->get_value( x, y );
	    ss << "J = " << val << " A/(m rad)\n";
	}
    }

    return( wstring_to_utf8( ss.str() ) );
}






