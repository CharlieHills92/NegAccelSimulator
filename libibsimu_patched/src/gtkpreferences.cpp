/*! \file gtkpreferences.cpp
 *  \brief Preferences for plot windows
 */

/* Copyright (c) 2005-2009,2012-2013,2022 Taneli Kalvas. All rights reserved.
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

#include "gtkpreferences.hpp"
#include "gtkframewindow.hpp"


GTKPreferences::GTKPreferences( GTKFrameWindow *gtkwindow, GtkWidget *window, Frame *frame )
    : _gtkwindow(gtkwindow), _window(window), _frame(frame)
{

}


GTKPreferences::~GTKPreferences()
{

}


void GTKPreferences::run( void )
{
    GtkWidget *dialog = gtk_dialog_new_with_buttons( "Plot preferences",
						     GTK_WINDOW(_window),
						     (GtkDialogFlags)(GTK_DIALOG_MODAL | 
								      GTK_DIALOG_DESTROY_WITH_PARENT), 
						     dgettext("gtk30","_Ok"), GTK_RESPONSE_ACCEPT,
						     dgettext("gtk30","_Cancel"), GTK_RESPONSE_REJECT,
						     NULL );
    gtk_window_set_resizable( GTK_WINDOW(dialog), FALSE );
    GtkWidget *mainbox = gtk_dialog_get_content_area( GTK_DIALOG(dialog) );

    // ****************************************************************************

    GtkWidget *grid = gtk_grid_new();
    gtk_grid_set_column_homogeneous( GTK_GRID(grid), TRUE );
    gtk_widget_set_hexpand( grid, TRUE );
    uint32_t yl = 0;

    // ****************************************************************************

    // Fontsize
    GtkWidget *label = gtk_label_new( "Fontsize" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    GtkWidget *fontsize_entry = gtk_entry_new();
    char s[128];
    snprintf( s, 128, "%g", _frame->get_font_size() );
    gtk_entry_set_text( GTK_ENTRY(fontsize_entry), s );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), fontsize_entry, 1, yl, 1, 1 );
    yl++;

    // Range xmin
    label = gtk_label_new( "Range xmin" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    GtkWidget *rxmin_entry = gtk_entry_new();
    double min, max;
    _frame->get_ranges( PLOT_AXIS_X1, min, max );
    snprintf( s, 128, "%g", min );
    gtk_entry_set_text( GTK_ENTRY(rxmin_entry), s );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), rxmin_entry, 1, yl, 1, 1 );
    yl++;

    // Range xmax
    label = gtk_label_new( "Range xmax" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    GtkWidget *rxmax_entry = gtk_entry_new();
    _frame->get_ranges( PLOT_AXIS_X1, min, max );
    snprintf( s, 128, "%g", max );
    gtk_entry_set_text( GTK_ENTRY(rxmax_entry), s );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), rxmax_entry, 1, yl, 1, 1 );
    yl++;

    // Range ymin
    label = gtk_label_new( "Range ymin" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    GtkWidget *rymin_entry = gtk_entry_new();
    _frame->get_ranges( PLOT_AXIS_Y1, min, max );
    snprintf( s, 128, "%g", min );
    gtk_entry_set_text( GTK_ENTRY(rymin_entry), s );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), rymin_entry, 1, yl, 1, 1 );
    yl++;
    
    // Range ymax
    label = gtk_label_new( "Range ymax" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    GtkWidget *rymax_entry = gtk_entry_new();
    _frame->get_ranges( PLOT_AXIS_Y1, min, max );
    snprintf( s, 128, "%g", max );
    gtk_entry_set_text( GTK_ENTRY(rymax_entry), s );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), rymax_entry, 1, yl, 1, 1 );
    yl++;

    // Notebook, page 1
    label = gtk_label_new( "Frame" );
    GtkWidget *notebook = gtk_notebook_new();
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), grid, label );

    // Notebook, additional pages
    void *pdata = _gtkwindow->build_preferences( notebook );

    // ****************************************************************************

    // Pack notebook
    gtk_box_pack_start( GTK_BOX(mainbox), notebook, FALSE, TRUE, 0 );

    // ****************************************************************************

    gtk_widget_show_all( dialog );
    if( gtk_dialog_run( GTK_DIALOG(dialog) ) == GTK_RESPONSE_ACCEPT ) {

	// Fontsize
	double fontsize = atof( gtk_entry_get_text( GTK_ENTRY(fontsize_entry) ) );
	_frame->set_font_size( fontsize );

	// Ranges
	double rxmin = atof( gtk_entry_get_text( GTK_ENTRY(rxmin_entry) ) );
	double rxmax = atof( gtk_entry_get_text( GTK_ENTRY(rxmax_entry) ) );
	double rymin = atof( gtk_entry_get_text( GTK_ENTRY(rymin_entry) ) );
	double rymax = atof( gtk_entry_get_text( GTK_ENTRY(rymax_entry) ) );
	_frame->set_ranges( PLOT_AXIS_X1, rxmin, rxmax );
	_frame->set_ranges( PLOT_AXIS_Y1, rymin, rymax );

	// Read additional pages of notebook
	_gtkwindow->read_preferences( notebook, pdata );

	// Refresh frame
	_gtkwindow->draw_and_expose();
    }

    gtk_widget_destroy( dialog );
}

