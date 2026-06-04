/*! \file gtkfielddiagdialog.cpp
 *  \brief Dialog for constructing field diagnostic windows
 */

/* Copyright (c) 2005-2013,2022 Taneli Kalvas. All rights reserved.
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

#include "gtkfielddiagdialog.hpp"
#include "meshvectorfield.hpp"
#include "multimeshvectorfield.hpp"
#include "ibsimu.hpp"


GTKFieldDiagDialog::GTKFieldDiagDialog( GtkWidget *window, GTKPlotter &plotter, 
					double x1[3], double x2[3] )
    : _window(window), _plotter(&plotter)
{
    _x1[0] = x1[0];
    _x1[1] = x1[1];
    _x1[2] = x1[2];
    _x2[0] = x2[0];
    _x2[1] = x2[1];
    _x2[2] = x2[2];

    _geom = _plotter->get_geometry();
}


GTKFieldDiagDialog::~GTKFieldDiagDialog()
{

}


void GTKFieldDiagDialog::run( void )
{
    GtkWidget *dialog = gtk_dialog_new_with_buttons( "Make field diagnostics",
                                                     GTK_WINDOW(_window),
                                                     (GtkDialogFlags)(GTK_DIALOG_MODAL | 
								      GTK_DIALOG_DESTROY_WITH_PARENT), 
                                                     dgettext("gtk30","_OK"), GTK_RESPONSE_ACCEPT,
                                                     dgettext("gtk30","_Cancel"), GTK_RESPONSE_REJECT,
                                                     NULL );
    gtk_window_set_resizable( GTK_WINDOW(dialog), FALSE );
    GtkWidget *mainbox = gtk_dialog_get_content_area( GTK_DIALOG(dialog) );

    // ****************************************************************************

    GtkWidget *grid = gtk_grid_new();

    // Labels
    GtkWidget *label = gtk_label_new( "" );
    gtk_label_set_markup( GTK_LABEL(label), "<span weight=\"bold\">Start</span>" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, 0, 2, 1 );
    label = gtk_label_new( "" );
    gtk_label_set_markup( GTK_LABEL(label), "<span weight=\"bold\">End</span>" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 2, 0, 2, 1 );

    // x1
    label = gtk_label_new( "x1" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, 1, 1, 1 );
    GtkWidget *entry_x1 = gtk_entry_new();
    gtk_entry_set_max_length( GTK_ENTRY(entry_x1), 30 );
    gtk_widget_set_margin_start( entry_x1, 5 );
    gtk_widget_set_margin_end( entry_x1, 5 );
    char s[128];
    snprintf( s, 128, "%g", _x1[0] );
    gtk_entry_set_text( GTK_ENTRY(entry_x1), s );
    gtk_grid_attach( GTK_GRID(grid), entry_x1, 1, 1, 1, 1 );

    // x2
    label = gtk_label_new( "x2" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 2, 1, 1, 1 );
    GtkWidget *entry_x2 = gtk_entry_new();
    gtk_entry_set_max_length( GTK_ENTRY(entry_x2), 30 );
    gtk_widget_set_margin_start( entry_x2, 5 );
    snprintf( s, 128, "%g", _x2[0] );
    gtk_entry_set_text( GTK_ENTRY(entry_x2), s );
    gtk_grid_attach( GTK_GRID(grid), entry_x2, 3, 1, 1, 1 );

    // y1
    if( _geom->geom_mode() == MODE_CYL )
	label = gtk_label_new( "r1" );
    else
	label = gtk_label_new( "y1" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, 2, 1, 1 );
    GtkWidget *entry_y1 = gtk_entry_new();
    gtk_entry_set_max_length( GTK_ENTRY(entry_y1), 30 );
    gtk_widget_set_margin_start( entry_y1, 5 );
    gtk_widget_set_margin_end( entry_y1, 5 );
    snprintf( s, 128, "%g", _x1[1] );
    gtk_entry_set_text( GTK_ENTRY(entry_y1), s );
    gtk_grid_attach( GTK_GRID(grid), entry_y1, 1, 2, 1, 1 );

    // y2
    if( _geom->geom_mode() == MODE_CYL )
	label = gtk_label_new( "r2" );
    else
	label = gtk_label_new( "y2" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 2, 2, 1, 1 );
    GtkWidget *entry_y2 = gtk_entry_new();
    gtk_entry_set_max_length( GTK_ENTRY(entry_y2), 30 );
    gtk_widget_set_margin_start( entry_y2, 5 );
    snprintf( s, 128, "%g", _x2[1] );
    gtk_entry_set_text( GTK_ENTRY(entry_y2), s );
    gtk_grid_attach( GTK_GRID(grid), entry_y2, 3, 2, 1, 1 );

    // z1
    label = gtk_label_new( "z1" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, 3, 1, 1 );
    GtkWidget *entry_z1 = gtk_entry_new();
    gtk_entry_set_max_length( GTK_ENTRY(entry_z1), 30 );
    gtk_widget_set_margin_start( entry_z1, 5 );
    gtk_widget_set_margin_end( entry_z1, 5 );
    if( _geom->geom_mode() == MODE_2D || _geom->geom_mode() == MODE_CYL )
	gtk_widget_set_sensitive( entry_z1, FALSE );
    snprintf( s, 128, "%g", _x1[2] );
    gtk_entry_set_text( GTK_ENTRY(entry_z1), s );
    gtk_grid_attach( GTK_GRID(grid), entry_z1, 1, 3, 1, 1 );

    // z2
    label = gtk_label_new( "z2" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 2, 3, 1, 1 );
    GtkWidget *entry_z2 = gtk_entry_new();
    gtk_entry_set_max_length( GTK_ENTRY(entry_z2), 30 );
    gtk_widget_set_margin_start( entry_z2, 5 );
    if( _geom->geom_mode() == MODE_2D || _geom->geom_mode() == MODE_CYL )
	gtk_widget_set_sensitive( entry_z2, FALSE );
    snprintf( s, 128, "%g", _x2[2] );
    gtk_entry_set_text( GTK_ENTRY(entry_z2), s );
    gtk_grid_attach( GTK_GRID(grid), entry_z2, 3, 3, 1, 1 );

    // N
    label = gtk_label_new( "Samples" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, 4, 2, 1 );
    GtkAdjustment *n_adj = gtk_adjustment_new( 100, 0, 10000, 10, 100, 0 );
    GtkWidget *n_spin = gtk_spin_button_new( n_adj, 10, 0 );
    gtk_grid_attach( GTK_GRID(grid), n_spin, 2, 4, 2, 1 );

    gtk_box_pack_start( GTK_BOX(mainbox), grid, FALSE, TRUE, 0 );

    // ****************************************************************************

    // Rest of the dialog is in one grid with three columns
    grid = gtk_grid_new();
    gtk_grid_set_column_homogeneous( GTK_GRID(grid), TRUE );
    gtk_widget_set_hexpand( grid, TRUE );
    uint32_t yl = 0;

    // ****************************************************************************

    // Separator
    GtkWidget *separator = gtk_separator_new( GTK_ORIENTATION_HORIZONTAL );
    gtk_widget_set_margin_start( separator, 5 );
    gtk_widget_set_margin_end( separator, 5 );
    gtk_widget_set_margin_top( separator, 5 );
    gtk_widget_set_margin_bottom( separator, 5 );
    gtk_grid_attach( GTK_GRID(grid), separator, 0, yl, 3, 1 );
    yl++;

    // ****************************************************************************
    // Axis selection titles

    label = gtk_label_new( "Axis x1" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 1, yl, 1, 1 );
    label = gtk_label_new( "Axis x2" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 2, yl, 1, 1 );
    yl++;

    // Distance
    label = gtk_label_new( "Distance" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_axis1_dist = gtk_radio_button_new_from_widget( NULL );
    GtkWidget *radio_axis2_dist = gtk_radio_button_new_from_widget( NULL );
    gtk_grid_attach( GTK_GRID(grid), radio_axis1_dist, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_axis2_dist, 2, yl, 1, 1 );
    yl++;

    // X
    label = gtk_label_new( "X" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_axis1_x = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_axis1_dist) );
    GtkWidget *radio_axis2_x = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_axis2_dist) );
    gtk_grid_attach( GTK_GRID(grid), radio_axis1_x, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_axis2_x, 2, yl, 1, 1 );
    yl++;

    // Y
    if(_geom->geom_mode() == MODE_CYL )
	label = gtk_label_new( "R" );
    else
	label = gtk_label_new( "Y" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_axis1_y = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_axis1_dist) );
    GtkWidget *radio_axis2_y = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_axis2_dist) );
    gtk_grid_attach( GTK_GRID(grid), radio_axis1_y, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_axis2_y, 2, yl, 1, 1 );
    yl++;

    // Z
    label = gtk_label_new( "Z" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_axis1_z = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_axis1_dist) );
    GtkWidget *radio_axis2_z = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_axis2_dist) );
    if( _geom->geom_mode() == MODE_2D || _geom->geom_mode() == MODE_CYL ) {
	gtk_widget_set_sensitive( radio_axis1_z, FALSE );
	gtk_widget_set_sensitive( radio_axis2_z, FALSE );
    }
    gtk_grid_attach( GTK_GRID(grid), radio_axis1_z, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_axis2_z, 2, yl, 1, 1 );
    yl++;

    // None
    label = gtk_label_new( "None" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_axis1_none = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_axis1_dist) );
    GtkWidget *radio_axis2_none = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_axis2_dist) );
    gtk_grid_attach( GTK_GRID(grid), radio_axis1_none, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_axis2_none, 2, yl, 1, 1 );
    yl++;

    // ****************************************************************************

    // Separator
    separator = gtk_separator_new( GTK_ORIENTATION_HORIZONTAL );
    gtk_widget_set_margin_start( separator, 5 );
    gtk_widget_set_margin_end( separator, 5 );
    gtk_widget_set_margin_top( separator, 5 );
    gtk_widget_set_margin_bottom( separator, 5 );
    gtk_grid_attach( GTK_GRID(grid), separator, 0, yl, 3, 1 );
    yl++;

    // ****************************************************************************

    // Field selection titles
    label = gtk_label_new( "Graph 1" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_widget_set_margin_start( label, 5 );
    gtk_widget_set_margin_end( label, 5 );
    gtk_grid_attach( GTK_GRID(grid), label, 1, yl, 1, 1 );
    label = gtk_label_new( "Graph 2" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_widget_set_margin_start( label, 5 );
    gtk_grid_attach( GTK_GRID(grid), label, 2, yl, 1, 1 );
    yl++;

    // Epot
    label = gtk_label_new( "Electric potential" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_g1_epot = gtk_radio_button_new_from_widget( NULL );
    GtkWidget *radio_g2_epot = gtk_radio_button_new_from_widget( NULL );
    if( !_plotter->get_epot() ) {
	gtk_widget_set_sensitive( radio_g1_epot, FALSE );
	gtk_widget_set_sensitive( radio_g2_epot, FALSE );
    }
    gtk_grid_attach( GTK_GRID(grid), radio_g1_epot, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_g2_epot, 2, yl, 1, 1 );
    yl++;

    // E
    label = gtk_label_new( "|E|" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_g1_e = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g1_epot) );
    GtkWidget *radio_g2_e = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g2_epot) );
    if( !_plotter->get_efield() ) {
	gtk_widget_set_sensitive( radio_g1_e, FALSE );
	gtk_widget_set_sensitive( radio_g2_e, FALSE );
    }
    gtk_grid_attach( GTK_GRID(grid), radio_g1_e, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_g2_e, 2, yl, 1, 1 );
    yl++;
    
    // Ex
    label = gtk_label_new( "Ex" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_g1_ex = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g1_epot) );
    GtkWidget *radio_g2_ex = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g2_epot) );
    if( !_plotter->get_efield() ) {
	gtk_widget_set_sensitive( radio_g1_ex, FALSE );
	gtk_widget_set_sensitive( radio_g2_ex, FALSE );
    }
    gtk_grid_attach( GTK_GRID(grid), radio_g1_ex, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_g2_ex, 2, yl, 1, 1 );
    yl++;

    // Ey
    if( _geom->geom_mode() == MODE_CYL )
	label = gtk_label_new( "Er" );
    else
	label = gtk_label_new( "Ey" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_g1_ey = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g1_epot) );
    GtkWidget *radio_g2_ey = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g2_epot) );
    if( !_plotter->get_efield() ) {
	gtk_widget_set_sensitive( radio_g1_ey, FALSE );
	gtk_widget_set_sensitive( radio_g2_ey, FALSE );
    }
    gtk_grid_attach( GTK_GRID(grid), radio_g1_ey, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_g2_ey, 2, yl, 1, 1 );
    yl++;

    // Ez
    label = gtk_label_new( "Ez" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_g1_ez = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g1_epot) );
    GtkWidget *radio_g2_ez = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g2_epot) );
    if( !_plotter->get_epot() || _geom->geom_mode() == MODE_2D || _geom->geom_mode() == MODE_CYL ||
	!_plotter->get_efield() ) {
	gtk_widget_set_sensitive( radio_g1_ez, FALSE );
	gtk_widget_set_sensitive( radio_g2_ez, FALSE );
    }
    gtk_grid_attach( GTK_GRID(grid), radio_g1_ez, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_g2_ez, 2, yl, 1, 1 );
    yl++;
    
    // Space charge
    label = gtk_label_new( "Space charge" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_g1_scharge = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g1_epot) );
    GtkWidget *radio_g2_scharge = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g2_epot) );
    if( !_plotter->get_scharge() ) {
	gtk_widget_set_sensitive( radio_g1_scharge, FALSE );
	gtk_widget_set_sensitive( radio_g2_scharge, FALSE );
    }
    gtk_grid_attach( GTK_GRID(grid), radio_g1_scharge, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_g2_scharge, 2, yl, 1, 1 );
    yl++;

    // Bfield components available
    bool bfield_fout[3] = {false, false, false};
    if( _plotter->get_bfield() ) {
	const VectorField *bfield = _plotter->get_bfield();
	const MeshVectorField *meshbfield = dynamic_cast<const MeshVectorField *>( bfield );
	const MultiMeshVectorField *multimeshbfield = dynamic_cast<const MultiMeshVectorField *>( bfield );
	if( meshbfield != NULL )
	    meshbfield->get_defined_components( bfield_fout );
	else if( multimeshbfield != NULL )
	    multimeshbfield->get_defined_components( bfield_fout );
	else {
	    bfield_fout[0] = bfield_fout[1] = bfield_fout[2] = true;
	}
    }

    // B
    label = gtk_label_new( "|B|" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_g1_b = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g1_epot) );
    if( !_plotter->get_bfield() )
	gtk_widget_set_sensitive( radio_g1_b, FALSE );
    GtkWidget *radio_g2_b = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g2_epot) );
    if( !_plotter->get_bfield() )
	gtk_widget_set_sensitive( radio_g2_b, FALSE );
    gtk_grid_attach( GTK_GRID(grid), radio_g1_b, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_g2_b, 2, yl, 1, 1 );
    yl++;

    // Bx
    label = gtk_label_new( "Bx" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_g1_bx = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g1_epot) );
    if( !bfield_fout[0] )
	gtk_widget_set_sensitive( radio_g1_bx, FALSE );
    GtkWidget *radio_g2_bx = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g2_epot) );
    if( !bfield_fout[0] )
	gtk_widget_set_sensitive( radio_g2_bx, FALSE );
    gtk_grid_attach( GTK_GRID(grid), radio_g1_bx, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_g2_bx, 2, yl, 1, 1 );
    yl++;

    // By
    if( _geom->geom_mode() == MODE_CYL )
	label = gtk_label_new( "Br" );
    else
	label = gtk_label_new( "By" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_g1_by = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g1_epot) );
    if( !bfield_fout[1] )
	gtk_widget_set_sensitive( radio_g1_by, FALSE );
    GtkWidget *radio_g2_by = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g2_epot) );
    if( !bfield_fout[1] )
	gtk_widget_set_sensitive( radio_g2_by, FALSE );
    gtk_grid_attach( GTK_GRID(grid), radio_g1_by, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_g2_by, 2, yl, 1, 1 );
    yl++;

    // Bz
    if( _geom->geom_mode() == MODE_CYL )
	label = gtk_label_new( "Btheta" );
    else
	label = gtk_label_new( "Bz" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_g1_bz = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g1_epot) );
    if( !bfield_fout[2] )
	gtk_widget_set_sensitive( radio_g1_bz, FALSE );
    GtkWidget *radio_g2_bz = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g2_epot) );
    if( !bfield_fout[2] )
	gtk_widget_set_sensitive( radio_g2_bz, FALSE );
    gtk_grid_attach( GTK_GRID(grid), radio_g1_bz, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_g2_bz, 2, yl, 1, 1 );
    yl++;

    // None
    label = gtk_label_new( "None" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    GtkWidget *radio_g1_none = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g1_epot) );
    GtkWidget *radio_g2_none = gtk_radio_button_new_from_widget( GTK_RADIO_BUTTON(radio_g2_epot) );
    if( !_plotter->get_epot() ) {
	// Default to none if epot not available
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(radio_g1_none), TRUE );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(radio_g2_none), TRUE );
    }
    gtk_grid_attach( GTK_GRID(grid), radio_g1_none, 1, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), radio_g2_none, 2, yl, 1, 1 );
    yl++;

    // ****************************************************************************

    // Separator
    separator = gtk_separator_new( GTK_ORIENTATION_HORIZONTAL );
    gtk_widget_set_margin_start( separator, 5 );
    gtk_widget_set_margin_end( separator, 5 );
    gtk_widget_set_margin_top( separator, 5 );
    gtk_widget_set_margin_bottom( separator, 5 );
    gtk_grid_attach( GTK_GRID(grid), separator, 0, yl, 3, 1 );
    yl++;

    // ****************************************************************************

    gtk_box_pack_start( GTK_BOX(mainbox), grid, TRUE, TRUE, 0 );

    // ****************************************************************************

    gtk_widget_show_all( dialog );

    if( gtk_dialog_run( GTK_DIALOG(dialog) ) == GTK_RESPONSE_ACCEPT ) {

	// Read in options and start particle diagnostics plot

	// Read coordinates
	Vec3D x1;
	Vec3D x2;

	const char *entry_text = gtk_entry_get_text( GTK_ENTRY(entry_x1) );
	x1[0] = atof(entry_text);
	entry_text = gtk_entry_get_text( GTK_ENTRY(entry_y1) );
	x1[1] = atof(entry_text);
	entry_text = gtk_entry_get_text( GTK_ENTRY(entry_z1) );
	x1[2] = atof(entry_text);

	entry_text = gtk_entry_get_text( GTK_ENTRY(entry_x2) );
	x2[0] = atof(entry_text);
	entry_text = gtk_entry_get_text( GTK_ENTRY(entry_y2) );
	x2[1] = atof(entry_text);
	entry_text = gtk_entry_get_text( GTK_ENTRY(entry_z2) );
	x2[2] = atof(entry_text);

	// Read N
	size_t N = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(n_spin) );
	
	// Read distance axis selection
	field_loc_type_e loc[2];
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_axis1_dist) ) )
	    loc[0] = FIELDD_LOC_DIST;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_axis1_x) ) )
	    loc[0] = FIELDD_LOC_X;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_axis1_y) ) )
	    loc[0] = FIELDD_LOC_Y;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_axis1_z) ) )
	    loc[0] = FIELDD_LOC_Z;
	else
	    loc[0] = FIELDD_LOC_NONE;

	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_axis2_dist) ) )
	    loc[1] = FIELDD_LOC_DIST;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_axis2_x) ) )
	    loc[1] = FIELDD_LOC_X;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_axis2_y) ) )
	    loc[1] = FIELDD_LOC_Y;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_axis2_z) ) )
	    loc[1] = FIELDD_LOC_Z;
	else
	    loc[1] = FIELDD_LOC_NONE;

	// Read plots to be made
	field_diag_type_e diag[2];
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g1_epot) ) )
	    diag[0] = FIELD_EPOT;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g1_e) ) )
	    diag[0] = FIELD_EFIELD;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g1_ex) ) )
	    diag[0] = FIELD_EFIELD_X;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g1_ey) ) )
	    diag[0] = FIELD_EFIELD_Y;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g1_ez) ) )
	    diag[0] = FIELD_EFIELD_Z;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g1_scharge) ) )
	    diag[0] = FIELD_SCHARGE;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g1_b) ) )
	    diag[0] = FIELD_BFIELD;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g1_bx) ) )
	    diag[0] = FIELD_BFIELD_X;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g1_by) ) )
	    diag[0] = FIELD_BFIELD_Y;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g1_bz) ) )
	    diag[0] = FIELD_BFIELD_Z;
	else
	    diag[0] = FIELD_NONE;
	
	if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g2_epot) ) )
	    diag[1] = FIELD_EPOT;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g2_e) ) )
	    diag[1] = FIELD_EFIELD;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g2_ex) ) )
	    diag[1] = FIELD_EFIELD_X;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g2_ey) ) )
	    diag[1] = FIELD_EFIELD_Y;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g2_ez) ) )
	    diag[1] = FIELD_EFIELD_Z;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g2_scharge) ) )
	    diag[1] = FIELD_SCHARGE;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g2_b) ) )
	    diag[1] = FIELD_BFIELD;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g2_bx) ) )
	    diag[1] = FIELD_BFIELD_X;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g2_by) ) )
	    diag[1] = FIELD_BFIELD_Y;
	else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(radio_g2_bz) ) )
	    diag[1] = FIELD_BFIELD_Z;
	else
	    diag[1] = FIELD_NONE;

	ibsimu.message(MSG_DEBUG_GENERAL,2) << "Making field diagnostic plot\n"
					    << "  x1 = " << x1[0] << "\n"
					    << "  y1 = " << x1[1] << "\n"
					    << "  z1 = " << x1[2] << "\n"
					    << "  x2 = " << x2[0] << "\n"
					    << "  y2 = " << x2[1] << "\n"
					    << "  z2 = " << x2[2] << "\n"
					    << "  N  = " << N << "\n"
					    << "  g1 = " << diag[0] << "\n"
					    << "  g2 = " << diag[1] << "\n"
					    << "  d1 = " << loc[0] << "\n"
					    << "  d2 = " << loc[1] << "\n";

	try {
	    _plotter->new_field_plot_window( N, x1, x2, diag, loc );
	} catch( Error e ) {
	    //GTKErrorDialog error( _window, e );
	    e.print_error_message( std::cout );
	    std::cout << "Trying to continue\n";
	}
    }

    gtk_widget_destroy( dialog );
}
