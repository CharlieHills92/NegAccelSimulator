/*! \file gtkgeomwindow.cpp
 *  \brief %Geometry view window
 */

/* Copyright (c) 2005-2009, 2011-2013, 2016, 2022 Taneli Kalvas. All rights reserved.
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
#include <locale>
#include "gtkgeomwindow.hpp"
#include "gtkparticlediagdialog.hpp"
#include "gtkfielddiagdialog.hpp"
#include "epot_efield.hpp"
#include "icons.hpp"
#include "ibsimu.hpp"
#include "fonts.hpp"


#define TOOL_UNKNOWN  -1
#define TOOL_PARTICLE_DIAG 4
#define TOOL_FIELD_DIAG 5


GTKGeomWindow::GTKGeomWindow( class GTKPlotter       &plotter, 
			      const Geometry         &geom,
			      const EpotField        *epot,
			      const EpotEfield       *efield,
			      const MeshScalarField  *scharge,
			      const MeshScalarField  *tdens,
			      const VectorField      *bfield,
			      const ParticleDataBase *pdb )
  : GTKFrameWindow(plotter), _geomplot(_frame,geom), 
    _geom(geom), _epot(epot), _efield(efield), _scharge(scharge), _tdens(tdens), 
    _bfield(bfield), _pdb(pdb), _prefdata(NULL)
{
    //std::cout << "GTKGeomWindow constructor\n";

    // Setup GeomPlot
    _geomplot.set_epot( epot );
    _geomplot.set_scharge( scharge );
    _geomplot.set_trajdens( tdens );
    _geomplot.set_bfield( bfield );
    _geomplot.set_efield( efield );
    _geomplot.set_particle_database( pdb );

    // Set window title
    gtk_window_set_title( GTK_WINDOW(_window), "Simulation geometry" );

    // Adding geometry window specific tools to toolbar
    // Creating separator
    GtkToolItem *toolitem = gtk_separator_tool_item_new();
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );

    // Creating "Particle diagnostics" button
    GdkPixbuf *pixbuf = gdk_pixbuf_new_from_inline( -1, icon_particle_diag_inline, FALSE, NULL );
    GtkWidget *icon = gtk_image_new_from_pixbuf( pixbuf );
    toolitem = gtk_radio_tool_button_new_from_widget( GTK_RADIO_TOOL_BUTTON(_radioitem) );
    gtk_tool_button_set_label( GTK_TOOL_BUTTON(toolitem), "Particle diagnostics" );
#if GTK_CHECK_VERSION(2,12,0)
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Particle diagnostics" );
#endif
    gtk_tool_button_set_icon_widget( GTK_TOOL_BUTTON(toolitem), icon );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "toggled",
		      G_CALLBACK(menuitem_tool_change_signal),
		      (gpointer)this );

    // Creating "Field diagnostics" button
    pixbuf = gdk_pixbuf_new_from_inline( -1, icon_field_diag_inline, FALSE, NULL );
    icon = gtk_image_new_from_pixbuf( pixbuf );
    toolitem = gtk_radio_tool_button_new_from_widget( GTK_RADIO_TOOL_BUTTON(_radioitem) );
    gtk_tool_button_set_label( GTK_TOOL_BUTTON(toolitem), "Field diagnostics" );
#if GTK_CHECK_VERSION(2,12,0)
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Field diagnostics" );
#endif
    gtk_tool_button_set_icon_widget( GTK_TOOL_BUTTON(toolitem), icon );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    g_signal_connect( G_OBJECT(toolitem), "toggled",
		      G_CALLBACK(menuitem_tool_change_signal),
		      (gpointer)this );

    // Creating "Geom 3D" button
    pixbuf = gdk_pixbuf_new_from_inline( -1, icon_geom3d_inline, FALSE, NULL );
    icon = gtk_image_new_from_pixbuf( pixbuf );
    toolitem = gtk_tool_button_new( icon, "3D geometry view" );
#if GTK_CHECK_VERSION(2,12,0)
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "3D geometry view" );
#endif
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
    if( _geom.geom_mode() != MODE_3D || !_geom.surface_built() )
	gtk_widget_set_sensitive( GTK_WIDGET(toolitem), FALSE );
    g_signal_connect( G_OBJECT(toolitem), "clicked",
		      G_CALLBACK(menuitem_geom3d_signal),
		      (gpointer)this );

    // Creating separator
    toolitem = gtk_separator_tool_item_new();
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );

    // Creating view combobox and level spinbutton
    //_combobox = gtk_combo_box_new_text();
    _combobox = gtk_combo_box_text_new();
    if( geom.geom_mode() == MODE_3D ) {
	gtk_combo_box_text_append( GTK_COMBO_BOX_TEXT(_combobox), NULL, "XY" );
	gtk_combo_box_text_append( GTK_COMBO_BOX_TEXT(_combobox), NULL, "XZ" );
        gtk_combo_box_text_append( GTK_COMBO_BOX_TEXT(_combobox), NULL, "YX" );
        gtk_combo_box_text_append( GTK_COMBO_BOX_TEXT(_combobox), NULL, "YZ" );
        gtk_combo_box_text_append( GTK_COMBO_BOX_TEXT(_combobox), NULL, "ZX" );
        gtk_combo_box_text_append( GTK_COMBO_BOX_TEXT(_combobox), NULL, "ZY" );
	/*
	gtk_combo_box_append_text( GTK_COMBO_BOX(_combobox), "XY" );
	gtk_combo_box_append_text( GTK_COMBO_BOX(_combobox), "XZ" );
        gtk_combo_box_append_text( GTK_COMBO_BOX(_combobox), "YX" );
        gtk_combo_box_append_text( GTK_COMBO_BOX(_combobox), "YZ" );
        gtk_combo_box_append_text( GTK_COMBO_BOX(_combobox), "ZX" );
        gtk_combo_box_append_text( GTK_COMBO_BOX(_combobox), "ZY" );
	*/
	gtk_combo_box_set_active( GTK_COMBO_BOX(_combobox), 0 );
	_spinbutton = gtk_spin_button_new_with_range( 0, geom.size(2)-1, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON(_spinbutton), geom.size(2)/2 );
	_geomplot.set_view( VIEW_XY, geom.size(2)/2 );
    } else {
	gtk_combo_box_text_append( GTK_COMBO_BOX_TEXT(_combobox), NULL, "XY" );
	gtk_combo_box_text_append( GTK_COMBO_BOX_TEXT(_combobox), NULL, "YX" );
        gtk_combo_box_set_active( GTK_COMBO_BOX(_combobox), 0 );
        _spinbutton = gtk_spin_button_new_with_range( 0, 0, 1 );
        gtk_spin_button_set_value( GTK_SPIN_BUTTON(_spinbutton), 0 );
    }
    toolitem = gtk_tool_item_new();
    gtk_container_add( GTK_CONTAINER(toolitem), _combobox );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
#if GTK_CHECK_VERSION(2,12,0)
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Select view direction" );
#endif

    gtk_spin_button_set_digits( GTK_SPIN_BUTTON(_spinbutton), 0 );
    toolitem = gtk_tool_item_new();
    gtk_container_add( GTK_CONTAINER(toolitem), _spinbutton );
    gtk_toolbar_insert( GTK_TOOLBAR(_toolbar), toolitem, -1 );
#if GTK_CHECK_VERSION(2,12,0)
    gtk_widget_set_tooltip_text( GTK_WIDGET(toolitem), "Select view level" );
#endif
    g_signal_connect( G_OBJECT(_combobox), "changed",
                      G_CALLBACK(combobox_signal),
                      (gpointer)this );
    g_signal_connect( G_OBJECT(_spinbutton), "value-changed",
                      G_CALLBACK(spinbutton_signal),
                      (gpointer)this );

    // Drawing area signals
    g_signal_connect( G_OBJECT(_darea), "draw",
		      G_CALLBACK(darea_draw_signal2), 
		      (gpointer)this );
    g_signal_connect( G_OBJECT(_darea), "button_press_event",
		      G_CALLBACK(darea_button_signal2),
		      (gpointer)this );
    g_signal_connect( G_OBJECT(_darea), "button_release_event",
		      G_CALLBACK(darea_button_signal2),
		      (gpointer)this );
    g_signal_connect( G_OBJECT(_darea), "motion_notify_event",
		      G_CALLBACK(darea_motion_signal2),
		      (gpointer)this );

    update_view();
    show();
}


void GTKGeomWindow::field_diag( int action, double x, double y )
{
    int x0, y0, width = 0, height = 0;

    if( action == 0 ) {
	// Start
	_start[0] = _end[0] = (int)floor(x+0.5);
	_start[1] = _end[1] = (int)floor(y+0.5);

    } else if( action == 1 || action == 2 ) {
	// Erase old line
	x0 = _start[0] < _end[0] ? _start[0] : _end[0];
	y0 = _start[1] < _end[1] ? _start[1] : _end[1];
	width = abs( _start[0] - _end[0] );
	height = abs( _start[1] - _end[1] );
	gtk_widget_queue_draw_area( _darea, x0-1, y0-1, width+3, height+3 );
    }

    if( action == 1 ) {
	// Draw new line
	_end[0] = (int)floor(x+0.5);
	_end[1] = (int)floor(y+0.5);
	_draw_action = 1;
	gtk_widget_queue_draw( _darea );

    } else if( action == 2 ) {
	_draw_action = 0;
	gtk_widget_queue_draw( _darea );
	
	// Done, start field diagnostics dialog
	Coordmapper cm = _frame.get_coordmapper( PLOT_AXIS_X1, PLOT_AXIS_Y1 );
	double x1[2] = {(double)_start[0], (double)_start[1]};
	double x2[2] = {(double)_end[0], (double)_end[1]};
	cm.inv_transform( x1[0], x1[1] );
	cm.inv_transform( x2[0], x2[1] );
	double p1[3];
	double p2[3];
	p1[_geomplot.vb(0)] = x1[0];
	p2[_geomplot.vb(0)] = x2[0];
	p1[_geomplot.vb(1)] = x1[1];
	p2[_geomplot.vb(1)] = x2[1];
	p1[_geomplot.vb(2)] = _geomplot.get_level_si();
	p2[_geomplot.vb(2)] = _geomplot.get_level_si();
	GTKFieldDiagDialog fielddiag( _window, _plotter, p1, p2 );
	fielddiag.run();
    }
}


void GTKGeomWindow::particle_diag( int action, double x, double y )
{
    int x0, y0, width = 0, height = 0;

    if( action == 0 ) {
	// Start
	_start[0] = _end[0] = (int)floor(x+0.5);
	_start[1] = _end[1] = (int)floor(y+0.5);
    } else if( action == 1 || action == 2 ) {
	// Erase old line
	x0 = _start[0] < _end[0] ? _start[0] : _end[0];
	y0 = _start[1] < _end[1] ? _start[1] : _end[1];
	width = abs( _start[0] - _end[0] );
	height = abs( _start[1] - _end[1] );
	if( width > height )
	    gtk_widget_queue_draw_area( _darea, x0-1, _start[1]-1, width+3, 3 );
	else
	    gtk_widget_queue_draw_area( _darea, _start[0]-1, y0-1, 3, height+3 );
    }

    if( action == 1 ) {
	// Draw new line
	_end[0] = (int)floor(x+0.5);
	_end[1] = (int)floor(y+0.5);
	x0 = _start[0] < _end[0] ? _start[0] : _end[0];
	y0 = _start[1] < _end[1] ? _start[1] : _end[1];
	width = abs( _start[0] - _end[0] );
	height = abs( _start[1] - _end[1] );
	_draw_action = 1;
	gtk_widget_queue_draw( _darea );
    }

    if( action == 2 ) {
	// If no particles
	if( _pdb == NULL ) {
	    GtkWidget *dialog = gtk_message_dialog_new( GTK_WINDOW(_window),
							GTK_DIALOG_DESTROY_WITH_PARENT,
							GTK_MESSAGE_ERROR,
							GTK_BUTTONS_CLOSE,
							"No particle database found." );
	    gtk_dialog_run( GTK_DIALOG(dialog) );
	    gtk_widget_destroy( dialog );
	    return;
	}

	_draw_action = 0;
	gtk_widget_queue_draw( _darea );
	
	// Done, start particle diagnostics dialog
	int crd;
	double val;
	Coordmapper cm = _frame.get_coordmapper( PLOT_AXIS_X1, PLOT_AXIS_Y1 );
	double c[2] = {(double)_start[0], (double)_start[1]};
	cm.inv_transform( c[0], c[1] );
	if( width > height ) {
	    crd = _geomplot.vb(1);
	    val = c[1];
	} else {
	    crd = _geomplot.vb(0);
	    val = c[0];
	}
	GTKParticleDiagDialog particlediag( _window, _plotter, crd, val );
	particlediag.run();
    }
}


void GTKGeomWindow::draw2( cairo_t *cairo )
{
    if( _tool == TOOL_PARTICLE_DIAG && _draw_action == 1 ) {
	int width = abs( _start[0] - _end[0] );
	int height = abs( _start[1] - _end[1] );
    	if( width > height ) {
	    cairo_move_to( cairo, _start[0], _start[1] );
	    cairo_line_to( cairo, _end[0], _start[1] );
	} else {
	    cairo_move_to( cairo, _start[0], _start[1] );
	    cairo_line_to( cairo, _start[0], _end[1] );
	}
	cairo_set_source_rgb( cairo, 0, 0, 0 );
	cairo_set_line_width( cairo, 1.0 );
	cairo_stroke( cairo );
    } else if( _tool == TOOL_FIELD_DIAG && _draw_action == 1 ) {
	cairo_move_to( cairo, _start[0], _start[1] );
	cairo_line_to( cairo, _end[0], _end[1] );
	cairo_set_source_rgb( cairo, 0, 0, 0 );
	cairo_set_line_width( cairo, 1.0 );
	cairo_stroke( cairo );
    }
}

    
void GTKGeomWindow::darea_motion2( GdkEventMotion *event )
{
    if( _tool == TOOL_PARTICLE_DIAG && (event->state & GDK_BUTTON1_MASK) ) {
	particle_diag( 1, event->x, event->y );
    } else if( _tool == TOOL_FIELD_DIAG && (event->state & GDK_BUTTON1_MASK) ) {
	field_diag( 1, event->x, event->y );
    }
}


void GTKGeomWindow::darea_button2( GdkEventButton *event )
{
    if( _tool == TOOL_PARTICLE_DIAG && event->type == GDK_BUTTON_PRESS && event->button == 1 ) {
	particle_diag( 0, event->x, event->y );
    } else if( _tool == TOOL_PARTICLE_DIAG && event->type == GDK_BUTTON_RELEASE && event->button == 1 ) {
	particle_diag( 2, event->x, event->y );
    } else if( _tool == TOOL_FIELD_DIAG && event->type == GDK_BUTTON_PRESS && event->button == 1 ) {
	field_diag( 0, event->x, event->y );
    } else if( _tool == TOOL_FIELD_DIAG && event->type == GDK_BUTTON_RELEASE && event->button == 1 ) {
	field_diag( 2, event->x, event->y );
    }
}


gboolean GTKGeomWindow::darea_draw_signal2( GtkWidget *widget, 
					    cairo_t *cairo,
					    gpointer object )
{
    GTKGeomWindow *plotter = (GTKGeomWindow *)object;
    plotter->draw2( cairo );
    return( FALSE );
}


gboolean GTKGeomWindow::darea_motion_signal2( GtkWidget *widget, 
					      GdkEventMotion *event,
					      gpointer object )
{
    GTKGeomWindow *plotter = (GTKGeomWindow *)object;
    plotter->darea_motion2( event );
    return( FALSE );
}


gboolean GTKGeomWindow::darea_button_signal2( GtkWidget *widget, 
					      GdkEventButton *event,
					      gpointer object )
{
    GTKGeomWindow *plotter = (GTKGeomWindow *)object;
    plotter->darea_button2( event );
    return( FALSE );
}


GTKGeomWindow::~GTKGeomWindow()
{
    if( _prefdata )
	delete _prefdata;
}


void GTKGeomWindow::update_view()
{
}


void GTKGeomWindow::zoom_fit( void )
{
    //std::cout << "Zoom fit\n";
    double min = -std::numeric_limits<double>::infinity();
    double max = std::numeric_limits<double>::infinity();
    _frame.set_ranges( PLOT_AXIS_X1, min, max );
    _frame.set_ranges( PLOT_AXIS_Y1, min, max );
    _frame.ruler_autorange_enable( PLOT_AXIS_X1, false, false );
    _frame.ruler_autorange_enable( PLOT_AXIS_Y1, false, false );

    draw_and_expose();
}


std::string GTKGeomWindow::track_text( double x, double y )
{
    double xc[3];
    int    i[3];
    std::wstringstream ss;
    std::locale mylocale("");
    ss.imbue( mylocale );

    xc[_geomplot.vb(0)] = x;
    xc[_geomplot.vb(1)] = y;
    xc[_geomplot.vb(2)] = _geomplot.get_level_si();
    
    i[_geomplot.vb(0)] = (int)floor( (x-_geom.origo(_geomplot.vb(0)))/_geom.h() );
    i[_geomplot.vb(1)] = (int)floor( (y-_geom.origo(_geomplot.vb(1)))/_geom.h() );
    i[_geomplot.vb(2)] = _geomplot.get_level();

    uint32_t smesh = _geom.mesh_check( i[0], i[1], i[2] );

    Vec3D loc(xc[0],xc[1],xc[2]);

    ss << "x = " << xc[0] << " m\n"
       << "y = " << xc[1] << " m\n"
       << "z = " << xc[2] << " m\n";
    ss << "i = " << i[0] << "\n"
       << "j = " << i[1] << "\n"
       << "k = " << i[2] << "\n";

    // Mesh id for debugging
    ss << "solid mesh = 0x" << std::hex << std::setfill(L'0') << std::setw(8)
       << smesh << std::dec << std::setfill(L' ') << ", ";
    if( (smesh & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEAR_SOLID ) {
	ss << "near solid: ";
	uint32_t ind = smesh & SMESH_NEAR_SOLID_INDEX_MASK;
	const uint8_t *ptr = _geom.nearsolid_ptr(ind);
	uint8_t sflag = ptr[0];
	uint8_t mask = 0x01;
	bool first = true;
	for( uint32_t a = 0; a < 6; a++ ) {
	    if( mask & sflag ) {
		if( !first )
		    ss << ", ";
		ss << "ndist[" << a << "] = " << (int)*ptr << " ";
		ptr++;
		first = false;
	    }
	    mask = mask << 1;
	}
	ss << "\n";
    } else if( (smesh & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_PURE_VACUUM )
	ss << "pure vacuum\n";
    else if( (smesh & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEUMANN )
	ss << "neumann\n";
    else if( (smesh & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET )
	ss << "dirichlet\n";

    if( _geom.have_solid_data() )
	ss << "solid = " << _geom.inside( loc ) << "\n";
    else
	ss << "solid data unavailable\n";

    if( _epot )
	ss << "epot = " << (*_epot)( loc ) << " V\n";
    if( _efield ) {
	Vec3D E = (*_efield)( loc );
	ss << "efield = " << E << " V/m\n";
	ss << "|efield| = " << E.norm2() << " V/m\n";
    }
    if( _bfield ) {
	Vec3D B = (*_bfield)( loc );
	ss << "bfield = " << B << " T\n";
	ss << "|bfield| = " << B.norm2() << " T\n";
    }
    if( _scharge )
	ss << "scharge = " << (*_scharge)( loc ) << " C/m3\n";
    if( _tdens )
	ss << "trajdens = " << (*_tdens)( loc ) << " A/m2\n";

    return( wstring_to_utf8( ss.str() ) );
}


void GTKGeomWindow::field_activate( void )
{
    field_type_e field_type = FIELD_NONE;
    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->field_none_radio) ) )
	field_type = FIELD_NONE;
    else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->field_J_radio) ) )
	field_type = FIELD_TRAJDENS;
    else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->field_rho_radio) ) )
	field_type = FIELD_SCHARGE;
    else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->field_phi_radio) ) )
	field_type = FIELD_EPOT;
    else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->field_E_radio) ) )
	field_type = FIELD_EFIELD;
    else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->field_Ex_radio) ) )
	field_type = FIELD_EFIELD_X;
    else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->field_Ey_radio) ) )
	field_type = FIELD_EFIELD_Y;
    else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->field_Ez_radio) ) )
	field_type = FIELD_EFIELD_Z;
    else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->field_B_radio) ) )
	field_type = FIELD_BFIELD;
    else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->field_Bx_radio) ) )
	field_type = FIELD_BFIELD_X;
    else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->field_By_radio) ) )
	field_type = FIELD_BFIELD_Y;
    else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->field_Bz_radio) ) )
	field_type = FIELD_BFIELD_Z;
    else
	throw( ErrorUnimplemented( ERROR_LOCATION ) );	
    _geomplot.set_fieldgraph_plot( field_type );

    double zmin, zmax;
    _geomplot.fieldgraph()->build_plot();
    _geomplot.fieldgraph()->get_zrange( zmin, zmax );
    //std::cout << "zmin = " << zmin << ", zmax = " << zmax << "\n";
    //std::string s = to_string( zmin );
    //gtk_entry_set_text( GTK_ENTRY(_prefdata->zmin_entry), s.c_str() );
    //s = to_string( zmax );
    //gtk_entry_set_text( GTK_ENTRY(_prefdata->zmax_entry), s.c_str() );
    char s[128];
    snprintf( s, 128, "%g", zmin );
    gtk_entry_set_text( GTK_ENTRY(_prefdata->zmin_entry), s );
    snprintf( s, 128, "%g", zmax );
    gtk_entry_set_text( GTK_ENTRY(_prefdata->zmax_entry), s );    
}


void GTKGeomWindow::field_toggled( GtkToggleButton *togglebutton,
				   gpointer         user_data )
{
    GTKGeomWindow *geomwindow = (GTKGeomWindow *)user_data;
    if( gtk_toggle_button_get_active( togglebutton ) )
	geomwindow->field_activate();
}


void *GTKGeomWindow::build_preferences( GtkWidget *notebook )
{
    if( _prefdata )
	delete _prefdata;
    _prefdata = new PreferencesData;

    // ****************************************************************************
    // Geometry notebook

    GtkWidget *grid = gtk_grid_new();
    gtk_grid_set_column_homogeneous( GTK_GRID(grid), TRUE );
    gtk_widget_set_hexpand( grid, TRUE );
    uint32_t yl = 0;

    // ****************************************************************************

    // Manual eqlines
    GtkWidget *label = gtk_label_new( "Manual eqlines" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    _prefdata->manual_eqlines_entry = gtk_entry_new();
    std::vector<double> eqlines_manual = _geomplot.get_eqlines_manual();
    std::string s;
    for( size_t a = 0; a < eqlines_manual.size(); a++ ) {
	char ss[128];
	snprintf( ss, 128, "%g", eqlines_manual[a] );
	s += ss;
	if( a != eqlines_manual.size()-1 )
	    s += " ";
    }
    gtk_entry_set_text( GTK_ENTRY(_prefdata->manual_eqlines_entry), s.c_str() );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->manual_eqlines_entry, 1, yl, 1, 1 );
    yl++;

    // Automatic eqlines
    label = gtk_label_new( "Automatic eqlines" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    size_t eqlines_auto = _geomplot.get_eqlines_auto();
    GtkAdjustment *automatic_eqlines_adj = gtk_adjustment_new( eqlines_auto, 0, 1000, 1, 10, 0 );
    _prefdata->automatic_eqlines_spin = gtk_spin_button_new( automatic_eqlines_adj, 1, 0 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->automatic_eqlines_spin, 1, yl, 1, 1 );
    yl++;

    // Particle division
    label = gtk_label_new( "Trajectory division" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    size_t particle_div = _geomplot.get_particle_div();
    GtkAdjustment *particle_div_adj = gtk_adjustment_new( particle_div, 0, 1000000, 1, 10, 0 );
    _prefdata->particle_div_spin = gtk_spin_button_new( particle_div_adj, 1, 0 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->particle_div_spin, 1, yl, 1, 1 );
    yl++;

    // Particle division, offset
    label = gtk_label_new( "Trajectory offset" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    size_t particle_offset = _geomplot.get_particle_offset();
    GtkAdjustment *particle_offset_adj = gtk_adjustment_new( particle_offset, 0, 1000000, 1, 10, 0 );
    _prefdata->particle_offset_spin = gtk_spin_button_new( particle_offset_adj, 1, 0 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->particle_offset_spin, 1, yl, 1, 1 );
    yl++;

    // QM discretation
    label = gtk_label_new( "Q/M discretation" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    _prefdata->qmdiscretation_check = gtk_check_button_new_with_label( "on/off" );
    bool qm_discretation = _geomplot.get_qm_discretation();
    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->qmdiscretation_check), qm_discretation );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->qmdiscretation_check, 1, yl, 1, 1 );
    yl++;

    // Mesh
    label = gtk_label_new( "Mesh" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    _prefdata->meshen_check = gtk_check_button_new_with_label( "on/off" );
    bool mesh = _geomplot.get_mesh();
    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->meshen_check), mesh );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 1, 1 );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->meshen_check, 1, yl, 1, 1 );
    yl++;

    // Add notebook page
    label = gtk_label_new( "Geometry" );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), grid, label );

    // ****************************************************************************
    // FieldGraph notebook, 4x4 grid

    grid = gtk_grid_new();
    gtk_grid_set_column_homogeneous( GTK_GRID(grid), TRUE );
    gtk_widget_set_hexpand( grid, TRUE );
    yl = 0;

    // ****************************************************************************

    // Field selection
    label = gtk_label_new( "Field selection" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 4, 1 );
    yl++;

    // First grid row
    _prefdata->field_none_radio = gtk_radio_button_new_with_label_from_widget( NULL, 
									       "None" );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->field_none_radio, 0, yl, 1, 1 );
    _prefdata->field_J_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->field_none_radio), 
									    "J" );
    if( !_tdens )
	gtk_widget_set_sensitive( _prefdata->field_J_radio, false );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->field_J_radio, 1, yl, 1, 1 );
    _prefdata->field_rho_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->field_none_radio), 
									      "rho" );
    if( !_scharge )
	gtk_widget_set_sensitive( _prefdata->field_rho_radio, false );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->field_rho_radio, 2, yl, 1, 1 );
    _prefdata->field_phi_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->field_none_radio), 
									      "phi" );
    if( !_epot )
	gtk_widget_set_sensitive( _prefdata->field_phi_radio, false );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->field_phi_radio, 3, yl, 1, 1 );
    yl++;

    // Second grid row
    _prefdata->field_E_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->field_none_radio), 
									    "|E|" );
    if( !_efield )
	gtk_widget_set_sensitive( _prefdata->field_E_radio, false );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->field_E_radio, 0, yl, 1, 1 );
    _prefdata->field_Ex_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->field_none_radio), 
									     "Ex" );
    if( !_efield )
	gtk_widget_set_sensitive( _prefdata->field_Ex_radio, false );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->field_Ex_radio, 1, yl, 1, 1 );
    _prefdata->field_Ey_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->field_none_radio), 
									     "Ey" );
    if( !_efield )
	gtk_widget_set_sensitive( _prefdata->field_Ey_radio, false );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->field_Ey_radio, 2, yl, 1, 1 );
    _prefdata->field_Ez_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->field_none_radio), 
									     "Ez" );
    if( !_efield )
	gtk_widget_set_sensitive( _prefdata->field_Ez_radio, false );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->field_Ez_radio, 3, yl, 1, 1 );
    yl++;

    // Third grid row
    _prefdata->field_B_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->field_none_radio), 
									    "|B|" );
    if( !_bfield )
	gtk_widget_set_sensitive( _prefdata->field_B_radio, false );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->field_B_radio, 0, yl, 1, 1 );
    _prefdata->field_Bx_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->field_none_radio), 
									     "Bx" );
    if( !_bfield )
	gtk_widget_set_sensitive( _prefdata->field_Bx_radio, false );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->field_Bx_radio, 1, yl, 1, 1 );
    _prefdata->field_By_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->field_none_radio), 
									     "By" );
    if( !_bfield )
	gtk_widget_set_sensitive( _prefdata->field_By_radio, false );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->field_By_radio, 2, yl, 1, 1 );
    _prefdata->field_Bz_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->field_none_radio), 
									     "Bz" );
    if( !_bfield )
	gtk_widget_set_sensitive( _prefdata->field_Bz_radio, false );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->field_Bz_radio, 3, yl, 1, 1 );
    yl++;

    field_type_e ftype = FIELD_NONE;
    if( _geomplot.fieldgraph() )
	ftype = _geomplot.fieldgraph()->field_type();
    switch( ftype ) {
    case FIELD_NONE:
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->field_none_radio), true );	
	break;
    case FIELD_TRAJDENS:
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->field_J_radio), true );	
	break;
    case FIELD_SCHARGE:
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->field_rho_radio), true );	
	break;
    case FIELD_EPOT:
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->field_phi_radio), true );	
	break;
    case FIELD_EFIELD:
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->field_E_radio), true );	
	break;
    case FIELD_EFIELD_X:
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->field_Ex_radio), true );	
	break;
    case FIELD_EFIELD_Y:
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->field_Ey_radio), true );	
	break;
    case FIELD_EFIELD_Z:
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->field_Ez_radio), true );	
	break;
    case FIELD_BFIELD:
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->field_B_radio), true );	
	break;
    case FIELD_BFIELD_X:
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->field_Bx_radio), true );	
	break;
    case FIELD_BFIELD_Y:
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->field_By_radio), true );	
	break;
    case FIELD_BFIELD_Z:
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->field_Bz_radio), true );	
	break;
    };
    g_signal_connect( _prefdata->field_none_radio, "toggled",
		      G_CALLBACK(field_toggled), (gpointer)this );
    g_signal_connect( _prefdata->field_J_radio, "toggled",
		      G_CALLBACK(field_toggled), (gpointer)this );
    g_signal_connect( _prefdata->field_rho_radio, "toggled",
		      G_CALLBACK(field_toggled), (gpointer)this );
    g_signal_connect( _prefdata->field_phi_radio, "toggled",
		      G_CALLBACK(field_toggled), (gpointer)this );
    g_signal_connect( _prefdata->field_E_radio, "toggled",
		      G_CALLBACK(field_toggled), (gpointer)this );
    g_signal_connect( _prefdata->field_Ex_radio, "toggled",
		      G_CALLBACK(field_toggled), (gpointer)this );
    g_signal_connect( _prefdata->field_Ey_radio, "toggled",
		      G_CALLBACK(field_toggled), (gpointer)this );
    g_signal_connect( _prefdata->field_Ez_radio, "toggled",
		      G_CALLBACK(field_toggled), (gpointer)this );
    g_signal_connect( _prefdata->field_B_radio, "toggled",
		      G_CALLBACK(field_toggled), (gpointer)this );
    g_signal_connect( _prefdata->field_Bx_radio, "toggled",
		      G_CALLBACK(field_toggled), (gpointer)this );
    g_signal_connect( _prefdata->field_By_radio, "toggled",
		      G_CALLBACK(field_toggled), (gpointer)this );
    g_signal_connect( _prefdata->field_Bz_radio, "toggled",
		      G_CALLBACK(field_toggled), (gpointer)this );

    // Zscaling
    zscale_e zscale = _geomplot.fieldgraph()->get_zscale();
    label = gtk_label_new( "Zscale" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 2, 1 );

    _prefdata->zscale_lin_radio = gtk_radio_button_new_with_label_from_widget( NULL, 
									   "Linear" );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->zscale_lin_radio, 2, yl, 2, 1 );
    yl++;
    _prefdata->zscale_log_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->zscale_lin_radio), 
									   "Logarithmic" );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->zscale_log_radio, 2, yl, 2, 1 );
    yl++;
    _prefdata->zscale_rellog_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->zscale_lin_radio), 
									      "Relative log" );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->zscale_rellog_radio, 2, yl, 2, 1 );
    yl++;

    if( zscale == ZSCALE_LINEAR )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->zscale_lin_radio), true );
    else if( zscale == ZSCALE_LOG )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->zscale_log_radio), true );
    else
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->zscale_rellog_radio), true );

    // Interpolation type
    interpolation_e interp = _geomplot.fieldgraph()->get_interpolation();    
    label = gtk_label_new( "Interpolation" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 2, 1 );

    _prefdata->int_closest_radio = gtk_radio_button_new_with_label_from_widget( NULL, 
									    "Closest" );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->int_closest_radio, 2, yl, 2, 1 );
    yl++;
    _prefdata->int_bilinear_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->int_closest_radio), 
									     "Bilinear" );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->int_bilinear_radio, 2, yl, 2, 1 );
    yl++;
    _prefdata->int_bicubic_radio = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(_prefdata->int_closest_radio), 
									    "Bicubic" );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->int_bicubic_radio, 2, yl, 2, 1 );
    yl++;

    if( interp == INTERPOLATION_CLOSEST )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->int_closest_radio), true );
    else if( interp == INTERPOLATION_BILINEAR )
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->int_bilinear_radio), true );
    else
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(_prefdata->int_bicubic_radio), true );

    // Ranges
    double zmin, zmax;
    _geomplot.fieldgraph()->get_zrange( zmin, zmax );

    // Range zmin
    label = gtk_label_new( "Range zmin" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 2, 1 );
    _prefdata->zmin_entry = gtk_entry_new();
    //s = to_string( zmin );
    //gtk_entry_set_text( GTK_ENTRY(_prefdata->zmin_entry), s.c_str() );
    char st[128];
    snprintf( st, 128, "%g", zmin );
    gtk_entry_set_text( GTK_ENTRY(_prefdata->zmin_entry), st );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->zmin_entry, 2, yl, 2, 1 );
    yl++;

    // Range zmax
    label = gtk_label_new( "Range zmax" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 2, 1 );
    _prefdata->zmax_entry = gtk_entry_new();
    //s = to_string( zmax );
    //gtk_entry_set_text( GTK_ENTRY(_prefdata->zmax_entry), s.c_str() );
    snprintf( st, 128, "%g", zmax );
    gtk_entry_set_text( GTK_ENTRY(_prefdata->zmax_entry), st );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->zmax_entry, 2, yl, 2, 1 );
    yl++;

    // Palette steps
    int steps = _geomplot.fieldgraph()->palette().get_steps();
    label = gtk_label_new( "Palette steps" );
    gtk_label_set_xalign( GTK_LABEL(label), 0 );
    gtk_label_set_yalign( GTK_LABEL(label), 0.5 );
    gtk_grid_attach( GTK_GRID(grid), label, 0, yl, 2, 1 );
    _prefdata->palette_steps_entry = gtk_entry_new();
    snprintf( st, 128, "%d", steps );
    gtk_entry_set_text( GTK_ENTRY(_prefdata->palette_steps_entry), st );
    gtk_grid_attach( GTK_GRID(grid), _prefdata->palette_steps_entry, 2, yl, 2, 1 );
    yl++;

    // ****************************************************************************

    // Add notebook page
    label = gtk_label_new( "FieldGraph" );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), grid, label );

    return( NULL );
}


void GTKGeomWindow::read_preferences( GtkWidget *notebook, void *pdata )
{
    // Eqlines
    size_t eqlines_auto = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(_prefdata->automatic_eqlines_spin) );
    _geomplot.set_eqlines_auto( eqlines_auto );
    std::vector<double> eqlines_manual;
    char *str = (char *)gtk_entry_get_text( GTK_ENTRY(_prefdata->manual_eqlines_entry) );
    while( *str != '\0' ) {
	if( !isdigit( *str ) && *str != '-' && *str != '.' && *str != '+' ) {
	    str++;
	    continue;
	}
	double val = strtod( str, &str );
	//std::cout << "val = " << val << "\n";
	eqlines_manual.push_back( val );
    }
    _geomplot.set_eqlines_manual( eqlines_manual );

    // Particle division and offset
    uint32_t particle_div = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(_prefdata->particle_div_spin) );
    uint32_t particle_offset = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(_prefdata->particle_offset_spin) );
    _geomplot.set_particle_div( particle_div, particle_offset );

    // QM discretation
    bool qm_discretation = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->qmdiscretation_check) );
    _geomplot.set_qm_discretation( qm_discretation );

    // Mesh
    bool mesh = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->meshen_check) );
    _geomplot.set_mesh( mesh );


    double zmin = atof( gtk_entry_get_text( GTK_ENTRY(_prefdata->zmin_entry) ) );
    double zmax = atof( gtk_entry_get_text( GTK_ENTRY(_prefdata->zmax_entry) ) );
    _geomplot.fieldgraph()->set_zrange( zmin, zmax );

    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->int_closest_radio) ) )
	_geomplot.fieldgraph()->set_interpolation( INTERPOLATION_CLOSEST );
    else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->int_bilinear_radio) ) )
	_geomplot.fieldgraph()->set_interpolation( INTERPOLATION_BILINEAR );
    else
	_geomplot.fieldgraph()->set_interpolation( INTERPOLATION_BICUBIC );

    if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->zscale_lin_radio) ) )
	_geomplot.fieldgraph()->set_zscale( ZSCALE_LINEAR );
    else if( gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(_prefdata->zscale_log_radio) ) )
	_geomplot.fieldgraph()->set_zscale( ZSCALE_LOG );
    else
	_geomplot.fieldgraph()->set_zscale( ZSCALE_RELLOG );

    int psteps = atoi( gtk_entry_get_text( GTK_ENTRY(_prefdata->palette_steps_entry) ) );
    _geomplot.fieldgraph()->palette().set_steps( psteps );
}


void GTKGeomWindow::combobox( GtkComboBox *combobox )
{
    //std::cout << "Combobox\n";

    int level;
    int levelmax;
    view_e view;
    if( _geom.geom_mode() == MODE_3D )
	view = (view_e)gtk_combo_box_get_active( GTK_COMBO_BOX(_combobox) );
    else {
	if( gtk_combo_box_get_active( GTK_COMBO_BOX(_combobox) ) == 0 )
	    view = VIEW_XY;
	else
	    view = VIEW_YX;
    }

    switch( view ) {
    case VIEW_XY:
	levelmax = _geom.size(2);
	level = levelmax/2;
        break;
    case VIEW_XZ:
	levelmax = _geom.size(1);
	level = levelmax/2;
        break;
    case VIEW_YX:
	levelmax = _geom.size(2);
	level = levelmax/2;
        break;
    case VIEW_YZ:
	levelmax = _geom.size(0);
	level = levelmax/2;
        break;
    case VIEW_ZX:
	levelmax = _geom.size(1);
	level = levelmax/2;
        break;
    case VIEW_ZY:
	levelmax = _geom.size(0);
	level = levelmax/2;
        break;
    default:
	throw( ErrorUnimplemented( ERROR_LOCATION ) );
	break;
    }

    // Set ranges to defaults on view change
    double min = -std::numeric_limits<double>::infinity();
    double max = std::numeric_limits<double>::infinity();
    _frame.set_ranges( PLOT_AXIS_X1, min, max );
    _frame.set_ranges( PLOT_AXIS_Y1, min, max );
    _frame.ruler_autorange_enable( PLOT_AXIS_X1, false, false );
    _frame.ruler_autorange_enable( PLOT_AXIS_Y1, false, false );

    // Update spinbutton
    /*
    g_signal_handlers_block_by_func( G_OBJECT(_spinbutton),
				     (void *)spinbutton_signal,
				     (gpointer)this );
    */
    g_signal_handlers_disconnect_by_func( G_OBJECT(_spinbutton),
					  (void *)spinbutton_signal,
					  (gpointer)this );
    gtk_spin_button_set_range( GTK_SPIN_BUTTON(_spinbutton), 0, levelmax-1 );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(_spinbutton), level );
    /*
    g_signal_handlers_unblock_by_func( G_OBJECT(_spinbutton),
				       (void *)spinbutton_signal,
				       (gpointer)this );
    */
    //g_signal_stop_emission_by_name( G_OBJECT(_spinbutton), "value-changed" );
    g_signal_connect( G_OBJECT(_spinbutton), "value-changed",
                      G_CALLBACK(spinbutton_signal),
                      (gpointer)this );


    // Update view
    _geomplot.set_view( view, level );
    update_view();
    draw_and_expose();
}


void GTKGeomWindow::spinbutton( GtkSpinButton *spinbutton )
{
    //std::cout << "Spinbutton\n";
    //g_signal_stop_emission_by_name( G_OBJECT(_spinbutton), "value-changed" );

    int level = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(_spinbutton) );

    std::wstringstream ss;
    std::locale mylocale("");
    ss.imbue( mylocale );
    switch( _geomplot.get_view() ) {
    case VIEW_XY:
	ss << "z = "<< _geom.origo(2)+level*_geom.h() << " m";
	break;
    case VIEW_XZ:
	ss << "y = "<< _geom.origo(1)+level*_geom.h() << " m";
	break;
    case VIEW_YZ:
	ss << "x = "<< _geom.origo(0)+level*_geom.h() << " m";
	break;
    case VIEW_ZX:
	ss << "y = "<< _geom.origo(1)+level*_geom.h() << " m";
	break;
    case VIEW_ZY:
	ss << "x = "<< _geom.origo(0)+level*_geom.h() << " m";
	break;
    case VIEW_YX:
	ss << "z = "<< _geom.origo(2)+level*_geom.h() << " m";
	break;
    default:
	throw( ErrorUnimplemented( ERROR_LOCATION ) );
	break;
    }
    gtk_statusbar_pop( GTK_STATUSBAR(_statusbar), 0 );
    std::string s8 = wstring_to_utf8( ss.str() );
    gtk_statusbar_push( GTK_STATUSBAR(_statusbar), 0, s8.c_str() );

    // Update view
    _geomplot.set_view( _geomplot.get_view(), level );
    update_view();
    draw_and_expose();
}


void GTKGeomWindow::menuitem_tool_change( GtkToolButton *button )
{
    //std::cout << "GTKGeomWindow: ";
    
    int tool;
    const char *label = gtk_tool_button_get_label( button );
    if( !strcmp( label, "Particle diagnostics" ) ) {
	//std::cout << "TOOL_PARTICLE_DIAG ";
	tool = TOOL_PARTICLE_DIAG;
    } else if( !strcmp( label, "Field diagnostics" ) ) {
	//std::cout << "TOOL_FIELD_DIAG ";
	tool = TOOL_FIELD_DIAG;
    } else {
	//std::cout << "TOOL_UNKNOWN\n";
	tool = TOOL_UNKNOWN;
	return;
    }

    if( !gtk_toggle_tool_button_get_active( GTK_TOGGLE_TOOL_BUTTON(button) ) ) {
	// Disable tool
	//std::cout << "disable\n";
	_tool = TOOL_UNKNOWN;
    } else {
	// Enable tool
	//std::cout << "enable\n";
	_tool = tool;
    }
}


void GTKGeomWindow::geom3d_launch( void )
{
    _plotter.new_geometry_3d_plot_window();
}


void GTKGeomWindow::combobox_signal( GtkComboBox *combobox,
				     gpointer object )
{
    GTKGeomWindow *window = (GTKGeomWindow *)object;
    window->combobox( combobox );
}


void GTKGeomWindow::spinbutton_signal( GtkSpinButton *spinbutton,
				       gpointer object )
{
    GTKGeomWindow *window = (GTKGeomWindow *)object;
    window->spinbutton( spinbutton );
}


void GTKGeomWindow::menuitem_tool_change_signal( GtkToolButton *button,
						 gpointer object )
{
    GTKGeomWindow *window = (GTKGeomWindow *)object;
    window->menuitem_tool_change( button );
}


void GTKGeomWindow::menuitem_geom3d_signal( GtkToolButton *button,
					    gpointer object )
{
    GTKGeomWindow *window = (GTKGeomWindow *)object;
    window->geom3d_launch();
}


