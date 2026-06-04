/*! \file gtkgeomwindow.hpp
 *  \brief %Geometry view window
 */

/* Copyright (c) 2005-2013, 2022 Taneli Kalvas. All rights reserved.
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

#ifndef GTKGEOMWINDOW_HPP
#define GTKGEOMWINDOW_HPP 1


#include <vector>

#include "gtkframewindow.hpp"
#include "geomplot.hpp"
#include "meshscalarfield.hpp"
#include "epot_field.hpp"
#include "epot_efield.hpp"
#include "vectorfield.hpp"



/*! \brief Interactive geometry plotter window.
 */
class GTKGeomWindow : public GTKFrameWindow {

    struct PreferencesData {
	GtkWidget *manual_eqlines_entry;
	GtkWidget *automatic_eqlines_spin;
	GtkWidget *particle_div_spin;
	GtkWidget *particle_offset_spin;

	GtkWidget *qmdiscretation_check;
	GtkWidget *meshen_check;

	GtkWidget *field_none_radio;
	GtkWidget *field_J_radio;
	GtkWidget *field_rho_radio;
	GtkWidget *field_phi_radio;
	
	GtkWidget *field_E_radio;
	GtkWidget *field_Ex_radio;
	GtkWidget *field_Ey_radio;
	GtkWidget *field_Ez_radio;
	
	GtkWidget *field_B_radio;
	GtkWidget *field_Bx_radio;
	GtkWidget *field_By_radio;
	GtkWidget *field_Bz_radio;
	
	GtkWidget *int_closest_radio;
	GtkWidget *int_bilinear_radio;
	GtkWidget *int_bicubic_radio;
	
	GtkWidget *zscale_lin_radio;
	GtkWidget *zscale_log_radio;
	GtkWidget *zscale_rellog_radio;
	
	GtkWidget *zmin_entry;
	GtkWidget *zmax_entry;
	
	GtkWidget *palette_steps_entry;
    };

    GeomPlot                 _geomplot;

    const Geometry          &_geom;
    const EpotField         *_epot;
    const EpotEfield        *_efield;
    const MeshScalarField   *_scharge;
    const MeshScalarField   *_tdens;
    const VectorField       *_bfield;
    const ParticleDataBase  *_pdb;

    GtkWidget               *_spinbutton;
    GtkWidget               *_combobox;

    PreferencesData         *_prefdata;

    void update_view();

    virtual void zoom_fit( void );
    virtual std::string track_text( double x, double y );
 
    virtual void *build_preferences( GtkWidget *notebook );
    virtual void read_preferences( GtkWidget *notebook, void *pdata );

    void combobox( GtkComboBox *combobox );
    void spinbutton( GtkSpinButton *spinbutton );
    void menuitem_tool_change( GtkToolButton *button );
    void field_diag( int action, double x, double y );
    void particle_diag( int action, double x, double y );
    void draw2( cairo_t *cairo );
    void darea_motion2( GdkEventMotion *event );
    void darea_button2( GdkEventButton *event );
    void field_activate( void );
    void geom3d_launch( void );

    static void combobox_signal( GtkComboBox *combobox,
				 gpointer object );
    static void spinbutton_signal( GtkSpinButton *spinbutton,
				   gpointer object );
    static void menuitem_tool_change_signal( GtkToolButton *button,
					     gpointer object );
    static void menuitem_geom3d_signal( GtkToolButton *button,
					gpointer object );
    static gboolean darea_draw_signal2( GtkWidget *widget, 
					cairo_t *cairo,
					gpointer object );
    static gboolean darea_motion_signal2( GtkWidget *widget, 
					  GdkEventMotion *event,
					  gpointer object );
    static gboolean darea_button_signal2( GtkWidget *widget, 
					  GdkEventButton *event,
					  gpointer object );
    static void field_toggled( GtkToggleButton *togglebutton,
			       gpointer         user_data );

public:

    /*! \brief Constructor.
     */
    GTKGeomWindow( class GTKPlotter &plotter,
		   const Geometry &geom,
		   const EpotField *epot,
		   const EpotEfield *efield,
		   const MeshScalarField *scharge,
		   const MeshScalarField *tdens,
		   const VectorField *bfield,
		   const ParticleDataBase *pdb );
    
    /*! \brief Destructor.
     */
    virtual ~GTKGeomWindow();
};


#endif
