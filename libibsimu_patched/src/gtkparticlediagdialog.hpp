/*! \file gtkparticlediagdialog.hpp
 *  \brief Dialog for constructing particle diagnostic windows
 */

/* Copyright (c) 2005-2012 Taneli Kalvas. All rights reserved.
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

#ifndef GTKPARTICLEDIAGDIALOG_HPP
#define GTKPARTICLEDIAGDIALOG_HPP 1


#include <gtk/gtk.h>
#include "gtkplotter.hpp"


/*! \brief Dialog window for starting interactive particle diagnostics.
 */
class GTKParticleDiagDialog
{

    GtkWidget      *_window;
    GTKPlotter     *_plotter;
    const Geometry *_geom;

    int             _plane;
    double          _val;

    GtkWidget      *_radio_plane_x;
    GtkWidget      *_radio_plane_y;
    GtkWidget      *_radio_plane_z;

    GtkWidget      *_radio_emit_xx;
    GtkWidget      *_radio_emit_yy;
    GtkWidget      *_radio_emit_ra; // Only for cyl
    GtkWidget      *_radio_emit_zz; // Or converted emittance in case of cyl

    GtkWidget      *_radio_prof_yz;
    GtkWidget      *_radio_prof_xz;
    GtkWidget      *_radio_prof_xy;

    GtkWidget      *_radio_plot_scatter;
    GtkWidget      *_radio_plot_colormap;

    GtkWidget      *_radio_prof_x;
    GtkWidget      *_radio_prof_y;
    GtkWidget      *_radio_prof_z;

    GtkWidget      *_radio_prof_xp;
    GtkWidget      *_radio_prof_yp;
    GtkWidget      *_radio_prof_zp;

    GtkWidget      *_radio_energy;
    GtkWidget      *_radio_qm;
    GtkWidget      *_radio_charge;
    GtkWidget      *_radio_mass;

    void plot1d_toggled2( GtkToggleButton *togglebutton );
    static void plot1d_toggled( GtkToggleButton *togglebutton,
				gpointer         user_data );

    void conversion_toggled2( GtkToggleButton *togglebutton );
    static void conversion_toggled( GtkToggleButton *togglebutton,
				    gpointer         user_data );

    void plane_activated( void );
    static void plane_toggled( GtkToggleButton *togglebutton,
			       gpointer         user_data );


public:

    /*! \brief Dialog window for starting interactive particle diagnostics.
     *
     *  The \a window is the parent of the dialog, the \a plotter is
     *  the main data container and \a plane and \a val define the
     *  plane of diagnostics.
     */
    GTKParticleDiagDialog( GtkWidget *window, GTKPlotter &plotter, int plane, double val );

    /*! \brief Destructor.
     */
    ~GTKParticleDiagDialog();

    /*! \brief Run the dialog.
     */
    void run( void );
};


#endif
