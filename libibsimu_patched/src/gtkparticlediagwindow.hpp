/*! \file gtkparticlediagwindow.hpp
 *  \brief %Particle diagnostic window.
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

#ifndef GTKPARTICLEDIAGWINDOW_HPP
#define GTKPARTICLEDIAGWINDOW_HPP 1


#include "gtkframewindow.hpp"
#include "particledatabase.hpp"
#include "particlediagplot.hpp"
#include "types.hpp"


/*! \brief Interactive particle diagnostics plotter.
 */
class GTKParticleDiagWindow : public GTKFrameWindow {

    ParticleDiagPlot _plot;

    struct PreferencesData {
	GtkWidget *plot_scatter_radio;
	GtkWidget *plot_colormap_radio;
	GtkWidget *histo_n_spin;
	GtkWidget *histo_m_spin;
	GtkWidget *histo_acc_closest_radio;
	GtkWidget *histo_acc_bilinear_radio;
	GtkWidget *int_closest_radio;
	GtkWidget *int_bilinear_radio;
	GtkWidget *int_bicubic_radio;
	GtkWidget *style_histo_radio;
	GtkWidget *style_line_radio;
	GtkWidget *dot_size_spin;
	GtkWidget *ellipse_check;
    };
    PreferencesData *_prefdata;

    virtual std::string track_text( double x, double y );

    void plot_2d_type_toggled_process( void );

    void build_preferences_1d( GtkWidget *notebook );
    void build_preferences_2d( GtkWidget *notebook );
    virtual void *build_preferences( GtkWidget *notebook );
    virtual void read_preferences( GtkWidget *notebook, void *pdata );

    void export_data( void );
    static void plot_2d_type_toggled( GtkToggleButton *togglebutton,
				      gpointer         user_data );
    static void menuitem_export_signal( GtkToolButton *button,
					gpointer object );

public:

    /*! \brief Constructor for diagnostics window.
     *
     * \a style is the style of plot with 0 being scatter plot and 1
     * being colormap (histogram) plot.
     */
    GTKParticleDiagWindow( GTKPlotter &plotter, const ParticleDataBase &pdb, 
			   const Geometry &geom,
			   coordinate_axis_e axis, double level, 
			   particle_diag_plot_type_e type,
			   trajectory_diagnostic_e diagx, 
			   trajectory_diagnostic_e diagy );

    /*! \brief Destructor.
     */
    virtual ~GTKParticleDiagWindow();

};


#endif
