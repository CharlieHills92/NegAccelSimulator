/*! \file gtkframewindow.hpp
 *  \brief Window for GTK plots with frames
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

#ifndef GTKFRAMEWINDOW_HPP
#define GTKFRAMEWINDOW_HPP 1


#include <cairo.h>
#include <gtk/gtk.h>
#include "gtkwindow.hpp"
#include "gtkplotter.hpp"
#include "frame.hpp"
#include "graph3d.hpp"



/*! \brief Base class for interactive plotters.
 *
 *  Provides all the common features for the interactive gtk based
 *  plotter windows.
 */
class GTKFrameWindow : public GTKWindow {

    size_t            _width;
    size_t            _height;

    GtkWidget        *_trackwindow;
    GtkWidget        *_tracklabel;

    cairo_t          *_cairo;
    cairo_surface_t  *_surface;

    void hardcopy( void );
    void move( int action, double x, double y );
    void track( int action, double x, double y );
    void zoom_out( double x, double y );
    void zoom_in( double x, double y );
    void zoom_window( int action, double x, double y );

    void darea_motion( GdkEventMotion *event );
    void darea_enter( GdkEventCrossing *event );
    void darea_leave( GdkEventCrossing *event );
    void darea_button( GdkEventButton *event );

    void configure( void );
    void delete_window( void );

    void menuitem_tool_change( GtkToolButton *button );
    void menuitem_preferences( GtkMenuItem *menuitem );

    void frame_draw( void );

    static gboolean window_delete_signal( GtkWidget *widget, 
					  GdkEventExpose *event, 
					  gpointer object );

    static void menuitem_quit_signal( GtkMenuItem *menuitem,
				      gpointer object );
    static void menuitem_tool_change_signal( GtkToolButton *button,
					     gpointer object );
    static void menuitem_hardcopy_signal( GtkToolButton *button,
					  gpointer object );
    static void menuitem_zoom_fit_signal( GtkToolButton *button,
					  gpointer object );
    static void menuitem_preferences_signal( GtkMenuItem *menuitem,
					     gpointer object );

    static gboolean darea_configure_signal( GtkWidget *widget, 
					    GdkEventConfigure *event, 
					    gpointer object );
    static gboolean darea_draw_signal( GtkWidget *widget, 
				       cairo_t *cairo,
				       gpointer object );
    static gboolean darea_button_signal( GtkWidget *widget, 
					 GdkEventButton *event,
					 gpointer object );
    static gboolean darea_motion_signal( GtkWidget *widget, 
					 GdkEventMotion *event,
					 gpointer object );
    static gboolean darea_enter_signal( GtkWidget *widget, 
					GdkEventCrossing *event,
					gpointer object );
    static gboolean darea_leave_signal( GtkWidget *widget, 
					GdkEventCrossing *event,
					gpointer object );

protected:

    GTKPlotter       &_plotter;

    Frame             _frame;
    GtkWidget        *_window;
    GtkWidget        *_darea;
    GtkWidget        *_menubar;
    GtkWidget        *_menu_file;
    GtkWidget        *_toolbar;
    GtkWidget        *_statusbar;

    GtkToolItem      *_radioitem;

    double            _track_px; // Track x coordinate in pixels
    double            _track_py; // Track y coordinate in pixels

    int               _start[2];     /*!< \brief Tool specific coordinates. */
    int               _end[2];       /*!< \brief Tool specific coordinates. */
    int               _tool;         /*!< \brief Selected tool. */
    int               _draw_action;  /*!< \brief Flag for toll specific draw actions. */

    GTKFrameWindow( GTKPlotter &plotter );

    virtual void zoom_fit( void );

    /*! \brief Function which providies tracking text at (x,y)
     *
     *  Virtual function. Provides empty string by default but can be
     *  overwritten by implementations of GTKWindow.
     */
    virtual std::string track_text( double x, double y );

    void draw( cairo_t *cairo );
    //void expose( int x, int y, int width, int height );

public:

    /*! \brief Destructor.
     */
    virtual ~GTKFrameWindow();

    virtual void *build_preferences( GtkWidget *notebook );
    virtual void read_preferences( GtkWidget *notebook, void *pdata );

    void draw_and_expose( void );
    void show( void );
};


#endif
