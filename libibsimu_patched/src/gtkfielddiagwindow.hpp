/*! \file gtkfielddiagwindow.hpp
 *  \brief Field diagnostic window
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

#ifndef GTKFIELDDIAGWINDOW_HPP
#define GTKFIELDDIAGWINDOW_HPP 1


#include "gtkframewindow.hpp"
#include "fielddiagplot.hpp"


/*! \brief Interactive field diagnostic plotter
 */
class GTKFieldDiagWindow : public GTKFrameWindow {

    const Geometry            &_geom;
    FieldDiagPlot              _plot;

    double                     _x1min;
    double                     _x1max;
    double                     _x2min;
    double                     _x2max;

    virtual void zoom_fit( void );
    virtual void *build_preferences( GtkWidget *notebook );
    virtual void read_preferences( GtkWidget *notebook, void *pdata );

    virtual std::string track_text( double x, double y );

    void export_data( void );
    static void menuitem_export_signal( GtkToolButton *button,
					gpointer object );

public:

    /*! \brief Make new field diagnostic window.
     *
     *  Field is diagnosed on a line from \a x1 to \a x2. On the first
     *  y-axis, \a g1type field is plotted. On the second \a g2type
     *  field is plotted. The data displayed on first x-axis is
     *  specified by \a dist1type and the second x-axis on \a
     *  dist2type.
     */
    GTKFieldDiagWindow( GTKPlotter &plotter, const Geometry &geom, size_t N, 
			const Vec3D &x1, const Vec3D &x2, 
			const field_diag_type_e diag[2], const field_loc_type_e loc[2] );

    /*! \brief Destructor.
     */
    virtual ~GTKFieldDiagWindow();

};


#endif
