/*! \file gtkparticlediagexportdialog.cpp
 *  \brief Dialog for exporting field diagnostic data.
 */

/* Copyright (c) 2005-2011,2013,2022 Taneli Kalvas. All rights reserved.
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


#include "config.h"
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "gtkparticlediagexportdialog.hpp"


GTKParticleDiagExportDialog::GTKParticleDiagExportDialog( GtkWidget *window, ParticleDiagPlot *plot )
    : _window(window), _plot(plot)
{

}


GTKParticleDiagExportDialog::~GTKParticleDiagExportDialog()
{

}


void GTKParticleDiagExportDialog::run( void )
{
   GtkWidget *dialog = gtk_file_chooser_dialog_new( "Export particle data",
						    GTK_WINDOW(_window),
						    GTK_FILE_CHOOSER_ACTION_SAVE,
						    dgettext("gtk30","_Cancel"), GTK_RESPONSE_CANCEL,
						    dgettext("gtk30","_Save"), GTK_RESPONSE_ACCEPT,
						    NULL );

#ifdef HAVE_GETCWD
    size_t size = 1024;
    char *buf;
    while( 1 ) {
        buf = new char[size];
        char *ret = getcwd( buf, size );
        if( ret )
            break;
        delete [] buf;
    }
    GFile *gfile = g_file_new_for_path( buf );
    gtk_file_chooser_set_current_folder_file( GTK_FILE_CHOOSER(dialog), gfile, NULL );
    g_object_unref( gfile );
#endif
   gtk_file_chooser_set_current_name( GTK_FILE_CHOOSER(dialog), "particles.txt" );
   gtk_file_chooser_set_show_hidden( GTK_FILE_CHOOSER(dialog), TRUE );
   gtk_file_chooser_set_local_only(GTK_FILE_CHOOSER(dialog), TRUE );
   
   gtk_widget_show_all( dialog );
   if( gtk_dialog_run( GTK_DIALOG(dialog) ) == GTK_RESPONSE_ACCEPT ) {

       char *filename = gtk_file_chooser_get_filename( GTK_FILE_CHOOSER(dialog) );

       // Write output to filename
       _plot->export_data( filename );

       g_free( filename );
   }

   gtk_widget_destroy( dialog );
}






