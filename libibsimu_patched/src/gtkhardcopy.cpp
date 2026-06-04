/*! \file gtkhardcopy.cpp
 *  \brief Dialog window for producing hard copies
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


#include "config.h"
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <cmath>
#include <cairo-svg.h>
#include <cairo-ps.h>
#include <cairo-pdf.h>
#include "gtkhardcopy.hpp"
#include "softwarerenderer.hpp"



GTKHardcopy::GTKHardcopy( GtkWidget *window, Frame *frame, size_t width, size_t height )
    : _window(window), _frame(frame), _geom3dplot(NULL), _width(width), _height(height)
{
    _aspect = (double)_height/_width;
}


GTKHardcopy::GTKHardcopy( GtkWidget *window, Geom3DPlot *geom3dplot, size_t width, size_t height )
    : _window(window), _frame(NULL), _geom3dplot(geom3dplot), _width(width), _height(height)
{
    _aspect = (double)_height/_width;
}


GTKHardcopy::~GTKHardcopy()
{

}


int GTKHardcopy::type_from_extension( const char *filename )
{
    int b;
    int len = strlen( filename );
    
    /* Search the beginning of extension (where '.' should be) */
    b = len-1;
    while( b > 0 && *(filename+b) != '.' )
	b--;
    if( b == 0 )
	return( 0 );
    
    if( !strcmp( filename+b, ".png" ) )
	return( 1 );
    else if( !strcmp( filename+b, ".svg" ) )
	return( 2 );
    else if( !strcmp( filename+b, ".eps" ) || !strcmp( filename+b, ".ps" ) )
	return( 3 );
    else if( !strcmp( filename+b, ".pdf" ) )
	return( 4 );

    return( 0 );
}


// Strip path and old extension and add the selected extension.
void GTKHardcopy::ensure_extension( std::string &filename, const std::string &extension )
{
    int a, b;
    int len = filename.length();

    /* Search the beginning of filename excluding path and last '/' */
    a = len-1;
    while( a > 0 && filename[a] != '/' )
	a--;
    if( filename[a] == '/' )
	a++;

    /* Search the beginning of extension (where '.' should be) */
    b = len-1;
    while( b > 0 && filename[b] != '.' )
	b--;
    if( b <= a )
	b = len;

    // Copy filename body to output
    filename = filename.substr( a, b-a );
    // Append extension
    filename += extension;
}


void GTKHardcopy::treeview_changed_signal( GtkTreeSelection *selection,
					   gpointer userdata )
{
    GTKHardcopy *hardcopy = (GTKHardcopy *)userdata;
    hardcopy->treeview_changed( selection );
}


void GTKHardcopy::treeview_changed( GtkTreeSelection *selection )
{
    char *filename_ptr = gtk_file_chooser_get_filename( GTK_FILE_CHOOSER(_dialog) );
    std::string filename( filename_ptr );
    GtkTreeModel *model;
    GtkTreeIter iter;
    int active;

    gtk_tree_selection_get_selected( selection, &model, &iter );
    gtk_tree_model_get( model, &iter, 2, &active, -1 );

    switch( active ) {
    case 0:
	gtk_expander_set_label( GTK_EXPANDER(_expander), "Select File _Type (By Extension)" );
	return;
	break;
    case 1:
	ensure_extension( filename, ".png" );
	gtk_expander_set_label( GTK_EXPANDER(_expander), "Select File _Type (PNG Image)" );
	break;
    case 2:
	ensure_extension( filename, ".svg" );
	gtk_expander_set_label( GTK_EXPANDER(_expander), "Select File _Type (SVG Image)" );
	break;
    case 3:
	ensure_extension( filename, ".eps" );
	gtk_expander_set_label( GTK_EXPANDER(_expander), "Select File _Type (Encapsulated Postscript)" );
	break;
#ifdef CAIRO_HAS_PDF_SURFACE
    case 4:
	ensure_extension( filename, ".pdf" );
	gtk_expander_set_label( GTK_EXPANDER(_expander), "Select File _Type (Portable Document Format)" );
	break;
#endif
    default:
	throw( Error( ERROR_LOCATION, "unsupported file type" ) );
    }

    gtk_file_chooser_set_current_name( GTK_FILE_CHOOSER(_dialog), filename.c_str() );
    g_free( filename_ptr );
}


void GTKHardcopy::spinx( void )
{
    int val = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(_spinx) );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(_spiny), floor(_aspect*val+0.5) );
}


void GTKHardcopy::spiny( void )
{
    int val = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(_spiny) );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(_spinx), floor(val/_aspect+0.5) );
}


void GTKHardcopy::spinx_signal( GtkSpinButton *spinbutton,
				gpointer object )
{
    GTKHardcopy *plotter = (GTKHardcopy *)object;
    plotter->spinx();
}


void GTKHardcopy::spiny_signal( GtkSpinButton *spinbutton,
				gpointer object )
{
    GTKHardcopy *plotter = (GTKHardcopy *)object;
    plotter->spiny();
}


void GTKHardcopy::run( void )
{
    _dialog = gtk_file_chooser_dialog_new( "Make a hardcopy",
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
    gtk_file_chooser_set_current_folder_file( GTK_FILE_CHOOSER(_dialog), gfile, NULL );
    g_object_unref( gfile );
#endif
    gtk_file_chooser_set_current_name( GTK_FILE_CHOOSER(_dialog), "hardcopy" );
    gtk_file_chooser_set_show_hidden( GTK_FILE_CHOOSER(_dialog), TRUE );
    gtk_file_chooser_set_local_only(GTK_FILE_CHOOSER(_dialog), TRUE );

    // Set filters
    GtkFileFilter *filter = gtk_file_filter_new();
    gtk_file_filter_set_name( filter, "All files" );
    gtk_file_filter_add_pattern( filter, "*" );
    gtk_file_chooser_add_filter( GTK_FILE_CHOOSER(_dialog), filter );
    filter = gtk_file_filter_new();
    gtk_file_filter_set_name( filter, "PNG image (*.png)" );
    gtk_file_filter_add_pattern( filter, "*.png" );
    gtk_file_chooser_add_filter( GTK_FILE_CHOOSER(_dialog), filter );
    filter = gtk_file_filter_new();
    gtk_file_filter_set_name( filter, "SVG image (*.svg)" );
    gtk_file_filter_add_pattern( filter, "*.svg" );
    gtk_file_chooser_add_filter( GTK_FILE_CHOOSER(_dialog), filter );
    filter = gtk_file_filter_new();
    gtk_file_filter_set_name( filter, "Postscript (*.ps,*.eps)" );
    gtk_file_filter_add_pattern( filter, "*.ps" );
    gtk_file_filter_add_pattern( filter, "*.eps" );
    gtk_file_chooser_add_filter( GTK_FILE_CHOOSER(_dialog), filter );
#ifdef CAIRO_HAS_PDF_SURFACE
    filter = gtk_file_filter_new();
    gtk_file_filter_set_name( filter, "PDF (*.pdf)" );
    gtk_file_filter_add_pattern( filter, "*.pdf" );
    gtk_file_chooser_add_filter( GTK_FILE_CHOOSER(_dialog), filter );
#endif 

    // File type
    GtkListStore *list_store = gtk_list_store_new( 3, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_INT );
    GtkTreeIter iter;
    gtk_list_store_append( list_store, &iter );
    gtk_list_store_set( list_store, &iter, 
			0, "By Extension", 
			1, "",    
			2, 0, -1);
    gtk_list_store_append( list_store, &iter );
    gtk_list_store_set( list_store, &iter, 
			0, "PNG Image",    
			1, "png", 
			2, 1, -1);
    if( _frame ) {
	gtk_list_store_append( list_store, &iter );
	gtk_list_store_set( list_store, &iter, 
			    0, "SVG Image",    
			    1, "svg", 
			    2, 2, -1);
	gtk_list_store_append( list_store, &iter );
	gtk_list_store_set( list_store, &iter, 
			    0, "Encapsulated Postscript", 
			    1, "eps",
			    2, 3, -1);
#ifdef CAIRO_HAS_PDF_SURFACE
	gtk_list_store_append( list_store, &iter );
	gtk_list_store_set( list_store, &iter, 
			    0, "Portable Document Format", 
			    1, "pdf",
			    2, 4, -1);
#endif
    }
    GtkWidget *treeview = gtk_tree_view_new_with_model( GTK_TREE_MODEL(list_store) );  
    GtkCellRenderer *renderer;
    GtkTreeViewColumn *column;
    renderer = gtk_cell_renderer_text_new();
    column = gtk_tree_view_column_new_with_attributes( "File Type", renderer,
						       "text", 0, NULL );
    gtk_tree_view_append_column( GTK_TREE_VIEW(treeview), column );
    renderer = gtk_cell_renderer_text_new();
    column = gtk_tree_view_column_new_with_attributes( "Extension", renderer,
						       "text", 1, NULL );
    gtk_tree_view_append_column( GTK_TREE_VIEW(treeview), column );
    GtkTreeSelection *selection = gtk_tree_view_get_selection( GTK_TREE_VIEW(treeview) );
    gtk_tree_selection_set_mode( selection, GTK_SELECTION_SINGLE );
    gtk_tree_model_get_iter_first( GTK_TREE_MODEL(list_store), &iter );
    gtk_tree_selection_select_iter( selection, &iter );
    g_signal_connect( G_OBJECT(selection), "changed",
		      G_CALLBACK(treeview_changed_signal),
		      (gpointer)this );

    // Expander for treeview object (file type selector)
    _expander = gtk_expander_new_with_mnemonic( "Select File _Type (By Extension)" );
    gtk_expander_set_use_underline( GTK_EXPANDER(_expander), TRUE );
    gtk_expander_set_expanded( GTK_EXPANDER(_expander), FALSE );
    gtk_container_add( GTK_CONTAINER(_expander), treeview );

    // Spinbutton for resolution
    GtkWidget *hbox = gtk_box_new( GTK_ORIENTATION_HORIZONTAL, 2 );
    GtkWidget *resolabel = gtk_label_new( "Image resolution:" );
    _spinx = gtk_spin_button_new_with_range( 1, 10000, 1 );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(_spinx), _width );
    gtk_widget_set_size_request( _spinx, 200, -1 );
    g_signal_connect( G_OBJECT(_spinx), "value-changed",
                      G_CALLBACK(spinx_signal),
                      (gpointer)this );
    _spiny = gtk_spin_button_new_with_range( 1, 10000, 1 );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(_spiny), _height );
    gtk_widget_set_size_request( _spiny, 200, -1 );
    g_signal_connect( G_OBJECT(_spiny), "value-changed",
                      G_CALLBACK(spiny_signal),
                      (gpointer)this );
    gtk_box_pack_start( GTK_BOX(hbox), resolabel, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), _spinx, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(hbox), _spiny, FALSE, TRUE, 0 );

    // Pack Extra widgets for file chooser in a vbox
    GtkWidget *vbox = gtk_box_new( GTK_ORIENTATION_VERTICAL, 2 );
    gtk_box_pack_start( GTK_BOX(vbox), _expander, FALSE, TRUE, 0 );
    gtk_box_pack_start( GTK_BOX(vbox), hbox, FALSE, TRUE, 0 );

    // Set extra widgets to file chooser
    gtk_file_chooser_set_extra_widget( GTK_FILE_CHOOSER(_dialog), vbox );

    gtk_widget_show_all( _dialog );
    if( gtk_dialog_run( GTK_DIALOG(_dialog) ) == GTK_RESPONSE_ACCEPT ) {

	// Process filename extension
	GtkTreeModel *model;
	GtkTreeIter iter;
	int active = 0;
	gtk_tree_selection_get_selected( selection, &model, &iter );
	gtk_tree_model_get( model, &iter, 2, &active, -1 );
	char *filename = gtk_file_chooser_get_filename( GTK_FILE_CHOOSER(_dialog) );

	if( active == 0 ) {
	    // Automatic file type from extension
	    active = type_from_extension( filename );
	    if( active == 0 ) {
		// Print error message
		GtkWidget *error = gtk_message_dialog_new( GTK_WINDOW(_dialog), GTK_DIALOG_MODAL, 
							   GTK_MESSAGE_ERROR, GTK_BUTTONS_OK,
							   "Unknown file name extension" );
		gtk_widget_show_all( error );
		gtk_dialog_run( GTK_DIALOG(error) );
		gtk_widget_destroy( error );

		gtk_widget_destroy( _dialog );
		return;
	    }
#ifndef CAIRO_HAS_PDF_SURFACE
	    if( active == 4 ) {
		// Print error message
		GtkWidget *error = gtk_message_dialog_new( GTK_WINDOW(_dialog), GTK_DIALOG_MODAL, 
							   GTK_MESSAGE_ERROR, GTK_BUTTONS_OK,
							   "No pdf output compiled in" );
		gtk_widget_show_all( error );
		gtk_dialog_run( GTK_DIALOG(error) );
		gtk_widget_destroy( error );

		gtk_widget_destroy( _dialog );
		return;		
	    }
#endif
	}

	if( _frame == NULL && active != 1 ) {
	    GtkWidget *error = gtk_message_dialog_new( GTK_WINDOW(_dialog), GTK_DIALOG_MODAL, 
						       GTK_MESSAGE_ERROR, GTK_BUTTONS_OK,
						       "Only png output possible from 3D viewer" );
	    gtk_widget_show_all( error );
	    gtk_dialog_run( GTK_DIALOG(error) );
	    gtk_widget_destroy( error );
	    
	    gtk_widget_destroy( _dialog );
	    return;
	}

	// Get resolution
	_width = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(_spinx) );
	_height = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(_spiny) );

	switch( active ) {
	case 1:
	    write_png( filename );
	    break;
	case 2:
	    write_svg( filename );
	    break;
	case 3:
	    write_eps( filename );
	    break;
#ifdef CAIRO_HAS_PDF_SURFACE
	case 4:
	    write_pdf( filename );
	    break;
#endif
	default:
	    throw( Error( ERROR_LOCATION, "unsupported file type" ) );
	}

	g_free( filename );
    }

    gtk_widget_destroy( _dialog );

}




/* ****************************************************** *
 * PNG
 * ****************************************************** */

void GTKHardcopy::get_image_size( cairo_surface_t *p_surface, 
				  int &width, int &height )
{
    int i, j, stride, w, h;
    unsigned char *buf;

    w = h = 0;
    buf = cairo_image_surface_get_data( p_surface );
    stride = cairo_image_surface_get_stride( p_surface );

    for( j = 0; j < cairo_image_surface_get_height( p_surface ); j++ ) {
        for( i = 0; i < cairo_image_surface_get_width( p_surface ); i++ ) {
            if( buf[4*i+j*stride+1] != 255 || /* R */
                buf[4*i+j*stride+2] != 255 || /* G */ 
                buf[4*i+j*stride+3] != 255 )  /* B */{
                if( i > w )
                    w = i;
                if( j > h )
                    h = j;
            }
        }
    }

    width = w+1;
    height = h+1;
}


void GTKHardcopy::unpremultiply_data( png_structp png, 
				      png_row_infop row_info, 
				      png_bytep data )
{
    unsigned int i;

    for( i = 0; i < row_info->rowbytes; i += 4 ) {
        unsigned char *b = &data[i];
        unsigned long pixel;
        unsigned char alpha;

        memcpy( &pixel, b, sizeof(unsigned long) );
        alpha = (pixel & 0xff000000) >> 24;
        if( alpha == 0 ) {
            b[0] = b[1] = b[2] = b[3] = 0;
        } else {
            b[0] = (((pixel & 0xff0000) >> 16) * 255 + alpha / 2) / alpha;
            b[1] = (((pixel & 0x00ff00) >>  8) * 255 + alpha / 2) / alpha;
            b[2] = (((pixel & 0x0000ff) >>  0) * 255 + alpha / 2) / alpha;
            b[3] = alpha;
        }
    }
}


void GTKHardcopy::write_to_png( cairo_surface_t *p_surface, 
				int width, int height, 
				const char *filename )
{
    int j, stride;
    int bpp = 8;
    FILE *fp;
    unsigned char *buf;
    png_structp png;
    png_infop   info;
    png_color_16 white;
    png_byte   *row_pointers[height];
    

    buf = cairo_image_surface_get_data( p_surface );
    stride = cairo_image_surface_get_stride( p_surface );

    if( !(fp = fopen( filename, "wb" )) )
	throw( Error( ERROR_LOCATION, "couldn\'t open file " + to_string(filename) + " for writing" ) );
    if( !(png = png_create_write_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL )) )
	throw( Error( ERROR_LOCATION, "libpng error" ) );
    if( !(info = png_create_info_struct( png )) )
	throw( Error( ERROR_LOCATION, "libpng error" ) );
    if( setjmp(png_jmpbuf(png)) )
	throw( Error( ERROR_LOCATION, "libpng error" ) );
	
    png_init_io( png, fp );
    png_set_IHDR( png, info, width, height, bpp, PNG_COLOR_TYPE_RGB_ALPHA,
		  PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE );
    white.red = 0xff;
    white.blue = 0xff;
    white.green = 0xff;
    png_set_bKGD( png, info, &white );

    png_write_info( png, info );
    png_set_write_user_transform_fn( png, unpremultiply_data );

    for( j = 0; j < height; j++ )
        row_pointers[j] = &buf[j*stride];
    png_write_image( png, row_pointers );
    png_write_end( png, info );
    png_destroy_write_struct( &png, &info );
    fclose( fp );
}


void GTKHardcopy::write_png( const char *filename )
{
    // Prepare new cairo memory surface
    cairo_surface_t *surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, _width, _height );
    if( cairo_surface_status( surface ) )
	throw( Error( ERROR_LOCATION, "error creating cairo surface" ) );

    if( _frame )
	write_png_frame( surface, filename );
    else
	write_png_geom3dplot( surface, filename );

    // Free surface
    cairo_surface_destroy( surface );
}


void GTKHardcopy::write_png_frame( cairo_surface_t *surface, const char *filename )
{
    cairo_t *cairo = cairo_create( surface );
    if( cairo_status( cairo ) )
	throw( Error( ERROR_LOCATION, "error creating cairo" ) );

    // Set cairoplot canvas size and draw
    _frame->set_geometry( _width, _height, 0, 0 );
    _frame->draw( cairo );

    // Write clipped png
    int rwidth, rheight;
    get_image_size( surface, rwidth, rheight );
    write_to_png( surface, rwidth, rheight, filename );

    cairo_destroy( cairo );
}


void GTKHardcopy::write_png_geom3dplot( cairo_surface_t *surface, const char *filename )
{
    cairo_t *cairo = cairo_create( surface );
    if( cairo_status( cairo ) )
	throw( Error( ERROR_LOCATION, "error creating cairo" ) );
    
    SoftwareRenderer r( surface );
    _geom3dplot->draw( cairo, &r );

    // Write png as is
    cairo_surface_write_to_png( surface, filename );

    cairo_destroy( cairo );
}



/* ****************************************************** *
 * EPS
 * ****************************************************** */


void GTKHardcopy::write_eps( const char *filename )
{
    // Prepare new cairo memory surface
    cairo_surface_t *surface = cairo_ps_surface_create( filename, _width, _height );
    if( cairo_surface_status( surface ) )
	throw( Error( ERROR_LOCATION, "error creating cairo surface" ) );
    cairo_t *cairo = cairo_create( surface );
    if( cairo_status( cairo ) )
	throw( Error( ERROR_LOCATION, "error creating cairo" ) );

    // Set cairoplot canvas size and draw
    _frame->set_geometry( _width, _height, 0, 0 );
    _frame->draw( cairo );

    // Free cairo and surface
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );    
}


/* ****************************************************** *
 * SVG
 * ****************************************************** */


void GTKHardcopy::write_svg( const char *filename )
{
    // Prepare new cairo memory surface
    cairo_surface_t *surface = cairo_svg_surface_create( filename, _width, _height );
    if( cairo_surface_status( surface ) )
	throw( Error( ERROR_LOCATION, "error creating cairo surface" ) );
    cairo_t *cairo = cairo_create( surface );
    if( cairo_status( cairo ) )
	throw( Error( ERROR_LOCATION, "error creating cairo" ) );

    // Set cairoplot canvas size and draw
    _frame->set_geometry( _width, _height, 0, 0 );
    _frame->draw( cairo );

    // Free cairo and surface
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );    
}


/* ****************************************************** *
 * PDF
 * ****************************************************** */


void GTKHardcopy::write_pdf( const char *filename )
{
#ifdef CAIRO_HAS_PDF_SURFACE
    // Prepare new cairo memory surface
    cairo_surface_t *surface = cairo_pdf_surface_create( filename, _width, _height );
    if( cairo_surface_status( surface ) )
	throw( Error( ERROR_LOCATION, "error creating cairo surface" ) );
    cairo_t *cairo = cairo_create( surface );
    if( cairo_status( cairo ) )
	throw( Error( ERROR_LOCATION, "error creating cairo" ) );

    // Set cairoplot canvas size and draw
    _frame->set_geometry( _width, _height, 0, 0 );
    _frame->draw( cairo );

    // Free cairo and surface
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );    
#endif
}

