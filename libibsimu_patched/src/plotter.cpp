/*! \file plotter.cpp
 *  \brief Basis for file output plotters
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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

#include "plotter.hpp"
#include "ibsimu.hpp"



Plotter::Plotter()
    : _width(640), _height(480)
{

}


Plotter::~Plotter()
{

}


void Plotter::set_font_size( size_t size )
{
    _frame.set_font_size( size );
}


void Plotter::set_ranges( double xmin, double ymin, double xmax, double ymax )
{
    _frame.set_ranges( PLOT_AXIS_X1, xmin, xmax );
    _frame.set_ranges( PLOT_AXIS_Y1, ymin, ymax );
}


void Plotter::get_ranges( double &xmin, double &ymin, double &xmax, double &ymax )
{
    _frame.get_ranges( PLOT_AXIS_X1, xmin, xmax );
    _frame.get_ranges( PLOT_AXIS_Y1, ymin, ymax );
}


/* *******************************************************************************
 * PNG
 * ******************************************************************************* */


#ifdef CAIRO_HAS_PNG_FUNCTIONS
void Plotter::png_get_image_size( cairo_surface_t *p_surface, 
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


void Plotter::png_unpremultiply_data( png_structp png, 
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


void Plotter::write_to_png( cairo_surface_t *p_surface, 
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
    png_set_write_user_transform_fn( png, png_unpremultiply_data );

    for( j = 0; j < height; j++ )
        row_pointers[j] = &buf[j*stride];
    png_write_image( png, row_pointers );
    png_write_end( png, info );
    png_destroy_write_struct( &png, &info );
    fclose( fp );
}


void Plotter::plot_png( const std::string &filename )
{
    ibsimu.message( 1 ) << "Plotting to PNG-file \"" << filename << "\"\n";    

    // Build plot
    build_plot();

    // Prepare new cairo memory surface
    cairo_surface_t *surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, _width, _height );
    if( cairo_surface_status( surface ) )
        throw( Error( ERROR_LOCATION, "error creating cairo surface" ) );
    cairo_t *cairo = cairo_create( surface );
    cairo_status_t status = cairo_status( cairo );
    if( status )
        throw( Error( ERROR_LOCATION, (std::string)"error creating cairo: " 
		      +  cairo_status_to_string(status) ) );

    // Clear canvas to white and draw frame
    cairo_rectangle( cairo, 0.0, 0.0, _width, _height );
    cairo_set_source_rgba( cairo, 1.0, 1.0, 1.0, 1.0 );
    cairo_fill( cairo );
    _frame.set_geometry( _width, _height, 0, 0 );
    _frame.draw( cairo );

    // Write png
    int rwidth, rheight;
    png_get_image_size( surface, rwidth, rheight );
    write_to_png( surface, rwidth, rheight, filename.c_str() );

    // Free cairo and surface
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );
}
#endif


/* *******************************************************************************
 * EPS
 * ******************************************************************************* */


#ifdef CAIRO_HAS_PS_SURFACE
void Plotter::plot_eps( const std::string &filename )
{
    ibsimu.message( 1 ) << "Plotting to EPS-file \"" << filename << "\"\n";    

    // Build plot
    build_plot();

    // Prepare new cairo memory surface
    cairo_surface_t *surface = cairo_ps_surface_create( filename.c_str(), _width, _height );
    if( cairo_surface_status( surface ) )
        throw( Error( ERROR_LOCATION, "error creating cairo surface" ) );
#if CAIRO_VERSION >= CAIRO_VERSION_ENCODE(1, 6, 0)
    cairo_ps_surface_set_eps( surface, true );
#endif
    cairo_t *cairo = cairo_create( surface );
    cairo_status_t status = cairo_status( cairo );
    if( status )
        throw( Error( ERROR_LOCATION, (std::string)"error creating cairo: " 
		      +  cairo_status_to_string(status) ) );

    // Clear canvas to white and draw frame
    _frame.set_geometry( _width, _height, 0, 0 );
    _frame.draw( cairo );

    // Free cairo and surface
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );
}
#endif


/* *******************************************************************************
 * PDF
 * ******************************************************************************* */



#ifdef CAIRO_HAS_PDF_SURFACE
void Plotter::plot_pdf( const std::string &filename )
{
    ibsimu.message( 1 ) << "Plotting to PDF-file \"" << filename << "\"\n";    

    // Build plot
    build_plot();

    // Prepare new cairo memory surface
    cairo_surface_t *surface = cairo_pdf_surface_create( filename.c_str(), _width, _height );
    if( cairo_surface_status( surface ) )
        throw( Error( ERROR_LOCATION, "error creating cairo surface" ) );
    cairo_t *cairo = cairo_create( surface );
    cairo_status_t status = cairo_status( cairo );
    if( status )
        throw( Error( ERROR_LOCATION, (std::string)"error creating cairo: " 
		      +  cairo_status_to_string(status) ) );

    _frame.set_geometry( _width, _height, 0, 0 );
    _frame.draw( cairo );

    // Free cairo and surface
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );
}
#endif


/* *******************************************************************************
 * SVG
 * ******************************************************************************* */


#ifdef CAIRO_HAS_SVG_SURFACE
void Plotter::plot_svg( const std::string &filename )
{
    ibsimu.message( 1 ) << "Plotting to SVG-file \"" << filename << "\"\n";    

    // Build plot
    build_plot();

    // Prepare new cairo memory surface
    cairo_surface_t *surface = cairo_svg_surface_create( filename.c_str(), _width, _height );
    if( cairo_surface_status( surface ) )
        throw( Error( ERROR_LOCATION, "error creating cairo surface" ) );
    cairo_t *cairo = cairo_create( surface );
    cairo_status_t status = cairo_status( cairo );
    if( status )
        throw( Error( ERROR_LOCATION, (std::string)"error creating cairo: " 
		      +  cairo_status_to_string(status) ) );

    // Set cairoplot canvas size and draw
    _frame.set_geometry( _width, _height, 0, 0 );
    _frame.draw( cairo );

    // Free cairo and surface
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );
}
#endif






