/*! \file plotter.hpp
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

#ifndef PLOTTER_HPP
#define PLOTTER_HPP 1


#include <png.h>
#include <cairo.h>
#ifdef CAIRO_HAS_SVG_SURFACE
#include <cairo-svg.h>
#endif
#if CAIRO_HAS_PS_SURFACE
#include <cairo-ps.h>
#endif
#ifdef CAIRO_HAS_PDF_SURFACE
#include <cairo-pdf.h>
#endif
#include "frame.hpp"


/*! \brief Non-interactive plotter.
 *
 *  Plotter skeleton for building plots non-interactively.
 */
class Plotter {

    size_t _width;
    size_t _height;

#ifdef CAIRO_HAS_PNG_FUNCTIONS
    static void png_get_image_size( cairo_surface_t *p_surface, 
				    int &width, int &height );
    static void png_unpremultiply_data( png_structp png, 
					png_row_infop row_info, 
					png_bytep data );
    static void write_to_png( cairo_surface_t *p_surface, 
			      int width, int height, 
			      const char *filename );
#endif

    virtual void build_plot( void ) = 0;

protected:

    Frame _frame;

    /*! \brief Constructor for plotter.
     *
     *  Not intended to be used on its own.
     */
    Plotter();

    /*! \brief Destructor for plotter.
     */
    virtual ~Plotter();

public:

    /*! \brief Set size of plot.
     */
    void set_size( size_t width, size_t height ) {
	_width = width;
	_height = height;
    } 

    /*! \brief Set font size for plot.
     */
    void set_font_size( size_t size );

    /*! \brief Set ranges of plot in x- and y-directions.
     */
    void set_ranges( double xmin, double ymin, double xmax, double ymax );

    /*! \brief Get ranges of last plot in x- and y-directions.
     */
    void get_ranges( double &xmin, double &ymin, double &xmax, double &ymax );

#ifdef CAIRO_HAS_PNG_FUNCTIONS
    /*! \brief Make a plot to a PNG-file.
     *
     *  Only defined if cairo supports PNG and CAIRO_HAS_PNG_FUNCTIONS is defined.
     */
    void plot_png( const std::string &filename );
#endif

#ifdef CAIRO_HAS_PS_SURFACE
    /*! \brief Make a plot to a EPS-file.
     *
     *  Only defined if cairo supports EPS and CAIRO_HAS_EPS_SURFACE is defined.
     */
    void plot_eps( const std::string &filename );
#endif

#ifdef CAIRO_HAS_PDF_SURFACE
    /*! \brief Make a plot to a PDF-file.
     *
     *  Only defined if cairo supports PDF and CAIRO_HAS_PDF_SURFACE is defined.
     */
    void plot_pdf( const std::string &filename );
#endif

#ifdef CAIRO_HAS_SVG_SURFACE
    /*! \brief Make a plot to a SVG-file.
     *
     *  Only defined if cairo supports SVG and CAIRO_HAS_SVG_SURFACE is defined.
     */
    void plot_svg( const std::string &filename );
#endif

};



#endif






