/*! \file lineclip.hpp
 *  \brief Floating point line clipping for cairo
 */

/* Copyright (c) 2005-2010,2012 Taneli Kalvas. All rights reserved.
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

#ifndef LINECLIP_HPP
#define LINECLIP_HPP 1

#include <cairo.h>


/*! \brief Line clipper.
 *
 *  Cairo graphics coordinates are internally handled using fixed
 *  point algebra for speed. This causes problems in clipping
 *  algorithm when there are large scale differences in
 *  coordinates. This class is provided for the user to overcome this
 *  problem by using floating point algebra for line clipping.
 */
class LineClip {
    cairo_t  *p_dc;          /* Cairo context */
    double    clip[4];       /* xmin, ymin, xmax, ymax  */
    
    double    first[2];      /* Last user given moveto point (first point of path) */
    
    double    last[2];       /* last user given point, nan if not available */
    int       last_outcode;  /* outcode of last user given point */
    int       last_op;       /* last operation 0=lineto, 1=moveto, 2=no op */
    
    double    drawn[2];      /* last drawn point */
    int       drawn_outcode; /* outcode of last drawn point */
    
    int       coord_alloc;   /* Allocated size of coordinate database */
    double   *coord;         /* Buffer for storing line coordinates in curve_to */

    int outcode( double x, double y );
    int exit_outcode( double x, double y );
    void get_point( double *coords, double t,
		    double x0, double y0,
		    double x1, double y1,
		    double x2, double y2,
		    double x3, double y3 );

public:

    /*! \brief Construct line clipper
     */
    LineClip( cairo_t *cairo );

    /*! \brief Destructor
     */
    ~LineClip();

    /*! \brief Set clipping area
     */
    void set( double xmin, double ymin, double xmax, double ymax );

    /*! \brief Reset clip
     *
     *  Sets clip limits to infinity and forgets last coordinate.
     */
    void reset();

    /*! \brief Move to (x,y)
     */
    void move_to( double x, double y );

    /*! \brief Line to (x,y)
     */
    void line_to( double x, double y );

    /*! \brief Curve to (x,y)
     *
     *  Curve is drawn with lines in cairo. Number is subdivisions is
     *  made large enough to contain the error at typically less than
     *  1 pixel.
     */
    void curve_to( double x1, double y1,
		   double x2, double y2,
		   double x3, double y3 );

    /*! \brief Close path
     */
    void close_path();

    /*! \brief Close path and fill enclosed area
     */
    void fill();
};


#endif
