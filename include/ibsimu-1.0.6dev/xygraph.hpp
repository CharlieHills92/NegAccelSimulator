/*! \file xygraph.hpp
 *  \brief XY-graph
 */

/* Copyright (c) 2005-2011,2013,2014 Taneli Kalvas. All rights reserved.
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

#ifndef XYGRAPH_HPP
#define XYGRAPH_HPP 1


#include <cairo.h>
#include <vector>
#include "vec3d.hpp"
#include "graph.hpp"
#include "coordmapper.hpp"


/*! \brief XYGraph line style
 */
enum line_style_e {
    XYGRAPH_LINE_DISABLE = 0,   /*!< \brief Disable line plotting */
    XYGRAPH_LINE_SOLID          /*!< \brief Regular line between data points */
};


/*! \brief XYGraph point style
 */
enum point_style_e {
    XYGRAPH_POINT_DISABLE = 0,   /*!< \brief Disable data point plotting */
    XYGRAPH_POINT_CIRCLE,        /*!< \brief Draw a circle */
    XYGRAPH_POINT_BOX            /*!< \brief Draw a box */
};


/*! \brief Class for XY-type simple graph plots. 
 *
 *  Implementation of %Graph. Used in Frame type plots.
 */
class XYGraph : public Graph {

    double                 _linewidth;
    Vec3D                  _color;
    line_style_e           _linestyle;
    point_style_e          _pointstyle;
    bool                   _point_filled;
    double                 _point_scale;
    bool                   _histogram;
    bool                   _extend_histogram;
    
    std::vector<double>    _xdata;
    std::vector<double>    _ydata;

    void plot_point( cairo_t *cairo, double x, double y );
    void plot_standard_lines( cairo_t *cairo, const Coordmapper *cm, const double range[4], size_t N );
    void plot_histogram_lines( cairo_t *cairo, const Coordmapper *cm, const double range[4], size_t N );

public:

    /*! \brief Default constructor for empty graph.
     */
    XYGraph();

    /*! \brief Constructor for basic graph with datapoints \a xdata
     *  and \a ydata.
     *
     *  Internal copies of the data from \a xdata and \a ydata are made.
     */
    XYGraph( const std::vector<double> &xdata, 
	     const std::vector<double> &ydata );

    /*! \brief Destructor.
     */
    virtual ~XYGraph() {}

    /*! \brief Plot graph with cairo.
     *
     *  Plot the graph using \a cairo and coordinate mapper \a cm. The
     *  visible range of plot is given in array \a range in order \a
     *  xmin, \a ymin, \a xmax, \a ymax. The graph should be able to
     *  handle any range values. Also \a min > \a max.
     *
     *  Called by Frame during drawing.
     */
    virtual void plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] );

    /*! \brief Plot sample for legend.
     *
     *  Plot graph sample for legend at cairo coordinates \a (x,y).
     */
    virtual void plot_sample( cairo_t *cairo, double x, double y, double width, double height );

    /*! \brief Get bounding box of graph.
     *
     *  Returns the bounding box of the graph in array \a bbox in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void get_bbox( double bbox[4] );

    /*! \brief Set new data arrays.
     */
    void set_data( const std::vector<double> &xdata, 
		   const std::vector<double> &ydata );

    /*! \brief Set line width.
     *
     *  Default to width 1.0.
     */
    void set_line_width( double linewidth );

    /*! \brief Set graph color.
     *
     *  Defaults to red (1,0,0).
     */
    void set_color( const Vec3D &color );

    /*! \brief Set line style.
     *
     *  Defaults to no lines drawn (XYGRAPH_LINE_DISABLE)
     */
    void set_line_style( line_style_e linestyle,
			 double linewidth = 1.0 );

    /*! \brief Set point style.
     *
     *  Defaults to filled XYGRAPH_POINT_CIRCLE with scale 3.0.
     */
    void set_point_style( point_style_e pointstyle, bool filled = true, double scale = 1.0 );

    /*! \brief Set histogram style.
     *
     *  Set to true for histogram style plots.
     */
    void set_histogram( bool histo );

    /*! \brief Extend histogram.
     *
     *  Set to true for histogram to be extended in x-direction to
     *  cover all of range.
     */
    void extend_histogram( bool extend );

};


#endif
