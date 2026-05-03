/*! \file frame.hpp
 *  \brief %Frame for plots
 */

/* Copyright (c) 2005-2013 Taneli Kalvas. All rights reserved.
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

#ifndef FRAME_HPP
#define FRAME_HPP 1


#include <string>
#include <vector>
#include <cairo.h>
#include "vec3d.hpp"
#include "ruler.hpp"
#include "coordmapper.hpp"
#include "graph.hpp"
#include "legend.hpp"



/*! \brief Fixed aspect ratio setting for frame
 *
 *  Default mode is PLOT_FIXED_ASPECT_DISABLED, which means that the
 *  frame axes are independent of each other. The mode
 *  PLOT_FIXED_ASPECT_EXTEND_RANGE has fixed aspect ratio, which is
 *  made by preserving the plot margins and extending the range of the
 *  smaller axis. The mode PLOT_FIXED_ASPECT_INCREASE_MARGIN fixes the
 *  aspect ratio by increasing the margin of the larger axis.
 */
enum PlotFixedMode {
    PLOT_FIXED_ASPECT_DISABLED = 0,
    PLOT_FIXED_ASPECT_EXTEND_RANGE,
    PLOT_FIXED_ASPECT_INCREASE_MARGIN
};


/*! \brief Axis specification for plots.
 *
 *  PLOT_AXIS_X1 is the bottom axis, PLOT_AXIS_Y1 is the left side
 *  axis, PLOT_AXIS_X2 is the top axis, PLOT_AXIS_Y2 is the right side
 *  axis and PLOT_AXIS_Z is the depth axis for three dimensional plots
 *  or colormap plots.
 */
enum PlotAxis {
    PLOT_AXIS_X1 = 0,
    PLOT_AXIS_Y1,
    PLOT_AXIS_X2,
    PLOT_AXIS_Y2,
    PLOT_AXIS_Z
};



/*! \brief General purpose plotter frame for cairographics surfaces.
 *
 *  Plotter frame provides a coordinate frame and simple linear
 *  independent coordinate transformations for axes from logical plot
 *  coordinates to canvas coordinates. Also the plotter can do
 *  numbered tics to the frame and labels for x- and y-axes. Frame x-
 *  and y-ranges can be set and coordinates can be forced to fixed
 *  aspect ratio. Content to the plot is drawn with Graph objects,
 *  which can be inserted into and removed from the plot.
 *
 *  Ruler autorange is different from frame/graph autorange:
 *  - Frame autorange: Frame asks graph of the range they need 
 *    and sets the ruler range according to graph requests
 *  - Ruler autorange: Rounds up the set range to the next tic mark.
 *
 */
class Frame {

    /*! \brief Drawable object in frame
     */
    struct DObj {
	PlotAxis        _xaxis;         /*!< \brief X-axis for graph. */
	PlotAxis        _yaxis;         /*!< \brief Y-axis for graph. */
	Graph          *_graph;         /*!< \brief Graph to plot in the frame. */

	DObj( PlotAxis xaxis, PlotAxis yaxis, Graph *graph )
	    : _xaxis(xaxis), _yaxis(yaxis), _graph(graph) {}
    };

    Ruler               _ruler[4];      /*!< \brief Rulers for x1, y1, x2, y2 axes. */
    Coordmapper1D       _cm[4];         /*!< \brief Coordinate mappers for x1, y1, x2, y2 axes. */
    bool                _enable[4];     /*!< \brief Ruler enabled. */
    bool                _fenable[4];    /*!< \brief Ruler force enable, false = auto, true = force enable. */
    bool                _autorange[8];  /*!< \brief Frame autorange enable for x1, y1, x2, y2 axes. */
    double              _range_min[4];  /*!< \brief Range (auto or manual) value for x1, y1, x2, y2 axes.  */
    double              _range_max[4];  /*!< \brief Range (auto or manual) value for x1, y1, x2, y2 axes.  */

    int                 _offx;          /*!< \brief Offset in x-direction. */
    int                 _offy;          /*!< \brief Offset in y-direction. */

    int                 _width;         /*!< \brief Width of plot. */
    int                 _height;        /*!< \brief Height of plot. */

    double              _fontsize;      /*!< \brief Font size for labels. */
    double              _titlespace;    /*!< \brief Space from frame title to next item down. */
    double              _cmlspace;      /*!< \brief Space from colormap legend to next item on the left. */
    Vec3D               _bg;            /*!< \brief Background color. */
    Vec3D               _fg;            /*!< \brief Frame color. */

    std::vector<DObj>   _dobj;          /*!< \brief Graphs to plot in the frame. */

    bool                _legend_enable; /*!< \brief Legend enable. */
    legend_position_e   _legend_pos;    /*!< \brief Legend position. */
    MultiEntryLegend    _legend;        /*!< \brief Graph legend. */
    ColormapLegend     *_cm_legend;     /*!< \brief Colormap legend. */
    bool                _cml_enable;    /*!< \brief Colormap legend enable. */

    Label               _title;         /*!< \brief Title for plot. */

    PlotFixedMode       _fixedaspect;   /*!< \brief Fixed aspect ration to 1:1, 0 = disabled, 
					 * 1 = extend range, 2 = increase margin. */
    bool                _automargin;    /*!< \brief Automatic margin calculation. */
    double              _tmargin[4];    /*!< \brief Frame total margins from edge of plot to 
					 * edge of frame.
					 *
					 * If automargin=true, the total margin is calculated
				 	 * by adding the amount of space needed for
				 	 * tics, labels, legends and titles (xmin, ymin, xmax, ymax). */

    void build_colormap_legend( void );
    void draw_colormap_legend( cairo_t *cairo );
    void draw_legend( cairo_t *cairo );
    void calculate_autoranging( void );
    void calculate_ruler_autoenable( void );
    void calculate_rulers( cairo_t *cairo, bool ruler_tic_bbox_test );
    void calculate_frame( cairo_t *cairo );
    void draw_frame( cairo_t *cairo );
    void mirror_rulers( void );
    void set_frame_clipping( cairo_t *cairo );
    void unset_frame_clipping( cairo_t *cairo );

public:

    /*! \brief Default constructor.
     */
    Frame();

    /*! \brief Copy constructor.
     */
    Frame( const Frame &frame );

    /*! \brief Destructor.
     */
    ~Frame();

    /*! \brief Set size of plot.
     *
     *  Sets the size of the plot to \a width by \a height and offsets
     *  \a offx and \a offy. Inside this size are the user defined
     *  margins and possibly automatic margin and the plot area
     *  itself. See set_automargin() for more information.
     */
    void set_geometry( int width, int height, int offx, int offy ) {
	_width = width;
	_height = height;
	_offx = offx;
	_offy = offy;
    }

    /*! \brief Set font size used for tics and labels.
     */
    void set_font_size( double size );

    /*! \brief Get font size used for tics and labels.
     */
    double get_font_size( void ) {
	return( _fontsize );
    }

    /*! \brief Set background color.
     */
    void set_background( Vec3D &bg ) {
	_bg = bg;
    }    

    /*! \brief Set foreground color.
     */
    void set_foreground( Vec3D &fg ) {
	_fg = fg;
    }

    /*! \brief Get a 2d coordinate mapper for selected axes.
     *
     *  The returned coordinate mapper is the last mapper used and
     *  built while drawing a plot.
     */
    Coordmapper get_coordmapper( PlotAxis xaxis, PlotAxis yaxis ) const;

    /*! \brief Get margins.
     */
    void get_margins( double margin[4] ) const;

    /*! \brief Get frame edge locations.
     *
     *  Frame edge location coordinates are returned in order
     *  (left(x), top(y), right(x), bottom(y)).
     */
    void get_frame_edges( double edge[4] ) const;

    /*! \brief Set title for plot.
     */
    void set_title( const std::string &title );

    /*! \brief Set axis label.
     */
    void set_axis_label( PlotAxis axis, const std::string &label );

    /*! \brief Force enable ruler for \a axis.
     *
     *  Normally rulers for axes that are not referenced by graphs
     *  are not shown and mirrored tics from the other ruler are shown
     *  instead. Ruler enable is automatic by default.
     */
    void force_enable_ruler( PlotAxis axis, bool force );

    /*! \brief Enable/disable ruler autorange.
     *
     *  Autoranging in rulers rounds up/down ("outwards") the ruler
     *  ranges to the next closest tic mark. Ruler autorange is
     *  enabled by default.
     */
    void ruler_autorange_enable( PlotAxis axis, bool min, bool max );

    /*! \brief Set coordinate ranges for axis
     *
     *  Use plus or minus infinite for full autoranging according to
     *  graph bounding boxes. Use finite numbers for manual range
     *  setting.
     */
    void set_ranges( PlotAxis axis, double min, double max );

    /*! \brief Get coordinate ranges for axis.
     *
     *  After draw() has been called, this function returns real
     *  ranges used.
     */
    void get_ranges( PlotAxis axis, double &min, double &max ) const;
	
    /*! \brief Set fixed aspect ratio mode
     *
     *  When fixed aspect ratio mode is enabled, the coordinates are
     *  guarateed to have fixed 1:1 aspect ratio. If extend is true,
     *  the ranges are adjusted for plot. The coordinate axis to be
     *  adjusted is chosen so that the image size is increased, never
     *  decreased. If extend is false, the ranges will stay fixed, but
     *  the margins of the plot are adjusted for correct aspect ratio.
     *
     *  Fixed aspect ratio mode is disabled by default.
     */
    void set_fixed_aspect( PlotFixedMode mode );

    /*! \brief Set automatic margin setting on or off.
     *
     *  When automargin is enabled the margins set by user are
     *  increased to fit in labels and tics. In other words the plot
     *  including all tics and labels are made to fit inside the area
     *  inside user set margins.
     *
     *  When automargin is disabled the plot window (plot frame) is
     *  drawn exactly at the plot size minus user set margins area.
     *
     *  Automargin is enabled by default.
     */
    void set_automargin( bool enable );

    /*! \brief Add graph to frame.
     *
     *  Add a \a graph to be plotted in the frame. The graph uses the
     *  axes \a xaxis and \a yaxis. A \a legend for the graph may be
     *  defined. Legend defaults to NULL (no legend). A pointers to
     *  graph and legend is added to frame. The graph/legend is not
     *  freed by frame.
     */
    void add_graph( PlotAxis xaxis, PlotAxis yaxis, 
		    Graph *graph, LegendEntry *legend = NULL );

    /*! \brief Clear all graphs from frame.
     *
     *  The list of graphs and legend entries is cleared. The
     *  graphs/legends are not freed, only the references are erased.
     */
    void clear_graphs( void );

    /*! \brief Set colormap legend enable/disable.
     */
    void enable_colormap_legend( bool enable );

    /*! \brief Set legend enable/disable.
     */
    void enable_legend( bool enable );

    /*! \brief Set legend position. 
     */
    void set_legend_position( legend_position_e pos );

    /*! \brief Set colormap legend.
     */
    //void set_colormap_legend( ColormapLegend *legend );

    /*! \brief Draw frame and plot contents.
     *
     *  The graphs inside the frame are asked for bounding boxes. The
     *  frame geometry is then calculated according to margin, size
     *  settings and the ranges. The coordinate transformation is
     *  initialized and the graphs are called to plot themselves.
     *
     *  The frame guarantees to call Graph::get_bbox() once and
     *  Graph::plot() once for each call to this function.
     *
     *  Graphs are plotted in the order of definition. Last definition
     *  is on the top.
     */
    void draw( cairo_t *cairo );
};


#endif

