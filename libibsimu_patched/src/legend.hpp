/*! \file legend.hpp
 *  \brief Plot legends
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

#ifndef LEGEND_HPP
#define LEGEND_HPP 1


#include <vector>
#include <string>
#include "graph.hpp"
#include "label.hpp"
#include "vec3d.hpp"
#include "colormap.hpp"


/*! \brief Legend position. 
 *
 *  Legend location enum/bitmask.
 */
enum legend_position_e {
    LEGEND_POS_BOTTOM_LEFT   = 0,
    LEGEND_POS_BOTTOM_CENTER = 1,
    LEGEND_POS_BOTTOM_RIGHT  = 2,

    LEGEND_POS_MIDDLE_LEFT   = 4+0,
    LEGEND_POS_MIDDLE_CENTER = 4+1,
    LEGEND_POS_MIDDLE_RIGHT  = 4+2,

    LEGEND_POS_TOP_LEFT      = 8+0,
    LEGEND_POS_TOP_CENTER    = 8+1,
    LEGEND_POS_TOP_RIGHT     = 8+2
};


#define LEGEND_POS_VERTICAL_MASK   12
#define LEGEND_POS_HORIZONTAL_MASK 3

#define LEGEND_POS_BOTTOM 0
#define LEGEND_POS_MIDDLE 4
#define LEGEND_POS_TOP    8

#define LEGEND_POS_LEFT   0
#define LEGEND_POS_CENTER 1
#define LEGEND_POS_RIGHT  2


/*! \brief Class for legend entry.
 *
 *  The legend entries contain a reference to the graph drawn so that
 *  if the style used in the graph is changed, the sample in legend
 *  is automatically changed.
 *
 *  Plotting of graph sample is done by plot_sample() in Graph.
 */
class LegendEntry {
    
    Graph        &_graph;      /*!< \brief Reference to graph drawn. */
    Label         _label;      /*!< \brief Label for legend entry. */

public:
    
    /*! \brief Contructor for legend entry.
     */
    LegendEntry( Graph &graph, const std::string &label ) 
	: _graph(graph), _label(label) {}

    /*! \brief Copy constructor.
     */
    LegendEntry( const LegendEntry &le ) 
	: _graph(le._graph), _label(le._label) {}

    /*! \brief Destructor.
     */
    ~LegendEntry() {}

    /*! \brief Assignment operator.
     */
    LegendEntry &operator=( const LegendEntry &le ) {
	_graph = le._graph;
	_label = le._label;
	return( *this );
    }

    /*! \brief Plot legend entry at \a (x,y).
     *
     *  The point \a (x,y) is the lower left point of the entry.
     */
    void plot( cairo_t *cairo, double x, double y );

    /*! \brief Get size of legend entry.
     */
    void get_size( cairo_t *cairo, double &width, double &height );

    /*! \brief Set font size for legend labels.
     */
    void set_font_size( double fontsize );
};


/*! \brief Base class for legend definition.
 *
 *  Legend is an object that contains a key to the plot styles used in
 *  graphs. The key contains a sample of the plot style used and a
 *  corresponding text label.
 *
 *  The Colormap legend is a special case because in addition to the
 *  plot style, the plot z-range is shown in the legend.
 *
 *  The size of legend can be queried and the location can be set.
 */
class Legend {

public:

    /*! \brief Default constructor for legend.
     */
    Legend() {}

    /*! \brief Virtual destructor.
     */
    virtual ~Legend() {}    

    /*! \brief Plot legend at \a (x,y).
     *
     *  The point \a (x,y) is the lower left point of the entry.
     */
    virtual void plot( cairo_t *cairo, double x, double y ) = 0;

    /*! \brief Get size of legend.
     */
    virtual void get_size( cairo_t *cairo, double &width, double &height ) = 0;
};


/*! \brief %Legend for presenting plot styles.
 */
class MultiEntryLegend : public Legend {

    double                     _fontsize; /*!< \brief Font size for labels. */
    std::vector<LegendEntry *> _entry;    /*!< \brief Legend entries. */

public:

    /*! \brief Default constructor for legend.
     */
    MultiEntryLegend();
    
    /*! \brief Virtual destructor.
     */
    virtual ~MultiEntryLegend() {}    

    /*! \brief Plot legend at \a (x,y).
     *
     *  The point \a (x,y) is the lower left point of the entry.
     */
    virtual void plot( cairo_t *cairo, double x, double y );

    /*! \brief Get size of legend.
     */
    virtual void get_size( cairo_t *cairo, double &width, double &height );

    /*! \brief Set font size for legend labels.
     */
    void set_font_size( double fontsize );

    /*! \brief Get font size
     */
    double get_font_size( void );

    /*! \brief Add entry to legend.
     */
    void add_entry( LegendEntry *entry );

    /*! \brief Clear legend entries.
     */
    void clear_entries( void );
};


/*! \brief %Legend for presenting colormap key.
 *
 *  Doesn't use Ruler,
 */
class ColormapLegend : public Legend {

    struct Tic {
	double  _loc;       /*!< \brief Tic height location. */
	Label   _label;     /*!< \brief Tic label. */

	/*! \brief Constructor for coordinate tic.
	 */
	Tic( double loc, const std::string &label ) : _loc(loc), _label(label) {}
    };

    double      _width;      /*!< \brief Width of legend box. */
    double      _height;     /*!< \brief Height of legend box. */

    double      _fontsize;   /*!< \brief Font size for labels. */
    Vec3D       _color;

    double      _ticlen_in;
    double      _ticlen_out;
    double      _ticspace;

    double           _range[2]; /*!< \brief z-range. */
    std::vector<Tic> _tic;

    Colormap   &_colormap;

    void build_legend( double x, double y );
    void plot_colormap_palette_to_image_surface( cairo_surface_t *surface, int plim[4] );
    void plot_colomap_palette( cairo_t *cairo, int plim[4] );
    
public:

    /*! \brief Default constructor for legend.
     */
    ColormapLegend( Colormap &colormap );

    /*! \brief Virtual destructor.
     */
    virtual ~ColormapLegend() {}

    /*! \brief Plot legend at \a (x,y).
     *
     *  The point \a (x,y) is the lower left point of the legend.
     */
    virtual void plot( cairo_t *cairo, double x, double y );

    /*! \brief Get size of legend.
     */
    virtual void get_size( cairo_t *cairo, double &width, double &height );

    /*! \brief Set font size for legend labels.
     */
    void set_font_size( double fontsize );

    /*! \brief Set ruler color.
     */
    void set_color( const Vec3D &color );

    /*! \brief Set height of legend.
     */
    void set_height( double height );
};


#endif

