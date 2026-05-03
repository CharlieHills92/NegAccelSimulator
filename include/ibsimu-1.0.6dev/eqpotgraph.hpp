/*! \file eqpotgraph.hpp
 *  \brief %Graph for plotting equipotential lines
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

#ifndef EQPOTGRAPH_HPP
#define EQPOTGRAPH_HPP 1


#include <vector>
#include "geometry.hpp"
#include "meshscalarfield.hpp"
#include "graph3d.hpp"
#include "vec3d.hpp"


/*! \brief Equipotential line plot.
 *
 *  Class for constructing and drawing equipotential line plots.
 */
class EqPotGraph : public Graph3D {

    struct Line {
	double x[4];              /* Coordinates of line segment (x1,y1,x2,y2) */
	
	Line( double _x1, double _y1, double _x2, double _y2 ) { 
	    x[0] = _x1; 
	    x[1] = _y1;
	    x[2] = _x2; 
	    x[3] = _y2;
	}

	double &operator[]( int i ) { return( x[i] ); }
	const double &operator[]( int i ) const { return( x[i] ); }
    };

    struct EqPotLines {
	double            pot;    /* Potential value for equipotential line. */
	std::vector<Line> x;      /* Line segment coordinates. */
	
	EqPotLines( double pot ) : pot(pot) {}
    };

    Vec3D                               _color;
    const MeshScalarField              &_epot;
    const Geometry                     &_geom;
    bool                                _data_built;

    std::vector<double>                 _eqlines_manual;
    size_t                              _eqlines_auto;
    std::vector<EqPotLines *>           _lines;

    view_e                              _oview;
    double                              _olevel;

    bool                                _cache;

    
    bool is_solid( uint32_t sol ) const;
    bool is_near( uint32_t sol ) const;
    bool eqline_exists( double pot1, uint32_t sol1, 
			double pot2, uint32_t sol2, 
			double pot ) const;
    void build_data( void );

public:

    /*! \brief Constructor for equipotential line plot.
     *
     *  Makes a plot object for plotting equipotential data from
     *  scalarfield \a field in geometry \a g.
     */
    EqPotGraph( const MeshScalarField &epot, const Geometry &geom );

    /*! \brief Destructor,
     */
    virtual ~EqPotGraph();

    /*! \brief Disable internal cache.
     *
     *  Makes equipotential lines to be calculated at every plot().
     */
    void disable_cache( void );

    /*! \brief Add manual equipotential lines to be plotted at
     *  specified potentials.
     */
    void set_eqlines_manual( const std::vector<double> &pot );

    /*! \brief Set \a N automatic equipotential lines to be plotted
     *  between minimum potential and maximum potentials.
     */
    void set_eqlines_auto( size_t N );

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

    /*! \brief Get bounding box of drawable.
     *
     *  Returns the bounding box of the drawable in array \a bbox in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void get_bbox( double bbox[4] );
};


#endif
