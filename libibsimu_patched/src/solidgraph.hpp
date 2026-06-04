/*! \file solidgraph.hpp
 *  \brief %Graph for plotting solids
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

#ifndef SOLIDGRAPH_HPP
#define SOLIDGRAPH_HPP 1


#include <vector>
#include "geometry.hpp"
#include "graph3d.hpp"
#include "vec3d.hpp"
#include "lineclip.hpp"


/*! \brief A 2D cut view of the geometry solids.
 *
 *  Class for constructing and plotting a view of the geometry
 *  solids. The view data is stored inside the object in a cache to
 *  speed up frequent use (in interactive plotter).
 */
class SolidGraph : public Graph3D {

    struct Point {
	double x[2];

	Point( double _x, double _y ) { x[0] = _x; x[1] = _y; }
	
	double &operator[]( int32_t i ) { return( x[i] ); }
	const double &operator[]( int32_t i ) const { return( x[i] ); }
    };

    class SolidPoints {
	uint32_t           _N;    /* Solid number of electrode, N>=7 */
	std::vector<Point> _p;    /* Coordinate points of electrode boundary */
    public:
	SolidPoints( uint32_t N ) : _N(N) {}
	
	void add_point( Point x );
	uint32_t N( void ) { return( _N ); }
	size_t size( void ) { return( _p.size() ); }
	Point &operator[]( int32_t i ) { return( _p[i] ); }
	const Point &operator[]( int32_t i ) const { return( _p[i] ); }
    };

    const Geometry                      &_geom;

    Vec3D                                _color;
    std::vector<SolidPoints *>           _solid;

    view_e                               _oview;
    int                                  _olevel;

    bool                                 _cache;

    static void step( int32_t &nextx, int32_t &nexty, int32_t dir );
    static bool is_edge( uint32_t N, uint32_t node1, uint32_t node2 );
    void loop( SolidPoints *solid, int32_t x, int32_t y,
	       char *done, uint32_t sizex, uint32_t sizey );
    void build_solid( uint32_t N, char *done, uint32_t sizex, uint32_t sizey );
    void build_data( void );
    void plot_data( cairo_t *cairo, LineClip &lc, const Coordmapper *cm );
    void clear_data( void );
    
public:

    /*! \brief Constructor for %SolidGraph drawable from geometry \a g.
     */
    SolidGraph( const Geometry &geom );

    /*! \brief Destructor.
     */
    virtual ~SolidGraph();

    /*! \brief Disable internal cache.
     *
     *  Makes solid boundaries to be calculated at every plot().
     */
    void disable_cache( void );

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
     *  Plot graph sample for legend at cairo coordinates \a x.
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
