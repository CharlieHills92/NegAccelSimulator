/*! \file particlegraph.hpp
 *  \brief %Graph for particle plots
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

#ifndef PARTICLEGRAPH_HPP
#define PARTICLEGRAPH_HPP 1


#include <cairo.h>
#include <vector>
#include "vec3d.hpp"
#include "graph3d.hpp"
#include "coordmapper.hpp"
#include "geometry.hpp"
#include "particledatabase.hpp"
#include "lineclip.hpp"


/*! \brief Class for drawing particle trajectories.
 *
 *  Draws particle trajectories. Using the \a particlediv variable
 *  only one per \a particlediv trajectories is plotted. The different
 *  q/m values are discriminated by colors if enabled with \a
 *  qm_dircr. The trajectories are drawn with straight lines if the
 *  particle database interpolation is set to linear or curves if
 *  polynomial interpolation is used.
 *
 *  Implementation of %Graph3D.
 */
class ParticleGraph : public Graph3D {
    
    const Geometry         &_geom;            /*!< \brief Reference to simulation geometry. */
    const ParticleDataBase &_pdb;             /*!< \brief Reference to particle database. */
    uint32_t                _particle_div;    /*!< \brief Particle plot divisor. */
    uint32_t                _particle_offset; /*!< \brief Particle plot offset. */

    std::vector<Vec3D>      _color;           /*!< \brief Colors for trajectories. */

    double                  _ox[5];           /*!< \brief Workspace for particleplot_draw_curve() */
    size_t                  _coordsize;       /*!< \brief Size of array _coord divided by two */
    double                 *_coord;           /*!< \brief Workspace for particleplot_draw_curve() */
    bool                    _qm_discr;        /*!< \brief q/m discriminator enable, default true */

    void get_point( const Coordmapper *cm, double *coord, double s, 
		    double Ax, double Bx, double Cx, double Dx, 
		    double Ay, double By, double Cy, double Dy ) const;

    void draw_linear( const Coordmapper *cm, LineClip &lc, 
		      double x[5], bool first ) const;

    void draw_curve( const Coordmapper *cm, LineClip &lc, 
		     double x[5], bool first );

public:

    /*! \brief Constructor for particle plotter.
     */
    ParticleGraph( const Geometry &g, const ParticleDataBase &pdb,
		   uint32_t particle_div = 11, uint32_t particle_offset = 0, 
		   bool qm_discr = true );
    
    /*! \brief Destructor.
     */
    virtual ~ParticleGraph();

    /*! \brief Set particle divisor and offset.
     *
     *  Set \a particle_div to zero for no plotting, one for plotting
     *  every particle, two for plotting every second particle, three
     *  for plotting every third particle, etc. Defaults to
     *  11. Plotter skips the first \a particle_offset particles.
     */
    void set_particle_div( uint32_t particle_div, uint32_t particle_offset );

    /*! \brief Enable q/m discretation
     *
     *  If q/m discretation is enabled, different q/m particles will
     *  be plotted with different colors. Default is enabled.
     */
    void set_qm_discretation( bool qm_discr );

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

    /*! \brief Get bounding box of graph.
     *
     *  Returns the bounding box of the graph in array \a bbox in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void get_bbox( double bbox[4] );

    /*! \brief Add a color to the list of trajectory colors.
     */
    void add_color( const Vec3D &color );

    /*! \brief Clear the list of trajectory colors.
     */
    void clear_colors( void );

};



#endif
