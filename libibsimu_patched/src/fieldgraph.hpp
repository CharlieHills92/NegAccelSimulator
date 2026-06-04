/*! \file fieldgraph.hpp
 *  \brief %Graph for plotting fields
 */

/* Copyright (c) 2005-2012,2014 Taneli Kalvas. All rights reserved.
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

#ifndef FIELDGRAPH_HPP
#define FIELDGRAPH_HPP 1


#include <vector>
#include "graph3d.hpp"
#include "meshcolormap.hpp"
#include "geometry.hpp"
#include "field.hpp"
#include "types.hpp"


/*! \brief Class for drawing fields with colormap
 *
 *  Implementation of %Graph3D. Used in Frame type plots.
 */
class FieldGraph : public Graph3D, public MeshColormap {

    const Geometry         &_geom;            /*!< \brief Geometry. */
    const Field            *_field;           /*!< \brief Field to be plotted. */
    field_type_e            _field_type;      /*!< \brief Field type used. */

    bool                    _first;
    view_e                  _oview;
    double                  _olevel;
    bool                    _enabled;         /*!< \brief Is plotting enabled */

    double                  _zmin;
    double                  _zmax;

public:

    /*! \brief Constructor for plotting \a field.
     */
    FieldGraph( const Geometry &geom, const Field *field, field_type_e field_type );

    /*! \brief Destructor.
     */
    virtual ~FieldGraph();

    /*! \brief Get field type.
     */
    field_type_e field_type( void ) const;
    
    /*! \brief Set field to be plotted.
     *
     *  The \a field_type can be FIELD_NONE and \a field NULL for no plotting.
     */
    //void set_field( field_type_e field_type, const ScalarField *field );

    /*! \brief Set field to be plotted.
     *
     *  The \a field_type can be FIELD_NONE and \a field NULL for no plotting.
     */
    //void set_field( field_type_e field_type, const VectorField *field );

    /*! \brief Enable/disable plot.
     */
    void enable( bool enable );

    /*! \brief Set zrange for plot.
     *
     *  The zrange defaults to automatically scaled range for the
     *  whole field for scalarfields and automatically scaled range
     *  for the view plane only for vectorfields.
     */
    void set_zrange( double zmin, double zmax );

    /*! \brief Build plot.
     *
     *  Reads in field data and produces the colormap for the
     *  plot. Automatically called by plot(), but can be called to
     *  produce automatic zranges.
     */
    void build_plot( void );
    
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
