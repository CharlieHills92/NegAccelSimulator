/*! \file meshcolormap.hpp
 *  \brief Mesh based colormap graph for plotting
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

#ifndef MESHCOLORMAP_HPP
#define MESHCOLORMAP_HPP 1


#include <cairo.h>
#include <vector>
#include "palette.hpp"
#include "graph.hpp"
#include "colormap.hpp"
#include "coordmapper.hpp"
#include "interpolation.hpp"


/*! \brief Interpolation type enum.
 *
 *  The interpolation can be either 0th order (\a
 *  INTERPOLATION_CLOSEST), 1st order (\a INTERPOLATION_BILINEAR) or
 *  3rd order (\a INTERPOLATION_BICUBIC).
 */
enum interpolation_e {
    INTERPOLATION_CLOSEST = 0, /*!< \brief Closest point interpolation */
    INTERPOLATION_BILINEAR,    /*!< \brief Bilinear interpolation */
    INTERPOLATION_BICUBIC      /*!< \brief Bicubic interpolation */
};

/*! \brief Z-scale enum
 *
 *  The zscale can be either linear scale (\a ZSCALE_LINEAR),
 *  logarithmic scale (\a ZSCALE_LOG) or relative logarithmic (\a
 *  ZSCALE_RELLOG). The relative logarithmic scaling follows the
 *  relation \f[ \frac{\log(0.001+x)-\log(0.001)}{\log(1.001)-\log(0.001)} \f],
 *  where \a x is prescaled to range [0,1]. The z-ranges completely
 *  contained on the negative side are inverted to positive and
 *  z-ranges both on negative and positive sides are scaled separately
 *  to provide magnification close to zero.
 */
enum zscale_e {
    ZSCALE_LINEAR = 0, /*!< \brief Linear scale */
    ZSCALE_LOG,        /*!< \brief Logarithmic scale */
    ZSCALE_RELLOG      /*!< \brief Relative logarithmic scale */
};


/*! \brief Class for colormap type plots on a regular rectangular mesh.
 *
 *  Implementation of %Colormap. Used in Frame type plots.
 */
class MeshColormap : public Colormap {

    Palette                _palette;       /*!< \brief Palette for plotting. */

    interpolation_e        _interpolation; /*!< \brief Interpolation mode. */
    zscale_e               _zscale;        /*!< \brief zscale mode. */
    
    double                 _zmin;          /*!< \brief Minimum z-value. */
    double                 _zmax;          /*!< \brief Maximum z-value. */

    double                 _datarange[4];  /*!< \brief Data ranges: xmin, ymin, xmax, ymax. */
    size_t                 _n;             /*!< \brief Size of data-array in x-direction. */
    size_t                 _m;             /*!< \brief Size of data-array in y-direction. */

    std::vector<double>    _f;             /*!< \brief Function value data, y major order. */

    Interpolation2D       *_intrp;         /*!< \brief Data interpolation. */

    bool                   _zscale_prepared;
    int                    _sign;
    double                 _scale_A;
    double                 _scale_B;
    double                 _scale_C;

    void prepare_zscaling( void );
    void prepare_data_interpolation( void );
    void plot_to_image_surface( cairo_surface_t *surface, const Coordmapper *cm, int plim[4] );

public:

    /*! \brief Default constructor for empty colormap graph.
     */
    MeshColormap();

    /*! \brief Copy constructor.
     */
    MeshColormap( const MeshColormap &colormap );

    /*! \brief Constructor for colormap from data.
     *
     *  Data is defined as \a n by \a m array of data, where x and y
     *  ranges are defined in datarange in order xmin, ymin, xmax,
     *  ymax. Z-values are defined in vector \a data in y major
     *  order. Internal copy of the data from data is made.
     */
    MeshColormap( const double datarange[4], size_t n, size_t m, 
		  const std::vector<double> &data );

    /*! \brief Destructor.
     */
    virtual ~MeshColormap();

    /*! \brief Clears colormap data.
     */
    void clear_data( void );
		   
    /*! \brief Define colormap from data.
     *
     *  Data is defined as \a n by \a m array of data, where x and y
     *  ranges are defined in datarange in order xmin, ymin, xmax,
     *  ymax. Z-values are defined in vector \a data in y major
     *  order. Internal copy of the data from data is made.
     *
     *  Overrides old data and resets z ranges.
     */
    void set_data( const double datarange[4], size_t n, size_t m, 
		   const std::vector<double> &data );

    /*! \brief Does colormap have data?
     */
    bool has_data( void ) const;

    /*! \brief Set interpolation mode.
     */
    interpolation_e get_interpolation( void ) const;

    /*! \brief Set interpolation mode.
     *
     *  Can be either \a INTERPOLATION_CLOSEST, \a
     *  INTERPOLATION_BILINEAR or \a INTERPOLATION_BICUBIC.
     */
    void set_interpolation( interpolation_e interpolation );

    /*! \brief Scale value \a val according to zscale mode.
     *  
     *  Scales value val, where \a zmin <= \a val <= \a zmax into the
     *  range between 0.0 and 1.0 according to zscale mode.
     */
    virtual double zscale( double val );

    /*! \brief Inverse scale value \a val according to zscale mode.
     *  
     *  Inverse function of zscale(). Scales value val, where \a 0 <= \a val <= \a 1 into the
     *  range between zmin and zmax according to zscale mode.
     */
    virtual double zscale_inv( double val );

    /*! \brief Get zscale mode.
     */
    zscale_e get_zscale( void ) const;

    /*! \brief Set zscale mode.
     *
     * Set a prescaling for z-axis. Defaults to \a ZSCALE_LINEAR,
     * which is a linear scaling from the z-axis to palette. Other
     * possibilities are \a ZSCALE_LOG, which is a standard
     * logarithmic scaling from z-axis to palette (providing
     * magnification close to zero) and \a ZSCALE_RELLOG, which is a
     * special relative logrithmic scaling following the relation \f[
     * \frac{\log(0.001+x)-\log(0.001)}{\log(1.001)-\log(0.001)} \f],
     * where \a x is prescaled to range [0,1]. The z-ranges completely
     * contained on the negative side are inverted to positive and
     * z-ranges both on negative and positive sides are scaled
     * separately to provide magnification close to zero.
     */
    void set_zscale( zscale_e zscale );

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

    /*! \brief Set colormap palette.
     */
    void set_palette( const Palette &palette );

    /*! \brief Get a reference to colormap palette.
     */
    virtual Palette &palette( void ) { return( _palette ); }

    /*! \brief Get a reference to colormap palette.
     */
    const Palette &palette( void ) const { return( _palette ); }

    /*! \brief Get zrange for colormap plot.
     */
    void get_zrange( double &min, double &max ) const;

    /*! \brief Set zrange for colormap plot.
     *
     *  The zrange defaults to automatically scaled ranging for
     *  colormap input data.
     */
    void set_zrange( double min, double max );

    /*! \brief Get value of interpolated colormap data.
     */
    double get_value( double x, double y ) const;
};


#endif
