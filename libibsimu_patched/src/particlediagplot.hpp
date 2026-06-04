/*! \file particlediagplot.hpp
 *  \brief %Particle diagnostic plot
 */

/* Copyright (c) 2005-2011,2013-2014 Taneli Kalvas. All rights reserved.
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

#ifndef PARTICLEDIAGPLOT_HPP
#define PARTICLEDIAGPLOT_HPP 1


#include "frame.hpp"
#include "geometry.hpp"
#include "particledatabase.hpp"
#include "types.hpp"
#include "histogram.hpp"
#include "trajectorydiagnostics.hpp"

#include "xygraph.hpp"
#include "meshcolormap.hpp"



/*! \brief %Particle diagnostic plot type.
 */
enum particle_diag_plot_type_e {
    PARTICLE_DIAG_PLOT_NONE = 0, /*!< \brief Dummy diagnostic type. Does nothing. */
    PARTICLE_DIAG_PLOT_SCATTER,  /*!< \brief Scatter plot. */
    PARTICLE_DIAG_PLOT_HISTO2D,  /*!< \brief 2D histogram. */
    PARTICLE_DIAG_PLOT_HISTO1D   /*!< \brief 1D histogram. */
};


/*! \brief %Particle diagnostic plot.
 *
 *  Two dimensional histograms have particle trajectory currents
 *  always taken in account.  Similarly profile plot
 *  (PARTICLE_DIAG_PLOT_HISTO1D) always takes in account the
 *  trajectory current. In cylindrical symmetry cases the output is
 *  scaled to have constant area per histogram bin.  One dimensional
 *  Emittance plots in (r,r') space are scaled to have constant area
 *  per histogram bin.
 *  
 */
class ParticleDiagPlot {
    
    Frame                     &_frame;

    const Geometry            &_geom;
    const ParticleDataBase    &_pdb;

    bool                       _free_plane;

    coordinate_axis_e          _axis;
    double                     _level;

    Vec3D                      _c;
    Vec3D                      _o;
    Vec3D                      _p;

    particle_diag_plot_type_e  _type;
    trajectory_diagnostic_e    _diagx;
    trajectory_diagnostic_e    _diagy;
    trajectory_diagnostic_e    _diagz;

    int                        _pdb_it_no;
    bool                       _update;
    double                     _Isum;
    TrajectoryDiagnosticData  *_tdata;  /*!< \brief Trajectory data (scatter) */
    Histogram                 *_histo;  /*!< \brief Histogram data */
    Emittance                 *_emit;   /*!< \brief Emittance data */

    XYGraph                   *_datagraph;
    XYGraph                   *_ellipse;
    bool                       _ellipse_enable;

    MeshColormap              *_colormap;
    std::vector<double>        _zdata;

    size_t                     _histogram_n;
    size_t                     _histogram_m;
    histogram_accumulation_e   _histogram_accumulation;
    bool                       _histogram_style;
    interpolation_e            _interpolation;
    double                     _dot_size;
    
    void build_data( void );
    void merge_bbox( double bbox[4], const double bb[4] );
    
public:

    /*! \brief Constructor for particle diagnostic plot.
     *
     *  Make a diagnostic plot using plotter \a frame for making a
     *  diagnostic at plane \a axis = \a level using particle data
     *  from \a pdb in geometry \a geom. The particle diagnostic is
     *  defined by diagnostic \a type and diagnostic axes \a diagx and
     *  \a diagy.
     *
     *  The EmittanceConv class is automatically used to convert
     *  particle data from cylindrically symmetric simulations if
     *  (y,y') or (z,z') plot is requested with plot type
     *  PARTICLE_DIAG_PLOT_HISTO2D.
     */
    ParticleDiagPlot( Frame &frame, const Geometry &geom, const ParticleDataBase &pdb, 
		      coordinate_axis_e axis, double level, 
		      particle_diag_plot_type_e type,
		      trajectory_diagnostic_e diagx, trajectory_diagnostic_e diagy = DIAG_NONE );

    /*! \brief Constructor for particle diagnostic plot.
     *
     *  Make a diagnostic plot using plotter \a frame for making a
     *  diagnostic at plane defined by center point \a c and two
     *  vectors defining the coordinate axes \a o and \a p. The
     *  diagnostics is made using particle data from \a pdb in
     *  geometry \a geom. The particle diagnostic is defined by
     *  diagnostic \a type and diagnostic axes \a diagx and \a diagy.
     */
    ParticleDiagPlot( Frame &frame, const Geometry &geom, const ParticleDataBase &pdb, 
		      const Vec3D &c, const Vec3D &o, const Vec3D &p,
		      particle_diag_plot_type_e type,
		      trajectory_diagnostic_e diagx, trajectory_diagnostic_e diagy = DIAG_NONE );

    /*! \brief Destructor.
     */
    ~ParticleDiagPlot();

    /*! \brief Enable/disable emittance fit for emittance plots?
     */
    void set_emittance_ellipse( bool enable ) {
	_ellipse_enable = enable;
    }

    /*! \brief Is emittance fit enabled for emittance plots?
     */
    bool get_emittance_ellipse( void ) {
	return( _ellipse_enable );
    }

    /*! \brief Set diagnostic plane.
     */
    void set_view( coordinate_axis_e axis, double level ) {	
	_update = true;
	_axis = axis;
	_level = level;
    }

    /*! \brief Get diagnostic plane definition if it an even coorinate plane.
     */
    void get_view( coordinate_axis_e &axis, double &level ) {
	axis = _axis;
	level = _level;
    }

    /*! \brief Set plot type.
     */
    void set_type( particle_diag_plot_type_e type ) {
	_update = true;
	_type = type;
    }

    /*! \brief Get plot type.
     */
    particle_diag_plot_type_e get_type( void ) {
	return( _type );
    }

    /*! \brief Set plot type and diagnostic axes of plot.
     */
    void set_plot( particle_diag_plot_type_e type,
		   trajectory_diagnostic_e diagx, trajectory_diagnostic_e diagy ) {
	_update = true;
	_type = type;
	_diagx = diagx;
	_diagy = diagy;
    }

    /*! \brief Get plot type and diagnostic axes of plot.
     */
    void get_plot( particle_diag_plot_type_e &type,
		   trajectory_diagnostic_e &diagx, trajectory_diagnostic_e &diagy ) {
	type = _type;
	diagx = _diagx;
	diagy = _diagy;
    }

    /*! \brief Set number of histogram bins in x-direction to use for colormap plot.
     */
    void set_histogram_n( size_t n ) {
	_update = true;
	_histogram_n = n;
    }

    /*! \brief Get number of histogram bins in x-direction to use for colormap plot.
     */
    size_t get_histogram_n( void ) {
	return( _histogram_n );
    }

    /*! \brief Set number of histogram bins in y-direction to use for colormap plot.
     */
    void set_histogram_m( size_t m ) {
	_update = true;
	_histogram_m = m;
    }

    /*! \brief Get number of histogram bins in y-direction to use for colormap plot.
     */
    size_t get_histogram_m( void ) {
	return( _histogram_m );
    }

    /*! \brief Set histogram accumulation type.
     */
    void set_histogram_accumulation( histogram_accumulation_e accumulation ) {
	_update = true;
	_histogram_accumulation = accumulation;
    }

    /*! \brief Get histogram accumulation type.
     */
    histogram_accumulation_e get_histogram_accumulation( void ) {
	return( _histogram_accumulation );
    }

    /*! \brief Set 1d histogram style.
     */
    void set_histogram_style( bool style ) {
	_histogram_style = style;
    }

    /*! \brief Get 1d histogram style.
     */
    bool get_histogram_style( void ) {
	return( _histogram_style );
    }

    /*! \brief Set the type of interpolation used in colormap plot.
     */
    void set_colormap_interpolation( interpolation_e interpolation ) {
	_interpolation = interpolation;
	if( _colormap )
	    _colormap->set_interpolation( interpolation );
    }

    /*! \brief Get the type of interpolation used in colormap plot.
     */
    interpolation_e get_colormap_interpolation( void ) {
	return( _interpolation );
    }

    /*! \brief Get a pointer to histogram in the plot.
     */
    const MeshColormap *get_colormap( void ) const {
	return( _colormap );
    }

    /*! \brief Set dot size for scatter plot.
     */
    void set_dot_size( double size ) {
	_dot_size = size;
	if( _datagraph )
	    _datagraph->set_point_style( XYGRAPH_POINT_CIRCLE, true, _dot_size );
    }

    /*! \brief Get dot size for scatter plot.
     */
    double get_dot_size( void ) {
	return( _dot_size );
    }

    /*! \brief Return a pointer to histogram.
     *
     *  Histogram might not exist in the plot object.
     */
    const Histogram *get_histogram( void ) {
	return( _histo );
    }

    /*! \brief Return total current in diagnostic.
     *
     *  No mirroring applied to returned current.
     */
    double get_isum( void ) {
	return( _Isum );
    }

    /*! \brief Calculate Emittance fit.
     */
    const Emittance &calculate_emittance( void );

    /*! \brief Export plotted data as ASCII.
     *
     *  If the plot is a scatter plot, each particle is exported with
     *  plotted coordinates and current it is carrying. If the plot is
     *  a histogram a grid is outputted with coordinates and current
     *  density in grid cells.
     */
    void export_data( const std::string &filename );

    /*! \brief Rebuild plot.
     */
    void build_plot( void );
};


#endif

