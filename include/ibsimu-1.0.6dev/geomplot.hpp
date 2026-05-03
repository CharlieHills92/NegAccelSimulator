/*! \file geomplot.hpp
 *  \brief %Geometry plotting
 */

/* Copyright (c) 2005-2012,2017 Taneli Kalvas. All rights reserved.
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

#ifndef GEOMPLOT_HPP
#define GEOMPLOT_HPP 1


#include "types.hpp"
#include "frame.hpp"
#include "geometry.hpp"
#include "meshscalarfield.hpp"
#include "vectorfield.hpp"
#include "epot_field.hpp"
#include "particledatabase.hpp"
#include "solidgraph.hpp"
#include "eqpotgraph.hpp"
#include "particlegraph.hpp"
#include "meshgraph.hpp"
#include "fieldgraph.hpp"




/*! \brief %Geometry plotter class.
 *
 *  Collection of graphs for building Geometry plots containing
 *  solids, equipotential lines, space charge field, particle
 *  trajectories and mesh lines. Uses Fieldgraph, SolidGraph,
 *  EqPotGraph, ParticleGraph and MeshGraph for plotting.
 */
class GeomPlot {

    Frame                   *_frame;

    const Geometry          &_geom;
    const EpotField         *_epot;
    const MeshScalarField   *_scharge;
    const MeshScalarField   *_tdens;
    const VectorField       *_bfield;
    const VectorField       *_efield;
    const ParticleDataBase  *_pdb;

    FieldGraph              *_fieldgraph;

    SolidGraph              *_solidgraph;
    EqPotGraph              *_eqpotgraph;
    ParticleGraph           *_particlegraph;
    MeshGraph               *_meshgraph;

    view_e                  _view;
    int                     _level;
    int                     _vb[3];

    uint32_t                _eqlines_auto;
    std::vector<double>     _eqlines_manual;
    uint32_t                _particle_div;
    uint32_t                _particle_offset;
    bool                    _scharge_field;
    bool                    _qm_discretation;
    bool                    _mesh;

    bool                    _cache;

    void reset_graphs( void );

public:

    /*! \brief Constructor for new geometry plot
     *
     *  Builds a new geometry plot in the plot frame. Default graph
     *  (SolidGraph) is added to the plot and view is set as
     *  XY-view. The default plane of view is the midplane for 3D
     *  geometries and 0 for others.
     */
    GeomPlot( Frame &frame, const Geometry &geom );

    /*! \brief Destructor for geometry plotter.
     */
    ~GeomPlot();

    /*! \brief Rebuild plot.
     */
    void build_plot( void );

    /*! \brief Disable plotting caches from use.
     *
     *  Used by the interactive plotter.
     */
    void disable_cache( void );

    /*! \brief Set electric potential.
     */
    void set_epot( const EpotField *epot );

    /*! \brief Set a vector of manual equipotential lines.
     */
    void set_eqlines_manual( const std::vector<double> &pot );

    /*! \brief Get a vector of manual equipotential lines.
     */
    std::vector<double> get_eqlines_manual( void ) const {
	return( _eqlines_manual );
    }

    /*! \brief Set the number of automatic equipotential lines.
     *
     *  The automatic lines are plotted at potentials pot = min +
     *  (i+0.5)*(max-min)/N ), where i is running from 0 to N and min
     *  is the minimum potential in the system and max is the maximum
     *  potential in the system.
     */
    void set_eqlines_auto( uint32_t N );

    /*! \brief Get the number of automatic equipotential lines.
     */
    uint32_t get_eqlines_auto( void ) const {
	return( _eqlines_auto );
    }

    /*! \brief Set magnetic field.
     */
    void set_bfield( const VectorField *bfield );

    /*! \brief Get magnetic field.
     */
    const VectorField *get_bfield( void ) const {
	return( _bfield );
    }

    /*! \brief Set electric field.
     */
    void set_efield( const VectorField *efield );

    /*! \brief Get electric field.
     */
    const VectorField *get_efield( void ) const {
	return( _efield );
    }

    /*! \brief Set trajectory density field.
     */
    void set_trajdens( const MeshScalarField *tdens );

    /*! \brief Get trajectory density field.
     */
    const MeshScalarField *get_trajdens( void ) const {
	return( _tdens );
    }

    /*! \brief Set space charge density field.
     */
    void set_scharge( const MeshScalarField *scharge );

    /*! \brief Get space charge density field.
     */
    const MeshScalarField *get_scharge( void ) const {
	return( _scharge );
    }

    /*! \brief Set field graph plotting type.
     */
    void set_fieldgraph_plot( field_type_e fieldplot );

    /*! \brief Get field graph object.
     */
    const FieldGraph *fieldgraph( void ) const {
	return( _fieldgraph );
    }

    /*! \brief Get field graph object.
     */
    FieldGraph *fieldgraph( void ) {
	return( _fieldgraph );
    }

    /*! \brief Set colormap legend enable/disable.
     */
    void enable_colormap_legend( bool enable ) {
	_frame->enable_colormap_legend( enable );
    }

    /*! \brief Set particle database used for particle plotting.
     */
    void set_particle_database( const ParticleDataBase *pdb ) {
	set_particledatabase( pdb );
    }

    /*! \brief Set particle database used for particle plotting.
     */
    void set_particledatabase( const ParticleDataBase *pdb );

    /*! \brief Set particle divisor and offset.
     *
     *  Set \a particle_div to zero for no plotting, one for plotting
     *  every particle, two for plotting every second particle, three
     *  for plotting every third particle, etc. Defaults to 11.
     *  %Plotter skips the first \a particle_offset particles.
     *
     *  This function controls plotting of particle trajectories.
     */
    void set_particle_div( uint32_t particle_div, uint32_t particle_offset = 0 );

    /*! \brief Get particle divisor.
     */
    uint32_t get_particle_div( void ) const {
	return( _particle_div );
    }

    /*! \brief Get particle offset.
     */
    uint32_t get_particle_offset( void ) const {
	return( _particle_offset );
    }

    /*! \brief Set q/m particle discretation.
     *
     *  If enabled, the different q/m values will be plotted with
     *  different colors. Otherwise, all particles are plotted with
     *  same color.
     */
    void set_qm_discretation( bool enable );

    /*! \brief Get q/m particle discretation.
     */
    bool get_qm_discretation( void ) const {
	return( _qm_discretation );
    }

    /*! \brief Set mesh plotting.
     *
     *  If enabled, the mesh squares are plotted. Mesh plotting is
     *  disabled by default.
     */
    void set_mesh( bool enable );

    /*! \brief Get mesh plotting.
     */
    bool get_mesh( void ) const {
	return( _mesh );
    }

    /*! \brief Set view.
     *
     *  Sets the viewplane to the geometry. The viewplane is set by
     *  direction \a view and depth \a level set as mesh level. Level
     *  is checked and limited to existing levels. Level -1 (default)
     *  means half the range (midplane).
     */
    void set_view( view_e view, int level = -1 );

    /*! \brief Set view in SI units.
     *
     *  Sets the viewplane to the geometry. The viewplane is set by
     *  direction \a view and depth \a level set as coordinates in SI
     *  units. The level is limited to existing values.
     */
    void set_view_si( view_e view, double level );

    /*! \brief Get view.
     */
    view_e get_view( void ) const {
	return( _view );
    }

    /*! \brief Get level of view in mesh squares.
     */
    int get_level( void ) const {
	return( _level );
    }

    /*! \brief Get level of view in SI units.
     */
    double get_level_si( void ) const {
	return( _geom.origo(_vb[2])+_level*_geom.h() );
    }

    /*! \brief Get component \a i of view base vector.
     */
    int vb( int i ) const {
	return( _vb[i] );
    }

    /*! \brief Get the view base vector.
     */
    void get_vb( int vb[3] ) const {
	vb[0] = _vb[0];
	vb[1] = _vb[1];
	vb[2] = _vb[2];
    }

    /*! \brief Get a reference to the frame
     */
    Frame &frame( void ) {
	return( *_frame );
    }
};


#endif

