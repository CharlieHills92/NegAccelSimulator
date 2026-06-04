/*! \file geom3dplot.hpp
 *  \brief %Geometry 3d plotter
 */

/* Copyright (c) 2012,2013,2022 Taneli Kalvas. All rights reserved.
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


#include <limits>
#include "geom3dplot.hpp"
#include "renderer.hpp"
#include "ibsimu.hpp"



Geom3DPlot::Geom3DPlot( const Geometry &geom,
			const ParticleDataBase *pdb )
    : _geom(geom), _pdb(pdb), _sdata_enable(false),
      _particle_div(100), _particle_offset(0), _bbox(true)

{
    _clevel[0] = 0;
    _clevel[1] = _geom.size(0)-1;
    _clevel[2] = 0;
    _clevel[3] = _geom.size(1)-1;
    _clevel[4] = 0;
    _clevel[5] = _geom.size(2)-1;

    // Default to plotting of all solids enabled
    for( uint32_t a = 0; a < geom.number_of_solids(); a++ )
	_senable.push_back( true );

    reset_camera_and_rotation();
    rebuild_model();
}


Geom3DPlot::~Geom3DPlot()
{

}


void Geom3DPlot::set_size( uint32_t width, uint32_t height )
{
    _width = width;
    _height = height;
}


void Geom3DPlot::set_model_transformation( const Transformation &modeltrans )
{
    _modeltrans = modeltrans;
}


Transformation Geom3DPlot::get_model_transformation( void ) const
{
    return( _modeltrans );
}


void Geom3DPlot::set_projection_frustum( double xnear, double xfar, double zoom )
{
    _near = xnear;
    _far  = xfar;
    _zoom = zoom;
}


void Geom3DPlot::get_projection_frustum( double &xnear, double &xfar, double &zoom ) const
{
    xnear = _near;
    xfar  = _far;
    zoom  = _zoom;
}


void Geom3DPlot::set_projection_zoom( double zoom )
{
    _zoom = zoom;
}


double Geom3DPlot::get_projection_zoom( void ) const
{
    return( _zoom );
}


void Geom3DPlot::set_view_look_at( const Vec3D &camera, 
				   const Vec3D &target,
				   const Vec3D &up )
{
    _camera = camera;
    _target = target;
    _up = up;
}


void Geom3DPlot::get_view_look_at( Vec3D &camera, 
				   Vec3D &target,
				   Vec3D &up ) const
{
    camera = _camera;
    target = _target;
    up = _up;
}


void Geom3DPlot::set_surface_triangle_color_range( double min, double max )
{
    _sdata_range[0] = min;
    _sdata_range[1] = max;
}


void Geom3DPlot::get_surface_triangle_color_range( double &min, double &max ) const
{
    min = _sdata_range[0];
    max = _sdata_range[1];
}


void Geom3DPlot::set_surface_triangle_data( const std::vector<double> *data )
{
    if( data == NULL ) {
	_sdata.clear();
	return;
    }
    _sdata = *data;

    // Autorange
    _sdata_range[0] = std::numeric_limits<double>::infinity();
    _sdata_range[1] = -std::numeric_limits<double>::infinity();
    for( uint32_t a = 0; a < _sdata.size(); a++ ) {
	if( _sdata[a] < _sdata_range[0] )
	    _sdata_range[0] = _sdata[a];
	if( _sdata[a] > _sdata_range[1] )
	    _sdata_range[1] = _sdata[a];
    }
    ibsimu.message(1) << "Surface data from " << _sdata_range[0] << " to " 
		      << _sdata_range[1] << "\n";

    _sdata_palette.clear();
    _sdata_palette.push_back( Vec3D(0,0,1), 0 );
    _sdata_palette.push_back( Vec3D(1,1,0), 1 );
    _sdata_palette.push_back( Vec3D(1,0,0), 2 );
    _sdata_palette.push_back( Vec3D(0,0,0), 3 );
    _sdata_palette.normalize();

    rebuild_model();
}


void Geom3DPlot::set_particle_div( uint32_t particle_div, uint32_t particle_offset )
{
    _particle_div = particle_div;
    _particle_offset = particle_offset;
}
    

uint32_t Geom3DPlot::get_particle_div( void ) const
{
    return( _particle_div );
}


uint32_t Geom3DPlot::get_particle_offset( void ) const
{
    return( _particle_offset );
}


void Geom3DPlot::set_bbox( bool bbox )
{
    _bbox = bbox;
}


bool Geom3DPlot::get_bbox( void ) const
{
    return( _bbox );
}


void Geom3DPlot::set_clevel( uint32_t direction, uint32_t level )
{
    if( direction > 5 )
	throw( ErrorRange( ERROR_LOCATION, "incorrect cut direction" ) );
    if( level >= _geom.size(direction/2) )
	throw( ErrorRange( ERROR_LOCATION, level, _geom.size(direction/2), "incorrect cut level" ) );
    _clevel[direction] = level;
}


uint32_t Geom3DPlot::get_clevel( uint32_t direction ) const
{
    if( direction > 5 )
	throw( ErrorRange( ERROR_LOCATION, "incorrect cut direction" ) );
    return( _clevel[direction] );
}


void Geom3DPlot::set_solid_plot( uint32_t a, bool enable )
{
    if( a < 7 )
	throw( ErrorRange( ERROR_LOCATION, "invalid solid number" ) );
    _senable[a-7] = enable;
}


bool Geom3DPlot::get_solid_plot( uint32_t a ) const
{
    if( a < 7 )
	throw( ErrorRange( ERROR_LOCATION, "invalid solid number" ) );
    return( _senable[a-7] );
}


void Geom3DPlot::set_surface_triangle_data_plot( bool enable )
{
    _sdata_enable = enable;
}


bool Geom3DPlot::get_surface_triangle_data_plot( void ) const
{
    return( _sdata_enable );
}


bool Geom3DPlot::node_disabled( uint32_t i, uint32_t j, uint32_t k ) const
{
    uint32_t smesh = _geom.mesh( i, j, k );
    if( (smesh & SMESH_NODE_ID_MASK) != SMESH_NODE_ID_DIRICHLET ||
	(smesh & SMESH_BOUNDARY_NUMBER_MASK) < 7 )
	return( false );
    uint32_t solid = smesh & SMESH_BOUNDARY_NUMBER_MASK;
    return( !_senable[solid-7] );
}


// Check if any of the cell nodes belongs to disabled solids
bool Geom3DPlot::cell_disabled( uint32_t i, uint32_t j, uint32_t k ) const
{
    if( node_disabled( i,   j,   k   ) ||
	node_disabled( i+1, j,   k   ) ||
	node_disabled( i,   j+1, k   ) ||
	node_disabled( i+1, j+1, k   ) ||
	node_disabled( i,   j,   k+1 ) ||
	node_disabled( i+1, j,   k+1 ) ||
	node_disabled( i,   j+1, k+1 ) ||
	node_disabled( i+1, j+1, k+1 ) )
	return( true );
    return( false );
}


bool Geom3DPlot::face_disabled( const int i[3], const int vb[3] ) const
{
    int32_t j[3] = { i[0], i[1], i[2] };

    // Node 1 (x,y)
    if( node_disabled( j[0], j[1], j[2] ) )
	return( true );

    // Node 2 (x+1,y)
    j[vb[0]]++;
    if( node_disabled( j[0], j[1], j[2] ) )
	return( true );

    // Node 3 (x+1,y+1)
    j[vb[1]]++;
    if( node_disabled( j[0], j[1], j[2] ) )
	return( true );

    // Node 4 (x,y+1)
    j[vb[0]]--;
    if( node_disabled( j[0], j[1], j[2] ) )
	return( true );

    return( false );
}


void Geom3DPlot::build_geometry_surface( void )
{
    // Reserve space
    _gsurface.reserve( 4*_geom.surface_trianglec() );
    if( _sdata.size() != 0 && _sdata_enable ) {
	if( _sdata.size() != _geom.surface_vertexc() )
	    throw( Error( ERROR_LOCATION, "surface data count does not match vertex count" ) );
	_gcolor.reserve( _geom.surface_vertexc() );    
    }

    // Go through all grid cells inside cut levels
    for( uint32_t k = _clevel[4]; k < _clevel[5]; k++ ) {
	for( uint32_t j = _clevel[2]; j < _clevel[3]; j++ ) {
	    for( uint32_t i = _clevel[0]; i < _clevel[1]; i++ ) {

		// Skip grid cell if it contains solid nodes, which are not enabled
		if( cell_disabled( i, j, k ) )
		    continue;

		// Copy surface triangles inside the cell
		int32_t ptr = _geom.surface_triangle_ptr( i, j, k );
		int32_t tric = _geom.surface_trianglec( i, j, k );
		for( int32_t a = 0; a < tric; a++ ) {
		    const VTriangle &tri = _geom.surface_triangle( ptr+a );
		    const Vec3D &x0 = _geom.surface_vertex( tri[0] );
		    const Vec3D &x1 = _geom.surface_vertex( tri[1] );
		    const Vec3D &x2 = _geom.surface_vertex( tri[2] );
		    Vec3D norm = cross( x1-x0, x2-x0 );
		    norm.normalize();
		    _gsurface.push_back( norm );
		    _gsurface.push_back( _scale*(x0-_center) );
		    _gsurface.push_back( _scale*(x1-_center) );
		    _gsurface.push_back( _scale*(x2-_center) );
		    if( _sdata.size() != 0 && _sdata_enable ) {
			_gcolor.push_back( sdata_palette( _sdata[tri[0]] ) );
			_gcolor.push_back( sdata_palette( _sdata[tri[1]] ) );
			_gcolor.push_back( sdata_palette( _sdata[tri[2]] ) );
		    }
		}
	    }
	}
    }
}


Vec3D Geom3DPlot::sdata_palette( double sval ) const
{
    double vscaled = (sval-_sdata_range[0]) / 
	(_sdata_range[1]-_sdata_range[0]);
    return( _sdata_palette(vscaled) );
}


void Geom3DPlot::rebuild_model( void )
{
    clear_surface_data();
    build_geometry_surface();
    build_cut_planes();
}


void Geom3DPlot::reset_camera_and_rotation( void )
{
    // Calculate scale
    _center = 0.5*(_geom.origo() + _geom.max());
    Vec3D size = _geom.max() - _geom.origo();
    _scale = 1.0/size.norm2();

    // Init camera
    _near = 1.0;
    _far = 3.0;
    _camera = Vec3D(0,0,2);
    _target = Vec3D(0,0,0);
    _up = Vec3D(0,1,0);
    _zoom = 0.25;

    // Init model transformation
    _modeltrans = Transformation::unity();
}


void Geom3DPlot::draw_bbox( Renderer *r )
{
    if( !_bbox )
	return;

    r->set_color( Vec3D(0,0,0) );
    r->disable_lighting();

    Vec3D pad = 0.01*_geom.h()*Vec3D(1,1,1);
    Vec3D min = _scale*(_geom.origo()-pad-_center);
    Vec3D max = _scale*(_geom.max()+pad-_center);

    Vec3D x0 = min;
    Vec3D x1 = Vec3D(  );
    r->line( Vec3D( min[0], min[1], min[2] ),
		     Vec3D( min[0], max[1], min[2] ) );
    r->line( Vec3D( min[0], max[1], min[2] ),
		     Vec3D( max[0], max[1], min[2] ) );
    r->line( Vec3D( max[0], max[1], min[2] ),
		     Vec3D( max[0], min[1], min[2] ) );
    r->line( Vec3D( max[0], min[1], min[2] ),
		     Vec3D( min[0], min[1], min[2] ) );

    r->line( Vec3D( min[0], min[1], max[2] ),
		     Vec3D( min[0], max[1], max[2] ) );
    r->line( Vec3D( min[0], max[1], max[2] ),
		     Vec3D( max[0], max[1], max[2] ) );
    r->line( Vec3D( max[0], max[1], max[2] ),
		     Vec3D( max[0], min[1], max[2] ) );
    r->line( Vec3D( max[0], min[1], max[2] ),
		     Vec3D( min[0], min[1], max[2] ) );

    r->line( Vec3D( min[0], min[1], min[2] ),
		     Vec3D( min[0], min[1], max[2] ) );
    r->line( Vec3D( min[0], max[1], min[2] ),
		     Vec3D( min[0], max[1], max[2] ) );
    r->line( Vec3D( max[0], max[1], min[2] ),
		     Vec3D( max[0], max[1], max[2] ) );
    r->line( Vec3D( max[0], min[1], min[2] ),
		     Vec3D( max[0], min[1], max[2] ) );
}


void Geom3DPlot::draw_model( Renderer *r )
{
    if( _sdata.size() != 0 && _sdata_enable ) {
	for( size_t a = 0, b = 0; a < _gsurface.size(); a+=4, b+=3 ) {
	    //r->set_material_ambient_color( 0.2*_gcolor[a] );
	    //r->set_material_diffuse_color( 0.8*_gcolor[a] );
	    //r->flat_triangle( _gsurface[b+1], _gsurface[b+2], _gsurface[b+3], _gsurface[b+0] );
	    r->shaded_triangle( _gsurface[a+1], _gcolor[b+0],
				_gsurface[a+2], _gcolor[b+1],
				_gsurface[a+3], _gcolor[b+2],
				_gsurface[a+0] );
	}
    } else {
	for( size_t a = 0; a < _gsurface.size(); a+=4 ) {
	    r->flat_triangle( _gsurface[a+1], _gsurface[a+2], _gsurface[a+3], _gsurface[a+0] );
	}
    }
}


void Geom3DPlot::draw_cut_planes( Renderer *r )
{
    for( int32_t p = 0; p < 6; p++ ) {
	
	Vec3D norm( 0.0, 0.0, 0.0 );
	if( p == 0 ) norm[0] = 1.0;
	else if( p == 1 ) norm[0] = -1.0;
	else if( p == 2 ) norm[1] = 1.0;
	else if( p == 3 ) norm[1] = -1.0;
	else if( p == 4 ) norm[2] = 1.0;
	else norm[2] = -1.0;

	for( uint32_t a = 0; a < _csurface[p].size(); a+=3 )
	    r->flat_triangle( _csurface[p][a+0], _csurface[p][a+1], _csurface[p][a+2], norm );
    }
}


int32_t Geom3DPlot::case2d( const int i[3], const int vb[3] )
{
    int res = 0;
    uint32_t node;
    int32_t j[3] = { i[0], i[1], i[2] };

    // Node 1 (x,y)
    node = _geom.mesh((j[2]*_geom.size(1) + j[1])*_geom.size(0) + j[0]);
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += 1;

    // Node 2 (x+1,y)
    j[vb[0]]++;
    node = _geom.mesh((j[2]*_geom.size(1) + j[1])*_geom.size(0) + j[0]);
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += 2;

    // Node 3 (x+1,y+1)
    j[vb[1]]++;
    node = _geom.mesh((j[2]*_geom.size(1) + j[1])*_geom.size(0) + j[0]);
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += 4;

    // Node 4 (x,y+1)
    j[vb[0]]--;
    node = _geom.mesh((j[2]*_geom.size(1) + j[1])*_geom.size(0) + j[0]);
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += 8;

    return( res );
}


void Geom3DPlot::cplane_add_vertex( int32_t p, const int32_t i[3], const int32_t vb[3], double dx, double dy )
{
    double u[3] = { (double)i[0], (double)i[1], (double)i[2] };
    u[vb[0]] += dx;
    u[vb[1]] += dy;

    Vec3D x( _geom.origo(0)+_geom.h()*u[0], _geom.origo(1)+_geom.h()*u[1], _geom.origo(2)+_geom.h()*u[2] );
    _csurface[p].push_back( _scale*(x-_center) );
}


double Geom3DPlot::cplane_dist( const int32_t i[3], const int32_t vb[3], int32_t dx, int32_t dy, int32_t dir )
{
    int32_t j[3] = { i[0], i[1], i[2] };
    j[vb[0]] += dx;
    j[vb[1]] += dy;
    return( _geom.solid_dist( j[0], j[1], j[2], dir )/255.0 );
}


void Geom3DPlot::build_cut_plane( int32_t p, int32_t const vb[3], int32_t level )
{
    int32_t i[3];
    i[vb[2]] = level;

    // Go through mesh squares in cut plane
    int32_t minx = _clevel[2*vb[0]];
    int32_t maxx = _clevel[2*vb[0]+1];
    int32_t miny = _clevel[2*vb[1]];
    int32_t maxy = _clevel[2*vb[1]+1];
    for( int32_t y = miny; y < maxy; y++ ) {
	i[vb[1]] = y;

	for( int32_t x = minx; x < maxx; x++ ) {
	    i[vb[0]] = x;

	    if( face_disabled( i, vb ) )
		continue;

	    // Construct case number
	    int32_t cn = case2d( i, vb );
	    //std::cout << cn << "\n";

	    double dist;
	    double dist2;
	    switch( cn ) {
	    case 0:
		// Fully outside
		break;
	    case 1:
		// (i,j) in
		dist = cplane_dist( i, vb, 1, 0, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 0.0 );
		dist = cplane_dist( i, vb, 0, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 0.0, 1.0-dist );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		break;
	    case 2:
		// (i+1,j) in
		dist = cplane_dist( i, vb, 0, 0, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 1.0, 1.0-dist );
		break;
	    case 3:
		// (i,j) and (i+1,j) in
		dist = cplane_dist( i, vb, 0, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 0.0, 1.0-dist );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0-dist );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 1.0, 1.0-dist );
		break;
	    case 4:
		// (i+1,j+1) in
		dist = cplane_dist( i, vb, 0, 1, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist, 1.0 );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 1.0, dist );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		break;
	    case 5:
		// (i,j) and (i+1,j+1) in
		dist = cplane_dist( i, vb, 0, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 0.0, 1.0-dist );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 0.0 );
		dist = cplane_dist( i, vb, 0, 1, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist, 1.0 );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 1.0, dist );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		break;
	    case 6:
		// (i+1,j) and (i+1,j+1) in
		dist = cplane_dist( i, vb, 0, 1, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist, 1.0 );
		dist = cplane_dist( i, vb, 0, 0, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, dist, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		break;
	    case 7:
		// (i,j), (i+1,j) and (i+1,j+1) in
		dist = cplane_dist( i, vb, 0, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 0.0, 1.0-dist );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		dist2 = cplane_dist( i, vb, 0, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, dist2, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0-dist );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, dist2, 1.0 );
		break;
	    case 8:
		dist = cplane_dist( i, vb, 0, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 0.0, dist );
		dist = cplane_dist( i, vb, 1, 1, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		break;
	    case 9:
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 1, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 1.0 );
		cplane_add_vertex( p, i, vb, 1.0-dist, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 0.0 );
		break;
	    case 10:
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		dist = cplane_dist( i, vb, 0, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 0.0, dist );
		dist = cplane_dist( i, vb, 1, 1, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 1.0 );
		dist = cplane_dist( i, vb, 0, 0, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 1.0, 1.0-dist );
		break;
	    case 11:
		dist = cplane_dist( i, vb, 1, 1, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		dist2 = cplane_dist( i, vb, 1, 1, 2*vb[1] );
		cplane_add_vertex( p, i, vb, 1.0, 1.0-dist2 );
		cplane_add_vertex( p, i, vb, 1.0-dist, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0-dist2 );
		break;
	    case 12:
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		dist = cplane_dist( i, vb, 0, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 0.0, dist );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 1.0, dist );
		cplane_add_vertex( p, i, vb, 1.0, dist );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		break;
	    case 13:
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		dist = cplane_dist( i, vb, 1, 0, 2*vb[0] );
		cplane_add_vertex( p, i, vb, 1.0-dist, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0-dist, 0.0 );
		dist2 = cplane_dist( i, vb, 1, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 1.0, dist2 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 1.0, dist2 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		break;
	    case 14:
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		dist = cplane_dist( i, vb, 0, 0, 2*vb[1]+1 );
		cplane_add_vertex( p, i, vb, 0.0, dist );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, dist );
		dist2 = cplane_dist( i, vb, 0, 0, 2*vb[0]+1 );
		cplane_add_vertex( p, i, vb, dist2, 0.0 );
		cplane_add_vertex( p, i, vb, dist2, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		break;
	    case 15:
		// Fully inside
		cplane_add_vertex( p, i, vb, 0.0, 0.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 1.0, 1.0 );
		cplane_add_vertex( p, i, vb, 0.0, 1.0 );
		cplane_add_vertex( p, i, vb, 1.0, 0.0 );
		break;
	    }
	    
	}
    }
}


void Geom3DPlot::build_cut_planes( void )
{
    int vb[3];

    // X=0
    vb[0] = 1;
    vb[1] = 2;
    vb[2] = 0;
    build_cut_plane( 0, vb, _clevel[0] );

    // X=size(0)-1
    vb[0] = 2;
    vb[1] = 1;
    vb[2] = 0;
    build_cut_plane( 1, vb, _clevel[1] );

    // Y=0
    vb[0] = 2;
    vb[1] = 0;
    vb[2] = 1;
    build_cut_plane( 2, vb, _clevel[2] );

    // Y=size(1)-1
    vb[0] = 0;
    vb[1] = 2;
    vb[2] = 1;
    build_cut_plane( 3, vb, _clevel[3] );

    // Z=0
    vb[0] = 0;
    vb[1] = 1;
    vb[2] = 2;
    build_cut_plane( 4, vb, _clevel[4] );

    // Z=size(2)-1
    vb[0] = 1;
    vb[1] = 0;
    vb[2] = 2;
    build_cut_plane( 5, vb, _clevel[5] );
}


void Geom3DPlot::draw_beam( Renderer *r )
{
    if( !_pdb || _particle_div == 0 )
	return;

    r->set_color( Vec3D(1,0,0) );
    r->disable_lighting();

    // Loop through all particles
    for( size_t a = _particle_offset; a < _pdb->size(); a += _particle_div ) {

	// pdiv plotting if one or less trajectory points
	if( _pdb->traj_size( a ) <= 1 )
	    continue;

	// Loop through all particle trajectory points
	double t;
	Vec3D ox, x, v;
	_pdb->trajectory_point( t, ox, v, a, 0 );
	for( size_t b = 1; b < _pdb->traj_size( a ); b++ ) {
	    
	    _pdb->trajectory_point( t, x, v, a, b );
	    r->line( _scale*(ox-_center), 
			     _scale*(x-_center) );
	    ox = x;
	}	
    }
}


void Geom3DPlot::clear_surface_data( void )
{
    _gsurface.clear();
    _gcolor.clear();
    for( int a = 0; a < 6; a++ )
	_csurface[a].clear();
}


void Geom3DPlot::draw( cairo_t *cairo, Renderer *r )
{
    r->start_rendering();

    // View
    r->set_projection_frustum( -_zoom*_width/_height, _zoom*_width/_height,
			       -_zoom, _zoom,
			       _near, _far );
    r->set_view_look_at( _camera, _target, _up );

    r->set_light_location( Vec3D(0,0,-1) );
    r->set_light_diffuse_color( Vec3D(1,1,1) );
    r->set_light_ambient_color( Vec3D(1,1,1) );

    r->set_material_ambient_color( Vec3D(0.0,0.0,0.2) );
    r->set_material_diffuse_color( Vec3D(0.0,0.0,0.8) );

    r->set_model_transformation( _modeltrans );

    r->enable_view_settings();

    // Draw model
    r->enable_lighting();
    draw_cut_planes( r );
    draw_model( r );
    draw_beam( r );
    draw_bbox( r );

    r->end_rendering( cairo );
}

