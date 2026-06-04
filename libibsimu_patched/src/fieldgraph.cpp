/*! \file fieldgraph.cpp
 *  \brief %Graph for plotting fields
 */

/* Copyright (c) 2005-2014 Taneli Kalvas. All rights reserved.
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
#include "fieldgraph.hpp"
#include "scalarfield.hpp"
#include "vectorfield.hpp"
#include "meshscalarfield.hpp"
#include "meshvectorfield.hpp"
#include "epot_efield.hpp"
#include "ibsimu.hpp"


//#define DEBUG_FIELDGRAPH 1


FieldGraph::FieldGraph( const Geometry &geom, const Field *field, field_type_e field_type )
    : Graph3D(geom), _geom(geom), _field(field), _field_type(field_type),
      _first(true), _oview(VIEW_XY), _olevel(0), _enabled(true), _zmin(0.0), _zmax(0.0)
{
}


FieldGraph::~FieldGraph()
{
}


field_type_e FieldGraph::field_type( void ) const
{
    return( _field_type );
}

/*
void FieldGraph::set_field( field_type_e field_type, const ScalarField *field )
{
    // Set parameters
    _first = true;
    _field_type = field_type;
    _scalarfield = field;
    _vectorfield = NULL;

    if( field == NULL ) {
	_field_type = FIELD_NONE;
	return;
    }

    _scalarfield->get_minmax( _zmin, _zmax );

    // Set default colormap palette
    std::vector<Palette::Entry> pentry;
    double zspan = _zmax - _zmin;
    if( _zmin >= -1.0e-6*zspan && _zmax >= 0.0 ) {
	// Red palette
	pentry.push_back( Palette::Entry( Vec3D(1,1,1), 0 ) );
	pentry.push_back( Palette::Entry( Vec3D(1,1,0), 1 ) );
	pentry.push_back( Palette::Entry( Vec3D(1,0,0), 2 ) );
	pentry.push_back( Palette::Entry( Vec3D(0,0,0), 3 ) );
    } else if( _zmax <= 1.0e-6*zspan && _zmin <= 0.0 ) {
	// Red palette
	pentry.push_back( Palette::Entry( Vec3D(1,1,1), 3 ) );
	pentry.push_back( Palette::Entry( Vec3D(1,1,0), 2 ) );
	pentry.push_back( Palette::Entry( Vec3D(1,0,0), 1 ) );
	pentry.push_back( Palette::Entry( Vec3D(0,0,0), 0 ) );
    } else {
	// Mixed palette, forcing zero at white
	pentry.push_back( Palette::Entry( Vec3D(0,0,0), _zmin ) );
	pentry.push_back( Palette::Entry( Vec3D(0,0,1), 0.67*_zmin ) );
	pentry.push_back( Palette::Entry( Vec3D(0,1,1), 0.33*_zmin ) );
	pentry.push_back( Palette::Entry( Vec3D(1,1,1), 0 ) );
	pentry.push_back( Palette::Entry( Vec3D(1,1,0), 0.33*_zmax ) );
	pentry.push_back( Palette::Entry( Vec3D(1,0,0), 0.67*_zmax ) );
	pentry.push_back( Palette::Entry( Vec3D(0,0,0), _zmax ) );
    }
    Palette p( pentry );
    set_palette( p );
    MeshColormap::set_zrange( _zmin, _zmax );
}


void FieldGraph::set_field( field_type_e field_type, const VectorField *field )
{
    // Set parameters
    _first = true;
    _field_type = field_type;
    _scalarfield = NULL;
    _vectorfield = field;

    if( field == NULL ) {
	_field_type = FIELD_NONE;
	return;
    }

    const MeshVectorField *mvf = dynamic_cast<const MeshVectorField *>( field );
    if( mvf ) {

	Vec3D min, max;
	switch( _field_type ) {
	case FIELD_EFIELD:
	case FIELD_BFIELD:
	    mvf->get_minmax( _zmin, _zmax );
	    break;
	case FIELD_EFIELD_X:
	case FIELD_BFIELD_X:
	    mvf->get_minmax( min, max );
	    _zmin = min[0];
	    _zmax = max[0];
	    break;
	case FIELD_EFIELD_Y:
	case FIELD_BFIELD_Y:
	    mvf->get_minmax( min, max );
	    _zmin = min[1];
	    _zmax = max[1];
	    break;
	case FIELD_EFIELD_Z:
	case FIELD_BFIELD_Z:
	    mvf->get_minmax( min, max );
	    _zmin = min[2];
	    _zmax = max[2];
	    break;
	default:
	    throw( Error( ERROR_LOCATION, "unknown field type" ) );
	    break;
	}

    } else {

	Vec3D min, max;
	switch( _field_type ) {
	case FIELD_EFIELD:
	case FIELD_BFIELD:
	    field->get_minmax( _geom, _zmin, _zmax );
	    break;
	case FIELD_EFIELD_X:
	case FIELD_BFIELD_X:
	    field->get_minmax( _geom, min, max );
	    _zmin = min[0];
	    _zmax = max[0];
	    break;
	case FIELD_EFIELD_Y:
	case FIELD_BFIELD_Y:
	    field->get_minmax( _geom, min, max );
	    _zmin = min[1];
	    _zmax = max[1];
	    break;
	case FIELD_EFIELD_Z:
	case FIELD_BFIELD_Z:
	    field->get_minmax( _geom, min, max );
	    _zmin = min[2];
	    _zmax = max[2];
	    break;
	default:
	    throw( Error( ERROR_LOCATION, "unknown field type" ) );
	    break;
	}

    }

    // Set default colormap palette
    std::vector<Palette::Entry> pentry;
    pentry.push_back( Palette::Entry( Vec3D(0,0,0), -1.00 ) );
    pentry.push_back( Palette::Entry( Vec3D(0,0,1), -0.67 ) );
    pentry.push_back( Palette::Entry( Vec3D(0,1,1), -0.33 ) );
    pentry.push_back( Palette::Entry( Vec3D(1,1,1),  0.00 ) );
    pentry.push_back( Palette::Entry( Vec3D(1,1,0),  0.33 ) );
    pentry.push_back( Palette::Entry( Vec3D(1,0,0),  0.67 ) );
    pentry.push_back( Palette::Entry( Vec3D(0,0,0),  1.00 ) );
    Palette p( pentry );
    set_palette( p );
    MeshColormap::set_zrange( _zmin, _zmax );
}
*/

void FieldGraph::enable( bool enable )
{
    _enabled = enable;
}


void FieldGraph::set_zrange( double zmin, double zmax )
{
    _zmin = zmin;
    _zmax = zmax;

    MeshColormap::set_zrange( _zmin, _zmax );
}


void FieldGraph::build_plot( void )
{
#ifdef DEBUG_FIELDGRAPH
    std::cout << "FieldGraph::build_plot()\n";
#endif

    if( !_field )
	return;

    const VectorField *vfield = dynamic_cast<const VectorField *>( _field );
    const ScalarField *sfield = dynamic_cast<const ScalarField *>( _field );
    const MeshScalarField *msfield = dynamic_cast<const MeshScalarField *>( _field );
    const MeshVectorField *mvfield = dynamic_cast<const MeshVectorField *>( _field );
    const EpotEfield *epotefield = dynamic_cast<const EpotEfield *>( _field );

    // Range on the plane of colormap
    double datamin = std::numeric_limits<double>::infinity();
    double datamax = -std::numeric_limits<double>::infinity();
    // Global range from Mesh-type field (if available)
    double globalmin = std::numeric_limits<double>::infinity();
    double globalmax = -std::numeric_limits<double>::infinity();

    // Check field type existence and get ranges
    switch( _field_type ) {
    case FIELD_NONE:
#ifdef DEBUG_FIELDGRAPH
	std::cout << "FieldGraph::build_plot(): none field\n";
#endif
	return;
    case FIELD_EPOT:
    case FIELD_SCHARGE:
    case FIELD_TRAJDENS:
#ifdef DEBUG_FIELDGRAPH
	std::cout << "FieldGraph::build_plot(): scalar field\n";
#endif
	if( sfield == NULL )
	    throw( Error( ERROR_LOCATION, "not a scalar field" ) );
	if( msfield )
	    msfield->get_minmax( globalmin, globalmax );
	break;
    case FIELD_EFIELD:
    case FIELD_BFIELD:
#ifdef DEBUG_FIELDGRAPH
	std::cout << "FieldGraph::build_plot(): vector field magnitude\n";
#endif
	if( vfield == NULL )
	    throw( Error( ERROR_LOCATION, "not a vector field" ) );
	if( mvfield )
	    mvfield->get_minmax( globalmin, globalmax );
	break;
    case FIELD_EFIELD_X:
    case FIELD_BFIELD_X:
#ifdef DEBUG_FIELDGRAPH
	std::cout << "FieldGraph::build_plot(): vector field x\n";
#endif
	if( vfield == NULL )
	    throw( Error( ERROR_LOCATION, "not a vector field" ) );
	if( mvfield ) {
	    Vec3D min, max;
	    mvfield->get_minmax( min, max );
	    globalmin = min[0];
	    globalmax = max[0];
	}
	break;
    case FIELD_EFIELD_Y:
    case FIELD_BFIELD_Y:
#ifdef DEBUG_FIELDGRAPH
	std::cout << "FieldGraph::build_plot(): vector field y\n";
#endif
	if( vfield == NULL )
	    throw( Error( ERROR_LOCATION, "not a vector field" ) );
	if( mvfield ) {
	    Vec3D min, max;
	    mvfield->get_minmax( min, max );
	    globalmin = min[1];
	    globalmax = max[1];
	}
	break;
    case FIELD_EFIELD_Z:
    case FIELD_BFIELD_Z:
#ifdef DEBUG_FIELDGRAPH
	std::cout << "FieldGraph::build_plot(): vector field z\n";
#endif
	if( vfield == NULL )
	    throw( Error( ERROR_LOCATION, "not a vector field" ) );
	if( mvfield ) {
	    Vec3D min, max;
	    mvfield->get_minmax( min, max );
	    globalmin = min[2];
	    globalmax = max[2];
	}
	break;
    default:
	throw( Error( ERROR_LOCATION, "unknown field type" ) );
	break;
    }

    Mesh mesh = _geom;
    if( epotefield ) {
	double h = 0.5*mesh.h();
	Int3D one(1,1,1);
	Int3D size = 2*(mesh.size()-one)+one;
	mesh.reset( mesh.geom_mode(), size, mesh.origo(), h );
    }

    // Build data for colormap based on geometry mesh
    // Could use mesh of field itself in case of mesh-based fields
    // Vectorfield can have a transformation!!
    double range[4] = { mesh.origo( _vb[0] ),
			mesh.origo( _vb[1] ),
			mesh.max  ( _vb[0] ),
			mesh.max  ( _vb[1] ) };
    size_t n = mesh.size( _vb[0] );
    size_t m = mesh.size( _vb[1] );

    std::vector<double> data;
    data.reserve( n*m );
    Vec3D x;
    x[_vb[2]] = _geom.origo(_vb[2]) + _level*_geom.h();
    for( size_t j = 0; j < m; j++ ) {
	x[_vb[1]] = mesh.origo(_vb[1]) + j*mesh.h();
	for( size_t i = 0; i < n; i++ ) {
	    x[_vb[0]] = mesh.origo(_vb[0]) + i*mesh.h();

	    Vec3D F;
	    double dataent;

	    switch( _field_type ) {
	    case FIELD_NONE:
		break;
	    case FIELD_EPOT:
	    case FIELD_SCHARGE:
	    case FIELD_TRAJDENS:
		dataent = (*sfield)( x );
		break;
	    case FIELD_EFIELD:
	    case FIELD_BFIELD:
		F = (*vfield)( x );
		dataent = F.norm2();
		break;
	    case FIELD_EFIELD_X:
	    case FIELD_BFIELD_X:
		F = (*vfield)( x );
		dataent = F[0];
		break;
	    case FIELD_EFIELD_Y:
	    case FIELD_BFIELD_Y:
		F = (*vfield)( x );
		dataent = F[1];
		break;
	    case FIELD_EFIELD_Z:
	    case FIELD_BFIELD_Z:
		F = (*vfield)( x );
		dataent = F[2];
		break;
	    default:
		throw( Error( ERROR_LOCATION, "unknown field type" ) );
		break;
	    }

	    data.push_back( dataent );
	    if( dataent < datamin )
		datamin = dataent;
	    if( dataent > datamax )
		datamax = dataent;
	}
    }

    // Use global minimum and maximum if possible, otherwise use
    // (datamin, datamax), which is the local extremes in this level.
    if( msfield || mvfield ) {
	_zmin = globalmin;
	_zmax = globalmax;
    } else {
	_zmin = datamin;
	_zmax = datamax;
    }

    // Set up colormap
#ifdef DEBUG_FIELDGRAPH
    std::cout << "FieldGraph::build_plot(): set up colormap\n";
#endif
    set_data( range, n, m, data );
    set_zrange( _zmin, _zmax );

    // Set palette
    std::vector<Palette::Entry> pentry;
    double zspan = _zmax - _zmin;
    if( _zmin >= -1.0e-6*zspan && _zmax >= 0.0 ) {
	// "Positive" palette
	pentry.push_back( Palette::Entry( Vec3D(1,1,1), 0 ) );
	pentry.push_back( Palette::Entry( Vec3D(1,1,0), 1 ) );
	pentry.push_back( Palette::Entry( Vec3D(1,0,0), 2 ) );
	pentry.push_back( Palette::Entry( Vec3D(0,0,0), 3 ) );
    } else if( _zmax <= 1.0e-6*zspan && _zmin <= 0.0 ) {
	// "Negative" palette
	pentry.push_back( Palette::Entry( Vec3D(1,1,1), 3 ) );
	pentry.push_back( Palette::Entry( Vec3D(1,1,0), 2 ) );
	pentry.push_back( Palette::Entry( Vec3D(1,0,0), 1 ) );
	pentry.push_back( Palette::Entry( Vec3D(0,0,0), 0 ) );
    } else {
	// "Mixed" palette - red palette for now
	pentry.push_back( Palette::Entry( Vec3D(1,1,1), 0 ) );
	pentry.push_back( Palette::Entry( Vec3D(1,1,0), 1 ) );
	pentry.push_back( Palette::Entry( Vec3D(1,0,0), 2 ) );
	pentry.push_back( Palette::Entry( Vec3D(0,0,0), 3 ) );
	/*
	  // Does not work in practice
	pentry.push_back( Palette::Entry( Vec3D(0,0,0), _zmin ) );
	pentry.push_back( Palette::Entry( Vec3D(0,0,1), 0.67*_zmin ) );
	pentry.push_back( Palette::Entry( Vec3D(0,1,1), 0.33*_zmin ) );
	pentry.push_back( Palette::Entry( Vec3D(1,1,1), 0 ) );
	pentry.push_back( Palette::Entry( Vec3D(1,1,0), 0.33*_zmax ) );
	pentry.push_back( Palette::Entry( Vec3D(1,0,0), 0.67*_zmax ) );
	pentry.push_back( Palette::Entry( Vec3D(0,0,0), _zmax ) );
	*/
    }

    Palette p( pentry );
    set_palette( p );

    _first = false;
}


void FieldGraph::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] )
{
    if( _first || _oview != _view || _olevel != _level ) {
	// First plot or changed view happened
	build_plot();
    }

    _oview = _view;
    _olevel = _level;

    if( _field_type != FIELD_NONE && _enabled )
	MeshColormap::plot( cairo, cm, range );
}


void FieldGraph::plot_sample( cairo_t *cairo, double x, double y, double width, double height )
{

}


void FieldGraph::get_bbox( double bbox[4] )
{
    bbox[0] = _geom.origo( _vb[0] );
    bbox[1] = _geom.origo( _vb[1] );
    bbox[2] = _geom.max( _vb[0] );
    bbox[3] = _geom.max( _vb[1] );
}


