/*! \file geometry.cpp
 *  \brief %Geometry definition.
 */

/* Copyright (c) 2005-2014,2019 Taneli Kalvas. All rights reserved.
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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string.h>
#include <limits>
#include "stlfile.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "dxf_solid.hpp"
#include "stl_solid.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "file.hpp"
#include "timer.hpp"


//#define DEBUG_BRACKET_SURFACE 1
//#define MC_DEBUG 1
#define DEBUG_SURFACE_TRIANGLE_NORMAL 1


Bound::Bound( bound_e type, double value ) 
    : _type(type), _value(value), _functor(NULL)
{

}


Bound::Bound( bound_e type, const CallbackFunctorD_V *functor )
    : _type(type), _value(0.0), _functor(functor)
{

}


bound_e Bound::type( void ) const
{
    return( _type );
}


void Bound::set_value( double value )
{
    _functor = NULL;
    _value = value;
}


double Bound::value( void ) const
{
    if( _functor )
	throw( Error( ERROR_LOCATION, "non-constant boundary value" ) );
    return( _value );
}


double Bound::value( const Vec3D &x ) const
{
    if( _functor )
	return( (*_functor)( x ) );
    return( _value );
}


bool Bound::is_constant() const
{
    return( !_functor );
}


Bound::Bound( std::istream &is ) 
{
    // Before version number 2, data was started by type 0 or 1.
    int32_t version = read_int32( is );
    if( version < 2 ) {
	_type = (bound_e)version;
	_value = read_double( is );
	_functor = NULL;
    } else {
	_type = (bound_e)read_int32( is );
	_value = read_double( is );
	if( read_int32( is ) )
	    ibsimu.message( MSG_WARNING, 1 ) << "Warning: functor boundary value was loaded.\n";
	_functor = NULL;
    }
}


void Bound::save( std::ostream &os ) const
{
    write_int32( os, 2 ); // version
    write_int32( os, _type );
    write_double( os, _value );
    if( _functor ) {
	ibsimu.message( MSG_WARNING, 1 ) << "Warning: functor boundary value was saved.\n";
	write_int32( os, 1 );
    } else {
	write_int32( os, 0 );
    }
}


std::ostream &operator<<( std::ostream &os, const Bound &b )
{
    if( b._type == BOUND_NEUMANN )
	os << "(BOUND_NEUMANN, ";
    else 
	os << "(BOUND_DIRICHLET, ";
    if( b._functor )
	os << "<functor>)";
    else
	os << b._value << ")";
    return( os );
}



Geometry::Geometry( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h )
    : Mesh(geom_mode,size,origo,h)
{
    check_definition();

    ibsimu.message( 1 ) << "Constructing geometry\n";
    ibsimu.inc_indent();
    ibsimu.message( 1 ) << "origo     = " << _origo << "\n";
    ibsimu.message( 1 ) << "size      = " << _size << "\n";
    ibsimu.message( 1 ) << "max       = " << _max << "\n";
    ibsimu.message( 1 ) << "h         = " << _h << "\n";
    ibsimu.message( 1 ) << "nodecount = " << nodecount() << "\n";

    _n = 0;
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );

    _built = false;
    _smesh = new uint32_t[_size[0]*_size[1]*_size[2]];

    pthread_mutex_init( &_mutex, NULL );
    pthread_cond_init( &_cond, NULL );

    ibsimu.message( 1 ) << "Done\n";
    ibsimu.dec_indent();
}


Geometry::Geometry( std::istream &is )
    : Mesh(is)
{
    check_definition();

    ibsimu.message( 1 ) << "Constructing Geometry from stream\n";
    ibsimu.inc_indent();

    _n = read_int32( is );
    for( uint32_t a = 0; a < _n; a++ ) {
	int32_t fileid = read_int32( is );
	if( fileid == FILEID_NULL )
	    _sdata.push_back( NULL );
	else if( fileid == FILEID_FUNCSOLID )
	    _sdata.push_back( new FuncSolid( is ) );
	else if( fileid == FILEID_DXFSOLID )
	    _sdata.push_back( new DXFSolid( is ) );
	else if( fileid == FILEID_STLSOLID )
	    _sdata.push_back( new STLSolid( is ) );
	else
	    throw( Error( ERROR_LOCATION, "unknown solid type" ) );
    }

    for( uint32_t a = 0; a < _n+6; a++ )
	_bound.push_back( Bound( is ) );
    
    _built = read_int8( is );
    _smesh = new uint32_t[_size[0]*_size[1]*_size[2]];
    read_compressed_block( is, sizeof(uint32_t)*_size[0]*_size[1]*_size[2], 
			   (int8_t *)_smesh );

    uint32_t nearsolidsize = read_int32( is );
    _nearsolid.resize( nearsolidsize );
    read_compressed_block( is, sizeof(uint8_t)*nearsolidsize, 
			   (int8_t *)&_nearsolid[0] );

    pthread_mutex_init( &_mutex, NULL );
    pthread_cond_init( &_cond, NULL );

    ibsimu.dec_indent();
}


Geometry::~Geometry()
{
    delete [] _smesh;
    pthread_mutex_destroy( &_mutex );
    pthread_cond_destroy( &_cond );
}


void Geometry::check_definition()
{
    // Geometry specific checks
    if( _geom_mode == MODE_3D ) {
	if( _size[0] < 3 || _size[1] < 3 || _size[2] < 3 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    } else if( _geom_mode == MODE_2D || _geom_mode == MODE_CYL ) {
	if( _geom_mode == MODE_CYL ) {
	    if( _origo(1) != 0.0 )
		throw( Error( ERROR_LOCATION, "j=0 node not on axis" ) );
	}
	if( _size[0] < 3 || _size[1] < 3 || _size[2] != 1 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    } else {
	if( _size[0] < 3 || _size[1] != 1 || _size[2] != 1 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    }
}


void Geometry::set_solid( uint32_t n, const Solid *s )
{
    if( n <= 6 || n > _n+7 )
	throw( Error( ERROR_LOCATION, "illegal solid number " + to_string(n) ) );

    if( n <= _n+6 ) {
	delete _sdata[n-7];
    } else {
	_sdata.push_back( 0 );
	_bound.push_back( Bound(BOUND_DIRICHLET,0.0) );
	_n++;
    }

    _sdata[n-7] = s;
}


uint32_t Geometry::number_of_solids() const
{
    return( _n );
}


uint32_t Geometry::number_of_boundaries() const
{
    return( _n+6 );
}


const Solid *Geometry::get_solid( uint32_t n ) const
{
    if( n <= 6 || n > _n+6 )
	throw( Error( ERROR_LOCATION, "illegal solid number " + to_string(n) ) );
    
    return( _sdata[n-7] );
}


void Geometry::set_boundary( uint32_t n, const Bound &b )
{
    if( n <= 0 || n > _n+6 )
	throw( Error( ERROR_LOCATION, "illegal solid number " + to_string(n) ) );

    if( n >= 7 && b.type() != BOUND_DIRICHLET )
	throw( Error( ERROR_LOCATION, "trying to set solid " + to_string(n) + " as Neumann boundary" ) );

    _bound[n-1] = b;
}


Bound Geometry::get_boundary( uint32_t n ) const
{
    if( n <= 0 || n > _n+6 )
	throw( Error( ERROR_LOCATION, "illegal solid number " + to_string(n) ) );

    return( _bound[n-1] );
}


bool Geometry::have_solid_data( void ) const
{
    for( ssize_t a = _sdata.size()-1; a >= 0 ; a-- ) {
	if( !_sdata[a] )
	    return( false );
    }

    return( true );
}


uint32_t Geometry::inside( const Vec3D &x ) const
{
    for( ssize_t a = _sdata.size()-1; a >= 0 ; a-- ) {
	if( !_sdata[a] )
	    throw( Error( ERROR_LOCATION, "solid " + to_string(a) + " not defined" ) );
	else if( _sdata[a]->inside( x ) )
	    return( a+7 );
    }
    double eps = 1.0e-6*h();
    if( x[2] > _max[2]+eps )
	return( 6 );
    else if( x[2] < _origo[2]-eps )
	return( 5 );
    else if( x[1] > _max[1]+eps )
	return( 4 );
    else if( x[1] < _origo[1]-eps )
	return( 3 );
    else if( x[0] > _max[0]+eps )
	return( 2 );
    else if( x[0] < _origo[0]-eps )
	return( 1 );

    return( 0 );
}


bool Geometry::inside( uint32_t n, const Vec3D &x ) const
{
    if( n == 0 || n > _n+6 )
	 throw( Error( ERROR_LOCATION, "illegal solid number n=" + to_string(n) ) );
    else if( n <= 6 ) {
	double eps = 1.0e-6*h();
	if( n == 6 && x[2] > _max[2]+eps )
	    return( true );
	else if( n == 5 && x[2] < _origo[2]-eps )
	    return( true );
	else if( n == 4 && x[1] > _max[1]+eps )
	    return( true );
	else if( n == 3 && x[1] < _origo[1]-eps )
	    return( true );
	else if( n == 2 && x[0] > _max[0]+eps )
	    return( true );
	else if( n == 1 && x[0] < _origo[0]-eps )
	    return( true );
	else
	    return( false );
    }

    return( _sdata[n-7]->inside( x ) );
}


uint32_t Geometry::mesh_check( int32_t i ) const
{
    if( i < 0 )
	return( SMESH_NODE_ID_DIRICHLET | 1 );
    else if( i >= _size[0] )
	return( SMESH_NODE_ID_DIRICHLET | 2 );

    return( _smesh[i] );
}


uint32_t Geometry::mesh_check( int32_t i, int32_t j ) const
{
    if( i < 0 )
	return( SMESH_NODE_ID_DIRICHLET | 1 );
    else if( i >= _size[0] )
	return( SMESH_NODE_ID_DIRICHLET | 2 );
    if( j < 0 )
	return( SMESH_NODE_ID_DIRICHLET | 3 );
    else if( j >= _size[1] )
	return( SMESH_NODE_ID_DIRICHLET | 4 );

    return( _smesh[i + j*_size[0]] );
}


uint32_t Geometry::mesh_check( int32_t i, int32_t j, int32_t k ) const
{
    if( i < 0 )
	return( SMESH_NODE_ID_DIRICHLET | 1 );
    else if( i >= _size[0] )
	return( SMESH_NODE_ID_DIRICHLET | 2 );
    if( j < 0 )
	return( SMESH_NODE_ID_DIRICHLET | 3 );
    else if( j >= _size[1] )
	return( SMESH_NODE_ID_DIRICHLET | 4 );
    if( k < 0 )
	return( SMESH_NODE_ID_DIRICHLET | 5 );
    else if( k >= _size[2] )
	return( SMESH_NODE_ID_DIRICHLET | 6 );

    return( _smesh[i + _size[0]*(j + k*_size[1])] );
}


double Geometry::bracket_surface( uint32_t n, const Vec3D &xin, const Vec3D &xout, Vec3D &xsurf ) const
{
    Vec3D xl = xin;
    Vec3D xh = xout;

#ifdef DEBUG_BRACKET_SURFACE
    std::cout << "    Bracketing surface between xl and xh\n";
#endif

    // Do iteration
    for( uint32_t a = 0; a < 8; a++ ) {

	xsurf = 0.5*(xl+xh);

#ifdef DEBUG_BRACKET_SURFACE
	std::cout << "      xl            = " << xl << "\n";
	std::cout << "      xh            = " << xh << "\n";
	std::cout << "      Testing xsurf = " << xsurf << "\n";
#endif

	if( inside( n, xsurf ) ) {
#ifdef DEBUG_BRACKET_SURFACE
	    std::cout << "        inside\n";
#endif
	    xl = xsurf;
	} else {
#ifdef DEBUG_BRACKET_SURFACE
	    std::cout << "        outside\n";
#endif
	    xh = xsurf;
	}
    }

    // Calculate best guess, the midpoint
    xsurf = 0.5*(xl+xh);

#ifdef DEBUG_BRACKET_SURFACE
	std::cout << "      Best guess    = " << xsurf << "\n";
#endif

    // Return parametric distance calculated using axis where
    // coordinate difference is largest
    int a;
    Vec3D dif = xout - xin;
    dif.abs();
    if( dif[0] > dif[1] ) {
 	if( dif[0] > dif[2] )
	    a = 0;
	else
	    a = 2;
    } else {
	if( dif[1] > dif[2] )
	    a = 1;
	else
	    a = 2;
    }

    if( xin[a] == xout[a]  )
	throw( Error( ERROR_LOCATION, "xin and xout are the same point" ) );

    return( (xsurf[a] - xin[a]) / (xout[a] - xin[a]) );
}


uint8_t Geometry::bracket_ndist( int32_t i, int32_t j, int32_t k, int32_t solid, int sign, int coord ) const
{
    Vec3D vout( _origo[0]+i*_h, _origo[1]+j*_h, _origo[2]+k*_h );
    uint32_t bp = 0x80;
    uint32_t surf = 0x80;

    // Do iteration
    double step = sign*_h;
    for( uint32_t a = 0; a < 8; a++ ) {
	Vec3D vtest( vout );
	vtest[coord] += step*surf/255.0;
	if( inside( solid, vtest ) )
	    surf -= bp;
	bp = bp >> 1;
	surf += bp;
    }

    // Don't allow zero distance to be returned
    if( surf == 0 )
	return( 1 );
    return( surf );
}


Vec3D Geometry::surface_normal( const Vec3D &x ) const
{
    switch( geom_mode() ) {
    case MODE_2D:
    case MODE_CYL:
	return( surface_normal_2d( x ) );
    case MODE_3D:
	return( surface_normal_3d( x ) );
    default:
	throw( ErrorUnimplemented( ERROR_LOCATION ) );
    }
}


Vec3D Geometry::surface_normal_3d( const Vec3D &x ) const
{
    double dx = 0.1*_h;

    Vec3D px[8];
    uint32_t p[8];
    int surc;
    bool sur[12];

    do {
	// Build corner coordinates
	px[0] = x+Vec3D( -dx, -dx, -dx );
	px[1] = x+Vec3D( -dx, +dx, -dx );
	px[2] = x+Vec3D( +dx, +dx, -dx );
	px[3] = x+Vec3D( +dx, -dx, -dx );
	px[4] = x+Vec3D( -dx, -dx, +dx );
	px[5] = x+Vec3D( -dx, +dx, +dx );
	px[6] = x+Vec3D( +dx, +dx, +dx );
	px[7] = x+Vec3D( +dx, -dx, +dx );

	// Check inside() in corners
	p[0] = inside( px[0] );
	p[1] = inside( px[1] );
	p[2] = inside( px[2] );
	p[3] = inside( px[3] );
	p[4] = inside( px[4] );
	p[5] = inside( px[5] );
	p[6] = inside( px[6] );
	p[7] = inside( px[7] );

	// Check if at least three intersections
	surc = 0;
	if( p[0] != p[1] ) {
	    surc++;
	    sur[0] = true;
	} else {
	    sur[0] = false;
	}
	if( p[1] != p[2] ) {
	    surc++;
	    sur[1] = true;
	} else {
	    sur[1] = false;
	}
	if( p[2] != p[3] ) {
	    surc++;
	    sur[2] = true;
	} else {
	    sur[2] = false;
	}
	if( p[3] != p[0] ) {
	    surc++;
	    sur[3] = true;
	} else {
	    sur[3] = false;
	}
	//
	if( p[4] != p[5] ) {
	    surc++;
	    sur[4] = true;
	} else {
	    sur[4] = false;
	}
	if( p[5] != p[6] ) {
	    surc++;
	    sur[5] = true;
	} else {
	    sur[5] = false;
	}
	if( p[6] != p[7] ) {
	    surc++;
	    sur[6] = true;
	} else {
	    sur[6] = false;
	}
	if( p[7] != p[4] ) {
	    surc++;
	    sur[7] = true;
	} else {
	    sur[7] = false;
	}
	//
	if( p[0] != p[4] ) {
	    surc++;
	    sur[8] = true;
	} else {
	    sur[8] = false;
	}
	if( p[1] != p[5] ) {
	    surc++;
	    sur[9] = true;
	} else {
	    sur[9] = false;
	}
	if( p[2] != p[6] ) {
	    surc++;
	    sur[10] = true;
	} else {
	    sur[10] = false;
	}
	if( p[3] != p[7] ) {
	    surc++;
	    sur[11] = true;
	} else {
	    sur[11] = false;
	}

	dx *= 2.0;
	if( dx >= _h )
	    return( Vec3D(0,0,0) );

    } while( surc < 3 );

    // Bracket surface locations
    int b = 0;
    Vec3D xs[12];
    if( sur[0] ) {
	if( p[0] )
	    bracket_surface( p[0], px[0], px[1], xs[b++] );
	else
	    bracket_surface( p[1], px[1], px[0], xs[b++] );
    }
    if( sur[1] ) {
	if( p[1] )
	    bracket_surface( p[1], px[1], px[2], xs[b++] );
	else
	    bracket_surface( p[2], px[2], px[1], xs[b++] );
    }
    if( sur[2] ) {
	if( p[2] )
	    bracket_surface( p[2], px[2], px[3], xs[b++] );
	else
	    bracket_surface( p[3], px[3], px[2], xs[b++] );
    }
    if( sur[3] ) {
	if( p[3] )
	    bracket_surface( p[3], px[3], px[0], xs[b++] );
	else
	    bracket_surface( p[0], px[0], px[3], xs[b++] );
    }
    //
    if( sur[4] ) {
	if( p[4] )
	    bracket_surface( p[4], px[4], px[5], xs[b++] );
	else
	    bracket_surface( p[5], px[5], px[4], xs[b++] );
    }
    if( sur[5] ) {
	if( p[5] )
	    bracket_surface( p[5], px[5], px[6], xs[b++] );
	else
	    bracket_surface( p[6], px[6], px[5], xs[b++] );
    }
    if( sur[6] ) {
	if( p[6] )
	    bracket_surface( p[6], px[6], px[7], xs[b++] );
	else
	    bracket_surface( p[7], px[7], px[6], xs[b++] );
    }
    if( sur[7] ) {
	if( p[7] )
	    bracket_surface( p[7], px[7], px[4], xs[b++] );
	else
	    bracket_surface( p[4], px[4], px[7], xs[b++] );
    }
    //
    if( sur[8] ) {
	if( p[0] )
	    bracket_surface( p[0], px[0], px[4], xs[b++] );
	else
	    bracket_surface( p[4], px[4], px[0], xs[b++] );
    }
    if( sur[9] ) {
	if( p[1] )
	    bracket_surface( p[1], px[1], px[5], xs[b++] );
	else
	    bracket_surface( p[5], px[5], px[1], xs[b++] );
    }
    if( sur[10] ) {
	if( p[2] )
	    bracket_surface( p[2], px[2], px[6], xs[b++] );
	else
	    bracket_surface( p[6], px[6], px[2], xs[b++] );
    }
    if( sur[11] ) {
	if( p[3] )
	    bracket_surface( p[3], px[3], px[7], xs[b++] );
	else
	    bracket_surface( p[7], px[7], px[3], xs[b++] );
    }

    // Choose for three first intersections for normal
    // Three points furthest from each other would be better, but more hassle to search
    Vec3D v1 = xs[1]-xs[0];
    Vec3D v2 = xs[2]-xs[0];
    Vec3D n = cross( v1, v2 );
    n.normalize();

    // Check polarity
    if( inside( x+0.1*h()*n ) )
	n *= -1.0;

    return( n );
}


Vec3D Geometry::surface_normal_2d( const Vec3D &x ) const
{
    double dx = 0.1*_h;

    Vec3D px[4];
    uint32_t p[4];
    int surc;
    bool sur[4];

    do {
	// Build corner coordinates
	px[0] = x+Vec3D(-dx,-dx,0);
	px[1] = x+Vec3D(-dx,+dx,0);
	px[2] = x+Vec3D(+dx,+dx,0);
	px[3] = x+Vec3D(+dx,-dx,0);

	// Check inside() in corners
	p[0] = inside( px[0] );
	p[1] = inside( px[1] );
	p[2] = inside( px[2] );
	p[3] = inside( px[3] );

	// Check if two intersections
	surc = 0;
	if( p[0] != p[1] ) {
	    surc++;
	    sur[0] = true;
	} else {
	    sur[0] = false;
	}
	if( p[1] != p[2] ) {
	    surc++;
	    sur[1] = true;
	} else {
	    sur[1] = false;
	}
	if( p[2] != p[3] ) {
	    surc++;
	    sur[2] = true;
	} else {
	    sur[2] = false;
	}
	if( p[3] != p[0] ) {
	    surc++;
	    sur[3] = true;
	} else {
	    sur[3] = false;
	}

	dx *= 2.0;
	if( dx >= _h )
	    return( Vec3D(0,0,0) );

    } while( surc != 2 );

    // Bracket surface locations
    int pol = 0;
    int b = 0;
    Vec3D xs[2];
    for( int a = 0; a < 4; a++ ) {
	if( sur[a] ) {
	    pol = p[a];
	    if( p[a] )
		bracket_surface( p[a], px[a], px[(a+1)%4], xs[b++] );
	    else
		bracket_surface( p[(a+1)%4], px[(a+1)%4], px[a], xs[b++] );
	}
    }

    Vec3D n( xs[1][1]-xs[0][1], xs[0][0]-xs[1][0], 0.0 );
    if( !pol )
	n *= -1.0;

    n.normalize();

    //std::cout << "normal = " << n << "\n";

    return( n );
}


uint32_t Geometry::is_solid( int32_t i, int32_t j, int32_t k ) const
{
    uint32_t a = mesh_check( i, j, k );
    if( (a & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(a & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	return( a & SMESH_BOUNDARY_NUMBER_MASK );
    return( 0 );
}


bool Geometry::is_near_solid( int32_t i, int32_t j, int32_t k ) const
{
    return( is_solid(i-1,j,  k  ) ||
	    is_solid(i+1,j,  k  ) ||
	    is_solid(i,  j-1,k  ) || 
	    is_solid(i,  j+1,k  ) ||
	    is_solid(i,  j,  k-1) || 
	    is_solid(i,  j,  k+1) );
}


void Geometry::add_near_solid_entry( uint32_t &near_solid_index, int32_t i, int32_t j, int32_t k )
{
    uint32_t solid = 0;         // Solid number of neighbour
    int nsolids = 0;            // Number of near neighbour solids nodes.
    uint8_t neighbours = 0;     // Bit flags for neighbours
    uint8_t ndist[7];           // Neighbour distances

    // X
    if( (solid = is_solid(i-1,j,k)) ) {
	neighbours += 0x20;
	nsolids++;
	ndist[nsolids] = bracket_ndist( i,j,k, solid, -1, 0 );
    }
    neighbours = neighbours >> 1;
    if( ( solid = is_solid(i+1,j,k)) ) {
	neighbours += 0x20;
	nsolids++;
	ndist[nsolids] = bracket_ndist( i,j,k, solid, +1, 0 );
    }
    neighbours = neighbours >> 1;

    // Y
    if( (solid = is_solid(i,j-1,k)) ) {
	neighbours += 0x20;
	nsolids++;
	ndist[nsolids] = bracket_ndist( i,j,k, solid, -1, 1 );
    }
    neighbours = neighbours >> 1;
    if( (solid = is_solid(i,j+1,k)) ) {
	neighbours += 0x20;
	nsolids++;
	ndist[nsolids] = bracket_ndist( i,j,k, solid, +1, 1 );
    }
    neighbours = neighbours >> 1;

    // Z
    if( (solid = is_solid(i,j,k-1)) ) {
	neighbours += 0x20;
	nsolids++;
	ndist[nsolids] = bracket_ndist( i,j,k, solid, -1, 2 );
    }
    neighbours = neighbours >> 1;
    if( (solid = is_solid(i,j,k+1)) ) {
	neighbours += 0x20;
	nsolids++;
	ndist[nsolids] = bracket_ndist( i,j,k, solid, +1, 2 );
    }

    uint32_t ind = near_solid_index;
    near_solid_index += nsolids+1;
    _nearsolid.resize( near_solid_index );
    ndist[0] = neighbours;
    memcpy( (void *)&_nearsolid[ind], (void *)ndist, nsolids+1 );
}


void Geometry::build_mesh_parallel_near_solid( uint32_t ind, int32_t i, int32_t j, int32_t k )
{
    uint32_t solid = 0;         // Solid number of neighbour
    int nsolids = 0;            // Number of near neighbour solids nodes.
    uint8_t neighbours = 0;     // Bit flags for neighbours
    uint8_t ndist[7];           // Neighbour distances
    
    // X
    if( (solid = is_solid(i-1,j,k)) ) {
	neighbours += 0x20;
	nsolids++;
	ndist[nsolids] = bracket_ndist( i,j,k, solid, -1, 0 );
    }
    neighbours = neighbours >> 1;
    if( ( solid = is_solid(i+1,j,k)) ) {
	neighbours += 0x20;
	nsolids++;
	ndist[nsolids] = bracket_ndist( i,j,k, solid, +1, 0 );
    }
    neighbours = neighbours >> 1;

    // Y
    if( (solid = is_solid(i,j-1,k)) ) {
	neighbours += 0x20;
	nsolids++;
	ndist[nsolids] = bracket_ndist( i,j,k, solid, -1, 1 );
    }
    neighbours = neighbours >> 1;
    if( (solid = is_solid(i,j+1,k)) ) {
	neighbours += 0x20;
	nsolids++;
	ndist[nsolids] = bracket_ndist( i,j,k, solid, +1, 1 );
    }
    neighbours = neighbours >> 1;

    // Z
    if( (solid = is_solid(i,j,k-1)) ) {
	neighbours += 0x20;
	nsolids++;
	ndist[nsolids] = bracket_ndist( i,j,k, solid, -1, 2 );
    }
    neighbours = neighbours >> 1;
    if( (solid = is_solid(i,j,k+1)) ) {
	neighbours += 0x20;
	nsolids++;
	ndist[nsolids] = bracket_ndist( i,j,k, solid, +1, 2 );
    }

    ndist[0] = neighbours;
    memcpy( (void *)&_nearsolid[ind], (void *)ndist, nsolids+1 );
}


void Geometry::build_mesh_parallel_prepare_near_solid( uint32_t &near_solid_index, 
						       int32_t i, int32_t j, int32_t k )
{
    int nsolids = 0;

    // Count number of near neighbour solid nodes
    if( is_solid(i-1,j,k) )
	nsolids++;
    if( is_solid(i+1,j,k) )
	nsolids++;
    if( is_solid(i,j-1,k) )
	nsolids++;
    if( is_solid(i,j+1,k) )
	nsolids++;
    if( is_solid(i,j,k-1) )
	nsolids++;
    if( is_solid(i,j,k+1) )
	nsolids++;

    near_solid_index += 1 + sizeof(uint8_t)*nsolids;
}


void Geometry::build_mesh_parallel_prepare_3d( void )
{
    uint32_t near_solid_index = 0;

    // Mark all nodes and count near solid data size
    for( int32_t k = 0; k < _size[2]; k++ ) {
	for( int32_t j = 0; j < _size[1]; j++ ) {
	    for( int32_t i = 0; i < _size[0]; i++ ) {

		if( mesh(i,j,k) != 0 )
		    continue;

		if( is_near_solid(i,j,k) ) {
		    // Near solid
		    mesh(i,j,k) = SMESH_NODE_ID_NEAR_SOLID | near_solid_index;
		    build_mesh_parallel_prepare_near_solid( near_solid_index, i, j, k );

		} else if( i == 0 && get_boundary(1).type() == BOUND_DIRICHLET ) {
		    // Xmin Dirichlet boundary
		    mesh(i,j,k) = SMESH_NODE_ID_DIRICHLET | 1;

		} else if( i == _size[0]-1 && get_boundary(2).type() == BOUND_DIRICHLET ) {
		    // Xmax Dirichlet boundary
		    mesh(i,j,k) = SMESH_NODE_ID_DIRICHLET | 2;

		} else if( j == 0 && get_boundary(3).type() == BOUND_DIRICHLET ) {
		    // Ymin Dirichlet boundary
		    mesh(i,j,k) = SMESH_NODE_ID_DIRICHLET | 3;

		} else if( j == _size[1]-1 && get_boundary(4).type() == BOUND_DIRICHLET ) {
		    // Ymax Dirichlet boundary
		    mesh(i,j,k) = SMESH_NODE_ID_DIRICHLET | 4;

		} else if( k == 0 && get_boundary(5).type() == BOUND_DIRICHLET ) {
		    // Zmin Dirichlet boundary
		    mesh(i,j,k) = SMESH_NODE_ID_DIRICHLET | 5;

		} else if( k == _size[2]-1 && get_boundary(6).type() == BOUND_DIRICHLET ) {
		    // Zmax Dirichlet boundary
		    mesh(i,j,k) = SMESH_NODE_ID_DIRICHLET | 6;

		} else if( i == 0 ) {
		    // Xmin Neumann boundary
		    mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | 1;

		} else if( i == _size[0]-1 ) {
		    // Xmax Neumann boundary
		    mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | 2;
		    
		} else if( j == 0 ) {
		    // Ymin Neumann boundary
		    mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | 3;

		} else if( j == _size[1]-1 ) {
		    // Ymax Neumann boundary
		    mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | 4;

		} else if( k == 0 ) {
		    // Zmin Neumann boundary
		    mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | 5;

		} else if( k == _size[2]-1 ) {
		    // Zmax Neumann boundary
		    mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | 6;

		} else {
		    // Pure vacuum
		    mesh(i,j,k) = SMESH_NODE_ID_PURE_VACUUM;
		}
	    }
	}
    }

    // Reserve space for near solid data
    _nearsolid.resize( near_solid_index );
}


void Geometry::build_mesh_parallel_prepare_2d( void )
{
    uint32_t near_solid_index = 0;

    // Mark all nodes and count near solid data size
    for( int32_t j = 0; j < _size[1]; j++ ) {
	for( int32_t i = 0; i < _size[0]; i++ ) {

	    if( mesh(i,j) != 0 )
		continue;

	    if( is_near_solid(i,j,0) ) {
		// Near solid
		mesh(i,j) = SMESH_NODE_ID_NEAR_SOLID | near_solid_index;
		build_mesh_parallel_prepare_near_solid( near_solid_index, i, j, 0 );

	    } else if( i == 0 && get_boundary(1).type() == BOUND_DIRICHLET ) {
		// Xmin Dirichlet boundary
		mesh(i,j) = SMESH_NODE_ID_DIRICHLET | 1;
		
	    } else if( i == _size[0]-1 && get_boundary(2).type() == BOUND_DIRICHLET ) {
		// Xmax Dirichlet boundary
		mesh(i,j) = SMESH_NODE_ID_DIRICHLET | 2;
		
	    } else if( j == 0 && get_boundary(3).type() == BOUND_DIRICHLET ) {
		// Ymin Dirichlet boundary
		mesh(i,j) = SMESH_NODE_ID_DIRICHLET | 3;
		
	    } else if( j == _size[1]-1 && get_boundary(4).type() == BOUND_DIRICHLET ) {
		// Ymax Dirichlet boundary
		mesh(i,j) = SMESH_NODE_ID_DIRICHLET | 4;
		
	    } else if( i == 0 ) {
		// Xmin Neumann boundary
		mesh(i,j) = SMESH_NODE_ID_NEUMANN | 1;
		
	    } else if( i == _size[0]-1 ) {
		// Xmax Neumann boundary
		mesh(i,j) = SMESH_NODE_ID_NEUMANN | 2;
		
	    } else if( j == 0 ) {
		// Ymin Neumann boundary
		mesh(i,j) = SMESH_NODE_ID_NEUMANN | 3;
		
	    } else if( j == _size[1]-1 ) {
		// Ymax Neumann boundary
		mesh(i,j) = SMESH_NODE_ID_NEUMANN | 4;
		
	    } else {
		// Pure vacuum
		mesh(i,j) = SMESH_NODE_ID_PURE_VACUUM;
	    }
	}
    }

    // Reserve space for near solid data
    _nearsolid.resize( near_solid_index );
}


void Geometry::build_mesh_parallel_prepare_1d( void )
{
    uint32_t near_solid_index = 0;

    // Mark all nodes and count near solid data size
    for( int32_t i = 0; i < _size[0]; i++ ) {

	if( mesh(i) != 0 )
	    continue;
	
	if( is_near_solid(i,0,0) ) {
	    // Near solid
	    mesh(i) = SMESH_NODE_ID_NEAR_SOLID | near_solid_index;
	    build_mesh_parallel_prepare_near_solid( near_solid_index, i, 0, 0 );
	    
	} else if( i == 0 && get_boundary(1).type() == BOUND_DIRICHLET ) {
	    // Xmin Dirichlet boundary
	    mesh(i) = SMESH_NODE_ID_DIRICHLET | 1;
	    
	} else if( i == _size[0]-1 && get_boundary(2).type() == BOUND_DIRICHLET ) {
	    // Xmax Dirichlet boundary
	    mesh(i) = SMESH_NODE_ID_DIRICHLET | 2;
	    
	} else if( i == 0 ) {
	    // Xmin Neumann boundary
	    mesh(i) = SMESH_NODE_ID_NEUMANN | 1;
	    
	} else if( i == _size[0]-1 ) {
	    // Xmax Neumann boundary
	    mesh(i) = SMESH_NODE_ID_NEUMANN | 2;
	    
	} else {
	    // Pure vacuum
	    mesh(i) = SMESH_NODE_ID_PURE_VACUUM;
	}
    }

    // Reserve space for near solid data
    _nearsolid.resize( near_solid_index );
}


void Geometry::build_mesh_parallel_thread_3d( BuildMeshData *bmd )
{
    // Parallel: Mark solid (Dirichlet) nodes. Others left to zero.
    for( int32_t k = bmd->index; k < _size[2]; k += ibsimu.get_thread_count() ) {
	double z = k*_h+_origo[2];
	for( int32_t j = 0; j < _size[1]; j++ ) {
	    double y = j*_h+_origo[1];
	    for( int32_t i = 0; i < _size[0]; i++ ) {
		double x = i*_h+_origo[0];
		uint32_t nid = inside( Vec3D(x,y,z) );
		if( nid )
		    mesh(i,j,k) = SMESH_NODE_ID_DIRICHLET | nid;
		else
		    mesh(i,j,k) = 0;
	    }
	}
    }

    // Serial part: edges and preparation for near solid nodes
    pthread_mutex_lock( &_mutex );
    _done++;
    if( _done == ibsimu.get_thread_count() ) {
	// Mark rest of mesh nodes and prepare for building near solid data
	build_mesh_parallel_prepare_3d();
	pthread_cond_broadcast( &_cond );
    } else {
	do {
	    pthread_cond_wait( &_cond, &_mutex );
	} while( _done != ibsimu.get_thread_count() );
    }
    pthread_mutex_unlock( &_mutex );

    // Parallel: Build near solid data
    for( int32_t k = bmd->index; k < _size[2]; k += ibsimu.get_thread_count() ) {
	for( int32_t j = 0; j < _size[1]; j++ ) {
	    for( int32_t i = 0; i < _size[0]; i++ ) {

		uint32_t node = mesh(i,j,k);
		uint32_t node_id = node & SMESH_NODE_ID_MASK;
		if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		    uint32_t index = node & SMESH_NEAR_SOLID_INDEX_MASK;
		    build_mesh_parallel_near_solid( index, i, j, k );
		}
	    }
	}
    }
}


void Geometry::build_mesh_parallel_thread_2d( BuildMeshData *bmd )
{
    // Parallel: Mark solid (Dirichlet) nodes. Others left to zero.
    for( int32_t j = bmd->index; j < _size[1]; j += ibsimu.get_thread_count() ) {
	double y = j*_h+_origo[1];
	for( int32_t i = 0; i < _size[0]; i++ ) {
	    double x = i*_h+_origo[0];
	    uint32_t nid = inside( Vec3D(x,y) );
	    if( nid )
		mesh(i,j) = SMESH_NODE_ID_DIRICHLET | nid;
	    else
		mesh(i,j) = 0;
	}
    }

    // Serial part: edges and preparation for near solid nodes
    pthread_mutex_lock( &_mutex );
    _done++;
    if( _done == ibsimu.get_thread_count() ) {
	// Mark rest of mesh nodes and prepare for building near solid data
	build_mesh_parallel_prepare_2d();
	pthread_cond_broadcast( &_cond );
    } else {
	do {
	    pthread_cond_wait( &_cond, &_mutex );
	} while( _done != ibsimu.get_thread_count() );
    }
    pthread_mutex_unlock( &_mutex );

    // Parallel: Build near solid data
    for( int32_t j = bmd->index; j < _size[1]; j += ibsimu.get_thread_count() ) {
	for( int32_t i = 0; i < _size[0]; i++ ) {
	    
	    uint32_t node = mesh(i,j);
	    uint32_t node_id = node & SMESH_NODE_ID_MASK;
	    if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
		uint32_t index = node & SMESH_NEAR_SOLID_INDEX_MASK;
		build_mesh_parallel_near_solid( index, i, j, 0 );
	    }
	}
    }
}


void Geometry::build_mesh_parallel_thread_1d( void )
{
    // Mark solid (Dirichlet) nodes. Others left to zero.
    for( int32_t i = 0; i < _size[0]; i++ ) {
	double x = i*_h+_origo[0];
	uint32_t nid = inside( Vec3D(x) );
	if( nid )
	    mesh(i) = SMESH_NODE_ID_DIRICHLET | nid;
	else
	    mesh(i) = 0;
    }

    // Mark rest of mesh nodes and prepare for building near solid data
    build_mesh_parallel_prepare_1d();

    // Build near solid data
    for( int32_t i = 0; i < _size[0]; i++ ) {	    
	uint32_t node = mesh(i);
	uint32_t node_id = node & SMESH_NODE_ID_MASK;
	if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
	    uint32_t index = node & SMESH_NEAR_SOLID_INDEX_MASK;
	    build_mesh_parallel_near_solid( index, i, 0, 0 );
	}
    }
}


void Geometry::build_mesh_parallel_thread( BuildMeshData *bmd )
{
    if( _geom_mode == MODE_3D )
	build_mesh_parallel_thread_3d( bmd );
    else if( _geom_mode == MODE_2D || _geom_mode == MODE_CYL )
	build_mesh_parallel_thread_2d( bmd );
    else
	throw( ErrorAssert( ERROR_LOCATION ) );
}


void *Geometry::build_mesh_parallel_entry( void *data )
{
    BuildMeshData *bmd = (BuildMeshData *)data;
    bmd->geom->build_mesh_parallel_thread( bmd );
    return( NULL );
}


void Geometry::build_mesh_parallel( void )
{
    // Mesh building in 1d is not parallelized
    if( _geom_mode == MODE_1D ) {
	build_mesh_parallel_thread_1d();
	return;
    }

    BuildMeshData bmd[ibsimu.get_thread_count()];

    // Prepare common data
    _done = 0;

    // Start threads
    for( uint32_t a = 0; a < ibsimu.get_thread_count(); a++ ) {
	bmd[a].index = a;
	bmd[a].geom = this;
	pthread_create( &bmd[a].thread, NULL, build_mesh_parallel_entry, (void *)&bmd[a] );
    }

    // Join threads
    for( uint32_t a = 0; a < ibsimu.get_thread_count(); a++ )
	pthread_join( bmd[a].thread, NULL );    
}


void Geometry::build_mesh( void )
{
    Timer t;
    ibsimu.message( 1 ) << "Building mesh\n";
    ibsimu.inc_indent();

    // Check the solid definitions
    for( ssize_t a = _sdata.size()-1; a >= 0 ; a-- ) {
	if( !_sdata[a] )
	    throw( Error( ERROR_LOCATION, "solid " + to_string(a) + " not defined" ) );
    }

    _built = true;
    build_mesh_parallel();

    // Report node counts
    int b;
    uint32_t ncount = nodecount();
    uint32_t nvacuum = 0;
    uint32_t nnearsolid = 0;
    uint32_t nneumann = 0;
    uint32_t ndirichlet = 0;
    int nsolid[_n];
    for( uint32_t a = 0; a < _n; a++ )
	nsolid[a] = 0;

    for( uint32_t a = 0; a < ncount; a++ ) {
	
	switch( mesh(a) & SMESH_NODE_ID_MASK ) {
	case SMESH_NODE_ID_NEAR_SOLID:
	    nnearsolid++;
	    break;
	case SMESH_NODE_ID_PURE_VACUUM:
	    nvacuum++;
	    break;
	case SMESH_NODE_ID_NEUMANN:
	    nneumann++;
	    if( (b = (mesh(a) & SMESH_BOUNDARY_NUMBER_MASK)) >= 7 )
		nsolid[b-7]++;
	    break;
	case SMESH_NODE_ID_DIRICHLET:
	    ndirichlet++;
	    if( (b = (mesh(a) & SMESH_BOUNDARY_NUMBER_MASK)) >= 7 )
		nsolid[b-7]++;
	    break;
	}
    }
    
    ibsimu.message( 1 ) << "Done. Built mesh with:\n";
    ibsimu.message( 1 ) << nvacuum << " pure vacuum nodes\n";
    ibsimu.message( 1 ) << nnearsolid << " near solid nodes\n";
    ibsimu.message( 1 ) << nneumann << " neumann nodes\n";
    ibsimu.message( 1 ) << ndirichlet << " dirichlet nodes\n";
    for( uint32_t a = 0; a < _n; a++ )
	ibsimu.message( 1 ) << nsolid[a] << " solid " << a+7 << " nodes\n";

    // End timer
    t.stop();

    ibsimu.message( 1 ) << "time used = " << t << "\n";
    ibsimu.dec_indent();
}


/* Faces array
 *
 * TrianglesA.lut.txt
 */
const int Geometry::mc_faces[15*256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    1, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 8, 3, 1, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    9, 2, 11, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 8, 3, 2, 11, 8, 11, 9, 8, -1, -1, -1, -1, -1, -1,
    3, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 10, 2, 8, 10, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    1, 9, 0, 2, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    1, 10, 2, 1, 9, 10, 9, 8, 10, -1, -1, -1, -1, -1, -1,
    3, 11, 1, 10, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 11, 1, 0, 8, 11, 8, 10, 11, -1, -1, -1, -1, -1, -1,
    3, 9, 0, 3, 10, 9, 10, 11, 9, -1, -1, -1, -1, -1, -1,
    9, 8, 11, 11, 8, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1,
    1, 2, 11, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    3, 4, 7, 3, 0, 4, 1, 2, 11, -1, -1, -1, -1, -1, -1,
    9, 2, 11, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1,
    2, 11, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1,
    8, 4, 7, 3, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    10, 4, 7, 10, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1,
    9, 0, 1, 8, 4, 7, 2, 3, 10, -1, -1, -1, -1, -1, -1,
    4, 7, 10, 9, 4, 10, 9, 10, 2, 9, 2, 1, -1, -1, -1,
    3, 11, 1, 3, 10, 11, 7, 8, 4, -1, -1, -1, -1, -1, -1,
    1, 10, 11, 1, 4, 10, 1, 0, 4, 7, 10, 4, -1, -1, -1,
    4, 7, 8, 9, 0, 10, 9, 10, 11, 10, 0, 3, -1, -1, -1,
    4, 7, 10, 4, 10, 9, 9, 10, 11, -1, -1, -1, -1, -1, -1,
    9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1,
    1, 2, 11, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    3, 0, 8, 1, 2, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1,
    5, 2, 11, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1,
    2, 11, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1,
    9, 5, 4, 2, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 10, 2, 0, 8, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1,
    0, 5, 4, 0, 1, 5, 2, 3, 10, -1, -1, -1, -1, -1, -1,
    2, 1, 5, 2, 5, 8, 2, 8, 10, 4, 8, 5, -1, -1, -1,
    11, 3, 10, 11, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1,
    4, 9, 5, 0, 8, 1, 8, 11, 1, 8, 10, 11, -1, -1, -1,
    5, 4, 0, 5, 0, 10, 5, 10, 11, 10, 0, 3, -1, -1, -1,
    5, 4, 8, 5, 8, 11, 11, 8, 10, -1, -1, -1, -1, -1, -1,
    9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1,
    0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1,
    1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    9, 7, 8, 9, 5, 7, 11, 1, 2, -1, -1, -1, -1, -1, -1,
    11, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1,
    8, 0, 2, 8, 2, 5, 8, 5, 7, 11, 5, 2, -1, -1, -1,
    2, 11, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1,
    7, 9, 5, 7, 8, 9, 3, 10, 2, -1, -1, -1, -1, -1, -1,
    9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 10, -1, -1, -1,
    2, 3, 10, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1,
    10, 2, 1, 10, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1,
    9, 5, 8, 8, 5, 7, 11, 1, 3, 11, 3, 10, -1, -1, -1,
    5, 7, 0, 5, 0, 9, 7, 10, 0, 1, 0, 11, 10, 11, 0,
    10, 11, 0, 10, 0, 3, 11, 5, 0, 8, 0, 7, 5, 7, 0,
    10, 11, 5, 7, 10, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    11, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 8, 3, 5, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    9, 0, 1, 5, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    1, 8, 3, 1, 9, 8, 5, 11, 6, -1, -1, -1, -1, -1, -1,
    1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1,
    9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1,
    5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1,
    2, 3, 10, 11, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    10, 0, 8, 10, 2, 0, 11, 6, 5, -1, -1, -1, -1, -1, -1,
    0, 1, 9, 2, 3, 10, 5, 11, 6, -1, -1, -1, -1, -1, -1,
    5, 11, 6, 1, 9, 2, 9, 10, 2, 9, 8, 10, -1, -1, -1,
    6, 3, 10, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1,
    0, 8, 10, 0, 10, 5, 0, 5, 1, 5, 10, 6, -1, -1, -1,
    3, 10, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1,
    6, 5, 9, 6, 9, 10, 10, 9, 8, -1, -1, -1, -1, -1, -1,
    5, 11, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 3, 0, 4, 7, 3, 6, 5, 11, -1, -1, -1, -1, -1, -1,
    1, 9, 0, 5, 11, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1,
    11, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1,
    6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1,
    1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1,
    8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1,
    7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9,
    3, 10, 2, 7, 8, 4, 11, 6, 5, -1, -1, -1, -1, -1, -1,
    5, 11, 6, 4, 7, 2, 4, 2, 0, 2, 7, 10, -1, -1, -1,
    0, 1, 9, 4, 7, 8, 2, 3, 10, 5, 11, 6, -1, -1, -1,
    9, 2, 1, 9, 10, 2, 9, 4, 10, 7, 10, 4, 5, 11, 6,
    8, 4, 7, 3, 10, 5, 3, 5, 1, 5, 10, 6, -1, -1, -1,
    5, 1, 10, 5, 10, 6, 1, 0, 10, 7, 10, 4, 0, 4, 10,
    0, 5, 9, 0, 6, 5, 0, 3, 6, 10, 6, 3, 8, 4, 7,
    6, 5, 9, 6, 9, 10, 4, 7, 9, 7, 10, 9, -1, -1, -1,
    11, 4, 9, 6, 4, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 11, 6, 4, 9, 11, 0, 8, 3, -1, -1, -1, -1, -1, -1,
    11, 0, 1, 11, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1,
    8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 11, -1, -1, -1,
    1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1,
    3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1,
    0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1,
    11, 4, 9, 11, 6, 4, 10, 2, 3, -1, -1, -1, -1, -1, -1,
    0, 8, 2, 2, 8, 10, 4, 9, 11, 4, 11, 6, -1, -1, -1,
    3, 10, 2, 0, 1, 6, 0, 6, 4, 6, 1, 11, -1, -1, -1,
    6, 4, 1, 6, 1, 11, 4, 8, 1, 2, 1, 10, 8, 10, 1,
    9, 6, 4, 9, 3, 6, 9, 1, 3, 10, 6, 3, -1, -1, -1,
    8, 10, 1, 8, 1, 0, 10, 6, 1, 9, 1, 4, 6, 4, 1,
    3, 10, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1,
    6, 4, 8, 10, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    7, 11, 6, 7, 8, 11, 8, 9, 11, -1, -1, -1, -1, -1, -1,
    0, 7, 3, 0, 11, 7, 0, 9, 11, 6, 7, 11, -1, -1, -1,
    11, 6, 7, 1, 11, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1,
    11, 6, 7, 11, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1,
    1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1,
    2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9,
    7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1,
    7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 10, 11, 6, 8, 11, 8, 9, 8, 6, 7, -1, -1, -1,
    2, 0, 7, 2, 7, 10, 0, 9, 7, 6, 7, 11, 9, 11, 7,
    1, 8, 0, 1, 7, 8, 1, 11, 7, 6, 7, 11, 2, 3, 10,
    10, 2, 1, 10, 1, 7, 11, 6, 1, 6, 7, 1, -1, -1, -1,
    8, 9, 6, 8, 6, 7, 9, 1, 6, 10, 6, 3, 1, 3, 6,
    0, 9, 1, 10, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    7, 8, 0, 7, 0, 6, 3, 10, 0, 10, 6, 0, -1, -1, -1,
    7, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    7, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    3, 0, 8, 10, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 9, 10, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    8, 1, 9, 8, 3, 1, 10, 7, 6, -1, -1, -1, -1, -1, -1,
    11, 1, 2, 6, 10, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    1, 2, 11, 3, 0, 8, 6, 10, 7, -1, -1, -1, -1, -1, -1,
    2, 9, 0, 2, 11, 9, 6, 10, 7, -1, -1, -1, -1, -1, -1,
    6, 10, 7, 2, 11, 3, 11, 8, 3, 11, 9, 8, -1, -1, -1,
    7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1,
    2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1,
    1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1,
    11, 7, 6, 11, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1,
    11, 7, 6, 1, 7, 11, 1, 8, 7, 1, 0, 8, -1, -1, -1,
    0, 3, 7, 0, 7, 11, 0, 11, 9, 6, 11, 7, -1, -1, -1,
    7, 6, 11, 7, 11, 8, 8, 11, 9, -1, -1, -1, -1, -1, -1,
    6, 8, 4, 10, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    3, 6, 10, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1,
    8, 6, 10, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1,
    9, 4, 6, 9, 6, 3, 9, 3, 1, 10, 3, 6, -1, -1, -1,
    6, 8, 4, 6, 10, 8, 2, 11, 1, -1, -1, -1, -1, -1, -1,
    1, 2, 11, 3, 0, 10, 0, 6, 10, 0, 4, 6, -1, -1, -1,
    4, 10, 8, 4, 6, 10, 0, 2, 9, 2, 11, 9, -1, -1, -1,
    11, 9, 3, 11, 3, 2, 9, 4, 3, 10, 3, 6, 4, 6, 3,
    8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1,
    0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1,
    1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1,
    8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 11, 1, -1, -1, -1,
    11, 1, 0, 11, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1,
    4, 6, 3, 4, 3, 8, 6, 11, 3, 0, 3, 9, 11, 9, 3,
    11, 9, 4, 6, 11, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 9, 5, 7, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 8, 3, 4, 9, 5, 10, 7, 6, -1, -1, -1, -1, -1, -1,
    5, 0, 1, 5, 4, 0, 7, 6, 10, -1, -1, -1, -1, -1, -1,
    10, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1,
    9, 5, 4, 11, 1, 2, 7, 6, 10, -1, -1, -1, -1, -1, -1,
    6, 10, 7, 1, 2, 11, 0, 8, 3, 4, 9, 5, -1, -1, -1,
    7, 6, 10, 5, 4, 11, 4, 2, 11, 4, 0, 2, -1, -1, -1,
    3, 4, 8, 3, 5, 4, 3, 2, 5, 11, 5, 2, 10, 7, 6,
    7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1,
    9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1,
    3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1,
    6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8,
    9, 5, 4, 11, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1,
    1, 6, 11, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4,
    4, 0, 11, 4, 11, 5, 0, 3, 11, 6, 11, 7, 3, 7, 11,
    7, 6, 11, 7, 11, 8, 5, 4, 11, 4, 8, 11, -1, -1, -1,
    6, 9, 5, 6, 10, 9, 10, 8, 9, -1, -1, -1, -1, -1, -1,
    3, 6, 10, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1,
    0, 10, 8, 0, 5, 10, 0, 1, 5, 5, 6, 10, -1, -1, -1,
    6, 10, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1,
    1, 2, 11, 9, 5, 10, 9, 10, 8, 10, 5, 6, -1, -1, -1,
    0, 10, 3, 0, 6, 10, 0, 9, 6, 5, 6, 9, 1, 2, 11,
    10, 8, 5, 10, 5, 6, 8, 0, 5, 11, 5, 2, 0, 2, 5,
    6, 10, 3, 6, 3, 5, 2, 11, 3, 11, 5, 3, -1, -1, -1,
    5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1,
    9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1,
    1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8,
    1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    1, 3, 6, 1, 6, 11, 3, 8, 6, 5, 6, 9, 8, 9, 6,
    11, 1, 0, 11, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1,
    0, 3, 8, 5, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    11, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    10, 5, 11, 7, 5, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    10, 5, 11, 10, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1,
    5, 10, 7, 5, 11, 10, 1, 9, 0, -1, -1, -1, -1, -1, -1,
    11, 7, 5, 11, 10, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1,
    10, 1, 2, 10, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1,
    0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 10, -1, -1, -1,
    9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 10, 7, -1, -1, -1,
    7, 5, 2, 7, 2, 10, 5, 9, 2, 3, 2, 8, 9, 8, 2,
    2, 5, 11, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1,
    8, 2, 0, 8, 5, 2, 8, 7, 5, 11, 2, 5, -1, -1, -1,
    9, 0, 1, 5, 11, 3, 5, 3, 7, 3, 11, 2, -1, -1, -1,
    9, 8, 2, 9, 2, 1, 8, 7, 2, 11, 2, 5, 7, 5, 2,
    1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1,
    9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1,
    9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    5, 8, 4, 5, 11, 8, 11, 10, 8, -1, -1, -1, -1, -1, -1,
    5, 0, 4, 5, 10, 0, 5, 11, 10, 10, 3, 0, -1, -1, -1,
    0, 1, 9, 8, 4, 11, 8, 11, 10, 11, 4, 5, -1, -1, -1,
    11, 10, 4, 11, 4, 5, 10, 3, 4, 9, 4, 1, 3, 1, 4,
    2, 5, 1, 2, 8, 5, 2, 10, 8, 4, 5, 8, -1, -1, -1,
    0, 4, 10, 0, 10, 3, 4, 5, 10, 2, 10, 1, 5, 1, 10,
    0, 2, 5, 0, 5, 9, 2, 10, 5, 4, 5, 8, 10, 8, 5,
    9, 4, 5, 2, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 5, 11, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1,
    5, 11, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1,
    3, 11, 2, 3, 5, 11, 3, 8, 5, 4, 5, 8, 0, 1, 9,
    5, 11, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1,
    8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1,
    0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1,
    9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 10, 7, 4, 9, 10, 9, 11, 10, -1, -1, -1, -1, -1, -1,
    0, 8, 3, 4, 9, 7, 9, 10, 7, 9, 11, 10, -1, -1, -1,
    1, 11, 10, 1, 10, 4, 1, 4, 0, 7, 4, 10, -1, -1, -1,
    3, 1, 4, 3, 4, 8, 1, 11, 4, 7, 4, 10, 11, 10, 4,
    4, 10, 7, 9, 10, 4, 9, 2, 10, 9, 1, 2, -1, -1, -1,
    9, 7, 4, 9, 10, 7, 9, 1, 10, 2, 10, 1, 0, 8, 3,
    10, 7, 4, 10, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1,
    10, 7, 4, 10, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1,
    2, 9, 11, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1,
    9, 11, 7, 9, 7, 4, 11, 2, 7, 8, 7, 0, 2, 0, 7,
    3, 7, 11, 3, 11, 2, 7, 4, 11, 1, 11, 0, 4, 0, 11,
    1, 11, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1,
    4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1,
    4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    9, 11, 8, 11, 10, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    3, 0, 9, 3, 9, 10, 10, 9, 11, -1, -1, -1, -1, -1, -1,
    0, 1, 11, 0, 11, 8, 8, 11, 10, -1, -1, -1, -1, -1, -1,
    3, 1, 11, 10, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    1, 2, 10, 1, 10, 9, 9, 10, 8, -1, -1, -1, -1, -1, -1,
    3, 0, 9, 3, 9, 10, 1, 2, 9, 2, 10, 9, -1, -1, -1,
    0, 2, 10, 8, 0, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    3, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 8, 2, 8, 11, 11, 8, 9, -1, -1, -1, -1, -1, -1,
    9, 11, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 8, 2, 8, 11, 0, 1, 8, 1, 11, 8, -1, -1, -1,
    1, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};


#define MC_V1 1
#define MC_V2 2
#define MC_V3 4
#define MC_V4 8
#define MC_V5 16
#define MC_V6 32
#define MC_V7 64
#define MC_V8 128


uint8_t Geometry::mc_case( int32_t i, int32_t j, int32_t k ) const
{
#ifdef MC_DEBUG
    std::cout << "mc_case( " << i << ", " << j << ", " << k << " ) ";
#endif

    // Go through mesh nodes surrounding cube (i,j,k)
    uint8_t res = 0;
    uint32_t node;
    uint32_t ptr = (k*_size[1] + j)*_size[0] + i;

    // Node 1 (i,j,k)
    node = _smesh[ptr];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += MC_V1;

    // Node 2 (i+1,j,k)
    node = _smesh[ptr+1];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += MC_V2;

    // Node 3 (i+1,j+1,k)
    node = _smesh[ptr+1+_size[0]];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += MC_V3;

    // Node 4 (i,j+1,k)
    node = _smesh[ptr+_size[0]];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += MC_V4;

    // Second z-level
    ptr += _size[0]*_size[1];

    // Node 5 (i,j,k+1)
    node = _smesh[ptr];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += MC_V5;

    // Node 6 (i+1,j,k+1)
    node = _smesh[ptr+1];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += MC_V6;

    // Node 7 (i+1,j+1,k+1)
    node = _smesh[ptr+1+_size[0]];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += MC_V7;

    // Node 8 (i,j+1,k+1)
    node = _smesh[ptr+_size[0]];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += MC_V8;

#ifdef MC_DEBUG
    std::cout << "res = " << res << "\n";
#endif

    return( res );
}


Vec3D Geometry::mc_surface( int32_t i, int32_t j, int32_t k, uint8_t cn, int32_t ei ) const
{
#ifdef MC_DEBUG
    std::cout << "mc_surface( " << i << ", " << j << ", " << k << ", " << cn << ", " << ei << " )\n";
#endif
    if( ei == 0 ) {
	if( cn & MC_V1 ) {
	    double dist = solid_dist( i+1, j, k, 0 )/255.0;
	    return( Vec3D( (i+1-dist)*_h+_origo[0], j*_h+_origo[1], k*_h+_origo[2] ) );
	} else {
	    double dist = solid_dist( i, j, k, 1 )/255.0;
	    return( Vec3D( (i+dist)*_h+_origo[0], j*_h+_origo[1], k*_h+_origo[2] ) );
	}
    } else if( ei == 1 ) {
	if( cn & MC_V2 ) {
	    double dist = solid_dist( i+1, j+1, k, 2 )/255.0;
	    return( Vec3D( (i+1)*_h+_origo[0], (j+1-dist)*_h+_origo[1], k*_h+_origo[2] ) );
	} else {
	    double dist = solid_dist( i+1, j, k, 3 )/255.0;
	    return( Vec3D( (i+1)*_h+_origo[0], (j+dist)*_h+_origo[1], k*_h+_origo[2] ) );
	}
    } else if( ei == 2 ) {
	if( cn & MC_V3 ) {
	    double dist = solid_dist( i, j+1, k, 1 )/255.0;
	    return( Vec3D( (i+dist)*_h+_origo[0], (j+1)*_h+_origo[1], k*_h+_origo[2] ) );
	} else {
	    double dist = solid_dist( i+1, j+1, k, 0 )/255.0;
	    return( Vec3D( (i+1-dist)*_h+_origo[0], (j+1)*_h+_origo[1], k*_h+_origo[2] ) );
	}
    } else if( ei == 3 ) {
	if( cn & MC_V1 ) {
	    double dist = solid_dist( i, j+1, k, 2 )/255.0;
	    return( Vec3D( i*_h+_origo[0], (j+1-dist)*_h+_origo[1], k*_h+_origo[2] ) );
	} else {
	    double dist = solid_dist( i, j, k, 3 )/255.0;
	    return( Vec3D( i*_h+_origo[0], (j+dist)*_h+_origo[1], k*_h+_origo[2] ) );
	}
    } else if( ei == 4 ) {
	if( cn & MC_V5 ) {
	    double dist = solid_dist( i+1, j, k+1, 0 )/255.0;
	    return( Vec3D( (i+1-dist)*_h+_origo[0], j*_h+_origo[1], (k+1)*_h+_origo[2] ) );
	} else {
	    double dist = solid_dist( i, j, k+1, 1 )/255.0;
	    return( Vec3D( (i+dist)*_h+_origo[0], j*_h+_origo[1], (k+1)*_h+_origo[2] ) );
	}
    } else if( ei == 5 ) {
	if( cn & MC_V6 ) {
	    double dist = solid_dist( i+1, j+1, k+1, 2 )/255.0;
	    return( Vec3D( (i+1)*_h+_origo[0], (j+1-dist)*_h+_origo[1], (k+1)*_h+_origo[2] ) );
	} else {
	    double dist = solid_dist( i+1, j, k+1, 3 )/255.0;
	    return( Vec3D( (i+1)*_h+_origo[0], (j+dist)*_h+_origo[1], (k+1)*_h+_origo[2] ) );
	}
    } else if( ei == 6 ) {
	if( cn & MC_V7 ) {
	    double dist = solid_dist( i, j+1, k+1, 1 )/255.0;
	    return( Vec3D( (i+dist)*_h+_origo[0], (j+1)*_h+_origo[1], (k+1)*_h+_origo[2] ) );
	} else {
	    double dist = solid_dist( i+1, j+1, k+1, 0 )/255.0;
	    return( Vec3D( (i+1-dist)*_h+_origo[0], (j+1)*_h+_origo[1], (k+1)*_h+_origo[2] ) );
	}
    } else if( ei == 7 ) {
	if( cn & MC_V5 ) {
	    double dist = solid_dist( i, j+1, k+1, 2 )/255.0;
	    return( Vec3D( i*_h+_origo[0], (j+1-dist)*_h+_origo[1], (k+1)*_h+_origo[2] ) );
	} else {
	    double dist = solid_dist( i, j, k+1, 3 )/255.0;
	    return( Vec3D( i*_h+_origo[0], (j+dist)*_h+_origo[1], (k+1)*_h+_origo[2] ) );
	}
    } else if( ei == 8 ) {
	if( cn & MC_V1 ) {
	    double dist = solid_dist( i, j, k+1, 4 )/255.0;
	    return( Vec3D( i*_h+_origo[0], j*_h+_origo[1], (k+1-dist)*_h+_origo[2] ) );
	} else {
	    double dist = solid_dist( i, j, k, 5 )/255.0;
	    return( Vec3D( i*_h+_origo[0], j*_h+_origo[1], (k+dist)*_h+_origo[2] ) );
	}
    } else if( ei == 9 ) {
	if( cn & MC_V2 ) {
	    double dist = solid_dist( i+1, j, k+1, 4 )/255.0;
	    return( Vec3D( (i+1)*_h+_origo[0], j*_h+_origo[1], (k+1-dist)*_h+_origo[2] ) );
	} else {
	    double dist = solid_dist( i+1, j, k, 5 )/255.0;
	    return( Vec3D( (i+1)*_h+_origo[0], j*_h+_origo[1], (k+dist)*_h+_origo[2] ) );
	}
    } else if( ei == 10 ) {
	if( cn & MC_V4 ) {
	    double dist = solid_dist( i, j+1, k+1, 4 )/255.0;
	    return( Vec3D( i*_h+_origo[0], (j+1)*_h+_origo[1], (k+1-dist)*_h+_origo[2] ) );
	} else {
	    double dist = solid_dist( i, j+1, k, 5 )/255.0;
	    return( Vec3D( i*_h+_origo[0], (j+1)*_h+_origo[1], (k+dist)*_h+_origo[2] ) );
	}
    } else if( ei == 11 ) {
	if( cn & MC_V3 ) {
	    double dist = solid_dist( i+1, j+1, k+1, 4 )/255.0;
	    return( Vec3D( (i+1)*_h+_origo[0], (j+1)*_h+_origo[1], (k+1-dist)*_h+_origo[2] ) );
	} else {
	    double dist = solid_dist( i+1, j+1, k, 5 )/255.0;
	    return( Vec3D( (i+1)*_h+_origo[0], (j+1)*_h+_origo[1], (k+dist)*_h+_origo[2] ) );
	}
    } else {
	throw( Error( ERROR_LOCATION, "algorithm error: unknown edge index" ) );
    }
}


uint32_t Geometry::mc_trianglec( uint8_t cn ) const
{
    int32_t offset = 15*cn;
    int32_t ti;
    for( ti = 0; ti < 5; ti++ ) {
	if( mc_faces[offset] == -1 )
	    break;
	offset += 3;
    }
    return( ti );
}


int32_t Geometry::mc_add_vertex_try( const Vec3D &x, int32_t i, int32_t j, int32_t k ) const
{
    if( i < 0 || j < 0 || k < 0 )
	return( -1 );

    uint8_t cn = mc_case( i, j, k );
    uint32_t tricount = mc_trianglec( cn );
    uint32_t triptr = surface_triangle_ptr( i, j, k );

    for( uint32_t a = 0; a < tricount; a++ ) {
	const VTriangle &tri = _surface.triangle(triptr+a);
	for( uint32_t b = 0; b < 3; b++ ) {
	    const Vec3D &y = _surface.vertex(tri[b]);
	    if( ssqr(x-y) < _surface_eps )
		return( tri[b] );
	}
    }

    return( -1 );
}


uint32_t Geometry::mc_add_vertex( const Vec3D &x, int32_t i, int32_t j, int32_t k, uint32_t firstv )
{
    // Check if vertex x exists in this cell
    for( uint32_t a = firstv; a < _surface.vertexc(); a++ ) {
	const Vec3D &y = _surface.vertex(a);
	if( ssqr(x-y) < _surface_eps )
	    return( a );
    }

    // Try if vertex x exists in neighbours
    int32_t v;
    if( (v=mc_add_vertex_try( x, i-1, j, k )) != -1 )
	return( v );
    if( (v=mc_add_vertex_try( x, i, j-1, k )) != -1 )
	return( v );
    if( (v=mc_add_vertex_try( x, i, j, k-1 )) != -1 )
	return( v );

    // New vertex
    return( _surface.add_vertex_no_check( x ) );
}


void Geometry::mc_triangulate( int32_t i, int32_t j, int32_t k )
{
#ifdef MC_DEBUG
    std::cout << "mc_triangulate( " << i << ", " << j << ", " << k << " )\n";
#endif
    uint8_t cn = mc_case( i, j, k );
    int offset = 15*cn;

#ifdef MC_DEBUG
    for( int ti = 0; ti < 5; ti++ ) {
	if( mc_faces[offset+3*ti] == -1 )
	    break;
	for( int a = 0; a < 3; a++ )
	    std::cout << "  triangle[" << ti << "] = " << mc_faces[offset+3*ti+a] << "\n";
    }
#endif

    // Go through (possible) five triangles
    uint32_t firstv = _surface.vertexc();
    for( int ti = 0; ti < 5; ti++ ) {

	// If all triangles processed
	if( mc_faces[offset] == -1 )
	    break;

	// Add vertices
	uint32_t v[3];
	for( int a = 0; a < 3; a++ ) {
	    int32_t ei = mc_faces[offset+a];
	    Vec3D x = mc_surface( i, j, k, cn, ei );
	    v[a] = mc_add_vertex( x, i, j, k, firstv );
	}
	
	// Add triangle
	_surface.add_triangle( v );
	offset += 3;
    }
}


void Geometry::build_surface( void )
{
    Timer t;
    ibsimu.message( 1 ) << "Building surface triangulation\n";
    ibsimu.inc_indent();

    if( _geom_mode != MODE_3D )
	throw( Error( ERROR_LOCATION, "incompatible geometry type" ) );
    if( !_built )
	throw( Error( ERROR_LOCATION, "solid mesh not built" ) );

    // Clear old data
    _surface.clear();
    _triptr.clear();

    // Reserve cube to triangle pointer array
    Int3D csize( _size[0]-1, _size[1]-1, _size[2]-1 );
    _triptr.reserve( csize[0]*csize[1]*csize[2] );

    // Precalculate vertex matching tolerance
    _surface_eps = (_h/512.0)*(_h/512.0);

    // Go though all cubes
    for( int32_t k = 0; k < csize[2]; k++ ) {
	for( int32_t j = 0; j < csize[1]; j++ ) {
	    for( int32_t i = 0; i < csize[0]; i++ ) {

		// Set triangle pointer
		_triptr.push_back( _surface.trianglec() );

		// Construct triangulation for (i,j,k)
		mc_triangulate( i, j, k );
	    }
	}
    }

    t.stop();
    ibsimu.message( 1 ) << _surface.trianglec() << " triangles\n";
    ibsimu.message( 1 ) << _surface.vertexc() << " vertices\n";
    ibsimu.message( 1 ) << "Done.\n";
    ibsimu.message( 1 ) << "time used = " << t << "\n";
    ibsimu.dec_indent();
}


uint8_t Geometry::surface_cell_face_case_2d( const int32_t i[3], const int32_t vb[3] ) const
{
    int res = 0;
    uint32_t node;
    int32_t j[3] = { i[0], i[1], i[2] };

    // Node 1 (x,y)
    node = mesh( (j[2]*_size[1] + j[1])*_size[0] + j[0] );
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += 1;

    // Node 2 (x+1,y)
    j[vb[0]]++;
    node = mesh( (j[2]*_size[1] + j[1])*_size[0] + j[0] );
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += 2;

    // Node 3 (x+1,y+1)
    j[vb[1]]++;
    node = mesh( (j[2]*_size[1] + j[1])*_size[0] + j[0] );
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += 4;

    // Node 4 (x,y+1)
    j[vb[0]]--;
    node = mesh( (j[2]*_size[1] + j[1])*_size[0] + j[0] );
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	res += 8;

    return( res );
}


double Geometry::surface_cell_face_dist( const int32_t i[3], const int32_t vb[3], 
					 int32_t dx, int32_t dy, 
					 int32_t dir ) const
{
    int32_t j[3] = { i[0], i[1], i[2] };
    j[vb[0]] += dx;
    j[vb[1]] += dy;
    return( solid_dist( j[0], j[1], j[2], dir )/255.0 );
}


uint32_t Geometry::surface_cell_face_add_vertex( VTriangleSurfaceSolid &solid,
						 const int32_t i[3], const int32_t vb[3], 
						 double dx, double dy ) const
{
    double u[3] = { (double)i[0], (double)i[1], (double)i[2] };
    u[vb[0]] += dx;
    u[vb[1]] += dy;

    Vec3D x( _origo[0]+_h*u[0], _origo[1]+_h*u[1], _origo[2]+_h*u[2] );
    return( solid.add_vertex( x ) );
}


void Geometry::surface_cell_face_add_triangles( VTriangleSurfaceSolid &solid,
						int32_t i0, int32_t i1, int32_t i2, 
						int32_t vb0, int32_t vb1, int32_t vb2 ) const 
{
    int32_t i[3] = { i0, i1, i2 };
    int32_t vb[3] = { vb0, vb1, vb2 };
    surface_cell_face_add_triangles( solid, i, vb );
}


void Geometry::surface_cell_face_add_triangles( VTriangleSurfaceSolid &solid,
					       const int32_t i[3], const int32_t vb[3] ) const
{
    uint8_t cn = surface_cell_face_case_2d( i, vb );
    double dist;
    double dist2;
    uint32_t v[3];
    switch( cn ) {
    case 0:
	break;
    case 1:
	dist = surface_cell_face_dist( i, vb, 1, 0, 2*vb[0] );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 1.0-dist, 0.0 );
	dist = surface_cell_face_dist( i, vb, 0, 1, 2*vb[1] );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0-dist );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 0.0 );
	solid.add_triangle( v );
	break;
    case 2:
	dist = surface_cell_face_dist( i, vb, 0, 0, 2*vb[0]+1 );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, dist, 0.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 0.0 );
	dist = surface_cell_face_dist( i, vb, 1, 1, 2*vb[1] );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0-dist );
	solid.add_triangle( v );
	break;
    case 3:
	dist = surface_cell_face_dist( i, vb, 0, 1, 2*vb[1] );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0-dist );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 0.0 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 0.0 );
	solid.add_triangle( v );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0-dist );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 0.0 );
	dist = surface_cell_face_dist( i, vb, 1, 1, 2*vb[1] );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0-dist );
	solid.add_triangle( v );
	break;
    case 4:
	dist = surface_cell_face_dist( i, vb, 0, 1, 2*vb[0]+1 );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, dist, 1.0 );
	dist = surface_cell_face_dist( i, vb, 1, 0, 2*vb[1]+1 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0, dist );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0 );
	solid.add_triangle( v );
	break;
    case 5:
	dist = surface_cell_face_dist( i, vb, 0, 1, 2*vb[1] );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0-dist );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 0.0 );
	dist = surface_cell_face_dist( i, vb, 1, 0, 2*vb[0] );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0-dist, 0.0 );
	solid.add_triangle( v );
	dist = surface_cell_face_dist( i, vb, 0, 1, 2*vb[0]+1 );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, dist, 1.0 );
	dist = surface_cell_face_dist( i, vb, 1, 0, 2*vb[1]+1 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0, dist );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0 );
	solid.add_triangle( v );
	break;
    case 6:
	dist = surface_cell_face_dist( i, vb, 0, 1, 2*vb[0]+1 );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, dist, 1.0 );
	dist = surface_cell_face_dist( i, vb, 0, 0, 2*vb[0]+1 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, dist, 0.0 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0 );
	solid.add_triangle( v );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, dist, 0.0 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 0.0 );
	solid.add_triangle( v );
	break;
    case 7:
	dist = surface_cell_face_dist( i, vb, 0, 1, 2*vb[1] );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0-dist );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 0.0 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 0.0 );
	solid.add_triangle( v );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 0.0 );
	dist2 = surface_cell_face_dist( i, vb, 0, 1, 2*vb[1] );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, dist2, 1.0 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0-dist );
	solid.add_triangle( v );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 0.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, dist2, 1.0 );
	solid.add_triangle( v );
	break;
    case 8:
	dist = surface_cell_face_dist( i, vb, 0, 0, 2*vb[1]+1 );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, dist );
	dist = surface_cell_face_dist( i, vb, 1, 1, 2*vb[0] );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0-dist, 1.0 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0 );
	solid.add_triangle( v );
	break;
    case 9:
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 0.0 );
	dist = surface_cell_face_dist( i, vb, 1, 1, 2*vb[0] );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0-dist, 1.0 );
	solid.add_triangle( v );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 1.0-dist, 1.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 0.0 );
	dist = surface_cell_face_dist( i, vb, 1, 0, 2*vb[0] );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0-dist, 0.0 );
	solid.add_triangle( v );
	break;
    case 10:
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0 );
	dist = surface_cell_face_dist( i, vb, 0, 0, 2*vb[1]+1 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 0.0, dist );
	dist = surface_cell_face_dist( i, vb, 1, 1, 2*vb[0] );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0-dist, 1.0 );
	solid.add_triangle( v );
	dist = surface_cell_face_dist( i, vb, 0, 0, 2*vb[0]+1 );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, dist, 0.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 0.0 );
	dist = surface_cell_face_dist( i, vb, 1, 1, 2*vb[1] );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0-dist );
	solid.add_triangle( v );
	break;
    case 11:
	dist = surface_cell_face_dist( i, vb, 1, 1, 2*vb[0] );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 1.0-dist, 1.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 0.0 );
	solid.add_triangle( v );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 0.0 );
	dist2 = surface_cell_face_dist( i, vb, 1, 1, 2*vb[1] );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0-dist2 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0-dist, 1.0 );
	solid.add_triangle( v );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 0.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 0.0 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0-dist2 );
	solid.add_triangle( v );
	break;
    case 12:
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0 );
	dist = surface_cell_face_dist( i, vb, 0, 0, 2*vb[1]+1 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 0.0, dist );
	dist = surface_cell_face_dist( i, vb, 1, 0, 2*vb[1]+1 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, dist );
	solid.add_triangle( v );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 1.0, dist );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0 );
	solid.add_triangle( v );
	break;
    case 13:
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 0.0 );
	dist = surface_cell_face_dist( i, vb, 1, 0, 2*vb[0] );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0-dist, 0.0 );
	solid.add_triangle( v );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 1.0-dist, 0.0 );
	dist2 = surface_cell_face_dist( i, vb, 1, 0, 2*vb[1]+1 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0, dist2 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0 );
	solid.add_triangle( v );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0, dist2 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0 );
	solid.add_triangle( v );
	break;
    case 14:
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0 );
	dist = surface_cell_face_dist( i, vb, 0, 0, 2*vb[1]+1 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 0.0, dist );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0 );
	solid.add_triangle( v );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 0.0, dist );
	dist2 = surface_cell_face_dist( i, vb, 0, 0, 2*vb[0]+1 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, dist2, 0.0 );
	solid.add_triangle( v );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, dist2, 0.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 0.0 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0 );
	solid.add_triangle( v );
	break;
    case 15:
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 0.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 0.0 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0 );
	solid.add_triangle( v );
	v[0] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 1.0 );
	v[1] = surface_cell_face_add_vertex( solid, i, vb, 0.0, 1.0 );
	v[2] = surface_cell_face_add_vertex( solid, i, vb, 1.0, 0.0 );
	solid.add_triangle( v );
	break;
    }
}


uint32_t Geometry::surface_inside( const Vec3D &x ) const
{
    // Grid cell
    int32_t i[3];
    for( int a = 0; a < 3; a++ ) {
	i[a] = floor( (x[a]-_origo[a])*_div_h );
	if( i[a] < 0 )
	    return( 2*a+1 );
	else if( i[a] >= _size[a]-1 )
	    return( 2*a+2 );
    }

    //ibsimu.message(1) << "surface_inside()\n";
    //ibsimu.message(1) << "x = " << x << "\n";
    //ibsimu.message(1) << "(i,j,k) = (" << i[0] << ", " << i[1] << ", " << i[2] << ")\n";

    // Check if point is trivially inside or outside according to
    // solid mesh (case numbers 0 and 255).
    uint8_t cn = mc_case( i[0], i[1], i[2] );
    uint32_t ptr = (i[2]*_size[1] + i[1])*_size[0] + i[0];
    uint32_t node = _smesh[ptr];

    //ibsimu.message(1) << "cn = " << (int)cn << "\n";
    
    if( cn == 0 )
	return( 0 );
    else if( cn == 255 )
	return( node & SMESH_BOUNDARY_NUMBER_MASK );

    // Build triangulated boundary representation of the solid in the
    // grid cell to resolve inclusion.
    VTriangleSurfaceSolid solid;
    double dx = _h/512.0;
    solid.set_vertex_matching_eps( dx );
    solid.set_signed_volume_eps( dx*dx*dx );

    // Go through (possible) five triangles
    int offset = 15*cn;
    for( int ti = 0; ti < 5; ti++ ) {

	// If all triangles processed
	if( mc_faces[offset] == -1 )
	    break;

	// Add vertices
	uint32_t v[3];
	for( int a = 0; a < 3; a++ ) {
	    int32_t ei = mc_faces[offset+a];
	    v[2-a] = solid.add_vertex( mc_surface( i[0], i[1], i[2], cn, ei ) );
	}
	solid.add_triangle( v );

	offset += 3;
    }

    // Add bounds of grid cell on six faces: -x, +x, -y, +y, -z, +z
    surface_cell_face_add_triangles( solid, i[0],   i[1],   i[2],   2, 1, 0 );
    surface_cell_face_add_triangles( solid, i[0]+1, i[1],   i[2],   1, 2, 0 );
    surface_cell_face_add_triangles( solid, i[0],   i[1],   i[2],   0, 2, 1 );
    surface_cell_face_add_triangles( solid, i[0],   i[1]+1, i[2],   2, 0, 1 );
    surface_cell_face_add_triangles( solid, i[0],   i[1],   i[2],   1, 0, 2 );
    surface_cell_face_add_triangles( solid, i[0],   i[1],   i[2]+1, 0, 1, 2 );

    solid.prepare_for_inside();

    //solid.check_data();
    //solid.debug_print( ibsimu.message(1) );
    
    // Shouldn't test close to edges to avoid errors on faces inside
    // solids. Isn't critical on faces at least partially outside.
    Vec3D y = x;
    double dd = _h/256;
    for( int a = 0; a < 3; a++ ) {
	double s = y[a]-(_origo[a]+_h*i[a]);
	if( fabs(s) < dd ) {
	    //ibsimu.message(1) << "close to face min, adjusting\n";
	    y[a] += dd;
	} else if( fabs(s-_h) < dd ) {
	    //ibsimu.message(1) << "close to face max, adjusting\n";
	    y[a] -= dd;
	}
    }

    if( solid.inside( y ) ) {
	//ibsimu.message(1) << "inside\n";
	return( surface_inside_solid_number(i[0],i[1],i[2]) );
    }

    //ibsimu.message(1) << "not inside\n";
    return( 0 );
}


uint32_t Geometry::surface_inside_solid_number( int32_t i, int32_t j,int32_t k ) const
{
    uint32_t ptr = (k*_size[1] + j)*_size[0] + i;
    uint32_t node = _smesh[ptr];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	return( node & SMESH_BOUNDARY_NUMBER_MASK );    
    node = _smesh[ptr+1];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	return( node & SMESH_BOUNDARY_NUMBER_MASK );
    node = _smesh[ptr+_size[0]];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	return( node & SMESH_BOUNDARY_NUMBER_MASK );
    node = _smesh[ptr+_size[0]+1];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	return( node & SMESH_BOUNDARY_NUMBER_MASK );

    ptr += _size[1]*_size[0];
    node = _smesh[ptr];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	return( node & SMESH_BOUNDARY_NUMBER_MASK );    
    node = _smesh[ptr+1];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	return( node & SMESH_BOUNDARY_NUMBER_MASK );
    node = _smesh[ptr+_size[0]];
    if( (node & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(node & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	return( node & SMESH_BOUNDARY_NUMBER_MASK );
    node = _smesh[ptr+_size[0]+1];
    return( node & SMESH_BOUNDARY_NUMBER_MASK );
}


Vec3D Geometry::surface_triangle_normal( const Vec3D &x ) const
{
    const double eps = 1e-6;

    // Calculate mesh cell indices i which encloses x and
    // Parametric position s.
    int32_t i[3];
    double s[3];
    for( int32_t k = 0; k < 3; k++ ) {
	i[k] = floor((x[k] - _origo[k])*_div_h);
	s[k] = x[k] - (_origo[k] + i[k]*_h);
    }

    // Test cel i and neighbouring cells if close to edge.
    // Start with dir 0, which is cell i, dir 1 is i-1 and dir
    // 2 is i+1. Same for all 3 coordinate directions
    for( int32_t k = 0; k < 3; k++ ) {
	for( int32_t dir = 0; dir < 3; dir++ ) {

	    if( dir == 0 && k != 0 ) {
		continue;
	    } else if( dir == 1 ) {
		if( s[k] > eps || i[k] == 0 )
		    continue;
		i[k] -= 1;
	    } else if( dir == 2 ) {
		if( s[k] < 1-eps || i[k] >= _size[k]-3 )
		    continue;
		i[k] += 1;
	    }

	    // Get pointer and count of triangles in the cell
	    uint32_t a = surface_triangle_ptr( i[0], i[1], i[2] );
	    uint32_t c = surface_trianglec( i[0], i[1], i[2] );

	    // Go through the triangles to find the correct one
	    for( uint32_t p = 0; p < c; p++ ) {
		const VTriangle &tri = surface_triangle(a+p);
		const Vec3D &v1 = _surface.vertex(tri[0]);
		const Vec3D &v2 = _surface.vertex(tri[1]);
		const Vec3D &v3 = _surface.vertex(tri[2]);
	
		Vec3D e1 = v2-v1;
		Vec3D e2 = v3-v1;
		Vec3D r = x-v1;
		Vec3D e3 = cross(e2,e1);
		e3.normalize();

		double dot11 = e1*e1;
		double dot12 = e1*e2;
		double dot1r = e1*r;
		double dot22 = e2*e2;
		double dot2r = e2*r;
		double dot3r = e3*r;

		double invdet = 1.0/(dot11*dot22 - dot12*dot12);
		double v = (dot22*dot1r - dot12*dot2r) * invdet;
		double w = (dot11*dot2r - dot12*dot1r) * invdet;

		// Check if point x is within triangle, allowing eps tolerance
		if( v >= -eps && v <= 1+eps && w >= -eps && w <= 1+eps &&
		    v+w <= 1+eps && fabs(dot3r) < eps*_h ) {
		    // Return normal of the triangle
		    return( e3 );
		}
	    }

	    if( dir == 1 ) i[k] += 1;
	    else if( dir == 2 ) i[k] -= 1;
	}    
    }

    // No triangle found at x, return NaN vector
    return( Vec3D(std::numeric_limits<double>::quiet_NaN()) );
}


Vec3D Geometry::surface_triangle_normal( int32_t a ) const
{
    const VTriangle &tri = _surface.triangle(a);
    const Vec3D &v1 = _surface.vertex(tri[0]);
    const Vec3D &v2 = _surface.vertex(tri[1]);
    const Vec3D &v3 = _surface.vertex(tri[2]);

    Vec3D e1 = v2-v1;
    Vec3D e2 = v3-v1;
    Vec3D e3 =  cross(e2,e1);
    e3.normalize();
    return( e3 );
}


int32_t Geometry::surface_trianglec( int32_t i, int32_t j, int32_t k ) const
{
    uint8_t cn = mc_case( i, j, k );
    int32_t offset = 15*cn;
    int32_t ti;
    for( ti = 0; ti < 5; ti++ ) {
	if( mc_faces[offset] == -1 )
	    break;
	offset += 3;
    }
    return( ti );
}


uint8_t Geometry::solid_dist( uint32_t i, uint32_t j, uint32_t k, uint32_t dir ) const
{
    uint32_t snode = _smesh[i + j*_size[0] + k*_size[0]*_size[1]];
    if( (snode & SMESH_NODE_ID_MASK) != SMESH_NODE_ID_NEAR_SOLID )
	throw( Error( ERROR_LOCATION, "not a near solid node" ) );

    const uint8_t *nptr = &_nearsolid[snode & SMESH_NEAR_SOLID_INDEX_MASK];
    uint8_t neighbours = nptr[0];
    nptr++;
    uint32_t a = 0;
    while( a < dir ) {
	if( neighbours & 0x01 )
	    nptr++;
	neighbours = neighbours >> 1;
	a++;
    }
    if( (neighbours & 0x01) == 0x00 )
	throw( Error( ERROR_LOCATION, (const std::string)"no near neighbour in selected direction"
		      ", dir = " + to_string(dir) 
		      + ", neighbours = " + to_string((int)_nearsolid[snode & SMESH_NEAR_SOLID_INDEX_MASK]) ) );

    return( *nptr );
}


uint8_t Geometry::solid_dist( uint32_t i, uint32_t dir ) const
{
    uint32_t snode = _smesh[i];
    if( (snode & SMESH_NODE_ID_MASK) != SMESH_NODE_ID_NEAR_SOLID )
	throw( Error( ERROR_LOCATION, "not a near solid node" ) );

    const uint8_t *nptr = &_nearsolid[snode & SMESH_NEAR_SOLID_INDEX_MASK];
    uint8_t neighbours = nptr[0];
    nptr++;
    uint32_t a = 0;
    while( a < dir ) {
	if( neighbours & 0x01 )
	    nptr++;
	neighbours = neighbours >> 1;
	a++;
    }
    if( (neighbours & 0x01) == 0x00 )
	throw( Error( ERROR_LOCATION, (const std::string)"no near neighbour in selected direction"
		      ", dir = " + to_string(dir) 
		      + ", neighbours = " + to_string((int)_nearsolid[snode & SMESH_NEAR_SOLID_INDEX_MASK]) ) );

    return( *nptr );
}


void Geometry::save( const std::string &filename, bool save_solids ) const
{
    ibsimu.message( 1 ) << "Saving Geometry to file \'" << filename << "\'.\n";

    std::ofstream os( filename.c_str(), std::ios_base::binary );
    if( !os.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\' for writing" ) );
    save( os, save_solids );
    os.close();
}


void Geometry::save( std::ostream &os, bool save_solids ) const
{
    Mesh::save( os );
    write_int32( os, _n );
    for( uint32_t a = 0; a < _n; a++ ) {
	if( save_solids )
	    _sdata[a]->save( os );
	else
	    write_int32( os, FILEID_NULL );
    }
    for( uint32_t a = 0; a < _n+6; a++ )
	_bound[a].save( os );

    write_int8( os, _built );
    write_compressed_block( os, _size[0]*_size[1]*_size[2]*sizeof(uint32_t), 
			    (int8_t *)_smesh );

    write_int32( os, _nearsolid.size() );
    write_compressed_block( os, _nearsolid.size()*sizeof(uint8_t), 
			    (int8_t *)&_nearsolid[0] );
}


void Geometry::debug_print( std::ostream &os ) const
{
    Mesh::debug_print( os );

    os << "**Geometry\n";

    os << "n = " << _n << "\n";
    if( _n == 0 )
	os << "no sdata\n";
    for( uint32_t a = 0; a < _n; a++ ) {
	os << "sdata[" << a << "]:\n";
	_sdata[a]->debug_print( os );
    }
    for( uint32_t a = 0; a < _n+6; a++ ) {
	os << "bound[" << a+1 << "] = " << _bound[a] << "\n";
    }
    os << "built = " << _built << "\n";
    if( (_geom_mode == MODE_1D || _geom_mode == MODE_2D || _geom_mode == MODE_CYL) && 
	_size[0] <= 20 && _size[1] <= 20 ) {
	os << "mesh visualization:\n";
	for( int32_t j = _size[1]-1; j >= 0; j-- ) {
	    for( int32_t i = 0; i < _size[0]; i++ ) {
		uint32_t ind = mesh(i,j);
		if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_PURE_VACUUM ) {
		    os << " ";
		} else if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEUMANN ) {
		    os << "N";
		} else if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET ) {
		    if( (ind & SMESH_BOUNDARY_NUMBER_MASK) <= 6 )
			os << "D";
		    else
			os << "S";
		} else if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEAR_SOLID ) {
		    os << ".";
		}
	    }
	    os << "\n";
	}
    }

    for( int32_t k = 0; k < _size[2]; k++ ) {
	for( int32_t j = 0; j < _size[1]; j++ ) {
	    for( int32_t i = 0; i < _size[0]; i++ ) {

		uint32_t ind = mesh(i,j,k);
		os << "smesh(" << i << ", " << j << ", " << k << ") = "
		   << "0x" << std::hex << std::setfill('0')
		   << std::setw(8) << ind << " ("
		   << std::dec << std::setfill(' ');

		if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_PURE_VACUUM ) {
		    os << "pure vacuum)\n";
		} else if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEUMANN ) {
		    os << "neumann, solid " << (ind & SMESH_BOUNDARY_NUMBER_MASK) << ")\n";
		} else if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET ) {
		    os << "dirichlet, solid " << (ind & SMESH_BOUNDARY_NUMBER_MASK) << ")\n";
		} else if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEAR_SOLID ) {		    
		    os << "near solid, index " << (ind & SMESH_NEAR_SOLID_INDEX_MASK) << ")\n";
		    ind = (ind & SMESH_NEAR_SOLID_INDEX_MASK);
		    uint8_t sflag = _nearsolid[ind];
		    os << std::hex << std::setfill('0')
		       << "  near solid flags = 0x" << std::setw(2) << (int)sflag << "\n"
		       << std::dec << std::setfill(' ');
		    uint8_t mask = 1;
		    uint8_t *ptr = (uint8_t *)&_nearsolid[ind+1];
		    for( uint32_t a = 0; a < 6; a++ ) {
			if( mask & sflag ) {
			    os << "  ndist[" << a << "] = " << (int)*ptr << "\n";
			    ptr++;
			}
			mask = mask << 1;
		    }
		}
	    }
	}
    }
}



