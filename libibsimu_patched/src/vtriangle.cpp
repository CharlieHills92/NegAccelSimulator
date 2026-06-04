/*! \file vtriangle.hpp
 *  \brief Vertex-based triangle representation
 */

/* Copyright (c) 2011-2013 Taneli Kalvas. All rights reserved.
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


#include <iomanip>
#include <limits>
#include "vtriangle.hpp"
#include "ibsimu.hpp"


VTriangle::VTriangle( uint32_t v1, uint32_t v2, uint32_t v3 )
{
    _v[0] = v1;
    _v[1] = v2;
    _v[2] = v3;
}


VTriangle::VTriangle( const uint32_t v[3] )
{
    _v[0] = v[0];
    _v[1] = v[1];
    _v[2] = v[2];
}


VTriangle::~VTriangle()
{
    
}


void VTriangle::debug_print( std::ostream &os ) const
{
    os << "**VTriangle\n";    
    os << "  v = "  
       << std::setw(6) << _v[0] << " "
       << std::setw(6) << _v[1] << " "
       << std::setw(6) << _v[2] << "\n";
}


/*
 * VTriangleSurface
 * ****************************************************** */


VTriangleSurface::VTriangleSurface( double vertex_matching_eps )
    : _vertex_matching_eps(vertex_matching_eps*vertex_matching_eps)
{

}


VTriangleSurface::~VTriangleSurface()
{

}


void VTriangleSurface::set_vertex_matching_eps( double vertex_matching_eps )
{
    _vertex_matching_eps = vertex_matching_eps*vertex_matching_eps;
}


void VTriangleSurface::clear( void )
{
    _vertex.clear();
    _triangle.clear();
}


uint32_t VTriangleSurface::add_vertex( const Vec3D &x )
{
    /* Search for vertex matching to x, allow equal to matching
     * epsilon to allow the use of zero as vertex_matching_eps.
    */
    uint32_t v;
    for( v = 0; v < _vertex.size(); v++ )
	if( ssqr(x-_vertex[v]) <= _vertex_matching_eps )
	    return( v ); // Vertex found

    // New vertex needed
    _vertex.push_back( x ); 
    return( _vertex.size()-1 );    
}


uint32_t VTriangleSurface::add_triangle( const Vec3D x[3] )
{
    uint32_t v[3];

    for( int a = 0; a < 3; a++ )
	v[a] = add_vertex( x[a] );

    _triangle.push_back( VTriangle( v ) );    
    return( _triangle.size()-1 );
}


uint32_t VTriangleSurface::add_triangle( const Vec3D &x1, const Vec3D &x2, const Vec3D &x3 )
{
    uint32_t v[3];

    v[0] = add_vertex( x1 );
    v[1] = add_vertex( x2 );
    v[2] = add_vertex( x3 );
    
    _triangle.push_back( VTriangle( v ) );
    return( _triangle.size()-1 );
}


void VTriangleSurface::debug_print( std::ostream &os ) const
{
    os << "**VTriangleSurface\n";

    os << "vertex_matching_eps = " << _vertex_matching_eps << "\n";
    os << "trianglec = " << _triangle.size() << "\n";
    os << "vertexc = " << _vertex.size() << "\n";

    for( uint32_t a = 0; a < _triangle.size(); a++ ) 
	os << "  triangle[" << a << "] = (" 
	   << _triangle[a][0] << ", "
	   << _triangle[a][1] << ", "
	   << _triangle[a][2] << ")\n";

    for( uint32_t a = 0; a < _vertex.size(); a++ ) 
	os << "  vertex[" << a << "] = " << _vertex[a] << "\n";
}


/*
 * VTriangleSurfaceSolid
 * ****************************************************** */


VTriangleSurfaceSolid::VTriangleSurfaceSolid( double vertex_matching_eps, 
					      double signed_volume_eps )
    : VTriangleSurface(vertex_matching_eps), 
      _signed_volume_eps(signed_volume_eps)
{

}


VTriangleSurfaceSolid::~VTriangleSurfaceSolid()
{

}


#define VTRI_EDGE1   0
#define VTRI_EDGE2   1
#define VTRI_EDGE3   2
#define VTRI_FACE    3
#define VTRI_OUTSIDE 4
#define VTRI_INSIDE  5


int VTriangleSurfaceSolid::signvol4( const Vec3D &q0, const Vec3D &q1, 
				     const Vec3D &q2, const Vec3D &q3 ) const
{
    double x = q0[0]*(-q1[1]*q2[2] + q1[1]*q3[2] + q2[1]*q1[2] - 
		       q2[1]*q3[2] - q3[1]*q1[2] + q3[1]*q2[2] ) +
	       q1[0]*( q0[1]*q2[2] - q0[1]*q3[2] - q2[1]*q0[2] + 
		       q2[1]*q3[2] + q3[1]*q0[2] - q3[1]*q2[2] ) +
	       q2[0]*(-q0[1]*q1[2] + q0[1]*q3[2] + q1[1]*q0[2] - 
		       q1[1]*q3[2] - q3[1]*q0[2] + q3[1]*q1[2] ) +
	       q3[0]*( q0[1]*q1[2] - q0[1]*q2[2] - q1[1]*q0[2] + 
		       q1[1]*q2[2] + q2[1]*q0[2] - q2[1]*q1[2] );
    if( x < _signed_volume_eps ) {
	if( x <= -_signed_volume_eps )
	    return( -1 );
	return( 0 );
    }
    return( 1 );
}


int VTriangleSurfaceSolid::signvol3( const Vec3D &q1, const Vec3D &q2, const Vec3D &q3 ) const
{
    double x = q1[0]*( q2[1]*q3[2] - q3[1]*q2[2] ) -
	       q1[1]*( q2[0]*q3[2] - q3[0]*q2[2] ) +
	       q1[2]*( q2[0]*q3[1] - q3[0]*q2[1] );
    if( x < _signed_volume_eps ) {
	if( x <= -_signed_volume_eps )
	    return( -1 );
	return( 0 );
    }
    return( 1 );
}


/* Classify inclusion of point p in tetrahedron (o,q1,q2,q3), when
 * The sense of tetrahedron (o,q1,q2,q3) is ss.
 */
int VTriangleSurfaceSolid::classify_original_tetrahedron( int ss, const Vec3D &p, 
							  const Vec3D &q1, const Vec3D &q2, const Vec3D &q3 ) const
{
    int s0, s1, s2, s3;

    if( ss > 0 ) {
	// Positive sense (o,q1,q2,q3)
	if( (s1 = signvol3( p, q2, q3 )) < 0 )
	    return( VTRI_OUTSIDE );
	if( (s2 = signvol3( p, q3, q1 )) < 0 )
	    return( VTRI_OUTSIDE );
	if( (s3 = signvol3( p, q1, q2 )) < 0 )
	    return( VTRI_OUTSIDE );
	if( (s0 = signvol4( p, q1, q2, q3 )) < 0 )
	    return( VTRI_OUTSIDE );
    } else {
	// Negative sense (o,q1,q2,q3)
	if( (s1 = signvol3( p, q2, q3 )) > 0 )
	    return( VTRI_OUTSIDE );
	if( (s2 = signvol3( p, q3, q1 )) > 0 )
	    return( VTRI_OUTSIDE );
	if( (s3 = signvol3( p, q1, q2 )) > 0 )
	    return( VTRI_OUTSIDE );
	if( (s0 = signvol4( p, q1, q2, q3 )) > 0 )
	    return( VTRI_OUTSIDE );
    }

    // Edges and faces
    if( s1 == 0 ) {
	if( s2 == 0 )
	    return( VTRI_EDGE3 );
	else if( s3 == 0 )
	    return( VTRI_EDGE2 );
	return( VTRI_FACE );
    } else if( s2 == 0 ) {
	if( s3 == 0 )
	    return( VTRI_EDGE1 );
	return( VTRI_FACE );
    } else if( s3 == 0 ) {
	return( VTRI_FACE );
    }

    return( VTRI_INSIDE );
}


bool VTriangleSurfaceSolid::inside( const Vec3D &p ) const
{
    // Fast test if outside bbox
    for( uint32_t a = 0; a < 3; a++ ) {
	if( p[a] < _bbox[0][a] )
	    return( false );
	if( p[a] > _bbox[1][a] )
	    return( false );
    }
    
    // Offset point
    Vec3D x = p+_offset;

    // Cleared positive and negative vertex arrays
    std::vector<bool> vpos( _vertex.size(), false );
    std::vector<bool> vneg( _vertex.size(), false );

    int incl = 0;
    for( uint32_t a = 0; a < _triangle.size(); a++ ) {

	int ss = signvol3( _vertex[_triangle[a][0]],
			   _vertex[_triangle[a][1]],
			   _vertex[_triangle[a][2]] );
	int stat = classify_original_tetrahedron( ss, x,
						  _vertex[_triangle[a][0]],
						  _vertex[_triangle[a][1]],
						  _vertex[_triangle[a][2]] );
	if( stat == VTRI_INSIDE ) {
	    incl += 2*ss;
	} else if( stat == VTRI_FACE ) {
	    incl += ss;
	} else if( stat != VTRI_OUTSIDE ) {
	    if( ss > 0 && !vpos[_triangle[a][stat]] ) {
		vpos[_triangle[a][stat]] = true;
		incl += 2*ss;
	    } else if( ss < 0 && !vneg[_triangle[a][stat]] ) {
		vneg[_triangle[a][stat]] = true;
		incl += 2*ss;
	    }
	}
    }

    if( incl > 0 )
	return( true ); 
    return( false );
}


void VTriangleSurfaceSolid::update_bbox( Vec3D &min, Vec3D &max, const Vec3D x ) const
{
    for( int a = 0; a < 3; a++ ) {
	if( x[a] < min[a] )
	    min[a] = x[a];
	if( x[a] > max[a] )
	    max[a] = x[a];
    }
}


void VTriangleSurfaceSolid::prepare_for_inside()
{
    // Calculate bbox
    Vec3D min = Vec3D( std::numeric_limits<double>::infinity(),
		       std::numeric_limits<double>::infinity(),
		       std::numeric_limits<double>::infinity() );
    Vec3D max = Vec3D( -std::numeric_limits<double>::infinity(),
		       -std::numeric_limits<double>::infinity(),
		       -std::numeric_limits<double>::infinity() );

    for( uint32_t a = 0; a < _vertex.size(); a++ )
	update_bbox( min, max, _vertex[a] );

    _bbox[0] = min;
    _bbox[1] = max;

    // Set offset
    _offset = -_bbox[0] + 1.0e-3*(_bbox[1]-_bbox[0]);

    // Apply offset to vertex data
    for( size_t a = 0; a < _vertex.size(); a++ )
	_vertex[a] += _offset;
}


void VTriangleSurfaceSolid::set_signed_volume_eps( double signed_volume_eps )
{
    _signed_volume_eps = signed_volume_eps;
}


void VTriangleSurfaceSolid::get_bbox( Vec3D &min, Vec3D &max ) const
{
    min = _bbox[0];
    max = _bbox[1];
}


void VTriangleSurfaceSolid::remove_duplicate_triangles( void )
{
    int removed = 0;

    /* Check all pairs
     */
    for( uint32_t a = 0; a < _triangle.size(); a++ ) {
	const VTriangle &tria = _triangle[a];
	int a1 = tria[0];
	int a2 = tria[1];
	int a3 = tria[2];
	for( uint32_t b = a+1; b < _triangle.size(); b++ ) {
	    const VTriangle &trib = _triangle[b];
	    int b1 = trib[0];
	    int b2 = trib[1];
	    int b3 = trib[2];
	    if( a1 == b1 && a2 == b2 && a3 == b3 ) {
		// Duplicate found, remove triangle b
		removed++;
		_triangle.erase(_triangle.begin()+b);
	    }
	}
    }

    if( removed ) {
	ibsimu.message(1) << "Removed " << removed << " duplicate triangles\n";
    }
}


void VTriangleSurfaceSolid::check_data( void ) const
{
    for( uint32_t a = 0; a < _triangle.size(); a++ ) {

	const VTriangle &tri = _triangle[a];
	Vec3D V1 = _vertex[tri[1]] - _vertex[tri[0]];
	Vec3D V2 = _vertex[tri[2]] - _vertex[tri[0]];
	double area = 0.5*norm2( cross( V1, V2 ) );
	if( fabs(area) == 0.0 )
	    throw( Error( ERROR_LOCATION, "Zero area triangle " + to_string(a) ) );
	for( uint32_t b = 0; b < 3; b++ ) {
	    
	    int e1 = tri[b];
	    int e2 = tri[(b+1)%3];

	    bool found = false;
	    uint32_t neighbour;
	    uint32_t c;
	    for( c = 0; c < _triangle.size(); c++ ) {

		if( c == a ) continue; // No self-checking
		const VTriangle &tric = _triangle[c];
		for( uint32_t d = 0; d < 3; d++ ) {
		    
		    int f1 = tric[d];
		    int f2 = tric[(d+1)%3];
		    if( f1 == e1 && f2 == e2 ) {
			// Incorrect orientation
			throw( Error( ERROR_LOCATION, "Incorrect orientation between neighbouring triangles " +
				      to_string(a) + " and " + to_string(c) ) );
		    } else if( f1 == e2 && f2 == e1 ) {
			// Correct orientation
			if( found == true ) {
			    const VTriangle &trin = _triangle[neighbour];
			    throw( Error( ERROR_LOCATION, "Double neighbours (" + to_string(neighbour) +
					  " and " + to_string(c) + ") found for triangle " +
					  to_string(a) + ".\n" + 
					  "Triangle " + to_string(a) + ": " + to_string(tri[0]) + ", " + to_string(tri[1]) + ", " + to_string(tri[2]) + ":\n" +
					  "  " + to_string(_vertex[tri[0]]-_offset) + "\n" +
					  "  " + to_string(_vertex[tri[1]]-_offset) + "\n" +
					  "  " + to_string(_vertex[tri[2]]-_offset) + "\n" +
					  "Triangle " + to_string(neighbour) + ": " + to_string(trin[0]) + ", " + to_string(trin[1]) + ", " + to_string(trin[2]) + ":\n" +
					  "  " + to_string(_vertex[trin[0]]-_offset) + "\n" +
					  "  " + to_string(_vertex[trin[1]]-_offset) + "\n" +
					  "  " + to_string(_vertex[trin[2]]-_offset) + "\n" +
					  "Triangle " + to_string(c) + ": " + to_string(tric[0]) + ", " + to_string(tric[1]) + ", " + to_string(tric[2]) + ":\n" +
					  "  " + to_string(_vertex[tric[0]]-_offset) + "\n" +
					  "  " + to_string(_vertex[tric[1]]-_offset) + "\n" +
					  "  " + to_string(_vertex[tric[2]]-_offset) ) );
			}
			found = true;
			neighbour = c;
		    }
		}
	    }
	    if( !found ) {
		throw( Error( ERROR_LOCATION, "Triangle " + to_string(a) + " neighbour not found."
			      "Triangle " + to_string(a) + ": " + to_string(tri[0]) + ", " + to_string(tri[1]) + ", " + to_string(tri[2]) + ":\n" +
			      "  " + to_string(_vertex[tri[0]]-_offset) + "\n" +
			      "  " + to_string(_vertex[tri[1]]-_offset) + "\n" +
			      "  " + to_string(_vertex[tri[2]]-_offset) + "\n" ) );
	    }
	}
    }
}


void VTriangleSurfaceSolid::clear( void )
{
    VTriangleSurface::clear();
    //_vpos.clear();
    //_vneg.clear();
}


void VTriangleSurfaceSolid::debug_print( std::ostream &os ) const
{
    VTriangleSurface::debug_print( os );

    os << "**VTriangleSurfaceSolid\n";

    os << "signed_volume_eps = " << _signed_volume_eps << "\n";

    os << "offset = " << _offset << "\n";
    os << "bbox[0] = " << _bbox[0] << "\n";
    os << "bbox[1] = " << _bbox[1] << "\n";
}

