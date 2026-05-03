/*! \file vtriangle.hpp
 *  \brief Vertex-based triangle representation
 */

/* Copyright (c) 2011-2013,2017 Taneli Kalvas. All rights reserved.
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


#ifndef VTRIANGLE_HPP
#define VTRIANGLE_HPP 1


#include <iostream>
#include <vector>
#include <vec3d.hpp>
#include <stdint.h>


/*! \brief Vertex-based triangle representation.
 */
class VTriangle {
    
    uint32_t _v[3];
    
public:
    
    VTriangle( uint32_t v1, uint32_t v2, uint32_t v3 );
    VTriangle( const uint32_t v[3] );
    ~VTriangle();
    
    const uint32_t &operator[]( uint32_t i ) const {
	return( _v[i] );
    }

    uint32_t &operator[]( uint32_t i ) {
	return( _v[i] );
    }
    
    void debug_print( std::ostream &os ) const;
};


/*! \brief VTriangle surface.
 *
 *  Surface mesh constructed of vertex triangles and vertices.
 *
 *  %VTriangleSurface can be constructed triangle-by-triangle with
 *  internal elimination of duplicated vertices. Construction of of
 *  surface with external (more intelligent) vertex handling is also
 *  possible.
 */
class VTriangleSurface {

    double                     _vertex_matching_eps; /*!< \brief Tolerance for vertex matching. */

protected:

    std::vector<Vec3D>         _vertex;    /*!< \brief List of vertices for surface triangles. */
    std::vector<VTriangle>     _triangle;  /*!< \brief List of surface triangles. */

public:

    /*! \brief Constructor for vertex triangle surface.
     *
     *  The vertex matching tolerance can be set with \a
     *  vertex_matching_eps (defaults to 1.0e-9).
     */
    VTriangleSurface( double vertex_matching_eps = 1.0e-9 );

    /*! \brief Destructor.
     */
    ~VTriangleSurface();

    /*! \brief Set vertex matching tolerance.
     *
     *  Defaults to 1.0e-9.
     */ 
   void set_vertex_matching_eps( double vertex_matching_eps );

    /*! \brief Return vertex count.
     */
    uint32_t vertexc( void ) const {
	return( _vertex.size() );
    }

    /*! \brief Return vertex \a i coordinates.
     */
    const Vec3D &vertex( uint32_t i ) const {
	return( _vertex[i] );
    }

    /*! \brief Return triangle count.
     */
    uint32_t trianglec( void ) const {
	return( _triangle.size() );
    }

    /*! \brief Return triangle \a i.
     */
    const VTriangle &triangle( uint32_t i ) const {
	return( _triangle[i] );
    }

    /*! \brief Add a vertex \a x without checking.
     *
     *  Returns the index of the new vertex.
     */
    uint32_t add_vertex_no_check( const Vec3D &x ) {
	_vertex.push_back( x );
	return( _vertex.size()-1 );
    }

    /*! \brief Add a vertex \a x with duplicate vertex elimination.
     *
     *  Returns the index of the (new or old) vertex.
     */
    uint32_t add_vertex( const Vec3D &x );

    /*! \brief Add a triangle consiting of vertices \a x1, \a x2 and
        \a x3 with duplicate vertex elimination.
     *
     *  The vertices are assumed to be defined with right-hand
     *  ordering for defining outside direction. Each vertex of the
     *  triangle is checked against all vertices in the surface for
     *  duplicate entries. A vertex pairs closer than \a
     *  vertex_matching_eps are considered to be equal.
     *
     *  Returns the index of the new triangle.
     */
    uint32_t add_triangle( const Vec3D &x1, const Vec3D &x2, const Vec3D &x3 );

    /*! \brief Add a triangle consiting of vertices \a x with
     *  duplicate vertex elimination.
     *
     *  Convenience function.
     */
    uint32_t add_triangle( const Vec3D x[3] );

    /*! \brief Add a triangle consiting of already defined vertices \a
     *  v1, \a v2 and \a v3.
     */
    uint32_t add_triangle( uint32_t v1, uint32_t v2, uint32_t v3 ) {
	_triangle.push_back( VTriangle( v1, v2, v3 ) );
	return( _triangle.size()-1 );
    }

    /*! \brief Add a triangle consiting of already defined vertices \a
     *  v.
     *
     *  Convenience function.
     */
    uint32_t add_triangle( const uint32_t v[3] ) {
	_triangle.push_back( VTriangle( v ) );
	return( _triangle.size()-1 );
    }

    /*! \brief Clear surface.
     */
    void clear( void );

    /*! \brief Debug print.
     */
    void debug_print( std::ostream &os ) const;
};


/*! \brief VTriangleSolid solid.
 *
 *  Closed surface mesh constructed of vertex triangles and vertices.
 *  The test inside() is provided for testing for point inclusion in
 *  the solid. The prepare_for_inside() function must be called after
 *  constructing surface, before using inside().
 */
class VTriangleSurfaceSolid : public VTriangleSurface {

    double                     _signed_volume_eps; /*!< \brief Tolerance for inside test. */

    Vec3D                      _offset;    /*!< \brief Offset applied on vertices. */
    Vec3D                      _bbox[2];   /*!< \brief Bounding box min, max. */

    int signvol4( const Vec3D &q0, const Vec3D &q1, 
		  const Vec3D &q2, const Vec3D &q3 ) const;
    int signvol3( const Vec3D &q1, const Vec3D &q2, const Vec3D &q3 ) const;
    int classify_original_tetrahedron( int ss, const Vec3D &p, 
				       const Vec3D &q1, const Vec3D &q2, const Vec3D &q3 ) const;

    void update_bbox( Vec3D &min, Vec3D &max, const Vec3D x ) const;

public:

    /*! \brief Constructor for vertex triangle surface solid.
     *
     *  The vertex matching tolerance can be set with \a
     *  vertex_matching_eps (defaults to 1.0e-9) and signed volume
     *  tolerance used in inside test with \a signed_volume_eps
     *  (defaults to 1.0e-15).
     */
    VTriangleSurfaceSolid( double vertex_matching_eps = 1.0e-9, 
			   double signed_volume_eps = 1.0e-15 );

    /*! \brief Destructor.
     */
    ~VTriangleSurfaceSolid();

    /*! \brief Set signed volume tolerance.
     *
     *  Defaults to 1.0e-15.
     */ 
    void set_signed_volume_eps( double signed_volume_eps );

    /*! \brief Return vertex \a i coordinates.
     */
    Vec3D vertex( uint32_t i ) const {
	return( VTriangleSurface::vertex(i) - _offset );
    }

    /*! \brief Prepare for inside tests.
     *
     *  The inside() test algorithm is only capable of working in the
     *  positive octant. Therefore an offset is applied to the data
     *  and a boundary box test is performed to avoid forbidden tests.
     *
     *  The bounding box and offset are computed with this function
     *  after constructing the surface data. Call exactly once after
     *  defining solid.
     */
    void prepare_for_inside();

    /*! \brief Return if point \a x is inside solid.
     *
     *  Uses algorithm from R. Segura, F. R. Feito, J. Ruiz de Miras,
     *  C. Ogayar and J. C. Torres, "An Efficient Point Classification
     *  Algorithm for Triangle Meshes", Journal of Graphics, GPU, and
     *  Game Tools 10, issue 3, 2005.
     *
     *  The prepare_for_inside() function must be called before using
     *  this function.
     */
    bool inside( const Vec3D &x ) const;

    /*! \brief Return bounding box in vectors \a min and \a max.
     *
     *  The prepare_for_inside() function must be called before using
     *  this function.
     */
    void get_bbox( Vec3D &min, Vec3D &max ) const;
    
    /*! \brief Remove duplicate triangles.
     */
    void remove_duplicate_triangles( void );

    /*! \brief Check data.
     *
     *  1. Triangle edge pairing. Every triangle must have exactly 
     *     one neighbouring triangle for each of the three edges.
     *  2. Neighbouring triangles must be defined in same direction.
     *  3. No zero area triangles are allowed.
     *
     *  Throws an error if not valid.
     */
    void check_data( void ) const;

    /*! \brief Clear surface.
     */
    void clear( void );

    /*! \brief Debug print.
     */
    void debug_print( std::ostream &os ) const;
};


#endif
