/*! \file mesh.hpp
 *  \brief Rectangular mesh definition.
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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

#ifndef MESH_HPP
#define MESH_HPP 1


#include <stdint.h>
#include <vector>
#include <iostream>
#include "file.hpp"
#include "vec3d.hpp"
#include "types.hpp"


/*! \brief %Mesh geometry definion.
 *
 *  Class contains mesh geometry definition. It stores geometry mode (\a
 *  geom_mode), number of mesh nodes in each direction (\a size), the
 *  mesh cell size (\a h) and the locations of mesh node \a (0,0,0) and
 *  \a (size[0]-1,size[1]-1,size[2]-1) (known as \a origo and \a max).
 *  The \a max point is internally calculated. Other parameters are
 *  given when %Mesh is constructed.
 *
 *  %Mesh is to be used as a base class in all classes, which store or
 *  process some kind of mesh data.
 */
class Mesh
{
protected:
    geom_mode_e                _geom_mode; /*!< \brief %Geometry mode */
    Int3D                      _size;      /*!< \brief Size of mesh */
    Vec3D                      _origo;     /*!< \brief Location of mesh point (0,0,0) [m] */
    Vec3D                      _max;       /*!< \brief Location of mesh point 
					    * (size[0]-1,size[1]-1,size[2]-1) [m] */
    double                     _h;         /*!< \brief Length of mesh step [m] */
    double                     _div_h;     /*!< \brief Reciprocal of length of mesh step [1/m] */

public:

    /*! \brief Default constructor for mesh definition.
     *
     *  Sets geometry mode to MODE3D, mesh cell size \a h to 1, mesh
     *  size \a size to (0,0,0) and origo \a origo to (0,0,0).
     */
    Mesh();

    /*! \brief Constructor for mesh definition.
     *
     *  Sets geometry mode, mesh cell size \a h, mesh size \a size and
     *  origo \a origo.
     */
    Mesh( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h );

    /*! \brief Constructoer for loading mesh from a stream \a is.
     */
    Mesh( std::istream &is );

    /*! \brief Destructor.
     */
    ~Mesh() {}

    /*! \brief Reset mesh definition.
     */
    void reset( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h );

    /*! \brief Returns geometry mode.
     */
    geom_mode_e geom_mode( void ) const { return( _geom_mode ); }

    /*! \brief Returns number of dimensions in geometry.
     */
    uint32_t dim( void ) const;

    /*! \brief Returns size array of geometry.
     */
    Int3D size( void ) const { return( _size ); }

    /*! \brief Returns size of solid mesh in direction \a i.
     */
    uint32_t size( int i ) const { return( _size[i] ); }
   
    /*! \brief Returns number of nodes in the mesh.
     */
    uint32_t nodecount( void ) const { return( _size[0]*_size[1]*_size[2] ); }

    /*! \brief Returns origo vector of geometry.
     */
    Vec3D origo( void ) const { return( _origo ); }

    /*! \brief Returns \a i-th component of vector origo.
     */
    double origo( int i ) const { return( _origo[i] ); }

    /*! \brief Returns vector pointing to the last mesh point opposite
     *  of origo.
     */
    Vec3D max( void ) const { return( _max ); }

    /*! \brief Returns \a i-th component of vector pointing to the
     *  last mesh point opposite of origo.
     */
    double max( int i ) const { return( _max[i] ); }

    /*! \brief Returns mesh cell size.
     */
    double h( void ) const { return( _h ); }

    /*! \brief Returns reciprocal of mesh cell size (1/h).
     */
    double div_h( void ) const { return( _div_h ); }

    /*! \brief Returns node closest to location \a x.
     *
     *  Calculated as \code floor((x[i]-origo(i))/h() + 0.5) \endcode
     *  for each component
     */
    Int3D closest_node( Vec3D x ) const;

    /*! \brief Returns node, which contains \a x.
     *
     *  Calculated as \code floor((x[i]-origo(i))/h()) \endcode
     *  for each component
     */
    Int3D mesh_number( Vec3D x ) const;

    /*! \brief Returns coordinates of node \a i.
     *
     *  Calculated as \code origo(i)+n[i]*h() \endcode
     *  for each component
     */
    Vec3D coord_of_node( Int3D n ) const;

    /*! \brief Saves geometry data to stream \a os.
     */
    void save( std::ostream &os ) const;

    /*! \brief Equality.
     *
     *  Allows small inequality.
     */
    bool operator==( const Mesh &m ) const;

    /*! \brief Non-equality.
     */
    bool operator!=( const Mesh &m ) const;

    /*! \brief Print debugging information to stream \a os.
     */
    void debug_print( std::ostream &os ) const;
};


#endif

