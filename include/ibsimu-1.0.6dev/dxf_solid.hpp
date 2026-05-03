/*! \file dxf_solid.hpp
 *  \brief %Solid definition using MyDXF
 */

/* Copyright (c) 2010-2012,2014 Taneli Kalvas. All rights reserved.
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

#ifndef DXFSOLID_HPP
#define DXFSOLID_HPP 1


#include <iostream>
#include "solid.hpp"


class MyDXFFile;
class MyDXFEntities;
class MyDXFEntitySelection;


/*! \brief %MyDXFFile solid class.
 *
 *  %DXFSolid is an implementation of %Solid using MyDXFFile
 *  entities. The solid is built from a two dimensional area defined
 *  by enclosing the area with dxf path objects in one layer. The
 *  solid volume in (three dimensional) simulation space is defined
 *  using a combination of two transformations. The first
 *  transformation is from simulation space to intermediate 3D space
 *  and it is made using the Transformation class. These intermediate
 *  3D space points are then mapped to two dimensional dxf space using
 *  an optional user defined function. The Transformation defaults to
 *  unity matrix and if the user defined function is left undefined it
 *  defaults to \f$ (x,y,z) \Rightarrow (x,y) \f$.
 */
class DXFSolid : public Solid {

    Vec3D                (*_func)(const Vec3D &);
    MyDXFEntities         *_entities;
    MyDXFEntitySelection  *_selection;
    
public:

    /*! \brief Constructor for making a solid from a DXF-file layer.
     *
     *  The entities from the DXF-file layer \a layername are copied
     *  to DXFSolid object. No dependency stays between dxffile and
     *  the object constructed. The transformations are initialized to
     *  unity.
     */
    DXFSolid( MyDXFFile *dxffile, const std::string &layername );

    /*! \brief Constructor for making a solid from entities in a
     *  MyDXFEntities object \a ent.
     *
     *  Entities from \a ent are copied to DXFSolid object. The
     *  entities are from \a dxffile. An internal copy of the entities
     *  object is made. The transformations are initialized to unity.
     */
    DXFSolid( MyDXFFile *dxffile, MyDXFEntities *ent );

    /*! \brief Constructor for loading solid data from stream \a is.
     */
    DXFSolid( std::istream &is );

    /*! \brief Destructor.
     */
    virtual ~DXFSolid();

    /*! \brief Return if 3D point \a x in simulation space is inside
     *  solid.
     */
    virtual bool inside( const Vec3D &x ) const;

    /*! \brief Print debugging information to stream \a os.
     */
    void debug_print( std::ostream &os ) const;

    /*! \brief Unity transformation.
     *
     *  Default tranformation: \f$ (x,y,z) \Rightarrow (x,y) \f$.
     */
    static Vec3D unity( const Vec3D &x );

    /*! \brief %Solid of revolution around x-axis.
     *
     *  Tranformation: \f$ (x,y,z) \Rightarrow (x,\sqrt(y^2+z^2)) \f$.
     */
    static Vec3D rotx( const Vec3D &x );

    /*! \brief %Solid of revolution around y-axis.
     *
     *  Tranformation: \f$ (x,y,z) \Rightarrow (y,\sqrt(x^2+z^2)) \f$.
     */
    static Vec3D roty( const Vec3D &x );

    /*! \brief %Solid of revolution around z-axis.
     *
     *  Tranformation: \f$ (x,y,z) \Rightarrow (z,\sqrt(x^2+y^2)) \f$.
     */
    static Vec3D rotz( const Vec3D &x );

    /*! \brief Define mapping from 3D space to 2D space.
     *
     *  The mapping function can be user defined or one of the
     *  predefined functions: unity(), rotx(), roty() or rotz(). The
     *  mapping function can return a vector with NaN components for
     *  guaranteed inside solid result. Similarly infinity is
     *  guaranteed to give free space result.
     */
    void define_2x3_mapping( Vec3D (*func)(const Vec3D &) );

    /*! \brief Saves solid data to stream \a os.
     */
    virtual void save( std::ostream &os ) const;
};


#endif

