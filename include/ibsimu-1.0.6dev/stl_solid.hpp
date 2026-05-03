/*! \file stl_solid.hpp
 *  \brief %Solid definition using Stereolithography CAD format
 */

/* Copyright (c) 2011,2012 Taneli Kalvas. All rights reserved.
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

#ifndef STLSOLID_HPP
#define STLSOLID_HPP 1


#include <iostream>
#include <vector>
#include "solid.hpp"
#include "transformation.hpp"


/*! \brief %STL solid
 *
 *  A solid object constructed from one or a union of several entities
 *  from STL-files.
 */
class STLSolid : public Solid {

    std::vector<class STLFile *> _stl;
    
public:

    /*! \brief Default constructor.
     *
     *  Create and empty object.
     */
    STLSolid();

    /*! \brief Constructor for making a solid reading a STL-file.
     */
    STLSolid( const std::string &filename );

    /*! \brief Constructor for loading solid data from stream \a is.
     */
    STLSolid( std::istream &is );

    /*! \brief Destructor.
     */
    virtual ~STLSolid();

    /*! \brief Return if 3D point \a x in simulation space is inside
     *  solid.
     */
    virtual bool inside( const Vec3D &x ) const;

    /*! \brief Add entity from STL-file to object.
     *
     *  The STLFile stl is owned by STLSolid after calling this
     *  function. It is destructed when STLSolid is destructed.
     */
    void add_stl_file( class STLFile *stl );

    /*! \brief Return a pointer to the STL-file \a i.
     */
    class STLFile *get_stl_file( uint32_t i = 0 ) const;

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;

    /*! \brief Saves solid data to stream.
     */
    virtual void save( std::ostream &s ) const;
};


#endif

