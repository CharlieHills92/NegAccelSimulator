/*! \file mydxfblocks.hpp
 *  \brief DXF blocks
 */

/* Copyright (c) 2010,2012 Taneli Kalvas. All rights reserved.
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

#ifndef MY_DXF_BLOCKS_HPP
#define MY_DXF_BLOCKS_HPP 1


#include <stdint.h>
#include <vector>
#include <string>
#include "vec3d.hpp"
#include "mydxffile.hpp"
#include "mydxfentities.hpp"
#include "transformation.hpp"


/*! \brief DXF block class.
 */
class MyDXFBlock
{

    std::string    _block_handle;
    std::string    _block_layer;

    std::string    _endblk_handle;
    std::string    _endblk_layer;

    std::string    _path;
    std::string    _owner_handle;
    std::string    _name;          // Name reference used by INSERT

    int16_t        _type;
    Vec3D          _p;             // Base point

    class MyDXFEntities *_entities;

public:

    MyDXFBlock( class MyDXFFile *dxf );
    ~MyDXFBlock();

    /*! \brief Write dxf file to stream.
     */
    void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Return name of block.
     */
    const std::string &name( void ) { return( _name ); }

    /*! \brief Get a pointer to the entities of block.
     */
    class MyDXFEntities *get_entities( void ) { return( _entities ); }

    /*! \brief Get a const pointer to the entities of block.
     */
    const class MyDXFEntities *get_entities( void ) const { return( _entities ); }

    /*! \brief Explode block into entities.
     *
     *  Tranform entities in the block with tranformation \a t and add
     *  the tranformed entities to the database \a ent.
     */
    void explode( class MyDXFEntities *ent, MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Plot block with cairo
     *
     *  Plot the entities withing the block using the transformation
     *  \a t from the object space to cairo coordinates. The
     *  visible range is specified by \a range (xmin,ymin,xmax,ymax)
     *  in cairo coordinates.
     */
    void plot( const class MyDXFFile *dxf, cairo_t *cairo, 
	       const Transformation *t, const double range[4] ) const;

    /*! \brief Return bounding box of entities within the block.
     */
    void get_bbox( Vec3D &min, Vec3D &max, 
		   const class MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Scale entities within block by factor \a s.
     */
    void scale( class MyDXFFile *dxf, double s );
    
    /*! \brief Translate entity by \a dx.
     */
    void translate( class MyDXFFile *dxf, const Vec3D &dx );

    /*! \brief Print debugging information to stream \a os.
     */
    void debug_print( std::ostream &os ) const;

    /*! \brief Print debugging information to stream \a os.
     */
    friend std::ostream &operator<<( std::ostream &os, const MyDXFBlock &blk );
};



/*! \brief DXF blocks class.
 *
 *  Container for data of a DXF file blocks.
 */
class MyDXFBlocks
{

    std::vector<MyDXFBlock *> _blocks;

public:

    MyDXFBlocks( class MyDXFFile *dxf );
    ~MyDXFBlocks();

    /*! \brief Write dxf file to stream.
     */
    void write( class MyDXFFile *dxf, std::ofstream &ostr );

    uint32_t size( void ) const { return( _blocks.size() ); }

    MyDXFBlock *get_by_name( const std::string &name );
    const MyDXFBlock *get_by_name( const std::string &name ) const;

    MyDXFBlock *operator()( int a ) { return( _blocks[a] ); }

    const MyDXFBlock *operator()( int a ) const { return( _blocks[a] ); }

    void clear( void );

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;
};





#endif




