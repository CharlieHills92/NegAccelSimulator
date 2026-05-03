/*! \file mydxfentities.hpp
 *  \brief DXF Entities
 */

/* Copyright (c) 2010-2011,2014 Taneli Kalvas. All rights reserved.
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

#ifndef MY_DXF_ENTITIES_HPP
#define MY_DXF_ENTITIES_HPP 1


#include <stdint.h>
#include <vector>
#include <cairo.h>
#include "mydxffile.hpp"
#include "vec3d.hpp"
#include "transformation.hpp"


#define MYDXF_PERT_EPS 1.1e-12


/*! \brief Entity type.
 */
enum EntityType {
    ENTITY_UNKNOWN = 0,
    ENTITY_LINE,
    ENTITY_LWPOLYLINE,
    ENTITY_SPLINE,
    ENTITY_ARC,
    ENTITY_CIRCLE,
    ENTITY_MTEXT,
    ENTITY_INSERT
};


/*! \brief DXF entity base class.
 *
 *  A general base class for all DXF entities. Contains data fields
 *  common to all entities.
 */
class MyDXFEntity
{

protected: 

    std::string _handle;
    std::string _layer;

    MyDXFEntity();

    //MyDXFEntity( const MyDXFEntity &ent );

    /*! \brief Propose a point to bounding box.
     *
     *  Updates bounding box value at min and max by including point p
     *  in the bounding box.
     */
    static void bbox_ppoint( Vec3D &min, Vec3D &max, const Vec3D &p );

    void write_common( class MyDXFFile *dxf, std::ofstream &ostr );
    void process_group( class MyDXFFile *dxf );
    void constructor_debug_print( void ) const;
    void debug_print_base( std::ostream &os ) const;

public:

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFEntity() {}

    /*! \brief Get a new copy of entity.
     */
    virtual MyDXFEntity *copy( void ) const = 0;

    /*! \brief Explode into entities.
     *
     *  Break entity into atomic entities and tranform the entities
     *  with tranformation \a t. Add the tranformed entities to the
     *  database \a ent.
     */
    virtual void explode( class MyDXFEntities *ent, MyDXFFile *dxf, const Transformation *t ) const = 0;

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr ) = 0;

    /*! \brief Scale entity by factor \a s.
     */
    virtual void scale( class MyDXFFile *dxf, double s ) = 0;

    /*! \brief Translate entity by \a dx.
     */
    virtual void translate( class MyDXFFile *dxf, const Vec3D &dx ) = 0;

    /*! \brief Rotate entity around origin
     *
     *  Rotate for \a a radians.
     */
    virtual void rotate_z( class MyDXFFile *dxf, double a ) = 0;

    /*! \brief Set layer.
     */
    void set_layer( const std::string &layer ) { _layer = layer; }

    /*! \brief Get layer.
     */
    std::string get_layer( void ) const { return( _layer ); }

    /*! \brief Get entity type.
     */
    virtual EntityType get_type( void ) const = 0;

    /*! \brief Set entity handle.
     */
    void set_handle( const std::string &handle ) { _handle = handle; }

    /*! \brief Get entity handle.
     */
    std::string get_handle( void ) const { return( _handle ); }

    /*! \brief Plot entity with cairo
     *
     *  Plot the entity using the transformation \a t from the
     *  object space to cairo coordinates. The visible range is
     *  specified by \a range (xmin,ymin,xmax,ymax) in cairo
     *  coordinates.
     */
    virtual void plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		       const Transformation *t, const double range[4] ) const = 0;

    /*! \brief Return bounding box of entity
     */
    virtual void get_bbox( Vec3D &min, Vec3D &max, 
			   const class MyDXFFile *dxf, const Transformation *t ) const = 0;

    /*! \brief Print debugging information to stream \a os.
     */
    virtual void debug_print( std::ostream &os ) const = 0;

    friend std::ostream &operator<<( std::ostream &os, const MyDXFEntity &ent );
};


/*! \brief DXF path entity base class.
 *
 *  A base class for two dimensional DXF entities, which can be part
 *  of a path. All path entities have a start point and an end point,
 *  that can be read and set.
 */
class MyDXFPathEntity : public MyDXFEntity
{

protected: 

    MyDXFPathEntity() {}

    MyDXFPathEntity( const MyDXFEntity &ent ) : MyDXFEntity(ent) {}

public:

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFPathEntity() {}

    /*! \brief Get start point of path entity.
     */
    virtual Vec3D start( void ) const = 0;

    /*! \brief Get end point of path entity.
     */
    virtual Vec3D end( void ) const = 0;

    /*! \brief Set start point of path entity.
     */
    virtual void set_start( const Vec3D &s ) = 0;

    /*! \brief Set end point of path entity.
     */
    virtual void set_end( const Vec3D &e ) = 0;

    /*! \brief Check for ray crossing.
     *
     *  Check if ray going from point (x,y) downwards (negative y
     *  direction) crosses the entity. Return 1 if crosses odd number
     *  of times and 0 if even number of times. Return 2 in case of
     *  exact crossing at boundaries. This function is used as a
     *  subroutine to inside_loop().
     */
    virtual int ray_cross( double x, double y ) const = 0;
};




/*! \brief DXF entity selection.
 *
 *  %MyDXFEntitySelection object is a list of indexes, which point to
 *  entities in the MyDXFEntities.
 */
class MyDXFEntitySelection
{

    std::vector<uint32_t> _selection;

public:

    /*! \brief Construct empty selection.
     */
    MyDXFEntitySelection() {}

    /*! \brief Destructor.
     */
    ~MyDXFEntitySelection() {}

    /*! \brief Return number of entities in selection.
     */
    uint32_t size() const { return( _selection.size() ); }

    /*! \brief Add entity number \a a in selection.
     */
    void add_entity( uint32_t a ) { _selection.push_back( a ); }

    /*! \brief Get a const reference to entity number in selection at location \a a.
     */
    const uint32_t &operator()( int a ) const { 
	if( a < 0 || a >= (int)_selection.size() )
	    throw( Error( ERROR_LOCATION, "index out of range" ) );
	return( _selection[a] ); 
    }

    /*! \brief Get reference to entity number in selection at location \a a.
     */
    uint32_t &operator()( int a ) { 
	if( a < 0 || a >= (int)_selection.size() )
	    throw( Error( ERROR_LOCATION, "index out of range" ) );
	return( _selection[a] ); 
    }

    friend std::ostream &operator<<( std::ostream &os, const MyDXFEntitySelection &sel );
};


/*! \brief DXF entity database.
 *
 *  A database of entities. Also responsible for reading entities from
 *  a DXF file. All supported entities are saved to database. All
 *  others are ignored.
 */
class MyDXFEntities
{

    MyDXFFile                  *_dxf;      /*!< \brief Pointer to parent dxffile */
    std::vector<MyDXFEntity *>  _entities; /*!< \brief Entities */

public:


    /*! \brief Construct empty entities database.
     */
    MyDXFEntities( class MyDXFFile *dxf );

    /*! \brief Construct new entities containing copies of selected intities in \a ent.
     */
    MyDXFEntities( class MyDXFFile *dxf, MyDXFEntities *ent, MyDXFEntitySelection *sel );

    /*! \brief Construct entities database by reading from DXF file.
     *
     *  Called with reading_block = true if called from inside BLOCKS
     *  section and reading_block = false if called inside ENTITIES section.
     */
    MyDXFEntities( class MyDXFFile *dxf, bool reading_blocks );

    /*! \brief Destructor.
     */
    ~MyDXFEntities();


    /*! \brief Write entities section of dxf file to stream.
     */
    void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Write a list of entities to stream.
     *
     *  This function only writes the entities withing the object. It
     *  can be called for writing the entities section or a block.
     */
    void write_entities( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Return number of entities.
     */
    uint32_t size() const { return( _entities.size() ); }

    /*! \brief Return const pointer to entity \a a.
     */
    const MyDXFEntity *get_entity( uint32_t a ) const { return( _entities[a] ); }

    /*! \brief Return pointer to entity \a a.
     */
    MyDXFEntity *get_entity( uint32_t a ) { return( _entities[a] ); }



    
    /*! \brief Add entity to list.
     *
     *  No copy of entity is made. The pointer is saved to the database.
     */
    void add_entity( MyDXFEntity *e ) { _entities.push_back( e ); }


    /*! \brief Make a new selection containg all entities from database.
     */
    MyDXFEntitySelection *selection_all( void ) const;

    /*! \brief Make a new selection containg entities from named layer.
     */
    MyDXFEntitySelection *selection_layer( const std::string &layername ) const;

    /*! \brief Make a new selection containg entities of given type.
     */
    MyDXFEntitySelection *selection_type( EntityType type ) const;

    /*! \brief Build complete loops.
     *
     * Make a new subselection of a selection. The new selection will
     * contain only objects, which make up one or several complete
     * loops. Ending point of one entity is ensured to be exactly the
     * starting point of another entity. Errors of the size eps are
     * accepted and fixed between end-to-end matching.
     */
    MyDXFEntitySelection *selection_path_loop( MyDXFEntitySelection *selection,
					       double eps = 1.0e-6 );





    /*! \brief Check if two entities are geometrically same.
     *
     *   Checks if entity \a a is the geometrically same as entity \a
     *   b within error limit \a eps.
     */
    bool geom_same( uint32_t a, uint32_t b, double eps = 1.0e-6 ) const;

    /* ! \brief Check if point is inside a loop defined by a selection
     *   of entities.
     *
     *   The check is done assuming a 2D drawing in xy-plane. The
     *   check is done using ray shooting algorithm. If exact crossing
     *   happens perturbation algorithm is used. New test is performed
     *   at eps distance from the first.
     */
    bool inside_loop( MyDXFEntitySelection *selection, double x, double y, double eps = 1.0e-6 );

    /*! \brief Plot selected entities with cairo
     *
     *  Plot the entities using the transformation \a t from the
     *  object space to cairo coordinates. The visible range is
     *  specified by \a range (xmin,ymin,xmax,ymax) in cairo
     *  coordinates.
     *
     *  Selection can be a NULL pointer to plot all entities.
     */
    void plot( const MyDXFEntitySelection *selection, const class MyDXFFile *dxf, 
	       cairo_t *cairo, const Transformation *t, const double range[4] ) const;

    /*! \brief Get bounding box containing all entities in selection.
     */
    void get_bbox( const MyDXFEntitySelection *selection, Vec3D &min, Vec3D &max, 
		   const class MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Scale selected entities by factor s.
     *
     *  Selection can be a NULL pointer to scale all entities.
     */
    void scale( MyDXFEntitySelection *selection, class MyDXFFile *dxf, double s );

    /*! \brief Translate selected entities by \a dx.
     */
    void translate( MyDXFEntitySelection *selection, class MyDXFFile *dxf, const Vec3D &dx );

    /*! \brief Rotate selected entities around origin.
     *
     *  Rotate for \a a radians.
     */
    void rotate_z( MyDXFEntitySelection *selection, double a );

    /*! \brief Remove selected entities.
     *
     *  Selection can be a NULL pointer to remove all entities. The
     *  selection is invalid after this operation and should not be
     *  used further. Also all other selections are invalidated by
     *  this operation because entity indices change.
     */
    void remove( MyDXFEntitySelection *selection );


    /*! \brief Explode selected insert entities.
     *
     *  The insert entities are expoded to contain just primitive
     *  entities with no dependencies to blocks. Selection can be a
     *  NULL pointer to explode all entities.
     */
    void explode( MyDXFEntitySelection *selection, class MyDXFFile *dxf );

    /*! \brief Explode all entities to \a ent.
     *
     *  Explode and add all entities into entities database \a ent
     *  using transformation \a t.
     */
    void explode( MyDXFEntities *ent, class MyDXFFile *dxf, const Transformation *t ) const;


    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;


};





#endif




