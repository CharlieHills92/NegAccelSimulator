/*! \file mydxfarc.chp
 *  \brief DXF arc entity
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


#ifndef MY_DXF_ARC_HPP
#define MY_DXF_ARC_HPP 1


#include "mydxfentities.hpp"


/*! \brief DXF arc entity class.
 *
 *  An arc entity is defined with a center point of the arc, radius,
 *  starting angle and ending angle.
 */
class MyDXFArc : public MyDXFPathEntity
{

    Vec3D    _pc;
    double   _r;
    double   _ang1; // Must be between 0 and 2 pi.
    double   _ang2; // Must be between 0 and 2 pi.

public:

    /*! \brief Default constructor.
     */
    MyDXFArc() : _r(1.0), _ang1(0.0), _ang2(2.0*M_PI) {};

    /*! \brief Construct arc entity by reading from DXF file.
     */
    MyDXFArc( class MyDXFFile *dxf );

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFArc() {}

    /*! \brief Get a new copy of entity.
     */
    virtual MyDXFArc *copy( void ) const { return( new MyDXFArc( *this ) ); }

    /*! \brief Explode into entities.
     *
     *  Break entity into atomic entities and tranform entities them
     *  with tranformation \a t. Add the tranformed entities to the
     *  database \a ent.
     */
    virtual void explode( class MyDXFEntities *ent, MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Get entity type.
     */
    virtual EntityType get_type( void ) const { return( ENTITY_ARC ); }

    /*! \brief Get center point of arc.
     */
    Vec3D center( void ) const { return( _pc ); }

    /*! \brief Get radius of arc.
     */
    double radius( void ) const { return( _r ); }

    /*! \brief Get start point of path entity.
     */
    virtual Vec3D start( void ) const {
	return( Vec3D(_pc[0] + _r*cos(_ang1), _pc[1] + _r*sin(_ang1), _pc[2] ) ); 
    }

    /*! \brief Get end point of path entity.
     */
    virtual Vec3D end( void ) const { 
	return( Vec3D(_pc[0] + _r*cos(_ang2), _pc[1] + _r*sin(_ang2), _pc[2] ) ); 
    }

    /*! \brief Set center point.
     */
    void set_pc( const Vec3D &pc ) { _pc = pc; }

    /*! \brief Set radius.
     */
    void set_r( double r ) { _r = r; }

    /*! \brief Set start angle.
     */
    void set_ang1( double ang1 );

    /*! \brief Set end angle.
     */
    void set_ang2( double ang2 );
    
    /*! \brief Get start angle.
     */
    double get_ang1( void ) const { return( _ang1 ); }

    /*! \brief Get end angle.
     */
    double get_ang2( void ) const { return( _ang2 ); }

    /*! \brief Set arc according to center and end points.
     *
     *  End point is only used for the angle.
     */
    void set_center_and_ends( const Vec3D &c, const Vec3D &s, const Vec3D &e );

    /*! \brief Set arc according to end points.
     *
     *  Resets the center point and angles according to starting point
     *  \a s and ending point \a e.  Arc is assumed to go in
     *  counter-clockwise direction from start to end. Radius of arc
     *  is kept constant if distance from \a s to \a e is less than
     *  double the radius. Otherwise the center point is set to the
     *  middle of start and end and the radius is set to half of the
     *  distance of start to end.
     */
    void set_center_point( const Vec3D &s, const Vec3D &e );

    /*! \brief Set start point of path entity.
     */
    virtual void set_start( const Vec3D &s );

    /*! \brief Set end point of path entity.
     */
    virtual void set_end( const Vec3D &e );

    /*! \brief Check for ray crossing.
     *
     *  Check if ray going from point (x,y) downwards (negative y
     *  direction) crosses the entity. Return 1 if crosses odd number
     *  of times and 0 if even number of times. Return 2 in case of
     *  exact crossing at boundaries. This function is used as a
     *  subroutine to inside_loop().
     */
    virtual int ray_cross( double x, double y ) const;

    /*! \brief Check if two entities are geometrically same.
     *
     *   Checks if entity \a a is the geometrically same as entity \a
     *   b within error limit \a eps.
     */
    bool geom_same( const MyDXFArc &arc, double eps = 1.0e-6 ) const;

    /*! \brief Plot entity with cairo
     *
     *  Plot the entity using the transformation \a from from the
     *  object space to cairo coordinates. The visible range is
     *  specified by \a range (xmin,ymin,xmax,ymax) in cairo
     *  coordinates.
     */
    virtual void plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		       const Transformation *t, const double range[4] ) const;

    /*! \brief Return bounding box of entity
     */
    virtual void get_bbox( Vec3D &min, Vec3D &max, 
			   const class MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Scale entity by factor \a s.
     */
    virtual void scale( class MyDXFFile *dxf, double s );

    /*! \brief Translate entity by \a dx.
     */
    virtual void translate( class MyDXFFile *dxf, const Vec3D &dx );

    /*! \brief Rotate entity around origin
     *
     *  Rotate for \a a radians.
     */
    virtual void rotate_z( class MyDXFFile *dxf, double a );

    /*! \brief Print debugging information to stream \a os.
     */
    virtual void debug_print( std::ostream &os ) const;
};


#endif

