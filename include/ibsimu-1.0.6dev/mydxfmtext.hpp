/*! \file mydxfmtext.hpp
 *  \brief DXF mtext entity
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


#ifndef MY_DXF_MTEXT_HPP
#define MY_DXF_MTEXT_HPP 1


#include "mydxfentities.hpp"


#define ATTACHMENT_POINT_TOP_LEFT       1
#define ATTACHMENT_POINT_TOP_CENTER     2
#define ATTACHMENT_POINT_TOP_RIGHT      3
#define ATTACHMENT_POINT_MIDDLE_LEFT    4
#define ATTACHMENT_POINT_MIDDLE_CENTER  5
#define ATTACHMENT_POINT_MIDDLE_RIGHT   6
#define ATTACHMENT_POINT_BOTTOM_LEFT    7
#define ATTACHMENT_POINT_BOTTOM_CENTER  8
#define ATTACHMENT_POINT_BOTTOM_RIGHT   9

#define DRAWING_DIRECTION_LEFT_TO_RIGHT 1
#define DRAWING_DIRECTION_TOP_TO_BOTTOM 3
#define DRAWING_DIRECTION_BY_STYLE      5


/*! \brief DXF text entity class.
 *
 *
 */
class MyDXFMText : public MyDXFEntity
{

    std::string _text;               /* 1, 3, Text string. If the text string is less than 250 characters, all 
				      * characters appear in group 1. If the text string is greater than 250
				      * characters, the string is divided into 250-character chunks, which
				      * appear in one or more group 3 codes. If group 3 codes are used, the
				      * last group is a group 1 and has fewer than 250 characters. */
    Vec3D       _p;                  /* 10, 20, 30, Insertion point. */
    double      _text_height;        /* 40, Nominal (initial) text height. */
    double      _rect_width;         /* 41, Reference rectangle width. */
    double      _text_width;         /* 42, Horizontal width of the characters that make up the mtext entity. This 
				      * value will always be equal to or less than the value of group code 41 
				      * (read-only, ignored if supplied). */
    double      _vert_height;        /* 43, Vertical height of the mtext entity (read-only, ignored if supplied). */
    double      _spacing_fac;        /* 44, Mtext line spacing factor (optional):
				      * Percentage of default (3-on-5) line spacing to be applied. Valid values
				      * range from 0.25 to 4.00 */
    double      _rotation;           /* 50, Rotation angle in radians. */
    int16_t     _attachment_point;   /* 71, Attachment point: 
				      * 1 = Top left; 2 = Top center; 3 = Top right
				      * 4 = Middle left; 5 = Middle center; 6 = Middle right
				      * 7 = Bottom left; 8 = Bottom center; 9 = Bottom right */
    int16_t     _drawing_direction;  /* 72, Drawing direction:
				      * 1 = Left to right
				      * 3 = Top to bottom
				      *	5 = By style (the flow direction is inherited from the associated text style) */
    int16_t     _line_spacing;       /* 73, Mtext line spacing style (optional):
				      *	1 = At least (taller characters will override)
				      * 2 = Exact (taller characters will not override) */
    std::string _style;              /* 7, Text style name (STANDARD if not provided) (optional). */
    Vec3D       _extrusion;          /* 210, 220, 230, Extrusion direction (optional; default = 0, 0, 1) */
    Vec3D       _xaxis;              /* 11, 21, 31, X-axis direction vector (in WCS).A group code 50 (rotation angle
				      * in radians) passed as DXF input is converted to the equivalent direction 
				      * vector (if both a code 50 and codes 11, 21, 31 are passed, the last one wins).
				      * This is provided as a convenience for conversions from text objects */
    
    
public:

    /*! \brief Default constructor.
     */
    MyDXFMText();

    /*! \brief Construct line entity by reading from DXF file.
     */
    MyDXFMText( class MyDXFFile *dxf );

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFMText() {}

    /*! \brief Get a new copy of entity.
     */
    virtual MyDXFMText *copy( void ) const { return( new MyDXFMText( *this ) ); }

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
    virtual EntityType get_type( void ) const { return( ENTITY_MTEXT ); }

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

