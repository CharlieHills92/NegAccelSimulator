/*! \file label.hpp
 *  \brief Plot labels
 */

/* Copyright (c) 2005-2010,2012,2013 Taneli Kalvas. All rights reserved.
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

#ifndef LABEL_HPP
#define LABEL_HPP 1


#include <string>
#include <cairo.h>
#include "vec3d.hpp"


/*! \brief Class for labels in plots.
 *
 *  Label can be used to draw text labels with latex formatting.
 */
class Label
{
    std::string          _text;      /*!< \brief Label content. */

    double               _size;      /*!< \brief Size of font. */
    std::string          _family;    /*!< \brief Family of font. */
    cairo_font_slant_t   _slant;     /*!< \brief Slant of font.
				      *
				      * Either CAIRO_FONT_SLANT_NORMAL, 
				      * CAIRO_FONT_SLANT_ITALIC or CAIRO_FONT_SLANT_OBLIQUE. */
    cairo_font_weight_t  _weight;    /*!< \brief Weight of font.
				      *
				      * Either CAIRO_FONT_WEIGHT_NORMAL or CAIRO_FONT_WEIGHT_BOLD.
				      */
    Vec3D                _color;     /*!< \brief Text color. */

    double               _xalign;    /*!< \brief Alignment of label in x-direction. */
    double               _yalign;    /*!< \brief Alignment of label in y-direction. */
    bool                 _yzeroext;  /*!< \brief Make align of y-direction using extents of zero. */
    double               _rotation;  /*!< \brief Rotation of label. */
    double               _xlocation; /*!< \brief Location of label in x-direction. */
    double               _ylocation; /*!< \brief Location of label in y-direction. */


    void process_parsed( cairo_t *cairo, const std::string &text, cairo_text_extents_t *extents0, 
			 double x0, double y0, double &x, double &y ) const;
    void parse_latex( cairo_t *cairo, const std::string &text, cairo_text_extents_t *extents0, 
		      double x0, double y0, double &x, double &y ) const;

public:

    Label();

    Label( const Label &label );

    Label( const std::string &text );

    ~Label();
    
    Label &operator=( const Label &label );

    /*! \brief Set label font size. 
     */
    void set_font_size( double size );

    /*! \brief Get label font size. 
     */
    double get_font_size( void ) const;

    /*! \brief Set label font family.
     */
    void set_font_family( const std::string &family );

    /*! \brief Set label font slant.
     */
    void set_font_slant(cairo_font_slant_t slant );
    
    /*! \brief Set label font weight.
     */
    void set_font_weight( cairo_font_weight_t weight );

    /*! \brief Set label color.
     */
    void set_color( const Vec3D &color );

    /*! \brief Set label location.
     */
    void set_location( double x, double y );

    /*! \brief Set label rotation.
     *
     *  Rotation is set in radians. The positive angle direction is
     *  counter-clockwise.
     */
    void set_rotation( double angle );

    /*! \brief Set label alignment.
     *
     *  Alignment of label relative to the set location
     *  point. Alignment (0,0) means the text is up and right from the
     *  point. (1,1) means the text is down and left from the
     *  point. Is either of aligment parameters is NaN or infinite,
     *  the text is laid so that the cursor starting position is at
     *  the set point.
     *
     *  If yzeroext is true, the label alignment in y-direction is
     *  made using the extents of character zero ("0") instead of the
     *  text itself.
     */
    void set_alignment( double x, double y, bool yzeroext = false );

    /*! \brief Set label text.
     */
    void set_text( const std::string &text );

    /*! \brief Get label text.
     */
    std::string get_text( void ) const;

    /*! \brief Draw label.
     */
    void draw( cairo_t *cairo );

    /*! \brief Get text extents of label.
     *
     *  The extents are independent of label rotation or alignment.
     */
    void get_extents( cairo_t *cairo, cairo_text_extents_t *extents );

    /*! \brief Get bounding box of label.
     *
     *  The bounding box takes in account the label rotation and
     *  aligment. The bounding box is always the size of the label
     *  in the direction of the axes. Bounding box is (xmin, ymin, xmax, ymax).
     */
    void get_bbox( cairo_t *cairo, double bbox[4] ) const;

    friend std::ostream &operator<<( std::ostream &os, const Label &label );
};


#endif
