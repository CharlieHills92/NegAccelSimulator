/*! \file palette.hpp
 *  \brief %Palette for colormaps
 */

/* Copyright (c) 2005-2013 Taneli Kalvas. All rights reserved.
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

#ifndef PALETTE_HPP
#define PALETTE_HPP 1


#include <vector>
#include "vec3d.hpp"


/*! \brief Class for palette definition.
 *
 *  %Palette is an object that contains a list of colors and
 *  corresponding values. The colors are interpolated linearly between
 *  the defined points for a smooth color palette. The palette values
 *  are normalized so that the end points of palette have values 0.0
 *  and 1.0.
 */
class Palette {    

public:

    /*! \brief %Palette defining entry.
     *
     *  An entry containing color-value pair, which are used for
     *  defining palettes.
     */
   struct Entry {
	Vec3D  _color; /*!< \brief Entry colors */
	double _val;   /*!< \brief Entry values. */

	Entry( const Vec3D &color, double val );

	bool operator<( const Entry &e ) const;
    };

private:

    std::vector<Entry> _entries;    /*!< \brief Palette entries. */
    int                _steps;      /*!< \brief Number of shades, less
				     *   than 1 for smooth palette,
				     *   larger for discrete steps. */

public:

    /*! \brief Default constructor for default black and white palette.
     */
    Palette();

    /*! \brief Constructor for defined palette.
     *
     *  The values are normed to range from 0.0 to 1.0.
     */
    Palette( const std::vector<Entry> &entries );

    /*! \brief Return the interpolated color value from palette.
     *
     *  Makes and interpolated color value at value \a x, where 0 <= x
     *  <= 1. If palette has no colors, black will be returned. If
     *  palette has one color, that color will be returned. With two
     *  or more colors the value returned is interpolated from the
     *  colors. Outside the defined range, the closest color value is
     *  returned (color 0 if x < 0 and color N-1 if x > 1). If stepped
     *  palette is enabled, the palette shades will be limited and
     *  hard limits will be shown on palette sweeps.
     */
    Vec3D operator()( double x ) const;

    /*! \brief Clear current palette.
     *
     *  Leaves palette with no colors. Used with push_back() to build
     *  new palettes on-line.
     */
    void clear( void );

    /*! \brief Pushes new entry to palette.
     *
     *  Palette entries are automatically sorted. The palette won't be
     *  normalized in range. This has to be manually done by calling
     *  normalize() after adding palette entries.
     */
    void push_back( const Vec3D &color, double val );

    /*! \brief Normalize palette entries.
     *
     *  Normalize palette to range from 0.0 to 1.0.
     */
    void normalize( void );

    /*! \brief Get palette steps.
     *
     *  If \a steps is less than or equal to 1 a regular interpolated
     *  palette will be used (default), otherwise \a steps is used as
     *  the number of separate shades in the palette.
     */
    int get_steps( void ) const;

    /*! \brief Set palette steps.
     *
     *  If \a steps is less than or equal to 1 a regular interpolated
     *  palette will be used (default), otherwise \a steps is used as
     *  the number of separate shades in the palette.
     *
     *  Defaults to 0.
     */
    void set_steps( int steps );

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;
};


#endif



















