/*! \file mydxfheader.hpp
 *  \brief DXF Header
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

#ifndef MY_DXF_HEADER_HPP
#define MY_DXF_HEADER_HPP 1


#include <stdint.h>
#include "vec3d.hpp"
#include "mydxffile.hpp"


/*! \brief DXF header class.
 *
 *  Container for data of a DXF file header. The data currently
 *  doesn't affect the functioning of the MyDXF library. All settings
 *  are currently ignored.
 */
class MyDXFHeader
{



public:

    std::string acadver;         /* 1 */
    double      angbase;         /* 50 */
    int16_t     angdir;          /* 70 */

    
    std::string handseed;        /* 5 */
    double      dimasz;          /* 40 */
    double      dimgap;          /* 40 */
    double      dimexo;          /* 40 */
    double      dimexe;          /* 40 */
    double      dimtxt;          /* 40 */
    int16_t     insunits;        /* 70 */
    Vec3D       plimmax;         /* 10, 20 */
    Vec3D       plimmin;         /* 10, 20 */


    int16_t     orthomode;       /* 70, on if nonzero */


    std::string pucsbase;        /* 2, Name of the UCS that defines the origin and orientation
                                  * of orthographic UCS settings (paper space only) */
    std::string pucsname;        /* 2, Current paper space UCS name */
    Vec3D       pucsorg;         /* 10, 20, 30, Current paper space UCS origin */
    Vec3D       pucsorgback;     /* 10, 20, 30, Point which becomes the new UCS origin after changing
				  * paper space UCS to BACK when PUCSBASE is set to WORLD */
    Vec3D       pucsorgbottom;   /* 10, 20, 30, Point which becomes the new UCS origin after changing
				  * paper space UCS to BOTTOM when PUCSBASE is set to WORLD */
    Vec3D       pucsorgfront;    /* 10, 20, 30, Point which becomes the new UCS origin after changing
				  * paper space UCS to FRONT when PUCSBASE is set to WORLD */
    Vec3D       pucsorgleft;     /* 10, 20, 30, Point which becomes the new UCS origin after changing
				  * paper space UCS to LEFT when PUCSBASE is set to WORLD */
    Vec3D       pucsorgright;    /* 10, 20, 30, Point which becomes the new UCS origin after changing
				  * paper space UCS to RIGHT when PUCSBASE is set to WORLD */
    Vec3D       pucsorgtop;      /* 10, 20, 30, Point which becomes the new UCS origin after changing
				  * paper space UCS to TOP when PUCSBASE is set to WORLD */
    std::string pucsorthoref;    /* 2, If paper space UCS is orthographic (PUCSORTHOVIEW not
				  * equal to 0), this is the name of the UCS that the orthographic
				  * UCS is relative to. If blank, UCS is relative to WORLD */
    int16_t     pucsorthoview;   /* 70, Orthographic view type of paper space UCS:
				  * 0 = UCS is not orthographic;
				  * 1 = Top; 2 = Bottom;
				  * 3 = Front; 4 = Back;
				  * 5 = Left; 6 = Right */
    Vec3D       pucsxdir;        /* 10, 20, 30, Current paper space UCS X axis */
    Vec3D       pucsydir;        /* 10, 20, 30, Current paper space UCS Y axis */


    std::string ucsbase;         /* 2, Name of the UCS that defines the origin and orientation
				  * of orthographic UCS settings */
    std::string ucsname;         /* 2, Name of current UCS */
    Vec3D       ucsorg;          /* 10, 20, 30, Origin of current UCS (in WCS) */
    Vec3D       ucsorgback;      /* 10, 20, 30, Point which becomes the new UCS origin after changing
				  * model space UCS to BACK when UCSBASE is set to WORLD */
    Vec3D       ucsorgbottom;    /* 10, 20, 30, Point which becomes the new UCS origin after changing
				  * model space UCS to BOTTOM when UCSBASE is set to WORLD */
    Vec3D       ucsorgfront;     /* 10, 20, 30, Point which becomes the new UCS origin after changing
				  * model space UCS to FRONT when UCSBASE is set to WORLD */
    Vec3D       ucsorgleft;      /* 10, 20, 30, Point which becomes the new UCS origin after changing
				  * model space UCS to LEFT when UCSBASE is set to WORLD */
    Vec3D       ucsorgright;     /* 10, 20, 30, Point which becomes the new UCS origin after changing
				  * model space UCS to RIGHT when UCSBASE is set to WORLD */
    Vec3D       ucsorgtop;       /* 10, 20, 30, Point which becomes the new UCS origin after changing
				  * model space UCS to TOP when UCSBASE is set to WORLD */
    std::string ucsorthoref;     /* 2, If model space UCS is orthographic (UCSORTHOVIEW not
				  * equal to 0), this is the name of the UCS that the orthographic
				  * UCS is relative to. If blank, UCS is relative to WORLD */
    int16_t     ucsorthoview;    /* 70, Orthographic view type of model space UCS:
				  * 0 = UCS is not orthographic;
				  * 1 = Top; 2 = Bottom;
				  * 3 = Front; 4 = Back;
				  * 5 = Left; 6 = Right */
    Vec3D       ucsxdir;         /* 10, 20, 30, Direction of the current UCS X axis (in WCS) */
    Vec3D       ucsydir;         /* 10, 20, 30, Direction of the current UCS Y axis (in WCS) */


    int16_t     worldview;       /* 70, 1 = Set UCS to WCS during DVIEW/VPOINT
				  * 0 = Don't change UCS */

    MyDXFHeader( class MyDXFFile *dxf );
    ~MyDXFHeader();

    /*! \brief Write dxf file to stream.
     */
    void write( class MyDXFFile *dxf, std::ofstream &_ostr );

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;
};





#endif



