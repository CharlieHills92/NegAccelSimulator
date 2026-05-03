/*! \file coordmapper.hpp
 *  \brief 1D and 2D coordinate transformations for plotter.
 */

/* Copyright (c) 2005-2010 Taneli Kalvas. All rights reserved.
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

#ifndef COORDMAPPER_HPP
#define COORDMAPPER_HPP 1


#include <iostream>


/*! \brief Linear 1D coordinate mapper.
 *
 *  Coordinate transformation is done like
 *
 *  x_new = xx * x + x0;
 *
 */
class Coordmapper1D {
    double _xx, _x0;   /*!< \brief Transformation coeffiecients for x-axis */

public:

    /*! \brief Default constructor for unitary transformation.
     */
    Coordmapper1D()
	: _xx(1.0), _x0(0.0) {}

    /*! \brief Constructor for fully defined transformation.
     */
    Coordmapper1D( double xx, double x0 )
	: _xx(xx), _x0(x0) {}

    /*! \brief Set transformation coefficients.
     */
    void set_transformation( double xx, double x0 ) {
	_xx = xx;
	_x0 = x0;
    }

    /*! \brief Make transformation from coordinates \a xin to
     *  coordinates \a xout.
     */
    void transform( double &xout, const double &xin ) const {
	xout = _xx * xin + _x0;
    }

    /*! \brief Make transformation from coordinates \a xin to
     *  coordinates \a xout.
     */
    void transform( double &x ) const {
	x = _xx * x + _x0;
    }

    /*! \brief Make inverse transformation for coordinate \a x.
     */
    void inv_transform( double &xout, const double &xin ) const {
	xout = (xin-_x0) / _xx;
    }

    /*! \brief Make inverse transformation for coordinate \a x.
     */
    void inv_transform( double &x ) const {
	x = (x-_x0) / _xx;
    }

    /*! \brief Debug print to stream.
     */ 
    void debug_print( std::ostream &os ) const {
	os << "**Coordmapper1D\n";
	os << "_xx = " << _xx << "\n";
	os << "_x0 = " << _x0 << "\n";
    }
};



/*! \brief Linear-linear 2D coordinate mapper.
 *
 */
class Coordmapper {
    Coordmapper1D _cmx;  /*!< \brief Transformation for x-axis */
    Coordmapper1D _cmy;  /*!< \brief Transformation for y-axis */

public:

    /*! \brief Default constructor for unitary transformation.
     */
    Coordmapper() {}

    /*! \brief Constructor for fully defined transformation using 1D
     *  coordinate mappers.
     */
    Coordmapper( Coordmapper1D cmx, Coordmapper1D cmy )
	: _cmx(cmx), _cmy(cmy) {}

    /*! \brief Constructor for fully defined transformation.
     */
    Coordmapper( double xx, double x0, double yy, double y0 )
	: _cmx(xx,x0), _cmy(yy,y0) {}

    /*! \brief Set transformation matrix coefficients.
     */
    void set_transformation( double xx, double x0, double yy, double y0 ) {
	_cmx.set_transformation( xx, x0 );
	_cmy.set_transformation( yy, y0 );
    }

    /*! \brief Make transformation for coordinates \a x, \a y.
     */
    void transform( double &x, double &y ) const {
	_cmx.transform( x );
	_cmy.transform( y );
    }

    /*! \brief Make transformation from coordinates \a xin to
     *  coordinates \a xout.
     */
    void transform( double *xout, const double *xin ) const {
	_cmx.transform( xout[0], xin[0] );
	_cmy.transform( xout[1], xin[1] );
    }

    /*! \brief Make inverse transformation for coordinates \a x, \a y.
     */
    void inv_transform( double &x, double &y ) const {
	_cmx.inv_transform( x );
	_cmy.inv_transform( y );
    }

    /*! \brief Make inverse transformation from coordinates \a xin to
     *  coordinates \a xout.
     */
    void inv_transform( double *xout, const double *xin ) const {
	_cmx.inv_transform( xout[0], xin[0] );
	_cmy.inv_transform( xout[1], xin[1] );
    }


    /*! \brief Debug print to stream.
     */ 
    void debug_print( std::ostream &os ) const {
	os << "**Coordmapper x:\n";
	_cmx.debug_print( os );
	os << "**Coordmapper y:\n";
	_cmy.debug_print( os );
    }
};


#endif





















