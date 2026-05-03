/*! \file axisymmetricvectorfield.hpp
 *  \brief Axisymmetric magnetic field
 */

/* Copyright (c) 2012,2014 Taneli Kalvas. All rights reserved.
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

#ifndef AXISYMMETRICVECTORFIELD_HPP
#define AXISYMMTERICVECTORFIELD_HPP 1


#include <vector>
#include <gsl/gsl_spline.h>
#include "types.hpp"
#include "vectorfield.hpp"


/*! \brief Axisymmetric magnetic field based on on-axis data.
 *
 *  The magnetic field is constructed using cubic spline
 *  interpolation of user given data on axis and the expansion of the
 *  field off-axis with
 *  \f{eqnarray*} B_z(r,z) &=& B_z(0,z) - \frac{r^2}{4} B_z''(0,z) \\
    B_r(r,z) &=& -\frac{r}{2} B_z'(0,z)
 *  \f}
 *  where the third and higher derivatives have been truncated (see J. Rodney, 
 *  M. Vaughan, "Representation of Axisymmetric Magnetic Fields in Computer 
 *  Programs", IEEE Transactions on electron devices 19, 144 (1972) for 
 *  more details).
 */
class AxisymmetricVectorField : public VectorField {

    geom_mode_e         _geom_mode;
    gsl_spline         *_spline;        /*!< \brief Spline used to interpolate on-axis data */
    gsl_interp_accel   *_accel;         /*!< \brief Interpolation accelerator */

    double              _origo;
    double              _h;
    std::vector<double> _Bz;
    uint32_t            _order;

    void get_fdm_derivatives( double Bder[6], double z ) const;

    const Vec3D eval_spline( double z, double r ) const;
    const Vec3D eval_fdm( double z, double r ) const;

public:

    /*! \brief Constructor.
     *
     *  Constructor for axisymmetric magnetic field based on an array
     *  of magnetic field data on the symmetry axis. If geom_mode is
     *  MODE_CYL the symmetry axis is x. If geom_mode is MODE_3D, the
     *  symmetry axis is z. Other geometry modes are not supported.
     *
     *  The magnetic field is constructed using cubic spline
     *  interpolation of user given data on axis and the expansion of the
     *  field off-axis with
     *  \f{eqnarray*} B_z(r,z) &=& B_z(0,z) - \frac{r^2}{4} B_z''(0,z) \\
                      B_r(r,z) &=& -\frac{r}{2} B_z'(0,z),
     *  \f}
     *  where the third and higher derivatives have been truncated.
     *
     *  The magnetic field outside the defined z is zero. No
     *  extrapolation is done.
     */
    AxisymmetricVectorField( geom_mode_e geom_mode, 
			     const std::vector<double> &z, 
			     const std::vector<double> &Bz );

    /*! \brief Constructor.
     *
     *  Constructor for axisymmetric magnetic field based on an array
     *  of magnetic field data on the symmetry axis. If geom_mode is
     *  MODE_CYL the symmetry axis is x. If geom_mode is MODE_3D, the
     *  symmetry axis is z. Other geometry modes are not supported.
     *
     *  The magnetic field is constructed using linear interpolation
     *  and finite difference derivatives up to 6th order for the Bz
     *  data on axis and the expansion of the field off-axis with
     *  \f{eqnarray*} B_z(r,z) &=& B_z(0,z) - \frac{r^2}{4} B_z''(0,z) + \frac{r^4}{64} B_z^{(iv)}(0,z) - \frac{r^6}{2304} B_z^{(vi)}(0,z) \\
                      B_r(r,z) &=& -\frac{r}{2} B_z'(0,z) + \frac{r^3}{16} B_z''' - \frac{r^5}{384} B_z^{(v)}, 
     *  \f}
     *  The order used in evaluation is defined by parameter
     *  \a order, which defaults to the maximum allowed order of
     *  6.
     *
     *  The magnetic field outside the defined z is extrapolated.
     */
    AxisymmetricVectorField( geom_mode_e geom_mode, 
			     double origo, double h,
			     const std::vector<double> &Bz,
			     uint32_t order = 6 );

    /*! \brief Copy constructor. 
     */
    AxisymmetricVectorField( const AxisymmetricVectorField &f );

    /*! \brief Virtual destructor.
     */
    virtual ~AxisymmetricVectorField();

    /*! \brief Operator for getting field value at \a x.
     */
    virtual const Vec3D operator()( const Vec3D &x ) const;

    /*! \brief Copy assignment. */
    AxisymmetricVectorField &operator=( const AxisymmetricVectorField &f );
};


#endif
