/*! \file epot_umfpacksolver.hpp
 *  \brief UMFPACK matrix solver for electric potential problem
 */

/* Copyright (c) 2011-2013 Taneli Kalvas. All rights reserved.
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


#ifndef EPOT_UMFPACKSOLVER_HPP
#define EPOT_UMFPACKSOLVER_HPP 1


#include "epot_matrixsolver.hpp"
#include "ccolmatrix.hpp"


/*! \brief UMFPACK matrix solver for Electric potential problem.
 */
class EpotUMFPACKSolver : public EpotMatrixSolver {

    void     *_numeric;         /*!< \brief Numeric data for LU decomposition. */

    double    _newton_res;      /*!< \brief Newton residual error. */
    double    _newton_step;     /*!< \brief Newton step size. */

    bool      _gnewton;         /*!< \brief Globally convergent version of Newton-Raphson. */
    double    _newton_r_eps;    /*!< \brief Accuracy request for Newton-Raphson residual. */
    double    _newton_step_eps; /*!< \brief Accuracy request for Newton-Raphson step. */
    uint32_t  _newton_imax;     /*!< \brief Maximum number of Newton-Raphson iterations. */


    /*! \brief Reset solver/problem settings.
     */
    virtual void reset_problem( void );

    /*! \brief Solve problem with given mesh based space charge.
     */
    virtual void subsolve( MeshScalarField &epot, const MeshScalarField &scharge );

    static void umfpack_error( const std::string func, int status );

    void umfpack_decompose( const CColMatrix &mat );
    void umfpack_solve( const CColMatrix &mat, const Vector &rhs, Vector &sol,
			bool force_decomposition = false );

public:

    /*! \brief Constructor.
     */
    EpotUMFPACKSolver( Geometry &geom, 
		       double newton_r_eps = 1.0e-5, 
		       double newton_step_eps = 1.0e-6, 
		       uint32_t newton_imax = 10,
		       bool gnewton = true );

    /*! \brief Construct from file.
     */
    EpotUMFPACKSolver( Geometry &geom, std::istream &s );

    /*! \brief Destructor.
     */
    virtual ~EpotUMFPACKSolver();

    /*! \brief Sets maximum iteration count for Newton-Raphson steps.
     */
    void set_newton_imax( uint32_t newton_imax );

    /*! \brief Enable/disable globally convergent Newton-Raphson.
     *
     *  Enabled by default.
     */
    void set_gnewton( bool enable );

    /*! \brief Sets the accuracy request for Newton-Raphson residual.
     */
    void set_newton_residual_eps( double newton_r_eps );

    /*! \brief Get last Newton-Raphson residual.
     */
    double get_newton_residual( void ) const;

    /*! \brief Sets the accuracy request for Newton-Raphson step size.
     */
    void set_newton_step_eps( double newton_step_eps );

    /*! \brief Get last Newton-Raphson step size.
     */
    double get_newton_step( void ) const;

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;

    /*! \brief Saves problem data to stream.
     */
    virtual void save( std::ostream &s ) const;
};


#endif
