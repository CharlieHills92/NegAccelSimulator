/*! \file epot_bicgstabsolver.hpp
 *  \brief BiCGSTAB matrix solver for electric potential problem
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


#ifndef EPOT_BICGSTABSOLVER_HPP
#define EPOT_BICGSTABSOLVER_HPP 1


#include "epot_matrixsolver.hpp"
#include "precond.hpp"


/*! \brief BiCGSTAB matrix solver for Electric potential problem.
 */
class EpotBiCGSTABSolver : public EpotMatrixSolver {

    double    _eps;             /*!< \brief Accuracy request. */
    uint32_t  _imax;            /*!< \brief Maximum iterations allowed. */

    uint32_t  _iter;            /*!< \brief Number of iteration rounds done. */
    double    _res;             /*!< \brief Scaled residual error. */
    double    _err;             /*!< \brief Error estimate. */

    double    _newton_res;      /*!< \brief Newton residual error max. */
    double    _newton_res_norm; /*!< \brief Newton residual error norm2. */
    double    _newton_step;     /*!< \brief Newton step size max. */
    double    _newton_step_norm;/*!< \brief Newton step size norm2. */

    bool      _gnewton;         /*!< \brief Globally convergent version of Newton-Raphson. */
    double    _newton_eps;      /*!< \brief Accuracy request for Newton-Raphson. */
    uint32_t  _newton_imax;     /*!< \brief Maximum number of Newton-Raphson iterations. */

    Precond  *_pc;              /*!< \brief Preconditioner. */

    MeshScalarField *_epot;     /*!< \brief Electric potential output. */
    void    (*_callback)(void); /*!< \brief Analysis callback. */
    void    (*_callback_nonlinear)(void); /*!< \brief Analysis callback. */

    /*! \brief Reset solver/problem settings.
     */
    virtual void reset_problem( void );

    /*! \brief Solve problem with given mesh based space charge.
     */
    virtual void subsolve( MeshScalarField &epot, const MeshScalarField &scharge );

    void bicgstab( const Matrix &mat, const Vector &rhs, Vector &sol,
		   const Precond &pc );

public:

    /*! \brief Constructor.
     */
    EpotBiCGSTABSolver( Geometry &geom,
			double eps = 1.0e-4, 
			uint32_t imax = 10000,
			double newton_eps = 1.0e-4, 
			uint32_t newton_imax = 10,
			bool gnewton = true );

    /*! \brief Destructor.
     */
    EpotBiCGSTABSolver( Geometry &geom, std::istream &s );

    /*! \brief Destructor.
     */
    virtual ~EpotBiCGSTABSolver();

    /*! \brief Set preconditioner to use.
     */
    void set_preconditioner( Precond &pc );

    /*! \brief Enable/disable globally convergent Newton-Raphson.
     *
     *  Enabled by default.
     */
    void set_gnewton( bool enable );

    /*! \brief Sets the accuracy request for BiCGSTAB solver.
     *
     *  Defaults to 1.0e-4.
     */
    void set_eps( double eps );

    /*! \brief Sets maximum iteration count for BiCGSTAB solver.
     *
     *  Defaults to 10000.
     */
    void set_imax( uint32_t imax );

    /*! \brief Sets maximum iteration count for Newton-Raphson steps.
     *
     *  Defaults to 10.
     */
    void set_newton_imax( uint32_t newton_imax );

    /*! \brief Sets the accuracy request for Newton-Raphson.
     *
     *  Based on error estimator. Defaults to 1.0e-4;
     */
    void set_newton_eps( double eps );

    /*! \brief Get last Newton-Raphson residual.
     */
    double get_newton_residual( void ) const;

    /*! \brief Get last Newton-Raphson residual norm.
     */
    double get_newton_residual_norm( void ) const;

    /*! \brief Get last Newton-Raphson step size.
     */
    double get_newton_step( void ) const;

    /*! \brief Get last Newton-Raphson step size norm.
     */
    double get_newton_step_norm( void ) const;

    /*! \brief Get scaled residual error.
     *
     *  Returns scaled 2-norm residual \f$ ||(A*x-b)|| / ||b|| \f$.
     *  Only available for linear solutions.
     */
    double get_scaled_residual( void ) const;

    /*! \brief Get estimate of relative solution error.
     *
     *  Returns an estimate for \f$ ||x-x^*|| / ||x^*|| \f$.  For
     *  linear solutions this is the scaled 2-norm residual 
     *  \f$ \max(I,J,K)^{-2} ||(A*x-b)|| / ||b|| \f$ and for nonlinear
     *  problem it is \f$ \max(I,J,K) ||R|| / 10^7 \f$. 
     */
    double get_error_estimate( void ) const;

    /*! \brief Get number of iteration rounds done with last solve().
     */
    uint32_t get_iter( void ) const;

    /*! \brief Set analysis callback.
     *
     *  If callback is set to non-NULL, the electric potential is
     *  constructed at each iteration and callback function is called.
     */
    void set_analysis_callback( void (*func)(void) );
    
    /*! \brief Set analysis callback for nonlinear iteration.
     *
     *  If callback is set to non-NULL, the electric potential is
     *  constructed at each newton step and callback function is called.
     */
    void set_analysis_callback_nonlinear( void (*func)(void) );
    
    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;

    /*! \brief Saves problem data to stream.
     */
    virtual void save( std::ostream &s ) const;
};


#endif
