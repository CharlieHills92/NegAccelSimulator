/*! \file epot_gssolver.hpp
 *  \brief Gauss-Seidel solver for electric potential problem
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

#ifndef EPOT_GSSOLVER_HPP
#define EPOT_GSSOLVER_HPP 1


#include "epot_solver.hpp"


/*! \brief Gauss-Seidel solver for Electric potential problem.
 */
class EpotGSSolver : public EpotSolver {

    MeshScalarField *_epot;
    MeshScalarField *_rhs;

    uint32_t         _iter;           /*!< \brief Number of iteration rounds done. */
    uint32_t         _imax;           /*!< \brief Maximum number of iteration rounds. */
    double           _eps;            /*!< \brief Accuracy request. */
    double           _step;           /*!< \brief Potential change norm. */
    double           _err;            /*!< \brief Error estimate. */
    double           _w;              /*!< \brief Relaxation coefficient. */
    double           _local_Ulim_fac; /*!< \brief Factor to calculate potential limit. */
    double           _local_Ulim;     /*!< \brief Potential limit for plasma calculation in local solver. */
    uint32_t         _local_imax;     /*!< \brief Maximum iterations for plasma calculation in local solver. */
    double           _local_eps;      /*!< \brief Convergence limit for plasma calculation in local solver. */
    
    double solve_shield_potential( double epf, double cof, double rhs, double p ) const;
    double solve_pexp_potential( double epf, double cof, double rhs, double p ) const;
    double solve_nsimp_potential( double epf, double cof, double rhs, double p ) const;
    void prepare_local_gnewton_settings( void );

    double gs_loop_3d( void ) const;
    double gs_process_near_solid_3d( const uint8_t *nearsolid_ptr, 
				     uint32_t a, uint32_t dj, uint32_t dk, uint8_t bindex ) const;
    double gs_process_pure_vacuum_3d( uint32_t a, uint32_t dj, uint32_t dk ) const;
    double gs_process_neumann_3d( uint32_t a, uint32_t dj, uint32_t dk, uint8_t bindex ) const;


    double gs_loop_cyl( void ) const;
    double gs_process_near_solid_cyl( const uint8_t *nearsolid_ptr, 
				      uint32_t i, uint32_t j, uint8_t bindex ) const;
    double gs_process_pure_vacuum_cyl( uint32_t i, uint32_t j ) const;
    double gs_process_neumann_cyl( uint32_t i, uint32_t j, uint8_t bindex ) const;


    double gs_loop_2d( void ) const;
    double gs_process_near_solid_2d( const uint8_t *nearsolid_ptr, 
				     uint32_t a, uint32_t dj,
				     uint8_t bindex ) const;
    double gs_process_pure_vacuum_2d( uint32_t a, uint32_t dj ) const;
    double gs_process_neumann_2d( uint32_t a, uint32_t dj, uint8_t bindex ) const;


    double gs_loop_1d( void ) const;
    double gs_process_near_solid_1d( const uint8_t *nearsolid_ptr, 
				     uint32_t i, uint8_t bindex ) const;
    double gs_process_pure_vacuum_1d( uint32_t i ) const;
    double gs_process_neumann_1d( uint32_t i, uint8_t bindex ) const;


    double near_solid_neumann_rhs_contribution( uint32_t i, uint32_t j, uint32_t k, 
						uint8_t bindex, const Vec3D &x ) const;
    void preprocess( const MeshScalarField &scharge );
    void postprocess( void );

    /*! \brief Return scaling coefficient.
     *
     *  Returns an estimate for converting potential change norm to
     *  relative error estimate norm.
     */
    double error_scale( void );

    /*! \brief Reset solver/problem settings.
     */
    virtual void reset_problem( void );

    /*! \brief Solve problem with given mesh based space charge.
     */
    virtual void subsolve( MeshScalarField &epot, const MeshScalarField &scharge );

public:

    /*! \brief Constructor.
     */
    EpotGSSolver( Geometry &geom );

    /*! \brief Construct from file.
     */
    EpotGSSolver( Geometry &geom, std::istream &s );

    /*! \brief Destructor.
     */
    virtual ~EpotGSSolver() {}

    /*! \brief Sets the accuracy request.
     *
     *  Defaults to 1.0e-4.
     */
    void set_eps( double eps );

    /*! \brief Sets maximum iteration count.
     *
     *  Defaults to 100000.
     */
    void set_imax( uint32_t imax );

    /*! \brief Sets relaxation parameter.
     *
     *  Defaults to 1.66.
     */
    void set_w( double w );

    /*! \brief Set node wise plasma solver parameters
     *
     *  The \a Ulim_fac gives the potential limit within which the
     *  nonlinear solver is used in multiples of \a Te or maximum
     *  compensating particle energy from the plasma potential
     *  (defaults to 10). The \a imax gives the maximum number of
     *  iterations done at each node (defaults to 1) and eps gives the
     *  convergence criterion (default to 1.0e-6 volts).
     */
    void set_plasma_solver_parameters( double Ulim_fac, uint32_t imax, double eps );

    /*! \brief Get potential change norm.
     *
     *  Returns 2-norm \f$ ||\Delta x|| \f$.
     */
    double get_potential_change_norm( void ) const;

    /*! \brief Get estimate of relative solution error.
     *
     *  Returns 2-norm \f$ 1e-3 \max(I,J,K) ||\Delta x|| \f$ in 2D and
     *  \f$ 1e-3 sqrt{\max(I,J,K)} ||\Delta x|| \f$ in 3D.
     */
    double get_error_estimate( void ) const;

    /*! \brief Get number of iteration rounds done with last solve().
     */
    uint32_t get_iter( void ) const;

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;

    /*! \brief Saves problem data to stream.
     */
    virtual void save( std::ostream &s ) const;
};

#endif
