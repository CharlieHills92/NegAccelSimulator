/*! \file epot_mgsolver.hpp
 *  \brief Multigrid solver for electric potential problem
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

#ifndef EPOT_MGSOLVER_HPP
#define EPOT_MGSOLVER_HPP 1


#include "epot_solver.hpp"
#include "epot_mgsubsolver.hpp"


/*! \brief Multigrid solver for Electric potential problem.
 */
class EpotMGSolver : public EpotSolver {

    std::vector<MeshScalarField *>   _epotv;
    std::vector<Geometry *>          _geomv;
    std::vector<EpotMGSubSolver *>   _epotsolverv;
    std::vector<MeshScalarField *>   _rhsv;
    std::vector<MeshScalarField *>   _workv;
    std::vector<MeshScalarField *>   _work2v;

    bool             _geom_prepared;  /*!< \brief Is geometry prepared? */
    uint32_t         _levels;         /*!< \brief Multigrid levels. */
    uint32_t         _npre;           /*!< \brief Pre cycle smoother rounds. */
    uint32_t         _npost;          /*!< \brief Post cycle smoother rounds. */
    uint32_t         _mgcycmax;       /*!< \brief Maximum number of MG cycles. */
    uint32_t         _mgcyc;          /*!< \brief Number of MG cycles taken. */
    double           _mgeps;          /*!< \brief Acceptable residual error from last MG cycle. */
    double           _mgerr;          /*!< \brief Error estimate from top level. */
    uint32_t         _gamma;          /*!< \brief Multigrid cycle coefficient, 1 for V-cycles, 2 for W-cycles. */
    double           _step;           /*!< \brief Potential change norm from top level. */
    double           _coarse_eps;     /*!< \brief Acceptable error for coarsest level. */
    double           _coarse_err;     /*!< \brief Acceptable error for coarsest level. */
    double           _w;              /*!< \brief Over-relaxation factor for coarsest level. */
    uint32_t         _imax;           /*!< \brief Maximum number of rounds for coarsest level. */
    uint32_t         _coarse_steps;   /*!< \brief Total number of rounds for coarsest level in last MG cycle. */
    double           _local_Ulim;     /*!< \brief Potential limit for plasma calculation in local solver. */
    uint32_t         _local_imax;     /*!< \brief Maximum iterations for plasma calculation in local solver. */
    double           _local_eps;      /*!< \brief Convergence limit for plasma calculation in local solver. */

    void prepare_local_gnewton_settings( void );
    void print_field( const MeshScalarField *F );

    void prepare_mg_geom( void );

    double near_solid_neumann_rhs_contribution( uint32_t i, uint32_t j, uint32_t k, 
						uint8_t bindex, const Vec3D &x ) const;
    void preprocess( MeshScalarField &epot, const MeshScalarField &scharge );
    void postprocess( void );

    void prolong_add_3d( MeshScalarField *out, int32_t i, int32_t j, int32_t k, double C );
    void prolong_add_cyl( MeshScalarField *out, int32_t i, int32_t j, double C );
    void prolong_add_2d( MeshScalarField *out, int32_t i, int32_t j, double C );
    void prolong_add_1d( MeshScalarField *out, int32_t i, double C );

    void correct( const Geometry *geom, MeshScalarField *sol, const MeshScalarField *corr );

    void restrict_3d( MeshScalarField *out, const MeshScalarField *in, bool defect );
    void restrict_cyl( MeshScalarField *out, const MeshScalarField *in, bool defect );
    void restrict_2d( MeshScalarField *out, const MeshScalarField *in, bool defect );
    void restrict_1d( MeshScalarField *out, const MeshScalarField *in, bool defect );

    /*! \brief Restrict field.
     *
     *  If restricting defect field it is known that odd points have zero defect.
     */
    void restrict( MeshScalarField *out, const MeshScalarField *in, bool defect );

    void prolong_3d( MeshScalarField *out, const MeshScalarField *in );
    void prolong_cyl( MeshScalarField *out, const MeshScalarField *in );
    void prolong_2d( MeshScalarField *out, const MeshScalarField *in );
    void prolong_1d( MeshScalarField *out, const MeshScalarField *in );
    void prolong( MeshScalarField *out, const MeshScalarField *in );

    void mg_recurse( uint32_t level );

    /*! \brief Reset solver/problem settings.
     */
    virtual void reset_problem( void );

    /*! \brief Solve problem with given mesh based space charge.
     */
    virtual void subsolve( MeshScalarField &epot, const MeshScalarField &scharge );

public:

    /*! \brief Constructor.
     */
    EpotMGSolver( Geometry &geom );

    /*! \brief Construct from file.
     */
    EpotMGSolver( Geometry &geom, std::istream &s );

    /*! \brief Destructor.
     */
    virtual ~EpotMGSolver();

    /*! \brief Sets the accuracy request for coarsest level SOR solver.
     *
     *  Defaults to 1e-10.
     */
    void set_eps( double eps );

    /*! \brief Sets the over-relaxation factor for coarsest level SOR solver.
     *
     *  Defaults to 1.7.
     */
    void set_w( double w );

    /*! \brief Sets maximum number of iteration rounds for coarsest level SOR solver.
     *
     *  Defaults to 10000.
     */
    void set_imax( uint32_t imax );

    /*! \brief Sets maximum number of local iterations to take for nonlinear problems.
     *
     *  Defaults to 1.
     */
    void set_local_imax( uint32_t local_imax );

    /*! \brief Sets multigrid levels.
     *
     *  Defaults to 1.
     */
    void set_levels( uint32_t levels );

    /*! \brief Sets maximum number of multigrid cycles to take.
     *
     *  Defaults to 100. Multigrid cycles are done until the
     *  residual error is less than \a mgeps or \a mgcyc cycles have
     *  been made. 
     */
    void set_mgcycmax( uint32_t mgcyc );

    /*! \brief Sets the accuracy request for finest level.
     *
     *  Defaults to 1.0e-4. Multigrid cycles are done until the
     *  residual error is less than \a mgeps or \a mgcyc cycles have
     *  been made.
     */
    void set_mgeps( double mgeps );

    /*! \brief Sets multigrid cycle coefficient.
     *
     *  Defaults to 1 (V-cycles). Use 2 for W-cycles.
     */
    void set_gamma( uint32_t gamma );

    /*! \brief Sets number of pre cycle smoother rounds.
     *
     *  Defaults to 5.
     */
    void set_npre( uint32_t npre );

    /*! \brief Sets number of post cycle smoother rounds.
     *
     *  Defaults to 5.
     */
    void set_npost( uint32_t npost );

    /*! \brief Get potential change norm.
     *
     *  Returns 2-norm \f$ ||\Delta x|| \f$.
     */
    double get_potential_change_norm( void ) const;

    /*! \brief Get estimate of relative solution error.
     *
     *  Returns 2-norm \f$ G(w) \max(I,J,K) ||\Delta x|| \f$ in 2D and
     *  \f$ F(w) sqrt{\max(I,J,K)} ||\Delta x|| \f$ in 3D, where G and
     *  F and third order polynomial fits to data from test problems.
     */
    double get_error_estimate( void ) const;

    /*! \brief Get number of multigrid cycles done.
     */
    uint32_t get_mgcyc( void ) const;

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;

    /*! \brief Saves problem data to stream.
     */
    virtual void save( std::ostream &s ) const;
};

#endif
