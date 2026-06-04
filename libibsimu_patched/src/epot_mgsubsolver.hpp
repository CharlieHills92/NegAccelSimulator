/*! \file epot_mgsolver.hpp
 *  \brief Multigrid solver for electric potential problem
 */

/* Copyright (c) 2011-2013,2021 Taneli Kalvas. All rights reserved.
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

#ifndef EPOT_MGSUBSOLVER_HPP
#define EPOT_MGSUBSOLVER_HPP 1


#include "epot_solver.hpp"


/*! \brief Subroutine class for EpotMGSolver.
 *
 *  Preprocesses the solid mesh and does relaxation rounds on one
 *  problem level.
 */
class EpotMGSubSolver : public EpotSolver {

    MeshScalarField        *_defect;
    MeshScalarField        *_epot;
    const MeshScalarField  *_rhs;

    double                  _Ulim;
    uint32_t                _imax;
    double                  _eps;

    virtual void reset_problem( void ) {}
    virtual void subsolve( MeshScalarField &epot, const MeshScalarField &scharge ) {}

    double solve_shield_potential( double epf, double cof, double rhs, double p ) const;
    double solve_nsimp_potential( double epf, double cof, double rhs, double p ) const;
    double solve_pexp_potential( double epf, double cof, double rhs, double p ) const;

    // 1D
    double rbgs_loop_1d( void ) const;
    double sor_loop_1d( double w ) const;
    double gs_process_near_solid_1d( const uint8_t *nearsolid_ptr, 
				     uint32_t i, uint8_t bindex ) const;
    double gs_process_pure_vacuum_1d( uint32_t i ) const;
    double gs_process_neumann_1d( uint32_t i, uint8_t bindex ) const;

    void   defect_1d( void ) const;
    double defect_near_solid_1d( const uint8_t *nearsolid_ptr, uint32_t i, uint8_t bindex ) const;
    double defect_pure_vacuum_1d( uint32_t i ) const;
    double defect_neumann_1d( uint32_t i, uint8_t bindex ) const;

    // 2D
    double rbgs_loop_2d( void ) const;
    double sor_loop_2d( double w ) const;
    double gs_process_near_solid_2d( const uint8_t *nearsolid_ptr, 
				     uint32_t a, uint32_t dj, uint8_t bindex ) const;
    double gs_process_pure_vacuum_2d( uint32_t a, uint32_t dj ) const;
    double gs_process_neumann_2d( uint32_t a, uint32_t dj, uint8_t bindex ) const;

    void   defect_2d( void ) const;
    double defect_near_solid_2d( const uint8_t *nearsolid_ptr, 
				 uint32_t a, uint32_t dj, uint8_t bindex ) const;
    double defect_pure_vacuum_2d( uint32_t a, uint32_t dj ) const;
    double defect_neumann_2d( uint32_t a, uint32_t dj, uint8_t bindex ) const;

    // CYL
    double rbgs_loop_cyl( void ) const;
    double sor_loop_cyl( double w ) const;
    double gs_process_near_solid_cyl( const uint8_t *nearsolid_ptr, 
				      uint32_t i, uint32_t j, uint8_t bindex ) const;
    double gs_process_pure_vacuum_cyl( uint32_t i, uint32_t j ) const;
    double gs_process_neumann_cyl( uint32_t i, uint32_t j, uint8_t bindex ) const;

    void   defect_cyl( void ) const;
    double defect_near_solid_cyl( const uint8_t *nearsolid_ptr, 
				  uint32_t i, uint32_t j, uint8_t bindex ) const;
    double defect_pure_vacuum_cyl( uint32_t i, uint32_t j ) const;
    double defect_neumann_cyl( uint32_t i, uint32_t j, uint8_t bindex ) const;

    // 3D
    double rbgs_loop_3d( void ) const;
    double sor_loop_3d( double w ) const;
    double gs_process_near_solid_3d( const uint8_t *nearsolid_ptr, 
				     uint32_t a, uint32_t dj, uint32_t dk, uint8_t bindex ) const;
    double gs_process_pure_vacuum_3d( uint32_t a, uint32_t dj, uint32_t dk ) const;
    double gs_process_neumann_3d( uint32_t a, uint32_t dj, uint32_t dk, uint8_t bindex ) const;

    void   defect_3d( bool after_smooth ) const;
    double defect_near_solid_3d( const uint8_t *nearsolid_ptr, 
				 uint32_t a, uint32_t dj, uint32_t dk, uint8_t bindex ) const;
    double defect_pure_vacuum_3d( uint32_t a, uint32_t dj, uint32_t dk ) const;
    double defect_neumann_3d( uint32_t a, uint32_t dj, uint32_t dk, uint8_t bindex ) const;

public:

    /*! \brief Constructor.
     *
     *  Construct subsolver for geometry \a geom. Use parameters from
     *  main level potential solver \a epsolver.
     */
    EpotMGSubSolver( const EpotSolver &epsolver, Geometry &geom,
		     double Ulim, uint32_t imax, double eps );

    /*! \brief Destructor.
     */
    virtual ~EpotMGSubSolver() {}

    /*! \brief Calculate defect
     *
     *  If calculating defect after RBGS smooth, the odd points are
     *  known to have zero defect.
     */
    void defect( MeshScalarField *defect, MeshScalarField *epot, const MeshScalarField *rhs,
		 bool after_smooth );

    /*! \brief Do a smoothing round with Red-Black Gauss-Seidel.
     */
    double mg_smooth( MeshScalarField *epot, const MeshScalarField *rhs );

    /*! \brief Do a solve round with SOR using over-relaxation factor \a w.
     */
    double mg_solve( MeshScalarField *epot, const MeshScalarField *rhs, double w );

    /*! \brief Return error scaling factor for SOR solver with over-relaxation factor w.
     */
    double error_scale( double w ) const;

    /*! \brief Return error scaling factor for MG smoother.
     */
    double error_scale_mg( void ) const;

    /*! \brief Preprocess.
     */
    void preprocess( MeshScalarField &epot );

    /*! \brief Postprocess.
     */
    void postprocess( void );

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const {}

    /*! \brief Saves problem data to stream.
     */
    virtual void save( std::ostream &s ) const {}
};


#endif
