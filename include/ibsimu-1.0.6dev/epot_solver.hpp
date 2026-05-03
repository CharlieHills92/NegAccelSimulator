/*! \file epot_solver.hpp
 *  \brief Poisson equation problem for solving electric potential.
 */

/* Copyright (c) 2005-2013,2017,2021 Taneli Kalvas. All rights reserved.
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

#ifndef EPOT_SOLVER_HPP
#define EPOT_SOLVER_HPP 1


#include <iostream>
#include <stdint.h>
#include "solver.hpp"
#include "callback.hpp"
#include "meshscalarfield.hpp"
#include "geometry.hpp"


/*! \brief Initial plasma volume definition.
 */
class InitialPlasma : public CallbackFunctorB_V {

    coordinate_axis_e _axis;
    double            _val;

public:

    /*! \brief Constructor setting initial plasma.
     *
     *  Initial plasma exists at coordinates less than \a val in \a
     *  axis direction.
     */
    InitialPlasma( coordinate_axis_e axis, double val ) 
        : _axis(axis), _val(val) {}

    /*! \brief Destructor.
     */
    ~InitialPlasma() {}

    /*! \brief Suppression function.
     */
    virtual bool operator()( const Vec3D &x ) const {
	switch( _axis ) {
	case AXIS_X:
	    if( x[0] < _val )
		return( true );
	    else 
		return( false );
	    break;
	case AXIS_Y:
	case AXIS_R:
	    if( x[1] < _val )
		return( true );
	    else 
		return( false );
	    break;
	case AXIS_Z:
	    if( x[2] < _val )
		return( true );
	    else 
		return( false );
	    break;
	}

	// Never here
	return( false );
    }
};


/*! \brief Plasma modes
 *
 *  Selection of modes for plasma calculation in electric potential
 *  problem. In a typical positive ion plasma extraction simulation
 *  the plasma mode is PLASMA_PEXP_INITIAL in the first iteration
 *  round to use initial guess for plasma meniscus
 *  location. Thereafter PLASMA_PEXP is used. For negative ion
 *  extraction the PLASMA_NSIMP_INITIAL is used for the first iteration
 *  and PLASMA_NSIMP thereafter. PLASMA_INITIAL is a macro, which
 *  equals to PLASMA_PEXP_INITIAL. It exists for backward compatibility.
 *  PLASMA_SHIELD is a plasma model implementing a shielding factor for
 *  space charge inside plasma.
 *
 */
enum plasma_mode_e {
    PLASMA_NONE = 0,
    PLASMA_PEXP_INITIAL,
    PLASMA_NSIMP_INITIAL,
    PLASMA_PEXP, 
    PLASMA_NSIMP,
    PLASMA_SHIELD
};


#define EPOT_SOLVER_BXMIN 1
#define EPOT_SOLVER_BXMAX 2
#define EPOT_SOLVER_BYMIN 4
#define EPOT_SOLVER_BYMAX 8
#define EPOT_SOLVER_BZMIN 16
#define EPOT_SOLVER_BZMAX 32


#define PLASMA_INITIAL PLASMA_PEXP_INITIAL

/*! \brief Class for constructing the linear/nonlinear problem for the
    solver.
 *
 *  %EpotProblem class constructs the Poisson (equation) problem in finite
 *  difference form from Geometry (mesh) and various parameters and it
 *  presents the problem to the Solver via matrix/vector
 *  representation. In case of linear problem (no plasma model), the
 *  Poisson equation 
 *  \f[ \nabla^2 \phi = -\frac{\rho}{\epsilon_0} \f]
 *  is
 *  \f[ \frac{\partial^2 \phi}{\partial x^2} = -\frac{\rho}{\epsilon_0} \f]
 *  in 1D, which is discretized into 
 *  \f[ \phi_{i-1} + \phi_{i+1} - 2\phi_{i} = -h^2 \frac{\rho_{i}}{\epsilon_0} \f]
 *  using finite differences. In 2D coordinates the discretized form is
 *  \f[ \phi_{i-1,j} + \phi_{i+1,j} + \phi_{i,j-1} + \phi_{i,j+1} - 4\phi_{i,j}
    = -h^2 \frac{\rho_{i,j}}{\epsilon_0} \f]
 *  and in 3D it is
 *  \f[ \phi_{i-1,j,k} + \phi_{i+1,j,k} + \phi_{i,j-1,k} + \phi_{i,j+1,k} 
    + \phi_{i,j,k-1} + \phi_{i,j,k+1} - 6\phi_{i,j} = -h^2 \frac{\rho_{i,j,k}}{\epsilon_0}. \f]
 *  In cylindrical coordinates the Poisson equation is
 *  \f[ \frac{\partial^2 \phi}{\partial r^2} + \frac{1}{r} \frac{\partial \phi}{\partial r}
    + \frac{1}{r^2} \frac{\partial^2 \phi}{\partial \theta^2} + \frac{\partial^2 \phi}{\partial z^2} 
    = -\frac{\rho}{\epsilon_0}, \f]
 *  where \f$ \frac{\partial^2 \phi}{\partial \theta^2} = 0 \f$ because of
 *  cylindrical symmetry of the simulations. Therefore the discretized form
 *  becomes
 *  \f[ \phi_{i-1,j} + \phi_{i+1,j} 
    + \left( 1 - \frac{h}{2r_j} \right) \phi_{i,j-1} 
    + \left( 1 + \frac{h}{2r_j} \right) \phi_{i,j+1} - 4\phi_{i,j}
    = -h^2 \frac{\rho_{i,j}}{\epsilon_0}, \f]
 *  where \f$ \frac{h}{2r_j} = \frac{h}{2jh} = \frac{1}{2j} \f$ because the 
 *  radius \f$ r_j = hj \f$. At the symmetry axis there is an exception because both 
 *  \f$ r \f$ and \f$ \frac{\partial \phi}{\partial r} \f$ approach zero. By 
 *  using Bernoulli-L'Hopital rule we can evaluate
 *  \f[ \lim_{r \rightarrow 0} \frac{1}{r} \frac{\partial \phi}{\partial r} 
    = \lim_{r \rightarrow 0} \frac{\partial^2 \phi}{\partial r^2}, \f]
 *  which turns the Poisson equation to
 *  \f[ 2 \frac{\partial^2 \phi}{\partial r^2} + \frac{\partial^2 \phi}{\partial z^2} 
    = -\frac{\rho}{\epsilon_0} \f]
 *  on axis. Discretation of the equation gives us
 *  \f[ \phi_{i-1,j} + \phi_{i+1,j} 
    + 2\phi_{i,j+1} + 2\phi_{i,j-1} - 6\phi_{i,j}
    = -h^2 \frac{\rho_{i,j}}{\epsilon_0}. \f]
 *  Here \f$ \phi_{i,j-1} = \phi_{i,j+1} \f$ so that the final form is
 *  \f[ \phi_{i-1,j} + \phi_{i+1,j} 
    + 4\phi_{i,j+1} - 6\phi_{i,j}
    = -h^2 \frac{\rho_{i,j}}{\epsilon_0}. \f]
 *
 *  In addition to the Poisson equation the problem matrix and vector
 *  also contain finite difference representations of the boundary
 *  conditions.  The Dirichlet boundary condition is defined by
 *  constant potential at boundary, i.e.  \f$ \phi_{i,j} =
 *  \phi_{\mathrm{const}} \f$. The Neumann boundary condition can be
 *  defined as first order discretation 
 *  \f[ \frac{\phi_{i+1}-\phi_{i}}{h} = N_{\mathrm{const}} \f]
 *  or second order discretation
 *  \f[ \frac{-\phi_{i+2}+4\phi_{i+1}-3\phi_{i}}{2h} = N_{\mathrm{const}} \f]
 *  selected by the user.
 *
 *  The plasma problems are described by using nonlinear models for
 *  screening charges in the plasma. For positive ion extraction for
 *  example, the screening charge is an electron population at the
 *  plasma potential \f$ \phi_p \f$ with a thermal energy distribution
 *  with temperature \f$ T_e \f$. The screening charge is therefore
 *  \f[ \rho_e = \rho_{e0} \exp \left( \frac{\phi-\phi_p}{kT_e/e} \right), \f]
 *  where electron charge density at plasma potential \f$ (\rho_{e0}) \f$ 
 *  is the same as the total positive beam space charge density for enabling
 *  plasma neutrality.
 */
class EpotSolver {

protected:

    Geometry           &_geom;             /*!< \brief Geometry reference. */

    plasma_mode_e       _plasma;           /*!< \brief Plasma simulation mode. */

    double              _rhoe;             /*!< \brief Electron charge density (C/m3), < 0. */
    double              _Te;               /*!< \brief Electron thermal energy, > 0. */
    double              _Up;               /*!< \brief Plasma potential, > 0. */

    std::vector<double> _rhoi;             /*!< \brief Charge density for positive ions, 
					    *   first fast protons, then thermal ions */
    std::vector<double> _Ei;               /*!< \brief Energy for positive ions,
					    *   first fast protons, then thermal ions */

    double              _force_pot;        /*!< \brief Potential to be forced. */
    CallbackFunctorB_V *_force_pot_func;   /*!< \brief Force region potential function. */
    CallbackFunctorD_V *_force_pot_func2;  /*!< \brief Force region potential function. */
    CallbackFunctorB_V *_init_plasma_func; /*!< \brief Initial plasma region function. */
    CallbackFunctorB_V *_plasma_calc_func; /*!< \brief Definition of plasma calculation region. */

    double              _plA;              /*!< \brief Plasma parameter.
				            *   For positive ion extraction: rho_th * h^2 / epsilon_0,
				            *   for negative ion extraction: rho_f * h^2 / epsilon_0 
					    *   for shield1: 1/Tm */
    double              _plB;              /*!< \brief Plasma parameter.
					    *   For positive ion extraction: 1/Te,
					    *   for negative ion extraction: E_f,i 
					    *   for shield1: phi_m/Tm */
    double              _plC;              /*!< \brief Plasma parameter for positive ion extraction. 
					    *    Up/Te */

    std::vector<double> _plD;              /*!< \brief Plasma parameter for negative ion extraction.
					    *   rho_th,i * h^2 / epsilon_0 */
    std::vector<double> _plE;              /*!< \brief Plasma parameter for negative ion extraction. 
					    *    1/Ti */
    

    /*! \brief Return non-linear right-hand-side and it's derivative
        for vacuum node in shielding model plasma.
     *
     *  Returns the non-linear part of the right-hand-side and it's
     *  derivative.
     */
    void shield_newton( double &rhs, double &drhs, double epot ) const;

    /*! \brief Return non-linear right-hand-side and it's derivative
        for vacuum node in positive ion plasma.
     *
     *  Returns the non-linear part of the right-hand-side and it's
     *  derivative resulting from the compensating space charge in
     *  positive ion extraction case.
     */
    void pexp_newton( double &rhs, double &drhs, double epot ) const;

    /*! \brief Return non-linear right-hand-side and it's derivative
        for vacuum node in negative ion plasma.
     *
     *  Returns the non-linear part of the right-hand-side and it's
     *  derivative resulting from the compensating space charge in
     *  negative ion extraction case.
     */
    void nsimp_newton( double &rhs, double &drhs, double epot ) const;

    /*! \brief Return bitmask indicating to which boundaries the node belongs to
     *
     *  1D version. The node i is checked. The return value has bits set
     *  accoding to, which boundaries the node belongs to. The lowest
     *  bit is xmin, second of xmax, third is ymin, fourth is ymax,
     *  fifth is zmin, sixth is zmax. If the node isn't a boundary
     *  node, a 0 is returned.
     */
    uint8_t boundary_index( uint32_t i ) const;

    /*! \brief Return bitmask indicating to which boundaries the node belongs to
     *
     *  2D version. The node (i,j) is checked. The return value has bits set
     *  accoding to, which boundaries the node belongs to. The lowest
     *  bit is xmin, second of xmax, third is ymin, fourth is ymax,
     *  fifth is zmin, sixth is zmax. If the node isn't a boundary
     *  node, a 0 is returned.
     */
    uint8_t boundary_index( uint32_t i, uint32_t j ) const;

    /*! \brief Return bitmask indicating to which boundaries the node belongs to
     *
     *  3D version. The node (i,j,k) is checked. The return value has bits set
     *  accoding to, which boundaries the node belongs to. The lowest
     *  bit is xmin, second of xmax, third is ymin, fourth is ymax,
     *  fifth is zmin, sixth is zmax. If the node isn't a boundary
     *  node, a 0 is returned.
     */
    uint8_t boundary_index( uint32_t i, uint32_t j, uint32_t k ) const;

    /*! \brief Return bitmask indicating to which boundaries the node belongs to
     *
     *  Version applicable to all dimensions. The node (i,j,k) is
     *  checked. The return value has bits set accoding to, which
     *  boundaries the node belongs to. The lowest bit is xmin, second
     *  of xmax, third is ymin, fourth is ymax, fifth is zmin, sixth
     *  is zmax. If the node isn't a boundary node, a 0 is returned.
     */
    uint8_t boundary_index_general( uint32_t i, uint32_t j, uint32_t k ) const;

    /*! \brief Do preprocessing action before solving.
     *
     *  Does precalculation of plasma paramters. Sets near solid
     *  neumann points as neumann and saves distance to a vector. Mark
     *  fixed and initial plasma nodes and set potentials.
     */
    void preprocess( MeshScalarField &epot );

    /*! \brief Do postprocessing action after solving.
     *
     *  Return near solid neumann points as near solid. Restore
     *  distance data. Remove fixed node tags.
     */
    void postprocess( void );

    /*! \brief Reset solver/problem settings.
     */
    virtual void reset_problem( void ) = 0;
    
    MeshScalarField *evaluate_scharge( const ScalarField &__scharge ) const;

    /*! \brief Solve problem with given mesh based space charge.
     */
    virtual void subsolve( MeshScalarField &epot, const MeshScalarField &scharge ) = 0;

public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Constructor for solver from \a geom.
     */
    EpotSolver( Geometry &geom );

    /*! \brief Constructor for solver from \a geom. Parameters from \a
     *  epsolver are copied to new solver.
     */
    EpotSolver( const EpotSolver &epsolver, Geometry &geom );

    /*! \brief Construct from file.
     */
    EpotSolver( Geometry &geom, std::istream &s );

    /*! \brief Destructor.
     */
    virtual ~EpotSolver() {}

/* ************************************** *
 * Problem constructing and solving       *
 * ************************************** */

    /*! \brief Copy parameters from solver \a epsolver.
     */
    void set_parameters( const EpotSolver &epsolver );

    /*! \brief Define forced potential volume.
     *
     *  Vacuum volume inside the volume defined by function \a
     *  force_pot_func will be forced to potential \a
     *  force_pot_func. This function is designed to be used with
     *  negative ion plasma extraction to stabilize plasma close
     *  non-physical boundaries.
     *
     *  This function may be removed in future.
     */
    void set_forced_potential_volume( double force_pot, 
				      CallbackFunctorB_V *force_pot_func );

    /*! \brief Define forced potential volume.
     *
     *  The \a force_pot_func functor will be called for every vacuum
     *  grid node at preprocessing stage of solve. If it returns a
     *  finite value, the node will be forced to the potential defined
     *  by the returned value. Otherwise the node will be treated as a
     *  regular vacuum node.
     */
    void set_forced_potential_volume( CallbackFunctorD_V *force_pot_func );

    /*! \brief Define initial plasma to the problem.
     *
     *  Initial plasma volume is defined in the area given by callback
     *  functor \a init_plasma_func.
     */
    void set_initial_plasma( double Up, 
			     CallbackFunctorB_V *init_plasma_func );

    /*! \brief Define plasma calculation region.
     *
     *  By default the plasma compensation function takes place in the
     *  whole geometry. With this function one can define a functor \a
     *  plasma_calc_func, which is called for all calculation nodes
     *  during Poisson solve. The function should return true if
     *  plasma calculation is enabled at that point and false if
     *  not.
     *
     *  Implemented only in matrix solvers.
     */
    void set_plasma_calc_region( CallbackFunctorB_V *plasma_calc_func );

    /*! \brief Enable plasma model for positive ion extraction problem.
     *
     *  Enable plasma model with background electron charge density of
     *  \a rhoe and electron temperature \a Te. The plasma potential is set to
     *  \a Up.
     *
     *  The boundary condition for the plasma should be
     *  BOUND_NEUMANN with gradient of zero V/m.
     */
    void set_pexp_plasma( double rhoe, double Te, double Up );

    /*! \brief Define initial plasma boundary location to negative ion
     *  extraction problem.
     *
     *  Initial plasma volume is defined in the area given by \a plasma_func.
     */
    void set_nsimp_initial_plasma( CallbackFunctorB_V *init_plasma_func );

    /*! \brief Enable plasma model for negative ion extraction problem.
     *
     *  The positive (analytic) space charges for the negative ion
     *  plasma extraction are set using this function. The positive
     *  ions consist of fast (directed) protons and any number of
     *  thermal positive ions trapped at the plasma boundary in the
     *  zero potential.
     *
     *  The parameters set are \a rhop, the space charge density of
     *  protons and \a Ep, the energy of protons at zero
     *  potential. Vectors \a rhoi and \a Ei are used to set the space
     *  charge densities and thermal energies of the trapped ions.
     *
     *  For a typical negative ion extraction case the compensating
     *  charges are positive (rhoi>0). This function allows also
     *  negative compensating species to be simulated with rhoi<0 and Ei<0.
     *
     *  The boundary condition for the plasma should be
     *  BOUND_DIRICHLET with zero volts.
     */
    void set_nsimp_plasma( double rhop, double Ep, 
			   std::vector<double> rhoi, std::vector<double> Ei );

    /*! \brief Enable shielding plasma model.
     *
     *  The Tm is the temperature of the compensating particles, Um is
     *  the potential on the surface where shielding factor is 1/2. If
     *  Tm is negative, the compensation is suitable for positive ion
     *  extraction and with positive Tm it is suitable for negative
     *  ion extraction.
     *
     *  The shielding function
     *  \f[ S=\frac{1}{2} \cdot \left[ 1 + \tanh\left( \frac{\phi-\phi_M}{T_M} \right)\right]
     *  multiplies the charge density from ray tracing to produce the 
     *  compensating effect of the plasma.
     */
    void set_shield_plasma( double Tm, double Um );

    /*! \brief Solve the problem.
     *
     *  The \a epot field is used as an initial guess for the
     *  solver. The space charge density field \a scharge is added to
     *  the problem vector before solving. The solution is returned in
     *  \a epot.
     */
    void solve( MeshScalarField &epot, const ScalarField &scharge );

/* ************************************** *
 * Solver interface                       *
 * ************************************** */

    /*! \brief Return true if problem is linear.
     */
    bool linear( void ) const;

/* ************************************** *
 * Misc                                   *
 * ************************************** */

    /*! \brief Get pointer to geometry.
     */
    const Geometry &geometry( void ) const;

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const = 0;

    /*! \brief Saves problem data to stream.
     */
    virtual void save( std::ostream &s ) const = 0;
};


#endif




















