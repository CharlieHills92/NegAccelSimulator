/*! \file convergence.hpp
 *  \brief Vlasov system convergence follower
 */

/* Copyright (c) 2011 Taneli Kalvas. All rights reserved.
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

#ifndef CONVERGENCE_HPP
#define CONVERGENCE_HPP 1


#include <iostream>
#include <stdint.h>
#include <vector>
#include "meshscalarfield.hpp"
#include "trajectorydiagnostics.hpp"


/*! \brief Vlasov system convergence follower.
 */
class Convergence {

    int32_t                 _iter;                 /*!< \brief Iteration counter. */

    struct ConvergenceField {
	MeshScalarField        *_field_old;        /*!< \brief Epot of last round. */
	const MeshScalarField  *_field;            /*!< \brief Epot pointer. */
	std::vector<double>     _field_diff_max;   /*!< \brief Convergence history of epot. */
	std::vector<double>     _field_diff_norm2; /*!< \brief Convergence history of epot. */

	ConvergenceField();
	~ConvergenceField();

	void evaluate_iteration( void );
	void clear( void );
    };

    ConvergenceField        _epot;                 /*!< \brief Convergence of epot. */
    ConvergenceField        _scharge;              /*!< \brief Convergence of scharge. */

    struct ConvergenceEmittance {
	std::vector<Emittance>  _emit_hist;        /*!< \brief Convergence history of emittance. */
	const Emittance        *_emit;             /*!< \brief Emittance pointer. */

	ConvergenceEmittance();
	ConvergenceEmittance( const Emittance *emit );
	~ConvergenceEmittance();

	void evaluate_iteration( void );
	void clear( void );
    };

    std::vector<ConvergenceEmittance> _emit;
   
public:

    /*! \brief Constructor for convergence class.
     */
    Convergence();

    /*! \brief Destructor for convergence class.
     */
    ~Convergence();

    /*! \brief Evaluate convergence of iteration round.
     */
    void evaluate_iteration( void );

    /*! \brief Print the history of convergence to stream \a os.
     */
    void print_history( std::ostream &os ) const;

    /*! \brief Add a reference to electric potential to be followed.
     */
    void add_epot( const MeshScalarField &epot );

    /*! \brief Add a reference to space charge density to be followed.
     */
    void add_scharge( const MeshScalarField &scharge );

    /*! \brief Add a reference to emittance to be followed.
     *
     *  Multiple emittances can be followed. Emittances must be
     *  referenced with running numbers \a i starting from
     *  0. Overwriting a definition is possible.
     */
    void add_emittance( uint32_t i, const Emittance &emit );

    /*! \brief Clear.
     *
     *  Clears references to followed objects, clears history and
     *  clears iteration counter.
     */
    void clear( void );
};


#endif

