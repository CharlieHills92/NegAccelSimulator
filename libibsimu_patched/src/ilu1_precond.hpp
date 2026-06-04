/*! \file ilu1_precond.hpp
 *  \brief ILU1 preconditioner for sparse matrices
 */

/* Copyright (c) 2012 Taneli Kalvas. All rights reserved.
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

#ifndef ILU1_PRECOND_HPP
#define ILU1_PRECOND_HPP 1


#include "matrix.hpp"
#include "crowmatrix.hpp"
#include "precond.hpp"


/*! \brief First level fill-in incomplete LU preconditioner class.
 */
class ILU1_Precond : public Precond {

    CRowMatrix *_LU;

    ILU1_Precond( const ILU1_Precond &pc ) {}
    const ILU1_Precond &operator=( const ILU1_Precond &pc ) { return( *this ); }

    static void add_element( int *ptr, int **col, int n, int &nz, int &asize, int c );

public:

    /*! \brief Constructor for an ILU1 preconditioner.
     */
    ILU1_Precond();

    /*! \brief Destructor.
     */
    ~ILU1_Precond();

    /*! \brief Get a new copy of preconditioner.
     *
     *  Does not copy matrix.
     */
    ILU1_Precond *copy( void ) const { return( new ILU1_Precond( *this ) ); }

    /*! \brief Prepare preconditioner for matrices with non-zero
     *  pattern equal to \a A.
     *
     *  Assumes ascending order matrix. Makes LU matrix a copy of \a A
     *  with elements where level is 1.
     */
    void prepare( const CRowMatrix &A );

    /*! \brief Construct preconditioner for matrix \a A.
     *
     *  Assumes ascending order matrix.
     */
    void construct( const CRowMatrix &A );

    /*! \brief Clear preconditioner.
     *
     *  Clears preconditioner. Both prepare() and construct()
     *  functions have to be called after this.
     */
    void clear( void );

    /*! \brief Return false if prepare is needed.
     *
     *  Returns true if prepare is not needed and false if it is.
     */
    bool is_prepared( void ) const;

    /*! \brief Return LU matrix.
     */
    const CRowMatrix *get_matrix( void ) const;

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;

    /*! \brief Return string indicating type of preconditioner.
     */
    std::string typestring( void ) const;
    
    /*! \brief Solve \a M* \a x = \a b and return \a x.
     */
    void solve( Vector &x, const Vector &b ) const;
};


#endif
