/*! \file hbio.hpp
 *  \brief Harwell Boeing sparse matrix file handling
 */

/* Copyright (c) 2005-2010,2012 Taneli Kalvas. All rights reserved.
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

#ifndef HBIO_HPP
#define HBIO_HPP 1


#include <cstdlib>
#include <iostream>
#include "error.hpp"
#include "mvector.hpp"
#include "ccolmatrix.hpp"


/*! \brief Harwell Boeing sparse matrix file format I/O class.
 *
 *  Class for reading and writing linear algebra problems (matrices
 *  and vectors) in the standard Harwell Boeing (HB) sparse matrix
 *  file format. The HBIO class supports a limited subset of HB file
 *  format features. The most basic features including title, key,
 *  selectable number accuracy, problem matrix, right hand side vector
 *  and the solution vector are supported. The HB format stores
 *  matrices in compressed column mode with fortran indexing (indices
 *  starting from 1), but the matrices are converted to 0 based
 *  indexing when read to memory.
 */
class HBIO {
    std::string filename;     //!< Filename
    std::string title;        //!< Title (max 72 chars)
    std::string key;          //!< Key (max 8 chars)
    int         valacc;       //!< Value accuracy
    int         rhsacc;       //!< Right hand side accuracy

    CColMatrix  mat;          //!< Matrix
    Vector      rhs;          //!< Right hand side vector
    Vector      sol;          //!< Solution vector

public:

    /*! \brief Constructor
     */
    HBIO();

    /*! \brief Destructor
     */
    ~HBIO();

    /*! \brief Write file \a filename
     */
    void write( const std::string filename ) const;

    /*! \brief Read file \a filename
     */
    void read( const std::string filename );

    /*! \brief Get file title (max 72 chars)
     */
    const std::string get_title( void ) const;

    /*! \brief Set file title (max 72 chars)
     */
    void set_title( const std::string ttitle );

    /*! \brief Get file key (max 8 chars)
     */
    const std::string get_key( void ) const;

    /*! \brief Set file key (max 8 chars)
     */
    void set_key( const std::string kkey );

    /*! \brief Get value accuracy (in chars)
     */
    int get_valacc( void ) const;

    /*! \brief Set value accuracy (in chars)
     */
    void set_valacc( int vvalacc );

    /*! \brief Get right hand side accuracy (in chars)
     */
    int get_rhsacc( void ) const;

    /*! \brief Set right hand side accuracy (in chars)
     */
    void set_rhsacc( int rrhsacc );

    /*! \brief Get matrix
     */
    void get_matrix( CColMatrix &mmat ) const;

    /*! \brief Set matrix
     */
    void set_matrix( const CColMatrix &mmat );

    /*! \brief Get right hand side vector
     */
    void get_rhs_vector( Vector &rrhs ) const;

    /*! \brief Set right hand side vector
     */
    void set_rhs_vector( const Vector &rrhs );

    /*! \brief Get solution vector
     */
    void get_solution_vector( Vector &ssol ) const;

    /*! \brief Set solution vector
     */
    void set_solution_vector( const Vector &ssol );
};


inline HBIO::HBIO()
{
    valacc = 8;
    rhsacc = 8;
}


inline HBIO::~HBIO()
{
}


inline const std::string HBIO::get_title( void ) const
{
    return( title );
}


inline void HBIO::set_title( const std::string ttitle )
{
    title = ttitle;
}


inline const std::string HBIO::get_key( void ) const
{
    return( key );
}


inline void HBIO::set_key( const std::string kkey )
{
    key = kkey;
}


inline int HBIO::get_valacc( void ) const
{
    return( valacc );
}


inline void HBIO::set_valacc( int vvalacc )
{
    valacc = vvalacc;
}


inline int HBIO::get_rhsacc( void ) const
{
    return( rhsacc );
}


inline void HBIO::set_rhsacc( int rrhsacc )
{
    rhsacc = rrhsacc;
}


#endif
