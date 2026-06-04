/*! \file hbio.cpp
 *  \brief Harwell Boeing sparse matrix file handling
 */

/* Copyright (c) 2005-2009,2012,2022 Taneli Kalvas. All rights reserved.
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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cmath>
#include "hbio.hpp"


/*
 * int  totcrd;       // Total number of data lines
 * int  ptrcrd;       // Number of pointer data lines
 * int  indcrd;       // Number of column index data lines
 * int  valcrd;       // Number of value data lines
 * int  rhscrd;       // Number of data lines for right hand side vectors
 * 
 * char mxtype[4];    // Matrix type, 0: R (real), C (complex), P (pattern)
 *                    //              1: U (unsymmetric), S (symmetric), H (Hermitian), 
 *                    //                 Z (skew symmetric), R (rectangular)
 *                    //              2: A (assembled), E (finite element matrix)
 * int  n;            // Number of rows
 * int  m;            // Number of columns
 * int  nz;           // Number of nonzero elements
 * int  nfelm;        // Number of finite element matrices, zero for assembled
 * 
 * int  ptrcount;     // Number of pointers per line
 * int  ptrlen;       // Length of pointer in chars
 * int  indcount;     // Number of column indices per line
 * int  indlen;       // Length of column index in chars
 * int  valcount;     // Number of values per line
 * int  vallen;       // Length of value in chars
 * int  valacc;       // Accuracy of value in chars
 * int  rhscount;     // Number of right hand sides per line
 * int  rhslen;       // Length of right hand side in chars
 * int  rhsacc;       // Accuracy of right hand side in chars
 * 
 * char rhstype[4];   // Right hand side type 0: F (full vector), M (unassembled vector)
 *                    //                      1: G (one or more guess vectors supplied), 
 *                    //                         " " (no guesses)
 *                    //                      2: X (one or more exact solutions supplied), 
 *                    //                         " " (no solutions)
 * int  nrhs;         // Number of right hand sides
 * int  nrhsix;       // Number of row indices
 */


void HBIO::write( const std::string filename ) const
{
    int  totcrd, ptrcrd, indcrd, valcrd, rhscrd;
    int  n, m, nz, nfelm, nrhs, nrhsix;
    int  ptrcount, ptrlen, indcount, indlen, valcount, vallen, rhscount, rhslen;
    char mxtype[4], rhstype[4];

    std::fstream fstr( filename.c_str(), std::fstream::out );   
    if( !fstr.is_open() )
	throw( Error( ERROR_LOCATION, "Could not open file " + filename ) );

    // Make checks
    if( rhs.size() != 0 && mat.columns() != rhs.size() )
	throw( Error( ERROR_LOCATION, "Matrix dimension does not match vector dimension" ) );

    // Set file parameters
    n      = mat.rows();
    m      = mat.columns();
    nz     = mat.nz_elements();
    nfelm  = 0;
    nrhs   = 0;
    nrhsix = 0;
    strcpy( mxtype, "RUA" );
    strcpy( rhstype, "F  " );
    if( sol.size() != 0 ) {
	rhstype[2] = 'X';
	nrhs++;
    }
    if( rhs.size() != 0 )
	nrhs++;

    ptrlen   = (int)log10( nz+1 ) + 2;
    ptrcount = 80 / ptrlen;
    ptrcrd   = (int)ceil( (m+1) / (double)ptrcount );

    indlen   = (int)log10( n+1 ) + 2;
    indcount = 80 / indlen;
    indcrd   = (int)ceil( nz / (double)indcount );

    vallen   = valacc+8;
    valcount = 80 / vallen;
    valcrd   = (int)ceil( nz / (double)valcount );

    if( rhs.size() != 0 || sol.size() != 0 ) {
	rhslen   = rhsacc+8;
	rhscount = 80 / rhslen;
	rhscrd   = (int)ceil( nrhs*n / (double)rhscount );
    } else {
	rhslen   = 0;
	rhscount = 0;
	rhscrd   = 0;
    }

    totcrd = ptrcrd + indcrd + valcrd + rhscrd;

    // Write line 1
    fstr.unsetf( std::ios::adjustfield );
    fstr.setf( std::ios::left );
    fstr << std::setw(72) << title.substr(0,72);
    fstr << std::setw(8) << key.substr(0,8) << "\n";

    // Write line 2
    fstr.unsetf( std::ios::adjustfield );
    fstr.setf( std::ios::right );
    fstr << std::setw(14) << totcrd;
    fstr << std::setw(14) << ptrcrd;
    fstr << std::setw(14) << indcrd;
    fstr << std::setw(14) << valcrd;
    fstr << std::setw(14) << rhscrd << "\n";

    // Write line 3
    fstr << std::setw(3) << mxtype;
    fstr << std::setw(11) << " ";
    fstr << std::setw(14) << n;
    fstr << std::setw(14) << m;
    fstr << std::setw(14) << nz;
    fstr << std::setw(14) << nfelm << "\n";

    // Write line 4
    fstr.unsetf( std::ios::adjustfield );
    fstr.setf( std::ios::left );
    std::ostringstream ss1, ss2, ss3, ss4;
    ss1 << '(' << ptrcount << 'I' << ptrlen << ')';
    fstr << std::setw(16) << ss1.str();
    ss2 << '(' << indcount << 'I' << indlen << ')';
    fstr << std::setw(16) << ss2.str();
    ss3 << '(' << valcount << 'E' << vallen << '.' << valacc << ')';
    fstr << std::setw(20) << ss3.str();
    if( rhs.size() != 0 || sol.size() != 0 ) {
	ss4 << '(' << rhscount << 'E' << rhslen << '.' << rhsacc << ')';
	fstr << std::setw(20) << ss4.str();
    }
    fstr << "\n";

    if( rhscrd > 0 ) {
	// Write line 5
	fstr << std::setw(3) << rhstype;
	fstr << std::setw(11) << " ";
	fstr.unsetf( std::ios::adjustfield );
	fstr.setf( std::ios::right );
	fstr << std::setw(14) << nrhs;
	fstr << std::setw(14) << nrhsix << "\n";
    }

    fstr.unsetf( std::ios::adjustfield );
    fstr.setf( std::ios::right );

    // Write matrix pointers
    for( int k = 0; k < m+1; ) {
	for( int i = 0; i < ptrcount && k < m+1; i++ )
	    fstr << std::setw(ptrlen) << mat.ptr(k++) + 1;
	fstr << "\n";
    }

    // Write matrix indices
    for( int k = 0; k < nz; ) {
	for( int i = 0; i < indcount && k < nz; i++ )
	    fstr << std::setw(indlen) << mat.row(k++) + 1;
	fstr << "\n";
    }

    fstr.unsetf( std::ios::floatfield );
    fstr.setf( std::ios::scientific );
    fstr.precision( valacc );

    // Write matrix values
    for( int k = 0; k < nz; ) {
	for( int i = 0; i < valcount && k < nz; i++ )
	    fstr << std::setw(vallen) << mat.val(k++);
	fstr << "\n";
    }

    fstr.precision( rhsacc );

    // Write right hand side vectors
    if( rhs.size() != 0 ) {
	for( int k = 0; k < n; ) {
	    for( int i = 0; i < rhscount && k < n; i++ )
		fstr << std::setw(rhslen) << rhs(k++);
	    fstr << "\n";
	}    
    }
    if( sol.size() != 0 ) {
	for( int k = 0; k < n; ) {
	    for( int i = 0; i < rhscount && k < n; i++ )
		fstr << std::setw(rhslen) << sol(k++);
	    fstr << "\n";
	}    
    }

    fstr.close();
}


void HBIO::read( const std::string filename )
{
    const char *s;
    char *s2;
    int line = 0;
    size_t pos;
    std::string buf;
    std::string buf2;

    std::fstream fstr( filename.c_str(), std::fstream::in );
    if( !fstr.is_open() )
	throw( Error( ERROR_LOCATION, "Could not open file " + filename ) );

    // Read line 1
    getline( fstr, buf );
    line++;
    buf2 = buf.substr( 0, 72 );
    if( (pos = buf2.find_last_not_of( ' ' )) != std::string::npos )
	buf2 = buf2.substr( 0, pos+1 );
    title = buf2;
    try {
	buf2 = buf.substr( 72, 8 );
	if( (pos = buf2.find_last_not_of( ' ' )) != std::string::npos )
	    buf2 = buf2.substr( 0, pos+1 );
	key = buf2;
    } catch( const std::out_of_range & ) {
	key = "";
    }

    // Read line 2
    getline( fstr, buf );
    line++;
    int totcrd = atoi( buf.substr(  0, 14 ).c_str() );
    int ptrcrd = atoi( buf.substr( 14, 14 ).c_str() );
    int indcrd = atoi( buf.substr( 28, 14 ).c_str() );
    int valcrd = atoi( buf.substr( 42, 14 ).c_str() );
    int rhscrd = atoi( buf.substr( 56, 14 ).c_str() );
    if( totcrd != ptrcrd + indcrd + valcrd + rhscrd )
	throw( Error( ERROR_LOCATION, "Invconsistent line counts on line 2" ) );

    // Read line 3
    getline( fstr, buf );
    line++;
    char mxtype[4];
    strcpy( mxtype, buf.substr( 0, 3 ).c_str() );
    int n      = atoi( buf.substr( 14, 14 ).c_str() );
    int m      = atoi( buf.substr( 28, 14 ).c_str() );
    int nz     = atoi( buf.substr( 42, 14 ).c_str() );
    int nfelm  = atoi( buf.substr( 56, 14 ).c_str() );
    if( n <= 0 || m <= 0 )
	throw( Error( ERROR_LOCATION, "Invalid matrix dimensions on line 3" ) );

    // Read line 4
    getline( fstr, buf );
    line++;
    // Read ptrfmt
    buf2 = buf.substr(  0, 16 );
    s = buf2.c_str();
    while( *s == ' ' ) s++;
    if( *s == '(' ) s++;
    int ptrcount = strtol( s, &s2, 10 );
    if( s == s2 || (*s2 != 'I' && *s2 != 'i') ) 
	throw( Error( ERROR_LOCATION, "Integer format expected on line 4" ) );
    s = s2 + 1;
    int ptrlen = strtol( s, &s2, 10 );
    if( s == s2 )
	throw( Error( ERROR_LOCATION, "Integer format expected on line 4" ) );

    // Read indfmt
    buf2 = buf.substr( 16, 16 );
    s = buf2.c_str();
    while( *s == ' ' ) s++;
    if( *s == '(' ) s++;
    int indcount = strtol( s, &s2, 10 );
    if( s == s2 || (*s2 != 'I' && *s2 != 'i') )
	throw( Error( ERROR_LOCATION, "Integer format expected in line 4" ) );
    s = s2 + 1;
    int indlen = strtol( s, &s2, 10 );
    if( s == s2 )
	throw( Error( ERROR_LOCATION, "Integer format expected in line 4" ) );

    // Read valfmt
    buf2 = buf.substr( 32, 20 );
    s = buf2.c_str();
    while( *s == ' ' ) s++;
    if( *s == '(' ) s++;
    int valcount = strtol( s, &s2, 10 );
    if( *s2 == 'P' || *s2 == 'p' ) {
	// Skip scaling parameter
	s = s2 + 1;
	valcount = strtol( s, &s2, 10 );
    }
    if( s == s2 || (*s2 != 'E' && *s2 != 'e' && 
		    *s2 != 'G' && *s2 != 'g' && 
		    *s2 != 'D' && *s2 != 'd' && 
		    *s2 != 'F' && *s2 != 'f') )
	throw( Error( ERROR_LOCATION, "Real format expected in line 4" ) );
    s = s2 + 1;
    int vallen = strtol( s, &s2, 10 );
    if( s == s2 || *s2 != '.' )
	throw( Error( ERROR_LOCATION, "Real format expected in line 4" ) );
    s = s2 + 1;
    valacc = strtol( s, &s2, 10 );
    if( s == s2 )
	throw( Error( ERROR_LOCATION, "Real format expected in line 4" ) );

    // Read rhsfmt
    int rhscount;
    int rhslen;
    if( rhscrd > 0 ) {
	buf2 = buf.substr( 52, 20 );
	s = buf2.c_str();
	while( *s == ' ' ) s++;
	if( *s == '(' ) s++;
	rhscount = strtol( s, &s2, 10 );
	if( *s2 == 'P' || *s2 == 'p' ) {
	    // Skip scaling parameter 
	    rhscount = strtol( s2, &s2, 10 );
	}
	if( s == s2 || (*s2 != 'E' && *s2 != 'e' && 
			*s2 != 'G' && *s2 != 'g' && 
			*s2 != 'D' && *s2 != 'd' && 
			*s2 != 'F' && *s2 != 'f') )
	throw( Error( ERROR_LOCATION, "Real format expected in line 4" ) );
	s = s2 + 1;
	rhslen = strtol( s, &s2, 10 );
	if( s == s2 || *s2 != '.' )
	    throw( Error( ERROR_LOCATION, "Real format expected in line 4" ) );
	s = s2 + 1;
	rhsacc = strtol( s, &s2, 10 );
	if( s == s2 )
	    throw( Error( ERROR_LOCATION, "Real format expected in line 4" ) );
    } else {
	rhscount = 0;
	rhslen   = 0;
	rhsacc   = 0;
    }

    char rhstype[4];
    int nrhs;
    int nrhsix;
    if( rhscrd > 0 ) {
	// Read line 5
	getline( fstr, buf );
	line++;
	strcpy( rhstype, buf.substr( 0, 3 ).c_str() );
	nrhs   = atoi( buf.substr( 14, 14 ).c_str() );
	nrhsix = atoi( buf.substr( 28, 14 ).c_str() );
	if( nrhs <= 0 )
	    throw( Error( ERROR_LOCATION, "Right hand side data not expected" ) );
    } else {
	strcpy( rhstype, "F  " );
	nrhs   = 0;
	nrhsix = 0;
    }

    // Keep compilers happy, un-unsed variable
    if( nrhsix ) ;

    // Make checks
    if( strcmp( mxtype, "RUA" ) )
	throw( Error( ERROR_LOCATION, "Only RUA matrices are supported" ) );
    if( rhstype[0] != 'F' )
	throw( Error( ERROR_LOCATION, "Only full vectors are supported" ) );

    if( ptrcrd*ptrcount < m+1 )
	throw( Error( ERROR_LOCATION, "Too few pointer data lines" ) );
    if( (ptrcrd-1)*ptrcount >= m+1 )
	throw( Error( ERROR_LOCATION, "Too many pointer data lines" ) );

    if( indcrd*indcount < nz )
	throw( Error( ERROR_LOCATION, "Too few index data lines" ) );
    if( (indcrd-1)*indcount >= nz )
	throw( Error( ERROR_LOCATION, "Too many index data lines" ) );

    if( valcrd*valcount < nz )
	throw( Error( ERROR_LOCATION, "Too few value data lines" ) );
    if( (valcrd-1)*valcount >= nz )
	throw( Error( ERROR_LOCATION, "Too many value data lines" ) );

    if( rhscrd*rhscount < n*nrhs )
	throw( Error( ERROR_LOCATION, "Too few right hand side data lines" ) );
    if( rhscrd != 0 && (rhscrd-1)*rhscount >= n*nrhs )
	throw( Error( ERROR_LOCATION, "Too many right hand side data lines" ) );

    if( nfelm > 0 )
	throw( Error( ERROR_LOCATION, "Finite element matrices not expected" ) );


    // Initialize 
    mat.clear();
    int *m_ptr = new int[m+1];
    int *m_row = new int[nz];
    double *m_val = new double[nz];
    rhs.resize( 0 );
    sol.resize( 0 );

    // Read matrix pointers
    int rloc;
    int k = 0;
    for( int i = 0; i < ptrcrd; i++ ) {
	getline( fstr, buf );
	line++;
	rloc = 0;
	for( int j = 0; j < ptrcount && k < m+1; j++, k++ ) {
	    buf2 = buf.substr( rloc, ptrlen );
	    m_ptr[k] = strtol( buf2.c_str(), &s2, 10 ) - 1;
	    if( m_ptr[k] < 0 || (!(*s2 >= '0' && *s2 <= '9') && *s2 != ' ' && *s2 != '\0') )
		throw( Error( ERROR_LOCATION, "Unexpected data in matrix pointers" ) );
	    rloc += ptrlen;
	}
    }

    // Read matrix indices
    k = 0;
    for( int i = 0; i < indcrd; i++ ) {
	getline( fstr, buf );
	line++;
	rloc = 0;
	for( int j = 0; j < indcount && k < nz; j++, k++ ) {
	    buf2 = buf.substr( rloc, indlen );
	    m_row[k] = strtol( buf2.c_str(), &s2, 10 ) - 1;
	    if( m_row[k] < 0 || (!(*s2 >= '0' && *s2 <= '9') && *s2 != ' ' && *s2 != '\0') )
		throw( Error( ERROR_LOCATION, "Unexpected data in matrix indices" ) );
	    rloc += indlen;
	}
    }

    // Read matrix values
    k = 0;
    for( int i = 0; i < valcrd; i++ ) {
	getline( fstr, buf );
	line++;
	rloc = 0;
	for( int j = 0; j < valcount && k < nz; j++, k++ ) {
	    buf2 = buf.substr( rloc, vallen );
	    m_val[k] = strtod( s = buf2.c_str(), &s2 );
	    if( *s2 == 'D' || *s2 == 'd' ) {
		buf2[s2-s] = 'E';
		m_val[k] = strtod( buf2.c_str(), &s2 );
	    }
	    if( !(*s2 >= '0' && *s2 <= '9') && *s2 != ' ' && *s2 != '\0' )
		throw( Error( ERROR_LOCATION, "Unexpected data in matrix values" ) );
	    rloc += vallen;
	}
    }

    // Set matrix
    mat = CColMatrix( n, m, nz, nz, m_ptr, m_row, m_val );

    // Read right hand side vectors
    if( rhscrd > 0 ) {

	// Read right hand side 
	rhs.resize( n );
	for( k = 0; k < n; ) {
	    getline( fstr, buf );
	    line++;
	    rloc = 0;
	    for( int j = 0; j < rhscount && k < n; j++, k++ ) {
		buf2 = buf.substr( rloc, rhslen );
		rhs[k] = strtod( s = buf2.c_str(), &s2 );
		if( *s2 == 'D' || *s2 == 'd' ) {
		    buf2[s2-s] = 'E';
		    rhs[k] = strtod( buf2.c_str(), &s2 );
		}
		if( !(*s2 >= '0' && *s2 <= '9') && *s2 != ' ' && *s2 != '\0' )
		    throw( Error( ERROR_LOCATION, "Unexpected data in right hand side vector values" ) );
		rloc += rhslen;
	    }
	}

	if( rhstype[2] == 'X' ) {
	    // Read solution vector
	    sol.resize( n );
	    for( k = 0; k < n; ) {
		getline( fstr, buf );
		line++;
		rloc = 0;
		for( int j = 0; j < rhscount && k < n; j++, k++ ) {
		    buf2 = buf.substr( rloc, rhslen );
		    sol[k] = strtod( s = buf2.c_str(), &s2 );
		    if( *s2 == 'D' || *s2 == 'd' ) {
			buf2[s2-s] = 'E';
			sol[k] = strtod( buf2.c_str(), &s2 );
		    }
		    if( !(*s2 >= '0' && *s2 <= '9') && *s2 != ' ' && *s2 != '\0' )
			throw( Error( ERROR_LOCATION, "Unexpected data in right hand side vector values" ) );
		    rloc += rhslen;
		}
	    }
	}
    }

    fstr.close();
}


void HBIO::get_matrix( CColMatrix &mmat ) const 
{
    mmat = mat;
}


void HBIO::set_matrix( const CColMatrix &mmat )
{
    mat = mmat;
}


void HBIO::get_rhs_vector( Vector &rrhs ) const 
{
    rrhs = rhs;
}


void HBIO::set_rhs_vector( const Vector &rrhs )
{
    rhs = rrhs;
}


void HBIO::get_solution_vector( Vector &ssol ) const
{
    ssol = sol;
}


void HBIO::set_solution_vector( const Vector &ssol )
{
    sol = ssol;
}

