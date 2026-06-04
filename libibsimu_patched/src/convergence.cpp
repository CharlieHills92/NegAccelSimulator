/*! \file convergence.cpp
 *  \brief Vlasov system convergence tester
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


#include "convergence.hpp"
#include "ibsimu.hpp"


Convergence::ConvergenceField::ConvergenceField()
    : _field_old(0), _field(0)
{
    
}


Convergence::ConvergenceField::~ConvergenceField()
{
    if( _field_old )
	delete( _field_old ); 
}


void Convergence::ConvergenceField::evaluate_iteration( void )
{
    // Calculate convergence estimate for epot
    if( _field ) {
	if( !_field_old ) {
	    // First round 
	    _field_diff_max.push_back( 0.0 );
	    _field_diff_norm2.push_back( 0.0 );
	    _field_old = new MeshScalarField( *_field );
	} else {
	    // Evaluate error estimate
	    double max = 0.0;
	    double norm2 = 0.0;
	    int nc = _field->nodecount();
	    for( int a = 0; a < nc; a++ ) {
		double diff = fabs( (*_field)(a) - (*_field_old)(a) );
		if( diff > max )
		    max = diff;
		norm2 += diff*diff;
	    }
	    _field_diff_max.push_back( max );
	    _field_diff_norm2.push_back( norm2 );
	    delete _field_old;
	    _field_old = new MeshScalarField( *_field );
	}
    }
}


void Convergence::ConvergenceField::clear( void )
{
    _field = NULL;
    if( _field_old )
	delete( _field_old ); 
    _field_old = NULL;
    _field_diff_max.clear();
    _field_diff_norm2.clear();
}


/* **************************************************************************** */


Convergence::ConvergenceEmittance::ConvergenceEmittance()
    : _emit(0)
{

}


Convergence::ConvergenceEmittance::ConvergenceEmittance( const Emittance *emit )
    : _emit(emit)
{

}


Convergence::ConvergenceEmittance::~ConvergenceEmittance()
{

}


void Convergence::ConvergenceEmittance::evaluate_iteration( void )
{
    if( _emit ) {
	_emit_hist.push_back( *_emit );
    }
}


void Convergence::ConvergenceEmittance::clear( void )
{
    _emit_hist.clear();
}


/* **************************************************************************** */


Convergence::Convergence()
    : _iter(0)
{
}


Convergence::~Convergence()
{
}


void Convergence::evaluate_iteration( void )
{
    ibsimu.message( 1 ) <<  "Iteration round " << _iter << "\n";
    ibsimu.inc_indent();

    _epot.evaluate_iteration();
    _scharge.evaluate_iteration();
    for( uint32_t a = 0; a < _emit.size(); a++ )
	_emit[a].evaluate_iteration();

    if( _epot._field )
	ibsimu.message( 1 ) << "Epot error = " << _epot._field_diff_max.back() << " V (max)\n"; 
    if( _scharge._field )
	ibsimu.message( 1 ) << "Space charge error = " << _scharge._field_diff_max.back() << " C/m3 (max)\n";
    for( uint32_t a = 0; a < _emit.size(); a++ ) {
	uint32_t n = _emit[a]._emit_hist.size();
	if( n >=  2 )
	    ibsimu.message( 1 ) << "Emittance " << a << " error = "
				<< fabs(_emit[a]._emit_hist[n-1].epsilon()-_emit[a]._emit_hist[n-2].epsilon())
				<< " mm mrad (max)\n";
			 
	else
	    ibsimu.message( 1 ) << "Emittance " << a << " error = " << 0.0 << " mm mrad (max)\n";
    }

    // Increase iteration counter
    _iter++;

    ibsimu.dec_indent();
}


void Convergence::print_history( std::ostream &os ) const
{
    os << std::setw(6) << "# iter" << " ";
    if( _epot._field )
	os << std::setw(14) << "epot diff" << " ";
    if( _scharge._field )
	os << std::setw(14) << "scharge diff" << " ";
    for( uint32_t a = 0; a < _emit.size(); a++ ) {
	std::string s = "emit" + to_string(a) + " epsilon";
	os << std::setw(14) << s << " ";
	os << std::setw(14) << "alpha" << " ";
	os << std::setw(14) << "beta" << " ";
    }
    os << "\n";
	    
    for( int a = 1; a < _iter; a++ ) {
	os << std::setw(6) << a << " ";
	if( _epot._field )
	    os << std::setw(14) << _epot._field_diff_max[a] << " ";
	if( _scharge._field )
	    os << std::setw(14) << _scharge._field_diff_max[a] << " ";
	for( uint32_t b = 0; b < _emit.size(); b++ ) {
	    os << std::setw(14) << _emit[b]._emit_hist[a].epsilon() << " "
	       << std::setw(14) << _emit[b]._emit_hist[a].alpha() << " "
	       << std::setw(14) << _emit[b]._emit_hist[a].beta() << " ";
	}
	os << "\n";
    }
}


void Convergence::add_epot( const MeshScalarField &epot )
{
    _epot._field = &epot;
}


void Convergence::add_scharge( const MeshScalarField &scharge )
{
    _scharge._field = &scharge;
}


void Convergence::add_emittance( uint32_t i, const Emittance &emit )
{
    if( i < _emit.size() )
	_emit[i]._emit = &emit;
    else if( i == _emit.size() )
	_emit.push_back( &emit );
    else
	throw( Error( ERROR_LOCATION, "invalid index i = " + to_string(i) ) );
}


void Convergence::clear( void )
{
    _epot.clear();
    _scharge.clear();
    _emit.clear();
    _iter = 0;
}
