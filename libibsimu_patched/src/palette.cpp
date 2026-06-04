/*! \file palette.cpp
 *  \brief Palette for colormaps
 */

/* Copyright (c) 2005-2009,2011-2013 Taneli Kalvas. All rights reserved.
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
#include <limits>
#include <math.h>
#include <algorithm>
#include "palette.hpp"



Palette::Entry::Entry( const Vec3D &color, double val )
{
    _color = color;
    _val = val;
}


bool Palette::Entry::operator<( const Entry &e ) const
{
    return( _val < e._val );
}






Palette::Palette()
    : _steps(0)
{
    _entries.push_back( Entry( Vec3D( 1.0, 1.0, 1.0 ), 0.0 ) );
    _entries.push_back( Entry( Vec3D( 0.0, 0.0, 0.0 ), 1.0 ) );
}


Palette::Palette( const std::vector<Entry> &entries )
    : _steps(0)
{
    // Search minimum and maximum
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();
    for( size_t a = 0; a < entries.size(); a++ ) {
	if( entries[a]._val < min )
	    min = entries[a]._val;
	if( entries[a]._val > max )
	    max = entries[a]._val;
    }
    
    // Calculate offset and coef
    double offset = -min;
    double coef = 1.0/(max-min);
    if( max-min == 0.0 )
	coef = 1.0;

    // Clear old palette
    _entries.clear();

    // Add new entries
    _entries.reserve( entries.size() );
    for( size_t a = 0; a < entries.size(); a++ ) {
	_entries.push_back( Entry( entries[a]._color, coef*(offset+entries[a]._val) ) );
    }

    // Sort entries
    sort( _entries.begin(), _entries.end() );
}


void Palette::clear( void )
{
    // Clear old palette
    _entries.clear();
}


void Palette::push_back( const Vec3D &color, double val )
{
    // Add new entry
    _entries.push_back( Entry( color, val ) );

    // Sort entries
    sort( _entries.begin(), _entries.end() );
}


int Palette::get_steps( void ) const
{
    return( _steps );
}


void Palette::set_steps( int steps )
{
    _steps = steps;
}


void Palette::normalize( void )
{
    // Search minimum and maximum
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();
    for( size_t a = 0; a < _entries.size(); a++ ) {
	if( _entries[a]._val < min )
	    min = _entries[a]._val;
	if( _entries[a]._val > max )
	    max = _entries[a]._val;
    }
    
    // Calculate offset and coef
    double offset = -min;
    double coef = 1.0/(max-min);
    if( max-min == 0.0 )
	coef = 1.0;

    // Renormalize entries
    for( size_t a = 0; a < _entries.size(); a++ )
	_entries[a]._val = coef*(offset+_entries[a]._val);
}


Vec3D Palette::operator()( double x ) const
{
    // If undefined
    if( _entries.size() == 0 )
	return( Vec3D( 0, 0, 0 ) );
    else if( _entries.size() == 1 )
	return( _entries[0]._color );

    // If out of limits
    if( x <= 0.0 )
	return( _entries[0]._color );
    else if( x >= 1.0 )
	return( _entries[_entries.size()-1]._color );

    // If stepped palette
    if( _steps > 1 ) {
	x *= _steps;
	x = floor( x );
	x = x / (_steps-1);
	if( x > 1.0 )
	    x = 1.0;
    }

    // Search correct index
    size_t a;
    for( a = 1; a < _entries.size(); a++ ) {
	if( x < _entries[a]._val )
	    break;
    }
    if( a == _entries.size() ) 
	a = _entries.size()-1;
    
    // Interpolate
    double t = (x-_entries[a-1]._val) / (_entries[a]._val-_entries[a-1]._val);
    Vec3D c = _entries[a-1]._color + 
	t*(_entries[a]._color-_entries[a-1]._color);

    return( c );
}


void Palette::debug_print( std::ostream &os ) const
{
    os << "**Palette\n";
    os << "steps = " << _steps << "\n";
    os << "size = " << _entries.size() << "\n";
    for( size_t a = 0; a < _entries.size(); a++ ) {
	os << "entries[" << a << "] = " 
	   << _entries[a]._val << " "
	   << _entries[a]._color[0] << " "
	   << _entries[a]._color[1] << " "
	   << _entries[a]._color[2] << " "
	   << _entries[a]._color[3] << "\n";
    }
}
