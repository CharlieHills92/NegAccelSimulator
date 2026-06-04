/*! \file mydxfentities.cpp
 *  \brief DXF entities
 */

/* Copyright (c) 2010-2012,2014 Taneli Kalvas. All rights reserved.
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
#include <math.h>
#include <limits>
#include "mydxfentities.hpp"
#include "mydxfline.hpp"
#include "mydxflwpolyline.hpp"
#include "mydxfspline.hpp"
#include "mydxfarc.hpp"
#include "mydxfcircle.hpp"
#include "mydxfmtext.hpp"
#include "mydxfinsert.hpp"
#include "polysolver.hpp"
#include "error.hpp"
#include "ibsimu.hpp"




//#define MYDXF_DEBUG_BBOX
//#define MYDXF_DEBUG 1
//#define DEBUG_INSIDE_LOOP 1




/* ************************************************************************** *
 * DXFEntity                                                                  *
 * ************************************************************************** */


MyDXFEntity::MyDXFEntity()
{

}


void MyDXFEntity::bbox_ppoint( Vec3D &min, Vec3D &max, const Vec3D &p )
{
    for( int a = 0; a < 3; a++ ) {
	if( p[a] < min[a] )
	    min[a] = p[a];
	if( p[a] > max[a] )
	    max[a] = p[a];
    }
}


void MyDXFEntity::debug_print_base( std::ostream &os ) const
{
    MyDXFFile::debug_print_format( os, "handle", _handle );
    MyDXFFile::debug_print_format( os, "layer", _layer );
}


std::ostream &operator<<( std::ostream &os, const MyDXFEntity &ent )
{
    ent.debug_print( os );
    ent.debug_print_base( os );
    return( os );
}


void MyDXFEntity::write_common( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 5, _handle.c_str() );
    dxf->write_group( 8, _layer.c_str() );
}


void MyDXFEntity::process_group( class MyDXFFile *dxf )
{
    if( dxf->group_get_code() == 5 ) {
	_handle = dxf->group_get_string();
    } else if( dxf->group_get_code() == 8 ) {
	_layer = dxf->group_get_string();
    }
}


/* ************************************************************************** *
 * DXFEntities                                                                *
 * ************************************************************************** */


MyDXFEntities::MyDXFEntities( class MyDXFFile *dxf )
    : _dxf(dxf)
{

}


MyDXFEntities::MyDXFEntities( class MyDXFFile *dxf, MyDXFEntities *ent, MyDXFEntitySelection *sel )
    : _dxf(dxf)
{
    for( size_t a = 0; a < sel->size(); a++ ) {
	MyDXFEntity *e = ent->get_entity( (*sel)(a) );
	_entities.push_back( e->copy() );
    }
}


MyDXFEntities::MyDXFEntities( class MyDXFFile *dxf, bool reading_blocks )
    : _dxf(dxf)
{
#ifdef MYDXF_DEBUG
    std::cout << "Reading ENTITIES\n";
#endif

    // Title of first entity already read if reading BLOCKS section
    if( !reading_blocks )
	dxf->read_group();

    while( dxf->group_get_code() != -1 ) {

	if( dxf->group_get_code() != 0 ) {
	    dxf->read_group();
	    continue; // Skip unknown input
	}
	else if( dxf->group_get_string() == "ENDSEC" ||
		 dxf->group_get_string() == "ENDBLK" ) {
	    break; // Done with entities
	}

	// Read entity type
	if( dxf->group_get_string() == "LINE" ) {
	    _entities.push_back( new MyDXFLine( dxf ) );
	} else if( dxf->group_get_string() == "CIRCLE" ) {
	    _entities.push_back( new MyDXFCircle( dxf ) );
	} else if( dxf->group_get_string() == "MTEXT" ) {
	    _entities.push_back( new MyDXFMText( dxf ) );
	} else if( dxf->group_get_string() == "INSERT" ) {
	    _entities.push_back( new MyDXFInsert( dxf ) );
	} else if( dxf->group_get_string() == "ARC" ) {
	    _entities.push_back( new MyDXFArc( dxf ) );
	} else if( dxf->group_get_string() == "LWPOLYLINE" ) {
	    _entities.push_back( new MyDXFLWPolyline( dxf ) );
	} else if( dxf->group_get_string() == "SPLINE" ) {
	    _entities.push_back( new MyDXFSpline( dxf ) );
	} else {
	    if( dxf->wlevel() >= 2 )
		std::cout << "Skipping unknown entity \'" << dxf->group_get_string() << "\'\n";
	    dxf->read_group();
	}
    }

#ifdef MYDXF_DEBUG
    std::cout << "Done with ENTITIES\n";
#endif
}


MyDXFEntities::~MyDXFEntities()
{
    // Free data
    for( size_t a = 0; a < _entities.size(); a++ )
	delete _entities[a];
}


void MyDXFEntities::write( class MyDXFFile *dxf, std::ofstream &ostr )
{
    dxf->write_group( 0, "SECTION" );
    dxf->write_group( 2, "ENTITIES" );

    write_entities( dxf, ostr );

    dxf->write_group( 0, "ENDSEC" );
}


void MyDXFEntities::write_entities( class MyDXFFile *dxf, std::ofstream &ostr )
{
    for( size_t a = 0; a < _entities.size(); a++ )
	_entities[a]->write( dxf, ostr );
}


bool MyDXFEntities::inside_loop( MyDXFEntitySelection *selection, double x, double y, double eps )
{
#ifdef DEBUG_INSIDE_LOOP
    std::cout << "inside_loop( x = " << x << ", y = " << y << ", eps = " << eps << " )\n";
#endif

    // Number of needed perturbations hard to estimate, depends on
    // entities. Guesstimate size+2.
    for( uint32_t b = 0; b < selection->size()+2; b++ ) {

#ifdef DEBUG_INSIDE_LOOP
	std::cout << "b = " << b << "\n";
#endif

	int stat = 0;
        int par = 0;
        for( uint32_t a = 0; a < selection->size(); a++ ) {
	    MyDXFEntity *e = _entities[(*selection)(a)];
	    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	    if( !pe )
		throw Error( ERROR_LOCATION, "Not a path entity" );

#ifdef DEBUG_INSIDE_LOOP
	    std::cout << "Entity " << a << "\n";
	    pe->debug_print( std::cout );
#endif
	    stat = pe->ray_cross( x, y );

#ifdef DEBUG_INSIDE_LOOP
	    std::cout << "stat = " << stat << "\n";
#endif

	    if( stat == 1 ) 
		par = !par;
	    else if( stat == 2 ) 
		break;
        }
        
        if( stat != 2 )
            return( par );

#ifdef DEBUG_INSIDE_LOOP
	std::cout << "Perturbation\n";
#endif

	// Perturbation is more than double the uncertainity, which is
	// used in ray_crossing routines.
        x += 2.1*MYDXF_PERT_EPS;
    }

    throw Error( ERROR_LOCATION, "Perturbation failed" );
}


MyDXFEntitySelection *MyDXFEntities::selection_layer( const std::string &layername ) const
{
    MyDXFEntitySelection *selection = new MyDXFEntitySelection();

    for( size_t a = 0; a < _entities.size(); a++ ) {
	if( _entities[a]->get_layer() == layername )
	    selection->add_entity( a );
    }

    return( selection );
}


MyDXFEntitySelection *MyDXFEntities::selection_type( EntityType type ) const
{
    MyDXFEntitySelection *selection = new MyDXFEntitySelection();

    for( size_t a = 0; a < _entities.size(); a++ ) {
	if( _entities[a]->get_type() == type )
	    selection->add_entity( a );
    }

    return( selection );
}


MyDXFEntitySelection *MyDXFEntities::selection_all( void ) const
{
    MyDXFEntitySelection *selection = new MyDXFEntitySelection();

    for( size_t a = 0; a < _entities.size(); a++ ) {
	selection->add_entity( a );
    }

    return( selection );
}


std::ostream &operator<<( std::ostream &os, const MyDXFEntitySelection &sel )
{
    if( sel._selection.size() > 1 ) {
	os << "{";
	for( size_t a = 0; a < sel._selection.size()-1; a++ )
	    os << sel._selection[a] << ", ";
	os << sel._selection[sel._selection.size()-1] << "}";
    } else {
	os << "{" << sel._selection[0] << "}";
    }
    
    return( os );
}


bool MyDXFEntities::geom_same( uint32_t a, uint32_t b, double eps ) const
{
    // Test lines
    const MyDXFLine *line1 = dynamic_cast<MyDXFLine *>( _entities[a] );
    const MyDXFLine *line2 = dynamic_cast<MyDXFLine *>( _entities[b] );
    if( line1 && line2 )
	return( line1->geom_same( *line2, eps ) );

    // Test arcs
    const MyDXFArc *arc1 = dynamic_cast<MyDXFArc *>( _entities[a] );
    const MyDXFArc *arc2 = dynamic_cast<MyDXFArc *>( _entities[b] );
    if( arc1 && arc2 )
	return( arc1->geom_same( *arc2, eps ) );
    
    // Test circles
    const MyDXFCircle *circle1 = dynamic_cast<MyDXFCircle *>( _entities[a] );
    const MyDXFCircle *circle2 = dynamic_cast<MyDXFCircle *>( _entities[b] );
    if( circle1 && circle2 )
	return( circle1->geom_same( *circle2, eps ) );

    // Test lwpolylines
    const MyDXFLWPolyline *lwpline1 = dynamic_cast<MyDXFLWPolyline *>( _entities[a] );
    const MyDXFLWPolyline *lwpline2 = dynamic_cast<MyDXFLWPolyline *>( _entities[b] );
    if( lwpline1 && lwpline2 )
	return( lwpline1->geom_same( *lwpline2, eps ) );

    return( false );
}


MyDXFEntitySelection *MyDXFEntities::selection_path_loop( MyDXFEntitySelection *selection,
							  double eps )
{
    uint32_t a;
    MyDXFEntitySelection *subsel = new MyDXFEntitySelection();

    std::vector<bool>     stdir;     // Direction of entities on stack
    std::vector<uint32_t> stack;     // Stack of entities for loop
    bool done[selection->size()];    // List of processed entities

    // Initialize
    for( a = 0; a < selection->size(); a++ ) {
	done[a] = false;
    }

    /*
    std::cout << "*************************************************************\n";
    std::cout << "selection = {";
    for( a = 0; a < selection->size()-1; a++ )
	std::cout << (*selection)(a) << ", ";
    std::cout << (*selection)(a) << "}";
    */
 
#ifdef MYDXF_DEBUG
    std::cout << "Removing non-path objects and checking for self-looped objects\n";
#endif
   // Remove non-path objects and process self-looped path objects
    for( a = 0; a < selection->size(); a++ ) {
	MyDXFEntity *e = _entities[(*selection)(a)];
	MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	if( !pe ) {
	    done[a] = true;
	    continue;
	}

#ifdef MYDXF_DEBUG
	std::cout << "Testing object " << (*selection)(a) << ":\n";
	std::cout << "  start = " << pe->start() << "\n";
	std::cout << "  end = " << pe->end() << "\n";
#endif

	if( pe->start() == pe->end() ) {
	    // Add entity to list and remove duplicates
#ifdef MYDXF_DEBUG
	    std::cout << "Adding self-looped object " << (*selection)(a) << " to selection";
#endif
	    subsel->add_entity( (*selection)(a) );
	    done[a] = true;
	    for( uint32_t b = a+1; b < selection->size(); b++ ) {
		MyDXFEntity *eb = _entities[(*selection)(b)];
		MyDXFPathEntity *peb = dynamic_cast<MyDXFPathEntity *>( eb );
		if( peb && peb->start() == peb->end() ) {
		    if( geom_same( (*selection)(a), (*selection)(b), eps ) ) {
			done[b] = true;
			continue;
		    }
		}
	    }
	}
    }

    // Loop until all entities done
    while( 1 ) {

#ifdef MYDXF_DEBUG
	std::cout << "\n\ndone = {";
	for( a = 0; a < selection->size()-1; a++ )
	    std::cout << done[a] << ", ";
	std::cout << done[a] << "}\n";
#endif

	if( stack.size() == 0 ) {
	    // Start with first unprocessed of the selected entities,
	    // which is a path object
	    for( a = 0; a < selection->size(); a++ ) {
		if( !done[a] )
		    break;
	    }
	    if( a == selection->size() )
		break; // No entities left
	    // Add to stack
	    done[a] = true;
	    stdir.push_back( true );
	    stack.push_back( (*selection)(a) );

#ifdef MYDXF_DEBUG
	    std::cout << "Starting stack with " << (*selection)(a) << "\n";
#endif
	}

	// Check if loop done, check if endpoint of last entity
	// on stack matches some starting point on stack
	Vec3D end, start;
	if( stdir.back() ) {
	    MyDXFEntity *e = _entities[ stack.back() ];
	    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	    end = pe->end();
	} else {
	    MyDXFEntity *e = _entities[ stack.back() ];
	    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	    end = pe->start();
	}
#ifdef MYDXF_DEBUG
	std::cout << "Check if we have loop m. "
		  << "end = " << end << "\n";
#endif

	for( a = 0; a < stack.size(); a++ ) {
	    if( stdir[a] ) {
		MyDXFEntity *e = _entities[ stack[a] ];
		MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
		start = pe->start();
	    } else {
		MyDXFEntity *e = _entities[ stack[a] ];
		MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
		start = pe->end();
	    }
#ifdef MYDXF_DEBUG
	    std::cout << "  Stack[" << a << "]: "
		      << "stdir = " << stdir[a] << ", "
		      << "start = " << start << "\n";
#endif
	    // If loop closed
	    if( norm2( start-end ) < eps ) {

#ifdef MYDXF_DEBUG
		std::cout << "  Match found\n";
#endif
		// Fix all ends from stack and add to subsel list
		for( int32_t b = stack.size()-1; b >= (int32_t)a; b-- ) {

		    uint32_t c = stack.back();
		    MyDXFEntity *e = _entities[c];
		    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );

		    if( stdir.back() )
			pe->set_end( start );
		    else
			pe->set_start( start );

		    if( stdir[b] )
			start = pe->start();
		    else
			start = pe->end();
		    
		    subsel->add_entity( c );
		    stdir.pop_back();
		    stack.pop_back();
		}
		break;
	    }
	}
	if( stack.size() == 0 )
	    continue;

	// Search for an unprocessed and selected entity, which is a
	// path object and start matches with the end of the last
	// entity on stack
	if( stdir.back() ) {
	    MyDXFEntity *e = _entities[ stack.back() ];
	    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	    end = pe->end();
	} else {
	    MyDXFEntity *e = _entities[ stack.back() ];
	    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	    end = pe->start();
	}
#ifdef MYDXF_DEBUG
	std::cout << "Search for ent. matching end = "
		  << end << "\n";
#endif

	for( a = 0; a < selection->size(); a++ ) {

	    if( done[a] || (*selection)(a) == stack.back() ) 
		continue;

	    MyDXFEntity *e = _entities[ (*selection)(a) ];
	    MyDXFPathEntity *pe = dynamic_cast<MyDXFPathEntity *>( e );
	    start = pe->start();

#ifdef MYDXF_DEBUG
	    std::cout << "  Entity " << std::setw(4) << (*selection)(a) << ": "
		      << "        start = " << start << "\n";
#endif
	    if( norm2( start-end ) < eps ) {
#ifdef MYDXF_DEBUG
		std::cout << "  Match\n";
#endif
		// Check if entities are same within eps
		if( geom_same( stack.back(), (*selection)(a), eps ) ) {
#ifdef MYDXF_DEBUG
		    std::cout << "  Removing duplicate\n";
#endif
		    // Remove duplicate
		    done[a] = true;
		    continue;
		}
		// Add to stack
		stdir.push_back( true );
		stack.push_back( (*selection)(a) );
		done[a] = true;
		break;
	    } else {

		start = pe->end();
#ifdef MYDXF_DEBUG
		std::cout << "  No match, other end: "
			  << "  end = " << start << "\n";
#endif

		if( norm2( start-end ) < eps ) {
#ifdef MYDXF_DEBUG
		    std::cout << "  Match\n";
#endif
		    // Check if entities are same within eps
		    if( geom_same( stack.back(), (*selection)(a), eps ) ) {
#ifdef MYDXF_DEBUG
			std::cout << "  Removing duplicate\n";
#endif
			// Remove duplicate
			done[a] = true;
			continue;
		    }
		    // Add to stack backwards
		    stdir.push_back( false );
		    stack.push_back( (*selection)(a) );
		    done[a] = true;
		    break;
		}
#ifdef MYDXF_DEBUG
		std::cout << "  No match\n";
#endif
	    }
	}
	if( a == selection->size() ) {
#ifdef MYDXF_DEBUG
	    std::cout << "  No matching entity found, removing last of stack\n";
#endif
	    if( _dxf->wlevel() ) {
		ibsimu.message(1) << "No match at " << end << ", removing entity\n";
	    }

	    // No matching entity found. Remove last entity of stack
	    // and mark it done
	    stdir.pop_back();
	    stack.pop_back();
	}

    }

    return( subsel );
}


void MyDXFEntities::plot( const MyDXFEntitySelection *selection, const class MyDXFFile *dxf, 
			  cairo_t *cairo, const Transformation *t, const double range[4] ) const
{
#ifdef MYDXF_DEBUG_PLOT
    std::cout << "MyDXFEntities::plot()\n";
#endif

    if( selection ) {
	// Go through selection
	for( uint32_t a = 0; a < selection->size(); a++ ) {
	    MyDXFEntity *e = _entities[(*selection)(a)];
	    e->plot( dxf, cairo, t, range );
	}
    } else {
	// Plot all entities 
	for( uint32_t a = 0; a < _entities.size(); a++ ) {
	    MyDXFEntity *e = _entities[a];
	    e->plot( dxf, cairo, t, range );
	}
    }
}


void MyDXFEntities::get_bbox( const MyDXFEntitySelection *selection, Vec3D &min, Vec3D &max, 
			      const class MyDXFFile *dxf, const Transformation *t ) const
{
    min = Vec3D( std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    max = Vec3D( -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity() );

    if( selection ) {
	for( size_t a = 0; a < selection->size(); a++ ) {
	    MyDXFEntity *e = _entities[(*selection)(a)];
	    Vec3D mi, ma;
	    e->get_bbox( mi, ma, dxf, t );

	    for( int b = 0; b < 3; b++ ) {
		if( mi[b] < min[b] )
		    min[b] = mi[b];
		if( ma[b] > max[b] )
		    max[b] = ma[b];
	    }
	}
    } else {
	for( size_t a = 0; a < _entities.size(); a++ ) {
	    MyDXFEntity *e = _entities[a];
	    Vec3D mi, ma;
	    e->get_bbox( mi, ma, dxf, t );
	    for( int b = 0; b < 3; b++ ) {
		if( mi[b] < min[b] )
		    min[b] = mi[b];
		if( ma[b] > max[b] )
		    max[b] = ma[b];
	    }
	}
    }

#ifdef MYDXF_DEBUG_BBOX
    std::cout << "Entities bbox\n";
    std::cout << "min = " << min << "\n";
    std::cout << "max = " << max << "\n";
#endif
}


void MyDXFEntities::scale( MyDXFEntitySelection *selection, class MyDXFFile *dxf, double s )
{
    if( selection ) {
	// Go through selection
	for( size_t a = 0; a < selection->size(); a++ ) {
	    MyDXFEntity *e = _entities[(*selection)(a)];
	    e->scale( dxf, s );
	}    
    } else {
	// Scale all entities 
	for( size_t a = 0; a < _entities.size(); a++ ) {
	    MyDXFEntity *e = _entities[a];
	    e->scale( dxf, s );
	}
    }
}


void MyDXFEntities::translate( MyDXFEntitySelection *selection, class MyDXFFile *dxf, const Vec3D &dx )
{
    if( selection ) {
	// Go through selection
	for( size_t a = 0; a < selection->size(); a++ ) {
	    MyDXFEntity *e = _entities[(*selection)(a)];
	    e->translate( dxf, dx );
	}    
    } else {
	// Scale all entities 
	for( size_t a = 0; a < _entities.size(); a++ ) {
	    MyDXFEntity *e = _entities[a];
	    e->translate( dxf, dx );
	}
    }
}



void MyDXFEntities::remove( MyDXFEntitySelection *selection )
{
    if( selection ) {
	// Go through selection, remove and mark as NULL
	for( size_t a = 0; a < selection->size(); a++ ) {
	    delete _entities[(*selection)(a)];
	    _entities[(*selection)(a)] = NULL;
	}
	// Move entities to fill gaps
	size_t a, b = 0;
	for( a = 0; a < _entities.size(); a++ ) {
	    if( !_entities[a] ) {
		// Search for entity to move to location a
		if( b < a ) b = a+1;
		for( ; b < _entities.size(); b++ )
		    if( _entities[b] )
			break;
		if( b == _entities.size() ) {
		    // No more entities to move
		    break;
		}
		_entities[a] = _entities[b];
		_entities[b] = NULL;
	    }
	}
	// Resize vector
	_entities.resize( a );
    } else {
	// Remove all entities
	for( size_t a = 0; a < _entities.size(); a++ )
	    delete _entities[a];
	_entities.resize( 0 );
    }
}


void MyDXFEntities::explode( MyDXFEntitySelection *selection, class MyDXFFile *dxf )
{
    Transformation t;

    if( selection ) {
	// Go through selection
	size_t size = selection->size();
	for( size_t a = 0; a < size; a++ ) {
	    MyDXFInsert *ei = dynamic_cast<MyDXFInsert *>( _entities[(*selection)(a)] );
	    if( ei )
		ei->explode( this, dxf, &t );
	}

    } else {
	// Explode all
	size_t size = _entities.size();
	for( size_t a = 0; a < size; a++ ) {
	    MyDXFInsert *ei = dynamic_cast<MyDXFInsert *>( _entities[a] );
	    if( ei )
		ei->explode( this, dxf, &t );
	}
    }
}


void MyDXFEntities::explode( MyDXFEntities *ent, class MyDXFFile *dxf, const Transformation *t ) const
{
    size_t size = _entities.size();
    for( size_t a = 0; a < size; a++ )
	_entities[a]->explode( ent, dxf, t );
}


void MyDXFEntities::debug_print( std::ostream &os ) const
{
    os << "*** Section ENTITIES **************************************\n";

    MyDXFFile::debug_print_format( os, "number_of_entities", (int)_entities.size() );
    os << "\n";
    
    for( size_t a = 0; a < _entities.size(); a++ ) {
	MyDXFEntity *e = _entities[a];
	e->debug_print( os );
	os << "\n";
    }
    
    os << "\n";
}




