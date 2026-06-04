/*! \file lineclip.cpp
 *  \brief Floating point line clipping for cairo
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

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <limits>
#include "lineclip.hpp"
#include "error.hpp"


#define OP_LINE 0
#define OP_MOVE 1
#define OP_NULL 2


//#define LINECLIP_DEBUG
//extern int ppdbg;

LineClip::LineClip( cairo_t *dc )
{
    p_dc           = dc;

    clip[0]        = -std::numeric_limits<double>::infinity();
    clip[1]        = -std::numeric_limits<double>::infinity();
    clip[2]        = std::numeric_limits<double>::infinity();
    clip[3]        = std::numeric_limits<double>::infinity();

    last[0]        = std::numeric_limits<double>::quiet_NaN();
    last[1]        = std::numeric_limits<double>::quiet_NaN();
    last_outcode   = 0;
    last_op        = OP_NULL;

    drawn[0]       = std::numeric_limits<double>::quiet_NaN();
    drawn[1]       = std::numeric_limits<double>::quiet_NaN();
    drawn_outcode  = 0;

    coord_alloc    = 0;
    coord          = 0;
}


LineClip::~LineClip()
{
    free( coord );
}


void LineClip::set( double xmin, double ymin, double xmax, double ymax )
{
    // Sort clip limits
    if( xmin < xmax ) {
	clip[0]        = xmin;
	clip[2]        = xmax;
    } else {
	clip[0]        = xmax;
	clip[2]        = xmin;
    }
    if( ymin < ymax ) {
	clip[1]        = ymin;
	clip[3]        = ymax;
    } else {
	clip[1]        = ymax;
	clip[3]        = ymin;
    }

    last[0]        = std::numeric_limits<double>::quiet_NaN();
    last[1]        = std::numeric_limits<double>::quiet_NaN();
    last_outcode   = 0;
    last_op        = OP_NULL;

    drawn[0]       = std::numeric_limits<double>::quiet_NaN();
    drawn[1]       = std::numeric_limits<double>::quiet_NaN();
    drawn_outcode  = 0;
}


void LineClip::reset()
{
    clip[0]        = -std::numeric_limits<double>::infinity();
    clip[1]        = -std::numeric_limits<double>::infinity();
    clip[2]        = std::numeric_limits<double>::infinity();
    clip[3]        = std::numeric_limits<double>::infinity();

    last[0]        = std::numeric_limits<double>::quiet_NaN();
    last[1]        = std::numeric_limits<double>::quiet_NaN();
    last_outcode   = 0;
    last_op        = OP_NULL;

    drawn[0]       = std::numeric_limits<double>::quiet_NaN();
    drawn[1]       = std::numeric_limits<double>::quiet_NaN();
    drawn_outcode  = 0;
}


int LineClip::outcode( double x, double y )
{
    int outcode = 0;

    if( y > clip[3] )
	outcode |= 8;
    else if( y < clip[1] )
	outcode |= 4;
    if( x > clip[2] )
	outcode |= 2;
    else if( x < clip[0] )
	outcode |= 1;

    return( outcode );
}


int LineClip::exit_outcode( double x, double y )
{
    int outcode = 0;

    if( y >= clip[3] )
	outcode = 8;
    else if( y <= clip[1] )
	outcode = 4;
    else if( x >= clip[2] )
	outcode = 2;
    else if( x <= clip[0] )
	outcode = 1;

    return( outcode );
}


void LineClip::move_to( double x, double y )
{
    drawn_outcode = 0;

    first[0] = x;
    first[1] = y;

    last_op = OP_MOVE;
    last[0] = x;
    last[1] = y;
    last_outcode = outcode( x, y );
}


#define AROUND_5(i,j)				\
    ( (((i)|(j)) == 5) ||			\
      ((i) ==  9 && (j) ==  4) ||		\
      ((i) ==  4 && (j) ==  9) ||		\
      ((i) ==  6 && (j) ==  1) ||		\
      ((i) ==  1 && (j) ==  6) )

#define AROUND_6(i,j)				\
    ( (((i)|(j)) == 6) ||			\
      ((i) ==  5 && (j) ==  2) ||		\
      ((i) ==  2 && (j) ==  5) ||		\
      ((i) == 10 && (j) ==  4) ||		\
      ((i) ==  4 && (j) == 10) )

#define AROUND_10(i,j)				\
    ( (((i)|(j)) == 10) ||			\
      ((i) ==  8 && (j) ==  6) ||		\
      ((i) ==  6 && (j) ==  8) ||		\
      ((i) ==  9 && (j) ==  2) ||		\
      ((i) ==  2 && (j) ==  9) )

#define AROUND_9(i,j)				\
    ( (((i)|(j)) == 9) ||			\
      ((i) == 10 && (j) ==  1) ||		\
      ((i) ==  1 && (j) == 10) ||		\
      ((i) ==  5 && (j) ==  8) ||		\
      ((i) ==  8 && (j) ==  5) )



void LineClip::line_to( double x, double y )
{
    int outcodeout;
    int outcode1_orig;
    int outcode2_orig;
    int outcode1, outcode2;
    double x1, y1, x2, y2;
    double xt, yt;

    x1 = last[0];
    y1 = last[1];
    outcode1_orig = outcode1 = last_outcode;
    x2 = x;
    y2 = y;
    outcode2_orig = outcode2 = outcode( x, y );

    //if( ppdbg ) {
    //	std::cout << "\n-------------------------------\nline_to() started\n";
    //	std::cout << "Line from " << outcode1_orig << ": (" << x1 << ", " << y1 << ")\n";
    //	std::cout << "Line to   " << outcode2_orig << ": (" << x2 << ", " << y2 << ")\n";
    //	std::cout << "clip[0] = " << clip[0] << ", " << "clip[1] = " << clip[1] << "\n";
    //	std::cout << "clip[2] = " << clip[2] << ", " << "clip[3] = " << clip[3] << "\n";
    //  }
  
    /* Do line clipping iteration */
    while( 1 ) {

	//if( ppdbg )
	//    std::cout << "\nClipping loop\n";

	if( (outcode1 | outcode2) == 0 ) {
	    /* Trivial inside */
	    //if( ppdbg )
	    //std::cout << "Trivial inside\n";
	    break;
	} else if( (outcode1 & outcode2) != 0 ) {
	    /* Trivial outside */
	    //if( ppdbg )
	    //std::cout << "Trivial outside\n";
	    break;
	} else {
	    /* Non-trivial case, at least one point outside clip region, do
	     * iteration */
	    //if( ppdbg )
	    //std::cout << "Nontrivial case\n";
	    
	    if( outcode1 != 0 ) outcodeout = outcode1;
	    else outcodeout = outcode2;

	    if( (outcodeout & 8) == 8 ) {
		/* Divide line at top of clip area */
		double tt = (clip[3] - y1)/(y2 - y1);
		xt = x2*tt + x1*(1.0-tt);
		yt = clip[3];
	    } else if( (outcodeout & 4) == 4 ) {
		/* Divide line at bottom of clip area */
		double tt = (clip[1] - y1)/(y2 - y1);
		xt = x2*tt + x1*(1.0-tt);
		yt = clip[1];
	    } else if( (outcodeout & 2) == 2 ) {
		/* Divide line at right of clip area */
		double tt = (clip[2] - x1)/(x2 - x1);
		yt = y2*tt + y1*(1.0-tt);
		xt = clip[2];
	    } else {
		/* Divide line at left of clip area */
		double tt = (clip[0] - x1)/(x2 - x1);
		yt = y2*tt + y1*(1.0-tt);
		xt = clip[0];
	    }

	    //if( ppdbg ) {
	    //std::cout << "xt = " << xt << "\n";
	    //std::cout << "yt = " << yt << "\n";
	    //}

	    /* Make outside point equal to the intersection point and redo
	     * clipping for this new point */
	    if( outcodeout == outcode1 ) {
		x1 = xt;
		y1 = yt;
		outcode1 = outcode( x1, y1 );
	    } else {
		x2 = xt;
		y2 = yt;
		outcode2 = outcode( x2, y2 );
	    }
	}
    }

    //if( ppdbg ) {
    //std::cout << "Processed data:\n";
    //std::cout << "Line from " << outcode1 << ": (" << x1 << ", " << y1 << ")\n";
    //std::cout << "Line to   " << outcode2 << ": (" << x2 << ", " << y2 << ")\n";
    //}

    /* Process different cases */
    if( (outcode1_orig | outcode2_orig) == 0 ) {
	/* From inside to inside *************************************************** */
	//if( ppdbg )
	//printf( "Inside to inside\n" );
	if( last_op == OP_MOVE )
	    cairo_move_to( p_dc, x1, y1 );
	cairo_line_to( p_dc, x2, y2 );

    } else if( outcode1_orig == 0 ) {
	/* From inside to outside ************************************************** */
	//if( ppdbg )
	//printf( "Inside to outside\n" );
	if( last_op == OP_MOVE )
	    cairo_move_to( p_dc, x1, y1 );
	cairo_line_to( p_dc, x2, y2 );    
	drawn[0] = x2;
	drawn[1] = y2;
	drawn_outcode = exit_outcode( x2, y2 );

    } else if( outcode2_orig == 0 || (outcode1 | outcode2) == 0 ) {
	/* From outside to inside or outside to outside, through inside ************ */
	//if( ppdbg )
	//printf( "Outside to inside\n" );
	if( last_op == OP_MOVE ) {
	    cairo_move_to( p_dc, x1, y1 );
	    cairo_line_to( p_dc, x2, y2 );
	} else {
	    //if( ppdbg ) {
	    //printf( "drawn_outcode = %d\n", drawn_outcode );
	    //printf( "exit_outcode = %d\n", exit_outcode( x1, y1 ) );
	    //}
	    if( drawn_outcode &&
		(drawn_outcode & exit_outcode( x1, y1 )) == 0 ) {
		/* No common edge with last drawn and first to be drawn ->
		 * needs two connecting lines */
		//if( ppdbg )
		//printf( "Going around a corner\n" );
		outcodeout = exit_outcode( x1, y1 );
		if( AROUND_5( drawn_outcode, outcodeout ) )
		    cairo_line_to( p_dc, clip[0], clip[1] );
		else if( AROUND_6( drawn_outcode, outcodeout ) )
		    cairo_line_to( p_dc, clip[2], clip[1] );
		else if( AROUND_9( drawn_outcode, outcodeout ) )
		    cairo_line_to( p_dc, clip[0], clip[3] );
		else if( AROUND_10( drawn_outcode, outcodeout ) )
		    cairo_line_to( p_dc, clip[2], clip[3] );
		else {
		    // TODO
		    throw( Error( ERROR_LOCATION, "illegal outcodes" ) ); 
		    //throw( int(1) );
		}
	    }
	    //if( ppdbg )
	    //printf( "Line: (%g, %g) to (%g, %g)\n", x1, y1, x2, y2 );
	    cairo_line_to( p_dc, x1, y1 );
	    cairo_line_to( p_dc, x2, y2 );
	}

	/* Process from inside to outside again */
	if( (outcode1 | outcode2) == 0 ) {
	    drawn[0] = x2;
	    drawn[1] = y2;
	    drawn_outcode = exit_outcode( x2, y2 );
	}

    } else {
	/* Outside only ************************************************************ */
	//if( ppdbg ) 
	//std::cout << "Outside only\n";

	if( last_op == OP_MOVE ) {
	    /* Start drawing at border */
	    if( outcode1_orig == 5 ) {
		xt = clip[0];
		yt = clip[1];
	    } else if( outcode1_orig == 1 ) {
		xt = clip[0];
		yt = last[1];
	    } else if( outcode1_orig == 9 ) {
		xt = clip[0];
		yt = clip[3];
	    } else if( outcode1_orig == 8 ) {
		xt = last[0];
		yt = clip[3];
	    } else if( outcode1_orig == 10 ) {
		xt = clip[2];
		yt = clip[3];
	    } else if( outcode1_orig == 2 ) {
		xt = clip[2];
		yt = last[1];
	    } else if( outcode1_orig == 6 ) {
		xt = clip[2];
		yt = clip[1];
	    } else {
		xt = last[0];
		yt = clip[1];
	    }

	    cairo_move_to( p_dc, xt, yt );
	    last[0] = xt;
	    last[1] = yt;
	    last_outcode = drawn_outcode = outcode1_orig;
	}

	//if( ppdbg ) 
	//printf( "%d to %d\n", drawn_outcode, outcode2_orig );

	if( (drawn_outcode & outcode2_orig) == 0 ) {
	    /* No common edge with last drawn and last point of this line ->
	     * draw line to common corner */

	    if( (drawn_outcode | outcode2_orig) == 15 ) {
		/* Going across */

		if( drawn_outcode != last_outcode ) {
		    /* Get a hint from last point */
		    if( outcode2_orig == 9 || outcode2_orig == 6 ) {
			if( last_outcode == 2 || last_outcode == 10 ||
			    last_outcode == 8 ) {
			    xt = clip[2];
			    yt = clip[3];
			    outcodeout = 10;
			} else {
			    xt = clip[0];
			    yt = clip[1];
			    outcodeout = 5;
			}
		    }
		    else /* outcode2_orig == 5 || outcode2_orig == 10 */ {
			if( last_outcode == 4 || last_outcode == 6 ||
			    last_outcode == 2 ) {
			    xt = clip[2];
			    yt = clip[1];
			    outcodeout = 6;
			} else {
			    xt = clip[0];
			    yt = clip[3];
			    outcodeout = 9;
			}
		    }

		} else {
		    /* Get a hint from clipping mid-points */

		    if( outcode2_orig == 9 || outcode2_orig == 6 ) {
			if( y1 >= clip[3] + (x1-clip[0])/(clip[2]-clip[0])*
			    (clip[1]-clip[3]) ) {
			    xt = clip[2];
			    yt = clip[3];
			    outcodeout = 10;
			} else {
			    xt = clip[0];
			    yt = clip[1];
			    outcodeout = 5;
			}
		    } else /* outcode2_orig == 5 || outcode2_orig == 10 */ {
			if( y1 >= clip[1] + (x1-clip[0])/(clip[2]-clip[0])*
			    (clip[3]-clip[1]) ) {
			    xt = clip[0];
			    yt = clip[3];
			    outcodeout = 9;
			} else {
			    xt = clip[2];
			    yt = clip[1];
			    outcodeout = 6;
			}
		    }
		}
		
	    } else {
		/* Going around a corner (not across) */

		if( AROUND_5( drawn_outcode, outcode2_orig ) ) {
		    xt = clip[0];
		    yt = clip[1];
		    outcodeout = 5;
		} else if( AROUND_6( drawn_outcode, outcode2_orig ) ) {
		    xt = clip[2];
		    yt = clip[1];
		    outcodeout = 6;
		} else if( AROUND_9( drawn_outcode, outcode2_orig ) ) {
		    xt = clip[0];
		    yt = clip[3];
		    outcodeout = 9;
		} else if( AROUND_10( drawn_outcode, outcode2_orig ) ) {
		    xt = clip[2];
		    yt = clip[3];
		    outcodeout = 10;
		} else {
		    // TODO
		    throw( Error( ERROR_LOCATION, "illegal outcodes" ) ); 
		    //throw( int(1) );
		}
	    }

	    cairo_line_to( p_dc, xt, yt );
	    drawn[0] = xt;
	    drawn[1] = yt;
	    drawn_outcode = outcodeout;

	    //if( ppdbg ) 
	    //printf( "Connected to corner at %d: (%g, %g)\n", outcodeout, xt, yt );
	}
    }

    /* Save new current point and outcode */
    last[0] = x;
    last[1] = y;
    last_outcode = outcode2_orig;
    last_op = OP_LINE;
}


/* Save a point at t (t=[0,1]) on cubic Bezier curve to array coords. */
    void LineClip::get_point( double *coords, double t,
  double x0, double y0,
  double x1, double y1,
  double x2, double y2,
  double x3, double y3 )
{
    double t2 = t*t;
    double t3 = t2*t;
    double it = (1.0-t);
    double it2 = it*it;
    double it3 = it2*it;

    coords[0] = it3*x0 + 3*t*it2*x1 + 3*t2*it*x2 + t3*x3;
    coords[1] = it3*y0 + 3*t*it2*y1 + 3*t2*it*y2 + t3*y3;
}


void LineClip::curve_to
( double x2, double y2,
  double x3, double y3,
  double x4, double y4 )
{
    double x1, y1;
    double xt, yt;
    int outcode1_orig, outcode2_orig;
    int outcodeout;

    /* Get last point */
    x1 = last[0];
    y1 = last[1];
    outcode1_orig = last_outcode;
    outcode2_orig = outcode( x4, y4 );

    /* Check trivially out cases */
    if( (x1 < clip[0] && x2 < clip[0] && x3 < clip[0] && x4 < clip[0]) ||
	(x1 > clip[2] && x2 > clip[2] && x3 > clip[2] && x4 > clip[2]) ||
	(y1 < clip[1] && y2 < clip[1] && y3 < clip[1] && y4 < clip[1]) ||
	(y1 > clip[3] && y2 > clip[3] && y3 > clip[3] && y4 > clip[3]) ) {

	//printf( "Trivially out\n" );

	if( last_op == OP_MOVE ) {
	    if( outcode1_orig == 5 ) {
		xt = clip[0];
		yt = clip[1];
	    } else if( outcode1_orig == 1 ) {
		xt = clip[0];
		yt = last[1];
	    } else if( outcode1_orig == 9 ) {
		xt = clip[0];
		yt = clip[3];
	    } else if( outcode1_orig == 8 ) {
		xt = last[0];
		yt = clip[3];
	    } else if( outcode1_orig == 10 ) {
		xt = clip[2];
		yt = clip[3];
	    } else if( outcode1_orig == 2 ) {
		xt = clip[2];
		yt = last[1];
	    } else if( outcode1_orig == 6 ) {
		xt = clip[2];
		yt = clip[1];
	    } else {
		xt = last[0];
		yt = clip[1];
	    }

	    cairo_move_to( p_dc, xt, yt );
	    last[0] = xt;
	    last[1] = yt;
	    last_outcode = drawn_outcode = outcode1_orig;
	}

	//printf( "%d to %d\n", drawn_outcode, outcode2_orig );
	if( (drawn_outcode & outcode2_orig) == 0 ) {
	    /* No common edge with last drawn and last point of this line ->
	     * draw line to common corner */

	    if( (drawn_outcode | outcode2_orig) == 15 ) {
		/* Going across */

		if( drawn_outcode != last_outcode ) {
		    /* Get a hint from last point */
		    if( outcode2_orig == 9 || outcode2_orig == 6 ) {
			if( last_outcode == 2 || last_outcode == 10 ||
			    last_outcode == 8 ) {
			    xt = clip[2];
			    yt = clip[3];
			    outcodeout = 10;
			} else {
			    xt = clip[0];
			    yt = clip[1];
			    outcodeout = 5;
			}
		    }
		    else /* outcode2_orig == 5 || outcode2_orig == 10 */ {
			if( last_outcode == 4 || last_outcode == 6 ||
			    last_outcode == 2 ) {
			    xt = clip[2];
			    yt = clip[1];
			    outcodeout = 6;
			} else {
			    xt = clip[0];
			    yt = clip[3];
			    outcodeout = 9;
			}
		    }

		} else {
		    /* Get a hint from clipping mid-points */

		    if( outcode2_orig == 9 || outcode2_orig == 6 ) {
			if( y1 >= clip[3] + (x1-clip[0])/(clip[2]-clip[0])*
			    (clip[1]-clip[3]) ) {
			    xt = clip[2];
			    yt = clip[3];
			    outcodeout = 10;
			} else {
			    xt = clip[0];
			    yt = clip[1];
			    outcodeout = 5;
			}
		    } else /* outcode2_orig == 5 || outcode2_orig == 10 */ {
			if( y1 >= clip[1] + (x1-clip[0])/(clip[2]-clip[0])*
			    (clip[3]-clip[1]) ) {
			    xt = clip[0];
			    yt = clip[3];
			    outcodeout = 9;
			} else {
			    xt = clip[2];
			    yt = clip[1];
			    outcodeout = 6;
			}
		    }
		}
	
	    } else {
		/* Going around a corner (not across) */

		if( AROUND_5( drawn_outcode, outcode2_orig ) ) {
		    xt = clip[0];
		    yt = clip[1];
		    outcodeout = 5;
		} else if( AROUND_6( drawn_outcode, outcode2_orig ) ) {
		    xt = clip[2];
		    yt = clip[1];
		    outcodeout = 6;
		} else if( AROUND_9( drawn_outcode, outcode2_orig ) ) {
		    xt = clip[0];
		    yt = clip[3];
		    outcodeout = 9;
		} else if( AROUND_10( drawn_outcode, outcode2_orig ) ) {
		    xt = clip[2];
		    yt = clip[3];
		    outcodeout = 10;
		} else {
		    // TODO
		    throw( Error( ERROR_LOCATION, "illegal outcodes" ) ); 
		    //throw( int(1) );
		}
	    }

	    cairo_line_to( p_dc, xt, yt );
	    drawn[0] = xt;
	    drawn[1] = yt;
	    drawn_outcode = outcodeout;
	    //printf( "Connected to corner at %d: (%g, %g)\n", outcodeout, xt, yt );
	}

	last_op = OP_LINE;
	last[0] = x4;
	last[1] = y4;
	last_outcode = outcode2_orig;
	return;
    }

    /* Chop curve into lines ********************************************************* */
    int pass;
    double err, maxerr;
    int coord_N, a;
  
    /* Allocate memory for initial chop to three points */
    if( coord_alloc < 3 ) {
	coord_alloc = 3;
	if( !(coord = (double *)realloc( coord, 2*coord_alloc*sizeof(double) )) ) {
	    // TODO
	    throw( Error( ERROR_LOCATION, "memory allocation error" ) ); 
	    //throw( int(1) );
	}
    }

    /* Fill database with three points */
    coord_N = 3;
    get_point( &coord[0], 0.0, x1, y1, x2, y2, x3, y3, x4, y4 );
    get_point( &coord[2], 0.5, x1, y1, x2, y2, x3, y3, x4, y4 );
    get_point( &coord[4], 1.0, x1, y1, x2, y2, x3, y3, x4, y4 );

    /* Loop while error is small enough */
    pass = 0;
    while( 1 ) {

	/* Calculate maximum error between the most recent points and
	 * the linear interpolation of points from the last round. */
	maxerr = 0.0;
	for( a = 1; a < coord_N; a += 2 ) {
	    xt = 0.5*(coord[2*(a-1)] + coord[2*(a+1)]);
	    yt = 0.5*(coord[2*(a-1)+1] + coord[2*(a+1)+1]);
	    xt -= coord[2*a];
	    yt -= coord[2*a+1];
	    err = sqrt( xt*xt + yt*yt );
	    if( err > maxerr )
		maxerr = err;
	}

	/* If error is less than 5 pixels two iterations in a row, it
	 * is good enough. Please note that the absolute error of the
	 * last iteration is typically less than 1 pixel at this level
	 * of iteration.
	 */
	if( maxerr < 5.0 ) {
	    if( pass == 1 )
		break;
	    pass = 1;
	} else
	    pass = 0;

	/* Add more points */
	coord_N += coord_N-1;
    
	/* Allocate space */
	if( coord_alloc < coord_N ) {
	    coord_alloc = coord_N;
	    if( !(coord = (double *)realloc( coord, 2*coord_alloc*sizeof(double) )) ) {
		// TODO
		throw( Error( ERROR_LOCATION, "memory allocation error" ) ); 
		//throw( int(1) );
	    }
	}
    
	/* Fill coordinates */
	for( a = coord_N-1; a > 0; a-- ) {
      
	    /* Even */
	    coord[2*a] = coord[a];
	    coord[2*a+1] = coord[a+1];
	    a--;
      
	    /* Odd */
	    get_point( &coord[2*a], (double)a/(coord_N-1), 
				x1, y1, x2, y2, x3, y3, x4, y4 );
	}
    }

    /* Error small enough -> draw lines */
    //printf( "Drawing with %d lines\n", coord_N-1 );
    for( a = 1; a < coord_N; a++ )
	line_to( coord[2*a], coord[2*a+1] );
}


void LineClip::close_path()
{
    /* Close path */
    line_to( first[0], first[1] );
}


void LineClip::fill()
{
    /* Close path and fill */
    line_to( first[0], first[1] );
    cairo_fill( p_dc );
}
