/*! \file mydxfheader.cpp
 *  \brief DXF Header
 */

/* Copyright (c) 2010,2012 Taneli Kalvas. All rights reserved.
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
#include "mydxfheader.hpp"
#include "error.hpp"


MyDXFHeader::MyDXFHeader( class MyDXFFile *dxf )
    : angbase(0.0), angdir(0), dimasz(0.0), dimgap(0.0), dimexo(0.0),
      dimexe(0.0), dimtxt(0.0), insunits(4), orthomode(0), pucsorthoview(0),
      ucsorthoview(0), worldview(0)
{
#ifdef MYDXF_DEBUG
    std::cout << "Reading section HEADER\n";
#endif

    // Read HEADER section
    //
    // Contains packages of groups:
    // A. Variable name: code=9, data=header name (string)
    // B. Variable data: code(s) and data according to field(s)
    //
    // Ends in ENDSEC

    std::string field;
    while( dxf->read_group() != -1 ) {
	
	if( dxf->group_get_code() == 0 && dxf->group_get_string() == "ENDSEC" )
	    break; // Done with header

	else if( dxf->group_get_code() == 9 )
	    field = dxf->group_get_string(); // Read header variable name

	else if( dxf->group_get_code() == 1 && field == "$ACADVER" )
	    acadver = dxf->group_get_string();
	else if( dxf->group_get_code() == 50 && field == "$ANGBASE" )
	    angbase = dxf->group_get_double();
	else if( dxf->group_get_code() == 70 && field == "$ANGDIR" )
	    angdir = dxf->group_get_int16();
	else if( dxf->group_get_code() == 70 && field == "$INSUNITS" )
	    insunits = dxf->group_get_int16();
	else if( dxf->group_get_code() == 5 && field == "$HANDSEED" )
	    handseed = dxf->group_get_string();
	else if( dxf->group_get_code() == 40 && field == "$DIMASZ" )
	    dimasz = dxf->group_get_double();
	else if( dxf->group_get_code() == 40 && field == "$DIMGAP" )
	    dimgap = dxf->group_get_double();
	else if( dxf->group_get_code() == 40 && field == "$DIMEXO" )
	    dimexo = dxf->group_get_double();
	else if( dxf->group_get_code() == 40 && field == "$DIMEXE" )
	    dimexe = dxf->group_get_double();
	else if( dxf->group_get_code() == 40 && field == "$DIMTXT" )
	    dimtxt = dxf->group_get_double();


	else if( dxf->group_get_code() == 10 && field == "$PLIMMAX" )
	    plimmax[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$PLIMMAX" )
	    plimmax[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 10 && field == "$PLIMMIN" )
	    plimmin[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 10 && field == "$PLIMMIN" )
	    plimmin[1] = dxf->group_get_double();


	else if( dxf->group_get_code() == 70 && field == "$ORTHOMODE" )
	    orthomode = dxf->group_get_int16();


	else if( dxf->group_get_code() == 2 && field == "$PUCSBASE" )
	    pucsbase = dxf->group_get_string();
	else if( dxf->group_get_code() == 2 && field == "$PUCSNAME" )
	    pucsname = dxf->group_get_string();

	else if( dxf->group_get_code() == 10 && field == "$PUCSORG" )
	    pucsorg[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$PUCSORG" )
	    pucsorg[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$PUCSORG" )
	    pucsorg[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 10 && field == "$PUCSORGBACK" )
	    pucsorgback[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$PUCSORGBACK" )
	    pucsorgback[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$PUCSORGBACK" )
	    pucsorgback[2] = dxf->group_get_double();
	
	else if( dxf->group_get_code() == 10 && field == "$PUCSORGBOTTOM" )
	    pucsorgbottom[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$PUCSORGBOTTOM" )
	    pucsorgbottom[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$PUCSORGBOTTOM" )
	    pucsorgbottom[2] = dxf->group_get_double();
	
	else if( dxf->group_get_code() == 10 && field == "$PUCSORGFRONT" )
	    pucsorgfront[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$PUCSORGFRONT" )
	    pucsorgfront[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$PUCSORGFRONT" )
	    pucsorgfront[2] = dxf->group_get_double();
	
	else if( dxf->group_get_code() == 10 && field == "$PUCSORGLEFT" )
	    pucsorgleft[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$PUCSORGLEFT" )
	    pucsorgleft[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$PUCSORGLEFT" )
	    pucsorgleft[2] = dxf->group_get_double();
	
	else if( dxf->group_get_code() == 10 && field == "$PUCSORGRIGHT" )
	    pucsorgright[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$PUCSORGRIGHT" )
	    pucsorgright[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$PUCSORGRIGHT" )
	    pucsorgright[2] = dxf->group_get_double();
	
	else if( dxf->group_get_code() == 10 && field == "$PUCSORGTOP" )
	    pucsorgtop[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$PUCSORGTOP" )
	    pucsorgtop[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$PUCSORGTOP" )
	    pucsorgtop[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 2 && field == "$PUCSORTHOREF" )
	    pucsorthoref = dxf->group_get_string();
	else if( dxf->group_get_code() == 70 && field == "$PUCSORTHOVIEW" )
	    pucsorthoview = dxf->group_get_int16();

	else if( dxf->group_get_code() == 10 && field == "$PUCSXDIR" )
	    pucsxdir[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$PUCSXDIR" )
	    pucsxdir[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$PUCSXDIR" )
	    pucsxdir[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 10 && field == "$PUCSYDIR" )
	    pucsydir[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$PUCSYDIR" )
	    pucsydir[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$PUCSYDIR" )
	    pucsydir[2] = dxf->group_get_double();



	else if( dxf->group_get_code() == 2 && field == "$UCSBASE" )
	    ucsbase = dxf->group_get_string();
	else if( dxf->group_get_code() == 2 && field == "$UCSNAME" )
	    ucsname = dxf->group_get_string();

	else if( dxf->group_get_code() == 10 && field == "$UCSORG" )
	    ucsorg[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$UCSORG" )
	    ucsorg[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$UCSORG" )
	    ucsorg[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 10 && field == "$UCSORGBACK" )
	    ucsorgback[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$UCSORGBACK" )
	    ucsorgback[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$UCSORGBACK" )
	    ucsorgback[2] = dxf->group_get_double();
	
	else if( dxf->group_get_code() == 10 && field == "$UCSORGBOTTOM" )
	    ucsorgbottom[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$UCSORGBOTTOM" )
	    ucsorgbottom[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$UCSORGBOTTOM" )
	    ucsorgbottom[2] = dxf->group_get_double();
	
	else if( dxf->group_get_code() == 10 && field == "$UCSORGFRONT" )
	    ucsorgfront[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$UCSORGFRONT" )
	    ucsorgfront[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$UCSORGFRONT" )
	    ucsorgfront[2] = dxf->group_get_double();
	
	else if( dxf->group_get_code() == 10 && field == "$UCSORGLEFT" )
	    ucsorgleft[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$UCSORGLEFT" )
	    ucsorgleft[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$UCSORGLEFT" )
	    ucsorgleft[2] = dxf->group_get_double();
	
	else if( dxf->group_get_code() == 10 && field == "$UCSORGRIGHT" )
	    ucsorgright[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$UCSORGRIGHT" )
	    ucsorgright[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$UCSORGRIGHT" )
	    ucsorgright[2] = dxf->group_get_double();
	
	else if( dxf->group_get_code() == 10 && field == "$UCSORGTOP" )
	    ucsorgtop[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$UCSORGTOP" )
	    ucsorgtop[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$UCSORGTOP" )
	    ucsorgtop[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 2 && field == "$UCSORTHOREF" )
	    ucsorthoref = dxf->group_get_string();
	else if( dxf->group_get_code() == 70 && field == "$UCSORTHOVIEW" )
	    ucsorthoview = dxf->group_get_int16();

	else if( dxf->group_get_code() == 10 && field == "$UCSXDIR" )
	    ucsxdir[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$UCSXDIR" )
	    ucsxdir[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$UCSXDIR" )
	    ucsxdir[2] = dxf->group_get_double();

	else if( dxf->group_get_code() == 10 && field == "$UCSYDIR" )
	    ucsydir[0] = dxf->group_get_double();
	else if( dxf->group_get_code() == 20 && field == "$UCSYDIR" )
	    ucsydir[1] = dxf->group_get_double();
	else if( dxf->group_get_code() == 30 && field == "$UCSYDIR" )
	    ucsydir[2] = dxf->group_get_double();



	else if( dxf->group_get_code() == 70 && field == "$WORLDVIEW" )
	    worldview = dxf->group_get_int16();
	


	else if( dxf->wlevel() >= 2 )
	    std::cout << "Skipping unknown header \'" << field << "\'\n";	
    }

}



MyDXFHeader::~MyDXFHeader()
{

}


void MyDXFHeader::write( class MyDXFFile *dxf, std::ofstream &_ostr )
{
    dxf->write_group(  0, "SECTION" );
    dxf->write_group(  2, "HEADER" );

    dxf->write_group(  9, "$ACADVER" );
    dxf->write_group(  1, acadver.c_str() );

    dxf->write_group(  9, "$ANGBASE" );
    dxf->write_group( 50, angbase );

    dxf->write_group(  9, "$ANGDIR" );
    dxf->write_group( 70, angdir );

    dxf->write_group(  9, "$HANDSEED" );
    dxf->write_group(  5, handseed.c_str() );

    dxf->write_group(  9, "$DIMASZ" );
    dxf->write_group( 40, dimasz );

    dxf->write_group(  9, "$DIMGAP" );
    dxf->write_group( 40, dimgap );

    dxf->write_group(  9, "$DIMEXO" );
    dxf->write_group( 40, dimexo );

    dxf->write_group(  9, "$DIMEXE" );
    dxf->write_group( 40, dimexe );

    dxf->write_group(  9, "$DIMTXT" );
    dxf->write_group( 40, dimtxt );

    dxf->write_group(  9, "$INSUNITS" );
    dxf->write_group( 70, insunits );

    dxf->write_group(  9, "$PLIMMAX" );
    dxf->write_group( 10, plimmax[0] );
    dxf->write_group( 20, plimmax[1] );

    dxf->write_group(  9, "$PLIMMIN" );
    dxf->write_group( 10, plimmin[0] );
    dxf->write_group( 20, plimmin[1] );

    dxf->write_group(  9, "$ORTHOMODE" );
    dxf->write_group( 70, orthomode );

    dxf->write_group(  9, "$PUCSBASE" );
    dxf->write_group(  2,  pucsbase.c_str() );

    dxf->write_group(  9, "$PUCSNAME" );
    dxf->write_group(  2,  pucsname.c_str() );

    dxf->write_group(  9, "$PUCSORG" );
    dxf->write_group( 10,  pucsorg[0] );
    dxf->write_group( 20,  pucsorg[1] );
    dxf->write_group( 30,  pucsorg[2] );

    dxf->write_group(  9, "$PUCSORGBACK" );
    dxf->write_group( 10,  pucsorgback[0] );
    dxf->write_group( 20,  pucsorgback[1] );
    dxf->write_group( 30,  pucsorgback[2] );

    dxf->write_group(  9, "$PUCSORGBOTTOM" );
    dxf->write_group( 10,  pucsorgbottom[0] );
    dxf->write_group( 20,  pucsorgbottom[1] );
    dxf->write_group( 30,  pucsorgbottom[2] );

    dxf->write_group(  9, "$PUCSORGFRONT" );
    dxf->write_group( 10,  pucsorgfront[0] );
    dxf->write_group( 20,  pucsorgfront[1] );
    dxf->write_group( 30,  pucsorgfront[2] );

    dxf->write_group(  9, "$PUCSORGLEFT" );
    dxf->write_group( 10,  pucsorgleft[0] );
    dxf->write_group( 20,  pucsorgleft[1] );
    dxf->write_group( 30,  pucsorgleft[2] );

    dxf->write_group(  9, "$PUCSORGRIGHT" );
    dxf->write_group( 10,  pucsorgright[0] );
    dxf->write_group( 20,  pucsorgright[1] );
    dxf->write_group( 30,  pucsorgright[2] );

    dxf->write_group(  9, "$PUCSORGTOP" );
    dxf->write_group( 10,  pucsorgtop[0] );
    dxf->write_group( 20,  pucsorgtop[1] );
    dxf->write_group( 30,  pucsorgtop[2] );

    dxf->write_group(  9, "$PUCSORTHOREF" );
    dxf->write_group(  2,  pucsorthoref.c_str() );

    dxf->write_group(  9, "$PUCSORTHOVIEW" );
    dxf->write_group( 70, pucsorthoview );

    dxf->write_group(  9, "$PUCSXDIR" );
    dxf->write_group( 10,  pucsxdir[0] );
    dxf->write_group( 20,  pucsxdir[1] );
    dxf->write_group( 30,  pucsxdir[2] );

    dxf->write_group(  9, "$PUCSYDIR" );
    dxf->write_group( 10,  pucsydir[0] );
    dxf->write_group( 20,  pucsydir[1] );
    dxf->write_group( 30,  pucsydir[2] );

    dxf->write_group(  9, "$UCSBASE" );
    dxf->write_group(  2,  ucsbase.c_str() );

    dxf->write_group(  9, "$UCSNAME" );
    dxf->write_group(  2,  ucsname.c_str() );

    dxf->write_group(  9, "$UCSORG" );
    dxf->write_group( 10,  ucsorg[0] );
    dxf->write_group( 20,  ucsorg[1] );
    dxf->write_group( 30,  ucsorg[2] );

    dxf->write_group(  9, "$UCSORGBACK" );
    dxf->write_group( 10,  ucsorgback[0] );
    dxf->write_group( 20,  ucsorgback[1] );
    dxf->write_group( 30,  ucsorgback[2] );

    dxf->write_group(  9, "$UCSORGBOTTOM" );
    dxf->write_group( 10,  ucsorgbottom[0] );
    dxf->write_group( 20,  ucsorgbottom[1] );
    dxf->write_group( 30,  ucsorgbottom[2] );

    dxf->write_group(  9, "$UCSORGFRONT" );
    dxf->write_group( 10,  ucsorgfront[0] );
    dxf->write_group( 20,  ucsorgfront[1] );
    dxf->write_group( 30,  ucsorgfront[2] );

    dxf->write_group(  9, "$UCSORGLEFT" );
    dxf->write_group( 10,  ucsorgleft[0] );
    dxf->write_group( 20,  ucsorgleft[1] );
    dxf->write_group( 30,  ucsorgleft[2] );

    dxf->write_group(  9, "$UCSORGRIGHT" );
    dxf->write_group( 10,  ucsorgright[0] );
    dxf->write_group( 20,  ucsorgright[1] );
    dxf->write_group( 30,  ucsorgright[2] );

    dxf->write_group(  9, "$UCSORGTOP" );
    dxf->write_group( 10,  ucsorgtop[0] );
    dxf->write_group( 20,  ucsorgtop[1] );
    dxf->write_group( 30,  ucsorgtop[2] );

    dxf->write_group(  9, "$UCSORTHOREF" );
    dxf->write_group(  2,  ucsorthoref.c_str() );

    dxf->write_group(  9, "$UCSORTHOVIEW" );
    dxf->write_group( 70, ucsorthoview );

    dxf->write_group(  9, "$UCSXDIR" );
    dxf->write_group( 10,  ucsxdir[0] );
    dxf->write_group( 20,  ucsxdir[1] );
    dxf->write_group( 30,  ucsxdir[2] );

    dxf->write_group(  9, "$UCSYDIR" );
    dxf->write_group( 10,  ucsydir[0] );
    dxf->write_group( 20,  ucsydir[1] );
    dxf->write_group( 30,  ucsydir[2] );

    dxf->write_group(  9, "$WORLDVIEW" );
    dxf->write_group( 70, worldview );

    dxf->write_group( 0, "ENDSEC" );
}



void MyDXFHeader::debug_print( std::ostream &os ) const
{
    os << "*** Section HEADER ****************************************\n";

    MyDXFFile::debug_print_format( os, "acadver", acadver );
    MyDXFFile::debug_print_format( os, "angbase", angbase );
    MyDXFFile::debug_print_format( os, "angdir", angdir );

    MyDXFFile::debug_print_format( os, "handseed", handseed );
    MyDXFFile::debug_print_format( os, "dimasz", dimasz );
    MyDXFFile::debug_print_format( os, "dimgap", dimgap );
    MyDXFFile::debug_print_format( os, "dimexo", dimexo );
    MyDXFFile::debug_print_format( os, "dimexe", dimexe );
    MyDXFFile::debug_print_format( os, "dimtxt", dimtxt );
    MyDXFFile::debug_print_format( os, "insunits", insunits );
    MyDXFFile::debug_print_format( os, "plimmax", plimmax );
    MyDXFFile::debug_print_format( os, "plimmin", plimmin );

    MyDXFFile::debug_print_format( os, "orthomode", orthomode );

    MyDXFFile::debug_print_format( os, "pucsbase", pucsbase );
    MyDXFFile::debug_print_format( os, "pucsname", pucsname );
    MyDXFFile::debug_print_format( os, "pucsorg", pucsorg );
    MyDXFFile::debug_print_format( os, "pucsorgback", pucsorgback );
    MyDXFFile::debug_print_format( os, "pucsorgbottom", pucsorgbottom );
    MyDXFFile::debug_print_format( os, "pucsorgfront", pucsorgfront );
    MyDXFFile::debug_print_format( os, "pucsorgleft", pucsorgleft );
    MyDXFFile::debug_print_format( os, "pucsorgright", pucsorgright );
    MyDXFFile::debug_print_format( os, "pucsorgtop", pucsorgtop );
    MyDXFFile::debug_print_format( os, "pucsorthoref", pucsorthoref );
    MyDXFFile::debug_print_format( os, "pucsorthoview", pucsorthoview );
    MyDXFFile::debug_print_format( os, "pucsxdir", pucsxdir );
    MyDXFFile::debug_print_format( os, "pucsydir", pucsydir );
    
    MyDXFFile::debug_print_format( os, "ucsbase ", ucsbase );
    MyDXFFile::debug_print_format( os, "ucsname", ucsname );
    MyDXFFile::debug_print_format( os, "ucsorg", ucsorg );
    MyDXFFile::debug_print_format( os, "ucsorgback", ucsorgback );
    MyDXFFile::debug_print_format( os, "ucsorgbottom", ucsorgbottom );
    MyDXFFile::debug_print_format( os, "ucsorgfront", ucsorgfront );
    MyDXFFile::debug_print_format( os, "ucsorgleft", ucsorgleft );
    MyDXFFile::debug_print_format( os, "ucsorgright", ucsorgright );
    MyDXFFile::debug_print_format( os, "ucsorgtop", ucsorgtop );
    MyDXFFile::debug_print_format( os, "ucsorthoref", ucsorthoref );
    MyDXFFile::debug_print_format( os, "ucsorthoview", ucsorthoview );
    MyDXFFile::debug_print_format( os, "ucsxdir", ucsxdir );
    MyDXFFile::debug_print_format( os, "ucsydir", ucsydir );

    MyDXFFile::debug_print_format( os, "worldview", worldview );

    os << "\n";
}
