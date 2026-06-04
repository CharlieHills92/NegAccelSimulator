/*! \file fonts.cpp
 *  \brief %Font handling
 */

/* Copyright (c) 2005-2009,2012,2016 Taneli Kalvas. All rights reserved.
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

#include "fonts.hpp"
#include <fontconfig/fcfreetype.h>
#include <cairo-ft.h>
#include <iostream>
#include "error.hpp"
#include "symbols.hpp"



// Global fontlib object
FontLib fontlib;



//#define CAIROPLOT_FONTS_DEBUG 1


std::string wstring_to_utf8( const std::wstring &in )
{
    std::string out;

    for( unsigned int a = 0; a < in.length(); a++ ) {
	
	unsigned long c = in[a];
	if( c < 0x0080 ) {
	    // 1 byte utf8 representation
	    out.append( 1, (char)c );
	} else if( c < 0x0800 ) {
	    // 2 byte utf8 representation
	    out.append( 1, (char)(0xC0 | ((c>>6) & 0x1F)) );
	    out.append( 1, (char)(0x80 | (c & 0x3F)) );
	} else if( c < 0x10000 ) {
	    // 3 byte utf8 representation
	    out.append( 1, (char)(0xE0 | ((c>>12) & 0x0F)) );
	    out.append( 1, (char)(0x80 | ((c>>6) & 0x3F)) );
	    out.append( 1, (char)(0x80 | (c & 0x3F)) );
	} else {
	    // 4 byte utf8 representation
	    out.append( 1, (char)(0xF0 | ((c>>18) & 0x07)) );
	    out.append( 1, (char)(0x80 | ((c>>12) & 0x3F)) );
	    out.append( 1, (char)(0x80 | ((c>>6) & 0x3F)) );
	    out.append( 1, (char)(0x80 | (c & 0x3F)) );

	}
    }

    return( out );
}


/*! \brief Returns length of utf8 character.
 */
int utf8_get_length( unsigned char c )
{
    if( (c & 0x80) == 0x00 )
	return( 1 );
    else if( (c & 0xE0) == 0xC0 )
	return( 2 );
    else if( (c & 0xF0) == 0xE0 )
	return( 3 );
    else if( (c & 0xF8) == 0xF0 )
	return( 4 );
    else if( (c & 0xC0) == 0x80 )
	throw( Error( ERROR_LOCATION, "unexpected character continuation byte" ) );
    else
	throw( Error( ERROR_LOCATION, "invalid utf-8 code" ) );
}


/*! \brief Returns next utf8 character from string pointed to by iterator.
 */
std::string utf8_next_character( std::string::const_iterator it )
{
    std::string character;

    int len = utf8_get_length( *it );
    while( len-- > 0 ) {
	character += *it;
	it++;
    }

    return( character );
}


/* *************************************************************** *
 * FONT
 * *************************************************************** */


Font::Font( FcPattern *pat )
    : _cff(NULL), _pmatch(pat), _slant(CAIRO_FONT_SLANT_NORMAL),
      _weight(CAIRO_FONT_WEIGHT_NORMAL)
{
    FcChar8 *family;
    FcPatternGetString( pat, "family", 0, &family );
    _family = to_string(family);
    _cff = cairo_ft_font_face_create_for_pattern( pat );
}


Font::Font( const std::string &family )
    : _cff(NULL), _pmatch(NULL), _family(family) ,_slant(CAIRO_FONT_SLANT_NORMAL),
      _weight(CAIRO_FONT_WEIGHT_NORMAL)
{

}


Font::Font( const std::string &family, cairo_font_slant_t slant, 
	    cairo_font_weight_t weight )
    : _cff(NULL), _pmatch(NULL), _family(family) ,_slant(slant), _weight(weight)
{

}

Font::Font( const Font &font )
    : _cff(font._cff), _pmatch(font._pmatch), _family(font._family) ,_slant(font._slant),
      _weight(font._weight)
{
    cairo_font_face_reference( _cff );
}


Font &Font::operator=( const Font &font )
{
    if( _cff )
	cairo_font_face_destroy( _cff );

    _cff = font._cff;
    cairo_font_face_reference( _cff );

    _pmatch = font._pmatch;
    _family = font._family;
    _slant = font._slant;
    _weight = font._weight;

    return( *this );
}


Font::~Font()
{
    if( _cff )
	cairo_font_face_destroy( _cff );
}


std::string Font::family( void ) const
{
    return( _family );
}


cairo_font_slant_t Font::slant( void ) const
{
    return( _slant );
}


cairo_font_weight_t Font::weight( void ) const
{
    return( _weight );
}


cairo_font_face_t *Font::font_face( void ) const
{
    return( _cff );
}


FcPattern *Font::fcpattern( void ) const
{
    return( _pmatch );
}


void Font::load( FontLib &fontlib )
{
    FcPattern *preq;
    FcPattern *pmatch;
    FcResult result = FcResultMatch;

    /* Create request pattern */
#ifdef CAIROPLOT_FONTS_DEBUG
    std::cout << "  loading " << _family;
#endif
    preq = FcPatternCreate();
    FcPatternAddString( preq, FC_FAMILY, (FcChar8 *)_family.c_str() );

    /* Weight */
    if( _weight == CAIRO_FONT_WEIGHT_NORMAL ) {
#ifdef CAIROPLOT_FONTS_DEBUG
	std::cout << "Normal ";
#endif
	FcPatternAddInteger( preq, FC_WEIGHT, FC_WEIGHT_NORMAL );
    } else if( _weight == CAIRO_FONT_WEIGHT_BOLD ) {
#ifdef CAIROPLOT_FONTS_DEBUG
	std::cout << "Bold ";
#endif
	FcPatternAddInteger( preq, FC_WEIGHT, FC_WEIGHT_BOLD );
    } else {
#ifdef CAIROPLOT_FONTS_DEBUG
	std::cout << "Unknown ";
#endif
    }
/* Slant */
    if( _slant == CAIRO_FONT_SLANT_NORMAL ) {
#ifdef CAIROPLOT_FONTS_DEBUG
	std::cout << "Roman ";
#endif
	FcPatternAddInteger( preq, FC_SLANT, FC_SLANT_ROMAN );
    } else if( _slant == CAIRO_FONT_SLANT_ITALIC ) {
#ifdef CAIROPLOT_FONTS_DEBUG
	std::cout << "Italic ";
#endif
	FcPatternAddInteger( preq, FC_SLANT, FC_SLANT_ITALIC );
    } else if( _slant == CAIRO_FONT_SLANT_OBLIQUE ) {
#ifdef CAIROPLOT_FONTS_DEBUG
	std::cout << "Oblique ";
#endif
	FcPatternAddInteger( preq, FC_SLANT, FC_SLANT_OBLIQUE );
    } else {
#ifdef CAIROPLOT_FONTS_DEBUG
	std::cout << "Unknown ";
#endif
    }
#ifdef CAIROPLOT_FONTS_DEBUG
    std::cout << "\n";
#endif

    FcPatternAddBool( preq, FC_OUTLINE, FcTrue );

    /* Find best matching font */
    FcConfigSubstitute( fontlib.fc(), preq, FcMatchPattern );
    FcDefaultSubstitute( preq );
    pmatch = FcFontMatch( fontlib.fc(), preq, &result );
    if( result != FcResultMatch )
	throw( Error( ERROR_LOCATION, "FcFontMatch failed" ) );

    FcPatternDestroy( preq );
    /* Embolden-tag gets processed automatically via cairo.
     * Matrix to italicize needs to be fed externally, pmatch is
     * saved to provide this information */
    _pmatch = pmatch;
    _cff = cairo_ft_font_face_create_for_pattern( pmatch );

}


/* *************************************************************** *
 * FONTLIB
 * *************************************************************** */


FontLib::FontLib()
{
    //std::cout << "FontLib constructor\n";
    
    // Initialize libraries
    if( FT_Init_FreeType( &_ft ) )
	throw( Error( ERROR_LOCATION, "couldn\'t initialize FreeType" ) );
    if( !(_fc = FcInitLoadConfigAndFonts()) )
	throw( Error( ERROR_LOCATION, "couldn\'t initialize FontConfig" ) );

    // Add automatic search fonts
    push_auto_search_font( "Symbol" );

    // Load default font
    push_font( "Times", CAIRO_FONT_SLANT_NORMAL,
	       CAIRO_FONT_WEIGHT_NORMAL );
}


FontLib::~FontLib()
{
    //std::cout << "FontLib destructor\n";

    // Free fonts
    while( pop_auto_search_font() );
    while( pop_font() );
    
    FT_Done_FreeType( _ft );
}


std::string FontLib::family( void ) const
{
    if( _loaded.size() == 0 )
	return( std::string() );
    return( _loaded[_loaded.size()-1].family() );
}


cairo_font_slant_t FontLib::slant( void ) const
{
    if( _loaded.size() == 0 )
	return( CAIRO_FONT_SLANT_NORMAL );
    return( _loaded[_loaded.size()-1].slant() );
}


cairo_font_weight_t FontLib::weight( void ) const
{
    if( _loaded.size() == 0 )
	return( CAIRO_FONT_WEIGHT_NORMAL );
    return( _loaded[_loaded.size()-1].weight() );
}


cairo_font_face_t *FontLib::font_face( void ) const
{
    if( _loaded.size() == 0 )
	return( NULL );
    return( _loaded[_loaded.size()-1].font_face() );
}


FcPattern *FontLib::fcpattern( void ) const
{
    if( _loaded.size() == 0 )
	return( NULL );
    return( _loaded[_loaded.size()-1].fcpattern() );
}


void FontLib::push_auto_search_font( const std::string &family )
{
    _search.push_back( Font( family ) );
}


int FontLib::pop_auto_search_font( void )
{
    _search.pop_back();
    return( _search.size() );
}


void FontLib::push_font( FcPattern *pat )
{
    Font font( pat );
    font.load( *this );
    _loaded.push_back( font );
}


void FontLib::push_font( const std::string &family, cairo_font_slant_t slant, 
			 cairo_font_weight_t weight )

{
    Font font( family, slant, weight );
    font.load( *this );
    _loaded.push_back( font );
}


int FontLib::pop_font( void )
{
    _loaded.pop_back();
    return( _loaded.size() );
}


/* Get matrix for drawing text in matrix for curretly loaded font. 
 * Returns 1 if matrix is modified. Original matrix is saved to orig_matrix.
 */
int FontLib::get_font_matrix( cairo_t *cairo, cairo_matrix_t *matrix, cairo_matrix_t *orig_matrix )
{
    cairo_matrix_t matrix1;
    cairo_matrix_t matrix2;
    FcMatrix *fcmatrix2 = NULL;
    //FcResult result = FcResultMatch;
    //result = FcPatternGetMatrix( fcpattern(), "matrix", 0, &fcmatrix2 );

    FcPatternGetMatrix( fcpattern(), "matrix", 0, &fcmatrix2 );
    if( fcmatrix2 != NULL ) {
	/* positive slanting goes to wrong direction -> xy and yx terms
	 * are inverted. Is my fontconfig wrong or is this really supposed
	 * to be like this? */
	matrix2.xx = fcmatrix2->xx;
	matrix2.xy = -fcmatrix2->xy;
	matrix2.yx = -fcmatrix2->yx;
	matrix2.yy = fcmatrix2->yy;
    
	cairo_get_font_matrix( cairo, &matrix1 );
	cairo_get_font_matrix( cairo, orig_matrix );
	cairo_matrix_multiply( matrix, &matrix1, &matrix2 );
	cairo_set_font_matrix( cairo, matrix );
	return( 1 );
    }
  
    cairo_get_font_matrix( cairo, matrix );
    return( 0 );
}


void FontLib::combine_extents( cairo_text_extents_t *extents1, double x1, double y1,
			       const cairo_text_extents_t *extents2, double x2, double y2 )
{
    double t;

    extents1->x_advance = x2+extents2->x_advance-x1;
    extents1->y_advance = y2+extents2->y_advance-y1;

    if( x2+extents2->x_bearing < x1+extents1->x_bearing ) {
	t = x2-x1+extents2->x_bearing;
	extents1->width += (extents1->x_bearing-t);
	extents1->x_bearing = t;
    }
    if( y2+extents2->y_bearing < y1+extents1->y_bearing ) {
	t = y2-y1+extents2->y_bearing;
	extents1->height += (extents1->y_bearing-t);
	extents1->y_bearing = t;
    }
  
    if( x2+extents2->x_bearing+extents2->width >
	x1+extents1->x_bearing+extents1->width )
	extents1->width = x2-x1+extents2->x_bearing
	    -extents1->x_bearing+extents2->width;
    if( y2+extents2->y_bearing+extents2->height >
	y1+extents1->y_bearing+extents1->height )
	extents1->height = y2-y1+extents2->y_bearing
	    -extents1->y_bearing+extents2->height;
}


// Return advance in bytes
int FontLib::process_substr( cairo_t *cairo, const std::string &str, cairo_text_extents_t *extents,
			     double x0, double y0, double &x, double &y )
{
    FT_Face face;
    cairo_scaled_font_t *scaled_font;
    cairo_matrix_t orig_matrix;
    cairo_matrix_t matrix;
    cairo_matrix_t ctm;
    cairo_font_options_t *options;
  
    int matrix_modified;
    int len = str.length();
    double dx, dy;

#ifdef CAIROPLOT_FONTS_DEBUG
    std::cout << "\nprocess_substring()\n";
    std::cout << "  str = \'" << str << "\'\n";
    std::cout << "  len = " << len << "\n";
    std::cout << "  str = ";
    for( int a = 0; a < len; a++ )
	std::cout << (int)str[a] << " ";
    std::cout << "\n";
#endif
  
    // Get font face
    cairo_set_font_face( cairo, font_face() );
    cairo_get_matrix( cairo, &ctm );
    matrix_modified = get_font_matrix( cairo, &matrix, &orig_matrix );
    options = cairo_font_options_create();
    scaled_font = cairo_scaled_font_create( font_face(), &matrix, &ctm, options );
    face = cairo_ft_scaled_font_lock_face( scaled_font );
  
    FcChar32 c;
    int charlen = FcUtf8ToUcs4( (const FcChar8 *)str.c_str(), &c, len );
    if( charlen <= 0 )
	throw( Error( ERROR_LOCATION, "couldn\'t convert character from utf-8 to ucs-4" ) );
    
#ifdef CAIROPLOT_FONTS_DEBUG
    std::cout << "  c = " << c << "\n";
#endif

    cairo_glyph_t gindex;
    gindex.index = FT_Get_Char_Index( face, c );
    if( gindex.index == 0 ) {
	// Couldn't find character in the current font
	if( matrix_modified )
	    cairo_set_font_matrix( cairo, &orig_matrix );

#ifdef CAIROPLOT_FONTS_DEBUG
	std::cout << "  printing failed\n";
#endif

	return( 0 );
    }
    
#ifdef CAIROPLOT_FONTS_DEBUG
    std::cout << "  printing succeeded\n";
#endif

    if( FT_Load_Glyph( face, gindex.index, FT_LOAD_DEFAULT ) )
	throw( Error( ERROR_LOCATION, "couldn\'t load requested glyph" ) );

    // Set glyph position
    gindex.x = x;
    gindex.y = y;

    // Advance to next position
    dx = face->glyph->advance.x / 64.0;
    dy = -face->glyph->advance.y / 64.0;
    cairo_device_to_user_distance( cairo, &dx, &dy );

    // Do extents
    if( extents != NULL ) {
	cairo_text_extents_t extents2;
	cairo_glyph_extents( cairo, &gindex, 1, &extents2 );
	combine_extents( extents, x0, y0, &extents2, x, y );
    }

    // Advance coordinates
    x += dx;
    y += dy;

    // Unlock face
    cairo_ft_scaled_font_unlock_face( scaled_font );
    cairo_scaled_font_destroy( scaled_font );
    cairo_font_options_destroy( options );

    // Show glyphs
    if( extents == NULL )
	cairo_show_glyphs( cairo, &gindex, 1 );

    return( charlen );
}



void FontLib::process( cairo_t *cairo, const std::string &str, cairo_text_extents_t *extents, double &x, double &y )
{
    int step = 0;
    std::string::const_iterator it = str.begin();
    double x0 = x;
    double y0 = y;

#ifdef CAIROPLOT_FONTS_DEBUG
    std::cout << "\nprocess()\n";
    std::cout << "  str = \'" << str << "\'\n";
#endif

    // Process the string
    while( it != str.end() ) {

	// Get the next character
	std::string character = utf8_next_character( it );

	// Go with the primary font
	if( (step = process_substr( cairo, character, extents, x0, y0, x, y )) ) {
	    it += step;
	    continue;
	}

	// Try fonts on search stack for the problem character
	int a;
	for( a = _search.size()-1; a >= 0; a-- ) {

	    push_font( _search[a].family(), slant(), weight() );
	    step = process_substr( cairo, character, extents, x0, y0, x, y );
	    pop_font();

	    if( step )
		break;
	}
	if( step ) {
	    it += step;
	    continue;
	}

	// Search from all fonts  
	FcChar32   c;
	FcUtf8ToUcs4( (const FcChar8 *)str.c_str(), &c, str.length() );
	FcFontSet *fs = FcConfigGetFonts( _fc, FcSetSystem );
	for( int i = 0; i < fs->nfont; i++ ) {
	    FcPattern *pat = *(fs->fonts+i);
	    FcCharSet *cs;
	    FcPatternGetCharSet( pat, FC_CHARSET, 0, &cs );
	    FcBool outline;
	    FcPatternGetBool( pat, FC_OUTLINE, 0, &outline );
	    if( FcCharSetHasChar( cs, c ) && outline ) {
		
		push_font( pat );
		step = process_substr( cairo, character, extents, x0, y0, x, y );
		pop_font();

		if( step )
		    break;
	    }
	}

	// No suitable font found, skip character
	it += utf8_get_length( character[0] );
    }
}


void FontLib::text_extents( cairo_t *cairo, const std::string &str, cairo_text_extents_t *extents )
{
    if( extents == NULL )
	return;
    extents->x_bearing = 0.0;
    extents->y_bearing = 0.0;
    extents->width     = 0.0;
    extents->height    = 0.0;
    extents->x_advance = 0.0;
    extents->y_advance = 0.0;
    double x = 0.0;
    double y = 0.0;
    process( cairo, str, extents, x, y );
}


void FontLib::draw_text( cairo_t *cairo, const std::string &str, double &x, double &y )
{
    process( cairo, str, NULL, x, y );
}

