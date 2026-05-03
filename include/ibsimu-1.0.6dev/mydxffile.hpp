/*! \file mydxffile.hpp
 *  \brief DXF File
 */

/* Copyright (c) 2010-2011 Taneli Kalvas. All rights reserved.
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

#ifndef MY_DXF_FILE_HPP
#define MY_DXF_FILE_HPP 1



//#define MYDXF_DEBUG 1
//#define MYDXF_DEBUG_PLOT 1



#include <fstream>
#include "mydxfheader.hpp"
#include "mydxftables.hpp"
#include "mydxfblocks.hpp"
#include "mydxfentities.hpp"




/*! \brief DXF file class.
 *
 *  This class is a memory representation of a dxf file read from the
 *  disc. The class can be used to read a dxf file. All supported
 *  features are saved to the hierarchy, all others are silently
 *  ignored.
 */
class MyDXFFile
{
    std::ifstream _istr;
    std::ofstream _ostr;
    bool          _ascii;
    int           _linec;

    int           _wlevel;

    int           _group_code;
    int           _group_type;

    std::string   _group_string;
    double        _group_double;
    bool          _group_bool;
    int8_t        _group_int8;
    int16_t       _group_int16;
    int32_t       _group_int32;
    int64_t       _group_int64;

    class MyDXFHeader    *_header;
    class MyDXFTables    *_tables;
    class MyDXFBlocks    *_blocks;
    class MyDXFEntities  *_entities;


public:
	
    /*! \brief Construct empty DXF file.
     */
    MyDXFFile();

    /*! \brief Construct by reading a DXF file.
     */
    MyDXFFile( const std::string &filename );

    /*! \brief Destructor.
     */
    ~MyDXFFile();

    /*! \brief Read DXF file.
     */
    void read( const std::string &filename );

    /*! \brief Write DXF file.
     */
    void write( const std::string &filename );

    /*! \brief Set the level of warning messages.
     *
     *  If \a wlevel is set to zero, no warnings will be printed. With
     *  increasing \a wlevel more warning messages are printed to
     *  standard output. With \a wlevel >= 1 problems handling
     *  entities are reported. With \a wlevel >= 2 all unsupported
     *  features are reported.
     */
    void set_warning_level( int wlevel ) { _wlevel = wlevel; }

    /*! \brief Get the level of warning messages.
     */
    int wlevel( void ) const { return( _wlevel ); }



    /*! \brief Write string group to output file.
     */
    void write_group( int code, const char *data );

    /*! \brief Write double group to output file.
     */
    void write_group( int code, double data );

    /*! \brief Write bool group to output file.
     */
    void write_group( int code, bool data );

    /*! \brief Write int8_t group to output file.
     */
    void write_group( int code, int8_t data );

    /*! \brief Write int16_t group to output file.
     */
    void write_group( int code, int16_t data );

    /*! \brief Write int32_t group to output file.
     */
    void write_group( int code, int32_t data );

    /*! \brief Write int64_t group to output file.
     */
    void write_group( int code, int64_t data );



    /*! \brief Read next group from open file and return group code.
     *
     *  Returns the group code read or -1 on EOF. An error is thrown
     *  on all other errors.
     */
    int read_group( void );

    /*! \brief Get code of the last group read.
     */
    int group_get_code( void ) const;

    /*! \brief Get the value of the last group read assuming it is a string.
     *
     *  An error is thrown if group type does not match.
     */
    std::string group_get_string( void ) const;

    /*! \brief Get the value of the last group read assuming it is a double.
     *
     *  An error is thrown if group type does not match.
     */
    double group_get_double( void ) const;

    /*! \brief Get the value of the last group read assuming it is a bool.
     *
     *  An error is thrown if group type does not match.
     */
    bool group_get_bool( void ) const;

    /*! \brief Get the value of the last group read assuming it is a int8.
     *
     *  An error is thrown if group type does not match.
     */
    int8_t group_get_int8( void ) const;

    /*! \brief Get the value of the last group read assuming it is a int16.
     *
     *  An error is thrown if group type does not match.
     */
    int16_t group_get_int16( void ) const;

    /*! \brief Get the value of the last group read assuming it is a int32.
     *
     *  An error is thrown if group type does not match.
     */
    int32_t group_get_int32( void ) const;

    /*! \brief Get the value of the last group read assuming it is a int64.
     *
     *  An error is thrown if group type does not match.
     */
    int64_t group_get_int64( void ) const;

    /*! \brief Get the current line number in DXF file during read.
     */
    int linec( void ) const { return( _linec ); }





    /*! \brief Get a pointer to the entities of DXF file.
     */
    class MyDXFEntities *get_entities( void ) { return( _entities ); };

    /*! \brief Get a const pointer to the entities of DXF file.
     */
    const class MyDXFEntities *get_entities( void ) const { return( _entities ); };



    /*! \brief Get a pointer to the blocks of DXF file.
     */
    class MyDXFBlocks *get_blocks( void ) { return( _blocks ); };

    /*! \brief Get a const pointer to the blocks of DXF file.
     */
    const class MyDXFBlocks *get_blocks( void ) const { return( _blocks ); };


    /*! \brief Get a pointer to the tables of DXF file.
     */
    class MyDXFTables *get_tables( void ) { return( _tables ); };

    /*! \brief Get a const pointer to the tables of DXF file.
     */
    const class MyDXFTables *get_tables( void ) const { return( _tables ); };




    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;

    /*! \brief Print debugging information to os with correct formatting
     */
    static void debug_print_format( std::ostream &os, 
				    const std::string &fieldname, 
				    const std::string &val );
    
    static void debug_print_format( std::ostream &os, 
				    const std::string &fieldname, 
				    double val );

    static void debug_print_format( std::ostream &os, 
				    const std::string &fieldname, 
				    int val );

    static void debug_print_format( std::ostream &os, 
				    const std::string &fieldname, 
				    const Vec3D &val );
};


#endif




