/*! \file mydxftables.hpp
 *  \brief DXF Tables
 */

/* Copyright (c) 2010-2012 Taneli Kalvas. All rights reserved.
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

#ifndef MY_DXF_TABLES_HPP
#define MY_DXF_TABLES_HPP 1


#include <vector>
#include <string>
#include <stdint.h>
#include "mydxffile.hpp"


/*! \brief DXF table entry
 */
class MyDXFTableEntry
{

    std::string _handle;           // 5 (for others) or 105 (for DIMSTYLE)
    std::string _handle_to_owner;  // 330

protected:

    /*! \brief Constructor.
     */
    MyDXFTableEntry();

    /*! \brief Process group not belonging to the child entry.
     */
    void process_group( class MyDXFFile *dxf );
	
    /*! \brief Write common groups.
     */
    void write_common( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Debug print common groups.
     */
    void debug_print_common( std::ostream &os ) const;
    
public:

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFTableEntry() {};

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr ) = 0;

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const = 0;
    
    friend std::ostream &operator<<( std::ostream &os, const MyDXFTableEntry &e );
};



/*! \brief DXF table entry for block record table.
 */
class MyDXFTableEntryBlockRecord : public MyDXFTableEntry
{

    std::string _name;             // 2
    int16_t     _units;            // 70
    int8_t      _explodability;    // 280
    int8_t      _scalability;      // 281
    std::string _handle_to_layout; // 340

public:

    /*! \brief Construct entry by reading from DXF file.
     */
    MyDXFTableEntryBlockRecord( class MyDXFFile *dxf );

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFTableEntryBlockRecord();

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;
};



/*! \brief DXF table entry for layer table.
 */
class MyDXFTableEntryLayer : public MyDXFTableEntry
{

    std::string _name;             // 2
    std::string _linetype;         // 6
    int16_t     _flags;            // 70
    int16_t     _color;            // 62 (negative if layer off)
    bool        _plotting;         // 290
    int8_t      _lineweight;       // 370

    std::string _handle_to_plot_style_name; // 390
    std::string _handle_to_material; // 347

public:

    /*! \brief Construct entry by reading from DXF file.
     */
    MyDXFTableEntryLayer( class MyDXFFile *dxf );

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFTableEntryLayer();

    /*! \brief Return layer name
     */
    const std::string &name( void ) const { return( _name ); }

    /*! \brief Return layer name
     */
    void set_name( const std::string &name ) { _name = name; }

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;
};



/*! \brief DXF table entry for vport table.
 */
class MyDXFTableEntryVport : public MyDXFTableEntry
{

    std::string _name;             /* 2, Viewport name */
    int16_t     _flags;            /* 70, Standard flag values (bit-coded values):
				    * 16 = If set, table entry is externally dependent on an xref
				    * 32 = If both this bit and bit 16 are set, the externally dependent xref has been successfully resolved
				    * 64 = If set, the table entry was referenced by at least one entity in the drawing the last time
				    * the drawing was edited. (This flag is for the benefit of AutoCAD commands. It can be ignored
				    * by most programs that read DXF files and does not need to be set by programs that write DXF
				    * files) */
    Vec3D       _vmin;             /* 10, 20, Lower-left corner of viewport */
    Vec3D       _vmax;             /* 11, 21, Upper-right corner of viewport */
    Vec3D       _vcenter;          /* 12, 22, View center point (in DCS) */
    Vec3D       _snap_base;        /* 13, 23, Snap base point (in DCS) */
    Vec3D       _snap_spacing;     /* 14, 24, Snap spacing X and Y */
    Vec3D       _grid_spacing;     /* 15, 25, Grid spacing X and Y */
    Vec3D       _view_direction;   /* 16, 26, 36, View direction from target point (in WCS) */
    Vec3D       _view_target;      /* 17, 27, 37, View target point (in WCS) */
    
    Vec3D       _ucs_origin;       /* 110, 120, 130, UCS origin */
    Vec3D       _ucs_x;            /* 111, 121, 131, UCS x-axis */
    Vec3D       _ucs_y;            /* 112, 122, 132, UCS y-axis */
    int16_t     _ucs_type;         /* 79, Orthographic type of UCS
				    * 0 = UCS is not orthographic
				    * 1 = Top; 2 = Bottom
				    * 3 = Front; 4 = Back
				    * 5 = Left; 6 = Right */
    double      _elevation;        /* 146, Elevation */

public:

    /*! \brief Construct entry by reading from DXF file.
     */
    MyDXFTableEntryVport( class MyDXFFile *dxf );

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFTableEntryVport();

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;
};



/*! \brief DXF table class.
 */
class MyDXFTable
{
    std::string               _name;            // 2
    std::string               _handle;          // 5
    std::string               _handle_to_owner; // 330

    std::vector<MyDXFTableEntry *> _entries;         // 70 (size)

public:

    MyDXFTable( const std::string &name, class MyDXFFile *dxf );
    ~MyDXFTable();

    /*! \brief Return number of table entries.
     */
    uint32_t size() const { return( _entries.size() ); }

    /*! \brief Return const pointer to entry \a a.
     */
    const MyDXFTableEntry *get_entry( uint32_t a ) const { return( _entries[a] ); }

    /*! \brief Return pointer to entry \a a.
     */
    MyDXFTableEntry *get_entry( uint32_t a ) { return( _entries[a] ); }

    /*! \brief Write dxf file to stream.
     */
    void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;
};



/*! \brief DXF tables class.
 *
 *  Container for data of a DXF file tables.
 */
class MyDXFTables
{

    MyDXFTable  *_blockrecord;
    MyDXFTable  *_layer;
    MyDXFTable  *_vport;

public:

    MyDXFTables( class MyDXFFile *dxf );
    ~MyDXFTables();

    /*! \brief Return const pointer to layers table.
     */
    const MyDXFTable *get_layers_table( void ) const { return( _layer ); }

    /*! \brief Return pointer to layers table.
     */
    MyDXFTable *get_layers_table( void ) { return( _layer ); }

    /*! \brief Write dxf file to stream.
     */
    void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;
};


#endif
