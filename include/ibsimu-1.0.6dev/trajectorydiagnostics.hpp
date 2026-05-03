/*! \file trajectorydiagnostics.hpp
 *  \brief Trajectory diagnostics
 */

/* Copyright (c) 2005-2012,2016,2021 Taneli Kalvas. All rights reserved.
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

#ifndef TRAJECTORYDIAGNOSTICS_HPP
#define TRAJECTORYDIAGNOSTICS_HPP 1


#include <vector>
#include <limits>
#include "error.hpp"
#include "histogram.hpp"
#include "types.hpp"


/*! \brief Class for trajectory diagnostic data column.
 *
 *  Contains one specified (trajectory_diagnostic_e) type of
 *  diagnostic data for N trajectories as a
 *  vector. %TrajectoryDiagnosticColumns are used in class
 *  TrajectoryDiagnosticData to store different types of diagnostic
 *  data and type information.
 */
class TrajectoryDiagnosticColumn 
{

    trajectory_diagnostic_e _diag; /*!< \brief Type of diagnostic in data. */
    std::vector<double>     _data; /*!< \brief Vector of diagnostic data. */

public:

    /*! \brief Create new diagnostic data column with type \a diag.
     */
    TrajectoryDiagnosticColumn( trajectory_diagnostic_e diag ) 
	: _diag(diag) {}

    /*! \brief Add mirrored trajectory diagnostic data to the column.
     *
     *  Mirror data along plane \a axis = \a level. The mirrored data
     *  is added to the data column. This function is used to get a
     *  complete data set of a particle beam, of which only half (or
     *  quarter) has been simulated.
     */
    void mirror( coordinate_axis_e axis, double level );

    /*! \brief Add diagnostic data point to column.
     */
    void add_data( double x ) {
	_data.push_back( x );
    }

    /*! \brief Get a reference to diagnostic data vector.
     */
    std::vector<double> &data( void ) { return( _data ); }

    /*! \brief Get a const reference to diagnostic data vector.
     */
    const std::vector<double> &data( void ) const { return( _data ); }

    /*! \brief Get size of diagnostic data vector.
     */
    size_t size( void ) const { return( _data.size() ); }

    /*! \brief Get diagnostic type.
     */
    trajectory_diagnostic_e diagnostic( void ) const { return( _diag ); }

    /*! \brief Get const reference to data element \a i.
     */
    const double &operator()( size_t i ) const { return( _data[i] ); }

    /*! \brief Get reference to data element \a i.
     */
    double &operator()( size_t i ) { return( _data[i] ); }

    /*! \brief Get const reference to data element \a i.
     */
    const double &operator[]( size_t i ) const { return( _data[i] ); }

    /*! \brief Get reference to data element \a i.
     */
    double &operator[]( size_t i ) { return( _data[i] ); }
};


/*! \brief Class for trajectory diagnostic data.
 *
 *  Contains a vector of diagnostic columns (TrajectoryDiagnosticColumn).
 */
class TrajectoryDiagnosticData 
{

    std::vector<TrajectoryDiagnosticColumn> _column; /*!< \brief Vector of diagnostic data columns. */

public:

    /*! \brief Create new empty diagnostic data object.
     */
    TrajectoryDiagnosticData() {}

    /*! \brief Create diagnostic data object with diagnostic types defined in vector \a diag.
     */
    TrajectoryDiagnosticData( std::vector<trajectory_diagnostic_e> diag ) {
	for( size_t a = 0; a < diag.size(); a++ )
	    _column.push_back( TrajectoryDiagnosticColumn( diag[a] ) );
    }

    /*! \brief Mirror data columns along plane at \a axis = \a level.
     */
    void mirror( coordinate_axis_e axis, double level ) {
	for( size_t a = 0; a < _column.size(); a++ )
	    _column[a].mirror( axis, level );
    }

    /*! \brief Clear all data and diagnostic types.
     */
    void clear() {
	_column.clear();
    }
    
    /*! \brief Add data column with type \a diag.
     */
    void add_data_column( trajectory_diagnostic_e diag ) {
	_column.push_back( TrajectoryDiagnosticColumn( diag ) );
    }

    /*! \brief Return number of data columns.
     */
    size_t diag_size() const {
	return( _column.size() );
    }

    /*! \brief Return number of trajectories in data.
     */
    size_t traj_size() const {
	if( _column.size() > 0 )
	    return( _column[0].size() );
	return( 0 );
    }

    /*! \brief Return \a i:th diagnostic type.
     */
    trajectory_diagnostic_e diagnostic( size_t i ) const {
	return( _column[i].diagnostic() ); 
    }

    /*! \brief Return \a i:th diagnostic type.
     */
    const TrajectoryDiagnosticColumn &operator()( size_t i ) const {
        if( i >= _column.size() )
            throw( ErrorRange( ERROR_LOCATION, i, _column.size() ) );
	return( _column[i] );
    }

    /*! \brief Return \a i:th diagnostic column.
     */
    TrajectoryDiagnosticColumn &operator()( size_t i ) {
        if( i >= _column.size() )
            throw( ErrorRange( ERROR_LOCATION, i, _column.size() ) );
    	return( _column[i] );
    }

    /*! \brief Return const reference to \a j:th trajectory data in \a
     *  i:th diagnostic column.
     */
    const double &operator()( size_t j, size_t i ) const {
        if( i >= _column.size() )
            throw( ErrorRange( ERROR_LOCATION, i, _column.size() ) );
    	return( _column[i](j) );
    }

    /*! \brief Return reference to \a j:th trajectory data in \a i:th
     *  diagnostic column.
     */
    double &operator()( size_t j, size_t i ) {
        if( i >= _column.size() )
            throw( ErrorRange( ERROR_LOCATION, i, _column.size() ) );
    	return( _column[i](j) );
    }

    /*! \brief Add data point to \a i:th diagnostic column.
     */
    void add_data( size_t i, double x ) {
	_column[i].add_data( x );
    }

    /*! \brief Export trajectory data as ASCII.
     */
    void export_data( const std::string &filename );
};


/*! \brief Class for emittance statistics.
 *
 *  %Emittance class does a statistical analysis on the particle
 *  distribution and it calculates averages \f$ <x> \f$ and 
 *  \f$ <x'> \f$ and the expectation values \f$ <x^2> \f$, 
 *  \f$ <x'^2> \f$ and \f$ <x x'> \f$. From these it calculates 
 *  the rms-emittance
 *  \f[ \epsilon = \sqrt{ <x^2><x'^2> - <x x'>^2 } \f]
 *  and the Twiss parameters
 *  \f[ \alpha = \frac{-<x x'>}{\epsilon}, \beta = \frac{<x^2>}{\epsilon}, \gamma = \frac{<x'^2>}{\epsilon} \f]
 *  In addition to these physical values, the class calculates the angle of the ellipse
 *  \f[ \theta = \frac{1}{2} \arctan2{\left( -2\alpha, \beta - \gamma \right)} \f]
 *  and the half-axis lengths
 *  \f[ r_1 = \sqrt{\frac{\epsilon}{2}} ( \sqrt{H+1} + \sqrt{H-1} ) \f]
 *  \f[ r_2 = \sqrt{\frac{\epsilon}{2}} ( \sqrt{H+1} - \sqrt{H-1} ), \f]
 *  where
 *  \f[ H = \frac{\beta + \gamma}{2} \f]
 *  
 */
class Emittance
{
protected:

    double _Isum;

    double _xave;
    double _xpave;

    double _x2;
    double _xp2;
    double _xxp;

    double _alpha;
    double _beta;
    double _gamma;
    double _epsilon;

    double _angle;
    double _rmajor;
    double _rminor;

public:

    /*! \brief Default constructor for emittance statistics
     */
    Emittance();

    /*! \brief Constructor for emittance statistics from trajectory
     *  diagnostic data columns \a x, \a xp and current \a I.
     */
    Emittance( const std::vector<double> &x,
	       const std::vector<double> &xp,
	       const std::vector<double> &I );

    /*! \brief Constructor for emittance statistics from trajectory
     *  diagnostic data columns \a x, \a xp, assuming even weights.
     */
    Emittance( const std::vector<double> &x,
	       const std::vector<double> &xp );

    /*! \brief Constructor for emittance statistics from trajectory
     *  diagnostic data in mesh form
     *
     *  The mesh has integer dimensions of (\a xsize, \a xpsize) and
     *  has the extents defined by \a range, where \a range = (\a
     *  xmin, \a xpmin, \a xmax, \a xpmax). Current data at each mesh
     *  node is given by vector I, where data is stored in x major
     *  order (I[xindex+xpindex*xsize]).
     */
    Emittance( size_t xsize, size_t xpsize, const double range[4],
	       const std::vector<double> &I );

    /*! \brief Return current of beam trajectories used to calculate the emittance.
     *
     *  The current is only defined in some cases. If emittance
     *  constructure Emittance( x, xp ) is used, this functions
     *  returns the number of trajectories.
     */
    double current( void ) const { return( _Isum ); }

    /*! \brief Return average position (center location) of emittance distribution.
     */
    double xave( void ) const { return( _xave ); }

    /*! \brief Return average angle (center location) of emittance distribution.
     */
    double xpave( void ) const { return( _xpave ); }

    /*! \brief Return \f$\alpha\f$ of emittance distribution.
     */
    double alpha( void ) const { return( _alpha ); }

    /*! \brief Return \f$\beta\f$ of emittance distribution.
     */
    double beta( void ) const { return( _beta ); }

    /*! \brief Return \f$\gamma\f$ of emittance distribution.
     */
    double gamma( void ) const { return( _gamma ); }

    /*! \brief Return rms emittance.
     */
    double epsilon( void ) const { return( _epsilon ); }

    /*! \brief Return angle of fitted rms ellipse.
     */
    double angle( void ) const { return( _angle ); }

    /*! \brief Return major radius of fitted rms ellipse.
     */
    double rmajor( void ) const { return( _rmajor ); }

    /*! \brief Return minor radius of fitted rms ellipse.
     */
    double rminor( void ) const { return( _rminor ); } 

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;
};



/*! \brief Class for emittance conversion from (r,r') to (x,x').
 */
class EmittanceConv : public Emittance
{

    Histogram2D *_grid;

public:

    /*! \brief Constructor for \a (x,x') emittance data and statistics
     *  from \a (r,r') data.
     *
     *  %EmittanceConv class reads particle diagnostic data arrays for
     *  \a r (radius), \a rp (radial angle), \a ap (skew angle) and \a
     *  I (current) and builds \a (x,x') data in a grid array of size
     *  \a n by \a m. Here the skew angle is \f$ r\omega/v_z \f$,
     *  where \f$ v_z \f$ is the velocity to the direction of beam
     *  propagation. The conversion is done by rotating each
     *  trajectory diagnostic point around the axis in \a rotn steps
     *  (defaults to 100). The output grid size can be forced by
     *  setting \a (xmin,xpmin,xmax,xpmax) variables, otherwise the
     *  grid is autotomatically sized to fit all data.
     *
     *  The emittance statistics is built using original data and not
     *  the gridded data for maximized precision.
     */
    EmittanceConv( uint32_t n, uint32_t m,
		   const std::vector<double> &r,
		   const std::vector<double> &rp,
		   const std::vector<double> &ap,
		   const std::vector<double> &I,
		   uint32_t rotn = 100,
		   double xmin = std::numeric_limits<double>::quiet_NaN(), 
		   double xpmin = std::numeric_limits<double>::quiet_NaN(), 
		   double xmax = std::numeric_limits<double>::quiet_NaN(), 
		   double xpmax = std::numeric_limits<double>::quiet_NaN() );

    /*! \brief Destructor for emittance converter.
     */
    ~EmittanceConv();

    /*! \brief Get a const reference to histogram built.
     */
    const Histogram2D &histogram( void ) const { return( *_grid ); }

    /*! \brief Free emittance histogram.
     */
    void free_histogram( void ) { delete _grid; _grid = NULL; }
};


#endif
