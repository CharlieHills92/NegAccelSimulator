/*
 * NIO1_geom.h
 *
 *  Created on: Oct 18, 2021
 *      Author: carlo
 */

#ifndef _NIO1_GEOM
#define _NIO1_GEOM

bool above_line( double r, double z, double R1, double R2, double z1, double z2 )
{
	return( r-R1 >= (R2-R1)/(z2-z1)*(z-z1) );
}

// PG: boolean definition
bool solid0_NIO1( double x, double y, double z )
{
    double r = sqrt(x*x+y*y);
    bool boolean;

    if( z < 3e-3 && z > 8.7e-3 )
    {
		boolean = false;
	}
	else if( z >= 3e-3 && z < 6.9e-3 )
	{
		boolean = above_line( r, z, 6.6e-3, 3.8e-3, 3e-3, 6.9e-3 );
	}
	else if( z >= 6.9e-3 && z <= 8.7e-3 )
	{
		boolean = above_line( r, z, 3.8e-3, 6.3e-3, 6.9e-3, 8.46e-3 );
	}
	//else if( z >= 8.46e-3 && z <= 8.7e-3 )
	//{
		//boolean = above_line( r, z, 6.3e-3, 6.4e-3, 8.46e-3, 8.7e-3 );
	//}
	else
	{
		boolean = false;
	}

    return boolean;
}

// EG: boolean definition
bool solid1_NIO1( double x, double y, double z )
{
    double r = sqrt(x*x+y*y);
    bool boolean;

    if( z < 13.7e-3 && z > 24e-3 )
    {
		boolean = false;
	}
	else if( z >= 13.7e-3 && z < 18.5e-3 )
	{
		boolean = above_line( r, z, 3.5e-3, 3.5e-3, 13.7e-3, 18.5e-3 );
	}
	else if( z >= 18.5e-3 && z <= 24e-3 )
	{
		boolean = above_line( r, z, 3.5e-3, 5e-3, 18.5e-3, 24e-3 );
	}
	else
	{
		boolean = false;
	}

    return boolean;
}

// GG: boolean definition
bool solid2_NIO1( double x, double y, double z )
{
    double r = sqrt(x*x+y*y);
    bool boolean;

    if( z < 48.6-3 && z > 59e-3 )
    {
		boolean = false;
	}
	else if( z >= 48.6e-3 && z < 57.5e-3 )
	{
		boolean = above_line( r, z, 3.5e-3, 3.5e-3, 48.6e-3, 57.5e-3 );
	}
	else if( z >= 57.5e-3 && z <= 59e-3 )
	{
		boolean = above_line( r, z, 3.5e-3, 4.4e-3, 57.5e-3, 59e-3 );
	}
	else
	{
		boolean = false;
	}

    return boolean;
}

// REP: boolean definition
bool solid3_NIO1( double x, double y, double z )
{
    double r = sqrt(x*x+y*y);
    bool boolean;

    if( z < 63-3 && z > 67e-3 )
    {
		boolean = false;
	}
	else if( z >= 63e-3 && z <= 67e-3 )
	{
		boolean = above_line( r, z, 4.4e-3, 4.4e-3, 63e-3, 67e-3 );
	}
	else
	{
		boolean = false;
	}

    return boolean;
}





#endif
