#ifndef _GEOM_FUN
#define _GEOM_FUN

#include <iostream>
#include <vector>
#include <cmath>
#include "../globals.h"

inline bool inSolidFromHoleProfile(double z_point, double r_point, std::vector<double>& z_profile, std::vector<double>& r_profile) {
    // Define the z-r profile data
    // z vectors should be in increasing order
    // std::vector<double> z_profile = {0, 2.00E-03, 8.00E-03, 1.00E-02};  // z values (height)
    // std::vector<double> r_profile = {9.00E-3, 7.00E-03, 7.00E-03, 1.00E-02};  // r values (radial distances)

    // Find the corresponding radial distance for the given z using linear interpolation
    for (size_t i = 0; i < z_profile.size()-1; ++i) {
        // if (z < z_profile[0] || z > z_profile[z_profile.size()-1]) {
        //     return false;
        // }

        // if ( interpMethod[i] == 0 ) { // Linear interpolation case
        if (z_profile[i] != z_profile[i+1]) {
            if (z_point >= z_profile[i] && z_point <= z_profile[i+1]) {
                // Linear interpolation between two neighboring points in the z-direction
                double z0 = z_profile[i], z1 = z_profile[i+1];
                double r0 = r_profile[i], r1 = r_profile[i+1];
                
                // Interpolate to find the corresponding r at the given z
                double r_interpolated = r0 + (z_point - z0) * (r1 - r0) / (z1 - z0);
                
                // Check if the radial distance of the point is within the solid boundary
                return r_point >= r_interpolated;
            }
        // }
        // else if ( interpMethod[i] != 0 ) { // Curvature radius among two points
        //     double cradius = abs(interpMethod[i]);
        //     double deltaZ = z_profile[i+1]-z_profile[i];
        //     double deltaR = r_profile[i+1]-r_profile[i];
        //     // double m = -deltaZ/deltaR;
        //     // double zc = 0.5*(z_profile[i]+z_profile[i+1])+interpMethod[i]/cradius*sqrt((cradius*cradius-(deltaZ*deltaZ+deltaR*deltaR)/4)/(1+m*m));
        //     // double rc = m*zc+0.5*(r_profile[i]+r_profile[i+1])-0.5*m*(z_profile[i]+z_profile[i+1]);

        }
    }

    return false;  // If no match is found
}


inline bool inTorus(double z_point, double r_point, double rc, double zc, double minorradius) {
    double newr = r_point-rc;
    double newz = z_point-zc;
    return newr*newr + newz*newz <= minorradius*minorradius;
}

// Overload function adding the rounding of the start and end points with given radii
// It assumes:
// * for curvature of the start point: rounding center has coordinate (z_profile[0]+rc_srtat,r_profile[0]) and radius rc_start
// * for curvature of the end point: rounding center has coordinate (z_profile.back()-rc_end,r_profile.back()) and radius rc_end
inline bool inSolidFromHoleProfile(double z_point, double r_point, std::vector<double>& z_profile, std::vector<double>& r_profile, double rc_start, double rc_end) {
    bool insolid = inSolidFromHoleProfile(z_point, r_point, z_profile, r_profile);
    if (rc_start>0) {insolid = insolid || inTorus(z_point,r_point,r_profile[0],z_profile[0]+rc_start,rc_start);}
    if (rc_end>0) {insolid = insolid || inTorus(z_point,r_point,r_profile.back(),z_profile.back()-rc_end,rc_end);}

    return insolid;
}

inline bool createGrid(double x, double y, double z, std::vector<double>& z_profile, std::vector<double>& r_profile, double rc_start=0, double rc_end=0) {
    // Check if the point's z value is within the profile range
    if (z < z_profile.front() || z > z_profile.back()) return false;  // Outside the profile range in the z-direction
    if (abs(x)>N_aperture_x*aperture_xdist/2+0.001 || abs(y)>N_aperture_y*aperture_ydist/2+0.001) return true; // Outside the profile range in the  x or y direction, so it is inside the solid

    // int additionX = (N_aperture_x+1) % 2;
    // int additionY = (N_aperture_y+1) % 2;
    int additionX=0;
    int additionY=0;

    if (abs(x)<N_aperture_x*aperture_xdist/2 && abs(y)<N_aperture_y*aperture_ydist/2) {
        x = x-(round(x/aperture_xdist)+additionX/2)*aperture_xdist;
        y = y-(round(y/aperture_ydist)+additionY/2)*aperture_ydist;
    }

    double r = sqrt(x*x+y*y);

    return inSolidFromHoleProfile(z,r,z_profile,r_profile,rc_start,rc_end);

}

#endif