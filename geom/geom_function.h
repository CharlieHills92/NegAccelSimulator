#ifndef _GEOM_FUN
#define _GEOM_FUN

#include <iostream>
#include <vector>
#include <cmath>
#include "../globals.h"

inline bool inSolidFromHoleProfile(
    double z_point,
    double r_point,
    const std::vector<double>& z_profile,
    const std::vector<double>& r_profile) {
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


inline bool inCircleZR(double z_point, double r_point, double z_center, double r_center, double radius) {
    double newr = r_point-r_center;
    double newz = z_point-z_center;
    return newr*newr + newz*newz <= radius*radius;
}


inline bool inTorus(double z_point, double r_point, double rc, double zc, double minorradius) {
    return inCircleZR(z_point, r_point, zc, rc, minorradius);
}

inline double normalizeAngle(double angle) {
    while (angle <= -M_PI) {
        angle += 2.0 * M_PI;
    }
    while (angle > M_PI) {
        angle -= 2.0 * M_PI;
    }

    return angle;
}

inline double positiveAngleDistance(double from, double to) {
    double delta = normalizeAngle(to - from);
    if (delta < 0.0) {
        delta += 2.0 * M_PI;
    }
    return delta;
}

inline bool angleOnMinorArc(double theta, double start_angle, double sweep) {
    const double tolerance = 1.0e-9;
    if (sweep >= 0.0) {
        return positiveAngleDistance(start_angle, theta) <= sweep + tolerance;
    }
    return positiveAngleDistance(theta, start_angle) <= -sweep + tolerance;
}

inline bool segmentRoundingBoundaryRadius(
    double z_point,
    double z0,
    double r0,
    double z1,
    double r1,
    double rounding_value,
    double& boundary_radius) {
    const double radius = std::abs(rounding_value);
    if (radius <= 0.0) {
        return false;
    }

    const double dz = z1 - z0;
    const double dr = r1 - r0;
    const double chord = std::sqrt(dz * dz + dr * dr);
    if (chord <= 1.0e-12 || chord > 2.0 * radius + 1.0e-12) {
        return false;
    }

    const double mid_z = 0.5 * (z0 + z1);
    const double mid_r = 0.5 * (r0 + r1);
    const double half_chord = 0.5 * chord;
    const double center_offset = std::sqrt(std::max(0.0, radius * radius - half_chord * half_chord));
    const double perp_z = -dr / chord;
    const double perp_r = dz / chord;

    double chosen_center_z = 0.0;
    double chosen_center_r = 0.0;
    double chosen_start_angle = 0.0;
    double chosen_sweep = 0.0;
    double chosen_mid_r = 0.0;
    bool found_center = false;

    for (int sign = -1; sign <= 1; sign += 2) {
        const double center_z = mid_z + sign * center_offset * perp_z;
        const double center_r = mid_r + sign * center_offset * perp_r;
        const double start_angle = std::atan2(r0 - center_r, z0 - center_z);
        const double end_angle = std::atan2(r1 - center_r, z1 - center_z);
        const double sweep = normalizeAngle(end_angle - start_angle);
        const double mid_angle = start_angle + 0.5 * sweep;
        const double mid_r_value = center_r + radius * std::sin(mid_angle);

        if (!found_center ||
            (rounding_value > 0.0 && mid_r_value < chosen_mid_r) ||
            (rounding_value < 0.0 && mid_r_value > chosen_mid_r)) {
            chosen_center_z = center_z;
            chosen_center_r = center_r;
            chosen_start_angle = start_angle;
            chosen_sweep = sweep;
            chosen_mid_r = mid_r_value;
            found_center = true;
        }
    }

    if (!found_center) {
        return false;
    }

    const double normalized_z = (z_point - chosen_center_z) / radius;
    if (normalized_z < -1.0 - 1.0e-12 || normalized_z > 1.0 + 1.0e-12) {
        return false;
    }

    const double clamped_z = std::max(-1.0, std::min(1.0, normalized_z));
    const double base_angle = std::acos(clamped_z);
    const double candidate_angles[2] = {base_angle, -base_angle};
    bool found_angle = false;
    double selected_radius = 0.0;

    for (int index = 0; index < 2; ++index) {
        const double theta = candidate_angles[index];
        if (!angleOnMinorArc(theta, chosen_start_angle, chosen_sweep)) {
            continue;
        }
        const double candidate_radius = chosen_center_r + radius * std::sin(theta);
        if (!found_angle ||
            (rounding_value > 0.0 && candidate_radius < selected_radius) ||
            (rounding_value < 0.0 && candidate_radius > selected_radius)) {
            selected_radius = candidate_radius;
            found_angle = true;
        }
    }

    if (!found_angle) {
        return false;
    }

    boundary_radius = selected_radius;
    return true;
}

inline bool segmentRoundedProfileContains(
    double z_point,
    double r_point,
    const std::vector<double>& z_profile,
    const std::vector<double>& r_profile,
    const std::vector<double>& rounding_radii) {
    if (z_profile.size() < 2 || z_profile.size() != r_profile.size()) {
        return false;
    }

    for (size_t index = 1; index < z_profile.size(); ++index) {
        const double z0 = z_profile[index - 1];
        const double z1 = z_profile[index];
        if (std::abs(z1 - z0) <= 1.0e-15) {
            continue;
        }
        if (z_point < z0 || z_point > z1) {
            continue;
        }

        double boundary_radius = 0.0;
        const double r0 = r_profile[index - 1];
        const double r1 = r_profile[index];
        const double rounding_value = (rounding_radii.size() == z_profile.size()) ? rounding_radii[index] : 0.0;
        if (!segmentRoundingBoundaryRadius(z_point, z0, r0, z1, r1, rounding_value, boundary_radius)) {
            boundary_radius = r0 + (z_point - z0) * (r1 - r0) / (z1 - z0);
        }
        return r_point >= boundary_radius;
    }

    return false;
}

inline bool inSolidFromHoleProfile(
    double z_point,
    double r_point,
    const std::vector<double>& z_profile,
    const std::vector<double>& r_profile,
    const std::vector<double>& rounding_radii) {
    if (rounding_radii.empty()) {
        return inSolidFromHoleProfile(z_point, r_point, z_profile, r_profile);
    }

    return segmentRoundedProfileContains(z_point, r_point, z_profile, r_profile, rounding_radii);
}

// Overload function adding the rounding of the start and end points with given radii
// It assumes:
// * for curvature of the start point: rounding center has coordinate (z_profile[0]+rc_srtat,r_profile[0]) and radius rc_start
// * for curvature of the end point: rounding center has coordinate (z_profile.back()-rc_end,r_profile.back()) and radius rc_end
inline bool inSolidFromHoleProfile(
    double z_point,
    double r_point,
    const std::vector<double>& z_profile,
    const std::vector<double>& r_profile,
    double rc_start,
    double rc_end) {
    std::vector<double> rounding_radii(z_profile.size(), 0.0);
    if (rounding_radii.size() > 1) {
        rounding_radii[1] = rc_start;
        rounding_radii.back() = rc_end;
    }
    return inSolidFromHoleProfile(z_point, r_point, z_profile, r_profile, rounding_radii);
}

inline bool createGrid(
    double x,
    double y,
    double z,
    const std::vector<double>& z_profile,
    const std::vector<double>& r_profile,
    const std::vector<double>& rounding_radii) {
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

    return inSolidFromHoleProfile(z, r, z_profile, r_profile, rounding_radii);

}

inline bool createGrid(
    double x,
    double y,
    double z,
    const std::vector<double>& z_profile,
    const std::vector<double>& r_profile,
    double rc_start=0,
    double rc_end=0) {
    std::vector<double> rounding_radii(z_profile.size(), 0.0);
    if (rounding_radii.size() > 1) {
        rounding_radii[1] = rc_start;
        rounding_radii.back() = rc_end;
    }
    return createGrid(x, y, z, z_profile, r_profile, rounding_radii);

}

#endif