#ifndef _MTF_GEOM
#define _MTF_GEOM

// #include "../param/TEST_parameters.h"
#include "../globals.h"
#include "geom_function.h"

////////////////// Geometry of the electrodes  /////////////////////////

/// MTF 2021 Accelerator, QST 
// C. Poggi, IO,  12/2024. File created from MITICA_geom.h file by P. Veltri

// int N_aperture_x=1;
// int N_aperture_y=1;

// PG
bool solid0_MTF( double x, double y, double z )
{
  std::vector<double> z_profile = {3e-3,7e-3,9e-3};  // z values (height)
  std::vector<double> r_profile = {10e-3,7e-3,8.5e-3};  // r values (radial distances)

  // Check if the point's z value is within the profile range
  if (z < z_profile.front() || z > z_profile.back()) return false;  // Outside the profile range in the z-direction
  if (abs(x)>N_aperture_x*aperture_xdist/2+0.001 || abs(y)>N_aperture_y*aperture_ydist/2+0.001) return true; // Outside the profile range in the  x or y direction, so it is inside the solid

  // int additionX = (N_aperture_x+1) % 2;
  // int additionY = (N_aperture_y+1) % 2;
  int additionX=0;
  int additionY=0;


  if (x*x+y*y>1e-4 && abs(z-0.006)<0.001) 
    return true;


  if (abs(x)<N_aperture_x*aperture_xdist/2 && abs(y)<N_aperture_y*aperture_ydist/2) {
    x = x-(round(x/aperture_xdist)+additionX/2)*aperture_xdist;
    y = y-(round(y/aperture_ydist)+additionY/2)*aperture_ydist;
  }

  double r = sqrt(x*x+y*y);

  bool insolid = inSolidFromHoleProfile(z,r,z_profile,r_profile);

  //bool stopexitfromotherholes = false;
  // additional boolean for rounding
  // insolid = insolid || inTorus(z,r,2e-3,9e-3,2e-3);

  return insolid;
}

// EG
bool solid1_MTF( double x, double y, double z )
{
  std::vector<double> z_profileEG = {15e-3,15.5e-3,22e-3,28e-3,28.5e-3};  // z values (height)
  std::vector<double> r_profileEG = {7e-3,6.5e-3,6.5e-3,7.5e-3,7.5e-3};  // r values (radial distances)
  std::vector<double> z_profileSG = {28.5e-3,30e-3,31.5e-3};  // z values (height)
  std::vector<double> r_profileSG = {8e-3,8e-3,9e-3};  // r values (radial distances)

  // Check if the point's z value is within the profile range
  if (z < z_profileEG.front() || z > z_profileSG.back()) return false;  // Outside the profile range in the z-direction
  //if (abs(x)>N_aperture_x*aperture_xdist/2+0.001 || abs(y)>N_aperture_y*aperture_ydist/2+0.001) return true; // Outside the profile range in the  x or y direction, so it is inside the solid

  // int additionX = (N_aperture_x+1) % 2;
  // int additionY = (N_aperture_y+1) % 2;
  int additionX=0;
  int additionY=0;

  double xORI = x;
  double yORI = y;

  bool insolid = false;
  double r = 0;

  if (z >= z_profileEG.front() && z <= z_profileEG.back()) {

    if (abs(xORI)<N_aperture_x*aperture_xdist/2 && abs(yORI)<N_aperture_y*aperture_ydist/2) {
      x = xORI-(round(xORI/aperture_xdist)+additionX/2)*aperture_xdist;
      y = yORI-(round(yORI/aperture_ydist)+additionY/2)*aperture_ydist;
    }

    r = sqrt(x*x+y*y);

    insolid = inSolidFromHoleProfile(z,r,z_profileEG,r_profileEG);
  
  }

  else if (z > z_profileSG.front() && z <= z_profileSG.back()) {
    // nrow determination - fix for proper alternating pattern with centered grid
    int nrow = int(floor((y + aperture_ydist/2) / aperture_ydist));
    
    y = yORI-ydisplacementSG; // vertical displacement of SG
    x = (nrow % 2 == 0) ? xORI+xdisplacementSG : xORI-xdisplacementSG; // horizontal displacement depending on the row
    if (abs(x)<N_aperture_x*aperture_xdist/2 && abs(y)<N_aperture_y*aperture_ydist/2) {
      x = x-(round(x/aperture_xdist)+additionX/2)*aperture_xdist;
      y = y-(round(y/aperture_ydist)+additionY/2)*aperture_ydist;
    }
    //logfile << "\tDisplacement grid created. Xdispl " << xdisplacementSG << " Ydispl " << ydisplacementSG << endl; 

    r = sqrt(x*x+y*y);

    insolid = inSolidFromHoleProfile(z,r,z_profileSG,r_profileSG);

    // additional boolean for rounding
    // insolid = insolid || inTorus(z,r,2e-3,9e-3,2e-3);
  }

  return insolid;
}

// AG1
bool solid2_MTF( double x, double y, double z )
{
  std::vector<double> z_profile = {119.5e-3,120.5e-3,126.5e-3,135e-3,135e-3,136.5e-3};  // z values (height)
  std::vector<double> r_profile = {8e-3,7e-3,7e-3,8.5e-3,7e-3,8.5e-3};  // r values (radial distances)

  // Check if the point's z value is within the profile range
  if (z < z_profile.front() || z > z_profile.back()) return false;  // Outside the profile range in the z-direction
  // if (abs(x)>N_aperture_x*aperture_xdist/2+0.001 || abs(y)>N_aperture_y*aperture_ydist/2+0.001) return true; // Outside the profile range in the  x or y direction, so it is inside the solid

  // int additionX = (N_aperture_x+1) % 2;
  // int additionY = (N_aperture_y+1) % 2;
  int additionX=0;
  int additionY=0;

  if (abs(x)<N_aperture_x*aperture_xdist/2 && abs(y)<N_aperture_y*aperture_ydist/2) {
    x = x-(round(x/aperture_xdist)+additionX/2)*aperture_xdist;
    y = y-(round(y/aperture_ydist)+additionY/2)*aperture_ydist;
  }

  double r = sqrt(x*x+y*y);

  bool insolid = inSolidFromHoleProfile(z,r,z_profile,r_profile);

  // additional boolean for rounding
  insolid = insolid || inTorus(z,r,z_profile[1],r_profile[1],1e-3);
  insolid = insolid || inTorus(z,r,z_profile[4],r_profile[5],1.5e-3);

  return insolid;
}


// AG2
bool solid3_MTF( double x, double y, double z )
{
  std::vector<double> z_profile = {224.5e-3,225.5e-3,231.5e-3,240e-3,240e-3,241.5e-3};  // z values (height)
  std::vector<double> r_profile = {8e-3,7e-3,7e-3,8.5e-3,7e-3,8.5e-3};  // r values (radial distances)

  // Check if the point's z value is within the profile range
  if (z < z_profile.front() || z > z_profile.back()) return false;  // Outside the profile range in the z-direction
  // if (abs(x)>N_aperture_x*aperture_xdist/2+0.001 || abs(y)>N_aperture_y*aperture_ydist/2+0.001) return true; // Outside the profile range in the  x or y direction, so it is inside the solid

  // int additionX = (N_aperture_x+1) % 2;
  // int additionY = (N_aperture_y+1) % 2;
  int additionX=0;
  int additionY=0;

  if (abs(x)<N_aperture_x*aperture_xdist/2 && abs(y)<N_aperture_y*aperture_ydist/2) {
    x = x-(round(x/aperture_xdist)+additionX/2)*aperture_xdist;
    y = y-(round(y/aperture_ydist)+additionY/2)*aperture_ydist;
  }

  double r = sqrt(x*x+y*y);

  bool insolid = inSolidFromHoleProfile(z,r,z_profile,r_profile);

  // additional boolean for rounding
  insolid = insolid || inTorus(z,r,z_profile[1],r_profile[1],1e-3);
  insolid = insolid || inTorus(z,r,z_profile[4],r_profile[5],1.5e-3);

  return insolid;
}

// AG3
bool solid4_MTF( double x, double y, double z )
{
  std::vector<double> z_profile = {329.5e-3,330.5e-3,336.5e-3,345e-3,345e-3,346.5e-3};  // z values (height)
  std::vector<double> r_profile = {9e-3,8e-3,8e-3,9e-3,8e-3,9.5e-3};  // r values (radial distances)

  // Check if the point's z value is within the profile range
  if (z < z_profile.front() || z > z_profile.back()) return false;  // Outside the profile range in the z-direction
  // if (abs(x)>N_aperture_x*aperture_xdist/2+0.001 || abs(y)>N_aperture_y*aperture_ydist/2+0.001) return true; // Outside the profile range in the  x or y direction, so it is inside the solid

  // int additionX = (N_aperture_x+1) % 2;
  // int additionY = (N_aperture_y+1) % 2;
  int additionX=0;
  int additionY=0;

  if (abs(x)<N_aperture_x*aperture_xdist/2 && abs(y)<N_aperture_y*aperture_ydist/2) {
    x = x-(round(x/aperture_xdist)+additionX/2)*aperture_xdist;
    y = y-(round(y/aperture_ydist)+additionY/2)*aperture_ydist;
  }

  double r = sqrt(x*x+y*y);

  bool insolid = inSolidFromHoleProfile(z,r,z_profile,r_profile);

  // additional boolean for rounding
  insolid = insolid || inTorus(z,r,z_profile[1],r_profile[1],1e-3);
  insolid = insolid || inTorus(z,r,z_profile[4],r_profile[5],1.5e-3);

  return insolid;
}

// AG4
bool solid5_MTF( double x, double y, double z )
{
  std::vector<double> z_profile = {434.5e-3,435.5e-3,441.5e-3,450e-3,450e-3,451.5e-3};  // z values (height)
  std::vector<double> r_profile = {9e-3,8e-3,8e-3,9e-3,8e-3,9.5e-3};  // r values (radial distances)

  // Check if the point's z value is within the profile range
  if (z < z_profile.front() || z > z_profile.back()) return false;  // Outside the profile range in the z-direction
  // if (abs(x)>N_aperture_x*aperture_xdist/2+0.001 || abs(y)>N_aperture_y*aperture_ydist/2+0.001) return true; // Outside the profile range in the  x or y direction, so it is inside the solid

  // int additionX = (N_aperture_x+1) % 2;
  // int additionY = (N_aperture_y+1) % 2;
  int additionX=0;
  int additionY=0;

  if (abs(x)<N_aperture_x*aperture_xdist/2 && abs(y)<N_aperture_y*aperture_ydist/2) {
    x = x-(round(x/aperture_xdist)+additionX/2)*aperture_xdist;
    y = y-(round(y/aperture_ydist)+additionY/2)*aperture_ydist;
  }

  double r = sqrt(x*x+y*y);

  bool insolid = inSolidFromHoleProfile(z,r,z_profile,r_profile);

  // additional boolean for rounding
  insolid = insolid || inTorus(z,r,z_profile[1],r_profile[1],1e-3);
  insolid = insolid || inTorus(z,r,z_profile[4],r_profile[5],1.5e-3);

  return insolid;
}

// GG
bool solid6_MTF( double x, double y, double z )
{
  std::vector<double> z_profile = {539.5e-3,540.5e-3,546.5e-3,555e-3,555e-3,556.5e-3};  // z values (height)
  std::vector<double> r_profile = {9e-3,8e-3,8e-3,9e-3,8e-3,9.5e-3};  // r values (radial distances)

  // Check if the point's z value is within the profile range
  if (z < z_profile.front() || z > z_profile.back()) return false;  // Outside the profile range in the z-direction
  // if (abs(x)>N_aperture_x*aperture_xdist/2+0.001 || abs(y)>N_aperture_y*aperture_ydist/2+0.001) return true; // Outside the profile range in the  x or y direction, so it is inside the solid

  // int additionX = (N_aperture_x+1) % 2;
  // int additionY = (N_aperture_y+1) % 2;
  int additionX=0;
  int additionY=0;

  if (abs(x)<N_aperture_x*aperture_xdist/2 && abs(y)<N_aperture_y*aperture_ydist/2) {
    x = x-(round(x/aperture_xdist)+additionX/2)*aperture_xdist;
    y = y-(round(y/aperture_ydist)+additionY/2)*aperture_ydist;
  }

  double r = sqrt(x*x+y*y);

  bool insolid = inSolidFromHoleProfile(z,r,z_profile,r_profile);

  // additional boolean for rounding
  insolid = insolid || inTorus(z,r,z_profile[1],r_profile[1],1e-3);
  insolid = insolid || inTorus(z,r,z_profile[4],r_profile[5],1.5e-3);

  return insolid;
}


////// End of Geometry 
////////////////////////////////////////////////////////////////////////////////

#endif


