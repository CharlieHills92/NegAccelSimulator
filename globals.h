#ifndef GLOBALS_H
#define GLOBALS_H

#include <fstream>

extern std::ofstream logfile;
extern bool debug;
extern int N_aperture_x;
extern int N_aperture_y;
extern double xdisplacementSG;
extern double ydisplacementSG;
extern double aperture_xdist;
extern double aperture_ydist;

void set_NX_NY(int NX, int NY);
void set_displacement_SG(double dx, double dy);
void set_delta_apertures(double dx, double dy);

#endif // GLOBALS_H
