// globals.cpp

#include "globals.h"

std::ofstream logfile;
bool debug = false;


int N_aperture_x=5;
int N_aperture_y=5;

double xdisplacementSG=-1e-3;
double ydisplacementSG=-0.7e-3;

double aperture_xdist=0.019;
double aperture_ydist=0.021;

void set_NX_NY(int NX, int NY) {
    N_aperture_x=NX;
    N_aperture_y=NY;
}

void set_delta_apertures(double dx, double dy) {
    aperture_xdist=dx;
    aperture_ydist=dy;
}

void set_displacement_SG(double dx, double dy) {
    xdisplacementSG=dx;
    ydisplacementSG=dy;
}