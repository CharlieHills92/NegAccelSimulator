#ifndef FUNCT_H_
#define FUNCT_H_

#include "ibsimu.hpp"
#include <vector>
// #include "particledatabase.hpp"
//~ #include "src/pier/matlabfun.h"

using namespace std;


enum particle_kind {
    PARTICLE_ALL=-100,
    PARTICLE_WRONG=-1,
    PARTICLE_HM = 0,
    PARTICLE_H0 = 1,
    PARTICLE_HP = 2,
    PARTICLE_H2P = 3,
    PARTICLE_H20 = 4,
    PARTICLE_E = 5
};

string get_particle_name(particle_kind pk);
int get_particle_int(particle_kind pk);

particle_kind int2kind(int num);

particle_kind identify_particle_species( double mass, double charge, double ION_MASS );
double density_at_z( double zc, double pressure, vector<double>& pos, vector<double>& dens, bool & ciaone );
void load_density_profile( string filename, vector<double>& pos, vector<double>& dens );
/// 1 D interpolation of a float in a interval of 2 float
double interp1(double x1,double x2, double y1, double y2, double xi);
/// 1 D interpolation of  a float wrt 2 vectors yi=interp1(x,y,xi);
double _interp1(vector<double> x, vector<double> y, double xi, bool & found);



#endif