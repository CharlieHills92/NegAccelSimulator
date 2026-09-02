#ifndef FUNCT_H_
#define FUNCT_H_

#include "ibsimu.hpp"
#include <vector>
// #include "particledatabase.hpp"
//~ #include "src/pier/matlabfun.h"

using namespace std;


enum particle_family {
    PARTICLE_FAMILY_H = 0,
    PARTICLE_FAMILY_D = 1
};


enum particle_kind {
    PARTICLE_ALL=-100,
    PARTICLE_WRONG=-1,
    PARTICLE_HM = 0,
    PARTICLE_H0 = 1,
    PARTICLE_HP = 2,
    PARTICLE_H2P = 3,
    PARTICLE_H20 = 4,
    PARTICLE_H3P = 5,
    PARTICLE_E = 6,
    PARTICLE_NEGATIVE_ION = PARTICLE_HM,
    PARTICLE_NEUTRAL_ATOM = PARTICLE_H0,
    PARTICLE_POSITIVE_ION = PARTICLE_HP,
    PARTICLE_MOLECULAR_POSITIVE_ION = PARTICLE_H2P,
    PARTICLE_MOLECULAR_NEUTRAL = PARTICLE_H20,
    PARTICLE_TRIATOMIC_POSITIVE_ION = PARTICLE_H3P,
    PARTICLE_ELECTRON = PARTICLE_E
};

string get_particle_name(particle_kind pk);
int get_particle_int(particle_kind pk);

particle_kind int2kind(int num);
size_t particle_kind_count();

particle_family infer_particle_family(double ion_mass_u);
void set_active_particle_family(particle_family family);
particle_family get_active_particle_family();

particle_kind identify_particle_species( double mass, double charge, double ION_MASS );
particle_kind particle_kind_from_config_name(const std::string& kind);
double particle_kind_charge_state(particle_kind kind);
double particle_kind_mass_u(particle_kind kind, double ion_mass_u);
bool particle_kind_is_electron(particle_kind kind);
bool particle_kind_is_positive_ion(particle_kind kind);

/*! \brief Is this an un-extracted primary negative ion heading back toward the source?
 *
 *  Primary beam particles carry gen == 0; anything produced by a stripping reaction or by
 *  surface backscatter has gen != 0 (surface-generated particles carry
 *  gen >= SURFACE_GENERATION_OFFSET = 101). A gen-0 negative ion still travelling in -z
 *  never got extracted, so the electrode it lands on is registering an artifact of the
 *  plasma/meniscus model rather than a beam load.
 *
 *  Backscattered or stripped negative ions travelling upstream ARE a real load and are
 *  deliberately not matched here -- that is the whole point of testing generation rather
 *  than just the sign of vz.
 */
bool is_unextracted_primary_negative_ion(double mass, double charge, double vz, int gen,
                                         double ION_MASS);
double density_at_z( double zc, double pressure, vector<double>& pos, vector<double>& dens, bool & ciaone );
void load_density_profile( string filename, vector<double>& pos, vector<double>& dens );
/// 1 D interpolation of a float in a interval of 2 float
double interp1(double x1,double x2, double y1, double y2, double xi);
/// 1 D interpolation of  a float wrt 2 vectors yi=interp1(x,y,xi);
double _interp1(vector<double> x, vector<double> y, double xi, bool & found);



#endif