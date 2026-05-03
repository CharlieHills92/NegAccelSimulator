#ifndef CROSS_SECTIONS_H_
#define CROSS_SECTIONS_H_

#include "particledatabase.hpp"
//~ #include "src/pier/matlabfun.h"

using namespace std;


// Energy in eV, sigma out in m2
double stripping_cross_at_E( double energy, double M_IONS, double & sigma_single, double & sigma_double);
double cs_HM_single_strip( double energy, double M_IONS );
double cs_HM_double_strip( double energy, double M_IONS );
// H-+H2->H-+H2++e OR H0+H2->H0+H2++e
double cs_bkg_ionization( double energy, double M_IONS );

// H0+H2->H++H2+e ionization of projectile H0
double  cs_proj_ionization_H0( double energy, double M_IONS );


// H++H2->H0+H2+ CX of projectile
double cs_CX_Hp( double energy, double M_IONS );

#endif