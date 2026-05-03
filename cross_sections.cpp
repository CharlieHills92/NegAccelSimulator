
#include "cross_sections.h"
#include <cmath>
//~ #include "src/pier/matlabfun.h"

using namespace std;


// Energy in eV, sigma out in m2
double stripping_cross_at_E( double energy, double M_IONS, double & sigma_single, double & sigma_double) {
	
	// Validate input parameters
	if (energy <= 0.0 || M_IONS <= 0.0 || !std::isfinite(energy) || !std::isfinite(M_IONS)) {
		sigma_single = 0.0;
		sigma_double = 0.0;
		return 0.0;
	}

	double single_strip=cs_HM_single_strip( energy, M_IONS );
	double double_strip=cs_HM_double_strip( energy, M_IONS );

	// Validate outputs
	if (!std::isfinite(single_strip) || single_strip < 0.0) single_strip = 0.0;
	if (!std::isfinite(double_strip) || double_strip < 0.0) double_strip = 0.0;

	sigma_single=single_strip;
	sigma_double=double_strip;

	
	return single_strip+double_strip;
}

double cs_HM_single_strip( double energy, double M_IONS ) {
	if (energy <= 0.0 || M_IONS <= 0.0 || !std::isfinite(energy) || !std::isfinite(M_IONS)) {
		return 0.0;
	}
	
	double E = energy*M_IONS;
	if (E <= 0.0 || !std::isfinite(E)) return 0.0;
	
	double lg10E=log10(E);
	if (!std::isfinite(lg10E)) return 0.0;
	
	double single_strip=0.;
	if( E>2.4 ) {
		single_strip=pow(10.,(-0.0013*pow(lg10E,6.)+0.0313*pow(lg10E,5.)-0.2912*pow(lg10E,4.)+1.2932*pow(lg10E,3.)-2.8823*pow(lg10E,2.)+3.276*lg10E-20.892));
		if (!std::isfinite(single_strip) || single_strip < 0.0) single_strip = 0.0;
	}
	return single_strip;
}

double cs_HM_double_strip( double energy, double M_IONS ) {
	if (energy <= 0.0 || M_IONS <= 0.0 || !std::isfinite(energy) || !std::isfinite(M_IONS)) {
		return 0.0;
	}
	
	double E = energy*M_IONS;
	if (E <= 0.0 || !std::isfinite(E)) return 0.0;

	double lg10E=log10(E);
	if (!std::isfinite(lg10E)) return 0.0;
	
	double double_strip=0.;
	
	if( E>1.e3 ) {
		double_strip=pow(10.,(-0.010114*pow(lg10E,6.)+0.303523*pow(lg10E,5.)-3.711695*pow(lg10E,4.)+23.674607*pow(lg10E,3.)-83.406121*pow(lg10E,2.)+154.867111*lg10E-138.733509-1));
		if (!std::isfinite(double_strip) || double_strip < 0.0) double_strip = 0.0;
	}
	
	return double_strip;
}

// H-+H2->H-+H2++e OR H0+H2->H0+H2++e
double cs_bkg_ionization( double energy, double M_IONS ) {
// From Fubiani thesis on EAMCC
	if (energy <= 0.0 || M_IONS <= 0.0 || !std::isfinite(energy) || !std::isfinite(M_IONS)) {
		return 0.0;
	}
	
	double E = energy*M_IONS;
	if (E <= 0.0 || !std::isfinite(E)) return 0.0;
	
	double cs=0.;
	double lg10E=log10(E);
	if (!std::isfinite(lg10E)) return 0.0;

    if( E>10 && E<1e7 ) {  // Fixed condition - was OR, should be AND
        cs=1e-4*pow(10.,3.644e-3*pow(lg10E,4.)-5.939e-2*pow(lg10E,3.)-1.002e-1*pow(lg10E,2.)+3.348*lg10E-25.);
        if (!std::isfinite(cs) || cs < 0.0) cs = 0.0;
    }

    return cs;
}

// H0+H2->H++H2+e ionization of projectile H0
double  cs_proj_ionization_H0( double energy, double M_IONS ) {
	double E = energy*M_IONS;
	double cs=0.;
	double lg10E=log10(E);
    cs=pow(10.,0.000186*pow(lg10E,6)-0.0004*pow(lg10E,5)-0.055208*pow(lg10E,4)+0.73408*pow(lg10E,3)-4.288931*pow(lg10E,2)+12.812774*lg10E-34.752437-1);

    return cs;
}


// H++H2->H0+H2+ CX of projectile
double cs_CX_Hp( double energy, double M_IONS ) {
	double E = energy*M_IONS;
	double cs=0.;
	double lg10E=log10(E);
    cs=pow(10.,1.26307E-02*pow(lg10E,6)-2.57671E-01*pow(lg10E,5)+ 2.04941E+00*pow(lg10E,4)-8.16387E+00*pow(lg10E,3)+1.67493E+01*pow(lg10E,2)-1.47594E+01*lg10E-1.80339E+01);
    return cs;
}