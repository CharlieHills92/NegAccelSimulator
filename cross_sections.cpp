
#include "cross_sections.h"
#include <atomic>
#include <cmath>
#include <iostream>
#include <map>

// Physical ceiling on a single cross section, as a backstop against a runaway polynomial
// fit. The largest legitimate value in this data set is negative_ion_single_stripping at
// ~1.1e-19 m^2, and resonant charge exchange reaches ~1e-18 m^2 at low energy, so 1e-17
// is roughly 100x anything real and cannot clip physics. It exists because a log10-sigma
// polynomial evaluated outside its fitted range does not fail, it diverges: the
// h0_process_2 fit reached 1.2e-7 m^2, a mean free path of 4e-13 m, and hence a collision
// on every integration step. The primary defence is a correct maximumEnergyEV per process
// (set from the fit's own data range); this catches the case where that is missing.
//~ #include "src/pier/matlabfun.h"

using namespace std;

namespace {

struct ConfiguredCrossSectionProcess {
	string source_path;
	vector<double> coefficients;
	unsigned int fit_degree;
	bool scale_energy_by_ion_mass;
	double minimum_energy_ev;
	double maximum_energy_ev;

	ConfiguredCrossSectionProcess()
		: fit_degree(0U),
		  scale_energy_by_ion_mass(true),
		  minimum_energy_ev(0.0),
		  maximum_energy_ev(-1.0) {}
};

map<string, ConfiguredCrossSectionProcess> g_cross_section_processes;

double evaluate_configured_cross_section( const string& reaction_id,
									  double energy,
									  double M_IONS ) {
	map<string, ConfiguredCrossSectionProcess>::const_iterator it =
		g_cross_section_processes.find(reaction_id);
	if( it == g_cross_section_processes.end() ) {
		return 0.0;
	}

	if( energy <= 0.0 || M_IONS <= 0.0 || !std::isfinite(energy) || !std::isfinite(M_IONS) ) {
		return 0.0;
	}

	const ConfiguredCrossSectionProcess& process = it->second;

	// Velocity-equivalence scaling. The .cross tables are tabulated against the LAB ENERGY
	// OF A HYDROGEN projectile. A projectile of mass A (in u) at energy E has the velocity
	// of a hydrogen projectile at E_H = E*(m_H/m_A) = E/A, so the table must be read at
	// E/M_IONS. Fubiani et al. 2008 states the same convention: "we use the same cross
	// sections as for the case of hydrogen for identical particle velocities (v_H = v_D)",
	// with Fig. 4 noting E_0^(D) = 2 E_0^(H).
	//
	// This previously MULTIPLIED, evaluating a 1 MeV deuteron's cross section at 2 MeV
	// instead of 500 keV -- a factor of 4 in the argument, on the falling side of sigma(E).
	// Benchmarked against the paper's MAMuG case: multiplying gives 14.7% stripping loss
	// against their quoted 29.0%, dividing gives 27.6%. It is a no-op for hydrogen
	// (M_IONS = 1), which is why it stayed latent. Note THCallback_surf_EAMCC's RN/RE/
	// eta_ion already divide (E0 = E0/mass_ratio); the two halves of the code disagreed.
	double evaluation_energy = process.scale_energy_by_ion_mass ? energy/M_IONS : energy;
	if( evaluation_energy <= 0.0 || !std::isfinite(evaluation_energy) ) {
		return 0.0;
	}
	if( evaluation_energy < process.minimum_energy_ev ) {
		return 0.0;
	}
	if( process.maximum_energy_ev > 0.0 && evaluation_energy > process.maximum_energy_ev ) {
		return 0.0;
	}

	double lg10E = log10(evaluation_energy);
	if( !std::isfinite(lg10E) ) {
		return 0.0;
	}

	double lg10Sigma = 0.0;
	for( size_t ii = 0; ii < process.coefficients.size(); ++ii ) {
		lg10Sigma = lg10Sigma*lg10E + process.coefficients[ii];
	}

	double sigma = pow(10.0, lg10Sigma);
	if( !std::isfinite(sigma) || sigma < 0.0 ) {
		return 0.0;
	}

	if( sigma > MAX_PHYSICAL_CROSS_SECTION_M2 ) {
		// Warn once, naming the first offender, then clamp. Clamping rather than zeroing
		// keeps the reaction present instead of silently deleting it, but the run is not
		// trustworthy until the process gets a proper energy bound -- so say so loudly.
		static std::atomic<bool> warned(false);
		bool expected = false;
		if( warned.compare_exchange_strong(expected, true) ) {
			std::cerr << "WARNING: cross section for process '" << reaction_id
			          << "' evaluated to " << sigma << " m^2 at " << evaluation_energy
			          << " eV, above the physical ceiling of "
			          << MAX_PHYSICAL_CROSS_SECTION_M2 << " m^2. This means its polynomial "
			          << "fit is being extrapolated outside the energy range it was fitted "
			          << "over. Set minimumEnergyEV/maximumEnergyEV for this process to the "
			          << "range of its .cross data. Clamping; further such warnings are "
			          << "suppressed." << std::endl;
		}
		sigma = MAX_PHYSICAL_CROSS_SECTION_M2;
	}

	return sigma;
}

} // namespace


void clear_cross_section_processes() {
	g_cross_section_processes.clear();
}


void configure_cross_section_process( const std::string& reaction_id,
									 const std::string& source_path,
									 const std::vector<double>& coefficients,
									 unsigned int fit_degree,
									 bool scale_energy_by_ion_mass,
									 double minimum_energy_ev,
									 double maximum_energy_ev ) {
	ConfiguredCrossSectionProcess process;
	process.source_path = source_path;
	process.coefficients = coefficients;
	process.fit_degree = fit_degree;
	process.scale_energy_by_ion_mass = scale_energy_by_ion_mass;
	process.minimum_energy_ev = minimum_energy_ev;
	process.maximum_energy_ev = maximum_energy_ev;
	g_cross_section_processes[reaction_id] = process;
}


double evaluate_cross_section_process( const std::string& process_id,
								   double energy,
								   double M_IONS ) {
	return evaluate_configured_cross_section( process_id, energy, M_IONS );
}


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
	return evaluate_cross_section_process( "negative_ion_single_stripping", energy, M_IONS );
}

double cs_HM_double_strip( double energy, double M_IONS ) {
	return evaluate_cross_section_process( "negative_ion_double_stripping", energy, M_IONS );
}

// H-+H2->H-+H2++e OR H0+H2->H0+H2++e
double cs_bkg_ionization( double energy, double M_IONS ) {
	return evaluate_cross_section_process( "background_gas_ionization", energy, M_IONS );
}

// H0+H2->H++H2+e ionization of projectile H0
double  cs_proj_ionization_H0( double energy, double M_IONS ) {
	return evaluate_cross_section_process( "neutral_projectile_ionization", energy, M_IONS );
}


// H++H2->H0+H2+ CX of projectile
double cs_CX_Hp( double energy, double M_IONS ) {
	return evaluate_cross_section_process( "positive_ion_charge_exchange", energy, M_IONS );
}