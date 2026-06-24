
#include "cross_sections.h"
#include <cmath>
#include <map>
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
	double evaluation_energy = process.scale_energy_by_ion_mass ? energy*M_IONS : energy;
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