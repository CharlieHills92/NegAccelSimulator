#ifndef THCALLBACK_H_
#define THCALLBACK_H_

#include "particledatabase.hpp"
#include "cross_sections.h"
#include "funct.h"
#include "SimulationParameters.h"
#include "random.hpp"
#include "globals.h"
#include "constants.hpp"
#include <cmath>
#include <algorithm>
//~ #include "src/pier/matlabfun.h"

using namespace std;

// Add safety constants
#define MAX_GENERATION_DEPTH 5
#define MAX_PARTICLE_COUNT 10000000
#define MIN_VELOCITY_THRESHOLD 1e3
#define MAX_VELOCITY_SCALE SPEED_C

//~ void get_dens_at_z( double z, double & nH2, double & nH0 ) {

    //~ double arr[] = { -0.05,0.,0.005,0.01,0.0155,0.026,0.03,0.035,0.05,0.055,0.062,0.064,0.074,0.08,0.084,0.1,0.15 };
    //~ vector<double> pos( arr, arr+sizeof(arr)/sizeof(arr[0]) );
    //~ double arr2[] = { 2.50E+19,2.50E+19,2.60E+19,2.65E+19,2.50E+19,1.70E+19,1.50E+19,1.41E+19,1.38E+19,1.34E+19,1.11E+19,1.00E+19,5.00E+18,2.50E+18,2.10E+18,2.00E+18,2.00E+18 };
    //~ vector<double> dens( arr2, arr2+sizeof(arr2)/sizeof(arr2[0]) );

//~ };


int check_periodicity( ParticleDataBase3D* _pdb, ParticleBase *particle, ParticleP3D *pcur, ParticleP3D *pend, const vector<double>& _periodicity );

namespace thcallback_detail {

enum ProductSpeedClass {
    PRODUCT_SPEED_WRONG = 0,
    PRODUCT_SPEED_FAST,
    PRODUCT_SPEED_SLOW
};

inline ProductSpeedClass product_speed_class_from_config_name(const std::string& speed_class) {
    if (speed_class == "fast") {
        return PRODUCT_SPEED_FAST;
    }
    if (speed_class == "slow") {
        return PRODUCT_SPEED_SLOW;
    }
    return PRODUCT_SPEED_WRONG;
}

struct GenericCrossSectionProcessProduct {
    particle_kind kind;
    ProductSpeedClass speed_class;
    uint count;

    GenericCrossSectionProcessProduct()
        : kind(PARTICLE_WRONG),
          speed_class(PRODUCT_SPEED_WRONG),
          count(1U) {}
};

struct GenericCrossSectionProcess {
    std::string process_id;
    particle_kind projectile_kind;
    bool consume_projectile;
    std::vector<GenericCrossSectionProcessProduct> products;

    GenericCrossSectionProcess()
        : projectile_kind(PARTICLE_WRONG),
          consume_projectile(false) {}
};

inline std::vector<GenericCrossSectionProcess> build_generic_processes(
    const std::vector<SimulationParameters::CrossSectionProcessDefinition>& definitions) {
    std::vector<GenericCrossSectionProcess> processes;
    for (std::vector<SimulationParameters::CrossSectionProcessDefinition>::const_iterator it = definitions.begin();
         it != definitions.end();
         ++it) {
        GenericCrossSectionProcess process;
        process.process_id = it->processId;
        process.projectile_kind = particle_kind_from_config_name(it->projectileKind);
        process.consume_projectile = it->projectileFate != "survive";
        if (process.projectile_kind == PARTICLE_WRONG) {
            continue;
        }

        for (std::vector<SimulationParameters::CrossSectionProcessProductDefinition>::const_iterator product_it = it->products.begin();
             product_it != it->products.end();
             ++product_it) {
            GenericCrossSectionProcessProduct product;
            product.kind = particle_kind_from_config_name(product_it->particleKind);
            product.speed_class = product_speed_class_from_config_name(product_it->speedClass);
            product.count = product_it->count;
            if (product.kind == PARTICLE_WRONG || product.speed_class == PRODUCT_SPEED_WRONG || product.count == 0U) {
                continue;
            }
            process.products.push_back(product);
        }

        processes.push_back(process);
    }

    return processes;
}

inline double slow_velocity_for_particle_kind(particle_kind kind, double ion_mass_u) {
    double mass_u = particle_kind_mass_u(kind, ion_mass_u);
    if (mass_u <= 0.0 || !std::isfinite(mass_u)) {
        return 0.0;
    }
    double energy_ev = particle_kind_is_electron(kind) ? 1.0 : 5.0;
    return std::sqrt(2.0*energy_ev*CHARGE_E/(mass_u*MASS_U));
}

inline bool sample_isotropic_velocity(Random* rng, double speed, double velocity[3]) {
    if (rng == NULL || speed <= 0.0 || !std::isfinite(speed)) {
        return false;
    }

    double rand_cos_theta[1];
    double rand_phi[1];
    rng->get(rand_cos_theta);
    rng->get(rand_phi);

    double cos_theta = 1.0 - 2.0*rand_cos_theta[0];
    if (!std::isfinite(cos_theta)) {
        return false;
    }
    double sin_theta_sq = std::max(0.0, 1.0 - cos_theta*cos_theta);
    double sin_theta = std::sqrt(sin_theta_sq);
    double phi = 2.0*std::acos(-1.0)*rand_phi[0];

    velocity[0] = speed*sin_theta*std::cos(phi);
    velocity[1] = speed*sin_theta*std::sin(phi);
    velocity[2] = speed*cos_theta;

    for (int ii = 0; ii < 3; ++ii) {
        if (!std::isfinite(velocity[ii]) || std::abs(velocity[ii]) > MAX_VELOCITY_SCALE) {
            return false;
        }
    }

    return true;
}

inline double child_particle_current(
    double parent_current,
    const GenericCrossSectionProcessProduct& product) {
    double sign = particle_kind_is_positive_ion(product.kind) ? -1.0 : 1.0;
    return sign*parent_current*static_cast<double>(product.count);
}

} // namespace thcallback_detail

class THCallback_test : public TrajectoryHandlerCallback {
public:

    THCallback_test() {}

    virtual ~THCallback_test() {}

    virtual void operator()( ParticleBase *particle, ParticlePBase *xcur, ParticlePBase *xend ) {
        ParticleP3D *pcur = (ParticleP3D *)( xcur );
        ParticleP3D *pend = (ParticleP3D *)( xend );
        
        //~ cout << "DEBUG -> v1=" << pcur->speed() << " v2=" << pend->speed() << " x1=" << (*pcur)[5] << " x2=" << (*pend)[5] << endl;
        
        // Kill particles with mass less than 0.5 atomic mass and x-coordinate more than 3 mm.
        if( particle->m() > 1.5e-27 && (*pcur)[5] >= 20.0e-3 ) {
            *pend = *pcur;
            particle->set_status( PARTICLE_STRIP );
        }
    }
};



class THCallback_strip : public TrajectoryHandlerCallback {

private:

    Random* _rng;
    ofstream* _debugstream;
    ParticleDataBase3D* _pdb;
    bool _debugprint;
    double _mass;
    vector<double> pos;
    vector<double> dens;
    vector<double> _periodicity;
    bool _periodic;
    vector<thcallback_detail::GenericCrossSectionProcess> _processes;
    
    int _counter;

    void handleConfiguredProcesses( ParticleBase *particle,
                                    double energy,
                                    double target_n,
                                    double delta ) {
        particle_kind projectile_kind = identify_particle_species( particle->m(), particle->q(), _mass );
        if( projectile_kind == PARTICLE_WRONG ) {
            return;
        }

        double sigma_total = 0.0;
        for( size_t ii = 0; ii < _processes.size(); ++ii ) {
            const thcallback_detail::GenericCrossSectionProcess& process = _processes[ii];
            if( !process.consume_projectile || process.projectile_kind != projectile_kind ) {
                continue;
            }

            double sigma = evaluate_cross_section_process( process.process_id, energy, _mass );
            if( std::isfinite(sigma) && sigma > 0.0 ) {
                sigma_total += sigma;
            }
        }

        if( sigma_total <= 0.0 || !std::isfinite(sigma_total) || target_n <= 0.0 || delta <= 0.0 ) {
            return;
        }

        double prob_collision = 1.0-exp(-sigma_total*target_n*delta);
        double rand1[1];
        _rng->get( rand1 );
        if( rand1[0] < prob_collision ) {
            particle->set_status( PARTICLE_STRIP );
        }
    }


public:

    THCallback_strip( bool debugprint, ParticleDataBase3D* pdb, double& mass,
                      const std::string& density_filename )
        : THCallback_strip( debugprint,
                            pdb,
                            mass,
                            vector<double>(),
                            std::vector<SimulationParameters::CrossSectionProcessDefinition>(),
                            density_filename ) {}

    THCallback_strip( bool debugprint, ParticleDataBase3D* pdb, double& mass,
                      const std::vector<SimulationParameters::CrossSectionProcessDefinition>& processes,
                      const std::string& density_filename )
        : THCallback_strip( debugprint, pdb, mass, vector<double>(), processes, density_filename ) {}

    THCallback_strip( bool debugprint, ParticleDataBase3D* pdb, double& mass,
                      const vector<double>& periodicity, const std::string& density_filename )
        : THCallback_strip( debugprint,
                            pdb,
                            mass,
                            periodicity,
                            std::vector<SimulationParameters::CrossSectionProcessDefinition>(),
                            density_filename ) {}

    THCallback_strip( bool debugprint, ParticleDataBase3D* pdb, double& mass,
                      const vector<double>& periodicity,
                      const std::vector<SimulationParameters::CrossSectionProcessDefinition>& processes,
                      const std::string& density_filename ) {
        //Random *rng; // Random number generator for two uniform random positions, and three Gaussian velocity components (to get a Maxwell-Boltzmann distribution for the total velocity).
        // if( ibsimu.get_rng_type() == RNG_SOBOL )  //In addition, two more numbers are needed for the cosine distribution and the angles.
        // 	_rng = new QRandom( 0 );
        // else
        // 	_rng = new MTRandom( 0 );
        _rng = new MTRandom( 1 );
        double qx[1];
        _rng->get( qx );
        //~ cout << "DEBUG: random generator initialized\n";
        // DEBUG
        _debugprint=debugprint;
        _pdb=pdb;
        _mass=mass;
        _periodicity=periodicity;
        _processes = thcallback_detail::build_generic_processes(processes);
        _periodic=false;
        if (_periodicity.size()>0) {
            _periodic=true;
            logfile << " THC_STRIPPING: THC PERIODICITY IS ON " << endl;
            logfile << " Periodic boundaries: ";
            for (int qq=0; qq<4; qq++) logfile << _periodicity[qq] << " ";
            logfile << endl << flush;
        }
        else {
            if (debug) {
                logfile << " THC_STRIPPING: THC PERIODICITY IS OFF " << endl;
            }
        }

        _counter=0;
        load_density_profile(density_filename, pos, dens);
        // cout << "Immhererad\n";
        // for(uint kk=0; kk<pos.size(); kk++) cout << pos[kk] << "\t" << dens[kk] << endl;
        // cout << "Immhererad\n";
        if( _debugprint ) {
            _debugstream=new ofstream("THCdebug.txt");
            *_debugstream << "zstart(m)" << " ";
            *_debugstream << "zend(m)" << " ";
            *_debugstream << "energy(eV)" << " ";
            *_debugstream << "sigma(m-2)" << " ";
            *_debugstream << "density(m-3)" << " ";
            *_debugstream << "delta(m)" << " ";
            *_debugstream << "mfp(m)" << " ";
            *_debugstream << "rand" << " ";
            *_debugstream << "prob_coll" << " ";
            *_debugstream << "collided" << "\n" << flush;
        }
    }

    virtual ~THCallback_strip() {
        if( _debugprint ) {
            _debugstream->close();
        }
        
    }

    virtual void operator()( ParticleBase *particle, ParticlePBase *xcur, ParticlePBase *xend ) {
        ParticleP3D *pcur = (ParticleP3D *)( xcur );
        ParticleP3D *pend = (ParticleP3D *)( xend );
        double vel = pcur->speed();
        // double energy = 0.5*particle->m()*vel*vel/1.6e-19;
		double energy = particle->m()*SPEED_C2*(1./sqrt(1.-(vel*vel)/(SPEED_C2))-1.)/CHARGE_E; // in eV
        // cout << "zLOC = " << (*pend)[5] << endl;
        if ( (*pend)[5]>0.007 && (*pend)[6]>0  ) {
            if( !_processes.empty() ) {
                if( particle->m() > 0.0 && energy >= 10.0 ) {
                    double zc=(*pcur)[5];
                    bool ciaone;
                    double target_n=density_at_z( zc, 0.3, pos, dens, ciaone );
                    double delta = sqrt( ((*pend)[1]-(*pcur)[1])*((*pend)[1]-(*pcur)[1])+((*pend)[3]-(*pcur)[3])*((*pend)[3]-(*pcur)[3])+((*pend)[5]-(*pcur)[5])*((*pend)[5]-(*pcur)[5]) );
                    handleConfiguredProcesses( particle, energy, target_n, delta );
                }
            }

            // PERIODIC DOMAIN
            if (_periodic) {
                // periodicity vector contains the dimensions of the two periodicities:
                // _periodicity[0]: positive x
                // _periodicity[1]: negative x
                // _periodicity[2]: positive y
                // _periodicity[3]: negative y
                check_periodicity(_pdb, particle, pcur, pend, _periodicity);
            }
        }

    }
        
};

class THCallback_secondaries : public TrajectoryHandlerCallback {

private:

    Random* _rng;
    ParticleDataBase3D* _pdb;
    double _mass;
    vector<double> pos;
    vector<double> dens;
    bool _periodic;
    vector<double> _periodicity;
    double _secondary_z_min;
    size_t _initial_particle_count;
    size_t _secondary_debug_count;
    vector<thcallback_detail::GenericCrossSectionProcess> _processes;

    bool shouldLogSecondaryDebug(size_t predicted_index) const {
         return predicted_index >= _initial_particle_count + 10300 &&
             predicted_index < _initial_particle_count + 10550;
    }

    void logSecondaryDebug(const char *branch,
                           ParticleBase *particle,
                           const ParticleP3D &parent_state,
                           const ParticleP3D &child_state,
                           double child_current,
                           double child_charge,
                           double child_mass,
                           int child_generation,
                           double energy,
                           double delta,
                           double target_n) {
        size_t predicted_index = _pdb->size();
        if (!shouldLogSecondaryDebug(predicted_index)) {
            return;
        }

        double child_speed = sqrt(child_state[2]*child_state[2] +
                                  child_state[4]*child_state[4] +
                                  child_state[6]*child_state[6]);

        logfile << "SECONDARY_DEBUG"
                << " event=" << _secondary_debug_count++
                << " predicted_index=" << predicted_index
                << " branch=" << branch
                << " parent_gen=" << particle->gen()
                << " parent_status=" << particle->get_status()
                << " parent_q=" << particle->q()
                << " parent_m=" << particle->m()
                << " parent_IQ=" << particle->IQ()
                << " child_gen=" << child_generation
                << " child_q=" << child_charge
                << " child_m=" << child_mass
                << " child_IQ=" << child_current
                << " energy=" << energy
                << " delta=" << delta
                << " target_n=" << target_n
                << " parent_t=" << parent_state[0]
                << " parent_x=" << parent_state[1]
                << " parent_vx=" << parent_state[2]
                << " parent_y=" << parent_state[3]
                << " parent_vy=" << parent_state[4]
                << " parent_z=" << parent_state[5]
                << " parent_vz=" << parent_state[6]
                << " child_t=" << child_state[0]
                << " child_x=" << child_state[1]
                << " child_vx=" << child_state[2]
                << " child_y=" << child_state[3]
                << " child_vy=" << child_state[4]
                << " child_z=" << child_state[5]
                << " child_vz=" << child_state[6]
                << " child_speed=" << child_speed
                << endl;
    }

    void addSecondaryParticleDebug(const char *branch,
                                   ParticleBase *particle,
                                   const ParticleP3D &parent_state,
                                   const ParticleP3D &child_state,
                                   double child_current,
                                   double child_charge,
                                   double child_mass,
                                   int child_generation,
                                   double energy,
                                   double delta,
                                   double target_n) {
        logSecondaryDebug(branch, particle, parent_state, child_state,
                          child_current, child_charge, child_mass,
                          child_generation, energy, delta, target_n);
        _pdb->add_particle(child_current, child_charge, child_mass, child_generation, child_state);
    }

    void handleConfiguredProcesses( ParticleBase *particle,
                                    ParticleP3D *pcur,
                                    ParticleP3D *pend,
                                    double energy,
                                    double target_n,
                                    double delta ) {
        particle_kind projectile_kind = identify_particle_species( particle->m(), particle->q(), _mass );
        if( projectile_kind == PARTICLE_WRONG ) {
            return;
        }

        std::vector<const thcallback_detail::GenericCrossSectionProcess*> matching_processes;
        std::vector<double> matching_sigmas;
        double sigma_total = 0.0;

        for( size_t ii = 0; ii < _processes.size(); ++ii ) {
            const thcallback_detail::GenericCrossSectionProcess& process = _processes[ii];
            if( process.projectile_kind != projectile_kind ) {
                continue;
            }

            double sigma = evaluate_cross_section_process( process.process_id, energy, _mass );
            if( !std::isfinite(sigma) || sigma <= 0.0 ) {
                continue;
            }

            matching_processes.push_back( &process );
            matching_sigmas.push_back( sigma );
            sigma_total += sigma;
        }

        if( matching_processes.empty() || sigma_total <= 0.0 || !std::isfinite(sigma_total) || target_n <= 0.0 || delta <= 0.0 ) {
            return;
        }

        double prob_collision = 1.0-exp(-sigma_total*target_n*delta);
        double rand1[1];
        _rng->get( rand1 );
        if( rand1[0] >= prob_collision ) {
            return;
        }

        *pend = *pcur;
        _rng->get( rand1 );
        double selector = rand1[0]*sigma_total;
        double sigma_cursor = 0.0;
        const thcallback_detail::GenericCrossSectionProcess* selected_process = matching_processes.back();
        for( size_t ii = 0; ii < matching_processes.size(); ++ii ) {
            sigma_cursor += matching_sigmas[ii];
            if( selector <= sigma_cursor ) {
                selected_process = matching_processes[ii];
                break;
            }
        }

        if( selected_process->consume_projectile ) {
            particle->set_status( PARTICLE_STRIP );
        }

        const ParticleP3D parent_state = *pend;
        int child_generation = particle->gen() + 1;

        for( size_t ii = 0; ii < selected_process->products.size(); ++ii ) {
            const thcallback_detail::GenericCrossSectionProcessProduct& product = selected_process->products[ii];
            double mass_u = particle_kind_mass_u( product.kind, _mass );
            double charge_state = particle_kind_charge_state( product.kind );
            if( mass_u <= 0.0 || !std::isfinite(mass_u) ) {
                continue;
            }

            ParticleP3D child_state = parent_state;
            if( product.speed_class == thcallback_detail::PRODUCT_SPEED_SLOW ) {
                double sampled_velocity[3];
                double slow_speed = thcallback_detail::slow_velocity_for_particle_kind( product.kind, _mass );
                if( !thcallback_detail::sample_isotropic_velocity( _rng, slow_speed, sampled_velocity ) ) {
                    continue;
                }

                child_state = ParticleP3D( parent_state[0], parent_state[1], sampled_velocity[0],
                                           parent_state[3], sampled_velocity[1], parent_state[5], sampled_velocity[2] );
            }

            try {
                addSecondaryParticleDebug( selected_process->process_id.c_str(),
                                           particle,
                                           *pcur,
                                           child_state,
                                           thcallback_detail::child_particle_current( particle->IQ(), product ),
                                           charge_state,
                                           mass_u,
                                           child_generation,
                                           energy,
                                           delta,
                                           target_n );
            } catch (Error& e) {
                logfile << "ERROR: Failed to add configured secondary particle for process '"
                        << selected_process->process_id << "': "
                        << e.get_error_message() << endl;
            } catch (const std::exception& e) {
                logfile << "ERROR: Failed to add configured secondary particle for process '"
                        << selected_process->process_id << "': "
                        << e.what() << endl;
            }
        }
    }

    // Position validation
    bool validatePosition(ParticleP3D *pend) {
        for (int i = 0; i < 7; i++) {
            if (!std::isfinite((*pend)[i])) {
                logfile << "ERROR: Invalid particle position component " << i << " in secondary generation: " << (*pend)[i] << endl;
                return false;
            }
        }
        return true;
    }

public:

    THCallback_secondaries( ParticleDataBase3D* pdb, double& mass,
                            const std::string& density_filename,
                            double secondary_z_min = 7.0e-3 )
        : THCallback_secondaries( pdb,
                                  mass,
                                  vector<double>(),
                                  std::vector<SimulationParameters::CrossSectionProcessDefinition>(),
                                  density_filename,
                                  secondary_z_min ) {}

    THCallback_secondaries( ParticleDataBase3D* pdb, double& mass,
                            const std::vector<SimulationParameters::CrossSectionProcessDefinition>& processes,
                            const std::string& density_filename,
                            double secondary_z_min = 7.0e-3 )
        : THCallback_secondaries( pdb, mass, vector<double>(), processes, density_filename, secondary_z_min ) {}

    THCallback_secondaries( ParticleDataBase3D* pdb, double& mass, const vector<double>& periodicity,
                            const std::string& density_filename, double secondary_z_min = 7.0e-3 )
        : THCallback_secondaries( pdb,
                                  mass,
                                  periodicity,
                                  std::vector<SimulationParameters::CrossSectionProcessDefinition>(),
                                  density_filename,
                                  secondary_z_min ) {}

    THCallback_secondaries( ParticleDataBase3D* pdb, double& mass, const vector<double>& periodicity,
                            const std::vector<SimulationParameters::CrossSectionProcessDefinition>& processes,
                            const std::string& density_filename, double secondary_z_min = 7.0e-3 ) {
        _pdb=pdb;
        _mass=mass;
        _processes = thcallback_detail::build_generic_processes(processes);
        _secondary_z_min = secondary_z_min;
        _initial_particle_count = pdb->size();
        _secondary_debug_count = 0;
        if (debug) cout << " MASS : " << _mass << endl;
        _rng = new MTRandom( 1 );
        double qx[1];
        _rng->get( qx );
        load_density_profile(density_filename,pos,dens);
        logfile << "SECONDARY_DEBUG initial_particle_count=" << _initial_particle_count
            << " debug_window_start=" << (_initial_particle_count + 10300)
            << " debug_window_end=" << (_initial_particle_count + 10550)
            << " secondary_z_min=" << _secondary_z_min
            << endl;
        _periodicity=periodicity;
        _periodic=false;
        if (_periodicity.size()>0) {
            _periodic=true;
            logfile << " THC_SECONDARIES: THC PERIODICITY IS ON " << endl;
            logfile << " Periodic boundaries: ";
            for (int qq=0; qq<4; qq++) logfile << _periodicity[qq] << " ";
            logfile << endl << flush;
        }
        else {
            if (debug) {
                logfile << " THC_SECONDARIES: THC PERIODICITY IS OFF " << endl;
            }
        }
    }

    virtual ~THCallback_secondaries() {}

    virtual void operator()( ParticleBase *particle, ParticlePBase *xcur, ParticlePBase *xend ) {
        ParticleP3D *pcur = (ParticleP3D *)( xcur );
        ParticleP3D *pend = (ParticleP3D *)( xend );
        
        // Input validation
        double vel = pcur->speed();
        if (!std::isfinite(vel)) {
            if (debug) logfile << "WARNING: Skipping secondary generation for invalid velocity: " << vel << endl;
            return;
        }
        
        // Validate particle position
        if (!validatePosition(pend)) {
			logfile << "ERROR: Invalid particle position detected. "
					<< "Mass=" << particle->m()
					<< " Charge=" << particle->q()
					<< " Current=" << particle->IQ()
					<< " status=" << particle->get_status()
					<< " x = " << (*pend)[1]
					<< " vx = " << (*pend)[2]
					<< " y = " << (*pend)[3]
					<< " vy = " << (*pend)[4]
					<< " z = " << (*pend)[5]
					<< " vz = " << (*pend)[6]
					<< " xs = " << (*pcur)[1]
					<< " vxs = " << (*pcur)[2]
					<< " ys = " << (*pcur)[3]
					<< " vys = " << (*pcur)[4]
					<< " zs = " << (*pcur)[5]
					<< " vzs = " << (*pcur)[6]
					<< endl << flush;
            return;
        }
        
        // Limit generation depth to prevent infinite chains
        int generation = particle->gen();
        if (generation > MAX_GENERATION_DEPTH) {
            if (debug) logfile << "WARNING: Maximum generation depth reached for particle, stopping secondaries" << endl;
            return;
        }
        
        // Check particle count to prevent memory overflow
        if (_pdb->size() > MAX_PARTICLE_COUNT) {
            logfile << "WARNING: Maximum particle count reached (" << _pdb->size() << "), stopping secondary generation" << endl;
            return;
        }
        
		double energy = particle->m()*SPEED_C2*(1./sqrt(1.-(vel*vel)/(SPEED_C2))-1.)/CHARGE_E; // in eV
                
        double zc=(*pcur)[5];
        bool ciaone = true;
        double target_n=density_at_z( zc, 0.3, pos, dens, ciaone );
        
        // Validate target density
        if (!std::isfinite(target_n) || target_n < 0) {
            if (debug) logfile << "WARNING: Invalid target density at z=" << zc << ": " << target_n << endl;
            return;
        }
        //double nu_scat=sigmav*target_n; //
        
        double delta = sqrt( ((*pend)[1]-(*pcur)[1])*((*pend)[1]-(*pcur)[1])+((*pend)[3]-(*pcur)[3])*((*pend)[3]-(*pcur)[3])+((*pend)[5]-(*pcur)[5])*((*pend)[5]-(*pcur)[5]) );
        
        // Validate step size
        if (!std::isfinite(delta) || delta < 0) {
            if (debug) logfile << "WARNING: Invalid step size in secondary generation: " << delta << endl;
            return;
        }

		if (delta == 0 ) {
			return;
		}
                    
		// Apply configured gas processes and optional periodic wrapping.
        if ( (*pend)[5] > _secondary_z_min && pend->speed() > 0.0 ) {
            if( !_processes.empty() ) {
                if( energy >= 10.0 && abs(particle->gen()%100)<5 ) {
                    handleConfiguredProcesses( particle, pcur, pend, energy, target_n, delta );
                }
            }
        
            // PERIODIC DOMAIN
            if (_periodic) {
                // periodicity vector contains the dimensions of the two periodicities:
                // _periodicity[0]: positive x
                // _periodicity[1]: negative x
                // _periodicity[2]: positive y
                // _periodicity[3]: negative y
                check_periodicity(_pdb, particle, pcur, pend, _periodicity);
            }
        }
    }

};



#endif /* THCALLBACK_H_ */
