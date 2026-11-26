#ifndef THCALLBACK_H_
#define THCALLBACK_H_

#include "particledatabase.hpp"
#include "cross_sections.h"
#include "funct.h"
#include "SimulationParameters.h"
#include "random.hpp"
#include "globals.h"
#include "constants.hpp"
#include "CallbackThreading.h"
#include <atomic>
#include <cmath>
#include <algorithm>
//~ #include "src/pier/matlabfun.h"

using namespace std;

// Add safety constants
#define MAX_GENERATION_DEPTH 5
#define MAX_PARTICLE_COUNT 10000000
#define MAX_SECONDARY_PARTICLES_PER_PASS 1000000
#define MIN_VELOCITY_THRESHOLD 1e3
#define MAX_VELOCITY_SCALE SPEED_C
// Upper bound on configured gas processes considered for one projectile species in a
// single step, so the per-step candidate list can live on the stack. Processes beyond this
// many that match the same projectile are ignored, and the constructor warns if the
// configuration could ever reach it.
#define MAX_CONFIGURED_PROCESSES 32

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

inline bool sample_isotropic_velocity(double speed, double velocity[3]) {
    if (speed <= 0.0 || !std::isfinite(speed)) {
        return false;
    }

    double rand_cos_theta = cbthread::rng().uniform();
    double rand_phi = cbthread::rng().uniform();

    double cos_theta = 1.0 - 2.0*rand_cos_theta;
    if (!std::isfinite(cos_theta)) {
        return false;
    }
    double sin_theta_sq = std::max(0.0, 1.0 - cos_theta*cos_theta);
    double sin_theta = std::sqrt(sin_theta_sq);
    double phi = 2.0*std::acos(-1.0)*rand_phi;

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
        if( cbthread::rng().uniform() < prob_collision ) {
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
        // Random numbers come from cbthread::rng(), a per-thread generator. See
        // CallbackThreading.h for why a shared ibsimu Random cannot be used once
        // iterate_trajectories() runs on more than one thread.
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

    ParticleDataBase3D* _pdb;
    double _mass;
    vector<double> pos;
    vector<double> dens;
    bool _periodic;
    vector<double> _periodicity;
    double _secondary_z_min;
    size_t _initial_particle_count;

    // Shared by every worker thread, so all mutable counters are atomic. The created
    // count also replaces _pdb->size() in the cap test: reading std::vector::size() while
    // another thread is inside add_particle()'s push_back() is a data race.
    std::atomic<size_t> _secondary_debug_count;
    std::atomic<size_t> _secondaries_created;
    std::atomic<size_t> _slow_neutrals_skipped;
    std::atomic<bool> _secondary_limit_warned;
    std::atomic<bool> _generation_limit_warned;
    vector<thcallback_detail::GenericCrossSectionProcess> _processes;

    // When true, slow neutral reaction products are accounted but not tracked. See
    // handleConfiguredProcesses.
    bool _skip_slow_neutral_products;

    /*! \brief Have we hit a particle cap?
     *
     *  Two distinct ceilings apply and they must not be conflated:
     *
     *   - MAX_SECONDARY_PARTICLES_PER_PASS bounds the secondaries CREATED in this pass, so
     *     that a runaway cascade cannot grow without limit. Primaries are not part of it:
     *     they are an input, not something the cascade produced.
     *   - MAX_PARTICLE_COUNT bounds the TOTAL population, primaries included.
     *
     *  The previous form added the primary count to the secondary counter and compared the
     *  sum against the SECONDARY cap. A run whose configured primary count already reached
     *  that cap therefore exhausted the budget on the very first call and produced no
     *  secondaries at all -- silently, and whatever the configured physics asked for. It
     *  also left MAX_PARTICLE_COUNT unreferenced, so the total ceiling was never enforced.
     */
    bool secondaryBudgetExhausted() const {
        const size_t created = _secondaries_created.load(std::memory_order_relaxed);
        return created >= MAX_SECONDARY_PARTICLES_PER_PASS ||
               _initial_particle_count + created >= MAX_PARTICLE_COUNT;
    }

    /*! \brief Emit a warn-once message, at most one thread at a time. */
    void warnOnce(std::atomic<bool>& flag, const char* message) {
        bool expected = false;
        if (flag.compare_exchange_strong(expected, true)) {
            CB_LOG_LOCK();
            logfile << "WARNING: " << message << endl;
        }
    }

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
        size_t predicted_index = _initial_particle_count +
                                 _secondaries_created.load(std::memory_order_relaxed);
        if (!shouldLogSecondaryDebug(predicted_index)) {
            return;
        }

        double child_speed = sqrt(child_state[2]*child_state[2] +
                                  child_state[4]*child_state[4] +
                                  child_state[6]*child_state[6]);

        CB_LOG_LOCK();
        logfile << "SECONDARY_DEBUG"
                << " event=" << _secondary_debug_count.fetch_add(1, std::memory_order_relaxed)
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
        if (secondaryBudgetExhausted()) {
            warnOnce(_secondary_limit_warned,
                     "Secondary particle pass limit reached; stopping further gas "
                     "secondary creation");
            return;
        }
        // Safe from a worker thread: add_particle() takes the scheduler's problem-vector
        // mutex (particledatabaseimp.hpp:626-634).
        _pdb->add_particle(child_current, child_charge, child_mass, child_generation, child_state);
        _secondaries_created.fetch_add(1, std::memory_order_relaxed);
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

        // Fixed-capacity scratch instead of two std::vectors. This runs on every
        // integration step of every particle, so the two heap allocations the vectors used
        // to make were pure per-step overhead -- and with eight worker threads they would
        // also have been eight threads contending on the allocator.
        const size_t process_count = _processes.size();
        if( process_count == 0 ) {
            return;
        }
        const thcallback_detail::GenericCrossSectionProcess* matching_processes[MAX_CONFIGURED_PROCESSES];
        double matching_sigmas[MAX_CONFIGURED_PROCESSES];
        size_t matching_count = 0;
        double sigma_total = 0.0;

        for( size_t ii = 0; ii < process_count && matching_count < MAX_CONFIGURED_PROCESSES; ++ii ) {
            const thcallback_detail::GenericCrossSectionProcess& process = _processes[ii];
            if( process.projectile_kind != projectile_kind ) {
                continue;
            }

            double sigma = evaluate_cross_section_process( process.process_id, energy, _mass );
            if( !std::isfinite(sigma) || sigma <= 0.0 ) {
                continue;
            }

            matching_processes[matching_count] = &process;
            matching_sigmas[matching_count] = sigma;
            ++matching_count;
            sigma_total += sigma;
        }

        if( matching_count == 0 || sigma_total <= 0.0 || !std::isfinite(sigma_total) || target_n <= 0.0 || delta <= 0.0 ) {
            return;
        }

        double prob_collision = 1.0-exp(-sigma_total*target_n*delta);
        if( cbthread::rng().uniform() >= prob_collision ) {
            return;
        }

        *pend = *pcur;
        double selector = cbthread::rng().uniform()*sigma_total;
        double sigma_cursor = 0.0;
        const thcallback_detail::GenericCrossSectionProcess* selected_process = matching_processes[matching_count-1];
        for( size_t ii = 0; ii < matching_count; ++ii ) {
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

            // Slow neutral products: created only if the caller asked for them.
            //
            // A 5 eV neutral is field-blind and moves at ~3e4 m/s, two orders of magnitude
            // below the beam, so within the 10000-step budget it travels a fraction of a
            // millimetre and terminates in vacuum by PARTICLE_NSTP. It therefore reaches no
            // electrode and contributes nothing to any grid power or current column, while
            // costing a full 10000-step integration each -- the most expensive particle in
            // the run, per unit of physics delivered.
            //
            // It is NOT free to drop, which is why this is opt-in: ibsimu deposits space
            // charge from IQ alone, without consulting q (particleiterator.hpp:923,932), so
            // these macroparticles currently deposit charge as though they were still
            // charged, and slowly, which weights them heavily. Skipping them changes the
            // space-charge field of the pass that creates them. That pass is always the
            // final one (secondaries are only enabled on the last iteration) and its space
            // charge is never fed back into a Poisson solve, so no trajectory or field the
            // simulation reports is affected -- but the saved scharge and any density
            // diagnostic drawn from it will change. See getStrippingSkipSlowNeutrals().
            //
            // Fast neutrals -- the stripped beam, the quantity of interest downstream --
            // are deliberately untouched: only PRODUCT_SPEED_SLOW is considered here.
            const bool is_slow = ( product.speed_class == thcallback_detail::PRODUCT_SPEED_SLOW );
            if( _skip_slow_neutral_products && is_slow && charge_state == 0.0 ) {
                _slow_neutrals_skipped.fetch_add( 1, std::memory_order_relaxed );
                continue;
            }

            ParticleP3D child_state = parent_state;
            if( is_slow ) {
                double sampled_velocity[3];
                double slow_speed = thcallback_detail::slow_velocity_for_particle_kind( product.kind, _mass );
                if( !thcallback_detail::sample_isotropic_velocity( slow_speed, sampled_velocity ) ) {
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
                CB_LOG_LOCK();
                logfile << "ERROR: Failed to add configured secondary particle for process '"
                        << selected_process->process_id << "': "
                        << e.get_error_message() << endl;
            } catch (const std::exception& e) {
                CB_LOG_LOCK();
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
                CB_LOG_LOCK();
                logfile << "ERROR: Invalid particle position component " << i << " in secondary generation: " << (*pend)[i] << endl;
                return false;
            }
        }
        return true;
    }

public:

    THCallback_secondaries( ParticleDataBase3D* pdb, double& mass,
                            const std::string& density_filename,
                            double secondary_z_min = 7.0e-3,
                            bool skip_slow_neutral_products = false )
        : THCallback_secondaries( pdb,
                                  mass,
                                  vector<double>(),
                                  std::vector<SimulationParameters::CrossSectionProcessDefinition>(),
                                  density_filename,
                                  secondary_z_min,
                                  skip_slow_neutral_products ) {}

    THCallback_secondaries( ParticleDataBase3D* pdb, double& mass,
                            const std::vector<SimulationParameters::CrossSectionProcessDefinition>& processes,
                            const std::string& density_filename,
                            double secondary_z_min = 7.0e-3,
                            bool skip_slow_neutral_products = false )
        : THCallback_secondaries( pdb, mass, vector<double>(), processes, density_filename,
                                  secondary_z_min, skip_slow_neutral_products ) {}

    THCallback_secondaries( ParticleDataBase3D* pdb, double& mass, const vector<double>& periodicity,
                            const std::string& density_filename, double secondary_z_min = 7.0e-3,
                            bool skip_slow_neutral_products = false )
        : THCallback_secondaries( pdb,
                                  mass,
                                  periodicity,
                                  std::vector<SimulationParameters::CrossSectionProcessDefinition>(),
                                  density_filename,
                                  secondary_z_min,
                                  skip_slow_neutral_products ) {}

    THCallback_secondaries( ParticleDataBase3D* pdb, double& mass, const vector<double>& periodicity,
                            const std::vector<SimulationParameters::CrossSectionProcessDefinition>& processes,
                            const std::string& density_filename, double secondary_z_min = 7.0e-3,
                            bool skip_slow_neutral_products = false )
        : _secondary_debug_count(0), _secondaries_created(0), _slow_neutrals_skipped(0),
          _secondary_limit_warned(false), _generation_limit_warned(false),
          _skip_slow_neutral_products(skip_slow_neutral_products) {
        _pdb=pdb;
        _mass=mass;
        _processes = thcallback_detail::build_generic_processes(processes);
        _secondary_z_min = secondary_z_min;
        _initial_particle_count = pdb->size();
        if (debug) cout << " MASS : " << _mass << endl;
        // Random numbers come from cbthread::rng(), a per-thread generator; see
        // CallbackThreading.h.
        load_density_profile(density_filename,pos,dens);
        if (_processes.size() > MAX_CONFIGURED_PROCESSES) {
            logfile << "WARNING: " << _processes.size() << " gas processes configured but only "
                    << MAX_CONFIGURED_PROCESSES << " can be considered per step; "
                    << "raise MAX_CONFIGURED_PROCESSES in THCallback.h" << endl;
        }
        logfile << "SECONDARY_DEBUG initial_particle_count=" << _initial_particle_count
            << " debug_window_start=" << (_initial_particle_count + 10300)
            << " debug_window_end=" << (_initial_particle_count + 10550)
            << " secondary_z_min=" << _secondary_z_min
            << " skip_slow_neutral_products=" << (_skip_slow_neutral_products ? 1 : 0)
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

    virtual ~THCallback_secondaries() {
        const size_t skipped = _slow_neutrals_skipped.load(std::memory_order_relaxed);
        if (skipped > 0) {
            logfile << "THC_SECONDARIES: " << skipped << " slow neutral product(s) accounted "
                    << "but not tracked (physics.stripping.skipSlowNeutralProducts); "
                    << _secondaries_created.load(std::memory_order_relaxed)
                    << " gas secondaries created" << endl;
        }
    }

    /*! \brief Gas secondary macroparticles created in the current pass. */
    size_t secondariesCreated() const {
        return _secondaries_created.load(std::memory_order_relaxed);
    }

    /*! \brief Slow neutral products accounted but not tracked. */
    size_t slowNeutralsSkipped() const {
        return _slow_neutrals_skipped.load(std::memory_order_relaxed);
    }

    virtual void operator()( ParticleBase *particle, ParticlePBase *xcur, ParticlePBase *xend ) {
        ParticleP3D *pcur = (ParticleP3D *)( xcur );
        ParticleP3D *pend = (ParticleP3D *)( xend );
        
        // Input validation
        double vel = pcur->speed();
        if (!std::isfinite(vel)) {
            if (debug) {
                CB_LOG_LOCK();
                logfile << "WARNING: Skipping secondary generation for invalid velocity: " << vel << endl;
            }
            return;
        }
        
        // Validate particle position
        if (!validatePosition(pend)) {
			CB_LOG_LOCK();
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
        if (generation >= MAX_GENERATION_DEPTH) {
            if (debug) {
                warnOnce(_generation_limit_warned,
                         "Maximum gas-secondary generation depth reached; stopping further "
                         "gas secondary creation");
            }
            return;
        }

        // Check particle count to prevent memory overflow
        if (secondaryBudgetExhausted()) {
            warnOnce(_secondary_limit_warned,
                     "Secondary particle pass limit reached; stopping further gas "
                     "secondary creation");
            return;
        }
        
		double energy = particle->m()*SPEED_C2*(1./sqrt(1.-(vel*vel)/(SPEED_C2))-1.)/CHARGE_E; // in eV
                
        double zc=(*pcur)[5];
        bool ciaone = true;
        double target_n=density_at_z( zc, 0.3, pos, dens, ciaone );
        
        // Validate target density
        if (!std::isfinite(target_n) || target_n < 0) {
            if (debug) {
                CB_LOG_LOCK();
                logfile << "WARNING: Invalid target density at z=" << zc << ": " << target_n << endl;
            }
            return;
        }
        //double nu_scat=sigmav*target_n; //
        
        double delta = sqrt( ((*pend)[1]-(*pcur)[1])*((*pend)[1]-(*pcur)[1])+((*pend)[3]-(*pcur)[3])*((*pend)[3]-(*pcur)[3])+((*pend)[5]-(*pcur)[5])*((*pend)[5]-(*pcur)[5]) );
        
        // Validate step size
        if (!std::isfinite(delta) || delta < 0) {
            if (debug) {
                CB_LOG_LOCK();
                logfile << "WARNING: Invalid step size in secondary generation: " << delta << endl;
            }
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
