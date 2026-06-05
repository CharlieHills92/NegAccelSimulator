#ifndef THCALLBACK_SURF_EAMCC_H_
#define THCALLBACK_SURF_EAMCC_H_

#include "particledatabase.hpp"
#include "random.hpp"
#include "globals.h"
#include "constants.hpp"
#include "geometry.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

// Safety constants (same as THCallback.h)
#define MAX_GENERATION_DEPTH 5
#define MAX_PARTICLE_COUNT 10000000
#define MIN_VELOCITY_THRESHOLD 1e3
#define MAX_VELOCITY_SCALE SPEED_C
#define MAX_SECONDARY_ELECTRONS 26
#define SECONDARY_ELECTRON_ENERGY 10.0  // eV

// Generation number offset to distinguish surface-generated from volume-generated particles
// Volume-generated: gen + 1 (e.g., 0->1, 1->2, 2->3, ...)
// Surface-generated: gen + SURFACE_GENERATION_OFFSET (e.g., 0->101, 1->102, 2->103, ...)
// 
// To distinguish in code:
//   - Check if gen >= SURFACE_GENERATION_OFFSET for surface-generated
//   - Use gen % 100 to get effective generation (0-4) for diagnostics
//   - Volume: gen % 100 = 1,2,3,4,5...
//   - Surface: gen % 100 = 1,2,3,4,5... (same effective gen, but gen >= 101)
#define SURFACE_GENERATION_OFFSET 101

/**
 * @class THCallback_surf_EAMCC
 * @brief Handles surface collisions with EAMCC secondary particle generation
 * 
 * This class implements the EAMCC (Electron and Ion Monte Carlo Collision) 
 * surface interaction model for secondary particle generation.
 * Based on the Fortran routines in EAMCCsecondaries folder.
 * 
 * Features:
 * - Electron backscattering and secondary emission
 * - Ion backscattering and secondary electron emission
 * - Angular-dependent emission probabilities
 * - Proper energy and charge conservation
 */
class THCallback_surf_EAMCC : public TrajectoryEndCallback {
private:
    Geometry &_geom;
    Random* _rng;
    ParticleDataBase3D* _pdb;
    double _mass;
    double _minimum_z;
    bool _debugprint;
    
    // Helper function for safe velocity scaling
    bool validateAndScaleVelocity(double vel, double minvel, const double originalVel[3], double scaledVel[3]) {
        if (vel < MIN_VELOCITY_THRESHOLD || !std::isfinite(vel)) {
            if (debug) logfile << "WARNING: Invalid velocity in secondary generation: " << vel << endl;
            return false;
        }
        
        double scale_factor = std::min(1.0, minvel/vel);
        
        for (int i = 0; i < 3; i++) {
            scaledVel[i] = scale_factor * originalVel[i];
            
            if (!std::isfinite(scaledVel[i]) || abs(scaledVel[i]) > MAX_VELOCITY_SCALE) {
                if (debug) logfile << "WARNING: Invalid scaled velocity component " << i << ": " << scaledVel[i] << endl;
                scaledVel[i] = 0.0;
            }
        }
        
        return true;
    }
    
    // Position validation
    bool validatePosition(const Vec3D& loc) {
        if (!std::isfinite(loc[0]) || !std::isfinite(loc[1]) || !std::isfinite(loc[2])) {
            logfile << "ERROR: Invalid particle position in secondary generation" << endl;
            return false;
        }
        return true;
    }

    /**
     * @brief Calculate backscattering probability for electrons (ETA0)
     * Based on ORNL Redbook data and experimental fits
     * @param E0 Incident electron energy in eV
     * @return Backscattering probability (0-1)
     */
    double ETA0(double E0) {
        // Data from ORNL Redbook Vol 1, D-4 for Cu
        // Energy points in eV
        const double E[] = {0., 100., 200., 500., 1000., 2000., 4000., 7000.,
                           10000., 20000., 40000., 70000., 100000.,
                           200000., 400000., 700000., 1.0e6, 4.0e6};
        // Backscattering probabilities
        const double eta[] = {0.0, 0.124, 0.161, 0.262, 0.266, 0.282, 0.315, 0.322,
                             0.323, 0.318, 0.309, 0.299, 0.285,
                             0.275, 0.260, 0.240, 0.210, 0.11};
        const int n = 18;
        
        if (E0 <= E[0]) return eta[0];
        if (E0 >= E[n-1]) return eta[n-1];
        
        // Linear interpolation
        for (int i = 1; i < n; i++) {
            if (E0 > E[i-1] && E0 <= E[i]) {
                return eta[i-1] + (E0 - E[i-1]) * (eta[i] - eta[i-1]) / (E[i] - E[i-1]);
            }
        }
        
        return eta[n-1];
    }

    /**
     * @brief Calculate CT factor for backscattered electrons
     * Based on P-F Staub et al. [J. of Phys. D 27, 1533 (1994)]
     * @param E0 Incident electron energy in eV
     * @param eta0 Backscattering probability at normal incidence
     * @return CT factor for angular dependence
     */
    double CTback(double E0, double eta0) {
        double kappa = 1.0 - exp(-1.83 * pow(E0 * 1e-3, 0.25));
        return kappa * log(1.0 / eta0);
    }

    /**
     * @brief Calculate true secondary electron emission yield (ETA1)
     * Based on CERN LHC Project Report 472 (fully conditioned copper surfaces)
     * @param Ep Incident particle energy in eV
     * @return Secondary electron yield
     */
    double ETA1(double Ep) {
        // Parameters from CERN report
        const double s = 1.35;
        const double Emax = 318.0;  // eV
        const double DeltaMax = 1.13;
        const double A0 = 20.69989;
        const double A1 = -7.07605;
        const double A2 = 0.483547;
        const double E0 = 56.914686;  // eV
        
        double EpEm = Ep / Emax;
        double deltaS = DeltaMax * s * EpEm / (s - 1.0 + pow(EpEm, s));
        
        if (EpEm >= 1.0) {
            return deltaS;
        }
        
        double LnEpE0 = log(Ep + E0);
        double f = A0 + A1 * LnEpE0 + A2 * LnEpE0 * LnEpE0;
        f = exp(f);
        double deltaT = deltaS / (1.0 - f);
        return deltaT;
    }

    /**
     * @brief Calculate CT factor for true secondary electron emission
     * Based on Koshikawa et al. [J. of Phys. D 6, 1360 (1973)]
     * @param E0 Incident electron energy in eV
     * @return CT factor for angular dependence
     */
    double CTsec(double E0) {
        const double E[] = {0.0, 500.0, 1000.0, 2000.0, 3000.0, 6000.0, 10000.0};
        const double CT[] = {0.692, 0.692, 0.731, 0.835, 0.956, 0.971, 1.005};
        const int n = 7;
        
        if (E0 <= E[0]) return CT[0];
        if (E0 >= E[n-1]) return CT[n-1];
        
        // Linear interpolation
        for (int i = 1; i < n; i++) {
            if (E0 > E[i-1] && E0 <= E[i]) {
                return CT[i-1] + (E0 - E[i-1]) * (CT[i] - CT[i-1]) / (E[i] - E[i-1]);
            }
        }
        
        return CT[n-1];
    }

    /**
     * @brief Calculate backscattering coefficient for ions (RN)
     * Based on ORNL Redbook Vol 3, E22
     * @param E0 Incident ion energy in keV
     * @param mass_ratio Mass ratio (for D, H2, He scaling)
     * @return Backscattering coefficient
     */
    double RN(double E0, double mass_ratio = 1.0) {
        // Adjust energy for mass scaling (velocity scaling)
        E0 = E0 / mass_ratio;
        
        // Fit parameters: log(RN) = sum a_i * log(E0)^i
        const double a[] = {-1.55626, -0.2385, -0.0673113, -0.00899714, 3.23125e-6};
        const int n = 5;
        
        if (E0 <= 0.0) return 0.0;
        
        double logRN = 0.0;
        double logE = log(E0);
        for (int i = 0; i < n; i++) {
            logRN += a[i] * pow(logE, i);
        }
        
        return exp(logRN);
    }

    /**
     * @brief Calculate energy ratio for backscattered ions (RE)
     * Based on ORNL Redbook Vol 3, E23
     * @param E0 Incident ion energy in keV
     * @param mass_ratio Mass ratio (for D, H2, He scaling)
     * @return Energy ratio RE/RN
     */
    double RE(double E0, double mass_ratio = 1.0) {
        // Adjust energy for mass scaling
        E0 = E0 / mass_ratio;
        
        // Fit parameters: log(RE) = sum a_i * log(E0)^i
        const double a[] = {-2.33586, -0.450527, -0.0732621, -0.00687502, -4.56457e-6};
        const int n = 5;
        
        if (E0 <= 0.0) return 0.0;
        
        double logRE = 0.0;
        double logE = log(E0);
        for (int i = 0; i < n; i++) {
            logRE += a[i] * pow(logE, i);
        }
        
        return exp(logRE);
    }

    /**
     * @brief Calculate secondary electron yield for ion impacts (eta_ion)
     * Based on ORNL Redbook Vol 3, C10 and Svensson et al PRB 25
     * @param E0 Incident ion energy in keV
     * @param mass_ratio Mass ratio (for D, H2, He scaling)
     * @return Secondary electron yield
     */
    double eta_ion(double E0, double mass_ratio = 1.0) {
        // Adjust energy for mass scaling
        E0 = E0 / mass_ratio;
        
        // Fit parameters: log(eta_ion) = sum a_i * log(E0)^i
        const double a[] = {-1.3854, 0.544935, -0.0124299, -0.00307087, -0.000558421};
        const int n = 5;
        
        if (E0 <= 0.0) return 0.0;
        
        double log_eta = 0.0;
        double logE = log(E0);
        for (int i = 0; i < n; i++) {
            log_eta += a[i] * pow(logE, i);
        }
        
        return exp(log_eta);
    }

    /**
     * @brief Handle electron surface collision
     * @param particle Incident particle
     * @param loc Impact location
     * @param vel Incident velocity
     * @param energy Incident energy in eV
     * @param theta Impact angle with surface normal (radians)
     * @param normal Surface normal vector
     * @param tang1 First tangent vector
     * @param tang2 Second tangent vector
     */
    void handleElectronSurfaceCollision(ParticleBase *particle, const Vec3D& loc, 
                                         const Vec3D& vel, double energy, double theta,
                                         const Vec3D& normal, const Vec3D& tang1, const Vec3D& tang2) {
        // Calculate backscattering probability
        double etaBACK = ETA0(energy);
        double CT = CTback(energy, etaBACK);
        etaBACK = etaBACK * exp(CT * (1.0 - cos(theta)));
        
        // Calculate secondary electron emission probability
        double etaSEC = ETA1(energy);
        CT = CTsec(energy);
        etaSEC = etaSEC * exp(CT * (1.0 - cos(theta)));
        
        if (_debugprint) {
            logfile << "Electron surface collision @ E=" << energy*1e-3 << " keV, "
                    << "theta=" << theta*180.0/M_PI << " deg" << endl;
            logfile << "  etaBACK=" << etaBACK << ", etaSEC=" << etaSEC << endl;
        }
        
        double Elost = energy;
        
        // Check for backscattering
        double rand[1];
        _rng->get(rand);
        
        if (rand[0] <= etaBACK && Elost >= 1.0) {
            // Backscattered electron
            // Use analytical energy distribution (simplified - could use ebackdist_anal)
            // For now, use average backscattered energy fraction
            double Eout = 0.5 * Elost;  // Simplified - should use proper distribution
            if (Eout > Elost) Eout = Elost;
            if (Eout < 0.0) Eout = 0.0;
            
            Elost -= Eout;
            
            if (_debugprint) {
                logfile << "  Electron backscattered with E=" << Eout*1e-3 << " keV" << endl;
            }
            
                // Create backscattered electron
            createSecondaryParticle(particle, loc, normal, tang1, tang2, Eout, 
                                    -1.0, MASS_E/MASS_U, particle->gen() + SURFACE_GENERATION_OFFSET);  // Electron
        }
        
        // Generate secondary electrons
        generateSecondaryElectrons(particle, loc, normal, tang1, tang2, etaSEC, Elost);
    }

    /**
     * @brief Handle ion surface collision
     * @param particle Incident particle
     * @param loc Impact location
     * @param vel Incident velocity
     * @param energy Incident energy in eV
     * @param theta Impact angle with surface normal (radians)
     * @param normal Surface normal vector
     * @param tang1 First tangent vector
     * @param tang2 Second tangent vector
     */
    void handleIonSurfaceCollision(ParticleBase *particle, const Vec3D& loc,
                                    const Vec3D& vel, double energy, double theta,
                                    const Vec3D& normal, const Vec3D& tang1, const Vec3D& tang2) {
        // Determine mass ratio for scaling
        double mass_ratio = 1.0;
        double particle_mass = particle->m() / MASS_U;  // In atomic mass units
        
        if (particle_mass > 1.5 && particle_mass < 2.0) {
            // H or D
            if (particle_mass > 1.7) mass_ratio = 2.0;  // D
        } else if (particle_mass > 3.0 && particle_mass < 4.0) {
            // H2 or D2
            if (particle_mass > 3.5) mass_ratio = 4.0;  // D2
            else mass_ratio = 2.0;  // H2
        } else if (particle_mass > 3.9 && particle_mass < 4.1) {
            // He
            mass_ratio = 4.0;
        }
        
        double E0_keV = energy * 1e-3;  // Convert to keV
        
        // Calculate backscattering probability
        double mu_i = 0.5;  // Free parameter
        double etaBACK = RN(E0_keV, mass_ratio);
        etaBACK = etaBACK / ((1.0 - mu_i) * cos(theta) + mu_i);
        if (etaBACK > 1.0) etaBACK = 1.0;
        
        // Calculate secondary electron emission probability
        double mu_e = 1.45;
        double etaSEC = eta_ion(E0_keV, mass_ratio);
        
        // Scale for molecular ions (H2, D2, He)
        if (particle_mass > 3.0 && particle_mass < 4.1) {
            etaSEC = etaSEC * 2.90 / 1.32;  // See Redbook Vol 1, C-10
        }
        
        etaSEC = etaSEC * exp(mu_e * (1.0 - cos(theta)));
        
        if (_debugprint) {
            logfile << "Ion surface collision @ E=" << E0_keV << " keV, "
                    << "theta=" << theta*180.0/M_PI << " deg" << endl;
            logfile << "  etaBACK=" << etaBACK << ", etaSEC=" << etaSEC << endl;
        }
        
        double Elost = energy;
        
        // Check for backscattering
        double rand[1];
        _rng->get(rand);
        
        if (rand[0] <= etaBACK && Elost >= 1.0) {
            // Backscattered ion
            double RE_val = RE(E0_keV, mass_ratio);
            double RN_val = RN(E0_keV, mass_ratio);
            double Eout = (RE_val / RN_val) * Elost;
            
            if (Eout >= 0.0 && Eout <= Elost) {
                Elost -= Eout;
                
                if (_debugprint) {
                    logfile << "  Ion backscattered with E=" << Eout*1e-3 << " keV" << endl;
                }
                
                // Create backscattered ion (same type as incident)
                createSecondaryParticle(particle, loc, normal, tang1, tang2, Eout,
                                        particle->q(), particle->m()/MASS_U, 
                                        particle->gen() + SURFACE_GENERATION_OFFSET);
            }
        }
        
        // Generate secondary electrons
        generateSecondaryElectrons(particle, loc, normal, tang1, tang2, etaSEC, Elost);
    }

    /**
     * @brief Generate secondary electrons after surface impact
     * @param particle Incident particle
     * @param loc Impact location
     * @param normal Surface normal
     * @param tang1 First tangent
     * @param tang2 Second tangent
     * @param etaSEC Secondary electron yield
     * @param Elost Energy available for secondary emission (in eV)
     */
    void generateSecondaryElectrons(ParticleBase *particle, const Vec3D& loc,
                                     const Vec3D& normal, const Vec3D& tang1, const Vec3D& tang2,
                                     double etaSEC, double& Elost) {
        if (Elost < SECONDARY_ELECTRON_ENERGY) return;
        
        double rand[1];
        _rng->get(rand);
        
        // Check for first secondary electron
        if (rand[0] <= etaSEC && Elost >= SECONDARY_ELECTRON_ENERGY) {
            Elost -= SECONDARY_ELECTRON_ENERGY;
            
            if (_debugprint) {
                logfile << "  Secondary electron emitted, E=" << SECONDARY_ELECTRON_ENERGY << " eV" << endl;
            }
            
            createSecondaryParticle(particle, loc, normal, tang1, tang2, 
                                    SECONDARY_ELECTRON_ENERGY,
                                    -1.0, MASS_E/MASS_U, particle->gen() + SURFACE_GENERATION_OFFSET);
            
            // Check for additional secondary electrons if yield > 1
            if (etaSEC > 1.0) {
                double eta_remaining = etaSEC - 1.0;
                int cnt_lsec = 1;
                
                while (eta_remaining > 0.0 && Elost >= SECONDARY_ELECTRON_ENERGY && 
                       cnt_lsec < MAX_SECONDARY_ELECTRONS) {
                    _rng->get(rand);
                    
                    if (rand[0] <= eta_remaining) {
                        Elost -= SECONDARY_ELECTRON_ENERGY;
                        cnt_lsec++;
                        
                        if (_debugprint) {
                            logfile << "  Additional secondary electron #" << cnt_lsec 
                                    << ", E=" << SECONDARY_ELECTRON_ENERGY << " eV" << endl;
                        }
                        
                        createSecondaryParticle(particle, loc, normal, tang1, tang2,
                                                SECONDARY_ELECTRON_ENERGY,
                                                -1.0, MASS_E/MASS_U, particle->gen() + SURFACE_GENERATION_OFFSET);
                    }
                    
                    eta_remaining -= 1.0;
                }
                
                if (cnt_lsec >= MAX_SECONDARY_ELECTRONS) {
                    logfile << "WARNING: Maximum secondary electron count reached" << endl;
                }
            }
        }
    }

    /**
     * @brief Create a secondary particle with proper angular distribution
     * Based on make_sec_wall_part Fortran routine
     * @param particle Incident particle (for current/charge info)
     * @param loc Impact location
     * @param normal Surface normal
     * @param tang1 First tangent
     * @param tang2 Second tangent
     * @param Ek Secondary particle energy in eV
     * @param charge Secondary particle charge
     * @param mass Secondary particle mass (in atomic mass units)
     * @param generation Generation number
     */
    void createSecondaryParticle(ParticleBase *particle, const Vec3D& loc,
                                 const Vec3D& normal, const Vec3D& tang1, const Vec3D& tang2,
                                 double Ek, double charge, double mass, int generation) {
        // Check generation depth using modulo 100 to get effective generation
        // This works for both volume-generated (gen+1) and surface-generated (gen+101) particles
        int effective_gen = generation % 100;
        if (effective_gen >= MAX_GENERATION_DEPTH) {
            if (_debugprint) logfile << "WARNING: Maximum generation depth reached (effective gen=" 
                              << effective_gen << ", full gen=" << generation << ")" << endl;
            return;
        }
        
        // Check particle count
        if (_pdb->size() > MAX_PARTICLE_COUNT) {
            logfile << "WARNING: Maximum particle count reached" << endl;
            return;
        }
        
        // Calculate velocity from energy (non-relativistic approximation for low energies)
        // For high energies, use relativistic formula
        double mc2 = mass * MASS_U * SPEED_C2;  // Rest energy in eV
        double gamma = 1.0 + Ek / mc2;
        double Vtot;
        
        if (gamma > 1.01) {
            // Relativistic case
            Vtot = sqrt(gamma * gamma - 1.0) / gamma * SPEED_C;
        } else {
            // Non-relativistic case
            Vtot = sqrt(2.0 * Ek * CHARGE_E / (mass * MASS_U));
        }
        
        // Validate velocity
        if (!std::isfinite(Vtot) || Vtot < 0.0 || Vtot > MAX_VELOCITY_SCALE) {
            if (debug) logfile << "WARNING: Invalid secondary particle velocity: " << Vtot << endl;
            return;
        }
        
        // Generate random emission angles
        // theta = asin(rnd) gives cosine distribution
        // phi = 2*pi*rnd gives uniform azimuthal distribution
        double rand[2];
        _rng->get(rand);
        double theta = asin(rand[0]);  // Emission angle from normal
        double phi = 2.0 * M_PI * rand[1];  // Azimuthal angle
        
        // Calculate velocity components in local coordinate system (u, t, n)
        double V11 = sin(theta) * cos(phi);
        double V22 = sin(theta) * sin(phi);
        double V33 = cos(theta);
        
        // Transform to global coordinates
        Vec3D vel_global = V11 * normal + V22 * tang1 + V33 * tang2;
        vel_global *= Vtot;
        
        // Adjust location slightly off the surface
        Vec3D loc_adj = loc + 0.01 * _geom.h() * normal;
        
        if (_debugprint) {
            double emission_angle = acos(V33) * 180.0 / M_PI;
            logfile << "  Secondary particle: E=" << Ek << " eV, "
                    << "angle=" << emission_angle << " deg, "
                    << "charge=" << charge << ", mass=" << mass << endl;
        }
        
        // Create ParticleP3D: (t, x, vx, y, vy, z, vz) - velocity components directly
        ParticleP3D sec_part(0.0,  // time
                            loc_adj[0], vel_global[0],
                            loc_adj[1], vel_global[1],
                            loc_adj[2], vel_global[2]);
        
        // Add to particle database
        try {
            _pdb->add_particle(particle->IQ(), charge, mass, generation, sec_part);
        } catch (const std::exception& e) {
            logfile << "ERROR: Failed to add secondary particle: " << e.what() << endl;
        }
    }

public:
    /**
     * @brief Constructor
     * @param geom Geometry object for surface normal calculation
     * @param pdb Particle database for adding secondary particles
     * @param mass Reference mass for particle type identification
     * @param debugprint Enable debug output
     */
    THCallback_surf_EAMCC(Geometry &geom,
                          ParticleDataBase3D* pdb,
                          double mass = 1.0,
                          bool debugprint = false,
                          double minimum_z = 7.0e-3)
        : _geom(geom), _pdb(pdb), _mass(mass), _minimum_z(minimum_z), _debugprint(debugprint) {
        _rng = new MTRandom(1);
        double qx[1];
        _rng->get(qx);  // Initialize RNG
        
        if (debug) {
            logfile << "THCallback_surf_EAMCC initialized" << endl;
        }
    }

    virtual ~THCallback_surf_EAMCC() {
        if (_rng) delete _rng;
    }

    /**
     * @brief Main callback function called when particle trajectory ends
     * @param particle Particle that hit the surface
     * @param pdb Particle database
     */
    virtual void operator()(ParticleBase *particle, class ParticleDataBase *pdb) {
        Particle3D *p3d = (Particle3D *)(particle);
        Vec3D loc = p3d->location();
        Vec3D vel = p3d->velocity();
        
        // Validate position
        if (!validatePosition(loc)) {
            return;
        }
        
        // Only generate secondary particles beyond the configured impact plane.
        if (loc[2] <= _minimum_z) {
            return;
        }
        
        // Limit generation depth to prevent infinite chains
        // Check effective generation (works for both volume and surface-generated particles)
        int generation = particle->gen();
        int effective_gen = generation % 100;
        if (effective_gen >= MAX_GENERATION_DEPTH) {
            if (_debugprint) {
                logfile << "WARNING: Maximum generation depth reached (gen=" << generation 
                        << ", effective=" << effective_gen << "), stopping surface secondary generation" << endl;
            }
            return;
        }
        
        // Check particle count to prevent memory overflow
        if (_pdb->size() > MAX_PARTICLE_COUNT) {
            logfile << "WARNING: Maximum particle count reached (" << _pdb->size() 
                    << "), stopping surface secondary generation" << endl;
            return;
        }
        
        // Calculate energy
        double vel_mag = vel.ssqr();
        if (vel_mag < 0.0 || !std::isfinite(vel_mag)) {
            if (debug) logfile << "WARNING: Invalid velocity magnitude" << endl;
            return;
        }
        vel_mag = sqrt(vel_mag);
        double energy = particle->m() * SPEED_C2 * (1.0 / sqrt(1.0 - (vel_mag * vel_mag) / (SPEED_C2)) - 1.0) / CHARGE_E;  // in eV
        
        if (!std::isfinite(energy) || energy < 0.0) {
            if (debug) logfile << "WARNING: Invalid particle energy: " << energy << " eV" << endl;
            return;
        }
        
        // Get surface normal
        Vec3D normal = _geom.surface_normal(loc);
        if (normal.ssqr() < 1e-10) {
            if (debug) logfile << "WARNING: Could not determine surface normal" << endl;
            return;
        }
        normal.normalize();
        
        // Calculate impact angle
        Vec3D vel_normalized = vel;
        vel_normalized.normalize();
        double cos_theta = -(vel_normalized[0]*normal[0] + vel_normalized[1]*normal[1] + vel_normalized[2]*normal[2]);  // Negative because particle is coming in
        if (cos_theta < 0.0) cos_theta = 0.0;  // Clamp to [0, pi/2]
        if (cos_theta > 1.0) cos_theta = 1.0;
        double theta = acos(cos_theta);
        
        // Create tangent vectors
        Vec3D tang1 = normal.arb_perpendicular();
        Vec3D tang2 = cross(normal, tang1);
        tang1.normalize();
        tang2.normalize();
        
        // Determine particle type and handle collision
        double particle_mass = particle->m() / MASS_U;
        double particle_charge = particle->q() / CHARGE_E;
        
        if (particle_mass < 1e-3 && particle_charge < 0.0) {
            // Electron
            handleElectronSurfaceCollision(particle, loc, vel, energy, theta, normal, tang1, tang2);
        } else if (particle_mass > 0.5) {
            // Ion
            handleIonSurfaceCollision(particle, loc, vel, energy, theta, normal, tang1, tang2);
        }
        
        // Mark original particle as absorbed (or could be backscattered, handled above)
        // The particle status will be set by the trajectory handler
    }
};

#endif /* THCALLBACK_SURF_EAMCC_H_ */

