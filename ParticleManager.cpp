/*
 * ParticleManager.cpp
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored from ManageSimulation)
 */

#include "ParticleManager.h"
#include "SimulationParameters.h"
#include "globals.h"
#include "my_diagnostics.h"

#include "particledatabase.hpp"
#include "geometry.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "funct.h"
#include "particles.hpp"

#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <map>

using namespace std;

ParticleManager::ParticleManager() : 
    particles(nullptr), collisions(nullptr), part_kind(nullptr), 
    part_gen(nullptr), rgn(nullptr) {
    initializeDiagnostics();
}

ParticleManager::~ParticleManager() {
    resetParticles();
}

void ParticleManager::resetParticles() {
    delete particles;
    particles = nullptr;
    
    for (auto* pdb : particles_species) {
        delete pdb;
    }
    particles_species.clear();
    particles_kinds.clear();
    part_number.clear();
    
    delete rgn;
    rgn = nullptr;
}

void ParticleManager::setParticles(ParticleDataBase3D* pdb) {
    if (particles != nullptr && particles != pdb) {
        delete particles;
    }
    particles = pdb;
}

void ParticleManager::initializeDiagnostics() {
    diagnostics = {
        DIAG_X, DIAG_Y, DIAG_Z,
        DIAG_VX, DIAG_VY, DIAG_VZ,
        DIAG_CURR, DIAG_MASS, DIAG_CHARGE, DIAG_NO
    };
}

void ParticleManager::defineParticleDatabase(const Geometry& geometry, const SimulationParameters& params) {
    ParticleDataBase3D* pdb = new ParticleDataBase3D(geometry);
    pdb->set_max_steps(3000);
    bool pmirror[6] = {false, false, false, false, false, false};
    pdb->set_mirror(pmirror);
    pdb->set_surface_collision(true);
    
    // Beam definition - add rectangular beam with energy
    pdb->clear();
    
    // Get geometry boundaries for beam positioning
    Vec3D ori = geometry.origo();
    Vec3D fin = geometry.max();
    double x_start = (ori[0] > -0.01) ? ori[0] : -0.01;
    double y_start = (ori[1] > -0.01) ? ori[1] : -0.01;
    double z_start = ori[2];
    double x_end = (fin[0] < 0.01) ? fin[0] : 0.01;
    double y_end = (fin[1] < 0.01) ? fin[1] : 0.01;

    // DEBUGGING: Uncomment to adjust beam width
    // x_start = -0.002; x_end = 0.002; // Adjusted for beam width
    // y_start = -0.002; y_end = 0.002; // Adjusted for beam width
    
    // Add rectangular beam with energy using parameters from scenario file
	pdb->clear();
    pdb->add_rectangular_beam_with_energy( 
        params.getNParticles(),        // N_PARTICLES
        params.getJIon(),             // J_ION [A/m^2]
        params.getQIons(),            // Q_IONS [e]
        params.getMIons(),            // M_IONS [amu]
        params.getE0Z(),              // E0_Z [eV] - axial energy
        params.getTPar(),             // TPAR [eV] - parallel temperature
        params.getTPerp(),            // TPERP [eV] - perpendicular temperature
        Vec3D(0.5*(x_start+x_end), 0.5*(y_start+y_end), z_start), // beam center
        Vec3D(1,0,0),                 // x direction vector
        Vec3D(0,1,0),                 // y direction vector
        0.5*(x_end-x_start),          // beam width in x
        0.5*(y_end-y_start)           // beam width in y
    );
    
    particles = pdb;
    
    ibsimu.message(1) << "Particle database created with " << params.getNParticles() 
                      << " particles, J=" << params.getJIon() << " A/m²" << endl;
}

void ParticleManager::defineParticleDatabase(const Geometry& geometry, const SimulationParameters& params,
                                            const string& emittername, bool backstream) {
    ParticleDataBase3D* pdb = new ParticleDataBase3D(geometry);
    pdb->set_max_steps(3000);
    bool pmirror[6] = {false, false, false, false, false, false};
    pdb->set_mirror(pmirror);
    pdb->set_surface_collision(true);
    
    // Beam definition moved outside the main iteration cycle
    // There is no need to change the starting particles at every iteration.
    pdb->clear();
    
    // Load emitter from file
    string emitter_file = emittername;
    loadEmitter(pdb, emitter_file, 1.0, 1.0, 0.0); // Forward particles with original charge and mass
    
    if (backstream) {
        // Load backstreaming particles to compensate space charge
        loadEmitter(pdb, emitter_file, -1.0, 1.0, 1.0); // Reversed velocities, different mass factor for visualization
    }

    particles = pdb;
    
    ibsimu.message(1) << "Particle database created with " << pdb->size() << " particles from emitter: " << emittername << endl;
}

void ParticleManager::fillParticleDatabases(const SimulationParameters& params) {
    // create pdb for each species
    if (params.getIncludeStripping() == 2) {
        ParticleDataBase3D* pdbHm = new ParticleDataBase3D(*particles);
        ParticleDataBase3D* pdbH0 = new ParticleDataBase3D(*particles);
        ParticleDataBase3D* pdbHp = new ParticleDataBase3D(*particles);
        ParticleDataBase3D* pdbH2p = new ParticleDataBase3D(*particles);
        ParticleDataBase3D* pdbH20 = new ParticleDataBase3D(*particles);
        ParticleDataBase3D* pdbe = new ParticleDataBase3D(*particles);

        *pdbHm = *particles;
        *pdbH0 = *particles;
        *pdbHp = *particles;
        *pdbH2p = *particles;
        *pdbH20 = *particles;
        *pdbe = *particles;

        if (debug) logfile << "\nDEBUG: SPECIES DATABASES CREATED. NOW CLEANING... \n" << flush;
        
        // Initialize counters for each species
        size_t nHm = 0, nH0 = 0, nHp = 0, nH2p = 0, nH20 = 0, ne = 0;
        
        // Filter particles by species - iterate through all particles and identify species
        for (size_t id = 0; id < particles->size(); id++) {
            Particle3D& pp = particles->particle(id);
            
            // Identify particle species based on mass and charge
            particle_kind species = identify_particle_species(pp.m(), pp.q(), params.getMIons());
            
            // Reset trajectory for particles that don't belong to each species database
            if (species != PARTICLE_HM) pdbHm->reset_trajectory(id);
            else nHm++;
            
            if (species != PARTICLE_H0) pdbH0->reset_trajectory(id);
            else nH0++;
            
            if (species != PARTICLE_HP) pdbHp->reset_trajectory(id);
            else nHp++;
            
            if (species != PARTICLE_H2P) pdbH2p->reset_trajectory(id);
            else nH2p++;
            
            if (species != PARTICLE_H20) pdbH20->reset_trajectory(id);
            else nH20++;
            
            if (species != PARTICLE_E) pdbe->reset_trajectory(id);
            else ne++;
        }
        
        // Store particle count for each species
        part_number.push_back(nHm);
        part_number.push_back(nH0);
        part_number.push_back(nHp);
        part_number.push_back(nH2p);
        part_number.push_back(nH20);
        part_number.push_back(ne);
        
        // Add species databases to collections
        particles_species.push_back(pdbHm); particles_kinds.push_back(PARTICLE_HM);
        particles_species.push_back(pdbH0); particles_kinds.push_back(PARTICLE_H0);
        particles_species.push_back(pdbHp); particles_kinds.push_back(PARTICLE_HP);
        particles_species.push_back(pdbH2p); particles_kinds.push_back(PARTICLE_H2P);
        particles_species.push_back(pdbH20); particles_kinds.push_back(PARTICLE_H20);
        particles_species.push_back(pdbe); particles_kinds.push_back(PARTICLE_E);
        
        logfile << "\nSPECIES DATABASES CREATED! \n" << flush;
        logfile << " #HM: " << nHm << endl << flush;
        logfile << " #H0: " << nH0 << endl << flush;
        logfile << " #HP: " << nHp << endl << flush;
        logfile << " #H2P: " << nH2p << endl << flush;
        logfile << " #H20: " << nH20 << endl << flush;
        logfile << " #EL: " << ne << endl << flush;
        logfile << " SUM: " << ne+nH20+nH2p+nHp+nH0+nHm << " \n\n" << flush;
    }
    else {
        ParticleDataBase3D* pdbHm = new ParticleDataBase3D(*particles);
        *pdbHm = *particles;
        particles_species.push_back(pdbHm); 
        particles_kinds.push_back(PARTICLE_HM);
        size_t nHm = pdbHm->size();
        part_number.push_back(nHm);
    }
}

void ParticleManager::analyzeParticleEndLocations(const SimulationParameters& params) {
    if (!particles) {
        ibsimu.message(1) << "No particle database available for end location analysis" << endl;
        return;
    }
    
    ibsimu.message(1) << "Analyzing particle end locations..." << endl;
    
    try {
        // Initialize counters for different particle fates
        size_t particles_extracted = 0;
        size_t particles_out = 0;
        size_t particles_stripped = 0;
        size_t particles_collided = 0;
        size_t particles_ok = 0;
        
        // Initialize generation counters
        map<int, size_t> generation_counts;
        
        // Initialize species counters
        map<particle_kind, size_t> species_counts;
        
        // Initialize collision location analysis
        vector<double> collision_z_positions;
        vector<double> extraction_z_positions;
        
        // Energy analysis
        vector<double> final_energies;
        vector<double> initial_energies;
        
        // Current analysis
        double total_initial_current = 0.0;
        double total_final_current = 0.0;
        double extracted_current = 0.0;
        
        // Analyze each particle
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D& particle = particles->particle(i);
            
            // Get particle properties
            double mass = particle.m();
            double charge = particle.q();
            double current = particle.IQ();
            int generation = particle.gen() % 100;  // Extract generation number (0-4+)
            particle_status_e status = particle.get_status();
            
            // Count by generation
            generation_counts[generation]++;
            
            // Identify particle species
            particle_kind species = identify_particle_species(mass, charge, params.getMIons());
            species_counts[species]++;
            
            // Add to current totals
            total_initial_current += abs(current);
            
            // Get trajectory information
            if (particle.traj_size() > 0) {
                // Get initial and final trajectory points
                ParticleP3D initial_point = particles->trajectory_point(i, 0);
                ParticleP3D final_point = particles->trajectory_point(i, particle.traj_size() - 1);
                
                // Calculate initial and final energies
                double initial_vx = initial_point(2);
                double initial_vy = initial_point(4);
                double initial_vz = initial_point(6);
                double initial_v_mag = sqrt(initial_vx*initial_vx + initial_vy*initial_vy + initial_vz*initial_vz);
                // Relativistic kinetic energy using constants c and c2 from constants.hpp
                // KE = (gamma - 1) * m * c^2, gamma = 1 / sqrt(1 - (v^2/c^2))
                double v2 = initial_v_mag * initial_v_mag;
                double beta2 = v2 / SPEED_C2;
                if (beta2 >= 1.0) beta2 = 0.999999999999; // clamp to avoid nan
                double gamma = 1.0 / sqrt(1.0 - beta2);
                double initial_energy_eV = (gamma - 1.0) * mass * SPEED_C2 / CHARGE_E;
                initial_energies.push_back(initial_energy_eV);
                
                double final_vx = final_point(2);
                double final_vy = final_point(4);
                double final_vz = final_point(6);
                double final_v_mag = sqrt(final_vx*final_vx + final_vy*final_vy + final_vz*final_vz);
                double final_energy_eV = 0.5 * mass * final_v_mag * final_v_mag / CHARGE_E;
                final_energies.push_back(final_energy_eV);
                
                // Get final position
                double final_z = final_point(5);
                
                // Analyze particle fate based on status
                switch (status) {
                    case PARTICLE_OK:
                        particles_ok++;
                        total_final_current += abs(current);
                        // Check if particle reached extraction plane
                        if (final_z >= 0.009) {  // EG exit at z=9mm
                            particles_extracted++;
                            extracted_current += abs(current);
                            extraction_z_positions.push_back(final_z);
                        }
                        break;
                        
                    case PARTICLE_COLL:
                        particles_collided++;
                        collision_z_positions.push_back(final_z);
                        break;
                        
                    case PARTICLE_STRIP:
                        particles_stripped++;
                        collision_z_positions.push_back(final_z);
                        break;
                        
                    case PARTICLE_OUT:
                        particles_out++;
                        break;
                        
                    default:
                        break;
                }
            }
        }
        
        // Calculate statistics
        size_t total_particles = particles->size();
        double extraction_efficiency = (total_particles > 0) ? 
            (double)particles_extracted / total_particles * 100.0 : 0.0;
        double stripping_rate = (total_particles > 0) ? 
            (double)particles_stripped / total_particles * 100.0 : 0.0;
        double collision_rate = (total_particles > 0) ? 
            (double)particles_collided / total_particles * 100.0 : 0.0;
        double out_rate = (total_particles > 0) ? 
            (double)particles_out / total_particles * 100.0 : 0.0;
        
        // Calculate current transmission
        double current_transmission = (total_initial_current > 0) ? 
            extracted_current / total_initial_current * 100.0 : 0.0;
        
        // Log detailed analysis results
        logfile << "\n=== PARTICLE END LOCATION ANALYSIS ===" << endl;
        logfile << "Total particles analyzed: " << total_particles << endl;
        logfile << "\nParticle Fates:" << endl;
        logfile << "  Extracted (OK at z>=9mm): " << particles_extracted 
                << " (" << extraction_efficiency << "%)" << endl;
        logfile << "  Still OK (not extracted): " << (particles_ok - particles_extracted) << endl;
        logfile << "  Collided with surfaces: " << particles_collided 
                << " (" << collision_rate << "%)" << endl;
        logfile << "  Stripped: " << particles_stripped 
                << " (" << stripping_rate << "%)" << endl;
        logfile << "  Out: " << particles_out 
                << " (" << out_rate << "%)" << endl;

        logfile << "\nCurrent Analysis:" << endl;
        logfile << "  Initial total current: " << total_initial_current << " A" << endl;
        logfile << "  Final total current: " << total_final_current << " A" << endl;
        logfile << "  Extracted current: " << extracted_current << " A" << endl;
        logfile << "  Current transmission: " << current_transmission << "%" << endl;
        
        // Generation analysis
        logfile << "\nGeneration Analysis:" << endl;
        for (const auto& gen_pair : generation_counts) {
            logfile << "  Generation " << gen_pair.first << ": " << gen_pair.second << " particles" << endl;
        }
        
        // Species analysis
        logfile << "\nSpecies Analysis:" << endl;
        for (const auto& species_pair : species_counts) {
            string species_name = get_particle_name(species_pair.first);
            logfile << "  " << species_name << ": " << species_pair.second << " particles" << endl;
        }
        
        // Collision location statistics
        if (!collision_z_positions.empty()) {
            sort(collision_z_positions.begin(), collision_z_positions.end());
            double min_collision_z = collision_z_positions.front();
            double max_collision_z = collision_z_positions.back();
            double median_collision_z = collision_z_positions[collision_z_positions.size() / 2];
            
            logfile << "\nCollision Location Statistics:" << endl;
            logfile << "  First collision at z = " << min_collision_z * 1000 << " mm" << endl;
            logfile << "  Last collision at z = " << max_collision_z * 1000 << " mm" << endl;
            logfile << "  Median collision at z = " << median_collision_z * 1000 << " mm" << endl;
        }
        
        // Extraction location statistics
        if (!extraction_z_positions.empty()) {
            sort(extraction_z_positions.begin(), extraction_z_positions.end());
            double min_extraction_z = extraction_z_positions.front();
            double max_extraction_z = extraction_z_positions.back();
            
            logfile << "\nExtraction Statistics:" << endl;
            logfile << "  First extraction at z = " << min_extraction_z * 1000 << " mm" << endl;
            logfile << "  Last extraction at z = " << max_extraction_z * 1000 << " mm" << endl;
        }
        
        // Energy statistics
        if (!initial_energies.empty() && !final_energies.empty()) {
            double avg_initial_energy = accumulate(initial_energies.begin(), initial_energies.end(), 0.0) / initial_energies.size();
            double avg_final_energy = accumulate(final_energies.begin(), final_energies.end(), 0.0) / final_energies.size();
            
            logfile << "\nEnergy Statistics:" << endl;
            logfile << "  Average initial energy: " << avg_initial_energy << " eV" << endl;
            logfile << "  Average final energy: " << avg_final_energy << " eV" << endl;
            logfile << "  Average energy gain: " << (avg_final_energy - avg_initial_energy) << " eV" << endl;
        }
        
        // Performance metrics for optimization
        logfile << "\nPerformance Metrics:" << endl;
        logfile << "  Extraction efficiency: " << extraction_efficiency << "%" << endl;
        logfile << "  Current transmission: " << current_transmission << "%" << endl;
        logfile << "  Beam losses: " << (collision_rate + out_rate) << "%" << endl;
        
        // IBSimu console output (summary)
        ibsimu.message(1) << "Analysis complete: " << particles_extracted << "/" << total_particles 
                          << " particles extracted (" << extraction_efficiency << "%)" << endl;
        ibsimu.message(1) << "Current transmission: " << current_transmission << "%" << endl;
        
        if (params.getIncludeStripping() > 0) {
            ibsimu.message(1) << "Stripping events: " << particles_stripped 
                              << " (" << stripping_rate << "%)" << endl;
        }
        
        logfile << "=== END PARTICLE ANALYSIS ===\n" << endl;
        
    } catch (const Error& e) {
        logfile << "ERROR in particle end location analysis: IBSimu Error occurred" << endl;
        ibsimu.message(1) << "Error during particle analysis" << endl;
    } catch (const std::exception& e) {
        logfile << "ERROR in particle end location analysis: " << e.what() << endl;
        ibsimu.message(1) << "Error during particle analysis: " << e.what() << endl;
    }
}

string ParticleManager::saveEmitter(const string& emitname, double zloc, const string& output_folder) {
    string filename = output_folder + emitname;
    
    if (!particles) {
        throw Error(ERROR_LOCATION, "No particle database available to save emitter");
    }
    
    ibsimu.message(1) << "Saving particles at z = " << zloc << " m to emitter file: " << filename << endl;
    
    try {
        // Create trajectory diagnostic data at the specified z location
        TrajectoryDiagnosticData tdata;
        
        // Get trajectories at the specified plane
        particles->trajectories_at_plane(tdata, AXIS_Z, zloc, diagnostics);
        
        // Write emitter file
        ofstream outfile(filename);
        if (!outfile.is_open()) {
            throw Error(ERROR_LOCATION, "Could not open file for writing: " + filename);
        }
        
        // Write header with detailed information
        outfile << "# Emitter file saved at z = " << zloc << " m" << endl;
        outfile << "# Format: x[m] y[m] z[m] vx[m/s] vy[m/s] vz[m/s] I[A] m[kg] q[C]" << endl;
        outfile << "# Number of particles: " << tdata.traj_size() << endl;
        
        if (tdata.traj_size() > 0) {
            for (uint32_t i = 0; i < tdata.traj_size(); i++) {
                // Write particle data: x, y, z, vx, vy, vz, current, mass, charge
                outfile << fixed << setprecision(6)
                        << tdata(i, 0) << " "  // x position
                        << tdata(i, 1) << " "  // y position  
                        << tdata(i, 2) << " "  // z position
                        << tdata(i, 3) << " "  // vx velocity
                        << tdata(i, 4) << " "  // vy velocity
                        << tdata(i, 5) << " "  // vz velocity
                        << tdata(i, 6) << " "  // current
                        << tdata(i, 7) << " "  // mass
                        << tdata(i, 8) << endl; // charge
            }
            
            ibsimu.message(1) << "Saved " << tdata.traj_size() << " particles to emitter file" << endl;
        } else {
            ibsimu.message(1) << "Warning: No particles found at z = " << zloc << " m" << endl;
        }
        
        outfile.close();
        
    } catch (const Error& e) {
        logfile << "ERROR in saveEmitter: IBSimu Error occurred" << endl;
        ibsimu.message(1) << "Error saving emitter file" << endl;
        throw;
    } catch (const std::exception& e) {
        logfile << "ERROR in saveEmitter: " << e.what() << endl;
        ibsimu.message(1) << "Error saving emitter file: " << e.what() << endl;
        throw Error(ERROR_LOCATION, "Error saving emitter: " + string(e.what()));
    }
    
    ibsimu.message(1) << "Emitter saved to: " << filename << endl;
    return filename;
}

bool ParticleManager::checkEGExtractedCurrent(double oriJ, double& extsimJ, const SimulationParameters& params) {
    // Implementation for checking extracted current density at EG (z=0.009m)
    // Based on the check_EGext function in ManageSimulation.cpp
    
    if (!particles) {
        extsimJ = 0.0;
        return false;
    }
    
    try {
        // Calculate current at EG exit location (z = 0.009m)
        double eg_exit_z = 0.009;  // EG exit location in meters
        
        // Create trajectory diagnostic data at EG exit
        TrajectoryDiagnosticData tdata;
        vector<trajectory_diagnostic_e> diagnostics = {
            DIAG_X, DIAG_Y, DIAG_Z, DIAG_VX, DIAG_VY, DIAG_VZ,
            DIAG_CURR, DIAG_MASS, DIAG_CHARGE, DIAG_NO
        };
        
        // Get trajectories at EG exit plane
        particles->trajectories_at_plane(tdata, AXIS_Z, eg_exit_z, diagnostics);
        
        // Calculate total current at EG exit
        double total_current = 0.0;
        size_t particle_count = tdata.traj_size();
        
        for (size_t i = 0; i < particle_count; ++i) {
            total_current += tdata(i, 6); // DIAG_CURR is column 6
        }
        
        // Calculate current density using PG aperture area (circle with 7mm radius)
        // Area = π × (7mm)² = π × (7e-3)²
        double pg_aperture_area = M_PI * 7e-3 * 7e-3;
        double current_density = abs(total_current) / pg_aperture_area;
        
        extsimJ = current_density;
        
        // Check if within tolerance
        double tolerance = params.getJTolerance();
        double error = abs(extsimJ - oriJ) / oriJ;
        bool within_tolerance = (error < tolerance);
        
        if (debug) {
            logfile << "DEBUG: EG extracted current calculation:" << endl;
            logfile << "  Particles at EG exit (z=" << eg_exit_z << "m): " << particle_count << endl;
            logfile << "  Total current: " << total_current << " A" << endl;
            logfile << "  PG aperture area: " << pg_aperture_area << " m²" << endl;
            logfile << "  Extracted current density: " << extsimJ << " A/m²" << endl;
            logfile << "  Target current density: " << oriJ << " A/m²" << endl;
            logfile << "  Error: " << error*100 << "%" << endl;
            logfile << "  Tolerance: " << tolerance*100 << "%" << endl;
            logfile << "  Within tolerance: " << (within_tolerance ? "YES" : "NO") << endl;
        }
        
        ibsimu.message(1) << "Current check: target=" << oriJ << " A/m², extracted=" << extsimJ 
                          << " A/m², error=" << error*100 << "%, tolerance=" << tolerance*100 << "%" << endl;
        
        return within_tolerance;
        
    } catch (const Error& e) {
        if (debug) logfile << "DEBUG: Error in EG extracted current calculation: IBSimu Error occurred" << endl;
        extsimJ = 0.0;
        return false;
    } catch (const std::exception& e) {
        if (debug) logfile << "DEBUG: Error in EG extracted current calculation: " << e.what() << endl;
        extsimJ = 0.0;
        return false;
    }
}

void ParticleManager::exportTrajectoriesToVTK(const string& filename) {
    if (!particles) {
        ibsimu.message(1) << "No particle database available for trajectories VTK export" << endl;
        return;
    }
    
    string vtkfile = filename + "_trajectories.vtk";
    
    try {
        ofstream file(vtkfile, ios::binary);
        if (!file.is_open()) {
            throw Error("Cannot open VTK file for writing: " + vtkfile);
        }
        
        // Count total trajectory points
        int total_points = 0;
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            total_points += particle.traj_size();
        }
        
        // Debug output
        ibsimu.message(1) << "Exporting " << particles->size() << " particle trajectories in binary format" << endl;
        ibsimu.message(1) << "Total trajectory points: " << total_points << endl;
        if (particles->size() > 0) {
            Particle3D &first_particle = particles->particle(0);
            if (first_particle.traj_size() > 0) {
                ibsimu.message(1) << "First particle: " << first_particle.traj_size() << " points" << endl;
                ParticleP3D start_point = particles->trajectory_point(0, 0);
                ParticleP3D end_point = particles->trajectory_point(0, first_particle.traj_size()-1);
                ibsimu.message(1) << "  Start: (" << start_point(1) << ", " << start_point(3) << ", " << start_point(5) << ")" << endl;
                ibsimu.message(1) << "  End:   (" << end_point(1) << ", " << end_point(3) << ", " << end_point(5) << ")" << endl;
            }
        }

        if (total_points == 0) {
            ibsimu.message(1) << "No trajectory points to export" << endl;
            return;
        }
        
        // Helper lambda to swap endianness for binary VTK (big-endian required)
        auto swap_endian_float = [](float value) -> float {
            union { float f; uint32_t i; } u;
            u.f = value;
            u.i = ((u.i >> 24) & 0x000000FF) | 
                  ((u.i >> 8)  & 0x0000FF00) | 
                  ((u.i << 8)  & 0x00FF0000) | 
                  ((u.i << 24) & 0xFF000000);
            return u.f;
        };
        
        auto swap_endian_int = [](int value) -> int {
            return ((value >> 24) & 0x000000FF) | 
                   ((value >> 8)  & 0x0000FF00) | 
                   ((value << 8)  & 0x00FF0000) | 
                   ((value << 24) & 0xFF000000);
        };
        
        // Write VTK header (ASCII)
        file << "# vtk DataFile Version 3.0\n";
        file << "IBSIMU Particle Trajectories Export\n";
        file << "BINARY\n";
        file << "DATASET POLYDATA\n";
        
        // Write points header (ASCII) then data (binary)
        file << "POINTS " << total_points << " float\n";
        
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            for (size_t j = 0; j < particle.traj_size(); j++) {
                ParticleP3D traj_point = particles->trajectory_point(i, j);
                // Write x, y, z coordinates as big-endian floats
                float x = static_cast<float>(traj_point(1));
                float y = static_cast<float>(traj_point(3));
                float z = static_cast<float>(traj_point(5));
                float x_be = swap_endian_float(x);
                float y_be = swap_endian_float(y);
                float z_be = swap_endian_float(z);
                file.write(reinterpret_cast<const char*>(&x_be), sizeof(float));
                file.write(reinterpret_cast<const char*>(&y_be), sizeof(float));
                file.write(reinterpret_cast<const char*>(&z_be), sizeof(float));
            }
        }
        
        // Write lines (trajectory polylines)
        int total_lines = particles->size();
        int total_line_data = 0;
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            total_line_data += particle.traj_size() + 1; // +1 for line length
        }
        
        file << "\nLINES " << total_lines << " " << total_line_data << "\n";
        
        int point_offset = 0;
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            int traj_size = particle.traj_size();
            // Write number of points in this line
            int traj_size_be = swap_endian_int(traj_size);
            file.write(reinterpret_cast<const char*>(&traj_size_be), sizeof(int));
            // Write point indices
            for (int j = 0; j < traj_size; j++) {
                int idx_be = swap_endian_int(point_offset + j);
                file.write(reinterpret_cast<const char*>(&idx_be), sizeof(int));
            }
            point_offset += traj_size;
        }
        
        // Write point data - particle velocities and energies
        file << "\nPOINT_DATA " << total_points << "\n";
        
        // Export velocity magnitude in km/s for better visualization
        file << "SCALARS velocity_magnitude_km_s float\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            for (size_t j = 0; j < particle.traj_size(); j++) {
                ParticleP3D traj_point = particles->trajectory_point(i, j);
                // Velocity components: vx=index 2, vy=index 4, vz=index 6 in m/s
                double vx = traj_point(2);
                double vy = traj_point(4);
                double vz = traj_point(6);
                double vel_mag = sqrt(vx*vx + vy*vy + vz*vz);
                // Convert to km/s for better scale
                float vel_km_s = static_cast<float>(vel_mag / 1000.0);
                float vel_km_s_be = swap_endian_float(vel_km_s);
                file.write(reinterpret_cast<const char*>(&vel_km_s_be), sizeof(float));
            }
        }
        
        // Export time in nanoseconds for better scale
        file << "\nSCALARS time_ns float\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            for (size_t j = 0; j < particle.traj_size(); j++) {
                ParticleP3D traj_point = particles->trajectory_point(i, j);
                // Time is at index 0 in seconds
                double time = traj_point(0);
                // Convert to nanoseconds for better scale
                float time_ns = static_cast<float>(time * 1e9);
                float time_ns_be = swap_endian_float(time_ns);
                file.write(reinterpret_cast<const char*>(&time_ns_be), sizeof(float));
            }
        }
        
        // Export relativistic kinetic energy in eV
        file << "\nSCALARS kinetic_energy_eV float\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            double mass = particle.m(); // mass in kg
            for (size_t j = 0; j < particle.traj_size(); j++) {
                ParticleP3D traj_point = particles->trajectory_point(i, j);
                // Velocity components: vx=index 2, vy=index 4, vz=index 6 in m/s
                double vx = traj_point(2);
                double vy = traj_point(4);
                double vz = traj_point(6);
                double v_mag = sqrt(vx*vx + vy*vy + vz*vz);
                
                // Relativistic kinetic energy calculation
                const double c = SPEED_C; // speed of light in m/s
                double beta = v_mag / SPEED_C;
                double beta2 = v_mag*v_mag / SPEED_C2;
                double gamma = 1.0 / sqrt(1.0 - beta*beta);
                
                // Relativistic kinetic energy: KE = (γ - 1)mc²
                double kinetic_energy_J = (gamma - 1.0) * mass * c * c;
                
                // Convert to eV (1 eV = 1.602176634e-19 J)
                float kinetic_energy_eV = static_cast<float>(kinetic_energy_J / CHARGE_E);
                float kinetic_energy_eV_be = swap_endian_float(kinetic_energy_eV);
                file.write(reinterpret_cast<const char*>(&kinetic_energy_eV_be), sizeof(float));
            }
        }
        
        // Export particle properties as cell data
        file << "\nCELL_DATA " << total_lines << "\n";
        
        // Export particle IDs
        file << "SCALARS particle_id int\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            int id_be = swap_endian_int(static_cast<int>(i));
            file.write(reinterpret_cast<const char*>(&id_be), sizeof(int));
        }
        
        // Export particle current in Amperes [A]
        file << "\nSCALARS current_A float\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            float current_A = static_cast<float>(particle.IQ());
            float current_A_be = swap_endian_float(current_A);
            file.write(reinterpret_cast<const char*>(&current_A_be), sizeof(float));
        }
        
        // Export particle charge in elementary charge units [e]
        file << "\nSCALARS charge_e float\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            // Convert from Coulombs to elementary charge units
            float charge_e = static_cast<float>(particle.q() / CHARGE_E);
            float charge_e_be = swap_endian_float(charge_e);
            file.write(reinterpret_cast<const char*>(&charge_e_be), sizeof(float));
        }
        
        // Export particle mass in atomic mass units [amu]
        file << "\nSCALARS mass_amu float\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            // Convert from kg to atomic mass units
            float mass_amu = static_cast<float>(particle.m() / MASS_U);
            float mass_amu_be = swap_endian_float(mass_amu);
            file.write(reinterpret_cast<const char*>(&mass_amu_be), sizeof(float));
        }

        // Export particle status as cell data
        file << "\nSCALARS particle_status int\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            int status_be = swap_endian_int(static_cast<int>(particle.get_status()));
            file.write(reinterpret_cast<const char*>(&status_be), sizeof(int));
        }

        // Export particle generation as cell data
        file << "\nSCALARS particle_generation int\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            int gen_be = swap_endian_int(particle.gen());
            file.write(reinterpret_cast<const char*>(&gen_be), sizeof(int));
        }
        
        file.close();
        ibsimu.message(1) << "Particle trajectories exported to binary VTK file: " << vtkfile << endl;
        ibsimu.message(1) << "  Format: Binary VTK with big-endian byte order" << endl;
        ibsimu.message(1) << "  Point data: velocity, time, kinetic energy" << endl;
        ibsimu.message(1) << "  Cell data: particle ID, current, charge, mass, status, generation" << endl;
        
    } catch (const Error& e) {
        ibsimu.message(1) << "Error exporting trajectories to VTK" << endl;
    }
}
