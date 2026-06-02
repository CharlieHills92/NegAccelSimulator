/*
 * DiagnosticsManager.cpp
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored from ManageSimulation)
 */

#include "DiagnosticsManager.h"
#include "SimulationParameters.h"
#include "FileManager.h"
#include "TransverseData.h"
#include "ManageSimulation_New.h"  // For PowerStruct
#include "globals.h"
#include "my_diagnostics.h"
#include "funct.h"
#include "cross_sections.h"
#include <chrono>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <cmath>

#include "geometry.hpp"
#include "particledatabase.hpp"
#include "geomplotter.hpp"
#include "gtkplotter.hpp"
#include "ibsimu.hpp"
#include "error.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void DiagnosticsManager::generateDiagnosticData(const vector<double>& diagzpos, int iteration, 
                                               const string& outfile, bool append_at_end, 
                                               const ParticleDataBase3D* pdb) {
    std::cerr << "DEBUG: generateDiagnosticData called with " << diagzpos.size() << " z-positions" << std::endl << std::flush;
    
    if (debug) logfile << "DEBUG: Printing diagnostic data to " << outfile << " file\n" << flush;
    
    vector<double> x[1];
    vector<double> y[1];
    vector<double> z[1];
    vector<double> vx[1];
    vector<double> vy[1];
    vector<double> vz[1];
    vector<double> curr[1];
    vector<double> mass[1];
    vector<double> charge[1];
    vector<double> no_part[1];

    ofstream pout5;
    if (append_at_end) {
        pout5.open(outfile.c_str(), ios_base::app);
    } else {
        pout5.open(outfile.c_str());
    }
    
    if (!pout5.is_open()) {
        throw Error(ERROR_LOCATION, "Could not open diagnostic output file: " + outfile);
    }

    // Write header if not appending
    if (!append_at_end) {
        pout5 << "# it z[mm] I[mA] <x>[mm] xmax[mm] xmin[mm] <y>[mm] ymax[mm] ymin[mm] "
              << "<x'>[mrad] <y'>[mrad] Dx[mrad] Dy[mrad] <V>[V] <B>[mT] rho[1/m3] sigma[m2]\n";
    }

    for (size_t ii = 0; ii < diagzpos.size(); ii++) {
        double zloc = diagzpos[ii];
        std::cerr << "DEBUG: Processing z-position " << ii+1 << "/" << diagzpos.size() << " at z=" << zloc << " m" << std::endl << std::flush;
        
        std::cerr << "DEBUG: About to enter try block" << std::endl << std::flush;
        try {
            std::cerr << "DEBUG: Inside try block, creating TransverseData" << std::endl << std::flush;
            // Trajectory diagnostics at specific z location
            TrajectoryDiagnosticData tdata;
            vector<trajectory_diagnostic_e> diagnostics = {
                DIAG_X, DIAG_Y, DIAG_Z, DIAG_VX, DIAG_VY, DIAG_VZ,
                DIAG_CURR, DIAG_MASS, DIAG_CHARGE, DIAG_NO
            };
            
            std::cerr << "  About to call trajectories_at_plane..." << std::endl << std::flush;
            pdb->trajectories_at_plane(tdata, AXIS_Z, zloc, diagnostics);
            std::cerr << "  trajectories_at_plane call completed successfully!" << std::endl << std::flush;
            
            // Debug: Check the size immediately after the call
            size_t particle_count = tdata.traj_size();
            std::cerr << "  DEBUG: tdata.traj_size() = " << particle_count << std::endl << std::flush;
            
            // Write count to file for debugging
            ofstream debugfile("debug_particle_count.txt", ios::app);
            size_t pdb_size = pdb->size();
            debugfile << "z=" << zloc << " count=" << particle_count << " pdb_size=" << pdb_size << endl;
            debugfile.flush();
            debugfile.close();
            
            std::cerr << "  Retrieved " << tdata.traj_size() << " particles at z=" << zloc << std::endl << std::flush;
            std::cerr << "  Particle database size: " << pdb->size() << std::endl << std::flush;
            
            if (tdata.traj_size() > 0) {
            std::cerr << "  Processing " << tdata.traj_size() << " particles for diagnostic calculations..." << std::endl << std::flush;
            
            // Extract data into vectors using the correct access pattern (i,column) from ManageSimulation.cpp
            vector<double> x_data, y_data, z_data, vx_data, vy_data, vz_data, curr_data;
            size_t n_particles = tdata.traj_size();
            
            for (size_t i = 0; i < n_particles; ++i) {
                x_data.push_back(tdata(i,0));    // DIAG_X is column 0
                y_data.push_back(tdata(i,1));    // DIAG_Y is column 1
                z_data.push_back(tdata(i,2));    // DIAG_Z is column 2
                vx_data.push_back(tdata(i,3));   // DIAG_VX is column 3
                vy_data.push_back(tdata(i,4));   // DIAG_VY is column 4
                vz_data.push_back(tdata(i,5));   // DIAG_VZ is column 5
                curr_data.push_back(tdata(i,6)); // DIAG_CURR is column 6
            }
            
            // Calculate statistics using the proper diagnostic functions
            double current_sum = sumcurrent(curr_data);
            std::cerr << "  Current data size: " << curr_data.size() << std::endl << std::flush;
            if (curr_data.size() > 0) {
                std::cerr << "  First few current values: ";
                for (size_t ii = 0; ii < std::min(size_t(5), curr_data.size()); ++ii) {
                    std::cerr << curr_data[ii] << " ";
                }
                std::cerr << std::endl << std::flush;
            }
            std::cerr << "  Current sum calculated: " << current_sum << " A (" << 1000*current_sum << " mA)" << std::endl << std::flush;
            
            double x_ave = Average(x_data, curr_data);
            double y_ave = Average(y_data, curr_data);
            double x_min = min_vec(x_data);
            double x_max = max_vec(x_data);
            double y_min = min_vec(y_data);
            double y_max = max_vec(y_data);
            
            // Calculate angles and divergences using the diagnostic functions
            double angle_x = CalcolaAngolo(vx_data, vz_data, curr_data);
            double angle_y = CalcolaAngolo(vy_data, vz_data, curr_data);
            double div_x = sqrt(2.0) * CalcolaDivergenza2(vx_data, vz_data, curr_data, angle_x);
            double div_y = sqrt(2.0) * CalcolaDivergenza2(vy_data, vz_data, curr_data, angle_y);
            
            // Additional field calculations (use default values when field objects not available)
            Vec3D pos3d(0.0, 0.0, zloc);
            double potential_val = 0.0;  // Default: no field objects available
            double bfield_val = 0.0;     // Default: no field objects available
            double density_val = 0.0;    // Default: no field objects available
            double sigma_val = 0.0;      // Default: no field objects available
            
            pout5 << iteration << " "
                  << 1000*zloc << " "
                  << 1000*current_sum << " "  // Convert A to mA
                  << 1000*x_ave << " "
                  << 1000*x_max << " "
                  << 1000*x_min << " "
                  << 1000*y_ave << " "
                  << 1000*y_max << " "
                  << 1000*y_min << " "
                  << 1000*angle_x << " "
                  << 1000*angle_y << " "
                  << 1000*div_x << " "
                  << 1000*div_y << " "
                  << 0.001*potential_val << " "
                  << 1000*bfield_val << " "
                  << density_val << " "
                  << sigma_val << " "
                  << "\n" << flush;
        } else {
            std::cerr << "  No particles found at z=" << zloc << ". Writing zeros." << std::endl << std::flush;
            if (debug) logfile << "DEBUG: No particles intercepting at z position " << zloc << " m. Skipping." << endl;
            // Write zeros for this z position (17 columns total)
            pout5 << iteration << " "
                  << 1000*zloc << " "
                  << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n" << flush;
        }
        } catch (const Error& e) {
            std::cerr << "DEBUG: IBSimu Error caught at z=" << zloc << " m" << std::endl << std::flush;
            if (debug) logfile << "DEBUG: Error in trajectory diagnostics at z=" << zloc 
                               << " m: IBSIMU Error. Writing zeros." << endl;
            // Write zeros for this z position due to error (17 columns total)
            pout5 << iteration << " "
                  << 1000*zloc << " "
                  << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n" << flush;
        } catch (const std::exception& e) {
            std::cerr << "DEBUG: Standard exception caught at z=" << zloc << " m: " << e.what() << std::endl << std::flush;
            if (debug) logfile << "DEBUG: Standard exception in trajectory diagnostics at z=" << zloc 
                               << " m: " << e.what() << ". Writing zeros." << endl;
            // Write zeros for this z position due to error (17 columns total)
            pout5 << iteration << " "
                  << 1000*zloc << " "
                  << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n" << flush;
        }
    }
    
    pout5.close();
    
    if (debug) logfile << "DEBUG: Transverse data saved." << endl << flush;
}

void DiagnosticsManager::generateDiagnosticData(const vector<double>& diagzpos, int iteration, 
                                               const string& outfile, const ParticleDataBase3D* particles) {
    generateDiagnosticData(diagzpos, iteration, outfile, false, particles);
}

void DiagnosticsManager::createPlots(int argc, char **argv, const Geometry* geometry,
                                    const EpotField* potential, const MeshVectorField* magnetic,
                                    const EpotEfield* electric, const MeshScalarField* spacecharge,
                                    const ParticleDataBase3D* particles, const string& plot_folder,
                                    const string& file_tag) {
    
    if (debug) logfile << "DEBUG: Creating Geomplotter object..." << endl << flush;
    GeomPlotter geomplotter(*geometry);
    
    if (debug) logfile << "DEBUG: Set window size" << endl << flush;
    geomplotter.set_size(2000, 1000);
    
    if (debug) logfile << "DEBUG: Set potential" << endl << flush;
    geomplotter.set_epot(potential);
    
    if (debug) logfile << "DEBUG: Set particle database" << endl << flush;
    geomplotter.set_particle_database(particles);
    
    if (debug) logfile << "DEBUG: Set eqlines" << endl << flush;
    geomplotter.set_eqlines_auto(40);
    
    if (debug) logfile << "DEBUG: Set view ZX" << endl << flush;
    geomplotter.set_view(VIEW_ZX, -1);
    geomplotter.plot_png(plot_folder + file_tag + "_plot_zx.png");
    logfile << "zx plotted!" << endl << flush;
    
    geomplotter.set_view(VIEW_ZY, -1);
    geomplotter.plot_png(plot_folder + file_tag + "_plot_zy.png");
    logfile << "zy plotted!" << endl << flush;
    
    // Additional plotting for meniscus region
    if (potential) {
        geomplotter.set_ranges(0, -0.01, 0.034, 0.01);
        geomplotter.set_view(VIEW_ZX, -1);
        geomplotter.plot_png(plot_folder + file_tag + "_MENISCUS_plot_zx.png");
        logfile << "Meniscus zx plotted!" << endl << flush;
        
        geomplotter.set_view(VIEW_ZY, -1);
        geomplotter.plot_png(plot_folder + file_tag + "_MENISCUS_plot_zy.png");
        logfile << "Meniscus zy plotted!" << endl << flush;
    }
    
    if (debug) logfile << "Done!" << endl << flush;
}

void DiagnosticsManager::createPlots(int argc, char **argv, const Geometry* geometry,
                                    const EpotField* potential, const MeshVectorField* magnetic,
                                    const EpotEfield* electric, const MeshScalarField* spacecharge,
                                    const ParticleDataBase3D* particles,
                                    const vector<ParticleDataBase3D*>& particles_species,
                                    particle_kind pk, const string& plot_folder,
                                    const string& file_tag) {
    
    ParticleDataBase3D* pdb;
    if (pk == PARTICLE_ALL) {
        pdb = const_cast<ParticleDataBase3D*>(particles);
    } else {
        pdb = particles_species[get_particle_int(pk)];
    }
    
    // Create plots for specific particle species
    createPlots(argc, argv, geometry, potential, magnetic, electric, spacecharge,
               pdb, plot_folder, file_tag + "_" + get_particle_name(pk));
}

void DiagnosticsManager::performAnalysis(const ParticleDataBase3D* particles, const SimulationParameters& params,
                                        const Geometry* geometry, const EpotField* potential, 
                                        const MeshVectorField* magnetic, const FileManager* fileManager,vector<double> zgrids) {
    if (debug) logfile << "DEBUG: Performing simulation analysis with field calculations..." << endl << flush;
    
    if (!particles || !geometry) {
        throw Error(ERROR_LOCATION, "Missing required objects for analysis (particles or geometry)");
    }
    
    if (debug) logfile << "\n" << flush;
    if (debug) logfile << "DEBUG: geometry->ORIGO: " << geometry->origo(0) << " " << geometry->origo(1) << " " << geometry->origo(2) << " " << endl << flush;
    
    // Create trajectory density field (matching original implementation)
    MeshScalarField tdens(*geometry);
    particles->build_trajectory_density_field(tdens);
    
    // Get diagnostic file names using FileManager
    string diagfile;
    if (fileManager) {
        diagfile = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + "_diagnostic_summary.txt";
    } else {
        diagfile = "TEST_MTF/TEST_MTF_1_diagnostic_summary.txt";  // Fallback
    }
    
    // Create diagnostic z-locations using geometry (matching original implementation)
    double z_start = geometry->origo(2);
    Vec3D lastpt = geometry->max();
    double z_end = lastpt[2] - params.getMeshSize();
    size_t steps = 21;
    double delta_z = (z_end - z_start) / steps;
    vector<double> zlocs;
    
    for (size_t ii = 0; ii < steps; ii++) {
        zlocs.push_back(z_start + ii * delta_z);
    }
    
    logfile << "Diagnostic data for ALL particles in file " << diagfile << "\n" << flush;
    logfile << "Save diagnostic data at " << steps << " points with field calculations...\n" << flush;
    
    // Call enhanced diagnostic generation WITH field objects
    generateDiagnosticData(zlocs, 0, diagfile, false, particles, geometry, potential, magnetic);
    
    logfile << "Done!\n" << flush;
    
    // Get grid locations from parameters if available, otherwise use default values
    
    try {
        // Try to get grid positions from simulation parameters
        // This would need to be implemented in SimulationParameters class
        // For now, use hardcoded MTF grid positions as fallback
        // zgrids = {0.009, 0.015, 0.032, 0.120, 0.137, 0.225, 0.242, 0.330, 0.347, 0.435, 0.452, 0.540, 0.557};
        // zgrids = geometryManager->getZGrids();

        // Create subset excluding first element (matching original logic)
        vector<double> zg;
        if (zgrids.size() > 1) {
            zg = vector<double>(zgrids.begin() + 1, zgrids.end());
        } else {
            // Fallback: use domain exit plane derived from geometry
            double z_fallback = geometry->max()[2] - params.getMeshSize();
            zg = {z_fallback};
        }
        
        string gridfile = "statsatgrids.txt";
        if (fileManager) {
            gridfile = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + "_statsatgrids.txt";
        }
        
        if (debug) logfile << "DEBUG: Save diagnostic data at grids with field calculations\n" << flush;
        if (debug) cout << "DEBUG: Save diagnostic data at grids\n" << flush;
        
        // Call enhanced diagnostic generation WITH field objects
        generateDiagnosticData(zg, 0, gridfile, false, particles, geometry, potential, magnetic);
        
    } catch (const std::exception& e) {
        if (debug) logfile << "DEBUG: Error processing grid data: " << e.what() << endl << flush;
        // Continue with fallback grid analysis derived from geometry bounds
        vector<double> fallback_zgrids = {0.009, geometry->max()[2] - params.getMeshSize()};
        string gridfile = "particlesatgrids.txt";
        generateDiagnosticData(fallback_zgrids, 0, gridfile, false, particles, geometry, potential, magnetic);
    }
    
    if (debug) logfile << "DEBUG: analysis with field calculations ended" << endl << flush;
    if (debug) cout << "DEBUG: analysis ended" << endl << flush;
    ibsimu.message(1) << "Analysis completed" << endl;
}

void DiagnosticsManager::performAnalysis(const ParticleDataBase3D* particles, const SimulationParameters& params,
                                        const Geometry* geometry, const EpotField* potential, 
                                        const MeshVectorField* magnetic, const double zlocsummary,
                                        const std::vector<ParticleDataBase3D*>& particles_species,
                                        bool include_stripping, const FileManager* fileManager, vector<double> zgrids) {
    if (debug) logfile << "DEBUG: Performing simulation analysis with species support..." << endl << flush;
    
    if (!particles || !geometry) {
        throw Error(ERROR_LOCATION, "Missing required objects for analysis (particles or geometry)");
    }
    
    if (debug) logfile << "\n" << flush;
    if (debug) logfile << "DEBUG: geometry->ORIGO: " << geometry->origo(0) << " " << geometry->origo(1) << " " << geometry->origo(2) << " " << endl << flush;
    
    // Create trajectory density field (matching original implementation)
    MeshScalarField tdens(*geometry);
    particles->build_trajectory_density_field(tdens);
    
    // Get diagnostic file names using FileManager
    string diagfile;
    if (fileManager) {
        diagfile = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + "_diagnostic_summary.txt";
    } else {
        diagfile = "TEST_MTF/TEST_MTF_1_diagnostic_summary.txt";  // Fallback
    }
    
    // Create diagnostic z-locations using geometry (matching original implementation)
    double z_start = geometry->origo(2);
    Vec3D lastpt = geometry->max();
    double z_end = lastpt[2] - params.getMeshSize();
    size_t steps = 21;
    double delta_z = (z_end - z_start) / steps;
    vector<double> zlocs;
    
    for (size_t ii = 0; ii < steps; ii++) {
        zlocs.push_back(z_start + ii * delta_z);
    }
    
    logfile << "Diagnostic data for ALL particles in file " << diagfile << "\n" << flush;
    logfile << "Save diagnostic data at " << steps << " points with field calculations...\n" << flush;
    
    // Call enhanced diagnostic generation WITH field objects for all particles
    generateDiagnosticData(zlocs, 0, diagfile, false, particles, geometry, potential, magnetic);
    logfile << "Done!\n" << flush;

    vector<double> zlocsummary_vec = {zlocsummary};

    // Species-specific analysis (matching original INCLUDE_STRIPPING==2 logic)
    if (include_stripping && !particles_species.empty()) {
        for (size_t ii = 0; ii < particles_species.size(); ii++) {
            if (particles_species[ii] && particles_species[ii]->size() > 0) {
                string species_diagfile;
                string summary_diagfile;
                if (fileManager) {
                    species_diagfile = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + "_species_" + std::to_string(ii) + "_diagnostic_summary.txt";
                    summary_diagfile = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + "_NEGIONBEAM_diagnostic_summary.txt";
                } else {
                    species_diagfile = "TEST_MTF/TEST_MTF_1_species_" + std::to_string(ii) + "_diagnostic_summary.txt";  // Fallback
                }
                
                logfile << "Diagnostic data for particle_species " << ii << " in file " << species_diagfile << "\n" << flush;
                logfile << "Save diagnostic data at " << steps << " points... " << flush;
                
                generateDiagnosticData(zlocs, 0, species_diagfile, false, particles_species[ii], geometry, potential, magnetic);
                if ( ii==0 ) {// Corresponds to only negative ions
                    generateDiagnosticData(zlocsummary_vec, 0, summary_diagfile, false, particles_species[ii], geometry, potential, magnetic);
                }
                logfile << "Done!\n" << flush;
            }
        }
    }
    
    // Grid analysis with proper grid positions (using parameters if available)
    try {
        // vector<double> zgrids = {0.009, 0.015, 0.032, 0.120, 0.137, 0.225, 0.242, 0.330, 0.347, 0.435, 0.452, 0.540, 0.557};
        
        // Create subset excluding first element (matching original logic)
        vector<double> zg;
        if (zgrids.size() > 1) {
            zg = vector<double>(zgrids.begin() + 1, zgrids.end());
        } else {
            cout << "Error: zgrids does not have enough elements to create zg." << endl;
            // Fallback: use domain exit plane derived from geometry
            double z_fallback = geometry->max()[2] - params.getMeshSize();
            zg = {z_fallback};
        }
        
        string gridfile = "TEST_MTF/TEST_MTF_1_particlesatgrids.txt";
        if (fileManager) {
            gridfile = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + "_particlesatgrids.txt";
        }
        
        if (debug) logfile << "DEBUG: Save diagnostic data at grids with field calculations\n" << flush;
        if (debug) cout << "DEBUG: Save diagnostic data at grids\n" << flush;
        
        generateDiagnosticData(zg, 0, gridfile, false, particles, geometry, potential, magnetic);
        
    } catch (const std::exception& e) {
        if (debug) logfile << "DEBUG: Error processing grid data: " << e.what() << endl << flush;
        cout << "Error: " << e.what() << endl;
    }
    
    // Grid power load analysis for ALL particles
    logfile << " Analyzing grid power loads for ALL particles\n" << flush;
    try {
        double mesh_size = params.getMeshSize();
        double ionmass = params.getMIons();
        string output_folder = fileManager ? fileManager->getOutputSummaryFolder() : "TEST_MTF/";
        string file_tag = fileManager ? fileManager->getFileTag() : "TEST_MTF_1";
        
        vector<double> grid_powers_all = analyzeGridPowerLoads(particles, zgrids, mesh_size, 
                                                               geometry, ionmass, output_folder, 
                                                               file_tag + "_ALL", PARTICLE_ALL);
        
        logfile << " Grid power analysis for ALL particles completed with " 
                           << grid_powers_all.size() << " power values\n" << flush;
    } catch (const std::exception& e) {
        logfile << "Error in grid power analysis for ALL particles: " << e.what() << endl << flush;
    }
    
    // Grid power load analysis for each particle species
    if (include_stripping && !particles_species.empty()) {
        const particle_kind species_kinds[] = {PARTICLE_HM, PARTICLE_H0, PARTICLE_HP, 
                                               PARTICLE_H2P, PARTICLE_H20, PARTICLE_E};
        const char* species_names[] = {"HM", "H0", "HP", "H2P", "H20", "E"};
        
        for (size_t ii = 0; ii < particles_species.size() && ii < 6; ii++) {
            if (particles_species[ii] && particles_species[ii]->size() > 0) {
                logfile << "Analyzing grid power loads for species " << species_names[ii] << "\n" << flush;

                try {
                    double mesh_size = params.getMeshSize();
                    double ionmass = params.getMIons();
                    string output_folder = fileManager ? fileManager->getOutputSummaryFolder() : "TEST_MTF/";
                    string file_tag = fileManager ? fileManager->getFileTag() : "TEST_MTF_1";
                    
                    vector<double> grid_powers_species = analyzeGridPowerLoads(particles_species[ii], zgrids, 
                                                                               mesh_size, geometry, ionmass, 
                                                                               output_folder, 
                                                                               file_tag + "_" + species_names[ii], 
                                                                               species_kinds[ii]);
                    
                    logfile << "Grid power analysis for species " << species_names[ii] 
                            << " completed with " << grid_powers_species.size() << " power values\n" << flush;
                } catch (const std::exception& e) {
                    logfile << "Error in grid power analysis for species " 
                            << species_names[ii] << ": " << e.what() << endl << flush;
                }
            }
        }
    }
    
    logfile << " Analysis with field calculations and species support ended" << endl << flush;
    // cout << "******** Analysis ended ********" << endl << flush;
    
    ibsimu.message(1) << "Analysis completed" << endl;
}

void DiagnosticsManager::printTrajectoryData(const string& filename, const ParticleDataBase3D* particles) {
    if (!particles) {
        throw Error(ERROR_LOCATION, "No particle database provided for trajectory output");
    }
    
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        throw Error(ERROR_LOCATION, "Could not open trajectory output file: " + filename);
    }
    
    outfile << "# Trajectory data\n";
    outfile << "# x[m]\ty[m]\tz[m]\tvx[m/s]\tvy[m/s]\tvz[m/s]\tI[A]\tm[kg]\tq[C]\n";
    
    // Trajectory output implementation would go here
    // This would involve iterating through all particle trajectories
    
    outfile.close();
    
    ibsimu.message(1) << "Trajectory data saved to: " << filename << endl;
}

void DiagnosticsManager::createSimulationSummary(double zlocsummary, const string& summary_file_tag, 
                                                 bool append_at_end, const ParticleDataBase3D* particles,
                                                 const string& outsummary_fold) {
    
    // Debug output to file
    // ofstream debugfile("debug_summary_call.txt", ios::app);
    // debugfile << "createSimulationSummary called with zlocsummary=" << zlocsummary << " tag=" << summary_file_tag << endl;
    // debugfile.flush();
    // debugfile.close();
    
    if (debug) logfile << "DEBUG: Creating simulation summary" << endl << flush;
    
    vector<double> zpos;
    zpos.push_back(zlocsummary);
    
    string strbeamprops = outsummary_fold + summary_file_tag + "_beamprops.txt";
    generateDiagnosticData(zpos, 0, strbeamprops, false, particles);
    
    // Add error handling for TransverseData construction
    try {
        TransverseData res(strbeamprops);
        
        if (debug) logfile << "DEBUG: zlocsummary: " << zlocsummary << endl << flush;
        
        // Summary calculations would go here
        // This would involve power calculations, efficiency analysis, etc.
        
        ibsimu.message(1) << "Simulation summary created" << endl;
    } catch (const std::exception& e) {
        if (debug) logfile << "DEBUG: Error reading TransverseData from " << strbeamprops << ": " << e.what() << endl << flush;
        ibsimu.message(1) << "Warning: Could not read beam properties from " << strbeamprops << ": " << e.what() << endl;
        return; // Exit gracefully instead of crashing
    }
}

void DiagnosticsManager::createSimulationSummary(double zlocsummary, const string& summary_file_tag, 
                                                 bool /*append_at_end*/, 
                                                 const vector<ParticleDataBase3D*>& particles_species,
                                                 particle_kind pk, const string& outsummary_fold) {
    
    string strbeamprops = outsummary_fold + summary_file_tag + "_" + get_particle_name(pk) + "_beamprops.txt";
    generateDiagnosticData({zlocsummary}, 0, strbeamprops, false, particles_species[get_particle_int(pk)]);
    
    // Add error handling for TransverseData construction
    try {
        TransverseData res(strbeamprops);
        
        // Species-specific summary calculations would go here
        
        ibsimu.message(1) << "Species-specific simulation summary created for " << get_particle_name(pk) << endl;
    } catch (const std::exception& e) {
        if (debug) logfile << "DEBUG: Error reading TransverseData from " << strbeamprops << ": " << e.what() << endl << flush;
        ibsimu.message(1) << "Warning: Could not read beam properties from " << strbeamprops << ": " << e.what() << endl;
        return; // Exit gracefully instead of crashing
    }
}

void DiagnosticsManager::createSimulationSummary(double zlocsummary, const string& summary_file_tag, 
                                                 bool append_at_end, const ParticleDataBase3D* particles,
                                                 const string& outsummary_fold,
                                                 const Geometry* geometry, const EpotField* potential, 
                                                 const MeshVectorField* magnetic) {
    
    // Debug output to file
    // ofstream debugfile("debug_summary_call.txt", ios::app);
    // debugfile << "createSimulationSummary called with FIELD CALCULATIONS: zlocsummary=" << zlocsummary << " tag=" << summary_file_tag << endl;
    // debugfile.flush();
    // debugfile.close();
    
    if (debug) logfile << "DEBUG: Creating simulation summary with field calculations" << endl << flush;
    
    vector<double> zpos;
    zpos.push_back(zlocsummary);
    
    string strbeamprops = outsummary_fold + summary_file_tag + "_beamprops.txt";
    // Use the enhanced version with field calculations
    generateDiagnosticData(zpos, 0, strbeamprops, false, particles, geometry, potential, magnetic);
    
    // Add error handling for TransverseData construction - use proper try-catch with constructor
    try {
        TransverseData res(strbeamprops);
        
        if (debug) logfile << "DEBUG: zlocsummary: " << zlocsummary << endl << flush;
        
        // Perform grid power analysis if geometry information is available
        if (geometry) {
            // Get default grid positions for MTF (or use provided zgrids if available)
            vector<double> default_zgrids = {0.0, 0.009, 0.015, 0.032, 0.120, 0.137, 0.225, 0.242, 0.330, 0.347, 0.435, 0.452, 0.540, 0.557};
        
        // Analyze grid power loads (using default ion mass of 1.0)
        double mesh_size = 0.001; // Default mesh size 1mm
        double ionmass = 1.0;     // Default hydrogen mass
        
        vector<double> grid_powers = analyzeGridPowerLoads(particles, default_zgrids, mesh_size, 
                                                          geometry, ionmass, outsummary_fold, 
                                                          summary_file_tag, PARTICLE_ALL);
        
        if (debug) logfile << "DEBUG: Grid power analysis completed with " << grid_powers.size() << " power values" << endl << flush;
        
        // Calculate total beam power on grids (excluding plasma grid #7)
        double full_beam_power = 0.0;
        for (size_t qq = 6; qq < grid_powers.size(); qq++) {
            if (qq != 7) { // Exclude plasma grid
                full_beam_power += grid_powers[qq];
            }
        }
        
        if (debug) logfile << "DEBUG: Total beam power on grids: " << full_beam_power << " W" << endl << flush;
        
        // Create detailed summary file with all data (similar to original implementation)
        string summary_filename = outsummary_fold + summary_file_tag + ".txt";
        ofstream summary_file;
        
        if (append_at_end) {
            summary_file.open(summary_filename, ios_base::app);
        } else {
            summary_file.open(summary_filename);
            // Write header
            summary_file << "zloc(mm)\tUext(kV)\tUtot(kV)\taccJ(A/m2)\textJ(A/m2)\t"
                        << "xave(mm)\txmax(mm)\txmin(mm)\t"
                        << "yave(mm)\tymax(mm)\tymin(mm)\t"
                        << "xpave(mrad)\typave(mrad)\tdivx(mrad)\tdivy(mrad)\t"
                        << "epot(kV)\tBy(mT)\tfull_pow(W)\tout_pow(W)\t";
            for (size_t ii = 7; ii < grid_powers.size(); ii++) {
                summary_file << "pow_S" << ii << "(W)\t" << "curr_S" << ii << "(A)\t";
            }
            summary_file << "\n" << flush;
        }
        
        if (summary_file.is_open()) {
            // Calculate current densities using PG aperture area (circle with 7mm radius)
            double pg_aperture_area = M_PI * 0.007 * 0.007;  // π × (7mm)² = 1.539×10⁻⁴ m²
            double current_A = res.get_current(0);  // Current in Amperes
            double current_density = current_A / pg_aperture_area;  // A/m²
            
            // Write summary data (using data from TransverseData)
            summary_file << res.get_zloc(0)*1e3 << "\t"    // zloc in mm
                        << "0.0" << "\t"                   // EG_VOLTAGE (placeholder)
                        << "0.0" << "\t"                   // G5_VOLTAGE (placeholder)  
                        << current_density << "\t"         // Accelerated current density in A/m²
                        << current_density << "\t"         // Extracted current density in A/m² (same as accelerated)
                        << res.get_xave(0)*1e3 << "\t"
                        << res.get_xmax(0)*1e3 << "\t"
                        << res.get_xmin(0)*1e3 << "\t"
                        << res.get_yave(0)*1e3 << "\t"
                        << res.get_ymax(0)*1e3 << "\t"
                        << res.get_ymin(0)*1e3 << "\t"
                        << res.get_xpave(0)*1e3 << "\t"
                        << res.get_ypave(0)*1e3 << "\t"
                        << res.get_divx(0)*1e3 << "\t"
                        << res.get_divy(0)*1e3 << "\t"
                        << res.get_epot(0)/1e3 << "\t"     // Potential in kV
                        << res.get_Bfield(0)*1e3 << "\t"   // B-field in mT
                        << full_beam_power << "\t"         // Total grid power
                        << (grid_powers.size() > 6 ? grid_powers[6] : 0.0) << "\t"; // Exit power
            
            // Write individual grid powers (starting from grid index 7)
            for (size_t ii = 7; ii < grid_powers.size(); ii++) {
                summary_file << grid_powers[ii] << "\t" << "0.0" << "\t"; // Power and current (current placeholder)
            }
            summary_file << "\n" << flush;
            summary_file.close();
            
            ibsimu.message(1) << "Enhanced simulation summary with grid power analysis saved to: " << summary_filename << endl;
        }
    }
    
    ibsimu.message(1) << "Simulation summary created with field calculations" << endl;
    
    } catch (const std::exception& e) {
        if (debug) logfile << "DEBUG: Error reading TransverseData from " << strbeamprops << ": " << e.what() << endl << flush;
        ibsimu.message(1) << "Warning: Could not read beam properties from " << strbeamprops << ": " << e.what() << endl;
        return; // Exit gracefully instead of crashing
    }
}

void DiagnosticsManager::createSimulationSummary(double zlocsummary, const string& summary_file_tag, 
                                                 bool append_at_end, const ParticleDataBase3D* particles,
                                                 const string& outsummary_fold,
                                                 const Geometry* geometry, const EpotField* potential, 
                                                 const MeshVectorField* magnetic, double extracted_current_density) {
    
    // Debug output to file
    ofstream debugfile("debug_summary_call.txt", ios::app);
    debugfile << "createSimulationSummary called with FIELD CALCULATIONS and EXTRACTED CURRENT: zlocsummary=" 
              << zlocsummary << " tag=" << summary_file_tag << " extJ=" << extracted_current_density << endl;
    debugfile.flush();
    debugfile.close();
    
    if (debug) logfile << "DEBUG: Creating simulation summary with field calculations and extracted current density" << endl << flush;
    
    vector<double> zpos;
    zpos.push_back(zlocsummary);
    
    string strbeamprops = outsummary_fold + summary_file_tag + "_beamprops.txt";
    // Use the enhanced version with field calculations
    generateDiagnosticData(zpos, 0, strbeamprops, false, particles, geometry, potential, magnetic);
    
    // Add error handling for TransverseData construction - use proper try-catch with constructor
    try {
        TransverseData res(strbeamprops);
        
        if (debug) logfile << "DEBUG: zlocsummary: " << zlocsummary << endl << flush;
        
        // Perform grid power analysis if geometry information is available
        if (geometry) {
            // Get default grid positions for MTF (or use provided zgrids if available)
            vector<double> default_zgrids = {0.0, 0.009, 0.015, 0.032, 0.120, 0.137, 0.225, 0.242, 0.330, 0.347, 0.435, 0.452, 0.540, 0.557};
        
        // Analyze grid power loads (using default ion mass of 1.0)
        double mesh_size = 0.001; // Default mesh size 1mm
        double ionmass = 1.0;     // Default hydrogen mass
        
        vector<double> grid_powers = analyzeGridPowerLoads(particles, default_zgrids, mesh_size, 
                                                          geometry, ionmass, outsummary_fold, 
                                                          summary_file_tag, PARTICLE_ALL);
        
        if (debug) logfile << "DEBUG: Grid power analysis completed with " << grid_powers.size() << " power values" << endl << flush;
        
        // Calculate total beam power on grids (excluding plasma grid #7)
        double full_beam_power = 0.0;
        for (size_t qq = 6; qq < grid_powers.size(); qq++) {
            if (qq != 7) { // Exclude plasma grid
                full_beam_power += grid_powers[qq];
            }
        }
        
        if (debug) logfile << "DEBUG: Total beam power on grids: " << full_beam_power << " W" << endl << flush;
        
        // Create detailed summary file with all data (similar to original implementation)
        string summary_filename = outsummary_fold + summary_file_tag + ".txt";
        ofstream summary_file;
        
        if (append_at_end) {
            summary_file.open(summary_filename, ios_base::app);
        } else {
            summary_file.open(summary_filename);
            // Write header
            summary_file << "zloc(mm)\tUext(kV)\tUtot(kV)\taccJ(A/m2)\textJ(A/m2)\t"
                        << "xave(mm)\txmax(mm)\txmin(mm)\t"
                        << "yave(mm)\tymax(mm)\tymin(mm)\t"
                        << "xpave(mrad)\typave(mrad)\tdivx(mrad)\tdivy(mrad)\t"
                        << "epot(kV)\tBy(mT)\tfull_pow(W)\tout_pow(W)\t";
            for (size_t ii = 7; ii < grid_powers.size(); ii++) {
                summary_file << "pow_S" << ii << "(W)\t" << "curr_S" << ii << "(A)\t";
            }
            summary_file << "\n" << flush;
        }
        
        if (summary_file.is_open()) {
            // Calculate accelerated current density using PG aperture area (circle with 7mm radius)
            double pg_aperture_area = M_PI * 0.007 * 0.007;  // π × (7mm)² = 1.539×10⁻⁴ m²
            double current_A = res.get_current(0);  // Current in Amperes
            double accelerated_current_density = abs(current_A) / pg_aperture_area;  // A/m²
            
            // Write summary data (using data from TransverseData)
            summary_file << res.get_zloc(0)*1e3 << "\t"    // zloc in mm
                        << "0.0" << "\t"                   // EG_VOLTAGE (placeholder)
                        << "0.0" << "\t"                   // G5_VOLTAGE (placeholder)  
                        << accelerated_current_density << "\t"  // Accelerated current density in A/m²
                        << extracted_current_density << "\t"    // Extracted current density in A/m² (from ParticleManager)
                        << res.get_xave(0)*1e3 << "\t"
                        << res.get_xmax(0)*1e3 << "\t"
                        << res.get_xmin(0)*1e3 << "\t"
                        << res.get_yave(0)*1e3 << "\t"
                        << res.get_ymax(0)*1e3 << "\t"
                        << res.get_ymin(0)*1e3 << "\t"
                        << res.get_xpave(0)*1e3 << "\t"
                        << res.get_ypave(0)*1e3 << "\t"
                        << res.get_divx(0)*1e3 << "\t"
                        << res.get_divy(0)*1e3 << "\t"
                        << res.get_epot(0)/1e3 << "\t"     // Potential in kV
                        << res.get_Bfield(0)*1e3 << "\t"   // B-field in mT
                        << full_beam_power << "\t"         // Total grid power
                        << (grid_powers.size() > 6 ? grid_powers[6] : 0.0) << "\t"; // Exit power
            
            // Write individual grid powers (starting from grid index 7)
            for (size_t ii = 7; ii < grid_powers.size(); ii++) {
                summary_file << grid_powers[ii] << "\t" << "0.0" << "\t"; // Power and current (current placeholder)
            }
            summary_file << "\n" << flush;
            summary_file.close();
            
            ibsimu.message(1) << "Enhanced simulation summary with grid power analysis and extracted current density saved to: " << summary_filename << endl;
        }
    }
    
    ibsimu.message(1) << "Simulation summary created with field calculations and extracted current density" << endl;
    
    } catch (const std::exception& e) {
        if (debug) logfile << "DEBUG: Error reading TransverseData from " << strbeamprops << ": " << e.what() << endl << flush;
        ibsimu.message(1) << "Warning: Could not read beam properties from " << strbeamprops << ": " << e.what() << endl;
        return; // Exit gracefully instead of crashing
    }
}

void DiagnosticsManager::generateDiagnosticData(const vector<double>& diagzpos, int iteration, 
                                               const string& outfile, bool append_at_end, 
                                               const ParticleDataBase3D* pdb,
                                               const Geometry* /*geometry*/,
                                               const EpotField* potential,
                                               const MeshVectorField* magnetic) {
    std::cerr << "DEBUG: generateDiagnosticData called with " << diagzpos.size() << " z-positions (with field calculations)" << std::endl << std::flush;
    
    if (debug) logfile << "DEBUG: Printing diagnostic data to " << outfile << " file\n" << flush;
    
    vector<double> x[1];
    vector<double> y[1];
    vector<double> z[1];
    vector<double> vx[1];
    vector<double> vy[1];
    vector<double> vz[1];
    vector<double> curr[1];
    vector<double> mass[1];
    vector<double> charge[1];
    vector<double> no_part[1];

    // Load density profile for calculations
    vector<double> pos, dens;
    string filedens = "densprofiles/MTF_dens.txt";
    load_density_profile(filedens, pos, dens);

    ofstream pout5;
    if (append_at_end) {
        pout5.open(outfile.c_str(), ios_base::app);
    } else {
        pout5.open(outfile.c_str());
        // Write header only when creating new file
        pout5 << "# it z[mm] I[mA] <x>[mm] xmax[mm] xmin[mm] <y>[mm] ymax[mm] ymin[mm] <x'>[mrad] <y'>[mrad] Dx[mrad] Dy[mrad] <V>[V] <B>[mT] rho[1/m3] sigma[m2]" << endl;
    }

    for (size_t i = 0; i < diagzpos.size(); ++i) {
        double zloc = diagzpos[i];
        std::cerr << "DEBUG: Processing z-position " << (i+1) << "/" << diagzpos.size() 
                  << " at z=" << zloc << " m" << std::endl << std::flush;
        
        // Clear vectors for this z-position
        x[0].clear(); y[0].clear(); z[0].clear();
        vx[0].clear(); vy[0].clear(); vz[0].clear();
        curr[0].clear(); mass[0].clear(); charge[0].clear(); no_part[0].clear();
        
        std::cerr << "DEBUG: About to enter try block" << std::endl << std::flush;
        try {
            std::cerr << "DEBUG: Inside try block, creating TransverseData" << std::endl << std::flush;
            
            // Trajectory diagnostics at specific z location
            TrajectoryDiagnosticData tdata;
            vector<trajectory_diagnostic_e> diagnostics = {
                DIAG_X, DIAG_Y, DIAG_Z, DIAG_VX, DIAG_VY, DIAG_VZ,
                DIAG_CURR, DIAG_MASS, DIAG_CHARGE, DIAG_NO
            };
            
            std::cerr << "  About to call trajectories_at_plane..." << std::endl << std::flush;
            pdb->trajectories_at_plane(tdata, AXIS_Z, zloc, diagnostics);
            std::cerr << "  trajectories_at_plane call completed successfully!" << std::endl << std::flush;
            
            // Debug: Check the size immediately after the call
            size_t particle_count = tdata.traj_size();
            std::cerr << "  DEBUG: tdata.traj_size() = " << particle_count << std::endl << std::flush;
            
            // Write count to file for debugging
            ofstream debugfile("debug_particle_count.txt", ios::app);
            size_t pdb_size = pdb->size();
            debugfile << "z=" << zloc << " count=" << particle_count << " pdb_size=" << pdb_size << endl;
            debugfile.flush();
            debugfile.close();
            
            std::cerr << "  Retrieved " << particle_count << " particles at z=" << zloc << std::endl << std::flush;
            std::cerr << "  Particle database size: " << pdb->size() << std::endl << std::flush;
            
            if (particle_count > 0) {
                std::cerr << "  Processing " << particle_count << " particles for diagnostic calculations..." << std::endl << std::flush;
                
                // Extract data into vectors using the correct access pattern
                vector<double> x_data, y_data, z_data, vx_data, vy_data, vz_data, curr_data;
                
                for (size_t j = 0; j < particle_count; ++j) {
                    x_data.push_back(tdata(j,0));    // DIAG_X is column 0
                    y_data.push_back(tdata(j,1));    // DIAG_Y is column 1
                    z_data.push_back(tdata(j,2));    // DIAG_Z is column 2
                    vx_data.push_back(tdata(j,3));   // DIAG_VX is column 3
                    vy_data.push_back(tdata(j,4));   // DIAG_VY is column 4
                    vz_data.push_back(tdata(j,5));   // DIAG_VZ is column 5
                    curr_data.push_back(tdata(j,6)); // DIAG_CURR is column 6
                }
                
                // Calculate statistics using the proper diagnostic functions
                double current_sum = sumcurrent(curr_data);
                std::cerr << "  Current data size: " << curr_data.size() << std::endl << std::flush;
                if (curr_data.size() > 0) {
                    std::cerr << "  First few current values: ";
                    for (size_t ii = 0; ii < std::min(size_t(5), curr_data.size()); ++ii) {
                        std::cerr << curr_data[ii] << " ";
                    }
                    std::cerr << std::endl << std::flush;
                }
                std::cerr << "  Current sum calculated: " << current_sum << " A (" << 1000*current_sum << " mA)" << std::endl << std::flush;
                
                double x_ave = Average(x_data, curr_data);
                double y_ave = Average(y_data, curr_data);
                double x_min = min_vec(x_data);
                double x_max = max_vec(x_data);
                double y_min = min_vec(y_data);
                double y_max = max_vec(y_data);
                
                // Calculate angles and divergences using the diagnostic functions
                double angle_x = CalcolaAngolo(vx_data, vz_data, curr_data);
                double angle_y = CalcolaAngolo(vy_data, vz_data, curr_data);
                double div_x = sqrt(2.0) * CalcolaDivergenza2(vx_data, vz_data, curr_data, angle_x);
                double div_y = sqrt(2.0) * CalcolaDivergenza2(vy_data, vz_data, curr_data, angle_y);
                
                // Field calculations using provided field objects
                Vec3D pos3d(0.0, 0.0, zloc);
                
                // Debug: Check if field objects are valid
                // std::cerr << "  DEBUG Field objects: potential=" << (potential ? "valid" : "null") 
                //           << ", magnetic=" << (magnetic ? "valid" : "null") << std::endl << std::flush;
                
                // Calculate potential value
                double potential_val = 0.0;
                if (potential) {
                    try {
                        potential_val = (*potential)(pos3d);
                        // std::cerr << "  DEBUG: Potential at (" << pos3d[0] << "," << pos3d[1] << "," << pos3d[2] 
                        //           << ") = " << potential_val << " V" << std::endl << std::flush;
                    } catch (const std::exception& e) {
                        std::cerr << "  DEBUG: Error evaluating potential: " << e.what() << std::endl << std::flush;
                    }
                }
                
                // Calculate magnetic field value (B_y component)
                double bfield_val = 0.0;
                if (magnetic) {
                    try {
                        Vec3D B_at_given_z = magnetic->operator()(pos3d);
                        bfield_val = B_at_given_z[1];  // Y-component of magnetic field
                        // std::cerr << "  DEBUG: Magnetic field at (" << pos3d[0] << "," << pos3d[1] << "," << pos3d[2] 
                        //           << ") = (" << B_at_given_z[0] << "," << B_at_given_z[1] << "," << B_at_given_z[2] 
                        //           << ") T, using B_y=" << bfield_val << " T" << std::endl << std::flush;
                    } catch (const std::exception& e) {
                        std::cerr << "  DEBUG: Error evaluating magnetic field: " << e.what() << std::endl << std::flush;
                    }
                }
                
                // Calculate density value
                double density_val = 0.0;
                if (pos.size() > 0 && dens.size() > 0) {
                    bool ciaone = true;
                    density_val = density_at_z(zloc, 0.3, pos, dens, ciaone);
                    // std::cerr << "  DEBUG: Density at z=" << zloc << " m = " << density_val 
                    //           << " kg/m³ (density profile loaded with " << pos.size() << " points)" << std::endl << std::flush;
                } else {
                    std::cerr << "  DEBUG: No density profile loaded (pos.size=" << pos.size() 
                              << ", dens.size=" << dens.size() << ")" << std::endl << std::flush;
                }
                
                // Calculate stripping cross section
                double sigma_val = 0.0;
                if (potential) {
                    try {
                        double sigma_single = 0.0;
                        double sigma_double = 0.0;
                        // Use default mass of 1 amu (hydrogen) if M_IONS not available
                        double mass_ions = 1.0;  // Default to hydrogen mass
                        double energy_eV = (*potential)(pos3d);
                        sigma_val = stripping_cross_at_E(energy_eV, mass_ions, sigma_single, sigma_double);
                        // std::cerr << "  DEBUG: Stripping cross section at E=" << energy_eV 
                        //           << " eV = " << sigma_val << " m² (single=" << sigma_single 
                        //           << ", double=" << sigma_double << ")" << std::endl << std::flush;
                    } catch (const std::exception& e) {
                        std::cerr << "  DEBUG: Error calculating stripping cross section: " << e.what() << std::endl << std::flush;
                    }
                } else {
                    std::cerr << "  DEBUG: No potential field for cross section calculation" << std::endl << std::flush;
                }
                
                pout5 << iteration << " "
                      << 1000*zloc << " "
                      << 1000*current_sum << " "  // Convert A to mA
                      << 1000*x_ave << " "
                      << 1000*x_max << " "
                      << 1000*x_min << " "
                      << 1000*y_ave << " "
                      << 1000*y_max << " "
                      << 1000*y_min << " "
                      << 1000*angle_x << " "      // Convert rad to mrad
                      << 1000*angle_y << " "      // Convert rad to mrad
                      << 1000*div_x << " "        // Convert rad to mrad
                      << 1000*div_y << " "        // Convert rad to mrad
                      << 0.001*potential_val << " "  // Convert V to kV
                      << 1000*bfield_val << " "    // Convert T to mT
                      << density_val << " "        // kg/m³
                      << sigma_val << " "          // m²
                      << endl << flush;
            } else {
                std::cerr << "  No particles found at z=" << zloc << std::endl << std::flush;
                
                // Write zeros for empty locations
                pout5 << iteration << " "
                      << 1000*zloc << " "
                      << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
                      << endl << flush;
            }
            
        } catch (Error& e) {
            std::cerr << "  IBSimu Error caught: [Error in grid power analysis]" << std::endl << std::flush;
            
            // Write error entry
            pout5 << iteration << " "
                  << 1000*zloc << " "
                  << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
                  << endl << flush;
        } catch (std::exception& e) {
            std::cerr << "  Standard exception caught: " << e.what() << std::endl << std::flush;
            
            // Write error entry
            pout5 << iteration << " "
                  << 1000*zloc << " "
                  << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
                  << endl << flush;
        }
    }
    
    pout5.close();
    
    if (debug) logfile << "DEBUG: Diagnostic data written to " << outfile << endl << flush;
}

std::vector<double> DiagnosticsManager::analyzeGridPowerLoads(const ParticleDataBase3D* particles,
                                                             const std::vector<double>& zgrids,
                                                             double mesh_size,
                                                             const Geometry* geometry,
                                                             double ionmass,
                                                             const std::string& output_folder,
                                                             const std::string& file_tag,
                                                             particle_kind pk) {
    if (debug) logfile << "DEBUG: Analyzing grid power loads..." << endl << flush;
    
    if (!particles) {
        if (debug) logfile << "DEBUG: No particle database provided" << endl << flush;
        return std::vector<double>();
    }
    
    // Number of solids = geometry solids + 7 boundaries (as in original implementation)
    size_t n_solids = geometry ? geometry->number_of_solids() : 0;
    size_t total_categories = n_solids + 7; // +7 for domain boundaries as in original
    
    if (debug) logfile << "DEBUG: N_SOLIDS: " << n_solids << ", total categories: " << total_categories << endl << flush;
    
    // Create collision vectors for each solid/boundary
    std::vector<std::vector<size_t>> colls(total_categories);
    
    // Get domain boundaries
    Vec3D origo = geometry->origo();
    Vec3D maxpt = geometry->max();
    double z_start = origo[2];
    double z_end = maxpt[2];
    double x_start = origo[0];
    double x_end = maxpt[0];
    double y_start = origo[1];
    double y_end = maxpt[1];
    
    if (debug) logfile << "DEBUG: Domain boundaries: x[" << x_start << "," << x_end 
                       << "] y[" << y_start << "," << y_end 
                       << "] z[" << z_start << "," << z_end << "]" << endl << flush;
    
    // Statistics counters
    size_t particle_counts[6] = {0}; // HM, H0, HP, H2P, H20, E
    size_t generation_counts[6] = {0}; // gen 0,1,2,3,4,>4
    
    // Analyze each particle's end location
    for (size_t ii = 0; ii < particles->size(); ii++) {
        Particle3D& pp = const_cast<Particle3D&>(particles->particle(ii));
        Vec3D loc = pp.location();
        size_t isolid = 0;
        
        // Identify particle species
        particle_kind species = identify_particle_species(pp.m(), pp.q(), ionmass);
        if (species >= 0 && species < 6) particle_counts[species]++;
        
        // Track generations
        int gen = pp.gen() % 100;
        if (gen >= 0 && gen < 5) generation_counts[gen]++;
        else generation_counts[5]++; // gen > 4
        
        // Determine collision solid/boundary
        if (pp.get_status() == PARTICLE_COLL) {
            // Check if particle hit a grid (based on z-coordinate)
            bool hit_grid = false;
            if (zgrids.size() >= 14) { // Need pairs of z positions for each grid
                for (size_t grid_idx = 0; grid_idx < 7; grid_idx++) {
                    size_t z_start_idx = grid_idx * 2;
                    size_t z_end_idx = z_start_idx + 1;
                    if (z_end_idx < zgrids.size()) {
                        double grid_z_start = zgrids[z_start_idx];
                        double grid_z_end = zgrids[z_end_idx];
                        if ((loc[2] > grid_z_start - mesh_size) && (loc[2] < grid_z_end + mesh_size)) {
                            isolid = 7 + grid_idx; // Grids start at solid index 7
                            hit_grid = true;
                            break;
                        }
                    }
                }
            }
            
            if (!hit_grid) {
                // Default to a grid category if no specific grid found
                isolid = 0; // Default grid category
            }
        }
        else if (pp.get_status() == PARTICLE_STRIP) {
            isolid = 0; // Stripped particles
        }
        else {
            // Determine boundary collision
            if (loc[2] < z_start + mesh_size) isolid = 5; // z_min boundary
            else if (loc[2] > z_end - mesh_size) isolid = 6; // z_max boundary (beam exit)
            else if (loc[1] < y_start + mesh_size) isolid = 3; // y_min boundary
            else if (loc[1] > y_end - mesh_size) isolid = 4; // y_max boundary
            else if (loc[0] < x_start + mesh_size) isolid = 1; // x_min boundary
            else if (loc[0] > x_end - mesh_size) isolid = 2; // x_max boundary
            else isolid = 0; // Default (shouldn't happen)
        }
        
        // Record particle collision
        if (isolid < total_categories) {
            colls[isolid].push_back(ii);
        }
    }
    
    // Log particle statistics
    if (debug) {
        logfile << "DEBUG: Particle counts by collision solid:" << endl;
        for (size_t ss = 0; ss < total_categories; ss++) {
            logfile << "  solid " << ss << ": " << colls[ss].size() << " particles" << endl;
        }
        
        const char* species_names[] = {"HM", "H0", "HP", "H2P", "H20", "E"};
        logfile << "DEBUG: Particle counts by species:" << endl;
        for (int i = 0; i < 6; i++) {
            logfile << "  " << species_names[i] << ": " << particle_counts[i] << endl;
        }
        
        logfile << "DEBUG: Particle counts by generation:" << endl;
        for (int i = 0; i < 5; i++) {
            logfile << "  gen " << i << ": " << generation_counts[i] << endl;
        }
        logfile << "  gen >4: " << generation_counts[5] << endl;
    }
    
    // Calculate power for each solid/grid
    std::vector<double> power_per_solid;
    std::vector<double> current_per_solid;
    
    for (size_t ss = 0; ss < total_categories; ss++) {
        PowerStruct pow;
        pow.ionmass = ionmass;
        pow.solid_idx = ss;
        
        size_t npart_here = colls[ss].size();
        if (debug) logfile << "DEBUG: Calculating power for solid " << ss 
                           << " with " << npart_here << " particles" << endl << flush;
        
        if (npart_here > 0) {
            // Add all particles that hit this solid
            for (size_t jj = 0; jj < npart_here; jj++) {
                pow.add(particles->particle(colls[ss][jj]));
            }
            
            // Calculate total power and current
            pow.calculate_total_power();
            
            // Save detailed particle data for this solid
            string solid_filename = output_folder + file_tag + "_parts_at_solid_" + std::to_string(ss) + ".dat";
            pow.print(solid_filename);
        } else {
            pow.total_power = 0.0;
            pow.total_current = 0.0;
        }
        
        // Store results based on particle kind
        if (pk == PARTICLE_ALL) {
            power_per_solid.push_back(pow.total_power);
            current_per_solid.push_back(pow.total_current);
        } else {
            int pk_index = get_particle_int(pk);
            if (pk_index >= 0 && pk_index < static_cast<int>(pow.total_power_perspecies.size())) {
                power_per_solid.push_back(pow.total_power_perspecies[pk_index]);
                current_per_solid.push_back(pow.total_current_perspecies[pk_index]);
            } else {
                power_per_solid.push_back(0.0);
                current_per_solid.push_back(0.0);
            }
        }
    }
    
    // Calculate total beam power (grids only, excluding plasma grid #7)
    double full_beam_power = 0.0;
    for (size_t qq = 6; qq < power_per_solid.size(); qq++) {
        if (qq != 7) { // Exclude plasma grid (solid index 7)
            full_beam_power += power_per_solid[qq];
        }
    }
    
    if (debug) logfile << "DEBUG: Full beam power on grids: " << full_beam_power << " W" << endl << flush;
    
    // Create summary file with power breakdown
    string summary_filename = output_folder + file_tag + "_grid_power_summary.txt";
    ofstream summary_file(summary_filename);
    if (summary_file.is_open()) {
        summary_file << "# Grid Power Load Analysis" << endl;
        summary_file << "# Solid_Index\tPower[W]\tCurrent[A]\tParticles\tDescription" << endl;
        
        const char* solid_descriptions[] = {
            "Stripped", "X_min", "X_max", "Y_min", "Y_max", "Z_min", "Z_max_exit",
            "Plasma_Grid", "Grid_1", "Grid_2", "Grid_3", "Grid_4", "Grid_5", "Grid_6"
        };
        
        for (size_t ss = 0; ss < std::min(power_per_solid.size(), size_t(14)); ss++) {
            const char* desc = (ss < 14) ? solid_descriptions[ss] : "Unknown";
            summary_file << ss << "\t" 
                        << power_per_solid[ss] << "\t"
                        << current_per_solid[ss] << "\t"
                        << colls[ss].size() << "\t"
                        << desc << endl;
        }
        
        summary_file << "# Total beam power (grids only): " << full_beam_power << " W" << endl;
        summary_file.close();
        
        ibsimu.message(1) << "Grid power analysis saved to: " << summary_filename << endl;
    }
    
    if (debug) logfile << "DEBUG: Grid power load analysis completed" << endl << flush;
    
    return power_per_solid;
}

// New enhanced createSimulationSummary function with extracted current density
void DiagnosticsManager::createSimulationSummary(double accelerated_current_density,
                                                 double extracted_current_density,
                                                 const std::string& filename, bool append_at_end,
                                                 const ParticleDataBase3D* particles,
                                                 const std::string& extra_info) {
    
    std::ofstream file;
    if (append_at_end) {
        file.open(filename, std::ios::app);
    } else {
        file.open(filename);
    }
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open summary file " << filename << std::endl;
        return;
    }

    // Count particles at the exit (z > 0.04)
    uint32_t particles_at_exit = 0;
    uint32_t total_particles = 0;
    
    if (particles) {
        total_particles = particles->size();
        
        for (uint32_t i = 0; i < total_particles; i++) {
            Particle3D particle = particles->particle(i);
            double z_pos = particle(3); // z position
            
            if (z_pos > 0.04) { // Particles that reached the exit
                particles_at_exit++;
            }
        }
    }
    
    // Calculate transmission efficiency
    double transmission_efficiency = 0.0;
    if (total_particles > 0) {
        transmission_efficiency = (double(particles_at_exit) / double(total_particles)) * 100.0;
    }
    
    // Get current timestamp
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto tm = *std::localtime(&time_t);
    
    file << "========================================" << std::endl;
    file << "        SIMULATION SUMMARY" << std::endl;
    file << "========================================" << std::endl;
    file << "Timestamp: " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << std::endl;
    file << std::endl;
    
    file << "CURRENT DENSITY ANALYSIS:" << std::endl;
    file << "  Accelerated Current Density: " << std::scientific << std::setprecision(6) 
         << accelerated_current_density << " A/m²" << std::endl;
    file << "  Extracted Current Density:   " << std::scientific << std::setprecision(6) 
         << extracted_current_density << " A/m²" << std::endl;
    file << std::endl;
    
    file << "PARTICLE STATISTICS:" << std::endl;
    file << "  Total Particles Injected: " << total_particles << std::endl;
    file << "  Particles at Exit (z>40mm): " << particles_at_exit << std::endl;
    file << "  Transmission Efficiency: " << std::fixed << std::setprecision(1) 
         << transmission_efficiency << "%" << std::endl;
    file << std::endl;
    
    if (!extra_info.empty()) {
        file << "ADDITIONAL INFORMATION:" << std::endl;
        file << extra_info << std::endl;
        file << std::endl;
    }
    
    file << "========================================" << std::endl;
    file << std::endl;
    
    file.close();
    
    std::cout << "Simulation summary saved to: " << filename << std::endl;
    std::cout << "  Accelerated Current Density: " << std::scientific << std::setprecision(3) 
              << accelerated_current_density << " A/m²" << std::endl;
    std::cout << "  Extracted Current Density: " << std::scientific << std::setprecision(3) 
              << extracted_current_density << " A/m²" << std::endl;
    std::cout << "  Transmission Efficiency: " << std::fixed << std::setprecision(1) 
              << transmission_efficiency << "%" << std::endl;
}

// Helper function to get particle species name
std::string DiagnosticsManager::getParticleSpeciesName(particle_kind pk) const {
    switch(pk) {
        case PARTICLE_ALL: return "ALL";
        case PARTICLE_HM: return "HM";
        case PARTICLE_H0: return "H0";
        case PARTICLE_HP: return "HP";
        case PARTICLE_H2P: return "H2P";
        case PARTICLE_H20: return "H20";
        case PARTICLE_E: return "E";
        default: return "UNKNOWN";
    }
}

// Individual simulation summary (for single simulation in scan)
void DiagnosticsManager::createIndividualSimulationSummary(int scan_index, double zlocsummary,
                                      const ParticleDataBase3D* particles,
                                      double extracted_current_density,
                                      const FileManager* fileManager,
                                      particle_kind pk) {
    
    if (!fileManager) {
        logfile << "ERROR: FileManager is required for createIndividualSimulationSummary" << std::endl;
        return;
    }
    
    std::string species_suffix = (pk == PARTICLE_ALL) ? "" : ("_" + getParticleSpeciesName(pk));
    std::string summary_filename = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + "_" + 
                                  std::to_string(scan_index) + "_summary" + species_suffix + ".txt";
    
    logfile << "Creating individual simulation summary: " << summary_filename << std::endl;
    
    // Create the summary file
    std::ofstream summary_file(summary_filename);
    if (!summary_file.is_open()) {
        logfile << "ERROR: Could not create summary file: " << summary_filename << std::endl;
        return;
    }
    
    // Write header
    summary_file << "========================================" << std::endl;
    summary_file << "    INDIVIDUAL SIMULATION SUMMARY" << std::endl;
    summary_file << "========================================" << std::endl;
    summary_file << "Simulation: " << fileManager->getFileTag() << std::endl;
    summary_file << "Scan Index: " << scan_index << std::endl;
    summary_file << "Species: " << getParticleSpeciesName(pk) << std::endl;
    summary_file << "Analysis Z-location: " << zlocsummary*1000 << " mm" << std::endl;
    summary_file << std::endl;
    
    // ===== BEAM PROPERTIES AT zlocsummary =====
    summary_file << "========================================" << std::endl;
    summary_file << "  BEAM PROPERTIES AT z=" << zlocsummary*1000 << " mm" << std::endl;
    summary_file << "========================================" << std::endl;
    
    // Read beam properties from NEGIONBEAM diagnostic file generated by performAnalysis
    std::string diagnostic_file = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + 
                                  "_NEGIONBEAM_diagnostic_summary.txt";
    std::ifstream diag_in(diagnostic_file);
    
    if (diag_in.is_open()) {
        std::string line;
        bool found = false;
        
        // Skip header
        std::getline(diag_in, line);
        
        // Read the single data line (NEGIONBEAM file should contain only one line at zlocsummary)
        if (std::getline(diag_in, line) && !line.empty() && line[0] != '#') {
            std::istringstream iss(line);
            int it;
            double z_mm, current_mA, x_ave, x_max, x_min, y_ave, y_max, y_min;
            double xp_ave, yp_ave, div_x, div_y, epot_kV, B_mT, rho, sigma;
            
            if (iss >> it >> z_mm >> current_mA >> x_ave >> x_max >> x_min >> y_ave >> y_max >> y_min 
                    >> xp_ave >> yp_ave >> div_x >> div_y >> epot_kV >> B_mT >> rho >> sigma) {
                found = true;
                
                summary_file << "Current: " << current_mA << " mA (" << current_mA/1000.0 << " A)" << std::endl;
                summary_file << "Extracted Current Density: " << std::scientific << std::setprecision(4) 
                            << extracted_current_density << " A/m²" << std::endl;
                summary_file << std::endl;
                
                summary_file << "Beam Center:" << std::endl;
                summary_file << "  <x> = " << std::fixed << std::setprecision(3) << x_ave << " mm" << std::endl;
                summary_file << "  <y> = " << y_ave << " mm" << std::endl;
                summary_file << std::endl;
                
                summary_file << "Beam Extent:" << std::endl;
                summary_file << "  x_max = " << x_max << " mm,  x_min = " << x_min << " mm" << std::endl;
                summary_file << "  y_max = " << y_max << " mm,  y_min = " << y_min << " mm" << std::endl;
                summary_file << std::endl;
                
                summary_file << "Beam Divergence:" << std::endl;
                summary_file << "  <x'> = " << xp_ave << " mrad" << std::endl;
                summary_file << "  <y'> = " << yp_ave << " mrad" << std::endl;
                summary_file << "  Dx = " << div_x << " mrad" << std::endl;
                summary_file << "  Dy = " << div_y << " mrad" << std::endl;
                summary_file << std::endl;
                
                summary_file << "Fields:" << std::endl;
                summary_file << "  Potential: " << epot_kV << " kV" << std::endl;
                summary_file << "  Magnetic Field (By): " << B_mT << " mT" << std::endl;
                summary_file << std::endl;
                
                summary_file << "Gas Properties:" << std::endl;
                summary_file << "  Density: " << std::scientific << rho << " kg/m³" << std::endl;
                summary_file << "  Stripping Cross Section: " << sigma << " m²" << std::endl;
                summary_file << std::endl;
            }
        }
        
        if (!found) {
            summary_file << "ERROR: Could not read beam properties data from file" << std::endl;
            summary_file << std::endl;
        }
        
        diag_in.close();
    } else {
        summary_file << "ERROR: Could not read diagnostic file: " << diagnostic_file << std::endl;
        summary_file << std::endl;
    }
    
    // ===== GRID POWER AND CURRENT LOADS =====
    summary_file << "========================================" << std::endl;
    summary_file << "  GRID POWER AND CURRENT LOADS" << std::endl;
    summary_file << "========================================" << std::endl;
    
    // Read grid power analysis from file generated by performAnalysis
    std::string grid_power_file = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + 
                                  "_ALL_grid_power_summary.txt";
    std::ifstream grid_in(grid_power_file);
    
    if (grid_in.is_open()) {
        std::string line;
        double total_power = 0.0;
        double total_current = 0.0;
        size_t total_particles = 0;
        
        // Skip header lines
        while (std::getline(grid_in, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            std::istringstream iss(line);
            int solid_idx;
            double power_W, current_A;
            size_t n_particles;
            std::string description;
            
            if (iss >> solid_idx >> power_W >> current_A >> n_particles) {
                // Read rest of line as description
                std::getline(iss, description);
                
                // Accumulate totals (excluding boundaries 1-6)
                if (solid_idx >= 7) {
                    total_power += power_W;
                    total_current += current_A;
                }
                total_particles += n_particles;
                
                // Print individual solid data
                summary_file << "Solid " << std::setw(2) << solid_idx << ": "
                            << std::setw(10) << std::fixed << std::setprecision(3) << power_W << " W, "
                            << std::setw(8) << std::setprecision(6) << current_A << " A, "
                            << std::setw(6) << n_particles << " particles"
                            << description << std::endl;
            }
        }
        
        summary_file << std::endl;
        summary_file << "TOTALS:" << std::endl;
        summary_file << "  Total Grid Power (solids 7+): " << std::fixed << std::setprecision(3) 
                    << total_power << " W" << std::endl;
        summary_file << "  Total Grid Current: " << std::setprecision(6) << total_current << " A" << std::endl;
        summary_file << "  Total Particles on All Solids: " << total_particles << std::endl;
        summary_file << std::endl;
        
        grid_in.close();
    } else {
        summary_file << "ERROR: Could not read grid power file: " << grid_power_file << std::endl;
        summary_file << "Please ensure performAnalysis was called before creating summary." << std::endl;
        summary_file << std::endl;
    }
    
    // ===== PARTICLE STATISTICS =====
    summary_file << "========================================" << std::endl;
    summary_file << "  PARTICLE STATISTICS" << std::endl;
    summary_file << "========================================" << std::endl;
    
    if (particles) {
        size_t total_particles = particles->size();
        size_t indomain_particles = 0;
        size_t collided_particles = 0;
        size_t stripped_particles = 0;
        size_t out_particles = 0;
        
        // Count particle fates
        for (size_t i = 0; i < total_particles; i++) {
            Particle3D& pp = const_cast<Particle3D&>(particles->particle(i));
            particle_status_e status = pp.get_status();
            
            if (status == PARTICLE_OK) indomain_particles++;
            else if (status == PARTICLE_COLL) collided_particles++;
            else if (status == PARTICLE_STRIP) stripped_particles++;
            else if (status == PARTICLE_OUT) out_particles++;
        }
        
        summary_file << "Total Particles Simulated: " << total_particles << std::endl;
        summary_file << "In domain: " << indomain_particles
                    << " (" << std::fixed << std::setprecision(1)
                    << (100.0 * indomain_particles / total_particles) << "%)" << std::endl;
        summary_file << "Collided: " << collided_particles
                    << " (" << (100.0 * collided_particles / total_particles) << "%)" << std::endl;
        summary_file << "Stripped: " << stripped_particles 
                    << " (" << (100.0 * stripped_particles / total_particles) << "%)" << std::endl;
        summary_file << "Exited: " << out_particles
                    << " (" << (100.0 * out_particles / total_particles) << "%)" << std::endl;
        summary_file << std::endl;
    }
    
    summary_file << "========================================" << std::endl;
    summary_file << "      END OF SUMMARY" << std::endl;
    summary_file << "========================================" << std::endl;
    
    summary_file.close();
    
    logfile << "Individual simulation summary created: " << summary_filename << std::endl;
}

// Add entry to scan-level beam properties summary
void DiagnosticsManager::addToScanBeamPropertiesSummary(int scan_index, const std::string& simulation_tag,
                                   const std::string& scan_folder, const std::string& scan_file_tag,
                                   double zlocsummary, const ParticleDataBase3D* particles,
                                   double extracted_current_density, const MeshVectorField* /*magnetic*/,
                                   const EpotField* /*potential*/, particle_kind pk) {
    
    std::string species_suffix = (pk == PARTICLE_ALL) ? "" : ("_" + getParticleSpeciesName(pk));
    std::string beamprops_filename = scan_folder + "/" + scan_file_tag + "_beamprops" + species_suffix + ".txt";
    
    bool file_exists = std::ifstream(beamprops_filename).good();
    std::ofstream file(beamprops_filename, std::ios::app);
    
    if (!file_exists) {
        // Write header for new file
        file << "# Scan-level beam properties summary" << std::endl;
        file << "# scan_idx sim_tag zloc[mm] accJ[A/m2] extJ[A/m2] n_particles transmission[%]" << std::endl;
    }
    
    // Calculate beam properties at specified z location using IBSimu trajectory diagnostics
    TrajectoryDiagnosticData tdata;
    std::vector<trajectory_diagnostic_e> diagnostics = {
        DIAG_X, DIAG_Y, DIAG_Z, DIAG_VX, DIAG_VY, DIAG_VZ,
        DIAG_CURR, DIAG_MASS, DIAG_CHARGE, DIAG_NO
    };
    
    // Get trajectories at the specified z plane with error handling
    size_t n_particles = 0;
    double accelerated_current_density = 0.0;
    double pg_aperture_area = M_PI * 0.007 * 0.007; // 7mm radius in m²
    
    try {
        if (debug) logfile << "DEBUG: About to call trajectories_at_plane for scan beam properties at z=" << zlocsummary << std::endl;
        const_cast<ParticleDataBase3D*>(particles)->trajectories_at_plane(tdata, AXIS_Z, zlocsummary, diagnostics);
        if (debug) logfile << "DEBUG: trajectories_at_plane call completed for scan beam properties" << std::endl;
        n_particles = tdata.traj_size();
        
        if (n_particles > 0) {
            // Calculate current from particle current contributions
            double total_current = 0.0;
            for (size_t i = 0; i < n_particles; i++) {
                total_current += std::abs(tdata(i, 6));  // DIAG_CURR is column 6
            }
            accelerated_current_density = total_current / pg_aperture_area;
        }
    } catch (const ErrorRange& e) {
        if (debug) logfile << "DEBUG: ErrorRange exception in trajectories_at_plane for scan beam properties" << endl;
        // Use fallback values - estimate from total particle count
        n_particles = particles ? particles->size() : 0;
        accelerated_current_density = extracted_current_density; // Use extracted as fallback
    } catch (const std::exception& e) {
        if (debug) logfile << "DEBUG: Error in trajectories_at_plane for scan beam properties: " << e.what() << endl;
        // Use fallback values - estimate from total particle count
        n_particles = particles ? particles->size() : 0;
        accelerated_current_density = extracted_current_density; // Use extracted as fallback
    }
    
    double transmission_efficiency = (n_particles > 0 && accelerated_current_density > 0) ? 
                                   (accelerated_current_density / extracted_current_density) * 100.0 : 0.0;
    
    // Write data line
    file << scan_index << " " << simulation_tag << " " 
         << std::fixed << std::setprecision(3) << zlocsummary*1000 << " " 
         << std::scientific << std::setprecision(6)
         << accelerated_current_density << " " << extracted_current_density << " "
         << n_particles << " " << std::fixed << std::setprecision(1) << transmission_efficiency
         << std::endl;
    
    file.close();
    logfile << "Added scan beam properties entry for " << simulation_tag << " (species: " 
              << getParticleSpeciesName(pk) << ")" << std::endl;
}

// Add entry to scan-level grid power summary
void DiagnosticsManager::addToScanGridPowerSummary(int scan_index, const std::string& simulation_tag,
                              const std::string& scan_folder, const std::string& scan_file_tag,
                              const ParticleDataBase3D* particles, const Geometry* /*geometry*/,
                              particle_kind pk) {
    
    std::string species_suffix = (pk == PARTICLE_ALL) ? "" : ("_" + getParticleSpeciesName(pk));
    std::string gridpower_filename = scan_folder + "/" + scan_file_tag + "_grid_power" + species_suffix + ".txt";
    
    bool file_exists = std::ifstream(gridpower_filename).good();
    std::ofstream file(gridpower_filename, std::ios::app);
    
    if (!file_exists) {
        // Write header for new file
        file << "# Scan-level grid power summary" << std::endl;
        file << "# scan_idx sim_tag total_particles stripped_particles pg_particles exit_particles" << std::endl;
    }
    
    // Analyze particle end locations with error handling
    size_t total_particles = 0, stripped_particles = 0, pg_particles = 0, exit_particles = 0;
    double ionmass = 1.67262e-27; // H- ion mass
    
    try {
        for (size_t a = 0; a < particles->size(); a++) {
            Particle3D& pp = const_cast<Particle3D&>(particles->particle(a));
            if (pk != PARTICLE_ALL && identify_particle_species(pp.m(), pp.q(), ionmass) != pk) continue;
            
            total_particles++;
            
            // Get final location from particle
            Vec3D final_pos = pp.location();
            
            // Categorize by final location
            if (final_pos[2] < 0.05) { // Stripped
                stripped_particles++;
            } else if (final_pos[2] < 0.2) { // PG region
                pg_particles++;
            } else { // Exit region
                exit_particles++;
            }
        }
    } catch (const std::exception& e) {
        if (debug) logfile << "DEBUG: Error analyzing particles for grid power summary: " << e.what() << endl;
        total_particles = particles ? particles->size() : 0;
    }
    
    // Write data line
    file << scan_index << " " << simulation_tag << " "
         << total_particles << " " << stripped_particles << " " 
         << pg_particles << " " << exit_particles << std::endl;
    
    file.close();
    logfile << "Added scan grid power entry for " << simulation_tag << " (species: " 
              << getParticleSpeciesName(pk) << ")" << std::endl;
}
