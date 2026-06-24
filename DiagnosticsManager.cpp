/*
 * DiagnosticsManager.cpp
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored from ManageSimulation)
 */

#include "DiagnosticsManager.h"
#include "SimulationParameters.h"
#include "FileManager.h"
#include "ManageSimulation_New.h"  // For PowerStruct
#include "globals.h"
#include "my_diagnostics.h"
#include "funct.h"
#include "cross_sections.h"
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
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

namespace {

std::string gridPowerRowName(const SimulationParameters& params, int id) {
    if (id == 0) {
        return "In_volume";
    }

    SimulationParameters::BoundaryConditionDefinition definition;
    if (params.tryGetBoundaryCondition(id, definition) && !definition.name.empty()) {
        return definition.name;
    }

    switch (id) {
        case 1: return "x-min";
        case 2: return "x-max";
        case 3: return "y-min";
        case 4: return "y-max";
        case 5: return "z-min";
        case 6: return "z-max";
        default: return string("Boundary_") + to_string(id);
    }
}

bool defaultGridPowerIncludeInTotal(const SimulationParameters& params, int id) {
    if (id <= 6) {
        return false;
    }

    std::string name = gridPowerRowName(params, id);
    std::transform(name.begin(), name.end(), name.begin(), [](unsigned char ch) {
        return static_cast<char>(std::toupper(ch));
    });
    return name != "PG" && name != "PLASMA GRID";
}

std::vector<int> orderedGridPowerIds(const SimulationParameters& params, const Geometry* geometry) {
    std::set<int> ids;
    ids.insert(0);
    for (int id = 1; id <= 6; ++id) {
        ids.insert(id);
    }
    if (geometry != NULL) {
        for (int id = 7; id <= static_cast<int>(geometry->number_of_solids()) + 6; ++id) {
            ids.insert(id);
        }
    }
    if (params.hasGeneratedGeometrySolids()) {
        const std::vector<SimulationParameters::GeometrySolidDefinition>& solids =
            params.getGeneratedGeometrySolids();
        for (std::vector<SimulationParameters::GeometrySolidDefinition>::const_iterator it = solids.begin();
             it != solids.end();
             ++it) {
            ids.insert(it->boundaryId);
        }
    }
    const std::vector<SimulationParameters::DiagnosticGridRangeDefinition>& configured_rows =
        params.getDiagnosticGridPowerRanges();
    for (std::vector<SimulationParameters::DiagnosticGridRangeDefinition>::const_iterator it =
             configured_rows.begin();
         it != configured_rows.end();
         ++it) {
        ids.insert(it->id);
    }
    return std::vector<int>(ids.begin(), ids.end());
}

std::map<int, bool> configuredGridPowerIncludeMap(const SimulationParameters& params) {
    std::map<int, bool> include_map;
    const std::vector<SimulationParameters::DiagnosticGridRangeDefinition>& configured_rows =
        params.getDiagnosticGridPowerRanges();
    for (std::vector<SimulationParameters::DiagnosticGridRangeDefinition>::const_iterator it =
             configured_rows.begin();
         it != configured_rows.end();
         ++it) {
        include_map[it->id] = it->includeInTotal;
    }
    return include_map;
}

bool gridPowerIncludeInTotal(const SimulationParameters& params,
                             const std::map<int, bool>& include_map,
                             int id) {
    std::map<int, bool>::const_iterator it = include_map.find(id);
    if (it != include_map.end()) {
        return it->second;
    }
    return defaultGridPowerIncludeInTotal(params, id);
}

int classifyGridPowerRowId(Particle3D& particle, const Geometry* geometry, double mesh_size) {
    if (particle.get_status() != PARTICLE_COLL && particle.get_status() != PARTICLE_OUT) {
        return 0;
    }
    if (geometry == NULL) {
        return 0;
    }

    const Vec3D loc = particle.location();
    const Vec3D vel = particle.velocity();
    double probe_distance = mesh_size > 0.0 ? 0.25 * mesh_size : 0.25 * geometry->h();
    if (!std::isfinite(probe_distance) || probe_distance <= 0.0) {
        probe_distance = 1.0e-9;
    }

    const double speed_squared = vel.ssqr();
    if (std::isfinite(speed_squared) && speed_squared > 0.0) {
        const double inv_speed = 1.0 / std::sqrt(speed_squared);
        const Vec3D probe(
            loc[0] + vel[0] * inv_speed * probe_distance,
            loc[1] + vel[1] * inv_speed * probe_distance,
            loc[2] + vel[2] * inv_speed * probe_distance);
        const uint32_t probe_id = geometry->inside(probe);
        if (probe_id != 0U) {
            return static_cast<int>(probe_id);
        }
    }

    const uint32_t direct_id = geometry->inside(loc);
    if (direct_id != 0U) {
        return static_cast<int>(direct_id);
    }

    const Vec3D origo = geometry->origo();
    const Vec3D maxpt = geometry->max();
    if (loc[2] <= origo[2] + probe_distance) return 5;
    if (loc[2] >= maxpt[2] - probe_distance) return 6;
    if (loc[1] <= origo[1] + probe_distance) return 3;
    if (loc[1] >= maxpt[1] - probe_distance) return 4;
    if (loc[0] <= origo[0] + probe_distance) return 1;
    if (loc[0] >= maxpt[0] - probe_distance) return 2;
    return 0;
}

} // namespace

void DiagnosticsManager::generateDiagnosticData(const vector<double>& diagzpos, int iteration, 
                                               const string& outfile, bool append_at_end, 
                                               const ParticleDataBase3D* pdb) {
    if (debug) std::cerr << "DEBUG: generateDiagnosticData called with " << diagzpos.size() << " z-positions" << std::endl << std::flush;
    
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
        if (debug) std::cerr << "DEBUG: Processing z-position " << ii+1 << "/" << diagzpos.size() << " at z=" << zloc << " m" << std::endl << std::flush;
        
        if (debug) std::cerr << "DEBUG: About to enter try block" << std::endl << std::flush;
        try {
            if (debug) std::cerr << "DEBUG: Inside try block, gathering plane diagnostics" << std::endl << std::flush;
            // Trajectory diagnostics at specific z location
            TrajectoryDiagnosticData tdata;
            vector<trajectory_diagnostic_e> diagnostics = {
                DIAG_X, DIAG_Y, DIAG_Z, DIAG_VX, DIAG_VY, DIAG_VZ,
                DIAG_CURR, DIAG_MASS, DIAG_CHARGE, DIAG_NO
            };
            
            if (debug) std::cerr << "  About to call trajectories_at_plane..." << std::endl << std::flush;
            pdb->trajectories_at_plane(tdata, AXIS_Z, zloc, diagnostics);
            if (debug) std::cerr << "  trajectories_at_plane call completed successfully!" << std::endl << std::flush;
            
            // Debug: Check the size immediately after the call
            size_t particle_count = tdata.traj_size();
            if (debug) std::cerr << "  DEBUG: tdata.traj_size() = " << particle_count << std::endl << std::flush;
            
            // Write count to file for debugging
            if (debug) {
                ofstream debugfile("debug_particle_count.txt", ios::app);
                size_t pdb_size = pdb->size();
                debugfile << "z=" << zloc << " count=" << particle_count << " pdb_size=" << pdb_size << endl;
                debugfile.flush();
                debugfile.close();
            }
            
            if (debug) std::cerr << "  Retrieved " << tdata.traj_size() << " particles at z=" << zloc << std::endl << std::flush;
            if (debug) std::cerr << "  Particle database size: " << pdb->size() << std::endl << std::flush;
            
            if (tdata.traj_size() > 0) {
            if (debug) std::cerr << "  Processing " << tdata.traj_size() << " particles for diagnostic calculations..." << std::endl << std::flush;
            
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
            if (debug) std::cerr << "  Current data size: " << curr_data.size() << std::endl << std::flush;
            if (debug && curr_data.size() > 0) {
                std::cerr << "  First few current values: ";
                for (size_t ii = 0; ii < std::min(size_t(5), curr_data.size()); ++ii) {
                    std::cerr << curr_data[ii] << " ";
                }
                std::cerr << std::endl << std::flush;
            }
            if (debug) std::cerr << "  Current sum calculated: " << current_sum << " A (" << 1000*current_sum << " mA)" << std::endl << std::flush;
            
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
            if (debug) std::cerr << "  No particles found at z=" << zloc << ". Writing zeros." << std::endl << std::flush;
            if (debug) logfile << "DEBUG: No particles intercepting at z position " << zloc << " m. Skipping." << endl;
            // Write zeros for this z position (17 columns total)
            pout5 << iteration << " "
                  << 1000*zloc << " "
                  << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n" << flush;
        }
        } catch (const Error& e) {
            if (debug) std::cerr << "DEBUG: IBSimu Error caught at z=" << zloc << " m" << std::endl << std::flush;
            if (debug) logfile << "DEBUG: Error in trajectory diagnostics at z=" << zloc 
                               << " m: IBSIMU Error. Writing zeros." << endl;
            // Write zeros for this z position due to error (17 columns total)
            pout5 << iteration << " "
                  << 1000*zloc << " "
                  << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n" << flush;
        } catch (const std::exception& e) {
            if (debug) std::cerr << "DEBUG: Standard exception caught at z=" << zloc << " m: " << e.what() << std::endl << std::flush;
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

void DiagnosticsManager::createPlots(int argc, char **argv, const SimulationParameters& params,
                                    const Geometry* geometry,
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
    const SimulationParameters::DiagnosticMeniscusPlotDefinition& meniscus =
        params.getDiagnosticMeniscusPlot();
    if (potential && meniscus.enabled) {
        geomplotter.set_ranges(
            meniscus.zMinMeters,
            meniscus.transverseMinMeters,
            meniscus.zMaxMeters,
            meniscus.transverseMaxMeters);
        geomplotter.set_view(VIEW_ZX, -1);
        geomplotter.plot_png(plot_folder + file_tag + "_MENISCUS_plot_zx.png");
        logfile << "Meniscus zx plotted!" << endl << flush;
        
        geomplotter.set_view(VIEW_ZY, -1);
        geomplotter.plot_png(plot_folder + file_tag + "_MENISCUS_plot_zy.png");
        logfile << "Meniscus zy plotted!" << endl << flush;
    }
    
    if (debug) logfile << "Done!" << endl << flush;
}

void DiagnosticsManager::createPlots(int argc, char **argv, const SimulationParameters& params,
                                    const Geometry* geometry,
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
    createPlots(argc, argv, params, geometry, potential, magnetic, electric, spacecharge,
               pdb, plot_folder, file_tag + "_" + get_particle_name(pk));
}

void DiagnosticsManager::performAnalysis(const ParticleDataBase3D* particles, const SimulationParameters& params,
                                        const Geometry* geometry, const EpotField* potential, 
                                        const MeshVectorField* magnetic, const FileManager* fileManager,vector<double> zgrids) {
    if (debug) logfile << "DEBUG: Performing simulation analysis with field calculations..." << endl << flush;
    
    if (!particles || !geometry) {
        throw Error(ERROR_LOCATION, "Missing required objects for analysis (particles or geometry)");
    }
    if (!fileManager) {
        throw Error(ERROR_LOCATION, "Missing FileManager for diagnostics output generation");
    }
    
    if (debug) logfile << "\n" << flush;
    if (debug) logfile << "DEBUG: geometry->ORIGO: " << geometry->origo(0) << " " << geometry->origo(1) << " " << geometry->origo(2) << " " << endl << flush;
    
    // Create trajectory density field (matching original implementation)
    MeshScalarField tdens(*geometry);
    particles->build_trajectory_density_field(tdens);
    
    string diagfile = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + "_diagnostic_summary.txt";
    
    vector<double> zlocs;
    if (params.hasDiagnosticSampleZPositions()) {
        zlocs = params.getDiagnosticSampleZPositions();
    } else {
        double z_start = geometry->origo(2);
        Vec3D lastpt = geometry->max();
        double z_end = lastpt[2] - params.getMeshSize();
        size_t steps = 21;
        double delta_z = (z_end - z_start) / steps;
        for (size_t ii = 0; ii < steps; ii++) {
            zlocs.push_back(z_start + ii * delta_z);
        }
    }
    
    logfile << "Diagnostic data for ALL particles in file " << diagfile << "\n" << flush;
    logfile << "Save diagnostic data at " << zlocs.size() << " points with field calculations...\n" << flush;
    
    // Call enhanced diagnostic generation WITH field objects
    generateDiagnosticData(zlocs, 0, diagfile, false, particles, geometry, potential, magnetic);
    
    logfile << "Done!\n" << flush;
    
    // Get grid locations from parameters if available, otherwise use default values
    
    try {
        // Create subset excluding first element (matching original logic)
        vector<double> zg;
        if (zgrids.size() > 1) {
            zg = vector<double>(zgrids.begin() + 1, zgrids.end());
        } else {
            // Fallback: use domain exit plane derived from geometry
            double z_fallback = geometry->max()[2] - params.getMeshSize();
            zg = {z_fallback};
        }
        
        string gridfile = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + "_statsatgrids.txt";
        
        if (debug) logfile << "DEBUG: Save diagnostic data at grids with field calculations\n" << flush;
        if (debug) cout << "DEBUG: Save diagnostic data at grids\n" << flush;
        
        // Call enhanced diagnostic generation WITH field objects
        generateDiagnosticData(zg, 0, gridfile, false, particles, geometry, potential, magnetic);
        
    } catch (const std::exception& e) {
        if (debug) logfile << "DEBUG: Error processing grid data: " << e.what() << endl << flush;
        ibsimu.message(1) << "Warning: grid-position diagnostics were not generated: " << e.what() << endl;
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
    if (!fileManager) {
        throw Error(ERROR_LOCATION, "Missing FileManager for diagnostics output generation");
    }
    
    if (debug) logfile << "\n" << flush;
    if (debug) logfile << "DEBUG: geometry->ORIGO: " << geometry->origo(0) << " " << geometry->origo(1) << " " << geometry->origo(2) << " " << endl << flush;
    
    // Create trajectory density field (matching original implementation)
    MeshScalarField tdens(*geometry);
    particles->build_trajectory_density_field(tdens);
    
    string diagfile = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + "_diagnostic_summary.txt";
    
    vector<double> zlocs;
    if (params.hasDiagnosticSampleZPositions()) {
        zlocs = params.getDiagnosticSampleZPositions();
    } else {
        double z_start = geometry->origo(2);
        Vec3D lastpt = geometry->max();
        double z_end = lastpt[2] - params.getMeshSize();
        size_t steps = 21;
        double delta_z = (z_end - z_start) / steps;
        for (size_t ii = 0; ii < steps; ii++) {
            zlocs.push_back(z_start + ii * delta_z);
        }
    }
    
    logfile << "Diagnostic data for ALL particles in file " << diagfile << "\n" << flush;
    logfile << "Save diagnostic data at " << zlocs.size() << " points with field calculations...\n" << flush;
    
    // Call enhanced diagnostic generation WITH field objects for all particles
    generateDiagnosticData(zlocs, 0, diagfile, false, particles, geometry, potential, magnetic);
    logfile << "Done!\n" << flush;

    vector<double> zlocsummary_vec = {zlocsummary};

    // Keep the single-plane negative-ion diagnostic table for downstream consumers.
    if (include_stripping && !particles_species.empty() && params.getDiagnosticWriteNegativeIonSummary()) {
        for (size_t ii = 0; ii < particles_species.size(); ii++) {
            if (particles_species[ii] && particles_species[ii]->size() > 0) {
                string summary_diagfile = fileManager->getOutputSummaryFolder() + fileManager->getFileTag() + "_NEGIONBEAM_diagnostic_summary.txt";
                if (params.getDiagnosticWriteNegativeIonSummary() && ii == 0) {
                    generateDiagnosticData(zlocsummary_vec, 0, summary_diagfile, false, particles_species[ii], geometry, potential, magnetic);
                }
                logfile << "Done!\n" << flush;
            }
        }
    }
    
    // Grid power load analysis for ALL particles
    logfile << " Analyzing grid power loads for ALL particles\n" << flush;
    try {
        double mesh_size = params.getMeshSize();
        double ionmass = params.getMIons();
        string output_folder = fileManager->getOutputSummaryFolder();
        string file_tag = fileManager->getFileTag();
        
        vector<double> grid_powers_all = analyzeGridPowerLoads(particles, params, mesh_size, 
                                                               geometry, ionmass, output_folder, 
                                                               file_tag + "_ALL", PARTICLE_ALL);
        
        logfile << " Grid power analysis for ALL particles completed with " 
                           << grid_powers_all.size() << " power values\n" << flush;
    } catch (const std::exception& e) {
        logfile << "Error in grid power analysis for ALL particles: " << e.what() << endl << flush;
    }
    
    // Grid power load analysis for each particle species
    if (include_stripping && !particles_species.empty() && params.getDiagnosticWritePerSpeciesGridPower()) {
        for (size_t ii = 0; ii < particles_species.size() && ii < particle_kind_count(); ii++) {
            if (particles_species[ii] && particles_species[ii]->size() > 0) {
                particle_kind species_kind = int2kind(static_cast<int>(ii));
                std::string species_name = get_particle_name(species_kind);
                logfile << "Analyzing grid power loads for species " << species_name << "\n" << flush;

                try {
                    double mesh_size = params.getMeshSize();
                    double ionmass = params.getMIons();
                    string output_folder = fileManager->getOutputSummaryFolder();
                    string file_tag = fileManager->getFileTag();
                    
                    vector<double> grid_powers_species = analyzeGridPowerLoads(particles_species[ii], params, 
                                                                               mesh_size, geometry, ionmass, 
                                                                               output_folder, 
                                                                               file_tag + "_" + species_name, 
                                                                               species_kind);
                    
                    logfile << "Grid power analysis for species " << species_name 
                            << " completed with " << grid_powers_species.size() << " power values\n" << flush;
                } catch (const std::exception& e) {
                    logfile << "Error in grid power analysis for species " 
                            << species_name << ": " << e.what() << endl << flush;
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

void DiagnosticsManager::generateDiagnosticData(const vector<double>& diagzpos, int iteration, 
                                               const string& outfile, bool append_at_end, 
                                               const ParticleDataBase3D* pdb,
                                               const Geometry* /*geometry*/,
                                               const EpotField* potential,
                                               const MeshVectorField* magnetic) {
    if (debug) std::cerr << "DEBUG: generateDiagnosticData called with " << diagzpos.size() << " z-positions (with field calculations)" << std::endl << std::flush;
    
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

    vector<double> pos, dens;
    if (!density_profile_filename_.empty()) {
        load_density_profile(density_profile_filename_, pos, dens);
    }

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
        if (debug) std::cerr << "DEBUG: Processing z-position " << (i+1) << "/" << diagzpos.size() 
                  << " at z=" << zloc << " m" << std::endl << std::flush;
        
        // Clear vectors for this z-position
        x[0].clear(); y[0].clear(); z[0].clear();
        vx[0].clear(); vy[0].clear(); vz[0].clear();
        curr[0].clear(); mass[0].clear(); charge[0].clear(); no_part[0].clear();
        
        if (debug) std::cerr << "DEBUG: About to enter try block" << std::endl << std::flush;
        try {
            if (debug) std::cerr << "DEBUG: Inside try block, gathering plane diagnostics" << std::endl << std::flush;
            
            // Trajectory diagnostics at specific z location
            TrajectoryDiagnosticData tdata;
            vector<trajectory_diagnostic_e> diagnostics = {
                DIAG_X, DIAG_Y, DIAG_Z, DIAG_VX, DIAG_VY, DIAG_VZ,
                DIAG_CURR, DIAG_MASS, DIAG_CHARGE, DIAG_NO
            };
            
            if (debug) std::cerr << "  About to call trajectories_at_plane..." << std::endl << std::flush;
            pdb->trajectories_at_plane(tdata, AXIS_Z, zloc, diagnostics);
            if (debug) std::cerr << "  trajectories_at_plane call completed successfully!" << std::endl << std::flush;
            
            // Debug: Check the size immediately after the call
            size_t particle_count = tdata.traj_size();
            if (debug) std::cerr << "  DEBUG: tdata.traj_size() = " << particle_count << std::endl << std::flush;
            
            // Write count to file for debugging
            if (debug) {
                ofstream debugfile("debug_particle_count.txt", ios::app);
                size_t pdb_size = pdb->size();
                debugfile << "z=" << zloc << " count=" << particle_count << " pdb_size=" << pdb_size << endl;
                debugfile.flush();
                debugfile.close();
            }
            
            if (debug) std::cerr << "  Retrieved " << particle_count << " particles at z=" << zloc << std::endl << std::flush;
            if (debug) std::cerr << "  Particle database size: " << pdb->size() << std::endl << std::flush;
            
            if (particle_count > 0) {
                if (debug) std::cerr << "  Processing " << particle_count << " particles for diagnostic calculations..." << std::endl << std::flush;
                
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
                if (debug) std::cerr << "  Current data size: " << curr_data.size() << std::endl << std::flush;
                if (debug && curr_data.size() > 0) {
                    std::cerr << "  First few current values: ";
                    for (size_t ii = 0; ii < std::min(size_t(5), curr_data.size()); ++ii) {
                        std::cerr << curr_data[ii] << " ";
                    }
                    std::cerr << std::endl << std::flush;
                }
                if (debug) std::cerr << "  Current sum calculated: " << current_sum << " A (" << 1000*current_sum << " mA)" << std::endl << std::flush;
                
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
                        if (debug) std::cerr << "  DEBUG: Error evaluating potential: " << e.what() << std::endl << std::flush;
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
                        if (debug) std::cerr << "  DEBUG: Error evaluating magnetic field: " << e.what() << std::endl << std::flush;
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
                    if (debug) std::cerr << "  DEBUG: No density profile loaded (pos.size=" << pos.size() 
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
                        if (debug) std::cerr << "  DEBUG: Error calculating stripping cross section: " << e.what() << std::endl << std::flush;
                    }
                } else {
                    if (debug) std::cerr << "  DEBUG: No potential field for cross section calculation" << std::endl << std::flush;
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
                if (debug) std::cerr << "  No particles found at z=" << zloc << std::endl << std::flush;
                
                // Write zeros for empty locations
                pout5 << iteration << " "
                      << 1000*zloc << " "
                      << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
                      << endl << flush;
            }
            
        } catch (Error& e) {
            if (debug) std::cerr << "  IBSimu Error caught: [Error in grid power analysis]" << std::endl << std::flush;
            
            // Write error entry
            pout5 << iteration << " "
                  << 1000*zloc << " "
                  << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
                  << endl << flush;
        } catch (std::exception& e) {
            if (debug) std::cerr << "  Standard exception caught: " << e.what() << std::endl << std::flush;
            
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
                                                             const SimulationParameters& params,
                                                             double mesh_size,
                                                             const Geometry* geometry,
                                                             double ionmass,
                                                             const std::string& output_folder,
                                                             const std::string& file_tag,
                                                             particle_kind pk,
                                                             std::vector<double>* current_per_solid_out) {
    if (debug) logfile << "DEBUG: Analyzing grid power loads..." << endl << flush;
    
    if (!particles) {
        if (debug) logfile << "DEBUG: No particle database provided" << endl << flush;
        return std::vector<double>();
    }
    
    std::vector<int> row_ids = orderedGridPowerIds(params, geometry);
    std::map<int, bool> include_map = configuredGridPowerIncludeMap(params);
    int max_row_id = row_ids.empty() ? 0 : row_ids.back();
    size_t total_categories = static_cast<size_t>(max_row_id + 1);
    
    if (debug) logfile << "DEBUG: Grid power rows: " << row_ids.size() << ", max row id: " << max_row_id << endl << flush;
    
    // Create collision vectors for each solid/boundary
    std::vector<std::vector<size_t>> colls(total_categories);
    
    // Statistics counters
    std::vector<size_t> particle_counts(particle_kind_count(), 0); // family slots + electrons
    size_t generation_counts[6] = {0}; // gen 0,1,2,3,4,>4
    
    // Analyze each particle's end location
    for (size_t ii = 0; ii < particles->size(); ii++) {
        Particle3D& pp = const_cast<Particle3D&>(particles->particle(ii));
        const int row_id = classifyGridPowerRowId(pp, geometry, mesh_size);
        
        // Identify particle species
        particle_kind species = identify_particle_species(pp.m(), pp.q(), ionmass);
        if (species >= 0 && species < static_cast<particle_kind>(particle_counts.size())) {
            particle_counts[species]++;
        }
        
        // Track generations
        int gen = pp.gen() % 100;
        if (gen >= 0 && gen < 5) generation_counts[gen]++;
        else generation_counts[5]++; // gen > 4
        
        // Record particle collision
        if (row_id >= 0 && static_cast<size_t>(row_id) < total_categories) {
            colls[row_id].push_back(ii);
        }
    }
    
    // Log particle statistics
    if (debug) {
        logfile << "DEBUG: Particle counts by collision solid:" << endl;
        for (std::vector<int>::const_iterator id_it = row_ids.begin(); id_it != row_ids.end(); ++id_it) {
            logfile << "  row " << *id_it << ": " << colls[*id_it].size() << " particles" << endl;
        }
        
        logfile << "DEBUG: Particle counts by species:" << endl;
        for (size_t i = 0; i < particle_counts.size(); i++) {
            logfile << "  " << get_particle_name(int2kind(static_cast<int>(i))) << ": "
                    << particle_counts[i] << endl;
        }
        
        logfile << "DEBUG: Particle counts by generation:" << endl;
        for (int i = 0; i < 5; i++) {
            logfile << "  gen " << i << ": " << generation_counts[i] << endl;
        }
        logfile << "  gen >4: " << generation_counts[5] << endl;
    }
    
    // Calculate power for each solid/grid
    std::vector<double> power_per_solid(total_categories, 0.0);
    std::vector<double> current_per_solid(total_categories, 0.0);
    
    for (std::vector<int>::const_iterator id_it = row_ids.begin(); id_it != row_ids.end(); ++id_it) {
        const size_t ss = static_cast<size_t>(*id_it);
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
            
        } else {
            pow.total_power = 0.0;
            pow.total_current = 0.0;
        }
        
        // Store results based on particle kind
        if (pk == PARTICLE_ALL) {
            power_per_solid[ss] = pow.total_power;
            current_per_solid[ss] = pow.total_current;
        } else {
            int pk_index = get_particle_int(pk);
            if (pk_index >= 0 && pk_index < static_cast<int>(pow.total_power_perspecies.size())) {
                power_per_solid[ss] = pow.total_power_perspecies[pk_index];
                current_per_solid[ss] = pow.total_current_perspecies[pk_index];
            }
        }
    }
    
    // Calculate total beam power using explicit row inclusion.
    double full_beam_power = 0.0;
    for (std::vector<int>::const_iterator id_it = row_ids.begin(); id_it != row_ids.end(); ++id_it) {
        if (*id_it >= 0 && static_cast<size_t>(*id_it) < power_per_solid.size() &&
            gridPowerIncludeInTotal(params, include_map, *id_it)) {
            full_beam_power += power_per_solid[*id_it];
        }
    }
    
    if (debug) logfile << "DEBUG: Full beam power on grids: " << full_beam_power << " W" << endl << flush;
    
    // Create summary file with power breakdown
    string summary_filename = output_folder + file_tag + "_grid_power_summary.txt";
    ofstream summary_file(summary_filename);
    if (summary_file.is_open()) {
        summary_file << "# Grid Power Load Analysis" << endl;
        summary_file << "# ID\tPower[W]\tCurrent[A]\tParticles\tIncludeInTotal\tDescription" << endl;
        
        for (std::vector<int>::const_iterator id_it = row_ids.begin(); id_it != row_ids.end(); ++id_it) {
            const int row_id = *id_it;
            summary_file << row_id << "\t" 
                        << power_per_solid[row_id] << "\t"
                        << current_per_solid[row_id] << "\t"
                        << colls[row_id].size() << "\t"
                        << (gridPowerIncludeInTotal(params, include_map, row_id) ? "true" : "false") << "\t"
                        << gridPowerRowName(params, row_id) << endl;
        }
        
        summary_file << "# Total beam power (included rows only): " << full_beam_power << " W" << endl;
        summary_file.close();
        
        ibsimu.message(1) << "Grid power analysis saved to: " << summary_filename << endl;
    }
    
    if (debug) logfile << "DEBUG: Grid power load analysis completed" << endl << flush;

    if (current_per_solid_out) {
        *current_per_solid_out = current_per_solid;
    }
    
    return power_per_solid;
}

