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
#include "SurfaceEventLedger.h"
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

namespace {

// ---------------------------------------------------------------------------------------
// Per-triangle power density map (Phase 6)
// ---------------------------------------------------------------------------------------

/// One impact or emission, ready to be binned onto the boundary triangulation.
struct SurfaceDeposit {
    Vec3D loc;
    double power_W;   ///< positive for an arrival, negative for an emission
    int row_id;
};

/// Lowest grid-power row id that corresponds to an actual defined solid.
///
/// Row 0 is In_volume (trajectories that ran out of steps or time in vacuum) and rows 1-6
/// are the six simulation-domain boundary planes -- Geometry::inside() returns 1..6 for
/// out-of-mesh coordinates (geometry.cpp:343-354). Neither is a physical surface, and
/// marching cubes only triangulates nodes with boundary number >= 7 (geometry.cpp:1669), so
/// no triangle exists there to bin onto. Deposits on those rows are reported separately
/// rather than counted as binning failures. Same threshold the EAMCC callback uses to
/// decide whether an impact hit a grid at all.
const int FIRST_SOLID_ROW_ID = 7;

/// Accumulated load on one boundary triangle.
struct TriangleLoad {
    double net_power_W;
    double gross_power_W;
    size_t impacts;
    int row_id;          ///< row of the FIRST contributing deposit; see mixed_row_triangles
    bool mixed_rows;

    TriangleLoad()
        : net_power_W(0.0), gross_power_W(0.0), impacts(0), row_id(-1), mixed_rows(false) {}
};

/// Squared distance from point p to triangle (v1,v2,v3).
///
/// Standard barycentric projection with the (v,w) coordinates clamped into the triangle, so
/// the nearest point may lie in the interior, on an edge, or at a vertex. Extends the
/// barycentric block written inline at libibsimu_patched/src/geometry.cpp:2416-2434, which
/// only tests containment and cannot answer "how far outside".
double pointTriangleDistanceSquared(const Vec3D& p, const Vec3D& v1,
                                    const Vec3D& v2, const Vec3D& v3) {
    const Vec3D e1 = v2 - v1;
    const Vec3D e2 = v3 - v1;
    const Vec3D d = p - v1;

    const double a = e1 * e1;
    const double b = e1 * e2;
    const double c = e2 * e2;
    const double dd = e1 * d;
    const double ee = e2 * d;

    double det = a * c - b * b;
    double v;
    double w;
    if (fabs(det) < 1e-30) {
        // Degenerate (zero-area) triangle: fall back to the nearest vertex.
        double best = (p - v1).ssqr();
        best = std::min(best, (p - v2).ssqr());
        best = std::min(best, (p - v3).ssqr());
        return best;
    }

    v = (c * dd - b * ee) / det;
    w = (a * ee - b * dd) / det;

    // Clamp (v, w) into the triangle v >= 0, w >= 0, v + w <= 1.
    if (v < 0.0) v = 0.0;
    if (w < 0.0) w = 0.0;
    if (v + w > 1.0) {
        const double s = v + w;
        v /= s;
        w /= s;
    }

    const Vec3D closest = v1 + v * e1 + w * e2;
    return (p - closest).ssqr();
}

/// Clamp a coordinate to a valid marching-cubes cell index along one axis.
///
/// Valid cell indices run 0 .. size(a)-2, because mc_case(i,j,k) reads nodes i and i+1.
/// Geometry::surface_triangle_ptr has no bounds check at all and surface_trianglec(i,j,k)
/// recomputes the MC case, so an unclamped index reads _smesh out of bounds at the mesh
/// maximum. Same clamping shape as Geometry::surface_inside
/// (libibsimu_patched/src/geometry.cpp:2247-2253), which returns a boundary index there.
int32_t clampedCellIndex(const Geometry& geom, const Vec3D& loc, int axis) {
    int32_t idx = static_cast<int32_t>(floor((loc[axis] - geom.origo(axis)) / geom.h()));
    const int32_t max_index = static_cast<int32_t>(geom.size(axis)) - 2;
    if (idx < 0) idx = 0;
    if (max_index >= 0 && idx > max_index) idx = max_index;
    return idx;
}

/// Nearest boundary triangle to loc, or -1 if none can be found.
///
/// Searches the containing cell first, then the 26-neighbourhood. The widening is not
/// optional: marching cubes only emits triangles for cells that straddle the boundary
/// (mc_case 0 and 255 produce none), and the impact point comes from bracket_surface
/// against the ANALYTIC solid rather than the MC surface, so it can sit up to ~h/2 from
/// any triangle -- and a mask thinner than a cell may produce no triangles at all nearby.
int32_t nearestSurfaceTriangle(const Geometry& geom, const Vec3D& loc) {
    const int32_t ci = clampedCellIndex(geom, loc, 0);
    const int32_t cj = clampedCellIndex(geom, loc, 1);
    const int32_t ck = clampedCellIndex(geom, loc, 2);

    int32_t best = -1;
    double best_dist2 = 0.0;

    // Two rounds: the containing cell alone, then the full 26-neighbourhood. Keeping them
    // separate means the common case costs one cell lookup.
    for (int round = 0; round < 2 && best < 0; round++) {
        const int span = (round == 0) ? 0 : 1;
        for (int di = -span; di <= span; di++) {
            for (int dj = -span; dj <= span; dj++) {
                for (int dk = -span; dk <= span; dk++) {
                    if (round == 1 && di == 0 && dj == 0 && dk == 0) continue;

                    const int32_t i = ci + di;
                    const int32_t j = cj + dj;
                    const int32_t k = ck + dk;
                    if (i < 0 || j < 0 || k < 0) continue;
                    if (i > static_cast<int32_t>(geom.size(0)) - 2) continue;
                    if (j > static_cast<int32_t>(geom.size(1)) - 2) continue;
                    if (k > static_cast<int32_t>(geom.size(2)) - 2) continue;

                    const int32_t count = geom.surface_trianglec(i, j, k);
                    if (count <= 0) continue;
                    const int32_t first = static_cast<int32_t>(geom.surface_triangle_ptr(i, j, k));

                    for (int32_t t = first; t < first + count; t++) {
                        if (t < 0 || t >= static_cast<int32_t>(geom.surface_trianglec())) continue;
                        const VTriangle& tri = geom.surface_triangle(t);
                        const double dist2 = pointTriangleDistanceSquared(
                            loc,
                            geom.surface_vertex(tri[0]),
                            geom.surface_vertex(tri[1]),
                            geom.surface_vertex(tri[2]));
                        if (best < 0 || dist2 < best_dist2) {
                            best = t;
                            best_dist2 = dist2;
                        }
                    }
                }
            }
        }
    }

    return best;
}

/// Area of boundary triangle a, in m^2.
double surfaceTriangleArea(const Geometry& geom, int32_t a) {
    const VTriangle& tri = geom.surface_triangle(a);
    const Vec3D& v1 = geom.surface_vertex(tri[0]);
    const Vec3D& v2 = geom.surface_vertex(tri[1]);
    const Vec3D& v3 = geom.surface_vertex(tri[2]);
    return 0.5 * cross(v2 - v1, v3 - v1).norm2();
}

float swapEndianFloat(float value) {
    union { float f; uint32_t i; } u;
    u.f = value;
    u.i = ((u.i >> 24) & 0x000000FFu) |
          ((u.i >> 8) & 0x0000FF00u) |
          ((u.i << 8) & 0x00FF0000u) |
          ((u.i << 24) & 0xFF000000u);
    return u.f;
}

int32_t swapEndianInt(int32_t value) {
    uint32_t u = static_cast<uint32_t>(value);
    u = ((u >> 24) & 0x000000FFu) |
        ((u >> 8) & 0x0000FF00u) |
        ((u << 8) & 0x00FF0000u) |
        ((u << 24) & 0xFF000000u);
    return static_cast<int32_t>(u);
}

void writeVTKFloatScalars(std::ofstream& file, const std::string& name,
                          const std::vector<float>& values) {
    file << "\nSCALARS " << name << " float 1\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < values.size(); i++) {
        float be = swapEndianFloat(values[i]);
        file.write(reinterpret_cast<const char*>(&be), sizeof(float));
    }
}

void writeVTKIntScalars(std::ofstream& file, const std::string& name,
                        const std::vector<int32_t>& values) {
    file << "\nSCALARS " << name << " int 1\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < values.size(); i++) {
        int32_t be = swapEndianInt(values[i]);
        file.write(reinterpret_cast<const char*>(&be), sizeof(int32_t));
    }
}

/// Human-readable species label for a breakdown key.
///
/// PARTICLE_WRONG (-1) and anything outside the known slots become "unclassified" rather
/// than being dropped, so that the per-species rows sum back to the row total.
std::string breakdownSpeciesName(int kind_int) {
    if (kind_int < 0 || kind_int >= static_cast<int>(particle_kind_count())) {
        return "unclassified";
    }
    return get_particle_name(int2kind(kind_int));
}

/// Where a particle of this RAW generation came from.
///
/// Volume secondaries carry gen+1 and surface secondaries gen+SURFACE_GENERATION_OFFSET
/// (101), so the raw value separates the two; gen 0 is an injected primary. Purely for
/// readability -- Gen is reported raw in its own column alongside this.
const char* breakdownOriginName(int gen) {
    if (gen == 0) return "primary";
    if (gen >= 100) return "surface";
    return "volume";
}

/**
 * @brief Write the per-solid x per-species x per-generation power/current breakdown.
 *
 * A new file with its own parser, so unlike <tag>_grid_power_summary.txt it is free to use
 * real named columns instead of comment lines -- the 5-or-6-tab-field contract that
 * parse_grid_power_summary_txt enforces does not apply here.
 *
 * Every cell is written, including zero-power ones and the "unclassified" species, so that
 * Sum(GrossPower over species and generations) reconciles exactly with the row's Power[W]
 * in the summary file.
 */
void writeGridPowerBreakdown(
        const std::string& output_folder,
        const std::string& file_tag,
        const SimulationParameters& params,
        const std::map<int, bool>& include_map,
        const std::vector<int>& row_ids,
        const std::vector<std::map<PowerBreakdownKey, PowerBreakdownEntry> >& arrivals,
        const std::vector<std::map<PowerBreakdownKey, PowerBreakdownEntry> >& emissions,
        const std::vector<double>& power_per_solid,
        const std::vector<double>& equivalent_current_per_solid,
        bool have_net) {
    const std::string filename = output_folder + file_tag + "_grid_power_breakdown.txt";
    std::ofstream out(filename.c_str());
    if (!out.is_open()) {
        ibsimu.message(1) << "WARNING: could not open " << filename << " for writing" << endl;
        return;
    }

    out << "# Grid Power Load Breakdown by species and generation" << endl;
    out << "# Origin: primary = injected (gen 0), volume = gas-phase secondary (gen 1-99),"
        << " surface = surface secondary (gen >= 101, i.e. gen + 101)." << endl;
    out << "# GrossPower is what ARRIVED and is always >= 0. EmittedPower/EmittedCurrent are"
        << " what the surface emitted, and are <= 0. Net = Gross + Emitted." << endl;
    out << "# EquivalentCurrent is the unsigned current this flux would carry if singly"
        << " charged. It is reported for EVERY species, neutrals included, and is NEVER"
        << " part of a current total -- Current[A] counts real charge only, so neutrals"
        << " contribute exactly 0 there." << endl;
    if (!have_net) {
        out << "# Net accounting unavailable for this file (no surface-emission ledger:"
            << " either surface collisions are disabled or this is a per-species pass)."
            << " Emitted and Net columns are 0 and equal to Gross respectively." << endl;
    }
    out << "# RowID\tSpecies\tGen\tOrigin\tGrossPower[W]\tEmittedPower[W]\tNetPower[W]"
        << "\tGrossCurrent[A]\tEmittedCurrent[A]\tNetCurrent[A]\tEquivalentCurrent[A]"
        << "\tParticles\tIncludeInTotal\tDescription" << endl;

    for (std::vector<int>::const_iterator id_it = row_ids.begin(); id_it != row_ids.end(); ++id_it) {
        const int row_id = *id_it;
        if (row_id < 0 || static_cast<size_t>(row_id) >= arrivals.size()) continue;

        // Union of the two key sets: a surface can emit a species that never arrives back
        // on it, and vice versa, so neither map alone is complete.
        std::set<PowerBreakdownKey> keys;
        for (std::map<PowerBreakdownKey, PowerBreakdownEntry>::const_iterator it =
                 arrivals[row_id].begin(); it != arrivals[row_id].end(); ++it) {
            keys.insert(it->first);
        }
        for (std::map<PowerBreakdownKey, PowerBreakdownEntry>::const_iterator it =
                 emissions[row_id].begin(); it != emissions[row_id].end(); ++it) {
            keys.insert(it->first);
        }

        const bool included = gridPowerIncludeInTotal(params, include_map, row_id);
        const std::string row_name = gridPowerRowName(params, row_id);

        for (std::set<PowerBreakdownKey>::const_iterator kit = keys.begin(); kit != keys.end(); ++kit) {
            PowerBreakdownEntry arr;
            std::map<PowerBreakdownKey, PowerBreakdownEntry>::const_iterator ait =
                arrivals[row_id].find(*kit);
            if (ait != arrivals[row_id].end()) arr = ait->second;

            PowerBreakdownEntry emi;
            std::map<PowerBreakdownKey, PowerBreakdownEntry>::const_iterator eit =
                emissions[row_id].find(*kit);
            if (eit != emissions[row_id].end()) emi = eit->second;

            out << row_id << "\t"
                << breakdownSpeciesName(kit->first) << "\t"
                << kit->second << "\t"
                << breakdownOriginName(kit->second) << "\t"
                << arr.power << "\t" << emi.power << "\t" << (arr.power + emi.power) << "\t"
                << arr.current << "\t" << emi.current << "\t" << (arr.current + emi.current) << "\t"
                << arr.equivalent_current << "\t"
                << arr.count << "\t"
                << (included ? "true" : "false") << "\t"
                << row_name << endl;
        }

        // Reconciliation check against the summary file's Power[W] for this row. Both are
        // built from the same PowerStruct, so a mismatch means a cell was dropped -- which
        // is exactly what the unclassified row exists to prevent.
        double row_gross = 0.0;
        double row_equivalent = 0.0;
        for (std::map<PowerBreakdownKey, PowerBreakdownEntry>::const_iterator it =
                 arrivals[row_id].begin(); it != arrivals[row_id].end(); ++it) {
            row_gross += it->second.power;
            row_equivalent += it->second.equivalent_current;
        }
        const double reference = power_per_solid[row_id];
        if (fabs(row_gross - reference) > 1e-9 * (1.0 + fabs(reference))) {
            ibsimu.message(1) << "WARNING: breakdown rows for grid-power row " << row_id
                              << " sum to " << row_gross << " W but the row total is "
                              << reference << " W." << endl;
        }
        const double eq_reference = equivalent_current_per_solid[row_id];
        if (fabs(row_equivalent - eq_reference) > 1e-9 * (1.0 + fabs(eq_reference))) {
            ibsimu.message(1) << "WARNING: breakdown equivalent current for grid-power row "
                              << row_id << " sums to " << row_equivalent << " A but the row"
                              << " total is " << eq_reference << " A." << endl;
        }
    }

    out.close();
    ibsimu.message(1) << "Grid power breakdown saved to: " << filename << endl;
}

/**
 * @brief Bin surface deposits onto the boundary triangulation and write the density map.
 *
 * Uses the marching-cubes boundary triangulation IBSimu already builds -- build_surface()
 * is called unconditionally in run_simulation (ManageSimulation_New.cpp:171) -- so no
 * additional surface reconstruction is needed.
 *
 * Two conventions worth stating explicitly, because both are easy to get wrong:
 *
 *  * The accumulator is keyed by TRIANGLE INDEX ALONE, never by (triangle, row_id).
 *    Vertices are shared between neighbouring cells (mc_add_vertex reuses from i-1/j-1/k-1,
 *    geometry.cpp:1868-1888) but each triangle is added exactly once, so iterating global
 *    triangle indices cannot double-count. A (triangle, row_id) key would emit one polygon
 *    per row_id and count the same area more than once.
 *  * Triangles are attributed to an electrode through the IMPACTING PARTICLE's row_id, not
 *    by re-deriving the solid from the triangle. surface_inside_solid_number inherits the
 *    documented limitation that the triangulation "does not separate solids which are
 *    touching ... higher than one grid cell" (geometry.hpp:505-515).
 */
void writePowerDensityMap(const std::string& vtk_folder,
                          const std::string& summary_folder,
                          const std::string& file_tag,
                          const SimulationParameters& params,
                          const std::map<int, bool>& include_map,
                          const Geometry* geometry,
                          const std::vector<int>& row_ids,
                          const std::vector<SurfaceDeposit>& deposits) {
    if (!geometry) {
        ibsimu.message(1) << "Power density map skipped: no geometry available" << endl;
        return;
    }
    if (!geometry->surface_built()) {
        ibsimu.message(1) << "Power density map skipped: boundary triangulation not built"
                          << endl;
        return;
    }
    if (deposits.empty()) {
        ibsimu.message(1) << "Power density map skipped: no surface deposits recorded" << endl;
        return;
    }

    const size_t total_categories = static_cast<size_t>(
        row_ids.empty() ? 0 : row_ids.back() + 1);

    std::map<int32_t, TriangleLoad> triangles;
    std::vector<double> unbinned_net_power(total_categories + 1, 0.0);
    std::vector<size_t> unbinned_count(total_categories + 1, 0);
    // Exact per-row net power, summed straight from the deposits. Deliberately NOT derived
    // from the triangle accumulator: because triangles are keyed by index alone, a triangle
    // that receives deposits from two different rows is attributed to just one of them, so
    // triangle-derived per-row power would not reconcile with the summary file. Triangle
    // count, wetted area and max density are necessarily triangle-attributed, and the
    // mixed-row count below says how often that attribution is ambiguous.
    std::map<int, double> row_deposit_net_power;
    size_t mixed_row_triangles = 0;

    // Deposits that are not on a surface at all: the domain boundary planes and In_volume.
    double nonsurface_net_power = 0.0;
    size_t nonsurface_count = 0;

    for (size_t dd = 0; dd < deposits.size(); dd++) {
        const SurfaceDeposit& dep = deposits[dd];

        if (dep.row_id < FIRST_SOLID_ROW_ID) {
            nonsurface_net_power += dep.power_W;
            nonsurface_count++;
            continue;
        }

        row_deposit_net_power[dep.row_id] += dep.power_W;

        const int32_t tri = nearestSurfaceTriangle(*geometry, dep.loc);

        if (tri < 0) {
            // Explicit bucket rather than a silent drop: a large unbinned fraction means
            // the neighbourhood search is failing, and the reader has to be able to see
            // that instead of trusting a map that quietly lost power.
            const size_t slot = (static_cast<size_t>(dep.row_id) < total_categories)
                                    ? static_cast<size_t>(dep.row_id)
                                    : total_categories;
            unbinned_net_power[slot] += dep.power_W;
            unbinned_count[slot]++;
            continue;
        }

        TriangleLoad& load = triangles[tri];
        load.net_power_W += dep.power_W;
        if (dep.power_W > 0.0) {
            load.gross_power_W += dep.power_W;
            load.impacts++;
        }
        if (load.row_id < 0) {
            load.row_id = dep.row_id;
        } else if (load.row_id != dep.row_id && !load.mixed_rows) {
            load.mixed_rows = true;
            mixed_row_triangles++;
        }
    }

    if (triangles.empty() && !row_deposit_net_power.empty()) {
        ibsimu.message(1) << "WARNING: power density map has no binned triangles; every "
                          << "solid deposit fell into the unbinned bucket." << endl;
    }

    // ---- Assemble the VTK payload -----------------------------------------------------
    // Only vertices referenced by a loaded triangle are written, with a remap, so the file
    // stays proportional to the wetted surface rather than to the whole boundary.
    std::vector<int32_t> tri_indices;
    std::map<uint32_t, int32_t> vertex_remap;
    std::vector<Vec3D> points;

    for (std::map<int32_t, TriangleLoad>::const_iterator it = triangles.begin();
         it != triangles.end(); ++it) {
        tri_indices.push_back(it->first);
        const VTriangle& tri = geometry->surface_triangle(it->first);
        for (int vv = 0; vv < 3; vv++) {
            const uint32_t global_v = tri[vv];
            if (vertex_remap.find(global_v) == vertex_remap.end()) {
                vertex_remap[global_v] = static_cast<int32_t>(points.size());
                points.push_back(geometry->surface_vertex(global_v));
            }
        }
    }

    std::vector<float> density_w_cm2(tri_indices.size(), 0.0f);
    std::vector<float> net_power_w(tri_indices.size(), 0.0f);
    std::vector<float> gross_power_w(tri_indices.size(), 0.0f);
    std::vector<float> area_cm2(tri_indices.size(), 0.0f);
    std::vector<int32_t> impacts(tri_indices.size(), 0);
    std::vector<int32_t> row_of_triangle(tri_indices.size(), -1);

    // Per-row aggregates. Wetted area only -- see the file header note.
    std::map<int, double> wetted_area_m2;
    std::map<int, double> row_max_density;
    std::map<int, size_t> row_triangles;

    for (size_t tt = 0; tt < tri_indices.size(); tt++) {
        const int32_t tri = tri_indices[tt];
        const TriangleLoad& load = triangles[tri];
        const double area = surfaceTriangleArea(*geometry, tri);
        const double area_cm2_value = area * 1.0e4;

        // A zero-area triangle would make the density infinite. Marching cubes can emit
        // degenerate triangles where the surface passes exactly through a node.
        const double density = (area_cm2_value > 1e-12)
                                   ? (load.net_power_W / area_cm2_value)
                                   : 0.0;

        density_w_cm2[tt] = static_cast<float>(density);
        net_power_w[tt] = static_cast<float>(load.net_power_W);
        gross_power_w[tt] = static_cast<float>(load.gross_power_W);
        area_cm2[tt] = static_cast<float>(area_cm2_value);
        impacts[tt] = static_cast<int32_t>(load.impacts);
        row_of_triangle[tt] = load.row_id;

        if (load.impacts > 0) {
            wetted_area_m2[load.row_id] += area;
        }
        row_triangles[load.row_id] += 1;
        std::map<int, double>::iterator mit = row_max_density.find(load.row_id);
        if (mit == row_max_density.end() || density > mit->second) {
            row_max_density[load.row_id] = density;
        }
    }

    // ---- VTK 3.0 BINARY POLYDATA ------------------------------------------------------
    // Hand-rolled, big-endian, following ParticleManager::exportTrajectoriesToVTK so that
    // no VTK library dependency is introduced.
    const std::string vtk_filename = vtk_folder + file_tag + "_power_density.vtk";
    std::ofstream vtk(vtk_filename.c_str(), std::ios::binary);
    if (!vtk.is_open()) {
        ibsimu.message(1) << "WARNING: could not open " << vtk_filename << " for writing" << endl;
    } else {
        vtk << "# vtk DataFile Version 3.0\n";
        vtk << "NegAccel surface power density map (net W/cm2)\n";
        vtk << "BINARY\n";
        vtk << "DATASET POLYDATA\n";

        vtk << "POINTS " << points.size() << " float\n";
        for (size_t pp = 0; pp < points.size(); pp++) {
            for (int a = 0; a < 3; a++) {
                float be = swapEndianFloat(static_cast<float>(points[pp][a]));
                vtk.write(reinterpret_cast<const char*>(&be), sizeof(float));
            }
        }

        // POLYGONS <count> <total ints>, each polygon prefixed by its vertex count.
        vtk << "\nPOLYGONS " << tri_indices.size() << " " << (4 * tri_indices.size()) << "\n";
        for (size_t tt = 0; tt < tri_indices.size(); tt++) {
            const VTriangle& tri = geometry->surface_triangle(tri_indices[tt]);
            int32_t header = swapEndianInt(3);
            vtk.write(reinterpret_cast<const char*>(&header), sizeof(int32_t));
            for (int vv = 0; vv < 3; vv++) {
                int32_t idx = swapEndianInt(vertex_remap[tri[vv]]);
                vtk.write(reinterpret_cast<const char*>(&idx), sizeof(int32_t));
            }
        }

        vtk << "\nCELL_DATA " << tri_indices.size() << "\n";
        writeVTKFloatScalars(vtk, "power_density_W_cm2", density_w_cm2);
        writeVTKFloatScalars(vtk, "net_power_W", net_power_w);
        writeVTKFloatScalars(vtk, "gross_power_W", gross_power_w);
        writeVTKFloatScalars(vtk, "area_cm2", area_cm2);
        writeVTKIntScalars(vtk, "impacts", impacts);
        writeVTKIntScalars(vtk, "row_id", row_of_triangle);
        vtk << "\n";
        vtk.close();
        ibsimu.message(1) << "Power density map saved to: " << vtk_filename << endl;
    }

    // ---- Text summary, so the result is checkable without ParaView ---------------------
    const std::string txt_filename = summary_folder + file_tag + "_power_density_summary.txt";
    std::ofstream txt(txt_filename.c_str());
    if (!txt.is_open()) {
        ibsimu.message(1) << "WARNING: could not open " << txt_filename << " for writing" << endl;
        return;
    }

    txt << "# Surface power density summary (companion to " << file_tag
        << "_power_density.vtk)" << endl;
    txt << "# Density is NET power per unit area: arrivals minus what the surface emitted "
        << "back." << endl;
    txt << "# WettedArea is the summed area of triangles that received at least one "
        << "impact. It is NOT the electrode area: the boundary triangulation cannot "
        << "reliably attribute a triangle to a solid, so triangles are attributed through "
        << "the impacting particle's row instead, and untouched area is unknown." << endl;
    txt << "# Impacts is written per triangle in the VTK file as well. With a few hundred "
        << "macroparticles on a grid the per-triangle values are Poisson noise, not a "
        << "resolved density -- read the counts before reading the density." << endl;
    txt << "# Only rows >= " << FIRST_SOLID_ROW_ID << " are mapped. Row 0 (In_volume) and "
        << "rows 1-6 (the simulation-domain boundary planes) are not physical surfaces and "
        << "carry no boundary triangulation, so their deposits are reported on the "
        << "NonSurface line below instead of being binned." << endl;
    txt << "# NetPower is summed from the deposits themselves and matches the summary "
        << "file's net power. Triangles, WettedArea and MaxDensity are attributed through "
        << "the row of each triangle's first deposit, which is the only attribution "
        << "available (see the WettedArea note above)." << endl;
    if (mixed_row_triangles > 0) {
        txt << "# NOTE: " << mixed_row_triangles << " triangle(s) received deposits "
            << "attributed to more than one row, so their area is credited to the row of "
            << "the first deposit only." << endl;
    }
    txt << "# RowID\tTriangles\tWettedArea[cm2]\tNetPower[W]\tMeanDensity[W/cm2]"
        << "\tMaxDensity[W/cm2]\tUnbinnedPower[W]\tUnbinnedDeposits\tIncludeInTotal"
        << "\tDescription" << endl;

    double total_row_power = 0.0;
    double total_unbinned = 0.0;
    size_t total_unbinned_count = 0;

    for (std::vector<int>::const_iterator id_it = row_ids.begin(); id_it != row_ids.end(); ++id_it) {
        const int row_id = *id_it;
        if (row_id < FIRST_SOLID_ROW_ID) continue;

        const double wetted = wetted_area_m2.count(row_id) ? wetted_area_m2[row_id] : 0.0;
        const double wetted_cm2 = wetted * 1.0e4;
        const double net_p = row_deposit_net_power.count(row_id)
                                 ? row_deposit_net_power[row_id] : 0.0;
        const double max_d = row_max_density.count(row_id) ? row_max_density[row_id] : 0.0;
        const size_t ntri = row_triangles.count(row_id) ? row_triangles[row_id] : 0;
        const double unb = (static_cast<size_t>(row_id) < unbinned_net_power.size())
                               ? unbinned_net_power[row_id] : 0.0;
        const size_t unb_n = (static_cast<size_t>(row_id) < unbinned_count.size())
                                 ? unbinned_count[row_id] : 0;

        total_row_power += net_p;
        total_unbinned += unb;
        total_unbinned_count += unb_n;

        txt << row_id << "\t" << ntri << "\t" << wetted_cm2 << "\t" << net_p << "\t"
            << (wetted_cm2 > 1e-12 ? net_p / wetted_cm2 : 0.0) << "\t" << max_d << "\t"
            << unb << "\t" << unb_n << "\t"
            << (gridPowerIncludeInTotal(params, include_map, row_id) ? "true" : "false")
            << "\t" << gridPowerRowName(params, row_id) << endl;
    }

    // Deposits whose row_id fell outside the configured rows entirely.
    const size_t overflow = unbinned_net_power.size() - 1;
    if (unbinned_count[overflow] > 0) {
        txt << "# UnattributedUnbinned\t" << unbinned_count[overflow] << "\t"
            << unbinned_net_power[overflow] << endl;
        total_unbinned += unbinned_net_power[overflow];
        total_unbinned_count += unbinned_count[overflow];
    }

    // Conservation statement over the solids only: binned + unbinned must reproduce the
    // net power the summary file reports for rows >= FIRST_SOLID_ROW_ID.
    txt << "# TotalSolidNetPower[W]\t" << total_row_power << endl;
    txt << "# TotalUnbinnedNetPower[W]\t" << total_unbinned << "\t(from "
        << total_unbinned_count << " deposit(s))" << endl;
    txt << "# NonSurfaceNetPower[W]\t" << nonsurface_net_power << "\t(from "
        << nonsurface_count << " deposit(s) on In_volume or a domain boundary plane; "
        << "not part of any electrode load)" << endl;

    txt.close();
    ibsimu.message(1) << "Power density summary saved to: " << txt_filename << endl;

    if (total_unbinned_count > 0) {
        const double denom = fabs(total_row_power);
        const double fraction = (denom > 0.0) ? fabs(total_unbinned) / denom : 1.0;
        if (fraction > 0.05) {
            ibsimu.message(1) << "WARNING: " << (100.0 * fraction) << "% of the solid "
                              << "surface power could not be binned onto a boundary "
                              << "triangle (" << total_unbinned_count << " deposits). The "
                              << "26-neighbourhood fallback is not finding triangles; the "
                              << "density map is incomplete." << endl;
        }
    }
}

} // namespace

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

    // Un-extracted primary negative ions are kept out of colls[] so that the Particles
    // column agrees with the Power and Current columns -- PowerStruct::add applies the same
    // predicate, so previously the count included particles the power/current sums did not.
    // They are tallied here instead and reported, rather than vanishing silently.
    std::vector<size_t> excluded_counts(total_categories, 0);
    std::vector<double> excluded_currents(total_categories, 0.0);
    
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
            if (is_unextracted_primary_negative_ion(pp.m(), pp.q(), pp.velocity()[2],
                                                    pp.gen(), ionmass)) {
                excluded_counts[row_id]++;
                excluded_currents[row_id] += pp.IQ();
            } else if (pk != PARTICLE_ALL && species != pk) {
                // Species-restricted pass: keep only the requested species, so that the
                // Particles column agrees with the Power and Current columns. Previously
                // colls[] held EVERY particle while only total_power_perspecies was
                // filtered, so the two columns described different populations -- H20
                // reported "58626 particles, 0 W" and H0 reported z-min as "50000
                // particles, 0 W". Filtering here makes the set pre-filtered, which is
                // what lets the power lookup below collapse to pow.total_power.
            } else {
                colls[row_id].push_back(ii);
            }
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
    std::vector<double> equivalent_current_per_solid(total_categories, 0.0);
    // Per-row (species, raw generation) breakdown, retained past the loop for the
    // breakdown file below.
    std::vector<std::map<PowerBreakdownKey, PowerBreakdownEntry> >
        breakdown_per_solid(total_categories);

    // Per-impact (location, power) pairs for the density map. Collected from PowerStruct's
    // own per-particle arrays rather than recomputed, so there is exactly one place in the
    // code that turns a particle into a wattage.
    const bool want_density_map = (pk == PARTICLE_ALL) &&
                                  params.getDiagnosticWritePowerDensityMap();
    std::vector<SurfaceDeposit> deposits;

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
            pow.total_equivalent_current = 0.0;
        }

        // colls[] is already restricted to the requested species (see the collection loop
        // above), so the plain totals are the per-species totals -- no second filter, and
        // no risk of the two disagreeing.
        power_per_solid[ss] = pow.total_power;
        current_per_solid[ss] = pow.total_current;
        equivalent_current_per_solid[ss] = pow.total_equivalent_current;
        breakdown_per_solid[ss] = pow.breakdown;

        if (want_density_map) {
            for (size_t jj = 0; jj < pow.wdata.size(); jj++) {
                SurfaceDeposit dep;
                dep.loc = Vec3D(pow.xdata[jj], pow.ydata[jj], pow.zdata[jj]);
                dep.power_W = pow.wdata[jj];
                dep.row_id = static_cast<int>(ss);
                deposits.push_back(dep);
            }
        }
    }
    
    // Net accounting: subtract what each surface EMITTED from what arrived at it.
    //
    // Gated on pk == PARTICLE_ALL. The per-species passes see a subset of arrivals, so
    // applying the full emission ledger to them would subtract all emissions from part of
    // the arrivals and drive the result meaninglessly negative. Per-species net belongs in
    // the breakdown table, where the emitted species is resolved.
    std::vector<double> emitted_power_per_solid(total_categories, 0.0);
    std::vector<double> emitted_current_per_solid(total_categories, 0.0);
    // Emission cells keyed the same way as the arrival breakdown, so the two can be joined
    // per (row, species, generation) in the breakdown file. Here the species/generation are
    // those of the EMITTED particle and the row is the surface it left.
    std::vector<std::map<PowerBreakdownKey, PowerBreakdownEntry> >
        emitted_breakdown_per_solid(total_categories);
    bool have_net = false;

    if (pk == PARTICLE_ALL && surface_event_ledger_ && !surface_event_ledger_->empty()) {
        have_net = true;
        const std::vector<SurfaceEvent>& events = surface_event_ledger_->events();
        size_t dropped = 0;
        for (size_t ee = 0; ee < events.size(); ee++) {
            const size_t ss = static_cast<size_t>(events[ee].solid);
            if (ss < total_categories) {
                // Events are already signed negative by the callback, so these accumulate
                // to <= 0 and are ADDED to the gross value to obtain net.
                emitted_power_per_solid[ss] += events[ee].power_W;
                emitted_current_per_solid[ss] += events[ee].current_A;

                PowerBreakdownEntry& cell = emitted_breakdown_per_solid[ss][
                    PowerBreakdownKey(events[ee].kind, events[ee].gen)];
                cell.power += events[ee].power_W;
                cell.current += events[ee].current_A;
                cell.count += 1;

                // Emissions are binned at their OWN location, which is what makes the
                // density map a net map rather than a gross one: a hot spot that
                // backscatters strongly shows a correspondingly reduced local load.
                if (want_density_map) {
                    SurfaceDeposit dep;
                    dep.loc = events[ee].loc;
                    dep.power_W = events[ee].power_W;   // already negative
                    dep.row_id = static_cast<int>(ss);
                    deposits.push_back(dep);
                }
            } else {
                dropped++;
            }
        }
        if (dropped > 0) {
            ibsimu.message(1) << "WARNING: " << dropped << " surface emission event(s) "
                              << "referenced a solid outside the grid-power row range and "
                              << "were not netted out." << endl;
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

        // Excluded tally. Emitted as comment lines so the 6-field data-row contract stays
        // intact for parse_grid_power_summary_txt(), which rejects any row that is not 5 or
        // 6 tab-separated fields. Rows above count only particles that contribute to Power
        // and Current; these are the un-extracted primary negative ions left out of both.
        size_t total_excluded_count = 0;
        double total_excluded_current = 0.0;
        for (std::vector<int>::const_iterator id_it = row_ids.begin(); id_it != row_ids.end(); ++id_it) {
            total_excluded_count += excluded_counts[*id_it];
            total_excluded_current += excluded_currents[*id_it];
        }

        summary_file << "# Excluded from the rows above: un-extracted primary negative ions "
                     << "(generation 0, moving upstream)" << endl;
        for (std::vector<int>::const_iterator id_it = row_ids.begin(); id_it != row_ids.end(); ++id_it) {
            if (excluded_counts[*id_it] > 0) {
                summary_file << "# ExcludedRow\t" << *id_it << "\t"
                             << excluded_counts[*id_it] << "\t"
                             << excluded_currents[*id_it] << endl;
            }
        }
        summary_file << "# TotalExcluded\t" << total_excluded_count << "\t"
                     << total_excluded_current << endl;

        // Net accounting. Same reason as above for using comment lines: Power[W] and
        // Current[A] keep their GROSS meaning (the GUI plots them, and gross is the right
        // quantity for a thermal-margin check), so net is reported alongside rather than
        // in place of them. parse_grid_power_summary_txt merges these back into the
        // matching row dicts by ID, so Python sees them as first-class per-row values.
        //
        // Sign convention: Emitted* are <= 0, and Net = Gross + Emitted.
        if (have_net) {
            double total_emitted_power = 0.0;
            double total_net_power = 0.0;
            double total_emitted_current = 0.0;
            double total_net_current = 0.0;

            summary_file << "# Net accounting: energy and charge carried back off each "
                         << "surface by the secondaries it emitted. Net = Gross + Emitted, "
                         << "Emitted <= 0." << endl;
            for (std::vector<int>::const_iterator id_it = row_ids.begin(); id_it != row_ids.end(); ++id_it) {
                const int row_id = *id_it;
                const double emitted_p = emitted_power_per_solid[row_id];
                const double emitted_i = emitted_current_per_solid[row_id];
                const double net_p = power_per_solid[row_id] + emitted_p;
                const double net_i = current_per_solid[row_id] + emitted_i;

                summary_file << "# NetRow\t" << row_id << "\t"
                             << emitted_p << "\t" << net_p << "\t"
                             << emitted_i << "\t" << net_i << endl;

                if (gridPowerIncludeInTotal(params, include_map, row_id)) {
                    total_emitted_power += emitted_p;
                    total_net_power += net_p;
                    total_emitted_current += emitted_i;
                    total_net_current += net_i;
                }

                // Invariant: net POWER cannot be negative. Emission is budget-capped in
                // THCallback_surf_EAMCC -- Eout is clamped to the running Elost and each
                // secondary electron is gated on Elost >= SECONDARY_ELECTRON_ENERGY -- so
                // the energy leaving a surface can never exceed the energy that arrived.
                // A negative value means the ledger and the arrival post-pass disagree,
                // most likely because an emission was attributed to a different solid than
                // the impact that produced it. Reported, not asserted, so that a
                // diagnostics discrepancy never aborts a completed simulation run.
                //
                // No analogous invariant exists for net CURRENT: one ion can liberate 2.7
                // electrons, so the charge balance legitimately changes sign.
                if (net_p < 0.0 && fabs(net_p) > 1e-12 * (1.0 + fabs(power_per_solid[row_id]))) {
                    ibsimu.message(1) << "WARNING: net power on row " << row_id << " is "
                                      << net_p << " W (gross " << power_per_solid[row_id]
                                      << " W, emitted " << emitted_p << " W). Emission "
                                      << "exceeds arrival, which the EAMCC energy budget "
                                      << "should forbid." << endl;
                }
            }
            summary_file << "# TotalNet\t" << total_emitted_power << "\t" << total_net_power
                         << "\t" << total_emitted_current << "\t" << total_net_current << endl;
        }

        summary_file.close();

        ibsimu.message(1) << "Grid power analysis saved to: " << summary_filename << endl;
    }

    writeGridPowerBreakdown(output_folder, file_tag, params, include_map, row_ids,
                            breakdown_per_solid, emitted_breakdown_per_solid,
                            power_per_solid, equivalent_current_per_solid, have_net);

    if (want_density_map) {
        writePowerDensityMap(vtk_folder_.empty() ? output_folder : vtk_folder_,
                             output_folder, file_tag, params, include_map, geometry,
                             row_ids, deposits);
    }


    if (debug) logfile << "DEBUG: Grid power load analysis completed" << endl << flush;

    if (current_per_solid_out) {
        *current_per_solid_out = current_per_solid;
    }
    
    return power_per_solid;
}

