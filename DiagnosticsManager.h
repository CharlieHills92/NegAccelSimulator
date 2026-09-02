/*
 * DiagnosticsManager.h
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored from ManageSimulation)
 */

#ifndef DIAGNOSTICSMANAGER_H_
#define DIAGNOSTICSMANAGER_H_

#include <string>
#include <vector>

// Project includes
#include "funct.h"

// Forward declarations
class ParticleDataBase3D;
class Geometry;
class EpotField;
class EpotEfield;
class MeshScalarField;
class SimulationParameters;
class FileManager;
class SurfaceEventLedger;
struct PowerStruct;

// IBSIMU includes
#include "particledatabase.hpp"
#include "meshvectorfield.hpp"

/**
 * @class DiagnosticsManager
 * @brief Handles analysis, plotting, and diagnostic output
 * 
 * This class is responsible for:
 * - Generating diagnostic data along z-axis
 * - Creating visualization plots
 * - Performing simulation analysis
 * - Creating simulation summaries
 * - Outputting trajectory data
 */
class DiagnosticsManager {
private:
    std::string density_profile_filename_;

    // Signed surface energy/charge balance for the tracking pass being analysed, or null
    // when surface collisions were disabled (in which case there is nothing to net out and
    // the gross columns are already correct). Not owned -- ManageSimulation owns it.
    const SurfaceEventLedger* surface_event_ledger_ = nullptr;

    // Destination folder for the power-density VTK map, injected rather than threaded
    // through analyzeGridPowerLoads (which already takes nine parameters). Empty means the
    // map falls back to the summary folder.
    std::string vtk_folder_;

public:
    /**
     * @brief Constructor
     */
    DiagnosticsManager() = default;

    /**
     * @brief Destructor
     */
    ~DiagnosticsManager() = default;

    void setDensityProfileFilename(const std::string& density_profile_filename) {
        density_profile_filename_ = density_profile_filename;
    }

    /**
     * @brief Attach the signed surface energy/charge ledger for net power/current output.
     *
     * Injected the same way as density_profile_filename_ so that analyzeGridPowerLoads and
     * performAnalysis keep their (already long) signatures. Pass nullptr to detach.
     */
    void setSurfaceEventLedger(const SurfaceEventLedger* ledger) {
        surface_event_ledger_ = ledger;
    }

    /// Where to write the power-density VTK map. Empty falls back to the summary folder.
    void setVTKFolder(const std::string& vtk_folder) {
        vtk_folder_ = vtk_folder;
    }

    /**
     * @brief Generate diagnostic data along z-axis
     * @param diagzpos Vector of z positions for diagnostics
     * @param iteration Current iteration number
     * @param outfile Output filename
     * @param append_at_end Whether to append to existing file
     * @param pdb Particle database to analyze
     */
    void generateDiagnosticData(const std::vector<double>& diagzpos, int iteration, 
                               const std::string& outfile, bool append_at_end, 
                               const ParticleDataBase3D* pdb);
    
    /**
     * @brief Generate diagnostic data with defaults
     * @param diagzpos Vector of z positions for diagnostics
     * @param iteration Current iteration number
     * @param outfile Output filename
     * @param particles Default particle database
     */
    void generateDiagnosticData(const std::vector<double>& diagzpos, int iteration, 
                               const std::string& outfile, const ParticleDataBase3D* particles);

    /**
     * @brief Generate diagnostic data with field calculations
     * @param diagzpos Vector of z positions for diagnostics
     * @param iteration Current iteration number
     * @param outfile Output filename
     * @param append_at_end Whether to append to existing file
     * @param pdb Particle database to analyze
     * @param geometry Simulation geometry
     * @param potential Electric potential field
     * @param magnetic Magnetic field
     */
    void generateDiagnosticData(const std::vector<double>& diagzpos, int iteration, 
                               const std::string& outfile, bool append_at_end, 
                               const ParticleDataBase3D* pdb,
                               const Geometry* geometry,
                               const EpotField* potential,
                               const MeshVectorField* magnetic);

    /**
     * @brief Create visualization plots
     * @param argc Command line argument count
     * @param argv Command line arguments
     * @param geometry Simulation geometry
     * @param potential Electric potential field
     * @param magnetic Magnetic field
     * @param electric Electric field
     * @param spacecharge Space charge field
     * @param particles Particle database
     * @param plot_folder Output folder for plots
     * @param file_tag File tag for naming
     */
    void createPlots(int argc, char **argv, const SimulationParameters& params,
                    const Geometry* geometry,
                    const EpotField* potential, const MeshVectorField* magnetic,
                    const EpotEfield* electric, const MeshScalarField* spacecharge,
                    const ParticleDataBase3D* particles, const std::string& plot_folder,
                    const std::string& file_tag);
    
    /**
     * @brief Create plots for specific particle species
     * @param argc Command line argument count
     * @param argv Command line arguments
     * @param geometry Simulation geometry
     * @param potential Electric potential field
     * @param magnetic Magnetic field
     * @param electric Electric field
     * @param spacecharge Space charge field
     * @param particles Particle database
     * @param particles_species Species-specific particle databases
     * @param pk Particle kind to plot
     * @param plot_folder Output folder for plots
     * @param file_tag File tag for naming
     */
    void createPlots(int argc, char **argv, const SimulationParameters& params,
                    const Geometry* geometry,
                    const EpotField* potential, const MeshVectorField* magnetic,
                    const EpotEfield* electric, const MeshScalarField* spacecharge,
                    const ParticleDataBase3D* particles,
                    const std::vector<ParticleDataBase3D*>& particles_species,
                    particle_kind pk, const std::string& plot_folder,
                    const std::string& file_tag);

    /**
     * @brief Perform comprehensive simulation analysis with field calculations
     * @param particles Particle database
     * @param params Simulation parameters
     * @param geometry Geometry object for field evaluation
     * @param potential Potential field for energy calculations
     * @param magnetic Magnetic field for B-field calculations
     * @param fileManager File manager for path generation
     */
    void performAnalysis(const ParticleDataBase3D* particles, const SimulationParameters& params,
                        const Geometry* geometry, const EpotField* potential, const MeshVectorField* magnetic,
                        const FileManager* fileManager, vector<double> zgrids);

    /**
     * @brief Perform comprehensive simulation analysis with species support
     * @param particles Main particle database
     * @param params Simulation parameters
     * @param geometry Geometry object for field evaluation
     * @param potential Potential field for energy calculations
     * @param magnetic Magnetic field for B-field calculations
     * @param particles_species Vector of species-specific particle databases
     * @param include_stripping Whether to include stripping analysis
     * @param fileManager File manager for path generation
     */
    void performAnalysis(const ParticleDataBase3D* particles, const SimulationParameters& params,
                        const Geometry* geometry, const EpotField* potential, const MeshVectorField* magnetic,
                        const double zlocsummary,const std::vector<ParticleDataBase3D*>& particles_species,
                        bool include_stripping, const FileManager* fileManager, vector<double> zgrids);

    /**
     * @brief Print trajectory data to text file
     * @param filename Output filename
     * @param particles Particle database
     */
    void printTrajectoryData(const std::string& filename, const ParticleDataBase3D* particles);

    /**
     * @brief Analyze particle end locations and calculate grid power loads
     * @param particles Particle database to analyze
     * @param mesh_size Mesh size for grid collision detection
     * @param geometry Geometry object for domain boundaries
     * @param ionmass Ion mass for power calculations
     * @param output_folder Output folder for grid power files
     * @param file_tag File tag for output naming
     * @param pk Particle kind (PARTICLE_ALL for all species)
     * @return Vector of power values per grid/solid
     */
    std::vector<double> analyzeGridPowerLoads(const ParticleDataBase3D* particles,
                                             const SimulationParameters& params,
                                             double mesh_size,
                                             const Geometry* geometry,
                                             double ionmass,
                                             const std::string& output_folder,
                                             const std::string& file_tag,
                                             particle_kind pk = PARTICLE_ALL,
                                             std::vector<double>* current_per_solid_out = nullptr);

private:
    /**
     * @brief Helper method to setup geometry plotter
     * @param geometry Simulation geometry
     * @param potential Electric potential field
     * @param magnetic Magnetic field
     * @param electric Electric field
     * @param spacecharge Space charge field
     * @param particles Particle database
     */
    void setupGeometryPlotter(const Geometry* geometry, const EpotField* potential,
                             const MeshVectorField* magnetic, const EpotEfield* electric,
                             const MeshScalarField* spacecharge, const ParticleDataBase3D* particles);

};

#endif /* DIAGNOSTICSMANAGER_H_ */
