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
public:
    /**
     * @brief Constructor
     */
    DiagnosticsManager() = default;

    /**
     * @brief Destructor
     */
    ~DiagnosticsManager() = default;

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
    void createPlots(int argc, char **argv, const Geometry* geometry,
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
    void createPlots(int argc, char **argv, const Geometry* geometry,
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
     * @brief Create simulation summary
     * @param zlocsummary Z location for summary
     * @param summary_file_tag File tag for summary
     * @param append_at_end Whether to append to existing file
     * @param particles Particle database
     * @param outsummary_fold Output summary folder
     */
    void createSimulationSummary(double zlocsummary, const std::string& summary_file_tag, 
                                bool append_at_end, const ParticleDataBase3D* particles,
                                const std::string& outsummary_fold);
    
    /**
     * @brief Create simulation summary for specific particle species
     * @param zlocsummary Z location for summary
     * @param summary_file_tag File tag for summary
     * @param append_at_end Whether to append to existing file
     * @param particles_species Species-specific particle databases
     * @param pk Particle kind
     * @param outsummary_fold Output summary folder
     */
    void createSimulationSummary(double zlocsummary, const std::string& summary_file_tag, 
                                bool append_at_end, 
                                const std::vector<ParticleDataBase3D*>& particles_species,
                                particle_kind pk, const std::string& outsummary_fold);

    /**
     * @brief Create simulation summary with field calculations
     * @param zlocsummary Z location for summary
     * @param summary_file_tag File tag for summary
     * @param append_at_end Whether to append to existing file
     * @param particles Particle database
     * @param outsummary_fold Output summary folder
     * @param geometry Simulation geometry
     * @param potential Electric potential field
     * @param magnetic Magnetic field
     */
    void createSimulationSummary(double zlocsummary, const std::string& summary_file_tag, 
                                bool append_at_end, const ParticleDataBase3D* particles,
                                const std::string& outsummary_fold,
                                const Geometry* geometry, const EpotField* potential, 
                                const MeshVectorField* magnetic);

    /**
     * @brief Create simulation summary with field calculations and extracted current density
     * @param zlocsummary Z location for summary analysis
     * @param summary_file_tag File tag for summary outputs
     * @param append_at_end Whether to append to existing files
     * @param particles Particle database
     * @param outsummary_fold Output folder for summary
     * @param geometry Geometry object
     * @param potential Electric potential field
     * @param magnetic Magnetic field
     * @param extracted_current_density Extracted current density at EG exit (A/m²)
     */
    void createSimulationSummary(double zlocsummary, const std::string& summary_file_tag, 
                                bool append_at_end, const ParticleDataBase3D* particles,
                                const std::string& outsummary_fold,
                                const Geometry* geometry, const EpotField* potential, 
                                const MeshVectorField* magnetic, double extracted_current_density);

    /**
     * @brief Analyze particle end locations and calculate grid power loads
     * @param particles Particle database to analyze
     * @param zgrids Vector of grid z-positions
     * @param mesh_size Mesh size for grid collision detection
     * @param geometry Geometry object for domain boundaries
     * @param ionmass Ion mass for power calculations
     * @param output_folder Output folder for grid power files
     * @param file_tag File tag for output naming
     * @param pk Particle kind (PARTICLE_ALL for all species)
     * @return Vector of power values per grid/solid
     */
    std::vector<double> analyzeGridPowerLoads(const ParticleDataBase3D* particles,
                                             const std::vector<double>& zgrids,
                                             double mesh_size,
                                             const Geometry* geometry,
                                             double ionmass,
                                             const std::string& output_folder,
                                             const std::string& file_tag,
                                             particle_kind pk = PARTICLE_ALL,
                                             std::vector<double>* current_per_solid_out = nullptr);

    /**
     * @brief Create enhanced simulation summary with current density analysis
     * @param accelerated_current_density Current density at accelerator exit (A/m²)
     * @param extracted_current_density Current density at EG exit (A/m²)
     * @param filename Output filename for summary
     * @param append_at_end Whether to append to existing file
     * @param particles Particle database for analysis
     * @param extra_info Additional information to include in summary
     */
    void createSimulationSummary(double accelerated_current_density,
                                double extracted_current_density,
                                const std::string& filename, bool append_at_end,
                                const ParticleDataBase3D* particles,
                                const std::string& extra_info);

    /**
     * @brief Create individual simulation summary (for single simulation in scan)
     * @param scan_index Index of simulation in scan (0-based)
     * @param simulation_tag Tag identifying the simulation (e.g. "MTF_STRIP_TEST_0")
     * @param zlocsummary Z location for summary analysis
     * @param particles Particle database
     * @param geometry Geometry object
     * @param potential Electric potential field
     * @param magnetic Magnetic field
     * @param extracted_current_density Extracted current density at EG exit (A/m²)
     * @param outsummary_fold Output folder for individual simulation summary
     * @param pk Particle kind (default: PARTICLE_ALL)
     */
    void createIndividualSimulationSummary(int scan_index, double zlocsummary,
                                          const ParticleDataBase3D* particles,
                                          double extracted_current_density,
                                          const FileManager* fileManager,
                                          particle_kind pk = PARTICLE_ALL);

    /**
     * @brief Add entry to scan-level beam properties summary
     * @param scan_index Index of simulation in scan (0-based)
     * @param simulation_tag Tag identifying the simulation
     * @param scan_folder Main scan folder path
     * @param scan_file_tag Base scan file tag (e.g. "MTF_STRIP_TEST")
     * @param zlocsummary Z location for analysis
     * @param particles Particle database
     * @param extracted_current_density Extracted current density at EG exit (A/m²)
     * @param magnetic Magnetic field for field analysis
     * @param potential Electric potential field
     * @param pk Particle kind (default: PARTICLE_ALL)
     */
    void addToScanBeamPropertiesSummary(int scan_index, const std::string& simulation_tag,
                                       const std::string& scan_folder, const std::string& scan_file_tag,
                                       double zlocsummary, const ParticleDataBase3D* particles,
                                       double extracted_current_density, const MeshVectorField* magnetic,
                                       const EpotField* potential, particle_kind pk = PARTICLE_ALL);

    /**
     * @brief Add entry to scan-level grid power summary
     * @param scan_index Index of simulation in scan (0-based)
     * @param simulation_tag Tag identifying the simulation
     * @param scan_folder Main scan folder path
     * @param scan_file_tag Base scan file tag (e.g. "MTF_STRIP_TEST")
     * @param particles Particle database
     * @param geometry Geometry object
     * @param pk Particle kind (default: PARTICLE_ALL)
     */
    void addToScanGridPowerSummary(int scan_index, const std::string& simulation_tag,
                                  const std::string& scan_folder, const std::string& scan_file_tag,
                                  const ParticleDataBase3D* particles, const Geometry* geometry,
                                  particle_kind pk = PARTICLE_ALL);

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

    /**
     * @brief Get particle species name string
     * @param pk Particle kind
     * @return String representation of particle species
     */
    std::string getParticleSpeciesName(particle_kind pk) const;
};

#endif /* DIAGNOSTICSMANAGER_H_ */
