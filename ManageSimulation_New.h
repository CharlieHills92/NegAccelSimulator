/*
 * ManageSimulation_New.h
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored ManageSimulation)
 */

#ifndef MANAGESIMULATION_NEW_H_
#define MANAGESIMULATION_NEW_H_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "funct.h"

// Component headers
#include "SimulationParameters.h"
#include "FileManager.h"
#include "GeometryManager.h"
#include "FieldManager.h"
#include "ParticleManager.h"
#include "DiagnosticsManager.h"
#include "SurfaceEventLedger.h"

// Surface-collision (EAMCC) trajectory-end callback. Only the definition in
// ManageSimulation_New.cpp needs the complete type, so a forward declaration
// keeps THCallback_surf_EAMCC.h out of every translation unit that includes us.
class THCallback_surf_EAMCC;

/**
 * @brief One (species, generation) cell of a surface's power/current breakdown.
 */
struct PowerBreakdownEntry {
    double power;                ///< W, always >= 0 (arrival side)
    double current;              ///< A, signed; exactly 0 for neutral species
    double equivalent_current;   ///< A, unsigned "if singly charged" flux; never a current total
    size_t count;                ///< macroparticles

    PowerBreakdownEntry()
        : power(0.0), current(0.0), equivalent_current(0.0), count(0) {}
};

/// Breakdown key: (particle_kind as int, RAW generation).
///
/// The generation is deliberately kept raw rather than reduced modulo 100: surface
/// secondaries carry gen >= SURFACE_GENERATION_OFFSET (101) and volume secondaries 1-4, and
/// telling those apart is the whole point of the breakdown. PARTICLE_WRONG (-1) is a valid
/// key and becomes the "unclassified" row, so that the per-species sums reconcile with the
/// row total instead of quietly falling short.
typedef std::pair<int, int> PowerBreakdownKey;

// PowerStruct (keeping for compatibility)
struct PowerStruct {
    double ionmass;
    size_t solid_idx;
    double total_power;
    double total_current;
    // Unsigned "if singly charged" flux current. Kept separate from total_current so that
    // neutral flux is reportable without ever entering an electrical current total.
    double total_equivalent_current;
    std::vector<double> total_power_perspecies;
    std::vector<double> total_current_perspecies;
    std::vector<double> total_equivalent_current_perspecies;
    // Replaces the former total_power_pergen / total_current_pergen fixed 6-element
    // vectors. Those were written but never read anywhere, and indexing them by the raw
    // generation overflowed the heap once surface secondaries (gen >= 101) existed.
    std::map<PowerBreakdownKey, PowerBreakdownEntry> breakdown;
    std::vector<double> xdata;
    std::vector<double> ydata;
    std::vector<double> zdata;
    std::vector<double> vxdata;
    std::vector<double> vydata;
    std::vector<double> vzdata;
    std::vector<double> curdata;
    std::vector<double> eqcurdata;
    std::vector<double> mdata;
    std::vector<double> wdata;
    std::vector<double> qdata;
    std::vector<int> gendata;
    std::vector<int> kinddata;

    void add(const Particle3D &pp);
    void calculate_total_power();
    void print(const std::string &filename);
};

/**
 * @class ManageSimulation
 * @brief Main orchestrator class for particle beam simulations
 * 
 * This refactored class coordinates between specialized manager components:
 * - SimulationParameters: Input parameter management
 * - FileManager: File operations and directory management  
 * - GeometryManager: Geometry creation for different accelerator types
 * - FieldManager: Magnetic and electric field management
 * - ParticleManager: Particle database and particle operations
 * - DiagnosticsManager: Analysis, plotting, and diagnostic output
 */
class ManageSimulation {
private:
    // Manager components
    std::unique_ptr<SimulationParameters> parameters;
    std::unique_ptr<FileManager> fileManager;
    std::unique_ptr<GeometryManager> geometryManager;
    std::unique_ptr<FieldManager> fieldManager;
    std::unique_ptr<ParticleManager> particleManager;
    std::unique_ptr<DiagnosticsManager> diagnosticsManager;

    // Surface-collision callback. It caches a Geometry& and a ParticleDataBase3D*,
    // both of which are destroyed and recreated by every create_geometry()/define_pdb()
    // pass (the JEXT current-matching loop calls run_simulation up to 5 times). Owning
    // it here and rebuilding it per run_simulation() call keeps those references valid.
    std::unique_ptr<THCallback_surf_EAMCC> surfaceCollisionCallback;

    // Signed surface energy/charge balance filled by the callback above during the final
    // tracking pass, then read by DiagnosticsManager to form the net power and current.
    // Owned here because it must survive iterate_trajectories() and be readable during
    // the analysis that follows, but must not survive into the next run_simulation() pass.
    SurfaceEventLedger surfaceEventLedger;

    // Helper methods
    void initializeIbsimu();
    void initializeComponents(const std::string& scan_name, const std::string& foldername);

public:
    /**
     * @brief Constructor with scan name and folder
     * @param scan_name Scan name identifier
     * @param foldername Output folder name
     */
    ManageSimulation(const std::string& scan_name, const std::string& foldername);
    
    /**
     * @brief Constructor with scan name only (uses same name for folder)
     * @param scan_name Scan name identifier
     */
    ManageSimulation(const std::string& scan_name) : ManageSimulation(scan_name, scan_name) {}
    
    /**
     * @brief Destructor
     */
    ~ManageSimulation();
    
    /**
     * @brief Reset simulation variables and free memory
     */
    void ResetSimulation();

    // Core simulation methods
    /**
     * @brief Create geometry based on accelerator type
     */
    void create_geometry();
    
    /**
     * @brief Create geometry with custom parameters
     * @param z_start Starting z position
     * @param z_end Ending z position
     * @param meshsize_multiplier Mesh size scaling factor
     */
    void create_geometry(double z_start, double z_end, double meshsize_multiplier);

    /**
     * @brief Add magnetic field to simulation
     */
    void add_Bfield();
    
    /**
     * @brief Define particle database
     */
    void define_pdb();
    
    /**
     * @brief Define particle database with custom emitter
     * @param emittername Name of emitter configuration
     * @param backstream Include backstreaming particles
     */
    void define_pdb(const std::string& emittername, bool backstream);

    /**
     * @brief Run the simulation
     */
    void run_simulation();
    
    /**
     * @brief Run simulation with optional particle database cycling
     * @param pdbincycle Whether to recreate particle database each iteration
     */
    void run_simulation(bool pdbincycle);
    
    /**
     * @brief Check extracted current density at EG
     * @param oriJ Original current density target
     * @param extsimJ Simulated extracted current density (output)
     * @return True if within tolerance
     */
    bool check_EGext(double oriJ, double& extsimJ);
    
    /**
     * @brief Load existing simulation data
     * @return True if successful
     */
    bool load_simulation();

    /**
     * @brief Trace particles using loaded fields (potential, magnetic, space charge)
     * @param use_stripping Whether to enable stripping callbacks during tracing
     * @return True if successful
     */
    bool trace_particles_with_loaded_fields(bool use_stripping = true);

    // Analysis and output methods
    /**
     * @brief Generate diagnostic data along z-axis
     */
    void diagnostic_data_alongZ(const std::vector<double>& diagzpos, int iteration, 
                               const std::string& outfile, bool append_at_end);
    void diagnostic_data_alongZ(const std::vector<double>& diagzpos, int iteration, 
                               const std::string& outfile);
    
    /**
     * @brief Create visualization plots
     */
    void plot_simulation(int argc, char **argv);
    void plot_simulation(int argc, char **argv, particle_kind pk);
    
    /**
     * @brief Perform simulation analysis
     */
    void analysis(double zlocsummary);

    /**
     * @brief Fill particle databases by species
     */
    void fill_particle_dbs();
    
    /**
     * @brief Print trajectory data to text file
     */
    void print_traj_to_txt(const std::string& filename);
    
    /**
     * @brief Save emitter configuration
     */
    std::string save_emitter(const std::string& emitname, double zloc);
    
    /**
     * @brief Analyze final particle locations
     */
    void particles_end_location();

    /**
     * @brief Analyze grid power loads from particle collisions
     * @param pk Particle kind to analyze (PARTICLE_ALL for all species)
     * @return Vector of power values per grid/solid
     */
    std::vector<double> analyze_grid_power_loads(particle_kind pk = PARTICLE_ALL);

    /**
     * @brief Export geometry and fields to VTK format for ParaView visualization
     * @param base_filename Base filename for VTK exports (without extension)
     */
    void export_for_paraview(const std::string& base_filename);

    /**
     * @brief Export geometry to VTK format immediately after creation
     * @param base_filename Base filename for VTK exports (without extension)
     */
    void export_geometry_to_vtk(const std::string& base_filename);
    
    /**
     * @brief Export complete simulation results to VTK format for ParaView visualization
     * @param base_filename Base filename for VTK exports (without extension)
     * @param scharge Space-charge field
     */
    void export_simulation_results_to_vtk(const std::string& base_filename,
                                          const MeshScalarField* scharge);
    
    void save_results_to_vtk();

    // Getters (delegating to appropriate managers)
    std::string get_file_tag() const { return fileManager->getFileTag(); }
    void set_file_tag(const std::string& filetag);
    
    uint get_iterations() const { return parameters->getIterations(); }
    void set_iterations(uint itn) { parameters->setIterations(itn); }
    
    double get_M_IONS() const { return parameters->getMIons(); }
    double get_Q_IONS() const { return parameters->getQIons(); }
    double get_J_ION() const { return parameters->getJIon(); }
    void set_J_ION(double j) { parameters->setJIon(j); }
    
    uint get_N_PARTICLES() const { return parameters->getNParticles(); }
    double get_EG_VOLTAGE() const { return parameters->getEGVoltage(); }
    double get_GG_VOLTAGE() const { return parameters->getGGVoltage(); }
    double get_REP_VOLTAGE() const { return parameters->getREPVoltage(); }
    double get_G1_VOLTAGE() const { return parameters->getG1Voltage(); }
    double get_G2_VOLTAGE() const { return parameters->getG2Voltage(); }
    double get_G3_VOLTAGE() const { return parameters->getG3Voltage(); }
    double get_G4_VOLTAGE() const { return parameters->getG4Voltage(); }
    double get_G5_VOLTAGE() const { return parameters->getG5Voltage(); }
    double get_MESH_SIZE() const { return parameters->getMeshSize(); }
    
    uint get_split_domain() const { return parameters->getSplitDomain(); }
    double get_tolerance() const { return parameters->getJTolerance(); }
    double get_ext_J() const { return parameters->getEGExtJ(); }
    void set_ext_J(double setJ) { parameters->setEGExtJ(setJ); }
    uint get_stripping() const { return parameters->getIncludeStripping(); }
    double get_domain_z_size() const { return parameters->getDomainZSizeOrDefault(); }
    double get_domain_z_start() const { return parameters->getDomainZStart(); }
    double get_domain_z_end() const { return get_domain_z_start() + get_domain_z_size(); }

    // Component access (for advanced usage)
    SimulationParameters* getParameters() const { return parameters.get(); }
    FileManager* getFileManager() const { return fileManager.get(); }
    GeometryManager* getGeometryManager() const { return geometryManager.get(); }
    FieldManager* getFieldManager() const { return fieldManager.get(); }
    ParticleManager* getParticleManager() const { return particleManager.get(); }
    DiagnosticsManager* getDiagnosticsManager() const { return diagnosticsManager.get(); }
};

#endif /* MANAGESIMULATION_NEW_H_ */
