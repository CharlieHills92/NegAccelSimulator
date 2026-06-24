/*
 * ParticleManager.h
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored from ManageSimulation)
 */

#ifndef PARTICLEMANAGER_H_
#define PARTICLEMANAGER_H_

#include <string>
#include <vector>

// Forward declarations
class ParticleDataBase3D;
class Geometry;
class Random;
class SimulationParameters;

// Project includes
#include "funct.h"

// IBSIMU includes
#include "particledatabase.hpp"
#include "trajectorydiagnostics.hpp"
#include "random.hpp"

/**
 * @class ParticleManager
 * @brief Handles particle database and particle-related operations
 * 
 * This class is responsible for:
 * - Creating and managing particle databases
 * - Handling different particle species
 * - Managing particle emitters
 * - Tracking particle generations and collisions
 */
class ParticleManager {
private:
    ParticleDataBase3D* particles;
    std::vector<ParticleDataBase3D*> particles_species;
    std::vector<particle_kind> particles_kinds;
    
    // Diagnostics and analysis
    std::vector<trajectory_diagnostic_e> diagnostics;
    TrajectoryDiagnosticData tdata;
    std::vector<std::vector<size_t>>* collisions;
    std::vector<particle_kind>* part_kind;
    std::vector<size_t> part_number;
    std::vector<int>* part_gen;
    
    Random* rgn;

public:
    /**
     * @brief Constructor
     */
    ParticleManager();

    /**
     * @brief Destructor
     */
    ~ParticleManager();

    /**
     * @brief Define particle database
     * @param geometry Simulation geometry
     * @param params Simulation parameters
     */
    void defineParticleDatabase(const Geometry& geometry, const SimulationParameters& params);
    
    /**
     * @brief Define particle database with custom emitter
     * @param geometry Simulation geometry
     * @param params Simulation parameters
     * @param emittername Name of emitter configuration
     * @param backstream Include backstreaming particles
     */
    void defineParticleDatabase(const Geometry& geometry, const SimulationParameters& params,
                               const std::string& emittername, bool backstream);

    /**
     * @brief Fill particle databases by species
     * @param params Simulation parameters
     */
    void fillParticleDatabases(const SimulationParameters& params);

    /**
     * @brief Analyze final particle locations
     * @param params Simulation parameters
     */
    void analyzeParticleEndLocations(const SimulationParameters& params);

    /**
     * @brief Save emitter configuration
     * @param emitname Emitter name
     * @param zloc Z location
     * @param output_folder Output folder path
     * @return Filename of saved emitter
     */
    std::string saveEmitter(const std::string& emitname, double zloc, const std::string& output_folder);

    /**
     * @brief Check extracted current density at EG
     * @param oriJ Original current density target
     * @param extsimJ Simulated extracted current density (output)
     * @param params Simulation parameters
     * @return True if within tolerance
     */
    bool checkEGExtractedCurrent(double oriJ, double& extsimJ, const SimulationParameters& params);

    /**
     * @brief Export particle trajectories to VTK format for ParaView visualization
     * @param filename Output VTK filename (without extension)
        * @param ion_mass_u Active ion mass in atomic mass units for family-aware species tagging
     */
        void exportTrajectoriesToVTK(const std::string& filename, double ion_mass_u);

    // Getters
    ParticleDataBase3D* getParticles() { return particles; }
    const std::vector<ParticleDataBase3D*>& getParticleSpecies() const { return particles_species; }
    const std::vector<particle_kind>& getParticleKinds() const { return particles_kinds; }
    const std::vector<trajectory_diagnostic_e>& getDiagnostics() const { return diagnostics; }
    const TrajectoryDiagnosticData& getTrajectoryData() const { return tdata; }
    std::vector<std::vector<size_t>>* getCollisions() const { return collisions; }
    std::vector<particle_kind>* getParticleKind() const { return part_kind; }
    const std::vector<size_t>& getParticleNumber() const { return part_number; }
    std::vector<int>* getParticleGen() const { return part_gen; }
    Random* getRandom() const { return rgn; }

    // Setters
    void setParticles(ParticleDataBase3D* pdb);

    // Cleanup
    void resetParticles();

private:
    /**
     * @brief Initialize diagnostics
     */
    void initializeDiagnostics();
};

#endif /* PARTICLEMANAGER_H_ */
