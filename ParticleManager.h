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

/**
 * @brief One beamlet of a patterned extraction, as the runtime sees it.
 *
 * Built from the particleSources of the runtime case: the multi-beamlet expansion stamps
 * every source it generates with the beamlet it feeds, and this is that stamp plus the
 * source centre, kept so that particles created during tracing (stripping secondaries,
 * surface products) can still be attributed to the beamlet they were born in.
 */
struct BeamletTag {
    int index;      // 1-based, in the pattern's own order; 1 for an unpatterned case
    int column;     // signed lattice node, centred on the pattern origin
    int row;
    double centerX; // m, the emitting source centre
    double centerY;

    BeamletTag() : index(1), column(0), row(0), centerX(0.0), centerY(0.0) {}
};

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

    // Beamlet bookkeeping. beamlet_tags holds one entry per distinct beamlet, in the order
    // the sources appear; particle_beamlet_slot maps a particle index in the database to a
    // slot in it, and is filled as the emitters run, so the attribution is exact for every
    // particle the case emits rather than inferred from where it happens to be.
    std::vector<BeamletTag> beamlet_tags;
    std::vector<int> particle_beamlet_slot;

    /** @brief Find, or create, the beamlet_tags slot for a beamlet index. */
    size_t beamletSlotFor(int beamlet_index, int column, int row, double center_x, double center_y);

    /** @brief Record that particles [first, last) in the database came from a beamlet slot. */
    void tagParticleRange(size_t first, size_t last, size_t slot);

public:
    /**
     * @brief Beamlet a particle belongs to, or NULL when the case has no beamlet tags.
     * @param index Particle index in the database
     * @param birth_x Transverse birth position, used only for particles that were not
     *                emitted by a source (secondaries created during tracing), which are
     *                attributed to the nearest beamlet axis.
     * @param birth_y As birth_x
     */
    const BeamletTag* beamletOfParticle(size_t index, double birth_x, double birth_y) const;

    /** @brief Number of distinct beamlets the current database was emitted from. */
    size_t beamletTagCount() const { return beamlet_tags.size(); }

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
     * @brief Append the beamlet halo to a database that already holds the core beam.
     *
     * Models the halo as in Pimazzoni et al., Fus. Eng. Des. 192 (2023) 113621, Sec. 2:
     * particles are born on the downstream chamfer of the plasma grid, aimed along the local
     * surface normal, and carry a fixed fraction of the core current.
     *
     * Does nothing unless beamHalo.enabled is set, so a case that does not ask for a halo is
     * byte-for-byte unaffected.
     *
     * One halo is emitted per particle source, around that source's own axis and carrying
     * that source's own share of the current, so multi-beamlet extraction gives every
     * beamlet its own halo rather than piling them all onto the first one.
     *
     * @param pdb            database already populated with the core sources
     * @param params         run configuration
     * @param core_current_A summed core beam current [A], reported in the run log
     * @return number of halo macroparticles added, 0 if disabled or not applicable
     */
    size_t addBeamHalo(ParticleDataBase3D* pdb, const SimulationParameters& params,
                       double core_current_A);

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
     *
     * Point data: velocity_magnitude_km_s, time_ns, kinetic_energy_eV.
     * Cell data (one value per trajectory): particle_id, current_A, charge_e, mass_amu,
     * particle_status, particle_generation, particle_kind, and -- when the case was
     * emitted from a beamlet pattern -- beamlet_id (1-based), beamlet_column, beamlet_row
     * and beamlet_center_x_mm / beamlet_center_y_mm. Colour or Threshold on beamlet_id in
     * ParaView to see or isolate the beamlets.
     *
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
