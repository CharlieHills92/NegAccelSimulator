/*
 * FieldManager.h
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored from ManageSimulation)
 */

#ifndef FIELDMANAGER_H_
#define FIELDMANAGER_H_

#include <string>
#include <vector>
#include "SimulationParameters.h"

// Forward declarations
class MeshVectorField;
class MeshScalarField;
class EpotField;
class EpotEfield;
class Geometry;

/**
 * @class FieldManager
 * @brief Handles magnetic and electric field management
 * 
 * This class is responsible for:
 * - Loading and combining magnetic fields
 * - Managing electric fields and potentials
 * - Setting field boundary conditions
 */
class FieldManager {
private:
    EpotField* potential;
    MeshScalarField* spacecharge;
    MeshVectorField* magnetic;
    MeshScalarField* gas_dens;
    EpotEfield* electric;

    // Helper methods
    MeshVectorField* loadMagneticFieldDefinition(const SimulationParameters::MagneticFieldDefinition& definition,
                                                 const std::string& base_directory,
                                                 const MeshVectorField* reference_mesh = NULL);

public:
    /**
     * @brief Constructor
     */
    FieldManager();

    /**
     * @brief Destructor
     */
    ~FieldManager();

    /**
     * @brief Add magnetic field to simulation
     * @param params Simulation parameters
    * @param bfield_fold Magnetic field source path
     */
    void addMagneticField(const SimulationParameters& params, const std::string& bfield_fold);

    /**
    * @brief Get the configured magnetic field source path
     * @param params Simulation parameters
    * @return Source directory or file path
     */
    static std::string getBFieldFolder(const SimulationParameters& params);

    // Getters
    EpotField* getPotential() const { return potential; }
    MeshScalarField* getSpacecharge() const { return spacecharge; }
    MeshVectorField* getMagnetic() const { return magnetic; }
    MeshScalarField* getGasDens() const { return gas_dens; }
    EpotEfield* getElectric() const { return electric; }

    // Setters
    void setPotential(EpotField* pot) { potential = pot; }
    void setSpacecharge(MeshScalarField* sc) { spacecharge = sc; }
    void setMagnetic(MeshVectorField* mag) { magnetic = mag; }
    void setGasDens(MeshScalarField* gd) { gas_dens = gd; }
    void setElectric(EpotEfield* elec) { electric = elec; }

    // Cleanup
    void resetFields();
};

#endif /* FIELDMANAGER_H_ */
