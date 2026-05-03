/*
 * GeometryManager.h
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored from ManageSimulation)
 */

#ifndef GEOMETRYMANAGER_H_
#define GEOMETRYMANAGER_H_

#include <string>
#include <vector>

// Forward declarations
class Geometry;
class MITICA_Accelerator;
class SimulationParameters;
class MeshScalarField;
class ParticleDataBase3D;

// IBSIMU includes
#include "geometry.hpp"
#include "func_solid.hpp"

/**
 * @class GeometryManager
 * @brief Handles geometry creation for different accelerator types
 * 
 * This class is responsible for:
 * - Creating geometries for SPIDER, MITICA, and MTF accelerators
 * - Managing accelerator-specific solid objects
 * - Setting boundary conditions
 */
class GeometryManager {
private:
    Geometry* geometry;
    MITICA_Accelerator* accelerator;
    std::vector<double> zgrids;

    // Helper methods for specific accelerator types
    void createSPIDERGeometry(const Int3D& meshsize, const Vec3D& origo, double mesh_size, 
                             const SimulationParameters& params, const std::string& geomfile);

public:
    /**
     * @brief Constructor
     */
    GeometryManager();

    /**
     * @brief Destructor
     */
    ~GeometryManager();

    /**
     * @brief Create geometry based on accelerator type from parameters
     * @param params Simulation parameters
     * @param geomfile Output geometry file path
     */
    void createGeometry(const SimulationParameters& params, const std::string& geomfile);
    
    /**
     * @brief Create geometry with custom parameters
     * @param params Simulation parameters
     * @param geomfile Output geometry file path
     * @param z_start Starting z position
     * @param z_end Ending z position
     * @param meshsize_multiplier Mesh size scaling factor
     */
    void createGeometry(const SimulationParameters& params, const std::string& geomfile,
                       double z_start, double z_end, double meshsize_multiplier);

    /**
     * @brief Create SPIDER geometry
     * @param params Simulation parameters
     * @param geomfile Output geometry file path
     */
    void createGeometrySPIDER(const SimulationParameters& params, const std::string& geomfile);
    
    /**
     * @brief Create SPIDER geometry with custom parameters
     * @param params Simulation parameters
     * @param geomfile Output geometry file path
     * @param z_start Starting z position
     * @param z_end Ending z position
     * @param meshsize_multiplier Mesh size scaling factor
     */
    void createGeometrySPIDER(const SimulationParameters& params, const std::string& geomfile,
                             double z_start, double z_end, double meshsize_multiplier);

    /**
     * @brief Create MITICA geometry
     * @param params Simulation parameters
     * @param geomfile Output geometry file path
     * @param z_start Starting z position
     * @param z_end Ending z position
     * @param meshsize_multiplier Mesh size scaling factor
     */
    void createGeometryMITICA(const SimulationParameters& params, const std::string& geomfile,
                             double z_start, double z_end, double meshsize_multiplier);

    /**
     * @brief Create MTF geometry
     * @param params Simulation parameters
     * @param geomfile Output geometry file path
     * @param z_start Starting z position
     * @param z_end Ending z position
     * @param meshsize_multiplier Mesh size scaling factor
     */
    void createGeometryMTF(const SimulationParameters& params, const std::string& geomfile,
                          double z_start, double z_end, double meshsize_multiplier);

    // Getters
    Geometry* getGeometry() const { return geometry; }
    MITICA_Accelerator* getAccelerator() const { return accelerator; }
    const std::vector<double>& getZGrids() const { return zgrids; }

    // ParaView/VTK export methods
    /**
     * @brief Export geometry to VTK format for ParaView visualization (simple solid/vacuum)
     * @param filename Output VTK filename (without extension)
     */
    void exportGeometryToVTK(const std::string& filename);
    
    /**
     * @brief Export detailed geometry to VTK format with individual solid IDs
     * @param filename Output VTK filename (without extension)
     */
    void exportDetailedGeometryToVTK(const std::string& filename);
    
    /**
     * @brief Export potential field to VTK format for ParaView visualization
     * @param epot Electric potential field
     * @param filename Output VTK filename (without extension)
     */
    void exportPotentialToVTK(const class EpotField& epot, const std::string& filename);
    
    /**
     * @brief Export solid boundaries to VTK format for ParaView visualization
     * @param filename Output VTK filename (without extension)
     */
    void exportSolidsToVTK(const std::string& filename);
    
    /**
     * @brief Export space-charge field to VTK format for ParaView visualization
     * @param scharge Space-charge field
     * @param filename Output VTK filename (without extension)
     */
    void exportSpacechargeToVTK(const class MeshScalarField& scharge, const std::string& filename);

    // Cleanup
    void resetGeometry();
};

#endif /* GEOMETRYMANAGER_H_ */
