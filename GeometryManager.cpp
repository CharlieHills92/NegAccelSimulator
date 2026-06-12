/*
 * GeometryManager.cpp
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored from ManageSimulation)
 */

#include "GeometryManager.h"
#include "SimulationParameters.h"

#include "geometry.hpp"
#include "func_solid.hpp"
#include "geom/geom_function.h"
#include "ibsimu.hpp"
#include "error.hpp"
#include "epot_field.hpp"
#include "particledatabase.hpp"
#include "particles.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>

using namespace std;

namespace {

std::string uppercaseCopy(const std::string& value) {
    std::string upper = value;
    std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char character) {
        return static_cast<char>(std::toupper(character));
    });
    return upper;
}

Bound resolveConfiguredBoundary(const SimulationParameters& params,
                               int boundary_id,
                               const Bound& fallback) {
    SimulationParameters::BoundaryConditionDefinition boundary_definition;
    if (!params.tryGetBoundaryCondition(boundary_id, boundary_definition)) {
        return fallback;
    }

    if (boundary_definition.conditionType == "neumann") {
        return Bound(BOUND_NEUMANN, boundary_definition.value);
    }

    return Bound(BOUND_DIRICHLET, boundary_definition.value);
}

void applyConfiguredBoundary(Geometry* geom,
                             const SimulationParameters& params,
                             int boundary_id,
                             const Bound& fallback) {
    geom->set_boundary(boundary_id, resolveConfiguredBoundary(params, boundary_id, fallback));
}

double resolveGeneratedSolidVoltage(const SimulationParameters& params,
                                   const SimulationParameters::GeometrySolidDefinition& solid) {
    if (solid.hasExplicitVoltage) {
        return solid.voltageVolts;
    }

    const std::string role = solid.role;
    const std::string name = uppercaseCopy(solid.name);

    if (role == "plasma" || role == "wall" || name == "PG") {
        return 0.0;
    }
    if (role == "extraction_grid" || name == "EG") {
        return params.getEGVoltage();
    }
    if (role == "ground_grid" || name == "GG" || name == "G1" || name == "AG1") {
        return params.getG1Voltage();
    }
    if (role == "repeller" || name == "REP" || name == "G2" || name == "AG2") {
        return params.getG2Voltage();
    }
    if (role == "accelerator_grid") {
        switch (solid.stage) {
            case 1:
                return params.getG1Voltage();
            case 2:
                return params.getG2Voltage();
            case 3:
                return params.getG3Voltage();
            case 4:
                return params.getG4Voltage();
            case 5:
                return params.getG5Voltage();
            default:
                break;
        }
    }
    if (name == "G3" || name == "AG3") {
        return params.getG3Voltage();
    }
    if (name == "G4" || name == "AG4") {
        return params.getG4Voltage();
    }
    if (name == "G5" || name == "AG5") {
        return params.getG5Voltage();
    }

    return 0.0;
}

double resolveGeneratedBoundaryVoltage(
    const SimulationParameters& params,
    const std::vector<SimulationParameters::GeometrySolidDefinition>& grouped_solids,
    uint32_t boundary_id) {
    const double tolerance = 1.0e-9;
    bool has_explicit_voltage = false;
    double explicit_voltage = 0.0;
    bool has_named_voltage = false;
    double named_voltage = 0.0;
    bool has_any_fallback = false;
    double fallback_voltage = 0.0;

    for (std::vector<SimulationParameters::GeometrySolidDefinition>::const_iterator solid_it =
             grouped_solids.begin();
         solid_it != grouped_solids.end();
         ++solid_it) {
        if (solid_it->hasExplicitVoltage) {
            if (!has_explicit_voltage) {
                explicit_voltage = solid_it->voltageVolts;
                has_explicit_voltage = true;
            } else if (std::abs(explicit_voltage - solid_it->voltageVolts) > tolerance) {
                throw Error(ERROR_LOCATION,
                            "Conflicting explicit voltages for generated solids sharing boundaryId " +
                                std::to_string(boundary_id));
            }
            continue;
        }

        const double resolved_voltage = resolveGeneratedSolidVoltage(params, *solid_it);
        if (!has_any_fallback) {
            fallback_voltage = resolved_voltage;
            has_any_fallback = true;
        }
        if (std::abs(resolved_voltage) <= tolerance) {
            continue;
        }
        if (!has_named_voltage) {
            named_voltage = resolved_voltage;
            has_named_voltage = true;
        } else if (std::abs(named_voltage - resolved_voltage) > tolerance) {
            throw Error(ERROR_LOCATION,
                        "Conflicting inferred voltages for generated solids sharing boundaryId " +
                            std::to_string(boundary_id));
        }
    }

    if (has_explicit_voltage) {
        return explicit_voltage;
    }
    if (has_named_voltage) {
        return named_voltage;
    }
    return has_any_fallback ? fallback_voltage : 0.0;
}

uint32_t resolveGeneratedSolidBoundaryId(const SimulationParameters::GeometrySolidDefinition& solid,
                                         size_t index) {
    if (solid.boundaryId >= 7) {
        return static_cast<uint32_t>(solid.boundaryId);
    }
    return static_cast<uint32_t>(7 + index);
}

bool evaluateGeneratedSolidProfile(double x,
                                   double y,
                                   double z,
                                   const std::vector<double>& z_profile,
                                   const std::vector<double>& r_profile,
                                   const std::vector<double>& rounding_radii,
                                   const SimulationParameters::GeometryAperturePattern& aperture_pattern) {
    if (z_profile.empty() || z < z_profile.front() || z > z_profile.back()) {
        return false;
    }

    double local_x = x - aperture_pattern.xOffsetMeters;
    double local_y = y - aperture_pattern.yOffsetMeters;
    const std::string& layout = aperture_pattern.layout;
    if (layout == "staggered-grid") {
        const double pitch_y = aperture_pattern.pitchYMeters;
        const int row = static_cast<int>(std::floor((local_y + pitch_y / 2.0) / pitch_y));
        local_x += (row % 2 == 0)
            ? aperture_pattern.rowShiftXMeters
            : -aperture_pattern.rowShiftXMeters;
    }

    if (layout != "single") {
        const double half_x = aperture_pattern.countX * aperture_pattern.pitchXMeters / 2.0;
        const double half_y = aperture_pattern.countY * aperture_pattern.pitchYMeters / 2.0;
        if (std::abs(local_x) > half_x + aperture_pattern.marginMeters ||
            std::abs(local_y) > half_y + aperture_pattern.marginMeters) {
            return aperture_pattern.outsidePatternIsSolid;
        }

        if (std::abs(local_x) < half_x && std::abs(local_y) < half_y) {
            local_x -= std::round(local_x / aperture_pattern.pitchXMeters) * aperture_pattern.pitchXMeters;
            local_y -= std::round(local_y / aperture_pattern.pitchYMeters) * aperture_pattern.pitchYMeters;
        }
    }

    const double radius = std::sqrt(local_x * local_x + local_y * local_y);
    return inSolidFromHoleProfile(z, radius, z_profile, r_profile, rounding_radii);
}

bool evaluateGeneratedSolid(double x,
                            double y,
                            double z,
                            const SimulationParameters::GeometrySolidDefinition& solid) {
    return evaluateGeneratedSolidProfile(
        x,
        y,
        z,
        solid.zProfileMeters,
        solid.rProfileMeters,
        solid.roundingRadiiMeters,
        solid.aperturePattern);
}

std::vector<double> buildGeneratedZGrids(
    const std::vector<SimulationParameters::GeometrySolidDefinition>& solids) {
    std::vector<double> z_grid_pairs;
    z_grid_pairs.reserve(solids.size() * 2);

    for (std::vector<SimulationParameters::GeometrySolidDefinition>::const_iterator solid_it = solids.begin();
         solid_it != solids.end(); ++solid_it) {
        z_grid_pairs.push_back(solid_it->zProfileMeters.front());
        z_grid_pairs.push_back(solid_it->zProfileMeters.back());
    }

    return z_grid_pairs;
}

void applyDefaultDomainBoundaries(Geometry* geom, const SimulationParameters& params) {
    applyConfiguredBoundary(geom, params, 1, Bound(BOUND_NEUMANN, 0.0));
    applyConfiguredBoundary(geom, params, 2, Bound(BOUND_NEUMANN, 0.0));
    applyConfiguredBoundary(geom, params, 3, Bound(BOUND_NEUMANN, 0.0));
    applyConfiguredBoundary(geom, params, 4, Bound(BOUND_NEUMANN, 0.0));
    applyConfiguredBoundary(geom, params, 5, Bound(BOUND_DIRICHLET, 0.0));
    applyConfiguredBoundary(geom, params, 6, Bound(BOUND_NEUMANN, 0.0));
}

} // namespace

GeometryManager::GeometryManager() : 
    geometry(nullptr) {
}

GeometryManager::~GeometryManager() {
    resetGeometry();
}

void GeometryManager::resetGeometry() {
    delete geometry;
    geometry = nullptr;
}

void GeometryManager::createGeometry(const SimulationParameters& params, const string& geomfile) {
    createGeometry(params, geomfile, params.getDomainZStart(), params.getDomainZSizeOrDefault(), 1.0);
}

void GeometryManager::createGeometry(const SimulationParameters& params, const string& geomfile,
                                    double z_start, double z_end, double meshsize_multiplier) {
    if (params.hasGeneratedGeometrySolids()) {
        createGeneratedGeometry(params, geomfile, z_start, z_end, meshsize_multiplier);
        return;
    }

    throw Error(ERROR_LOCATION,
                "Only explicit geometry.source.mode='generated-data' is supported by the current runtime");
}

void GeometryManager::createGeneratedGeometry(const SimulationParameters& params, const string& geomfile,
                                             double z_start, double z_end, double meshsize_multiplier) {
    resetGeometry();

    const std::vector<SimulationParameters::GeometrySolidDefinition>& solids =
        params.getGeneratedGeometrySolids();
    if (solids.empty()) {
        throw Error(ERROR_LOCATION, "Generated geometry requested without any solid definitions");
    }

    const double mesh_size = params.getMeshSize() * meshsize_multiplier;
    const double x_size = params.getDomainXSizeOrDefault();
    const double y_size = params.getDomainYSizeOrDefault();
    const double z_size = (z_end - z_start > 0.0) ? (z_end - z_start) : params.getDomainZSizeOrDefault();

    Int3D meshsize(
        static_cast<int>(floor(x_size / mesh_size)) + 1,
        static_cast<int>(floor(y_size / mesh_size)) + 1,
        static_cast<int>(floor(z_size / mesh_size)) + 1);
    Vec3D origo(-x_size / 2.0, -y_size / 2.0, z_start);

    Geometry* geom = new Geometry(MODE_3D, meshsize, origo, mesh_size);
    applyDefaultDomainBoundaries(geom, params);

    std::map<uint32_t, std::vector<SimulationParameters::GeometrySolidDefinition> > grouped_solids;
    for (size_t index = 0; index < solids.size(); ++index) {
        const SimulationParameters::GeometrySolidDefinition solid = solids[index];
        const uint32_t boundary_id = resolveGeneratedSolidBoundaryId(solid, index);
        grouped_solids[boundary_id].push_back(solid);
    }

    for (std::map<uint32_t, std::vector<SimulationParameters::GeometrySolidDefinition> >::const_iterator
             boundary_it = grouped_solids.begin();
         boundary_it != grouped_solids.end();
         ++boundary_it) {
        const uint32_t boundary_id = boundary_it->first;
        const std::vector<SimulationParameters::GeometrySolidDefinition> boundary_solids =
            boundary_it->second;

        Solid* generated_solid = new FuncSolid(
            [boundary_solids](double x, double y, double z) {
                for (std::vector<SimulationParameters::GeometrySolidDefinition>::const_iterator solid_it =
                         boundary_solids.begin();
                     solid_it != boundary_solids.end();
                     ++solid_it) {
                    if (evaluateGeneratedSolid(x, y, z, *solid_it)) {
                        return true;
                    }
                }
                return false;
            });
        geom->set_solid(boundary_id, generated_solid);
        applyConfiguredBoundary(
            geom,
            params,
            static_cast<int>(boundary_id),
            Bound(BOUND_DIRICHLET, resolveGeneratedBoundaryVoltage(params, boundary_solids, boundary_id)));

        std::string grouped_names;
        for (std::vector<SimulationParameters::GeometrySolidDefinition>::const_iterator solid_it =
                 boundary_solids.begin();
             solid_it != boundary_solids.end();
             ++solid_it) {
            if (!grouped_names.empty()) {
                grouped_names += ", ";
            }
            grouped_names += solid_it->name;
        }

        ibsimu.message(1) << "\tGenerated solid group [" << grouped_names << "] created (ID: "
                          << boundary_id << ")" << endl;
    }

    geom->build_mesh();
    geom->build_surface();
    geom->save(geomfile);

    geometry = geom;
    zgrids = buildGeneratedZGrids(solids);

    ibsimu.message(1) << "Generated geometry created and saved to: " << geomfile << endl;
}

void GeometryManager::exportGeometryToVTK(const string& filename) {
    if (!geometry) {
        ibsimu.message(1) << "No geometry available for VTK export" << endl;
        return;
    }
    
    string vtkfile = filename + "_geometry.vtk";
    
    try {
        // Get mesh information
        Int3D size = geometry->size();
        Vec3D origo = geometry->origo();
        double h = geometry->h();
        
        ofstream file(vtkfile, ios::binary);
        if (!file.is_open()) {
            throw Error("Cannot open VTK file for writing: " + vtkfile);
        }
        
        // Helper lambda to swap endianness for binary VTK (big-endian required)
        auto swap_endian_int = [](int value) -> int {
            return ((value >> 24) & 0x000000FF) | 
                   ((value >> 8)  & 0x0000FF00) | 
                   ((value << 8)  & 0x00FF0000) | 
                   ((value << 24) & 0xFF000000);
        };
        
        auto swap_endian_float = [](float value) -> float {
            union { float f; uint32_t i; } u;
            u.f = value;
            u.i = ((u.i >> 24) & 0x000000FF) | 
                  ((u.i >> 8)  & 0x0000FF00) | 
                  ((u.i << 8)  & 0x00FF0000) | 
                  ((u.i << 24) & 0xFF000000);
            return u.f;
        };
        
        // Write VTK header (ASCII)
        file << "# vtk DataFile Version 3.0" << endl;
        file << "IBSIMU Geometry Export (Solid/Vacuum)" << endl;
        file << "BINARY" << endl;
        file << "DATASET STRUCTURED_POINTS" << endl;
        
        // Write dimensions (ASCII)
        file << "DIMENSIONS " << size[0] << " " << size[1] << " " << size[2] << endl;
        float origo_f[3] = {static_cast<float>(origo[0]), static_cast<float>(origo[1]), static_cast<float>(origo[2])};
        file << "ORIGIN " << origo_f[0] << " " << origo_f[1] << " " << origo_f[2] << endl;
        float spacing = static_cast<float>(h);
        file << "SPACING " << spacing << " " << spacing << " " << spacing << endl;
        
        // Write point data header (ASCII)
        file << "POINT_DATA " << size[0] * size[1] * size[2] << endl;
        file << "SCALARS solid_vacuum int" << endl;
        file << "LOOKUP_TABLE default" << endl;
        
        // Export solid/vacuum distinction using IBSIMU's inside() function (binary)
        for (int k = 0; k < size[2]; k++) {
            for (int j = 0; j < size[1]; j++) {
                for (int i = 0; i < size[0]; i++) {
                    // Calculate world coordinates from mesh indices
                    Vec3D point(origo[0] + i * h, origo[1] + j * h, origo[2] + k * h);
                    
                    // Use IBSIMU's inside() function to get solid ID
                    uint32_t solid_id = geometry->inside(point);
                    // Convert to binary: 0 = vacuum/boundary, 1 = electrode solid
                    int solid_indicator = (solid_id > 6) ? 1 : 0;
                    int solid_indicator_be = swap_endian_int(solid_indicator);
                    file.write(reinterpret_cast<const char*>(&solid_indicator_be), sizeof(int));
                }
            }
        }
        
        file.close();
        ibsimu.message(1) << "Geometry exported to binary VTK file: " << vtkfile << endl;
        
    } catch (const Error& e) {
        ibsimu.message(1) << "Error exporting geometry to VTK" << endl;
    }
}

void GeometryManager::exportDetailedGeometryToVTK(const string& filename) {
    if (!geometry) {
        ibsimu.message(1) << "No geometry available for detailed VTK export" << endl;
        return;
    }
    
    string vtkfile = filename + "_geometry_detailed.vtk";
    
    try {
        // Get mesh information
        Int3D size = geometry->size();
        Vec3D origo = geometry->origo();
        double h = geometry->h();
        
        ofstream file(vtkfile, ios::binary);
        if (!file.is_open()) {
            throw Error("Cannot open VTK file for writing: " + vtkfile);
        }
        
        // Helper lambda to swap endianness
        auto swap_endian_int = [](int value) -> int {
            return ((value >> 24) & 0x000000FF) | 
                   ((value >> 8)  & 0x0000FF00) | 
                   ((value << 8)  & 0x00FF0000) | 
                   ((value << 24) & 0xFF000000);
        };
        
        // Write VTK header (ASCII)
        file << "# vtk DataFile Version 3.0" << endl;
        file << "IBSIMU Detailed Geometry Export (Individual Solids)" << endl;
        file << "BINARY" << endl;
        file << "DATASET STRUCTURED_POINTS" << endl;
        
        // Write dimensions (ASCII)
        file << "DIMENSIONS " << size[0] << " " << size[1] << " " << size[2] << endl;
        float origo_f[3] = {static_cast<float>(origo[0]), static_cast<float>(origo[1]), static_cast<float>(origo[2])};
        file << "ORIGIN " << origo_f[0] << " " << origo_f[1] << " " << origo_f[2] << endl;
        float spacing = static_cast<float>(h);
        file << "SPACING " << spacing << " " << spacing << " " << spacing << endl;
        
        // Write multiple scalar fields for comprehensive analysis
        file << "POINT_DATA " << size[0] * size[1] * size[2] << endl;
        
        // First field: Solid IDs using IBSIMU's inside() function (binary)
        file << "SCALARS solid_id int" << endl;
        file << "LOOKUP_TABLE default" << endl;
        
        for (int k = 0; k < size[2]; k++) {
            for (int j = 0; j < size[1]; j++) {
                for (int i = 0; i < size[0]; i++) {
                    // Calculate world coordinates from mesh indices
                    Vec3D point(origo[0] + i * h, origo[1] + j * h, origo[2] + k * h);
                    
                    // Use IBSIMU's inside() function to get solid ID
                    uint32_t solid_id = geometry->inside(point);
                    int solid_id_int = static_cast<int>(solid_id);
                    int solid_id_int_be = swap_endian_int(solid_id_int);
                    file.write(reinterpret_cast<const char*>(&solid_id_int_be), sizeof(int));
                }
            }
        }
        
        // Second field: Binary solid/vacuum for easy thresholding
        file << "\nSCALARS is_solid int" << endl;
        file << "LOOKUP_TABLE default" << endl;
        
        for (int k = 0; k < size[2]; k++) {
            for (int j = 0; j < size[1]; j++) {
                for (int i = 0; i < size[0]; i++) {
                    // Calculate world coordinates from mesh indices
                    Vec3D point(origo[0] + i * h, origo[1] + j * h, origo[2] + k * h);
                    
                    // Use IBSIMU's inside() function and convert to binary
                    uint32_t solid_id = geometry->inside(point);
                    int is_solid = (solid_id > 6) ? 1 : 0; // Electrodes are > 6, boundaries/vacuum are 0-6
                    int is_solid_be = swap_endian_int(is_solid);
                    file.write(reinterpret_cast<const char*>(&is_solid_be), sizeof(int));
                }
            }
        }
        
        file.close();
        ibsimu.message(1) << "Detailed geometry exported to binary VTK file: " << vtkfile << endl;
        
    } catch (const Error& e) {
        ibsimu.message(1) << "Error exporting detailed geometry to VTK" << endl;
    }
}

void GeometryManager::exportPotentialToVTK(const EpotField& epot, const string& filename) {
    if (!geometry) {
        ibsimu.message(1) << "No geometry available for potential VTK export" << endl;
        return;
    }
    
    string vtkfile = filename + "_potential.vtk";
    
    try {
        // Get mesh information
        Int3D size = geometry->size();
        Vec3D origo = geometry->origo();
        double h = geometry->h();
        
        ofstream file(vtkfile, ios::binary);
        if (!file.is_open()) {
            throw Error("Cannot open VTK file for writing: " + vtkfile);
        }
        
        // Helper lambda to swap endianness
        auto swap_endian_float = [](float value) -> float {
            union { float f; uint32_t i; } u;
            u.f = value;
            u.i = ((u.i >> 24) & 0x000000FF) | 
                  ((u.i >> 8)  & 0x0000FF00) | 
                  ((u.i << 8)  & 0x00FF0000) | 
                  ((u.i << 24) & 0xFF000000);
            return u.f;
        };
        
        // Write VTK header (ASCII)
        file << "# vtk DataFile Version 3.0" << endl;
        file << "IBSIMU Potential Field Export" << endl;
        file << "BINARY" << endl;
        file << "DATASET STRUCTURED_POINTS" << endl;
        
        // Write dimensions (ASCII)
        file << "DIMENSIONS " << size[0] << " " << size[1] << " " << size[2] << endl;
        float origo_f[3] = {static_cast<float>(origo[0]), static_cast<float>(origo[1]), static_cast<float>(origo[2])};
        file << "ORIGIN " << origo_f[0] << " " << origo_f[1] << " " << origo_f[2] << endl;
        float spacing = static_cast<float>(h);
        file << "SPACING " << spacing << " " << spacing << " " << spacing << endl;
        
        // Write point data header (ASCII)
        file << "POINT_DATA " << size[0] * size[1] * size[2] << endl;
        file << "SCALARS potential float" << endl;
        file << "LOOKUP_TABLE default" << endl;
        
        // Export potential values at each mesh point (binary)
        for (int k = 0; k < size[2]; k++) {
            for (int j = 0; j < size[1]; j++) {
                for (int i = 0; i < size[0]; i++) {
                    double potential = epot(i, j, k);
                    float potential_f = static_cast<float>(potential);
                    float potential_f_be = swap_endian_float(potential_f);
                    file.write(reinterpret_cast<const char*>(&potential_f_be), sizeof(float));
                }
            }
        }
        
        file.close();
        ibsimu.message(1) << "Potential field exported to binary VTK file: " << vtkfile << endl;
        
    } catch (const Error& e) {
        ibsimu.message(1) << "Error exporting potential to VTK" << endl;
    }
}

void GeometryManager::exportSolidsToVTK(const string& filename) {
    if (!geometry) {
        ibsimu.message(1) << "No geometry available for solids VTK export" << endl;
        return;
    }
    
    string vtkfile = filename + "_solids.vtk";
    
    try {
        // Get mesh information
        Int3D size = geometry->size();
        Vec3D origo = geometry->origo();
        double h = geometry->h();
        
        ofstream file(vtkfile, ios::binary);
        if (!file.is_open()) {
            throw Error("Cannot open VTK file for writing: " + vtkfile);
        }
        
        // Helper lambda to swap endianness
        auto swap_endian_int = [](int value) -> int {
            return ((value >> 24) & 0x000000FF) | 
                   ((value >> 8)  & 0x0000FF00) | 
                   ((value << 8)  & 0x00FF0000) | 
                   ((value << 24) & 0xFF000000);
        };
        
        // Write VTK header (ASCII)
        file << "# vtk DataFile Version 3.0" << endl;
        file << "IBSIMU Solid Boundaries Export" << endl;
        file << "BINARY" << endl;
        file << "DATASET STRUCTURED_POINTS" << endl;
        
        // Write dimensions (ASCII)
        file << "DIMENSIONS " << size[0] << " " << size[1] << " " << size[2] << endl;
        float origo_f[3] = {static_cast<float>(origo[0]), static_cast<float>(origo[1]), static_cast<float>(origo[2])};
        file << "ORIGIN " << origo_f[0] << " " << origo_f[1] << " " << origo_f[2] << endl;
        float spacing = static_cast<float>(h);
        file << "SPACING " << spacing << " " << spacing << " " << spacing << endl;
        
        // Write point data header (ASCII)
        file << "POINT_DATA " << size[0] * size[1] * size[2] << endl;
        file << "SCALARS material_id int" << endl;
        file << "LOOKUP_TABLE default" << endl;
        
        // Export solid IDs using IBSIMU's inside() function (binary)
        for (int k = 0; k < size[2]; k++) {
            for (int j = 0; j < size[1]; j++) {
                for (int i = 0; i < size[0]; i++) {
                    // Calculate world coordinates from mesh indices
                    Vec3D point(origo[0] + i * h, origo[1] + j * h, origo[2] + k * h);
                    
                    // Use IBSIMU's inside() function to get solid ID
                    // Returns 0 for vacuum, 1-6 for boundaries, solid number for electrodes
                    uint32_t solid_id = geometry->inside(point);
                    int solid_id_int = static_cast<int>(solid_id);
                    int solid_id_int_be = swap_endian_int(solid_id_int);
                    file.write(reinterpret_cast<const char*>(&solid_id_int_be), sizeof(int));
                }
            }
        }
        
        file.close();
        ibsimu.message(1) << "Solid boundaries exported to binary VTK file: " << vtkfile << endl;
        
    } catch (const Error& e) {
        ibsimu.message(1) << "Error exporting solids to VTK" << endl;
    }
}

void GeometryManager::exportSpacechargeToVTK(const MeshScalarField& scharge, const string& filename) {
    if (!geometry) {
        ibsimu.message(1) << "No geometry available for space-charge VTK export" << endl;
        return;
    }
    
    string vtkfile = filename + "_scharge.vtk";
    
    try {
        // Get mesh information
        Int3D size = geometry->size();
        Vec3D origo = geometry->origo();
        double h = geometry->h();
        
        ofstream file(vtkfile, ios::binary);
        if (!file.is_open()) {
            throw Error("Cannot open VTK file for writing: " + vtkfile);
        }
        
        // Helper lambda to swap endianness
        auto swap_endian_float = [](float value) -> float {
            union { float f; uint32_t i; } u;
            u.f = value;
            u.i = ((u.i >> 24) & 0x000000FF) | 
                  ((u.i >> 8)  & 0x0000FF00) | 
                  ((u.i << 8)  & 0x00FF0000) | 
                  ((u.i << 24) & 0xFF000000);
            return u.f;
        };
        
        // Write VTK header (ASCII)
        file << "# vtk DataFile Version 3.0" << endl;
        file << "IBSIMU Space Charge Field Export" << endl;
        file << "BINARY" << endl;
        file << "DATASET STRUCTURED_POINTS" << endl;
        
        // Write dimensions
        file << "DIMENSIONS " << size[0] << " " << size[1] << " " << size[2] << endl;
        file << "ORIGIN " << origo[0] << " " << origo[1] << " " << origo[2] << endl;
        file << "SPACING " << h << " " << h << " " << h << endl;
        
        // Write point data header (ASCII)
        file << "POINT_DATA " << size[0] * size[1] * size[2] << endl;
        file << "SCALARS space_charge float 1" << endl;
        file << "LOOKUP_TABLE default" << endl;
        
        // Export space-charge values at each mesh point (binary)
        for (int k = 0; k < size[2]; k++) {
            for (int j = 0; j < size[1]; j++) {
                for (int i = 0; i < size[0]; i++) {
                    float scharge_val = static_cast<float>(scharge(i, j, k));
                    float swapped_val = swap_endian_float(scharge_val);
                    file.write(reinterpret_cast<const char*>(&swapped_val), sizeof(float));
                }
            }
        }
        
        file.close();
        ibsimu.message(1) << "Space charge field exported to binary VTK file: " << vtkfile << endl;
        
    } catch (const Error& e) {
        ibsimu.message(1) << "Error exporting space charge to VTK" << endl;
    }
}
