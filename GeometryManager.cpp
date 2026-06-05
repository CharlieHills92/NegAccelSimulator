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
#include "geom/SPIDER_geom.h"
#include "geom/MITICA_geom2.h"
#include "geom/MTF_geom.h"
#include "ibsimu.hpp"
#include "error.hpp"
#include "epot_field.hpp"
#include "particledatabase.hpp"
#include "particles.hpp"

#include <map>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

GeometryManager::GeometryManager() : 
    geometry(nullptr), accelerator(nullptr) {
}

GeometryManager::~GeometryManager() {
    resetGeometry();
}

void GeometryManager::resetGeometry() {
    delete geometry;
    geometry = nullptr;
    delete accelerator;
    accelerator = nullptr;
}

void GeometryManager::createGeometry(const SimulationParameters& params, const string& geomfile) {
    createGeometry(params, geomfile, 0.0, params.getDomainZSizeOrDefault(), 1.0);
}

void GeometryManager::createGeometry(const SimulationParameters& params, const string& geomfile,
                                    double z_start, double z_end, double meshsize_multiplier) {
    const string& geometry_template = params.getGeometryTemplate();

    if (geometry_template == "SPIDER") {
            createGeometrySPIDER(params, geomfile, z_start, z_end, meshsize_multiplier);
            return;
    }
    if (geometry_template == "MITICA") {
            createGeometryMITICA(params, geomfile, z_start, z_end, meshsize_multiplier);
            return;
    }
    if (geometry_template == "MTF") {
            createGeometryMTF(params, geomfile, z_start, z_end, meshsize_multiplier);
            return;
    }

    const int accel_type = static_cast<int>(params.getAcceleratorIdx());
    switch(accel_type) {
        case 1:
            createGeometrySPIDER(params, geomfile, z_start, z_end, meshsize_multiplier);
            break;
        case 2:
            createGeometryMITICA(params, geomfile, z_start, z_end, meshsize_multiplier);
            break;
        case 3:
            createGeometryMTF(params, geomfile, z_start, z_end, meshsize_multiplier);
            break;
        default:
            throw Error(ERROR_LOCATION,
                        "No supported geometry selection found. geometry.source.template='" +
                            geometry_template + "', accelerator index=" + std::to_string(accel_type));
    }
}

void GeometryManager::createGeometrySPIDER(const SimulationParameters& params, const string& geomfile) {
    constexpr double DEFAULT_Z_START = 0.0;
    constexpr double DEFAULT_Z_END = 90.0e-3;
    createGeometrySPIDER(params, geomfile, DEFAULT_Z_START, DEFAULT_Z_END, 1.0);
}

void GeometryManager::createGeometrySPIDER(const SimulationParameters& params, const string& geomfile,
                                          double z_start, double z_end, double meshsize_multiplier) {
    // Apply mesh size multiplier
    double modified_mesh_size = params.getMeshSize() * meshsize_multiplier;
    
    // Use configurable domain sizes with accelerator-specific defaults
    double x_size = params.getDomainXSizeOrDefault();
    double y_size = params.getDomainYSizeOrDefault();
    double z_size = (z_end - z_start > 0) ? (z_end - z_start) : params.getDomainZSizeOrDefault();
    
    // Calculate mesh dimensions
    Int3D meshsize( 
        static_cast<int>(floor(x_size / modified_mesh_size)) + 1,
        static_cast<int>(floor(y_size / modified_mesh_size)) + 1,
        static_cast<int>(floor(z_size / modified_mesh_size)) + 1
    );
    
    Vec3D origo(-x_size/2, -y_size/2, z_start);
    
    ibsimu.message(1) << "SPIDER geometry:" << endl
                      << "  Domain size: " << x_size*1000 << "mm x " << y_size*1000 << "mm x " << z_size*1000 << "mm" << endl
                      << "\tMesh size: " << modified_mesh_size << " m" << endl
                      << "\tMesh points: " << meshsize[0] << "x" << meshsize[1] << "x" << meshsize[2] << endl;
    
    zgrids = {0.003, 0.009, 0.015, 0.032, 0.120, 0.137, 0.225, 0.242, 0.330, 0.347, 0.435, 0.452, 0.540, 0.557};

    createSPIDERGeometry(meshsize, origo, modified_mesh_size, params, geomfile);
}

void GeometryManager::createSPIDERGeometry(const Int3D& meshsize, const Vec3D& origo, 
                                          double mesh_size, const SimulationParameters& params,
                                          const string& geomfile) {
    Geometry* geom = new Geometry(MODE_3D, meshsize, origo, mesh_size);
    
    // Create solid objects
    struct SolidInfo {
        uint32_t id;
        bool (*solid_func)(double, double, double);  // Function pointer type
        const char* name;
        double voltage;
    };
    
    const vector<SolidInfo> solids = {
        {7U, solid0_SPIDER, "PG", 0.0},
        {8U, solid1_SPIDER, "EG", params.getEGVoltage()},
        {9U, solid2_SPIDER, "GG", params.getGGVoltage()}
    };
    
    for (const auto& solid_info : solids) {
        Solid* solid = new FuncSolid(solid_info.solid_func);
        geom->set_solid(solid_info.id, solid);
        geom->set_boundary(solid_info.id, Bound(BOUND_DIRICHLET, solid_info.voltage));
        
        ibsimu.message(1) << "\t" << solid_info.name << " solid created (ID: " 
                          << solid_info.id << ", V: " << solid_info.voltage << ")" << endl;
    }
    
    // Set boundary conditions for walls
    const vector<Bound> wall_boundaries = {
        Bound(BOUND_NEUMANN, 0.0),  // 1: x-min
        Bound(BOUND_NEUMANN, 0.0),  // 2: x-max  
        Bound(BOUND_NEUMANN, 0.0),  // 3: y-min
        Bound(BOUND_NEUMANN, 0.0),  // 4: y-max
        Bound(BOUND_DIRICHLET, 0.0), // 5: z-min
        Bound(BOUND_NEUMANN, 0.0),  // 6: z-max
    };
    
    for (size_t i = 0; i < wall_boundaries.size(); ++i) {
        geom->set_boundary(i + 1, wall_boundaries[i]);
    }
    
    geom->build_mesh();
    geom->build_surface();
    geom->save(geomfile);
    
    geometry = geom;
    
    ibsimu.message(1) << "SPIDER geometry created and saved to: " << geomfile << endl;
}

void GeometryManager::createGeometryMITICA(const SimulationParameters& params, const string& geomfile,
                                          double z_start, double z_end, double meshsize_multiplier) {
    double start = z_start;
    double h = params.getMeshSize() * meshsize_multiplier;
    
    // Use configurable domain sizes with accelerator-specific defaults
    double x_size = params.getDomainXSizeOrDefault();
    double y_size = params.getDomainYSizeOrDefault();
    double z_size = (z_end - z_start > 0.0) ? (z_end - z_start) : params.getDomainZSizeOrDefault();
    
    if (debug) cout << "DEBUG MITICA Geometry: Using domain sizes X=" << x_size << "m, Y=" << y_size << "m, Z=" << z_size << "m" << endl;
    if (debug) cout << "DEBUG MITICA Geometry: z_start=" << start << "m, z_end=" << z_end << "m, z_end-start=" << (z_end - start) << "m" << endl;
    
    double sizereq[3] = { x_size, y_size, z_size };
    
    Int3D meshsize((int)floor(sizereq[0]/h) + 1,
                   (int)floor(sizereq[1]/h) + 1,
                   (int)floor(sizereq[2]/h) + 1);
                   
    Vec3D origo(-sizereq[0]/2, -sizereq[1]/2, start);
    Geometry* geom = new Geometry(MODE_3D, meshsize, origo, h);
    
    ibsimu.message(1) << "MITICA geometry mesh size " << h << "\n";
    ibsimu.message(1) << "MITICA domain size: " << x_size*1000 << "mm x " << y_size*1000 << "mm x " << z_size*1000 << "mm" << endl;

    MITICA_Accelerator* acc = new MITICA_Accelerator(params.getExtGap(), params.getAccGap());

    zgrids = {0.003, 0.009, 0.015, 0.032, 0.120, 0.137, 0.225, 0.242, 0.330, 0.347, 0.435, 0.452, 0.540, 0.557};

    Solid *s0 = new FuncSolid(acc->create_PGsolid());
    geom->set_solid(7U, s0); // PG
    Solid *s1 = new FuncSolid(acc->create_EGsolid());
    geom->set_solid(8U, s1); // EG
    Solid *s2 = new FuncSolid(acc->create_AG1solid());
    geom->set_solid(9U, s2); // G1
    Solid *s3 = new FuncSolid(acc->create_AG2solid());
    geom->set_solid(10U, s3); // G2
    Solid *s4 = new FuncSolid(acc->create_AG3solid());
    geom->set_solid(11U, s4); // G3
    Solid *s5 = new FuncSolid(acc->create_AG4solid());
    geom->set_solid(12U, s5); // G4
    Solid *s6 = new FuncSolid(acc->create_AG5solid());
    geom->set_solid(13U, s6); // G5
    
    cout << " GEOMETRY CREATED \n";
    
    geom->set_boundary(1, Bound(BOUND_NEUMANN, 0.0));
    geom->set_boundary(2, Bound(BOUND_NEUMANN, 0.0));
    geom->set_boundary(3, Bound(BOUND_NEUMANN, 0.0));
    geom->set_boundary(4, Bound(BOUND_NEUMANN, 0.0));
    geom->set_boundary(5, Bound(BOUND_DIRICHLET, 0.0));
    geom->set_boundary(6, Bound(BOUND_NEUMANN, 0.0));
    
    geom->set_boundary(7, Bound(BOUND_DIRICHLET, 0.0));
    geom->set_boundary(8, Bound(BOUND_DIRICHLET, params.getEGVoltage()));
    geom->set_boundary(9, Bound(BOUND_DIRICHLET, params.getG1Voltage()));
    geom->set_boundary(10, Bound(BOUND_DIRICHLET, params.getG2Voltage()));
    geom->set_boundary(11, Bound(BOUND_DIRICHLET, params.getG3Voltage()));
    geom->set_boundary(12, Bound(BOUND_DIRICHLET, params.getG4Voltage()));
    geom->set_boundary(13, Bound(BOUND_DIRICHLET, params.getG5Voltage()));

    cout << " MITICA GEOMETRY CREATED \n";
    
    geom->build_mesh();
    geom->build_surface();
    geom->save(geomfile);

    geometry = geom;
    accelerator = acc;
    
    ibsimu.message(1) << "MITICA geometry created and saved to: " << geomfile << endl;
}

void GeometryManager::createGeometryMTF(const SimulationParameters& params, const string& geomfile,
                                       double z_start, double z_end, double meshsize_multiplier) {
    double start = z_start;
    double h = params.getMeshSize() * meshsize_multiplier;
    
    // Use configurable domain sizes with accelerator-specific defaults
    double x_size = params.getDomainXSizeOrDefault();
    double y_size = params.getDomainYSizeOrDefault();
    double z_size = (z_end - z_start > 0.0) ? (z_end - z_start) : params.getDomainZSizeOrDefault();
    
    if (debug) cout << "DEBUG MTF Geometry: Using domain sizes X=" << x_size << "m, Y=" << y_size << "m, Z=" << z_size << "m" << endl;
    if (debug) cout << "DEBUG MTF Geometry: z_start=" << start << "m, z_end=" << z_end << "m, z_end-start=" << (z_end - start) << "m" << endl;
    
    double sizereq[3] = { x_size, y_size, z_size };
    
    Int3D meshsize((int)floor(sizereq[0]/h) + 1,
                   (int)floor(sizereq[1]/h) + 1,
                   (int)floor(sizereq[2]/h) + 1);
                   
    Vec3D origo(-sizereq[0]/2, -sizereq[1]/2, start);
    Geometry* geom = new Geometry(MODE_3D, meshsize, origo, h);
    
    ibsimu.message(1) << "MTF geometry mesh size " << h << "\n";
    ibsimu.message(1) << "MTF domain size: " << x_size*1000 << "mm x " << y_size*1000 << "mm x " << z_size*1000 << "mm" << endl;

    // MTF specific z-grids for potential calculation
    zgrids = {0.003, 0.009, 0.015, 0.0315, 0.1195, 0.1365, 0.2245, 0.2415, 0.3295, 0.3465, 0.4345, 0.4515, 0.5395, 0.5565};

    // Create MTF solids using the MTF geometry functions
    Solid *s0 = new FuncSolid(solid0_MTF);
    geom->set_solid(7U, s0); // PG
    Solid *s1 = new FuncSolid(solid1_MTF);
    geom->set_solid(8U, s1); // EG
    Solid *s2 = new FuncSolid(solid2_MTF);
    geom->set_solid(9U, s2); // G1
    Solid *s3 = new FuncSolid(solid3_MTF);
    geom->set_solid(10U, s3); // G2
    Solid *s4 = new FuncSolid(solid4_MTF);
    geom->set_solid(11U, s4); // G3
    Solid *s5 = new FuncSolid(solid5_MTF);
    geom->set_solid(12U, s5); // G4
    Solid *s6 = new FuncSolid(solid6_MTF);
    geom->set_solid(13U, s6); // G5

    // Set boundary conditions for walls
    geom->set_boundary(1, Bound(BOUND_NEUMANN, 0.0));
    geom->set_boundary(2, Bound(BOUND_NEUMANN, 0.0));
    geom->set_boundary(3, Bound(BOUND_NEUMANN, 0.0));
    geom->set_boundary(4, Bound(BOUND_NEUMANN, 0.0));
    
    // Set boundary condition for z-min (depends on start position)
    if (z_start < 0.02) {
        geom->set_boundary(5, Bound(BOUND_DIRICHLET, 0.0));
    } else {
        geom->set_boundary(5, Bound(BOUND_DIRICHLET, params.getEGVoltage()));
    }
    
    // Set boundary conditions for electrodes
    geom->set_boundary(7, Bound(BOUND_DIRICHLET, 0.0));                    // PG
    geom->set_boundary(8, Bound(BOUND_DIRICHLET, params.getEGVoltage()));  // EG
    geom->set_boundary(9, Bound(BOUND_DIRICHLET, params.getG1Voltage()));  // G1
    geom->set_boundary(10, Bound(BOUND_DIRICHLET, params.getG2Voltage())); // G2
    geom->set_boundary(11, Bound(BOUND_DIRICHLET, params.getG3Voltage())); // G3
    geom->set_boundary(12, Bound(BOUND_DIRICHLET, params.getG4Voltage())); // G4
    geom->set_boundary(13, Bound(BOUND_DIRICHLET, params.getG5Voltage())); // G5

    cout << " MTF GEOMETRY CREATED \n";
    
    geom->build_mesh();
    geom->build_surface();
    geom->save(geomfile);

    geometry = geom;
    
    ibsimu.message(1) << "MTF geometry created and saved to: " << geomfile << endl;
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
