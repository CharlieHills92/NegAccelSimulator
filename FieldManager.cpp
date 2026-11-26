/*
 * FieldManager.cpp
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored from ManageSimulation)
 */

#include "FieldManager.h"
#include "SimulationParameters.h"
#include "globals.h"

#include "meshvectorfield.hpp"
#include "epot_efield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "funct.h"

#include <map>
#include <vector>

using namespace std;

namespace {

string ensureTrailingSlash(const string& path) {
    if (path.empty() || path[path.size() - 1] == '/') {
        return path;
    }
    return path + "/";
}

} // namespace

FieldManager::FieldManager() : 
    potential(nullptr), spacecharge(nullptr), magnetic(nullptr),
    gas_dens(nullptr), electric(nullptr) {
}

FieldManager::~FieldManager() {
    resetFields();
}

void FieldManager::resetFields() {
    delete potential;
    potential = nullptr;
    delete spacecharge;
    spacecharge = nullptr;
    delete magnetic;
    magnetic = nullptr;
    delete gas_dens;
    gas_dens = nullptr;
    delete electric;
    electric = nullptr;
}

string FieldManager::getBFieldFolder(const SimulationParameters& params) {
    if (params.getBIsOn() == 0U) {
        return string();
    }
    return params.getMagneticFieldDirectory();
}

MeshVectorField* FieldManager::loadMagneticFieldDefinition(
    const SimulationParameters::MagneticFieldDefinition& definition,
    const string& base_directory,
    const MeshVectorField* reference_mesh) {
    constexpr bool fout[3] = {true, true, true};
    constexpr double FIELD_SCALE = 1.0e-3;

    if (definition.sourceType == "constant") {
        MeshVectorField* constant_field = NULL;
        if (reference_mesh != NULL) {
            constant_field = new MeshVectorField(static_cast<const Mesh&>(*reference_mesh), fout);
        } else {
            constant_field = new MeshVectorField(MODE_3D, fout, Int3D(1, 1, 1), Vec3D(0.0, 0.0, 0.0), 1.0);
        }

        Vec3D field_value(definition.constantValue[0] * definition.scale,
                          definition.constantValue[1] * definition.scale,
                          definition.constantValue[2] * definition.scale);

        Int3D size = constant_field->size();
        for (int32_t i = 0; i < size[0]; ++i) {
            for (int32_t j = 0; j < size[1]; ++j) {
                for (int32_t k = 0; k < size[2]; ++k) {
                    constant_field->set(i, j, k, field_value);
                }
            }
        }
        return constant_field;
    }

    if (definition.sourceType == "file") {
        string filename = definition.filePath;
        if (!filename.empty() && filename[0] != '/' && !base_directory.empty()) {
            filename = ensureTrailingSlash(base_directory) + filename;
        }
        return new MeshVectorField(MODE_3D, fout, FIELD_SCALE, definition.scale, filename);
    }

    throw Error(ERROR_LOCATION, "Unsupported magnetic field sourceType: " + definition.sourceType);
}

void FieldManager::addMagneticField(const SimulationParameters& params, const string& bfield_fold) {
    const string bfield_folder = ensureTrailingSlash(bfield_fold);
    
    if (!bfield_folder.empty()) {
        ibsimu.message(1) << "Loading magnetic field from: " << bfield_folder << endl;
    }
    
    MeshVectorField* bfield = nullptr;
    
    if(params.getBIsOn() == 0) {
        bfield = new MeshVectorField();
        ibsimu.message(1) << "\t*** NO MAGNETIC FIELD SET ***" << endl;
        logfile << "\t*** NO MAGNETIC FIELD SET ***" << endl << flush;
    }
    else {
        vector<MeshVectorField*> tempBfield;
        
        try {
            const vector<SimulationParameters::MagneticFieldDefinition>& definitions = params.getMagneticFields();

            // Load file-backed fields first so constants can reuse a compatible mesh for summation.
            for (size_t ii = 0; ii < definitions.size(); ++ii) {
                const SimulationParameters::MagneticFieldDefinition& definition = definitions[ii];
                if (definition.sourceType != "file") {
                    continue;
                }
                tempBfield.push_back(loadMagneticFieldDefinition(definition, bfield_folder, NULL));
                ibsimu.message(1) << "\tMagnetic field loaded: " << definition.name
                                  << " (type: " << definition.sourceType
                                  << ", scale: " << definition.scale << ")" << endl;
                logfile << "\tMagnetic field loaded: " << definition.name
                        << " (type: " << definition.sourceType
                        << ", scale: " << definition.scale << ")" << endl << flush;
            }

            const MeshVectorField* reference_mesh = tempBfield.empty() ? NULL : tempBfield[0];
            for (size_t ii = 0; ii < definitions.size(); ++ii) {
                const SimulationParameters::MagneticFieldDefinition& definition = definitions[ii];
                if (definition.sourceType != "constant") {
                    continue;
                }
                tempBfield.push_back(loadMagneticFieldDefinition(definition, bfield_folder, reference_mesh));
                if (reference_mesh == NULL && !tempBfield.empty()) {
                    reference_mesh = tempBfield[0];
                }
                ibsimu.message(1) << "\tMagnetic field loaded: " << definition.name
                                  << " (type: " << definition.sourceType
                                  << ", scale: " << definition.scale << ")" << endl;
                logfile << "\tMagnetic field loaded: " << definition.name
                        << " (type: " << definition.sourceType
                        << ", scale: " << definition.scale << ")" << endl << flush;
            }
            
            // Combine all loaded fields
            if (!tempBfield.empty()) {
                bfield = tempBfield[0];  // Take ownership of first field
                
                // Add remaining fields and clean up
                for (size_t ii = 1; ii < tempBfield.size(); ++ii) {
                    *bfield += *tempBfield[ii];
                    delete tempBfield[ii];
                }
                
                ibsimu.message(1) << "\t*** MAGNETIC FIELD COMBINED ***" << endl;
            } else {
                bfield = new MeshVectorField();
                ibsimu.message(1) << "\t*** NO MAGNETIC FIELDS LOADED ***" << endl;
            }
            
        } catch (const exception& e) {
            // Clean up on error
            for (auto* field : tempBfield) {
                delete field;
            }
            throw Error(ERROR_LOCATION, "Error loading magnetic fields: " + string(e.what()));
        }
        
        logfile << "\t*** MAGNETIC FIELD PROCESSING COMPLETE ***" << endl << flush;
    }
    
    // Set field extrapolation
    constexpr field_extrpl_e bfldextrpl[6] = { 
        FIELD_ZERO, FIELD_ZERO, FIELD_ZERO, 
        FIELD_ZERO, FIELD_ZERO, FIELD_ZERO 
    };
    bfield->set_extrapolation(bfldextrpl);
    
    magnetic = bfield;
}
