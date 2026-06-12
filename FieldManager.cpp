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

    if (params.getMagneticFieldSourceMode() == "directory") {
        return params.getMagneticFieldDirectory();
    }
    if (params.getMagneticFieldSourceMode() == "file") {
        return params.getMagneticFieldFile();
    }

    throw Error(ERROR_LOCATION,
                "Magnetic field is enabled but no explicit source path was configured");
}

void FieldManager::addMagneticField(const SimulationParameters& params, const string& bfield_fold) {
    const string& source_mode = params.getMagneticFieldSourceMode();
    const string bfield_folder = ensureTrailingSlash(bfield_fold);

    constexpr bool fout[3] = {true, true, true};
    constexpr double FIELD_SCALE = 1.0e-3;
    
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
            if (source_mode == "file") {
                if (params.getPGFilterScale() != 0.0 || params.getCESMADCMScale() != 0.0) {
                    throw Error(ERROR_LOCATION,
                                "PG filter and CESM+ADCM companion fields require magneticField.sourceMode='directory'");
                }

                ibsimu.message(1) << "\tLoading magnetic field from file " << bfield_fold << "..." << endl;
                logfile << "\tLoading magnetic field from file " << bfield_fold << "..." << endl << flush;
                tempBfield.push_back(new MeshVectorField(MODE_3D, fout, FIELD_SCALE,
                                                        params.getExtFieldScale(), bfield_fold));
                ibsimu.message(1) << "\tMagnetic field loaded: " << bfield_fold
                                  << " (scale: " << params.getExtFieldScale() << ")" << endl;
                logfile << "\tMagnetic field loaded with scaling " << params.getExtFieldScale() << endl << flush;
            } else {
                ibsimu.message(1) << "\tLoading fields from folder " << bfield_folder << "..." << endl;
                logfile << "\tLoading fields from folder " << bfield_folder << "..." << endl << flush;

                if (params.getPGFilterScale() != 0.0) {
                    const string pg_filename = bfield_folder + "PGfilter.fld";
                    tempBfield.push_back(new MeshVectorField(MODE_3D, fout, FIELD_SCALE,
                                                            params.getPGFilterScale(), pg_filename));
                    ibsimu.message(1) << "\tPG filter field loaded: " << pg_filename
                                      << " (scale: " << params.getPGFilterScale() << ")" << endl;
                    logfile << "\tPG filter field loaded with scaling " << params.getPGFilterScale() << endl << flush;
                }

                if (params.getCESMADCMScale() != 0.0) {
                    const string cesm_filename = bfield_folder + "CESMADCMfield.fld";
                    tempBfield.push_back(new MeshVectorField(MODE_3D, fout, FIELD_SCALE,
                                                            params.getCESMADCMScale(), cesm_filename));
                    ibsimu.message(1) << "\tCESM+ADCM field loaded: " << cesm_filename
                                      << " (scale: " << params.getCESMADCMScale() << ")" << endl;
                    logfile << "\tCESM+ADCM field loaded with scaling " << params.getCESMADCMScale() << endl << flush;
                }

                if (params.getExtFieldCase() > 0) {
                    loadExternalField(tempBfield, fout, FIELD_SCALE, params, bfield_folder);
                }
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

void FieldManager::loadExternalField(vector<MeshVectorField*>& tempBfield, 
                                    const bool fout[3], double fieldScale,
                                    const SimulationParameters& params, const string& bfield_folder) {
    
    switch (params.getExtFieldCase()) {
        case 1: {
            const string ext_filename = bfield_folder + "EXTfield.fld";
            tempBfield.push_back(new MeshVectorField(MODE_3D, fout, fieldScale, 
                                                    params.getExtFieldScale(), ext_filename));
            ibsimu.message(1) << "\tExternal field loaded: " << ext_filename 
                              << " (scale: " << params.getExtFieldScale() << ")" << endl;
            break;
        }
        case 2: {
            // Load both compensated and uncompensated fields
            const string nocomp_filename = bfield_folder + "EXTfield_nocomp_9p.fld";
            const string comp_filename = bfield_folder + "EXTfield_comp.fld";
            
            tempBfield.push_back(new MeshVectorField(MODE_3D, fout, fieldScale, 
                                                    1.0 - params.getExtFieldScale(), nocomp_filename));
            tempBfield.push_back(new MeshVectorField(MODE_3D, fout, fieldScale, 
                                                    params.getExtFieldScale(), comp_filename));
            
            ibsimu.message(1) << "\tExternal fields loaded:" << endl
                              << "\t\t" << nocomp_filename << " (scale: " << (1.0 - params.getExtFieldScale()) << ")" << endl
                              << "\t\t" << comp_filename << " (scale: " << params.getExtFieldScale() << ")" << endl;
            break;
        }
        default:
            ibsimu.message(1) << "\tUnknown external field case: " << params.getExtFieldCase() << endl;
            break;
    }
    
    logfile << "\tExternal field loaded with scaling " << params.getExtFieldScale() << endl << flush;
}
