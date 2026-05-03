/*
 * SimulationParameters.cpp
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored from ManageSimulation)
 */

#include "SimulationParameters.h"
#include "error.hpp"
#include "ibsimu.hpp"

#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;

SimulationParameters::SimulationParameters() :
    ACCELERATOR_IDX(0), B_ISON(0), INCLUDE_STRIPPING(0), INCLUDE_SURFACE_COLLISIONS(0), ELECTRONS(0),
    M_IONS(0), Q_IONS(0), J_ION(0), TPERP(0), TPAR(0), E0_Z(0), U_PLASMA(0), N_PARTICLES(0),
    EG_VOLTAGE(0), GG_VOLTAGE(0), REP_VOLTAGE(0), G1_VOLTAGE(0), G2_VOLTAGE(0), 
    G3_VOLTAGE(0), G4_VOLTAGE(0), G5_VOLTAGE(0),
    MESH_SIZE(0), ITERATIONS(0), PGFILTER_SCALE(0), CESMADCM_SCALE(0), MGSOLVER(0),
    EXTFIELD_CASE(0), EXTFIELD_SCALE(0), SPLIT_DOMAIN(0), JTOLERANCE(0), 
    ALPHA_COEFF(0), T_POSITIVE(0), SHIELD_MODEL(0), N_SOLIDS(0),
    EXT_GAP(0), ACC_GAP(0), DOMAIN_X_SIZE(-1), DOMAIN_Y_SIZE(-1), DOMAIN_Z_SIZE(-1),
    EGEXTJ(0), domain_ii(0) {
}

void SimulationParameters::readParametersFromFile(const string& input) {
    // Check if input is a scenario file (.scn) or input file (.inp)
    size_t ext_pos = input.find_last_of('.');
    string extension = (ext_pos != string::npos) ? input.substr(ext_pos) : "";
    
    if (extension == ".scn") {
        // Parse scenario file directly
        parseScenarioFile(input);
        
        cout << "Using scenario-based configuration from: " << input << endl;
    } else {
        // Parse legacy .inp file
        parseInputFile(input);
        
        cout << "Using legacy input file configuration from: " << input << endl;
    }
}

void SimulationParameters::parseInputFile(const string& inputFile) {
    cout << "Parsing legacy input file: " << inputFile << endl;
    
    // First set all default values
    setDefaultValues();
    
    ifstream file(inputFile);
    if (!file.is_open()) {
        throw Error(ERROR_LOCATION, "Could not open input file: " + inputFile);
    }

    string line;
    int line_number = 0;
    
    while (getline(file, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        line_number++;
        
        // Extract value from "VALUE // DESCRIPTION" format
        string value;
        size_t comment_pos = line.find("//");
        if (comment_pos != string::npos) {
            value = line.substr(0, comment_pos);
        } else {
            value = line;
        }
        
        // Trim whitespace
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);
        
        if (value.empty()) continue;
        
        // Handle special KEY=VALUE format for domain sizes
        if (value.find("DOMAIN_") != string::npos && value.find("=") != string::npos) {
            size_t equals_pos = value.find('=');
            string key = value.substr(0, equals_pos);
            string val = value.substr(equals_pos + 1);
            
            // Trim whitespace
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            val.erase(0, val.find_first_not_of(" \t"));
            val.erase(val.find_last_not_of(" \t") + 1);
            
            if (key == "DOMAIN_X_SIZE") {
                DOMAIN_X_SIZE = stod(val) * 1e-3; // Convert mm to m
                cout << "Set " << key << " = " << val << "mm (" << DOMAIN_X_SIZE << "m)" << endl;
            } else if (key == "DOMAIN_Y_SIZE") {
                DOMAIN_Y_SIZE = stod(val) * 1e-3; // Convert mm to m
                cout << "Set " << key << " = " << val << "mm (" << DOMAIN_Y_SIZE << "m)" << endl;
            } else if (key == "DOMAIN_Z_SIZE") {
                DOMAIN_Z_SIZE = stod(val) * 1e-3; // Convert mm to m
                cout << "Set " << key << " = " << val << "mm (" << DOMAIN_Z_SIZE << "m)" << endl;
            }
            continue;
        }
        
        // Map line position to parameter assignment (based on .inp file format)
        try {
            if (line_number == 1) {
                ACCELERATOR_IDX = static_cast<uint>(stod(value));
                cout << "Line " << line_number << ": ACCELERATOR_IDX = " << value << endl;
            } else if (line_number == 2) {
                B_ISON = static_cast<uint>(stod(value));
                cout << "Line " << line_number << ": B_ISON = " << value << endl;
            } else if (line_number == 3) {
                INCLUDE_STRIPPING = static_cast<uint>(stod(value));
                cout << "Line " << line_number << ": INCLUDE_STRIPPING = " << value << endl;
            } else if (line_number == 4) {
                INCLUDE_SURFACE_COLLISIONS = static_cast<uint>(stod(value));
                cout << "Line " << line_number << ": INCLUDE_SURFACE_COLLISIONS = " << value << endl;
            } else if (line_number == 5) {
                ELECTRONS = stod(value);
                cout << "Line " << line_number << ": ELECTRONS = " << value << endl;
            } else if (line_number == 6) {
                M_IONS = stod(value);
                cout << "Line " << line_number << ": M_IONS = " << value << endl;
            } else if (line_number == 7) {
                Q_IONS = stod(value);
                cout << "Line " << line_number << ": Q_IONS = " << value << endl;
            } else if (line_number == 8) {
                J_ION = stod(value);
                cout << "Line " << line_number << ": J_ION = " << value << endl;
            } else if (line_number == 9) {
                TPERP = stod(value);
                cout << "Line " << line_number << ": TPERP = " << value << endl;
            } else if (line_number == 10) {
                TPAR = stod(value);
                cout << "Line " << line_number << ": TPAR = " << value << endl;
            } else if (line_number == 11) {
                E0_Z = stod(value);
                cout << "Line " << line_number << ": E0_Z = " << value << endl;
            } else if (line_number == 12) {
                U_PLASMA = stod(value);
                cout << "Line " << line_number << ": U_PLASMA = " << value << endl;
            } else if (line_number == 13) {
                N_PARTICLES = static_cast<uint>(stod(value));
                cout << "Line " << line_number << ": N_PARTICLES = " << value << endl;
            } else if (line_number == 14) {
                EG_VOLTAGE = stod(value);
                cout << "Line " << line_number << ": EG_VOLTAGE = " << value << endl;
            } else if (line_number == 15) {
                GG_VOLTAGE = stod(value);
                cout << "Line " << line_number << ": GG_VOLTAGE = " << value << endl;
            } else if (line_number == 16) {
                REP_VOLTAGE = stod(value);
                cout << "Line " << line_number << ": REP_VOLTAGE = " << value << endl;
            } else if (line_number == 17) {
                G3_VOLTAGE = stod(value);
                cout << "Line " << line_number << ": G3_VOLTAGE = " << value << endl;
            } else if (line_number == 18) {
                G4_VOLTAGE = stod(value);
                cout << "Line " << line_number << ": G4_VOLTAGE = " << value << endl;
            } else if (line_number == 19) {
                G5_VOLTAGE = stod(value);
                cout << "Line " << line_number << ": G5_VOLTAGE = " << value << endl;
            } else if (line_number == 20) {
                MESH_SIZE = stod(value);
                cout << "Line " << line_number << ": MESH_SIZE = " << value << endl;
            } else if (line_number == 21) {
                ITERATIONS = static_cast<uint>(stod(value));
                cout << "Line " << line_number << ": ITERATIONS = " << value << endl;
            } else if (line_number == 22) {
                PGFILTER_SCALE = stod(value);
                cout << "Line " << line_number << ": PGFILTER_SCALE = " << value << endl;
            } else if (line_number == 23) {
                CESMADCM_SCALE = stod(value);
                cout << "Line " << line_number << ": CESMADCM_SCALE = " << value << endl;
            } else if (line_number == 24) {
                EXTFIELD_CASE = static_cast<uint>(stod(value));
                cout << "Line " << line_number << ": EXTFIELD_CASE = " << value << endl;
            } else if (line_number == 25) {
                EXTFIELD_SCALE = stod(value);
                cout << "Line " << line_number << ": EXTFIELD_SCALE = " << value << endl;
            } else if (line_number == 26) {
                SPLIT_DOMAIN = static_cast<uint>(stod(value));
                cout << "Line " << line_number << ": SPLIT_DOMAIN = " << value << endl;
            } else if (line_number == 27) {
                JTOLERANCE = stod(value);
                cout << "Line " << line_number << ": JTOLERANCE = " << value << endl;
            } else if (line_number == 28) {
                ALPHA_COEFF = stod(value);
                cout << "Line " << line_number << ": ALPHA_COEFF = " << value << endl;
            } else if (line_number == 29) {
                T_POSITIVE = stod(value);
                cout << "Line " << line_number << ": T_POSITIVE = " << value << endl;
            } else if (line_number == 30) {
                EXT_GAP = stod(value);
                cout << "Line " << line_number << ": EXT_GAP = " << value << endl;
            } else if (line_number == 31) {
                ACC_GAP = stod(value);
                cout << "Line " << line_number << ": ACC_GAP = " << value << endl;
            } else if (line_number == 32) {
                DOMAIN_X_SIZE = stod(value) * 1e-3; // Convert mm to m
                cout << "Line " << line_number << ": DOMAIN_X_SIZE = " << value << "mm (" << DOMAIN_X_SIZE << "m)" << endl;
            } else if (line_number == 33) {
                DOMAIN_Y_SIZE = stod(value) * 1e-3; // Convert mm to m
                cout << "Line " << line_number << ": DOMAIN_Y_SIZE = " << value << "mm (" << DOMAIN_Y_SIZE << "m)" << endl;
            } else if (line_number == 34) {
                DOMAIN_Z_SIZE = stod(value) * 1e-3; // Convert mm to m
                cout << "Line " << line_number << ": DOMAIN_Z_SIZE = " << value << "mm (" << DOMAIN_Z_SIZE << "m)" << endl;
            } else if (line_number == 35) {
                MGSOLVER = static_cast<uint>(stod(value));
                cout << "Line " << line_number << ": MGSOLVER = " << value << endl;
            } else if (line_number == 36) {
                SHIELD_MODEL = static_cast<uint>(stod(value));
                cout << "Line " << line_number << ": SHIELD_MODEL = " << value << endl;
            } else if (line_number > 36) {
                // Handle additional lines beyond the standard 36 parameters
                // These might include domain size specifications or other extensions
                cout << "Line " << line_number << ": " << line << " (additional parameter)" << endl;
            }
        } catch (const exception& e) {
            cout << "Warning: Could not parse line " << line_number << " value '" << value << "': " << e.what() << endl;
        }
    }
    
    file.close();
    
    // Set dependent variables
    G1_VOLTAGE = GG_VOLTAGE;
    G2_VOLTAGE = REP_VOLTAGE;
    
    cout << "Input file parsed successfully." << endl;
}

double SimulationParameters::getDomainXSizeOrDefault() const {
    if (DOMAIN_X_SIZE > 0) {
        return DOMAIN_X_SIZE;
    }
    
    // Return accelerator-specific defaults
    switch (static_cast<uint>(ACCELERATOR_IDX)) {
        case 1: // SPIDER
        case 2: // MITICA  
        case 3: // MTF
        default:
            return 30.0e-3; // 30mm default for all accelerators
    }
}

double SimulationParameters::getDomainYSizeOrDefault() const {
    if (DOMAIN_Y_SIZE > 0) {
        return DOMAIN_Y_SIZE;
    }
    
    // Return accelerator-specific defaults
    switch (static_cast<uint>(ACCELERATOR_IDX)) {
        case 1: // SPIDER
        case 2: // MITICA
        case 3: // MTF
        default:
            return 30.0e-3; // 30mm default for all accelerators
    }
}

double SimulationParameters::getDomainZSizeOrDefault() const {
    if (DOMAIN_Z_SIZE > 0) {
        return DOMAIN_Z_SIZE;
    }
    
    // Return accelerator-specific defaults
    switch (static_cast<uint>(ACCELERATOR_IDX)) {
        case 1: // SPIDER
            return 80.0e-3;  // 80mm
        case 2: // MITICA
            return 567.0e-3; // 567mm
        case 3: // MTF
        default:
            return 567.0e-3; // 567mm default
    }
}

void SimulationParameters::parseScenarioFile(const string& scenarioFile) {
    cout << "Parsing scenario file: " << scenarioFile << endl;
    
    // First set all default values
    setDefaultValues();
    
    ifstream file(scenarioFile);
    if (!file.is_open()) {
        throw Error(ERROR_LOCATION, "Could not open scenario file: " + scenarioFile);
    }
    
    string line;
    while (getline(file, line)) {
        // Remove comments and trim whitespace
        size_t comment_pos = line.find('#');
        if (comment_pos != string::npos) {
            line = line.substr(0, comment_pos);
        }
        
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" 	"));
        line.erase(line.find_last_not_of(" 	") + 1);
        
        if (line.empty()) continue;
        
        // Parse key-value pairs (handle both "KEY=VALUE" and "KEY VALUE" formats)
        string key, value;
        size_t equals_pos = line.find('=');
        
        if (equals_pos != string::npos) {
            // Handle "KEY=VALUE" format
            key = line.substr(0, equals_pos);
            value = line.substr(equals_pos + 1);
        } else {
            // Handle "KEY VALUE" format (typical for .scn files)
            istringstream iss(line);
            iss >> key >> value;
        }
        
        // Trim whitespace from key and value
        key.erase(0, key.find_first_not_of(" 	"));
        key.erase(key.find_last_not_of(" 	") + 1);
        value.erase(0, value.find_first_not_of(" 	"));
        value.erase(value.find_last_not_of(" 	") + 1);
        
        if (key.empty() || value.empty()) continue;
        
        // Parse user-defined parameters
        try {
            if (key == "ACCELERATOR_IDX") {
                ACCELERATOR_IDX = static_cast<uint>(stod(value));
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "B_ISON") {
                B_ISON = static_cast<uint>(stod(value));
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "INCLUDE_STRIPPING" || key == "STRIPPING") {
                INCLUDE_STRIPPING = static_cast<uint>(stod(value));
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "INCLUDE_SURFACE_COLLISIONS" || key == "SURFACE_COLLISIONS" || key == "SURFACE") {
                INCLUDE_SURFACE_COLLISIONS = static_cast<uint>(stod(value));
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "ELECTRONS") {
                ELECTRONS = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "M_IONS") {
                M_IONS = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "Q_IONS") {
                Q_IONS = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "J" || key == "J_ION") {
                J_ION = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "N_PARTICLES") {
                N_PARTICLES = static_cast<uint>(stod(value));
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "MESH_SIZE") {
                MESH_SIZE = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "ITERATIONS") {
                ITERATIONS = static_cast<uint>(stod(value));
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "JTOLERANCE") {
                JTOLERANCE = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "PGFILTER_SCALE") {
                PGFILTER_SCALE = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "CESMADCM_SCALE") {
                CESMADCM_SCALE = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "EXTFIELD_SCALE") {
                EXTFIELD_SCALE = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "EXTFIELD_CASE") {
                EXTFIELD_CASE = static_cast<uint>(stod(value));
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "T_PERP") {
                TPERP = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "E0_Z") {
                E0_Z = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "SPLIT_DOMAIN") {
                SPLIT_DOMAIN = static_cast<uint>(stod(value));
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "ACC_GAP") {
                ACC_GAP = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "EG_V") {
                EG_VOLTAGE = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "AG1_V") {
                GG_VOLTAGE = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "AG2_V") {
                REP_VOLTAGE = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "AG3_V") {
                G3_VOLTAGE = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "AG4_V") {
                G4_VOLTAGE = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "AG5_V") {
                G5_VOLTAGE = stod(value);
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "DOMAIN_X_SIZE") {
                DOMAIN_X_SIZE = stod(value) * 1e-3; // Convert mm to m
                cout << "Set " << key << " = " << value << "mm (" << DOMAIN_X_SIZE << "m)" << endl;
            } else if (key == "DOMAIN_Y_SIZE") {
                DOMAIN_Y_SIZE = stod(value) * 1e-3; // Convert mm to m
                cout << "Set " << key << " = " << value << "mm (" << DOMAIN_Y_SIZE << "m)" << endl;
            } else if (key == "DOMAIN_Z_SIZE") {
                DOMAIN_Z_SIZE = stod(value) * 1e-3; // Convert mm to m
                cout << "Set " << key << " = " << value << "mm (" << DOMAIN_Z_SIZE << "m)" << endl;
            } else if (key == "MGSOLVER") {
                MGSOLVER = static_cast<uint>(stod(value));
                cout << "Set " << key << " = " << value << endl;
            } else if (key == "SHIELD_MODEL") {
                SHIELD_MODEL = static_cast<uint>(stod(value));
                cout << "Set " << key << " = " << value << endl;
            }
            // Add more parameters as needed
        } catch (const exception& e) {
            cout << "Warning: Could not parse parameter " << key << " = " << value << endl;
        }
    }
    
    file.close();
    
    // Set dependent variables
    G1_VOLTAGE = GG_VOLTAGE;
    G2_VOLTAGE = REP_VOLTAGE;
    
    cout << "Scenario file parsed successfully." << endl;
}

void SimulationParameters::setDefaultValues() {
    // Set defaults based on accelerator type
    uint accel_type = static_cast<uint>(ACCELERATOR_IDX);
    
    // Default values common to all accelerators
    B_ISON = 1;
    INCLUDE_STRIPPING = 0;
    INCLUDE_SURFACE_COLLISIONS = 0;  // Disabled by default
    ELECTRONS = 0;
    M_IONS = 1.0;
    Q_IONS = -1.0;
    J_ION = 330.0;
    TPERP = 0.0;
    TPAR = 0.0;
    E0_Z = 3.0;
    U_PLASMA = 3.0;
    N_PARTICLES = 150000;
    
    // Default voltages for MTF
    if (accel_type == 3) {
        EG_VOLTAGE = 8000.0;
        GG_VOLTAGE = 182000.0;
        REP_VOLTAGE = 356000.0;
        G3_VOLTAGE = 530000.0;
        G4_VOLTAGE = 704000.0;
        G5_VOLTAGE = 878000.0;
    } else {
        // Default voltages for other accelerators
        EG_VOLTAGE = 8000.0;
        GG_VOLTAGE = 182000.0;
        REP_VOLTAGE = 356000.0;
        G3_VOLTAGE = 530000.0;
        G4_VOLTAGE = 704000.0;
        G5_VOLTAGE = 878000.0;
    }
    
    // Numerical defaults
    MESH_SIZE = 0.0003;
    ITERATIONS = 5;
    PGFILTER_SCALE = 0.0;
    CESMADCM_SCALE = 0.0;
    EXTFIELD_CASE = 0;
    EXTFIELD_SCALE = 1.0;
    SPLIT_DOMAIN = 0;
    JTOLERANCE = 1.0;
    ALPHA_COEFF = 0.3;
    T_POSITIVE = 0.8;
    MGSOLVER = 0; // Default to no multigrid solver
    SHIELD_MODEL = 0; // Default to no shield model (use nsimp)

    // Geometry defaults
    EXT_GAP = 0.006;
    ACC_GAP = 0.088;
    
    // Domain size defaults (if not already set)
    if (DOMAIN_X_SIZE < 0) DOMAIN_X_SIZE = 30.0e-3; // 30mm
    if (DOMAIN_Y_SIZE < 0) DOMAIN_Y_SIZE = 30.0e-3; // 30mm
    if (DOMAIN_Z_SIZE < 0) {
        switch (accel_type) {
            case 1: DOMAIN_Z_SIZE = 80.0e-3; break;   // SPIDER: 80mm
            case 2: DOMAIN_Z_SIZE = 567.0e-3; break;  // MITICA: 567mm
            case 3: DOMAIN_Z_SIZE = 567.0e-3; break;  // MTF: 567mm
            default: DOMAIN_Z_SIZE = 567.0e-3; break;
        }
    }
    
    // Other defaults
    EGEXTJ = 0.0;
    domain_ii = 0;
    N_SOLIDS = 0;
}

void SimulationParameters::generateInputFile(const string& inputFile) {
    cout << "Generating complete input file: " << inputFile << endl;
    
    ofstream file(inputFile);
    if (!file.is_open()) {
        throw Error(ERROR_LOCATION, "Could not create input file: " + inputFile);
    }
    
    // Write all parameters in the expected order
    file << ACCELERATOR_IDX << "\t// ACCELERATOR_IDX" << endl;
    file << B_ISON << "\t// B_ISON" << endl;
    file << INCLUDE_STRIPPING << "\t// INCLUDE_STRIPPING" << endl;
    file << INCLUDE_SURFACE_COLLISIONS << "\t// INCLUDE_SURFACE_COLLISIONS" << endl;
    file << ELECTRONS << "\t// ELECTRONS" << endl;
    file << M_IONS << "\t// M_IONS [amu]" << endl;
    file << Q_IONS << "\t// Q_IONS [e]" << endl;
    file << J_ION << "\t// J [A/m^2]" << endl;
    file << TPERP << "\t// TPERP [eV]" << endl;
    file << TPAR << "\t// TPAR [eV]" << endl;
    file << E0_Z << "\t// AXIAL ENERGY [eV]" << endl;
    file << U_PLASMA << "\t// PLASMA POTENTIAL [V]" << endl;
    file << N_PARTICLES << "\t// N_PARTICLES" << endl;
    file << EG_VOLTAGE << "\t// EG VOLTAGE [V]" << endl;
    file << GG_VOLTAGE << "\t// GG/G1 VOLTAGE [V]" << endl;
    file << REP_VOLTAGE << "\t// REP/G2 VOLTAGE [V]" << endl;
    file << G3_VOLTAGE << "\t// G3 VOLTAGE [V]" << endl;
    file << G4_VOLTAGE << "\t// G4 VOLTAGE [V]" << endl;
    file << G5_VOLTAGE << "\t// G5 VOLTAGE [V]" << endl;
    file << MESH_SIZE << "\t// MESH SIZE [m]" << endl;
    file << ITERATIONS << "\t// ITERATIONS" << endl;
    file << PGFILTER_SCALE << "\t// PG FILTER SCALE" << endl;
    file << CESMADCM_SCALE << "\t// CESM+ADCM SCALE" << endl;
    file << EXTFIELD_CASE << "\t// EXT FIELD CASE" << endl;
    file << EXTFIELD_SCALE << "\t// EXT FIELD SCALE" << endl;
    file << SPLIT_DOMAIN << "\t// SPLIT DOMAIN" << endl;
    file << JTOLERANCE << "\t// J TOLERANCE" << endl;
    file << ALPHA_COEFF << "\t// SC AVERAGING COEFFICIENT" << endl;
    file << T_POSITIVE << "\t// POSITIVE ION TEMPERATURE [eV]" << endl;
    file << EXT_GAP << "\t// EXTRACTION GAP LENGTH [m]" << endl;
    file << ACC_GAP << "\t// ACCELERATOR GAP LENGTH [m]" << endl;
    file << (DOMAIN_X_SIZE * 1e3) << "\t// DOMAIN X SIZE IN mm" << endl;
    file << (DOMAIN_Y_SIZE * 1e3) << "\t// DOMAIN Y SIZE IN mm" << endl;
    file << (DOMAIN_Z_SIZE * 1e3) << "\t// DOMAIN Z SIZE IN mm" << endl;
    file << MGSOLVER << "\t// MGSOLVER (0 BICGSTAB, 1 MG)" << endl;
    file << SHIELD_MODEL << "\t// MENISCUS MODEL (0 NSIMP, 1 SHIELD)" << endl;

    file.close();
    cout << "Input file generated successfully." << endl;
}

SimulationParameters::ScanParameters SimulationParameters::parseScanFile(const std::string& scnFile) {
    ScanParameters scanParams;
    ifstream file(scnFile);
    
    if (!file.is_open()) {
        throw runtime_error("Could not open scenario file: " + scnFile);
    }
    
    string line;
    while (getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '/' || line[0] == '#') continue;
        
        // Parse key-value pairs
        string key, value_str;
        size_t equals_pos = line.find('=');
        
        if (equals_pos != string::npos) {
            key = line.substr(0, equals_pos);
            value_str = line.substr(equals_pos + 1);
        } else {
            istringstream iss(line);
            iss >> key >> value_str;
        }
        
        // Trim whitespace
        key.erase(0, key.find_first_not_of(" 	"));
        key.erase(key.find_last_not_of(" 	") + 1);
        value_str.erase(0, value_str.find_first_not_of(" 	"));
        value_str.erase(value_str.find_last_not_of(" 	") + 1);
        
        if (key.empty()) continue;
        
        // Handle scan length
        if (key == "SCAN_LENGTH") {
            scanParams.scan_length = stoi(value_str);
            continue;
        }
        
        // // Parse multiple values for scan parameters
        // if (key == "STRIPPING" || key == "INCLUDE_STRIPPING") {
        //     scanParams.scan_parameter_name = "STRIPPING";
        //     istringstream values_stream(value_str);
        //     string val;
        //     while (values_stream >> val) {
        //         scanParams.stripping_values.push_back(stod(val));
        //     }
        // }
        // Add other scan parameters as needed (voltage scans, current scans, etc.)
    }
    
    file.close();
    
    // // If no scan parameters found but scan_length > 1, create default values
    // if (scanParams.stripping_values.empty() && scanParams.scan_length > 1) {
    //     // Default stripping scan: 0 and 1
    //     scanParams.stripping_values = {0, 1};
    //     scanParams.scan_parameter_name = "STRIPPING";
    // }
    
    return scanParams;
}

SimulationParameters SimulationParameters::createScanCase(const ScanParameters& scanParams, int caseIndex) {
    SimulationParameters caseParams = *this; // Copy current parameters
    
    // // Modify the varying parameter for this case
    // if (scanParams.scan_parameter_name == "STRIPPING" && 
    //     caseIndex < static_cast<int>(scanParams.stripping_values.size())) {
    //     caseParams.setIncludeStripping(static_cast<uint>(scanParams.stripping_values[caseIndex]));
    // }
    
    return caseParams;
}
