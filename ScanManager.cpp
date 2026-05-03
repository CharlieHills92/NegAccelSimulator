#include "ScanManager.h"
#include <sys/stat.h>

ScanManager::ScanManager() : scanLength(1) {
}

ScanManager::~ScanManager() {
}

std::unordered_map<std::string, std::vector<double>> ScanManager::readValuesFromFile(const std::string& filename) {
    std::unordered_map<std::string, std::vector<double>> data;
    std::ifstream file(filename);

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            // Skip empty lines and comments
            if (line.empty() || line[0] == '#' || line[0] == '/') {
                continue;
            }
            
            std::istringstream iss(line);
            std::string var_name, equals_sign;
            iss >> var_name >> equals_sign;

            // Skip lines that don't have the expected format
            if (equals_sign != "=") {
                continue;
            }

            double value;
            std::vector<double> values;
            while (iss >> value) {
                values.push_back(value);
            }

            if (!values.empty()) {
                data[var_name] = values;
                std::cout << "Loaded parameter " << var_name << " with " << values.size() << " values: ";
                for (double v : values) {
                    std::cout << v << " ";
                }
                std::cout << std::endl;
            }
        }
        file.close();
        std::cout << "File " << filename << " read and data loaded!" << std::endl;
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }

    return data;
}

double ScanManager::getParameterValue(const std::string& paramName, int caseIndex) {
    // First try to find the parameter directly as it appears in the .scn file
    if (scanData.find(paramName) != scanData.end()) {
        const std::vector<double>& values = scanData[paramName];
        if (values.size() <= caseIndex) {
            return values[0]; // Use first value if not enough values
        } else {
            return values[caseIndex];
        }
    }
    
    // If not found, try with parameter name mappings for backwards compatibility
    std::string actualParamName = paramName;
    if (paramName == "J_E") actualParamName = "ELECTRONS";
    else if (paramName == "EG_V") actualParamName = "EG_VOLTAGE";
    else if (paramName == "AG1_V") actualParamName = "GG_VOLTAGE";
    else if (paramName == "AG2_V") actualParamName = "REP_VOLTAGE";
    else if (paramName == "AG3_V") actualParamName = "G3_VOLTAGE";
    else if (paramName == "AG4_V") actualParamName = "G4_VOLTAGE";
    else if (paramName == "AG5_V") actualParamName = "G5_VOLTAGE";
    else if (paramName == "T_PERP") actualParamName = "TPERP";
    else if (paramName == "T_PAR") actualParamName = "TPAR";
    
    if (actualParamName != paramName && scanData.find(actualParamName) != scanData.end()) {
        const std::vector<double>& values = scanData[actualParamName];
        if (values.size() <= caseIndex) {
            return values[0]; // Use first value if not enough values
        } else {
            return values[caseIndex];
        }
    }
    
    // If still not found, return default value
    auto it = defaultValues.find(paramName);
    if (it != defaultValues.end()) {
        return it->second;
    } else {
        std::cerr << "Warning: No default value for parameter " << paramName << std::endl;
        return 0.0;
    }
}

bool ScanManager::loadScanFile(const std::string& scanFile) {
    std::cout << "Loading scan file: " << scanFile << std::endl;
    
    // Extract file tag from filename (C++11 compatible)
    size_t lastSlash = scanFile.find_last_of("/\\");
    size_t lastDot = scanFile.find_last_of('.');
    if (lastSlash == std::string::npos) lastSlash = 0;
    else lastSlash++;
    
    if (lastDot != std::string::npos && lastDot > lastSlash) {
        scanFileTag = scanFile.substr(lastSlash, lastDot - lastSlash);
    } else {
        scanFileTag = scanFile.substr(lastSlash);
    }
    
    // Read the scan data
    scanData = readValuesFromFile(scanFile);
    
    // Determine scan length
    if (scanData.find("SCAN_LENGTH") != scanData.end()) {
        scanLength = static_cast<int>(scanData["SCAN_LENGTH"][0]);
    } else {
        // Find maximum number of values for any parameter
        scanLength = 1;
        for (const auto& pair : scanData) {
            if (pair.second.size() > scanLength) {
                scanLength = pair.second.size();
            }
        }
    }
    
    std::cout << "Scan length determined: " << scanLength << std::endl;
    
    // Create scan cases
    scanCases.clear();
    for (int i = 0; i < scanLength; i++) {
        ScanCase scanCase;
        
        scanCase.accelerator_idx = static_cast<int>(getParameterValue("ACCELERATOR_IDX", i));
        scanCase.b_ison = static_cast<int>(getParameterValue("B_ISON", i));
        scanCase.stripping = static_cast<int>(getParameterValue("STRIPPING", i));
        scanCase.include_surface_collisions = static_cast<int>(getParameterValue("INCLUDE_SURFACE_COLLISIONS", i));
        scanCase.j_e = getParameterValue("J_E", i);
        scanCase.m_ions = getParameterValue("M_IONS", i);
        scanCase.q_ions = static_cast<int>(getParameterValue("Q_IONS", i));
        scanCase.j_ion = getParameterValue("J_ION", i);
        scanCase.t_perp = getParameterValue("T_PERP", i);
        scanCase.t_par = getParameterValue("T_PAR", i);
        scanCase.e0_z = getParameterValue("E0_Z", i);
        scanCase.u_plasma = getParameterValue("U_PLASMA", i);
        scanCase.n_particles = static_cast<int>(getParameterValue("N_PARTICLES", i));
        scanCase.eg_v = getParameterValue("EG_V", i);
        scanCase.ag1_v = getParameterValue("AG1_V", i);
        scanCase.ag2_v = getParameterValue("AG2_V", i);
        scanCase.ag3_v = getParameterValue("AG3_V", i);
        scanCase.ag4_v = getParameterValue("AG4_V", i);
        scanCase.ag5_v = getParameterValue("AG5_V", i);
        scanCase.mesh_size = getParameterValue("MESH_SIZE", i);
        scanCase.iterations = static_cast<int>(getParameterValue("ITERATIONS", i));
        scanCase.pgfilter_scale = getParameterValue("PGFILTER_SCALE", i);
        scanCase.cesmadcm_scale = getParameterValue("CESMADCM_SCALE", i);
        scanCase.extfield_case = static_cast<int>(getParameterValue("EXTFIELD_CASE", i));
        scanCase.extfield_scale = getParameterValue("EXTFIELD_SCALE", i);
        scanCase.split_domain = static_cast<int>(getParameterValue("SPLIT_DOMAIN", i));
        scanCase.jtolerance = getParameterValue("JTOLERANCE", i);
        scanCase.alpha_coeff = getParameterValue("ALPHA_COEFF", i);
        scanCase.t_positive = getParameterValue("T_POSITIVE", i);
        scanCase.ext_gap = getParameterValue("EXT_GAP", i);
        scanCase.acc_gap = getParameterValue("ACC_GAP", i);
        scanCase.domain_x_size = getParameterValue("DOMAIN_X_SIZE", i);
        scanCase.domain_y_size = getParameterValue("DOMAIN_Y_SIZE", i);
        scanCase.domain_z_size = getParameterValue("DOMAIN_Z_SIZE", i);
        scanCase.mgsolver = static_cast<uint>(getParameterValue("MGSOLVER", i));
        scanCase.shield_model = static_cast<uint>(getParameterValue("SHIELD_MODEL", i));

        scanCase.input_file_tag = scanFileTag + "_" + std::to_string(i);
        
        scanCases.push_back(scanCase);
        
        std::cout << "Created scan case " << i << " with STRIPPING=" << scanCase.stripping << std::endl;
    }
    
    return !scanCases.empty();
}

std::vector<std::string> ScanManager::createInputFileFromScan() {
    std::vector<std::string> inputFileTags;
    
    // Create output directory if it doesn't exist (C++11 compatible)
    if (mkdir(scanFileTag.c_str(), 0755) == -1) {
        // Directory might already exist, that's okay
    }
    
    for (int i = 0; i < scanLength; i++) {
        const ScanCase& scanCase = scanCases[i];
        
        std::string outfile = scanFileTag + "/" + scanCase.input_file_tag + ".inp";
        std::cout << "Creating file... " << outfile << "... ";
        
        std::ofstream pout(outfile);
        if (!pout.is_open()) {
            std::cerr << "Failed to create file " << outfile << std::endl;
            continue;
        }
        
        pout << scanCase.accelerator_idx << "\t// ACCELERATOR INDEX (MITICA=2)" << std::endl
             << scanCase.b_ison << "\t// B_ISON (TRUE=1)" << std::endl
             << scanCase.stripping << "\t// STRIPPING (TRUE=1)" << std::endl
             << scanCase.include_surface_collisions << "\t// INCLUDE_SURFACE_COLLISIONS (TRUE=1)" << std::endl
             << scanCase.j_e << "\t// J ELECTRONS" << std::endl
             << scanCase.m_ions << "\t// M_IONS" << std::endl
             << scanCase.q_ions << "\t// Q_IONS" << std::endl
             << scanCase.j_ion << "\t// J_IONS" << std::endl
             << scanCase.t_perp << "\t// T_PERP" << std::endl
             << scanCase.t_par << "\t// T_PAR" << std::endl
             << scanCase.e0_z << "\t// BEAM STARTING AXIAL ENERGY" << std::endl
             << scanCase.u_plasma << "\t// PLASMA POTENTIAL" << std::endl
             << scanCase.n_particles << "\t// NUMBER OF PARTICLES" << std::endl
             << scanCase.eg_v << "\t// EG VOLTAGE" << std::endl
             << scanCase.ag1_v << "\t// AG1 VOLTAGE" << std::endl
             << scanCase.ag2_v << "\t// AG2 VOLTAGE" << std::endl
             << scanCase.ag3_v << "\t// AG3 VOLTAGE" << std::endl
             << scanCase.ag4_v << "\t// AG4 VOLTAGE" << std::endl
             << scanCase.ag5_v << "\t// AG5 VOLTAGE" << std::endl
             << scanCase.mesh_size << "\t// MESH SIZE" << std::endl
             << scanCase.iterations << "\t// ITERATIONS" << std::endl
             << scanCase.pgfilter_scale << "\t// PG FILTER FIELD SCALE" << std::endl
             << scanCase.cesmadcm_scale << "\t// CESM+ADCM FIELD SCALE" << std::endl
             << scanCase.extfield_case << "\t// EXTERNAL FIELD CASE" << std::endl
             << scanCase.extfield_scale << "\t// EXTERNAL FIELD SCALE" << std::endl
             << scanCase.split_domain << "\t// SPLIT DOMAIN (TRUE=1)" << std::endl
             << scanCase.jtolerance << "\t// TOLERANCE ON EXTRACTED CURRENT" << std::endl
             << scanCase.alpha_coeff << "\t// ALPHA COEFFICIENT FOR SC AVERAGING" << std::endl
             << scanCase.t_positive << "\t// POSITIVE ION TEMPERATURE" << std::endl
             << scanCase.ext_gap << "\t// EXTRACTION GAP LENGTH" << std::endl
             << scanCase.acc_gap << "\t// ACCELERATION GAP LENGTH" << std::endl
             << scanCase.domain_x_size << "\t// DOMAIN X SIZE [mm]" << std::endl
             << scanCase.domain_y_size << "\t// DOMAIN Y SIZE [mm]" << std::endl
             << scanCase.domain_z_size << "\t// DOMAIN Z SIZE [mm]" << std::endl
             << scanCase.mgsolver << "\t// MGSOLVER" << std::endl
             << scanCase.shield_model << "\t// MENISCUS MODEL" << std::endl;

        pout.close();
        inputFileTags.push_back(scanCase.input_file_tag);
        std::cout << "Done!" << std::endl;
    }
    
    std::cout << "Input files created!" << std::endl;
    return inputFileTags;
}
