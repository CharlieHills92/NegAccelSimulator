#ifndef SCANMANAGER_H
#define SCANMANAGER_H

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <iostream>

class ScanManager {
public:
    struct ScanCase {
        int accelerator_idx;
        int b_ison;
        int stripping;
        int include_surface_collisions;
        double j_e;
        double m_ions;
        int q_ions;
        double j_ion;
        double t_perp;
        double t_par;
        double e0_z;
        double u_plasma;
        int n_particles;
        double eg_v;
        double ag1_v;
        double ag2_v;
        double ag3_v;
        double ag4_v;
        double ag5_v;
        double mesh_size;
        int iterations;
        double pgfilter_scale;
        double cesmadcm_scale;
        int extfield_case;
        double extfield_scale;
        int split_domain;
        double jtolerance;
        double alpha_coeff;
        double t_positive;
        double ext_gap;
        double acc_gap;
        double domain_x_size;
        double domain_y_size;
        double domain_z_size;
        uint mgsolver;
        uint shield_model;
        
        std::string input_file_tag;
    };

private:
    std::unordered_map<std::string, std::vector<double>> scanData;
    std::vector<ScanCase> scanCases;
    std::string scanFileTag;
    int scanLength;

    // Default values
    const std::unordered_map<std::string, double> defaultValues = {
        {"ACCELERATOR_IDX", 2},
        {"B_ISON", 1},
        {"STRIPPING", 0},
        {"INCLUDE_SURFACE_COLLISIONS", 0},
        {"J_E", 0},
        {"M_IONS", 1},
        {"Q_IONS", -1},
        {"J_ION", 330},
        {"T_PERP", 1},
        {"T_PAR", 0},
        {"E0_Z", 4},
        {"U_PLASMA", 3},
        {"N_PARTICLES", 120000},
        {"EG_V", 8e3},
        {"AG1_V", 182e3},
        {"AG2_V", 356e3},
        {"AG3_V", 530e3},
        {"AG4_V", 704e3},
        {"AG5_V", 878e3},
        {"MESH_SIZE", 0.5e-3},
        {"ITERATIONS", 15},
        {"PGFILTER_SCALE", 0},
        {"CESMADCM_SCALE", 1},
        {"EXTFIELD_CASE", 0},
        {"EXTFIELD_SCALE", 0},
        {"SPLIT_DOMAIN", 1},
        {"JTOLERANCE", 0.1},
        {"ALPHA_COEFF", 0.3},
        {"T_POSITIVE", 0.8},
        {"EXT_GAP", 0.006},
        {"ACC_GAP", 0.088},
        {"DOMAIN_X_SIZE", 30},
        {"DOMAIN_Y_SIZE", 30},
        {"DOMAIN_Z_SIZE", 567},
        {"MGSOLVER", 0U},
        {"MENISCUS_MODEL", 0U},
        {"SHIELD_MODEL", 0U}
    };

    std::unordered_map<std::string, std::vector<double>> readValuesFromFile(const std::string& filename);
    double getParameterValue(const std::string& paramName, int caseIndex);
    void createInputFiles();

public:
    ScanManager();
    ~ScanManager();

    bool loadScanFile(const std::string& scanFile);
    std::vector<ScanCase> getScanCases() const { return scanCases; }
    int getScanLength() const { return scanLength; }
    std::string getScanFileTag() const { return scanFileTag; }
    
    // Create individual .inp files for each scan case
    std::vector<std::string> createInputFileFromScan();
};

#endif // SCANMANAGER_H
