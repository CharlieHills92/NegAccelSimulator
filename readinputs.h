#ifndef READINPUTS_H_
#define READINPUTS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>

// Function to read values from file
std::unordered_map<std::string, std::vector<double>> read_values_from_file(const std::string& filename) {
    std::unordered_map<std::string, std::vector<double>> data;
    std::ifstream file(filename);

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string var_name;
            iss >> var_name;

            double value;
            std::vector<double> values;
            while (iss >> value) {
                values.push_back(value);
            }

            data[var_name] = values;
        }
        file.close();
        std::cout << "File read and data loaded!\n";
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }

    return data;
};

std::vector<std::string> create_input_file_from_scan( std::string scan_file_tag ) {
    //std::string filename = "your_file.txt"; // Change this to your file's name
    std::string filename;
    // Check if scan_file_tag already has .scn extension
    if (scan_file_tag.size() >= 4 && scan_file_tag.substr(scan_file_tag.size() - 4) == ".scn") {
        filename = scan_file_tag;
    } else {
        filename = scan_file_tag + ".scn";
    }

	std::cout << " Scan file name " << filename << "\n";
    auto data = read_values_from_file(filename);

    std::string variable_name = "SCAN_LENGTH";
    uint scan_length = 1; // Default value to prevent infinite loop
    if (data.find(variable_name) != data.end()) {
        const std::vector<double>& values = data[variable_name];
        scan_length = uint(values[0]);
    }
    else {
        std::cerr << "Variable " << variable_name << " not defined! Using default value 1" << std::endl;
        scan_length = 1; // Explicitly set to prevent infinite loop
    }

    std::vector<std::string> input_file_tags;
    
    vector<double> J_ION;
    vector<double> J_E;
    vector<double> EG_V;
    vector<double> AG1_V;
    vector<double> AG2_V;
    vector<double> AG3_V;
    vector<double> AG4_V;
    vector<double> AG5_V;
    vector<double> PGFILTER_SCALE;
    vector<double> CESMADCM_SCALE;
    vector<double> EXTFIELD_SCALE;
    vector<uint> EXTFIELD_CASE;
    vector<double> T_PERP;
    vector<double> T_PAR;
    vector<double> E0_Z;
    vector<double> U_PLASMA;
    vector<uint> B_ISON;
    vector<uint> ACCELERATOR_IDX;
    vector<uint> STRIPPING;
    vector<double> M_IONS;
    vector<int> Q_IONS;
    vector<uint> N_PARTICLES;
    vector<double> MESH_SIZE;
    vector<uint> ITERATIONS;
    vector<uint> SPLIT_DOMAIN;
    vector<double> JTOLERANCE;
    vector<double> ALPHA_COEFF;
    vector<double> T_POSITIVE;
    vector<double> EXT_GAP;
    vector<double> ACC_GAP;
    vector<double> DOMAIN_X_SIZE;
    vector<double> DOMAIN_Y_SIZE;
    vector<double> DOMAIN_Z_SIZE;
    vector<uint> MGSOLVER;

    double temp;

    for (uint ii=0; ii<scan_length; ii++) {

        std::string inpfile = scan_file_tag+"_"+to_string(ii);
        input_file_tags.push_back(inpfile);

        std::cout <<  "Expected " << ii << ": " << inpfile << "\n";
        std::cout <<  "Input file tag " << ii << ": " << input_file_tags[ii] << "\n";

        // J_ION
        variable_name = "J_ION";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=330.; // default value
        J_ION.push_back(temp);

        std::cout <<  variable_name+" caricata: " << J_ION.back() << "\n";

        // J_E
        variable_name = "J_E";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=0.; // default value
        J_E.push_back(temp);

        // EG_V
        variable_name = "EG_V";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=8e3; // default value
        EG_V.push_back(temp);

        // AG1_V
        variable_name = "AG1_V";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=182e3; // default value
        AG1_V.push_back(temp);

        // AG2_V
        variable_name = "AG2_V";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=356e3; // default value
        AG2_V.push_back(temp);

        // AG3_V
        variable_name = "AG3_V";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=530e3; // default value
        AG3_V.push_back(temp);

        // AG4_V
        variable_name = "AG4_V";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=704e3; // default value
        AG4_V.push_back(temp);

        // AG5_V
        variable_name = "AG5_V";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=878e3; // default value
        AG5_V.push_back(temp);

        // PGFILTER_SCALE
        variable_name = "PGFILTER_SCALE";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=0.; // default value
        PGFILTER_SCALE.push_back(temp);

        // CESMADCM_SCALE
        variable_name = "CESMADCM_SCALE";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=1.; // default value
        CESMADCM_SCALE.push_back(temp);

        // EXTFIELD_SCALE
        variable_name = "EXTFIELD_SCALE";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=0.; // default value
        EXTFIELD_SCALE.push_back(temp);

        // EXTFIELD_CASE
        variable_name = "EXTFIELD_CASE";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=0; // default value
        EXTFIELD_CASE.push_back(uint(temp));

        // T_PERP
        variable_name = "T_PERP";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=1; // default value
        T_PERP.push_back(temp);

        // T_PAR
        variable_name = "T_PAR";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=0; // default value
        T_PAR.push_back(temp);

        // E0_Z
        variable_name = "E0_Z";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=4.; // default value
        E0_Z.push_back(temp);

        // U_PLASMA
        variable_name = "U_PLASMA";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=3.; // default value
        U_PLASMA.push_back(temp);

        // B_ISON
        variable_name = "B_ISON";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=1; // default value
        B_ISON.push_back(uint(temp));

        // ACCELERATOR_IDX
        variable_name = "ACCELERATOR_IDX";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=2; // default value
        ACCELERATOR_IDX.push_back(uint(temp));

        // STRIPPING
        variable_name = "STRIPPING";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=0; // default value
        STRIPPING.push_back(uint(temp));

        // M_IONS
        variable_name = "M_IONS";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=1; // default value
        M_IONS.push_back(temp);

        // Q_IONS
        variable_name = "Q_IONS";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=-1; // default value
        Q_IONS.push_back(int(temp));

        // N_PARTICLES
        variable_name = "N_PARTICLES";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=120000; // default value
        N_PARTICLES.push_back(uint(temp));

        // MESH_SIZE
        variable_name = "MESH_SIZE";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=0.5e-3; // default value
        MESH_SIZE.push_back(temp);

        // ITERATIONS
        variable_name = "ITERATIONS";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=15; // default value
        ITERATIONS.push_back(uint(temp));

        // SPLIT DOMAIN
        variable_name = "SPLIT_DOMAIN";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=1; // default value
        SPLIT_DOMAIN.push_back(uint(temp));

        // J TOLERANCE
        variable_name = "JTOLERANCE";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=0.1; // default value
        JTOLERANCE.push_back(double(temp));

        // ALPHA COEFFICIENT
        variable_name = "ALPHA_COEFF";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=0.3; // default value
        ALPHA_COEFF.push_back(double(temp));

        // POSITIVE ION TEMPERATURE
        variable_name = "T_POSITIVE";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=0.8; // default value
        T_POSITIVE.push_back(double(temp));

        // EXTRACTION GAP LENGTH
        variable_name = "EXT_GAP";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=0.006; // default value
        EXT_GAP.push_back(double(temp));

        // ACCELERATION GAP LENGTH
        variable_name = "ACC_GAP";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=0.088; // default value
        ACC_GAP.push_back(double(temp));

        // DOMAIN X SIZE
        variable_name = "DOMAIN_X_SIZE";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=30.; // default value in mm
        DOMAIN_X_SIZE.push_back(double(temp));

        // DOMAIN Y SIZE
        variable_name = "DOMAIN_Y_SIZE";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=30.; // default value in mm
        DOMAIN_Y_SIZE.push_back(double(temp));

        // DOMAIN Z SIZE
        variable_name = "DOMAIN_Z_SIZE";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=567.; // default value in mm (MTF default)
        DOMAIN_Z_SIZE.push_back(double(temp));

        // DOMAIN Z SIZE
        variable_name = "MGSOLVER";
        if (data.find(variable_name) != data.end()){
            const std::vector<double>& values = data[variable_name];
            if (values.size()<scan_length) temp=values[0];
            else temp=values[ii];
        }
        else temp=0; // default value in mm (MTF default)
        MGSOLVER.push_back(uint(temp));

    }

    std::cout << "Input values loaded!\n"; 

    for ( uint ii=0; ii<scan_length; ii++ ) {

        std::string outfile = scan_file_tag+"/"+input_file_tags[ii]+".inp";
        std::cout << "Creating file... " << outfile << "... "; 
        std::ofstream pout;
        pout.open( outfile.c_str() );

        pout    << ACCELERATOR_IDX[ii] << "\t// ACCELERATOR INDEX (MITICA=2) \n"
                << B_ISON[ii] << "\t// B_ISON (TRUE=1) \n"
                << STRIPPING[ii] << "\t// STRIPPING (TRUE=1) \n"
                << J_E[ii] << "\t// J ELECTRONS\n"
                << M_IONS[ii] << "\t// M_IONS \n"
                << Q_IONS[ii] << "\t// Q_IONS \n"
                << J_ION[ii] << "\t// J_IONS \n"
                << T_PERP[ii] << "\t// T_PERP \n"
                << T_PAR[ii] << "\t// T_PAR \n"
                << E0_Z[ii] << "\t// BEAM STARTING AXIAL ENERGY \n"
                << U_PLASMA[ii] << "\t// PLASMA POTENTIAL \n"
                << N_PARTICLES[ii] << "\t// NUMBER OF PARTICLES \n"
                << EG_V[ii] << "\t// EG VOLTAGE \n"
                << AG1_V[ii] << "\t// AG1 VOLTAGE \n"
                << AG2_V[ii] << "\t// AG2 VOLTAGE \n"
                << AG3_V[ii] << "\t// AG3 VOLTAGE \n"
                << AG4_V[ii] << "\t// AG4 VOLTAGE \n"
                << AG5_V[ii] << "\t// AG5 VOLTAGE \n"
                << MESH_SIZE[ii] << "\t// MESH SIZE \n"
                << ITERATIONS[ii] << "\t// ITERATIONS \n"
                << PGFILTER_SCALE[ii] << "\t// PG FILTER FIELD SCALE \n"
                << CESMADCM_SCALE[ii] << "\t// CESM+ADCM FIELD SCALE \n"
                << EXTFIELD_CASE[ii] << "\t// EXTERNAL FIELD CASE \n"
                << EXTFIELD_SCALE[ii] << "\t// EXTERNAL FIELD SCALE \n"
                << SPLIT_DOMAIN[ii] <<  "\t// SPLIT DOMAIN (TRUE=1) \n"
                << JTOLERANCE[ii] <<  "\t// TOLERANCE ON EXTRACTED CURRENT \n"
                << ALPHA_COEFF[ii] <<  "\t// ALPHA COEFFICIENT FOR SC AVERAGING \n"
                << T_POSITIVE[ii] <<  "\t// POSITIVE ION TEMPERATURE \n"
                << EXT_GAP[ii] <<  "\t// EXTRACTION GAP LENGTH \n"
                << ACC_GAP[ii] <<  "\t// ACCELERATION GAP LENGTH \n"
                << DOMAIN_X_SIZE[ii] << "\t// DOMAIN X SIZE [mm] \n"
                << DOMAIN_Y_SIZE[ii] << "\t// DOMAIN Y SIZE [mm] \n"
                << DOMAIN_Z_SIZE[ii] << "\t// DOMAIN Z SIZE [mm] \n"
                << MGSOLVER[ii] << "\t// MGSOLVER \n";

        pout.close();
        std::cout << "Done!\n"; 
    }


    std::cout << "Input files created!\n";

    return input_file_tags;




    // // Retrieving values of variable "J_ION"
    // std::string variable_name = "J_ION";
    // if (data.find(variable_name) != data.end()) {
    //     std::cout << "Values of variable " << variable_name << ": ";
    //     for (double value : data[variable_name]) {
    //         std::cout << value << " ";
    //     }
    //     std::cout << std::endl;
    // } else {
    //     std::cerr << "Variable " << variable_name << " not found!" << std::endl;
    // }

    // return 0;
};



#endif /* READINPUTS_H_ */

