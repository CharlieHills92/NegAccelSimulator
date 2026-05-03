#include "SimulationParameters.h"
#include <iostream>
#include <string>

using namespace std;

int main() {
    try {
        SimulationParameters params;
        string scenario_file = "MTF_FULL_80x80.scn";
        
        cout << "Generating input file from scenario: " << scenario_file << endl;
        params.parseScenarioFile(scenario_file);
        
        string input_file = "MTF_FULL_80x80.inp";
        params.generateInputFile(input_file);
        
        cout << "Generated: " << input_file << endl;
        cout << "Domain: " << params.getDomainXSizeOrDefault()*1000 << "x" 
             << params.getDomainYSizeOrDefault()*1000 << "mm" << endl;
        cout << "Mesh: " << params.getMeshSize() << "m" << endl;
        
        return 0;
    } catch (const exception& e) {
        cout << "Error: " << e.what() << endl;
        return 1;
    }
}
