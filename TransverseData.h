/*
 * TransverseData.h
 *
 *  Created on: Mar 05, 2024
 *      Author: Carlo Poggi
 */

#ifndef TRANSVERSEDATA_H_
#define TRANSVERSEDATA_H_




#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>

using namespace std;

class TransverseData {
private:
    vector<double> it;
    vector<double> zloc;
    vector<double> current;
    vector<double> xave;
    vector<double> xmax;
    vector<double> xmin;
    vector<double> yave;
    vector<double> ymax;
    vector<double> ymin;
    vector<double> xpave;
    vector<double> ypave;
    vector<double> divx;
    vector<double> divy;
    vector<double> epot;
    vector<double> Bfield;

public:

	TransverseData( vector<double>& iteridx, vector<double>& zz, vector<double>& cc,
                    vector<double>& xav, vector<double>& xmi, vector<double>& xma,
                    vector<double>& yav, vector<double>& ymi, vector<double>& yma,
                    vector<double>& xpav, vector<double>& ypav, vector<double>& dvx, vector<double>& dvy,
                    vector<double>& ep, vector<double>& Bf ) {
                        it=iteridx; zloc=zz; current=cc; xave=xav; xmin=xmi; xmax=xma; yave=yav; ymin=ymi; ymax=yma;
                        xpave=xpav; ypave=ypav; divx=dvx; divy=dvy; epot=ep; Bfield=Bf;
                    };
	TransverseData( const string& filename ) {readTrasverseDataFromFile( filename);};
	~TransverseData() {clear();};

    vector<double> get_it() {return it;};
    vector<double> get_zloc() {return zloc;};
    vector<double> get_current() {return current;};
    vector<double> get_xave() {return xave;};
    vector<double> get_xmax() {return xmax;};
    vector<double> get_xmin() {return xmin;};
    vector<double> get_yave() {return yave;};
    vector<double> get_ymax() {return ymax;};
    vector<double> get_ymin() {return ymin;};
    vector<double> get_xpave() {return xpave;};
    vector<double> get_ypave() {return ypave;};
    vector<double> get_divx() {return divx;};
    vector<double> get_divy() {return divy;};
    vector<double> get_epot() {return epot;};
    vector<double> get_Bfield() {return Bfield;};

    double get_it(size_t idx) {return it[idx];};
    double get_zloc(size_t idx) {return zloc[idx];};
    double get_current(size_t idx) {return current[idx];};
    double get_xave(size_t idx) {return xave[idx];};
    double get_xmax(size_t idx) {return xmax[idx];};
    double get_xmin(size_t idx) {return xmin[idx];};
    double get_yave(size_t idx) {return yave[idx];};
    double get_ymax(size_t idx) {return ymax[idx];};
    double get_ymin(size_t idx) {return ymin[idx];};
    double get_xpave(size_t idx) {return xpave[idx];};
    double get_ypave(size_t idx) {return ypave[idx];};
    double get_divx(size_t idx) {return divx[idx];};
    double get_divy(size_t idx) {return divy[idx];};
    double get_epot(size_t idx) {return epot[idx];};
    double get_Bfield(size_t idx) {return Bfield[idx];};

    void clear() {
        it.clear();
        zloc.clear();
        current.clear();
        xave.clear();
        xmax.clear();
        xmin.clear();
        yave.clear();
        ymax.clear();
        ymin.clear();
        xpave.clear();
        ypave.clear();
        divx.clear();
        divy.clear();
        epot.clear();
        Bfield.clear();
    };

    // Function to read data from file and return it as a struct
    void readTrasverseDataFromFile(const string& filename) {

        // Read data from the file
        ifstream file(filename);
        if (file.is_open()) {
            string line;
            // Skip the header line
            getline(file, line);
            
            // Read each line in the file
            while (getline(file, line)) {
                istringstream iss(line);
                double value;

                // Read values from the line and store them in variables
                iss >> value;
                it.push_back(value/1e3);
                iss >> value;
                zloc.push_back(value/1e3);
                iss >> value;
                current.push_back(value);
                iss >> value;
                xave.push_back(value/1e3);
                iss >> value;
                xmax.push_back(value/1e3);
                iss >> value;
                xmin.push_back(value/1e3);
                iss >> value;
                yave.push_back(value/1e3);
                iss >> value;
                ymax.push_back(value/1e3);
                iss >> value;
                ymin.push_back(value/1e3);
                iss >> value;
                xpave.push_back(value/1e3);
                iss >> value;
                ypave.push_back(value/1e3);
                iss >> value;
                divx.push_back(value/1e3);
                iss >> value;
                divy.push_back(value/1e3);
                iss >> value;
                epot.push_back(value*1e3);
                iss >> value;
                Bfield.push_back(value/1e3);
            }

            file.close();
        } else {
            cerr << "Unable to open file." << endl;
        }

    };



};





#endif



