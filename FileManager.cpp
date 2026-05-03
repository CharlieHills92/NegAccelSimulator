/*
 * FileManager.cpp
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored from ManageSimulation)
 */

#include "FileManager.h"
#include "error.hpp"
#include "funct.h"
#include "globals.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>
#include <cstring>
#include <iostream>

using namespace std;

FileManager::FileManager(const string& foldername) {
    createDirectoryStructure(foldername);
}

void FileManager::createDirectoryStructure(const string& foldername) {
    // Create main directory
    if (mkdir(foldername.c_str(), 0777) == -1) {
        if (errno != EEXIST) {
            throw Error(ERROR_LOCATION, "Could not create directory " + foldername + ": " + strerror(errno));
        }
    } else {
        logfile << "\t" << foldername << " directory created" << endl;
    }

    // Set up directory structure
    ref_fold = foldername + "/";
    outsummary_fold = ref_fold + "Summary/";
    plot_fold = ref_fold + "Plots/";
    data_fold = ref_fold + "Data/";
    
    // Create subdirectories
    const vector<string> subdirs = {outsummary_fold, plot_fold, data_fold};
    
    for (const string& dir : subdirs) {
        if (mkdir(dir.c_str(), 0777) == -1 && errno != EEXIST) {
            throw Error(ERROR_LOCATION, "Could not create directory " + dir + ": " + strerror(errno));
        } else if (errno != EEXIST) {
            logfile << "\t" << dir << " directory created" << endl;
        }
    }
}

void FileManager::setFileTag(const string& filetag) {
    file_tag = filetag;
    updateFilePaths();
}

void FileManager::updateFilePaths() {
    inputfile = ref_fold + file_tag + ".inp";
    geomfile = data_fold + file_tag + "_geom.dat";
    pdbfile = data_fold + file_tag + "_pdb.dat";
    epotfile = data_fold + file_tag + "_epot.dat";
    schargefile = data_fold + file_tag + "_scharge.dat";
    bfieldfile = data_fold + file_tag + "_bfield.dat";
    diagfile = outsummary_fold + file_tag + "_diag.txt";
    
    // Clear existing species files and rebuild
    diagfile_species.clear();
    for (int i = 0; i < 6; ++i) {
        diagfile_species.push_back(outsummary_fold + file_tag + "_diag_" + 
                                    get_particle_name(int2kind(i)) + ".txt");
    }
    
    iterationsfileEG = outsummary_fold + file_tag + "_it_EG.txt";
    iterationsfileOUT = outsummary_fold + file_tag + "_it_OUT.txt";
}
