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

namespace {

string trimTrailingSlashes(const string& path) {
    string normalized = path;
    while (!normalized.empty() && normalized[normalized.size() - 1] == '/') {
        normalized.erase(normalized.size() - 1);
    }
    return normalized;
}

string normalizeSubdirectory(const string& directory, const string& fallback) {
    const string trimmed = trimTrailingSlashes(directory);
    return trimmed.empty() ? fallback : trimmed;
}

void ensureDirectoryExists(const string& path) {
    if (path.empty()) {
        return;
    }

    string current;
    for (size_t index = 0; index < path.size(); ++index) {
        const char ch = path[index];
        current += ch;
        if (ch != '/' && index + 1 != path.size()) {
            continue;
        }

        const string directory = trimTrailingSlashes(current);
        if (directory.empty()) {
            continue;
        }

        if (mkdir(directory.c_str(), 0777) == -1 && errno != EEXIST) {
            throw Error(ERROR_LOCATION, "Could not create directory " + directory + ": " + strerror(errno));
        }
    }
}

string joinDirectory(const string& root, const string& child) {
    return trimTrailingSlashes(root) + "/" + trimTrailingSlashes(child) + "/";
}

} // namespace

FileManager::FileManager(const string& foldername,
                         const string& summary_dir,
                         const string& plot_dir,
                         const string& data_dir,
                         const string& vtk_dir) {
    createDirectoryStructure(foldername, summary_dir, plot_dir, data_dir, vtk_dir);
}

void FileManager::createDirectoryStructure(const string& foldername,
                                           const string& summary_dir,
                                           const string& plot_dir,
                                           const string& data_dir,
                                           const string& vtk_dir) {
    ensureDirectoryExists(foldername);

    // Set up directory structure
    ref_fold = trimTrailingSlashes(foldername) + "/";
    outsummary_fold = joinDirectory(ref_fold, normalizeSubdirectory(summary_dir, "Summary"));
    plot_fold = joinDirectory(ref_fold, normalizeSubdirectory(plot_dir, "Plots"));
    data_fold = joinDirectory(ref_fold, normalizeSubdirectory(data_dir, "Data"));
    vtk_fold = joinDirectory(ref_fold, normalizeSubdirectory(vtk_dir, "VTK"));
    
    // Create subdirectories
    const vector<string> subdirs = {outsummary_fold, plot_fold, data_fold, vtk_fold};
    
    for (const string& dir : subdirs) {
        ensureDirectoryExists(dir);
    }
}

void FileManager::setFileTag(const string& filetag) {
    file_tag = filetag;
    updateFilePaths();
}

void FileManager::updateFilePaths() {
    inputfile = ref_fold + file_tag + ".json";
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
