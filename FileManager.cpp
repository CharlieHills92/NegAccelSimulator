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

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>

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

bool startsWith(const string& value, const string& prefix) {
    return value.size() >= prefix.size() && value.compare(0, prefix.size(), prefix) == 0;
}

bool endsWith(const string& value, const string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool isAllDigits(const string& value) {
    if (value.empty()) {
        return false;
    }

    for (size_t index = 0; index < value.size(); ++index) {
        if (!std::isdigit(static_cast<unsigned char>(value[index]))) {
            return false;
        }
    }

    return true;
}

bool isLegacySpeciesDiagnostic(const string& file_tag, const string& filename) {
    const string prefix = file_tag + "_species_";
    const string suffix = "_diagnostic_summary.txt";
    if (!startsWith(filename, prefix) || !endsWith(filename, suffix)) {
        return false;
    }

    const string species_index = filename.substr(
        prefix.size(), filename.size() - prefix.size() - suffix.size());
    return isAllDigits(species_index);
}

bool isLegacyRunSummary(const string& file_tag, const string& filename) {
    const string prefix = file_tag + "_";
    const string suffix = "_summary.txt";
    if (!startsWith(filename, prefix) || !endsWith(filename, suffix)) {
        return false;
    }

    const string summary_index = filename.substr(
        prefix.size(), filename.size() - prefix.size() - suffix.size());
    return isAllDigits(summary_index);
}

bool isLegacyPartsAtSolidDump(const string& file_tag, const string& filename) {
    const string prefix = file_tag + "_";
    const string marker = "_parts_at_solid_";
    const string suffix = ".dat";
    if (!startsWith(filename, prefix) || !endsWith(filename, suffix)) {
        return false;
    }

    const size_t marker_pos = filename.find(marker, prefix.size());
    if (marker_pos == string::npos) {
        return false;
    }

    const size_t digits_start = marker_pos + marker.size();
    const string solid_index = filename.substr(
        digits_start, filename.size() - digits_start - suffix.size());
    return isAllDigits(solid_index);
}

bool isLegacySummaryArtifact(const string& file_tag, const string& filename) {
    return filename == file_tag + "_particlesatgrids.txt" ||
           isLegacySpeciesDiagnostic(file_tag, filename) ||
           isLegacyRunSummary(file_tag, filename) ||
           isLegacyPartsAtSolidDump(file_tag, filename);
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

string FileManager::buildIterationDiagnosticFile(const string& label) const {
    return outsummary_fold + file_tag + "_it_" + label + ".txt";
}

string FileManager::buildIterationDiagnosticFileForZPosition(double z_position_meters) const {
    const long long z_position_microns = static_cast<long long>(std::llround(z_position_meters * 1.0e6));
    const long long magnitude_microns = z_position_microns < 0 ? -z_position_microns : z_position_microns;
    ostringstream label;
    if (z_position_microns < 0) {
        label << "zm" << magnitude_microns << "um";
    } else {
        label << "z" << magnitude_microns << "um";
    }
    return buildIterationDiagnosticFile(label.str());
}

string FileManager::buildIterationVTKBase(unsigned int iteration) const {
    ostringstream builder;
    builder << vtk_fold << file_tag << "_it" << setw(4) << setfill('0') << iteration;
    return builder.str();
}

vector<string> FileManager::removeLegacySummaryArtifacts() const {
    vector<string> removed_files;

    if (file_tag.empty() || outsummary_fold.empty()) {
        return removed_files;
    }

    DIR* directory = opendir(outsummary_fold.c_str());
    if (!directory) {
        if (errno == ENOENT) {
            return removed_files;
        }
        throw Error(ERROR_LOCATION,
                    "Could not open summary directory " + outsummary_fold + ": " + strerror(errno));
    }

    errno = 0;
    dirent* entry = NULL;
    while ((entry = readdir(directory)) != NULL) {
        const string filename = entry->d_name;
        if (filename == "." || filename == "..") {
            continue;
        }

        if (!isLegacySummaryArtifact(file_tag, filename)) {
            continue;
        }

        const string filepath = outsummary_fold + filename;
        if (std::remove(filepath.c_str()) == 0) {
            removed_files.push_back(filename);
            continue;
        }

        if (errno != ENOENT) {
            const int remove_errno = errno;
            closedir(directory);
            throw Error(ERROR_LOCATION,
                        "Could not remove legacy summary artifact " + filepath + ": " +
                        strerror(remove_errno));
        }
    }

    const int read_errno = errno;
    closedir(directory);
    if (read_errno != 0) {
        throw Error(ERROR_LOCATION,
                    "Error while reading summary directory " + outsummary_fold + ": " +
                    strerror(read_errno));
    }

    return removed_files;
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
    for (size_t i = 0; i < particle_kind_count(); ++i) {
        diagfile_species.push_back(outsummary_fold + file_tag + "_diag_" + 
                                    get_particle_name(int2kind(static_cast<int>(i))) + ".txt");
    }
    
    iterationsfileEG = buildIterationDiagnosticFile("EG");
    iterationsfileOUT = buildIterationDiagnosticFile("OUT");
}
