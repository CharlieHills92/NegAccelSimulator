/*
 * FileManager.h
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored from ManageSimulation)
 */

#ifndef FILEMANAGER_H_
#define FILEMANAGER_H_

#include <string>
#include <vector>

/**
 * @class FileManager
 * @brief Handles file operations and directory management for simulations
 * 
 * This class is responsible for:
 * - Creating and managing directory structures
 * - Generating and managing file paths
 * - File naming conventions
 */
class FileManager {
private:
    // File management
    std::string file_tag;
    std::string inputfile;
    std::string geomfile;
    std::string pdbfile;
    std::string epotfile;
    std::string schargefile;
    std::string bfieldfn;
    std::string bfieldfile;
    std::string diagfile;
    std::vector<std::string> diagfile_species;
    std::string iterationsfileEG;
    std::string iterationsfileOUT;
    
    // Directory structure
    std::string ref_fold;
    std::string outsummary_fold;
    std::string plot_fold;
    std::string bfield_fold;
    std::string data_fold;
    std::string vtk_fold;

public:
    /**
     * @brief Constructor
     * @param foldername Base folder name for the simulation
     */
    FileManager(const std::string& foldername,
                const std::string& summary_dir,
                const std::string& plot_dir,
                const std::string& data_dir,
                const std::string& vtk_dir);

    /**
     * @brief Destructor
     */
    ~FileManager() = default;

    /**
     * @brief Create directory structure for simulation
     * @param foldername Base folder name
     */
    void createDirectoryStructure(const std::string& foldername,
                                  const std::string& summary_dir,
                                  const std::string& plot_dir,
                                  const std::string& data_dir,
                                  const std::string& vtk_dir);

    /**
     * @brief Set file tag and update all file paths
     * @param filetag File tag to use for naming
     */
    void setFileTag(const std::string& filetag);

    /**
     * @brief Remove legacy summary artifacts for the active case tag.
     * @return Filenames removed from the summary directory.
     */
    std::vector<std::string> removeLegacySummaryArtifacts() const;

    // Directory getters
    const std::string& getRefFolder() const { return ref_fold; }
    const std::string& getOutputSummaryFolder() const { return outsummary_fold; }
    const std::string& getPlotFolder() const { return plot_fold; }
    const std::string& getBFieldFolder() const { return bfield_fold; }
    const std::string& getDataFolder() const { return data_fold; }
    const std::string& getVTKFolder() const { return vtk_fold; }

    // File path getters
    const std::string& getFileTag() const { return file_tag; }
    const std::string& getInputFile() const { return inputfile; }
    const std::string& getGeomFile() const { return geomfile; }
    const std::string& getPdbFile() const { return pdbfile; }
    const std::string& getEpotFile() const { return epotfile; }
    const std::string& getSchargeFile() const { return schargefile; }
    const std::string& getBFieldFn() const { return bfieldfn; }
    const std::string& getBFieldFile() const { return bfieldfile; }
    const std::string& getDiagFile() const { return diagfile; }
    const std::vector<std::string>& getDiagFileSpecies() const { return diagfile_species; }
    const std::string& getIterationsFileEG() const { return iterationsfileEG; }
    const std::string& getIterationsFileOUT() const { return iterationsfileOUT; }
    std::string buildIterationDiagnosticFile(const std::string& label) const;
    std::string buildIterationDiagnosticFileForZPosition(double z_position_meters) const;
    std::string buildIterationVTKBase(unsigned int iteration) const;

    // Setters for specific fields
    void setBFieldFolder(const std::string& folder) { bfield_fold = folder; }
    void setBFieldFn(const std::string& fn) { bfieldfn = fn; }

private:
    /**
     * @brief Update all file paths based on current file_tag
     */
    void updateFilePaths();
};

#endif /* FILEMANAGER_H_ */
