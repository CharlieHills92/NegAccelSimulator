#ifndef SIMULATIONPARAMETERS_H_
#define SIMULATIONPARAMETERS_H_

#include <sys/types.h>

#include <string>
#include <vector>

class SimulationParameters {
private:
    uint ACCELERATOR_IDX;
    uint B_ISON;
    uint INCLUDE_STRIPPING;
    uint INCLUDE_SURFACE_COLLISIONS;
    double ELECTRONS;

    double M_IONS;
    double Q_IONS;
    double J_ION;
    double TPERP;
    double TPAR;
    double E0_Z;
    double U_PLASMA;
    uint N_PARTICLES;

    double EG_VOLTAGE;
    double GG_VOLTAGE;
    double REP_VOLTAGE;
    double G1_VOLTAGE;
    double G2_VOLTAGE;
    double G3_VOLTAGE;
    double G4_VOLTAGE;
    double G5_VOLTAGE;

    double MESH_SIZE;
    uint ITERATIONS;
    double PGFILTER_SCALE;
    double CESMADCM_SCALE;
    uint EXTFIELD_CASE;
    double EXTFIELD_SCALE;
    uint SPLIT_DOMAIN;
    double JTOLERANCE;
    double ALPHA_COEFF;
    double T_POSITIVE;
    uint N_SOLIDS;
    uint MGSOLVER;
    uint SHIELD_MODEL;

    double EXT_GAP;
    double ACC_GAP;

    double DOMAIN_X_SIZE;
    double DOMAIN_Y_SIZE;
    double DOMAIN_Z_SIZE;

    double EGEXTJ;
    uint domain_ii;
    std::string GEOMETRY_TEMPLATE;
    std::string STRIPPING_DENSITY_PROFILE;
    double STRIPPING_MIN_Z;
    double SURFACE_COLLISIONS_MIN_Z;
    uint SURFACE_COLLISIONS_DEBUG;
    uint PERIODIC_BOUNDARIES_ENABLED;
    double PERIODIC_X_MIN;
    double PERIODIC_X_MAX;
    double PERIODIC_Y_MIN;
    double PERIODIC_Y_MAX;
    std::string OUTPUT_SUMMARY_DIRECTORY;
    std::string OUTPUT_PLOTS_DIRECTORY;
    std::string OUTPUT_DATA_DIRECTORY;
    std::string OUTPUT_VTK_DIRECTORY;
    uint OUTPUT_SUMMARY_ENABLED;
    uint OUTPUT_PLOTS_ENABLED;
    uint OUTPUT_DATA_ENABLED;
    uint OUTPUT_VTK_ENABLED;
    uint OUTPUT_VTK_EXPORT_GEOMETRY;
    uint OUTPUT_VTK_EXPORT_SIMULATION_STATE;
    uint OUTPUT_VTK_EXPORT_TRACED_PARTICLES;
    std::string OUTPUT_LOGGING_CONSOLE_LEVEL;
    std::string OUTPUT_LOGGING_FILE_LEVEL;
    uint OUTPUT_LOGGING_CAPTURE_STDOUT;
    uint OUTPUT_LOGGING_WRITE_DEBUG_ARTIFACTS;
    std::string OUTPUT_LOGGING_STRUCTURED_LOG_FILE;

public:
    SimulationParameters();
    ~SimulationParameters() = default;

    void readParametersFromFile(const std::string& input);
    void setDefaultValues();

    uint getAcceleratorIdx() const { return ACCELERATOR_IDX; }
    uint getBIsOn() const { return B_ISON; }
    uint getIncludeStripping() const { return INCLUDE_STRIPPING; }
    uint getIncludeSurfaceCollisions() const { return INCLUDE_SURFACE_COLLISIONS; }
    double getElectrons() const { return ELECTRONS; }

    double getMIons() const { return M_IONS; }
    double getQIons() const { return Q_IONS; }
    double getJIon() const { return J_ION; }
    double getTPerp() const { return TPERP; }
    double getTPar() const { return TPAR; }
    double getE0Z() const { return E0_Z; }
    double getUPlasma() const { return U_PLASMA; }
    uint getNParticles() const { return N_PARTICLES; }

    double getEGVoltage() const { return EG_VOLTAGE; }
    double getGGVoltage() const { return GG_VOLTAGE; }
    double getREPVoltage() const { return REP_VOLTAGE; }
    double getG1Voltage() const { return G1_VOLTAGE; }
    double getG2Voltage() const { return G2_VOLTAGE; }
    double getG3Voltage() const { return G3_VOLTAGE; }
    double getG4Voltage() const { return G4_VOLTAGE; }
    double getG5Voltage() const { return G5_VOLTAGE; }

    double getMeshSize() const { return MESH_SIZE; }
    uint getIterations() const { return ITERATIONS; }
    double getPGFilterScale() const { return PGFILTER_SCALE; }
    double getCESMADCMScale() const { return CESMADCM_SCALE; }
    uint getExtFieldCase() const { return EXTFIELD_CASE; }
    double getExtFieldScale() const { return EXTFIELD_SCALE; }
    uint getSplitDomain() const { return SPLIT_DOMAIN; }
    double getJTolerance() const { return JTOLERANCE; }
    double getAlphaCoeff() const { return ALPHA_COEFF; }
    double getTPositive() const { return T_POSITIVE; }
    uint getNSolids() const { return N_SOLIDS; }
    uint getMGSolver() const { return MGSOLVER; }
    uint getShieldModel() const { return SHIELD_MODEL; }

    double getExtGap() const { return EXT_GAP; }
    double getAccGap() const { return ACC_GAP; }
    double getEGExtJ() const { return EGEXTJ; }
    uint getDomainII() const { return domain_ii; }
    const std::string& getGeometryTemplate() const { return GEOMETRY_TEMPLATE; }
    const std::string& getStrippingDensityProfile() const { return STRIPPING_DENSITY_PROFILE; }
    double getStrippingMinimumZ() const;
    double getSurfaceCollisionsMinimumZ() const { return SURFACE_COLLISIONS_MIN_Z; }
    bool getSurfaceCollisionsDebug() const { return SURFACE_COLLISIONS_DEBUG != 0U; }
    bool getPeriodicBoundariesEnabled() const { return PERIODIC_BOUNDARIES_ENABLED != 0U; }
    std::vector<double> getPeriodicityBounds() const;
    const std::string& getOutputSummaryDirectory() const { return OUTPUT_SUMMARY_DIRECTORY; }
    const std::string& getOutputPlotsDirectory() const { return OUTPUT_PLOTS_DIRECTORY; }
    const std::string& getOutputDataDirectory() const { return OUTPUT_DATA_DIRECTORY; }
    const std::string& getOutputVTKDirectory() const { return OUTPUT_VTK_DIRECTORY; }
    bool getOutputSummaryEnabled() const { return OUTPUT_SUMMARY_ENABLED != 0U; }
    bool getOutputPlotsEnabled() const { return OUTPUT_PLOTS_ENABLED != 0U; }
    bool getOutputDataEnabled() const { return OUTPUT_DATA_ENABLED != 0U; }
    bool getOutputVTKEnabled() const { return OUTPUT_VTK_ENABLED != 0U; }
    bool getOutputVTKExportGeometry() const { return OUTPUT_VTK_EXPORT_GEOMETRY != 0U; }
    bool getOutputVTKExportSimulationState() const { return OUTPUT_VTK_EXPORT_SIMULATION_STATE != 0U; }
    bool getOutputVTKExportTracedParticles() const { return OUTPUT_VTK_EXPORT_TRACED_PARTICLES != 0U; }
    const std::string& getOutputLoggingConsoleLevel() const { return OUTPUT_LOGGING_CONSOLE_LEVEL; }
    const std::string& getOutputLoggingFileLevel() const { return OUTPUT_LOGGING_FILE_LEVEL; }
    bool getOutputLoggingCaptureStdout() const { return OUTPUT_LOGGING_CAPTURE_STDOUT != 0U; }
    bool getOutputLoggingWriteDebugArtifacts() const { return OUTPUT_LOGGING_WRITE_DEBUG_ARTIFACTS != 0U; }
    const std::string& getOutputLoggingStructuredLogFile() const { return OUTPUT_LOGGING_STRUCTURED_LOG_FILE; }

    double getDomainXSize() const { return DOMAIN_X_SIZE; }
    double getDomainYSize() const { return DOMAIN_Y_SIZE; }
    double getDomainZSize() const { return DOMAIN_Z_SIZE; }

    double getDomainXSizeOrDefault() const;
    double getDomainYSizeOrDefault() const;
    double getDomainZSizeOrDefault() const;

    void setJIon(double j) { J_ION = j; }
    void setIterations(uint itn) { ITERATIONS = itn; }
    void setEGExtJ(double setJ) { EGEXTJ = setJ; }
    void setIncludeStripping(uint strip) { INCLUDE_STRIPPING = strip; }
    void setIncludeSurfaceCollisions(uint surf) { INCLUDE_SURFACE_COLLISIONS = surf; }

private:
    void parseJsonFile(const std::string& configFile);
    void finalizeDerivedParameters();
};

#endif /* SIMULATIONPARAMETERS_H_ */
