#ifndef SIMULATIONPARAMETERS_H_
#define SIMULATIONPARAMETERS_H_

#include <sys/types.h>

#include <map>
#include <string>
#include <vector>

class SimulationParameters {
public:
    struct GeometryAperturePattern {
        std::string layout;
        uint countX;
        uint countY;
        double pitchXMeters;
        double pitchYMeters;
        double marginMeters;
        double xOffsetMeters;
        double yOffsetMeters;
        double rowShiftXMeters;
        bool outsidePatternIsSolid;

        GeometryAperturePattern()
            : layout("single"),
              countX(1U),
              countY(1U),
              pitchXMeters(0.0),
              pitchYMeters(0.0),
              marginMeters(0.0),
              xOffsetMeters(0.0),
              yOffsetMeters(0.0),
              rowShiftXMeters(0.0),
              outsidePatternIsSolid(true) {}
    };

    struct GeometrySolidDefinition {
        std::string name;
        std::string kind;
        std::string role;
        int stage;
        int boundaryId;
        bool hasExplicitVoltage;
        double voltageVolts;
        std::vector<double> zProfileMeters;
        std::vector<double> rProfileMeters;
        std::vector<double> roundingRadiiMeters;
        GeometryAperturePattern aperturePattern;

        GeometrySolidDefinition()
            : stage(-1),
              boundaryId(-1),
              hasExplicitVoltage(false),
              voltageVolts(0.0) {}
    };

    struct BoundaryConditionDefinition {
        std::string name;
        std::string conditionType;
        double value;

        BoundaryConditionDefinition()
            : conditionType("dirichlet"),
              value(0.0) {}
    };

    struct ParticleTypeDefinition {
        std::string id;
        std::string name;
        std::string kind;
        double chargeState;
        double massU;
        bool sourceable;

        ParticleTypeDefinition()
            : chargeState(0.0),
              massU(0.0),
              sourceable(false) {}
    };

    struct ParticleSourceDefinition {
        std::string id;
        std::string name;
        std::string particleTypeId;
        std::string kind;
        std::string sourceModel;
        double chargeState;
        double massU;
        uint particleCount;
        double currentDensityAm2;
        double perpendicularTemperatureEV;
        double parallelTemperatureEV;
        double axialEnergyEV;
        std::vector<double> centerMeters;
        std::vector<double> mainDirection;
        std::vector<double> inPlaneReferenceDirection;
        double widthMeters;
        double heightMeters;

        ParticleSourceDefinition()
            : chargeState(0.0),
              massU(0.0),
              particleCount(0U),
              currentDensityAm2(0.0),
              perpendicularTemperatureEV(0.0),
              parallelTemperatureEV(0.0),
              axialEnergyEV(0.0),
              centerMeters(3, 0.0),
              mainDirection(3, 0.0),
              inPlaneReferenceDirection(3, 0.0),
              widthMeters(0.0),
              heightMeters(0.0) {}
    };

    struct MagneticFieldDefinition {
        std::string name;
        std::string sourceType;
        double scale;
        std::vector<double> constantValue;
        std::string filePath;

        MagneticFieldDefinition()
            : scale(1.0),
              constantValue(3, 0.0) {}
    };

        struct CrossSectionProcessProductDefinition {
            std::string particleKind;
            std::string speedClass;
            uint count;

            CrossSectionProcessProductDefinition()
                : count(1U) {}
        };

        struct CrossSectionProcessDefinition {
            std::string processId;
            std::string name;
            std::string projectileKind;
            std::string projectileFate;
                std::string sourcePath;
                uint fitDegree;
                std::vector<double> coefficients;
                bool scaleEnergyByIonMass;
                double minimumEnergyEV;
                double maximumEnergyEV;
            std::vector<CrossSectionProcessProductDefinition> products;

                CrossSectionProcessDefinition()
                        : fitDegree(0U),
                            scaleEnergyByIonMass(true),
                            minimumEnergyEV(0.0),
                            maximumEnergyEV(-1.0) {}
        };

    struct DiagnosticGridRangeDefinition {
                int id;
                bool includeInTotal;

        DiagnosticGridRangeDefinition()
                        : id(0),
                            includeInTotal(false) {}
    };

    struct DiagnosticMeniscusPlotDefinition {
        bool enabled;
        double zMinMeters;
        double zMaxMeters;
        double transverseMinMeters;
        double transverseMaxMeters;

        DiagnosticMeniscusPlotDefinition()
            : enabled(true),
              zMinMeters(0.0),
              zMaxMeters(0.034),
              transverseMinMeters(-0.01),
              transverseMaxMeters(0.01) {}
    };

private:
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
    double INITIAL_PLASMA_MAX_Z;
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
    uint SPLIT_DOMAIN;
    double JTOLERANCE;
    double ALPHA_COEFF;
    double T_POSITIVE;
    uint N_SOLIDS;
    uint MGSOLVER;
    uint SHIELD_MODEL;
    double BICGSTAB_EPS;
    uint BICGSTAB_MAX_ITERATIONS;
    double BICGSTAB_NEWTON_EPS;
    uint BICGSTAB_NEWTON_MAX_ITERATIONS;
    uint BICGSTAB_GLOBALLY_CONVERGENT_NEWTON;
    uint MG_LEVELS;
    double MG_TOLERANCE;
    uint MG_MAX_CYCLES;
    uint MG_GAMMA;
    uint MG_PRE_SMOOTH;
    uint MG_POST_SMOOTH;
    double MG_COARSE_RELAXATION;
    uint MG_COARSE_MAX_ITERATIONS;
    uint MG_LOCAL_PLASMA_MAX_ITERATIONS;

    double EXT_GAP;

    double DOMAIN_X_SIZE;
    double DOMAIN_Y_SIZE;
    double DOMAIN_Z_SIZE;
    double DOMAIN_Z_START;

    double EGEXTJ;
    uint domain_ii;
    std::string GEOMETRY_SOURCE_MODE;
    std::vector<GeometrySolidDefinition> GENERATED_GEOMETRY_SOLIDS;
    std::vector<ParticleTypeDefinition> PARTICLE_TYPES;
    std::vector<ParticleSourceDefinition> PARTICLE_SOURCES;
    std::map<int, BoundaryConditionDefinition> EXPLICIT_BOUNDARY_CONDITIONS;
    std::string MAGNETIC_FIELD_SOURCE_MODE;
    std::string MAGNETIC_FIELD_DIRECTORY;
    std::string MAGNETIC_FIELD_FILE;
    std::vector<MagneticFieldDefinition> MAGNETIC_FIELDS;
    std::vector<CrossSectionProcessDefinition> CROSS_SECTION_PROCESSES;
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
    uint OUTPUT_ITERATION_ENABLED;
    uint OUTPUT_ITERATION_EVERY_N_ITERATIONS;
    uint OUTPUT_ITERATION_EXPORT_PLANE_DIAGNOSTICS;
    uint OUTPUT_ITERATION_EXPORT_SIMULATION_STATE;
    uint OUTPUT_ITERATION_EXPORT_TRACED_PARTICLES;
    std::vector<double> OUTPUT_ITERATION_PLANE_Z_POSITIONS;
    std::vector<double> DIAGNOSTIC_SAMPLE_Z_POSITIONS;
    double DIAGNOSTIC_SUMMARY_Z_POSITION;
    double DIAGNOSTIC_EMITTER_EXPORT_Z_POSITION;
    double DIAGNOSTIC_TRANSMISSION_PLANE_Z_POSITION;
    double DIAGNOSTIC_APERTURE_RADIUS;
    uint DIAGNOSTIC_WRITE_PER_SPECIES_DIAGNOSTICS;
    uint DIAGNOSTIC_WRITE_PER_SPECIES_GRID_POWER;
    uint DIAGNOSTIC_WRITE_PER_SPECIES_PLOTS;
    uint DIAGNOSTIC_WRITE_NEGATIVE_ION_SUMMARY;
    uint DIAGNOSTIC_WRITE_POWER_DENSITY_MAP;
    std::vector<DiagnosticGridRangeDefinition> DIAGNOSTIC_GRID_POWER_RANGES;
    DiagnosticMeniscusPlotDefinition DIAGNOSTIC_MENISCUS_PLOT;

public:
    SimulationParameters();
    ~SimulationParameters() = default;

    void readParametersFromFile(const std::string& input);
    void setDefaultValues();

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
    double getInitialPlasmaMaxZ() const { return INITIAL_PLASMA_MAX_Z; }
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
    uint getSplitDomain() const { return SPLIT_DOMAIN; }
    double getJTolerance() const { return JTOLERANCE; }
    double getAlphaCoeff() const { return ALPHA_COEFF; }
    double getTPositive() const { return T_POSITIVE; }
    uint getNSolids() const { return N_SOLIDS; }
    uint getMGSolver() const { return MGSOLVER; }
    uint getShieldModel() const { return SHIELD_MODEL; }
    double getBiCGSTABEps() const { return BICGSTAB_EPS; }
    uint getBiCGSTABMaxIterations() const { return BICGSTAB_MAX_ITERATIONS; }
    double getBiCGSTABNewtonEps() const { return BICGSTAB_NEWTON_EPS; }
    uint getBiCGSTABNewtonMaxIterations() const { return BICGSTAB_NEWTON_MAX_ITERATIONS; }
    bool getBiCGSTABGloballyConvergentNewton() const { return BICGSTAB_GLOBALLY_CONVERGENT_NEWTON != 0U; }
    uint getMGLevels() const { return MG_LEVELS; }
    double getMGTolerance() const { return MG_TOLERANCE; }
    uint getMGMaxCycles() const { return MG_MAX_CYCLES; }
    uint getMGGamma() const { return MG_GAMMA; }
    uint getMGPreSmooth() const { return MG_PRE_SMOOTH; }
    uint getMGPostSmooth() const { return MG_POST_SMOOTH; }
    double getMGCoarseRelaxation() const { return MG_COARSE_RELAXATION; }
    uint getMGCoarseMaxIterations() const { return MG_COARSE_MAX_ITERATIONS; }
    uint getMGLocalPlasmaMaxIterations() const { return MG_LOCAL_PLASMA_MAX_ITERATIONS; }

    double getExtGap() const { return EXT_GAP; }
    double getEGExtJ() const { return EGEXTJ; }
    uint getDomainII() const { return domain_ii; }
    const std::string& getGeometrySourceMode() const { return GEOMETRY_SOURCE_MODE; }
    bool hasGeneratedGeometrySolids() const {
        return GEOMETRY_SOURCE_MODE == "generated-data" && !GENERATED_GEOMETRY_SOLIDS.empty();
    }
    const std::vector<GeometrySolidDefinition>& getGeneratedGeometrySolids() const {
        return GENERATED_GEOMETRY_SOLIDS;
    }
    bool hasParticleSources() const {
        return !PARTICLE_SOURCES.empty();
    }
    const std::vector<ParticleTypeDefinition>& getParticleTypes() const {
        return PARTICLE_TYPES;
    }
    const std::vector<ParticleSourceDefinition>& getParticleSources() const {
        return PARTICLE_SOURCES;
    }
    bool tryGetBoundaryCondition(int boundaryId, BoundaryConditionDefinition& definition) const {
        std::map<int, BoundaryConditionDefinition>::const_iterator it =
            EXPLICIT_BOUNDARY_CONDITIONS.find(boundaryId);
        if (it == EXPLICIT_BOUNDARY_CONDITIONS.end()) {
            return false;
        }
        definition = it->second;
        return true;
    }
    const std::string& getMagneticFieldSourceMode() const { return MAGNETIC_FIELD_SOURCE_MODE; }
    const std::string& getMagneticFieldDirectory() const { return MAGNETIC_FIELD_DIRECTORY; }
    const std::string& getMagneticFieldFile() const { return MAGNETIC_FIELD_FILE; }
    const std::vector<MagneticFieldDefinition>& getMagneticFields() const { return MAGNETIC_FIELDS; }
    const std::vector<CrossSectionProcessDefinition>& getCrossSectionProcesses() const {
        return CROSS_SECTION_PROCESSES;
    }
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
    bool getOutputIterationEnabled() const { return OUTPUT_ITERATION_ENABLED != 0U; }
    uint getOutputIterationEveryNIterations() const { return OUTPUT_ITERATION_EVERY_N_ITERATIONS; }
    bool getOutputIterationExportPlaneDiagnostics() const {
        return OUTPUT_ITERATION_EXPORT_PLANE_DIAGNOSTICS != 0U;
    }
    bool getOutputIterationExportSimulationState() const {
        return OUTPUT_ITERATION_EXPORT_SIMULATION_STATE != 0U;
    }
    bool getOutputIterationExportTracedParticles() const {
        return OUTPUT_ITERATION_EXPORT_TRACED_PARTICLES != 0U;
    }
    bool hasOutputIterationPlaneZPositions() const {
        return !OUTPUT_ITERATION_PLANE_Z_POSITIONS.empty();
    }
    const std::vector<double>& getOutputIterationPlaneZPositions() const {
        return OUTPUT_ITERATION_PLANE_Z_POSITIONS;
    }
    bool hasDiagnosticSampleZPositions() const { return !DIAGNOSTIC_SAMPLE_Z_POSITIONS.empty(); }
    const std::vector<double>& getDiagnosticSampleZPositions() const { return DIAGNOSTIC_SAMPLE_Z_POSITIONS; }
    bool hasDiagnosticSummaryZPosition() const { return DIAGNOSTIC_SUMMARY_Z_POSITION >= 0.0; }
    double getDiagnosticSummaryZPosition() const { return DIAGNOSTIC_SUMMARY_Z_POSITION; }
    bool hasDiagnosticEmitterExportZPosition() const { return DIAGNOSTIC_EMITTER_EXPORT_Z_POSITION >= 0.0; }
    double getDiagnosticEmitterExportZPosition() const { return DIAGNOSTIC_EMITTER_EXPORT_Z_POSITION; }
    bool hasDiagnosticTransmissionPlaneZPosition() const { return DIAGNOSTIC_TRANSMISSION_PLANE_Z_POSITION >= 0.0; }
    double getDiagnosticTransmissionPlaneZPosition() const { return DIAGNOSTIC_TRANSMISSION_PLANE_Z_POSITION; }
    double getDiagnosticApertureRadiusMeters() const { return DIAGNOSTIC_APERTURE_RADIUS; }
    bool getDiagnosticWritePerSpeciesDiagnostics() const { return DIAGNOSTIC_WRITE_PER_SPECIES_DIAGNOSTICS != 0U; }
    bool getDiagnosticWritePerSpeciesGridPower() const { return DIAGNOSTIC_WRITE_PER_SPECIES_GRID_POWER != 0U; }
    bool getDiagnosticWritePerSpeciesPlots() const { return DIAGNOSTIC_WRITE_PER_SPECIES_PLOTS != 0U; }
    bool getDiagnosticWriteNegativeIonSummary() const { return DIAGNOSTIC_WRITE_NEGATIVE_ION_SUMMARY != 0U; }
    bool getDiagnosticWritePowerDensityMap() const { return DIAGNOSTIC_WRITE_POWER_DENSITY_MAP != 0U; }
    bool hasDiagnosticGridPowerRanges() const { return !DIAGNOSTIC_GRID_POWER_RANGES.empty(); }
    const std::vector<DiagnosticGridRangeDefinition>& getDiagnosticGridPowerRanges() const {
        return DIAGNOSTIC_GRID_POWER_RANGES;
    }
    const DiagnosticMeniscusPlotDefinition& getDiagnosticMeniscusPlot() const {
        return DIAGNOSTIC_MENISCUS_PLOT;
    }

    double getDomainXSize() const { return DOMAIN_X_SIZE; }
    double getDomainYSize() const { return DOMAIN_Y_SIZE; }
    double getDomainZSize() const { return DOMAIN_Z_SIZE; }
    double getDomainZStart() const { return DOMAIN_Z_START; }

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
