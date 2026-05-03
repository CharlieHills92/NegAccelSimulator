#ifndef SIMULATIONPARAMETERS_H_
#define SIMULATIONPARAMETERS_H_

#include <string>
#include <vector>
#include <fstream>

/**
 * @class SimulationParameters
 * @brief Handles reading, parsing, and managing simulation input parameters
 * 
 * This class is responsible for:
 * - Reading input files
 * - Parsing and validating parameters
 * - Providing access to simulation parameters
 */
class SimulationParameters {
public:
    // Parameter structure for metadata-driven parsing
    struct ParamInfo {
        void* variable;
        const char* name;
        const char* unit;
        enum Type { DOUBLE, UINT } type;
    };

private:
    // Simulation parameters
    uint ACCELERATOR_IDX;
    uint B_ISON;
    uint INCLUDE_STRIPPING;
    uint INCLUDE_SURFACE_COLLISIONS;  ///< Enable surface collision callbacks (EAMCC model): 0=disabled, 1=enabled
    double ELECTRONS;

    // Ion beam parameters
    double M_IONS;          ///< Ion mass [u]
    double Q_IONS;          ///< Ion charge [e]
    double J_ION;           ///< Ion current density [A/m^2]
    double TPERP;           ///< Perpendicular temperature [eV]
    double TPAR;            ///< Parallel temperature [eV]
    double E0_Z;            ///< Axial energy [eV]
    double U_PLASMA;        ///< Plasma potential for nsimp model, or meniscus voltage for shield model [V]
    uint N_PARTICLES;       ///< Number of simulation particles
    
    // Grid/voltage parameters
    double EG_VOLTAGE;      ///< Extraction grid voltage [V]
    double GG_VOLTAGE;      ///< Ground grid voltage [V]
    double REP_VOLTAGE;     ///< Repeller voltage [V]
    double G1_VOLTAGE;      ///< Grid 1 voltage [V]
    double G2_VOLTAGE;      ///< Grid 2 voltage [V]
    double G3_VOLTAGE;      ///< Grid 3 voltage [V]
    double G4_VOLTAGE;      ///< Grid 4 voltage [V]
    double G5_VOLTAGE;      ///< Grid 5 voltage [V]
    
    // Numerical parameters
    double MESH_SIZE;       ///< Mesh size [m]
    uint ITERATIONS;        ///< Number of iterations
    double PGFILTER_SCALE;  ///< PG filter scale factor
    double CESMADCM_SCALE;  ///< CESM+ADCM scale factor
    uint EXTFIELD_CASE;     ///< External field case
    double EXTFIELD_SCALE;  ///< External field scale factor
    uint SPLIT_DOMAIN;      ///< Domain splitting flag
    double JTOLERANCE;      ///< Current density tolerance
    double ALPHA_COEFF;     ///< Space charge averaging coefficient
    double T_POSITIVE;      ///< Positive ion temperature for nsimp model or meniscus temperature for shield model [eV]
    uint N_SOLIDS;          ///< Number of solid objects
    uint MGSOLVER;          ///< Multigrid solver type
    uint SHIELD_MODEL;      ///< Select SHIELD model for meniscus calculation. If 0, nsimp model is used.

    // Geometry parameters
    double EXT_GAP;         ///< Extraction gap length
    double ACC_GAP;         ///< Accelerator gap length
    
    // Domain size parameters
    double DOMAIN_X_SIZE;   ///< X domain size [m] (optional, uses accelerator defaults if not set)
    double DOMAIN_Y_SIZE;   ///< Y domain size [m] (optional, uses accelerator defaults if not set)
    double DOMAIN_Z_SIZE;   ///< Z domain size [m] (optional, uses accelerator defaults if not set)

    double EGEXTJ;          ///< Simulated extracted current density [A/m^2]
    uint domain_ii;         ///< Current domain number

    // Helper methods (no longer needed with new parsing approach)

public:
    /**
     * @brief Constructor
     */
    SimulationParameters();

    /**
     * @brief Destructor
     */
    ~SimulationParameters() = default;

    /**
     * @brief Read simulation parameters from input file
     * @param input Input filename
     */
    void readParametersFromFile(const std::string& input);

    // Getters
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
    
    // Domain size getters with accelerator-specific defaults
    double getDomainXSize() const { return DOMAIN_X_SIZE; }
    double getDomainYSize() const { return DOMAIN_Y_SIZE; }
    double getDomainZSize() const { return DOMAIN_Z_SIZE; }
    
    /**
     * @brief Get domain X size with fallback to accelerator defaults
     * @return X domain size in meters
     */
    double getDomainXSizeOrDefault() const;
    
    /**
     * @brief Get domain Y size with fallback to accelerator defaults
     * @return Y domain size in meters
     */
    double getDomainYSizeOrDefault() const;
    
    /**
     * @brief Get domain Z size with fallback to accelerator defaults
     * @return Z domain size in meters
     */
    double getDomainZSizeOrDefault() const;

    /**
     * @brief Parse scenario file (.scn) with user-defined parameters
     * @param scenarioFile Path to the .scn file
     */
    void parseScenarioFile(const std::string& scenarioFile);
    
    /**
     * @brief Parse input file (.inp) with legacy position-based format
     * @param inputFile Path to the .inp file
     */
    void parseInputFile(const std::string& inputFile);
    
    /**
     * @brief Generate complete input file (.inp) with all parameters
     * @param inputFile Path to the .inp file to generate
     */
    void generateInputFile(const std::string& inputFile);
    
    /**
     * @brief Set default values for all parameters based on accelerator type
     */
    void setDefaultValues();

    // Scan functionality
    struct ScanParameters {
        std::vector<double> stripping_values;
        std::vector<double> voltage_values;
        std::vector<double> current_values;
        std::string scan_parameter_name;
        int scan_length = 1;
    };
    
    /**
     * @brief Parse scenario file (.scn) and extract scan parameters
     * @param scnFile Path to the .scn file
     * @return ScanParameters structure with scan information
     */
    ScanParameters parseScanFile(const std::string& scnFile);
    
    /**
     * @brief Create individual parameter set for a specific scan case
     * @param scanParams Scan parameters from .scn file
     * @param caseIndex Index of the scan case (0-based)
     * @return SimulationParameters object configured for the specific case
     */
    SimulationParameters createScanCase(const ScanParameters& scanParams, int caseIndex);

    // Setters
    void setJIon(double j) { J_ION = j; }
    void setIterations(uint itn) { ITERATIONS = itn; }
    void setEGExtJ(double setJ) { EGEXTJ = setJ; }
    void setIncludeStripping(uint strip) { INCLUDE_STRIPPING = strip; }
    void setIncludeSurfaceCollisions(uint surf) { INCLUDE_SURFACE_COLLISIONS = surf; }
};

#endif /* SIMULATIONPARAMETERS_H_ */
