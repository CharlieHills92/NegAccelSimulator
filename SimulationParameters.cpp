#include "SimulationParameters.h"

#include "error.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>

using json = nlohmann::json;
using namespace std;

namespace {

uint shieldModelFromType(const string& shieldType) {
    return shieldType == "shield" ? 1U : 0U;
}

uint solverFromType(const string& solverType) {
    return solverType == "multigrid" ? 1U : 0U;
}

uint strippingModeFromString(const string& mode) {
    if (mode == "disabled") {
        return 0U;
    }
    if (mode == "primaryOnly") {
        return 1U;
    }
    if (mode == "withSecondaries") {
        return 2U;
    }

    throw runtime_error("Unsupported stripping mode: " + mode);
}

bool hasObject(const json& root, const char* key) {
    return root.contains(key) && root.at(key).is_object();
}

std::vector<double> requireNumberArray(const json& document,
                                       const char* key,
                                       const std::string& context) {
    if (!document.contains(key) || !document.at(key).is_array()) {
        throw runtime_error(context + "." + key + " must be an array");
    }

    std::vector<double> values;
    for (json::const_iterator it = document.at(key).begin(); it != document.at(key).end(); ++it) {
        if (!it->is_number()) {
            throw runtime_error(context + "." + key + " must contain only numeric values");
        }
        values.push_back(it->get<double>());
    }

    if (values.size() < 2) {
        throw runtime_error(context + "." + key + " must contain at least two values");
    }

    return values;
}

SimulationParameters::GeometryAperturePattern parseGeometryAperturePattern(
    const json& aperture_pattern,
    const std::string& context) {
    SimulationParameters::GeometryAperturePattern parsed;
    parsed.layout = aperture_pattern.value("layout", std::string("single"));
    parsed.countX = aperture_pattern.value("countX", 1U);
    parsed.countY = aperture_pattern.value("countY", 1U);
    parsed.pitchXMeters = aperture_pattern.value("pitchXMeters", 0.0);
    parsed.pitchYMeters = aperture_pattern.value("pitchYMeters", 0.0);
    parsed.marginMeters = aperture_pattern.value("marginMeters", 0.0);
    parsed.xOffsetMeters = aperture_pattern.value("xOffsetMeters", 0.0);
    parsed.yOffsetMeters = aperture_pattern.value("yOffsetMeters", 0.0);
    parsed.rowShiftXMeters = aperture_pattern.value("rowShiftXMeters", 0.0);
    parsed.outsidePatternIsSolid = aperture_pattern.value("outsidePatternIsSolid", true);

    if (parsed.layout != "single" &&
        parsed.layout != "rectangular-grid" &&
        parsed.layout != "staggered-grid") {
        throw runtime_error(context + ".layout must be single, rectangular-grid, or staggered-grid");
    }

    if (parsed.layout != "single") {
        if (parsed.countX == 0U || parsed.countY == 0U) {
            throw runtime_error(context + ".countX and countY must be positive for patterned solids");
        }
        if (parsed.pitchXMeters <= 0.0 || parsed.pitchYMeters <= 0.0) {
            throw runtime_error(context + ".pitchXMeters and pitchYMeters must be positive for patterned solids");
        }
    }

    return parsed;
}

std::string normalizeGeneratedSolidKind(const std::string& kind,
                                        const std::string& context) {
    if (kind == "solid" || kind == "diagnosticPlane") {
        return kind;
    }
    if (kind == "grid" ||
        kind == "aperture" ||
        kind == "wall" ||
        kind == "plasma" ||
        kind == "custom") {
        return "solid";
    }
    throw runtime_error(context + ".kind must be either solid or diagnosticPlane");
}

SimulationParameters::GeometrySolidDefinition parseGeneratedGeometrySolid(
    const json& solid,
    size_t index) {
    const std::string context = "geometry.solids[" + std::to_string(index) + "]";
    if (!solid.is_object()) {
        throw runtime_error(context + " must be an object");
    }

    if (!solid.contains("name") || !solid.at("name").is_string()) {
        throw runtime_error(context + ".name must be a string");
    }
    if (!solid.contains("kind") || !solid.at("kind").is_string()) {
        throw runtime_error(context + ".kind must be a string");
    }
    SimulationParameters::GeometrySolidDefinition parsed;
    parsed.name = solid.at("name").get<std::string>();
    parsed.kind = normalizeGeneratedSolidKind(solid.at("kind").get<std::string>(), context);
    parsed.role = solid.value("role", std::string());
    parsed.stage = solid.value("stage", -1);
    parsed.boundaryId = solid.value("boundaryId", -1);
    if (solid.contains("voltageVolts")) {
        parsed.hasExplicitVoltage = true;
        parsed.voltageVolts = solid.at("voltageVolts").get<double>();
    }

    if (!solid.contains("zProfileMeters") || !solid.contains("rProfileMeters")) {
        throw runtime_error(context + ".zProfileMeters and .rProfileMeters are required");
    }

    parsed.zProfileMeters = requireNumberArray(solid, "zProfileMeters", context);
    parsed.rProfileMeters = requireNumberArray(solid, "rProfileMeters", context);

    if (parsed.zProfileMeters.size() != parsed.rProfileMeters.size()) {
        throw runtime_error(context + ".zProfileMeters and rProfileMeters must have the same length");
    }
    bool has_positive_span = false;
    for (size_t point_index = 1; point_index < parsed.zProfileMeters.size(); ++point_index) {
        if (parsed.zProfileMeters[point_index] < parsed.zProfileMeters[point_index - 1]) {
            throw runtime_error(context + ".zProfileMeters must be sorted in non-decreasing order");
        }
        if (parsed.zProfileMeters[point_index] > parsed.zProfileMeters[point_index - 1]) {
            has_positive_span = true;
        }
    }
    if (!has_positive_span) {
        throw runtime_error(context + ".zProfileMeters must span a non-zero z range");
    }

    parsed.roundingRadiiMeters.assign(parsed.zProfileMeters.size(), 0.0);
    if (solid.contains("roundingRadiiMeters")) {
        parsed.roundingRadiiMeters = requireNumberArray(solid, "roundingRadiiMeters", context);
        if (parsed.roundingRadiiMeters.size() != parsed.zProfileMeters.size()) {
            throw runtime_error(context + ".roundingRadiiMeters must have the same length as zProfileMeters");
        }
    }
    if (hasObject(solid, "aperturePattern")) {
        parsed.aperturePattern = parseGeometryAperturePattern(
            solid.at("aperturePattern"), context + ".aperturePattern");
    }

    return parsed;
}

SimulationParameters::BoundaryConditionDefinition parseBoundaryConditionDefinition(
    const json& boundary,
    size_t index,
    int& boundary_id) {
    const std::string context = "boundaryConditions.boundaries[" + std::to_string(index) + "]";
    if (!boundary.is_object()) {
        throw runtime_error(context + " must be an object");
    }
    if (!boundary.contains("boundaryId") || !boundary.at("boundaryId").is_number_integer()) {
        throw runtime_error(context + ".boundaryId must be an integer");
    }
    if (!boundary.contains("value") || !boundary.at("value").is_number()) {
        throw runtime_error(context + ".value must be numeric");
    }

    boundary_id = boundary.at("boundaryId").get<int>();
    if (boundary_id < 1) {
        throw runtime_error(context + ".boundaryId must be >= 1");
    }

    SimulationParameters::BoundaryConditionDefinition parsed;
    parsed.name = boundary.value("name", std::string());
    parsed.conditionType = boundary.value("conditionType", std::string("dirichlet"));
    std::transform(parsed.conditionType.begin(),
                   parsed.conditionType.end(),
                   parsed.conditionType.begin(),
                   [](unsigned char character) {
                       return static_cast<char>(std::tolower(character));
                   });
    if (parsed.conditionType != "dirichlet" && parsed.conditionType != "neumann") {
        throw runtime_error(context + ".conditionType must be dirichlet or neumann");
    }
    parsed.value = boundary.at("value").get<double>();
    return parsed;
}

const json* findDensityProfile(const json& profiles, const string& profile_name) {
    if (!profiles.is_array()) {
        return nullptr;
    }

    for (json::const_iterator it = profiles.begin(); it != profiles.end(); ++it) {
        if (!it->is_object()) {
            continue;
        }
        if (it->value("name", string()) == profile_name) {
            return &(*it);
        }
    }

    return nullptr;
}

} // namespace

SimulationParameters::SimulationParameters()
    : B_ISON(0),
      INCLUDE_STRIPPING(0),
      INCLUDE_SURFACE_COLLISIONS(0),
      ELECTRONS(0.0),
      M_IONS(0.0),
      Q_IONS(0.0),
      J_ION(0.0),
      TPERP(0.0),
      TPAR(0.0),
      E0_Z(0.0),
      U_PLASMA(0.0),
      N_PARTICLES(0U),
      EG_VOLTAGE(0.0),
      GG_VOLTAGE(0.0),
      REP_VOLTAGE(0.0),
      G1_VOLTAGE(0.0),
      G2_VOLTAGE(0.0),
      G3_VOLTAGE(0.0),
      G4_VOLTAGE(0.0),
      G5_VOLTAGE(0.0),
      MESH_SIZE(0.0),
      ITERATIONS(0U),
      PGFILTER_SCALE(0.0),
      CESMADCM_SCALE(0.0),
      EXTFIELD_CASE(0U),
      EXTFIELD_SCALE(0.0),
      SPLIT_DOMAIN(0U),
      JTOLERANCE(0.0),
      ALPHA_COEFF(0.0),
      T_POSITIVE(0.0),
      N_SOLIDS(0U),
      MGSOLVER(0U),
      SHIELD_MODEL(0U),
      EXT_GAP(0.0),
      DOMAIN_X_SIZE(-1.0),
      DOMAIN_Y_SIZE(-1.0),
      DOMAIN_Z_SIZE(-1.0),
      DOMAIN_Z_START(0.0),
      EGEXTJ(0.0),
      domain_ii(0U),
      GEOMETRY_SOURCE_MODE(),
      GENERATED_GEOMETRY_SOLIDS(),
      EXPLICIT_BOUNDARY_CONDITIONS(),
      MAGNETIC_FIELD_SOURCE_MODE("none"),
      MAGNETIC_FIELD_DIRECTORY(),
      MAGNETIC_FIELD_FILE(),
      STRIPPING_DENSITY_PROFILE(),
      STRIPPING_MIN_Z(-1.0),
      SURFACE_COLLISIONS_MIN_Z(7.0e-3),
      SURFACE_COLLISIONS_DEBUG(0U),
      PERIODIC_BOUNDARIES_ENABLED(0U),
      PERIODIC_X_MIN(0.0),
      PERIODIC_X_MAX(0.0),
      PERIODIC_Y_MIN(0.0),
      PERIODIC_Y_MAX(0.0),
      OUTPUT_SUMMARY_DIRECTORY("Summary"),
      OUTPUT_PLOTS_DIRECTORY("Plots"),
      OUTPUT_DATA_DIRECTORY("Data"),
      OUTPUT_VTK_DIRECTORY("VTK"),
      OUTPUT_SUMMARY_ENABLED(1U),
      OUTPUT_PLOTS_ENABLED(1U),
      OUTPUT_DATA_ENABLED(1U),
      OUTPUT_VTK_ENABLED(1U),
      OUTPUT_VTK_EXPORT_GEOMETRY(1U),
      OUTPUT_VTK_EXPORT_SIMULATION_STATE(1U),
      OUTPUT_VTK_EXPORT_TRACED_PARTICLES(1U),
      OUTPUT_LOGGING_CONSOLE_LEVEL("info"),
      OUTPUT_LOGGING_FILE_LEVEL("debug"),
      OUTPUT_LOGGING_CAPTURE_STDOUT(1U),
      OUTPUT_LOGGING_WRITE_DEBUG_ARTIFACTS(0U),
      OUTPUT_LOGGING_STRUCTURED_LOG_FILE("run.log") {
}

void SimulationParameters::readParametersFromFile(const string& input) {
    const size_t ext_pos = input.find_last_of('.');
    const string extension = (ext_pos != string::npos) ? input.substr(ext_pos) : "";

    if (extension != ".json") {
        throw Error(ERROR_LOCATION, "Only JSON configuration files are supported: " + input);
    }

    parseJsonFile(input);
    cout << "Using JSON configuration from: " << input << endl;
}

void SimulationParameters::parseJsonFile(const string& configFile) {
    ifstream file(configFile.c_str());
    if (!file.is_open()) {
        throw Error(ERROR_LOCATION, "Could not open JSON configuration: " + configFile);
    }

    json root;
    try {
        file >> root;
    } catch (const json::exception& e) {
        throw Error(ERROR_LOCATION, "Invalid JSON configuration " + configFile + ": " + string(e.what()));
    }

    if (!root.is_object()) {
        throw Error(ERROR_LOCATION, "Top-level JSON configuration must be an object: " + configFile);
    }

    string selected_density_profile_name;
    *this = SimulationParameters();
    setDefaultValues();

    if (hasObject(root, "gasDensity")) {
        const json& gas_density = root.at("gasDensity");
        if (gas_density.contains("defaultProfile")) {
            selected_density_profile_name = gas_density.at("defaultProfile").get<string>();
        }
    }

    if (!hasObject(root, "geometry")) {
        throw Error(ERROR_LOCATION, "geometry object is required in the runtime JSON configuration");
    }

    const json& geometry = root.at("geometry");
    if (!hasObject(geometry, "source")) {
        throw Error(ERROR_LOCATION, "geometry.source object is required in the runtime JSON configuration");
    }

    const json& source = geometry.at("source");
    GEOMETRY_SOURCE_MODE = source.value("mode", std::string());
    if (GEOMETRY_SOURCE_MODE != "generated-data") {
        throw Error(ERROR_LOCATION,
                    "Only geometry.source.mode='generated-data' is supported by the current C++ runtime");
    }
    if (!geometry.contains("solids") || !geometry.at("solids").is_array()) {
        throw Error(ERROR_LOCATION,
                    "geometry.solids must be an array when geometry.source.mode='generated-data'");
    }

    GENERATED_GEOMETRY_SOLIDS.clear();
    for (json::const_iterator it = geometry.at("solids").begin();
         it != geometry.at("solids").end(); ++it) {
        try {
            GENERATED_GEOMETRY_SOLIDS.push_back(
                parseGeneratedGeometrySolid(*it, GENERATED_GEOMETRY_SOLIDS.size()));
        } catch (const std::runtime_error& e) {
            throw Error(ERROR_LOCATION, e.what());
        }
    }
    N_SOLIDS = static_cast<uint>(GENERATED_GEOMETRY_SOLIDS.size());

    if (!hasObject(geometry, "mesh") || !geometry.at("mesh").contains("sizeMeters")) {
        throw Error(ERROR_LOCATION, "geometry.mesh.sizeMeters is required in the runtime JSON configuration");
    }
    MESH_SIZE = geometry.at("mesh").at("sizeMeters").get<double>();

    if (!hasObject(geometry, "domain")) {
        throw Error(ERROR_LOCATION, "geometry.domain object is required in the runtime JSON configuration");
    }
    const json& domain = geometry.at("domain");
    if (!domain.contains("xSizeMeters") || !domain.contains("ySizeMeters") ||
        !domain.contains("zSizeMeters")) {
        throw Error(ERROR_LOCATION,
                    "geometry.domain.xSizeMeters, ySizeMeters, and zSizeMeters are required");
    }
    DOMAIN_X_SIZE = domain.at("xSizeMeters").get<double>();
    DOMAIN_Y_SIZE = domain.at("ySizeMeters").get<double>();
    DOMAIN_Z_SIZE = domain.at("zSizeMeters").get<double>();
    if (domain.contains("zStartMeters")) {
        DOMAIN_Z_START = domain.at("zStartMeters").get<double>();
    }

    if (hasObject(geometry, "gaps")) {
        const json& gaps = geometry.at("gaps");
        if (gaps.contains("extractionGapMeters")) {
            EXT_GAP = gaps.at("extractionGapMeters").get<double>();
        }
        if (gaps.contains("accelerationGapMeters")) {
            throw Error(
                ERROR_LOCATION,
                "geometry.gaps.accelerationGapMeters is no longer supported by the runtime JSON contract");
        }
    }

    if (hasObject(root, "particleSources") && hasObject(root.at("particleSources"), "negativeIonBeam")) {
        const json& beam = root.at("particleSources").at("negativeIonBeam");
        if (beam.contains("massU")) {
            M_IONS = beam.at("massU").get<double>();
        }
        if (beam.contains("chargeState")) {
            Q_IONS = beam.at("chargeState").get<double>();
        }
        if (beam.contains("currentDensityAm2")) {
            J_ION = beam.at("currentDensityAm2").get<double>();
        }
        if (beam.contains("perpendicularTemperatureEV")) {
            TPERP = beam.at("perpendicularTemperatureEV").get<double>();
        }
        if (beam.contains("parallelTemperatureEV")) {
            TPAR = beam.at("parallelTemperatureEV").get<double>();
        }
        if (beam.contains("axialEnergyEV")) {
            E0_Z = beam.at("axialEnergyEV").get<double>();
        }
        if (beam.contains("plasmaPotentialVolts")) {
            U_PLASMA = beam.at("plasmaPotentialVolts").get<double>();
        }
        if (beam.contains("electronsModelWeight")) {
            ELECTRONS = beam.at("electronsModelWeight").get<double>();
        }
    }

    if (hasObject(root, "simulation")) {
        const json& simulation = root.at("simulation");
        if (simulation.contains("particleCount")) {
            N_PARTICLES = simulation.at("particleCount").get<uint>();
        }
        if (simulation.contains("iterations")) {
            ITERATIONS = simulation.at("iterations").get<uint>();
        }

        if (hasObject(simulation, "domainDecomposition")) {
            const json& decomposition = simulation.at("domainDecomposition");
            if (decomposition.contains("splitDomain")) {
                SPLIT_DOMAIN = decomposition.at("splitDomain").get<bool>() ? 1U : 0U;
            }
        }

        if (hasObject(simulation, "solver")) {
            const json& solver = simulation.at("solver");
            if (solver.contains("type")) {
                MGSOLVER = solverFromType(solver.at("type").get<string>());
            }
            if (solver.contains("shieldModel")) {
                SHIELD_MODEL = shieldModelFromType(solver.at("shieldModel").get<string>());
            }
            if (solver.contains("positiveIonTemperatureEV")) {
                T_POSITIVE = solver.at("positiveIonTemperatureEV").get<double>();
            }
            if (solver.contains("plasmaPotentialVolts")) {
                U_PLASMA = solver.at("plasmaPotentialVolts").get<double>();
            }
        }

        if (hasObject(simulation, "spaceCharge")) {
            const json& space_charge = simulation.at("spaceCharge");
            if (space_charge.contains("alphaCoeff")) {
                ALPHA_COEFF = space_charge.at("alphaCoeff").get<double>();
            }
            if (space_charge.contains("pgFilterScale")) {
                PGFILTER_SCALE = space_charge.at("pgFilterScale").get<double>();
            }
            if (space_charge.contains("cesmadcmScale")) {
                CESMADCM_SCALE = space_charge.at("cesmadcmScale").get<double>();
            }
        }

        if (hasObject(simulation, "convergence")) {
            const json& convergence = simulation.at("convergence");
            if (convergence.contains("currentDensityTolerance")) {
                JTOLERANCE = convergence.at("currentDensityTolerance").get<double>();
            }
        }
    }

    if (hasObject(root, "boundaryConditions")) {
        const json& boundary_conditions = root.at("boundaryConditions");

        if (hasObject(boundary_conditions, "plasma")) {
            const json& plasma = boundary_conditions.at("plasma");
            if (plasma.contains("potentialV")) {
                U_PLASMA = plasma.at("potentialV").get<double>();
            }
            if (plasma.contains("positiveIonTemperatureEv")) {
                T_POSITIVE = plasma.at("positiveIonTemperatureEv").get<double>();
            }
        }

        if (boundary_conditions.contains("electrodes") && boundary_conditions.at("electrodes").is_array()) {
            const json& electrodes = boundary_conditions.at("electrodes");
            for (json::const_iterator it = electrodes.begin(); it != electrodes.end(); ++it) {
                if (!it->is_object() || !it->contains("voltageVolts")) {
                    continue;
                }

                const json& electrode = *it;
                const string name = electrode.value("name", string());
                const string role = electrode.value("role", string());
                const int stage = electrode.value("stage", -1);
                const double voltage = electrode.at("voltageVolts").get<double>();

                if (role == "extraction_grid" || name == "EG") {
                    EG_VOLTAGE = voltage;
                } else if (stage == 1 || role == "ground_grid" || name == "AG1" || name == "GG" || name == "G1") {
                    GG_VOLTAGE = voltage;
                } else if (stage == 2 || role == "repeller" || name == "AG2" || name == "REP" || name == "G2") {
                    REP_VOLTAGE = voltage;
                } else if (stage == 3 || name == "AG3" || name == "G3") {
                    G3_VOLTAGE = voltage;
                } else if (stage == 4 || name == "AG4" || name == "G4") {
                    G4_VOLTAGE = voltage;
                } else if (stage == 5 || name == "AG5" || name == "G5") {
                    G5_VOLTAGE = voltage;
                }
            }
        }

        if (boundary_conditions.contains("boundaries") && boundary_conditions.at("boundaries").is_array()) {
            const json& boundaries = boundary_conditions.at("boundaries");
            EXPLICIT_BOUNDARY_CONDITIONS.clear();
            for (json::const_iterator it = boundaries.begin(); it != boundaries.end(); ++it) {
                int boundary_id = -1;
                BoundaryConditionDefinition parsed = parseBoundaryConditionDefinition(
                    *it, EXPLICIT_BOUNDARY_CONDITIONS.size(), boundary_id);
                if (!EXPLICIT_BOUNDARY_CONDITIONS.insert(std::make_pair(boundary_id, parsed)).second) {
                    throw runtime_error("Duplicate boundaryConditions.boundaries boundaryId: " +
                                        std::to_string(boundary_id));
                }
            }
        }

        if (hasObject(boundary_conditions, "periodicBoundaries")) {
            const json& periodic_boundaries = boundary_conditions.at("periodicBoundaries");
            if (periodic_boundaries.contains("enabled")) {
                PERIODIC_BOUNDARIES_ENABLED = periodic_boundaries.at("enabled").get<bool>() ? 1U : 0U;
            }
            if (periodic_boundaries.contains("xMinMeters")) {
                PERIODIC_X_MIN = periodic_boundaries.at("xMinMeters").get<double>();
            }
            if (periodic_boundaries.contains("xMaxMeters")) {
                PERIODIC_X_MAX = periodic_boundaries.at("xMaxMeters").get<double>();
            }
            if (periodic_boundaries.contains("yMinMeters")) {
                PERIODIC_Y_MIN = periodic_boundaries.at("yMinMeters").get<double>();
            }
            if (periodic_boundaries.contains("yMaxMeters")) {
                PERIODIC_Y_MAX = periodic_boundaries.at("yMaxMeters").get<double>();
            }
        }
    }

    if (hasObject(root, "externalMagneticField")) {
        const json& magnetic_field = root.at("externalMagneticField");
        if (magnetic_field.contains("enabled")) {
            B_ISON = magnetic_field.at("enabled").get<bool>() ? 1U : 0U;
        }
        if (magnetic_field.contains("scale")) {
            EXTFIELD_SCALE = magnetic_field.at("scale").get<double>();
        }
        if (B_ISON != 0U) {
            if (!magnetic_field.contains("sourceMode")) {
                throw Error(ERROR_LOCATION,
                            "externalMagneticField.sourceMode is required when the magnetic field is enabled");
            }

            MAGNETIC_FIELD_SOURCE_MODE = magnetic_field.at("sourceMode").get<string>();
            if (MAGNETIC_FIELD_SOURCE_MODE == "directory") {
                if (!magnetic_field.contains("directory")) {
                    throw Error(ERROR_LOCATION,
                                "externalMagneticField.directory is required when sourceMode='directory'");
                }
                MAGNETIC_FIELD_DIRECTORY = magnetic_field.at("directory").get<string>();
                if (magnetic_field.contains("case")) {
                    EXTFIELD_CASE = magnetic_field.at("case").get<uint>();
                }
            } else if (MAGNETIC_FIELD_SOURCE_MODE == "file") {
                if (!magnetic_field.contains("file")) {
                    throw Error(ERROR_LOCATION,
                                "externalMagneticField.file is required when sourceMode='file'");
                }
                MAGNETIC_FIELD_FILE = magnetic_field.at("file").get<string>();
            } else {
                throw Error(ERROR_LOCATION,
                            "externalMagneticField.sourceMode must be either 'directory' or 'file'");
            }
        } else {
            MAGNETIC_FIELD_SOURCE_MODE = "none";
        }
    }

    if (hasObject(root, "outputs")) {
        const json& outputs = root.at("outputs");

        if (hasObject(outputs, "summary")) {
            const json& summary = outputs.at("summary");
            if (summary.contains("enabled")) {
                OUTPUT_SUMMARY_ENABLED = summary.at("enabled").get<bool>() ? 1U : 0U;
            }
            if (summary.contains("directory")) {
                OUTPUT_SUMMARY_DIRECTORY = summary.at("directory").get<string>();
            }
        }
        if (hasObject(outputs, "plots")) {
            const json& plots = outputs.at("plots");
            if (plots.contains("enabled")) {
                OUTPUT_PLOTS_ENABLED = plots.at("enabled").get<bool>() ? 1U : 0U;
            }
            if (plots.contains("directory")) {
                OUTPUT_PLOTS_DIRECTORY = plots.at("directory").get<string>();
            }
        }
        if (hasObject(outputs, "data")) {
            const json& data = outputs.at("data");
            if (data.contains("enabled")) {
                OUTPUT_DATA_ENABLED = data.at("enabled").get<bool>() ? 1U : 0U;
            }
            if (data.contains("directory")) {
                OUTPUT_DATA_DIRECTORY = data.at("directory").get<string>();
            }
        }
        if (hasObject(outputs, "vtk")) {
            const json& vtk = outputs.at("vtk");
            if (vtk.contains("enabled")) {
                OUTPUT_VTK_ENABLED = vtk.at("enabled").get<bool>() ? 1U : 0U;
            }
            if (vtk.contains("directory")) {
                OUTPUT_VTK_DIRECTORY = vtk.at("directory").get<string>();
            }
            if (vtk.contains("exportGeometry")) {
                OUTPUT_VTK_EXPORT_GEOMETRY = vtk.at("exportGeometry").get<bool>() ? 1U : 0U;
            }
            if (vtk.contains("exportSimulationState")) {
                OUTPUT_VTK_EXPORT_SIMULATION_STATE = vtk.at("exportSimulationState").get<bool>() ? 1U : 0U;
            }
            if (vtk.contains("exportTracedParticles")) {
                OUTPUT_VTK_EXPORT_TRACED_PARTICLES = vtk.at("exportTracedParticles").get<bool>() ? 1U : 0U;
            }
        }
        if (hasObject(outputs, "logging")) {
            const json& logging = outputs.at("logging");
            if (logging.contains("consoleLevel")) {
                OUTPUT_LOGGING_CONSOLE_LEVEL = logging.at("consoleLevel").get<string>();
            }
            if (logging.contains("fileLevel")) {
                OUTPUT_LOGGING_FILE_LEVEL = logging.at("fileLevel").get<string>();
            }
            if (logging.contains("captureStdout")) {
                OUTPUT_LOGGING_CAPTURE_STDOUT = logging.at("captureStdout").get<bool>() ? 1U : 0U;
            }
            if (logging.contains("writeDebugArtifacts")) {
                OUTPUT_LOGGING_WRITE_DEBUG_ARTIFACTS = logging.at("writeDebugArtifacts").get<bool>() ? 1U : 0U;
            }
            if (logging.contains("structuredLogFile")) {
                OUTPUT_LOGGING_STRUCTURED_LOG_FILE = logging.at("structuredLogFile").get<string>();
            }
        }
    }

    if (hasObject(root, "physics")) {
        const json& physics = root.at("physics");
        if (hasObject(physics, "stripping")) {
            const json& stripping = physics.at("stripping");
            if (stripping.contains("mode")) {
                INCLUDE_STRIPPING = strippingModeFromString(stripping.at("mode").get<string>());
            }
            if (stripping.contains("minimumZMeters")) {
                STRIPPING_MIN_Z = stripping.at("minimumZMeters").get<double>();
            }
            if (stripping.contains("densityProfile")) {
                selected_density_profile_name = stripping.at("densityProfile").get<string>();
            }
        }

        if (hasObject(physics, "surfaceCollisions")) {
            const json& surface_collisions = physics.at("surfaceCollisions");
            if (surface_collisions.contains("enabled")) {
                INCLUDE_SURFACE_COLLISIONS = surface_collisions.at("enabled").get<bool>() ? 1U : 0U;
            }
            if (surface_collisions.contains("debug")) {
                SURFACE_COLLISIONS_DEBUG = surface_collisions.at("debug").get<bool>() ? 1U : 0U;
            }
            if (surface_collisions.contains("minimumImpactZMeters")) {
                SURFACE_COLLISIONS_MIN_Z = surface_collisions.at("minimumImpactZMeters").get<double>();
            }
        }
    }

    if (INCLUDE_STRIPPING > 0) {
        if (selected_density_profile_name.empty()) {
            throw Error(ERROR_LOCATION,
                        "physics.stripping.densityProfile is required when stripping is enabled");
        }
        if (!hasObject(root, "gasDensity")) {
            throw Error(ERROR_LOCATION,
                        "gasDensity object is required when stripping is enabled");
        }

        const json& gas_density = root.at("gasDensity");
        if (!gas_density.contains("profiles")) {
            throw Error(ERROR_LOCATION, "gasDensity.profiles is required when selecting a density profile");
        }

        const json* selected_density_profile = findDensityProfile(
            gas_density.at("profiles"), selected_density_profile_name);
        if (!selected_density_profile) {
            throw Error(ERROR_LOCATION,
                        "Gas density profile not found: " + selected_density_profile_name);
        }

        if (!hasObject(*selected_density_profile, "source")) {
            throw Error(ERROR_LOCATION,
                        "Gas density profile '" + selected_density_profile_name + "' is missing a source object");
        }

        const json& density_source = selected_density_profile->at("source");
        const string mode = density_source.value("mode", string());
        if (mode != "file" || !density_source.contains("path")) {
            throw Error(ERROR_LOCATION,
                        "Gas density profile '" + selected_density_profile_name +
                            "' must use source.mode='file' with a path for the current C++ runtime");
        }

        STRIPPING_DENSITY_PROFILE = density_source.at("path").get<string>();
    }

    finalizeDerivedParameters();
}

double SimulationParameters::getStrippingMinimumZ() const {
    if (STRIPPING_MIN_Z >= 0.0) {
        return STRIPPING_MIN_Z;
    }

    return 7.0e-3 + MESH_SIZE;
}

std::vector<double> SimulationParameters::getPeriodicityBounds() const {
    if (!getPeriodicBoundariesEnabled()) {
        return std::vector<double>();
    }

    return std::vector<double>{PERIODIC_X_MAX, PERIODIC_X_MIN, PERIODIC_Y_MAX, PERIODIC_Y_MIN};
}

double SimulationParameters::getDomainXSizeOrDefault() const {
    return DOMAIN_X_SIZE;
}

double SimulationParameters::getDomainYSizeOrDefault() const {
    return DOMAIN_Y_SIZE;
}

double SimulationParameters::getDomainZSizeOrDefault() const {
    return DOMAIN_Z_SIZE;
}

void SimulationParameters::setDefaultValues() {
    GEOMETRY_SOURCE_MODE.clear();
    GENERATED_GEOMETRY_SOLIDS.clear();
    EXPLICIT_BOUNDARY_CONDITIONS.clear();
    MAGNETIC_FIELD_SOURCE_MODE = "none";
    MAGNETIC_FIELD_DIRECTORY.clear();
    MAGNETIC_FIELD_FILE.clear();

    B_ISON = 0U;
    INCLUDE_STRIPPING = 0U;
    INCLUDE_SURFACE_COLLISIONS = 0U;
    ELECTRONS = 0.0;
    M_IONS = 1.0;
    Q_IONS = -1.0;
    J_ION = 330.0;
    TPERP = 0.0;
    TPAR = 0.0;
    E0_Z = 3.0;
    U_PLASMA = 3.0;
    N_PARTICLES = 150000U;

    EG_VOLTAGE = 8000.0;
    GG_VOLTAGE = 182000.0;
    REP_VOLTAGE = 356000.0;
    G3_VOLTAGE = 530000.0;
    G4_VOLTAGE = 704000.0;
    G5_VOLTAGE = 878000.0;

    MESH_SIZE = 0.0003;
    ITERATIONS = 5U;
    PGFILTER_SCALE = 0.0;
    CESMADCM_SCALE = 0.0;
    EXTFIELD_CASE = 0U;
    EXTFIELD_SCALE = 1.0;
    SPLIT_DOMAIN = 0U;
    JTOLERANCE = 1.0;
    ALPHA_COEFF = 0.3;
    T_POSITIVE = 0.8;
    MGSOLVER = 0U;
    SHIELD_MODEL = 0U;

    EXT_GAP = 0.006;
    DOMAIN_X_SIZE = -1.0;
    DOMAIN_Y_SIZE = -1.0;
    DOMAIN_Z_SIZE = -1.0;
    DOMAIN_Z_START = 0.0;

    EGEXTJ = 0.0;
    domain_ii = 0U;
    STRIPPING_DENSITY_PROFILE.clear();
    STRIPPING_MIN_Z = -1.0;
    SURFACE_COLLISIONS_MIN_Z = 7.0e-3;
    SURFACE_COLLISIONS_DEBUG = 0U;
    PERIODIC_BOUNDARIES_ENABLED = 0U;
    PERIODIC_X_MIN = 0.0;
    PERIODIC_X_MAX = 0.0;
    PERIODIC_Y_MIN = 0.0;
    PERIODIC_Y_MAX = 0.0;
    OUTPUT_SUMMARY_DIRECTORY = "Summary";
    OUTPUT_PLOTS_DIRECTORY = "Plots";
    OUTPUT_DATA_DIRECTORY = "Data";
    OUTPUT_VTK_DIRECTORY = "VTK";
    OUTPUT_SUMMARY_ENABLED = 1U;
    OUTPUT_PLOTS_ENABLED = 1U;
    OUTPUT_DATA_ENABLED = 1U;
    OUTPUT_VTK_ENABLED = 1U;
    OUTPUT_VTK_EXPORT_GEOMETRY = 1U;
    OUTPUT_VTK_EXPORT_SIMULATION_STATE = 1U;
    OUTPUT_VTK_EXPORT_TRACED_PARTICLES = 1U;
    OUTPUT_LOGGING_CONSOLE_LEVEL = "info";
    OUTPUT_LOGGING_FILE_LEVEL = "debug";
    OUTPUT_LOGGING_CAPTURE_STDOUT = 1U;
    OUTPUT_LOGGING_WRITE_DEBUG_ARTIFACTS = 0U;
    OUTPUT_LOGGING_STRUCTURED_LOG_FILE = "run.log";
    N_SOLIDS = 0U;
    finalizeDerivedParameters();
}

void SimulationParameters::finalizeDerivedParameters() {
    G1_VOLTAGE = GG_VOLTAGE;
    G2_VOLTAGE = REP_VOLTAGE;
}
