#include "SimulationParameters.h"

#include "StrippingUtils.h"
#include "error.hpp"

#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>
#include <stdexcept>

using json = nlohmann::json;
using namespace std;

namespace {

const char* geometryTemplateFromAcceleratorIndex(uint accelerator_index) {
    switch (accelerator_index) {
        case 1U:
            return "SPIDER";
        case 2U:
            return "MITICA";
        case 3U:
            return "MTF";
        default:
            return "";
    }
}

bool isKnownGeometryTemplate(const string& geometry_template) {
    return geometry_template == "SPIDER" || geometry_template == "MITICA" || geometry_template == "MTF";
}

string defaultGeometryTemplate(uint accelerator_index) {
    return geometryTemplateFromAcceleratorIndex(accelerator_index);
}

double defaultDomainXSizeMeters(const string& geometry_template) {
    if (geometry_template == "SPIDER") {
        return 26.0e-3;
    }
    if (geometry_template == "MITICA") {
        return 30.0e-3;
    }
    if (geometry_template == "MTF") {
        return 80.0e-3;
    }

    return 30.0e-3;
}

double defaultDomainYSizeMeters(const string& geometry_template) {
    if (geometry_template == "SPIDER") {
        return 28.0e-3;
    }
    if (geometry_template == "MITICA") {
        return 30.0e-3;
    }
    if (geometry_template == "MTF") {
        return 80.0e-3;
    }

    return 30.0e-3;
}

double defaultDomainZSizeMeters(const string& geometry_template) {
    if (geometry_template == "SPIDER") {
        return 80.0e-3;
    }
    if (geometry_template == "MITICA") {
        return 554.0e-3;
    }
    if (geometry_template == "MTF") {
        return 567.0e-3;
    }

    return 567.0e-3;
}

void applyDefaultPeriodicity(const string& geometry_template,
                             double& x_min,
                             double& x_max,
                             double& y_min,
                             double& y_max,
                             uint& enabled) {
    if (geometry_template == "SPIDER") {
        x_min = -10.0e-3;
        x_max = 10.0e-3;
        y_min = -11.0e-3;
        y_max = 11.0e-3;
        enabled = 1U;
        return;
    }
    if (geometry_template == "MITICA") {
        x_min = -15.0e-3;
        x_max = 15.0e-3;
        y_min = -15.0e-3;
        y_max = 15.0e-3;
        enabled = 1U;
        return;
    }
    if (geometry_template == "MTF") {
        x_min = -40.0e-3;
        x_max = 40.0e-3;
        y_min = -40.0e-3;
        y_max = 40.0e-3;
        enabled = 1U;
        return;
    }

    x_min = 0.0;
    x_max = 0.0;
    y_min = 0.0;
    y_max = 0.0;
    enabled = 0U;
}

uint acceleratorIndexFromType(const string& acceleratorType) {
    if (acceleratorType == "SPIDER") {
        return 1U;
    }
    if (acceleratorType == "MITICA") {
        return 2U;
    }
    if (acceleratorType == "MTF") {
        return 3U;
    }
    if (acceleratorType == "ELISE") {
        return 4U;
    }
    if (acceleratorType == "NIO1") {
        return 5U;
    }
    if (acceleratorType == "BUG" || acceleratorType == "CUSTOM") {
        return 0U;
    }

    throw runtime_error("Unsupported accelerator type: " + acceleratorType);
}

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
    : ACCELERATOR_IDX(0),
      B_ISON(0),
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
      ACC_GAP(0.0),
      DOMAIN_X_SIZE(-1.0),
      DOMAIN_Y_SIZE(-1.0),
      DOMAIN_Z_SIZE(-1.0),
      EGEXTJ(0.0),
            domain_ii(0U),
            GEOMETRY_TEMPLATE(),
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
    uint accelerator_index = 0U;
    if (hasObject(root, "accelerator")) {
        const json& accelerator = root.at("accelerator");
        if (accelerator.contains("legacyIndex")) {
            accelerator_index = accelerator.at("legacyIndex").get<uint>();
        } else if (accelerator.contains("type")) {
            accelerator_index = acceleratorIndexFromType(accelerator.at("type").get<string>());
        }
    }

    *this = SimulationParameters();
    ACCELERATOR_IDX = accelerator_index;
    setDefaultValues();

    if (hasObject(root, "gasDensity")) {
        const json& gas_density = root.at("gasDensity");
        if (gas_density.contains("defaultProfile")) {
            selected_density_profile_name = gas_density.at("defaultProfile").get<string>();
        }
    }

    if (hasObject(root, "accelerator")) {
        const json& accelerator = root.at("accelerator");
        if (accelerator.contains("legacyIndex")) {
            ACCELERATOR_IDX = accelerator.at("legacyIndex").get<uint>();
        } else if (accelerator.contains("type")) {
            ACCELERATOR_IDX = acceleratorIndexFromType(accelerator.at("type").get<string>());
        }
    }

    if (hasObject(root, "geometry")) {
        const json& geometry = root.at("geometry");

        if (hasObject(geometry, "source")) {
            const json& source = geometry.at("source");
            const string mode = source.value("mode", string());
            if (mode == "builtin-generator" && source.contains("template")) {
                GEOMETRY_TEMPLATE = source.at("template").get<string>();
                if (!isKnownGeometryTemplate(GEOMETRY_TEMPLATE)) {
                    throw Error(ERROR_LOCATION, "Unsupported built-in geometry template: " + GEOMETRY_TEMPLATE);
                }

                if (DOMAIN_X_SIZE <= 0.0 || DOMAIN_X_SIZE == defaultDomainXSizeMeters(defaultGeometryTemplate(ACCELERATOR_IDX))) {
                    DOMAIN_X_SIZE = defaultDomainXSizeMeters(GEOMETRY_TEMPLATE);
                }
                if (DOMAIN_Y_SIZE <= 0.0 || DOMAIN_Y_SIZE == defaultDomainYSizeMeters(defaultGeometryTemplate(ACCELERATOR_IDX))) {
                    DOMAIN_Y_SIZE = defaultDomainYSizeMeters(GEOMETRY_TEMPLATE);
                }
                if (DOMAIN_Z_SIZE <= 0.0 || DOMAIN_Z_SIZE == defaultDomainZSizeMeters(defaultGeometryTemplate(ACCELERATOR_IDX))) {
                    DOMAIN_Z_SIZE = defaultDomainZSizeMeters(GEOMETRY_TEMPLATE);
                }

                applyDefaultPeriodicity(GEOMETRY_TEMPLATE,
                                        PERIODIC_X_MIN,
                                        PERIODIC_X_MAX,
                                        PERIODIC_Y_MIN,
                                        PERIODIC_Y_MAX,
                                        PERIODIC_BOUNDARIES_ENABLED);
            }
        }

        if (hasObject(geometry, "mesh")) {
            const json& mesh = geometry.at("mesh");
            if (mesh.contains("sizeMeters")) {
                MESH_SIZE = mesh.at("sizeMeters").get<double>();
            }
        }

        if (hasObject(geometry, "domain")) {
            const json& domain = geometry.at("domain");
            if (domain.contains("xSizeMeters")) {
                DOMAIN_X_SIZE = domain.at("xSizeMeters").get<double>();
            }
            if (domain.contains("ySizeMeters")) {
                DOMAIN_Y_SIZE = domain.at("ySizeMeters").get<double>();
            }
            if (domain.contains("zSizeMeters")) {
                DOMAIN_Z_SIZE = domain.at("zSizeMeters").get<double>();
            }
        }

        if (hasObject(geometry, "gaps")) {
            const json& gaps = geometry.at("gaps");
            if (gaps.contains("extractionGapMeters")) {
                EXT_GAP = gaps.at("extractionGapMeters").get<double>();
            }
            if (gaps.contains("accelerationGapMeters")) {
                ACC_GAP = gaps.at("accelerationGapMeters").get<double>();
            }
        }

        if (geometry.contains("solids") && geometry.at("solids").is_array()) {
            N_SOLIDS = static_cast<uint>(geometry.at("solids").size());
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
        if (magnetic_field.contains("case")) {
            EXTFIELD_CASE = magnetic_field.at("case").get<uint>();
        }
        if (magnetic_field.contains("scale")) {
            EXTFIELD_SCALE = magnetic_field.at("scale").get<double>();
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

    if (INCLUDE_STRIPPING > 0 && !selected_density_profile_name.empty() && hasObject(root, "gasDensity")) {
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
    return DOMAIN_X_SIZE > 0.0 ? DOMAIN_X_SIZE : defaultDomainXSizeMeters(GEOMETRY_TEMPLATE);
}

double SimulationParameters::getDomainYSizeOrDefault() const {
    return DOMAIN_Y_SIZE > 0.0 ? DOMAIN_Y_SIZE : defaultDomainYSizeMeters(GEOMETRY_TEMPLATE);
}

double SimulationParameters::getDomainZSizeOrDefault() const {
    if (DOMAIN_Z_SIZE > 0.0) {
        return DOMAIN_Z_SIZE;
    }

    return defaultDomainZSizeMeters(GEOMETRY_TEMPLATE);
}

void SimulationParameters::setDefaultValues() {
    const uint accel_type = ACCELERATOR_IDX;

    GEOMETRY_TEMPLATE = defaultGeometryTemplate(accel_type);

    B_ISON = 1U;
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
    ACC_GAP = 0.088;

    if (DOMAIN_X_SIZE < 0.0) {
        DOMAIN_X_SIZE = defaultDomainXSizeMeters(GEOMETRY_TEMPLATE);
    }
    if (DOMAIN_Y_SIZE < 0.0) {
        DOMAIN_Y_SIZE = defaultDomainYSizeMeters(GEOMETRY_TEMPLATE);
    }
    if (DOMAIN_Z_SIZE < 0.0) {
        DOMAIN_Z_SIZE = defaultDomainZSizeMeters(GEOMETRY_TEMPLATE);
    }

    EGEXTJ = 0.0;
    domain_ii = 0U;
    STRIPPING_DENSITY_PROFILE = getDensityProfileFilename(ACCELERATOR_IDX);
    STRIPPING_MIN_Z = -1.0;
    SURFACE_COLLISIONS_MIN_Z = 7.0e-3;
    SURFACE_COLLISIONS_DEBUG = 0U;
    applyDefaultPeriodicity(GEOMETRY_TEMPLATE,
                            PERIODIC_X_MIN,
                            PERIODIC_X_MAX,
                            PERIODIC_Y_MIN,
                            PERIODIC_Y_MAX,
                            PERIODIC_BOUNDARIES_ENABLED);
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
