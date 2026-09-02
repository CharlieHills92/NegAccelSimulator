#include "SimulationParameters.h"

#include "cross_sections.h"
#include "error.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <set>
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

SimulationParameters::DiagnosticGridRangeDefinition parseDiagnosticGridRangeDefinition(
    const json& range,
    size_t range_index) {
    const std::string context = "diagnostics.gridPower.ranges[" + std::to_string(range_index) + "]";
    if (!range.is_object()) {
        throw runtime_error(context + " must be an object");
    }
    if (!range.contains("id") || !range.at("id").is_number_integer()) {
        throw runtime_error(context + ".id must be an integer");
    }
    if (!range.contains("includeInTotal") || !range.at("includeInTotal").is_boolean()) {
        throw runtime_error(context + ".includeInTotal must be boolean");
    }

    SimulationParameters::DiagnosticGridRangeDefinition definition;
    definition.id = range.at("id").get<int>();
    definition.includeInTotal = range.at("includeInTotal").get<bool>();
    if (definition.id < 0) {
        throw runtime_error(context + ".id must be >= 0");
    }
    return definition;
}

std::vector<double> requireFixedNumberArray(const json& document,
                                            const char* key,
                                            size_t expected_size,
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

    if (values.size() != expected_size) {
        throw runtime_error(
            context + "." + key + " must contain exactly " + std::to_string(expected_size) + " values");
    }

    return values;
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

bool isSupportedParticleKind(const std::string& kind) {
    static const char* SUPPORTED_PARTICLE_KINDS[] = {
        "H-", "H0", "H+", "H2+", "H20", "H3+",
        "D-", "D0", "D+", "D2+", "D20", "D3+",
        "e-"
    };
    for (size_t ii = 0; ii < sizeof(SUPPORTED_PARTICLE_KINDS) / sizeof(SUPPORTED_PARTICLE_KINDS[0]); ++ii) {
        if (kind == SUPPORTED_PARTICLE_KINDS[ii]) {
            return true;
        }
    }
    return false;
}

std::string particleFamilyPrefixFromIonMass(double ion_mass_u) {
    return ion_mass_u >= 1.5 ? "D" : "H";
}

std::string normalizedProcessIdSuffix(const std::string& projectile_kind) {
    std::string suffix;
    for (size_t ii = 0; ii < projectile_kind.size(); ++ii) {
        const char character = projectile_kind[ii];
        if (std::isalnum(static_cast<unsigned char>(character))) {
            suffix.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(character))));
            continue;
        }
        if (character == '+') {
            suffix += "plus";
            continue;
        }
        if (character == '-') {
            suffix += "minus";
        }
    }
    return suffix;
}

SimulationParameters::CrossSectionProcessProductDefinition makeLegacyProductDefinition(
    const std::string& particle_kind,
    uint count,
    const std::string& speed_class) {
    SimulationParameters::CrossSectionProcessProductDefinition product;
    product.particleKind = particle_kind;
    product.speedClass = speed_class;
    product.count = count;
    return product;
}

bool isSupportedProductSpeedClass(const std::string& speed_class) {
    return speed_class == "fast" || speed_class == "slow";
}

std::string inferLegacyProjectileKind(const std::string& reaction_id,
                                      double ion_mass_u,
                                      const std::string& context) {
    const std::string family = particleFamilyPrefixFromIonMass(ion_mass_u);
    if (reaction_id == "negative_ion_single_stripping" ||
        reaction_id == "negative_ion_double_stripping") {
        return family + "-";
    }
    if (reaction_id == "neutral_projectile_ionization") {
        return family + "0";
    }
    if (reaction_id == "positive_ion_charge_exchange") {
        return family + "+";
    }
    if (reaction_id == "background_gas_ionization") {
        throw runtime_error(
            context +
            ".projectileKind is required when migrating background_gas_ionization because the legacy process applies to multiple projectile species");
    }
    throw runtime_error(context + ".processId is not supported: " + reaction_id);
}

std::string inferLegacyProjectileFate(const std::string& reaction_id,
                                      const std::string& context) {
    if (reaction_id == "background_gas_ionization") {
        return "survive";
    }
    if (reaction_id == "negative_ion_single_stripping" ||
        reaction_id == "negative_ion_double_stripping" ||
        reaction_id == "neutral_projectile_ionization" ||
        reaction_id == "positive_ion_charge_exchange") {
        return "consume";
    }
    throw runtime_error(context + ".processId is not supported: " + reaction_id);
}

std::vector<SimulationParameters::CrossSectionProcessProductDefinition> inferLegacyProducts(
    const std::string& reaction_id,
    double ion_mass_u,
    const std::string& context) {
    const std::string family = particleFamilyPrefixFromIonMass(ion_mass_u);
    if (reaction_id == "negative_ion_single_stripping") {
        return {
            makeLegacyProductDefinition(family + "0", 1U, "fast"),
            makeLegacyProductDefinition("e-", 1U, "fast"),
        };
    }
    if (reaction_id == "negative_ion_double_stripping") {
        return {
            makeLegacyProductDefinition(family + "+", 1U, "fast"),
            makeLegacyProductDefinition("e-", 2U, "fast"),
        };
    }
    if (reaction_id == "background_gas_ionization") {
        return {
            makeLegacyProductDefinition(family + "2+", 1U, "slow"),
            makeLegacyProductDefinition("e-", 1U, "slow"),
        };
    }
    if (reaction_id == "neutral_projectile_ionization") {
        return {
            makeLegacyProductDefinition(family + "+", 1U, "fast"),
            makeLegacyProductDefinition("e-", 1U, "fast"),
        };
    }
    if (reaction_id == "positive_ion_charge_exchange") {
        return {
            makeLegacyProductDefinition(family + "2+", 1U, "slow"),
            makeLegacyProductDefinition(family + "0", 1U, "fast"),
        };
    }
    throw runtime_error(context + ".processId is not supported: " + reaction_id);
}

SimulationParameters::CrossSectionProcessProductDefinition parseCrossSectionProcessProductDefinition(
    const json& product,
    size_t process_index,
    size_t product_index) {
    const std::string context = "physics.reactions[" + std::to_string(process_index) + "].products[" +
                                std::to_string(product_index) + "]";
    if (!product.is_object()) {
        throw runtime_error(context + " must be an object");
    }
    if (!product.contains("particleKind") || !product.at("particleKind").is_string()) {
        throw runtime_error(context + ".particleKind must be a string");
    }

    SimulationParameters::CrossSectionProcessProductDefinition parsed;
    parsed.particleKind = product.at("particleKind").get<std::string>();
    if (!isSupportedParticleKind(parsed.particleKind)) {
        throw runtime_error(context + ".particleKind is not supported: " + parsed.particleKind);
    }

    if (!product.contains("speedClass") || !product.at("speedClass").is_string()) {
        throw runtime_error(context + ".speedClass must be a string");
    }
    parsed.speedClass = product.at("speedClass").get<std::string>();
    if (!isSupportedProductSpeedClass(parsed.speedClass)) {
        throw runtime_error(context + ".speedClass must be either fast or slow");
    }

    if (product.contains("count")) {
        if (!product.at("count").is_number_unsigned()) {
            throw runtime_error(context + ".count must be an unsigned integer");
        }
        parsed.count = product.at("count").get<uint>();
        if (parsed.count == 0U) {
            throw runtime_error(context + ".count must be >= 1");
        }
    }

    if (product.contains("chargeState") && !product.at("chargeState").is_number()) {
        throw runtime_error(context + ".chargeState must be numeric when provided");
    }
    if (product.contains("massU") && !product.at("massU").is_number()) {
        throw runtime_error(context + ".massU must be numeric when provided");
    }

    return parsed;
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

SimulationParameters::ParticleTypeDefinition parseParticleTypeDefinition(
    const json& particle_type,
    size_t index) {
    const std::string context = "particleTypes[" + std::to_string(index) + "]";
    if (!particle_type.is_object()) {
        throw runtime_error(context + " must be an object");
    }
    if (!particle_type.contains("id") || !particle_type.at("id").is_string()) {
        throw runtime_error(context + ".id must be a string");
    }
    if (!particle_type.contains("kind") || !particle_type.at("kind").is_string()) {
        throw runtime_error(context + ".kind must be a string");
    }
    if (!particle_type.contains("chargeState") || !particle_type.at("chargeState").is_number()) {
        throw runtime_error(context + ".chargeState must be numeric");
    }
    if (!particle_type.contains("massU") || !particle_type.at("massU").is_number()) {
        throw runtime_error(context + ".massU must be numeric");
    }

    SimulationParameters::ParticleTypeDefinition parsed;
    parsed.id = particle_type.at("id").get<std::string>();
    parsed.name = particle_type.value("name", parsed.id);
    parsed.kind = particle_type.at("kind").get<std::string>();
    parsed.chargeState = particle_type.at("chargeState").get<double>();
    parsed.massU = particle_type.at("massU").get<double>();
    parsed.sourceable = particle_type.value("sourceable", false);
    return parsed;
}

SimulationParameters::ParticleSourceDefinition parseParticleSourceDefinition(
    const json& particle_source,
    size_t index) {
    const std::string context = "particleSources[" + std::to_string(index) + "]";
    if (!particle_source.is_object()) {
        throw runtime_error(context + " must be an object");
    }
    if (!particle_source.contains("id") || !particle_source.at("id").is_string()) {
        throw runtime_error(context + ".id must be a string");
    }
    if (!particle_source.contains("particleTypeId") || !particle_source.at("particleTypeId").is_string()) {
        throw runtime_error(context + ".particleTypeId must be a string");
    }
    if (!particle_source.contains("kind") || !particle_source.at("kind").is_string()) {
        throw runtime_error(context + ".kind must be a string");
    }
    if (!particle_source.contains("sourceModel") || !particle_source.at("sourceModel").is_string()) {
        throw runtime_error(context + ".sourceModel must be a string");
    }
    if (particle_source.at("sourceModel").get<std::string>() != "uniform") {
        throw runtime_error(context + ".sourceModel must be uniform in the current C++ runtime");
    }
    if (!particle_source.contains("chargeState") || !particle_source.at("chargeState").is_number()) {
        throw runtime_error(context + ".chargeState must be numeric");
    }
    if (!particle_source.contains("massU") || !particle_source.at("massU").is_number()) {
        throw runtime_error(context + ".massU must be numeric");
    }
    if (!particle_source.contains("particleCount") || !particle_source.at("particleCount").is_number_unsigned()) {
        throw runtime_error(context + ".particleCount must be an unsigned integer");
    }
    if (!particle_source.contains("currentDensityAm2") || !particle_source.at("currentDensityAm2").is_number()) {
        throw runtime_error(context + ".currentDensityAm2 must be numeric");
    }
    if (!particle_source.contains("axialEnergyEV") || !particle_source.at("axialEnergyEV").is_number()) {
        throw runtime_error(context + ".axialEnergyEV must be numeric");
    }
    if (!hasObject(particle_source, "uniform")) {
        throw runtime_error(context + ".uniform object is required");
    }

    const json& uniform = particle_source.at("uniform");
    SimulationParameters::ParticleSourceDefinition parsed;
    parsed.id = particle_source.at("id").get<std::string>();
    parsed.name = particle_source.value("name", parsed.id);
    parsed.particleTypeId = particle_source.at("particleTypeId").get<std::string>();
    parsed.kind = particle_source.at("kind").get<std::string>();
    parsed.sourceModel = particle_source.at("sourceModel").get<std::string>();
    parsed.chargeState = particle_source.at("chargeState").get<double>();
    parsed.massU = particle_source.at("massU").get<double>();
    parsed.particleCount = particle_source.at("particleCount").get<uint>();
    parsed.currentDensityAm2 = particle_source.at("currentDensityAm2").get<double>();
    parsed.perpendicularTemperatureEV = particle_source.value("perpendicularTemperatureEV", 0.0);
    parsed.parallelTemperatureEV = particle_source.value("parallelTemperatureEV", 0.0);
    parsed.axialEnergyEV = particle_source.at("axialEnergyEV").get<double>();
    parsed.centerMeters = requireFixedNumberArray(uniform, "centerMeters", 3, context + ".uniform");
    parsed.mainDirection = requireFixedNumberArray(uniform, "mainDirection", 3, context + ".uniform");
    parsed.inPlaneReferenceDirection = requireFixedNumberArray(
        uniform,
        "inPlaneReferenceDirection",
        3,
        context + ".uniform");
    if (!uniform.contains("widthMeters") || !uniform.at("widthMeters").is_number()) {
        throw runtime_error(context + ".uniform.widthMeters must be numeric");
    }
    if (!uniform.contains("heightMeters") || !uniform.at("heightMeters").is_number()) {
        throw runtime_error(context + ".uniform.heightMeters must be numeric");
    }
    parsed.widthMeters = uniform.at("widthMeters").get<double>();
    parsed.heightMeters = uniform.at("heightMeters").get<double>();
    if (parsed.widthMeters <= 0.0 || parsed.heightMeters <= 0.0) {
        throw runtime_error(context + ".uniform.widthMeters and heightMeters must be > 0");
    }
    return parsed;
}

SimulationParameters::MagneticFieldDefinition parseMagneticFieldDefinition(
    const json& field,
    size_t index) {
    const std::string context = "externalMagneticField.fields[" + std::to_string(index) + "]";
    if (!field.is_object()) {
        throw runtime_error(context + " must be an object");
    }
    if (!field.contains("name") || !field.at("name").is_string()) {
        throw runtime_error(context + ".name must be a string");
    }
    if (!field.contains("sourceType") || !field.at("sourceType").is_string()) {
        throw runtime_error(context + ".sourceType must be a string");
    }

    SimulationParameters::MagneticFieldDefinition parsed;
    parsed.name = field.at("name").get<std::string>();
    parsed.sourceType = field.at("sourceType").get<std::string>();
    parsed.scale = field.value("scale", 1.0);
    if (parsed.scale < 0.0) {
        throw runtime_error(context + ".scale must be >= 0");
    }

    if (parsed.sourceType == "constant") {
        if (!hasObject(field, "constantValue")) {
            throw runtime_error(context + ".constantValue object is required for sourceType='constant'");
        }
        const json& constant = field.at("constantValue");
        if (!constant.contains("bx") || !constant.at("bx").is_number() ||
            !constant.contains("by") || !constant.at("by").is_number() ||
            !constant.contains("bz") || !constant.at("bz").is_number()) {
            throw runtime_error(context + ".constantValue must contain numeric bx, by, bz");
        }
        parsed.constantValue[0] = constant.at("bx").get<double>();
        parsed.constantValue[1] = constant.at("by").get<double>();
        parsed.constantValue[2] = constant.at("bz").get<double>();
    } else if (parsed.sourceType == "file") {
        if (!field.contains("filePath") || !field.at("filePath").is_string()) {
            throw runtime_error(context + ".filePath must be a string for sourceType='file'");
        }
        parsed.filePath = field.at("filePath").get<std::string>();
    } else {
        throw runtime_error(context + ".sourceType must be one of constant or file");
    }

    return parsed;
}

SimulationParameters::CrossSectionProcessDefinition parseCrossSectionProcessDefinition(
    const json& process,
    size_t index,
    double ion_mass_u) {
    const std::string context = "physics.reactions[" + std::to_string(index) + "]";
    if (!process.is_object()) {
        throw runtime_error(context + " must be an object");
    }
    if (!process.contains("processId") && !process.contains("reactionId")) {
        throw runtime_error(context + ".processId must be a string");
    }
    if (!process.contains("sourcePath") || !process.at("sourcePath").is_string()) {
        throw runtime_error(context + ".sourcePath must be a string");
    }
    if (!process.contains("fitDegree") || !process.at("fitDegree").is_number_unsigned()) {
        throw runtime_error(context + ".fitDegree must be an unsigned integer");
    }
    if (!process.contains("coefficients") || !process.at("coefficients").is_array()) {
        throw runtime_error(context + ".coefficients must be an array");
    }

    SimulationParameters::CrossSectionProcessDefinition parsed;
    parsed.processId = process.contains("processId")
                           ? process.at("processId").get<std::string>()
                           : process.at("reactionId").get<std::string>();
    if (process.contains("reactionId") && process.at("reactionId").is_string() &&
        process.contains("processId") && process.at("processId").is_string() &&
        process.at("reactionId").get<std::string>() != parsed.processId) {
        throw runtime_error(context + ".reactionId must match processId when both are provided");
    }
    parsed.name = process.value("name", parsed.processId);
    parsed.sourcePath = process.at("sourcePath").get<std::string>();
    parsed.fitDegree = process.at("fitDegree").get<uint>();
    if (parsed.fitDegree > 6U) {
        throw runtime_error(context + ".fitDegree must be <= 6");
    }
    parsed.scaleEnergyByIonMass = process.value("scaleEnergyByIonMass", true);
    parsed.minimumEnergyEV = process.value("minimumEnergyEV", 0.0);
    parsed.maximumEnergyEV = process.value("maximumEnergyEV", -1.0);
    parsed.coefficients = requireFixedNumberArray(
        process,
        "coefficients",
        static_cast<size_t>(parsed.fitDegree + 1U),
        context);

    if (process.contains("projectileKind")) {
        if (!process.at("projectileKind").is_string()) {
            throw runtime_error(context + ".projectileKind must be a string");
        }
        parsed.projectileKind = process.at("projectileKind").get<std::string>();
        if (!isSupportedParticleKind(parsed.projectileKind) || parsed.projectileKind == "e-") {
            throw runtime_error(context + ".projectileKind must be a supported heavy-particle kind");
        }
    } else {
        parsed.projectileKind = inferLegacyProjectileKind(parsed.processId, ion_mass_u, context);
    }

    if (process.contains("projectileFate")) {
        if (!process.at("projectileFate").is_string()) {
            throw runtime_error(context + ".projectileFate must be a string");
        }
        parsed.projectileFate = process.at("projectileFate").get<std::string>();
    } else {
        parsed.projectileFate = inferLegacyProjectileFate(parsed.processId, context);
    }
    if (parsed.projectileFate != "consume" && parsed.projectileFate != "survive") {
        throw runtime_error(context + ".projectileFate must be either consume or survive");
    }

    if (process.contains("products")) {
        if (!process.at("products").is_array()) {
            throw runtime_error(context + ".products must be an array");
        }
        for (json::const_iterator it = process.at("products").begin();
             it != process.at("products").end();
             ++it) {
            parsed.products.push_back(
                parseCrossSectionProcessProductDefinition(*it, index, parsed.products.size()));
        }
    } else {
        parsed.products = inferLegacyProducts(parsed.processId, ion_mass_u, context);
    }

    if (parsed.minimumEnergyEV < 0.0) {
        throw runtime_error(context + ".minimumEnergyEV must be >= 0");
    }
    if (parsed.maximumEnergyEV > 0.0 && parsed.maximumEnergyEV < parsed.minimumEnergyEV) {
        throw runtime_error(context + ".maximumEnergyEV must be >= minimumEnergyEV when provided");
    }

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
    INITIAL_PLASMA_MAX_Z(0.0),
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
      SPLIT_DOMAIN(0U),
      JTOLERANCE(0.0),
      ALPHA_COEFF(0.0),
      T_POSITIVE(0.0),
      N_SOLIDS(0U),
      MGSOLVER(0U),
      SHIELD_MODEL(0U),
    BICGSTAB_EPS(1.0e-4),
    BICGSTAB_MAX_ITERATIONS(10000U),
    BICGSTAB_NEWTON_EPS(1.0e-4),
    BICGSTAB_NEWTON_MAX_ITERATIONS(10U),
    BICGSTAB_GLOBALLY_CONVERGENT_NEWTON(1U),
    MG_LEVELS(1U),
    MG_TOLERANCE(1.0e-4),
    MG_MAX_CYCLES(100U),
    MG_GAMMA(1U),
    MG_PRE_SMOOTH(5U),
    MG_POST_SMOOTH(5U),
    MG_COARSE_RELAXATION(1.7),
    MG_COARSE_MAX_ITERATIONS(10000U),
    MG_LOCAL_PLASMA_MAX_ITERATIONS(1U),
      EXT_GAP(0.0),
      DOMAIN_X_SIZE(-1.0),
      DOMAIN_Y_SIZE(-1.0),
      DOMAIN_Z_SIZE(-1.0),
      DOMAIN_Z_START(0.0),
      EGEXTJ(0.0),
      domain_ii(0U),
      GEOMETRY_SOURCE_MODE(),
      GENERATED_GEOMETRY_SOLIDS(),
    PARTICLE_TYPES(),
    PARTICLE_SOURCES(),
      EXPLICIT_BOUNDARY_CONDITIONS(),
      MAGNETIC_FIELD_SOURCE_MODE("none"),
      MAGNETIC_FIELD_DIRECTORY(),
      MAGNETIC_FIELD_FILE(),
            MAGNETIC_FIELDS(),
            CROSS_SECTION_PROCESSES(),
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
        OUTPUT_LOGGING_STRUCTURED_LOG_FILE("run.log"),
        OUTPUT_ITERATION_ENABLED(1U),
        OUTPUT_ITERATION_EVERY_N_ITERATIONS(1U),
        OUTPUT_ITERATION_EXPORT_PLANE_DIAGNOSTICS(1U),
        OUTPUT_ITERATION_EXPORT_SIMULATION_STATE(0U),
        OUTPUT_ITERATION_EXPORT_TRACED_PARTICLES(0U),
            OUTPUT_ITERATION_PLANE_Z_POSITIONS(),
    DIAGNOSTIC_SAMPLE_Z_POSITIONS(),
    DIAGNOSTIC_SUMMARY_Z_POSITION(-1.0),
    DIAGNOSTIC_EMITTER_EXPORT_Z_POSITION(-1.0),
    DIAGNOSTIC_TRANSMISSION_PLANE_Z_POSITION(-1.0),
    DIAGNOSTIC_APERTURE_RADIUS(7.0e-3),
    DIAGNOSTIC_WRITE_PER_SPECIES_DIAGNOSTICS(1U),
    // Off by default: the per-species grid-power passes cost 7 extra full analyses of the
    // particle database, and the <tag>_grid_power_breakdown.txt file now resolves species
    // and generation from the single PARTICLE_ALL pass. Kept as an opt-in cross-check.
    DIAGNOSTIC_WRITE_PER_SPECIES_GRID_POWER(0U),
    DIAGNOSTIC_WRITE_PER_SPECIES_PLOTS(1U),
    DIAGNOSTIC_WRITE_NEGATIVE_ION_SUMMARY(1U),
    // Off by default: writing a per-triangle map costs an extra pass over every impact
    // plus a VTK file, and is only meaningful once the impact statistics support it (a few
    // hundred macroparticles on a grid is Poisson noise, not a density map).
    DIAGNOSTIC_WRITE_POWER_DENSITY_MAP(0U),
    DIAGNOSTIC_GRID_POWER_RANGES(),
    DIAGNOSTIC_MENISCUS_PLOT() {
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

    PARTICLE_TYPES.clear();
    PARTICLE_SOURCES.clear();

    if (root.contains("particleTypes")) {
        if (!root.at("particleTypes").is_array()) {
            throw Error(ERROR_LOCATION, "particleTypes must be an array in the runtime JSON configuration");
        }
        for (json::const_iterator it = root.at("particleTypes").begin();
             it != root.at("particleTypes").end();
             ++it) {
            try {
                PARTICLE_TYPES.push_back(parseParticleTypeDefinition(*it, PARTICLE_TYPES.size()));
            } catch (const std::runtime_error& e) {
                throw Error(ERROR_LOCATION, e.what());
            }
        }
    }

    if (root.contains("particleSources") && root.at("particleSources").is_array()) {
        for (json::const_iterator it = root.at("particleSources").begin();
             it != root.at("particleSources").end();
             ++it) {
            try {
                PARTICLE_SOURCES.push_back(parseParticleSourceDefinition(*it, PARTICLE_SOURCES.size()));
            } catch (const std::runtime_error& e) {
                throw Error(ERROR_LOCATION, e.what());
            }
        }
    } else if (hasObject(root, "particleSources") && hasObject(root.at("particleSources"), "negativeIonBeam")) {
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
    } else if (root.contains("particleSources")) {
        throw Error(
            ERROR_LOCATION,
            "particleSources must be either an array of explicit sources or the legacy negativeIonBeam object");
    }

    if (!PARTICLE_SOURCES.empty()) {
        const SimulationParameters::ParticleSourceDefinition* primary_source = nullptr;
        double aggregate_negative_current_density = 0.0;
        uint aggregate_particle_count = 0U;
        std::string heavy_negative_kind;
        for (std::vector<SimulationParameters::ParticleSourceDefinition>::const_iterator it = PARTICLE_SOURCES.begin();
             it != PARTICLE_SOURCES.end();
             ++it) {
            if (primary_source == nullptr && it->kind != "e-") {
                primary_source = &(*it);
            }
            if (it->chargeState < 0.0 && it->kind != "e-") {
                if (!heavy_negative_kind.empty() && heavy_negative_kind != it->kind) {
                    throw Error(
                        ERROR_LOCATION,
                        "Mixed heavy-ion source kinds are not yet supported by the current C++ runtime: " +
                            heavy_negative_kind + ", " + it->kind);
                }
                heavy_negative_kind = it->kind;
                aggregate_negative_current_density += it->currentDensityAm2;
            }
            aggregate_particle_count += it->particleCount;
        }
        if (primary_source == nullptr) {
            primary_source = &PARTICLE_SOURCES.front();
        }

        M_IONS = primary_source->massU;
        Q_IONS = primary_source->chargeState;
        J_ION = (aggregate_negative_current_density > 0.0)
                    ? aggregate_negative_current_density
                    : primary_source->currentDensityAm2;
        TPERP = primary_source->perpendicularTemperatureEV;
        TPAR = primary_source->parallelTemperatureEV;
        E0_Z = primary_source->axialEnergyEV;
        N_PARTICLES = aggregate_particle_count;
    }

    if (hasObject(root, "simulation")) {
        const json& simulation = root.at("simulation");
        const bool has_explicit_particle_count = simulation.contains("particleCount");
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
            if (hasObject(solver, "bicgstab")) {
                const json& bicgstab = solver.at("bicgstab");
                if (bicgstab.contains("eps")) {
                    BICGSTAB_EPS = bicgstab.at("eps").get<double>();
                }
                if (bicgstab.contains("maxIterations")) {
                    BICGSTAB_MAX_ITERATIONS = bicgstab.at("maxIterations").get<uint>();
                }
                if (bicgstab.contains("newtonEps")) {
                    BICGSTAB_NEWTON_EPS = bicgstab.at("newtonEps").get<double>();
                }
                if (bicgstab.contains("newtonMaxIterations")) {
                    BICGSTAB_NEWTON_MAX_ITERATIONS = bicgstab.at("newtonMaxIterations").get<uint>();
                }
                if (bicgstab.contains("globallyConvergentNewton")) {
                    BICGSTAB_GLOBALLY_CONVERGENT_NEWTON =
                        bicgstab.at("globallyConvergentNewton").get<bool>() ? 1U : 0U;
                }
            }
            if (hasObject(solver, "multigrid")) {
                const json& multigrid = solver.at("multigrid");
                if (multigrid.contains("levels")) {
                    MG_LEVELS = multigrid.at("levels").get<uint>();
                }
                if (multigrid.contains("mgTolerance")) {
                    MG_TOLERANCE = multigrid.at("mgTolerance").get<double>();
                }
                if (multigrid.contains("maxCycles")) {
                    MG_MAX_CYCLES = multigrid.at("maxCycles").get<uint>();
                }
                if (multigrid.contains("gamma")) {
                    MG_GAMMA = multigrid.at("gamma").get<uint>();
                }
                if (multigrid.contains("preSmooth")) {
                    MG_PRE_SMOOTH = multigrid.at("preSmooth").get<uint>();
                }
                if (multigrid.contains("postSmooth")) {
                    MG_POST_SMOOTH = multigrid.at("postSmooth").get<uint>();
                }
                if (multigrid.contains("coarseRelaxation")) {
                    MG_COARSE_RELAXATION = multigrid.at("coarseRelaxation").get<double>();
                }
                if (multigrid.contains("coarseMaxIterations")) {
                    MG_COARSE_MAX_ITERATIONS = multigrid.at("coarseMaxIterations").get<uint>();
                }
                if (multigrid.contains("localPlasmaMaxIterations")) {
                    MG_LOCAL_PLASMA_MAX_ITERATIONS = multigrid.at("localPlasmaMaxIterations").get<uint>();
                }
            }
            if (solver.contains("plasmaModel")) {
                SHIELD_MODEL = shieldModelFromType(solver.at("plasmaModel").get<string>());
            } else if (solver.contains("shieldModel")) {
                SHIELD_MODEL = shieldModelFromType(solver.at("shieldModel").get<string>());
            }
            if (solver.contains("initialPlasmaMaxZMeters")) {
                INITIAL_PLASMA_MAX_Z = solver.at("initialPlasmaMaxZMeters").get<double>();
            }
            if (SHIELD_MODEL == 0U) {
                if (solver.contains("positiveIonTemperatureEV")) {
                    T_POSITIVE = solver.at("positiveIonTemperatureEV").get<double>();
                }
                if (solver.contains("plasmaPotentialVolts")) {
                    U_PLASMA = solver.at("plasmaPotentialVolts").get<double>();
                }
            } else {
                if (solver.contains("tanhWidthEV")) {
                    T_POSITIVE = solver.at("tanhWidthEV").get<double>();
                } else if (solver.contains("positiveIonTemperatureEV")) {
                    T_POSITIVE = solver.at("positiveIonTemperatureEV").get<double>();
                }
                if (solver.contains("meniscusVoltageVolts")) {
                    U_PLASMA = solver.at("meniscusVoltageVolts").get<double>();
                } else if (solver.contains("plasmaPotentialVolts")) {
                    U_PLASMA = solver.at("plasmaPotentialVolts").get<double>();
                }
            }
        }

        if (hasObject(simulation, "spaceCharge")) {
            const json& space_charge = simulation.at("spaceCharge");
            if (space_charge.contains("alphaCoeff")) {
                ALPHA_COEFF = space_charge.at("alphaCoeff").get<double>();
            }
        }

        if (hasObject(simulation, "convergence")) {
            const json& convergence = simulation.at("convergence");
            if (convergence.contains("currentDensityTolerance")) {
                JTOLERANCE = convergence.at("currentDensityTolerance").get<double>();
            }
        }

        if (!PARTICLE_SOURCES.empty() && has_explicit_particle_count) {
            uint derived_particle_count = 0U;
            for (std::vector<SimulationParameters::ParticleSourceDefinition>::const_iterator it = PARTICLE_SOURCES.begin();
                 it != PARTICLE_SOURCES.end();
                 ++it) {
                derived_particle_count += it->particleCount;
            }
            if (N_PARTICLES != derived_particle_count) {
                throw Error(
                    ERROR_LOCATION,
                    "simulation.particleCount must match the sum of particleSources[*].particleCount");
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
        MAGNETIC_FIELDS.clear();
        if (magnetic_field.contains("enabled")) {
            B_ISON = magnetic_field.at("enabled").get<bool>() ? 1U : 0U;
        }
        MAGNETIC_FIELD_DIRECTORY = magnetic_field.value("directory", string());
        if (B_ISON != 0U) {
            if (!magnetic_field.contains("fields") || !magnetic_field.at("fields").is_array()) {
                throw Error(ERROR_LOCATION,
                            "externalMagneticField.fields array is required when the magnetic field is enabled");
            }

            for (json::const_iterator it = magnetic_field.at("fields").begin();
                 it != magnetic_field.at("fields").end();
                 ++it) {
                try {
                    MAGNETIC_FIELDS.push_back(parseMagneticFieldDefinition(*it, MAGNETIC_FIELDS.size()));
                } catch (const std::runtime_error& e) {
                    throw Error(ERROR_LOCATION, e.what());
                }
            }

            MAGNETIC_FIELD_SOURCE_MODE = MAGNETIC_FIELDS.empty() ? "none" : "list";
        } else {
            MAGNETIC_FIELD_SOURCE_MODE = "none";
            MAGNETIC_FIELDS.clear();
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
        if (hasObject(outputs, "iteration")) {
            const json& iteration = outputs.at("iteration");
            if (iteration.contains("enabled")) {
                OUTPUT_ITERATION_ENABLED = iteration.at("enabled").get<bool>() ? 1U : 0U;
            }
            if (iteration.contains("everyNIterations")) {
                OUTPUT_ITERATION_EVERY_N_ITERATIONS = iteration.at("everyNIterations").get<uint>();
                if (OUTPUT_ITERATION_EVERY_N_ITERATIONS == 0U) {
                    throw Error(ERROR_LOCATION, "outputs.iteration.everyNIterations must be >= 1");
                }
            }
            if (iteration.contains("exportPlaneDiagnostics")) {
                OUTPUT_ITERATION_EXPORT_PLANE_DIAGNOSTICS =
                    iteration.at("exportPlaneDiagnostics").get<bool>() ? 1U : 0U;
            }
            if (iteration.contains("exportSimulationState")) {
                OUTPUT_ITERATION_EXPORT_SIMULATION_STATE =
                    iteration.at("exportSimulationState").get<bool>() ? 1U : 0U;
            }
            if (iteration.contains("exportTracedParticles")) {
                OUTPUT_ITERATION_EXPORT_TRACED_PARTICLES =
                    iteration.at("exportTracedParticles").get<bool>() ? 1U : 0U;
            }
            if (iteration.contains("planeZPositionsMeters")) {
                if (!iteration.at("planeZPositionsMeters").is_array()) {
                    throw Error(ERROR_LOCATION, "outputs.iteration.planeZPositionsMeters must be an array");
                }
                OUTPUT_ITERATION_PLANE_Z_POSITIONS.clear();
                const double domain_z_min = DOMAIN_Z_START;
                const double domain_z_max = DOMAIN_Z_START + DOMAIN_Z_SIZE;
                for (json::const_iterator it = iteration.at("planeZPositionsMeters").begin();
                     it != iteration.at("planeZPositionsMeters").end();
                     ++it) {
                    if (!it->is_number()) {
                        throw Error(
                            ERROR_LOCATION,
                            "outputs.iteration.planeZPositionsMeters entries must be numeric"
                        );
                    }
                    const double plane = it->get<double>();
                    if (plane < domain_z_min || plane > domain_z_max) {
                        throw Error(
                            ERROR_LOCATION,
                            "outputs.iteration.planeZPositionsMeters value " + std::to_string(plane) +
                                " is outside geometry.domain z range"
                        );
                    }
                    OUTPUT_ITERATION_PLANE_Z_POSITIONS.push_back(plane);
                }
            }
        }
    }

    if (hasObject(root, "diagnostics")) {
        const json& diagnostics = root.at("diagnostics");

        if (hasObject(diagnostics, "planes")) {
            const json& planes = diagnostics.at("planes");
            if (planes.contains("sampleZPositionsMeters")) {
                if (!planes.at("sampleZPositionsMeters").is_array()) {
                    throw Error(ERROR_LOCATION, "diagnostics.planes.sampleZPositionsMeters must be an array");
                }
                DIAGNOSTIC_SAMPLE_Z_POSITIONS.clear();
                for (json::const_iterator it = planes.at("sampleZPositionsMeters").begin();
                     it != planes.at("sampleZPositionsMeters").end();
                     ++it) {
                    if (!it->is_number()) {
                        throw Error(ERROR_LOCATION,
                                    "diagnostics.planes.sampleZPositionsMeters must contain only numeric values");
                    }
                    DIAGNOSTIC_SAMPLE_Z_POSITIONS.push_back(it->get<double>());
                }
            }
            if (planes.contains("summaryZPositionMeters")) {
                DIAGNOSTIC_SUMMARY_Z_POSITION = planes.at("summaryZPositionMeters").get<double>();
            }
            if (planes.contains("emitterExportZPositionMeters")) {
                DIAGNOSTIC_EMITTER_EXPORT_Z_POSITION =
                    planes.at("emitterExportZPositionMeters").get<double>();
            }
        }

        if (hasObject(diagnostics, "species")) {
            const json& species = diagnostics.at("species");
            if (species.contains("writePerSpeciesDiagnostics")) {
                DIAGNOSTIC_WRITE_PER_SPECIES_DIAGNOSTICS =
                    species.at("writePerSpeciesDiagnostics").get<bool>() ? 1U : 0U;
            }
            if (species.contains("writePerSpeciesGridPower")) {
                DIAGNOSTIC_WRITE_PER_SPECIES_GRID_POWER =
                    species.at("writePerSpeciesGridPower").get<bool>() ? 1U : 0U;
            }
            if (species.contains("writePerSpeciesPlots")) {
                DIAGNOSTIC_WRITE_PER_SPECIES_PLOTS =
                    species.at("writePerSpeciesPlots").get<bool>() ? 1U : 0U;
            }
            if (species.contains("writeNegativeIonSummary")) {
                DIAGNOSTIC_WRITE_NEGATIVE_ION_SUMMARY =
                    species.at("writeNegativeIonSummary").get<bool>() ? 1U : 0U;
            }
        }

        if (hasObject(diagnostics, "gridPower")) {
            const json& grid_power = diagnostics.at("gridPower");
            if (grid_power.contains("writePowerDensityMap")) {
                DIAGNOSTIC_WRITE_POWER_DENSITY_MAP =
                    grid_power.at("writePowerDensityMap").get<bool>() ? 1U : 0U;
            }
            if (grid_power.contains("ranges")) {
                if (!grid_power.at("ranges").is_array()) {
                    throw Error(ERROR_LOCATION, "diagnostics.gridPower.ranges must be an array");
                }
                DIAGNOSTIC_GRID_POWER_RANGES.clear();
                std::set<int> seen_grid_power_ids;
                for (json::const_iterator it = grid_power.at("ranges").begin();
                     it != grid_power.at("ranges").end();
                     ++it) {
                    try {
                        DiagnosticGridRangeDefinition definition =
                            parseDiagnosticGridRangeDefinition(*it, DIAGNOSTIC_GRID_POWER_RANGES.size());
                        if (!seen_grid_power_ids.insert(definition.id).second) {
                            throw runtime_error(
                                string("diagnostics.gridPower.ranges contains duplicate id: ") +
                                std::to_string(definition.id));
                        }
                        DIAGNOSTIC_GRID_POWER_RANGES.push_back(definition);
                    } catch (const std::runtime_error& e) {
                        throw Error(ERROR_LOCATION, e.what());
                    }
                }
            }
        }

        if (hasObject(diagnostics, "summary")) {
            const json& summary = diagnostics.at("summary");
            if (summary.contains("apertureRadiusMeters")) {
                DIAGNOSTIC_APERTURE_RADIUS = summary.at("apertureRadiusMeters").get<double>();
            }
            if (summary.contains("transmissionPlaneZPositionMeters")) {
                DIAGNOSTIC_TRANSMISSION_PLANE_Z_POSITION =
                    summary.at("transmissionPlaneZPositionMeters").get<double>();
            }
        }

        if (hasObject(diagnostics, "plots")) {
            const json& plots = diagnostics.at("plots");
            if (hasObject(plots, "meniscus")) {
                const json& meniscus = plots.at("meniscus");
                if (meniscus.contains("enabled")) {
                    DIAGNOSTIC_MENISCUS_PLOT.enabled = meniscus.at("enabled").get<bool>();
                }
                if (meniscus.contains("zMinMeters")) {
                    DIAGNOSTIC_MENISCUS_PLOT.zMinMeters = meniscus.at("zMinMeters").get<double>();
                }
                if (meniscus.contains("zMaxMeters")) {
                    DIAGNOSTIC_MENISCUS_PLOT.zMaxMeters = meniscus.at("zMaxMeters").get<double>();
                }
                if (meniscus.contains("transverseMinMeters")) {
                    DIAGNOSTIC_MENISCUS_PLOT.transverseMinMeters =
                        meniscus.at("transverseMinMeters").get<double>();
                }
                if (meniscus.contains("transverseMaxMeters")) {
                    DIAGNOSTIC_MENISCUS_PLOT.transverseMaxMeters =
                        meniscus.at("transverseMaxMeters").get<double>();
                }
            }
        }
    }

    if (hasObject(root, "physics")) {
        const json& physics = root.at("physics");
        if (physics.contains("reactions")) {
            if (!physics.at("reactions").is_array()) {
                throw Error(ERROR_LOCATION, "physics.reactions must be an array when provided");
            }

            CROSS_SECTION_PROCESSES.clear();
            std::set<std::string> seen_process_ids;
            for (json::const_iterator it = physics.at("reactions").begin();
                 it != physics.at("reactions").end();
                 ++it) {
                try {
                    if (it->is_object() &&
                        it->value("reactionId", string()) == "background_gas_ionization" &&
                        !it->contains("projectileKind")) {
                        const std::string family = particleFamilyPrefixFromIonMass(M_IONS);
                        const std::string projectile_kinds[] = { family + "-", family + "0" };
                        for (size_t projectile_index = 0; projectile_index < 2U; ++projectile_index) {
                            json expanded_process = *it;
                            expanded_process["projectileKind"] = projectile_kinds[projectile_index];
                            CrossSectionProcessDefinition process_definition =
                                parseCrossSectionProcessDefinition(
                                    expanded_process,
                                    CROSS_SECTION_PROCESSES.size(),
                                    M_IONS);
                            process_definition.processId =
                                string("background_gas_ionization_") +
                                normalizedProcessIdSuffix(projectile_kinds[projectile_index]);
                            process_definition.name = process_definition.name +
                                                      " (" + projectile_kinds[projectile_index] + ")";
                            if (!seen_process_ids.insert(process_definition.processId).second) {
                                throw runtime_error(
                                    string("physics.reactions contains duplicate processId: ") +
                                    process_definition.processId);
                            }
                            CROSS_SECTION_PROCESSES.push_back(process_definition);
                        }
                        continue;
                    }

                    CrossSectionProcessDefinition process_definition =
                        parseCrossSectionProcessDefinition(*it, CROSS_SECTION_PROCESSES.size(), M_IONS);
                    if (!seen_process_ids.insert(process_definition.processId).second) {
                        throw runtime_error(
                            string("physics.reactions contains duplicate processId: ") +
                            process_definition.processId);
                    }
                    CROSS_SECTION_PROCESSES.push_back(process_definition);
                } catch (const std::runtime_error& e) {
                    throw Error(ERROR_LOCATION, e.what());
                }
            }
        }

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
        if (CROSS_SECTION_PROCESSES.empty()) {
            throw Error(ERROR_LOCATION,
                        "physics.reactions must contain at least one process when stripping is enabled");
        }
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

    clear_cross_section_processes();
    for (std::vector<CrossSectionProcessDefinition>::const_iterator it = CROSS_SECTION_PROCESSES.begin();
         it != CROSS_SECTION_PROCESSES.end();
         ++it) {
        configure_cross_section_process(it->processId,
                                        it->sourcePath,
                                        it->coefficients,
                                        it->fitDegree,
                                        it->scaleEnergyByIonMass,
                                        it->minimumEnergyEV,
                                        it->maximumEnergyEV);
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
    PARTICLE_TYPES.clear();
    PARTICLE_SOURCES.clear();
    EXPLICIT_BOUNDARY_CONDITIONS.clear();
    MAGNETIC_FIELD_SOURCE_MODE = "none";
    MAGNETIC_FIELD_DIRECTORY.clear();
    MAGNETIC_FIELD_FILE.clear();
    MAGNETIC_FIELDS.clear();
    CROSS_SECTION_PROCESSES.clear();
    clear_cross_section_processes();

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
    INITIAL_PLASMA_MAX_Z = 7.0e-3;
    N_PARTICLES = 150000U;

    EG_VOLTAGE = 8000.0;
    GG_VOLTAGE = 182000.0;
    REP_VOLTAGE = 356000.0;
    G3_VOLTAGE = 530000.0;
    G4_VOLTAGE = 704000.0;
    G5_VOLTAGE = 878000.0;

    MESH_SIZE = 0.0003;
    ITERATIONS = 5U;
    SPLIT_DOMAIN = 0U;
    JTOLERANCE = 1.0;
    ALPHA_COEFF = 0.3;
    T_POSITIVE = 0.8;
    MGSOLVER = 0U;
    SHIELD_MODEL = 0U;
    BICGSTAB_EPS = 1.0e-4;
    BICGSTAB_MAX_ITERATIONS = 10000U;
    BICGSTAB_NEWTON_EPS = 1.0e-4;
    BICGSTAB_NEWTON_MAX_ITERATIONS = 10U;
    BICGSTAB_GLOBALLY_CONVERGENT_NEWTON = 1U;
    MG_LEVELS = 1U;
    MG_TOLERANCE = 1.0e-4;
    MG_MAX_CYCLES = 100U;
    MG_GAMMA = 1U;
    MG_PRE_SMOOTH = 5U;
    MG_POST_SMOOTH = 5U;
    MG_COARSE_RELAXATION = 1.7;
    MG_COARSE_MAX_ITERATIONS = 10000U;
    MG_LOCAL_PLASMA_MAX_ITERATIONS = 1U;

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
    DIAGNOSTIC_SAMPLE_Z_POSITIONS.clear();
    DIAGNOSTIC_SUMMARY_Z_POSITION = -1.0;
    DIAGNOSTIC_EMITTER_EXPORT_Z_POSITION = -1.0;
    DIAGNOSTIC_TRANSMISSION_PLANE_Z_POSITION = -1.0;
    DIAGNOSTIC_APERTURE_RADIUS = 7.0e-3;
    DIAGNOSTIC_WRITE_PER_SPECIES_DIAGNOSTICS = 1U;
    DIAGNOSTIC_WRITE_PER_SPECIES_GRID_POWER = 0U;   // see the constructor for why
    DIAGNOSTIC_WRITE_PER_SPECIES_PLOTS = 1U;
    DIAGNOSTIC_WRITE_NEGATIVE_ION_SUMMARY = 1U;
    DIAGNOSTIC_WRITE_POWER_DENSITY_MAP = 0U;
    DIAGNOSTIC_GRID_POWER_RANGES.clear();
    DIAGNOSTIC_MENISCUS_PLOT = DiagnosticMeniscusPlotDefinition();
    N_SOLIDS = 0U;
    finalizeDerivedParameters();
}

void SimulationParameters::finalizeDerivedParameters() {
    G1_VOLTAGE = GG_VOLTAGE;
    G2_VOLTAGE = REP_VOLTAGE;
}
