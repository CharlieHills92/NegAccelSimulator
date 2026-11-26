/*
 * ParticleManager.cpp
 *
 *  Created on: Aug 04, 2025
 *      Author: GitHub Copilot (Refactored from ManageSimulation)
 */

#include "ParticleManager.h"
#include "SimulationParameters.h"
#include "globals.h"
#include "my_diagnostics.h"

#include <random>
#include <cstdlib>

#include "particledatabase.hpp"
#include "geometry.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "funct.h"
#include "particles.hpp"

#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <map>

using namespace std;

ParticleManager::ParticleManager() : 
    particles(nullptr), collisions(nullptr), part_kind(nullptr), 
    part_gen(nullptr), rgn(nullptr) {
    initializeDiagnostics();
}

ParticleManager::~ParticleManager() {
    resetParticles();
}

size_t ParticleManager::beamletSlotFor(int beamlet_index,
                                       int column,
                                       int row,
                                       double center_x,
                                       double center_y) {
    // Several authored sources (different species, say) can feed the same beamlet, so the
    // slot is keyed on the beamlet index and not on the source.
    for (size_t s = 0; s < beamlet_tags.size(); s++) {
        if (beamlet_tags[s].index == beamlet_index) return s;
    }
    BeamletTag tag;
    tag.index = beamlet_index;
    tag.column = column;
    tag.row = row;
    tag.centerX = center_x;
    tag.centerY = center_y;
    beamlet_tags.push_back(tag);
    return beamlet_tags.size() - 1;
}

void ParticleManager::tagParticleRange(size_t first, size_t last, size_t slot) {
    if (last <= first) return;
    if (particle_beamlet_slot.size() < last) {
        particle_beamlet_slot.resize(last, -1);
    }
    for (size_t i = first; i < last; i++) {
        particle_beamlet_slot[i] = static_cast<int>(slot);
    }
}

const BeamletTag* ParticleManager::beamletOfParticle(size_t index,
                                                     double birth_x,
                                                     double birth_y) const {
    if (beamlet_tags.empty()) return NULL;

    if (index < particle_beamlet_slot.size()) {
        const int slot = particle_beamlet_slot[index];
        if (slot >= 0 && static_cast<size_t>(slot) < beamlet_tags.size()) {
            return &beamlet_tags[slot];
        }
    }

    // Not emitted by a source: a secondary created during tracing. Attribute it to the
    // beamlet whose axis it was born closest to. That is the beamlet whose beam produced
    // it, except in the gap between two apertures, where either answer is defensible.
    size_t best = 0;
    double best_d2 = -1.0;
    for (size_t s = 0; s < beamlet_tags.size(); s++) {
        const double dx = birth_x - beamlet_tags[s].centerX;
        const double dy = birth_y - beamlet_tags[s].centerY;
        const double d2 = dx * dx + dy * dy;
        if (best_d2 < 0.0 || d2 < best_d2) {
            best_d2 = d2;
            best = s;
        }
    }
    return &beamlet_tags[best];
}

void ParticleManager::resetParticles() {
    delete particles;
    particles = nullptr;
    
    for (auto* pdb : particles_species) {
        delete pdb;
    }
    particles_species.clear();
    particles_kinds.clear();
    part_number.clear();
    
    delete rgn;
    rgn = nullptr;
}

void ParticleManager::setParticles(ParticleDataBase3D* pdb) {
    if (particles != nullptr && particles != pdb) {
        delete particles;
    }
    particles = pdb;
}

void ParticleManager::initializeDiagnostics() {
    diagnostics = {
        DIAG_X, DIAG_Y, DIAG_Z,
        DIAG_VX, DIAG_VY, DIAG_VZ,
        DIAG_CURR, DIAG_MASS, DIAG_CHARGE, DIAG_NO
    };
}

namespace {

/// The plasma grid, identified the way the configuration actually labels it.
///
/// Cases in this project carry no "role" field, so name is the primary key and role is
/// accepted when present. The final fallback is the solid that starts furthest upstream,
/// which is the plasma grid in every accelerator layout this code builds.
const SimulationParameters::GeometrySolidDefinition* findPlasmaGrid(
        const std::vector<SimulationParameters::GeometrySolidDefinition>& solids) {
    const SimulationParameters::GeometrySolidDefinition* upstream = NULL;
    for (size_t i = 0; i < solids.size(); i++) {
        const SimulationParameters::GeometrySolidDefinition& s = solids[i];
        if (s.zProfileMeters.size() < 2 || s.rProfileMeters.size() != s.zProfileMeters.size())
            continue;

        std::string name = s.name;
        std::transform(name.begin(), name.end(), name.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
        if (name == "PG" || name == "PLASMA GRID" || name == "PLASMAGRID" ||
            s.role == "plasmaGrid") {
            return &s;
        }
        if (upstream == NULL || s.zProfileMeters.front() < upstream->zProfileMeters.front()) {
            upstream = &s;
        }
    }
    return upstream;
}

}  // namespace

size_t ParticleManager::addBeamHalo(ParticleDataBase3D* pdb,
                                    const SimulationParameters& params,
                                    double core_current_A) {
    if (!params.getBeamHaloEnabled()) return 0U;

    if (!params.hasGeneratedGeometrySolids()) {
        ibsimu.message(1) << "WARNING: beam halo requested but the case has no generated "
                          << "geometry solids; no halo emitted." << endl;
        return 0U;
    }
    const SimulationParameters::GeometrySolidDefinition* pg =
        findPlasmaGrid(params.getGeneratedGeometrySolids());
    if (pg == NULL) {
        ibsimu.message(1) << "WARNING: beam halo requested but the plasma grid could not be "
                          << "identified; no halo emitted." << endl;
        return 0U;
    }

    // The downstream chamfer is the part of the bore profile after the throat: from the
    // narrowest radius to the downstream face. Halo particles are born there and aimed
    // along the local surface normal, following Pimazzoni et al., Fus. Eng. Des. 192 (2023)
    // 113621, Sec. 2.
    const std::vector<double>& pz = pg->zProfileMeters;
    const std::vector<double>& pr = pg->rProfileMeters;
    size_t throat = 0;
    for (size_t i = 1; i < pr.size(); i++)
        if (pr[i] < pr[throat]) throat = i;
    if (throat + 1 >= pz.size()) {
        ibsimu.message(1) << "WARNING: beam halo requested but the plasma-grid profile has "
                          << "no downstream chamfer after its throat; no halo emitted." << endl;
        return 0U;
    }

    // Lateral surface of each revolved segment, 2*pi*r_mid*ds, as the sampling weight so
    // that particles are spread uniformly over the real emitting area rather than uniformly
    // along the profile, which would over-weight the narrow end.
    std::vector<double> cum;
    double area = 0.0;
    for (size_t i = throat; i + 1 < pz.size(); i++) {
        const double dz = pz[i + 1] - pz[i];
        const double dr = pr[i + 1] - pr[i];
        area += 2.0 * M_PI * 0.5 * (pr[i] + pr[i + 1]) * sqrt(dz * dz + dr * dr);
        cum.push_back(area);
    }
    if (!(area > 0.0)) {
        ibsimu.message(1) << "WARNING: beam halo chamfer has zero area; no halo emitted." << endl;
        return 0U;
    }

    // Macroparticle count scales with the requested current fraction, so each halo particle
    // carries the same charge as a core one. Statistics stay uniform across the beam and the
    // added cost is exactly the halo fraction -- 8 % more particles for an 8 % halo.
    const double fraction = 0.01 * params.getBeamHaloCurrentFractionPercent();

    const std::vector<SimulationParameters::ParticleSourceDefinition>& srcs =
        params.getParticleSources();
    if (srcs.empty()) return 0U;

    // Start a quarter cell off the surface: a particle born exactly on the boundary is
    // reported as an immediate collision and never enters the beam.
    const double standoff = 0.25 * params.getMeshSize();

    std::mt19937_64 rng(20261005ULL);
    std::uniform_real_distribution<double> U(0.0, 1.0);

    // One halo per source, around that source's own axis. The plasma-grid chamfer is the
    // same revolved profile at every aperture of a patterned grid, so with multi-beamlet
    // extraction each beamlet gets its own ring carrying its own share of the current. Doing
    // this from srcs.front() alone -- as this did before multi-beamlet extraction existed --
    // put the halo of every beamlet onto the first one, which is harmless with a single
    // source and badly wrong with several.
    size_t n_halo_total = 0U;
    double halo_current_A = 0.0;
    for (std::vector<SimulationParameters::ParticleSourceDefinition>::const_iterator it = srcs.begin();
         it != srcs.end();
         ++it) {
        const SimulationParameters::ParticleSourceDefinition& ref = *it;

        // Macroparticle count scales with the requested current fraction, so each halo
        // particle carries the same charge as a core one. Statistics stay uniform across the
        // beam and the added cost is exactly the halo fraction -- 8 % more particles for an
        // 8 % halo.
        const size_t n_halo =
            static_cast<size_t>(llround(fraction * static_cast<double>(ref.particleCount)));
        if (n_halo == 0U) continue;

        const double src_current_A = ref.currentDensityAm2 * ref.widthMeters * ref.heightMeters;
        const double cx = ref.centerMeters[0];
        const double cy = ref.centerMeters[1];

        const double m = ref.massU * MASS_U;
        const double v = sqrt(2.0 * params.getBeamHaloAxialEnergyEV() * CHARGE_E / m);
        const double IQ = fraction * src_current_A / static_cast<double>(n_halo);

        // The halo of a beamlet belongs to that beamlet: tag it the same way the core is
        // tagged, so a halo trajectory is not mistaken for an unattributed secondary.
        const size_t halo_first = pdb->size();
        const size_t halo_slot = beamletSlotFor(
            ref.beamletIndex, ref.beamletColumn, ref.beamletRow, cx, cy);

        for (size_t a = 0; a < n_halo; a++) {
            // Pick a segment by area, then a point within it.
            const double target = U(rng) * area;
            size_t seg = 0;
            while (seg + 1 < cum.size() && cum[seg] < target) seg++;
            const size_t i = throat + seg;
            // Place the point within the segment by AREA, not uniformly along the profile.
            // Revolving gives dA = 2*pi*r ds with r linear in t, so the density is p(t) ~ r(t)
            // and sampling t uniformly would over-weight the narrow end of the chamfer. The
            // inverse of that CDF is the positive root of (dr/2) t^2 + r1 t - U (r1 + dr/2) = 0.
            const double r1 = pr[i];
            const double dr_seg = pr[i + 1] - pr[i];
            const double u = U(rng);
            double t;
            if (fabs(dr_seg) < 1e-15 * (fabs(r1) + 1e-30)) {
                t = u;                                // cylindrical land: area is linear in t
            } else {
                const double disc = r1 * r1 + 2.0 * dr_seg * u * (r1 + 0.5 * dr_seg);
                t = (-r1 + sqrt(disc > 0.0 ? disc : 0.0)) / dr_seg;
                if (t < 0.0) t = 0.0;
                if (t > 1.0) t = 1.0;
            }
            const double z = pz[i] + t * (pz[i + 1] - pz[i]);
            const double r = pr[i] + t * dr_seg;

            const double dz = pz[i + 1] - pz[i];
            const double dr = pr[i + 1] - pr[i];
            const double len = sqrt(dz * dz + dr * dr);
            // Surface normal pointing out of the metal into the bore: radially inward (-dz) and,
            // for a widening chamfer, downstream (+dr).
            const double n_r = -dz / len;
            const double n_z = dr / len;

            const double phi = 2.0 * M_PI * U(rng);
            const double cphi = cos(phi), sphi = sin(phi);

            ParticleP3D x;
            x[0] = 0.0;
            x[1] = cx + (r + standoff * n_r) * cphi;  x[2] = v * n_r * cphi;
            x[3] = cy + (r + standoff * n_r) * sphi;  x[4] = v * n_r * sphi;
            x[5] = z + standoff * n_z;                x[6] = v * n_z;
            pdb->add_particle(IQ, ref.chargeState, ref.massU, 0, x);
        }

        tagParticleRange(halo_first, pdb->size(), halo_slot);

        n_halo_total += n_halo;
        halo_current_A += fraction * src_current_A;
    }

    if (n_halo_total == 0U) {
        ibsimu.message(1) << "WARNING: beam halo fraction yields zero macroparticles; "
                          << "no halo emitted." << endl;
        return 0U;
    }

    // add_particle deliberately leaves _rhosum untouched, unlike add_rectangular_beam_*.
    // That is what keeps the plasma reference density (ManageSimulation_New.cpp, which reads
    // get_rhosum() to configure set_nsimp_plasma) equal to the core beam alone, so enabling
    // the halo cannot move the meniscus.
    ibsimu.message(1) << "Beam halo: " << n_halo_total << " particles over " << srcs.size()
                      << " beamlet(s), " << params.getBeamHaloCurrentFractionPercent() << " % of "
                      << core_current_A * 1e3 << " mA core = " << halo_current_A * 1e3
                      << " mA, from the plasma-grid chamfer of solid '" << pg->name
                      << "' (z " << pz[throat] << " to " << pz.back() << " m), at "
                      << params.getBeamHaloAxialEnergyEV() << " eV" << endl;
    return n_halo_total;
}

void ParticleManager::defineParticleDatabase(const Geometry& geometry, const SimulationParameters& params) {
    ParticleDataBase3D* pdb = new ParticleDataBase3D(geometry);
    pdb->set_max_steps(3000);
    bool pmirror[6] = {false, false, false, false, false, false};
    pdb->set_mirror(pmirror);
    pdb->set_surface_collision(true);

    pdb->clear();

    // The database is rebuilt from scratch here (and once per iteration when pdbincycle is
    // set), so the beamlet map has to be rebuilt with it or it would describe the previous
    // database.
    beamlet_tags.clear();
    particle_beamlet_slot.clear();

    if (params.hasParticleSources()) {
        const std::vector<SimulationParameters::ParticleSourceDefinition>& sources =
            params.getParticleSources();
        uint total_particles = 0U;
        double total_current_density = 0.0;
        double total_current = 0.0;              // A, summed over the core sources

        for (std::vector<SimulationParameters::ParticleSourceDefinition>::const_iterator it = sources.begin();
             it != sources.end();
             ++it) {
            const SimulationParameters::ParticleSourceDefinition& source = *it;

            if (source.sourceModel != "uniform") {
                throw Error(
                    ERROR_LOCATION,
                    "Only uniform particleSources are supported by the current C++ runtime");
            }

            Vec3D center(source.centerMeters[0], source.centerMeters[1], source.centerMeters[2]);
            Vec3D main_direction(source.mainDirection[0], source.mainDirection[1], source.mainDirection[2]);
            Vec3D reference_direction(
                source.inPlaneReferenceDirection[0],
                source.inPlaneReferenceDirection[1],
                source.inPlaneReferenceDirection[2]);

            Vec3D dir2 = cross(main_direction, reference_direction);
            Vec3D dir1 = cross(dir2, main_direction);
            if (main_direction.ssqr() == 0.0 || dir1.ssqr() == 0.0 || dir2.ssqr() == 0.0) {
                throw Error(
                    ERROR_LOCATION,
                    "Particle source '" + source.id + "' has invalid mainDirection or inPlaneReferenceDirection");
            }

            // Bracket the emitter: every particle the database gains across this call came
            // out of this source, hence out of this beamlet. Taking the indices rather than
            // trusting particleCount keeps the map right even if IBSimu ever emits a
            // different number than it was asked for.
            const size_t first_particle = pdb->size();

            pdb->add_rectangular_beam_with_energy(
                source.particleCount,
                source.currentDensityAm2,
                source.chargeState,
                source.massU,
                source.axialEnergyEV,
                source.parallelTemperatureEV,
                source.perpendicularTemperatureEV,
                center,
                dir1,
                dir2,
                0.5 * source.widthMeters,
                0.5 * source.heightMeters);

            tagParticleRange(
                first_particle,
                pdb->size(),
                beamletSlotFor(
                    source.beamletIndex,
                    source.beamletColumn,
                    source.beamletRow,
                    source.centerMeters[0],
                    source.centerMeters[1]));

            total_particles += source.particleCount;
            total_current_density += source.currentDensityAm2;
            total_current += source.currentDensityAm2 * source.widthMeters * source.heightMeters;
        }

        const size_t halo_particles = addBeamHalo(pdb, params, total_current);

        particles = pdb;
        ibsimu.message(1) << "Particle database created with " << total_particles
                          << " particles from " << sources.size() << " configured source(s), total J="
                          << total_current_density << " A/m²" << endl;
        if (halo_particles > 0U) {
            ibsimu.message(1) << "  plus " << halo_particles << " beam-halo particles" << endl;
        }
        return;
    }

    // Legacy single-source fallback - beam center inferred from current geometry bounds.
    Vec3D ori = geometry.origo();
    Vec3D fin = geometry.max();
    double x_start = (ori[0] > -0.01) ? ori[0] : -0.01;
    double y_start = (ori[1] > -0.01) ? ori[1] : -0.01;
    double z_start = ori[2];
    double x_end = (fin[0] < 0.01) ? fin[0] : 0.01;
    double y_end = (fin[1] < 0.01) ? fin[1] : 0.01;

    pdb->add_rectangular_beam_with_energy(
        params.getNParticles(),
        params.getJIon(),
        params.getQIons(),
        params.getMIons(),
        params.getE0Z(),
        params.getTPar(),
        params.getTPerp(),
        Vec3D(0.5 * (x_start + x_end), 0.5 * (y_start + y_end), z_start),
        Vec3D(1, 0, 0),
        Vec3D(0, 1, 0),
        0.5 * (x_end - x_start),
        0.5 * (y_end - y_start));

    particles = pdb;

    ibsimu.message(1) << "Particle database created with " << params.getNParticles()
                      << " particles, J=" << params.getJIon() << " A/m²" << endl;
}

void ParticleManager::defineParticleDatabase(const Geometry& geometry, const SimulationParameters& params,
                                            const string& emittername, bool backstream) {
    ParticleDataBase3D* pdb = new ParticleDataBase3D(geometry);
    pdb->set_max_steps(3000);
    bool pmirror[6] = {false, false, false, false, false, false};
    pdb->set_mirror(pmirror);
    pdb->set_surface_collision(true);
    
    // Beam definition moved outside the main iteration cycle
    // There is no need to change the starting particles at every iteration.
    pdb->clear();

    // An emitter file carries no beamlet pattern, so leave the map empty: the VTK export
    // then omits the beamlet arrays rather than labelling every trajectory "beamlet 1".
    beamlet_tags.clear();
    particle_beamlet_slot.clear();

    // Load emitter from file
    string emitter_file = emittername;
    loadEmitter(pdb, emitter_file, 1.0, 1.0, 0.0); // Forward particles with original charge and mass
    
    if (backstream) {
        // Load backstreaming particles to compensate space charge
        loadEmitter(pdb, emitter_file, -1.0, 1.0, 1.0); // Reversed velocities, different mass factor for visualization
    }

    particles = pdb;
    
    ibsimu.message(1) << "Particle database created with " << pdb->size() << " particles from emitter: " << emittername << endl;
}

void ParticleManager::fillParticleDatabases(const SimulationParameters& params) {
    // create pdb for each species
    if (params.getIncludeStripping() == 2) {
        const particle_kind species_kinds[] = {
            PARTICLE_NEGATIVE_ION,
            PARTICLE_NEUTRAL_ATOM,
            PARTICLE_POSITIVE_ION,
            PARTICLE_MOLECULAR_POSITIVE_ION,
            PARTICLE_MOLECULAR_NEUTRAL,
            PARTICLE_TRIATOMIC_POSITIVE_ION,
            PARTICLE_ELECTRON
        };

        ParticleDataBase3D* pdbHm = new ParticleDataBase3D(*particles);
        ParticleDataBase3D* pdbH0 = new ParticleDataBase3D(*particles);
        ParticleDataBase3D* pdbHp = new ParticleDataBase3D(*particles);
        ParticleDataBase3D* pdbH2p = new ParticleDataBase3D(*particles);
        ParticleDataBase3D* pdbH20 = new ParticleDataBase3D(*particles);
        ParticleDataBase3D* pdbH3p = new ParticleDataBase3D(*particles);
        ParticleDataBase3D* pdbe = new ParticleDataBase3D(*particles);

        *pdbHm = *particles;
        *pdbH0 = *particles;
        *pdbHp = *particles;
        *pdbH2p = *particles;
        *pdbH20 = *particles;
        *pdbH3p = *particles;
        *pdbe = *particles;

        if (debug) logfile << "\nDEBUG: SPECIES DATABASES CREATED. NOW CLEANING... \n" << flush;
        
        // Initialize counters for each species
        size_t nHm = 0, nH0 = 0, nHp = 0, nH2p = 0, nH20 = 0, nH3p = 0, ne = 0;
        
        // Filter particles by species - iterate through all particles and identify species
        for (size_t id = 0; id < particles->size(); id++) {
            Particle3D& pp = particles->particle(id);
            
            // Identify particle species based on mass and charge
            particle_kind species = identify_particle_species(pp.m(), pp.q(), params.getMIons());
            
            // Reset trajectory for particles that don't belong to each species database
            if (species != PARTICLE_NEGATIVE_ION) pdbHm->reset_trajectory(id);
            else nHm++;
            
            if (species != PARTICLE_NEUTRAL_ATOM) pdbH0->reset_trajectory(id);
            else nH0++;
            
            if (species != PARTICLE_POSITIVE_ION) pdbHp->reset_trajectory(id);
            else nHp++;
            
            if (species != PARTICLE_MOLECULAR_POSITIVE_ION) pdbH2p->reset_trajectory(id);
            else nH2p++;
            
            if (species != PARTICLE_MOLECULAR_NEUTRAL) pdbH20->reset_trajectory(id);
            else nH20++;

            if (species != PARTICLE_TRIATOMIC_POSITIVE_ION) pdbH3p->reset_trajectory(id);
            else nH3p++;
            
            if (species != PARTICLE_ELECTRON) pdbe->reset_trajectory(id);
            else ne++;
        }

        ParticleDataBase3D* species_databases[] = {pdbHm, pdbH0, pdbHp, pdbH2p, pdbH20, pdbH3p, pdbe};
        size_t species_counts[] = {nHm, nH0, nHp, nH2p, nH20, nH3p, ne};
        
        // Store particle count for each species
        for (size_t index = 0; index < particle_kind_count(); ++index) {
            part_number.push_back(species_counts[index]);
            particles_species.push_back(species_databases[index]);
            particles_kinds.push_back(species_kinds[index]);
        }
        
        logfile << "\nSPECIES DATABASES CREATED! \n" << flush;
        size_t species_total = 0;
        for (size_t index = 0; index < particle_kind_count(); ++index) {
            logfile << " #" << get_particle_name(species_kinds[index]) << ": "
                    << species_counts[index] << endl << flush;
            species_total += species_counts[index];
        }
        logfile << " SUM: " << species_total << " \n\n" << flush;
    }
    else {
        ParticleDataBase3D* pdbHm = new ParticleDataBase3D(*particles);
        *pdbHm = *particles;
        particles_species.push_back(pdbHm); 
        particles_kinds.push_back(PARTICLE_NEGATIVE_ION);
        size_t nHm = pdbHm->size();
        part_number.push_back(nHm);
    }
}

void ParticleManager::analyzeParticleEndLocations(const SimulationParameters& params) {
    if (!particles) {
        ibsimu.message(1) << "No particle database available for end location analysis" << endl;
        return;
    }
    
    ibsimu.message(1) << "Analyzing particle end locations..." << endl;
    
    try {
        // Initialize counters for different particle fates
        size_t particles_extracted = 0;
        size_t particles_out = 0;
        size_t particles_stripped = 0;
        size_t particles_collided = 0;
        size_t particles_ok = 0;
        
        // Initialize generation counters
        map<int, size_t> generation_counts;
        
        // Initialize species counters
        map<particle_kind, size_t> species_counts;
        
        // Initialize collision location analysis
        vector<double> collision_z_positions;
        vector<double> extraction_z_positions;
        
        // Energy analysis
        vector<double> final_energies;
        vector<double> initial_energies;
        
        // Current analysis
        double total_initial_current = 0.0;
        double total_final_current = 0.0;
        double extracted_current = 0.0;
        double extraction_plane_z = -1.0;
        if (params.hasDiagnosticTransmissionPlaneZPosition()) {
            extraction_plane_z = params.getDiagnosticTransmissionPlaneZPosition();
        } else if (params.hasGeneratedGeometrySolids()) {
            const std::vector<SimulationParameters::GeometrySolidDefinition>& solids =
                params.getGeneratedGeometrySolids();
            for (std::vector<SimulationParameters::GeometrySolidDefinition>::const_iterator solid_it = solids.begin();
                 solid_it != solids.end();
                 ++solid_it) {
                if (solid_it->zProfileMeters.empty()) {
                    continue;
                }
                const double solid_z_min = *std::min_element(solid_it->zProfileMeters.begin(),
                                                             solid_it->zProfileMeters.end());
                if (extraction_plane_z < 0.0 || solid_z_min < extraction_plane_z) {
                    extraction_plane_z = solid_z_min;
                }
            }
        }
        
        // Analyze each particle
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D& particle = particles->particle(i);
            
            // Get particle properties
            double mass = particle.m();
            double charge = particle.q();
            double current = particle.IQ();
            int generation = particle.gen() % 100;  // Extract generation number (0-4+)
            particle_status_e status = particle.get_status();
            
            // Count by generation
            generation_counts[generation]++;
            
            // Identify particle species
            particle_kind species = identify_particle_species(mass, charge, params.getMIons());
            species_counts[species]++;
            
            // Add to current totals
            total_initial_current += abs(current);
            
            // Get trajectory information
            if (particle.traj_size() > 0) {
                // Get initial and final trajectory points
                ParticleP3D initial_point = particles->trajectory_point(i, 0);
                ParticleP3D final_point = particles->trajectory_point(i, particle.traj_size() - 1);
                
                // Calculate initial and final energies
                double initial_vx = initial_point(2);
                double initial_vy = initial_point(4);
                double initial_vz = initial_point(6);
                double initial_v_mag = sqrt(initial_vx*initial_vx + initial_vy*initial_vy + initial_vz*initial_vz);
                // Relativistic kinetic energy using constants c and c2 from constants.hpp
                // KE = (gamma - 1) * m * c^2, gamma = 1 / sqrt(1 - (v^2/c^2))
                double v2 = initial_v_mag * initial_v_mag;
                double beta2 = v2 / SPEED_C2;
                if (beta2 >= 1.0) beta2 = 0.999999999999; // clamp to avoid nan
                double gamma = 1.0 / sqrt(1.0 - beta2);
                double initial_energy_eV = (gamma - 1.0) * mass * SPEED_C2 / CHARGE_E;
                initial_energies.push_back(initial_energy_eV);
                
                double final_vx = final_point(2);
                double final_vy = final_point(4);
                double final_vz = final_point(6);
                double final_v_mag = sqrt(final_vx*final_vx + final_vy*final_vy + final_vz*final_vz);
                double final_energy_eV = 0.5 * mass * final_v_mag * final_v_mag / CHARGE_E;
                final_energies.push_back(final_energy_eV);
                
                // Get final position
                double final_z = final_point(5);

                // Analyze particle fate based on status
                switch (status) {
                    case PARTICLE_OK:
                        particles_ok++;
                        total_final_current += abs(current);
                        // Check if particle reached extraction plane
                        if (extraction_plane_z >= 0.0 && final_z >= extraction_plane_z) {
                            particles_extracted++;
                            extracted_current += abs(current);
                            extraction_z_positions.push_back(final_z);
                        }
                        break;
                        
                    case PARTICLE_COLL:
                        particles_collided++;
                        collision_z_positions.push_back(final_z);
                        break;
                        
                    case PARTICLE_STRIP:
                        particles_stripped++;
                        collision_z_positions.push_back(final_z);
                        break;
                        
                    case PARTICLE_OUT:
                        particles_out++;
                        break;
                        
                    default:
                        break;
                }
            }
        }
        
        // Calculate statistics
        size_t total_particles = particles->size();
        double extraction_efficiency = (total_particles > 0) ? 
            (double)particles_extracted / total_particles * 100.0 : 0.0;
        double stripping_rate = (total_particles > 0) ? 
            (double)particles_stripped / total_particles * 100.0 : 0.0;
        double collision_rate = (total_particles > 0) ? 
            (double)particles_collided / total_particles * 100.0 : 0.0;
        double out_rate = (total_particles > 0) ? 
            (double)particles_out / total_particles * 100.0 : 0.0;
        
        // Calculate current transmission
        double current_transmission = (total_initial_current > 0) ? 
            extracted_current / total_initial_current * 100.0 : 0.0;
        
        // Log detailed analysis results
        logfile << "\n=== PARTICLE END LOCATION ANALYSIS ===" << endl;
        logfile << "Total particles analyzed: " << total_particles << endl;
        logfile << "\nParticle Fates:" << endl;
        logfile << "  Extracted (OK at z>=9mm): " << particles_extracted 
                << " (" << extraction_efficiency << "%)" << endl;
        logfile << "  Still OK (not extracted): " << (particles_ok - particles_extracted) << endl;
        logfile << "  Collided with surfaces: " << particles_collided 
                << " (" << collision_rate << "%)" << endl;
        logfile << "  Stripped: " << particles_stripped 
                << " (" << stripping_rate << "%)" << endl;
        logfile << "  Out: " << particles_out 
                << " (" << out_rate << "%)" << endl;

        logfile << "\nCurrent Analysis:" << endl;
        logfile << "  Initial total current: " << total_initial_current << " A" << endl;
        logfile << "  Final total current: " << total_final_current << " A" << endl;
        logfile << "  Extracted current: " << extracted_current << " A" << endl;
        logfile << "  Current transmission: " << current_transmission << "%" << endl;
        
        // Generation analysis
        logfile << "\nGeneration Analysis:" << endl;
        for (const auto& gen_pair : generation_counts) {
            logfile << "  Generation " << gen_pair.first << ": " << gen_pair.second << " particles" << endl;
        }
        
        // Species analysis
        logfile << "\nSpecies Analysis:" << endl;
        for (const auto& species_pair : species_counts) {
            string species_name = get_particle_name(species_pair.first);
            logfile << "  " << species_name << ": " << species_pair.second << " particles" << endl;
        }
        
        // Collision location statistics
        if (!collision_z_positions.empty()) {
            sort(collision_z_positions.begin(), collision_z_positions.end());
            double min_collision_z = collision_z_positions.front();
            double max_collision_z = collision_z_positions.back();
            double median_collision_z = collision_z_positions[collision_z_positions.size() / 2];
            
            logfile << "\nCollision Location Statistics:" << endl;
            logfile << "  First collision at z = " << min_collision_z * 1000 << " mm" << endl;
            logfile << "  Last collision at z = " << max_collision_z * 1000 << " mm" << endl;
            logfile << "  Median collision at z = " << median_collision_z * 1000 << " mm" << endl;
        }
        
        // Extraction location statistics
        if (!extraction_z_positions.empty()) {
            sort(extraction_z_positions.begin(), extraction_z_positions.end());
            double min_extraction_z = extraction_z_positions.front();
            double max_extraction_z = extraction_z_positions.back();
            
            logfile << "\nExtraction Statistics:" << endl;
            logfile << "  First extraction at z = " << min_extraction_z * 1000 << " mm" << endl;
            logfile << "  Last extraction at z = " << max_extraction_z * 1000 << " mm" << endl;
        }
        
        // Energy statistics
        if (!initial_energies.empty() && !final_energies.empty()) {
            double avg_initial_energy = accumulate(initial_energies.begin(), initial_energies.end(), 0.0) / initial_energies.size();
            double avg_final_energy = accumulate(final_energies.begin(), final_energies.end(), 0.0) / final_energies.size();
            
            logfile << "\nEnergy Statistics:" << endl;
            logfile << "  Average initial energy: " << avg_initial_energy << " eV" << endl;
            logfile << "  Average final energy: " << avg_final_energy << " eV" << endl;
            logfile << "  Average energy gain: " << (avg_final_energy - avg_initial_energy) << " eV" << endl;
        }
        
        // Performance metrics for optimization
        logfile << "\nPerformance Metrics:" << endl;
        logfile << "  Extraction efficiency: " << extraction_efficiency << "%" << endl;
        logfile << "  Current transmission: " << current_transmission << "%" << endl;
        logfile << "  Beam losses: " << (collision_rate + out_rate) << "%" << endl;
        
        // IBSimu console output (summary)
        ibsimu.message(1) << "Analysis complete: " << particles_extracted << "/" << total_particles 
                          << " particles extracted (" << extraction_efficiency << "%)" << endl;
        ibsimu.message(1) << "Current transmission: " << current_transmission << "%" << endl;
        
        if (params.getIncludeStripping() > 0) {
            ibsimu.message(1) << "Stripping events: " << particles_stripped 
                              << " (" << stripping_rate << "%)" << endl;
        }
        
        logfile << "=== END PARTICLE ANALYSIS ===\n" << endl;
        
    } catch (const Error& e) {
        logfile << "ERROR in particle end location analysis: IBSimu Error occurred" << endl;
        ibsimu.message(1) << "Error during particle analysis" << endl;
    } catch (const std::exception& e) {
        logfile << "ERROR in particle end location analysis: " << e.what() << endl;
        ibsimu.message(1) << "Error during particle analysis: " << e.what() << endl;
    }
}

string ParticleManager::saveEmitter(const string& emitname, double zloc, const string& output_folder) {
    string filename = output_folder + emitname;
    
    if (!particles) {
        throw Error(ERROR_LOCATION, "No particle database available to save emitter");
    }
    
    ibsimu.message(1) << "Saving particles at z = " << zloc << " m to emitter file: " << filename << endl;
    
    try {
        // Create trajectory diagnostic data at the specified z location
        TrajectoryDiagnosticData tdata;
        
        // Get trajectories at the specified plane
        particles->trajectories_at_plane(tdata, AXIS_Z, zloc, diagnostics);
        
        // Write emitter file
        ofstream outfile(filename);
        if (!outfile.is_open()) {
            throw Error(ERROR_LOCATION, "Could not open file for writing: " + filename);
        }
        
        // Write header with detailed information
        outfile << "# Emitter file saved at z = " << zloc << " m" << endl;
        outfile << "# Format: x[m] y[m] z[m] vx[m/s] vy[m/s] vz[m/s] I[A] m[kg] q[C]" << endl;
        outfile << "# Number of particles: " << tdata.traj_size() << endl;
        
        if (tdata.traj_size() > 0) {
            for (uint32_t i = 0; i < tdata.traj_size(); i++) {
                // Write particle data: x, y, z, vx, vy, vz, current, mass, charge
                outfile << fixed << setprecision(6)
                        << tdata(i, 0) << " "  // x position
                        << tdata(i, 1) << " "  // y position  
                        << tdata(i, 2) << " "  // z position
                        << tdata(i, 3) << " "  // vx velocity
                        << tdata(i, 4) << " "  // vy velocity
                        << tdata(i, 5) << " "  // vz velocity
                        << tdata(i, 6) << " "  // current
                        << tdata(i, 7) << " "  // mass
                        << tdata(i, 8) << endl; // charge
            }
            
            ibsimu.message(1) << "Saved " << tdata.traj_size() << " particles to emitter file" << endl;
        } else {
            ibsimu.message(1) << "Warning: No particles found at z = " << zloc << " m" << endl;
        }
        
        outfile.close();
        
    } catch (const Error& e) {
        logfile << "ERROR in saveEmitter: IBSimu Error occurred" << endl;
        ibsimu.message(1) << "Error saving emitter file" << endl;
        throw;
    } catch (const std::exception& e) {
        logfile << "ERROR in saveEmitter: " << e.what() << endl;
        ibsimu.message(1) << "Error saving emitter file: " << e.what() << endl;
        throw Error(ERROR_LOCATION, "Error saving emitter: " + string(e.what()));
    }
    
    ibsimu.message(1) << "Emitter saved to: " << filename << endl;
    return filename;
}

bool ParticleManager::checkEGExtractedCurrent(double oriJ, double& extsimJ, const SimulationParameters& params) {
    ParticleDataBase3D* extracted_particles = particles;
    int hm_index = get_particle_int(PARTICLE_NEGATIVE_ION);
    if (hm_index >= 0 && hm_index < static_cast<int>(particles_species.size()) && particles_species[hm_index]) {
        extracted_particles = particles_species[hm_index];
    }

    if (!extracted_particles) {
        extsimJ = 0.0;
        return false;
    }
    
    try {
        double eg_exit_z = -1.0;
        if (params.hasDiagnosticTransmissionPlaneZPosition()) {
            eg_exit_z = params.getDiagnosticTransmissionPlaneZPosition();
        } else if (params.hasGeneratedGeometrySolids()) {
            const std::vector<SimulationParameters::GeometrySolidDefinition>& solids =
                params.getGeneratedGeometrySolids();
            for (std::vector<SimulationParameters::GeometrySolidDefinition>::const_iterator solid_it = solids.begin();
                 solid_it != solids.end();
                 ++solid_it) {
                if (solid_it->zProfileMeters.empty()) {
                    continue;
                }
                const double solid_z_min = *std::min_element(solid_it->zProfileMeters.begin(),
                                                             solid_it->zProfileMeters.end());
                if (eg_exit_z < 0.0 || solid_z_min < eg_exit_z) {
                    eg_exit_z = solid_z_min;
                }
            }
        }
        if (eg_exit_z < 0.0) {
            extsimJ = 0.0;
            if (debug) logfile << "DEBUG: EG extracted current calculation skipped because no transmission plane is configured" << endl;
            return false;
        }
        
        // Create trajectory diagnostic data at EG exit
        TrajectoryDiagnosticData tdata;
        vector<trajectory_diagnostic_e> diagnostics = {
            DIAG_X, DIAG_Y, DIAG_Z, DIAG_VX, DIAG_VY, DIAG_VZ,
            DIAG_CURR, DIAG_MASS, DIAG_CHARGE, DIAG_NO
        };
        
        // Get trajectories at EG exit plane
        extracted_particles->trajectories_at_plane(tdata, AXIS_Z, eg_exit_z, diagnostics);
        
        // Calculate total current at EG exit
        double total_current = 0.0;
        size_t particle_count = tdata.traj_size();
        
        for (size_t i = 0; i < particle_count; ++i) {
            total_current += tdata(i, 6); // DIAG_CURR is column 6
        }
        
        double aperture_radius = params.getDiagnosticApertureRadiusMeters();
        double pg_aperture_area = M_PI * aperture_radius * aperture_radius;
        double current_density = abs(total_current) / pg_aperture_area;
        
        extsimJ = current_density;
        
        // Check if within tolerance
        double tolerance = params.getJTolerance();
        double error = abs(extsimJ - oriJ) / oriJ;
        bool within_tolerance = (error < tolerance);
        
        if (debug) {
            logfile << "DEBUG: EG extracted current calculation:" << endl;
            logfile << "  Current source: "
                    << ((extracted_particles == particles) ? "all-particle database (fallback)" : "family negative-ion database")
                    << endl;
            logfile << "  Particles at EG exit (z=" << eg_exit_z << "m): " << particle_count << endl;
            logfile << "  Total current: " << total_current << " A" << endl;
            logfile << "  PG aperture area: " << pg_aperture_area << " m²" << endl;
            logfile << "  Extracted current density: " << extsimJ << " A/m²" << endl;
            logfile << "  Target current density: " << oriJ << " A/m²" << endl;
            logfile << "  Error: " << error*100 << "%" << endl;
            logfile << "  Tolerance: " << tolerance*100 << "%" << endl;
            logfile << "  Within tolerance: " << (within_tolerance ? "YES" : "NO") << endl;
        }
        
        ibsimu.message(1) << "Current check: target=" << oriJ << " A/m², extracted=" << extsimJ 
                          << " A/m², error=" << error*100 << "%, tolerance=" << tolerance*100 << "%" << endl;
        
        return within_tolerance;
        
    } catch (const Error& e) {
        if (debug) logfile << "DEBUG: Error in EG extracted current calculation: IBSimu Error occurred" << endl;
        extsimJ = 0.0;
        return false;
    } catch (const std::exception& e) {
        if (debug) logfile << "DEBUG: Error in EG extracted current calculation: " << e.what() << endl;
        extsimJ = 0.0;
        return false;
    }
}

void ParticleManager::exportTrajectoriesToVTK(const string& filename, double ion_mass_u) {
    if (!particles) {
        ibsimu.message(1) << "No particle database available for trajectories VTK export" << endl;
        return;
    }
    
    string vtkfile = filename + "_trajectories.vtk";
    
    try {
        ofstream file(vtkfile, ios::binary);
        if (!file.is_open()) {
            throw Error("Cannot open VTK file for writing: " + vtkfile);
        }
        
        // Count total trajectory points
        int total_points = 0;
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            total_points += particle.traj_size();
        }
        
        // Debug output
        ibsimu.message(1) << "Exporting " << particles->size() << " particle trajectories in binary format" << endl;
        ibsimu.message(1) << "Total trajectory points: " << total_points << endl;
        if (particles->size() > 0) {
            Particle3D &first_particle = particles->particle(0);
            if (first_particle.traj_size() > 0) {
                ibsimu.message(1) << "First particle: " << first_particle.traj_size() << " points" << endl;
                ParticleP3D start_point = particles->trajectory_point(0, 0);
                ParticleP3D end_point = particles->trajectory_point(0, first_particle.traj_size()-1);
                ibsimu.message(1) << "  Start: (" << start_point(1) << ", " << start_point(3) << ", " << start_point(5) << ")" << endl;
                ibsimu.message(1) << "  End:   (" << end_point(1) << ", " << end_point(3) << ", " << end_point(5) << ")" << endl;
            }
        }

        if (total_points == 0) {
            ibsimu.message(1) << "No trajectory points to export" << endl;
            return;
        }
        
        // Helper lambda to swap endianness for binary VTK (big-endian required)
        auto swap_endian_float = [](float value) -> float {
            union { float f; uint32_t i; } u;
            u.f = value;
            u.i = ((u.i >> 24) & 0x000000FF) | 
                  ((u.i >> 8)  & 0x0000FF00) | 
                  ((u.i << 8)  & 0x00FF0000) | 
                  ((u.i << 24) & 0xFF000000);
            return u.f;
        };
        
        auto swap_endian_int = [](int value) -> int {
            return ((value >> 24) & 0x000000FF) | 
                   ((value >> 8)  & 0x0000FF00) | 
                   ((value << 8)  & 0x00FF0000) | 
                   ((value << 24) & 0xFF000000);
        };
        
        // Write VTK header (ASCII)
        file << "# vtk DataFile Version 3.0\n";
        file << "IBSIMU Particle Trajectories Export\n";
        file << "BINARY\n";
        file << "DATASET POLYDATA\n";
        
        // Write points header (ASCII) then data (binary)
        file << "POINTS " << total_points << " float\n";
        
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            for (size_t j = 0; j < particle.traj_size(); j++) {
                ParticleP3D traj_point = particles->trajectory_point(i, j);
                // Write x, y, z coordinates as big-endian floats
                float x = static_cast<float>(traj_point(1));
                float y = static_cast<float>(traj_point(3));
                float z = static_cast<float>(traj_point(5));
                float x_be = swap_endian_float(x);
                float y_be = swap_endian_float(y);
                float z_be = swap_endian_float(z);
                file.write(reinterpret_cast<const char*>(&x_be), sizeof(float));
                file.write(reinterpret_cast<const char*>(&y_be), sizeof(float));
                file.write(reinterpret_cast<const char*>(&z_be), sizeof(float));
            }
        }
        
        // Write lines (trajectory polylines)
        int total_lines = particles->size();
        int total_line_data = 0;
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            total_line_data += particle.traj_size() + 1; // +1 for line length
        }
        
        file << "\nLINES " << total_lines << " " << total_line_data << "\n";
        
        int point_offset = 0;
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            int traj_size = particle.traj_size();
            // Write number of points in this line
            int traj_size_be = swap_endian_int(traj_size);
            file.write(reinterpret_cast<const char*>(&traj_size_be), sizeof(int));
            // Write point indices
            for (int j = 0; j < traj_size; j++) {
                int idx_be = swap_endian_int(point_offset + j);
                file.write(reinterpret_cast<const char*>(&idx_be), sizeof(int));
            }
            point_offset += traj_size;
        }
        
        // Write point data - particle velocities and energies
        file << "\nPOINT_DATA " << total_points << "\n";
        
        // Export velocity magnitude in km/s for better visualization
        file << "SCALARS velocity_magnitude_km_s float\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            for (size_t j = 0; j < particle.traj_size(); j++) {
                ParticleP3D traj_point = particles->trajectory_point(i, j);
                // Velocity components: vx=index 2, vy=index 4, vz=index 6 in m/s
                double vx = traj_point(2);
                double vy = traj_point(4);
                double vz = traj_point(6);
                double vel_mag = sqrt(vx*vx + vy*vy + vz*vz);
                // Convert to km/s for better scale
                float vel_km_s = static_cast<float>(vel_mag / 1000.0);
                float vel_km_s_be = swap_endian_float(vel_km_s);
                file.write(reinterpret_cast<const char*>(&vel_km_s_be), sizeof(float));
            }
        }
        
        // Export time in nanoseconds for better scale
        file << "\nSCALARS time_ns float\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            for (size_t j = 0; j < particle.traj_size(); j++) {
                ParticleP3D traj_point = particles->trajectory_point(i, j);
                // Time is at index 0 in seconds
                double time = traj_point(0);
                // Convert to nanoseconds for better scale
                float time_ns = static_cast<float>(time * 1e9);
                float time_ns_be = swap_endian_float(time_ns);
                file.write(reinterpret_cast<const char*>(&time_ns_be), sizeof(float));
            }
        }
        
        // Export relativistic kinetic energy in eV
        file << "\nSCALARS kinetic_energy_eV float\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            double mass = particle.m(); // mass in kg
            for (size_t j = 0; j < particle.traj_size(); j++) {
                ParticleP3D traj_point = particles->trajectory_point(i, j);
                // Velocity components: vx=index 2, vy=index 4, vz=index 6 in m/s
                double vx = traj_point(2);
                double vy = traj_point(4);
                double vz = traj_point(6);
                double v_mag = sqrt(vx*vx + vy*vy + vz*vz);
                
                // Relativistic kinetic energy calculation
                const double c = SPEED_C; // speed of light in m/s
                double beta = v_mag / SPEED_C;
                double beta2 = v_mag*v_mag / SPEED_C2;
                double gamma = 1.0 / sqrt(1.0 - beta*beta);
                
                // Relativistic kinetic energy: KE = (γ - 1)mc²
                double kinetic_energy_J = (gamma - 1.0) * mass * c * c;
                
                // Convert to eV (1 eV = 1.602176634e-19 J)
                float kinetic_energy_eV = static_cast<float>(kinetic_energy_J / CHARGE_E);
                float kinetic_energy_eV_be = swap_endian_float(kinetic_energy_eV);
                file.write(reinterpret_cast<const char*>(&kinetic_energy_eV_be), sizeof(float));
            }
        }
        
        // Export particle properties as cell data
        file << "\nCELL_DATA " << total_lines << "\n";
        
        // Export particle IDs
        file << "SCALARS particle_id int\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            int id_be = swap_endian_int(static_cast<int>(i));
            file.write(reinterpret_cast<const char*>(&id_be), sizeof(int));
        }
        
        // Export particle current in Amperes [A]
        file << "\nSCALARS current_A float\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            float current_A = static_cast<float>(particle.IQ());
            float current_A_be = swap_endian_float(current_A);
            file.write(reinterpret_cast<const char*>(&current_A_be), sizeof(float));
        }
        
        // Export particle charge in elementary charge units [e]
        file << "\nSCALARS charge_e float\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            // Convert from Coulombs to elementary charge units
            float charge_e = static_cast<float>(particle.q() / CHARGE_E);
            float charge_e_be = swap_endian_float(charge_e);
            file.write(reinterpret_cast<const char*>(&charge_e_be), sizeof(float));
        }
        
        // Export particle mass in atomic mass units [amu]
        file << "\nSCALARS mass_amu float\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            // Convert from kg to atomic mass units
            float mass_amu = static_cast<float>(particle.m() / MASS_U);
            float mass_amu_be = swap_endian_float(mass_amu);
            file.write(reinterpret_cast<const char*>(&mass_amu_be), sizeof(float));
        }

        // Export particle status as cell data
        file << "\nSCALARS particle_status int\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            int status_be = swap_endian_int(static_cast<int>(particle.get_status()));
            file.write(reinterpret_cast<const char*>(&status_be), sizeof(int));
        }

        // Export particle generation as cell data
        file << "\nSCALARS particle_generation int\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            int gen_be = swap_endian_int(particle.gen());
            file.write(reinterpret_cast<const char*>(&gen_be), sizeof(int));
        }

        // Export particle species as family-neutral particle_kind ids
        file << "\nSCALARS particle_kind int\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < particles->size(); i++) {
            Particle3D &particle = particles->particle(i);
            int kind_be = swap_endian_int(
                static_cast<int>(identify_particle_species(particle.m(), particle.q(), ion_mass_u))
            );
            file.write(reinterpret_cast<const char*>(&kind_be), sizeof(int));
        }

        // Beamlet of origin, as cell data: one value per trajectory, which is the natural
        // association because a particle belongs to one beamlet for its whole life. In
        // ParaView: colour by beamlet_id for the whole pattern at once, or Threshold on it
        // (or on beamlet_row / beamlet_column) to isolate one beamlet, one row or one
        // column. The arrays are written only when the case actually has beamlet tags, so
        // a file from a run without them is byte for byte what it was before.
        if (!beamlet_tags.empty()) {
            // Birth point of each particle, needed to place secondaries created during
            // tracing, which no emitter tagged.
            std::vector<int> beamlet_id(particles->size(), 0);
            std::vector<int> beamlet_column(particles->size(), 0);
            std::vector<int> beamlet_row(particles->size(), 0);
            for (size_t i = 0; i < particles->size(); i++) {
                Particle3D &particle = particles->particle(i);
                double birth_x = 0.0;
                double birth_y = 0.0;
                if (particle.traj_size() > 0) {
                    ParticleP3D birth = particles->trajectory_point(i, 0);
                    birth_x = birth(1);
                    birth_y = birth(3);
                }
                const BeamletTag* tag = beamletOfParticle(i, birth_x, birth_y);
                if (tag != NULL) {
                    beamlet_id[i] = tag->index;
                    beamlet_column[i] = tag->column;
                    beamlet_row[i] = tag->row;
                }
            }

            file << "\nSCALARS beamlet_id int\n";
            file << "LOOKUP_TABLE default\n";
            for (size_t i = 0; i < particles->size(); i++) {
                int be = swap_endian_int(beamlet_id[i]);
                file.write(reinterpret_cast<const char*>(&be), sizeof(int));
            }

            file << "\nSCALARS beamlet_column int\n";
            file << "LOOKUP_TABLE default\n";
            for (size_t i = 0; i < particles->size(); i++) {
                int be = swap_endian_int(beamlet_column[i]);
                file.write(reinterpret_cast<const char*>(&be), sizeof(int));
            }

            file << "\nSCALARS beamlet_row int\n";
            file << "LOOKUP_TABLE default\n";
            for (size_t i = 0; i < particles->size(); i++) {
                int be = swap_endian_int(beamlet_row[i]);
                file.write(reinterpret_cast<const char*>(&be), sizeof(int));
            }

            // Where the beamlet's own axis is, so a trajectory can be read as a
            // displacement from its aperture without joining the run's JSON to the VTK.
            file << "\nSCALARS beamlet_center_x_mm float\n";
            file << "LOOKUP_TABLE default\n";
            for (size_t i = 0; i < particles->size(); i++) {
                const int slot = beamlet_id[i];
                float value = 0.0f;
                for (size_t s = 0; s < beamlet_tags.size(); s++) {
                    if (beamlet_tags[s].index == slot) {
                        value = static_cast<float>(1.0e3 * beamlet_tags[s].centerX);
                        break;
                    }
                }
                float be = swap_endian_float(value);
                file.write(reinterpret_cast<const char*>(&be), sizeof(float));
            }

            file << "\nSCALARS beamlet_center_y_mm float\n";
            file << "LOOKUP_TABLE default\n";
            for (size_t i = 0; i < particles->size(); i++) {
                const int slot = beamlet_id[i];
                float value = 0.0f;
                for (size_t s = 0; s < beamlet_tags.size(); s++) {
                    if (beamlet_tags[s].index == slot) {
                        value = static_cast<float>(1.0e3 * beamlet_tags[s].centerY);
                        break;
                    }
                }
                float be = swap_endian_float(value);
                file.write(reinterpret_cast<const char*>(&be), sizeof(float));
            }
        }

        file.close();
        ibsimu.message(1) << "Particle trajectories exported to binary VTK file: " << vtkfile << endl;
        ibsimu.message(1) << "  Format: Binary VTK with big-endian byte order" << endl;
        ibsimu.message(1) << "  Point data: velocity, time, kinetic energy" << endl;
        ibsimu.message(1) << "  Cell data: particle ID, current, charge, mass, status, generation, species" << endl;
        if (!beamlet_tags.empty()) {
            ibsimu.message(1) << "  Cell data: beamlet_id, beamlet_column, beamlet_row, "
                              << "beamlet_center_x_mm, beamlet_center_y_mm over "
                              << beamlet_tags.size() << " beamlet(s)" << endl;
        }
        
    } catch (const Error& e) {
        ibsimu.message(1) << "Error exporting trajectories to VTK" << endl;
    }
}
