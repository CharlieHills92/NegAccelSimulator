/*
 * SurfaceEventLedger.h
 *
 * Signed surface energy/charge balance for the grid power diagnostics.
 *
 * The grid power/current columns are built by a post-pass over the particle database
 * (DiagnosticsManager::analyzeGridPowerLoads), which sees only where trajectories ENDED.
 * That post-pass is blind to two things:
 *
 *   1. Energy carried AWAY from an electrode by the secondaries it emitted. The incident
 *      particle deposits its full kinetic energy at the impact point, and each secondary
 *      then deposits its own energy wherever it lands -- so the energy the secondaries
 *      took with them is counted twice, once at the source electrode and once at the
 *      destination. For electron impacts this is of order 20% of the incident energy.
 *   2. Charge carried away from an electrode by those same secondaries. Without it the
 *      Current column is an "arriving current", not the net drain current a power supply
 *      would actually measure.
 *
 * Both are the same conservation statement applied to two conserved quantities, so both
 * are fixed by one ledger of signed events: arrivals positive, emissions negative. The
 * emitting side is only visible from inside the surface-collision callback, which is why
 * the callback records here rather than the diagnostics recomputing it.
 *
 * THREADING INVARIANT: the surface-collision path is single-threaded today --
 * ManageSimulation::initializeIbsimu() forces ibsimu.set_thread_count(1) whenever
 * getIncludeSurfaceCollisions() > 0 (ManageSimulation_New.cpp:101-103). The mutex below
 * is future-proofing only. If that thread-count heuristic is ever relaxed, note that
 * THCallback_surf_EAMCC's single shared _rng (MTRandom) would race as well, and the
 * callback's _current_hit_solid scratch member would need to become thread-local.
 */

#ifndef SURFACEEVENTLEDGER_H_
#define SURFACEEVENTLEDGER_H_

#include <cstddef>
#include <mutex>
#include <vector>

#include "vec3d.hpp"

/**
 * @brief One signed contribution to a surface's energy/charge balance.
 *
 * Sign convention: positive means "delivered to the surface", negative means "taken away
 * from the surface". The ledger currently only holds emissions (all negative), because
 * arrivals are already enumerated from the particle database; storing the same convention
 * for both keeps a future unification straightforward.
 */
struct SurfaceEvent {
    uint32_t solid;     ///< Geometry::inside() index of the emitting/receiving solid (>= 7)
    Vec3D loc;          ///< Emission point. Required at this granularity for the per-triangle map.
    int kind;           ///< particle_kind of the emitted particle (as int, to avoid a funct.h dep)
    int gen;            ///< RAW generation of the emitted particle (surface secondaries carry >= 101)
    double power_W;     ///< Energy rate: negative for emission
    double current_A;   ///< Charge rate, signed by the particle's own charge: 0 for neutrals
};

/**
 * @brief Append-only collection of SurfaceEvents for one trajectory-tracking pass.
 *
 * Written during pdb->iterate_trajectories() by THCallback_surf_EAMCC, then read by
 * DiagnosticsManager after tracking has finished. Reads and writes never overlap, which
 * is why events() may hand out a bare reference.
 */
class SurfaceEventLedger {
public:
    SurfaceEventLedger() {}

    /**
     * @brief Record one emission.
     *
     * Callers must pass power_W <= 0 and the correspondingly signed current_A; the ledger
     * does not flip signs, so that a future arrival-side user can share the same call.
     */
    void record(uint32_t solid, const Vec3D& loc, int kind, int gen,
                double power_W, double current_A);

    /**
     * @brief Drop every recorded event.
     *
     * Must be called before each tracking pass. The particle database and geometry are
     * rebuilt between JEXT current-matching passes, so events from an earlier pass refer
     * to a geometry that no longer exists.
     */
    void clear();

    /// Events recorded so far. Only valid to call once tracking has stopped.
    const std::vector<SurfaceEvent>& events() const { return _events; }

    size_t size() const;
    bool empty() const;

private:
    // Non-copyable: DiagnosticsManager holds a bare pointer to the owner's instance.
    SurfaceEventLedger(const SurfaceEventLedger&);
    SurfaceEventLedger& operator=(const SurfaceEventLedger&);

    std::vector<SurfaceEvent> _events;
    mutable std::mutex _mutex;
};

#endif /* SURFACEEVENTLEDGER_H_ */
