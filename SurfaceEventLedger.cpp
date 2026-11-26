/*
 * SurfaceEventLedger.cpp
 */

#include "SurfaceEventLedger.h"

void SurfaceEventLedger::record(uint32_t solid, const Vec3D& loc, int kind, int gen,
                                double power_W, double current_A) {
    SurfaceEvent ev;
    ev.solid = solid;
    ev.loc = loc;
    ev.kind = kind;
    ev.gen = gen;
    ev.power_W = power_W;
    ev.current_A = current_A;

    std::lock_guard<std::mutex> guard(_mutex);
    _events.push_back(ev);
}

void SurfaceEventLedger::clear() {
    std::lock_guard<std::mutex> guard(_mutex);
    _events.clear();
}

size_t SurfaceEventLedger::size() const {
    std::lock_guard<std::mutex> guard(_mutex);
    return _events.size();
}

bool SurfaceEventLedger::empty() const {
    std::lock_guard<std::mutex> guard(_mutex);
    return _events.empty();
}
