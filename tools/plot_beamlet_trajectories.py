#!/usr/bin/env python3
"""Trajectory figure for a multi-beamlet run.

Reads the binary VTK polydata written by ParticleManager::exportTrajectoriesToVTK and draws
the traced trajectories in both longitudinal planes plus the transverse pattern, so a
multi-beamlet case can be checked by eye: are all the beamlets there, are they where the
apertures are, and does each one behave like the single-beamlet reference.

    python3 plot_beamlet_trajectories.py QUAD_RUN_DIR CASE_JSON --out fig.png
                                         [--reference SINGLE_RUN_DIR]
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

C_BEAMLET = ["#1f6feb", "#e8590c", "#2f9e44", "#9c36b5"]
C_METAL = "#9aa0a6"
C_METAL_E = "#5f6368"
C_REF = "#868e96"


# ── VTK ──────────────────────────────────────────────────────────────────────

def read_trajectory_vtk(path: Path, with_cell_data: bool = False):
    """Return a list of (n, 3) arrays in metres, one per traced particle.

    The writer emits ASCII headers with big-endian binary blocks between them, so the file
    is walked by locating each header and slicing the bytes that follow it.

    With ``with_cell_data``, also return the per-trajectory integer cell arrays, which is
    where ``beamlet_id`` lives. They are returned in step with the tracks: a one-point
    trajectory is dropped from both, so index i of the list and index i of each array are
    still the same particle.
    """
    raw = path.read_bytes()

    m = re.search(rb"POINTS (\d+) float\n", raw)
    if m is None:
        raise ValueError(f"{path}: no POINTS header")
    n_points = int(m.group(1))
    start = m.end()
    points = np.frombuffer(raw, dtype=">f4", count=3 * n_points, offset=start)
    points = points.reshape(n_points, 3).astype(float)

    tail_at = start + 12 * n_points
    m = re.compile(rb"\nLINES (\d+) (\d+)\n").search(raw, tail_at)
    if m is None:
        raise ValueError(f"{path}: no LINES header")
    n_lines, n_line_data = int(m.group(1)), int(m.group(2))
    line_data = np.frombuffer(raw, dtype=">i4", count=n_line_data, offset=m.end())

    tracks, kept, cursor = [], [], 0
    for line_index in range(n_lines):
        count = int(line_data[cursor])
        idx = line_data[cursor + 1: cursor + 1 + count]
        cursor += 1 + count
        if count >= 2:
            tracks.append(points[idx])
            kept.append(line_index)
    if not with_cell_data:
        return tracks

    # CELL_DATA blocks: "SCALARS <name> <type>", "LOOKUP_TABLE default", then n_lines
    # big-endian values. Walk them by header so an array added later is picked up without
    # changing this reader.
    cell: dict[str, np.ndarray] = {}
    keep = np.asarray(kept, dtype=int)
    at = raw.find(b"\nCELL_DATA ")
    if at >= 0:
        pattern = re.compile(rb"SCALARS (\S+) (int|float)\nLOOKUP_TABLE default\n")
        cursor = at
        while True:
            m = pattern.search(raw, cursor)
            if m is None:
                break
            name = m.group(1).decode()
            dtype = ">i4" if m.group(2) == b"int" else ">f4"
            values = np.frombuffer(raw, dtype=dtype, count=n_lines, offset=m.end())
            cell[name] = values[keep] if keep.size else values[:0]
            cursor = m.end() + 4 * n_lines
    return tracks, cell


def find_trajectory_file(target: Path) -> Path:
    """Accept the .vtk itself, or a directory to search for the largest one."""
    if target.is_file():
        return target
    hits = sorted(target.rglob("*_trajectories.vtk"))
    if not hits:
        raise SystemExit(f"no *_trajectories.vtk under {target}")
    return max(hits, key=lambda p: p.stat().st_size)


# ── geometry ─────────────────────────────────────────────────────────────────

def load_case(path: Path):
    case = json.loads(path.read_text())
    geometry = case["geometry"]
    solids = [s for s in geometry.get("solids", []) if s.get("kind", "solid") == "solid"]
    centres = []
    for source in case["particleSources"]:
        c = source["uniform"]["centerMeters"]
        centres.append((c[0], c[1]))
    return case, geometry, solids, centres


def aperture_nodes(pattern: dict, axis: str) -> list[float]:
    """Transverse positions of the whole apertures the geometry opens along one axis."""
    if not isinstance(pattern, dict) or pattern.get("layout", "single") == "single":
        return [float(pattern.get(f"{axis}OffsetMeters", 0.0)) if pattern else 0.0]
    count = int(pattern.get("countX" if axis == "x" else "countY", 1))
    pitch = float(pattern.get("pitchXMeters" if axis == "x" else "pitchYMeters", 0.0))
    offset = float(pattern.get(f"{axis}OffsetMeters", 0.0))
    # Exactly `count` nodes, floor-centred for an even count -- the single lattice rule
    # shared with axis_node_indices() in the beamlet expansion and with the patterned
    # region in evaluateGeneratedSolidProfile.
    start = -(count // 2)
    return [offset + k * pitch for k in range(start, start + count)]


def draw_grids(ax, solids, axis: str, lim: float):
    """Electrodes in the (z, x) or (z, y) plane, as metal between the apertures."""
    for solid in solids:
        z = np.asarray(solid["zProfileMeters"], dtype=float) * 1e3
        r = np.asarray(solid["rProfileMeters"], dtype=float) * 1e3
        z0, z1, r_ap = z.min(), z.max(), r.min()
        nodes = [n * 1e3 for n in aperture_nodes(solid.get("aperturePattern", {}), axis)
                 if abs(n * 1e3) < lim + r_ap]
        edges = sorted([-lim] + [e for n in nodes for e in (n - r_ap, n + r_ap)] + [lim])
        for lo, hi in zip(edges[0::2], edges[1::2]):
            lo, hi = max(lo, -lim), min(hi, lim)
            if hi > lo:
                ax.add_patch(Rectangle((z0, lo), z1 - z0, hi - lo, facecolor=C_METAL,
                                       edgecolor=C_METAL_E, lw=0.4, zorder=4, clip_on=True))
        ax.text(0.5 * (z0 + z1), lim * 0.93, solid.get("name", ""), ha="center", va="top",
                fontsize=6.5, rotation=90, color="#333", zorder=6, clip_on=True)


def beamlet_of(track, centres, axis_pair=(0, 1)) -> int:
    """Nearest source centre to the trajectory's first point.

    The fallback for a file written before the runtime stamped ``beamlet_id`` on each
    trajectory. It guesses from the emission point, which is right for a particle the
    source emitted and wrong for one created during tracing.
    """
    x0, y0 = track[0][axis_pair[0]], track[0][axis_pair[1]]
    d = [(x0 - cx) ** 2 + (y0 - cy) ** 2 for cx, cy in centres]
    return int(np.argmin(d))


def group_by_beamlet(tracks, cell_data, centres) -> dict[int, list]:
    """Tracks grouped by beamlet, 0-based, using the file's own tag where it has one."""
    groups: dict[int, list] = {i: [] for i in range(len(centres))}
    tags = cell_data.get("beamlet_id") if cell_data else None
    if tags is not None and len(tags) == len(tracks):
        for track, tag in zip(tracks, tags):
            # beamlet_id is 1-based and is 0 only when the runtime had no tag to give.
            index = int(tag) - 1
            groups.setdefault(index if index >= 0 else 0, []).append(track)
        return groups
    for track in tracks:
        groups[beamlet_of(track, centres)].append(track)
    return groups


# ── figure ───────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("case_json")
    ap.add_argument("--out", default="beamlet_trajectories.png")
    ap.add_argument("--reference", help="single-beamlet run directory to overlay")
    ap.add_argument("--max-tracks", type=int, default=260, help="per beamlet, for legibility")
    args = ap.parse_args()

    case, geometry, solids, centres = load_case(Path(args.case_json))
    tracks, cell_data = read_trajectory_vtk(
        find_trajectory_file(Path(args.run_dir)), with_cell_data=True)
    ref_tracks = (read_trajectory_vtk(find_trajectory_file(Path(args.reference)))
                  if args.reference else [])

    lim = 0.5e3 * float(geometry["domain"]["xSizeMeters"])
    z_max = 1e3 * (float(geometry["domain"]["zStartMeters"])
                   + float(geometry["domain"]["zSizeMeters"]))

    groups = group_by_beamlet(tracks, cell_data, centres)
    print("beamlet assignment: %s"
          % ("beamlet_id from the VTK" if "beamlet_id" in cell_data
             else "nearest source centre (file carries no beamlet_id)"))

    rng = np.random.default_rng(7)

    def thinned(items):
        if len(items) <= args.max_tracks:
            return items
        return [items[i] for i in rng.choice(len(items), args.max_tracks, replace=False)]

    fig = plt.figure(figsize=(15.0, 13.0))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.0, 1.45], hspace=0.34, wspace=0.28)

    # (a)/(b) longitudinal views, full column
    for row, (comp, label) in enumerate(((0, "x"), (1, "y"))):
        ax = fig.add_subplot(gs[row, :])
        draw_grids(ax, solids, label, lim)
        for i, items in groups.items():
            for t in thinned(items):
                ax.plot(t[:, 2] * 1e3, t[:, comp] * 1e3, lw=0.35, alpha=0.55,
                        color=C_BEAMLET[i % len(C_BEAMLET)], zorder=3)
        for i, c in enumerate(centres):
            ax.plot([], [], color=C_BEAMLET[i % len(C_BEAMLET)], lw=1.6,
                    label="beamlet %d  (%.1f, %.1f) mm" % (i + 1, c[0] * 1e3, c[1] * 1e3))
        ax.set_xlim(-5, z_max + 5)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("z  [mm]")
        ax.set_ylabel(f"{label}  [mm]")
        ax.set_title("(%s) trajectories, %s-z plane" % ("ab"[row], label), loc="left",
                     fontsize=11, fontweight="bold")
        if row == 0:
            ax.legend(loc="upper right", fontsize=7.5, ncol=4, framealpha=0.9)

    # (c) extraction zoom
    ax = fig.add_subplot(gs[2, 0])
    draw_grids(ax, solids, "x", lim)
    for i, items in groups.items():
        for t in thinned(items):
            ax.plot(t[:, 2] * 1e3, t[:, 0] * 1e3, lw=0.4, alpha=0.6,
                    color=C_BEAMLET[i % len(C_BEAMLET)], zorder=3)
    ax.set_xlim(-1, 40)
    ax.set_ylim(-26, 26)
    ax.set_xlabel("z  [mm]")
    ax.set_ylabel("x  [mm]")
    ax.set_title("(c) extraction region", loc="left", fontsize=11, fontweight="bold")

    # (d) transverse pattern at the exit
    ax = fig.add_subplot(gs[2, 1])
    for i, items in groups.items():
        pts = np.array([t[-1] for t in items])
        if pts.size:
            ax.scatter(pts[:, 0] * 1e3, pts[:, 1] * 1e3, s=1.6, alpha=0.45,
                       color=C_BEAMLET[i % len(C_BEAMLET)], linewidths=0)
    pg = next((s for s in solids if s["name"] == "PG"), None)
    if pg is not None:
        r_ap = 1e3 * min(pg["rProfileMeters"])
        for nx in aperture_nodes(pg.get("aperturePattern", {}), "x"):
            for ny in aperture_nodes(pg.get("aperturePattern", {}), "y"):
                if abs(nx * 1e3) < lim + r_ap and abs(ny * 1e3) < lim + r_ap:
                    ax.add_patch(plt.Circle((nx * 1e3, ny * 1e3), r_ap, fill=False,
                                            ec="#adb5bd", lw=0.6, ls="--", zorder=1,
                                            clip_on=True))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x  [mm]")
    ax.set_ylabel("y  [mm]")
    ax.set_title("(d) final positions (dashed = PG apertures)", loc="left",
                 fontsize=11, fontweight="bold")

    # (e) the single-beamlet reference against beamlet 1
    ax = fig.add_subplot(gs[2, 2])
    draw_grids(ax, solids, "x", lim)
    if ref_tracks:
        for t in thinned(ref_tracks):
            ax.plot(t[:, 2] * 1e3, t[:, 0] * 1e3, lw=0.4, alpha=0.5, color=C_REF, zorder=2)
        ax.plot([], [], color=C_REF, lw=1.6, label="1 beamlet alone")
    for t in thinned(groups[0]):
        ax.plot(t[:, 2] * 1e3, t[:, 0] * 1e3, lw=0.4, alpha=0.5, color=C_BEAMLET[0], zorder=3)
    ax.plot([], [], color=C_BEAMLET[0], lw=1.6, label="same beamlet, in the 2x2")
    ax.set_xlim(-1, 40)
    ax.set_ylim(-26, 26)
    ax.set_xlabel("z  [mm]")
    ax.set_ylabel("x  [mm]")
    ax.set_title("(e) effect of the neighbours", loc="left", fontsize=11, fontweight="bold")
    ax.legend(loc="upper right", fontsize=7.5, framealpha=0.9)

    fig.suptitle("2 x 2 beamlet extraction -- %s" % case["metadata"].get("caseTag", ""),
                 fontsize=13, fontweight="bold", y=0.995)
    fig.savefig(args.out, dpi=170, facecolor="white")
    print("wrote", args.out)
    for i, items in groups.items():
        print("  beamlet %d: %5d traced trajectories" % (i + 1, len(items)))
    if ref_tracks:
        print("  reference: %5d traced trajectories" % len(ref_tracks))


if __name__ == "__main__":
    main()
