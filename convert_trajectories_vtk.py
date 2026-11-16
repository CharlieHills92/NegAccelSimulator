#!/usr/bin/env python3
"""Convert an IBSimu trajectory VTK file to a final-particle text table.

Usage:
    python3 convert_trajectories_vtk.py simulation_trajectories.vtk
    python3 convert_trajectories_vtk.py input.vtk output.txt

The VTK exporter stores one trajectory as one POLYDATA line and stores the
particle ID and status as CELL_DATA.  The converter writes the last point of
each line, so the output contains one row per particle.
"""

import argparse
import struct
from pathlib import Path


STATUS_NAMES = {
    0: "OK",
    1: "OUT",
    2: "COLLIDED",
    6: "STRIPPED",
}


class VtkReader:
    def __init__(self, path):
        self.handle = path.open("rb")

    def line(self):
        value = self.handle.readline()
        if not value:
            return ""
        return value.decode("ascii").strip()

    def next_nonempty_line(self):
        while True:
            raw_value = self.handle.readline()
            if not raw_value:
                return ""
            value = raw_value.decode("ascii").strip()
            if value:
                return value

    def values(self, count, scalar_type, binary):
        scalar_type = scalar_type.lower()
        if scalar_type in ("int", "integer"):
            format_code = "i"
        elif scalar_type in ("float", "double"):
            format_code = "f" if scalar_type == "float" else "d"
        else:
            raise ValueError(f"Unsupported VTK scalar type: {scalar_type}")

        if binary:
            size = struct.calcsize(">" + format_code) * count
            raw = self.handle.read(size)
            if len(raw) != size:
                raise ValueError("Unexpected end of file in VTK scalar data")
            return list(struct.unpack(">" + format_code * count, raw))

        result = []
        while len(result) < count:
            text = self.line()
            if not text:
                raise ValueError("Unexpected end of file in VTK scalar data")
            result.extend(
                int(token) if format_code == "i" else float(token)
                for token in text.split()
            )
        return result[:count]

    def read(self):
        header = [self.line() for _ in range(4)]
        if header[3] != "DATASET POLYDATA":
            raise ValueError("Input must be a VTK DATASET POLYDATA trajectory file")
        binary = header[2].upper() == "BINARY"

        line = self.next_nonempty_line()
        while line and not line.startswith("POINTS "):
            line = self.line()
        if not line:
            raise ValueError("VTK file has no POINTS section")
        point_parts = line.split()
        point_count = int(point_parts[1])
        point_values = self.values(point_count * 3, point_parts[2], binary)
        points = [tuple(point_values[index:index + 3])
                  for index in range(0, len(point_values), 3)]

        line = self.next_nonempty_line()
        while line and not line.startswith("LINES "):
            line = self.line()
        if not line:
            raise ValueError("VTK file has no LINES section")
        line_parts = line.split()
        line_count = int(line_parts[1])
        line_values = self.values(int(line_parts[2]), "int", binary)
        trajectories = []
        offset = 0
        for _ in range(line_count):
            size = int(line_values[offset])
            offset += 1
            trajectories.append([int(value) for value in line_values[offset:offset + size]])
            offset += size

        point_data = {}
        cell_data = {}
        line = self.next_nonempty_line()
        while line and not line.startswith("CELL_DATA "):
            if line.startswith("POINT_DATA "):
                point_data_count = int(line.split()[1])
                line = self.line()
                while line and not line.startswith("CELL_DATA "):
                    if not line.startswith("SCALARS "):
                        line = self.line()
                        continue
                    scalar_parts = line.split()
                    line = self.line()  # LOOKUP_TABLE
                    point_data[scalar_parts[1]] = self.values(
                        point_data_count, scalar_parts[2], binary
                    )
                    line = self.next_nonempty_line()
                break
            line = self.line()

        if line.startswith("CELL_DATA "):
            cell_count = int(line.split()[1])
            line = self.line()
            while line:
                if not line.startswith("SCALARS "):
                    line = self.line()
                    continue
                scalar_parts = line.split()
                scalar_name = scalar_parts[1]
                line = self.line()  # LOOKUP_TABLE
                cell_data[scalar_name] = self.values(cell_count, scalar_parts[2], binary)
                line = self.next_nonempty_line()

        self.handle.close()
        return points, trajectories, point_data, cell_data


def convert(input_path, output_path):
    points, trajectories, point_data, cell_data = VtkReader(input_path).read()
    particle_ids = cell_data.get("particle_id")
    statuses = cell_data.get("particle_status", cell_data.get("status"))
    if particle_ids is None or statuses is None:
        raise ValueError("VTK CELL_DATA must contain particle_id and particle_status")
    if len(particle_ids) != len(trajectories) or len(statuses) != len(trajectories):
        raise ValueError("Particle metadata count does not match trajectory count")

    def point_scalar(name, point_index):
        values = point_data.get(name)
        return "NA" if values is None else f"{values[point_index]:.17g}"

    def cell_scalar(name, particle_index):
        values = cell_data.get(name)
        return "NA" if values is None else f"{values[particle_index]:.17g}"

    with output_path.open("w", encoding="utf-8") as output:
        output.write(
            "particle_id\tx\ty\tz\tstatus\tstatus_id\t"
            "kinetic_energy_eV\tpower_W\tcurrent_A\tcharge_e\tmass_amu\t"
            "velocity_km_s\ttime_ns\tgeneration\tparticle_kind\ttrajectory_points\n"
        )
        for index, trajectory in enumerate(trajectories):
            if not trajectory:
                raise ValueError(f"Trajectory {index} has no points")
            final_point_index = trajectory[-1]
            x, y, z = points[final_point_index]
            status_id = int(statuses[index])
            status = STATUS_NAMES.get(status_id, f"UNKNOWN({status_id})")
            energy = point_data.get("kinetic_energy_eV")
            current = cell_data.get("current_A")
            power = (
                "NA" if energy is None or current is None
                else f"{abs(current[index]) * energy[final_point_index]:.17g}"
            )
            output.write(
                f"{int(particle_ids[index])}\t{x:.17g}\t{y:.17g}\t{z:.17g}\t"
                f"{status}\t{status_id}\t{point_scalar('kinetic_energy_eV', final_point_index)}\t"
                f"{power}\t{cell_scalar('current_A', index)}\t{cell_scalar('charge_e', index)}\t"
                f"{cell_scalar('mass_amu', index)}\t{point_scalar('velocity_magnitude_km_s', final_point_index)}\t"
                f"{point_scalar('time_ns', final_point_index)}\t{cell_scalar('particle_generation', index)}\t"
                f"{cell_scalar('particle_kind', index)}\t{len(trajectory)}\n"
            )

    print(f"Wrote {len(trajectories)} particles to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Write final position and status for every particle in a trajectory VTK file."
    )
    parser.add_argument("input_vtk", type=Path, help="Input trajectory .vtk file")
    parser.add_argument(
        "output_txt", type=Path, nargs="?",
        help="Output .txt file (default: input name with .txt extension)",
    )
    args = parser.parse_args()
    output_path = args.output_txt or args.input_vtk.with_suffix(".txt")
    convert(args.input_vtk, output_path)


if __name__ == "__main__":
    main()