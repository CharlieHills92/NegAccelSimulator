#!/usr/bin/env python3
"""
Plot particle trajectories from VTK file at a specific z-plane.

This script reads the VTK trajectory file exported by ParticleManager::exportTrajectoriesToVTK
and plots the particle positions in the x-y plane at a specified z location.

Usage:
    python plot_trajectories_vtk.py <vtk_file> [z_plane]

Example:
    python plot_trajectories_vtk.py output/simulation_trajectories.vtk 0.565
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os


def read_vtk_polydata(filename):
    """
    Read VTK POLYDATA file containing particle trajectories.
    
    Returns:
        points: numpy array of shape (n_points, 3) with x, y, z coordinates
        lines: list of arrays, each containing point indices for one trajectory
        point_data: dictionary of point data arrays (velocity, time, energy)
        cell_data: dictionary of cell data arrays (particle properties)
    """
    print(f"Reading VTK file: {filename}")
    
    with open(filename, 'r') as f:
        lines_text = f.readlines()
    
    # Parse VTK file
    i = 0
    points = []
    lines = []
    point_data = {}
    cell_data = {}
    
    # Skip header (first 4 lines)
    while i < len(lines_text) and not lines_text[i].startswith('POINTS'):
        i += 1
    
    # Read POINTS
    if i < len(lines_text):
        parts = lines_text[i].split()
        n_points = int(parts[1])
        print(f"Reading {n_points} points...")
        i += 1
        
        for j in range(n_points):
            x, y, z = map(float, lines_text[i].split())
            points.append([x, y, z])
            i += 1
    
    points = np.array(points)
    
    # Read LINES
    while i < len(lines_text) and not lines_text[i].startswith('LINES'):
        i += 1
    
    if i < len(lines_text):
        parts = lines_text[i].split()
        n_lines = int(parts[1])
        print(f"Reading {n_lines} trajectory lines...")
        i += 1
        
        for j in range(n_lines):
            line_data = list(map(int, lines_text[i].split()))
            n_points_in_line = line_data[0]
            point_indices = line_data[1:]
            lines.append(np.array(point_indices))
            i += 1
    
    # Read POINT_DATA
    while i < len(lines_text) and not lines_text[i].startswith('POINT_DATA'):
        i += 1
    
    if i < len(lines_text):
        n_point_data = int(lines_text[i].split()[1])
        i += 1
        
        # Read scalar fields
        while i < len(lines_text) and lines_text[i].startswith('SCALARS'):
            scalar_name = lines_text[i].split()[1]
            i += 1  # Skip LOOKUP_TABLE line
            i += 1
            
            scalar_data = []
            for j in range(n_point_data):
                if i >= len(lines_text) or lines_text[i].startswith('SCALARS') or lines_text[i].startswith('CELL_DATA'):
                    break
                scalar_data.append(float(lines_text[i].strip()))
                i += 1
            
            point_data[scalar_name] = np.array(scalar_data)
            print(f"  Read point data: {scalar_name} ({len(scalar_data)} values)")
    
    # Read CELL_DATA
    while i < len(lines_text) and not lines_text[i].startswith('CELL_DATA'):
        i += 1
    
    if i < len(lines_text):
        n_cell_data = int(lines_text[i].split()[1])
        i += 1
        
        # Read scalar fields
        while i < len(lines_text) and lines_text[i].startswith('SCALARS'):
            scalar_name = lines_text[i].split()[1]
            scalar_type = lines_text[i].split()[2]
            i += 1  # Skip LOOKUP_TABLE line
            i += 1
            
            scalar_data = []
            for j in range(n_cell_data):
                if i >= len(lines_text):
                    break
                if scalar_type == 'int':
                    scalar_data.append(int(lines_text[i].strip()))
                else:
                    scalar_data.append(float(lines_text[i].strip()))
                i += 1
            
            cell_data[scalar_name] = np.array(scalar_data)
            print(f"  Read cell data: {scalar_name} ({len(scalar_data)} values)")
    
    print(f"Successfully read {len(points)} points and {len(lines)} trajectories")
    return points, lines, point_data, cell_data


def interpolate_trajectory_at_z(points, line_indices, z_plane, tolerance=1e-6):
    """
    Find the x, y position where a trajectory crosses a specific z-plane.
    
    Args:
        points: numpy array of all points (n_points, 3)
        line_indices: array of point indices for one trajectory
        z_plane: z coordinate of the plane
        tolerance: tolerance for z matching (in meters)
    
    Returns:
        (x, y) tuple if crossing found, None otherwise
    """
    traj_points = points[line_indices]
    
    # Find segments that cross the z_plane
    for i in range(len(traj_points) - 1):
        p1 = traj_points[i]
        p2 = traj_points[i + 1]
        
        z1, z2 = p1[2], p2[2]
        
        # Check if z_plane is between z1 and z2
        if (z1 <= z_plane <= z2) or (z2 <= z_plane <= z1):
            # Linear interpolation
            if abs(z2 - z1) < 1e-12:  # Avoid division by zero
                # Segment is parallel to z_plane, take midpoint
                x = (p1[0] + p2[0]) / 2
                y = (p1[1] + p2[1]) / 2
            else:
                # Interpolation factor
                t = (z_plane - z1) / (z2 - z1)
                x = p1[0] + t * (p2[0] - p1[0])
                y = p1[1] + t * (p2[1] - p1[1])
            
            return (x, y)
    
    return None


def plot_trajectories_at_z_plane(vtk_file, z_plane=0.565, output_file=None):
    """
    Plot particle trajectories at a specific z-plane in the x-y plane.
    
    Args:
        vtk_file: path to VTK trajectory file
        z_plane: z coordinate of the plane (in meters)
        output_file: optional output file for saving the plot
    """
    # Read VTK file
    points, lines, point_data, cell_data = read_vtk_polydata(vtk_file)
    
    # Extract trajectory crossings at z_plane
    crossings = []
    particle_ids = []
    
    print(f"\nFinding trajectory crossings at z = {z_plane} m ({z_plane*1000} mm)...")
    
    for i, line_indices in enumerate(lines):
        crossing = interpolate_trajectory_at_z(points, line_indices, z_plane)
        if crossing is not None:
            crossings.append(crossing)
            particle_ids.append(i)
    
    print(f"Found {len(crossings)} particles crossing z = {z_plane} m")
    
    if len(crossings) == 0:
        print("No trajectories found at the specified z-plane!")
        print(f"z-range in data: [{points[:, 2].min():.6f}, {points[:, 2].max():.6f}] m")
        return
    
    # Convert to arrays
    crossings = np.array(crossings)
    x_vals = crossings[:, 0]
    y_vals = crossings[:, 1]
    
    # Create figure with two subplots
    fig = plt.figure(figsize=(16, 7))
    
    # Subplot 1: Scatter plot of particle positions
    ax1 = fig.add_subplot(121)
    
    # Get particle status if available for coloring
    if 'particle_status' in cell_data and len(particle_ids) > 0:
        status_values = cell_data['particle_status'][particle_ids]
        scatter = ax1.scatter(x_vals * 1000, y_vals * 1000, c=status_values, 
                             cmap='viridis', s=20, alpha=0.6)
        cbar = plt.colorbar(scatter, ax=ax1)
        cbar.set_label('Particle Status', rotation=270, labelpad=20)
    else:
        ax1.scatter(x_vals * 1000, y_vals * 1000, s=20, alpha=0.6, c='blue')
    
    ax1.set_xlabel('x [mm]', fontsize=12)
    ax1.set_ylabel('y [mm]', fontsize=12)
    ax1.set_title(f'Particle Positions at z = {z_plane*1000:.1f} mm\n({len(crossings)} particles)', 
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal', adjustable='box')
    
    # Add statistics to the plot
    x_mean, y_mean = x_vals.mean(), y_vals.mean()
    x_std, y_std = x_vals.std(), y_vals.std()
    ax1.axvline(x_mean * 1000, color='red', linestyle='--', alpha=0.5, label=f'x_mean = {x_mean*1000:.2f} mm')
    ax1.axhline(y_mean * 1000, color='red', linestyle='--', alpha=0.5, label=f'y_mean = {y_mean*1000:.2f} mm')
    ax1.legend(loc='upper right', fontsize=10)
    
    # Subplot 2: 2D histogram / density plot
    ax2 = fig.add_subplot(122)
    
    hist, xedges, yedges = np.histogram2d(x_vals * 1000, y_vals * 1000, bins=50)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    
    im = ax2.imshow(hist.T, origin='lower', extent=extent, aspect='auto', cmap='hot', interpolation='bilinear')
    cbar2 = plt.colorbar(im, ax=ax2)
    cbar2.set_label('Particle Count', rotation=270, labelpad=20)
    
    ax2.set_xlabel('x [mm]', fontsize=12)
    ax2.set_ylabel('y [mm]', fontsize=12)
    ax2.set_title(f'Beam Density at z = {z_plane*1000:.1f} mm', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, color='white')
    
    # Add beam statistics as text
    stats_text = f'Statistics:\n'
    stats_text += f'N = {len(crossings)}\n'
    stats_text += f'<x> = {x_mean*1000:.3f} mm\n'
    stats_text += f'<y> = {y_mean*1000:.3f} mm\n'
    stats_text += f'σ_x = {x_std*1000:.3f} mm\n'
    stats_text += f'σ_y = {y_std*1000:.3f} mm\n'
    stats_text += f'x: [{x_vals.min()*1000:.2f}, {x_vals.max()*1000:.2f}] mm\n'
    stats_text += f'y: [{y_vals.min()*1000:.2f}, {y_vals.max()*1000:.2f}] mm'
    
    ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes, 
            fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    # Save or show
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\nPlot saved to: {output_file}")
    else:
        plt.show()
    
    # Print statistics
    print(f"\nBeam Statistics at z = {z_plane*1000:.1f} mm:")
    print(f"  Number of particles: {len(crossings)}")
    print(f"  <x> = {x_mean*1000:.3f} mm,  σ_x = {x_std*1000:.3f} mm")
    print(f"  <y> = {y_mean*1000:.3f} mm,  σ_y = {y_std*1000:.3f} mm")
    print(f"  x range: [{x_vals.min()*1000:.3f}, {x_vals.max()*1000:.3f}] mm")
    print(f"  y range: [{y_vals.min()*1000:.3f}, {y_vals.max()*1000:.3f}] mm")
    print(f"  Beam width (2σ): {2*x_std*1000:.3f} mm × {2*y_std*1000:.3f} mm")


def plot_3d_trajectories(vtk_file, max_trajectories=100, output_file=None):
    """
    Create a 3D plot of particle trajectories.
    
    Args:
        vtk_file: path to VTK trajectory file
        max_trajectories: maximum number of trajectories to plot (for performance)
        output_file: optional output file for saving the plot
    """
    # Read VTK file
    points, lines, point_data, cell_data = read_vtk_polydata(vtk_file)
    
    # Create 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot subset of trajectories
    n_plot = min(len(lines), max_trajectories)
    print(f"\nPlotting {n_plot} trajectories in 3D...")
    
    for i in range(n_plot):
        line_indices = lines[i]
        traj_points = points[line_indices]
        
        # Color by particle status if available
        if 'particle_status' in cell_data:
            status = cell_data['particle_status'][i]
            color = plt.cm.viridis(status / 4.0)  # Normalize to 0-1
        else:
            color = 'blue'
        
        ax.plot(traj_points[:, 0] * 1000, 
               traj_points[:, 1] * 1000, 
               traj_points[:, 2] * 1000,
               alpha=0.5, linewidth=0.5, color=color)
    
    ax.set_xlabel('x [mm]', fontsize=12)
    ax.set_ylabel('y [mm]', fontsize=12)
    ax.set_zlabel('z [mm]', fontsize=12)
    ax.set_title(f'3D Particle Trajectories\n({n_plot} particles shown)', 
                fontsize=14, fontweight='bold')
    
    # Set equal aspect ratio
    max_range = np.array([points[:, 0].max() - points[:, 0].min(),
                         points[:, 1].max() - points[:, 1].min(),
                         points[:, 2].max() - points[:, 2].min()]).max() / 2.0
    
    mid_x = (points[:, 0].max() + points[:, 0].min()) * 0.5
    mid_y = (points[:, 1].max() + points[:, 1].min()) * 0.5
    mid_z = (points[:, 2].max() + points[:, 2].min()) * 0.5
    
    ax.set_xlim((mid_x - max_range) * 1000, (mid_x + max_range) * 1000)
    ax.set_ylim((mid_y - max_range) * 1000, (mid_y + max_range) * 1000)
    ax.set_zlim((mid_z - max_range) * 1000, (mid_z + max_range) * 1000)
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"3D plot saved to: {output_file}")
    else:
        plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    
    vtk_file = sys.argv[1]
    
    if not os.path.exists(vtk_file):
        print(f"Error: File not found: {vtk_file}")
        sys.exit(1)
    
    # Get z_plane from command line or use default
    z_plane = float(sys.argv[2]) if len(sys.argv) > 2 else 0.565
    
    # Generate output filename
    base_name = os.path.splitext(vtk_file)[0]
    output_2d = f"{base_name}_z{z_plane*1000:.0f}mm.png"
    output_3d = f"{base_name}_3d.png"
    
    # Plot trajectories at z-plane
    plot_trajectories_at_z_plane(vtk_file, z_plane, output_file=output_2d)
    
    # Also create 3D visualization
    response = input("\nCreate 3D trajectory plot? (y/n): ")
    if response.lower() == 'y':
        plot_3d_trajectories(vtk_file, max_trajectories=200, output_file=output_3d)
    
    print("\nDone!")
