#!/usr/bin/env python3
"""
Analysis script for Poiseuille flow molecular dynamics simulation results.
This script processes velocity profile data to calculate viscosity and create plots.
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os

def read_velocity_profile(filename):
    """Read velocity profile data from file."""
    try:
        data = np.loadtxt(filename, comments='#')
        y_positions = data[:, 0]
        velocities = data[:, 1]
        return y_positions, velocities
    except:
        print(f"Error reading file: {filename}")
        return None, None

def calculate_viscosity_from_profile(y_positions, velocities, applied_force=0.1):
    """Calculate viscosity from velocity profile using finite differences."""
    # Remove zero velocities for gradient calculation
    non_zero_mask = velocities != 0
    if np.sum(non_zero_mask) < 3:
        return 0.0
    
    y_nz = y_positions[non_zero_mask]
    v_nz = velocities[non_zero_mask]
    
    # Calculate velocity gradient using central differences
    if len(y_nz) < 3:
        return 0.0
    
    # Sort by y position
    sorted_indices = np.argsort(y_nz)
    y_sorted = y_nz[sorted_indices]
    v_sorted = v_nz[sorted_indices]
    
    # Calculate gradient
    dy = np.diff(y_sorted)
    dv = np.diff(v_sorted)
    
    # Avoid division by zero
    valid_gradients = dy > 1e-6
    if np.sum(valid_gradients) == 0:
        return 0.0
    
    gradients = dv[valid_gradients] / dy[valid_gradients]
    avg_gradient = np.mean(np.abs(gradients))
    
    # Viscosity = applied_force / velocity_gradient
    viscosity = applied_force / avg_gradient if avg_gradient > 1e-6 else 0.0
    
    return viscosity

def plot_velocity_profile(y_positions, velocities, step, output_dir='.'):
    """Create a plot of the velocity profile."""
    plt.figure(figsize=(10, 6))
    plt.plot(y_positions, velocities, 'bo-', linewidth=2, markersize=4)
    plt.xlabel('Y Position')
    plt.ylabel('Velocity (Z-direction)')
    plt.title(f'Velocity Profile - Step {step}')
    plt.grid(True, alpha=0.3)
    
    # Add theoretical Poiseuille profile for comparison (parabolic)
    y_theory = np.linspace(min(y_positions), max(y_positions), 100)
    # Simple parabolic profile (normalized)
    y_center = (max(y_positions) + min(y_positions)) / 2
    width = max(y_positions) - min(y_positions)
    v_max = max(velocities) if max(velocities) > 0 else 1.0
    v_theory = v_max * (1 - 4 * ((y_theory - y_center) / width)**2)
    v_theory = np.maximum(v_theory, 0)  # Ensure non-negative
    
    plt.plot(y_theory, v_theory, 'r--', linewidth=2, alpha=0.7, 
             label='Theoretical Poiseuille')
    plt.legend()
    
    filename = os.path.join(output_dir, f'velocity_profile_plot_{step}.png')
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    
    return filename

def analyze_all_profiles():
    """Analyze all velocity profile files in the current directory."""
    profile_files = sorted(glob.glob('velocity_profile_*.dat'))
    
    if not profile_files:
        print("No velocity profile files found!")
        return
    
    print(f"Found {len(profile_files)} velocity profile files")
    print("\nAnalysis Results:")
    print("=" * 50)
    print(f"{'Step':<10} {'Viscosity':<15} {'Max Velocity':<15}")
    print("-" * 50)
    
    steps = []
    viscosities = []
    max_velocities = []
    
    for filename in profile_files:
        # Extract step number from filename
        step_str = filename.replace('velocity_profile_', '').replace('.dat', '')
        try:
            step = int(step_str)
        except:
            step = 0
        
        y_pos, velocities = read_velocity_profile(filename)
        if y_pos is not None and velocities is not None:
            viscosity = calculate_viscosity_from_profile(y_pos, velocities)
            max_vel = np.max(velocities)
            
            print(f"{step:<10} {viscosity:<15.6f} {max_vel:<15.6f}")
            
            steps.append(step)
            viscosities.append(viscosity)
            max_velocities.append(max_vel)
            
            # Create plot for this profile
            plot_velocity_profile(y_pos, velocities, step)
    
    if len(steps) > 1:
        # Plot viscosity evolution
        plt.figure(figsize=(12, 5))
        
        plt.subplot(1, 2, 1)
        plt.plot(steps, viscosities, 'bo-', linewidth=2)
        plt.xlabel('Simulation Step')
        plt.ylabel('Viscosity')
        plt.title('Viscosity Evolution')
        plt.grid(True, alpha=0.3)
        
        plt.subplot(1, 2, 2)
        plt.plot(steps, max_velocities, 'ro-', linewidth=2)
        plt.xlabel('Simulation Step')
        plt.ylabel('Maximum Velocity')
        plt.title('Maximum Velocity Evolution')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('viscosity_evolution.png', dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"\nSummary:")
        print(f"Average viscosity: {np.mean(viscosities):.6f}")
        print(f"Final viscosity: {viscosities[-1]:.6f}")
        print(f"Final max velocity: {max_velocities[-1]:.6f}")
        
        print(f"\nPlots generated:")
        print(f"- velocity_profile_plot_*.png (individual profiles)")
        print(f"- viscosity_evolution.png (time evolution)")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "--help" or sys.argv[1] == "-h":
            print(__doc__)
            print("\nUsage:")
            print("  python3 analyze.py        # Analyze all velocity profiles")
            print("  python3 analyze.py --help # Show this help")
        else:
            print("Unknown argument. Use --help for usage information.")
    else:
        analyze_all_profiles()