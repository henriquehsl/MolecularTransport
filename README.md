# MolecularTransport

Molecular dynamics simulations of viscosity in Poiseuille flow.

## Overview

This repository contains C++ code for simulating viscosity in Poiseuille flow using molecular dynamics. The simulation models fluid particles flowing between parallel walls under the influence of an external driving force, allowing for the measurement of viscosity through velocity profile analysis.

## Features

- **Molecular Dynamics Engine**: Complete MD simulation with Lennard-Jones interactions
- **Poiseuille Flow Setup**: Parallel wall geometry with reflecting boundary conditions
- **Velocity Profile Analysis**: Automatic calculation of velocity profiles across the channel
- **Viscosity Measurement**: Estimation of fluid viscosity from velocity gradients
- **Output Formats**: XYZ files for particle positions and data files for velocity profiles

## Building and Running

### Prerequisites
- C++ compiler with C++11 support (g++, clang++)
- Make utility
- Python 3 (optional, for analysis tools)
- numpy and matplotlib (optional, for plotting)

### Compilation
```bash
make
```

### Running the simulation
```bash
make run
```

Or run directly:
```bash
./poiseuille_flow
```

### Other useful commands
```bash
make debug     # Compile with debug flags
make clean     # Remove build artifacts and output files
make analyze   # Analyze velocity profiles (requires Python, numpy, matplotlib)
make help      # Show available targets
```

**Note**: Analysis features require Python dependencies:
```bash
pip install numpy matplotlib
```

## Simulation Parameters

The simulation uses the following default parameters (see `parameters.txt` for details):
- 500 particles in a 10×10×20 simulation box
- Lennard-Jones interactions with σ=1.0, ε=1.0
- External driving force of 0.1 in the flow direction
- 10,000 simulation steps with output every 100 steps

## Output Files

The simulation generates:
- `output_*.xyz`: Particle positions in XYZ format for visualization
- `velocity_profile_*.dat`: Velocity profiles for viscosity analysis

## Theory

Poiseuille flow describes laminar flow of a viscous fluid through a channel. In this simulation:
- Particles interact via Lennard-Jones potential
- Walls provide boundary conditions in the y-direction
- An external force drives flow in the z-direction
- Viscosity is estimated from the velocity gradient: η ∝ F/∂v/∂y

## Visualization

The XYZ output files can be visualized using molecular visualization software such as:
- VMD (Visual Molecular Dynamics)
- OVITO
- PyMOL

## Units

The simulation uses reduced units where:
- Length: σ (Lennard-Jones sigma parameter)
- Energy: ε (Lennard-Jones epsilon parameter)  
- Mass: particle mass
- Time: σ√(m/ε)
