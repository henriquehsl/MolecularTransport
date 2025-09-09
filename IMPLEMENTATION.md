# Molecular Dynamics Implementation Details

## Physics Background

### Poiseuille Flow
Poiseuille flow is the laminar flow of a viscous fluid through a pipe or between parallel plates. The velocity profile in steady state is parabolic for flow between parallel plates:

```
v(y) = (dP/dx) * (h²/2μ) * (1 - (2y/h)²)
```

Where:
- `dP/dx` is the pressure gradient (or applied force per unit volume)
- `h` is the channel height
- `μ` is the dynamic viscosity
- `y` is the distance from the channel center

### Viscosity Calculation
From the velocity profile, viscosity can be calculated using:

```
μ = F / (dv/dy)
```

Where `F` is the applied force and `dv/dy` is the velocity gradient.

## Implementation Details

### Molecular Dynamics Components

1. **Particle Structure**
   - Position (x, y, z)
   - Velocity (vx, vy, vz)  
   - Force (fx, fy, fz)

2. **Force Calculation**
   - Lennard-Jones potential: `U(r) = 4ε[(σ/r)¹² - (σ/r)⁶]`
   - Force: `F = -dU/dr`
   - Wall repulsion forces
   - External driving force

3. **Integration Algorithm**
   - Velocity Verlet integration for stability
   - Time step: dt = 0.005 (reduced units)

4. **Boundary Conditions**
   - Periodic boundaries in x and z directions
   - Reflecting walls in y direction
   - External force applied in z direction

### Key Parameters

- **System Size**: 10×10×20 (length×height×width)
- **Particle Count**: 500 (adjustable)
- **Temperature**: 1.0 (reduced units)
- **Time Step**: 0.005
- **Simulation Steps**: 10,000

### Output Files

1. **output_*.xyz**: Particle positions in XYZ format
   - Compatible with VMD, OVITO, PyMOL
   - Contains atomic coordinates for each time step

2. **velocity_profile_*.dat**: Velocity profile data
   - Column 1: y-position
   - Column 2: average velocity in z-direction
   - Used for viscosity calculation

### Analysis Tools

The `analyze.py` script provides:
- Viscosity calculation from velocity profiles
- Velocity profile visualization
- Time evolution plots
- Comparison with theoretical Poiseuille profile

### Reduced Units

The simulation uses reduced units where:
- Length unit: σ (Lennard-Jones sigma)
- Energy unit: ε (Lennard-Jones epsilon)
- Mass unit: m (particle mass)
- Time unit: σ√(m/ε)

### Typical Results

Expected viscosity values range from 0.05 to 0.2 in reduced units, depending on:
- Applied force strength
- Temperature
- Density
- System size

## Usage Examples

### Basic Simulation
```bash
make run
```

### Analysis
```bash
make analyze
```

### Custom Parameters
Modify the `SimParams` structure in `main.cpp` and recompile.

## Validation

The simulation can be validated by:
1. Checking parabolic velocity profile development
2. Comparing with theoretical Poiseuille flow
3. Verifying steady-state viscosity values
4. Testing different system sizes and forces

## Extensions

Possible extensions include:
- Parameter input from file
- Multiple fluid components
- Different wall materials
- Temperature gradients
- Non-Newtonian fluids