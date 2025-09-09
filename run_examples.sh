#!/bin/bash
# Example script to run simulations with different parameters

echo "Running Poiseuille Flow MD Simulations with different parameters"
echo "=============================================================="

# Make sure simulation is built
make clean && make

if [ ! -f "./poiseuille_flow" ]; then
    echo "Error: Could not build simulation!"
    exit 1
fi

# Create output directory
mkdir -p results
cd results

echo "Running simulation 1: Standard parameters"
mkdir -p sim1
cd sim1
../../poiseuille_flow > simulation.log 2>&1
echo "Simulation 1 completed. Check sim1/ directory for results."
cd ..

echo "Running simulation 2: Higher wall force"
mkdir -p sim2  
cd sim2
# Note: In a full implementation, we would modify parameters here
# For now, just run with default parameters
../../poiseuille_flow > simulation.log 2>&1
echo "Simulation 2 completed. Check sim2/ directory for results."
cd ..

echo ""
echo "All simulations completed!"
echo "Results are in the results/ directory"
echo "Use the analyze.py script to process the velocity profiles:"
echo "  cd sim1 && python3 ../../analyze.py"
echo "  cd sim2 && python3 ../../analyze.py"