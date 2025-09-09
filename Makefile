# Makefile for Molecular Dynamics Poiseuille Flow Simulation

CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra
LDFLAGS = -lm
TARGET = poiseuille_flow
SOURCES = main.cpp

# Default target
all: $(TARGET)

# Build the main executable
$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(LDFLAGS)

# Debug build
debug: CXXFLAGS = -std=c++11 -g -Wall -Wextra -DDEBUG
debug: $(TARGET)

# Clean build artifacts
clean:
	rm -f $(TARGET)
	rm -f output_*.xyz
	rm -f velocity_profile_*.dat

# Run the simulation
run: $(TARGET)
	./$(TARGET)

# Run with different parameters (example)
run-small: $(TARGET)
	./$(TARGET)

# Analyze results (requires Python with matplotlib)
analyze:
	@if [ -f "velocity_profile_0.dat" ]; then \
		python3 analyze.py; \
	else \
		echo "No velocity profile files found. Run simulation first."; \
	fi

# Help target
help:
	@echo "Available targets:"
	@echo "  all       - Build the simulation (default)"
	@echo "  debug     - Build with debug flags"
	@echo "  clean     - Remove build artifacts and output files"
	@echo "  run       - Build and run the simulation"
	@echo "  analyze   - Analyze velocity profiles (requires Python)"
	@echo "  help      - Show this help message"

.PHONY: all debug clean run run-small analyze help