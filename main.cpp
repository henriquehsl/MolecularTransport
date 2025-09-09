#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <random>
#include <iomanip>

// Constants
const double PI = 3.14159265359;
const double KB = 1.0; // Boltzmann constant (reduced units)

// Particle structure
struct Particle {
    double x, y, z;    // position
    double vx, vy, vz; // velocity
    double fx, fy, fz; // force
    
    Particle(double x = 0, double y = 0, double z = 0, 
             double vx = 0, double vy = 0, double vz = 0) 
        : x(x), y(y), z(z), vx(vx), vy(vy), vz(vz), fx(0), fy(0), fz(0) {}
};

// Simulation parameters
struct SimParams {
    int nparticles = 500;
    double box_x = 10.0;
    double box_y = 10.0;
    double box_z = 20.0;
    double wall_force = 0.1;  // External force driving flow
    double temperature = 1.0;
    double density = 0.8;
    double dt = 0.005;
    int nsteps = 10000;
    int output_freq = 100;
    double sigma = 1.0;       // Lennard-Jones sigma
    double epsilon = 1.0;     // Lennard-Jones epsilon
    double cutoff = 2.5;      // LJ cutoff distance
};

class PoiseuilleFlowMD {
private:
    std::vector<Particle> particles;
    SimParams params;
    std::mt19937 rng;
    
    // Performance tracking
    std::vector<double> velocity_profile;
    std::vector<double> stress_tensor;
    int profile_bins = 50;
    
public:
    PoiseuilleFlowMD(const SimParams& p) : params(p), rng(std::random_device{}()) {
        velocity_profile.resize(profile_bins, 0.0);
        stress_tensor.resize(profile_bins, 0.0);
    }
    
    void initialize() {
        particles.clear();
        particles.reserve(params.nparticles);
        
        // Initialize particles in a simple cubic lattice
        int nx = static_cast<int>(std::cbrt(params.nparticles)) + 1;
        double spacing = std::min(std::min(params.box_x, params.box_y), params.box_z) / nx;
        
        std::uniform_real_distribution<double> vel_dist(-0.5, 0.5);
        
        for (int i = 0; i < params.nparticles; ++i) {
            double x = (i % nx) * spacing + 0.5;
            double y = ((i / nx) % nx) * spacing + 0.5;
            double z = (i / (nx * nx)) * spacing + 0.5;
            
            // Ensure particles are within walls (z-direction is flow direction)
            if (y >= params.box_y - 1.0) y = params.box_y - 1.1;
            if (y <= 1.0) y = 1.1;
            
            // Random initial velocities
            double vx = vel_dist(rng) * std::sqrt(params.temperature);
            double vy = vel_dist(rng) * std::sqrt(params.temperature);
            double vz = vel_dist(rng) * std::sqrt(params.temperature);
            
            particles.emplace_back(x, y, z, vx, vy, vz);
        }
        
        std::cout << "Initialized " << particles.size() << " particles" << std::endl;
    }
    
    void compute_forces() {
        // Reset forces
        for (auto& p : particles) {
            p.fx = p.fy = p.fz = 0.0;
        }
        
        // Lennard-Jones forces between particles
        for (size_t i = 0; i < particles.size(); ++i) {
            for (size_t j = i + 1; j < particles.size(); ++j) {
                double dx = particles[i].x - particles[j].x;
                double dy = particles[i].y - particles[j].y;
                double dz = particles[i].z - particles[j].z;
                
                // Periodic boundary conditions in x and z
                if (dx > params.box_x / 2) dx -= params.box_x;
                if (dx < -params.box_x / 2) dx += params.box_x;
                if (dz > params.box_z / 2) dz -= params.box_z;
                if (dz < -params.box_z / 2) dz += params.box_z;
                
                double r2 = dx*dx + dy*dy + dz*dz;
                
                if (r2 < params.cutoff * params.cutoff && r2 > 0.01) {
                    double r6 = std::pow(params.sigma * params.sigma / r2, 3);
                    double r12 = r6 * r6;
                    double force_mag = 24 * params.epsilon / r2 * (2 * r12 - r6);
                    
                    double fx = force_mag * dx;
                    double fy = force_mag * dy;
                    double fz = force_mag * dz;
                    
                    particles[i].fx += fx;
                    particles[i].fy += fy;
                    particles[i].fz += fz;
                    particles[j].fx -= fx;
                    particles[j].fy -= fy;
                    particles[j].fz -= fz;
                }
            }
        }
        
        // Wall forces and external driving force
        for (auto& p : particles) {
            // Wall repulsion (y-direction walls)
            if (p.y < 1.0) {
                double wall_dist = 1.0 - p.y;
                p.fy += 100.0 * wall_dist * wall_dist;
            }
            if (p.y > params.box_y - 1.0) {
                double wall_dist = p.y - (params.box_y - 1.0);
                p.fy -= 100.0 * wall_dist * wall_dist;
            }
            
            // External force driving flow in z-direction
            p.fz += params.wall_force;
        }
    }
    
    void integrate() {
        // Velocity Verlet integration
        for (auto& p : particles) {
            // Update positions
            p.x += p.vx * params.dt + 0.5 * p.fx * params.dt * params.dt;
            p.y += p.vy * params.dt + 0.5 * p.fy * params.dt * params.dt;
            p.z += p.vz * params.dt + 0.5 * p.fz * params.dt * params.dt;
            
            // Apply periodic boundary conditions
            if (p.x < 0) p.x += params.box_x;
            if (p.x >= params.box_x) p.x -= params.box_x;
            if (p.z < 0) p.z += params.box_z;
            if (p.z >= params.box_z) p.z -= params.box_z;
            
            // Reflect at walls in y-direction
            if (p.y < 0.5) {
                p.y = 0.5;
                p.vy = -p.vy;
            }
            if (p.y > params.box_y - 0.5) {
                p.y = params.box_y - 0.5;
                p.vy = -p.vy;
            }
        }
        
        // Store old forces
        std::vector<double> old_fx(particles.size()), old_fy(particles.size()), old_fz(particles.size());
        for (size_t i = 0; i < particles.size(); ++i) {
            old_fx[i] = particles[i].fx;
            old_fy[i] = particles[i].fy;
            old_fz[i] = particles[i].fz;
        }
        
        // Compute new forces
        compute_forces();
        
        // Update velocities
        for (size_t i = 0; i < particles.size(); ++i) {
            particles[i].vx += 0.5 * (old_fx[i] + particles[i].fx) * params.dt;
            particles[i].vy += 0.5 * (old_fy[i] + particles[i].fy) * params.dt;
            particles[i].vz += 0.5 * (old_fz[i] + particles[i].fz) * params.dt;
        }
    }
    
    void analyze_flow() {
        // Reset analysis arrays
        std::fill(velocity_profile.begin(), velocity_profile.end(), 0.0);
        std::fill(stress_tensor.begin(), stress_tensor.end(), 0.0);
        std::vector<int> bin_count(profile_bins, 0);
        
        double bin_width = params.box_y / profile_bins;
        
        for (const auto& p : particles) {
            int bin = static_cast<int>(p.y / bin_width);
            if (bin >= 0 && bin < profile_bins) {
                velocity_profile[bin] += p.vz;
                bin_count[bin]++;
            }
        }
        
        // Average velocities in each bin
        for (int i = 0; i < profile_bins; ++i) {
            if (bin_count[i] > 0) {
                velocity_profile[i] /= bin_count[i];
            }
        }
    }
    
    double calculate_viscosity() {
        analyze_flow();
        
        // Calculate velocity gradient (dv/dy) for viscosity estimation
        double gradient_sum = 0.0;
        int gradient_count = 0;
        double bin_width = params.box_y / profile_bins;
        
        for (int i = 1; i < profile_bins - 1; ++i) {
            double gradient = (velocity_profile[i+1] - velocity_profile[i-1]) / (2 * bin_width);
            gradient_sum += std::abs(gradient);
            gradient_count++;
        }
        
        double avg_gradient = gradient_count > 0 ? gradient_sum / gradient_count : 0.0;
        
        // Viscosity = applied_force / velocity_gradient (simplified)
        double viscosity = avg_gradient > 1e-6 ? params.wall_force / avg_gradient : 0.0;
        
        return viscosity;
    }
    
    void output_state(int step) {
        std::string filename = "output_" + std::to_string(step) + ".xyz";
        std::ofstream file(filename);
        
        file << particles.size() << std::endl;
        file << "Step " << step << std::endl;
        
        for (const auto& p : particles) {
            file << "Ar " << std::fixed << std::setprecision(6) 
                 << p.x << " " << p.y << " " << p.z << std::endl;
        }
        
        file.close();
    }
    
    void output_velocity_profile(int step) {
        std::string filename = "velocity_profile_" + std::to_string(step) + ".dat";
        std::ofstream file(filename);
        
        double bin_width = params.box_y / profile_bins;
        
        file << "# y_position velocity_z" << std::endl;
        for (int i = 0; i < profile_bins; ++i) {
            double y_pos = (i + 0.5) * bin_width;
            file << std::fixed << std::setprecision(6) 
                 << y_pos << " " << velocity_profile[i] << std::endl;
        }
        
        file.close();
    }
    
    void run() {
        std::cout << "Starting Poiseuille flow simulation..." << std::endl;
        std::cout << "Parameters:" << std::endl;
        std::cout << "  Particles: " << params.nparticles << std::endl;
        std::cout << "  Box size: " << params.box_x << " x " << params.box_y << " x " << params.box_z << std::endl;
        std::cout << "  Steps: " << params.nsteps << std::endl;
        std::cout << "  dt: " << params.dt << std::endl;
        
        initialize();
        compute_forces();
        
        for (int step = 0; step < params.nsteps; ++step) {
            integrate();
            
            if (step % params.output_freq == 0) {
                double viscosity = calculate_viscosity();
                std::cout << "Step " << step << ", Viscosity: " << std::fixed << std::setprecision(6) 
                         << viscosity << std::endl;
                
                if (step % (params.output_freq * 10) == 0) {
                    output_state(step);
                    output_velocity_profile(step);
                }
            }
        }
        
        // Final analysis
        double final_viscosity = calculate_viscosity();
        std::cout << std::endl << "Final viscosity estimate: " << final_viscosity << std::endl;
        
        output_state(params.nsteps);
        output_velocity_profile(params.nsteps);
    }
};

int main() {
    SimParams params;
    
    std::cout << "Molecular Dynamics Simulation of Poiseuille Flow" << std::endl;
    std::cout << "================================================" << std::endl;
    
    PoiseuilleFlowMD simulation(params);
    simulation.run();
    
    std::cout << std::endl << "Simulation completed!" << std::endl;
    std::cout << "Output files: output_*.xyz (particle positions)" << std::endl;
    std::cout << "              velocity_profile_*.dat (velocity profiles)" << std::endl;
    
    return 0;
}