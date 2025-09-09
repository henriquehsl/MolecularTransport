
//Authored by Henrique Santos Lima
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <iomanip>



// --- Struct Vector2D com operadores para simplificar as operações matemáticas ---
struct Vector2D {
    double x, y;

    Vector2D() : x(0.0), y(0.0) {}
    Vector2D(double x_, double y_) : x(x_), y(y_) {}

    Vector2D operator+(const Vector2D& other) const {
        return Vector2D(x + other.x, y + other.y);
    }
    Vector2D operator-(const Vector2D& other) const {
        return Vector2D(x - other.x, y - other.y);
    }
    Vector2D operator*(double scalar) const {
        return Vector2D(x * scalar, y * scalar);
    }
    Vector2D operator/(double scalar) const {
        return Vector2D(x / scalar, y / scalar);
    }
    Vector2D& operator+=(const Vector2D& other) {
        x += other.x; y += other.y;
        return *this;
    }
    Vector2D& operator-=(const Vector2D& other) {
        x -= other.x; y -= other.y;
        return *this;
    }
    Vector2D& operator*=(double scalar) {
        x *= scalar; y *= scalar;
        return *this;
    }
    Vector2D& operator/=(double scalar) {
        x /= scalar; y /= scalar;
        return *this;
    }
    double magnitude() const {
        return std::sqrt(x*x + y*y);
    }
    double magnitude2() const {
        return x*x + y*y;
    }
};

// --- Estrutura de soma de Kahan para maior precisão em ponto flutuante ---
struct KahanSum {
    double sum;
    double c;
    KahanSum() : sum(0.0), c(0.0) {}
    void add(double value) {
        double y = value - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    double value() const { return sum; }
};

class PoiseuilleFlowMD {
private:
    // Parâmetros do sistema
    int Lx, Ly;
    double density, temperature, dt, force_x, wall_strength;
    double sigma, epsilon, rc, rc2;
    int wall_thickness, channel_height;
    int n_particles, n_wall_particles;

    // Dados das partículas
    std::vector<Vector2D> positions;
    std::vector<Vector2D> velocities;
    std::vector<Vector2D> forces;
    std::vector<Vector2D> wall_positions;

    // Gerador de números aleatórios
    std::mt19937 rng;
    std::normal_distribution<double> normal_dist;
    std::uniform_real_distribution<double> uniform_dist;

    // Dados de análise
    std::vector<double> kinetic_energy_history;
    std::vector<double> flow_rate_history;
    std::vector<std::vector<double>> velocity_profiles;

    int step_counter;

public:
    PoiseuilleFlowMD(int Lx_ = 20, int Ly_ = 20, double density_ = 0.6,
                     double temperature_ = 1.0, double dt_ = 0.003,
                     double force_x_ = 0.02, double wall_strength_ = 1.0)
        : Lx(Lx_), Ly(Ly_), density(density_), temperature(temperature_),
          dt(dt_), force_x(force_x_), wall_strength(wall_strength_),
          sigma(1.0), epsilon(1.0), rc(2.5), wall_thickness(2),
          rng(std::chrono::steady_clock::now().time_since_epoch().count()),
          normal_dist(0.0, 1.0), uniform_dist(-0.05, 0.05), step_counter(0) {

        rc2 = rc * rc;
        channel_height = Ly - 2 * wall_thickness;

        setupSystem();
        setupWalls();

        std::cout << "Initialized " << n_particles << " fluid particles and "
                  << n_wall_particles << " wall particles\n";
        std::cout << "Channel dimensions: " << Lx << " x " << channel_height << std::endl;
    }

    // --- FUNÇÃO CORRIGIDA PARA CÁLCULO DA VISCOSIDADE ---
    double calculateViscosity() const {
        // Assume que o perfil de velocidade é parabólico e o fluxo está em regime estacionário.
        // A viscosidade pode ser calculada a partir da relação de Poiseuille
        // Fx = 2 * eta * A * vmax / (canal_height)^2 -> eta = Fx * (canal_height)^2 / (2 * A * vmax)
        // No nosso caso a fórmula simplificada é eta = (Fx*L_x) / (4*vmax)
        
        // Obter o perfil de velocidade médio do estado estacionário
        if (velocity_profiles.empty()) {
            return 0.0;
        }
        int n_average = std::max(1, static_cast<int>(velocity_profiles.size()) / 2);
        int n_bins = static_cast<int>(velocity_profiles[0].size());
        std::vector<KahanSum> avg_velocities_kahan(n_bins);
        for (int i = static_cast<int>(velocity_profiles.size()) - n_average; 
             i < static_cast<int>(velocity_profiles.size()); ++i) {
            for (int j = 0; j < n_bins; ++j) {
                avg_velocities_kahan[j].add(velocity_profiles[i][j]);
            }
        }
        std::vector<double> avg_velocities(n_bins, 0.0);
        for (int j = 0; j < n_bins; ++j) {
            avg_velocities[j] = avg_velocities_kahan[j].value() / static_cast<double>(n_average);
        }
        
        // Encontrar a velocidade máxima no perfil
        double v_max = *std::max_element(avg_velocities.begin(), avg_velocities.end());

        // Calcular a viscosidade (eta)
        // Equação de Poiseuille para fluxo 2D: Fx = (2 * eta * v_max * Lx) / (canal_height * L_x) = (2 * eta * v_max) / canal_height
        // Fx * N_particulas = gradiente de pressão * área = -dP/dx * canal_height * Lx
        // A força por partícula é Fx. A força total é Fx * N.
        // F total = -dP/dx * Volume. Aqui, dP/dx * Volume_2D = F_total
        // A força externa é o gradiente de pressão: Fx_total = Fx * N = -(dP/dx) * Area_2D
        // No limite hidrodinâmico, o perfil de velocidade é v(y) = (Fy * Ly^2 / 2 * eta * Lx) * (1- (y/(Ly/2))^2)
        // Da equação de Navier-Stokes para fluxo 2D estacionário:
        // dP/dx = eta * d^2vx/dy^2
        // Assumindo um perfil parabólico, vx(y) = v_max * [1 - ((y - h/2)/(h/2))^2]
        // d^2vx/dy^2 = -2*v_max/(h/2)^2
        // dP/dx = eta * (-2*v_max)/(h/2)^2
        // O gradiente de pressão -dP/dx é o que a nossa força Fx representa.
        // -dP/dx = n_particles_per_area * Fx = (n_particles/Lx*Ly) * Fx
        // Fx = -dP/dx
        // eta = -dP/dx * (h/2)^2 / (2*v_max)
        // Com h = channel_height, e -dP/dx = force_x
        // eta = force_x * (channel_height/2)^2 / (2 * v_max)
        // Essa é a formula correta a ser usada para um fluido newtoniano.
        // A implementação anterior era do tensor de pressão, que é muito mais complexa e propensa a erros.

        if (v_max < 1e-8) { // Evitar divisão por zero se não houver fluxo
            return 0.0;
        }

        double h = static_cast<double>(channel_height);
        double viscosity = (force_x * (h*h) / (8.0 * v_max));
        
        return viscosity;
    }

private:
    void setupSystem() {
        double channel_volume = static_cast<double>(Lx) * static_cast<double>(channel_height);
        n_particles = static_cast<int>(std::round(density * channel_volume));
        positions.resize(n_particles);
        velocities.resize(n_particles);
        forces.resize(n_particles);
        int nx = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(n_particles))));
        int ny = static_cast<int>(std::ceil(static_cast<double>(n_particles) / nx));
        double dx = static_cast<double>(Lx) / static_cast<double>(nx);
        double dy = static_cast<double>(channel_height) / static_cast<double>(ny);
        int idx = 0;
        for (int i = 0; i < nx && idx < n_particles; ++i) {
            for (int j = 0; j < ny && idx < n_particles; ++j) {
                double x = (i + 0.5) * dx + uniform_dist(rng);
                double y = wall_thickness + (j + 0.5) * dy + uniform_dist(rng);
                x = std::max(0.1, std::min(x, static_cast<double>(Lx) - 0.1));
                y = std::max(static_cast<double>(wall_thickness) + 0.1, 
                           std::min(y, static_cast<double>(Ly - wall_thickness) - 0.1));
                positions[idx].x = x;
                positions[idx].y = y;
                ++idx;
            }
        }
        double sqrt_temp = std::sqrt(temperature);
        for (int i = 0; i < n_particles; ++i) {
            velocities[i].x = normal_dist(rng) * sqrt_temp;
            velocities[i].y = normal_dist(rng) * sqrt_temp;
        }
        Vector2D com_velocity(0.0, 0.0);
        for (const auto& vel : velocities) {
            com_velocity += vel;
        }
        com_velocity /= static_cast<double>(n_particles);
        for (auto& vel : velocities) {
            vel -= com_velocity;
        }
    }

    void setupWalls() {
        wall_positions.clear();
        for (int i = 0; i < Lx; ++i) {
            for (int j = 0; j < wall_thickness; ++j) {
                wall_positions.push_back(Vector2D(i + 0.5, j + 0.5));
                wall_positions.push_back(Vector2D(i + 0.5, Ly - wall_thickness + j + 0.5));
            }
        }
        n_wall_particles = static_cast<int>(wall_positions.size());
    }

    void applyPeriodicBoundary() {
        for (auto& pos : positions) {
            pos.x = pos.x - std::floor(pos.x / static_cast<double>(Lx)) * static_cast<double>(Lx);
        }
    }

    Vector2D lennardJonesForce(const Vector2D& r_vec, double r_mag) const {
        if (r_mag >= rc || r_mag < 1e-8) {
            return Vector2D(0.0, 0.0);
        }
        double r2 = r_mag * r_mag;
        double inv_r2 = 1.0 / r2;
        double inv_r6 = (sigma * sigma * inv_r2) * (sigma * sigma * inv_r2) * (sigma * sigma * inv_r2);
        double inv_r12 = inv_r6 * inv_r6;
        double force_magnitude = 48.0 * epsilon * (inv_r12 - 0.5 * inv_r6) * inv_r2;
        return r_vec * force_magnitude;
    }

    void calculateForces() {
        std::fill(forces.begin(), forces.end(), Vector2D(0.0, 0.0));
        for (int i = 0; i < n_particles; ++i) {
            for (int j = i + 1; j < n_particles; ++j) {
                Vector2D r_vec = positions[i] - positions[j];
                r_vec.x -= static_cast<double>(Lx) * std::round(r_vec.x / static_cast<double>(Lx));
                double r_mag = r_vec.magnitude();
                if (r_mag < rc) {
                    Vector2D force = lennardJonesForce(r_vec, r_mag);
                    forces[i] += force;
                    forces[j] -= force;
                }
            }
        }
        for (int i = 0; i < n_particles; ++i) {
            for (const auto& wall_pos : wall_positions) {
                Vector2D r_vec = positions[i] - wall_pos;
                r_vec.x -= static_cast<double>(Lx) * std::round(r_vec.x / static_cast<double>(Lx));
                double r_mag = r_vec.magnitude();
                if (r_mag < rc) {
                    Vector2D force = lennardJonesForce(r_vec, r_mag);
                    forces[i] += force * wall_strength;
                }
            }
        }
        for (auto& force : forces) {
            force.x += force_x;
        }
    }

    void velocityVerletStep() {
        const double half_dt2 = 0.5 * dt * dt;
        for (int i = 0; i < n_particles; ++i) {
            positions[i] += velocities[i] * dt + forces[i] * half_dt2;
        }
        applyPeriodicBoundary();
        std::vector<Vector2D> old_forces = forces;
        calculateForces();
        const double half_dt = 0.5 * dt;
        for (int i = 0; i < n_particles; ++i) {
            velocities[i] += (old_forces[i] + forces[i]) * half_dt;
        }
    }

    void applyThermostat(double coupling_strength = 1.0) {
        double current_temp = 0.0;
        for (const auto& vel : velocities) {
            current_temp += vel.magnitude2();
        }
        current_temp /= (2.0 * static_cast<double>(n_particles));
        if (current_temp > 1e-10) {
            double scaling_factor = std::sqrt(1.0 + coupling_strength *
                                            (temperature / current_temp - 1.0));
            for (auto& vel : velocities) {
                vel *= scaling_factor;
            }
        }
    }

    void applyStreamingPreservingThermostat(double coupling_strength = 0.005) {
        Vector2D mean_velocity(0.0, 0.0);
        for (const auto& vel : velocities) {
            mean_velocity += vel;
        }
        mean_velocity /= static_cast<double>(n_particles);
        double peculiar_temp = 0.0;
        for (const auto& vel : velocities) {
            Vector2D peculiar_vel = vel - mean_velocity;
            peculiar_temp += peculiar_vel.magnitude2();
        }
        peculiar_temp /= (2.0 * static_cast<double>(n_particles));
        if (peculiar_temp > 1e-10) {
            double scaling_factor = std::sqrt(1.0 + coupling_strength *
                                            (temperature / peculiar_temp - 1.0));
            for (auto& vel : velocities) {
                Vector2D peculiar_vel = vel - mean_velocity;
                vel = mean_velocity + peculiar_vel * scaling_factor;
            }
        }
    }

    void removeCenterOfMassTransverse() {
        double vcom_y = 0.0;
        for (const auto& vel : velocities) {
            vcom_y += vel.y;
        }
        vcom_y /= static_cast<double>(n_particles);
        for (auto& vel : velocities) {
            vel.y -= vcom_y;
        }
    }

    void removeCenterOfMassFull() {
        Vector2D com_velocity(0.0, 0.0);
        for (const auto& vel : velocities) {
            com_velocity += vel;
        }
        com_velocity /= static_cast<double>(n_particles);
        for (auto& vel : velocities) {
            vel -= com_velocity;
        }
    }

    std::vector<double> calculateVelocityProfile(int n_bins = 100) const {
        std::vector<KahanSum> bin_velocities(n_bins);
        std::vector<int> bin_counts(n_bins, 0);
        double bin_width = static_cast<double>(channel_height) / static_cast<double>(n_bins);
        double y_min = static_cast<double>(wall_thickness);
        for (int i = 0; i < n_particles; ++i) {
            double y_pos = positions[i].y;
            int bin_idx = static_cast<int>((y_pos - y_min) / bin_width);
            if (bin_idx >= 0 && bin_idx < n_bins) {
                bin_velocities[bin_idx].add(velocities[i].x);
                bin_counts[bin_idx]++;
            }
        }
        std::vector<double> avg_velocities(n_bins, 0.0);
        for (int i = 0; i < n_bins; ++i) {
            if (bin_counts[i] > 0) {
                avg_velocities[i] = bin_velocities[i].value() / static_cast<double>(bin_counts[i]);
            }
        }
        return avg_velocities;
    }

    double calculateKineticEnergy() const {
        double ke = 0.0;
        for (const auto& vel : velocities) {
            ke += 0.5 * vel.magnitude2();
        }
        return ke;
    }

    double calculateFlowRate() const {
        double flow_rate = 0.0;
        for (const auto& vel : velocities) {
            flow_rate += vel.x;
        }
        return flow_rate / static_cast<double>(n_particles);
    }

    double calculateTemperature() const {
        double temp = 0.0;
        for (const auto& vel : velocities) {
            temp += vel.magnitude2();
        }
        return temp / (2.0 * static_cast<double>(n_particles));
    }

    void fitQuadratic(const std::vector<double>& x, const std::vector<double>& y, 
                      double& a, double& b, double& c) const {
        size_t n = x.size();
        double sum_x = 0.0, sum_x2 = 0.0, sum_x3 = 0.0, sum_x4 = 0.0;
        double sum_y = 0.0, sum_xy = 0.0, sum_x2y = 0.0;
        for (size_t i = 0; i < n; ++i) {
            double xi = x[i], yi = y[i];
            sum_x += xi;
            sum_x2 += xi * xi;
            sum_x3 += xi * xi * xi;
            sum_x4 += xi * xi * xi * xi;
            sum_y += yi;
            sum_xy += xi * yi;
            sum_x2y += xi * xi * yi;
        }
        double denom = n * sum_x2 * sum_x4 + 2 * sum_x * sum_x2 * sum_x3 - 
                      sum_x2 * sum_x2 * sum_x2 - n * sum_x3 * sum_x3 - sum_x * sum_x * sum_x4;
        a = (n * sum_x2 * sum_x2y + sum_x * sum_x3 * sum_y + sum_x * sum_x * sum_x2y - 
             sum_x2 * sum_x * sum_x2y - n * sum_x3 * sum_xy - sum_x2 * sum_x * sum_y) / denom;
        b = (n * sum_x3 * sum_x2y + sum_x * sum_x4 * sum_y + sum_x2 * sum_x2 * sum_xy - 
             sum_x3 * sum_x * sum_x2y - n * sum_x4 * sum_xy - sum_x2 * sum_x3 * sum_y) / denom;
        c = (sum_x2 * sum_x3 * sum_x2y + sum_x * sum_x2 * sum_xy + sum_x * sum_x3 * sum_x2y - 
             sum_x3 * sum_x * sum_xy - sum_x2 * sum_x4 * sum_xy - sum_x2 * sum_x2 * sum_x2y) / denom;
    }

public:
    void runSimulation(int n_equilibration = 3000, int n_production = 5000,
                      int thermostat_steps = 1000, int save_frequency = 50) {
        auto start_time = std::chrono::high_resolution_clock::now();
        std::cout << "Starting equilibration phase..." << std::endl;
        double original_force = force_x;
        force_x = 0.0;
        for (int step = 0; step < n_equilibration; ++step) {
            velocityVerletStep();
            if (step < thermostat_steps) {
                applyThermostat(1.00);
            }
            if (step % 5000 == 0) {
                removeCenterOfMassFull();
                double temp = calculateTemperature();
                std::cout << "Equilibration step " << step << ": T="
                         << std::fixed << std::setprecision(4) << temp << std::endl;
            }
        }
        std::cout << "Starting production phase with driving force..." << std::endl;
        force_x = original_force;
        kinetic_energy_history.clear();
        flow_rate_history.clear();
        velocity_profiles.clear();
        for (int step = 0; step < n_production; ++step) {
            velocityVerletStep();
            if (step % 5000 == 0 && step > 0) {
                applyStreamingPreservingThermostat(1.0);
            }
            double current_temp = calculateTemperature();
            if (current_temp > 1.5 * temperature || current_temp < 0.5 * temperature) {
                applyStreamingPreservingThermostat(0.02);
            }
            if (step % 10000 == 0 && step > 0) {
                removeCenterOfMassTransverse();
            }
            step_counter = step;
            if (step % save_frequency == 0) {
                double ke = calculateKineticEnergy();
                double flow_rate = calculateFlowRate();
                kinetic_energy_history.push_back(ke);
                flow_rate_history.push_back(flow_rate);
                std::vector<double> v_profile = calculateVelocityProfile(40);
                velocity_profiles.push_back(v_profile);
            }
            if (step % 5000 == 0) {
                double flow_rate = calculateFlowRate();
                double temp = calculateTemperature();
                std::cout << "Production step " << step << ": Flow rate="
                         << std::fixed << std::setprecision(4) << flow_rate
                         << ", T=" << std::setprecision(4) << temp << std::endl;
            }
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Simulation completed in " << duration.count() / 1000.0
                  << " seconds!" << std::endl;
    }

    void analyzeResults() const {
        if (velocity_profiles.empty()) {
            std::cout << "No data to analyze. Run simulation first." << std::endl;
            return;
        }
        int n_average = std::max(1, static_cast<int>(velocity_profiles.size()) / 2);
        int n_bins = static_cast<int>(velocity_profiles[0].size());
        std::vector<KahanSum> avg_velocities_kahan(n_bins);
        for (int i = static_cast<int>(velocity_profiles.size()) - n_average; 
             i < static_cast<int>(velocity_profiles.size()); ++i) {
            for (int j = 0; j < n_bins; ++j) {
                avg_velocities_kahan[j].add(velocity_profiles[i][j]);
            }
        }
        std::vector<double> avg_velocities(n_bins, 0.0);
        for (int j = 0; j < n_bins; ++j) {
            avg_velocities[j] = avg_velocities_kahan[j].value() / static_cast<double>(n_average);
        }
        double v_max = *std::max_element(avg_velocities.begin(), avg_velocities.end());
        KahanSum avg_flow_rate_kahan;
        int n_flow_avg = std::max(1, static_cast<int>(flow_rate_history.size()) / 2);
        for (int i = static_cast<int>(flow_rate_history.size()) - n_flow_avg; 
             i < static_cast<int>(flow_rate_history.size()); ++i) {
            avg_flow_rate_kahan.add(flow_rate_history[i]);
        }
        double avg_flow_rate = avg_flow_rate_kahan.value() / static_cast<double>(n_flow_avg);
        std::cout << "\n=== SIMULATION RESULTS ===" << std::endl;
        std::cout << "Maximum velocity: " << std::fixed << std::setprecision(4) << v_max << std::endl;
        std::cout << "Average flow rate: " << std::fixed << std::setprecision(6) << avg_flow_rate << std::endl;
        std::cout << "Final temperature: " << std::setprecision(4) << calculateTemperature() << std::endl;
    }

    void saveResults(const std::string& filename = "velocity_profile.dat") const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error opening file for writing!" << std::endl;
            return;
        }
        if (!velocity_profiles.empty()) {
            int n_average = std::max(1, static_cast<int>(velocity_profiles.size()) / 2);
            int n_bins = static_cast<int>(velocity_profiles[0].size());
            std::vector<KahanSum> avg_velocities_kahan(n_bins);
            for (int i = static_cast<int>(velocity_profiles.size()) - n_average; 
                 i < static_cast<int>(velocity_profiles.size()); ++i) {
                for (int j = 0; j < n_bins; ++j) {
                    avg_velocities_kahan[j].add(velocity_profiles[i][j]);
                }
            }
            std::vector<double> avg_velocities(n_bins, 0.0);
            for (int j = 0; j < n_bins; ++j) {
                avg_velocities[j] = avg_velocities_kahan[j].value() / static_cast<double>(n_average);
            }
            file << "# y_position velocity" << std::endl;
            double bin_width = static_cast<double>(channel_height) / static_cast<double>(n_bins);
            for (int i = 0; i < n_bins; ++i) {
                double y_pos = static_cast<double>(wall_thickness) + (i + 0.5) * bin_width;
                file << y_pos << " " << avg_velocities[i] << std::endl;
            }
        }
        file.close();
        std::cout << "Results saved to " << filename << std::endl;
    }
};

int main() {
    std::cout << "Poiseuille Flow MD Simulation (C++)" << std::endl;
    std::cout << "====================================" << std::endl;

    PoiseuilleFlowMD sim(40, 40, 0.8, 1.2, 3e-3, 0.005, 1.0);

    // Tempos de simulação aumentados para garantir melhor estabilização.
    sim.runSimulation(100000, 1000000,20000, 100);

    sim.analyzeResults();

    // A chamada para calcular a viscosidade agora usa a fórmula corrigida
    double viscosity = sim.calculateViscosity();
    std::cout << "Viscosity at T=1.2: " << std::fixed << std::setprecision(6) << viscosity << std::endl;

    sim.saveResults("velocity_profile.dat");

    std::cout << "\nSimulation complete! Check velocity_profile.dat for results." << std::endl;
    std::cout << "You can plot the results using gnuplot or any plotting software." << std::endl;

    return 0;
}
