// simulation_corrected_timeavg.cpp
// Single-thread 2D LJ deformation test (NPT equilibration -> NVT + axial strain)
// Time-averaged temperature profile and output filenames with T, N, rho
// Authored by Henrique Santos Lima

#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <string>

using Vec2 = std::array<double,2>;
using Mat2 = std::array<double,4>;

constexpr int   N_side       = 20;
constexpr int   N            = N_side * N_side;
constexpr double rho         = 0.8;
constexpr double T0          = 0.8;
constexpr double dt          = 1e-7;
constexpr int   steps        = 1000000;
constexpr int   sample_every = 500;
constexpr int   deform_every = 2000;
constexpr double d_eps       = 1e-3;
constexpr int   equil_steps  = 50000;

constexpr double rc          = 2.5;
constexpr double rc2         = rc*rc;
constexpr double eps_LJ      = 1.0;
constexpr double sig_LJ      = 1.0;
constexpr double kB          = 1.0;
constexpr double m           = 1.0;
constexpr int    dim         = 2;

// Thermostat / barostat (Berendsen)
constexpr double tau_T       = 50.0 * dt;
constexpr double T_target    = T0;
constexpr double P_target    = 0.0;
constexpr double tau_P       = 500.0 * dt;
constexpr double beta        = 1.0;

// LJ precompute
constexpr double sig2 = sig_LJ * sig_LJ;
constexpr double sig6 = sig2 * sig2 * sig2;
constexpr double sig12 = sig6 * sig6;

// Shift LJ potential at cutoff for continuity in U
constexpr bool   use_shift = true;
constexpr double inv_rc2   = 1.0 / rc2;
constexpr double inv_rc6   = sig6 * inv_rc2 * inv_rc2 * inv_rc2;
constexpr double inv_rc12  = inv_rc6 * inv_rc6;
constexpr double U_shift   = 4.0 * eps_LJ * (inv_rc12 - inv_rc6); // so that U(rc)=0 if use_shift

constexpr double MIN_DR2   = 1e-12;

inline void wrap(double &x, double L) {
    x -= L * std::floor(x / L);
}

inline void minimum_image(double &dx, double &dy, double Lx, double Ly) {
    dx -= Lx * std::round(dx / Lx);
    dy -= Ly * std::round(dy / Ly);
}

struct ForceResults {
    Mat2 virial;    // (1/2) Σ_{pairs} r_ij ⊗ f_ij
    double U;       // potential energy
    double P_pair;  // Σ_{pairs} r_ij · f_ij
};

void compute_forces(
    const std::vector<Vec2> &r,
    double Lx, double Ly,
    std::vector<Vec2> &F,
    ForceResults &res)
{
    std::fill(F.begin(), F.end(), Vec2{0.0, 0.0});
    res.virial = {0.0,0.0,0.0,0.0};
    res.U = 0.0;
    res.P_pair = 0.0;

    for (int i = 0; i < N-1; ++i) {
        const double xi = r[i][0];
        const double yi = r[i][1];
        for (int j = i+1; j < N; ++j) {
            double dx = r[j][0] - xi;
            double dy = r[j][1] - yi;
            minimum_image(dx, dy, Lx, Ly);
            double dr2 = dx*dx + dy*dy;
            if (dr2 >= rc2 || dr2 < MIN_DR2) continue;

            double inv_r2  = 1.0 / dr2;
            double inv_r6  = sig6 * inv_r2 * inv_r2 * inv_r2;
            double inv_r12 = inv_r6 * inv_r6;

            double fac = 24.0 * eps_LJ * (2.0*inv_r12 - inv_r6) * inv_r2;

            double fx = dx * fac;
            double fy = dy * fac;
            F[i][0] += fx; F[i][1] += fy;
            F[j][0] -= fx; F[j][1] -= fy;

            double pot = 4.0 * eps_LJ * (inv_r12 - inv_r6);
            if (use_shift) pot -= U_shift;
            res.U += pot;

            res.virial[0] += 0.5 * dx * fx;
            res.virial[1] += 0.5 * dx * fy;
            res.virial[2] += 0.5 * dy * fx;
            res.virial[3] += 0.5 * dy * fy;

            res.P_pair += dx * fx + dy * fy;
        }
    }
}

double global_temperature(const std::vector<Vec2> &v) {
    double vx_sum = 0.0, vy_sum = 0.0;
    for (const auto &vi : v) {
        vx_sum += vi[0];
        vy_sum += vi[1];
    }
    vx_sum /= N;
    vy_sum /= N;

    double sum_sq = 0.0;
    for (const auto &vi : v) {
        double vx = vi[0] - vx_sum;
        double vy = vi[1] - vy_sum;
        sum_sq += vx*vx + vy*vy;
    }
    return (m * sum_sq) / (dim * N * kB);
}

Mat2 compute_cauchy_stress(const std::vector<Vec2> &v, const ForceResults &fr, double V) {
    double vx_sum = 0.0, vy_sum = 0.0;
    for (const auto &vi : v) {
        vx_sum += vi[0];
        vy_sum += vi[1];
    }
    vx_sum /= N;
    vy_sum /= N;

    Mat2 kin = {0.0,0.0,0.0,0.0};
    for (const auto &vi : v) {
        double vx = vi[0] - vx_sum;
        double vy = vi[1] - vy_sum;
        kin[0] += m * vx * vx;
        kin[1] += m * vx * vy;
        kin[3] += m * vy * vy;
    }
    kin[2] = kin[1];

    Mat2 sigma;
    sigma[0] = (kin[0] + fr.virial[0]) / V;
    sigma[1] = (kin[1] + fr.virial[1]) / V;
    sigma[2] = (kin[2] + fr.virial[2]) / V;
    sigma[3] = (kin[3] + fr.virial[3]) / V;
    return sigma;
}

std::vector<double> temperature_profile(
    const std::vector<Vec2> &v,
    const std::vector<Vec2> &r,
    double Ly,
    int nbins)
{
    std::vector<int> counts(nbins,0);
    std::vector<double> vx_sum(nbins,0.0), vy_sum(nbins,0.0), v2_sum(nbins,0.0);

    for (int i=0; i<N; ++i) {
        int idx = static_cast<int>(std::floor(r[i][1]/Ly * nbins));
        if (idx < 0) idx = 0;
        if (idx >= nbins) idx = nbins - 1;
        double vx = v[i][0];
        double vy = v[i][1];
        counts[idx]++;
        vx_sum[idx] += vx;
        vy_sum[idx] += vy;
        v2_sum[idx] += vx*vx + vy*vy;
    }

    std::vector<double> Tbins(nbins,0.0);
    for (int b=0; b<nbins; ++b) {
        if (counts[b] > 1) {
            double invc = 1.0 / counts[b];
            double vxm = vx_sum[b]*invc;
            double vym = vy_sum[b]*invc;
            double var = v2_sum[b]*invc - (vxm*vxm + vym*vym);
            if (var > 0)
                Tbins[b] = (m * var) / (dim * kB);
            else
                Tbins[b] = T0;
        } else {
            Tbins[b] = T0;
        }
    }
    return Tbins;
}

// Moving average smoothing for vector<double>
std::vector<double> moving_average(const std::vector<double>& data, int window) {
    std::vector<double> smoothed(data.size(), 0.0);
    int N = data.size();
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        int count = 0;
        for (int j = std::max(0, i-window); j <= std::min(N-1, i+window); ++j) {
            sum += data[j];
            count++;
        }
        smoothed[i] = sum / count;
    }
    return smoothed;
}

int main() {
    auto t0 = std::chrono::high_resolution_clock::now();

    double A0 = N / rho;
    double Ly = std::sqrt(A0);
    double Lx = A0 / Ly;

    std::vector<Vec2> r(N);
    for (int i=0; i<N_side; ++i) {
        for (int j=0; j<N_side; ++j) {
            int idx = i*N_side + j;
            r[idx][0] = (i + 0.5)/N_side * Lx;
            r[idx][1] = (j + 0.5)/N_side * Ly;
        }
    }

    std::mt19937_64 rng(42);
    std::normal_distribution<double> dist(0.0, std::sqrt(kB*T0/m));
    std::vector<Vec2> v(N);
    for (auto &vi : v) {
        vi[0] = dist(rng);
        vi[1] = dist(rng);
    }

    double vxm=0.0, vym=0.0;
    for (auto &vi : v) { vxm += vi[0]; vym += vi[1]; }
    vxm /= N; vym /= N;
    for (auto &vi : v) { vi[0] -= vxm; vi[1] -= vym; }

    {
        double Tinit = global_temperature(v);
        if (Tinit > 0) {
            double s = std::sqrt(T0 / Tinit);
            for (auto &vi : v) {
                vi[0] *= s;
                vi[1] *= s;
            }
        }
    }

    std::vector<Vec2> F(N);
    ForceResults fr;
    compute_forces(r, Lx, Ly, F, fr);

    std::vector<Mat2> strain_history;
    std::vector<Mat2> stress_history;
    std::vector<double> temp_history;
    std::vector<double> time_history;
    strain_history.reserve(steps / sample_every + 2);
    stress_history.reserve(steps / sample_every + 2);
    temp_history.reserve(steps / sample_every + 2);
    time_history.reserve(steps / sample_every + 2);

    std::cout << "Equilibrating (NPT Berendsen) ..." << std::endl;

    for (int step=0; step<equil_steps; ++step) {
        for (int i=0; i<N; ++i) {
            v[i][0] += 0.5 * dt * F[i][0] / m;
            v[i][1] += 0.5 * dt * F[i][1] / m;
        }
        for (int i=0; i<N; ++i) {
            r[i][0] += dt * v[i][0];
            r[i][1] += dt * v[i][1];
            wrap(r[i][0], Lx);
            wrap(r[i][1], Ly);
        }
        compute_forces(r, Lx, Ly, F, fr);
        for (int i=0; i<N; ++i) {
            v[i][0] += 0.5 * dt * F[i][0] / m;
            v[i][1] += 0.5 * dt * F[i][1] / m;
        }

        if (step % 10 == 0) {
            double Tcurr = global_temperature(v);
            if (Tcurr > 1e-14) {
                double lambda = std::sqrt(1.0 + (dt / tau_T) * (T_target / Tcurr - 1.0));
                lambda = std::clamp(lambda, 0.95, 1.05);
                for (auto &vi : v) {
                    vi[0] *= lambda;
                    vi[1] *= lambda;
                }
            }
        }

        if (step % 20 == 0) {
            double Tcurr = global_temperature(v);
            double V = Lx * Ly;
            double Pcurr = (N * kB * Tcurr + (fr.P_pair / dim)) / V;
            double mu = 1.0 - beta * (dt / tau_P) * (Pcurr - P_target);
            mu = std::clamp(mu, 0.995, 1.005);
            double scale = std::sqrt(mu);
            Lx *= scale;
            Ly *= scale;
            for (int i=0; i<N; ++i) {
                r[i][0] *= scale;
                r[i][1] *= scale;
                v[i][0] *= scale;
                v[i][1] *= scale;
                wrap(r[i][0], Lx);
                wrap(r[i][1], Ly);
            }
        }

        if ((step & 127) == 0) {
            double sx=0, sy=0;
            for (auto &vi : v) { sx += vi[0]; sy += vi[1]; }
            sx /= N; sy /= N;
            for (auto &vi : v) { vi[0] -= sx; vi[1] -= sy; }
        }

        if ((step % 5000) == 0) {
            double Tcurr = global_temperature(v);
            if (std::isnan(Tcurr) || std::isinf(Tcurr) || Tcurr > 10.0 * T0) {
                std::cerr << "Instability during equilibration at step " << step << "\n";
                return 1;
            }
            std::cout << "Equil step " << step << " T=" << Tcurr << std::endl;
        }
    }

    double Lx0 = Lx;
    double Ly0 = Ly;

    std::cout << "Deformation (NVT Berendsen) ..." << std::endl;

    for (int step=0; step<steps; ++step) {
        for (int i=0; i<N; ++i) {
            v[i][0] += 0.5 * dt * F[i][0] / m;
            v[i][1] += 0.5 * dt * F[i][1] / m;
        }
        for (int i=0; i<N; ++i) {
            r[i][0] += dt * v[i][0];
            r[i][1] += dt * v[i][1];
            wrap(r[i][0], Lx);
            wrap(r[i][1], Ly);
        }

        if ((step+1) % deform_every == 0) {
            double Lx_new = Lx * (1.0 + d_eps);
            double scale = Lx_new / Lx;
            for (int i=0; i<N; ++i) {
                double xnew = r[i][0] * scale;
                wrap(xnew, Lx_new);
                r[i][0] = xnew;
            }
            Lx = Lx_new;
        }

        compute_forces(r, Lx, Ly, F, fr);

        for (int i=0; i<N; ++i) {
            v[i][0] += 0.5 * dt * F[i][0] / m;
            v[i][1] += 0.5 * dt * F[i][1] / m;
        }

        if (step % 10 == 0) {
            double Tcurr = global_temperature(v);
            if (Tcurr > 1e-14) {
                double lambda = std::sqrt(1.0 + (dt / tau_T) * (T_target / Tcurr - 1.0));
                lambda = std::clamp(lambda, 0.98, 1.02);
                for (auto &vi : v) {
                    vi[0] *= lambda;
                    vi[1] *= lambda;
                }
            }
        }

        if ((step & 255) == 0) {
            double sx=0, sy=0;
            for (auto &vi : v) { sx += vi[0]; sy += vi[1]; }
            sx /= N; sy /= N;
            for (auto &vi : v) { vi[0] -= sx; vi[1] -= sy; }
        }

        if (step % sample_every == 0) {
            double F11 = Lx / Lx0;
            double F22 = Ly / Ly0;
            Mat2 E = {
                0.5 * (F11*F11 - 1.0), 0.0,
                0.0, 0.5 * (F22*F22 - 1.0)
            };
            double V = Lx * Ly;
            Mat2 sigma = compute_cauchy_stress(v, fr, V);
            double Tcurr = global_temperature(v);

            if (std::isnan(Tcurr) || std::isinf(Tcurr) || Tcurr > 5.0 * T0) {
                std::cerr << "Instability at deformation step " << step << "\n";
                break;
            }

            strain_history.push_back(E);
            stress_history.push_back(sigma);
            temp_history.push_back(Tcurr);
            time_history.push_back(step * dt);

            if (step % 5000 == 0) {
                std::cout << "Step " << step
                          << " T=" << Tcurr
                          << " Exx=" << E[0]
                          << " sigma_xx=" << sigma[0]
                          << std::endl;
            }
        }
    }

    int nbins = 20;
    int n_avg_steps = 5000; // Number of steps to average over for profile
    std::vector<double> Tbin_sum(nbins, 0.0);
    std::vector<int> Tbin_count(nbins, 0);
    std::vector<double> ycent(nbins);

    for (int b=0; b<nbins; ++b)
        ycent[b] = (b + 0.5) * Ly / nbins;

    // Time-averaged temperature profile over last n_avg_steps
    for (int avg_step = 0; avg_step < n_avg_steps; ++avg_step) {
        auto Tbins = temperature_profile(v, r, Ly, nbins);
        for (int b = 0; b < nbins; ++b) {
            Tbin_sum[b] += Tbins[b];
            Tbin_count[b] += 1;
        }
        // Advance MD (one velocity-verlet step)
        for (int i=0; i<N; ++i) {
            v[i][0] += 0.5 * dt * F[i][0] / m;
            v[i][1] += 0.5 * dt * F[i][1] / m;
        }
        for (int i=0; i<N; ++i) {
            r[i][0] += dt * v[i][0];
            r[i][1] += dt * v[i][1];
            wrap(r[i][0], Lx);
            wrap(r[i][1], Ly);
        }
        compute_forces(r, Lx, Ly, F, fr);
        for (int i=0; i<N; ++i) {
            v[i][0] += 0.5 * dt * F[i][0] / m;
            v[i][1] += 0.5 * dt * F[i][1] / m;
        }
    }

    std::vector<double> Tbin_avg(nbins, 0.0);
    for (int b = 0; b < nbins; ++b) {
        Tbin_avg[b] = Tbin_sum[b] / Tbin_count[b];
    }

    // Optionally apply moving average smoothing with window=1
    std::vector<double> Tbin_smooth = moving_average(Tbin_avg, 1);

    // Filenames with T, N, rho
    std::string suffix = "_T_" + std::to_string(T0)
                       + "_N_" + std::to_string(N)
                       + "_rho_" + std::to_string(rho);

    {
        std::ofstream f("temperature_profile" + suffix + ".csv");
        f << "y_center,temperature\n";
        for (int b=0; b<nbins; ++b)
            f << ycent[b] << "," << Tbin_smooth[b] << "\n";
    }
    {
        std::ofstream f("stress_strain_xx" + suffix + ".csv");
        f << "E_xx,sigma_xx\n";
        for (size_t i=0; i<strain_history.size(); ++i)
            f << strain_history[i][0] << "," << stress_history[i][0] << "\n";
    }
    {
        std::ofstream f("stress_strain_xy" + suffix + ".csv");
        f << "E_xy,sigma_xy\n";
        for (size_t i=0; i<strain_history.size(); ++i)
            f << strain_history[i][1] << "," << stress_history[i][1] << "\n";
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t1 - t0;
    std::cout << "Completed in " << elapsed.count() << " s\n";
    std::cout << "Final axial Green-Lagrange strain Exx: "
              << (strain_history.empty() ? 0.0 : strain_history.back()[0]) << "\n";
    std::cout << "Data saved with suffix " << suffix << std::endl;
    return 0;
}