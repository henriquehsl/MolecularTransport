#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <ctime>
using namespace std;

default_random_engine gen;
normal_distribution<double> d(0.0, 1.0);

#define deltaT 0.125
#define eps 2.0

class alphaXY {
public:
    double E_c(const vector<double>& v, int N);
    void initial_cond(vector<double>& x, vector<double>& v, int N);
    void a_func(const vector<double>& x, const vector<double>& v, vector<double>& a, const vector<double>& W, double Const, int N, double nu_l, double nu_r, double T, double dt);
    void verlet(vector<double>& x, vector<double>& v, vector<double>& a, const vector<double>& W, int N, double nu_l, double nu_r, double Const, double T, double dt);
    long double J_flux(const vector<double>& x, const vector<double>& v, const vector<double>& W, int N, int l);
    void transient(vector<double>& x, vector<double>& v, vector<double>& a, const vector<double>& W, double alpha, long long int tran, int N, double Const, double T, double dt);
    void experiment(vector<double>& x, vector<double>& v, vector<double>& a, const vector<double>& W, double alpha, vector<long double>& kappa, char* info, int n, ofstream& out, long long int tran, long long int tsteps, int N, double Const, double T, double dt);
    void simulation(vector<double>& x, vector<double>& v, vector<double>& a, const vector<double>& W, double alpha, vector<long double>& kappa, char* info, ofstream& out, int Nexp, long long int tran, long long int tsteps, int N, double Const, double T, double dt1, double dt2);
};

double alphaXY::E_c(const vector<double>& v, int N) {
    double Ek = 0.0;
    for (int i = 0; i < N; ++i)
        Ek += 0.5 * v[i] * v[i];
    return Ek;
}

void alphaXY::initial_cond(vector<double>& x, vector<double>& v, int N) {
    for (int i = 0; i < N; ++i) {
        x[i] = 0.0;
        v[i] = 0.0;
    }
}

void alphaXY::a_func(const vector<double>& x, const vector<double>& v, vector<double>& a, const vector<double>& W, double Const, int N, double nu_l, double nu_r, double T, double dt) {
    double T_l = T * (1 + deltaT), T_r = T * (1 - deltaT);
    double sqrt2Tldt = sqrt(2.0 * T_l / dt);
    double sqrt2Trdt = sqrt(2.0 * T_r / dt);

    for (int i = 0; i < N; ++i) {
        a[i] = 0.0;
        for (int j = 0; j < N; ++j) {
            if (j != i)
                a[i] += Const * sin(x[j] - x[i]) * W[N * i + j];
        }
    }
    a[0] += -v[0] + nu_l * sqrt2Tldt;
    a[N - 1] += -v[N - 1] + nu_r * sqrt2Trdt;
}

void alphaXY::verlet(vector<double>& x, vector<double>& v, vector<double>& a, const vector<double>& W, int N, double nu_l, double nu_r, double Const, double T, double dt) {
    static vector<double> a_now;
    if (a_now.size() != N) a_now.resize(N);

    a_func(x, v, a, W, Const, N, nu_l, nu_r, T, dt);
    for (int i = 0; i < N; ++i) {
        a_now[i] = a[i];
        x[i] += v[i] * dt + 0.5 * a_now[i] * dt * dt;
    }
    a_func(x, v, a, W, Const, N, nu_l, nu_r, T, dt);
    for (int i = 0; i < N; ++i) {
        v[i] += 0.5 * dt * (a_now[i] + a[i]);
    }
}

long double alphaXY::J_flux(const vector<double>& x, const vector<double>& v, const vector<double>& W, int N, int l) {
    double psi = 0.0;
    for (int k = 0; k < l; ++k) {
        psi += (v[l] + v[k]) * sin(x[k] - x[l]) * W[N * l + k];
    }
    return psi;
}

void alphaXY::transient(vector<double>& x, vector<double>& v, vector<double>& a, const vector<double>& W, double alpha, long long int tran, int N, double Const, double T, double dt) {
    initial_cond(x, v, N);
    double nu_r, nu_l;
    for (long long int i = 0; i < tran; ++i) {
        nu_l = d(gen);
        nu_r = d(gen);
        verlet(x, v, a, W, N, nu_l, nu_r, Const, T, dt);
    }
}

void alphaXY::experiment(vector<double>& x, vector<double>& v, vector<double>& a, const vector<double>& W, double alpha, vector<long double>& kappa, char* info, int n, ofstream& out, long long int tran, long long int tsteps, int N, double Const, double T, double dt) {
    int p = static_cast<int>(std::round(0.15 * N));
    for (int ii = 0; ii < N; ++ii) kappa[ii] = 0.0;
    double nu_r, nu_l;
    for (long long int i = 0; i < tsteps; ++i) {
        nu_l = d(gen);
        nu_r = d(gen);
        verlet(x, v, a, W, N, nu_l, nu_r, Const, T, dt);
        for (int ii = p; ii < N - p; ++ii) {
            kappa[ii] += J_flux(x, v, W, N, ii);
        }
    }
    double DT = 2.0 * deltaT * T;
    for (int i = p; i < N - p; ++i) {
        kappa[i] = (kappa[i] * N * Const) / (2.0 * DT * tsteps);
    }
}

void alphaXY::simulation(vector<double>& x, vector<double>& v, vector<double>& a, const vector<double>& W, double alpha, vector<long double>& kappa, char* info, ofstream& out, int Nexp, long long int tran, long long int tsteps, int N, double Const, double T, double dt1, double dt2) {
    int p = static_cast<int>(std::round(0.15 * N));
    sprintf(info, "exp_%d_N_%d_alpha_%0.2f_T_%0.4f_trans_%dM_time_%dM_BULK_1dalphaXY.txt", Nexp, N, alpha, T, int(tran / 1e6), int(tsteps / 1e6));
    out.open(info);
    transient(x, v, a, W, alpha, tran, N, Const, T, dt1);
    double kappa_av = 0.0;
    for (int n = 1; n <= Nexp; ++n) {
        experiment(x, v, a, W, alpha, kappa, info, n, out, tran, tsteps, N, Const, T, dt2);
        for (int i = p; i < N - p; ++i) {
            kappa_av += kappa[i];
        }
    }
    kappa_av = kappa_av / (Nexp * (N - 2 * p));
    out << kappa_av / N << "\t" << endl;
    out.close();
}

int main() {
    clock_t begin = clock();
    srand(time(NULL));

    int N = 20;
    long long int tran = 165000;
    long long int tsteps = 40000;
    const double T = 0.03;
    double alpha = 0.0;
    const int Nexp = 20;

    cout << "T= " << T << "\t" << N << "\t" << tran << "\t" << tsteps << endl;

    char info[90];
    ofstream out;
    vector<long double> kappa(N);
    vector<double> W(N * N);

    // Fill weight matrix
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double val = pow(abs(i - j), -1.0 * alpha);
            W[N * i + j] = val;
            W[N * j + i] = val;
        }
        W[N * i + i] = 0.0;
    }

    // Compute normalization
    double N_tilde = 0.0, Const;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j)
            N_tilde += 2 * W[N * i + j] / N;
    }
    Const = eps / N_tilde;

    vector<double> x(N), v(N), a(N);

    alphaXY o;
    o.simulation(x, v, a, W, alpha, kappa, info, out, Nexp, tran, tsteps, N, Const, T, 0.01, 0.01);

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "time elapsed: " << elapsed_secs << " seconds" << endl;
}