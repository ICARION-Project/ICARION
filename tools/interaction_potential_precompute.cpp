// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include <cstdlib>
#include <array>
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <fstream>
#include <limits>
#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <H5Cpp.h>
#include <json/json.h>

#include "utils/constants.h"
#include "core/types/Vec3.h"
#include "core/io/moleculeLoader.h"
#include "core/physics/collisions/InteractionPotentialOfflineSampleSet.h"

namespace {

struct Options {
    std::string input;
    std::string output;
    std::string species_id;
    std::string gas = "He";
    std::string gas_model = "mono";
    std::string mixing_rule = "lb";
    std::string polarization = "partial";
    std::string potential = "lj1264";
    std::string param_model = "sigma_epsilon";
    std::string format = "hdf5";
    std::string orient_grid = "lebedev";
    int n_orientations = 50;
    int n_trials = 20000;
    int threads = 0;
    int v_bins = 32;
    double v_min = 0.0;
    double v_max = 0.0;
    double temperature_K = 0.0;
    double epsilon_deflection = 1e-4;
    double b_guess_m = 20.0 * ANGSTROM_TO_M;
    double b_growth = 2.0;
    double b_rel_tol = 1e-3;
    double eta_dt = 0.02;
    double sigma_scale = 1.0;
    double epsilon_scale = 1.0;
    double pol_damp_A = 0.0;
    double mmff_energy_scale = 1.0;
    double mmff_distance_scale = 1.0;
    double n2_bond_A = 1.0976;
    bool n2_quadrupole = false;
    double n2_q_site_e = -0.4825;
    double n2_q_center_e = 0.965;
    bool n2_average_lj = true;
    bool n2_aniso_pol = false;
    double n2_alpha_par_A3 = 0.0;
    double n2_alpha_perp_A3 = 0.0;
    bool store_full_cdf = false;
    int checkpoint_cells = 0;
    bool resume = false;
    uint64_t seed = 0;
    double gas_sigma_m = 0.0;
    double gas_epsilon_J = 0.0;
    std::string gas_params_file = "data/forcefields/gas_lj_params.json";
    std::string element_params_file;
    std::string lebedev_file;
    bool scan_deflection = false;
    double scan_b_min_m = 0.0;
    double scan_b_max_m = 0.0;
    int scan_b_steps = 200;
    double scan_v_rel = 0.0;
    int scan_orient_index = 0;
    std::string scan_paths_output;
    int scan_path_stride = 10;
};

void print_usage() {
    std::cout
        << "Usage: interaction_potential_precompute --input <species_database.json> --species <id> --output <file.h5> [options]\\n"
        << "Options:\\n"
        << "  --gas <He|N2|Ar|...>        Gas species (default: He)\\n"
        << "  --gas-model <mono|diatomic> Gas model (default: mono; diatomic supported for N2)\\n"
        << "  --mixing-rule <lb|mmff|pair> Mixing rule for sigma_ij/epsilon_ij (default: lb)\\n"
        << "  --polarization <partial|total|pairwise|none>  Polarization model (default: partial)\\n"
        << "  --potential <lj1264|exp6>   Potential model (default: lj1264)\\n"
        << "  --param-model <sigma_epsilon|mmff_reij>  Element parameterization (default: sigma_epsilon)\\n"
        << "  --format <hdf5|json>        Output format (default: hdf5)\\n"
        << "  --orient-grid <lebedev|qmc|random>  Orientation grid (default: lebedev)\\n"
        << "  --n-orientations <int>      Number of orientations (default: 50)\\n"
        << "  --n-trials <int>            Monte Carlo trials per (E,Ω) (default: 20000)\\n"
        << "  --threads <int>             OpenMP thread count (0 = runtime default)\\n"
        << "  --v-bins <int>              log(v_rel) bins (default: 32)\\n"
        << "  --v-min <float>             Minimum v_rel (m/s), default from temperature\\n"
        << "  --v-max <float>             Maximum v_rel (m/s), default from temperature\\n"
        << "  --temperature-K <float>     Temperature for default v range (default: species ref or 300K)\\n"
        << "  --epsilon <float>           Deflection tolerance ε (default: 1e-4)\\n"
        << "  --b-guess-A <float>         Initial b_guess in Angstrom (default: 20)\\n"
        << "  --b-growth <float>          b growth factor (default: 2.0)\\n"
        << "  --b-rel-tol <float>         b_max binary search tolerance (default: 1e-3)\\n"
        << "  --eta-dt <float>            dt control factor η (default: 0.02)\\n"
        << "  --sigma-scale <float>       Scale factor for LJ sigma (default: 1.0)\\n"
        << "  --epsilon-scale <float>     Scale factor for LJ epsilon (default: 1.0)\\n"
        << "  --pol-damp-A <float>        Polarization damping radius in Angstrom (default: 0, disabled)\\n"
        << "  --mmff-ener-scale <float>   MMFF e_ij energy scale (default: 1.0)\\n"
        << "  --mmff-dist-scale <float>   MMFF R* distance scale (default: 1.0)\\n"
        << "  --n2-bond-A <float>         N2 bond length in Angstrom (default: 1.0976)\\n"
        << "  --n2-quadrupole             Enable N2 quadrupole (5-charge) model (default: off)\\n"
        << "  --n2-q-site-e <float>       N2 site charge in e (default: -0.4825)\\n"
        << "  --n2-q-center-e <float>     N2 center charge in e (default: +0.965)\\n"
        << "  --n2-average-lj <0|1>       Average LJ over two N sites (default: 1)\\n"
        << "  --n2-aniso-pol              Enable anisotropic N2 polarizability (pairwise model)\\n"
        << "  --n2-alpha-par-A3 <float>   N2 parallel polarizability in A^3 (optional)\\n"
        << "  --n2-alpha-perp-A3 <float>  N2 perpendicular polarizability in A^3 (optional)\\n"
        << "  --gas-sigma-A <float>       Override gas sigma in Angstrom\\n"
        << "  --gas-epsilon-eV <float>    Override gas epsilon in eV\\n"
        << "  --gas-params <file>         Gas LJ params JSON (default: data/forcefields/gas_lj_params.json)\\n"
        << "  --element-params <file>     Element LJ params JSON (default: gas params file)\\n"
        << "  --lebedev-file <file>       Lebedev grid file (default: data/forcefields/lebedev_011.txt)\\n"
        << "  --store-full-cdf            Store full Δp CDF (default: off)\\n"
        << "  --checkpoint-cells <int>    Write HDF5 checkpoint every N completed (Ω,v) cells (default: 0=final only)\\n"
        << "  --resume                    Resume from existing HDF5 checkpoint file\\n"
        << "  --seed <uint64>             RNG seed (default: 0 = hash(species,gas,version))\\n"
        << "  --scan-deflection           Output deflection scan (θ(b), 1−cosθ) to --output CSV\\n"
        << "  --scan-v-rel <float>        Relative speed for deflection scan (m/s)\\n"
        << "  --scan-b-min-A <float>      Minimum impact parameter for scan (Angstrom)\\n"
        << "  --scan-b-max-A <float>      Maximum impact parameter for scan (Angstrom, default: b_max(E))\\n"
        << "  --scan-b-steps <int>        Number of b samples for scan (default: 200)\\n"
        << "  --scan-orient-idx <int>     Orientation index for scan (default: 0)\\n"
        << "  --scan-paths-out <file>     Optional CSV output for actual trajectory points\\n"
        << "  --scan-path-stride <int>    Store every N accepted steps in scan-paths-out (default: 10)\\n";
}

bool parse_args(int argc, char** argv, Options& opt) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--input" && i + 1 < argc) {
            opt.input = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            opt.output = argv[++i];
        } else if (arg == "--species" && i + 1 < argc) {
            opt.species_id = argv[++i];
        } else if (arg == "--gas" && i + 1 < argc) {
            opt.gas = argv[++i];
        } else if (arg == "--gas-model" && i + 1 < argc) {
            opt.gas_model = argv[++i];
        } else if (arg == "--mixing-rule" && i + 1 < argc) {
            opt.mixing_rule = argv[++i];
        } else if (arg == "--polarization" && i + 1 < argc) {
            opt.polarization = argv[++i];
        } else if (arg == "--potential" && i + 1 < argc) {
            opt.potential = argv[++i];
        } else if (arg == "--param-model" && i + 1 < argc) {
            opt.param_model = argv[++i];
        } else if (arg == "--format" && i + 1 < argc) {
            opt.format = argv[++i];
        } else if (arg == "--orient-grid" && i + 1 < argc) {
            opt.orient_grid = argv[++i];
        } else if (arg == "--n-orientations" && i + 1 < argc) {
            opt.n_orientations = std::stoi(argv[++i]);
        } else if (arg == "--n-trials" && i + 1 < argc) {
            opt.n_trials = std::stoi(argv[++i]);
        } else if (arg == "--threads" && i + 1 < argc) {
            opt.threads = std::stoi(argv[++i]);
        } else if (arg == "--v-bins" && i + 1 < argc) {
            opt.v_bins = std::stoi(argv[++i]);
        } else if (arg == "--v-min" && i + 1 < argc) {
            opt.v_min = std::stod(argv[++i]);
        } else if (arg == "--v-max" && i + 1 < argc) {
            opt.v_max = std::stod(argv[++i]);
        } else if (arg == "--temperature-K" && i + 1 < argc) {
            opt.temperature_K = std::stod(argv[++i]);
        } else if (arg == "--epsilon" && i + 1 < argc) {
            opt.epsilon_deflection = std::stod(argv[++i]);
        } else if (arg == "--b-guess-A" && i + 1 < argc) {
            opt.b_guess_m = std::stod(argv[++i]) * ANGSTROM_TO_M;
        } else if (arg == "--b-growth" && i + 1 < argc) {
            opt.b_growth = std::stod(argv[++i]);
        } else if (arg == "--b-rel-tol" && i + 1 < argc) {
            opt.b_rel_tol = std::stod(argv[++i]);
        } else if (arg == "--eta-dt" && i + 1 < argc) {
            opt.eta_dt = std::stod(argv[++i]);
        } else if (arg == "--sigma-scale" && i + 1 < argc) {
            opt.sigma_scale = std::stod(argv[++i]);
        } else if (arg == "--epsilon-scale" && i + 1 < argc) {
            opt.epsilon_scale = std::stod(argv[++i]);
        } else if (arg == "--pol-damp-A" && i + 1 < argc) {
            opt.pol_damp_A = std::stod(argv[++i]);
        } else if (arg == "--mmff-ener-scale" && i + 1 < argc) {
            opt.mmff_energy_scale = std::stod(argv[++i]);
        } else if (arg == "--mmff-dist-scale" && i + 1 < argc) {
            opt.mmff_distance_scale = std::stod(argv[++i]);
        } else if (arg == "--n2-bond-A" && i + 1 < argc) {
            opt.n2_bond_A = std::stod(argv[++i]);
        } else if (arg == "--n2-quadrupole") {
            opt.n2_quadrupole = true;
        } else if (arg == "--n2-q-site-e" && i + 1 < argc) {
            opt.n2_q_site_e = std::stod(argv[++i]);
        } else if (arg == "--n2-q-center-e" && i + 1 < argc) {
            opt.n2_q_center_e = std::stod(argv[++i]);
        } else if (arg == "--n2-average-lj" && i + 1 < argc) {
            opt.n2_average_lj = std::stoi(argv[++i]) != 0;
        } else if (arg == "--n2-aniso-pol") {
            opt.n2_aniso_pol = true;
        } else if (arg == "--n2-alpha-par-A3" && i + 1 < argc) {
            opt.n2_alpha_par_A3 = std::stod(argv[++i]);
        } else if (arg == "--n2-alpha-perp-A3" && i + 1 < argc) {
            opt.n2_alpha_perp_A3 = std::stod(argv[++i]);
        } else if (arg == "--gas-sigma-A" && i + 1 < argc) {
            opt.gas_sigma_m = std::stod(argv[++i]) * ANGSTROM_TO_M;
        } else if (arg == "--gas-epsilon-eV" && i + 1 < argc) {
            opt.gas_epsilon_J = std::stod(argv[++i]) * ELEM_CHARGE_C;
        } else if (arg == "--gas-params" && i + 1 < argc) {
            opt.gas_params_file = argv[++i];
        } else if (arg == "--element-params" && i + 1 < argc) {
            opt.element_params_file = argv[++i];
        } else if (arg == "--lebedev-file" && i + 1 < argc) {
            opt.lebedev_file = argv[++i];
        } else if (arg == "--store-full-cdf") {
            opt.store_full_cdf = true;
        } else if (arg == "--checkpoint-cells" && i + 1 < argc) {
            opt.checkpoint_cells = std::stoi(argv[++i]);
        } else if (arg == "--resume") {
            opt.resume = true;
        } else if (arg == "--seed" && i + 1 < argc) {
            opt.seed = static_cast<uint64_t>(std::stoull(argv[++i]));
        } else if (arg == "--scan-deflection") {
            opt.scan_deflection = true;
        } else if (arg == "--scan-v-rel" && i + 1 < argc) {
            opt.scan_v_rel = std::stod(argv[++i]);
        } else if (arg == "--scan-b-min-A" && i + 1 < argc) {
            opt.scan_b_min_m = std::stod(argv[++i]) * ANGSTROM_TO_M;
        } else if (arg == "--scan-b-max-A" && i + 1 < argc) {
            opt.scan_b_max_m = std::stod(argv[++i]) * ANGSTROM_TO_M;
        } else if (arg == "--scan-b-steps" && i + 1 < argc) {
            opt.scan_b_steps = std::stoi(argv[++i]);
        } else if (arg == "--scan-orient-idx" && i + 1 < argc) {
            opt.scan_orient_index = std::stoi(argv[++i]);
        } else if (arg == "--scan-paths-out" && i + 1 < argc) {
            opt.scan_paths_output = argv[++i];
        } else if (arg == "--scan-path-stride" && i + 1 < argc) {
            opt.scan_path_stride = std::stoi(argv[++i]);
        } else if (arg == "--help" || arg == "-h") {
            return false;
        }
    }

    if (opt.input.empty() || opt.output.empty() || opt.species_id.empty()) {
        return false;
    }
    if (opt.n_orientations <= 0) {
        return false;
    }
    if (opt.scan_deflection) {
        if (opt.scan_v_rel <= 0.0 || opt.scan_b_steps < 2 || opt.scan_path_stride <= 0) {
            return false;
        }
        return true;
    }
    if (opt.n_trials <= 0 || opt.v_bins <= 0) {
        return false;
    }
    if (opt.threads < 0) {
        return false;
    }
    if (opt.checkpoint_cells < 0) {
        return false;
    }
    return true;
}

double radical_inverse(uint64_t n, int base) {
    double inv_base = 1.0 / static_cast<double>(base);
    double fraction = 0.0;
    double factor = inv_base;
    while (n > 0) {
        const int digit = static_cast<int>(n % static_cast<uint64_t>(base));
        fraction += digit * factor;
        n /= static_cast<uint64_t>(base);
        factor *= inv_base;
    }
    return fraction;
}

inline double dot_local(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double norm2(const Vec3& v) {
    return dot_local(v, v);
}

inline double norm(const Vec3& v) {
    return std::sqrt(norm2(v));
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

Vec3 rotate_vec_by_quat(const std::array<double, 4>& q_wxyz, const Vec3& v) {
    const double w = q_wxyz[0];
    const Vec3 u(q_wxyz[1], q_wxyz[2], q_wxyz[3]);
    const double u2 = dot_local(u, u);
    const Vec3 term1 = u * (2.0 * dot_local(u, v));
    const Vec3 term2 = v * (w * w - u2);
    const Vec3 term3 = cross(u, v) * (2.0 * w);
    return term1 + term2 + term3;
}

std::array<double, 4> sample_shoemake(double u1, double u2, double u3) {
    const double sqrt1 = std::sqrt(1.0 - u1);
    const double sqrt2 = std::sqrt(u1);
    const double theta1 = 2.0 * M_PI * u2;
    const double theta2 = 2.0 * M_PI * u3;
    const double x = sqrt1 * std::sin(theta1);
    const double y = sqrt1 * std::cos(theta1);
    const double z = sqrt2 * std::sin(theta2);
    const double w = sqrt2 * std::cos(theta2);
    return {w, x, y, z};
}

Vec3 sample_unit_vector(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    const double u = uni(rng);
    const double v = uni(rng);
    const double z = 2.0 * u - 1.0;
    const double t = 2.0 * M_PI * v;
    const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
    return Vec3(r * std::cos(t), r * std::sin(t), z);
}

std::vector<std::array<double, 4>> generate_orientations_random(int n, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::vector<std::array<double, 4>> out;
    out.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        out.push_back(sample_shoemake(uni(rng), uni(rng), uni(rng)));
    }
    return out;
}

std::vector<std::array<double, 4>> generate_orientations_qmc(int n) {
    std::vector<std::array<double, 4>> out;
    out.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double u1 = radical_inverse(static_cast<uint64_t>(i + 1), 2);
        const double u2 = radical_inverse(static_cast<uint64_t>(i + 1), 3);
        const double u3 = radical_inverse(static_cast<uint64_t>(i + 1), 5);
        out.push_back(sample_shoemake(u1, u2, u3));
    }
    return out;
}

std::array<double, 4> quat_from_z_to_dir(const Vec3& dir) {
    const Vec3 z(0.0, 0.0, 1.0);
    const double cos_theta = dot_local(z, dir);
    if (cos_theta > 0.999999) {
        return {1.0, 0.0, 0.0, 0.0};
    }
    if (cos_theta < -0.999999) {
        return {0.0, 1.0, 0.0, 0.0};
    }
    const Vec3 axis = cross(z, dir);
    const double axis_norm = norm(axis);
    if (axis_norm <= 0.0) {
        return {1.0, 0.0, 0.0, 0.0};
    }
    const Vec3 axis_n = axis / axis_norm;
    const double angle = std::acos(std::max(-1.0, std::min(1.0, cos_theta)));
    const double half = 0.5 * angle;
    const double s = std::sin(half);
    return {std::cos(half), axis_n.x * s, axis_n.y * s, axis_n.z * s};
}

std::vector<std::array<double, 4>> load_lebedev_orientations(const std::string& path, int expected_points = -1) {
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("Failed to open Lebedev grid file: " + path);
    }
    std::vector<std::array<double, 4>> out;
    if (expected_points > 0) {
        out.reserve(static_cast<size_t>(expected_points));
    }
    double phi_deg = 0.0;
    double theta_deg = 0.0;
    double weight = 0.0;
    while (ifs >> phi_deg >> theta_deg >> weight) {
        const double phi = phi_deg * M_PI / 180.0;
        const double theta = theta_deg * M_PI / 180.0;
        const double sin_t = std::sin(theta);
        const Vec3 dir(sin_t * std::cos(phi), sin_t * std::sin(phi), std::cos(theta));
        out.push_back(quat_from_z_to_dir(dir));
    }
    if (expected_points > 0 && static_cast<int>(out.size()) != expected_points) {
        throw std::runtime_error("Lebedev grid size mismatch: expected " +
                                 std::to_string(expected_points) + ", got " +
                                 std::to_string(out.size()));
    }
    return out;
}

std::string lebedev_file_for_n(int n) {
    switch (n) {
        case 50:  return "data/forcefields/lebedev_011.txt";
        case 110: return "data/forcefields/lebedev_017.txt";
        case 194: return "data/forcefields/lebedev_023.txt";
        case 302: return "data/forcefields/lebedev_029.txt";
        case 590: return "data/forcefields/lebedev_041.txt";
        default:  return {};
    }
}

uint64_t default_seed_for(const std::string& species, const std::string& gas) {
    const std::string key = species + "|" + gas + "|ipm_v1";
    return static_cast<uint64_t>(std::hash<std::string>{}(key));
}

uint64_t splitmix64(uint64_t x) {
    x += 0x9E3779B97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

uint64_t seed_for_cell(uint64_t base, uint64_t oi, uint64_t ki) {
    uint64_t x = base;
    x ^= 0x9E3779B97f4a7c15ULL * (oi + 1);
    x ^= 0xBF58476D1CE4E5B9ULL * (ki + 1);
    return splitmix64(x);
}

struct AtomParam {
    Vec3 pos;
    double sigma;
    double epsilon;
    double inv_sigma;
    double lj_c6;
    double lj_c12;
    double exp6_pref_exp;
    double exp6_pref_disp;
    double c4;
    double q_c;
};

struct TrajectoryConfig {
    double eta_dt = 0.02;
    double energy_rel_tol_high = 1e-6;
    double energy_rel_tol_low = 1e-8;
    int relax_steps = 100;
    double r_cut_m = 0.0;
    int max_steps = 200000;
    double max_step_m = 0.0; // optional cap on |v|*dt
    double pol_damp_m = 0.0;
    bool store_path = false;
    int path_stride = 1;
};

struct TrajectoryResult {
    Vec3 v_out;
    double max_rel_energy_drift = 0.0;
    bool reached_asymptotic = false;
    std::vector<Vec3> path_points;
};

struct Lj1264Samples {
    std::vector<std::array<double, 4>> orientations_quat;
    std::vector<double> logv_bins;
    std::vector<double> sigma_mt;   // size N_orient * K
    std::vector<double> b_max;      // size N_orient * K
    std::vector<long long> cdf_offsets; // size N_orient * K
    std::vector<long long> cdf_counts;  // size N_orient * K
    std::vector<double> cdf_values;
    std::vector<double> dp_samples; // size total_samples * 3
    std::vector<double> dp_stats;   // size N_orient * K * 4
};

struct LocalCdfSamples {
    std::vector<double> cdf;
    std::vector<double> dp;
};

struct GasParams {
    double sigma_m = 0.0;
    double epsilon_J = 0.0;
    double alpha_m3 = 0.0;
    double alpha_par_m3 = 0.0;
    double alpha_perp_m3 = 0.0;
};

struct MmffGasParams {
    double alpha_A3 = 0.0;
    double N = 0.0;
    double A = 0.0;
    double G = 0.0;
};

enum class GasGeometry {
    Monoatomic,
    DiatomicN2
};

struct GasModelConfig {
    GasGeometry geometry = GasGeometry::Monoatomic;
    double bond_m = 0.0;
    bool quadrupole = false;
    double q_site_c = 0.0;
    double q_center_c = 0.0;
    bool average_lj = true;
    bool anisotropic_pol = false;
    double alpha_par_m3 = 0.0;
    double alpha_perp_m3 = 0.0;
};

struct ElementParams {
    double sigma_m = 0.0;
    double epsilon_eV = 0.0;
    bool has_sigma_eps = false;
    double alpha_A3 = 0.0;
    double N = 0.0;
    double A = 0.0;
    double G = 0.0;
    bool has_mmff = false;
};

// Tool-local enum types mirror CLI choices. They stay here because they are
// implementation details of the offline precompute algorithm, not runtime
// configuration enums.
enum class MixingRule {
    LB,
    MMFF,
    PAIR
};

enum class PotentialType {
    LJ1264,
    EXP6
};

enum class ParamModel {
    SIGMA_EPSILON,
    MMFF_REIJ
};

enum class PolarizationModel {
    PARTIAL,
    TOTAL,
    PAIRWISE,
    NONE
};

MixingRule mixing_rule_or_throw(const std::string& name) {
    if (name == "lb") return MixingRule::LB;
    if (name == "mmff") return MixingRule::MMFF;
    if (name == "pair" || name == "none") return MixingRule::PAIR;
    throw std::runtime_error("Unsupported mixing rule: " + name + " (use lb, mmff, or pair)");
}

PotentialType potential_type_or_throw(const std::string& name) {
    if (name == "lj1264") return PotentialType::LJ1264;
    if (name == "exp6") return PotentialType::EXP6;
    throw std::runtime_error("Unsupported potential: " + name + " (use lj1264 or exp6)");
}

ParamModel param_model_or_throw(const std::string& name) {
    if (name == "sigma_epsilon" || name == "sigma" || name == "lj") {
        return ParamModel::SIGMA_EPSILON;
    }
    if (name == "mmff_reij" || name == "mmff") {
        return ParamModel::MMFF_REIJ;
    }
    throw std::runtime_error("Unsupported param model: " + name + " (use sigma_epsilon or mmff_reij)");
}

PolarizationModel polarization_model_or_throw(const std::string& name) {
    if (name == "partial") return PolarizationModel::PARTIAL;
    if (name == "total" || name == "center") return PolarizationModel::TOTAL;
    if (name == "pairwise") return PolarizationModel::PAIRWISE;
    if (name == "none") return PolarizationModel::NONE;
    throw std::runtime_error("Unsupported polarization model: " + name + " (use partial, total, pairwise, or none)");
}

double mix_sigma(MixingRule rule, double sigma_i, double sigma_g) {
    if (rule == MixingRule::PAIR) {
        return sigma_i;
    }
    if (rule == MixingRule::LB) {
        return 0.5 * (sigma_i + sigma_g);
    }
    const double rsum = sigma_i + sigma_g;
    if (rsum <= 0.0) {
        return 0.0;
    }
    const double gamma = (sigma_i - sigma_g) / rsum;
    const double coeff1 = -MMFF_MIXING_BETA * gamma * gamma;
    return 0.5 * rsum * (1.0 + MMFF_MIXING_B * (1.0 - std::exp(coeff1)));
}

double mix_epsilon(MixingRule rule, double epsilon_i_J, double epsilon_g_J) {
    if (rule == MixingRule::PAIR) {
        return epsilon_i_J;
    }
    if (epsilon_i_J <= 0.0 || epsilon_g_J <= 0.0) {
        return 0.0;
    }
    return std::sqrt(epsilon_i_J * epsilon_g_J);
}

GasParams load_gas_lj_params_or_throw(const std::string& path, const std::string& gas) {
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("Failed to open gas params file: " + path);
    }
    Json::CharReaderBuilder builder;
    Json::Value root;
    std::string errs;
    if (!Json::parseFromStream(builder, ifs, &root, &errs)) {
        throw std::runtime_error("Failed to parse gas params file: " + errs);
    }
    if (!root.isMember("gases") || !root["gases"].isObject()) {
        throw std::runtime_error("Gas params file must contain a 'gases' object");
    }
    const auto& gases = root["gases"];
    if (!gases.isMember(gas)) {
        throw std::runtime_error("Gas params file has no entry for: " + gas);
    }
    const auto& entry = gases[gas];
    if (!entry.isMember("sigma_A") && !entry.isMember("sigma")) {
        throw std::runtime_error("Gas params entry missing sigma for: " + gas);
    }
    if (!entry.isMember("epsilon_eV") && !entry.isMember("epsilon")) {
        throw std::runtime_error("Gas params entry missing epsilon for: " + gas);
    }
    GasParams gp;
    gp.sigma_m = (entry.isMember("sigma_A") ? entry["sigma_A"].asDouble() : entry["sigma"].asDouble()) * ANGSTROM_TO_M;
    gp.epsilon_J = (entry.isMember("epsilon_eV") ? entry["epsilon_eV"].asDouble() : entry["epsilon"].asDouble()) * ELEM_CHARGE_C;
    if (entry.isMember("alpha")) {
        gp.alpha_m3 = entry["alpha"].asDouble() * ANGSTROM3_TO_M3;
    } else if (entry.isMember("alpha_A3")) {
        gp.alpha_m3 = entry["alpha_A3"].asDouble() * ANGSTROM3_TO_M3;
    }
    if (entry.isMember("alpha_parallel_A3")) {
        gp.alpha_par_m3 = entry["alpha_parallel_A3"].asDouble() * ANGSTROM3_TO_M3;
    }
    if (entry.isMember("alpha_perp_A3")) {
        gp.alpha_perp_m3 = entry["alpha_perp_A3"].asDouble() * ANGSTROM3_TO_M3;
    }
    return gp;
}

std::unordered_map<std::string, ElementParams> load_element_params_or_throw(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("Failed to open element params file: " + path);
    }
    Json::CharReaderBuilder builder;
    Json::Value root;
    std::string errs;
    if (!Json::parseFromStream(builder, ifs, &root, &errs)) {
        throw std::runtime_error("Failed to parse element params file: " + errs);
    }
    if (!root.isMember("elements") || !root["elements"].isObject()) {
        throw std::runtime_error("Element params file must contain an 'elements' object");
    }
    std::unordered_map<std::string, ElementParams> out;
    const auto& elems = root["elements"];
    for (const auto& key : elems.getMemberNames()) {
        const auto& entry = elems[key];
        ElementParams ep;
        if (entry.isMember("sigma_A") || entry.isMember("sigma")) {
            ep.sigma_m = (entry.isMember("sigma_A") ? entry["sigma_A"].asDouble() : entry["sigma"].asDouble()) * ANGSTROM_TO_M;
            ep.epsilon_eV = entry.isMember("epsilon_eV") ? entry["epsilon_eV"].asDouble() : entry["epsilon"].asDouble();
            ep.has_sigma_eps = true;
        }
        if (entry.isMember("alpha_A3") || entry.isMember("alpha")) {
            ep.alpha_A3 = entry.isMember("alpha_A3") ? entry["alpha_A3"].asDouble() : entry["alpha"].asDouble();
        }
        if (entry.isMember("N")) {
            ep.N = entry["N"].asDouble();
        } else if (entry.isMember("Ni")) {
            ep.N = entry["Ni"].asDouble();
        }
        if (entry.isMember("A")) {
            ep.A = entry["A"].asDouble();
        } else if (entry.isMember("Ai")) {
            ep.A = entry["Ai"].asDouble();
        }
        if (entry.isMember("G")) {
            ep.G = entry["G"].asDouble();
        } else if (entry.isMember("Gi")) {
            ep.G = entry["Gi"].asDouble();
        }
        if (ep.alpha_A3 > 0.0 && ep.N > 0.0 && ep.A > 0.0 && ep.G > 0.0) {
            ep.has_mmff = true;
        }
        out.emplace(key, ep);
    }
    return out;
}

GasParams gas_params_or_throw(const std::string& gas) {
    GasParams gp;
    if (gas == "He") {
        gp.alpha_m3 = POLARIZABILITY_HE_SI;
    } else if (gas == "Ar") {
        gp.alpha_m3 = POLARIZABILITY_AR_SI;
    } else if (gas == "CO2") {
        gp.alpha_m3 = POLARIZABILITY_CO2_SI;
    } else if (gas == "Ne") {
        gp.alpha_m3 = POLARIZABILITY_NE_SI;
    } else if (gas == "N2") {
        gp.alpha_m3 = POLARIZABILITY_N2_SI;
    } else if (gas == "O2") {
        gp.alpha_m3 = POLARIZABILITY_O2_SI;
    } else if (gas == "H2O") {
        gp.alpha_m3 = POLARIZABILITY_H2O_SI;
    } else {
        throw std::runtime_error("Unsupported gas for polarizability: " + gas);
    }

    return gp;
}

MmffGasParams mmff_gas_params_or_throw(const std::string& gas, double alpha_m3) {
    MmffGasParams gp;
    gp.alpha_A3 = alpha_m3 / ANGSTROM3_TO_M3;
    if (gas == "He") {
        gp.N = 1.42;
        gp.A = 4.40;
        gp.G = 1.209;
        return gp;
    }
    if (gas == "N2") {
        gp.N = 5.918;
        gp.A = 3.1781544;
        gp.G = 1.16175013;
        return gp;
    }
    throw std::runtime_error("MMFF R*/e_ij parameters only implemented for He and N2 gases.");
}

double gas_mass_kg_or_throw(const std::string& gas) {
    if (gas == "He") return MOLAR_MASS_HE_KG;
    if (gas == "Ar") return MOLAR_MASS_AR_KG;
    if (gas == "CO2") return MOLAR_MASS_CO2_KG;
    if (gas == "Ne") return MOLAR_MASS_NE_KG;
    if (gas == "N2") return MOLAR_MASS_N2_KG;
    if (gas == "O2") return MOLAR_MASS_O2_KG;
    if (gas == "H2O") return MOLAR_MASS_H2O_KG;
    throw std::runtime_error("Unsupported gas mass for: " + gas);
}

std::vector<AtomParam> build_atom_params_or_throw(
    const ICARION::io::Molecule& molecule,
    const std::string& gas_name,
    const GasParams& gas,
    double ion_charge_e,
    double sigma_scale,
    double epsilon_scale,
    double mmff_energy_scale,
    double mmff_distance_scale,
    const std::unordered_map<std::string, ElementParams>& element_params,
    MixingRule mixing_rule,
    PolarizationModel polarization_model,
    ParamModel param_model
) {
    std::vector<AtomParam> atoms;
    atoms.reserve(molecule.atoms.size());
    bool has_partial = false;
    MmffGasParams mmff_gas;
    if (param_model == ParamModel::MMFF_REIJ) {
        mmff_gas = mmff_gas_params_or_throw(gas_name, gas.alpha_m3);
    }
    for (const auto& atom : molecule.atoms) {
        auto it_elem = element_params.find(atom.element);
        if (it_elem == element_params.end()) {
            throw std::runtime_error("Missing element LJ parameters for: " + atom.element);
        }
        if (std::abs(atom.partial_charge_e) > 0.0) {
            has_partial = true;
        }
        AtomParam p;
        p.pos = atom.pos_m - molecule.center_of_mass_m;
        if (param_model == ParamModel::SIGMA_EPSILON) {
            // Standard path: combine per-element sigma/epsilon with the selected
            // gas parameters to create pairwise Lennard-Jones coefficients.
            if (!it_elem->second.has_sigma_eps) {
                throw std::runtime_error("Element params missing sigma/epsilon for: " + atom.element);
            }
            const double sigma_i = it_elem->second.sigma_m * sigma_scale;
            const double eps_eV = it_elem->second.epsilon_eV * epsilon_scale;
            const double sigma_ij = mix_sigma(mixing_rule, sigma_i, gas.sigma_m);
            const double epsilon_ij = mix_epsilon(
                mixing_rule,
                eps_eV * ELEM_CHARGE_C,
                gas.epsilon_J
            );
            p.sigma = sigma_ij;
            p.epsilon = epsilon_ij;
        } else {
            // MMFF path: reconstruct R*_ij and e_ij from MMFF-style alpha/N/A/G
            // parameters. The resulting pair length/energy are converted to SI
            // before the trajectory integrator sees them.
            if (!it_elem->second.has_mmff) {
                throw std::runtime_error("Element params missing MMFF (alpha/N/A/G) for: " + atom.element);
            }
            const double alpha_i = it_elem->second.alpha_A3;
            const double Ni = it_elem->second.N;
            const double Ai = it_elem->second.A;
            const double Gi = it_elem->second.G;
            const double RiiStar = Ai * std::pow(alpha_i, 0.25);
            const double RgasStar = mmff_gas.A * std::pow(mmff_gas.alpha_A3, 0.25);
            const double rsum = RiiStar + RgasStar;
            if (rsum <= 0.0) {
                throw std::runtime_error("Invalid MMFF R* for element: " + atom.element);
            }
            const double gamma = (RiiStar - RgasStar) / rsum;
            const double coeff1 = -MMFF_MIXING_BETA * gamma * gamma;
            double RijStar = 0.5 * rsum * (1.0 + MMFF_MIXING_B * (1.0 - std::exp(coeff1)));
            const double coeff2 = 181.16 * Gi * mmff_gas.G * alpha_i * mmff_gas.alpha_A3;
            const double coeff3 = std::sqrt(alpha_i / Ni) + std::sqrt(mmff_gas.alpha_A3 / mmff_gas.N);
            double eij = coeff2 / (coeff3 * std::pow(RijStar, 6.0));
            RijStar = RijStar * 1.0e-10 * mmff_distance_scale;
            eij = eij * (KCAL_TO_J / AVOGADRO_CONSTANT) * mmff_energy_scale;
            eij = eij / 1.1195;
            p.sigma = RijStar;
            p.epsilon = eij;
        }
        if (p.sigma > 0.0 && p.epsilon != 0.0) {
            p.inv_sigma = 1.0 / p.sigma;
            const double sigma2 = p.sigma * p.sigma;
            const double sigma6 = sigma2 * sigma2 * sigma2;
            const double sigma12 = sigma6 * sigma6;
            // Cache coefficients used in the hot force loop. This keeps the
            // trajectory integrator from repeatedly rebuilding powers of sigma.
            p.lj_c6 = 4.0 * p.epsilon * sigma6;
            p.lj_c12 = 4.0 * p.epsilon * sigma12;
            p.exp6_pref_exp = p.epsilon * EXP6_PREEXP_FACTOR;
            p.exp6_pref_disp = p.epsilon * EXP6_DISPERSION_FACTOR;
        } else {
            p.inv_sigma = 0.0;
            p.lj_c6 = 0.0;
            p.lj_c12 = 0.0;
            p.exp6_pref_exp = 0.0;
            p.exp6_pref_disp = 0.0;
        }
        const double q_c = atom.partial_charge_e * ELEM_CHARGE_C;
        p.q_c = q_c;
        if (polarization_model == PolarizationModel::PARTIAL) {
            // Partial-charge polarization is represented as a per-charge C4/r^4
            // attraction against the neutral gas polarizability.
            p.c4 = 0.5 * gas.alpha_m3 * COULOMB_CONST * q_c * q_c;
        } else {
            p.c4 = 0.0;
        }
        atoms.push_back(p);
    }

    if (polarization_model == PolarizationModel::TOTAL ||
        ((polarization_model == PolarizationModel::PARTIAL ||
          polarization_model == PolarizationModel::PAIRWISE) && !has_partial)) {
        const double q_c = ion_charge_e * ELEM_CHARGE_C;
        if (std::abs(q_c) > 0.0) {
            // If no partial-charge representation is available, use a single
            // ion-center charge so polarization is still present and explicit.
            AtomParam center;
            center.pos = Vec3(0.0, 0.0, 0.0);
            center.sigma = 0.0;
            center.epsilon = 0.0;
            center.inv_sigma = 0.0;
            center.lj_c6 = 0.0;
            center.lj_c12 = 0.0;
            center.exp6_pref_exp = 0.0;
            center.exp6_pref_disp = 0.0;
            center.q_c = q_c;
            center.c4 = 0.5 * gas.alpha_m3 * COULOMB_CONST * q_c * q_c;
            atoms.push_back(center);
        }
    }

    return atoms;
}

[[maybe_unused]] double find_b_max(
    double b_guess,
    double epsilon_deflection,
    double growth_factor,
    double rel_tol,
    const std::function<double(double)>& deflection_weight
);

[[maybe_unused]] static void write_hdf5(
    const std::string& output,
    const std::string& species_id,
    const std::string& gas,
    const Options& opt,
    const Lj1264Samples& samples,
    const std::vector<unsigned char>* cell_done = nullptr,
    long long completed_cells = -1
) {
    const hsize_t n_orient = static_cast<hsize_t>(samples.orientations_quat.size());
    const hsize_t n_bins = static_cast<hsize_t>(samples.logv_bins.size());

    H5::H5File file(output, H5F_ACC_TRUNC);

    auto write_attr_str = [&](const std::string& name, const std::string& value) {
        H5::StrType str_type(H5::PredType::C_S1, value.size());
        H5::DataSpace space(H5S_SCALAR);
        H5::Attribute attr = file.createAttribute(name, str_type, space);
        attr.write(str_type, value);
    };

    auto write_attr_double = [&](const std::string& name, double value) {
        H5::DataSpace space(H5S_SCALAR);
        H5::Attribute attr = file.createAttribute(name, H5::PredType::NATIVE_DOUBLE, space);
        attr.write(H5::PredType::NATIVE_DOUBLE, &value);
    };

    auto write_attr_int = [&](const std::string& name, long long value) {
        H5::DataSpace space(H5S_SCALAR);
        H5::Attribute attr = file.createAttribute(name, H5::PredType::NATIVE_LLONG, space);
        attr.write(H5::PredType::NATIVE_LLONG, &value);
    };

    write_attr_str("format", ICARION::physics::INTERACTION_POTENTIAL_OFFLINE_SAMPLE_SET_FORMAT);
    write_attr_str("species_id", species_id);
    write_attr_str("gas", gas);
    write_attr_str("gas_model", opt.gas_model);
    write_attr_str("mixing_rule", opt.mixing_rule);
    write_attr_str("polarization", opt.polarization);
    write_attr_str("potential", opt.potential);
    write_attr_str("param_model", opt.param_model);
    write_attr_str("units", ICARION::physics::INTERACTION_POTENTIAL_OFFLINE_SAMPLE_SET_UNITS);
    write_attr_int("version", 1);
    write_attr_int("n_orientations", static_cast<long long>(n_orient));
    write_attr_int("v_bins", static_cast<long long>(n_bins));
    write_attr_int("seed", static_cast<long long>(opt.seed));
    write_attr_double("epsilon_deflection", opt.epsilon_deflection);
    write_attr_double("sigma_scale", opt.sigma_scale);
    write_attr_double("epsilon_scale", opt.epsilon_scale);
    write_attr_double("pol_damp_A", opt.pol_damp_A);
    write_attr_double("mmff_energy_scale", opt.mmff_energy_scale);
    write_attr_double("mmff_distance_scale", opt.mmff_distance_scale);
    if (completed_cells >= 0) {
        write_attr_int("completed_cells", completed_cells);
    }
    if (opt.gas_model != "mono") {
        write_attr_double("n2_bond_m", opt.n2_bond_A * ANGSTROM_TO_M);
        write_attr_int("n2_quadrupole", opt.n2_quadrupole ? 1 : 0);
        write_attr_double("n2_q_site_e", opt.n2_q_site_e);
        write_attr_double("n2_q_center_e", opt.n2_q_center_e);
        write_attr_int("n2_average_lj", opt.n2_average_lj ? 1 : 0);
        write_attr_int("n2_aniso_pol", opt.n2_aniso_pol ? 1 : 0);
        if (opt.n2_aniso_pol) {
            write_attr_double("n2_alpha_par_A3", opt.n2_alpha_par_A3);
            write_attr_double("n2_alpha_perp_A3", opt.n2_alpha_perp_A3);
        }
    }

    // Runtime lookup axis 1: log(v_rel). The handler interpolates/samples by
    // relative speed, so this grid is stored explicitly rather than inferred.
    if (!samples.logv_bins.empty()) {
        hsize_t dims[1] = {n_bins};
        H5::DataSpace space(1, dims);
        H5::DataSet dset = file.createDataSet("logv_bins", H5::PredType::NATIVE_DOUBLE, space);
        dset.write(samples.logv_bins.data(), H5::PredType::NATIVE_DOUBLE);
    }

    // Runtime lookup axis 0: molecular orientation, stored as quaternions for
    // reproducibility and post-hoc inspection.
    if (!samples.orientations_quat.empty()) {
        hsize_t dims[2] = {n_orient, 4};
        H5::DataSpace space(2, dims);
        H5::DataSet dset = file.createDataSet("orientations_quat", H5::PredType::NATIVE_DOUBLE, space);
        std::vector<double> flat;
        flat.reserve(samples.orientations_quat.size() * 4);
        for (const auto& q : samples.orientations_quat) {
            flat.push_back(q[0]);
            flat.push_back(q[1]);
            flat.push_back(q[2]);
            flat.push_back(q[3]);
        }
        dset.write(flat.data(), H5::PredType::NATIVE_DOUBLE);
    }

    // Main runtime rate table: momentum-transfer cross section per
    // orientation/speed cell.
    if (!samples.sigma_mt.empty()) {
        hsize_t dims[2] = {n_orient, n_bins};
        H5::DataSpace space(2, dims);
        H5::DataSet dset = file.createDataSet("sigma_mt_m2", H5::PredType::NATIVE_DOUBLE, space);
        dset.write(samples.sigma_mt.data(), H5::PredType::NATIVE_DOUBLE);
    }

    // Impact-parameter cutoff used to generate each cell. This is not needed
    // for the current runtime rate calculation, but it is critical provenance
    // for auditing and reproducing a table.
    if (!samples.b_max.empty()) {
        hsize_t dims[2] = {n_orient, n_bins};
        H5::DataSpace space(2, dims);
        H5::DataSet dset = file.createDataSet("b_max_m", H5::PredType::NATIVE_DOUBLE, space);
        dset.write(samples.b_max.data(), H5::PredType::NATIVE_DOUBLE);
    }

    // Optional high-fidelity momentum-kick payload. Per-cell CDF rows are
    // flattened into cdf_values/dp_samples, while offsets/counts preserve the
    // (orientation, speed) cell boundaries.
    if (!samples.cdf_offsets.empty()) {
        hsize_t dims[2] = {n_orient, n_bins};
        H5::DataSpace space(2, dims);
        H5::DataSet dset = file.createDataSet("cdf_offsets", H5::PredType::NATIVE_LLONG, space);
        dset.write(samples.cdf_offsets.data(), H5::PredType::NATIVE_LLONG);
    }

    if (!samples.cdf_counts.empty()) {
        hsize_t dims[2] = {n_orient, n_bins};
        H5::DataSpace space(2, dims);
        H5::DataSet dset = file.createDataSet("cdf_counts", H5::PredType::NATIVE_LLONG, space);
        dset.write(samples.cdf_counts.data(), H5::PredType::NATIVE_LLONG);
    }

    if (!samples.cdf_values.empty()) {
        hsize_t dims[1] = {static_cast<hsize_t>(samples.cdf_values.size())};
        H5::DataSpace space(1, dims);
        H5::DataSet dset = file.createDataSet("cdf_values", H5::PredType::NATIVE_DOUBLE, space);
        dset.write(samples.cdf_values.data(), H5::PredType::NATIVE_DOUBLE);
    }

    if (!samples.dp_samples.empty()) {
        hsize_t dims[2] = {static_cast<hsize_t>(samples.dp_samples.size() / 3), 3};
        H5::DataSpace space(2, dims);
        H5::DataSet dset = file.createDataSet("dp_samples", H5::PredType::NATIVE_DOUBLE, space);
        dset.write(samples.dp_samples.data(), H5::PredType::NATIVE_DOUBLE);
    }

    if (!samples.dp_stats.empty()) {
        // Always write the compact fallback moments. Runtime can still sample a
        // statistically consistent kick if the full CDF payload is absent.
        hsize_t dims[3] = {n_orient, n_bins, 4};
        H5::DataSpace space(3, dims);
        H5::DataSet dset = file.createDataSet("dp_stats", H5::PredType::NATIVE_DOUBLE, space);
        dset.write(samples.dp_stats.data(), H5::PredType::NATIVE_DOUBLE);
    }

    if (cell_done && !cell_done->empty()) {
        hsize_t dims[2] = {n_orient, n_bins};
        H5::DataSpace space(2, dims);
        H5::DataSet dset = file.createDataSet("cell_done", H5::PredType::NATIVE_UCHAR, space);
        dset.write(cell_done->data(), H5::PredType::NATIVE_UCHAR);
    }
}

static bool load_hdf5_resume_checkpoint(
    const std::string& input_h5,
    const std::vector<double>& expected_logv,
    Lj1264Samples& samples,
    std::vector<unsigned char>& cell_done,
    size_t& completed_cells,
    std::string& error
) {
    completed_cells = 0;
    if (!std::filesystem::exists(input_h5)) {
        return false;
    }

    try {
        H5::H5File file(input_h5, H5F_ACC_RDONLY);

        const hsize_t n_orient = static_cast<hsize_t>(samples.orientations_quat.size());
        const hsize_t n_bins = static_cast<hsize_t>(samples.logv_bins.size());
        const size_t n_cells = static_cast<size_t>(n_orient * n_bins);

        {
            H5::DataSet logv_dset = file.openDataSet("logv_bins");
            H5::DataSpace logv_space = logv_dset.getSpace();
            if (logv_space.getSimpleExtentNdims() != 1) {
                error = "Resume file has invalid logv_bins rank.";
                return false;
            }
            hsize_t dims[1] = {0};
            logv_space.getSimpleExtentDims(dims);
            if (dims[0] != n_bins) {
                error = "Resume file v_bins mismatch.";
                return false;
            }
            std::vector<double> file_logv(static_cast<size_t>(dims[0]), 0.0);
            logv_dset.read(file_logv.data(), H5::PredType::NATIVE_DOUBLE);
            const double tol = 1e-12;
            for (size_t i = 0; i < file_logv.size(); ++i) {
                if (std::abs(file_logv[i] - expected_logv[i]) > tol) {
                    error = "Resume file logv_bins mismatch.";
                    return false;
                }
            }
        }

        auto read_2d_double = [&](const std::string& name, std::vector<double>& target) {
            H5::DataSet ds = file.openDataSet(name);
            H5::DataSpace sp = ds.getSpace();
            if (sp.getSimpleExtentNdims() != 2) {
                throw std::runtime_error("Dataset '" + name + "' rank mismatch");
            }
            hsize_t dims[2] = {0, 0};
            sp.getSimpleExtentDims(dims);
            if (dims[0] != n_orient || dims[1] != n_bins) {
                throw std::runtime_error("Dataset '" + name + "' shape mismatch");
            }
            ds.read(target.data(), H5::PredType::NATIVE_DOUBLE);
        };

        read_2d_double("sigma_mt_m2", samples.sigma_mt);
        read_2d_double("b_max_m", samples.b_max);

        {
            H5::DataSet ds = file.openDataSet("dp_stats");
            H5::DataSpace sp = ds.getSpace();
            if (sp.getSimpleExtentNdims() != 3) {
                throw std::runtime_error("Dataset 'dp_stats' rank mismatch");
            }
            hsize_t dims[3] = {0, 0, 0};
            sp.getSimpleExtentDims(dims);
            if (dims[0] != n_orient || dims[1] != n_bins || dims[2] != 4) {
                throw std::runtime_error("Dataset 'dp_stats' shape mismatch");
            }
            ds.read(samples.dp_stats.data(), H5::PredType::NATIVE_DOUBLE);
        }

        cell_done.assign(n_cells, 0);
        bool has_cell_done = false;
        try {
            H5::DataSet ds = file.openDataSet("cell_done");
            H5::DataSpace sp = ds.getSpace();
            if (sp.getSimpleExtentNdims() != 2) {
                throw std::runtime_error("Dataset 'cell_done' rank mismatch");
            }
            hsize_t dims[2] = {0, 0};
            sp.getSimpleExtentDims(dims);
            if (dims[0] != n_orient || dims[1] != n_bins) {
                throw std::runtime_error("Dataset 'cell_done' shape mismatch");
            }
            ds.read(cell_done.data(), H5::PredType::NATIVE_UCHAR);
            has_cell_done = true;
        } catch (const H5::Exception&) {
            has_cell_done = false;
        }

        if (has_cell_done) {
            for (unsigned char flag : cell_done) {
                if (flag != 0) {
                    ++completed_cells;
                }
            }
        } else {
            for (size_t idx = 0; idx < n_cells; ++idx) {
                const bool done = (samples.b_max[idx] > 0.0);
                cell_done[idx] = done ? 1 : 0;
                if (done) {
                    ++completed_cells;
                }
            }
        }

        return true;
    } catch (const H5::Exception& e) {
        error = std::string("HDF5 resume read failed: ") + e.getDetailMsg();
        return false;
    } catch (const std::exception& e) {
        error = std::string("Resume read failed: ") + e.what();
        return false;
    }
}

[[maybe_unused]] static void compute_vdw_force_and_potential(
    const std::vector<AtomParam>& atoms,
    const Vec3& r,
    Vec3& force,
    double& potential,
    PotentialType potential_type
) {
    force = Vec3(0.0, 0.0, 0.0);
    potential = 0.0;
    if (potential_type == PotentialType::LJ1264) {
        // 12-6 pair potential. The 4*epsilon*sigma^n factors are precomputed in
        // AtomParam; only inverse powers of r are evaluated here.
        for (const auto& atom : atoms) {
            if (atom.lj_c6 == 0.0 && atom.lj_c12 == 0.0) {
                continue;
            }
            const Vec3 dr = Vec3(r.x - atom.pos.x, r.y - atom.pos.y, r.z - atom.pos.z);
            const double r2 = norm2(dr);
            if (r2 <= 0.0) {
                continue;
            }
            const double rnorm = std::sqrt(r2);
            const double inv_r = 1.0 / rnorm;
            const double inv_r2 = inv_r * inv_r;
            const double inv_r4 = inv_r2 * inv_r2;
            const double inv_r6 = inv_r4 * inv_r2;
            const double inv_r12 = inv_r6 * inv_r6;
            const double inv_r7 = inv_r6 * inv_r;
            const double inv_r13 = inv_r12 * inv_r;

            potential += atom.lj_c12 * inv_r12 - atom.lj_c6 * inv_r6;
            const double dVdr = -12.0 * atom.lj_c12 * inv_r13 + 6.0 * atom.lj_c6 * inv_r7;

            const double f_mag = -dVdr;
            const double scale = f_mag * inv_r;
            force += dr * scale;
        }
        return;
    }

    // Exp-6 pair potential. R is dimensionless r/sigma; prefactors were cached
    // from epsilon so this loop mirrors the runtime trajectory equations.
    for (const auto& atom : atoms) {
        if (atom.inv_sigma <= 0.0 || atom.epsilon == 0.0) {
            continue;
        }
        const Vec3 dr = Vec3(r.x - atom.pos.x, r.y - atom.pos.y, r.z - atom.pos.z);
        const double r2 = norm2(dr);
        if (r2 <= 0.0) {
            continue;
        }
        const double rnorm = std::sqrt(r2);
        const double inv_r = 1.0 / rnorm;
        const double R = rnorm * atom.inv_sigma;
        if (R <= 0.0) {
            continue;
        }
        const double R2 = R * R;
        const double R6 = R2 * R2 * R2;
        const double inv_R6 = 1.0 / R6;
        const double expterm = atom.exp6_pref_exp * std::exp(-EXP6_BETA * R);
        const double dispterm = atom.exp6_pref_disp * inv_R6;
        potential += expterm - dispterm;
        const double dVdr = (-EXP6_BETA * expterm + 6.0 * dispterm / R) * atom.inv_sigma;

        const double f_mag = -dVdr;
        const double scale = f_mag * inv_r;
        force += dr * scale;
    }
}

[[maybe_unused]] static void compute_polarization_force_and_potential(
    const std::vector<AtomParam>& atoms,
    const Vec3& r,
    Vec3& force,
    double& potential,
    PolarizationModel polarization_model,
    double alpha_m3,
    double pol_damp_m,
    const GasModelConfig* gas_model,
    const Vec3& gas_axis
) {
    if (polarization_model == PolarizationModel::NONE) {
        return;
    }

    if (polarization_model == PolarizationModel::PARTIAL ||
        polarization_model == PolarizationModel::TOTAL) {
        // Isotropic induced-dipole attraction. Optional damping replaces r^2 by
        // r^2 + a^2 to avoid an unphysical singularity at very short range.
        const bool damped = pol_damp_m > 0.0;
        const double damp2 = pol_damp_m * pol_damp_m;
        for (const auto& atom : atoms) {
            if (atom.c4 == 0.0) {
                continue;
            }
            const Vec3 dr = Vec3(r.x - atom.pos.x, r.y - atom.pos.y, r.z - atom.pos.z);
            const double r2 = norm2(dr);
            if (r2 <= 0.0) {
                continue;
            }
            const double rnorm = std::sqrt(r2);
            const double inv_r = 1.0 / rnorm;

            if (damped) {
                const double r2d = r2 + damp2;
                const double inv_r2d = 1.0 / r2d;
                const double inv_r4d = inv_r2d * inv_r2d;
                const double inv_r6d = inv_r4d * inv_r2d;
                potential += -atom.c4 * inv_r4d;
                const double dVdr = 4.0 * atom.c4 * rnorm * inv_r6d;
                force += dr * (-dVdr * inv_r);
            } else {
                const double inv_r2 = inv_r * inv_r;
                const double inv_r4 = inv_r2 * inv_r2;
                const double inv_r5 = inv_r4 * inv_r;
                potential += -atom.c4 * inv_r4;
                const double dVdr = 4.0 * atom.c4 * inv_r5;
                force += dr * (-dVdr * inv_r);
            }
        }
        return;
    }

    double rx = 0.0;
    double ry = 0.0;
    double rz = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    double sum4 = 0.0;
    double sum5 = 0.0;
    double sum6 = 0.0;

    // Pairwise polarization first reconstructs the electric field and field
    // gradient from explicit ion charges at the gas probe position.
    for (const auto& atom : atoms) {
        if (atom.q_c == 0.0) {
            continue;
        }
        const Vec3 dr = Vec3(r.x - atom.pos.x, r.y - atom.pos.y, r.z - atom.pos.z);
        const double r2 = norm2(dr);
        if (r2 <= 0.0) {
            continue;
        }
        const double rnorm = std::sqrt(r2);

        const double r3 = r2 * rnorm;
        const double r5 = r3 * r2;
        const double rxyz3i = atom.q_c / r3;
        const double rxyz5i = -3.0 * atom.q_c / r5;
        rx += dr.x * rxyz3i;
        ry += dr.y * rxyz3i;
        rz += dr.z * rxyz3i;
        sum1 += rxyz3i + (dr.x * dr.x) * rxyz5i;
        sum2 += (dr.x * dr.y) * rxyz5i;
        sum3 += (dr.x * dr.z) * rxyz5i;
        sum4 += rxyz3i + (dr.y * dr.y) * rxyz5i;
        sum5 += (dr.y * dr.z) * rxyz5i;
        sum6 += rxyz3i + (dr.z * dr.z) * rxyz5i;
    }

    const double ex = rx;
    const double ey = ry;
    const double ez = rz;
    const double e2 = ex * ex + ey * ey + ez * ez;

    if (gas_model && gas_model->anisotropic_pol) {
        // N2 can use a tensor polarizability: alpha_parallel along the sampled
        // gas axis and alpha_perp in the transverse plane.
        const double alpha_par = gas_model->alpha_par_m3;
        const double alpha_perp = gas_model->alpha_perp_m3;
        const double alpha_delta = alpha_par - alpha_perp;
        const double en = ex * gas_axis.x + ey * gas_axis.y + ez * gas_axis.z;

        const double scale_aniso = 0.5 * COULOMB_CONST;
        potential += -scale_aniso * (alpha_perp * e2 + alpha_delta * en * en);

        force.x += scale_aniso * (2.0 * alpha_perp * (ex * sum1 + ey * sum2 + ez * sum3)
                 + 2.0 * alpha_delta * en * (gas_axis.x * sum1 + gas_axis.y * sum2 + gas_axis.z * sum3));
        force.y += scale_aniso * (2.0 * alpha_perp * (ex * sum2 + ey * sum4 + ez * sum5)
                 + 2.0 * alpha_delta * en * (gas_axis.x * sum2 + gas_axis.y * sum4 + gas_axis.z * sum5));
        force.z += scale_aniso * (2.0 * alpha_perp * (ex * sum3 + ey * sum5 + ez * sum6)
                 + 2.0 * alpha_delta * en * (gas_axis.x * sum3 + gas_axis.y * sum5 + gas_axis.z * sum6));
    } else {
        // Isotropic tensor fallback: U = -1/2 alpha |E|^2.
        const double dipol = 0.5 * alpha_m3 * COULOMB_CONST;
        potential += -dipol * e2;
        force.x += dipol * (2.0 * ex * sum1 + 2.0 * ey * sum2 + 2.0 * ez * sum3);
        force.y += dipol * (2.0 * ex * sum2 + 2.0 * ey * sum4 + 2.0 * ez * sum5);
        force.z += dipol * (2.0 * ex * sum3 + 2.0 * ey * sum5 + 2.0 * ez * sum6);
    }
}

[[maybe_unused]] static void compute_quadrupole_force_and_potential(
    const std::vector<AtomParam>& atoms,
    const Vec3& r,
    const GasModelConfig& gas_model,
    const Vec3& gas_axis,
    Vec3& force,
    double& potential
) {
    if (!gas_model.quadrupole) {
        return;
    }
    // Optional N2 quadrupole: two off-center charges plus a center charge, all
    // evaluated against explicit ion partial charges.
    const Vec3 offset = gas_axis * (0.5 * gas_model.bond_m);
    const Vec3 site_pos[2] = {r + offset, r - offset};
    const double q_site = gas_model.q_site_c;
    const double q_center = gas_model.q_center_c;

    for (const auto& atom : atoms) {
        if (atom.q_c == 0.0) {
            continue;
        }
        const double qi = atom.q_c;
        for (int s = 0; s < 2; ++s) {
            const Vec3 dr = Vec3(site_pos[s].x - atom.pos.x, site_pos[s].y - atom.pos.y, site_pos[s].z - atom.pos.z);
            const double r2 = norm2(dr);
            if (r2 <= 0.0) {
                continue;
            }
            const double inv_r = 1.0 / std::sqrt(r2);
            const double inv_r3 = inv_r * inv_r * inv_r;
            const double pref = COULOMB_CONST * qi * q_site;
            potential += pref * inv_r;
            force += dr * (pref * inv_r3);
        }
        if (q_center != 0.0) {
            const Vec3 drc = Vec3(r.x - atom.pos.x, r.y - atom.pos.y, r.z - atom.pos.z);
            const double r2c = norm2(drc);
            if (r2c > 0.0) {
                const double inv_r = 1.0 / std::sqrt(r2c);
                const double inv_r3 = inv_r * inv_r * inv_r;
                const double pref = COULOMB_CONST * qi * q_center;
                potential += pref * inv_r;
                force += drc * (pref * inv_r3);
            }
        }
    }
}

[[maybe_unused]] static void compute_force_and_potential(
    const std::vector<AtomParam>& atoms,
    const Vec3& r,
    Vec3& force,
    double& potential,
    PotentialType potential_type,
    PolarizationModel polarization_model,
    double alpha_m3,
    double pol_damp_m,
    const GasModelConfig* gas_model,
    const Vec3& gas_axis
) {
    force = Vec3(0.0, 0.0, 0.0);
    potential = 0.0;

    if (gas_model && gas_model->geometry == GasGeometry::DiatomicN2) {
        // For diatomic gas, evaluate the van-der-Waals interaction at both gas
        // sites and average if requested. Polarization/quadrupole terms are then
        // added at the molecular probe position using the sampled gas axis.
        Vec3 f1, f2;
        double v1 = 0.0;
        double v2 = 0.0;
        const Vec3 offset = gas_axis * (0.5 * gas_model->bond_m);
        compute_vdw_force_and_potential(atoms, r + offset, f1, v1, potential_type);
        compute_vdw_force_and_potential(atoms, r - offset, f2, v2, potential_type);
        const double scale = gas_model->average_lj ? 0.5 : 1.0;
        force += (f1 + f2) * scale;
        potential += (v1 + v2) * scale;
        compute_polarization_force_and_potential(atoms, r, force, potential,
                                                 polarization_model, alpha_m3, pol_damp_m,
                                                 gas_model, gas_axis);
        compute_quadrupole_force_and_potential(atoms, r, *gas_model, gas_axis, force, potential);
        return;
    }

    // Monoatomic/default path: one probe point carries all pair and polarization
    // interactions.
    compute_vdw_force_and_potential(atoms, r, force, potential, potential_type);
    compute_polarization_force_and_potential(atoms, r, force, potential,
                                             polarization_model, alpha_m3, pol_damp_m,
                                             nullptr, Vec3(0.0, 0.0, 1.0));
}

[[maybe_unused]] static TrajectoryResult integrate_trajectory(
    const std::vector<AtomParam>& atoms,
    const TrajectoryConfig& cfg,
    const Vec3& r0,
    const Vec3& v0,
    double m_red,
    PotentialType potential_type,
    PolarizationModel polarization_model,
    double alpha_m3,
    double pol_damp_m,
    const GasModelConfig* gas_model,
    const Vec3& gas_axis
) {
    TrajectoryResult result;
    Vec3 r = r0;
    Vec3 v = v0;
    Vec3 force;
    double potential = 0.0;
    compute_force_and_potential(atoms, r, force, potential, potential_type, polarization_model,
                                alpha_m3, pol_damp_m, gas_model, gas_axis);

    const double inv_m = 1.0 / m_red;
    double energy_prev = 0.5 * m_red * norm2(v) + potential;
    double dt = 0.0;
    int good_steps = 0;

    if (cfg.store_path) {
        result.path_points.push_back(r);
    }

    for (int step = 0; step < cfg.max_steps; ++step) {
        const double r_norm = norm(r);
        const double f_norm = norm(force);
        if (f_norm > 0.0 && r_norm > 0.0) {
            const double dt_force = cfg.eta_dt * std::sqrt(m_red * r_norm / f_norm);
            dt = (dt <= 0.0) ? dt_force : std::min(dt, dt_force);
        }
        if (cfg.max_step_m > 0.0) {
            const double v_norm = norm(v);
            if (v_norm > 1e-12) {
                const double dt_step = cfg.max_step_m / v_norm;
                dt = (dt <= 0.0) ? dt_step : std::min(dt, dt_step);
            }
        }
        if (dt <= 0.0) {
            dt = 1e-15;
        }

        Vec3 v_half = v + force * (0.5 * dt * inv_m);
        Vec3 r_trial = r + v_half * dt;
        Vec3 force_new;
        double potential_new = 0.0;
        compute_force_and_potential(atoms, r_trial, force_new, potential_new, potential_type, polarization_model,
                                    alpha_m3, pol_damp_m, gas_model, gas_axis);
        Vec3 v_trial = v_half + force_new * (0.5 * dt * inv_m);

        const double energy_new = 0.5 * m_red * norm2(v_trial) + potential_new;
        const double rel_drift = std::abs(energy_new - energy_prev) / std::max(1e-12, std::abs(energy_prev));
        result.max_rel_energy_drift = std::max(result.max_rel_energy_drift, rel_drift);

        if (rel_drift > cfg.energy_rel_tol_high) {
            dt *= 0.5;
            continue;
        }

        r = r_trial;
        v = v_trial;
        force = force_new;
        potential = potential_new;
        energy_prev = energy_new;

        if (cfg.store_path && (step % std::max(1, cfg.path_stride) == 0)) {
            result.path_points.push_back(r);
        }

        if (rel_drift < cfg.energy_rel_tol_low) {
            ++good_steps;
            if (good_steps >= cfg.relax_steps) {
                dt *= 1.1;
                good_steps = 0;
            }
        } else {
            good_steps = 0;
        }

        const double kin = 0.5 * m_red * norm2(v);
        const double r_norm_new = norm(r);
        const double f_norm_new = norm(force);
        if (r_norm_new > cfg.r_cut_m && std::abs(potential) / std::max(1e-12, kin) < 1e-6 &&
            (f_norm_new * r_norm_new) / std::max(1e-12, kin) < 1e-6) {
            result.reached_asymptotic = true;
            break;
        }
    }

    if (cfg.store_path) {
        if (result.path_points.empty() ||
            norm(result.path_points.back() - r) > 1e-15) {
            result.path_points.push_back(r);
        }
    }

    result.v_out = v;
    return result;
}

[[maybe_unused]] static double deflection_weight_for_b(
    const std::vector<AtomParam>& atoms,
    const TrajectoryConfig& base_cfg,
    double b,
    double v_rel,
    double m_red,
    double z0_margin_m,
    double r_cut_margin_m,
    PotentialType potential_type,
    PolarizationModel polarization_model,
    double alpha_m3,
    const GasModelConfig* gas_model,
    const Vec3& gas_axis
) {
    TrajectoryConfig cfg = base_cfg;
    cfg.r_cut_m = b + r_cut_margin_m;

    const Vec3 r0(b, 0.0, -(b + z0_margin_m));
    const Vec3 v0(0.0, 0.0, v_rel);
    const TrajectoryResult result =
        integrate_trajectory(atoms, cfg, r0, v0, m_red, potential_type, polarization_model,
                             alpha_m3, cfg.pol_damp_m, gas_model, gas_axis);

    const double v2 = norm2(v0);
    if (v2 <= 0.0) {
        return 0.0;
    }
    const double cos_theta = dot_local(v0, result.v_out) / std::max(1e-12, v2);
    const double clamped = std::max(-1.0, std::min(1.0, cos_theta));
    return 1.0 - clamped;
}

[[maybe_unused]] static double determine_b_max(
    const std::vector<AtomParam>& atoms,
    const TrajectoryConfig& base_cfg,
    double v_rel,
    double m_red,
    double b_guess,
    double epsilon_deflection,
    double growth_factor,
    double rel_tol,
    PotentialType potential_type,
    PolarizationModel polarization_model,
    double alpha_m3,
    const GasModelConfig* gas_model,
    const Vec3& gas_axis
) {
    const double z0_margin_m = 20.0 * ANGSTROM_TO_M;
    const double r_cut_margin_m = 50.0 * ANGSTROM_TO_M;
    auto weight = [&](double b) {
        return deflection_weight_for_b(
            atoms, base_cfg, b, v_rel, m_red, z0_margin_m, r_cut_margin_m,
            potential_type, polarization_model, alpha_m3, gas_model, gas_axis);
    };
    return find_b_max(b_guess, epsilon_deflection, growth_factor, rel_tol, weight);
}

[[maybe_unused]] static double determine_b_max_diatomic_robust(
    const std::vector<AtomParam>& atoms,
    const TrajectoryConfig& base_cfg,
    double v_rel,
    double m_red,
    double b_guess,
    double epsilon_deflection,
    double growth_factor,
    double rel_tol,
    PotentialType potential_type,
    PolarizationModel polarization_model,
    double alpha_m3,
    const GasModelConfig* gas_model,
    std::mt19937_64& rng
) {
    const std::array<Vec3, 6> canonical_axes = {
        Vec3(1.0, 0.0, 0.0),
        Vec3(-1.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.0, -1.0, 0.0),
        Vec3(0.0, 0.0, 1.0),
        Vec3(0.0, 0.0, -1.0)
    };

    double b_max = 0.0;
    for (const Vec3& axis : canonical_axes) {
        b_max = std::max(
            b_max,
            determine_b_max(
                atoms,
                base_cfg,
                v_rel,
                m_red,
                b_guess,
                epsilon_deflection,
                growth_factor,
                rel_tol,
                potential_type,
                polarization_model,
                alpha_m3,
                gas_model,
                axis
            )
        );
    }

    constexpr int EXTRA_RANDOM_AXES = 4;
    for (int i = 0; i < EXTRA_RANDOM_AXES; ++i) {
        const Vec3 axis = sample_unit_vector(rng);
        b_max = std::max(
            b_max,
            determine_b_max(
                atoms,
                base_cfg,
                v_rel,
                m_red,
                b_guess,
                epsilon_deflection,
                growth_factor,
                rel_tol,
                potential_type,
                polarization_model,
                alpha_m3,
                gas_model,
                axis
            )
        );
    }

    return b_max;
}

[[maybe_unused]] double find_b_max(
    double b_guess,
    double epsilon_deflection,
    double growth_factor,
    double rel_tol,
    const std::function<double(double)>& deflection_weight
) {
    auto w = [&](double b) { return deflection_weight(b); };
    double b_low = b_guess;
    double w_low = w(b_low);
    if (w_low < epsilon_deflection) {
        return b_low;
    }

    double b_high = b_low;
    double w_high = w_low;
    int stable_hits = 0;
    while (stable_hits < 2) {
        b_high *= growth_factor;
        w_high = w(b_high);
        if (w_high < epsilon_deflection) {
            ++stable_hits;
        } else {
            stable_hits = 0;
        }
        if (b_high > 1e-6) {
            break;
        }
    }

    if (w_high >= epsilon_deflection) {
        return b_high;
    }

    for (int iter = 0; iter < 64; ++iter) {
        const double b_mid = 0.5 * (b_low + b_high);
        const double w_mid = w(b_mid);
        if (w_mid < epsilon_deflection) {
            b_high = b_mid;
        } else {
            b_low = b_mid;
        }
        if (std::abs(b_high - b_low) / b_high < rel_tol) {
            break;
        }
    }
    return b_high;
}

}  // namespace

int main(int argc, char** argv) {
    Options opt;
    if (!parse_args(argc, argv, opt)) {
        print_usage();
        return EXIT_FAILURE;
    }
    if (opt.gas_model == "n2") {
        opt.gas_model = "diatomic";
    }
    if (opt.gas_model != "mono" && opt.gas_model != "diatomic") {
        std::cerr << "Unsupported --gas-model '" << opt.gas_model << "' (use mono or diatomic).\\n";
        return EXIT_FAILURE;
    }

    int resolved_threads = 1;
#ifdef _OPENMP
    if (opt.threads > 0) {
        omp_set_num_threads(opt.threads);
    }
    resolved_threads = omp_get_max_threads();
#else
    if (opt.threads > 0) {
        std::cerr << "[interaction_potential_precompute] Warning: --threads was set but OpenMP is not enabled in this build; ignoring.\\n";
    }
#endif

    std::cout << "[interaction_potential_precompute] Inputs:\\n";
    std::cout << "  input:      " << opt.input << "\\n";
    std::cout << "  output:     " << opt.output << "\\n";
    std::cout << "  species:    " << opt.species_id << "\\n";
    std::cout << "  gas:        " << opt.gas << "\\n";
    std::cout << "  gas_model:  " << opt.gas_model << "\\n";
    std::cout << "  mixing:     " << opt.mixing_rule << "\\n";
    std::cout << "  polarization: " << opt.polarization << "\\n";
    std::cout << "  potential:  " << opt.potential << "\\n";
    std::cout << "  format:     " << opt.format << "\\n";
    std::cout << "  orient grid:" << opt.orient_grid << " (" << opt.n_orientations << ")\\n";
    std::cout << "  trials:     " << opt.n_trials << "\\n";
    std::cout << "  threads:    " << resolved_threads;
#ifdef _OPENMP
    if (opt.threads == 0) {
        std::cout << " (OpenMP default)";
    }
#else
    std::cout << " (OpenMP disabled)";
#endif
    std::cout << "\\n";
    std::cout << "  v bins:     " << opt.v_bins << "\\n";
    std::cout << "  epsilon:    " << opt.epsilon_deflection << "\\n";
    std::cout << "  b_guess:    " << (opt.b_guess_m / ANGSTROM_TO_M) << " Å\\n";
    std::cout << "  b_growth:   " << opt.b_growth << "\\n";
    std::cout << "  b_rel_tol:  " << opt.b_rel_tol << "\\n";
    std::cout << "  eta_dt:     " << opt.eta_dt << "\\n";
    std::cout << "  sigma_scale:" << opt.sigma_scale << "\\n";
    std::cout << "  eps_scale:  " << opt.epsilon_scale << "\\n";
    std::cout << "  pol_damp_A: " << opt.pol_damp_A << "\\n";
    std::cout << "  full CDF:   " << (opt.store_full_cdf ? "yes" : "no") << "\\n";
    std::cout << "  resume:     " << (opt.resume ? "yes" : "no") << "\\n";
    std::cout << "  checkpoint: " << opt.checkpoint_cells << " cells\\n";
    if (opt.seed == 0) {
        opt.seed = default_seed_for(opt.species_id, opt.gas);
    }
    std::cout << "  seed:       " << opt.seed << "\\n";
    if (opt.gas_model != "mono") {
        std::cout << "  n2_bond_A:  " << opt.n2_bond_A << "\\n";
        std::cout << "  n2_quadrupole: " << (opt.n2_quadrupole ? "on" : "off") << "\\n";
        std::cout << "  n2_q_site_e: " << opt.n2_q_site_e << "\\n";
        std::cout << "  n2_q_center_e: " << opt.n2_q_center_e << "\\n";
        std::cout << "  n2_average_lj: " << (opt.n2_average_lj ? "yes" : "no") << "\\n";
        std::cout << "  n2_aniso_pol: " << (opt.n2_aniso_pol ? "on" : "off") << "\\n";
        if (opt.n2_aniso_pol) {
            if (opt.n2_alpha_par_A3 > 0.0) {
                std::cout << "  n2_alpha_par_A3: " << opt.n2_alpha_par_A3 << "\\n";
            }
            if (opt.n2_alpha_perp_A3 > 0.0) {
                std::cout << "  n2_alpha_perp_A3: " << opt.n2_alpha_perp_A3 << "\\n";
            }
        }
    }
    if (opt.scan_deflection) {
        std::cout << "  scan:       deflection (v_rel=" << opt.scan_v_rel
                  << " m/s, b_min=" << (opt.scan_b_min_m / ANGSTROM_TO_M)
                  << " Å, b_max=" << (opt.scan_b_max_m / ANGSTROM_TO_M)
                  << " Å, steps=" << opt.scan_b_steps
                  << ", orient=" << opt.scan_orient_index << ")\\n";
    }

    std::mt19937_64 rng(opt.seed);
    std::vector<std::array<double, 4>> orientations;
    if (opt.orient_grid == "random") {
        orientations = generate_orientations_random(opt.n_orientations, rng);
    } else if (opt.orient_grid == "lebedev") {
        std::string lebedev_path = opt.lebedev_file;
        // Allow --lebedev-file to accept a point count (e.g. 110) as well as a file path
        if (!lebedev_path.empty() &&
            lebedev_path.find_first_not_of("0123456789") == std::string::npos) {
            lebedev_path = lebedev_file_for_n(std::stoi(lebedev_path));
        }
        if (lebedev_path.empty()) {
            lebedev_path = lebedev_file_for_n(opt.n_orientations);
        }
        if (lebedev_path.empty()) {
            throw std::runtime_error(
                "Lebedev grid size not supported. Use --lebedev-file or one of: "
                "50, 110, 194, 302, 590"
            );
        }
        orientations = load_lebedev_orientations(lebedev_path);
        if (opt.n_orientations != static_cast<int>(orientations.size())) {
            std::cerr << "[interaction_potential_precompute] Lebedev grid size (" << orientations.size()
                      << ") overrides --n-orientations.\\n";
        }
    } else {
        orientations = generate_orientations_qmc(opt.n_orientations);
    }
    std::cout << "  orientations generated: " << orientations.size() << "\\n";

    Json::Value root;
    {
        std::ifstream ifs(opt.input);
        if (!ifs) {
            std::cerr << "Failed to open input JSON: " << opt.input << "\\n";
            return EXIT_FAILURE;
        }
        Json::CharReaderBuilder builder;
        std::string errs;
        if (!Json::parseFromStream(builder, ifs, &root, &errs)) {
            std::cerr << "Failed to parse input JSON: " << errs << "\\n";
            return EXIT_FAILURE;
        }
    }

    if (!root.isMember("species") || !root["species"].isObject()) {
        std::cerr << "Input JSON must contain a 'species' object\\n";
        return EXIT_FAILURE;
    }
    const Json::Value& species = root["species"];
    if (!species.isMember(opt.species_id)) {
        std::cerr << "Species '" << opt.species_id << "' not found in input\\n";
        return EXIT_FAILURE;
    }
    const Json::Value& props = species[opt.species_id];
    if (!props.isMember("geometry_file") || !props["geometry_file"].isString()) {
        std::cerr << "Species '" << opt.species_id << "' has no geometry_file\\n";
        return EXIT_FAILURE;
    }
    if (!props.isMember("charge") || !props["charge"].isNumeric()) {
        std::cerr << "Species '" << opt.species_id << "' has no charge\\n";
        return EXIT_FAILURE;
    }

    std::filesystem::path base = std::filesystem::path(opt.input).parent_path();
    std::filesystem::path geom_path = base / props["geometry_file"].asString();
    if (!std::filesystem::exists(geom_path)) {
        std::cerr << "Geometry file not found: " << geom_path << "\\n";
        return EXIT_FAILURE;
    }

    ICARION::io::Molecule molecule;
    try {
        molecule = ICARION::io::load_molecule(geom_path.string());
    } catch (const std::exception& e) {
        std::cerr << "Failed to load geometry: " << e.what() << "\\n";
        return EXIT_FAILURE;
    }

    if (opt.element_params_file.empty()) {
        opt.element_params_file = opt.gas_params_file;
    }

    GasParams gas_params;
    try {
        gas_params = gas_params_or_throw(opt.gas);
        if (gas_params.sigma_m <= 0.0 || gas_params.epsilon_J <= 0.0) {
            const auto file_params = load_gas_lj_params_or_throw(opt.gas_params_file, opt.gas);
            gas_params.sigma_m = file_params.sigma_m;
            gas_params.epsilon_J = file_params.epsilon_J;
            if (file_params.alpha_m3 > 0.0) {
                gas_params.alpha_m3 = file_params.alpha_m3;
            }
            if (file_params.alpha_par_m3 > 0.0) {
                gas_params.alpha_par_m3 = file_params.alpha_par_m3;
            }
            if (file_params.alpha_perp_m3 > 0.0) {
                gas_params.alpha_perp_m3 = file_params.alpha_perp_m3;
            }
        }
        // Explicit CLI overrides take precedence over built-in/database values.
        if (opt.gas_sigma_m > 0.0) {
            gas_params.sigma_m = opt.gas_sigma_m;
        }
        if (opt.gas_epsilon_J > 0.0) {
            gas_params.epsilon_J = opt.gas_epsilon_J;
        }
        if (gas_params.sigma_m <= 0.0 || gas_params.epsilon_J <= 0.0) {
            throw std::runtime_error("Missing gas LJ parameters (sigma/epsilon)");
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\\n";
        std::cerr << "Provide --gas-sigma-A/--gas-epsilon-eV or fill " << opt.gas_params_file << "\\n";
        return EXIT_FAILURE;
    }

    std::unordered_map<std::string, ElementParams> element_params;
    try {
        element_params = load_element_params_or_throw(opt.element_params_file);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\\n";
        return EXIT_FAILURE;
    }

    MixingRule mixing_rule;
    try {
        mixing_rule = mixing_rule_or_throw(opt.mixing_rule);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\\n";
        return EXIT_FAILURE;
    }

    PolarizationModel polarization_model;
    try {
        polarization_model = polarization_model_or_throw(opt.polarization);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\\n";
        return EXIT_FAILURE;
    }

    PotentialType potential_type;
    try {
        potential_type = potential_type_or_throw(opt.potential);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\\n";
        return EXIT_FAILURE;
    }

    ParamModel param_model;
    try {
        param_model = param_model_or_throw(opt.param_model);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\\n";
        return EXIT_FAILURE;
    }
    if (param_model == ParamModel::MMFF_REIJ && potential_type != PotentialType::EXP6) {
        std::cerr << "MMFF R*/e_ij parameterization requires exp6 potential.\\n";
        return EXIT_FAILURE;
    }

    GasModelConfig gas_model;
    const GasModelConfig* gas_model_ptr = nullptr;
    if (opt.gas_model != "mono") {
        if (opt.gas != "N2") {
            std::cerr << "Diatomic gas model is only supported for N2.\\n";
            return EXIT_FAILURE;
        }
        gas_model.geometry = GasGeometry::DiatomicN2;
        gas_model.bond_m = opt.n2_bond_A * ANGSTROM_TO_M;
        if (gas_model.bond_m <= 0.0) {
            std::cerr << "Invalid N2 bond length.\\n";
            return EXIT_FAILURE;
        }
        gas_model.quadrupole = opt.n2_quadrupole;
        gas_model.q_site_c = opt.n2_q_site_e * ELEM_CHARGE_C;
        gas_model.q_center_c = opt.n2_q_center_e * ELEM_CHARGE_C;
        gas_model.average_lj = opt.n2_average_lj;
        if (opt.n2_aniso_pol) {
            if (polarization_model != PolarizationModel::PAIRWISE) {
                std::cerr << "Warning: --n2-aniso-pol requires polarization=pairwise for effect.\\n";
            }
            gas_model.anisotropic_pol = true;
            const double alpha_par_m3 = (opt.n2_alpha_par_A3 > 0.0)
                ? opt.n2_alpha_par_A3 * ANGSTROM3_TO_M3
                : (gas_params.alpha_par_m3 > 0.0 ? gas_params.alpha_par_m3 : gas_params.alpha_m3);
            const double alpha_perp_m3 = (opt.n2_alpha_perp_A3 > 0.0)
                ? opt.n2_alpha_perp_A3 * ANGSTROM3_TO_M3
                : (gas_params.alpha_perp_m3 > 0.0 ? gas_params.alpha_perp_m3 : gas_params.alpha_m3);
            gas_model.alpha_par_m3 = alpha_par_m3;
            gas_model.alpha_perp_m3 = alpha_perp_m3;
            if (opt.n2_alpha_par_A3 <= 0.0) {
                opt.n2_alpha_par_A3 = alpha_par_m3 / ANGSTROM3_TO_M3;
            }
            if (opt.n2_alpha_perp_A3 <= 0.0) {
                opt.n2_alpha_perp_A3 = alpha_perp_m3 / ANGSTROM3_TO_M3;
            }
        }
        gas_model_ptr = &gas_model;
    }

    std::vector<AtomParam> atoms;
    try {
        atoms = build_atom_params_or_throw(
            molecule,
            opt.gas,
            gas_params,
            props["charge"].asDouble(),
            opt.sigma_scale,
            opt.epsilon_scale,
            opt.mmff_energy_scale,
            opt.mmff_distance_scale,
            element_params,
            mixing_rule,
            polarization_model,
            param_model
        );
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\\n";
        return EXIT_FAILURE;
    }
    std::cout << "  atoms loaded: " << atoms.size() << "\\n";
    double min_sigma = 0.0;
    for (const auto& atom : atoms) {
        if (atom.sigma > 0.0) {
            min_sigma = (min_sigma <= 0.0) ? atom.sigma : std::min(min_sigma, atom.sigma);
        }
    }

    double m_ion = 0.0;
    if (props.isMember("mass_amu") && props["mass_amu"].isNumeric()) {
        m_ion = props["mass_amu"].asDouble() * AMU_TO_KG;
    } else if (molecule.total_mass_u > 0.0) {
        m_ion = molecule.total_mass_u * AMU_TO_KG;
    } else {
        std::cerr << "Species mass not found (mass_amu)\\n";
        return EXIT_FAILURE;
    }

    double m_gas = 0.0;
    try {
        m_gas = gas_mass_kg_or_throw(opt.gas);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\\n";
        return EXIT_FAILURE;
    }
    const double m_red = (m_ion * m_gas) / (m_ion + m_gas);

    double temperature_K = opt.temperature_K;
    if (temperature_K <= 0.0 && props.isMember("reference_temperature_K")) {
        temperature_K = props["reference_temperature_K"].asDouble();
    }
    if (temperature_K <= 0.0) {
        temperature_K = 300.0;
    }

    const double v_th = std::sqrt(8.0 * BOLTZMANN_CONSTANT * temperature_K / (M_PI * m_red));
    double v_min = (opt.v_min > 0.0) ? opt.v_min : 0.2 * v_th;
    double v_max = (opt.v_max > 0.0) ? opt.v_max : 5.0 * v_th;
    if (v_min <= 0.0 || v_max <= v_min) {
        std::cerr << "Invalid v_min/v_max for log grid\\n";
        return EXIT_FAILURE;
    }

    std::vector<double> logv_bins;
    logv_bins.reserve(static_cast<size_t>(opt.v_bins));
    const double logv_min = std::log(v_min);
    const double logv_max = std::log(v_max);
    // Use a log-spaced velocity grid because low relative speeds dominate many
    // thermal regimes while high-speed tails still need bounded coverage.
    for (int i = 0; i < opt.v_bins; ++i) {
        const double t = (opt.v_bins == 1) ? 0.0 : static_cast<double>(i) / (opt.v_bins - 1);
        logv_bins.push_back(logv_min + t * (logv_max - logv_min));
    }
    std::cout << "  v_rel range: [" << v_min << ", " << v_max << "] m/s\\n";

    Lj1264Samples samples;
    samples.orientations_quat = orientations;
    samples.logv_bins = logv_bins;
    const size_t n_orient = orientations.size();
    const size_t n_bins = logv_bins.size();
    samples.sigma_mt.assign(n_orient * n_bins, 0.0);
    samples.b_max.assign(n_orient * n_bins, 0.0);
    samples.cdf_offsets.assign(n_orient * n_bins, -1);
    samples.cdf_counts.assign(n_orient * n_bins, 0);
    samples.dp_stats.assign(n_orient * n_bins * 4, 0.0);

    const double z0_margin = 20.0 * ANGSTROM_TO_M;
    const double rcut_margin = 50.0 * ANGSTROM_TO_M;
    TrajectoryConfig base_traj_cfg;
    base_traj_cfg.eta_dt = opt.eta_dt;
    base_traj_cfg.pol_damp_m = opt.pol_damp_A * ANGSTROM_TO_M;
    if (min_sigma > 0.0) {
        base_traj_cfg.max_step_m = 0.2 * min_sigma;
    }

    std::vector<std::vector<AtomParam>> atoms_by_orientation(n_orient);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (size_t oi = 0; oi < n_orient; ++oi) {
        // Pre-rotate the ion once per orientation. Cell workers then only vary
        // v_rel, impact parameter, azimuth, and optional gas-axis orientation.
        auto rotated = atoms;
        for (auto& atom : rotated) {
            atom.pos = rotate_vec_by_quat(orientations[oi], atom.pos);
        }
        atoms_by_orientation[oi] = std::move(rotated);
    }

    if (opt.scan_deflection) {
        if (orientations.empty()) {
            std::cerr << "No orientations available for scan\\n";
            return EXIT_FAILURE;
        }
        if (opt.scan_orient_index < 0 || opt.scan_orient_index >= static_cast<int>(orientations.size())) {
            std::cerr << "scan-orient-idx out of range (0.." << (orientations.size() - 1) << ")\\n";
            return EXIT_FAILURE;
        }
        const std::vector<AtomParam>& atoms_rot = atoms_by_orientation[static_cast<size_t>(opt.scan_orient_index)];

        const double v_rel = opt.scan_v_rel;
        const Vec3 gas_axis_bmax(0.0, 0.0, 1.0);
        double b_max = opt.scan_b_max_m;
        if (b_max <= 0.0) {
            b_max = determine_b_max(atoms_rot, base_traj_cfg, v_rel, m_red,
                                    opt.b_guess_m, opt.epsilon_deflection,
                                    opt.b_growth, opt.b_rel_tol, potential_type,
                                    polarization_model, gas_params.alpha_m3,
                                    gas_model_ptr, gas_axis_bmax);
        }
        const double b_min = std::max(0.0, opt.scan_b_min_m);
        if (b_max <= b_min) {
            std::cerr << "scan b_max must be > b_min\\n";
            return EXIT_FAILURE;
        }

        std::ofstream ofs(opt.output);
        if (!ofs) {
            std::cerr << "Failed to open scan output: " << opt.output << "\\n";
            return EXIT_FAILURE;
        }
        std::ofstream path_ofs;
        if (!opt.scan_paths_output.empty()) {
            path_ofs.open(opt.scan_paths_output);
            if (!path_ofs) {
                std::cerr << "Failed to open scan paths output: " << opt.scan_paths_output << "\\n";
                return EXIT_FAILURE;
            }
            path_ofs << "# species=" << opt.species_id << "\n";
            path_ofs << "# gas=" << opt.gas << "\n";
            path_ofs << "# v_rel_mps=" << v_rel << "\n";
            path_ofs << "# orient_idx=" << opt.scan_orient_index << "\n";
            path_ofs << "# path_stride=" << opt.scan_path_stride << "\n";
            path_ofs << "traj_idx,b_A,point_idx,z_A,x_A,y_A\n";
        }
        ofs << "# species=" << opt.species_id << "\n";
        ofs << "# gas=" << opt.gas << "\n";
        ofs << "# gas_model=" << opt.gas_model << "\n";
        ofs << "# potential=" << opt.potential << "\n";
        ofs << "# mixing=" << opt.mixing_rule << "\n";
        ofs << "# polarization=" << opt.polarization << "\n";
        ofs << "# orient_idx=" << opt.scan_orient_index << "\n";
        if (opt.gas_model != "mono") {
            ofs << "# n2_bond_A=" << opt.n2_bond_A << "\n";
            ofs << "# n2_quadrupole=" << (opt.n2_quadrupole ? 1 : 0) << "\n";
            ofs << "# n2_q_site_e=" << opt.n2_q_site_e << "\n";
            ofs << "# n2_q_center_e=" << opt.n2_q_center_e << "\n";
            ofs << "# n2_average_lj=" << (opt.n2_average_lj ? 1 : 0) << "\n";
            ofs << "# n2_aniso_pol=" << (opt.n2_aniso_pol ? 1 : 0) << "\n";
            if (opt.n2_aniso_pol) {
                ofs << "# n2_alpha_par_A3=" << opt.n2_alpha_par_A3 << "\n";
                ofs << "# n2_alpha_perp_A3=" << opt.n2_alpha_perp_A3 << "\n";
            }
            ofs << "# n2_axis=[0,0,1]\n";
        }
        ofs << "# v_rel_mps=" << v_rel << "\n";
        ofs << "# b_min_A=" << (b_min / ANGSTROM_TO_M) << "\n";
        ofs << "# b_max_A=" << (b_max / ANGSTROM_TO_M) << "\n";
        ofs << "b_A,theta_rad,theta_deg,one_minus_cos,max_rel_energy_drift,reached_asymptotic,phi_xz_rad,phi_xz_deg,vout_x_norm,vout_z_norm\n";

        for (int i = 0; i < opt.scan_b_steps; ++i) {
            const double t = (opt.scan_b_steps == 1) ? 0.0 : static_cast<double>(i) / (opt.scan_b_steps - 1);
            const double b = b_min + t * (b_max - b_min);
            TrajectoryConfig cfg = base_traj_cfg;
            cfg.r_cut_m = b + rcut_margin;
            cfg.store_path = !opt.scan_paths_output.empty();
            cfg.path_stride = std::max(1, opt.scan_path_stride);

            const Vec3 r0(b, 0.0, -(b + z0_margin));
            const Vec3 v0(0.0, 0.0, v_rel);
            const TrajectoryResult res =
                integrate_trajectory(atoms_rot, cfg, r0, v0, m_red, potential_type,
                                    polarization_model, gas_params.alpha_m3, cfg.pol_damp_m,
                                    gas_model_ptr, gas_axis_bmax);

            const double v2 = std::max(1e-12, norm2(v0));
            const double cos_theta = dot_local(v0, res.v_out) / v2;
            const double clamped = std::max(-1.0, std::min(1.0, cos_theta));
            const double theta = std::acos(clamped);
            const double one_minus_cos = 1.0 - clamped;
            const double vout_norm = std::sqrt(std::max(1e-24, norm2(res.v_out)));
            const double vout_x_norm = res.v_out.x / vout_norm;
            const double vout_z_norm = res.v_out.z / vout_norm;
            const double phi_xz = std::atan2(vout_x_norm, vout_z_norm);

            ofs << std::setprecision(10)
                << (b / ANGSTROM_TO_M) << "," << theta << "," << (theta * 180.0 / M_PI) << ","
                << one_minus_cos << "," << res.max_rel_energy_drift << "," << (res.reached_asymptotic ? 1 : 0)
                << "," << phi_xz << "," << (phi_xz * 180.0 / M_PI) << "," << vout_x_norm << "," << vout_z_norm
                << "\n";

            if (path_ofs) {
                for (size_t pi = 0; pi < res.path_points.size(); ++pi) {
                    const Vec3& rp = res.path_points[pi];
                    path_ofs << i << "," << (b / ANGSTROM_TO_M) << "," << pi << ","
                             << (rp.z / ANGSTROM_TO_M) << "," << (rp.x / ANGSTROM_TO_M) << ","
                             << (rp.y / ANGSTROM_TO_M) << "\n";
                }
            }
        }
        std::cout << "  deflection scan written: " << opt.output << "\\n";
        if (path_ofs) {
            std::cout << "  trajectory paths written: " << opt.scan_paths_output << "\\n";
        }
        return EXIT_SUCCESS;
    }

    std::vector<LocalCdfSamples> local_cdfs;
    if (opt.store_full_cdf) {
        if (opt.resume || opt.checkpoint_cells > 0) {
            std::cerr << "Checkpoint/resume is currently supported only without --store-full-cdf.\\n";
            return EXIT_FAILURE;
        }
        local_cdfs.resize(n_orient * n_bins);
    }

    auto compute_cell = [&](size_t oi, size_t ki, LocalCdfSamples* cdf_out) {
        const std::vector<AtomParam>& atoms_rot = atoms_by_orientation[oi];

        // One runtime lookup cell: fixed molecular orientation and fixed relative
        // speed bin. All trajectories in this cell share b_max and contribute to
        // the momentum-transfer cross section used by the collision handler.
        std::mt19937_64 local_rng(seed_for_cell(opt.seed, static_cast<uint64_t>(oi), static_cast<uint64_t>(ki)));
        std::uniform_real_distribution<double> uni(0.0, 1.0);
        const double v_rel = std::exp(logv_bins[ki]);
        const Vec3 v0(0.0, 0.0, v_rel);
        const double inv_v2 = 1.0 / std::max(1e-12, v_rel * v_rel);
        const bool sample_gas_axis = (gas_model_ptr && gas_model_ptr->geometry == GasGeometry::DiatomicN2);
        const Vec3 gas_axis_bmax(0.0, 0.0, 1.0);
        double b_max = 0.0;
        TrajectoryConfig traj_cfg = base_traj_cfg;

        double sum_w = 0.0;
        double sum_par = 0.0;
        double sum_par2 = 0.0;
        double sum_perp = 0.0;
        double sum_perp2 = 0.0;

        std::vector<double> local_weights;
        std::vector<double> local_dp;
        if (cdf_out) {
            local_weights.reserve(static_cast<size_t>(opt.n_trials));
            local_dp.reserve(static_cast<size_t>(opt.n_trials * 3));
        }

        size_t attempted_trials = 0;
        size_t accepted_trials = 0;
        size_t rejected_non_asymptotic = 0;
        const size_t target_trials = static_cast<size_t>(opt.n_trials);
        const size_t max_attempts = std::max<size_t>(target_trials, 1) * 5;

        auto accumulate_result = [&](const TrajectoryResult& res) {
            attempted_trials += 1;
            if (!res.reached_asymptotic) {
                rejected_non_asymptotic += 1;
                return;
            }
            accepted_trials += 1;

            // Momentum-transfer weighting: sigma_mt integrates (1 - cos theta)
            // over impact-parameter space. The same weight also defines the
            // optional runtime CDF over sampled momentum kicks.
            const double cos_theta = dot_local(v0, res.v_out) * inv_v2;
            const double clamped = std::max(-1.0, std::min(1.0, cos_theta));
            const double w = 1.0 - clamped;
            sum_w += w;

            const Vec3 dp_rel = (res.v_out - v0) * m_red;
            const Vec3 dp_ion = dp_rel;
            const double dp_par = dp_ion.z;
            const double dp_perp = std::sqrt(dp_ion.x * dp_ion.x + dp_ion.y * dp_ion.y);
            sum_par += w * dp_par;
            sum_par2 += w * dp_par * dp_par;
            sum_perp += w * dp_perp;
            sum_perp2 += w * dp_perp * dp_perp;

            if (cdf_out && w > 0.0) {
                local_weights.push_back(w);
                local_dp.push_back(dp_ion.x);
                local_dp.push_back(dp_ion.y);
                local_dp.push_back(dp_ion.z);
            }
        };

        if (sample_gas_axis) {
            // Diatomic gas orientation changes the potential, so use a robust
            // b_max search that samples gas-axis orientations while expanding
            // the impact-parameter cutoff.
            b_max = determine_b_max_diatomic_robust(
                atoms_rot,
                base_traj_cfg,
                v_rel,
                m_red,
                opt.b_guess_m,
                opt.epsilon_deflection,
                opt.b_growth,
                opt.b_rel_tol,
                potential_type,
                polarization_model,
                gas_params.alpha_m3,
                gas_model_ptr,
                local_rng
            );
        } else {
            // Monoatomic gas has no orientation degree of freedom; a single
            // asymptotic-deflection search is enough for this cell.
            b_max = determine_b_max(atoms_rot, base_traj_cfg, v_rel, m_red,
                                    opt.b_guess_m, opt.epsilon_deflection,
                                    opt.b_growth, opt.b_rel_tol, potential_type,
                                    polarization_model, gas_params.alpha_m3,
                                    gas_model_ptr, gas_axis_bmax);
        }
        samples.b_max[oi * n_bins + ki] = b_max;
        traj_cfg.r_cut_m = b_max + rcut_margin;

        while (accepted_trials < target_trials && attempted_trials < max_attempts) {
            const size_t remaining = target_trials - accepted_trials;
            const size_t trials_this_round = remaining;
            for (size_t t = 0; t < trials_this_round; ++t) {
                // Uniform disk sampling: sqrt(u) gives a uniform distribution
                // in area over the impact-parameter disk.
                const double u1 = uni(local_rng);
                const double u2 = uni(local_rng);
                const double b = b_max * std::sqrt(u1);
                const double phi = 2.0 * M_PI * u2;
                const Vec3 r0(b * std::cos(phi), b * std::sin(phi), -(b_max + z0_margin));
                const Vec3 gas_axis = sample_gas_axis
                    ? sample_unit_vector(local_rng)
                    : Vec3(0.0, 0.0, 1.0);
                const auto res = integrate_trajectory(atoms_rot, traj_cfg, r0, v0, m_red, potential_type,
                                                      polarization_model, gas_params.alpha_m3, traj_cfg.pol_damp_m,
                                                      gas_model_ptr, gas_axis);
                accumulate_result(res);
            }
        }

        if (accepted_trials == 0) {
            throw std::runtime_error(
                "No asymptotic trajectories converged for cell (orientation="
                + std::to_string(oi) + ", v-bin=" + std::to_string(ki) + ")"
            );
        }
        if (accepted_trials < target_trials) {
#ifdef _OPENMP
#pragma omp critical(ipm_asymptotic_quality_log)
#endif
            {
                std::cerr << "[interaction_potential_precompute] Warning: asymptotic quality gate accepted "
                          << accepted_trials << "/" << target_trials
                          << " trajectories for cell (oi=" << oi << ", ki=" << ki << ")"
                          << " after " << attempted_trials << " attempts"
                          << " (rejected_non_asymptotic=" << rejected_non_asymptotic << ")\n";
            }
        }

        // Convert the Monte Carlo mean over the sampled disk into the physical
        // momentum-transfer cross section for this orientation/speed cell.
        const double sigma_mt = (M_PI * b_max * b_max / static_cast<double>(accepted_trials)) * sum_w;
        samples.sigma_mt[oi * n_bins + ki] = sigma_mt;

        const double inv_w = (sum_w > 0.0) ? 1.0 / sum_w : 0.0;
        const double mean_par = sum_par * inv_w;
        const double mean_perp = sum_perp * inv_w;
        const double var_par = std::max(0.0, sum_par2 * inv_w - mean_par * mean_par);
        const double var_perp = std::max(0.0, sum_perp2 * inv_w - mean_perp * mean_perp);
        const size_t stats_idx = (oi * n_bins + ki) * 4;
        samples.dp_stats[stats_idx + 0] = mean_par;
        samples.dp_stats[stats_idx + 1] = mean_perp;
        samples.dp_stats[stats_idx + 2] = var_par;
        samples.dp_stats[stats_idx + 3] = var_perp;

        if (cdf_out && !local_weights.empty()) {
            // Optional high-fidelity runtime path: persist an inverse-CDF table
            // over momentum kicks so the handler can sample full dp vectors,
            // not just the mean/variance fallback in dp_stats.
            double total_w = 0.0;
            for (double w : local_weights) {
                total_w += w;
            }
            double acc = 0.0;
            const size_t n_local = local_weights.size();
            cdf_out->cdf.resize(n_local);
            for (size_t idx = 0; idx < n_local; ++idx) {
                acc += local_weights[idx];
                cdf_out->cdf[idx] = acc / total_w;
            }
            cdf_out->dp = std::move(local_dp);
        }
    };

    if (opt.store_full_cdf) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (size_t oi = 0; oi < n_orient; ++oi) {
            for (size_t ki = 0; ki < n_bins; ++ki) {
                const size_t cell_idx = oi * n_bins + ki;
                compute_cell(oi, ki, &local_cdfs[cell_idx]);
            }
        }
    } else {
        const size_t total_cells = n_orient * n_bins;
        std::vector<unsigned char> cell_done(total_cells, 0);
        size_t completed_cells = 0;
        int next_progress_pct = 10;

        if (opt.resume) {
            std::string resume_error;
            // Resume is intentionally checkpoint-level only: completed cells are
            // loaded from the same HDF5 format that final output uses.
            const bool resume_loaded = load_hdf5_resume_checkpoint(
                opt.output,
                logv_bins,
                samples,
                cell_done,
                completed_cells,
                resume_error
            );
            if (resume_loaded) {
                std::cout << "  resume loaded: " << completed_cells << "/" << total_cells << " cells" << std::endl;
                const int completed_pct = (total_cells > 0)
                    ? static_cast<int>((100.0 * static_cast<double>(completed_cells)) / static_cast<double>(total_cells))
                    : 100;
                while (next_progress_pct <= completed_pct && next_progress_pct <= 100) {
                    next_progress_pct += 10;
                }
            } else if (std::filesystem::exists(opt.output)) {
                std::cerr << resume_error << "\n";
                return EXIT_FAILURE;
            } else {
                std::cout << "  resume requested but no existing file found; starting fresh.\n";
            }
        }

        std::vector<size_t> todo_indices;
        todo_indices.reserve(total_cells - completed_cells);
        for (size_t idx = 0; idx < total_cells; ++idx) {
            if (cell_done[idx] == 0) {
                todo_indices.push_back(idx);
            }
        }

        const size_t batch_cells = (opt.checkpoint_cells > 0)
            ? static_cast<size_t>(opt.checkpoint_cells)
            : todo_indices.size();

        auto flush_checkpoint = [&](bool final_flush) {
            if (!(opt.format == "hdf5" || opt.format == "h5")) {
                return;
            }
            if (!final_flush && opt.checkpoint_cells <= 0) {
                return;
            }
            write_hdf5(opt.output, opt.species_id, opt.gas, opt, samples, &cell_done,
                       static_cast<long long>(completed_cells));
            std::cout << "  checkpoint: " << completed_cells << "/" << total_cells << " cells" << std::endl;
        };

        auto log_progress_deciles = [&]() {
            if (total_cells == 0) {
                return;
            }
            const int completed_pct = static_cast<int>(
                (100.0 * static_cast<double>(completed_cells)) / static_cast<double>(total_cells)
            );
            while (next_progress_pct <= completed_pct && next_progress_pct <= 100) {
                std::cout << "  progress: " << next_progress_pct << "% ("
                          << completed_cells << "/" << total_cells << " cells)" << std::endl;
                next_progress_pct += 10;
            }
        };

        for (size_t offset = 0; offset < todo_indices.size(); offset += std::max<size_t>(1, batch_cells)) {
            const size_t end = std::min(todo_indices.size(), offset + std::max<size_t>(1, batch_cells));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for (size_t pos = offset; pos < end; ++pos) {
                const size_t idx = todo_indices[pos];
                const size_t oi = idx / n_bins;
                const size_t ki = idx % n_bins;
                compute_cell(oi, ki, nullptr);
                cell_done[idx] = 1;
            }

            completed_cells += (end - offset);
            log_progress_deciles();
            flush_checkpoint(false);
        }

        if (completed_cells == total_cells && next_progress_pct <= 100) {
            std::cout << "  progress: 100% (" << completed_cells << "/" << total_cells << " cells)" << std::endl;
        }
        flush_checkpoint(true);
    }

    if (opt.store_full_cdf) {
        size_t running_offset = 0;
        for (size_t oi = 0; oi < n_orient; ++oi) {
            for (size_t ki = 0; ki < n_bins; ++ki) {
                const size_t idx = oi * n_bins + ki;
                const auto& local = local_cdfs[idx];
                samples.cdf_counts[idx] = static_cast<long long>(local.cdf.size());
                if (!local.cdf.empty()) {
                    // Flatten per-cell CDF/momentum samples into one contiguous
                    // payload; cdf_offsets/cdf_counts keep the cell boundaries.
                    samples.cdf_offsets[idx] = static_cast<long long>(running_offset);
                    samples.cdf_values.insert(samples.cdf_values.end(),
                                              local.cdf.begin(),
                                              local.cdf.end());
                    samples.dp_samples.insert(samples.dp_samples.end(),
                                              local.dp.begin(),
                                              local.dp.end());
                    running_offset += local.cdf.size();
                }
            }
        }
    }

    if (opt.format == "hdf5" || opt.format == "h5") {
        if (opt.store_full_cdf) {
            write_hdf5(opt.output, opt.species_id, opt.gas, opt, samples);
        }
    } else {
        std::cerr << "JSON output not implemented for interaction_potential_precompute\\n";
        return EXIT_FAILURE;
    }

    std::cout << "[interaction_potential_precompute] Done.\\n";
    return EXIT_SUCCESS;
}
