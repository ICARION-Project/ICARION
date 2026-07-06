// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <random>
#include <limits>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <filesystem>

#include <H5Cpp.h>
#include <json/json.h>

#include "utils/constants.h"
#include "core/io/moleculeLoader.h"
#include "core/types/Vec3.h"
#include "core/physics/collisions/EHSSOfflineSampleSet.h"
#include "core/physics/collisions/geometryUtils.h"
#include "core/physics/collisions/core/CollisionGeometry.h"

namespace {

enum class OutputFormat { Json, Hdf5 };

struct Options {
    std::string input;
    std::string output;
    std::string species_id;
    std::string gas = "He";
    std::string format;
    int n_orientations = 1000;
    int n_sigma_samples = 10000;
    int n_mu_samples = 2000;
    int max_tries = 0;
    uint64_t seed = 12345;
    bool center_on_com = true;
};

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
        } else if (arg == "--format" && i + 1 < argc) {
            opt.format = argv[++i];
        } else if (arg == "--n-orientations" && i + 1 < argc) {
            opt.n_orientations = std::stoi(argv[++i]);
        } else if (arg == "--n-sigma-samples" && i + 1 < argc) {
            opt.n_sigma_samples = std::stoi(argv[++i]);
        } else if (arg == "--n-mu-samples" && i + 1 < argc) {
            opt.n_mu_samples = std::stoi(argv[++i]);
        } else if (arg == "--max-tries" && i + 1 < argc) {
            opt.max_tries = std::stoi(argv[++i]);
        } else if (arg == "--seed" && i + 1 < argc) {
            opt.seed = static_cast<uint64_t>(std::stoull(argv[++i]));
        } else if (arg == "--no-center-on-com") {
            opt.center_on_com = false;
        } else if (arg == "--help" || arg == "-h") {
            return false;
        }
    }

    if (opt.input.empty() || opt.output.empty() || opt.species_id.empty()) {
        return false;
    }
    if (opt.n_orientations <= 0 || opt.n_sigma_samples <= 0 || opt.n_mu_samples <= 0) {
        return false;
    }
    return true;
}

OutputFormat resolve_format(const Options& opt) {
    if (!opt.format.empty()) {
        if (opt.format == "json") {
            return OutputFormat::Json;
        }
        if (opt.format == "hdf5" || opt.format == "h5") {
            return OutputFormat::Hdf5;
        }
        throw std::runtime_error("Unknown --format (use json or hdf5)");
    }

    const std::filesystem::path out_path(opt.output);
    const auto ext = out_path.extension().string();
    if (ext == ".h5" || ext == ".hdf5") {
        return OutputFormat::Hdf5;
    }
    return OutputFormat::Json;
}

std::unordered_map<std::string, double> gas_radii() {
    return {
        {"He", RADIUS_HE_M},
        {"Ar", RADIUS_AR_M},
        {"CO2", RADIUS_CO2_M},
        {"Ne", RADIUS_NE_M},
        {"N2", RADIUS_N2_M},
        {"O2", RADIUS_O2_M},
        {"H2O", RADIUS_H2O_M}
    };
}

void quaternion_to_rotation(double w, double x, double y, double z, double R[3][3]) {
    const double norm = std::sqrt(w * w + x * x + y * y + z * z);
    if (norm <= 0.0) {
        R[0][0] = 1.0; R[0][1] = 0.0; R[0][2] = 0.0;
        R[1][0] = 0.0; R[1][1] = 1.0; R[1][2] = 0.0;
        R[2][0] = 0.0; R[2][1] = 0.0; R[2][2] = 1.0;
        return;
    }

    w /= norm;
    x /= norm;
    y /= norm;
    z /= norm;

    R[0][0] = 1.0 - 2.0 * (y * y + z * z);
    R[0][1] = 2.0 * (x * y - z * w);
    R[0][2] = 2.0 * (x * z + y * w);
    R[1][0] = 2.0 * (x * y + z * w);
    R[1][1] = 1.0 - 2.0 * (x * x + z * z);
    R[1][2] = 2.0 * (y * z - x * w);
    R[2][0] = 2.0 * (x * z - y * w);
    R[2][1] = 2.0 * (y * z + x * w);
    R[2][2] = 1.0 - 2.0 * (x * x + y * y);
}

void sample_random_quaternion(std::mt19937& rng, double& w, double& x, double& y, double& z) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double u1 = uniform(rng);
    double u2 = uniform(rng);
    double u3 = uniform(rng);

    double q1 = std::sqrt(1.0 - u1) * std::sin(2.0 * M_PI * u2);
    double q2 = std::sqrt(1.0 - u1) * std::cos(2.0 * M_PI * u2);
    double q3 = std::sqrt(u1) * std::sin(2.0 * M_PI * u3);
    double q4 = std::sqrt(u1) * std::cos(2.0 * M_PI * u3);

    x = q1;
    y = q2;
    z = q3;
    w = q4;
}

struct CollisionResult {
    bool hit = false;
    Vec3 contact_normal{0.0, 0.0, 0.0};
};

CollisionResult detect_atom_collision(
    const std::vector<Vec3>& rotated_atoms,
    const std::vector<double>& atom_radii,
    const Vec3& neutral_offset,
    double gas_radius
) {
    constexpr double min_contact_dist_sq = 1e-24;
    const Vec3 collision_axis{0.0, 0.0, 1.0};
    CollisionResult best{false, Vec3{0.0, 0.0, 0.0}};
    double best_s_hit = 0.0;

    // Runtime EHSS uses the first hard-sphere contact along the relative
    // trajectory. The offline table must preserve that rule, so for each
    // projected neutral ray we test every atom sphere and keep the earliest
    // intersection along +z.
    for (size_t j = 0; j < rotated_atoms.size(); ++j) {
        const Vec3& ra = rotated_atoms[j];
        double Rsum = atom_radii[j] + gas_radius;
        double Rsum2 = Rsum * Rsum;

        Vec3 rel = neutral_offset - ra;
        double sstar = rel.z;  // dot(collision_axis, rel)
        Vec3 dvec{rel.x, rel.y, 0.0};
        double dmin2 = dvec.x * dvec.x + dvec.y * dvec.y;

        if (dmin2 <= Rsum2) {
            double h = std::sqrt(std::max(0.0, Rsum2 - dmin2));
            double s_hit = sstar - h;
            Vec3 p_hit_rel{rel.x, rel.y, rel.z - s_hit};

            double phr2 = p_hit_rel.x * p_hit_rel.x + p_hit_rel.y * p_hit_rel.y + p_hit_rel.z * p_hit_rel.z;
            Vec3 n_contact;
            if (phr2 <= min_contact_dist_sq) {
                n_contact = Vec3{0.0, 0.0, -1.0};
            } else {
                double inv = 1.0 / std::sqrt(phr2);
                n_contact = Vec3{p_hit_rel.x * inv, p_hit_rel.y * inv, p_hit_rel.z * inv};
                if (n_contact.z > 0.0) {
                    n_contact = Vec3{-n_contact.x, -n_contact.y, -n_contact.z};
                }
            }

            if (!best.hit || s_hit < best_s_hit) {
                best = CollisionResult{true, n_contact};
                best_s_hit = s_hit;
            }
        }
    }

    return best;
}

void write_json(
    const Options& opt,
    const std::vector<std::array<double, 4>>& orientations,
    const std::vector<double>& sigma_eff,
    const std::vector<std::vector<double>>& mu_samples
) {
    Json::Value out;
    // Keep JSON and HDF5 semantically identical: the runtime loader validates
    // the same format/unit markers and reads the same logical arrays.
    out["version"] = 1;
    out["species_id"] = opt.species_id;
    out["gas"] = opt.gas;
    out["format"] = ICARION::physics::EHSS_OFFLINE_SAMPLE_SET_FORMAT;
    out["units"] = ICARION::physics::EHSS_OFFLINE_SAMPLE_SET_UNITS;
    out["n_orientations"] = static_cast<int>(orientations.size());
    out["n_sigma_samples"] = opt.n_sigma_samples;
    out["n_mu_samples"] = opt.n_mu_samples;
    out["seed"] = static_cast<Json::UInt64>(opt.seed);

    Json::Value orient_json(Json::arrayValue);
    for (const auto& q : orientations) {
        Json::Value entry(Json::arrayValue);
        entry.append(q[0]);
        entry.append(q[1]);
        entry.append(q[2]);
        entry.append(q[3]);
        orient_json.append(entry);
    }
    out["orientations_quat"] = orient_json;

    Json::Value sigma_json(Json::arrayValue);
    for (double v : sigma_eff) {
        sigma_json.append(v);
    }
    out["sigma_eff_m2"] = sigma_json;

    Json::Value mu_json(Json::arrayValue);
    for (const auto& row : mu_samples) {
        Json::Value row_json(Json::arrayValue);
        for (double v : row) {
            row_json.append(v);
        }
        mu_json.append(row_json);
    }
    out["mu_samples"] = mu_json;

    std::ofstream ofs(opt.output);
    if (!ofs) {
        throw std::runtime_error("Failed to open output file: " + opt.output);
    }
    ofs << out.toStyledString();
}

void write_hdf5(
    const Options& opt,
    const std::vector<std::array<double, 4>>& orientations,
    const std::vector<double>& sigma_eff,
    const std::vector<std::vector<double>>& mu_samples
) {
    const hsize_t n_orient = static_cast<hsize_t>(orientations.size());
    const hsize_t n_mu = static_cast<hsize_t>(opt.n_mu_samples);

    H5::H5File file(opt.output, H5F_ACC_TRUNC);

    auto write_attr_str = [&](const std::string& name, const std::string& value) {
        H5::StrType str_type(H5::PredType::C_S1, value.size());
        H5::DataSpace space(H5S_SCALAR);
        H5::Attribute attr = file.createAttribute(name, str_type, space);
        attr.write(str_type, value);
    };

    auto write_attr_int = [&](const std::string& name, long long value) {
        H5::DataSpace space(H5S_SCALAR);
        H5::Attribute attr = file.createAttribute(name, H5::PredType::NATIVE_LLONG, space);
        attr.write(H5::PredType::NATIVE_LLONG, &value);
    };

    write_attr_str("format", ICARION::physics::EHSS_OFFLINE_SAMPLE_SET_FORMAT);
    write_attr_str("species_id", opt.species_id);
    write_attr_str("gas", opt.gas);
    write_attr_str("units", ICARION::physics::EHSS_OFFLINE_SAMPLE_SET_UNITS);
    write_attr_int("version", 1);
    write_attr_int("n_orientations", static_cast<long long>(n_orient));
    write_attr_int("n_sigma_samples", opt.n_sigma_samples);
    write_attr_int("n_mu_samples", opt.n_mu_samples);
    write_attr_int("seed", static_cast<long long>(opt.seed));

    // orientations: (N,4), stored as quaternions so a generated file remains
    // auditable against the molecular orientations used during precompute.
    {
        hsize_t dims[2] = {n_orient, 4};
        H5::DataSpace space(2, dims);
        H5::DataSet dset = file.createDataSet("orientations_quat", H5::PredType::NATIVE_DOUBLE, space);
        std::vector<double> flat;
        flat.reserve(orientations.size() * 4);
        for (const auto& q : orientations) {
            flat.push_back(q[0]);
            flat.push_back(q[1]);
            flat.push_back(q[2]);
            flat.push_back(q[3]);
        }
        dset.write(flat.data(), H5::PredType::NATIVE_DOUBLE);
    }

    // sigma_eff: (N,), one projected collision cross section per orientation.
    {
        hsize_t dims[1] = {n_orient};
        H5::DataSpace space(1, dims);
        H5::DataSet dset = file.createDataSet("sigma_eff_m2", H5::PredType::NATIVE_DOUBLE, space);
        dset.write(sigma_eff.data(), H5::PredType::NATIVE_DOUBLE);
    }

    // mu_samples: (N, M), conditional first-contact scattering samples. Runtime
    // samples a row after selecting an orientation and reconstructs the elastic
    // post-collision relative velocity from mu plus a random azimuth.
    {
        hsize_t dims[2] = {n_orient, n_mu};
        H5::DataSpace space(2, dims);
        H5::DataSet dset = file.createDataSet("mu_samples", H5::PredType::NATIVE_DOUBLE, space);
        std::vector<double> flat;
        flat.reserve(static_cast<size_t>(n_orient * n_mu));
        for (const auto& row : mu_samples) {
            flat.insert(flat.end(), row.begin(), row.end());
        }
        dset.write(flat.data(), H5::PredType::NATIVE_DOUBLE);
    }
}

}  // namespace

int main(int argc, char** argv) {
    Options opt;
    if (!parse_args(argc, argv, opt)) {
        std::cerr << "Usage: ehss_offline_precompute --input species.json --output out.json --species H3O+ "
                  << "[--gas He] [--n-orientations 1000] [--n-sigma-samples 10000] [--n-mu-samples 2000] "
                  << "[--seed 12345] [--format json|hdf5]\n";
        return 1;
    }

    const auto gases = gas_radii();
    auto it_gas = gases.find(opt.gas);
    if (it_gas == gases.end()) {
        std::cerr << "Unsupported gas '" << opt.gas << "'\n";
        return 1;
    }
    const double gas_radius = it_gas->second;

    OutputFormat fmt;
    try {
        fmt = resolve_format(opt);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    Json::Value root;
    {
        std::ifstream ifs(opt.input);
        if (!ifs) {
            std::cerr << "Failed to open input JSON: " << opt.input << "\n";
            return 1;
        }
        Json::CharReaderBuilder builder;
        std::string errs;
        if (!Json::parseFromStream(builder, ifs, &root, &errs)) {
            std::cerr << "Failed to parse input JSON: " << errs << "\n";
            return 1;
        }
    }

    if (!root.isMember("species") || !root["species"].isObject()) {
        std::cerr << "Input JSON must contain a 'species' object\n";
        return 1;
    }

    const Json::Value& species = root["species"];
    if (!species.isMember(opt.species_id)) {
        std::cerr << "Species '" << opt.species_id << "' not found in input\n";
        return 1;
    }

    const Json::Value& props = species[opt.species_id];
    if (!props.isMember("geometry_file") || !props["geometry_file"].isString()) {
        std::cerr << "Species '" << opt.species_id << "' has no geometry_file\n";
        return 1;
    }

    std::filesystem::path base = std::filesystem::path(opt.input).parent_path();
    std::filesystem::path geom_path = base / props["geometry_file"].asString();
    if (!std::filesystem::exists(geom_path)) {
        std::cerr << "Geometry file not found: " << geom_path << "\n";
        return 1;
    }

    ICARION::io::Molecule molecule;
    try {
        molecule = ICARION::io::load_molecule(geom_path.string());
    } catch (const std::exception& e) {
        std::cerr << "Failed to load geometry: " << e.what() << "\n";
        return 1;
    }

    auto geometry = ICARION::physics::convert_molecule_to_geometry(molecule, opt.center_on_com);
    const auto& centers = geometry.first;
    const auto& radii = geometry.second;
    if (centers.empty()) {
        std::cerr << "Geometry file contains no atoms\n";
        return 1;
    }

    std::mt19937 rng(static_cast<uint32_t>(opt.seed));
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    std::vector<std::array<double, 4>> orientations_quat;
    orientations_quat.reserve(opt.n_orientations);
    std::vector<double> sigma_eff_m2;
    sigma_eff_m2.reserve(opt.n_orientations);
    std::vector<std::vector<double>> mu_samples;
    mu_samples.reserve(opt.n_orientations);

    const auto t_start = std::chrono::steady_clock::now();
    const int progress_step = std::max(1, opt.n_orientations / 10);

    for (int k = 0; k < opt.n_orientations; ++k) {
        // Each row in the output represents one random molecular orientation.
        // All atoms are rotated once, then reused for the sigma and mu Monte
        // Carlo loops below.
        double w = 1.0, x = 0.0, y = 0.0, z = 0.0;
        sample_random_quaternion(rng, w, x, y, z);
        orientations_quat.push_back({w, x, y, z});

        double Rm[3][3];
        quaternion_to_rotation(w, x, y, z, Rm);

        std::vector<Vec3> rotated_centers;
        rotated_centers.reserve(centers.size());
        double b_max_base = 0.0;
        for (size_t i = 0; i < centers.size(); ++i) {
            Vec3 ra = ICARION::physics::collision_core::CollisionGeometry::rotate_vector(centers[i], Rm);
            rotated_centers.push_back(ra);
            double r_xy = std::sqrt(ra.x * ra.x + ra.y * ra.y) + radii[i];
            b_max_base = std::max(b_max_base, r_xy);
        }

        const double b_max = b_max_base + gas_radius;
        const double area_sample = M_PI * b_max * b_max;

        // Estimate the orientation-specific collision cross section by sampling
        // impact points in the projected disk and counting rays that hit any
        // atom-centered hard sphere.
        int hit_count = 0;
        for (int s = 0; s < opt.n_sigma_samples; ++s) {
            double u1 = uniform(rng);
            double u2 = uniform(rng);
            double r = b_max * std::sqrt(u1);
            double phi = 2.0 * M_PI * u2;
            Vec3 neutral_offset{r * std::cos(phi), r * std::sin(phi), 0.0};

            auto collision = detect_atom_collision(rotated_centers, radii, neutral_offset, gas_radius);
            if (collision.hit) {
                ++hit_count;
            }
        }

        double sigma_eff = (static_cast<double>(hit_count) / static_cast<double>(opt.n_sigma_samples)) * area_sample;
        sigma_eff_m2.push_back(sigma_eff);

        const double p_hit = (area_sample > 0.0) ? (sigma_eff / area_sample) : 0.0;
        int max_tries = opt.max_tries;
        if (max_tries <= 0) {
            const double p_safe = std::max(p_hit, 1e-6);
            max_tries = static_cast<int>(std::ceil(opt.n_mu_samples / p_safe * 4.0));
            max_tries = std::max(max_tries, opt.n_mu_samples * 20);
        }

        // Store the conditional scattering distribution seen by the runtime:
        // resample impact points until a hit occurs, take the first atom contact,
        // and keep mu = cos(theta) from that contact normal.
        std::vector<double> mu_row;
        mu_row.reserve(opt.n_mu_samples);
        int tries = 0;
        while (static_cast<int>(mu_row.size()) < opt.n_mu_samples) {
            if (tries++ > max_tries) {
                std::cerr << "Failed to collect enough mu samples at orientation " << k
                          << " (hit probability ~" << p_hit << "). "
                          << "Increase --max-tries or reduce --n-mu-samples.\n";
                return 1;
            }
            double u1 = uniform(rng);
            double u2 = uniform(rng);
            double r = b_max * std::sqrt(u1);
            double phi = 2.0 * M_PI * u2;
            Vec3 neutral_offset{r * std::cos(phi), r * std::sin(phi), 0.0};

            auto collision = detect_atom_collision(rotated_centers, radii, neutral_offset, gas_radius);
            if (!collision.hit) {
                continue;
            }
            double mu = -collision.contact_normal.z;
            mu = std::min(1.0, std::max(-1.0, mu));
            mu_row.push_back(mu);
        }
        mu_samples.push_back(std::move(mu_row));

        if ((k + 1) % progress_step == 0 || (k + 1) == opt.n_orientations) {
            const auto t_now = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = t_now - t_start;
            double pct = 100.0 * (static_cast<double>(k + 1) / static_cast<double>(opt.n_orientations));
            std::cout << "[ehss_offline_precompute] orientation "
                      << (k + 1) << "/" << opt.n_orientations
                      << " (" << std::fixed << std::setprecision(1) << pct << "%)"
                      << " elapsed " << std::setprecision(2) << elapsed.count() << "s\n";
        }
    }

    try {
        if (fmt == OutputFormat::Json) {
            write_json(opt, orientations_quat, sigma_eff_m2, mu_samples);
        } else {
            write_hdf5(opt, orientations_quat, sigma_eff_m2, mu_samples);
        }
    } catch (const std::exception& e) {
        std::cerr << "Failed to write output: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Wrote EHSS offline samples to " << opt.output << "\n";
    return 0;
}
