// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

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

#include "utils/constants.h"
#include "core/io/moleculeLoader.h"
#include "core/physics/collisions/geometryUtils.h"
#include "core/physics/collisions/core/CollisionGeometry.h"
#include <json/json.h>

namespace {

struct Options {
    std::string input;
    std::string output;
    std::string species_id;
    int n_orientations = 300;
    int n_samples = 8000;
    uint64_t seed = 12345;
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
        } else if (arg == "--n-orientations" && i + 1 < argc) {
            opt.n_orientations = std::stoi(argv[++i]);
        } else if (arg == "--n-samples" && i + 1 < argc) {
            opt.n_samples = std::stoi(argv[++i]);
        } else if (arg == "--seed" && i + 1 < argc) {
            opt.seed = static_cast<uint64_t>(std::stoull(argv[++i]));
        } else if (arg == "--help" || arg == "-h") {
            return false;
        }
    }

    if (opt.input.empty() || opt.output.empty() || opt.species_id.empty()) {
        return false;
    }
    if (opt.n_orientations <= 0 || opt.n_samples <= 0) {
        return false;
    }

    return true;
}

std::vector<std::pair<std::string, double>> gas_radii() {
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

}  // namespace

int main(int argc, char** argv) {
    Options opt;
    if (!parse_args(argc, argv, opt)) {
        std::cerr << "Usage: ehss_samples_precompute --input species.json --output out.json --species H3O+ "
                  << "[--n-orientations 300] [--n-samples 8000] [--seed 12345]\n";
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

    auto geometry = ICARION::physics::convert_molecule_to_geometry(molecule, true);
    const auto& centers = geometry.first;
    const auto& radii = geometry.second;
    if (centers.empty()) {
        std::cerr << "Geometry file contains no atoms\n";
        return 1;
    }

    const auto gases = gas_radii();
    double max_gas_radius = 0.0;
    for (const auto& g : gases) {
        max_gas_radius = std::max(max_gas_radius, g.second);
    }

    std::mt19937 rng(static_cast<uint32_t>(opt.seed));
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    std::vector<std::array<double, 4>> orientations_quat;
    orientations_quat.reserve(opt.n_orientations);

    std::unordered_map<std::string, std::vector<double>> areas_by_gas_m2;
    for (const auto& g : gases) {
        areas_by_gas_m2[g.first].reserve(opt.n_orientations);
    }

    const auto t_start = std::chrono::steady_clock::now();
    const int progress_step = std::max(1, opt.n_orientations / 10);
    for (int k = 0; k < opt.n_orientations; ++k) {
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

        const double b_max = b_max_base + max_gas_radius;
        const double area_sample = M_PI * b_max * b_max;

        std::vector<double> d_min;
        d_min.reserve(opt.n_samples);
        for (int s = 0; s < opt.n_samples; ++s) {
            double u1 = uniform(rng);
            double u2 = uniform(rng);
            double r = b_max * std::sqrt(u1);
            double phi = 2.0 * M_PI * u2;
            double px = r * std::cos(phi);
            double py = r * std::sin(phi);

            double min_val = std::numeric_limits<double>::max();
            for (size_t i = 0; i < rotated_centers.size(); ++i) {
                double dx = px - rotated_centers[i].x;
                double dy = py - rotated_centers[i].y;
                double dist = std::sqrt(dx * dx + dy * dy);
                double val = dist - radii[i];
                if (val < min_val) {
                    min_val = val;
                }
            }
            d_min.push_back(min_val);
        }

        for (const auto& g : gases) {
            const double Rn = g.second;
            int count = 0;
            for (double val : d_min) {
                if (val <= Rn) {
                    ++count;
                }
            }
            double area = (static_cast<double>(count) / static_cast<double>(opt.n_samples)) * area_sample;
            areas_by_gas_m2[g.first].push_back(area);
        }

        if ((k + 1) % progress_step == 0 || (k + 1) == opt.n_orientations) {
            const auto t_now = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = t_now - t_start;
            double pct = 100.0 * (static_cast<double>(k + 1) / static_cast<double>(opt.n_orientations));
            std::cout << "[ehss_samples_precompute] orientation "
                      << (k + 1) << "/" << opt.n_orientations
                      << " (" << std::fixed << std::setprecision(1) << pct << "%)"
                      << " elapsed " << std::setprecision(2) << elapsed.count() << "s\n";
        }
    }

    Json::Value out;
    out["version"] = 1;
    out["species_id"] = opt.species_id;
    out["units"] = "m,m2";
    out["n_orientations"] = opt.n_orientations;
    out["n_samples_per_orientation"] = opt.n_samples;
    out["method"] = "projection_union_monte_carlo";

    Json::Value orient_json(Json::arrayValue);
    for (const auto& q : orientations_quat) {
        Json::Value entry(Json::arrayValue);
        entry.append(q[0]);
        entry.append(q[1]);
        entry.append(q[2]);
        entry.append(q[3]);
        orient_json.append(entry);
    }
    out["orientations_quat"] = orient_json;

    Json::Value areas_json(Json::objectValue);
    for (const auto& g : gases) {
        const auto& vec = areas_by_gas_m2[g.first];
        Json::Value arr(Json::arrayValue);
        for (double v : vec) {
            arr.append(v);
        }
        areas_json[g.first] = arr;
    }
    out["areas_by_gas_m2"] = areas_json;

    std::ofstream ofs(opt.output);
    if (!ofs) {
        std::cerr << "Failed to open output file: " << opt.output << "\n";
        return 1;
    }
    ofs << out.toStyledString();
    std::cout << "Wrote EHSS samples to " << opt.output << "\n";
    return 0;
}
