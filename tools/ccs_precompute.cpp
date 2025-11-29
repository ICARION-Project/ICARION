// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file ccs_precompute.cpp
 * @brief CLI to derive gas-specific CCS maps from reference CCS (HSS/EHSS)
 *
 * Usage:
 *   ccs_precompute --input species.json --output out.json --ref-gas He --ref-ccs-A2 110.0 --model HSS --override
 *
 * Formula:
 *   CCS_gas = π (r_ion + r_gas)^2
 *   r_ion   = max(0, sqrt(CCS_ref/π) - r_ref)
 *
 * Supported gases: He, Ar, CO2, Ne, N2, O2, H2O (radii from utils/constants.h)
 */

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <filesystem>

#include "utils/constants.h"
#include "core/types/CollisionTypes.h"
#include "core/physics/collisions/core/CollisionGeometry.h"
#include "core/io/moleculeLoader.h"
#include "core/physics/collisions/geometryUtils.h"
#include <json/json.h>

namespace {

// Gas radii (matching utils/constants.h)
std::unordered_map<std::string, double> gas_radius_m = {
    {"He", RADIUS_HE_M},
    {"Ar", RADIUS_AR_M},
    {"CO2", RADIUS_CO2_M},
    {"Ne", RADIUS_NE_M},
    {"N2", RADIUS_N2_M},
    {"O2", RADIUS_O2_M},
    {"H2O", RADIUS_H2O_M},
};

double derive_sigma_for_gas(double sigma_ref_m2, const std::string& gas_ref, const std::string& gas_target) {
    auto it_ref = gas_radius_m.find(gas_ref);
    auto it_tgt = gas_radius_m.find(gas_target);
    if (it_ref == gas_radius_m.end() || it_tgt == gas_radius_m.end()) {
        return 0.0;
    }
    double r_ref = it_ref->second;
    double r_ion = std::max(0.0, std::sqrt(std::max(sigma_ref_m2, 0.0) / M_PI) - r_ref);
    double r_tgt = it_tgt->second;
    return M_PI * (r_ion + r_tgt) * (r_ion + r_tgt);
}

bool parse_args(int argc, char** argv,
                std::string& input,
                std::string& output,
                std::string& ref_gas,
                double& ref_ccs_a2,
                std::string& model,
                bool& override,
                int& n_orientations,
                std::unordered_set<std::string>& target_species) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--input" && i + 1 < argc) {
            input = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            output = argv[++i];
        } else if (arg == "--ref-gas" && i + 1 < argc) {
            ref_gas = argv[++i];
        } else if (arg == "--ref-ccs-A2" && i + 1 < argc) {
            ref_ccs_a2 = std::stod(argv[++i]);
        } else if (arg == "--model" && i + 1 < argc) {
            model = argv[++i];
        } else if (arg == "--override") {
            override = true;
        } else if (arg == "--n-orientations" && i + 1 < argc) {
            n_orientations = std::stoi(argv[++i]);
        } else if (arg == "--species" && i + 1 < argc) {
            std::string list = argv[++i];
            size_t start = 0;
            while (start < list.size()) {
                size_t pos = list.find(',', start);
                std::string token = (pos == std::string::npos) ? list.substr(start) : list.substr(start, pos - start);
                if (!token.empty()) {
                    target_species.insert(token);
                }
                if (pos == std::string::npos) break;
                start = pos + 1;
            }
        } else if (arg == "--help" || arg == "-h") {
            return false;
        }
    }
    if (input.empty() || output.empty() || ref_gas.empty() || ref_ccs_a2 <= 0.0 || target_species.empty()) {
        return false;
    }
    if (model != "HSS" && model != "EHSS") {
        model = "HSS";
    }
    return true;
}

double compute_oapa_ccs(const std::vector<Vec3>& centers,
                        const std::vector<double>& radii,
                        double neutral_radius,
                        int n_orientations) {
    ICARION::physics::EhssRng rng(12345);
    double A_sum = 0.0;
    double R[3][3];

    for (int k = 0; k < n_orientations; ++k) {
        ICARION::physics::collision_core::CollisionGeometry::generate_random_rotation(rng, R);

        double Rmax = 0.0;
        for (size_t i = 0; i < centers.size(); ++i) {
            Vec3 c{
                R[0][0] * centers[i].x + R[0][1] * centers[i].y + R[0][2] * centers[i].z,
                R[1][0] * centers[i].x + R[1][1] * centers[i].y + R[1][2] * centers[i].z,
                R[2][0] * centers[i].x + R[2][1] * centers[i].y + R[2][2] * centers[i].z
            };
            double r_eff = std::sqrt(c.x * c.x + c.y * c.y) + radii[i] + neutral_radius;
            Rmax = std::max(Rmax, r_eff);
        }

        double A_k = M_PI * Rmax * Rmax;
        A_sum += A_k;
    }

    return A_sum / static_cast<double>(n_orientations);
}

}  // namespace

int main(int argc, char** argv) {
    std::string input;
    std::string output;
    std::string ref_gas;
    double ref_ccs_a2 = 0.0;
    std::string model = "HSS";
    bool override = false;
    int n_orientations = 300;

    std::unordered_set<std::string> target_species;

    if (!parse_args(argc, argv, input, output, ref_gas, ref_ccs_a2, model, override, n_orientations, target_species)) {
        std::cerr << "Usage: ccs_precompute --input species.json --output out.json "
                  << "--species H3O+,NH4+ --ref-gas He --ref-ccs-A2 110.0 "
                  << "[--model HSS|EHSS] [--override] [--n-orientations 300]\n";
        return 1;
    }

    Json::Value root;
    {
        std::ifstream ifs(input);
        if (!ifs) {
            std::cerr << "Failed to open input JSON: " << input << "\n";
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

    constexpr double A2_TO_M2 = 1e-20;
    constexpr double M2_TO_A2 = 1e20;

    double sigma_ref_m2 = ref_ccs_a2 * A2_TO_M2;
    std::string key_model = "CCS_" + model;

    auto species_names = root["species"].getMemberNames();
    for (const auto& sid : species_names) {
        if (!target_species.count(sid)) {
            continue;  // do not touch other species
        }
        auto& props = root["species"][sid];
        // Skip neutrals
        if (props.get("charge", 0).asInt() == 0) continue;
        if (props.isMember(key_model) && !override) continue;

        Json::Value ccs_map(Json::objectValue);
        Json::Value ccs_ehss_map(Json::objectValue);

        // HSS approximation: kinetic diameter mapping
        for (const auto& g : gas_radius_m) {
            double sigma = derive_sigma_for_gas(sigma_ref_m2, ref_gas, g.first);
            if (sigma > 0.0) {
                ccs_map[g.first] = sigma * M2_TO_A2;
            }
        }

        // EHSS approximation: OAPA projection using geometry_file if present
        if (model == "EHSS") {
            if (!props.isMember("geometry_file") || !props["geometry_file"].isString()) {
                std::cerr << "[ccs_precompute] Warning: species '" << sid
                          << "' has no geometry_file; skipping EHSS map\n";
            } else {
                std::string geom_path = props["geometry_file"].asString();
                // resolve relative to input file directory
                std::filesystem::path base = std::filesystem::path(input).parent_path();
                std::filesystem::path abs_geom = base / geom_path;
                if (std::filesystem::exists(abs_geom)) {
                try {
                    auto molecule = ICARION::io::load_molecule(abs_geom.string());
                    auto geom = ICARION::physics::convert_molecule_to_geometry(molecule);
                    for (const auto& g : gas_radius_m) {
                        double area = compute_oapa_ccs(geom.first, geom.second, g.second, n_orientations);
                        ccs_ehss_map[g.first] = area * M2_TO_A2;
                    }
                } catch (const std::exception& e) {
                    std::cerr << "Warning: failed to load geometry for " << sid << ": " << e.what() << "\n";
                }
            }
        }

        props["CCS_reference_gas"] = ref_gas;
        props["CCS_model"] = model;
        props[key_model] = ccs_map;
        if (ccs_ehss_map.size() > 0) {
            props["CCS_EHSS"] = ccs_ehss_map;
        }
    }

    std::ofstream ofs(output);
    if (!ofs) {
        std::cerr << "Failed to open output file: " << output << "\n";
        return 1;
    }
    ofs << root.toStyledString();
    std::cout << "Wrote updated species database to " << output << "\n";
    return 0;
}
}
 
