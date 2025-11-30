// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "SolverConfig.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

using json = nlohmann::json;
/**
 * @brief Load solver configuration from a JSON file.
 * @param filename Path to the JSON configuration file.
 * @return SolverConfig object populated with the configuration data.
 */
SolverConfig load_config_json(const std::string& filename)
{
    SolverConfig cfg;
    std::ifstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("Cannot open config: " + filename);

    json j;
    f >> j;

    // ---------------------------
    // Solver section
    // ---------------------------
    if (j.contains("solver")) {
        auto s = j["solver"];
        if (s.contains("tol"))        cfg.tol       = s["tol"].get<double>();
        if (s.contains("max_iters"))  cfg.max_iters = s["max_iters"].get<int>();
        if (s.contains("restart"))    cfg.restart   = s["restart"].get<int>();
        if (s.contains("fmm")) {
            auto fmmp = s["fmm"];
            if (fmmp.contains("use_fmm")) cfg.fmm.use_fmm = fmmp["use_fmm"].get<bool>();
            if (fmmp.contains("fmm_order")) cfg.fmm.fmm_order = fmmp["fmm_order"].get<int>();
            if (fmmp.contains("theta")) cfg.fmm.theta = fmmp["theta"].get<double>();
            if (fmmp.contains("max_depth")) cfg.fmm.max_depth = fmmp["max_depth"].get<int>();
            if (fmmp.contains("leaf_max")) cfg.fmm.leaf_max = fmmp["leaf_max"].get<int>();
        }
    }

    // ---------------------------
    // Grid section
    // ---------------------------
    if (j.contains("grid")) {
        auto g = j["grid"];
        if (g.contains("origin"))     cfg.origin = { g["origin"][0], g["origin"][1], g["origin"][2] };
        if (g.contains("size"))       cfg.size   = { g["size"][0],   g["size"][1],   g["size"][2]   };
        if (g.contains("resolution")) {
            cfg.Nx = g["resolution"][0];
            cfg.Ny = g["resolution"][1];
            cfg.Nz = g["resolution"][2];
        }
        if (g.contains("output"))     cfg.output  = g["output"].get<std::string>();
    }

    // ---------------------------
    // Geometry section
    // ---------------------------
    if (j.contains("geometry")) {
        for (auto& g : j["geometry"]) {
            GeometryEntry entry;

            if (g.contains("type")) entry.type = g["type"].get<std::string>();
            if (g.contains("file")) entry.file = g["file"].get<std::string>();
            if (g.contains("V"))    entry.V    = g["V"].get<double>();
            if (g.contains("id"))   entry.id   = g["id"].get<int>();
            if (g.contains("size")) entry.size = { g["size"][0], g["size"][1], g["size"][2] };
            if (g.contains("center")) entry.center = { g["center"][0], g["center"][1], g["center"][2] };
            if (g.contains("z"))    entry.z    = g["z"].get<double>();
            if (g.contains("unit")) entry.unit = g["unit"].get<std::string>();

            // --- extra parameters ---
            if (g.contains("r_inner")) entry.r_inner = g["r_inner"].get<double>();
            if (g.contains("r_outer")) entry.r_outer = g["r_outer"].get<double>();
            if (g.contains("radius"))  entry.radius  = g["radius"].get<double>();
            if (g.contains("length"))  entry.length  = g["length"].get<double>();
            if (g.contains("width"))   entry.width   = g["width"].get<double>();
            if (g.contains("thickness")) entry.thickness = g["thickness"].get<double>();
            if (g.contains("count"))   entry.count   = g["count"].get<int>();
            if (g.contains("spacing"))     entry.spacing = { g["spacing"][0], g["spacing"][1], g["spacing"][2] };
            if (g.contains("start_z")) entry.start_z = g["start_z"].get<double>();
            if (g.contains("id_start")) entry.id_start = g["id_start"].get<int>();
            if (g.contains("n_theta")) entry.n_theta = g["n_theta"].get<int>();
            if (g.contains("n_phi"))   entry.n_phi   = g["n_phi"].get<int>();

            if (g.contains("translate"))   entry.translate = { g["translate"][0], g["translate"][1], g["translate"][2] };
            if (g.contains("start_V"))     entry.start_V = g["start_V"].get<double>();
            if (g.contains("delta_V"))     entry.delta_V = g["delta_V"].get<double>();

            // --- Boundary condition ---
            if (g.contains("boundary")) {
                auto b = g["boundary"];
                std::string type = b.value("type", "Dirichlet");

                if (type == "dirichlet") {
                    entry.boundary.type = BoundaryType::Dirichlet;
                    entry.boundary.potential = b.value("potential", 0.0);
                } else if (type == "neumann") {
                    entry.boundary.type = BoundaryType::Neumann;
                    entry.insulator = true;
                } else {
                    std::cerr << "Warning: Unknown boundary type '" << type << "'. Defaulting to Dirichlet." << std::endl;
                    entry.boundary.type = BoundaryType::Dirichlet;
                }
            } else if (g.contains("V")) {
                // Legacy support: if only "V" is given, assume Dirichlet
                entry.boundary.type = BoundaryType::Dirichlet;
                entry.boundary.potential = g["V"].get<double>();
            } else if (g.value("insulator", false)) {
                // If marked as insulator, set Neumann BC
                entry.boundary.type = BoundaryType::Neumann;
                entry.insulator = true;
            }
            // --- Store electrode-level BC in map ---
            BoundaryCondition bc_entry = entry.boundary;
            if (entry.boundary.type == BoundaryType::Neumann && g.contains("boundary")) {
                auto b = g["boundary"];
                bc_entry.flux = b.value("flux", 0.0);
            }
            
            // --- Handle ring_stack voltage ramps ---
            if (entry.type == "ring_stack" && entry.count > 0) {
                // Create a voltage ramp for each ring in the stack
                for (int i = 0; i < entry.count; ++i) {
                    BoundaryCondition bc_ring;
                    bc_ring.type = BoundaryType::Dirichlet;
                    bc_ring.potential = entry.start_V + i * entry.delta_V;
                    cfg.bc_by_electrode[entry.id_start + i] = bc_ring;
                }
            } else {
                // Single electrode
                cfg.bc_by_electrode[entry.id] = bc_entry;
            }

            std::string file = g.value("file", "");
            if (entry.boundary.type == BoundaryType::Dirichlet && (file.find("spacer") != std::string::npos || file.find("insulator") != std::string::npos)) {
                entry.boundary.type = BoundaryType::Neumann;
                entry.insulator = true;
            }

            cfg.geometry.push_back(entry);
        }
    }
    
    // ---------------------------
    // Adaptive integration parameters
    // ---------------------------
    if (j.contains("adaptive")) {
        auto g = j["adaptive"];
        if (g.contains("near_factor")) cfg.adaptive.near_factor = g["near_factor"].get<double>();
        if (g.contains("gauss_u"))     cfg.adaptive.gauss_u     = g["gauss_u"].get<int>();
        if (g.contains("gauss_v"))     cfg.adaptive.gauss_v     = g["gauss_v"].get<int>();
    }

    return cfg;
}
