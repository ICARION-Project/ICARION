// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#pragma once
#include "core/utils/mathUtils.h"
#include <string>
#include <unordered_map>
#include <vector>

enum class BoundaryType {
    Dirichlet,
    Neumann
};

struct BoundaryCondition {
    BoundaryType type = BoundaryType::Dirichlet;
    double potential = 0.0; // for Dirichlet [V]
    double flux = 0.0;      // for Neumann [V/m]
};

using ElectrodeBCMap = std::unordered_map<int, BoundaryCondition>; // electrode_id -> BC

// Struct for geometry entries
struct GeometryEntry {
    std::string type;       // "obj", "stl", "plate", "ring", etc.
    std::string file;       // for obj/stl
    double V = 0.0;         // potential in volts
    int id = 0;             // electrode ID

    // Common geometric parameters
    Vec3 size{0,0,0};       // plate/box size
    Vec3 center{0,0,0};     // optional center position
    double z = 0.0;         // for planar primitives
    std::string unit = "auto"; // "m", "mm", "um", or "auto"

    // Ring and cylinder parameters
    double r_inner = 0.0;
    double r_outer = 0.0;
    // Per-ring taper: change in radius per ring (m). Positive = increasing radius along +z
    double r_taper = 0.0;
    double radius  = 0.0;
    double length  = 0.0;
    double width   = 0.0;
    double thickness = 0.0;

    // Stack / multi-element
    int count = 1;
    Vec3 spacing = {0.0, 0.0, 0.0};
    double start_z = 0.0;
    int id_start = 0;
    Vec3 translate = {0.0,0.0,0.0};
    double start_V = 0.0;
    double delta_V = 0.0;

    // Sphere / disk resolution
    int n_theta = 64;
    int n_phi   = 32;

    // Boundary condition for this geometry
    BoundaryCondition boundary; // Boundary condition for this geometry
    bool insulator = false;    // Whether this geometry is an insulator (Neumann BC)
};

struct AdaptiveParams {
    double near_factor = 2.0;   // multiplier for near/far threshold
    int gauss_u = 4;            // Gauss-Legendre nodes in u direction
    int gauss_v = 4;            // Gauss-Legendre nodes in v direction
};

struct FMMParams {
    bool use_fmm = false;      // Whether to use FMM acceleration
    int fmm_order = 3;         // Multipole expansion order (p_max if adaptive)
    double theta = 0.7;        // MAC threshold: (rA + rB) / dist < theta
    int max_depth = 10;        // Maximum tree depth
    int leaf_max = 100;        // Maximum panels per leaf
    
    // Adaptive expansion order
    bool adaptive_p = false;        // Enable variable expansion order
    int adaptive_order_min = 3;     // Minimum order (p_min) for far interactions
    double adaptive_factor = 2.0;   // Separation threshold: if sep_ratio > factor, use p_min
};

struct SolverConfig {
    // solver
    double tol = 1e-6;
    int max_iters = 500;
    int restart = 50;

    // grid
    Vec3 origin{-5,-5,-1};
    Vec3 size{10,10,3};
    int Nx = 50, Ny = 50, Nz = 50;
    std::string output = "field_bem.h5";

    // geometry
    std::vector<GeometryEntry> geometry;

    // adaptive integration parameters
    AdaptiveParams adaptive;
    // boundary conditions
    ElectrodeBCMap bc_by_electrode; // Map of electrode ID to boundary condition

    // FMM parameters
    FMMParams fmm;

};

SolverConfig load_config_json(const std::string& filename);
