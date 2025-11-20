// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        field_update_api.h
 *   @brief       Thread-safe field update API for live field computation
 *
 *   @details
 *   Provides a thread-safe interface for managing live updates of electric
 *   field configurations in a simulation environment. Supports geometry
 *   updates, grid reconfiguration, and on-demand field recomputation.
 *
 *   @date        2025-11-09
 *   @version     1.0.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */
// --------------------------------------------------------------------
#pragma once

#include "core/physics/fields/fieldSampling.h"
#include <mutex>
#include <atomic>
#include <string>

// Re-export core types for backward compatibility
using ICARION::fields::Vec3d;
using ICARION::fields::Grid3DSnapshot;
using ICARION::fields::FieldSnapshot;
using ICARION::fields::sample_field;

struct GeometryConfig {
    // Example parameters (extend as needed)
    double radius_m{0.0};
    double length_m{0.0};
    double quad_dc_V{0.0};
    double rf_V{0.0};
    double rf_freq_Hz{0.0};
    
    // Origin offset (default: 0,0,0)
    double origin_x_m{0.0};
    double origin_y_m{0.0};
    double origin_z_m{0.0};
    // ... electrode offsets, tilts, etc.
};

// Thread-safe field server for live field computation
class FieldServer {
public:
    FieldServer();
    // Marks internal state as dirty so next "get" ensures recalc or fetch
    void set_geometry(const GeometryConfig& g);
    void set_grid(int nx, int ny, int nz, double dx, double dy, double dz, Vec3d origin);
    // Force recompute (blocking) and return snapshot
    FieldSnapshot recompute_now();
    // Non-blocking get; if dirty, returns last snapshot while background job may update
    FieldSnapshot get_snapshot() const;
    // Hook for live updates to ask server if a newer field is available
    bool has_newer(int known_version) const;

private:
    mutable std::mutex mtx_;
    GeometryConfig geo_;
    FieldSnapshot snapshot_;
    std::atomic<bool> dirty_{true};
    std::atomic<int> version_{0};

    void compute_locked_(); // called with mtx_ held
};
