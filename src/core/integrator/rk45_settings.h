// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
/**
 * RK45 settings header - small POD struct that can be included by both
 * the integrator implementation and global parameter definitions without
 * creating include cycles.
 */
#pragma once

struct RK45Settings {
    double absTol     = 1e-14;   ///< Absolute error tolerance
    double relTol     = 1e-12;   ///< Relative error tolerance
    double safety     = 0.84;    ///< Safety factor for step size adaptation
    double min_factor = 0.2;     ///< Minimum step size reduction factor
    double max_factor = 2.0;     ///< Maximum step size increase factor
    double dt_min     = 1e-12;   ///< Minimum allowed timestep [s]
    int max_rejects   = 1000;    ///< Maximum consecutive rejected steps before abort
};
