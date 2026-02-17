// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
/**
 * @file CollisionTypes.h
 * @brief Common types for collision physics (RNG, parameters)
 * 
 * Defines shared types used across collision modules:
 * - PhysicsRng: Random number generator wrapper
 * - EhssRng: Legacy alias kept for backward compatibility
 * - EHSSParams: Legacy parameter struct (for backwards compatibility)
 */

#pragma once

#include <cstdint>
#include <random>

namespace ICARION::physics {

/**
 * @brief Random number generator for physics simulations
 * 
 * Wraps std::mt19937_64 with uniform and normal distributions.
 * Provides reproducible random numbers for collision, reaction, and stochastic physics.
 * 
 * @note Previously named EhssRng (legacy EHSS-specific name).
 *       Now used globally for all stochastic physics processes.
 *       Prefer PhysicsRng in new code; EhssRng remains as an alias only.
 */
class PhysicsRng {
public:
    explicit PhysicsRng(uint64_t seed = 0xDEADBEEFCAFEBABEULL);
    
    /**
     * @brief Generate uniform random number in [0, 1)
     */
    double uniform01();
    
    /**
     * @brief Generate standard normal random number (mean=0, stddev=1)
     */
    double normal();

private:
    std::mt19937_64 eng_;
    std::uniform_real_distribution<double> uni_;
    std::normal_distribution<double> norm_;
};

// Legacy alias for backwards compatibility (prefer PhysicsRng in new code)
using EhssRng = PhysicsRng;

/**
 * @brief Parameters for EHSS collision models (LEGACY)
 * 
 * @deprecated This struct is kept for backwards compatibility with old test code.
 *             New code should pass parameters directly to collision functions.
 */
struct EHSSParams {
    double n, dt, mi, mn, kB, Tn, sigma_eff;
    double ubx, uby, ubz;
    double Rn;
    int num_atoms = 0;
    double max_extent = 0.0;
    int max_attempts = 256;
};

} // namespace ICARION::physics
