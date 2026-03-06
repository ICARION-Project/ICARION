// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include <string>
#include <vector>
#include <optional>
#include "core/utils/mathUtils.h"

namespace ICARION::config {

/**
 * @brief Position sampling distribution types
 */
enum class PositionDistribution {
    Point,           ///< Single point
    Gaussian,        ///< Normal distribution (3D independent)
    UniformSphere,   ///< Uniform in sphere
    UniformCylinder, ///< Uniform in cylinder (z-axis)
    UniformBox       ///< Uniform in box
};

/**
 * @brief Velocity sampling distribution types
 */
enum class VelocityDistribution {
    Fixed,      ///< All ions same velocity
    Thermal,    ///< Maxwell-Boltzmann (random directions)
    Kinetic,    ///< Fixed energy, directed
    Gaussian    ///< Normal distribution (3D independent)
};

/**
 * @brief Position distribution configuration
 * 
 * Each species can have its own position boundaries/distribution.
 */
struct PositionConfig {
    PositionDistribution type = PositionDistribution::Point;
    
    // Point / Gaussian mean / Sphere center / Cylinder center
    Vec3 center = {0, 0, 0};  // [m]
    
    // Gaussian: std per axis [m]
    Vec3 std_dev = {0, 0, 0};
    
    // Sphere: radius [m]
    double radius = 0.0;
    
    // Cylinder: radius [m], length [m] (along z-axis)
    double cylinder_radius = 0.0;
    double cylinder_length = 0.0;
    
    // Box: min/max corners [m]
    Vec3 box_min = {0, 0, 0};
    Vec3 box_max = {0, 0, 0};
};

/**
 * @brief Velocity distribution configuration
 */
struct VelocityConfig {
    VelocityDistribution type = VelocityDistribution::Fixed;
    
    // Fixed: velocity [m/s]
    Vec3 value = {0, 0, 0};
    
    // Thermal: temperature [K] (directions sampled randomly from MB distribution)
    double temperature_K = 0.0;
    
    // Kinetic: energy [eV], direction (unit vector), spread angle [deg]
    double energy_eV = 0.0;
    Vec3 direction = {0, 0, 1};
    double spread_angle_deg = 0.0;
    
    // Gaussian: mean [m/s], std [m/s]
    Vec3 mean = {0, 0, 0};
    Vec3 std_dev = {0, 0, 0};
};

/**
 * @brief Single ion species configuration
 * 
 * Each species has its own:
 * - Position boundaries (where ions are created)
 * - Velocity distribution
 * - Count and birth time
 */
struct IonSpeciesConfig {
    std::string species_id;           ///< Species ID from database (e.g., "H3O+")
    size_t count = 0;                 ///< Number of ions to generate
    PositionConfig position;          ///< Position sampling (boundaries for this species)
    VelocityConfig velocity;          ///< Velocity sampling
    double birth_time_s = 0.0;        ///< Birth time [s] (0 = born at start)
    bool use_birth_time_distribution = false;  ///< If true, sample birth time per ion
    double birth_time_min_s = 0.0;    ///< Truncation lower bound for sampled birth time [s]
    double birth_time_max_s = 0.0;    ///< Truncation upper bound for sampled birth time [s]
    double birth_time_mean_s = 0.0;   ///< Mean for Gaussian birth-time sampling [s]
    double birth_time_std_s = 0.0;    ///< Std for Gaussian birth-time sampling [s]
};

/**
 * @brief Complete ion configuration
 */
struct IonConfig {
    // Option 1: Load from file
    std::optional<std::string> from_file;  ///< Path to ion cloud JSON
    
    // Option 2: Generate from species list (each with own boundaries)
    std::vector<IonSpeciesConfig> species;
    
    /**
     * @brief Check if valid
     */
    bool is_valid() const {
        return from_file.has_value() || !species.empty();
    }
};

/**
 * @brief Convert enum to string for logging
 */
inline const char* to_string(PositionDistribution dist) {
    switch (dist) {
        case PositionDistribution::Point: return "Point";
        case PositionDistribution::Gaussian: return "Gaussian";
        case PositionDistribution::UniformSphere: return "UniformSphere";
        case PositionDistribution::UniformCylinder: return "UniformCylinder";
        case PositionDistribution::UniformBox: return "UniformBox";
        default: return "Unknown";
    }
}

inline const char* to_string(VelocityDistribution dist) {
    switch (dist) {
        case VelocityDistribution::Fixed: return "Fixed";
        case VelocityDistribution::Thermal: return "Thermal";
        case VelocityDistribution::Kinetic: return "Kinetic";
        case VelocityDistribution::Gaussian: return "Gaussian";
        default: return "Unknown";
    }
}

} // namespace ICARION::config
