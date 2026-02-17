// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <iostream>
#include "core/types/IonState.h"

namespace ICARION {
namespace safety {

// Bring Vec3 and IonState into scope
using ICARION::core::Vec3;
using ICARION::core::IonState;

/**
 * @brief Configuration for numerical safety checks (CPU-only helper)
 */
struct NumericalSafetyConfig {
    bool enable_nan_checks = true;         ///< Enable NaN/Inf detection
    bool enable_bounds_checks = true;      ///< Enable position/velocity bounds checking
    bool enable_logging = false;           ///< Log safety violations (performance impact; avoid in tight loops)
    
    // Bounds for reasonable values
    double max_position_m = 1.0;           ///< Maximum position magnitude (m)
    double max_velocity_ms = 1e6;          ///< Maximum velocity magnitude (m/s)
    double max_acceleration_ms2 = 1e12;    ///< Maximum acceleration magnitude (m/s²)
    
    // Error handling behavior
    bool throw_on_violation = true;        ///< Throw exception on safety violation
    bool attempt_recovery = false;         ///< Attempt to recover invalid values (best-effort, coarse)
};

/**
 * @brief Check if a single value is finite (not NaN or infinite)
 */
inline bool is_finite_value(double value) {
    return std::isfinite(value) && !std::isnan(value);
}

/**
 * @brief Check if a Vec3 is finite (no NaN/Inf components)
 */
inline bool is_finite_vec3(const Vec3& v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

// Alias for compatibility
inline bool is_finite(const Vec3& v) {
    return is_finite_vec3(v);
}

/**
 * @brief Check if a scalar value is finite
 */
inline bool is_finite(double val) {
    return std::isfinite(val);
}

/**
 * @brief Check if an IonState has finite position and velocity
 */
inline bool is_finite(const IonState& ion) {
    return is_finite_vec3(ion.pos) && is_finite_vec3(ion.vel);
}

/**
 * @brief Check if a value is within reasonable bounds
 */
inline bool is_within_bounds(double value, double max_magnitude) {
    return std::abs(value) <= max_magnitude;
}

/**
 * @brief Check if a Vec3 is within reasonable bounds
 */
template<typename Vec3Type>
inline bool is_vec3_within_bounds(const Vec3Type& vec, double max_magnitude) {
    double magnitude = std::sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
    return magnitude <= max_magnitude;
}

/**
 * @brief Comprehensive safety check for ion state (CPU helper; not on hot paths)
 * 
 * @param ion Ion state to check
 * @param config Safety configuration
 * @param context Context string for error messages
 * @param step Current integration step
 * @param ion_index Ion index for identification
 * @return true if ion state is safe, false if violations detected
 */
template<typename IonStateType>
bool check_ion_safety(const IonStateType& ion, 
                      const NumericalSafetyConfig& config,
                      const std::string& context = "integration",
                      int step = -1,
                      int ion_index = -1) {
    
    std::string error_prefix = "Numerical safety violation";
    if (ion_index >= 0) {
        error_prefix += " for ion " + std::to_string(ion_index);
    }
    if (step >= 0) {
        error_prefix += " at step " + std::to_string(step);
    }
    error_prefix += " in " + context + ": ";
    
    // Check for NaN/Inf values
    if (config.enable_nan_checks) {
        if (!is_finite_vec3(ion.pos)) {
            std::string msg = error_prefix + "Non-finite position detected (x=" + 
                            std::to_string(ion.pos.x) + ", y=" + std::to_string(ion.pos.y) + 
                            ", z=" + std::to_string(ion.pos.z) + ")";
            
            if (config.enable_logging) {
                std::cerr << "WARNING: " << msg << std::endl;
            }
            
            if (config.throw_on_violation) {
                throw std::runtime_error(msg);
            }
            return false;
        }
        
        if (!is_finite_vec3(ion.vel)) {
            std::string msg = error_prefix + "Non-finite velocity detected (vx=" + 
                            std::to_string(ion.vel.x) + ", vy=" + std::to_string(ion.vel.y) + 
                            ", vz=" + std::to_string(ion.vel.z) + ")";
            
            if (config.enable_logging) {
                std::cerr << "WARNING: " << msg << std::endl;
            }
            
            if (config.throw_on_violation) {
                throw std::runtime_error(msg);
            }
            return false;
        }
        
        // Check time value
        if (!is_finite_value(ion.t)) {
            std::string msg = error_prefix + "Non-finite time detected (t=" + std::to_string(ion.t) + ")";
            
            if (config.enable_logging) {
                std::cerr << "WARNING: " << msg << std::endl;
            }
            
            if (config.throw_on_violation) {
                throw std::runtime_error(msg);
            }
            return false;
        }
    }
    
    // Check bounds
    if (config.enable_bounds_checks) {
        if (!is_vec3_within_bounds(ion.pos, config.max_position_m)) {
            double pos_magnitude = std::sqrt(ion.pos.x*ion.pos.x + ion.pos.y*ion.pos.y + ion.pos.z*ion.pos.z);
            std::string msg = error_prefix + "Position magnitude exceeds bounds (" + 
                            std::to_string(pos_magnitude) + " m > " + 
                            std::to_string(config.max_position_m) + " m)";
            
            if (config.enable_logging) {
                std::cerr << "WARNING: " << msg << std::endl;
            }
            
            // Position bounds violation is not treated as fatal here; caller may still continue.\n+            return true;
        }
        
        if (!is_vec3_within_bounds(ion.vel, config.max_velocity_ms)) {
            double vel_magnitude = std::sqrt(ion.vel.x*ion.vel.x + ion.vel.y*ion.vel.y + ion.vel.z*ion.vel.z);
            std::string msg = error_prefix + "Velocity magnitude exceeds bounds (" + 
                            std::to_string(vel_magnitude) + " m/s > " + 
                            std::to_string(config.max_velocity_ms) + " m/s)";
            
            if (config.enable_logging) {
                std::cerr << "WARNING: " << msg << std::endl;
            }
            
            if (config.throw_on_violation) {
                throw std::runtime_error(msg);
            }
            return false;
        }
    }
    
    return true;
}

/**
 * @brief Attempt to recover an ion with invalid values
 * 
 * @param ion Ion to recover
 * @param backup_ion Previous valid ion state for recovery
 * @return true if recovery successful
 */
template<typename IonStateType>
bool attempt_ion_recovery(IonStateType& ion, const IonStateType& backup_ion) {
    bool recovered = false;
    
    // Recover non-finite positions
    if (!is_finite_vec3(ion.pos)) {
        ion.pos = backup_ion.pos;
        recovered = true;
    }
    
    // Recover non-finite velocities
    if (!is_finite_vec3(ion.vel)) {
        ion.vel = backup_ion.vel * 0.1; // Reduce velocity to prevent immediate re-occurrence
        recovered = true;
    }
    
    // Recover non-finite time
    if (!is_finite_value(ion.t)) {
        ion.t = backup_ion.t;
        recovered = true;
    }
    
    return recovered;
}

/**
 * @brief Debug macro for finite value checking
 * 
 * Usage: ICARION_CHECK_FINITE(ion.pos, "After RK4 step");
 */
#ifdef ICARION_DEBUG
#define ICARION_CHECK_FINITE(vec, context) \
    do { \
        if (!ICARION::safety::is_finite_vec3(vec)) { \
            std::cerr << "DEBUG: Non-finite value detected in " << context \
                      << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        } \
    } while(0)
#else
#define ICARION_CHECK_FINITE(vec, context) do {} while(0)
#endif

} // namespace safety
} // namespace ICARION
