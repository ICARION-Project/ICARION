// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_GEOMETRY_CONFIG_H
#define ICARION_CONFIG_GEOMETRY_CONFIG_H

#include "core/utils/mathUtils.h"
#include "../validation/ValidationResult.h"
#include <stdexcept>

namespace ICARION::config {

/**
 * @brief Geometry configuration (input parameters)
 * 
 * Contains user-specified geometric dimensions.
 * Derived values (bounds, transforms) are computed at initialization.
 * Clean separation: input data only in this struct.
 */
struct GeometryConfig {
    // === Basic cylindrical dimensions ===
    double length_m = 0.0;              ///< Axial length [m]
    double radius_m = 0.0;              ///< Radial extent [m]
    
    // === Orbitrap-specific ===
    double radius_in_m = 0.0;           ///< Inner electrode radius [m]
    double radius_out_m = 0.0;          ///< Outer electrode radius [m]
    double radius_char_m = 0.0;         ///< Characteristic radius [m] (hyperlogarithmic geometry)
    
    // === TOF-specific ===
    double acc_length_m = 0.0;          ///< Acceleration region length [m] (TOF)

    // === Multi-domain specific ===
    double end_aperture_m = 0.0;        ///< Exit aperture radius [m] (multi-domain transitions)
    
    // === Spatial transform ===
    Vec3 origin_m = {0.0, 0.0, 0.0};    ///< Domain origin in global coordinates
    
    // === Derived values (computed after loading) ===
    Vec3 min_bound = {0.0, 0.0, 0.0};   ///< Bounding box minimum
    Vec3 max_bound = {0.0, 0.0, 0.0};   ///< Bounding box maximum
    
    /**
     * @brief Compute bounding box after loading geometry
     * 
     * Called by loader after parsing JSON.
     * Simple cylindrical approximation for now.
     */
    void compute_bounds() {
        // Basic cylindrical bounding box
        // For more complex geometries, this can be overridden
        // must be changed for Orbitrap! #todo
        double r = (radius_m > 0.0) ? radius_m : 
                   (radius_out_m > 0.0) ? radius_out_m : 0.01;
        double l = (length_m > 0.0) ? length_m : 0.01;
        
        min_bound = origin_m - Vec3{r, r, l/2.0};
        max_bound = origin_m + Vec3{r, r, l/2.0};
    }
    
    /**
     * @brief Validate geometry parameters
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        if (length_m < 0.0) {
            result.add_error("Geometry length_m cannot be negative");
        }
        if (radius_m < 0.0) {
            result.add_error("Geometry radius_m cannot be negative");
        }
        if (radius_in_m < 0.0) {
            result.add_error("Geometry radius_in_m cannot be negative");
        }
        if (radius_out_m < 0.0) {
            result.add_error("Geometry radius_out_m cannot be negative");
        }
        
        // Orbitrap consistency check
        if (radius_in_m > 0.0 && radius_out_m > 0.0) {
            if (radius_in_m >= radius_out_m) {
                result.add_error("Orbitrap: radius_in_m must be < radius_out_m");
            }
        }
        
        return result;
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_GEOMETRY_CONFIG_H
