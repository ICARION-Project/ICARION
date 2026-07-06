// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#ifndef ICARION_CONFIG_DOMAIN_CONFIG_H
#define ICARION_CONFIG_DOMAIN_CONFIG_H

#include "GeometryConfig.h"
#include "EnvironmentConfig.h"
#include "FieldsConfig.h"
#include "BoundaryConfig.h"
#include "SolverEnums.h"
#include "core/config/types/InstrumentTypes.h"
#include "core/utils/mathUtils.h"
#include "../validation/ValidationResult.h"
#include <string>
#include <stdexcept>
#include <iostream>

namespace ICARION::config {

// Use canonical InstrumentType from instrument namespace
using Instrument = ICARION::instrument::InstrumentType;

/**
 * @brief Configuration for a single instrument domain
 * 
 * Replaces InstrumentDomain with clean separation of input vs. derived data.
 * All input parameters are loaded from JSON, derived quantities are computed
 * via finalize() after loading.
 */
struct DomainConfig {
    // === Identification ===
    std::string name = "domain";                        ///< User-defined domain name
    Instrument instrument = Instrument::UnknownInstrument;  ///< Instrument type
    
    // === Physical configuration ===
    GeometryConfig geometry;
    EnvironmentConfig environment;
    FieldsConfig fields;
    BoundaryConfig boundary;                            ///< Boundary action configuration
    
    // === Solver ===
    SolverType solver = SolverType::RK4;                ///< Integrator for this domain
    
    // === Coordinate transforms (for multi-domain) ===
    Mat3 rotation_global_to_local = Mat3::identity();
    Mat3 rotation_local_to_global = Mat3::identity();
    
    // === Runtime state (set after loading) ===
    int domain_index = -1;                              ///< Assigned by loader (0, 1, 2, ...)
    
    // === Live field computation (Future feature) ===
    bool use_live_field_server = false;                 ///< Use FieldServer instead of analytic
    
    /**
     * @brief Finalize domain configuration
     * 
     * Computes all derived quantities after loading from JSON:
     * - Geometry bounding boxes
     * - Environment thermodynamic properties
     * - Field angular frequencies
     * 
     * Must be called before using the domain in simulation.
     */
    void finalize() {
        geometry.compute_bounds();
        environment.compute_derived_properties();
        fields.compute_derived();
    }
    
    /**
     * @brief Validate complete domain configuration
     * 
     * @throws std::runtime_error if critically invalid (empty name, unknown instrument)
     * 
     * Prints warnings for unusual but potentially valid configurations
     * (e.g., LQIT without RF, IMS without drift field for SIFT-like simulations).
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        // Hard errors (blocking)
        if (name.empty()) {
            result.add_error("Domain name cannot be empty");
        }
        
        if (instrument == Instrument::UnknownInstrument) {
            result.add_error("Domain '" + name + "' has unknown instrument type");
        }
        
        // Validate sub-components
        result.merge(geometry.validate());
        result.merge(environment.validate());
        result.merge(fields.validate());
        result.merge(boundary.validate());
        
        // Instrument-specific validation (warnings only for missing typical fields)
        switch (instrument) {
            case Instrument::LQIT:
                // Check if RF voltage and frequency are both zero (only for static values)
                if (fields.rf.voltage_V.constant_value.has_value() && 
                    fields.rf.frequency_Hz.constant_value.has_value() &&
                    fields.rf.voltage_V.constant_value.value() == 0.0 && 
                    fields.rf.frequency_Hz.constant_value.value() == 0.0) {
                    result.add_warning("LQIT domain '" + name + "' has no RF field. "
                                      "This is unusual but may be intentional for testing.");
                }
                break;
                
            case Instrument::Orbitrap:
                if (geometry.radius_in_m == 0.0 || geometry.radius_out_m == 0.0) {
                    result.add_error("Orbitrap domain '" + name + "' missing inner/outer radius");
                }
                break;
                
            case Instrument::FTICR:
                if (!fields.magnetic.enabled || 
                    (fields.magnetic.field_strength_T.x == 0.0 && 
                     fields.magnetic.field_strength_T.y == 0.0 && 
                     fields.magnetic.field_strength_T.z == 0.0)) {
                    result.add_warning("FT-ICR domain '" + name + "' has no/zero magnetic field. "
                                      "Cyclotron motion will not occur.");
                }
                break;
                
            case Instrument::IMS:
            case Instrument::TIMS:
                // Check if EN_Td and axial_V are both zero (only for static values)
                if (fields.dc.EN_Td.constant_value.has_value() && 
                    fields.dc.axial_V.constant_value.has_value() &&
                    fields.dc.EN_Td.constant_value.value() == 0.0 && 
                    fields.dc.axial_V.constant_value.value() == 0.0 &&
                    !fields.tims.enabled) {
                    result.add_warning("IMS/TIMS domain '" + name + "' has no drift field. "
                                      "Ion transport may rely solely on gas flow.");
                }
                break;
                
            default:
                // Other instruments: no specific requirements
                break;
        }
        
        return result;
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_DOMAIN_CONFIG_H
