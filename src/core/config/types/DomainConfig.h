// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_DOMAIN_CONFIG_H
#define ICARION_CONFIG_DOMAIN_CONFIG_H

#include "GeometryConfig.h"
#include "EnvironmentConfig.h"
#include "FieldsConfig.h"
#include "SolverEnums.h"
#include "instrument/InstrumentTypes.h"
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
    
    // === Solver ===
    // FUTURE: Allow per-domain solver selection via SolverModule
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
    void validate() const {
        // Hard errors (blocking)
        if (name.empty()) {
            throw std::runtime_error("Domain name cannot be empty");
        }
        
        if (instrument == Instrument::UnknownInstrument) {
            throw std::runtime_error("Domain '" + name + "' has unknown instrument type");
        }
        
        // Validate sub-components
        geometry.validate();
        environment.validate();
        fields.validate();
        
        // Instrument-specific validation (warnings only for missing typical fields)
        switch (instrument) {
            case Instrument::LQIT:
                if (fields.rf.voltage_V == 0.0 && fields.rf.frequency_Hz == 0.0) {
                    std::cout << "[WARN] LQIT domain '" << name << "' has no RF field. "
                              << "This is unusual but may be intentional for testing." << std::endl;
                }
                break;
                
            case Instrument::Orbitrap:
                if (geometry.radius_in_m == 0.0 || geometry.radius_out_m == 0.0) {
                    throw std::runtime_error("Orbitrap domain '" + name + "' missing inner/outer radius");
                }
                break;
                
            case Instrument::FTICR:
                if (!fields.magnetic.enabled || 
                    (fields.magnetic.field_strength_T.x == 0.0 && 
                     fields.magnetic.field_strength_T.y == 0.0 && 
                     fields.magnetic.field_strength_T.z == 0.0)) {
                    std::cout << "[WARN] FT-ICR domain '" << name << "' has no/zero magnetic field. "
                              << "Cyclotron motion will not occur." << std::endl;
                }
                break;
                
            case Instrument::IMS:
                if (fields.dc.EN_Td == 0.0 && fields.dc.axial_V == 0.0) {
                    std::cout << "[WARN] IMS domain '" << name << "' has no drift field. "
                              << "Ion transport may rely solely on gas flow (e.g., SIFT mode)." << std::endl;
                }
                break;
                
            default:
                // Other instruments: no specific requirements
                break;
        }
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_DOMAIN_CONFIG_H
