// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
/**
 * @file CompositeFieldProvider.h
 * @brief Field provider supporting superposition of multiple field arrays with time-varying scaling
 */
#pragma once

#include "IFieldProvider.h"
#include "core/io/fieldArrayLoader.h"
#include "core/config/types/FieldsConfig.h"
#include "core/config/types/DomainConfig.h"
#include <vector>
#include <memory>
#include <cmath>

/**
 * @class CompositeFieldProvider
 * @brief Combines multiple field arrays with time-dependent scaling factors
 * 
 * Implements field superposition:
 *   E_total(r, t) = Σ_i scale_i(t) · E_i(r)
 * 
 * where scale_i(t) depends on the term type:
 * - Constant:  scale = constant_V
 * - DC_Axial:  scale = DC.axial_V(t)
 * - DC_Quad:   scale = DC.quad_V(t)
 * - DC_Radial: scale = DC.radial_V(t)
 * - RF:        scale = RF.voltage_V(t) * cos(2π*f*t + φ)
 * 
 * Each field array typically represents a unit-voltage solution (1V applied)
 * from a BEM/FEM solver, which is then scaled by the actual voltage.
 */
class CompositeFieldProvider : public IFieldProvider {
public:
    /**
     * @brief Field array term with scaling information
     */
    struct FieldTerm {
        const FieldArray* field_array;           ///< Pointer to field data (must remain valid)
        ICARION::config::FieldsConfig::FieldArrayTerm::ScaleKind kind;  ///< Scaling type
        double constant_scale = 1.0;             ///< Constant scaling factor
        double frequency_Hz = 0.0;               ///< RF frequency (0 = use domain RF)
        double phase_rad = 0.0;                  ///< RF phase offset
    };
    
    /**
     * @brief Construct from field terms and domain configuration
     * 
     * @param terms Field array terms with scaling information
     * @param domain_config Domain configuration (for DC/RF voltages)
     */
    CompositeFieldProvider(
        std::vector<FieldTerm> terms,
        const ICARION::config::DomainConfig* domain_config
    ) : terms_(terms), domain_(domain_config) {
        if (!domain_config) {
            throw std::invalid_argument("CompositeFieldProvider: domain_config cannot be null");
        }
    }
    
    /**
     * @brief Evaluate electric field at position (time-independent interface)
     * 
     * Uses t=0 for time-varying terms. For accurate time-varying fields,
     * use the time-dependent get_E(pos, t) interface.
     */
    Vec3 get_E(const Vec3& pos) const override {
        return get_E(pos, 0.0);
    }
    
    /**
     * @brief Evaluate electric field at position and time
     * 
     * @param pos Position [m] in simulation domain
     * @param t Current simulation time [s]
     * @return Electric field E [V/m] from superposition
     */
    Vec3 get_E(const Vec3& pos, double t) const override {
        Vec3 E_total{0.0, 0.0, 0.0};
        
        for (const auto& term : terms_) {
            if (!term.field_array) continue;
            
            // Compute time-dependent scaling factor
            double scale = compute_scale_factor(term, t);
            
            // Interpolate field and add scaled contribution
            Vec3 E_i = interpolate_field(*term.field_array, pos);
            E_total.x += scale * E_i.x;
            E_total.y += scale * E_i.y;
            E_total.z += scale * E_i.z;
        }
        
        return E_total;
    }
    
    /**
     * @brief Get potential at position
     * @return 0.0 (potential reconstruction not implemented)
     */
    double get_phi(const Vec3& pos) const override {
        (void)pos;
        return 0.0;
    }
    
    /**
     * @brief Get number of field terms
     */
    size_t num_terms() const { return terms_.size(); }
    
private:
    /**
     * @brief Compute time-dependent scaling factor for a field term
     */
    double compute_scale_factor(const FieldTerm& term, double t) const {
        using ScaleKind = ICARION::config::FieldsConfig::FieldArrayTerm::ScaleKind;
        
        switch (term.kind) {
            case ScaleKind::Constant:
                return term.constant_scale;
            
            case ScaleKind::DC_Axial:
                return eval_voltage(domain_->fields.dc.axial_V, t);
            
            case ScaleKind::DC_Quad:
                return eval_voltage(domain_->fields.dc.quad_V, t);
            
            case ScaleKind::DC_Radial:
                return eval_voltage(domain_->fields.dc.radial_V, t);
            
            case ScaleKind::RF: {
                // RF scaling: V(t) * cos(ωt + φ)
                double V = eval_voltage(domain_->fields.rf.voltage_V, t);
                double f = (term.frequency_Hz > 0.0) ? term.frequency_Hz 
                         : eval_voltage(domain_->fields.rf.frequency_Hz, t);
                double omega = 2.0 * M_PI * f;
                return V * std::cos(omega * t + term.phase_rad);
            }
            
            default:
                return 0.0;
        }
    }
    
    /**
     * @brief Evaluate voltage from ValueOrWaveform at time t
     */
    double eval_voltage(const ICARION::config::ValueOrWaveform& value, double t) const {
        // If constant value is set, use it
        if (value.constant_value.has_value()) {
            return value.constant_value.value();
        }
        
        // If waveform reference exists, evaluate it
        if (value.waveform_ref.has_value()) {
            const auto& waveform_name = value.waveform_ref.value();
            
            // Look up waveform in domain's library
            const auto& lib = domain_->fields.waveform_library;
            auto it = lib.find(waveform_name);
            if (it != lib.end()) {
                return it->second.evaluate(t);
            }
            
            // Waveform not found - return 0
            return 0.0;
        }
        
        // No value set
        return 0.0;
    }
    
    std::vector<FieldTerm> terms_;                    ///< Field array terms
    const ICARION::config::DomainConfig* domain_;     ///< Domain config (for DC/RF voltages)
};
