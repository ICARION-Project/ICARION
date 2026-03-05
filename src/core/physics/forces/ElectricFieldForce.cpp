// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "ElectricFieldForce.h"
#include "fieldsolver/utils/IFieldProvider.h"
#include <stdexcept>

namespace ICARION {
namespace physics {

// ============================================================================
// Constructors
// ============================================================================

ElectricFieldForce::ElectricFieldForce(const ICARION::config::DomainConfig& domain)
    : use_field_provider_(false)
    , domain_(&domain)  // SSOT: store config reference
{
    // Build analytical fallback model from config for legacy contexts
    fallback_model_ = std::make_unique<ICARION::config::AnalyticalFieldModel>(domain);
}

ElectricFieldForce::ElectricFieldForce(std::shared_ptr<::IFieldProvider> field_provider)
    : use_field_provider_(true)
    , field_provider_(std::move(field_provider))
    , domain_(nullptr)  // Field provider mode: no config needed
{
    if (!field_provider_) {
        throw std::invalid_argument(
            "ElectricFieldForce: Field provider is null. "
            "Use analytical constructor if field provider unavailable."
        );
    }
}

// ============================================================================
// IForce Interface Implementation
// ============================================================================

Vec3 ElectricFieldForce::compute(const core::IonEnsemble& ensemble, size_t ion_idx, double t,
                                 const ForceContext& ctx) const {
    ForceState state{};
    state.pos = ensemble.get_pos(ion_idx);
    state.vel = ensemble.get_vel(ion_idx);
    state.mass_kg = ensemble.mass_data()[ion_idx];
    state.ion_charge_C = ensemble.charge_data()[ion_idx];
    state.CCS_m2 = ensemble.CCS(ion_idx);
    state.reduced_mobility_cm2_Vs = ensemble.mobility(ion_idx);
    state.species_id = ensemble.species_id(ion_idx);
    state.active = ensemble.active_data()[ion_idx] != 0;
    state.born = ensemble.born_data()[ion_idx] != 0;
    state.current_domain_index = ensemble.domain_index(ion_idx);
    state.birth_time_s = ensemble.birth_time(ion_idx);
    state.ensemble_index = ion_idx;

    return compute_soa(state, t, ctx);
}

Vec3 ElectricFieldForce::compute_soa(const ForceState& state, double t,
                                     const ForceContext& ctx) const {
    Vec3 E_field{0.0, 0.0, 0.0};

    // Prefer SSOT field model; fall back to providers, then internal model.
    if (ctx.field_model) {
        E_field = ctx.field_model->E(state.pos, t);
    } else if (ctx.field_provider) {
        E_field = ctx.field_provider->get_E(state.pos, t);
    } else if (use_field_provider_ && field_provider_) {
        E_field = field_provider_->get_E(state.pos, t);
    } else if (fallback_model_) {
        E_field = fallback_model_->E(state.pos, t);
    }

    // F = q * E
    return E_field * state.ion_charge_C;
}

std::string ElectricFieldForce::name() const {
    if (use_field_provider_) {
        return "ElectricField(FieldProvider)";
    }
    
    if (!domain_) {
        return "ElectricField(NoConfig)";
    }
    
    using Inst = ICARION::config::Instrument;
    switch (domain_->instrument) {
        case Inst::LQIT:              return "ElectricField(LQIT)";
        case Inst::IMS:               return "ElectricField(IMS)";
        case Inst::TOF:               return "ElectricField(TOF)";
        case Inst::FTICR:             return "ElectricField(FTICR)";
        case Inst::Orbitrap:          return "ElectricField(Orbitrap)";
        case Inst::QuadrupoleRF:      return "ElectricField(QuadrupoleRF)";
        case Inst::NoFixedInstrument: return "ElectricField(NoFixedInstrument)";
        case Inst::UnknownInstrument: return "ElectricField(Unknown)";
        default:                      return "ElectricField(Unknown)";
    }
}

} // namespace physics
} // namespace ICARION
