// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

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

Vec3 ElectricFieldForce::compute(const IonState& ion, double t, const ForceContext& ctx) const {
    Vec3 E_field{0.0, 0.0, 0.0};
    
    // Prefer SSOT field model; fall back to providers, then internal model.
    if (ctx.field_model) {
        E_field = ctx.field_model->E(ion.pos, t);
    } else if (ctx.field_provider) {
        E_field = ctx.field_provider->get_E(ion.pos, t);
    } else if (use_field_provider_ && field_provider_) {
        E_field = field_provider_->get_E(ion.pos, t);
    } else if (fallback_model_) {
        E_field = fallback_model_->E(ion.pos, t);
    }
    
    // F = q * E
    return E_field * ion.ion_charge_C;
}

Vec3 ElectricFieldForce::compute_soa(const core::IonEnsemble& ensemble, size_t ion_idx, double t,
                                     const ForceContext& ctx) const {
    IonState ion;
    ion.pos = ensemble.get_pos(ion_idx);
    ion.vel = ensemble.get_vel(ion_idx);
    ion.mass_kg = ensemble.mass_data()[ion_idx];
    ion.ion_charge_C = ensemble.charge_data()[ion_idx];
    ion.active = ensemble.active_data()[ion_idx] != 0;
    ion.born = ensemble.born_data()[ion_idx] != 0;
    ion.current_domain_index = ensemble.domain_index(ion_idx);
    ion.CCS_m2 = ensemble.CCS(ion_idx);
    ion.reduced_mobility_cm2_Vs = ensemble.mobility(ion_idx);
    ion.species_id = ensemble.species_id(ion_idx);
    ion.birth_time_s = ensemble.birth_time(ion_idx);
    return compute(ion, t, ctx);
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
