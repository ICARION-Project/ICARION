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

Vec3 ElectricFieldForce::compute(const core::IonEnsemble& ensemble, size_t ion_idx, double t,
                                 const ForceContext& ctx) const {
    Vec3 E_field{0.0, 0.0, 0.0};

    const Vec3 pos = ensemble.get_pos(ion_idx);
    const double q = ensemble.charge_data()[ion_idx];

    // Prefer SSOT field model; fall back to providers, then internal model.
    if (ctx.field_model) {
        E_field = ctx.field_model->E(pos, t);
    } else if (ctx.field_provider) {
        E_field = ctx.field_provider->get_E(pos, t);
    } else if (use_field_provider_ && field_provider_) {
        E_field = field_provider_->get_E(pos, t);
    } else if (fallback_model_) {
        E_field = fallback_model_->E(pos, t);
    }

    // F = q * E
    return E_field * q;
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
