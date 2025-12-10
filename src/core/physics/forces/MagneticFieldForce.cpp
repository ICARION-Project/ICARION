// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "MagneticFieldForce.h"
#include "core/utils/mathUtils.h"
#include "fieldsolver/utils/IFieldProvider.h"
#include "core/config/types/FieldsConfig.h"

#include <stdexcept>

namespace ICARION {
namespace physics {

// ============================================================================
// Constructors
// ============================================================================

MagneticFieldForce::MagneticFieldForce(const ICARION::config::MagneticFieldConfig& magnetic_config)
    : use_field_provider_(false)
    , magnetic_config_(&magnetic_config)  // SSOT: store config reference
{
    // No validation needed - config can be all zeros (disabled force)
}

MagneticFieldForce::MagneticFieldForce(std::shared_ptr<::IFieldProvider> field_provider)
    : use_field_provider_(true)
    , field_provider_(std::move(field_provider))
    , magnetic_config_(nullptr)  // Field provider mode: no config needed
{
    if (!field_provider_) {
        throw std::invalid_argument(
            "MagneticFieldForce: Field provider is null. "
            "Use analytical constructor if field provider unavailable."
        );
    }
}

// ============================================================================
// IForce Interface Implementation
// ============================================================================

Vec3 MagneticFieldForce::compute(const IonState& ion, double t, const ForceContext& ctx) const {
    (void)t;  // Magnetic field is time-independent (for now)
    
    // Check if force is enabled (SSOT: read from config)
    if (magnetic_config_ && !magnetic_config_->enabled) {
        return Vec3{0.0, 0.0, 0.0};
    }
    
    // Get B-field
    Vec3 B_field{0.0, 0.0, 0.0};
    
    // Priority: context field provider > constructor field provider > analytical
    if (ctx.field_provider) {
        // TODO: FieldProvider doesn't support magnetic fields yet
        // B_field = ctx.field_provider->get_B(ion.pos);
    } else if (use_field_provider_ && field_provider_) {
        // TODO: FieldProvider doesn't support magnetic fields yet
        // B_field = field_provider_->get_B(ion.pos);
    } else {
        B_field = compute_analytical_field(ion.pos);
    }
    
    // F = q * (v × B) - Lorentz force
    Vec3 v_cross_B = cross(ion.vel, B_field);
    return v_cross_B * ion.ion_charge_C;
}

Vec3 MagneticFieldForce::compute_batch(const core::IonEnsemble& ensemble, size_t ion_idx, double t,
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

std::string MagneticFieldForce::name() const {
    if (use_field_provider_) {
        return "MagneticField(FieldProvider)";
    }
    return "MagneticField";
}

// ============================================================================
// Analytical Field Calculation
// ============================================================================

Vec3 MagneticFieldForce::compute_analytical_field(const Vec3& pos) const {
    if (!magnetic_config_) {
        return Vec3{0.0, 0.0, 0.0};  // No config available
    }
    
    // SSOT: Read directly from config!
    // Uniform field component (primary)
    Vec3 B_field = magnetic_config_->field_strength_T;
    
    // Gradient component (optional): B(r) = B₀ + ∇B·r
    const Vec3& grad = magnetic_config_->field_gradient_T_m;
    if (grad.x != 0.0 || grad.y != 0.0 || grad.z != 0.0) {
        B_field.x += grad.x * pos.x;
        B_field.y += grad.y * pos.y;
        B_field.z += grad.z * pos.z;
    }
    
    return B_field;
}

} // namespace physics
} // namespace ICARION
