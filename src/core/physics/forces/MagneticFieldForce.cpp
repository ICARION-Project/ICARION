// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

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

MagneticFieldForce::MagneticFieldForce(const config::MagneticFieldConfig& magnetic_config)
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
        B_field = ctx.field_provider->get_E(ion.pos);  // get_E() used for B-field
    } else if (use_field_provider_ && field_provider_) {
        B_field = field_provider_->get_E(ion.pos);
    } else {
        B_field = compute_analytical_field(ion.pos);
    }
    
    // F = q * (v × B) - Lorentz force
    Vec3 v_cross_B = cross(ion.vel, B_field);
    return v_cross_B * ion.ion_charge_C;
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
    Vec3 B_field{0.0, 0.0, magnetic_config_->field_strength_T};
    
    // Gradient component (optional): B(r) = B₀ + ∇B·r
    if (magnetic_config_->gradient_T_m != 0.0) {
        B_field.z += magnetic_config_->gradient_T_m * pos.z;
    }
    
    return B_field;
}

} // namespace physics
} // namespace ICARION
