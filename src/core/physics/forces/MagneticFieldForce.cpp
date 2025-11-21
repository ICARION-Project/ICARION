// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "MagneticFieldForce.h"
#include "core/utils/mathUtils.h"
#include "fieldsolver/utils/IFieldProvider.h"

#include <stdexcept>

namespace ICARION {
namespace physics {

// ============================================================================
// Constructors
// ============================================================================

MagneticFieldForce::MagneticFieldForce(const MagneticFieldParams& params)
    : use_field_provider_(false)
    , analytical_params_(params)
{
    // No validation needed - params can be all zeros (disabled force)
}

MagneticFieldForce::MagneticFieldForce(std::shared_ptr<::IFieldProvider> field_provider)
    : use_field_provider_(true)
    , field_provider_(std::move(field_provider))
{
    if (!field_provider_) {
        throw std::invalid_argument(
            "MagneticFieldForce: Field provider is null. "
            "Use analytical constructor if field provider unavailable."
        );
    }
    
    // Enable by default when using field provider
    analytical_params_.enabled = true;
}

// ============================================================================
// IForce Interface Implementation
// ============================================================================

Vec3 MagneticFieldForce::compute(const IonState& ion, double t, const ForceContext& ctx) const {
    (void)t;  // Magnetic field is time-independent (for now)
    
    // Check if force is enabled
    if (!analytical_params_.enabled) {
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
    // Linear model: B(r) = B₀ + ∇B·r
    return analytical_params_.uniform_field_T + Vec3{
        analytical_params_.gradient_T_per_m.x * pos.x,
        analytical_params_.gradient_T_per_m.y * pos.y,
        analytical_params_.gradient_T_per_m.z * pos.z
    };
}

} // namespace physics
} // namespace ICARION
