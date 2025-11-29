// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "ForceRegistry.h"

namespace ICARION::physics {

void ForceRegistry::add_force(std::unique_ptr<IForce> force) {
    if (force) {  // Null-check (defensive programming)
        forces_.push_back(std::move(force));
    }
}

Vec3 ForceRegistry::compute_total_force(
    const IonState& ion,
    double t,
    const ForceContext& context
) const {
    // Early exit if no forces
    if (forces_.empty()) {
        return Vec3{0, 0, 0};
    }
    
    // Accumulate forces
    Vec3 total_force{0, 0, 0};
    
    for (const auto& force : forces_) {
        // Check if force applies (allows conditional forces)
        if (force->applies_to(ion)) {
            total_force += force->compute(ion, t, context);
        }
    }
    
    return total_force;
}

} // namespace ICARION::physics
