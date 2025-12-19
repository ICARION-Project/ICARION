// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "ForceRegistry.h"

namespace ICARION::physics {

void ForceRegistry::add_force(std::unique_ptr<IForce> force) {
    if (force) {  // Null-check (defensive programming)
        forces_.push_back(std::move(force));
    }
}

void ForceRegistry::set_space_charge_model(SpaceChargeModelPtr model) {
    space_charge_model_ = std::move(model);
}

Vec3 ForceRegistry::compute_total_force(
    const core::IonEnsemble& ensemble,
    size_t ion_idx,
    double t,
    const ForceContext& context
) const {
    Vec3 total_force{0, 0, 0};

    if (!forces_.empty()) {
        // Build a minimal IonState once for applies_to checks
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

        for (const auto& force : forces_) {
            if (force->applies_to(ion)) {
                total_force += force->compute(ensemble, ion_idx, t, context);
            }
        }
    }

    if (space_charge_model_) {
        Vec3 E = space_charge_model_->sample_electric_field(ion_idx);
        total_force += E * ensemble.charge_data()[ion_idx];
    }
    return total_force;
}

Vec3 ForceRegistry::compute_total_force(
    const IonState& ion,
    double t,
    const ForceContext& context
) const {
    Vec3 total_force{0, 0, 0};

    if (!forces_.empty()) {
        for (const auto& force : forces_) {
            if (force->applies_to(ion)) {
                // For AoS path, we pass through the same context (no ensemble)
                ForceContext ctx = context;
                ctx.all_ions = nullptr;
                ctx.ion_ensemble = nullptr;
                total_force += force->compute(ion, t, ctx);
            }
        }
    }

    if (space_charge_model_) {
        // AoS path has no per-ion index; space-charge not supported here.
    }
    return total_force;
}

} // namespace ICARION::physics
