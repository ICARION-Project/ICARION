// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "ForceRegistry.h"

namespace ICARION::physics {

void ForceRegistry::add_force(std::unique_ptr<IForce> force) {
  if (force) { // Null-check (defensive programming)
    forces_.push_back(std::move(force));
  }
}

void ForceRegistry::set_space_charge_model(SpaceChargeModelPtr model) {
  space_charge_model_ = std::move(model);
}

Vec3 ForceRegistry::compute_total_force(const core::IonEnsemble &ensemble,
                                        size_t ion_idx, double t,
                                        const ForceContext &context) const {
  ForceState state;
  state.pos = ensemble.get_pos(ion_idx);
  state.vel = ensemble.get_vel(ion_idx);
  state.mass_kg = ensemble.mass_data()[ion_idx];
  state.ion_charge_C = ensemble.charge_data()[ion_idx];
  state.active = ensemble.active_data()[ion_idx] != 0;
  state.born = ensemble.born_data()[ion_idx] != 0;
  state.current_domain_index = ensemble.domain_index(ion_idx);
  state.CCS_m2 = ensemble.CCS(ion_idx);
  state.reduced_mobility_cm2_Vs = ensemble.mobility(ion_idx);
  state.species_id = ensemble.species_id(ion_idx);
  state.birth_time_s = ensemble.birth_time(ion_idx);
  state.ensemble_index = ion_idx;

  return compute_total_force_soa(state, t, context);
}

Vec3 ForceRegistry::compute_total_force_soa(const ForceState &state, double t,
                                            const ForceContext &context) const {
  Vec3 total_force{0, 0, 0};

  if (!forces_.empty()) {
    IonState ion = state.to_ion_state();

    for (const auto &force : forces_) {
      if (force->applies_to(ion)) {
        total_force += force->compute_soa(state, t, context);
      }
    }
  }

  if (space_charge_model_ && state.ensemble_index.has_value()) {
    Vec3 E = space_charge_model_->sample_electric_field(*state.ensemble_index);
    total_force += E * state.ion_charge_C;
  }
  return total_force;
}

} // namespace ICARION::physics
