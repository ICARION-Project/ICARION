// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include "core/types/IonEnsemble.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"
#include <optional>
#include <string>

namespace ICARION::physics {

// Forward declarations
struct ForceContext;

/**
 * @brief Lightweight SoA snapshot used for force evaluation
 */
struct ForceState {
  Vec3 pos{};
  Vec3 vel{};
  double mass_kg = 0.0;
  double ion_charge_C = 0.0;
  double CCS_m2 = 0.0;
  double reduced_mobility_cm2_Vs = 0.0;
  std::string species_id;
  bool active = true;
  bool born = true;
  uint32_t current_domain_index = 0;
  double birth_time_s = 0.0;
  std::optional<size_t> ensemble_index;

  IonState to_ion_state() const {
    IonState ion;
    ion.pos = pos;
    ion.vel = vel;
    ion.mass_kg = mass_kg;
    ion.ion_charge_C = ion_charge_C;
    ion.CCS_m2 = CCS_m2;
    ion.reduced_mobility_cm2_Vs = reduced_mobility_cm2_Vs;
    ion.species_id = species_id;
    ion.active = active;
    ion.born = born;
    ion.current_domain_index = current_domain_index;
    ion.birth_time_s = birth_time_s;
    return ion;
  }
};

/**
 * @brief Force contribution interface
 *
 * All forces (electric, magnetic, damping, space charge) implement this
 * interface. Enables modular force composition via ForceRegistry.
 *
 * Design Principles:
 * - Single Responsibility: Each force computes ONE physical contribution
 * - Const-correctness: Forces don't modify ion state (pure computation)
 * - Composability: Forces can be combined via ForceRegistry
 * - Testability: Each force can be unit tested in isolation
 *
 * Example Usage:
 * @code
 * // Create electric field force
 * auto e_force = std::make_unique<ElectricFieldForce>(domain, field_provider);
 *
 * // Compute force on ion
 * ForceContext ctx{};
 * Vec3 force = e_force->compute(ion, t, ctx);
 * @endcode
 *
 * @see ForceRegistry for force composition
 * @see ForceContext for shared context data
 */
class IForce {
public:
  virtual ~IForce() = default;

  /**
   * @brief Compute force contribution (SoA-only)
   *
   * @param ensemble Ion ensemble (SoA view)
   * @param ion_idx Index of ion in ensemble
   * @param t Current simulation time [s]
   * @param context Optional context for force computation (field provider,
   * domain, space charge, etc.)
   */
  virtual Vec3 compute(const core::IonEnsemble &ensemble, size_t ion_idx,
                       double t, const ForceContext &context) const = 0;

  /**
   * @brief Compute force from a SoA state snapshot (no AoS scratch)
   */
  virtual Vec3 compute_soa(const ForceState &state, double t,
                           const ForceContext &context) const = 0;

  /**
   * @brief Check if this force applies to given ion (AoS view for filtering
   * only)
   */
  virtual bool applies_to(const IonState &ion) const {
    (void)ion; // Unused parameter
    return true;
  }

  /**
   * @brief Indicates whether applies_to() depends on AoS state
   *
   * If false, ForceRegistry may skip IonState construction and applicability checks.
   */
  virtual bool requires_aos_state() const { return false; }

  /**
   * @brief Get force name for logging/debugging
   */
  virtual std::string name() const = 0;
};

} // namespace ICARION::physics
