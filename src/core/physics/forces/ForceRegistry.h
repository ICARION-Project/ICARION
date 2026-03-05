// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include "ForceContext.h"
#include "IForce.h"
#include "core/config/types/DomainConfig.h"
#include "core/physics/spacecharge/ISpaceChargeModel.h"
#include "core/types/IonEnsemble.h"
#include <memory>
#include <vector>

namespace ICARION::config {
class IFieldModel;
}

namespace ICARION::physics {

/**
 * @brief Manages and aggregates all active forces for simulation
 *
 * Manages a set of force contributors and aggregates their outputs.
 *
 * Typical Usage:
 * @code
 * // Setup forces
 * ForceRegistry registry;
 * registry.add_force(std::make_unique<ElectricFieldForce>(domain));
 * registry.add_force(std::make_unique<MagneticFieldForce>(mag_config));
 * registry.add_force(std::make_unique<DampingForce>(gamma));
 *
 * // In integration loop
 * ForceContext ctx{field_provider, domain, all_ions};
 * Vec3 total_force = registry.compute_total_force(ion, t, ctx);
 * Vec3 acceleration = total_force / ion.mass_kg;
 * @endcode
 *
 * Thread Safety:
 * - compute_total_force() is const and thread-safe (read-only)
 * - add_force() / clear() are NOT thread-safe (call during setup only)
 *
 * @see IForce for force interface
 * @see ForceContext for shared context data
 */
class ForceRegistry {
public:
  /**
   * @brief Construct registry with domain context (RECOMMENDED)
   *
   * @param domain Domain configuration (geometry, fields, environment)
   */
  explicit ForceRegistry(const config::DomainConfig &domain)
      : domain_(&domain) {}

  /**
   * @brief Add force to registry
   *
   * Transfers ownership of the force to the registry.
   * Forces are evaluated in the order they are added.
   *
   * @param force Force implementation (ownership transferred)
   *
   * Example:
   * @code
   * registry.add_force(std::make_unique<ElectricFieldForce>(domain));
   * @endcode
   *
   * @note NOT thread-safe (call during setup phase only)
   */
  void add_force(std::unique_ptr<IForce> force);

  /**
   * @brief Compute total force using SoA ensemble data (primary)
   */
  Vec3 compute_total_force(const core::IonEnsemble &ensemble, size_t ion_idx,
                           double t, const ForceContext &context = {}) const;

  /**
   * @brief Compute total force from a SoA snapshot (no AoS scratch)
   */
  Vec3 compute_total_force_soa(const ForceState &state, double t,
                               const ForceContext &context = {}) const;

  /**
   * @brief Get all registered forces (for inspection/debugging)
   *
   * @return Vector of force pointers (const reference, no ownership transfer)
   */
  const std::vector<std::unique_ptr<IForce>> &forces() const { return forces_; }

  /**
   * @brief Get number of registered forces
   *
   * @return Number of active forces
   */
  size_t size() const { return forces_.size(); }

  /**
   * @brief Check if registry is empty
   *
   * @return true if no forces registered
   */
  bool empty() const { return forces_.empty(); }

  /**
   * @brief Clear all forces (removes all registered forces)
   *
   * Useful for resetting simulation or switching force configurations.
   *
   * @note NOT thread-safe (call during setup phase only)
   */
  void clear() {
    forces_.clear();
    uses_aos_filter_ = false;
  }

  /**
   * @brief Get domain configuration (if available)
   *
   * @return Pointer to domain config, or nullptr if not set
   *
   * Allows forces to access domain context (geometry, fields, environment).
   */
  const config::DomainConfig *domain() const { return domain_; }

  /**
   * @brief Rebind domain pointer (non-owning)
   *
   * Used when a copied configuration becomes the active SSOT.
   */
  void set_domain(const config::DomainConfig *domain) { domain_ = domain; }

  /**
   * @brief Set optional field model (non-owning)
   *
   * Allows integrators to pass a constructed IFieldModel via ForceContext.
   * Caller must ensure lifetime exceeds registry usage.
   */
  void set_field_model(const config::IFieldModel *model) {
    field_model_ = model;
  }

  /**
   * @brief Access configured field model (if any)
   */
  const config::IFieldModel *field_model() const { return field_model_; }

  /**
   * @brief Assign shared space-charge model (optional).
   *
   * Registries can share the same model (e.g., global direct Coulomb solver).
   */
  void set_space_charge_model(SpaceChargeModelPtr model);

  /**
   * @brief Access space-charge model (may be shared).
   */
  ISpaceChargeModel *space_charge_model() const {
    return space_charge_model_.get();
  }

private:
  /**
   * @brief Vector of registered forces (owned by registry)
   */
  std::vector<std::unique_ptr<IForce>> forces_;

  /**
   * @brief Domain configuration (optional, nullptr if not set)
   *
   * Non-owning pointer (registry does not own domain config).
   * Must remain valid for lifetime of registry.
   */
  const config::DomainConfig *domain_ = nullptr;

  /**
   * @brief Optional field model (non-owning)
   *
   * Set by setup layer once field models are constructed.
   */
  const config::IFieldModel *field_model_ = nullptr;

  /**
   * @brief Optional shared space-charge model.
   */
  SpaceChargeModelPtr space_charge_model_;

  /**
   * @brief Track whether any force requires AoS applicability checks.
   */
  bool uses_aos_filter_ = false;
};

} // namespace ICARION::physics
