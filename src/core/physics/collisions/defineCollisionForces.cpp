/**
 * =====================================================================
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        defineCollisionForces.cpp
 *   @brief       Computes collisional damping forces acting on ions.
 *
 *   @details
 *   Implements multiple collision models including:
 *    - Friction (mobility-based damping)
 *    - Langevin (ion-neutral long-range interactions)
 *    - Hard-sphere collisions
 *    - EHSS / HSS stochastic models (no explicit damping)
 *
 *   Each function returns a force vector divided by the ion mass (acceleration).
 *
 *   ⚠️ **DEPRECATED (Phase 2B - November 2025):**
 *   This file duplicates logic already implemented in `DampingForce`.
 *   
 *   **Reason for duplication:**
 *   - This file is used by legacy `compute_accelerations()` (computeAccelerations.cpp)
 *   - DampingForce is used by new `ForceRegistry` architecture
 *   
 *   **Migration plan (Phase 4):**
 *   - Replace `compute_accelerations()` with `ForceRegistry::compute_total_force()`
 *   - Delete this file after integrator refactoring is complete
 *   
 *   **Affected models:**
 *   - HardSphereCollision() → DampingForce (DampingModel::HardSphere)
 *   - LangevinCollision()   → DampingForce (DampingModel::Langevin)
 *   - FrictionCollision()   → DampingForce (DampingModel::Friction)
 *   - EHSS/HSS              → No force (stochastic collisions only)
 *
 *   @date        2025-10-08
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#include "defineCollisionForces.h"
#include <cmath>
#include <algorithm>

namespace ICARION {
namespace physics {

/**
 * @brief Dispatches the collision model to the appropriate function.
 *
 * @param[in] ion     Current ion state.
 * @param[in] gParams Global simulation parameters including selected collision model.
 * @param[in] dom     Instrument domain properties.
 * @return Vec3 Acceleration vector due to collisional damping [m/s²].
 */
Vec3 CollisionForce(const IonState& ion, const GlobalParams& gParams, const InstrumentDomain& dom) {
    switch(gParams.collisionModel) {
        case CollisionModel::HSD: return HardSphereCollision(ion, dom);
        case CollisionModel::Langevin:   return LangevinCollision(ion, dom);
        case CollisionModel::Friction:   return FrictionCollision(ion, dom);
        case CollisionModel::EHSS:       return EHSSCollision(ion, dom);
        case CollisionModel::HSS:        return HSSCollision(ion, dom);  // Note: function name stays for now
        case CollisionModel::NoCollisions: return Vec3{0.0, 0.0, 0.0}; // No collisions
    }
    return Vec3{0.0, 0.0, 0.0}; // fallback
}

/**
 * @brief Hard-sphere collision model: damping from elastic collisions.
 *
 * Damping proportional to collision frequency derived from ion CCS,
 * neutral particle density, and thermal velocity.
 *
 * ⚠️ **DEPRECATED:** Duplicate of DampingForce::calculate_gamma(DampingModel::HardSphere)
 * Will be removed in Phase 4 (Integrator Refactor).
 *
 * @param[in] ion Current ion state.
 * @param[in] dom Instrument domain properties.
 * @return Vec3 Acceleration vector [m/s²].
 */
Vec3 HardSphereCollision(const IonState& ion, const InstrumentDomain& dom) {
    double momentum_transfer_rate = dom.env.particle_density_m_3 * ion.CCS_m2 * dom.env.mean_thermal_velocity_m_s
                  * ion.mass_kg / (ion.domain_neutral_mass_kg + ion.mass_kg);
    double damping = -ion.mass_kg * momentum_transfer_rate;
    return ion.vel * damping;
}

/**
 * @brief Langevin collision model: damping from long-range ion-neutral polarization.
 *
 * Computes velocity-dependent damping proportional to the Langevin collision frequency.
 *
 * ⚠️ **DEPRECATED:** Duplicate of DampingForce::calculate_gamma(DampingModel::Langevin)
 * Will be removed in Phase 4 (Integrator Refactor).
 *
 * @param[in] ion Current ion state.
 * @param[in] dom Instrument domain properties.
 * @return Vec3 Acceleration vector [m/s²].
 */
Vec3 LangevinCollision(const IonState& ion, const InstrumentDomain& dom) {
    double v_mag = std::max(1e-6, norm(ion.vel));
    double cs = M_PI * ion.ion_charge_C * sqrt(ion.domain_neutral_polarizability_m3 /
                (pow(4*M_PI*EPSILON_0,1.0) * (ion.mass_kg * ion.domain_neutral_mass_kg)/(ion.mass_kg + ion.domain_neutral_mass_kg)))
                / v_mag;
    double freq = dom.env.particle_density_m_3 * cs * dom.env.mean_thermal_velocity_m_s * ion.domain_neutral_mass_kg
                  / (ion.domain_neutral_mass_kg + ion.mass_kg);
    double damping = -ion.mass_kg * freq;
    return ion.vel * damping;
}

/**
 * @brief Friction collision model.
 *
 * Computes damping based on ion mobility. The force is opposite to ion velocity.
 *
 * ⚠️ **DEPRECATED:** Duplicate of DampingForce::calculate_gamma(DampingModel::Friction)
 * Will be removed in Phase 4 (Integrator Refactor).
 *
 * @param[in] ion  Current ion state.
 * @param[in] dom  Instrument domain parameters.
 * @return Vec3    Damping force vector [N].
 */
Vec3 FrictionCollision(const IonState& ion, const InstrumentDomain& dom) {
    double ion_mobility = (ion.reduced_mobility_cm2_Vs * 1e-4) * LOSCHMIDT_CONSTANT / dom.env.particle_density_m_3;
    double damping = -ion.ion_charge_C / ion_mobility;
    return ion.vel * damping;
}

/**
 * @brief EHSS collision model: returns zero damping force.
 * 
 * EHSS (Explicit Hard-Sphere Scattering) is a **stochastic** collision model.
 * Collisions are handled discretely by ICollisionHandler (Phase 2C), not as continuous force.
 * 
 * ⚠️ **NOTE:** This function will be removed in Phase 4 when EHSS is fully migrated
 * to ICollisionHandler architecture.
 */
Vec3 EHSSCollision(const IonState& ion, const InstrumentDomain& dom) {
    return Vec3{0.0, 0.0, 0.0};
}

/**
 * @brief HSS collision model: returns zero damping force.
 * 
 * HSS (Hard-Sphere Stochastic) is a **stochastic** collision model.
 * Collisions are handled discretely by ICollisionHandler (Phase 2C), not as continuous force.
 * 
 * ⚠️ **NOTE:** This function will be removed in Phase 4 when HSS is fully migrated
 * to ICollisionHandler architecture.
 */
Vec3 HSSCollision(const IonState& ion, const InstrumentDomain& dom) {
    return Vec3{0.0, 0.0, 0.0};
}

} // namespace physics
} // namespace ICARION

