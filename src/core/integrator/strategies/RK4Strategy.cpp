// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   ------------------------------------------------
 *   Modular framework for simulating ion trajectories in custom
 *   electric fields and background gas environments.
 *
 *   @file       RK4Strategy.cpp
 *   @brief      RK4 integration implementation
 *
 *   @date       2025-11-22
 *   @version    1.0.0
 *   @authors    ICARION Development Team
 *
 * =====================================================================
 */

#include "RK4Strategy.h"
#include "core/physics/forces/ForceContext.h"

namespace ICARION {
namespace integrator {

void RK4Strategy::step(
    IonState& ion,
    double t,
    double dt,
    const physics::ForceRegistry& force_registry,
    const std::vector<IonState>& all_ions
) {
    // Create context (SSOT: domain from ForceRegistry)
    physics::ForceContext ctx;
    ctx.domain = force_registry.domain();  // Get domain from registry
    ctx.all_ions = &all_ions;
    ctx.field_provider = nullptr;  // Not using field provider for now
    
    // =========================================================================
    // STAGE 1: k1 = f(t, y)
    // =========================================================================
    Vec3 F1 = force_registry.compute_total_force(ion, t, ctx);
    Vec3 a1 = F1 / ion.mass_kg;  // acceleration = F/m
    
    Vec3 k1_vel = ion.vel;       // k1 for position: dx/dt = v
    Vec3 k1_acc = a1;            // k1 for velocity: dv/dt = a
    
    // =========================================================================
    // STAGE 2: k2 = f(t + dt/2, y + k1*dt/2)
    // =========================================================================
    IonState ion_temp2 = ion;
    ion_temp2.pos = ion.pos + k1_vel * (dt * 0.5);  // x + k1_vel * dt/2
    ion_temp2.vel = ion.vel + k1_acc * (dt * 0.5);  // v + k1_acc * dt/2
    
    Vec3 F2 = force_registry.compute_total_force(ion_temp2, t + dt * 0.5, ctx);
    Vec3 a2 = F2 / ion_temp2.mass_kg;
    
    Vec3 k2_vel = ion_temp2.vel;  // k2 for position: dx/dt at midpoint
    Vec3 k2_acc = a2;             // k2 for velocity: dv/dt at midpoint
    
    // =========================================================================
    // STAGE 3: k3 = f(t + dt/2, y + k2*dt/2)
    // =========================================================================
    IonState ion_temp3 = ion;
    ion_temp3.pos = ion.pos + k2_vel * (dt * 0.5);  // x + k2_vel * dt/2
    ion_temp3.vel = ion.vel + k2_acc * (dt * 0.5);  // v + k2_acc * dt/2
    
    Vec3 F3 = force_registry.compute_total_force(ion_temp3, t + dt * 0.5, ctx);
    Vec3 a3 = F3 / ion_temp3.mass_kg;
    
    Vec3 k3_vel = ion_temp3.vel;  // k3 for position: dx/dt at midpoint
    Vec3 k3_acc = a3;             // k3 for velocity: dv/dt at midpoint
    
    // =========================================================================
    // STAGE 4: k4 = f(t + dt, y + k3*dt)
    // =========================================================================
    IonState ion_temp4 = ion;
    ion_temp4.pos = ion.pos + k3_vel * dt;  // x + k3_vel * dt
    ion_temp4.vel = ion.vel + k3_acc * dt;  // v + k3_acc * dt
    
    Vec3 F4 = force_registry.compute_total_force(ion_temp4, t + dt, ctx);
    Vec3 a4 = F4 / ion_temp4.mass_kg;
    
    Vec3 k4_vel = ion_temp4.vel;  // k4 for position: dx/dt at endpoint
    Vec3 k4_acc = a4;             // k4 for velocity: dv/dt at endpoint;
    
    // =========================================================================
    // FINAL UPDATE: y_new = y + (k1 + 2*k2 + 2*k3 + k4) * dt/6
    // =========================================================================
    // Position update: x_new = x + (k1_vel + 2*k2_vel + 2*k3_vel + k4_vel) * dt/6
    ion.pos = ion.pos + (k1_vel + k2_vel * 2.0 + k3_vel * 2.0 + k4_vel) * (dt / 6.0);
    
    // Velocity update: v_new = v + (k1_acc + 2*k2_acc + 2*k3_acc + k4_acc) * dt/6
    ion.vel = ion.vel + (k1_acc + k2_acc * 2.0 + k3_acc * 2.0 + k4_acc) * (dt / 6.0);
    
    // Note: mass and charge remain constant during integration
}

} // namespace integrator
} // namespace ICARION
