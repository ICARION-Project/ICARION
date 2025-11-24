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
    
    Vec3 vel1 = ion.vel;         // dx/dt at t
    Vec3 acc1 = a1;              // dv/dt at t
    
    // =========================================================================
    // STAGE 2: k2 = f(t + dt/2, y + k1*dt/2)
    // =========================================================================
    IonState ion_temp2 = ion;
    ion_temp2.pos += vel1 * (dt * 0.5);  // x + dx/dt * dt/2
    ion_temp2.vel += acc1 * (dt * 0.5);  // v + dv/dt * dt/2
    
    Vec3 F2 = force_registry.compute_total_force(ion_temp2, t + dt * 0.5, ctx);
    Vec3 a2 = F2 / ion_temp2.mass_kg;
    
    Vec3 vel2 = ion_temp2.vel;
    Vec3 acc2 = a2;
    
    // =========================================================================
    // STAGE 3: k3 = f(t + dt/2, y + k2*dt/2)
    // =========================================================================
    IonState ion_temp3 = ion;
    ion_temp3.pos += vel2 * (dt * 0.5);  // x + dx/dt(midpoint) * dt/2
    ion_temp3.vel += acc2 * (dt * 0.5);  // v + dv/dt(midpoint) * dt/2
    
    Vec3 F3 = force_registry.compute_total_force(ion_temp3, t + dt * 0.5, ctx);
    Vec3 a3 = F3 / ion_temp3.mass_kg;
    
    Vec3 vel3 = ion_temp3.vel;
    Vec3 acc3 = a3;
    
    // =========================================================================
    // STAGE 4: k4 = f(t + dt, y + k3*dt)
    // =========================================================================
    IonState ion_temp4 = ion;
    ion_temp4.pos += vel3 * dt;  // x + dx/dt(midpoint) * dt
    ion_temp4.vel += acc3 * dt;  // v + dv/dt(midpoint) * dt
    
    Vec3 F4 = force_registry.compute_total_force(ion_temp4, t + dt, ctx);
    Vec3 a4 = F4 / ion_temp4.mass_kg;
    
    Vec3 vel4 = ion_temp4.vel;
    Vec3 acc4 = a4;
    
    // =========================================================================
    // FINAL UPDATE: y_new = y + (k1 + 2*k2 + 2*k3 + k4) * dt/6
    // =========================================================================
    // Position update: x_new = x + (v1 + 2*v2 + 2*v3 + v4) * dt/6
    ion.pos = ion.pos + (vel1 + vel2 * 2.0 + vel3 * 2.0 + vel4) * (dt / 6.0);
    
    // Velocity update: v_new = v + (a1 + 2*a2 + 2*a3 + a4) * dt/6
    ion.vel = ion.vel + (acc1 + acc2 * 2.0 + acc3 * 2.0 + acc4) * (dt / 6.0);
    
    // Note: mass and charge remain constant during integration
}

} // namespace integrator
} // namespace ICARION
