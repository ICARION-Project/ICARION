// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "RK4Strategy.h"
#include "core/physics/forces/ForceContext.h"

namespace ICARION {
namespace integrator {

void RK4Strategy::step(
    core::IonEnsemble& ensemble,
    size_t ion_idx,
    double t,
    double dt,
    const physics::ForceRegistry& force_registry
) {
    // Get direct array pointers (cache-friendly)
    auto* pos_x = ensemble.pos_x_data();
    auto* pos_y = ensemble.pos_y_data();
    auto* pos_z = ensemble.pos_z_data();
    auto* vel_x = ensemble.vel_x_data();
    auto* vel_y = ensemble.vel_y_data();
    auto* vel_z = ensemble.vel_z_data();
    auto* mass = ensemble.mass_data();
    auto* charge = ensemble.charge_data();
    
    const size_t i = ion_idx;
    const double inv_mass = 1.0 / mass[i];
    
    // Create force context (domain from ForceRegistry)
    physics::ForceContext ctx;
    ctx.domain = force_registry.domain();
    ctx.all_ions = nullptr;  // AoS-only fallback
    ctx.field_provider = nullptr;  // No field provider hookup in RK4
    ctx.field_model = force_registry.field_model();
    ctx.ion_ensemble = &ensemble;
    ctx.ion_index = i;
    
    // Temporary IonState for force evaluation (minimal overhead)
    IonState ion_temp;
    ion_temp.mass_kg = mass[i];
    ion_temp.ion_charge_C = charge[i];
    
    // =========================================================================
    // STAGE 1: k1 = f(t, y)
    // =========================================================================
    ion_temp.pos = Vec3(pos_x[i], pos_y[i], pos_z[i]);
    ion_temp.vel = Vec3(vel_x[i], vel_y[i], vel_z[i]);
    
    Vec3 F1 = force_registry.compute_total_force(ion_temp, t, ctx);
    Vec3 a1 = F1 * inv_mass;
    
    Vec3 k1_vel = ion_temp.vel;
    Vec3 k1_acc = a1;
    
    // =========================================================================
    // STAGE 2: k2 = f(t + dt/2, y + k1*dt/2)
    // =========================================================================
    ion_temp.pos = Vec3(pos_x[i], pos_y[i], pos_z[i]) + k1_vel * (dt * 0.5);
    ion_temp.vel = Vec3(vel_x[i], vel_y[i], vel_z[i]) + k1_acc * (dt * 0.5);
    
    Vec3 F2 = force_registry.compute_total_force(ion_temp, t + dt * 0.5, ctx);
    Vec3 a2 = F2 * inv_mass;
    
    Vec3 k2_vel = ion_temp.vel;
    Vec3 k2_acc = a2;
    
    // =========================================================================
    // STAGE 3: k3 = f(t + dt/2, y + k2*dt/2)
    // =========================================================================
    ion_temp.pos = Vec3(pos_x[i], pos_y[i], pos_z[i]) + k2_vel * (dt * 0.5);
    ion_temp.vel = Vec3(vel_x[i], vel_y[i], vel_z[i]) + k2_acc * (dt * 0.5);
    
    Vec3 F3 = force_registry.compute_total_force(ion_temp, t + dt * 0.5, ctx);
    Vec3 a3 = F3 * inv_mass;
    
    Vec3 k3_vel = ion_temp.vel;
    Vec3 k3_acc = a3;
    
    // =========================================================================
    // STAGE 4: k4 = f(t + dt, y + k3*dt)
    // =========================================================================
    ion_temp.pos = Vec3(pos_x[i], pos_y[i], pos_z[i]) + k3_vel * dt;
    ion_temp.vel = Vec3(vel_x[i], vel_y[i], vel_z[i]) + k3_acc * dt;
    
    Vec3 F4 = force_registry.compute_total_force(ion_temp, t + dt, ctx);
    Vec3 a4 = F4 * inv_mass;
    
    Vec3 k4_vel = ion_temp.vel;
    Vec3 k4_acc = a4;
    
    // =========================================================================
    // FINAL UPDATE: Direct array write (cache-friendly!)
    // =========================================================================
    Vec3 pos_new = Vec3(pos_x[i], pos_y[i], pos_z[i]) + 
                   (k1_vel + k2_vel * 2.0 + k3_vel * 2.0 + k4_vel) * (dt / 6.0);
    Vec3 vel_new = Vec3(vel_x[i], vel_y[i], vel_z[i]) + 
                   (k1_acc + k2_acc * 2.0 + k3_acc * 2.0 + k4_acc) * (dt / 6.0);
    
    // Write back to SoA arrays
    pos_x[i] = pos_new.x;
    pos_y[i] = pos_new.y;
    pos_z[i] = pos_new.z;
    vel_x[i] = vel_new.x;
    vel_y[i] = vel_new.y;
    vel_z[i] = vel_new.z;
}

} // namespace integrator
} // namespace ICARION
