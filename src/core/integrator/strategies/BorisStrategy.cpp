// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "BorisStrategy.h"
#include "core/physics/forces/ForceContext.h"
#include <cmath>

namespace ICARION {
namespace integrator {

namespace {
    // Boris algorithm constants
    constexpr double HALF_TIMESTEP = 0.5;  // For half-step electric acceleration
    constexpr double S_NUMERATOR = 2.0;    // Numerator in s = 2t/(1 + t²) formula
}

void BorisStrategy::step(
    IonState& ion,
    double t,
    double dt,
    const physics::ForceRegistry& force_registry,
    const std::vector<IonState>& all_ions
) {
    // Get domain from ForceRegistry (SSOT!) for magnetic field access
    const config::DomainConfig* domain = force_registry.domain();
    if (!domain) {
        throw std::runtime_error("BorisStrategy: ForceRegistry has no domain configured");
    }
    
    // Create context
    physics::ForceContext ctx;
    ctx.domain = domain;
    ctx.all_ions = &all_ions;
    ctx.field_provider = nullptr;
    ctx.field_model = force_registry.field_model();
    
    // =========================================================================
    // STEP 1: Compute electric forces and get magnetic field
    // =========================================================================
    // Note: Boris method treats E and B separately.
    // The magnetic force is NOT computed explicitly - it's handled via rotation.
    // We compute only electric-like forces (E-field, space charge, drag, etc.)
    
    Vec3 F_electric = force_registry.compute_total_force(ion, t, ctx);
    Vec3 a_electric = F_electric / ion.mass_kg;
    
    // Get magnetic field at current position (uniform field only; gradients ignored)
    Vec3 B{0, 0, 0};
    if (domain->fields.magnetic.enabled) {
        // Uniform magnetic field (most common case)
        B = domain->fields.magnetic.field_strength_T;
        // Note: field_gradient is ignored for now (uniform field approximation)
    }
    
    // =========================================================================
    // STEP 2: Half-step electric acceleration (v^n → v^-)
    // =========================================================================
    Vec3 v_minus = ion.vel + a_electric * (dt * 0.5);
    
    // =========================================================================
    // STEP 3: Magnetic rotation (v^- → v^+)
    // =========================================================================
    // Boris rotation parameters
    // t = (q/m) * B * (dt/2)
    double qm = ion.ion_charge_C / ion.mass_kg;
    Vec3 t_vec = B * (qm * dt * HALF_TIMESTEP);
    
    // s = 2*t / (1 + |t|²)
    double t_mag_sq = t_vec.x * t_vec.x + t_vec.y * t_vec.y + t_vec.z * t_vec.z;
    double s_factor = S_NUMERATOR / (1.0 + t_mag_sq);
    Vec3 s_vec = t_vec * s_factor;
    
    // v' = v^- + v^- × t
    Vec3 v_minus_cross_t{
        v_minus.y * t_vec.z - v_minus.z * t_vec.y,
        v_minus.z * t_vec.x - v_minus.x * t_vec.z,
        v_minus.x * t_vec.y - v_minus.y * t_vec.x
    };
    Vec3 v_prime = v_minus + v_minus_cross_t;
    
    // v^+ = v^- + v' × s
    Vec3 v_prime_cross_s{
        v_prime.y * s_vec.z - v_prime.z * s_vec.y,
        v_prime.z * s_vec.x - v_prime.x * s_vec.z,
        v_prime.x * s_vec.y - v_prime.y * s_vec.x
    };
    Vec3 v_plus = v_minus + v_prime_cross_s;
    
    // =========================================================================
    // STEP 4: Position update (x^n → x^(n+1))
    // =========================================================================
    // Use v^+ (velocity at half-step) for leapfrog position update
    Vec3 pos_new = ion.pos + v_plus * dt;
    
    // =========================================================================
    // STEP 5: Half-step electric acceleration (v^+ → v^(n+1))
    // =========================================================================
    Vec3 v_new = v_plus + a_electric * (dt * HALF_TIMESTEP);
    
    // Update ion state
    ion.vel = v_new;
    ion.pos = pos_new;
    
    // Note: mass and charge remain constant during integration
}

void BorisStrategy::step_soa(
    core::IonEnsemble& ensemble,
    size_t ion_idx,
    double t,
    double dt,
    const physics::ForceRegistry& force_registry
) {
    const config::DomainConfig* domain = force_registry.domain();
    if (!domain) {
        throw std::runtime_error("BorisStrategy: ForceRegistry has no domain configured");
    }

    // Direct SoA pointers
    auto* pos_x = ensemble.pos_x_data();
    auto* pos_y = ensemble.pos_y_data();
    auto* pos_z = ensemble.pos_z_data();
    auto* vel_x = ensemble.vel_x_data();
    auto* vel_y = ensemble.vel_y_data();
    auto* vel_z = ensemble.vel_z_data();
    auto* mass = ensemble.mass_data();
    auto* charge = ensemble.charge_data();

    // Prepare minimal IonState for force evaluation (applies_to etc.)
    IonState ion;
    ion.pos = {pos_x[ion_idx], pos_y[ion_idx], pos_z[ion_idx]};
    ion.vel = {vel_x[ion_idx], vel_y[ion_idx], vel_z[ion_idx]};
    ion.mass_kg = mass[ion_idx];
    ion.ion_charge_C = charge[ion_idx];
    ion.active = ensemble.active_data()[ion_idx] != 0;
    ion.born = ensemble.born_data()[ion_idx] != 0;
    ion.current_domain_index = ensemble.domain_index(ion_idx);
    ion.CCS_m2 = ensemble.CCS(ion_idx);
    ion.reduced_mobility_cm2_Vs = ensemble.mobility(ion_idx);
    ion.species_id = ensemble.species_id(ion_idx);
    ion.birth_time_s = ensemble.birth_time(ion_idx);

    const double qm = charge[ion_idx] / mass[ion_idx];

    physics::ForceContext ctx;
    ctx.domain = domain;
    ctx.all_ions = nullptr;  // SoA path uses ion_ensemble
    ctx.field_provider = nullptr;
    ctx.field_model = force_registry.field_model();
    ctx.ion_ensemble = &ensemble;
    ctx.ion_index = ion_idx;

    Vec3 F_electric = force_registry.compute_total_force_soa(ensemble, ion_idx, t, ctx);
    Vec3 a_electric = F_electric / ion.mass_kg;

    Vec3 B{0, 0, 0};
    if (domain->fields.magnetic.enabled) {
        B = domain->fields.magnetic.field_strength_T;
    }

    Vec3 v_minus = ion.vel + a_electric * (dt * 0.5);

    Vec3 t_vec = B * (qm * dt * HALF_TIMESTEP);
    double t_mag_sq = t_vec.x * t_vec.x + t_vec.y * t_vec.y + t_vec.z * t_vec.z;
    double s_factor = S_NUMERATOR / (1.0 + t_mag_sq);
    Vec3 s_vec = t_vec * s_factor;

    Vec3 v_minus_cross_t{
        v_minus.y * t_vec.z - v_minus.z * t_vec.y,
        v_minus.z * t_vec.x - v_minus.x * t_vec.z,
        v_minus.x * t_vec.y - v_minus.y * t_vec.x
    };
    Vec3 v_prime = v_minus + v_minus_cross_t;

    Vec3 v_prime_cross_s{
        v_prime.y * s_vec.z - v_prime.z * s_vec.y,
        v_prime.z * s_vec.x - v_prime.x * s_vec.z,
        v_prime.x * s_vec.y - v_prime.y * s_vec.x
    };
    Vec3 v_plus = v_minus + v_prime_cross_s;

    Vec3 pos_new = ion.pos + v_plus * dt;
    Vec3 v_new = v_plus + a_electric * (dt * 0.5);

    pos_x[ion_idx] = pos_new.x;
    pos_y[ion_idx] = pos_new.y;
    pos_z[ion_idx] = pos_new.z;
    vel_x[ion_idx] = v_new.x;
    vel_y[ion_idx] = v_new.y;
    vel_z[ion_idx] = v_new.z;
}

} // namespace integrator
} // namespace ICARION
