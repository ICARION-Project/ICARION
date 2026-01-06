// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "RK4Strategy.h"
#include "core/physics/forces/ForceContext.h"
#include "core/types/IonEnsemble.h"
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ICARION {
namespace integrator {

namespace {
IonState make_state_from_ensemble(const core::IonEnsemble& ensemble, size_t i) {
    IonState s;
    s.pos = ensemble.get_pos(i);
    s.vel = ensemble.get_vel(i);
    s.mass_kg = ensemble.mass_data()[i];
    s.ion_charge_C = ensemble.charge_data()[i];
    s.active = ensemble.active_data()[i] != 0;
    s.born = ensemble.born_data()[i] != 0;
    s.birth_time_s = ensemble.birth_time(i);
    s.death_time_s = ensemble.death_time(i);
    s.t = ensemble.time(i);
    s.CCS_m2 = ensemble.CCS(i);
    s.reduced_mobility_cm2_Vs = ensemble.mobility(i);
    s.current_domain_index = ensemble.domain_index(i);
    s.species_id = ensemble.species_id(i);
    return s;
}

void write_state_to_scratch(const IonState& state,
                            double temperature,
                            double gas_density,
                            double neutral_mass,
                            core::IonEnsemble& scratch) {
    scratch.resize(1);
    scratch.set_pos(0, state.pos);
    scratch.set_vel(0, state.vel);
    scratch.mass_data()[0] = state.mass_kg;
    scratch.charge_data()[0] = state.ion_charge_C;
    scratch.active_data()[0] = state.active ? 1 : 0;
    scratch.born_data()[0] = state.born ? 1 : 0;
    scratch.birth_time_data()[0] = state.birth_time_s;
    scratch.death_time_data()[0] = state.death_time_s;
    scratch.time_data()[0] = state.t;
    scratch.CCS_data()[0] = state.CCS_m2;
    scratch.mobility_data()[0] = state.reduced_mobility_cm2_Vs;
    scratch.update_species(0, state.species_id, state.mass_kg, state.ion_charge_C,
                           state.CCS_m2, state.reduced_mobility_cm2_Vs);
    scratch.domain_index_data()[0] = state.current_domain_index;
    scratch.temperature_data()[0] = temperature;
    scratch.gas_density_data()[0] = gas_density;
    scratch.neutral_mass_data()[0] = neutral_mass;
}
}  // namespace

void RK4Strategy::step(core::IonEnsemble& ensemble,
                       size_t ion_idx,
                       double t,
                       double dt,
                       const physics::ForceRegistry& force_registry) {
    // Get direct array pointers (cache-friendly)
    auto* pos_x = ensemble.pos_x_data();
    auto* pos_y = ensemble.pos_y_data();
    auto* pos_z = ensemble.pos_z_data();
    auto* vel_x = ensemble.vel_x_data();
    auto* vel_y = ensemble.vel_y_data();
    auto* vel_z = ensemble.vel_z_data();
    auto* mass = ensemble.mass_data();
    
    const size_t i = ion_idx;
    const double inv_mass = 1.0 / mass[i];

    IonState state = make_state_from_ensemble(ensemble, i);
    core::IonEnsemble scratch;
    scratch.reserve(1);
    const double temp = ensemble.temperature(i);
    const double gas_density = ensemble.gas_density(i);
    const double neutral_mass = ensemble.neutral_mass(i);
    
    // Create force context (domain from ForceRegistry)
    physics::ForceContext ctx;
    ctx.domain = force_registry.domain();
    ctx.all_ions = nullptr;
    ctx.field_provider = nullptr;  // No field provider hookup in RK4
    ctx.field_model = force_registry.field_model();
    ctx.ion_ensemble = &scratch;
    ctx.ion_index = 0;

    auto eval_acc = [&](const Vec3& pos, const Vec3& vel, double t_eval) {
        state.pos = pos;
        state.vel = vel;
        write_state_to_scratch(state, temp, gas_density, neutral_mass, scratch);
        Vec3 F = force_registry.compute_total_force(scratch, 0, t_eval, ctx);
        return F * inv_mass;
    };
    
    // =========================================================================
    // STAGE 1: k1 = f(t, y)
    // =========================================================================
    Vec3 pos0{pos_x[i], pos_y[i], pos_z[i]};
    Vec3 vel0{vel_x[i], vel_y[i], vel_z[i]};

    Vec3 k1_vel = vel0;
    Vec3 k1_acc = eval_acc(pos0, vel0, t);
    
    // =========================================================================
    // STAGE 2: k2 = f(t + dt/2, y + k1*dt/2)
    // =========================================================================
    Vec3 pos_stage2 = pos0 + k1_vel * (dt * 0.5);
    Vec3 vel_stage2 = vel0 + k1_acc * (dt * 0.5);
    Vec3 k2_vel = vel_stage2;
    Vec3 k2_acc = eval_acc(pos_stage2, vel_stage2, t + dt * 0.5);
    
    // =========================================================================
    // STAGE 3: k3 = f(t + dt/2, y + k2*dt/2)
    // =========================================================================
    Vec3 pos_stage3 = pos0 + k2_vel * (dt * 0.5);
    Vec3 vel_stage3 = vel0 + k2_acc * (dt * 0.5);
    Vec3 k3_vel = vel_stage3;
    Vec3 k3_acc = eval_acc(pos_stage3, vel_stage3, t + dt * 0.5);
    
    // =========================================================================
    // STAGE 4: k4 = f(t + dt, y + k3*dt)
    // =========================================================================
    Vec3 pos_stage4 = pos0 + k3_vel * dt;
    Vec3 vel_stage4 = vel0 + k3_acc * dt;
    Vec3 k4_vel = vel_stage4;
    Vec3 k4_acc = eval_acc(pos_stage4, vel_stage4, t + dt);
    
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

bool RK4Strategy::step_batch(
    core::IonEnsemble& ensemble,
    double t,
    double dt,
    const std::vector<std::shared_ptr<physics::ForceRegistry>>& registries,
    const std::vector<int>& domain_indices) {
    const size_t n = ensemble.size();
    if (n == 0 || domain_indices.size() != n) {
        return false;
    }

    for (const auto& reg : registries) {
        if (reg && reg->space_charge_model()) {
            return false;
        }
    }

    auto* pos_x = ensemble.pos_x_data();
    auto* pos_y = ensemble.pos_y_data();
    auto* pos_z = ensemble.pos_z_data();
    auto* vel_x = ensemble.vel_x_data();
    auto* vel_y = ensemble.vel_y_data();
    auto* vel_z = ensemble.vel_z_data();
    auto* mass = ensemble.mass_data();
    auto* active = ensemble.active_data();

    std::vector<double> base_px(pos_x, pos_x + n);
    std::vector<double> base_py(pos_y, pos_y + n);
    std::vector<double> base_pz(pos_z, pos_z + n);
    std::vector<double> base_vx(vel_x, vel_x + n);
    std::vector<double> base_vy(vel_y, vel_y + n);
    std::vector<double> base_vz(vel_z, vel_z + n);

    std::vector<Vec3> k1_v(n), k2_v(n), k3_v(n), k4_v(n);
    std::vector<Vec3> k1_a(n), k2_a(n), k3_a(n), k4_a(n);

    constexpr int kOmpChunk = 128;
    const bool use_omp = parallel_enabled_;
#ifndef _OPENMP
    (void)use_omp;
#endif

    auto compute_accels = [&](double t_stage, std::vector<Vec3>& acc_out) {
        #pragma omp parallel if(use_omp)
        {
            #pragma omp for schedule(guided, kOmpChunk)
            for (int i = 0; i < static_cast<int>(n); ++i) {
                if (domain_indices[static_cast<size_t>(i)] < 0 || !active[i]) {
                    continue;
                }
                const int dom = domain_indices[static_cast<size_t>(i)];
                if (dom < 0 || dom >= static_cast<int>(registries.size())) {
                    continue;
                }
                const auto& reg = registries[static_cast<size_t>(dom)];
                if (!reg) {
                    continue;
                }
                physics::ForceContext ctx;
                ctx.domain = reg->domain();
                ctx.field_model = reg->field_model();
                ctx.ion_ensemble = &ensemble;
                ctx.ion_index = static_cast<size_t>(i);
                Vec3 F = reg->compute_total_force(ensemble, static_cast<size_t>(i), t_stage, ctx);
                acc_out[static_cast<size_t>(i)] = F * (1.0 / mass[i]);
            }
        }
    };

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk)
        for (int i = 0; i < static_cast<int>(n); ++i) {
            if (domain_indices[static_cast<size_t>(i)] < 0 || !active[i]) {
                continue;
            }
            k1_v[static_cast<size_t>(i)] = Vec3{base_vx[i], base_vy[i], base_vz[i]};
        }
    }
    compute_accels(t, k1_a);

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk)
        for (int i = 0; i < static_cast<int>(n); ++i) {
            if (domain_indices[static_cast<size_t>(i)] < 0 || !active[i]) {
                continue;
            }
            pos_x[i] = base_px[i] + k1_v[static_cast<size_t>(i)].x * (dt * 0.5);
            pos_y[i] = base_py[i] + k1_v[static_cast<size_t>(i)].y * (dt * 0.5);
            pos_z[i] = base_pz[i] + k1_v[static_cast<size_t>(i)].z * (dt * 0.5);
            vel_x[i] = base_vx[i] + k1_a[static_cast<size_t>(i)].x * (dt * 0.5);
            vel_y[i] = base_vy[i] + k1_a[static_cast<size_t>(i)].y * (dt * 0.5);
            vel_z[i] = base_vz[i] + k1_a[static_cast<size_t>(i)].z * (dt * 0.5);
            k2_v[static_cast<size_t>(i)] = Vec3{vel_x[i], vel_y[i], vel_z[i]};
        }
    }
    compute_accels(t + dt * 0.5, k2_a);

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk)
        for (int i = 0; i < static_cast<int>(n); ++i) {
            if (domain_indices[static_cast<size_t>(i)] < 0 || !active[i]) {
                continue;
            }
            pos_x[i] = base_px[i] + k2_v[static_cast<size_t>(i)].x * (dt * 0.5);
            pos_y[i] = base_py[i] + k2_v[static_cast<size_t>(i)].y * (dt * 0.5);
            pos_z[i] = base_pz[i] + k2_v[static_cast<size_t>(i)].z * (dt * 0.5);
            vel_x[i] = base_vx[i] + k2_a[static_cast<size_t>(i)].x * (dt * 0.5);
            vel_y[i] = base_vy[i] + k2_a[static_cast<size_t>(i)].y * (dt * 0.5);
            vel_z[i] = base_vz[i] + k2_a[static_cast<size_t>(i)].z * (dt * 0.5);
            k3_v[static_cast<size_t>(i)] = Vec3{vel_x[i], vel_y[i], vel_z[i]};
        }
    }
    compute_accels(t + dt * 0.5, k3_a);

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk)
        for (int i = 0; i < static_cast<int>(n); ++i) {
            if (domain_indices[static_cast<size_t>(i)] < 0 || !active[i]) {
                continue;
            }
            pos_x[i] = base_px[i] + k3_v[static_cast<size_t>(i)].x * dt;
            pos_y[i] = base_py[i] + k3_v[static_cast<size_t>(i)].y * dt;
            pos_z[i] = base_pz[i] + k3_v[static_cast<size_t>(i)].z * dt;
            vel_x[i] = base_vx[i] + k3_a[static_cast<size_t>(i)].x * dt;
            vel_y[i] = base_vy[i] + k3_a[static_cast<size_t>(i)].y * dt;
            vel_z[i] = base_vz[i] + k3_a[static_cast<size_t>(i)].z * dt;
            k4_v[static_cast<size_t>(i)] = Vec3{vel_x[i], vel_y[i], vel_z[i]};
        }
    }
    compute_accels(t + dt, k4_a);

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk)
        for (int i = 0; i < static_cast<int>(n); ++i) {
            if (domain_indices[static_cast<size_t>(i)] < 0 || !active[i]) {
                continue;
            }
            Vec3 pos_new = Vec3{base_px[i], base_py[i], base_pz[i]} +
                (k1_v[static_cast<size_t>(i)] + k2_v[static_cast<size_t>(i)] * 2.0 +
                 k3_v[static_cast<size_t>(i)] * 2.0 + k4_v[static_cast<size_t>(i)]) * (dt / 6.0);
            Vec3 vel_new = Vec3{base_vx[i], base_vy[i], base_vz[i]} +
                (k1_a[static_cast<size_t>(i)] + k2_a[static_cast<size_t>(i)] * 2.0 +
                 k3_a[static_cast<size_t>(i)] * 2.0 + k4_a[static_cast<size_t>(i)]) * (dt / 6.0);

            pos_x[i] = pos_new.x;
            pos_y[i] = pos_new.y;
            pos_z[i] = pos_new.z;
            vel_x[i] = vel_new.x;
            vel_y[i] = vel_new.y;
            vel_z[i] = vel_new.z;
        }
    }

    return true;
}

} // namespace integrator
} // namespace ICARION
