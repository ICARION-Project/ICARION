// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include "RK45Strategy.h"
#include "core/config/types/IFieldModel.h"
#include "core/physics/forces/ForceContext.h"
#include "core/physics/forces/IForce.h"
#include "core/types/IonEnsemble.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <optional>
#include <spdlog/spdlog.h>
#include <stdexcept>

namespace ICARION {
namespace integrator {

namespace {

bool rk45_force_debug_enabled() {
  static const bool enabled = [] {
    if (const char *env = std::getenv("ICARION_DEBUG_RK45_FORCE")) {
      switch (env[0]) {
      case '1':
      case 't':
      case 'T':
      case 'y':
      case 'Y':
        return true;
      default:
        break;
      }
    }
    return false;
  }();
  return enabled;
}

bool rk45_state_debug_enabled() {
  static const bool enabled = [] {
    if (const char *env = std::getenv("ICARION_DEBUG_RK45_STATE")) {
      switch (env[0]) {
      case '1':
      case 't':
      case 'T':
      case 'y':
      case 'Y':
        return true;
      default:
        break;
      }
    }
    return false;
  }();
  return enabled;
}

IonState make_state_from_ensemble(const core::IonEnsemble &ensemble, size_t i) {
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

physics::ForceState make_force_state(const IonState &state,
                                     std::optional<size_t> ion_index) {
  physics::ForceState fs;
  fs.pos = state.pos;
  fs.vel = state.vel;
  fs.mass_kg = state.mass_kg;
  fs.ion_charge_C = state.ion_charge_C;
  fs.CCS_m2 = state.CCS_m2;
  fs.reduced_mobility_cm2_Vs = state.reduced_mobility_cm2_Vs;
  fs.species_id = state.species_id;
  fs.active = state.active;
  fs.born = state.born;
  fs.current_domain_index = state.current_domain_index;
  fs.birth_time_s = state.birth_time_s;
  fs.ensemble_index = ion_index;
  return fs;
}

// Dormand-Prince 5(4) Butcher tableau coefficients
// Reference: Hairer, Norsett, Wanner (1993) "Solving Ordinary Differential
// Equations I"

// Time nodes for stages
constexpr double c2 = 1.0 / 5.0;
constexpr double c3 = 3.0 / 10.0;
constexpr double c4 = 4.0 / 5.0;
constexpr double c5 = 8.0 / 9.0;
constexpr double c6 = 1.0;
constexpr double c7 = 1.0;

// Runge-Kutta matrix (a_ij coefficients)
// Stage 2
constexpr double a21 = 1.0 / 5.0;

// Stage 3
constexpr double a31 = 3.0 / 40.0;
constexpr double a32 = 9.0 / 40.0;

// Stage 4
constexpr double a41 = 44.0 / 45.0;
constexpr double a42 = -56.0 / 15.0;
constexpr double a43 = 32.0 / 9.0;

// Stage 5
constexpr double a51 = 19372.0 / 6561.0;
constexpr double a52 = -25360.0 / 2187.0;
constexpr double a53 = 64448.0 / 6561.0;
constexpr double a54 = -212.0 / 729.0;

// Stage 6
constexpr double a61 = 9017.0 / 3168.0;
constexpr double a62 = -355.0 / 33.0;
constexpr double a63 = 46732.0 / 5247.0;
constexpr double a64 = 49.0 / 176.0;
constexpr double a65 = -5103.0 / 18656.0;

// Stage 7 (FSAL: used for 4th order solution and next step's k1)
constexpr double a71 = 35.0 / 384.0;
constexpr double a72 = 0.0;
constexpr double a73 = 500.0 / 1113.0;
constexpr double a74 = 125.0 / 192.0;
constexpr double a75 = -2187.0 / 6784.0;
constexpr double a76 = 11.0 / 84.0;

// 4th order solution weights (b_i coefficients) - SAME as stage 7
constexpr double b41 = 35.0 / 384.0;
constexpr double b42 = 0.0;
constexpr double b43 = 500.0 / 1113.0;
constexpr double b44 = 125.0 / 192.0;
constexpr double b45 = -2187.0 / 6784.0;
constexpr double b46 = 11.0 / 84.0;
constexpr double b47 = 0.0; // FSAL property

// 5th order solution weights (b*_i coefficients) - for error estimation
constexpr double b51 = 5179.0 / 57600.0;
constexpr double b52 = 0.0;
constexpr double b53 = 7571.0 / 16695.0;
constexpr double b54 = 393.0 / 640.0;
constexpr double b55 = -92097.0 / 339200.0;
constexpr double b56 = 187.0 / 2100.0;
constexpr double b57 = 1.0 / 40.0;

// PI controller coefficients (optimized for order 5)
constexpr double PI_BETA = 0.04;                  // Stabilization parameter
constexpr double PI_ALPHA = 0.2 - PI_BETA * 0.75; // Proportional gain
constexpr double REJECTION_EXPONENT =
    0.2; // 1/order for rejected steps (1/5 for order 5)

// Step control parameters
constexpr int MAX_REJECT_ATTEMPTS = 10;
constexpr double SAFETY_MARGIN = 0.9;
constexpr double MIN_ERROR_THRESHOLD =
    1e-10; // Below this, error estimate unreliable
constexpr double ERROR_ACCEPTANCE_THRESHOLD =
    1.0; // Error <= 1.0 means step accepted
constexpr double DT_MIN_TOLERANCE =
    1.001; // Accept step if dt near dt_min (1.001x factor)

// Spatial dimensions
constexpr int NUM_SPATIAL_DIMS = 3; // x, y, z components
constexpr int X_INDEX = 0;
constexpr int Y_INDEX = 1;
constexpr int Z_INDEX = 2;
} // namespace

RK45Strategy::RK45Strategy() : config_{}, stats_{} {
  // Use default configuration
}

RK45Strategy::RK45Strategy(const AdaptiveConfig &config)
    : config_(config), stats_{} {
  // Validate configuration
  if (config_.atol <= 0.0 || config_.rtol <= 0.0) {
    throw std::invalid_argument("RK45Strategy: tolerances must be positive");
  }
  if (config_.safety_factor <= 0.0 || config_.safety_factor >= 1.0) {
    throw std::invalid_argument(
        "RK45Strategy: safety factor must be in (0, 1)");
  }
  if (config_.absolute_min_step_s < 0.0) {
    throw std::invalid_argument(
        "RK45Strategy: absolute_min_step_s must be non-negative");
  }
}

RK45Strategy::RK45State &RK45Strategy::state_for(size_t ion_idx) {
  if (per_ion_state_.size() <= ion_idx) {
    per_ion_state_.resize(ion_idx + 1);
  }
  return per_ion_state_[ion_idx];
}

void RK45Strategy::compute_acceleration_batch(
    double &ax, double &ay, double &az, const core::IonEnsemble &ensemble,
    size_t ion_idx, double t, const physics::ForceRegistry &force_registry) {
  physics::ForceContext ctx;
  ctx.domain = force_registry.domain();
  ctx.all_ions = nullptr;
  ctx.field_provider = nullptr;
  ctx.field_model = force_registry.field_model();
  ctx.ion_ensemble = &ensemble;
  ctx.ion_index = ion_idx;

  Vec3 F = force_registry.compute_total_force(ensemble, ion_idx, t, ctx);
  const double inv_mass = 1.0 / ensemble.mass_data()[ion_idx];
  Vec3 a = F * inv_mass;

  ax = a.x;
  ay = a.y;
  az = a.z;
}

void RK45Strategy::compute_acceleration_state(
    double &ax, double &ay, double &az, const IonState &state, double t,
    const physics::ForceRegistry &force_registry,
    std::optional<size_t> ion_index) {
  physics::ForceContext ctx;
  ctx.domain = force_registry.domain();
  ctx.all_ions = nullptr;
  ctx.field_provider = nullptr;
  ctx.field_model = force_registry.field_model();
  ctx.ion_ensemble = nullptr;
  ctx.ion_index = ion_index.value_or(static_cast<size_t>(-1));

  physics::ForceState force_state = make_force_state(state, ion_index);
  Vec3 F = force_registry.compute_total_force_soa(force_state, t, ctx);
  const double inv_mass = 1.0 / state.mass_kg;
  Vec3 a = F * inv_mass;

  static thread_local int rk45_force_debug_counter = 0;
  if (rk45_force_debug_enabled() && rk45_force_debug_counter < 20 &&
      state.mass_kg > 0.0) {
    Vec3 E_ref{0.0, 0.0, 0.0};
    if (ctx.field_model) {
      E_ref = ctx.field_model->E(state.pos, t);
    }
    const double q_over_m = state.ion_charge_C / state.mass_kg;
    Vec3 a_from_field = E_ref * q_over_m;
    Vec3 delta = a - a_from_field;
    SPDLOG_INFO(
        "[RK45Force] species={}, domain={}, t={:.6e}, pos=({:.6e}, {:.6e}, "
        "{:.6e}), F=({:.6e}, {:.6e}, {:.6e}), a=({:.6e}, {:.6e}, {:.6e}), "
        "a_E=({:.6e}, {:.6e}, {:.6e}), delta=({:.6e}, {:.6e}, {:.6e})",
        state.species_id, state.current_domain_index, t, state.pos.x,
        state.pos.y, state.pos.z, F.x, F.y, F.z, a.x, a.y, a.z, a_from_field.x,
        a_from_field.y, a_from_field.z, delta.x, delta.y, delta.z);
    ++rk45_force_debug_counter;
  }

  ax = a.x;
  ay = a.y;
  az = a.z;
}

double RK45Strategy::estimate_error(const IonState &y4, const IonState &y5,
                                    const IonState &y_current) const {
  // Compute scaled error for position and velocity components
  // error = max_i |e_i| / (atol + rtol * |y_i|)

  double max_error = 0.0;

  // Position error (x, y, z)
  for (int i = 0; i < NUM_SPATIAL_DIMS; ++i) {
    double y5_val = (i == X_INDEX)   ? y5.pos.x
                    : (i == Y_INDEX) ? y5.pos.y
                                     : y5.pos.z;
    double y4_val = (i == X_INDEX)   ? y4.pos.x
                    : (i == Y_INDEX) ? y4.pos.y
                                     : y4.pos.z;
    double y_val = (i == X_INDEX)   ? y_current.pos.x
                   : (i == Y_INDEX) ? y_current.pos.y
                                    : y_current.pos.z;

    double err_abs = std::fabs(y5_val - y4_val);
    double scale = config_.atol + config_.rtol * std::fabs(y_val);
    double err_scaled = err_abs / scale;

    max_error = std::max(max_error, err_scaled);
  }

  // Velocity error (vx, vy, vz)
  for (int i = 0; i < NUM_SPATIAL_DIMS; ++i) {
    double y5_val = (i == X_INDEX)   ? y5.vel.x
                    : (i == Y_INDEX) ? y5.vel.y
                                     : y5.vel.z;
    double y4_val = (i == X_INDEX)   ? y4.vel.x
                    : (i == Y_INDEX) ? y4.vel.y
                                     : y4.vel.z;
    double y_val = (i == X_INDEX)   ? y_current.vel.x
                   : (i == Y_INDEX) ? y_current.vel.y
                                    : y_current.vel.z;

    double err_abs = std::fabs(y5_val - y4_val);
    double scale = config_.atol + config_.rtol * std::fabs(y_val);
    double err_scaled = err_abs / scale;

    max_error = std::max(max_error, err_scaled);
  }

  // Avoid division by zero in step control
  return std::max(max_error, MIN_ERROR_THRESHOLD);
}

double RK45Strategy::compute_new_step(double current_dt, double error,
                                      double dt_min, double dt_max,
                                      double prev_error) {
  // PI controller for step size adjustment
  // Gustafsson, K. (1991): "Control theoretic techniques for stepsize
  // selection" dt_new = dt * (error_prev / error)^(beta/order) * (1 /
  // error)^(alpha/order)

  double factor;

  if (error > 1.0) {
    // Step rejected: use conservative factor
    factor = config_.safety_factor * std::pow(1.0 / error, REJECTION_EXPONENT);
  } else {
    // Step accepted: use PI controller
    factor = config_.safety_factor * std::pow(prev_error / error, PI_BETA) *
             std::pow(1.0 / error, PI_ALPHA);
  }

  // Clamp growth/shrinkage per step
  factor = std::max(factor, config_.max_step_decrease);
  factor = std::min(factor, config_.max_step_increase);

  double dt_new = current_dt * factor;

  // Clamp to absolute limits
  dt_new = std::max(dt_new, dt_min);
  dt_new = std::min(dt_new, dt_max);

  return dt_new;
}

void RK45Strategy::step(core::IonEnsemble &ensemble, size_t ion_idx, double t,
                        double dt,
                        const physics::ForceRegistry &force_registry) {
  double dt_variable = dt;
  double dt_used = dt;

  IonState ion = make_state_from_ensemble(ensemble, ion_idx);
  auto &st = state_for(ion_idx);

  // Run adaptive step using SoA acceleration path; ignore dt update (fixed-step
  // interface)
  const double dt_initial = dt_variable;
  const config::DomainConfig *domain = force_registry.domain();
  if (!domain) {
    throw std::runtime_error(
        "RK45Strategy: ForceRegistry has no domain configured");
  }
  const double dt_min = std::max(dt_initial * config_.min_step_factor,
                                 config_.absolute_min_step_s);
  const double dt_max = dt_initial * config_.max_step_factor;

  double dt_work = dt_variable;
  bool step_accepted = false;
  int attempts = 0;
  IonState y0 = ion;

  while (!step_accepted && attempts < MAX_REJECT_ATTEMPTS) {
    Vec3 k1_v, k1_a;
    Vec3 k2_v, k2_a;
    Vec3 k3_v, k3_a;
    Vec3 k4_v, k4_a;
    Vec3 k5_v, k5_a;
    Vec3 k6_v, k6_a;
    Vec3 k7_v, k7_a;

    k1_v = y0.vel;
    if (st.fsal_available) {
      k1_a = Vec3{st.k1_ax, st.k1_ay, st.k1_az};
    } else {
      double ax, ay, az;
      compute_acceleration_state(ax, ay, az, y0, t, force_registry, ion_idx);
      k1_a = Vec3{ax, ay, az};
    }

    IonState y2 = y0;
    y2.pos += k1_v * (dt_work * a21);
    y2.vel += k1_a * (dt_work * a21);
    k2_v = y2.vel;
    double ax2, ay2, az2;
    compute_acceleration_state(ax2, ay2, az2, y2, t + c2 * dt_work,
                               force_registry, ion_idx);
    k2_a = Vec3{ax2, ay2, az2};

    IonState y3 = y0;
    y3.pos += (k1_v * a31 + k2_v * a32) * dt_work;
    y3.vel += (k1_a * a31 + k2_a * a32) * dt_work;
    k3_v = y3.vel;
    double ax3, ay3, az3;
    compute_acceleration_state(ax3, ay3, az3, y3, t + c3 * dt_work,
                               force_registry, ion_idx);
    k3_a = Vec3{ax3, ay3, az3};

    IonState y4_temp = y0;
    y4_temp.pos += (k1_v * a41 + k2_v * a42 + k3_v * a43) * dt_work;
    y4_temp.vel += (k1_a * a41 + k2_a * a42 + k3_a * a43) * dt_work;
    k4_v = y4_temp.vel;
    double ax4, ay4, az4;
    compute_acceleration_state(ax4, ay4, az4, y4_temp, t + c4 * dt_work,
                               force_registry, ion_idx);
    k4_a = Vec3{ax4, ay4, az4};

    IonState y5_temp = y0;
    y5_temp.pos +=
        (k1_v * a51 + k2_v * a52 + k3_v * a53 + k4_v * a54) * dt_work;
    y5_temp.vel +=
        (k1_a * a51 + k2_a * a52 + k3_a * a53 + k4_a * a54) * dt_work;
    k5_v = y5_temp.vel;
    double ax5, ay5, az5;
    compute_acceleration_state(ax5, ay5, az5, y5_temp, t + c5 * dt_work,
                               force_registry, ion_idx);
    k5_a = Vec3{ax5, ay5, az5};

    IonState y6 = y0;
    y6.pos += (k1_v * a61 + k2_v * a62 + k3_v * a63 + k4_v * a64 + k5_v * a65) *
              dt_work;
    y6.vel += (k1_a * a61 + k2_a * a62 + k3_a * a63 + k4_a * a64 + k5_a * a65) *
              dt_work;
    k6_v = y6.vel;
    double ax6, ay6, az6;
    compute_acceleration_state(ax6, ay6, az6, y6, t + c6 * dt_work,
                               force_registry, ion_idx);
    k6_a = Vec3{ax6, ay6, az6};

    IonState y7 = y0;
    y7.pos += (k1_v * a71 + k2_v * a72 + k3_v * a73 + k4_v * a74 + k5_v * a75 +
               k6_v * a76) *
              dt_work;
    y7.vel += (k1_a * a71 + k2_a * a72 + k3_a * a73 + k4_a * a74 + k5_a * a75 +
               k6_a * a76) *
              dt_work;
    k7_v = y7.vel;
    double ax7, ay7, az7;
    compute_acceleration_state(ax7, ay7, az7, y7, t + c7 * dt_work,
                               force_registry, ion_idx);
    k7_a = Vec3{ax7, ay7, az7};

    IonState y4 = y7; // FSAL property

    IonState y5 = y0;
    y5.pos += (k1_v * b51 + k2_v * b52 + k3_v * b53 + k4_v * b54 + k5_v * b55 +
               k6_v * b56 + k7_v * b57) *
              dt_work;
    y5.vel += (k1_a * b51 + k2_a * b52 + k3_a * b53 + k4_a * b54 + k5_a * b55 +
               k6_a * b56 + k7_a * b57) *
              dt_work;

    double error = estimate_error(y4, y5, y0);

    const bool force_accept = config_.accept_at_dt_min &&
        (dt_work <= dt_min * DT_MIN_TOLERANCE);
    if (error <= ERROR_ACCEPTANCE_THRESHOLD || force_accept) {
      ion = y4;
      step_accepted = true;

      st.k1_ax = k7_a.x;
      st.k1_ay = k7_a.y;
      st.k1_az = k7_a.z;
      st.fsal_available = true;

      if (stats_enabled_) {
        stats_.accepted_steps++;
        stats_.min_step_used = std::min(stats_.min_step_used, dt_work);
        stats_.max_step_used = std::max(stats_.max_step_used, dt_work);
        stats_.avg_error =
            (stats_.avg_error * (stats_.accepted_steps - 1) + error) /
            stats_.accepted_steps;
        stats_.sum_step_used += dt_work;
      }

      double dt_next =
          compute_new_step(dt_work, error, dt_min, dt_max, st.last_error);
      dt_variable = dt_next;
      dt_used = dt_work;
      st.last_error = error;
      last_dt_used_ = dt_used;
      last_dt_suggested_ = dt_variable;
    } else {
      if (stats_enabled_) {
        stats_.rejected_steps++;
      }
      double dt_next =
          compute_new_step(dt_work, error, dt_min, dt_max, st.last_error);
      dt_work = dt_next;
      attempts++;
      st.fsal_available = false;
    }
  }

  if (!step_accepted) {
    throw std::runtime_error("RK45Strategy: Failed to converge after " +
                             std::to_string(MAX_REJECT_ATTEMPTS) + " attempts");
  }

  // Write back to SoA
  const bool state_logging_enabled = rk45_state_debug_enabled();
  static thread_local int rk45_state_debug_counter = 0;
  bool captured_state = false;
  Vec3 pos_before;
  Vec3 vel_before;
  if (state_logging_enabled && rk45_state_debug_counter < 20) {
    pos_before = ensemble.get_pos(ion_idx);
    vel_before = ensemble.get_vel(ion_idx);
    captured_state = true;
  }

  ensemble.set_pos(ion_idx, ion.pos);
  ensemble.set_vel(ion_idx, ion.vel);

  if (captured_state) {
    Vec3 pos_after = ensemble.get_pos(ion_idx);
    Vec3 vel_after = ensemble.get_vel(ion_idx);
    SPDLOG_INFO("[RK45State] ion_idx={} t={:.6e} dt_used={:.6e} "
                "pos_before=({:.6e}, {:.6e}, {:.6e}) vel_before=({:.6e}, "
                "{:.6e}, {:.6e}) solver_pos=({:.6e}, {:.6e}, {:.6e}) "
                "solver_vel=({:.6e}, {:.6e}, {:.6e}) pos_after=({:.6e}, "
                "{:.6e}, {:.6e}) vel_after=({:.6e}, {:.6e}, {:.6e})",
                ion_idx, t, dt_used, pos_before.x, pos_before.y, pos_before.z,
                vel_before.x, vel_before.y, vel_before.z, ion.pos.x, ion.pos.y,
                ion.pos.z, ion.vel.x, ion.vel.y, ion.vel.z, pos_after.x,
                pos_after.y, pos_after.z, vel_after.x, vel_after.y,
                vel_after.z);
    rk45_state_debug_counter++;
  }
}

void RK45Strategy::step_adaptive(IonState &ion, double t, double &dt_inout,
                                 const physics::ForceRegistry &force_registry,
                                 const std::vector<IonState> &all_ions) {
  (void)all_ions;
  const double dt_initial = dt_inout;

  // Domain pointer is retrieved for consistency; current implementation
  // does not use domain fields directly.
  const config::DomainConfig *domain = force_registry.domain();
  if (!domain) {
    throw std::runtime_error(
        "RK45Strategy: ForceRegistry has no domain configured");
  }
  const double dt_min = std::max(dt_initial * config_.min_step_factor,
                                 config_.absolute_min_step_s);
  const double dt_max = dt_initial * config_.max_step_factor;

  double dt = dt_inout;
  bool step_accepted = false;
  int attempts = 0;

  IonState y0 = ion; // Save initial state
  auto &st = state_for(0);

  while (!step_accepted && attempts < MAX_REJECT_ATTEMPTS) {
    // =====================================================================
    // Compute 7 RK stages (k1...k7)
    // For 2nd order ODE: x' = v, v' = a = f(x, v, t) / m
    // Each stage k_i stores (velocity, acceleration) pair
    // =====================================================================

    Vec3 k1_v, k1_a;
    Vec3 k2_v, k2_a;
    Vec3 k3_v, k3_a;
    Vec3 k4_v, k4_a;
    Vec3 k5_v, k5_a;
    Vec3 k6_v, k6_a;
    Vec3 k7_v, k7_a;

    // Stage 1: k1 = f(t, y)
    k1_v = y0.vel; // dx/dt = v
    if (st.fsal_available) {
      // Reuse k7 from previous step (FSAL property)
      k1_a.x = st.k1_ax;
      k1_a.y = st.k1_ay;
      k1_a.z = st.k1_az;
    } else {
      double ax, ay, az;
      compute_acceleration_state(ax, ay, az, y0, t, force_registry);
      k1_a = Vec3{ax, ay, az}; // dv/dt = a
    }

    // Stage 2: k2 = f(t + c2*dt, y + dt*(a21*k1))
    IonState y2 = y0;
    y2.pos += k1_v * (dt * a21); // x + v*dt*a21
    y2.vel += k1_a * (dt * a21); // v + a*dt*a21
    k2_v = y2.vel;
    double ax2, ay2, az2;
    compute_acceleration_state(ax2, ay2, az2, y2, t + c2 * dt, force_registry);
    k2_a = Vec3{ax2, ay2, az2};

    // Stage 3: k3 = f(t + c3*dt, y + dt*(a31*k1 + a32*k2))
    IonState y3 = y0;
    y3.pos += (k1_v * a31 + k2_v * a32) * dt;
    y3.vel += (k1_a * a31 + k2_a * a32) * dt;
    k3_v = y3.vel;
    double ax3, ay3, az3;
    compute_acceleration_state(ax3, ay3, az3, y3, t + c3 * dt, force_registry);
    k3_a = Vec3{ax3, ay3, az3};

    // Stage 4: k4 = f(t + c4*dt, y + dt*(a41*k1 + a42*k2 + a43*k3))
    IonState y4_temp = y0;
    y4_temp.pos += (k1_v * a41 + k2_v * a42 + k3_v * a43) * dt;
    y4_temp.vel += (k1_a * a41 + k2_a * a42 + k3_a * a43) * dt;
    k4_v = y4_temp.vel;
    double ax4, ay4, az4;
    compute_acceleration_state(ax4, ay4, az4, y4_temp, t + c4 * dt,
                               force_registry);
    k4_a = Vec3{ax4, ay4, az4};

    // Stage 5: k5 = f(t + c5*dt, y + dt*(a51*k1 + ... + a54*k4))
    IonState y5_temp = y0;
    y5_temp.pos += (k1_v * a51 + k2_v * a52 + k3_v * a53 + k4_v * a54) * dt;
    y5_temp.vel += (k1_a * a51 + k2_a * a52 + k3_a * a53 + k4_a * a54) * dt;
    k5_v = y5_temp.vel;
    double ax5, ay5, az5;
    compute_acceleration_state(ax5, ay5, az5, y5_temp, t + c5 * dt,
                               force_registry);
    k5_a = Vec3{ax5, ay5, az5};

    // Stage 6: k6 = f(t + c6*dt, y + dt*(a61*k1 + ... + a65*k5))
    IonState y6 = y0;
    y6.pos +=
        (k1_v * a61 + k2_v * a62 + k3_v * a63 + k4_v * a64 + k5_v * a65) * dt;
    y6.vel +=
        (k1_a * a61 + k2_a * a62 + k3_a * a63 + k4_a * a64 + k5_a * a65) * dt;
    k6_v = y6.vel;
    double ax6, ay6, az6;
    compute_acceleration_state(ax6, ay6, az6, y6, t + c6 * dt, force_registry);
    k6_a = Vec3{ax6, ay6, az6};

    // Stage 7: k7 = f(t + dt, y + dt*(a71*k1 + ... + a76*k6))
    // This becomes the 4th-order solution AND next step's k1 (FSAL)
    IonState y7 = y0;
    y7.pos += (k1_v * a71 + k2_v * a72 + k3_v * a73 + k4_v * a74 + k5_v * a75 +
               k6_v * a76) *
              dt;
    y7.vel += (k1_a * a71 + k2_a * a72 + k3_a * a73 + k4_a * a74 + k5_a * a75 +
               k6_a * a76) *
              dt;
    k7_v = y7.vel;
    double ax7, ay7, az7;
    compute_acceleration_state(ax7, ay7, az7, y7, t + c7 * dt, force_registry);
    k7_a = Vec3{ax7, ay7, az7};

    // =====================================================================
    // Construct 4th and 5th order solutions
    // =====================================================================

    // 4th order solution (b4 coefficients) - SAME as y7 due to FSAL
    IonState y4 = y7;

    // 5th order solution (b5 coefficients)
    IonState y5 = y0;
    y5.pos += (k1_v * b51 + k2_v * b52 + k3_v * b53 + k4_v * b54 + k5_v * b55 +
               k6_v * b56 + k7_v * b57) *
              dt;
    y5.vel += (k1_a * b51 + k2_a * b52 + k3_a * b53 + k4_a * b54 + k5_a * b55 +
               k6_a * b56 + k7_a * b57) *
              dt;

    // =====================================================================
    // Error estimation and step acceptance
    // =====================================================================

    double error = estimate_error(y4, y5, y0);
    const bool force_accept = config_.accept_at_dt_min &&
        (dt <= dt_min * DT_MIN_TOLERANCE);

    if (error <= ERROR_ACCEPTANCE_THRESHOLD || force_accept) {
      // Accept step
      ion = y4; // Use 4th-order solution
      step_accepted = true;

      // Store k7 for next step's k1 (FSAL)
      st.k1_ax = k7_a.x;
      st.k1_ay = k7_a.y;
      st.k1_az = k7_a.z;
      st.fsal_available = true;

      // Update statistics
      if (stats_enabled_) {
        stats_.accepted_steps++;
        stats_.min_step_used = std::min(stats_.min_step_used, dt);
        stats_.max_step_used = std::max(stats_.max_step_used, dt);
        stats_.avg_error =
            (stats_.avg_error * (stats_.accepted_steps - 1) + error) /
            stats_.accepted_steps;
        stats_.sum_step_used += dt;
      }

      // Compute new step size for next step
      double dt_next =
          compute_new_step(dt, error, dt_min, dt_max, st.last_error);
      // Expose used dt to caller; next-step hint can be derived if needed
      dt_inout = dt;
      last_dt_used_ = dt;
      last_dt_suggested_ = dt_next;
      (void)dt_next;

      st.last_error = error;

    } else {
      // Reject step, reduce dt, retry
      if (stats_enabled_) {
        stats_.rejected_steps++;
      }

      double dt_next =
          compute_new_step(dt, error, dt_min, dt_max, st.last_error);
      dt = dt_next;

      attempts++;

      // Reset FSAL flag (must recompute k1)
      st.fsal_available = false;
    }
  }

  if (!step_accepted) {
    throw std::runtime_error("RK45Strategy: Failed to converge after " +
                             std::to_string(MAX_REJECT_ATTEMPTS) + " attempts");
  }
}

} // namespace integrator
} // namespace ICARION
