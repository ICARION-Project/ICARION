#include "core/integrator/integrator_helpers.h"
#include "core/log/Logger.h"
#include "fieldsolver/utils/field_update_api.h"
#include <cmath>
#include <random>
#include "H5Cpp.h"
#include "utils/vdw_radii.h"
#include <fstream>
#include <mutex>
#include <thread>
#include <atomic>
#include "core/utils/safety/numericalSafetyGuards.h"
#include "core/utils/safety/numericalSafetyLogger.h"
#include <unistd.h>
// Field provider wrappers
#include "fieldsolver/utils/GridFieldProvider.h"
// Numerical safety guards
#include "core/utils/safety/numericalSafetyGuards.h"

namespace integrator_helpers {

void log_step(double t, double t_end, int& step_count, double& next_write_time, double write_dt,
             const std::vector<IonState>& ions, std::vector<double>& times_buffer,
             std::vector<std::vector<IonState>>& trajectory_buffer, size_t buffer_max,
             const std::string& filename) {
    if (step_count == 0 || t >= next_write_time || t >= t_end) {
        times_buffer.push_back(t);
        trajectory_buffer.push_back(ions);
        if (write_dt > 0.0) {
            while (next_write_time <= t) next_write_time += write_dt;
        } else {
            next_write_time = t + 1.0;
        }
    }
    if (times_buffer.size() >= buffer_max) {
        ICARION::io::append_to_HDF5(filename, times_buffer, trajectory_buffer);
        times_buffer.clear();
        trajectory_buffer.clear();
    }
    step_count++;
}

void integrate_one_step(
    IonState& ion,
    const GlobalParams& gParams,
    const InstrumentDomain& dom,
    const std::unordered_map<std::string, Species>& speciesDB,
    const std::vector<ReactionEntry>& reaction_list,
    EhssRng& local_rng,
    const RK45Settings& rk45,
    const std::vector<IonState>& current_ions
    , const IFieldProvider* field_provider
) {
    double t_local = ion.t;
    double t_target = ion.t + ion.dt;
    double dt_local = ion.dt;

    // --- Prepare optional field provider if a snapshot or field array exists ---
    // If a provider was explicitly passed, prefer it. Otherwise, attempt to
    // construct a GridFieldProvider from a FieldServer snapshot or FieldArray.
    std::unique_ptr<IFieldProvider> provider_owner;
    const IFieldProvider* provider_ptr = nullptr;
    FieldSnapshot snapshot_local;
    if (field_provider != nullptr) {
        provider_ptr = field_provider;
    } else {
        if (dom.use_grid_field && dom.fieldServer != nullptr) {
            snapshot_local = dom.fieldServer->get_snapshot();
            provider_owner.reset(new GridFieldProvider(&snapshot_local));
        } else if (dom.fieldArrayLoaded) {
            provider_owner.reset(new GridFieldProvider(&dom.fieldArray));
        }
        provider_ptr = provider_owner ? provider_owner.get() : nullptr;
    }

    // Create safety configuration for numerical stability
    static ICARION::safety::NumericalSafetyConfig safety_config;
    safety_config.enable_nan_checks = true;
    safety_config.enable_bounds_checks = false; // Position bounds check can be too restrictive
    safety_config.throw_on_violation = true;
    safety_config.enable_logging = false; // Disable for performance unless debugging
    
    // Store initial ion state for potential recovery
    IonState ion_backup = ion;
    static int global_step_counter = 0;
    static int global_ion_counter = 0;
    global_step_counter++;

    switch (dom.solver_type) {
        case SolverType::Boris: {
            Vec3 E_field = (dom.instrument == Instrument::FTICR) ?
                FTICRField(ion, dom.DC.radial_V, std::sqrt((dom.geom.length_m*dom.geom.length_m/8.0)+(dom.geom.radius_m*dom.geom.radius_m/4.0)), dom.geom.length_m)
                : ElectricFieldVec(ion, dom, t_local);
            Vec3 B_field = MagneticFieldVec(ion, dom);
            
            // NUMERICAL SAFETY: Check field values
            ICARION_CHECK_FINITE(E_field, "Electric field calculation");
            ICARION_CHECK_FINITE(B_field, "Magnetic field calculation");
            
            const double q = ion.ion_charge_C;
            const double m = ion.mass_kg;
            Vec3 v_minus = ion.vel + (E_field * q / m) * (0.5 * dt_local);
            Vec3 t_vec = (B_field * q / m) * (0.5 * dt_local);
            double t_mag2 = dot(t_vec, t_vec);
            Vec3 s_vec = t_vec * (2.0 / (1.0 + t_mag2));
            Vec3 v_prime = v_minus + cross(v_minus, t_vec);
            Vec3 v_plus = v_minus + cross(v_prime, s_vec);
            // NUMERICAL SAFETY: Comprehensive safety check with logging
            ICARION_SAFETY_TIMER();
            
            if (safety_config.enable_nan_checks) {
                if (std::isnan(v_minus.x) || std::isnan(v_minus.y) || std::isnan(v_minus.z)) {
                    ICARION_SAFETY_LOG_NAN(global_ion_counter, global_step_counter, t_local,
                                           ion.pos, v_minus, "Boris v_minus calculation");
                    ICARION_SAFETY_MARK_VIOLATION();
                }
                if (std::isnan(v_plus.x) || std::isnan(v_plus.y) || std::isnan(v_plus.z)) {
                    ICARION_SAFETY_LOG_NAN(global_ion_counter, global_step_counter, t_local,
                                           ion.pos, v_plus, "Boris v_plus calculation");
                    ICARION_SAFETY_MARK_VIOLATION();
                }
            }
            
            ion.vel = v_plus + (E_field * q / m) * (0.5 * dt_local);
            ion.pos += ion.vel * dt_local;
            ion.t = t_target;
            
            // Final safety check with comprehensive logging
            if (!ICARION::safety::check_ion_safety(ion, safety_config, "Boris integration", 
                                                   global_step_counter, global_ion_counter)) {
                ICARION_SAFETY_MARK_VIOLATION();
                ICARION::log::Logger::integrator()->warn("Invalid ion state detected after Boris integration");
            }
            
            break;
        }
        case SolverType::RK4: {
            ICARION_SAFETY_TIMER();
            
            IonState k1 = compute_accelerations(t_local, ion, gParams, dom, provider_ptr, current_ions);
            if (safety_config.enable_nan_checks) {
                if (std::isnan(k1.pos.x) || std::isnan(k1.pos.y) || std::isnan(k1.pos.z) ||
                    std::isnan(k1.vel.x) || std::isnan(k1.vel.y) || std::isnan(k1.vel.z)) {
                    ICARION_SAFETY_LOG_NAN(global_ion_counter, global_step_counter, t_local,
                                           ion.pos, ion.vel, "RK4 k1 calculation");
                    ICARION_SAFETY_MARK_VIOLATION();
                }
            }
            
            IonState k2 = compute_accelerations(t_local + 0.5 * dt_local, ion + k1 * 0.5 * dt_local, gParams, dom, provider_ptr, current_ions);
            if (safety_config.enable_nan_checks) {
                if (std::isnan(k2.pos.x) || std::isnan(k2.pos.y) || std::isnan(k2.pos.z) ||
                    std::isnan(k2.vel.x) || std::isnan(k2.vel.y) || std::isnan(k2.vel.z)) {
                    ICARION_SAFETY_LOG_NAN(global_ion_counter, global_step_counter, t_local,
                                           ion.pos, ion.vel, "RK4 k2 calculation");
                    ICARION_SAFETY_MARK_VIOLATION();
                }
            }
            
            IonState k3 = compute_accelerations(t_local + 0.5 * dt_local, ion + k2 * 0.5 * dt_local, gParams, dom, provider_ptr, current_ions);
            if (safety_config.enable_nan_checks) {
                if (std::isnan(k3.pos.x) || std::isnan(k3.pos.y) || std::isnan(k3.pos.z) ||
                    std::isnan(k3.vel.x) || std::isnan(k3.vel.y) || std::isnan(k3.vel.z)) {
                    ICARION_SAFETY_LOG_NAN(global_ion_counter, global_step_counter, t_local,
                                           ion.pos, ion.vel, "RK4 k3 calculation");
                    ICARION_SAFETY_MARK_VIOLATION();
                }
            }
            
            IonState k4 = compute_accelerations(t_local + dt_local, ion + k3 * dt_local, gParams, dom, provider_ptr, current_ions);
            if (safety_config.enable_nan_checks) {
                if (std::isnan(k4.pos.x) || std::isnan(k4.pos.y) || std::isnan(k4.pos.z) ||
                    std::isnan(k4.vel.x) || std::isnan(k4.vel.y) || std::isnan(k4.vel.z)) {
                    ICARION_SAFETY_LOG_NAN(global_ion_counter, global_step_counter, t_local,
                                           ion.pos, ion.vel, "RK4 k4 calculation");
                    ICARION_SAFETY_MARK_VIOLATION();
                }
            }
            
            ion += (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt_local / 6.0);
            ion.t = t_target;
            
            // NUMERICAL SAFETY: Final safety check with comprehensive logging
            if (!ICARION::safety::check_ion_safety(ion, safety_config, "RK4 integration", 
                                                   global_step_counter, global_ion_counter)) {
                ICARION_SAFETY_MARK_VIOLATION();
                ICARION::log::Logger::integrator()->warn("Invalid ion state detected after RK4 integration");
            }
            if (!ICARION::safety::check_ion_safety(ion, safety_config, "RK4 integration", 
                                                   global_step_counter, global_ion_counter)) {
                // Ion state is invalid - this should have thrown an exception if configured
                ICARION::log::Logger::integrator()->warn("Invalid ion state detected after RK4 integration");
            }
            
            break;
        }
        case SolverType::RK45: {
            ICARION_SAFETY_TIMER();
            int rk45_substeps = 0;
            int rejected_steps = 0;
            
            while (t_local < t_target) {
                if (t_local + dt_local > t_target) dt_local = t_target - t_local;
                
                IonState y5, y4;
                ::integrator_helpers::rk45_dp_step(t_local, ion, dt_local, gParams, dom, y5, y4, current_ions, provider_ptr);
                
                // NUMERICAL SAFETY: Check intermediate RK45 results with logging
                if (safety_config.enable_nan_checks) {
                    if (!ICARION::safety::is_finite_vec3(y5.pos) || !ICARION::safety::is_finite_vec3(y5.vel)) {
                        ICARION_SAFETY_LOG_NAN(global_ion_counter, global_step_counter, t_local,
                                               y5.pos, y5.vel, "RK45 y5 calculation substep " + std::to_string(rk45_substeps));
                        ICARION_SAFETY_MARK_VIOLATION();
                        throw std::runtime_error("NaN/Inf in RK45 y5 at substep " + std::to_string(rk45_substeps));
                    }
                    if (!ICARION::safety::is_finite_vec3(y4.pos) || !ICARION::safety::is_finite_vec3(y4.vel)) {
                        ICARION_SAFETY_LOG_NAN(global_ion_counter, global_step_counter, t_local,
                                               y4.pos, y4.vel, "RK45 y4 calculation substep " + std::to_string(rk45_substeps));
                        ICARION_SAFETY_MARK_VIOLATION();
                        throw std::runtime_error("NaN/Inf in RK45 y4 at substep " + std::to_string(rk45_substeps));
                    }
                }
                
                double err = ::integrator_helpers::rk_error_norm(y5, y4, ion, rk45);
                
                // NUMERICAL SAFETY: Check error norm with logging
                if (!std::isfinite(err) || std::isnan(err)) {
                    ICARION_SAFETY_LOG_NAN(global_ion_counter, global_step_counter, t_local,
                                           ion.pos, ion.vel, "RK45 error norm calculation, err=" + std::to_string(err));
                    ICARION_SAFETY_MARK_VIOLATION();
                    throw std::runtime_error("Invalid error norm in RK45 at substep " + std::to_string(rk45_substeps));
                }
                
                if (err <= 1.0) { 
                    ion = y5; 
                    t_local += dt_local; 
                    
                    // NUMERICAL SAFETY: Check accepted step with logging
                    if (!ICARION::safety::check_ion_safety(ion, safety_config, 
                                                          "RK45 accepted step " + std::to_string(rk45_substeps), 
                                                          global_step_counter, global_ion_counter)) {
                        ICARION_SAFETY_MARK_VIOLATION();
                        ICARION::log::Logger::integrator()->warn("Invalid state after RK45 accepted step");
                    }
                } else {
                    // Step rejected - log for analysis
                    rejected_steps++;
                    double factor = rk45.safety * std::pow(std::max(1e-16, 1.0 / err), 0.25);
                    double new_dt = std::max(rk45.dt_min, std::min(dt_local * factor, dt_local * rk45.max_factor));
                    
                    ICARION_SAFETY_LOG_REJECTED_STEP(global_ion_counter, global_step_counter, t_local,
                                                     dt_local, new_dt, "Error too large: " + std::to_string(err));
                    dt_local = new_dt;
                }
                
                // Update timestep for next iteration if step was accepted
                if (err <= 1.0) {
                    double factor = rk45.safety * std::pow(std::max(1e-16, 1.0 / err), 0.25);
                    factor = std::max(rk45.min_factor, std::min(factor, rk45.max_factor));
                    dt_local = std::max(rk45.dt_min, dt_local * factor);
                }
                
                rk45_substeps++;
                
                // NUMERICAL SAFETY: Prevent infinite loops with detailed logging
                if (rk45_substeps > rk45.max_rejects) {
                    std::string error_msg = "RK45 exceeded maximum substeps (" + std::to_string(rk45.max_rejects) + 
                                          ") at global step " + std::to_string(global_step_counter) + 
                                          " for ion " + std::to_string(global_ion_counter) + 
                                          ". dt_local=" + std::to_string(dt_local) +
                                          ", rejected_steps=" + std::to_string(rejected_steps);
                    ICARION_SAFETY_LOG_REJECTED_STEP(global_ion_counter, global_step_counter, t_local,
                                                     dt_local, 0.0, "Maximum substeps exceeded");
                    ICARION_SAFETY_MARK_VIOLATION();
                    throw std::runtime_error(error_msg);
                }
            }
            ion.t = t_target;
            
            // NUMERICAL SAFETY: Final check after RK45 completion
            if (!ICARION::safety::check_ion_safety(ion, safety_config, "RK45 completion", 
                                                   global_step_counter, global_ion_counter)) {
                ICARION_SAFETY_MARK_VIOLATION();
                ICARION::log::Logger::integrator()->warn("Invalid ion state after RK45 completion");
            }
            
            break;
        }
    }
    
    global_ion_counter++;
}

// NOTE: load_geometry() and compute_CCS_from_geometry() removed in Phase 2D refactor
// Geometry loading now in src/core/physics/collisions/utils.{h,cpp}
// See load_geometry_from_file() and MolecularGeometry struct

// NOTE: handle_collision() removed in Phase 2D refactor
// Collision handling now delegated to ICollisionHandler implementations
// See CollisionHandlerFactory and handler->handle_collision() in integrate_trajectory()

void handle_reaction(IonState& y, EhssRng& rng, double dt, const GlobalParams& gParams,
                     const std::unordered_map<std::string, Species>& speciesDB,
                     const std::vector<ReactionEntry>& reaction_list) {
    if (!gParams.enable_reactions) return;
    const std::string current_species = y.species_id;
    for (const auto& rxn : reaction_list) {
        if (rxn.reactant != current_species) continue;
        double k_eff = rxn.rate_constant;
        for (const auto& term : rxn.order) {
            double conc_density = (rxn.neutral_concentration > 1e6) ? rxn.neutral_concentration : rxn.neutral_concentration * y.domain_particle_density_m3;
            if (conc_density < 0.0) conc_density = 0.0;
            k_eff *= std::pow(conc_density, term.exponent);
        }
        double P = 1.0 - std::exp(-k_eff * dt);
        if (rng.uniform01() < P) {
            const auto& prod = speciesDB.at(rxn.product);
            y.species_id = rxn.product; y.mass_kg = prod.mass_kg; y.reduced_mobility_cm2_Vs = prod.mobility; y.ion_charge_C = prod.charge;
            y.CCS_m2 = 3.0 / 16.0 / LOSCHMIDT_CONSTANT * y.ion_charge_C * std::sqrt(2.0 * M_PI / (BOLTZMANN_CONSTANT * y.domain_temperature_K * (y.mass_kg * y.domain_neutral_mass_kg) / (y.mass_kg + y.domain_neutral_mass_kg))) * 1.0 / (y.reduced_mobility_cm2_Vs * 1e-4);
        }
    }
}

double rk_error_norm(const IonState& y5, const IonState& y4, const IonState& y, const RK45Settings& rk45) {
    double s = 0.0;
    s += std::pow((y5.pos.x - y4.pos.x) / (rk45.absTol + rk45.relTol * std::max(std::fabs(y.pos.x), std::fabs(y5.pos.x))), 2);
    s += std::pow((y5.pos.y - y4.pos.y) / (rk45.absTol + rk45.relTol * std::max(std::fabs(y.pos.y), std::fabs(y5.pos.y))), 2);
    s += std::pow((y5.pos.z - y4.pos.z) / (rk45.absTol + rk45.relTol * std::max(std::fabs(y.pos.z), std::fabs(y5.pos.z))), 2);
    s += std::pow((y5.vel.x - y4.vel.x) / (rk45.absTol + rk45.relTol * std::max(std::fabs(y.vel.x), std::fabs(y5.vel.x))), 2);
    s += std::pow((y5.vel.y - y4.vel.y) / (rk45.absTol + rk45.relTol * std::max(std::fabs(y.vel.y), std::fabs(y5.vel.y))), 2);
    s += std::pow((y5.vel.z - y4.vel.z) / (rk45.absTol + rk45.relTol * std::max(std::fabs(y.vel.z), std::fabs(y5.vel.z))), 2);
    return std::sqrt(s / 6.0);
}

void rk45_dp_step(double t, const IonState& y, double h, const GlobalParams& gParams,
                  const InstrumentDomain& dom, IonState& y5, IonState& y4, const std::vector<IonState>& current_ions,
                  const IFieldProvider* field_provider) {
    const double c2 = 1.0/5.0, c3 = 3.0/10.0, c4 = 4.0/5.0, c5 = 8.0/9.0, c6 = 1.0, c7 = 1.0;
    const double a21 = 1.0/5.0;
    const double a31 = 3.0/40.0, a32 = 9.0/40.0;
    const double a41 = 44.0/45.0, a42 = -56.0/15.0, a43 = 32.0/9.0;
    const double a51 = 19372.0/6561.0, a52 = -25360.0/2187.0, a53 = 64448.0/6561.0, a54 = -212.0/729.0;
    const double a61 = 9017.0/3168.0, a62 = -355.0/33.0, a63 = 46732.0/5247.0, a64 = 49.0/176.0, a65 = -5103.0/18656.0;
    const double a71 = 35.0/384.0, a72 = 0.0, a73 = 500.0/1113.0, a74 = 125.0/192.0, a75 = -2187.0/6784.0, a76 = 11.0/84.0;
    const double b1 = 35.0/384.0, b2 = 0.0, b3 = 500.0/1113.0, b4 = 125.0/192.0, b5 = -2187.0/6784.0, b6 = 11.0/84.0, b7 = 0.0;
    const double bs1 = 5179.0/57600.0, bs2 = 0.0, bs3 = 7571.0/16695.0, bs4 = 393.0/640.0, bs5 = -92097.0/339200.0, bs6 = 187.0/2100.0, bs7 = 1.0/40.0;

    IonState k1 = compute_accelerations(t, y, gParams, dom, field_provider, current_ions) * h;
    IonState k2 = compute_accelerations(t + c2*h, y + k1 * a21, gParams, dom, field_provider, current_ions) * h;
    IonState k3 = compute_accelerations(t + c3*h, y + k1 * a31 + k2 * a32, gParams, dom, field_provider, current_ions) * h;
    IonState k4 = compute_accelerations(t + c4*h, y + k1 * a41 + k2 * a42 + k3 * a43, gParams, dom, field_provider, current_ions) * h;
    IonState k5 = compute_accelerations(t + c5*h, y + k1 * a51 + k2 * a52 + k3 * a53 + k4 * a54, gParams, dom, field_provider, current_ions) * h;
    IonState k6 = compute_accelerations(t + c6*h, y + k1 * a61 + k2 * a62 + k3 * a63 + k4 * a64 + k5 * a65, gParams, dom, field_provider, current_ions) * h;
    IonState k7 = compute_accelerations(t + c7*h, y + k1 * a71 + k2 * a72 + k3 * a73 + k4 * a74 + k5 * a75 + k6 * a76, gParams, dom, field_provider, current_ions) * h;

    y5 = y + k1 * b1 + k2 * b2 + k3 * b3 + k4 * b4 + k5 * b5 + k6 * b6 + k7 * b7;
    y4 = y + k1 * bs1 + k2 * bs2 + k3 * bs3 + k4 * bs4 + k5 * bs5 + k6 * bs6 + k7 * bs7;
}

} // namespace integrator_helpers
