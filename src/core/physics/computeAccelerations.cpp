/**
 * =====================================================================
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        computeAccelerations.cpp
 *   @brief       Computes ion accelerations and updated velocities.
 *
 *   @details
 *   Provides the function to calculate the time derivatives of an ion’s state
 *   vector for use in numerical integration (e.g., Runge-Kutta).
 *
 *   Acceleration contributions include:
 *   - Electric fields (DC, RF, AC, Orbitrap)
 *   - Collision damping (HardSphere, Langevin, Friction, EHSS, HSMC)
 *   - Background gas flow (adds to velocity)
 *
 *   Supports multiple instruments: LQIT, IMS, SIFDT-MS, Orbitrap, QuadrupoleRF, TOF.
 *
 *   @note
 *   Returned IonState contains updated velocity and acceleration.
 *   All units are SI; fields in V/m. Collision events and boundary checks
 *   are handled elsewhere.
 *
 *
 *   @date        2025-10-06
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#include "computeAccelerations.h"
//#include "core/physics/spacecharge/spaceChargeForces.h"  // Space charge force calculations
#include "core/physics/fields/fieldSampling.h"  // Core field sampling utilities

// Field server: Core-only uses shim, full build uses fieldsolver
#include "fieldsolver/utils/field_update_api.h"

#include "core/physics/fields/physics_fields_shared.h"
#include "core/physics/fields/physics_math_shared.h"

#include <cmath>
#include <functional>

namespace ICARION {
namespace physics {

/**
 * @brief Compute updated velocity and acceleration for a single ion.
 *
 * This function sums all relevant forces acting on the ion, converts them
 * to acceleration, and corrects velocity for background gas flow.
 *
 * @param[in] t      Current simulation time [s].
 * @param[in] y      Current ion state (position, velocity, mass, charge, etc.).
 * @param[in] gParams Global simulation parameters (collision model, etc.).
 * @param[in] dom    Instrument domain parameters (geometry, field settings, etc.).
 *
 * @return IonState  Updated ion state containing:
 *                   - vel: corrected velocity including gas flow [m/s]
 *                   - pos: current velocity (used as derivative in integration)
 *                   - additional properties unchanged
 */
IonState compute_accelerations(double t, const IonState& y, const GlobalParams& gParams, const InstrumentDomain& dom, 
                               const IFieldProvider* field_provider, const std::vector<IonState>& current_ions) {
    // --- Initialize total acceleration vector ---
    Vec3 acc_total{0.0, 0.0, 0.0};
    
    // --- Run guard checks for domain and global parameters ---
    run_guard_check_domain(dom, gParams.collisionModel);
    run_guard_check_global(gParams);

    // --- Validate current ion state ---
    y.validate();

    // --- Configure field contributions ---
    std::vector<std::function<Vec3(const IonState&, double)>> fieldModules;

    // Prefer external field_provider if supplied
    if (field_provider != nullptr) {
        fieldModules.push_back([field_provider](const IonState& ion, double t) {
            return field_provider->get_E(ion.pos);
        });
    }
    // Check if GridField mode is active
    else if (dom.use_grid_field && dom.fieldServer != nullptr) {
        // Use live field from FieldServer
        fieldModules.push_back([&dom](const IonState& ion, double t) {
            using ICARION::fields::FieldSnapshot;
            using ICARION::fields::Vec3d;
            using ICARION::fields::sample_field;
            
            FieldSnapshot snapshot = dom.fieldServer->get_snapshot();
            // Convert Vec3 to Vec3d for sampling
            Vec3d pos_d{ion.pos.x, ion.pos.y, ion.pos.z};
            Vec3d E_d = sample_field(snapshot, pos_d);
            // Convert back to Vec3
            return Vec3{E_d.x, E_d.y, E_d.z};
        });
    } 
    else if (!dom.FA_terms.empty() && dom.fieldArrayLoaded) {
            fieldModules.push_back([&dom](const IonState& ion, double t) {
                Vec3 E{0.0, 0.0, 0.0};
                for (const auto& term : dom.FA_terms) {
                    if (!term.loaded || !term.field.is_valid()) continue;
                    double scale = 0.0;
                    switch (term.kind) {
                        case InstrumentDomain::FAScaleKind::Constant:
                            scale = term.constant; // in Volts
                            break;
                        case InstrumentDomain::FAScaleKind::DC_Axial:
                            scale = dom.DC.axial_V;
                            break;
                        case InstrumentDomain::FAScaleKind::DC_Quad:
                            scale = dom.DC.quad_V;
                            break;
                        case InstrumentDomain::FAScaleKind::DC_Radial:
                            scale = dom.DC.radial_V;
                            break;
                        case InstrumentDomain::FAScaleKind::RF: {
                            double amp = dom.RF.voltage_V * term.constant; // allow additional per-term gain
                            double omega = (term.frequency_Hz > 0.0) ? (2.0 * M_PI * term.frequency_Hz) : dom.RF.angular_frequency_rad_s;
                            double phase = dom.RF.phase_rad + term.phase_rad;
                            scale = amp * std::sin(omega * t + phase);
                            break;
                        }
                    }
                    if (scale != 0.0) {
                        E += interpolate_field(term.field, ion.pos) * scale; // fields are for 1 V
                    }
                }
                return E;
            });
        } else if (!dom.FA_file.empty() && dom.fieldArrayLoaded) {
        fieldModules.push_back([&dom](const IonState& ion, double t) {
            return interpolate_field(dom.fieldArray, ion.pos);
        });
    } else {
        switch(dom.instrument) {
            case Instrument::LQIT:
                // (1) RF quadrupole radial field
                fieldModules.push_back([&dom](const IonState& ion, double t) {
                    return RFField(ion, dom.RF.voltage_V, dom.DC.quad_V,
                                            dom.RF.angular_frequency_rad_s, dom.geom.radius_m, t);
                });
                // (2) AC excitation field along X
                fieldModules.push_back([&dom](const IonState& ion, double t) {
                    double voltage = dom.AC.voltage_V; // default constant

                    // If a voltage_time_table is provided, use linear interpolation
                    // between samples to obtain the instantaneous AC voltage.
                    if (!dom.AC.voltage_time_table.empty()) {
                        const auto &tab = dom.AC.voltage_time_table;
                        if (t <= tab.front().first) {
                            voltage = tab.front().second;
                        } else if (t >= tab.back().first) {
                            voltage = tab.back().second;
                        } else {
                            auto it = std::upper_bound(tab.begin(), tab.end(), t, [](double tt, const std::pair<double,double>& p){ return tt < p.first; });
                            if (it != tab.begin() && it != tab.end()) {
                                auto hi = it;
                                auto lo = it - 1;
                                double t0 = lo->first;
                                double v0 = lo->second;
                                double t1 = hi->first;
                                double v1 = hi->second;
                                double alpha = (t - t0) / (t1 - t0);
                                voltage = v0 + alpha * (v1 - v0);
                            }
                        }
                    } else if (dom.AC.enable_voltage_sweep) {
                        double t_rel = t - dom.AC.start_time_s;
                        if (t_rel <= 0.0) {
                            voltage = dom.AC.voltage_V; // before ramp starts
                        } else if (t_rel < dom.AC.rise_time_s) {
                            voltage = dom.AC.voltage_V 
                                    + t_rel * dom.AC.amplitude_slope_V_s; // ramping
                        } else {
                            voltage = dom.AC.voltage_V 
                                    + dom.AC.amplitude_slope_V_s * dom.AC.rise_time_s; // max reached
                        }
                    }

                    // Frequency sweep (linear in Hz) support. If enabled, compute
                    // the current frequency as start_freq + slope * t_rel and
                    // convert to angular frequency. Fall back to the static
                    // angular_frequency_rad_s when sweep not enabled or params
                    // missing.
                    double omega = dom.AC.angular_frequency_rad_s;
                    if (dom.AC.enable_frequency_sweep) {
                        double t_rel = t - dom.AC.start_time_s;
                        // use shared helper for host/device parity
                        omega = compute_ac_omega_from_hz(dom.AC.ac_start_freq_Hz,
                                                         dom.AC.ac_sweep_slope_Hz_per_s,
                                                         t_rel,
                                                         dom.AC.enable_frequency_sweep ? 1 : 0);
                    }

                    // LQIT-lock semantics on host: match device behavior
                    double t_ac = t;
                    if (dom.AC.lqit_lock_enable) {
                        double phase = dom.AC.lqit_lock_phase_rad;
                        if (fabs(omega) > 1e-30) t_ac += phase / omega;
                        double bw = dom.AC.lqit_lock_bandwidth_Hz;
                        if (bw > 0.0) {
                            const double TWO_PI = 2.0 * 3.14159265358979323846;
                            double freq_Hz = omega / TWO_PI;
                            double rf_freq = dom.RF.angular_frequency_rad_s / TWO_PI;
                            double sigma = bw * 0.5;
                            double delta = fabs(freq_Hz - rf_freq);
                            double amp_scale = 1.0;
                            if (sigma > 0.0) amp_scale = exp_shared(- (delta*delta) / (2.0 * sigma * sigma));
                            voltage *= amp_scale;
                        }
                    }

                    return ACField(ion, voltage,
                                            omega, dom.geom.radius_m, t_ac, Vec3(1, 0, 0));
                });
                // (3) Axial DC fields at boundaries
                if (y.pos.z > dom.geom.length_m * 0.9){
                    fieldModules.push_back([&dom](const IonState& ion, double t) {
                    return DCField(ion, -dom.DC.axial_V, dom.geom.length_m * 0.1);
                });
                } else if (y.pos.z < dom.geom.length_m * 0.1) {
                    fieldModules.push_back([&dom](const IonState& ion, double t) {
                    return DCField(ion, dom.DC.axial_V, dom.geom.length_m * 0.1);
                });
                }
                break;

            case Instrument::IMS:
            case Instrument::SIFDT_MS:
                // Axial DC field
                fieldModules.push_back([&dom](const IonState& ion, double t) {
                    return DCField(ion, dom.DC.axial_V, dom.geom.length_m);
                });
                // Optional radial RF field
                fieldModules.push_back([&dom](const IonState& ion, double t) {
                    return RFField(ion, dom.RF.voltage_V, dom.DC.quad_V,
                                            dom.RF.angular_frequency_rad_s, dom.geom.radius_m, t);
                });
                break;

            case Instrument::Orbitrap:
                // Radial (X/Y) trapping field and axial (Z) field 
                fieldModules.push_back([&dom](const IonState& ion, double t){
                    double voltage = dom.DC.radial_V; 
                    if (dom.DC.enable_radial_voltage_sweep) {
                        double t_rel = t - dom.DC.radial_start_time_s;
                        if (t_rel <= 0.0) {
                            voltage = dom.DC.radial_V; // before ramp starts
                        } else if (t_rel < dom.DC.radial_rise_time_s) {
                            voltage = dom.DC.radial_V 
                                    + t_rel * dom.DC.radial_slope_V_s; // ramping
                        } else {
                            voltage = dom.DC.radial_V 
                                    + dom.DC.radial_slope_V_s * dom.DC.radial_rise_time_s; // max reached
                        }
                    } 

                    double k = 2.0 * voltage / (dom.geom.radius_char_m * dom.geom.radius_char_m
                                    * log(dom.geom.radius_out_m / dom.geom.radius_in_m)
                                    - 0.5*(dom.geom.radius_out_m * dom.geom.radius_out_m - dom.geom.radius_in_m * dom.geom.radius_in_m));
                    return OrbitrapField(ion, k, dom.geom.radius_char_m, dom.geom.length_m);
                });
                break;
            
            case Instrument::QuadrupoleRF:
                // Radial RF 
                fieldModules.push_back([&dom](const IonState& ion, double t) {
                    return RFField(ion, dom.RF.voltage_V, dom.DC.quad_V,
                                            dom.RF.angular_frequency_rad_s, dom.geom.radius_m, t);
                });
                // Axial DC
                fieldModules.push_back([&dom](const IonState& ion, double t) {
                    return DCField(ion, dom.DC.axial_V, dom.geom.length_m);
                });
                break;

            case Instrument::TOF:
                // axial (Z) DC field (in acceleration source only)
                if (y.pos.z < dom.geom.acc_length_m){
                    fieldModules.push_back([&dom](const IonState& ion, double t) {
                        return DCField(ion, dom.DC.axial_V, dom.geom.acc_length_m);
                    });
                } else {
                    fieldModules.push_back([&dom](const IonState& ion, double t) {
                        //No field
                        return Vec3{0.0, 0.0, 0.0};
                    });
                }
                break;
            
            case Instrument::FT_ICR:

                fieldModules.push_back([&dom](const IonState& ion, double t) {
                    // calculate characteristic distance d
                    double d = std::sqrt((dom.geom.length_m*dom.geom.length_m / 8.0) + 
                        (dom.geom.radius_m*dom.geom.radius_m / 4.0));
                    return FTICRField(ion, dom.DC.radial_V, d, dom.geom.length_m); // in V/m
                });
                break;

            case Instrument::NoFixedInstrument:
                // User-defined fields only; no preset instrument
                // Add DC, RF, AC fields if voltages are non-zero
                if (std::fabs(dom.DC.axial_V) > 1e-12) {
                    fieldModules.push_back([&dom](const IonState& ion, double t) {
                        return DCField(ion, dom.DC.axial_V, dom.geom.length_m); // in V/m
                    });
                }
                if (std::fabs(dom.RF.voltage_V) > 1e-12) {
                    fieldModules.push_back([&dom](const IonState& ion, double t) {
                        return RFField(ion, dom.RF.voltage_V, dom.DC.quad_V,
                                                dom.RF.angular_frequency_rad_s, dom.geom.radius_m, t);
                    });
                }
                if (std::fabs(dom.AC.voltage_V) > 1e-12) {
                    fieldModules.push_back([&dom](const IonState& ion, double t) {
                        double voltage = dom.AC.voltage_V; // default constant

                        // Use voltage_time_table if present (arbitrary waveform)
                        if (!dom.AC.voltage_time_table.empty()) {
                            const auto &tab = dom.AC.voltage_time_table;
                            if (t <= tab.front().first) {
                                voltage = tab.front().second;
                            } else if (t >= tab.back().first) {
                                voltage = tab.back().second;
                            } else {
                                auto it = std::upper_bound(tab.begin(), tab.end(), t, [](double tt, const std::pair<double,double>& p){ return tt < p.first; });
                                if (it != tab.begin() && it != tab.end()) {
                                    auto hi = it;
                                    auto lo = it - 1;
                                    double t0 = lo->first;
                                    double v0 = lo->second;
                                    double t1 = hi->first;
                                    double v1 = hi->second;
                                    double alpha = (t - t0) / (t1 - t0);
                                    voltage = v0 + alpha * (v1 - v0);
                                }
                            }
                        } else if (dom.AC.enable_voltage_sweep) {
                            double t_rel = t - dom.AC.start_time_s;
                            if (t_rel <= 0.0) {
                                voltage = dom.AC.voltage_V; // before ramp starts
                            } else if (t_rel < dom.AC.rise_time_s) {
                                voltage = dom.AC.voltage_V  
                                        + t_rel * dom.AC.amplitude_slope_V_s; // ramping
                            } else {
                                voltage = dom.AC.voltage_V 
                                        + dom.AC.amplitude_slope_V_s * dom.AC.rise_time_s; // max reached
                            }
                        }
                        // Support linear frequency sweep for non-fixed instruments
                        double omega = dom.AC.angular_frequency_rad_s;
                        if (dom.AC.enable_frequency_sweep) {
                            double t_rel = t - dom.AC.start_time_s;
                            omega = compute_ac_omega_from_hz(dom.AC.ac_start_freq_Hz,
                                                             dom.AC.ac_sweep_slope_Hz_per_s,
                                                             t_rel,
                                                             dom.AC.enable_frequency_sweep ? 1 : 0);
                        }
                        return ACField(ion, voltage, omega, dom.geom.radius_m, t, Vec3(1, 0, 0));
                    });
                }
        }
    }
    // --- Sum all electric field contributions ---
    for (auto& f : fieldModules) {
        acc_total += f(y, t) / y.mass_kg * y.ion_charge_C;
    }

    // --- Add acceleration from magnetic field, if enabled ---
    if (dom.B.enabled) {
        Vec3 B_field = MagneticFieldVec(y, dom);
        Vec3 magnetic_force = cross(y.vel, B_field);
        Vec3 F_mag = magnetic_force * y.ion_charge_C;
        acc_total += F_mag / y.mass_kg;
    }

    // --- Add acceleration from collision forces if OU thermalization is not used ---
    if (!gParams.enable_ou_thermalization) {
        acc_total += CollisionForce(y, gParams, dom) / y.mass_kg;
    }   

    // --- Add space charge field, if enabled ---
    if (gParams.enable_space_charge && gParams.spaceChargeSolver != nullptr) {
        Vec3 E_sc = gParams.spaceChargeSolver->fieldAt(y.pos);
        acc_total += E_sc * y.ion_charge_C / y.mass_kg;
    }

    
    // --- Correct velocity for background gas flow ---
    Vec3 vel_total = y.vel + y.domain_gas_velocity_m_s;
    
    // --- Return updated state ---
    return {vel_total, acc_total};
}

}  // namespace physics
}  // namespace ICARION
