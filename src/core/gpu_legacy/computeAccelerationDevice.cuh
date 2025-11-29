
#pragma once
#include "core/types/Vec3.h"
#include "core/param/GPUParams.h"
#include "core/physics/fields/physics_fields_shared.h"
#include "core/physics/fields/physics_math_shared.h"
#include "core/types/IonState_GPU.h"   // GPU-Variante von IonState
#include "core/types/Vec3.h"
#include <cuda_runtime.h>

// =============================================================
// Simple field helpers
// =============================================================

__device__ inline Vec3 DCField_device(const Vec3& pos, double axial_V, double L) {
    if (L <= 0.0) return {0.0, 0.0, 0.0};
    return {0.0, 0.0, axial_V / L}; // uniform axial field
}

// Helper: Transform global position to local domain coordinates
__device__ inline Vec3 global_to_local_device(const Vec3& global_pos, const Vec3& origin,
                                              const Vec3& rot_row0, const Vec3& rot_row1, const Vec3& rot_row2) {
    Vec3 diff;
    diff.x = global_pos.x - origin.x;
    diff.y = global_pos.y - origin.y;
    diff.z = global_pos.z - origin.z;
    
    Vec3 local_pos;
    local_pos.x = rot_row0.x * diff.x + rot_row0.y * diff.y + rot_row0.z * diff.z;
    local_pos.y = rot_row1.x * diff.x + rot_row1.y * diff.y + rot_row1.z * diff.z;
    local_pos.z = rot_row2.x * diff.x + rot_row2.y * diff.y + rot_row2.z * diff.z;
    
    return local_pos;
}

__device__ inline Vec3 RFField_device(const Vec3& pos, double Vrf, double Vquad,
                                      double omega, double r0, double t) {
    if (r0 <= 0.0) return {0.0, 0.0, 0.0};
    // Match CPU implementation: Ex = 2*x*(DC + RF*cos(omega*t)) / r0^2
    double c = cos_shared(omega * t);
    double fac = 2.0 * (Vquad + Vrf * c) / (r0 * r0);
    return { fac * pos.x, -fac * pos.y, 0.0 };
}

__device__ inline Vec3 ACField_device(const Vec3& pos, double Vac, double omega,
                                      double radius, double t, Vec3 dir) {
    // Match CPU implementation in defineFields.cpp:
    // dir is normalized; compute magnitude = -1/radius * (voltage * cos(omega*t))
    double len = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
    Vec3 dir_unit = (len > 0.0) ? Vec3{dir.x/len, dir.y/len, dir.z/len} : Vec3{0.0,0.0,0.0};
    double c = cos_shared(omega * t);
    double mag = 0.0;
    if (radius != 0.0) mag = -1.0 / radius * (Vac * c);
    return { dir_unit.x * mag, dir_unit.y * mag, dir_unit.z * mag };
}

// Orbitrap field matching CPU implementation in defineFields.cpp
// k is the force constant, r_char is characteristic radius, length_m is trap length
__device__ inline Vec3 OrbitrapField_device(const Vec3& pos, double k, double r_char, double length_m) {
    double r2 = pos.x*pos.x + pos.y*pos.y;
    if (r2 < 1e-18) r2 = 1e-18;
    double C = 1.0 - (r_char * r_char) / r2;
    double z_center = pos.z - 0.5 * length_m;

    double Ex = 0.5 * k * pos.x * C;
    double Ey = 0.5 * k * pos.y * C;
    double Ez = - k * z_center;

    return Vec3{Ex, Ey, Ez};
}

__device__ inline Vec3 MagneticField_device(const Vec3& pos, const Vec3& B_field_T) {
    // Placeholder: uniform magnetic field
    return B_field_T;
}

// simple FTICR device field (mirrors CPU FTICRField)
__device__ inline Vec3 FTICRField_device(const Vec3& pos, double voltage, double characteristic_length, double instrument_length_axial_m) {
    double factor = voltage / (characteristic_length * characteristic_length);
    double Ex = pos.x * factor;
    double Ey = pos.y * factor;
    double Ez = -2.0 * (pos.z - 0.5 * instrument_length_axial_m) * factor;
    return Vec3{Ex, Ey, Ez};
}

// =============================================================
// Core compute_accelerations_device
// =============================================================

__device__ inline void compute_accelerations_device(double t,
                                                    const IonStateGPU& y,
                                                    const GlobalParamsGPU& g,
                                                    const DomainGPU& dom,
                                                    Vec3& out_vel_total,
                                                    Vec3& out_acc_total)
{
    Vec3 acc_total = {0.0, 0.0, 0.0};

    // --- Field arrays (not yet implemented)
    if (dom.FA.loaded) {
        // placeholder: field interpolation
    } else {
        switch (dom.instrument) {

            // =====================================================
            // LQIT
            // =====================================================
            case LQIT: {
                // Transform global position to local coordinates for field calculations
                Vec3 local_pos = global_to_local_device(y.pos, dom.geom.origin_m, 
                                                        dom.geom.rot_row0, dom.geom.rot_row1, dom.geom.rot_row2);
                
                // (1) RF quadrupole
                acc_total += RFField_device(local_pos, dom.RF.voltage_V, dom.DC.quad_V,
                                            dom.RF.omega_rad_s, dom.geom.radius_m, t)
                             * (y.ion_charge_C / y.mass_kg);

                // (2) AC excitation
                double Vac = dom.AC.voltage_V;
                if (dom.AC.enable_voltage_sweep) {
                    double t_rel = t - dom.AC.start_time_s;
                    if (t_rel <= 0.0) Vac = dom.AC.voltage_V;
                    else if (t_rel < dom.AC.rise_time_s)
                        Vac = dom.AC.voltage_V + t_rel * dom.AC.amplitude_slope_V_s;
                    else
                        Vac = dom.AC.voltage_V + dom.AC.amplitude_slope_V_s * dom.AC.rise_time_s;
                }
                // Support frequency sweep on device (linear in Hz) to match CPU
                // Use shared helper to ensure identical formula on host and device.
                double omega_ac = dom.AC.omega_rad_s;
                if (dom.AC.enable_frequency_sweep) {
                    double t_rel = t - dom.AC.start_time_s;
                    omega_ac = compute_ac_omega_from_hz(dom.AC.ac_start_freq_Hz,
                                                        dom.AC.ac_sweep_slope_Hz_per_s,
                                                        t_rel,
                                                        dom.AC.enable_frequency_sweep);
                }

                // LQIT-lock semantics on device:
                // - If enabled, apply a phase offset to the AC term (phase_rad)
                //   by shifting the effective time: sin(omega*t + phase) == sin(omega*(t + phase/omega)).
                // - Optionally gate amplitude using a Gaussian around the RF frequency when
                //   lqit_lock_bandwidth_Hz > 0. This gives a smooth attenuation when the
                //   instantaneous AC frequency drifts away from the RF lock frequency.
                double t_ac = t;
                if (dom.AC.lqit_lock_enable) {
                    double phase = dom.AC.lqit_lock_phase_rad;
                    if (fabs(omega_ac) > 1e-30) t_ac += phase / omega_ac;

                    double bw = dom.AC.lqit_lock_bandwidth_Hz;
                    if (bw > 0.0) {
                        const double TWO_PI = 2.0 * 3.14159265358979323846;
                        double freq_Hz = omega_ac / TWO_PI;
                        double rf_freq = dom.RF.omega_rad_s / TWO_PI; // dom.RF.frequency_Hz may not be present on GPU
                        // bandwidth interpreted as full-width; convert to sigma
                        double sigma = bw * 0.5;
                        double delta = fabs(freq_Hz - rf_freq);
                        double amp_scale = 1.0;
                        if (sigma > 0.0) amp_scale = exp_shared(- (delta*delta) / (2.0 * sigma * sigma));
                        Vac *= amp_scale;
                    }
                }

                acc_total += ACField_device(local_pos, Vac, omega_ac, dom.geom.radius_m, t_ac, {1,0,0})
                             * (y.ion_charge_C / y.mass_kg);

                // (3) Axial DC near boundaries (already using local_pos from above)
                if (local_pos.z > dom.geom.length_m * 0.9) {
                    acc_total += DCField_device(local_pos, -dom.DC.axial_V, dom.geom.length_m * 0.1)
                                 * (y.ion_charge_C / y.mass_kg);
                } else if (local_pos.z < dom.geom.length_m * 0.1) {
                    acc_total += DCField_device(local_pos, dom.DC.axial_V, dom.geom.length_m * 0.1)
                                 * (y.ion_charge_C / y.mass_kg);
                }
            } break;

            // =====================================================
            // IMS (Ion Mobility Spectrometry)
            // =====================================================
            case IMS: {
                // Axial DC
                acc_total += DCField_device(y.pos, dom.DC.axial_V, dom.geom.length_m)
                             * (y.ion_charge_C / y.mass_kg);
                // Optional RF focusing
                acc_total += RFField_device(y.pos, dom.RF.voltage_V, dom.DC.quad_V,
                                            dom.RF.omega_rad_s, dom.geom.radius_m, t)
                             * (y.ion_charge_C / y.mass_kg);
            } break;

            // =====================================================
            // Orbitrap
            // =====================================================
            case Orbitrap: {
                // Map radial DC (radial_V) and DC sweep if present
                double radial_voltage = dom.DC.radial_V;
                if (dom.DC.enable_radial_voltage_sweep) {
                    double t_rel = t - dom.DC.radial_start_time_s;
                    if (t_rel <= 0.0) radial_voltage = dom.DC.radial_V;
                    else if (t_rel < dom.DC.radial_rise_time_s)
                        radial_voltage = dom.DC.radial_V + t_rel * dom.DC.radial_slope_V_s;
                    else
                        radial_voltage = dom.DC.radial_V + dom.DC.radial_slope_V_s * dom.DC.radial_rise_time_s;
                }
                // compute k same as CPU implementation
                double denom = dom.geom.radius_char_m * dom.geom.radius_char_m * log(dom.geom.radius_out_m / dom.geom.radius_in_m)
                               - 0.5*(dom.geom.radius_out_m * dom.geom.radius_out_m - dom.geom.radius_in_m * dom.geom.radius_in_m);
                double k = 0.0;
                if (fabs(denom) > 1e-18) k = 2.0 * radial_voltage / denom;

                acc_total += OrbitrapField_device(y.pos, k, dom.geom.radius_char_m, dom.geom.length_m) * (y.ion_charge_C / y.mass_kg);
            } break;

            // =====================================================
            // QuadrupoleRF (radial RF + axial DC)
            // =====================================================
            case QuadrupoleRF: {
                acc_total += RFField_device(y.pos, dom.RF.voltage_V, dom.DC.quad_V,
                                            dom.RF.omega_rad_s, dom.geom.radius_m, t) * (y.ion_charge_C / y.mass_kg);
                acc_total += DCField_device(y.pos, dom.DC.axial_V, dom.geom.length_m) * (y.ion_charge_C / y.mass_kg);
            } break;

            // =====================================================
            // TOF
            // =====================================================
            case TOF: {
                if (y.pos.z < dom.geom.acc_length_m) {
                    acc_total += DCField_device(y.pos, dom.DC.axial_V, dom.geom.acc_length_m) * (y.ion_charge_C / y.mass_kg);
                }
            } break;

            // =====================================================
            // FT-ICR
            // =====================================================
            case FT_ICR: {
                double d = sqrt((dom.geom.length_m*dom.geom.length_m / 8.0) + (dom.geom.radius_m*dom.geom.radius_m / 4.0));
                acc_total += FTICRField_device(y.pos, dom.DC.radial_V, d, dom.geom.length_m) * (y.ion_charge_C / y.mass_kg);
            } break;

            // =====================================================
            // NoFixedInstrument: apply any non-zero DC/RF/AC
            // =====================================================
            case NoFixedInstrument: {
                if (fabs(dom.DC.axial_V) > 1e-12) {
                    acc_total += DCField_device(y.pos, dom.DC.axial_V, dom.geom.length_m) * (y.ion_charge_C / y.mass_kg);
                }
                if (fabs(dom.RF.voltage_V) > 1e-12) {
                    acc_total += RFField_device(y.pos, dom.RF.voltage_V, dom.DC.quad_V,
                                                dom.RF.omega_rad_s, dom.geom.radius_m, t) * (y.ion_charge_C / y.mass_kg);
                }
                if (fabs(dom.AC.voltage_V) > 1e-12) {
                    // simple AC along x
                    double Vac = dom.AC.voltage_V;
                    if (dom.AC.enable_voltage_sweep) {
                        double t_rel = t - dom.AC.start_time_s;
                        if (t_rel <= 0.0) Vac = dom.AC.voltage_V;
                        else if (t_rel < dom.AC.rise_time_s)
                            Vac = dom.AC.voltage_V + t_rel * dom.AC.amplitude_slope_V_s;
                        else
                            Vac = dom.AC.voltage_V + dom.AC.amplitude_slope_V_s * dom.AC.rise_time_s;
                    }
                    acc_total += ACField_device(y.pos, Vac, dom.AC.omega_rad_s, dom.geom.radius_m, t, {1,0,0}) * (y.ion_charge_C / y.mass_kg);
                }
            } break;

            default:
                acc_total = {0.0, 0.0, 0.0};
                break;
        }
    }

    // --- Collision forces ---
    Vec3 collision_acc = {0.0, 0.0, 0.0};
    
    if (g.collisionModel == 2) {  // CollisionModel::Friction = 2
        // Friction collision: damping based on ion mobility
        // F = -q/K * v, where K is mobility
        const double LOSCHMIDT_CONSTANT = 2.6867811e25;  // m^-3
        double ion_mobility = (y.reduced_mobility_cm2_Vs * 1e-4) * LOSCHMIDT_CONSTANT / y.domain_particle_density_m3;
        double damping = -y.ion_charge_C / ion_mobility;
        Vec3 friction_force = y.vel * damping;
        collision_acc = friction_force / y.mass_kg;
    }
    // TODO: Add Langevin, EHSS, HSS collision models
    
    acc_total += collision_acc;

    // --- Gas flow correction
    Vec3 vel_total = y.vel + y.domain_gas_velocity_m_s;

    out_vel_total = vel_total;
    out_acc_total = acc_total;
}
