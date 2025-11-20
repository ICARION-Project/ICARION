// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        TransferDomainData.h
 *   @brief       Transfers instrument domain and global parameters to GPU-compatible structs.
 *
 * @details
 * Provides functions to convert CPU-side parameter structures into GPU-compatible versions.
 *
 *
 *   @date        2025-10-17
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */
#pragma once
#include "paramUtils.h"
#include "GPUParams.h"
#include "InstrumentEnums.h"

/**
 * @brief Convert a CPU InstrumentDomain into a GPU-compatible DomainGPU.
 */
inline DomainGPU fromCPU(const InstrumentDomain& dom)
{
    DomainGPU d{};

    // --- Instrument type ---
    d.instrument = static_cast<InstrumentGPU>(dom.instrument);

    // --- RF parameters ---
    d.RF.voltage_V        = dom.RF.voltage_V;
    d.RF.omega_rad_s      = dom.RF.angular_frequency_rad_s;
    d.RF.phase_rad        = dom.RF.phase_rad;

    // --- DC parameters ---
    d.DC.axial_V          = dom.DC.axial_V;
    d.DC.quad_V           = dom.DC.quad_V;
    d.DC.radial_V         = dom.DC.radial_V;
    d.DC.enable_radial_voltage_sweep = dom.DC.enable_radial_voltage_sweep ? 1 : 0;
    d.DC.radial_start_time_s = dom.DC.radial_start_time_s;
    d.DC.radial_rise_time_s  = dom.DC.radial_rise_time_s;
    d.DC.radial_slope_V_s    = dom.DC.radial_slope_V_s;

    // --- AC parameters ---
    d.AC.voltage_V        = dom.AC.voltage_V;
    d.AC.omega_rad_s      = dom.AC.angular_frequency_rad_s;
    d.AC.enable_voltage_sweep = dom.AC.enable_voltage_sweep ? 1 : 0;
    d.AC.enable_frequency_sweep = dom.AC.enable_frequency_sweep ? 1 : 0;
    d.AC.ac_start_freq_Hz = dom.AC.ac_start_freq_Hz;
    d.AC.ac_sweep_slope_Hz_per_s = dom.AC.ac_sweep_slope_Hz_per_s;
    d.AC.start_time_s     = dom.AC.start_time_s;
    d.AC.rise_time_s      = dom.AC.rise_time_s;
    d.AC.amplitude_slope_V_s = dom.AC.amplitude_slope_V_s;
    // LQIT locking parameters (ICARION)
    d.AC.lqit_lock_enable = dom.AC.lqit_lock_enable ? 1 : 0;
    d.AC.lqit_lock_phase_rad = dom.AC.lqit_lock_phase_rad;
    d.AC.lqit_lock_bandwidth_Hz = dom.AC.lqit_lock_bandwidth_Hz;

    // --- Geometry ---
    d.geom.radius_m       = dom.geom.radius_m;
    d.geom.length_m       = dom.geom.length_m;
    d.geom.radius_char_m  = dom.geom.radius_char_m;
    d.geom.radius_out_m   = dom.geom.radius_out_m;
    d.geom.radius_in_m    = dom.geom.radius_in_m;
    d.geom.acc_length_m   = dom.geom.acc_length_m;
    d.geom.end_aperture_m = dom.geom.end_aperture_m;
    d.geom.origin_m = dom.geom.origin_m;
    // Copy precomputed rotation (global -> local) matrix rows into GPU struct
    d.geom.rot_row0.x = dom.rotation_global_to_local.m[0][0];
    d.geom.rot_row0.y = dom.rotation_global_to_local.m[0][1];
    d.geom.rot_row0.z = dom.rotation_global_to_local.m[0][2];

    d.geom.rot_row1.x = dom.rotation_global_to_local.m[1][0];
    d.geom.rot_row1.y = dom.rotation_global_to_local.m[1][1];
    d.geom.rot_row1.z = dom.rotation_global_to_local.m[1][2];

    d.geom.rot_row2.x = dom.rotation_global_to_local.m[2][0];
    d.geom.rot_row2.y = dom.rotation_global_to_local.m[2][1];
    d.geom.rot_row2.z = dom.rotation_global_to_local.m[2][2];

    // --- Magnetic field ---
    d.B.enabled = dom.B.enabled ? 1 : 0;
    d.B.Bxyz = dom.B.field_strength_T;       // map CPU field strength -> GPU Bxyz
    d.B.B_gradient = dom.B.field_gradient_T_m; // map CPU gradient -> GPU gradient (T/m)

    // --- Environment parameters (needed for collision models) ---
    d.env.pressure_Pa = dom.env.pressure_Pa;
    d.env.temperature_K = dom.env.temperature_K;
    d.env.particle_density_m_3 = dom.env.particle_density_m_3;
    d.env.mean_thermal_velocity_m_s = dom.env.mean_thermal_velocity_m_s;
    d.env.neutral_mass_kg = dom.env.neutral_mass_kg;
    d.env.neutral_polarizability_m3 = dom.env.neutral_polarizability_m3;
    d.env.neutral_radius_m = dom.env.neutral_radius_m;
    d.env.flow_velocity_m_s = dom.env.gas_velocity_m_s;

    // --- FieldArray ---
    d.FA.loaded = dom.fieldArrayLoaded ? 1 : 0;
    // TODO: später: Pointer, dimensions, strides, texture object

    return d;
}

/**
 * @brief Convert CPU global simulation parameters to GPU version.
 */
inline GlobalParamsGPU fromCPU(const GlobalParams& g)
{
    GlobalParamsGPU gGPU{};
    // store collision model as integer on the GPU side
    gGPU.collisionModel = static_cast<int>(g.collisionModel);
    // Add more global flags here when needed
    return gGPU;
}

/**
 * @brief Convert GPU DomainGPU back to a CPU InstrumentDomain.
 */
inline InstrumentDomain toCPU(const DomainGPU& src)
{
    InstrumentDomain dom{};

    // --- Instrument type ---
    dom.instrument = static_cast<Instrument>(src.instrument);

    // --- RF parameters ---
    dom.RF.voltage_V               = src.RF.voltage_V;
    dom.RF.angular_frequency_rad_s = src.RF.omega_rad_s;
    dom.RF.phase_rad               = src.RF.phase_rad;

    // --- DC parameters ---
    dom.DC.axial_V                 = src.DC.axial_V;
    dom.DC.quad_V                  = src.DC.quad_V;
    dom.DC.radial_V                = src.DC.radial_V;
    dom.DC.enable_radial_voltage_sweep = (src.DC.enable_radial_voltage_sweep != 0);
    dom.DC.radial_start_time_s     = src.DC.radial_start_time_s;
    dom.DC.radial_rise_time_s      = src.DC.radial_rise_time_s;
    dom.DC.radial_slope_V_s        = src.DC.radial_slope_V_s;

    // --- AC parameters ---
    dom.AC.voltage_V               = src.AC.voltage_V;
    dom.AC.angular_frequency_rad_s = src.AC.omega_rad_s;
    dom.AC.enable_voltage_sweep    = (src.AC.enable_voltage_sweep != 0);
    dom.AC.enable_frequency_sweep  = (src.AC.enable_frequency_sweep != 0);
    dom.AC.ac_start_freq_Hz        = src.AC.ac_start_freq_Hz;
    dom.AC.ac_sweep_slope_Hz_per_s = src.AC.ac_sweep_slope_Hz_per_s;
    dom.AC.start_time_s            = src.AC.start_time_s;
    dom.AC.rise_time_s             = src.AC.rise_time_s;
    dom.AC.amplitude_slope_V_s     = src.AC.amplitude_slope_V_s;
    // LQIT locking parameters
    dom.AC.lqit_lock_enable = (src.AC.lqit_lock_enable != 0);
    dom.AC.lqit_lock_phase_rad = src.AC.lqit_lock_phase_rad;
    dom.AC.lqit_lock_bandwidth_Hz = src.AC.lqit_lock_bandwidth_Hz;

    // --- Geometry ---
    dom.geom.radius_m        = src.geom.radius_m;
    dom.geom.length_m        = src.geom.length_m;
    dom.geom.radius_char_m   = src.geom.radius_char_m;
    dom.geom.radius_out_m    = src.geom.radius_out_m;
    dom.geom.radius_in_m     = src.geom.radius_in_m;
    dom.geom.acc_length_m    = src.geom.acc_length_m;
    dom.geom.end_aperture_m  = src.geom.end_aperture_m;

    // --- Magnetic field ---
    dom.B.enabled        = (src.B.enabled != 0);
    // map GPU B back to CPU MagneticField members
    dom.B.field_strength_T   = src.B.Bxyz;
    dom.B.field_gradient_T_m = src.B.B_gradient;

    // --- FieldArray flags ---
    dom.fieldArrayLoaded = (src.FA.loaded != 0);

    return dom;
}

/**
 * @brief Convert GPU GlobalParamsGPU back to CPU GlobalParams.
 */
inline GlobalParams toCPU(const GlobalParamsGPU& src)
{
    GlobalParams g{};
    // convert integer back to CollisionModel enum
    g.collisionModel = static_cast<CollisionModel>(src.collisionModel);
    // add other globals if they exist
    return g;
}