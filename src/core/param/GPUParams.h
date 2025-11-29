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
 *   @file        GPUParams.h
 *   @brief       Defines GPU-compatible parameter structures.
 *
 *   @details
 *   Defines simplified versions of parameter structs for transfer to GPU memory.
 *
 *   @date        2025-10-17
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */
#pragma once
#include "core/types/Vec3.h"   
#include <cstdint>
#include "core/config/types/InstrumentTypes.h"

// ============================================================================
// GPU-compatible Instrument Enum (SSOT: InstrumentType is canonical)
// ============================================================================
/**
 * @brief GPU-compatible instrument type enum
 * 
 * Plain enum (not enum class) for CUDA kernel compatibility.
 * Values must match ICARION::instrument::InstrumentType exactly.
 * 
 * SSOT: InstrumentType (core/config/types/InstrumentTypes.h) is the single source of truth.
 * This enum is automatically validated at compile-time via static_assert.
 */
enum InstrumentGPU : int { 
    LQIT = 0,              ///< Linear Quadrupole Ion Trap
    IMS = 1,               ///< Ion Mobility Spectrometry
    Orbitrap = 2,          ///< Orbitrap Mass Analyzer
    QuadrupoleRF = 3,      ///< Quadrupole RF (includes SLIM)
    TOF = 4,               ///< Time-of-Flight
    FT_ICR = 5,            ///< Fourier Transform Ion Cyclotron Resonance
    NoFixedInstrument = 6, ///< Generic/custom instrument
    UnknownInstrument = 7  ///< Unrecognized instrument type
};

// SSOT Validation: Ensure GPU enum matches CPU InstrumentType
using InstrumentType = ICARION::instrument::InstrumentType;

static_assert(static_cast<int>(InstrumentType::LQIT) == static_cast<int>(InstrumentGPU::LQIT), 
    "SSOT violation: InstrumentGPU::LQIT must match InstrumentType::LQIT");
static_assert(static_cast<int>(InstrumentType::IMS) == static_cast<int>(InstrumentGPU::IMS), 
    "SSOT violation: InstrumentGPU::IMS must match InstrumentType::IMS");
static_assert(static_cast<int>(InstrumentType::Orbitrap) == static_cast<int>(InstrumentGPU::Orbitrap), 
    "SSOT violation: InstrumentGPU::Orbitrap must match InstrumentType::Orbitrap");
static_assert(static_cast<int>(InstrumentType::QuadrupoleRF) == static_cast<int>(InstrumentGPU::QuadrupoleRF), 
    "SSOT violation: InstrumentGPU::QuadrupoleRF must match InstrumentType::QuadrupoleRF");
static_assert(static_cast<int>(InstrumentType::TOF) == static_cast<int>(InstrumentGPU::TOF), 
    "SSOT violation: InstrumentGPU::TOF must match InstrumentType::TOF");
static_assert(static_cast<int>(InstrumentType::FTICR) == static_cast<int>(InstrumentGPU::FT_ICR), 
    "SSOT violation: InstrumentGPU::FT_ICR must match InstrumentType::FTICR");
static_assert(static_cast<int>(InstrumentType::NoFixedInstrument) == static_cast<int>(InstrumentGPU::NoFixedInstrument), 
    "SSOT violation: InstrumentGPU::NoFixedInstrument must match InstrumentType::NoFixedInstrument");
static_assert(static_cast<int>(InstrumentType::UnknownInstrument) == static_cast<int>(InstrumentGPU::UnknownInstrument), 
    "SSOT violation: InstrumentGPU::UnknownInstrument must match InstrumentType::UnknownInstrument");

// ============================================================================

/* RF parameters */
struct RFParamsGPU {
    double voltage_V;
    double omega_rad_s;
    double phase_rad; // optional RF phase (radians)
};

/* DC parameters */
struct DCParamsGPU {
    double axial_V;
    double quad_V;
    double radial_V;
    // sweeps
    int enable_radial_voltage_sweep;
    double radial_start_time_s, radial_rise_time_s, radial_slope_V_s;
};

/* AC parameters */
struct ACParamsGPU {
    double voltage_V;
    double omega_rad_s;
    int enable_voltage_sweep;
    // Frequency sweep support (linear in Hz)
    int enable_frequency_sweep;
    double ac_start_freq_Hz;
    double ac_sweep_slope_Hz_per_s;
    double start_time_s, rise_time_s, amplitude_slope_V_s;
    // ICARION / LQIT locking parameters
    int lqit_lock_enable;        // 0/1
    double lqit_lock_phase_rad;  // radians
    double lqit_lock_bandwidth_Hz;
};

/* Geometry */
struct GeomGPU {
    double radius_m;
    double length_m;
    double radius_char_m;
    double radius_out_m;
    double radius_in_m;
    double acc_length_m;
    double end_aperture_m;
    Vec3 origin_m; // global origin of the domain (meters)
    // Rotation matrix rows that map GLOBAL -> LOCAL coordinates on the device.
    // local = R_global_to_local * (global - origin)
    Vec3 rot_row0;
    Vec3 rot_row1;
    Vec3 rot_row2;
};

/* Environmental parameters */
struct EnvParamsGPU {
    double pressure_Pa;
    double temperature_K;
    double particle_density_m_3;
    double mean_thermal_velocity_m_s;
    double neutral_mass_kg;
    double neutral_polarizability_m3;
    double neutral_radius_m;  // Hard-sphere radius for EHSS collisions
    Vec3 flow_velocity_m_s;
};

/* Magnetic field parameters */
struct BParamsGPU {
    int enabled;
    Vec3 Bxyz;
    Vec3 B_gradient;
};

/* Optional field array structure */
struct FieldArrayGPU {
    int loaded;
    /* + dims, strides, pointer/tex */
};

/* Domain parameters */
struct DomainGPU {
    InstrumentGPU instrument;
    RFParamsGPU RF;
    DCParamsGPU DC;
    ACParamsGPU AC;
    GeomGPU     geom;
    BParamsGPU  B;
    EnvParamsGPU env;

    FieldArrayGPU FA; // vorerst FA.loaded=0 -> ignorieren
};

struct GlobalParamsGPU {
    int collisionModel;   
    // ... weitere globale Schalter
};

// ------------------------------------------------------------------
// ABI Contract: Alignment & Size Assertions
// Ensure GPU parameter structs have consistent layout between host and device.
// All structs should be at least 8-byte aligned (for double fields).
// Prefer 16-byte alignment for structs with Vec3 or for better cache/memory performance.
// ------------------------------------------------------------------
static_assert(sizeof(double) == 8, "expected 8-byte double");
static_assert(sizeof(int) == 4, "expected 4-byte int");
static_assert(sizeof(Vec3) == 24, "expected 24-byte Vec3 (3 doubles)");

// Individual parameter struct checks
static_assert(alignof(RFParamsGPU) >= 8, "RFParamsGPU alignment < 8 bytes");
static_assert(sizeof(RFParamsGPU) % 8 == 0, "RFParamsGPU size not multiple of 8 bytes");

static_assert(alignof(DCParamsGPU) >= 8, "DCParamsGPU alignment < 8 bytes");
static_assert(sizeof(DCParamsGPU) % 8 == 0, "DCParamsGPU size not multiple of 8 bytes");

static_assert(alignof(ACParamsGPU) >= 8, "ACParamsGPU alignment < 8 bytes");
static_assert(sizeof(ACParamsGPU) % 8 == 0, "ACParamsGPU size not multiple of 8 bytes");

static_assert(alignof(GeomGPU) >= 8, "GeomGPU alignment < 8 bytes");
static_assert(sizeof(GeomGPU) % 8 == 0, "GeomGPU size not multiple of 8 bytes");

static_assert(alignof(EnvParamsGPU) >= 8, "EnvParamsGPU alignment < 8 bytes");
static_assert(sizeof(EnvParamsGPU) % 8 == 0, "EnvParamsGPU size not multiple of 8 bytes");

static_assert(alignof(BParamsGPU) >= 8, "BParamsGPU alignment < 8 bytes");
static_assert(sizeof(BParamsGPU) % 8 == 0, "BParamsGPU size not multiple of 8 bytes");

// Composite struct checks
static_assert(alignof(DomainGPU) >= 8, "DomainGPU alignment < 8 bytes");
static_assert(sizeof(DomainGPU) % 8 == 0, "DomainGPU size not multiple of 8 bytes; check struct packing/padding");

static_assert(alignof(GlobalParamsGPU) >= 4, "GlobalParamsGPU alignment < 4 bytes");

// Size documentation (helpful for debugging transfer sizes)
// These will show up in compiler output if -Wno-error is not set
// Uncomment for development/debugging:
// static_assert(sizeof(ACParamsGPU) == 88, "ACParamsGPU size changed - update transfer code if needed");
// static_assert(sizeof(DomainGPU) == 448, "DomainGPU size changed - update transfer code if needed");
