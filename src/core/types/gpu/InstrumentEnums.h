// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file InstrumentEnums.h
 * @brief GPU-compatible instrument type enum (SSOT compliant)
 * 
 * @details
 * Defines InstrumentGPU enum for GPU kernels. This is synchronized with
 * the canonical CPU-side InstrumentType enum from core/config/types/InstrumentTypes.h
 * 
 * SSOT: Single source of truth is InstrumentType in InstrumentTypes.h
 * This header lives in param/ for CUDA-friendly plain enums; kept in sync via static_assert.
 * 
 * @version 2.0 - SSOT Migration (2025-11-23)
 * @author Christoph Schaefer
 */

#pragma once
#include "instrument/InstrumentTypes.h"

/**
 * @brief GPU-compatible instrument type enum
 * 
 * Plain enum (not enum class) for CUDA kernel compatibility.
 * Values must match ICARION::instrument::InstrumentType exactly.
 * 
 * SSOT: InstrumentType (instrument/InstrumentTypes.h) is the single source of truth.
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
    TIMS = 7,              ///< Trapped Ion Mobility Spectrometry
    UnknownInstrument = 8  ///< Unrecognized instrument type
};

// ============================================================================
// SSOT Validation: Ensure GPU enum matches CPU InstrumentType
// ============================================================================
// These static_asserts enforce synchronization with the canonical InstrumentType.
// If these fail, someone changed InstrumentType without updating InstrumentGPU.

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
static_assert(static_cast<int>(InstrumentType::TIMS) == static_cast<int>(InstrumentGPU::TIMS),
    "SSOT violation: InstrumentGPU::TIMS must match InstrumentType::TIMS");
static_assert(static_cast<int>(InstrumentType::UnknownInstrument) == static_cast<int>(InstrumentGPU::UnknownInstrument), 
    "SSOT violation: InstrumentGPU::UnknownInstrument must match InstrumentType::UnknownInstrument");
