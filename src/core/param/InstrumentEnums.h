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
 *   @file        InstrumentEnums.h
 *   @brief       Defines enumeration for supported instrument types.
 *
 *   @details
 *   Defines the `InstrumentGPU` enum listing supported instrument types
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
// Keep ordering identical to CPU-side `enum class Instrument` in paramUtils.h
// CPU order: LQIT=0, SIFDT_MS=1, IMS=2, Orbitrap=3, QuadrupoleRF=4, TOF=5, FT_ICR=6, NoFixedInstrument=7
enum InstrumentGPU : int { LQIT=0, SIFDT_MS=1, IMS=2, Orbitrap=3, QuadrupoleRF=4, TOF=5, FT_ICR=6, NoFixedInstrument=7 };

// Safety checks: ensure CPU-side `enum class Instrument` ordering remains identical
// to the GPU-side `InstrumentGPU` enum. If these static_asserts fail the build,
// someone has reordered the CPU enum without updating the GPU mapping.
static_assert(static_cast<int>(Instrument::LQIT) == static_cast<int>(InstrumentGPU::LQIT), "Instrument ordering mismatch: LQIT");
static_assert(static_cast<int>(Instrument::SIFDT_MS) == static_cast<int>(InstrumentGPU::SIFDT_MS), "Instrument ordering mismatch: SIFDT_MS");
static_assert(static_cast<int>(Instrument::IMS) == static_cast<int>(InstrumentGPU::IMS), "Instrument ordering mismatch: IMS");
static_assert(static_cast<int>(Instrument::Orbitrap) == static_cast<int>(InstrumentGPU::Orbitrap), "Instrument ordering mismatch: Orbitrap");
static_assert(static_cast<int>(Instrument::QuadrupoleRF) == static_cast<int>(InstrumentGPU::QuadrupoleRF), "Instrument ordering mismatch: QuadrupoleRF");
static_assert(static_cast<int>(Instrument::TOF) == static_cast<int>(InstrumentGPU::TOF), "Instrument ordering mismatch: TOF");
static_assert(static_cast<int>(Instrument::FT_ICR) == static_cast<int>(InstrumentGPU::FT_ICR), "Instrument ordering mismatch: FT_ICR");
static_assert(static_cast<int>(Instrument::NoFixedInstrument) == static_cast<int>(InstrumentGPU::NoFixedInstrument), "Instrument ordering mismatch: NoFixedInstrument");

