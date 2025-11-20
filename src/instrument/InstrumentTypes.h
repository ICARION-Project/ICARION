// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   ------------------------------------------------
 *   Modular framework for simulating ion trajectories in custom
 *   electric fields and background gas environments.
 *
 *   @file       InstrumentTypes.h
 *   @brief      Instrument type enumeration and conversion functions
 *
 *   @details
 *   Defines the canonical InstrumentType enum used throughout ICARION,
 *   along with utilities for string conversion and multi-domain support queries.
 *
 *   @date       2025-11-20
 *   @version    1.0.0
 *   @authors    ICARION Development Team
 *
 * =====================================================================
 */

#pragma once

#include <string>
#include <string_view>

namespace ICARION {
namespace instrument {

/**
 * @brief Supported instrument types in ICARION simulations
 * 
 * This is the canonical enumeration for all instrument types.
 * Note: SIFDT_MS has been consolidated into IMS (Ion Mobility Spectrometry).
 */
enum class InstrumentType {
    LQIT = 0,              ///< Linear Quadrupole Ion Trap
    IMS = 1,               ///< Ion Mobility Spectrometry (includes SIFDT-MS)
    Orbitrap = 2,          ///< Orbitrap Mass Analyzer
    QuadrupoleRF = 3,      ///< Quadrupole RF (includes SLIM)
    TOF = 4,               ///< Time-of-Flight
    FTICR = 5,             ///< Fourier Transform Ion Cyclotron Resonance
    NoFixedInstrument = 6, ///< Generic/custom instrument
    UnknownInstrument = 7  ///< Unrecognized instrument type
};

/**
 * @brief Convert string to InstrumentType
 * 
 * Accepts various naming conventions (case-insensitive):
 * - "IMS", "SIFDT-MS", "SIFDT_MS", "SIFDT" → InstrumentType::IMS
 * - "LQIT" → InstrumentType::LQIT
 * - "Orbitrap" → InstrumentType::Orbitrap
 * - etc.
 * 
 * @param type String representation of instrument type
 * @return Corresponding InstrumentType enum value
 */
InstrumentType instrumentTypeFromString(std::string_view type);

/**
 * @brief Convert InstrumentType to canonical string representation
 * 
 * @param type Instrument type enum
 * @return String representation (e.g., "IMS", "LQIT", "Orbitrap")
 */
std::string instrumentTypeToString(InstrumentType type);

/**
 * @brief Check if instrument type supports multiple simulation domains
 * 
 * @param type Instrument type enum
 * @return true if multi-domain simulations are supported
 */
bool instrumentSupportsMultipleDomains(InstrumentType type);

}  // namespace instrument
}  // namespace ICARION
