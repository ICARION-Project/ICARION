// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <string>
#include <string_view>

namespace ICARION {
namespace instrument {

/**
 * @brief Supported instrument types in ICARION simulations
 * 
 * Canonical enumeration used across config/physics/integrator and mirrored
 * by GPU enums (see param/InstrumentEnums.h).
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
