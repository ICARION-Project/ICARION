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
 *   @file       InstrumentTypes.cpp
 *   @brief      Implementation of instrument types and functions
 *
 *   @details
 *   Implements functions for mapping between instrument type enums and their string representations.
 *
 *   @date       2025-11-14
 *   @version    1.0.0
 *   @authors    ICARION Development Team
 *
 * =====================================================================
 */
 
#include "core/config/types/InstrumentTypes.h"

#include <algorithm>
#include <cctype>

namespace ICARION {
namespace instrument {

namespace {

std::string toUpper(std::string_view input) {
    std::string normalized(input);
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return normalized;
}

}  // namespace

InstrumentType instrumentTypeFromString(std::string_view type) {
    const std::string key = toUpper(type);
    if (key == "LQIT") return InstrumentType::LQIT;
    if (key == "IMS" || key == "SIFDT-MS" || key == "SIFDT_MS" || key == "SIFDT") return InstrumentType::IMS;
    if (key == "ORBITRAP") return InstrumentType::Orbitrap;
    if (key == "QUADRUPOLE" || key == "QUADRUPOLE_RF" || key == "SLIM") return InstrumentType::QuadrupoleRF;
    if (key == "TOF" || key == "TIME-OF-FLIGHT") return InstrumentType::TOF;
    if (key == "FT-ICR" || key == "FT_ICR" || key == "FTICR") return InstrumentType::FTICR;
    if (key == "NONE" || key == "NOFIXEDINSTRUMENT") return InstrumentType::NoFixedInstrument;
    if (key.empty()) return InstrumentType::UnknownInstrument;
    return InstrumentType::UnknownInstrument;
}

std::string instrumentTypeToString(InstrumentType type) {
    switch (type) {
        case InstrumentType::LQIT: return "LQIT";
        case InstrumentType::IMS: return "IMS";
        case InstrumentType::Orbitrap: return "Orbitrap";
        case InstrumentType::QuadrupoleRF: return "QuadrupoleRF";
        case InstrumentType::TOF: return "TOF";
        case InstrumentType::FTICR: return "FTICR";
        case InstrumentType::NoFixedInstrument: return "NoFixedInstrument";
        case InstrumentType::UnknownInstrument:
        default:
            return "UnknownInstrument";
    }
}

bool instrumentSupportsMultipleDomains(InstrumentType type) {
    switch (type) {
        case InstrumentType::IMS:
        case InstrumentType::QuadrupoleRF:
            return true;
        default:
            return false;
    }
}

}  // namespace instrument
}  // namespace ICARION
