// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#ifndef ICARION_CONFIG_ENUM_MAPPER_H
#define ICARION_CONFIG_ENUM_MAPPER_H

#include "../types/PhysicsEnums.h"
#include "../types/SolverEnums.h"
#include "core/config/types/InstrumentTypes.h"
#include <string>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>

namespace ICARION::config {

/**
 * @brief String to enum mapping utilities
 * 
 * Handles case-insensitive parsing with aliases for backward compatibility.
 */
class EnumMapper {
public:
    /**
     * @brief Parse instrument type from string
     * 
     * @param str Instrument name (case-insensitive)
     * @return InstrumentType Enum value
     * @throws std::runtime_error if unknown instrument
     */
    static ICARION::instrument::InstrumentType parse_instrument(const std::string& str) {
        using Instrument = ICARION::instrument::InstrumentType;
        
        static const std::unordered_map<std::string, Instrument> map = {
            {"lqit", Instrument::LQIT},
            {"ims", Instrument::IMS},
            {"sifdt-ms", Instrument::IMS},
            {"sifdt_ms", Instrument::IMS},
            {"sifdt", Instrument::IMS},
            {"tof", Instrument::TOF},
            {"quadrupole", Instrument::QuadrupoleRF},
            {"quadrupolerf", Instrument::QuadrupoleRF},
            {"quad", Instrument::QuadrupoleRF},
            {"orbitrap", Instrument::Orbitrap},
            {"ft-icr", Instrument::FTICR},
            {"ft_icr", Instrument::FTICR},
            {"fticr", Instrument::FTICR},
            {"icr", Instrument::FTICR},
            {"nofixedinstrument", Instrument::NoFixedInstrument},
            {"custom", Instrument::NoFixedInstrument}
        };
        
        std::string lower = to_lower(str);
        auto it = map.find(lower);
        if (it == map.end()) {
            throw std::runtime_error("Unknown instrument type: '" + str + "'");
        }
        return it->second;
    }
    
    /**
     * @brief Parse collision model from string
     * 
     * @param str Collision model name (case-insensitive)
     * @return CollisionModel Enum value
     * @throws std::runtime_error if unknown model
     */
    static CollisionModel parse_collision_model(const std::string& str) {
        static const std::unordered_map<std::string, CollisionModel> map = {
            {"none", CollisionModel::NoCollisions},
            {"nocollisions", CollisionModel::NoCollisions},
            {"vacuum", CollisionModel::NoCollisions},
            {"ehss", CollisionModel::EHSS},
            {"hss", CollisionModel::HSS},
            {"hsmc", CollisionModel::HSS},  // Legacy alias
            {"langevin", CollisionModel::Langevin},
            {"friction", CollisionModel::Friction},
            {"hsd", CollisionModel::HSD},
            {"hardsphere", CollisionModel::HSD},  // Legacy alias
            {"hard-sphere", CollisionModel::HSD},  // Legacy alias
            {"hs", CollisionModel::HSD}  // Legacy alias
        };
        
        std::string lower = to_lower(str);
        auto it = map.find(lower);
        if (it == map.end()) {
            throw std::runtime_error("Unknown collision model: '" + str + "'");
        }
        return it->second;
    }
    
    /**
     * @brief Parse solver type from string
     * 
     * @param str Solver name (case-insensitive)
     * @return SolverType Enum value
     * @throws std::runtime_error if unknown solver
     */
    static SolverType parse_solver(const std::string& str) {
        static const std::unordered_map<std::string, SolverType> map = {
            {"rk4", SolverType::RK4},
            {"rk45", SolverType::RK45},
            {"rkf45", SolverType::RK45},
            {"runge-kutta-4", SolverType::RK4},
            {"runge-kutta-fehlberg", SolverType::RK45},
            {"boris", SolverType::Boris}
        };
        
        std::string lower = to_lower(str);
        auto it = map.find(lower);
        if (it == map.end()) {
            throw std::runtime_error("Unknown solver type: '" + str + "'");
        }
        return it->second;
    }
    
    /**
     * @brief Convert instrument enum to string (for output/logging)
     * 
     * @param instrument Instrument enum value
     * @return String representation
     */
    static std::string instrument_to_string(ICARION::instrument::InstrumentType instrument) {
        using Instrument = ICARION::instrument::InstrumentType;
        
        switch (instrument) {
            case Instrument::LQIT: return "LQIT";
            case Instrument::IMS: return "IMS";
            case Instrument::TOF: return "TOF";
            case Instrument::QuadrupoleRF: return "Quadrupole";
            case Instrument::Orbitrap: return "Orbitrap";
            case Instrument::FTICR: return "FT-ICR";
            case Instrument::NoFixedInstrument: return "Custom";
            case Instrument::UnknownInstrument: return "Unknown";
            default: return "Unknown";
        }
    }
    
    /**
     * @brief Convert collision model enum to string
     * 
     * @param model Collision model enum value
     * @return String representation
     */
    static std::string collision_model_to_string(CollisionModel model) {
        switch (model) {
            case CollisionModel::NoCollisions: return "NoCollisions";
            case CollisionModel::HSD: return "HSD";
            case CollisionModel::Langevin: return "Langevin";
            case CollisionModel::Friction: return "Friction";
            case CollisionModel::EHSS: return "EHSS";
            case CollisionModel::HSS: return "HSS";
            case CollisionModel::UnknownCollisionModel: return "Unknown";
            default: return "Unknown";
        }
    }
    
    /**
     * @brief Convert solver enum to string
     * 
     * @param solver Solver enum value
     * @return String representation
     */
    static std::string solver_to_string(SolverType solver) {
        switch (solver) {
            case SolverType::RK4: return "RK4";
            case SolverType::RK45: return "RK45";
            case SolverType::Boris: return "Boris";
            default: return "Unknown";
        }
    }
    
private:
    /**
     * @brief Convert string to lowercase (helper)
     * 
     * @param str Input string
     * @return Lowercase string
     */
    static std::string to_lower(std::string str) {
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        return str;
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_ENUM_MAPPER_H
