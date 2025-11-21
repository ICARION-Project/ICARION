// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "CollisionHandlerFactory.h"
#include "core/io/logger.h"
#include <sstream>

namespace ICARION::physics {

std::unique_ptr<ICollisionHandler> CollisionHandlerFactory::create(
    const config::PhysicsConfig& config,
    const GeometryMap* geometry_map,
    double gamma_for_ou,
    bool enable_logging
) {
    using config::CollisionModel;
    
    // Handle stochastic collision models
    switch (config.collision_model) {
        
        case CollisionModel::EHSS: {
            // EHSS requires geometry map
            if (!geometry_map) {
                throw std::invalid_argument(
                    "[CollisionHandlerFactory] EHSS collision model requires "
                    "molecular geometry map, but geometry_map is nullptr"
                );
            }
            
            if (enable_logging) {
                ICARION::io::debug_log(
                    "[CollisionHandlerFactory] Creating EHSSCollisionHandler "
                    "with geometry map"
                );
            }
            
            return std::make_unique<EHSSCollisionHandler>(
                *geometry_map,
                enable_logging
            );
        }
        
        case CollisionModel::HSS: {
            if (enable_logging) {
                ICARION::io::debug_log(
                    "[CollisionHandlerFactory] Creating HSSCollisionHandler "
                    "(isotropic scattering)"
                );
            }
            
            return std::make_unique<HSSCollisionHandler>(enable_logging);
        }
        
        // Deterministic models - use DampingForce for damping
        case CollisionModel::Friction:
        case CollisionModel::Langevin:
        case CollisionModel::HSD: {
            // Check if OU thermalization should be added
            if (config.enable_ou_thermalization) {
                if (gamma_for_ou <= 0.0) {
                    throw std::invalid_argument(
                        "[CollisionHandlerFactory] OU thermalization enabled but "
                        "gamma_for_ou <= 0. Must provide valid damping coefficient."
                    );
                }
                
                if (enable_logging) {
                    std::ostringstream msg;
                    msg << "[CollisionHandlerFactory] Creating OUCollisionHandler "
                        << "for deterministic model (gamma=" << gamma_for_ou << ")";
                    ICARION::io::debug_log(msg.str());
                }
                
                return std::make_unique<OUCollisionHandler>(gamma_for_ou);
            }
            
            // No OU thermalization - deterministic damping only via DampingForce
            if (enable_logging) {
                ICARION::io::debug_log(
                    "[CollisionHandlerFactory] Deterministic model without OU - "
                    "returning nullptr (use DampingForce only)"
                );
            }
            return nullptr;
        }
        
        case CollisionModel::NoCollisions: {
            if (enable_logging) {
                ICARION::io::debug_log(
                    "[CollisionHandlerFactory] NoCollisions model - returning nullptr"
                );
            }
            return nullptr;
        }
        
        case CollisionModel::UnknownCollisionModel:
        default: {
            throw std::invalid_argument(
                "[CollisionHandlerFactory] Unknown or unsupported collision model"
            );
        }
    }
}

} // namespace ICARION::physics
