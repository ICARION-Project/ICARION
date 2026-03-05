// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "CollisionHandlerFactory.h"
#include "core/log/Logger.h"
#include <sstream>

namespace ICARION::physics {

std::unique_ptr<ICollisionHandler> CollisionHandlerFactory::create(
    const config::PhysicsConfig& config,
    const GeometryMap* geometry_map,
    double gamma_for_ou,
    bool enable_logging,
    const config::SpeciesDatabase* species_db,
    bool enable_gpu,
    unsigned long long gpu_seed,
    size_t gpu_threshold) {
    using config::CollisionModel;

#ifndef ICARION_USE_GPU
    (void)enable_gpu;
    (void)gpu_seed;
    (void)gpu_threshold;
#endif

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
                ICARION::log::debug_log(
                    "[CollisionHandlerFactory] Creating EHSSCollisionHandler "
                    "with geometry map"
                );
            }
            
            GeometryMap geometry_copy = *geometry_map;
            auto cpu_handler = std::make_unique<EHSSCollisionHandler>(
                std::move(geometry_copy),
                enable_logging,
                species_db
            );

#ifdef ICARION_USE_GPU
            if (enable_gpu) {
                return std::make_unique<GPUCollisionHandler>(
                    std::move(cpu_handler),
                    "EHSS",
                    geometry_map,
                    enable_logging,
                    species_db,
                    gpu_seed,
                    gpu_threshold
                );
            }
#endif
            return cpu_handler;
        }
        
        case CollisionModel::HSS: {
            if (enable_logging) {
                ICARION::log::debug_log(
                    "[CollisionHandlerFactory] Creating HSSCollisionHandler "
                    "(mixture-aware HSS)"
                );
            }
            
            auto cpu_handler = std::make_unique<HSSCollisionHandler>(enable_logging, species_db);
#ifdef ICARION_USE_GPU
            if (enable_gpu) {
                return std::make_unique<GPUCollisionHandler>(
                    std::move(cpu_handler),
                    "HSS",
                    nullptr,
                    enable_logging,
                    species_db,
                    gpu_seed,
                    gpu_threshold
                );
            }
#endif
            return cpu_handler;
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
                        << "for deterministic model (gamma=" << gamma_for_ou << ", apply_damping=false)";
                    ICARION::log::debug_log(msg.str());
                }
                
                // apply_damping=false: Only thermal kicks, DampingForce provides friction
                return std::make_unique<OUCollisionHandler>(gamma_for_ou, false);
            }
            
            // No OU thermalization - deterministic damping only via DampingForce
            if (enable_logging) {
                ICARION::log::debug_log(
                    "[CollisionHandlerFactory] Deterministic model without OU - "
                    "returning nullptr (use DampingForce only)"
                );
            }
            return nullptr;
        }
        
        case CollisionModel::NoCollisions: {
            if (enable_logging) {
                ICARION::log::debug_log(
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
