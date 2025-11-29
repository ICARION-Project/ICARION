// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "PhysicsSetup.h"
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/physics/forces/MagneticFieldForce.h"
#include "core/physics/forces/DampingForce.h"
#include "core/physics/forces/SpaceChargeDirect.h"
#include "core/physics/forces/SpaceChargeGrid.h"
#include "core/physics/spacecharge/spaceChargeSolver.h"
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/integrator/strategies/RK45Strategy.h"
#include "core/integrator/strategies/BorisStrategy.h"
#include "core/physics/collisions/CollisionHandlerFactory.h"
#include "core/physics/collisions/geometryUtils.h"
#include "core/physics/reactions/ReactionHandlerFactory.h"
#include "core/log/Logger.h"
#include <algorithm>
#include <unordered_set>

namespace ICARION::setup {

PhysicsModules PhysicsSetup::initialize(
    const config::FullConfig& config,
    const std::vector<core::IonState>& ions
) {
    log::Logger::main()->info("Initializing physics modules");
    
    PhysicsModules modules;
    
    // Create force registries for each domain
    modules.force_registries = create_force_registries(config);
    
    // Add space charge forces if enabled
    if (config.physics.enable_space_charge) {
        add_space_charge_forces(modules.force_registries, config, ions);
    }
    
    // Create integration strategy
    modules.integrator = create_integrator(config);
    
    // Create collision handler
    modules.collision_handler = create_collision_handler(config);
    
    // Create reaction handler
    modules.reaction_handler = create_reaction_handler(config);
    
    return modules;
}

std::vector<std::shared_ptr<physics::ForceRegistry>> PhysicsSetup::create_force_registries(
    const config::FullConfig& config
) {
    std::vector<std::shared_ptr<physics::ForceRegistry>> registries;
    
    for (const auto& domain : config.domains) {
        auto registry = std::make_shared<physics::ForceRegistry>(domain);
        
        // Add electric field force (always)
        registry->add_force(std::make_unique<physics::ElectricFieldForce>(domain));
        
        // Add magnetic field force if configured
        if (domain.fields.magnetic.enabled) {
            registry->add_force(std::make_unique<physics::MagneticFieldForce>(domain.fields.magnetic));
        }
        
        // Add damping force for Friction model
        if (config.physics.collision_model == config::CollisionModel::Friction) {
            registry->add_force(std::make_unique<physics::DampingForce>(
                domain.environment,
                physics::DampingModel::Friction,
                nullptr  // Species DB not available here, will use ion CCS
            ));
        }
        
        registries.push_back(registry);
    }
    
    // Log summary
    log::Logger::main()->info("Created {} ForceRegistry instances (one per domain)", registries.size());
    log::Logger::main()->info("  ✓ ElectricFieldForce added to all registries");
    
    size_t mag_count = 0;
    for (const auto& domain : config.domains) {
        if (domain.fields.magnetic.enabled) mag_count++;
    }
    if (mag_count > 0) {
        log::Logger::main()->info("  ✓ MagneticFieldForce added to {} registries", mag_count);
    }
    
    if (config.physics.collision_model == config::CollisionModel::Friction) {
        log::Logger::main()->info("  ✓ DampingForce added to all registries (Friction model)");
    }
    
    return registries;
}

void PhysicsSetup::add_space_charge_forces(
    std::vector<std::shared_ptr<physics::ForceRegistry>>& registries,
    const config::FullConfig& config,
    const std::vector<core::IonState>& ions
) {
    const size_t N = ions.size();
    constexpr size_t SPACE_CHARGE_THRESHOLD = 1000;
    
    if (N < SPACE_CHARGE_THRESHOLD) {
        // Use direct N-body Coulomb (exact, but O(N²))
        log::Logger::main()->info("Space charge: Using SpaceChargeDirect (N={} < {})",
                                  N, SPACE_CHARGE_THRESHOLD);
        log::Logger::main()->info("  → Direct N-body Coulomb (exact, O(N²))");
        
        constexpr double SOFTENING_LENGTH = 1e-10;  // 0.1 nm (prevents 1/r² divergence)
        for (auto& registry : registries) {
            registry->add_force(std::make_unique<physics::SpaceChargeDirect>(SOFTENING_LENGTH));
        }
        log::Logger::main()->info("  ✓ SpaceChargeDirect added to {} registries (ε={:.2e} m)",
                                  registries.size(), SOFTENING_LENGTH);
    } else {
        // Use grid-based Poisson solver (fast, but approximate)
        log::Logger::main()->info("Space charge: Using SpaceChargeGrid (N={} >= {})",
                                  N, SPACE_CHARGE_THRESHOLD);
        log::Logger::main()->info("  → Grid-based Poisson solver (fast, O(N log N))");
        
        // Estimate domain size from ion initial positions
        Vec3 min_pos = ions[0].pos;
        Vec3 max_pos = ions[0].pos;
        for (const auto& ion : ions) {
            min_pos.x = std::min(min_pos.x, ion.pos.x);
            min_pos.y = std::min(min_pos.y, ion.pos.y);
            min_pos.z = std::min(min_pos.z, ion.pos.z);
            max_pos.x = std::max(max_pos.x, ion.pos.x);
            max_pos.y = std::max(max_pos.y, ion.pos.y);
            max_pos.z = std::max(max_pos.z, ion.pos.z);
        }
        
        Vec3 domain_size = {max_pos.x - min_pos.x, max_pos.y - min_pos.y, max_pos.z - min_pos.z};
        Vec3 domain_center = {(min_pos.x + max_pos.x) / 2, (min_pos.y + max_pos.y) / 2, (min_pos.z + max_pos.z) / 2};
        
        // Add 50% margin to domain size (ions will move)
        domain_size = domain_size * 1.5;
        
        // Grid resolution: Aim for ~1mm cells (adjust based on domain)
        constexpr int TARGET_GRID_SIZE = 64;  // 64³ = 262k cells (good balance)
        double cell_size_x = domain_size.x / TARGET_GRID_SIZE;
        double cell_size_y = domain_size.y / TARGET_GRID_SIZE;
        double cell_size_z = domain_size.z / TARGET_GRID_SIZE;
        
        // Use uniform cell size (max of xyz)
        double cell_size = std::max({cell_size_x, cell_size_y, cell_size_z, 1e-4});  // Min 0.1mm
        
        Vec3 grid_origin = {
            domain_center.x - (TARGET_GRID_SIZE * cell_size) / 2,
            domain_center.y - (TARGET_GRID_SIZE * cell_size) / 2,
            domain_center.z - (TARGET_GRID_SIZE * cell_size) / 2
        };
        
        log::Logger::main()->info("  Grid: {}³ cells, {:.2e} m cell size", TARGET_GRID_SIZE, cell_size);
        log::Logger::main()->info("  Domain: [{:.3f}, {:.3f}] x [{:.3f}, {:.3f}] x [{:.3f}, {:.3f}] mm",
                                  grid_origin.x * 1e3, (grid_origin.x + TARGET_GRID_SIZE * cell_size) * 1e3,
                                  grid_origin.y * 1e3, (grid_origin.y + TARGET_GRID_SIZE * cell_size) * 1e3,
                                  grid_origin.z * 1e3, (grid_origin.z + TARGET_GRID_SIZE * cell_size) * 1e3);
        
        // Create solver
        auto sc_solver = std::make_shared<SpaceChargeSolver>(
            TARGET_GRID_SIZE, TARGET_GRID_SIZE, TARGET_GRID_SIZE,
            cell_size, cell_size, cell_size,
            grid_origin
        );
        
        // Wrap solver in IForce interface and add to registries
        for (auto& registry : registries) {
            registry->add_force(std::make_unique<physics::SpaceChargeGrid>(sc_solver));
        }
        log::Logger::main()->info("  ✓ SpaceChargeGrid added to {} registries", registries.size());
    }
}

std::shared_ptr<integrator::IIntegrationStrategy> PhysicsSetup::create_integrator(
    const config::FullConfig& config
) {
    std::shared_ptr<integrator::IIntegrationStrategy> integrator;
    
    if (config.simulation.integrator == "RK4" || config.simulation.integrator == "rk4") {
        integrator = std::make_shared<integrator::RK4Strategy>();
        log::Logger::main()->info("Using RK4 integrator");
    } else if (config.simulation.integrator == "RK45" || config.simulation.integrator == "rk45") {
        integrator = std::make_shared<integrator::RK45Strategy>();
        log::Logger::main()->info("Using RK45 integrator");
    } else if (config.simulation.integrator == "Boris" || config.simulation.integrator == "boris") {
        integrator = std::make_shared<integrator::BorisStrategy>();
        log::Logger::main()->info("Using Boris integrator");
    } else {
        log::Logger::main()->warn("Unknown integrator '{}', defaulting to RK45", config.simulation.integrator);
        integrator = std::make_shared<integrator::RK45Strategy>();
    }
    
    return integrator;
}

std::shared_ptr<physics::ICollisionHandler> PhysicsSetup::create_collision_handler(
    const config::FullConfig& config
) {
    // Load geometry map for EHSS (if needed)
    std::unique_ptr<physics::GeometryMap> geometry_map_ptr = nullptr;
    const physics::GeometryMap* geometry_map = nullptr;
    
    if (config.physics.collision_model == config::CollisionModel::EHSS) {
        // Collect all ion species from config
        std::unordered_set<std::string> species_ids;
        for (const auto& species : config.ions.species) {
            species_ids.insert(species.species_id);
        }
        
        try {
            log::Logger::main()->info("Loading molecular geometries for EHSS collision model");
            geometry_map_ptr = std::make_unique<physics::GeometryMap>(
                physics::load_geometry_map(species_ids, "/home/chsch95/ICARION/data/molecules/", false)
            );
            geometry_map = geometry_map_ptr.get();
            log::Logger::main()->info("Loaded {} molecular geometries", geometry_map->size());
        } catch (const std::exception& e) {
            log::Logger::main()->error("Failed to load molecular geometries: {}", e.what());
            log::Logger::main()->warn("Falling back to HSS collision model");
            // Don't exit, let CollisionHandlerFactory handle the fallback
        }
    }
    
    constexpr double gamma_for_ou = 0.0;  // OU damping coefficient not used for stochastic models
    
    return physics::CollisionHandlerFactory::create(
        config.physics,
        geometry_map,
        gamma_for_ou,
        false,  // enable_logging
        &config.species_db
    );
}

std::shared_ptr<physics::IReactionHandler> PhysicsSetup::create_reaction_handler(
    const config::FullConfig& config
) {
    return physics::ReactionHandlerFactory::create(config.physics);
}

}  // namespace ICARION::setup
