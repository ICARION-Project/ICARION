// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "PhysicsSetup.h"
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/physics/forces/MagneticFieldForce.h"
#include "core/physics/forces/DampingForce.h"
#include "core/physics/spacecharge/SpaceChargeModelFactory.h"
#include "core/physics/forces/SpaceChargeDirect.h"
#include "core/physics/forces/SpaceChargeGrid.h"
#include "core/physics/spacecharge/spaceChargeSolver.h"

#ifdef ICARION_USE_GPU
#include "core/physics/forces/SpaceChargeGPU.h"
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/spacecharge/GPUSpaceChargeP3M.h"
#endif
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/integrator/strategies/RK45Strategy.h"
#include "core/integrator/strategies/BorisStrategy.h"
#include "core/physics/collisions/CollisionHandlerFactory.h"
#include "core/physics/collisions/geometryUtils.h"
#include "core/physics/reactions/ReactionHandlerFactory.h"
#include "core/io/fieldArrayLoader.h"
#include "fieldsolver/utils/GridFieldProvider.h"
#include "fieldsolver/utils/CompositeFieldProvider.h"
#include "core/log/Logger.h"
#include "core/types/IonEnsemble.h"
#include <algorithm>
#include <unordered_set>

namespace ICARION::setup {

PhysicsModules PhysicsSetup::initialize(
    const config::FullConfig& config,
    const core::IonEnsemble& ions
) {
    log::Logger::main()->info("Initializing physics modules");
    
    PhysicsModules modules;
    
    // Create force registries for each domain
    // (This will populate field_arrays_storage_ if field arrays are present)
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
        
        // Check for field arrays and create field provider if present
        std::shared_ptr<::IFieldProvider> electric_field_provider = nullptr;
        
        if (!domain.fields.field_array_terms.empty()) {
            log::Logger::main()->info("Loading field arrays for domain '{}'", domain.name);
            
            // Load all field arrays for this domain
            std::vector<CompositeFieldProvider::FieldTerm> composite_terms;
            
            for (const auto& term : domain.fields.field_array_terms) {
                try {
                    FieldArray field = load_field_array(term.file);
                    
                    std::string scale_desc;
                    using ScaleKind = config::FieldsConfig::FieldArrayTerm::ScaleKind;
                    switch (term.kind) {
                        case ScaleKind::Constant:
                            scale_desc = fmt::format("constant={:.2f}V", term.constant);
                            break;
                        case ScaleKind::DC_Axial:
                            scale_desc = "scale=DC.axial_V";
                            break;
                        case ScaleKind::DC_Quad:
                            scale_desc = "scale=DC.quad_V";
                            break;
                        case ScaleKind::DC_Radial:
                            scale_desc = "scale=DC.radial_V";
                            break;
                        case ScaleKind::RF:
                            scale_desc = fmt::format("scale=RF (f={:.0f}Hz, φ={:.2f}rad)", 
                                                    term.frequency_Hz, term.phase_rad);
                            break;
                    }
                    
                    log::Logger::main()->info("  ✓ Loaded: {} (grid: {}×{}×{}, {})",
                                              term.file,
                                              field.nx,
                                              field.ny,
                                              field.nz,
                                              scale_desc);
                    
                    // Store field array on heap to get stable pointer
                    auto field_ptr = std::make_unique<FieldArray>(std::move(field));
                    const FieldArray* raw_ptr = field_ptr.get();
                    field_arrays_storage_.push_back(std::move(field_ptr));
                    
                    // Create composite term
                    CompositeFieldProvider::FieldTerm composite_term;
                    composite_term.field_array = raw_ptr;
                    composite_term.kind = term.kind;
                    composite_term.constant_scale = term.constant;
                    composite_term.frequency_Hz = term.frequency_Hz;
                    composite_term.phase_rad = term.phase_rad;
                    composite_terms.push_back(composite_term);
                    
                } catch (const std::exception& e) {
                    log::Logger::main()->error("Failed to load field array: {} - {}",
                                               term.file, e.what());
                    throw;
                }
            }
            
            // Create field provider with superposition support
            if (!composite_terms.empty()) {
                if (composite_terms.size() == 1) {
                    // Single field array: use simple GridFieldProvider for efficiency
                    electric_field_provider = std::make_shared<GridFieldProvider>(composite_terms[0].field_array);
                    log::Logger::main()->info("  ✓ Created GridFieldProvider for domain '{}'", domain.name);
                } else {
                    // Multiple field arrays: use CompositeFieldProvider for superposition
                    electric_field_provider = std::make_shared<CompositeFieldProvider>(composite_terms, &domain);
                    log::Logger::main()->info("  ✓ Created CompositeFieldProvider with {} field terms for domain '{}'", 
                                             composite_terms.size(), domain.name);
                }
            }
        }
        
        // Add electric field force
        if (electric_field_provider) {
            // Use field provider mode (interpolated fields from FieldArray) and set SSOT model
            registry->add_force(std::make_unique<physics::ElectricFieldForce>(electric_field_provider));
            auto model = std::make_unique<config::FieldProviderModel>(electric_field_provider);
            registry->set_field_model(model.get());
            field_models_storage_.push_back(std::move(model));
        } else {
            // Use analytical mode (instrument-specific fields) and set SSOT model
            registry->add_force(std::make_unique<physics::ElectricFieldForce>(domain));
            auto model = std::make_unique<config::AnalyticalFieldModel>(domain);
            registry->set_field_model(model.get());
            field_models_storage_.push_back(std::move(model));
        }
        
        // Add magnetic field force if configured (but NOT for Boris integrator)
        // NOTE: Boris integrator handles magnetic fields internally via rotation,
        //       so it doesn't need/want MagneticFieldForce in the registry
        if (domain.fields.magnetic.enabled && config.simulation.integrator != "Boris") {
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
        if (domain.fields.magnetic.enabled && config.simulation.integrator != "Boris") {
            mag_count++;
        }
    }
    if (mag_count > 0) {
        log::Logger::main()->info("  ✓ MagneticFieldForce added to {} registries (not Boris)", mag_count);
    } else if (config.simulation.integrator == "Boris") {
        size_t boris_magnetic_domains = 0;
        for (const auto& domain : config.domains) {
            if (domain.fields.magnetic.enabled) boris_magnetic_domains++;
        }
        if (boris_magnetic_domains > 0) {
            log::Logger::main()->info("  ℹ MagneticFieldForce skipped for {} domains (Boris handles B-field internally)", boris_magnetic_domains);
        }
    }
    
    if (config.physics.collision_model == config::CollisionModel::Friction) {
        log::Logger::main()->info("  ✓ DampingForce added to all registries (Friction model)");
    }
    
    return registries;
}

void PhysicsSetup::add_space_charge_forces(
    std::vector<std::shared_ptr<physics::ForceRegistry>>& registries,
    const config::FullConfig& config,
    const core::IonEnsemble& ions
) {
    const size_t N = ions.size();
    const size_t domain_count = registries.size();

    bool model_configured = false;
    for (size_t idx = 0; idx < domain_count; ++idx) {
        const auto& domain = config.domains[idx];
        auto model = physics::SpaceChargeModelFactory::create(config, domain, N);
        if (model) {
            registries[idx]->set_space_charge_model(model);
            log::Logger::main()->info("Space charge: {} uses {}", domain.name, model->name());
            model_configured = true;
        }
    }

    if (model_configured) {
        return;  // vNext path enabled
    }

    constexpr size_t SPACE_CHARGE_THRESHOLD = 1000;
    
    // Priority: GPU > Grid (CPU) > Direct (CPU)
    // - GPU: N >= 1000 + GPU available → SpaceChargeGPU (10-40× faster than CPU Grid)
    // - Grid: N >= 1000 + no GPU → SpaceChargeGrid (10-100× faster than Direct)
    // - Direct: N < 1000 → SpaceChargeDirect (exact, competitive for small N)
    
#ifdef ICARION_USE_GPU
    // Try GPU first (if N >= threshold and GPU available)
    if (N >= SPACE_CHARGE_THRESHOLD) {
        try {
            // Attempt to create GPU context
            auto gpu_ctx = icarion::gpu::GPUContext::create(0);  // Device 0
            
            if (gpu_ctx && gpu_ctx->is_valid()) {
                log::Logger::main()->info("Space charge: Using SpaceChargeGPU (N={} >= {}, GPU available)",
                                          N, SPACE_CHARGE_THRESHOLD);
                log::Logger::main()->info("  → GPU P³M algorithm (O(N log N), cuFFT-based Poisson solver)");
                
                // Estimate domain size from ion initial positions (SoA)
                const auto* px = ions.pos_x_data();
                const auto* py = ions.pos_y_data();
                const auto* pz = ions.pos_z_data();
                Vec3 min_pos{px[0], py[0], pz[0]};
                Vec3 max_pos{px[0], py[0], pz[0]};
                for (size_t i = 1; i < N; ++i) {
                    min_pos.x = std::min(min_pos.x, px[i]);
                    min_pos.y = std::min(min_pos.y, py[i]);
                    min_pos.z = std::min(min_pos.z, pz[i]);
                    max_pos.x = std::max(max_pos.x, px[i]);
                    max_pos.y = std::max(max_pos.y, py[i]);
                    max_pos.z = std::max(max_pos.z, pz[i]);
                }
                
                Vec3 domain_size = {max_pos.x - min_pos.x, max_pos.y - min_pos.y, max_pos.z - min_pos.z};
                
                // Add 50% margin (ions will move)
                domain_size = domain_size * 1.5;
                
                // Ensure minimum domain size (1cm cube)
                domain_size.x = std::max(domain_size.x, 0.01);
                domain_size.y = std::max(domain_size.y, 0.01);
                domain_size.z = std::max(domain_size.z, 0.01);
                
                Vec3 domain_center = {
                    (min_pos.x + max_pos.x) / 2,
                    (min_pos.y + max_pos.y) / 2,
                    (min_pos.z + max_pos.z) / 2
                };
                
                // Configure P³M solver: Target 30µm cells, 32-256 grid resolution
                constexpr double TARGET_CELL_SIZE = 30e-6;  // 30 µm
                int nx = std::clamp(static_cast<int>(domain_size.x / TARGET_CELL_SIZE), 32, 256);
                int ny = std::clamp(static_cast<int>(domain_size.y / TARGET_CELL_SIZE), 32, 256);
                int nz = std::clamp(static_cast<int>(domain_size.z / TARGET_CELL_SIZE), 32, 256);
                
                icarion::gpu::GPUSpaceChargeP3M::Config p3m_config;
                p3m_config.grid_nx = nx;
                p3m_config.grid_ny = ny;
                p3m_config.grid_nz = nz;
                p3m_config.domain_min = {
                    domain_center.x - domain_size.x / 2,
                    domain_center.y - domain_size.y / 2,
                    domain_center.z - domain_size.z / 2
                };
                p3m_config.domain_max = {
                    domain_center.x + domain_size.x / 2,
                    domain_center.y + domain_size.y / 2,
                    domain_center.z + domain_size.z / 2
                };
                
                log::Logger::main()->info("  Grid: {}×{}×{} cells, {:.1f}×{:.1f}×{:.1f} µm cell size",
                                          nx, ny, nz,
                                          domain_size.x / nx * 1e6,
                                          domain_size.y / ny * 1e6,
                                          domain_size.z / nz * 1e6);
                log::Logger::main()->info("  Domain: [{:.3f}, {:.3f}] x [{:.3f}, {:.3f}] x [{:.3f}, {:.3f}] mm",
                                          p3m_config.domain_min.x * 1e3, p3m_config.domain_max.x * 1e3,
                                          p3m_config.domain_min.y * 1e3, p3m_config.domain_max.y * 1e3,
                                          p3m_config.domain_min.z * 1e3, p3m_config.domain_max.z * 1e3);
                
                // Create GPU P³M solver
                auto gpu_solver_unique = icarion::gpu::GPUSpaceChargeP3M::create(*gpu_ctx, p3m_config);
                
                if (gpu_solver_unique) {
                    // Convert unique_ptr to shared_ptr (multiple registries need same solver)
                    auto gpu_solver = std::shared_ptr<icarion::gpu::GPUSpaceChargeP3M>(std::move(gpu_solver_unique));
                    
                    // Wrap in IForce interface
                    for (auto& registry : registries) {
                        registry->add_force(std::make_unique<physics::SpaceChargeGPU>(gpu_solver));
                    }
                    log::Logger::main()->info("  ✓ SpaceChargeGPU added to {} registries", registries.size());
                    return;  // GPU setup successful - exit function
                } else {
                    log::Logger::main()->warn("  ⚠ GPU P³M solver creation failed - falling back to CPU Grid");
                }
            } else {
                log::Logger::main()->info("  ℹ GPU context unavailable - using CPU Grid");
            }
        } catch (const std::exception& e) {
            log::Logger::main()->warn("  ⚠ GPU initialization failed: {} - falling back to CPU Grid", e.what());
        }
    }
#endif
    
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
        
        // Estimate domain size from ion initial positions (SoA)
        const auto* px = ions.pos_x_data();
        const auto* py = ions.pos_y_data();
        const auto* pz = ions.pos_z_data();
        Vec3 min_pos{px[0], py[0], pz[0]};
        Vec3 max_pos{px[0], py[0], pz[0]};
        for (size_t i = 1; i < N; ++i) {
            min_pos.x = std::min(min_pos.x, px[i]);
            min_pos.y = std::min(min_pos.y, py[i]);
            min_pos.z = std::min(min_pos.z, pz[i]);
            max_pos.x = std::max(max_pos.x, px[i]);
            max_pos.y = std::max(max_pos.y, py[i]);
            max_pos.z = std::max(max_pos.z, pz[i]);
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
