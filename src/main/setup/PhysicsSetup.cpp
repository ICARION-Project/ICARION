// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "PhysicsSetup.h"
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/physics/forces/MagneticFieldForce.h"
#include "core/physics/forces/DampingForce.h"
#include "core/physics/spacecharge/SpaceChargeModelFactory.h"
#include "core/integrator/strategies/IntegrationStrategyFactory.h"
#include "core/physics/collisions/CollisionHandlerFactory.h"
#include "core/physics/collisions/geometryUtils.h"
#include "core/physics/reactions/ReactionHandlerFactory.h"
#include "core/io/fieldArrayLoader.h"
#include "fieldsolver/utils/GridFieldProvider.h"
#include "fieldsolver/utils/CompositeFieldProvider.h"
#include "core/log/Logger.h"
#include "core/types/IonEnsemble.h"
#include <algorithm>
#include <cctype>
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
    if (registries.empty() || N == 0) {
        return;
    }

    const size_t domain_count = registries.size();
    auto create_model_for_domain = [&](size_t idx) {
        return physics::SpaceChargeModelFactory::create(
            config, config.domains[idx], N);
    };

    auto first_model = create_model_for_domain(0);
    if (!first_model) {
        log::Logger::main()->warn(
            "Space charge enabled, but no model could be created for domain '{}'. "
            "Skipping space-charge coupling.",
            config.domains.front().name);
        return;
    }

    if (first_model->name() == "SpaceChargeDirectModel") {
        for (auto& registry : registries) {
            registry->set_space_charge_model(first_model);
        }
        log::Logger::main()->info(
            "Space charge: Using {} shared across {} domains",
            first_model->name(),
            registries.size());
        return;
    }

    registries[0]->set_space_charge_model(first_model);
    log::Logger::main()->info("Space charge: {} uses {}",
                              config.domains[0].name,
                              first_model->name());

    size_t assigned_models = first_model ? 1 : 0;
    for (size_t idx = 1; idx < domain_count; ++idx) {
        auto model = create_model_for_domain(idx);
        if (!model) {
            log::Logger::main()->warn(
                "Space charge: No model available for domain '{}'; "
                "space charge disabled for this domain.",
                config.domains[idx].name);
            continue;
        }

        registries[idx]->set_space_charge_model(model);
        assigned_models++;
        log::Logger::main()->info("Space charge: {} uses {}",
                                  config.domains[idx].name,
                                  model->name());
    }

    if (assigned_models == 0) {
        log::Logger::main()->warn("Space charge: No domains received a model; "
                                  "space charge remains disabled.");
    }
}

std::shared_ptr<integrator::IIntegrationStrategy> PhysicsSetup::create_integrator(
    const config::FullConfig& config
) {
    std::string name = config.simulation.integrator;
    if (name.empty()) {
        name = "RK45";
    }
    std::transform(name.begin(), name.end(), name.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    std::shared_ptr<integrator::IIntegrationStrategy> integrator;
    try {
        auto ptr = integrator::IntegrationStrategyFactory::create(name, &config);
        integrator = std::shared_ptr<integrator::IIntegrationStrategy>(std::move(ptr));
    } catch (const std::exception& e) {
        log::Logger::main()->warn("Integrator '{}' invalid ({}); falling back to RK45",
                                  config.simulation.integrator, e.what());
        auto fallback = integrator::IntegrationStrategyFactory::create("RK45", &config);
        integrator = std::shared_ptr<integrator::IIntegrationStrategy>(std::move(fallback));
    }

    log::Logger::main()->info("Using {} integrator", integrator->name());
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
        &config.species_db,
        config.simulation.enable_gpu,
        static_cast<unsigned long long>(config.simulation.rng_seed),
        5000
    );
}

std::shared_ptr<physics::IReactionHandler> PhysicsSetup::create_reaction_handler(
    const config::FullConfig& config
) {
    return physics::ReactionHandlerFactory::create(config.physics);
}

}  // namespace ICARION::setup
