// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "SimulationEngine.h"

#include "core/physics/forces/DampingForce.h"
#include "core/utils/Profiler.h"

#include <sstream>

#ifndef ICARION_VERSION
#define ICARION_VERSION "unknown"
#endif

namespace ICARION {
namespace integrator {

#ifdef ICARION_USE_GPU
void SimulationEngine::initialize_gpu(bool enable_gpu) {
    PROFILE_SCOPE_IF_ENABLED("GPU Initialization");

    // Check if GPU is enabled in config
    if (!enable_gpu) {
        output_manager_->log_progress("GPU: Disabled in configuration, using CPU-only");
        return;
    }

    output_manager_->log_progress(
        std::string("GPU: Disabled for v") + ICARION_VERSION +
        " (experimental path remains built but not used). Running CPU-only.");
    return;

    output_manager_->log_progress("GPU: Experimental path (E/B-only). Space-charge, damping, magnetic forces, and multi-domain batches are NOT supported; falling back to CPU in those cases.");
    bool friction_damping_present = false;

    // Validate force setup: only a single ElectricFieldForce per domain, no space charge
    auto gpu_forces_supported = [&]() -> bool {
        if (config_.physics.enable_space_charge) return false;
        for (const auto& reg : force_registries_) {
            if (!reg) continue;
            if (reg->space_charge_model()) return false;
            const auto& forces = reg->forces();
            for (const auto& f : forces) {
                if (dynamic_cast<physics::DampingForce*>(f.get())) {
                    const auto name = f->name();
                    if (name.find("Friction") != std::string::npos) {
                        friction_damping_present = true;
                    }
                }
            }
            if (forces.size() != 1) return false;
            if (forces.front()->name().find("ElectricField") != 0) return false;
        }
        return true;
    };

    if (!gpu_forces_supported()) {
        output_manager_->log_progress("GPU: Disabled — unsupported force combination (requires exactly one ElectricFieldForce per domain, no space charge/damping/magnetic).");
        return;
    }

    try {
        // Check if CUDA is available
        if (!icarion::gpu::GPUContext::is_cuda_available()) {
            output_manager_->log_progress("GPU: CUDA not available, using CPU-only");
            return;
        }

        // Create GPU context (device 0)
        gpu_context_ = icarion::gpu::GPUContext::create(0);
        if (!gpu_context_) {
            output_manager_->log_progress("GPU: Failed to create context, using CPU-only");
            return;
        }

        // Log GPU properties
        const auto& props = gpu_context_->get_properties();
        std::ostringstream gpu_msg;
        gpu_msg << "GPU: " << props.name
                << " (Compute " << props.compute_capability_major << "." << props.compute_capability_minor
                << ", " << props.total_memory / (1024*1024) << " MB)";
        output_manager_->log_progress(gpu_msg.str());

        // Create GPU collision helper (if collision handler exists)
        if (collision_handler_) {
            // Determine collision model from config
            auto collision_model_enum = config_.physics.collision_model;
            std::string collision_model_str;

            // Only create GPU helper for supported models (HSS, EHSS)
            bool gpu_supported = false;
            if (collision_model_enum == config::CollisionModel::HSS) {
                collision_model_str = "HSS";
                gpu_supported = true;
            } else if (collision_model_enum == config::CollisionModel::EHSS) {
                collision_model_str = "EHSS";
                gpu_supported = true;
            }

            if (gpu_supported) {
                gpu_collision_helper_ = icarion::gpu::GPUCollisionHelper::create(
                    *gpu_context_,
                    gpu_collision_threshold_,
                    collision_model_str,
                    config_.simulation.rng_seed
                );

                if (gpu_collision_helper_) {
                    std::ostringstream collision_msg;
                    collision_msg << "GPU: Collision processing enabled for N >= "
                                  << gpu_collision_threshold_ << " ions (" << collision_model_str << ")";
                    output_manager_->log_progress(collision_msg.str());

                    // For EHSS, upload geometry data
                    if (collision_model_enum == config::CollisionModel::EHSS) {
                        // TODO(v1.1): Extract geometry from species database and upload to GPU
                        // This requires access to species_db molecular geometries
                        output_manager_->log_progress("GPU: EHSS geometry upload deferred (using HSS fallback for now)");
                    }
                }
            }
        }

        // NOTE: GPU Space Charge (P³M) initialization deferred to first use
        // via lazy initialization in try_gpu_space_charge().
        // This avoids coupling to SpaceChargeConfig which may not exist yet.
        // Full integration with DomainConfig.space_charge pending.

        (void)friction_damping_present; // GPU runtime disabled
    }
    catch (const std::exception& e) {
        output_manager_->log_progress(std::string("GPU: Initialization failed: ") + e.what());
        gpu_space_charge_.reset();
        gpu_collision_helper_.reset();
        gpu_context_.reset();
    }
}
#endif

} // namespace integrator
} // namespace ICARION
