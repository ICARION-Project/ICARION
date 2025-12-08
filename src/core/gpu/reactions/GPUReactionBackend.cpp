// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#ifdef ICARION_USE_GPU

#include "GPUReactionBackend.h"
#include "core/log/Logger.h"

namespace icarion::gpu {

GPUReactionBackend::GPUReactionBackend(std::shared_ptr<GPUContext> context,
                                       Config cfg)
    : context_(std::move(context)), config_(cfg) {}

bool GPUReactionBackend::process_batch(ICARION::core::IonEnsemble& ensemble,
                                       const std::vector<int>& domain_indices,
                                       double dt,
                                       const ICARION::config::ReactionDatabase& reaction_db,
                                       const ICARION::config::SpeciesDatabase& species_db,
                                       const std::vector<ICARION::config::DomainConfig>& domains) {
    (void)ensemble;
    (void)domain_indices;
    (void)dt;
    (void)reaction_db;
    (void)species_db;
    (void)domains;

    if (!warned_once_) {
        ICARION::log::Logger::main()->warn(
            "GPUReactionBackend: GPU reaction path is not implemented yet – falling back to CPU.");
        warned_once_ = true;
    }
    return false;
}

}  // namespace icarion::gpu

#endif  // ICARION_USE_GPU
