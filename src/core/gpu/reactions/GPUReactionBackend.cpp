// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#ifdef ICARION_USE_GPU

#include "GPUReactionBackend.h"
#include "core/log/Logger.h"
#include "core/physics/reactions/StochasticReactionHandler.h"
#include "core/gpu/reactions/reaction_kernels.cuh"
#include "utils/constants.h"
#include <unordered_map>
#include <cuda_runtime.h>
#include <curand_kernel.h>

namespace icarion::gpu {

GPUReactionBackend::GPUReactionBackend(std::shared_ptr<GPUContext> context,
                                       Config cfg)
    : context_(std::move(context)), config_(cfg) {
    cpu_fallback_ = std::make_unique<ICARION::physics::StochasticReactionHandler>(cfg.enable_logging);
}

GPUReactionBackend::~GPUReactionBackend() {
    if (d_rng_states_) cudaFree(d_rng_states_);
    if (d_species_idx_) cudaFree(d_species_idx_);
    if (d_domain_idx_) cudaFree(d_domain_idx_);
    if (d_active_) cudaFree(d_active_);
    if (d_born_) cudaFree(d_born_);
    if (d_reactions_) cudaFree(d_reactions_);
    if (d_offsets_) cudaFree(d_offsets_);
}

void GPUReactionBackend::ensure_rng_states(size_t n) {
    if (rng_capacity_ >= n) {
        return;
    }
    if (d_rng_states_) cudaFree(d_rng_states_);
    cudaMalloc(&d_rng_states_, n * sizeof(curandStateXORWOW));
    init_rng_states_xorwow(d_rng_states_, n, config_.rng_seed, context_->get_stream());
    rng_capacity_ = n;
}

void GPUReactionBackend::ensure_ion_buffers(size_t n) {
    if (ion_capacity_ >= n) return;
    if (d_species_idx_) cudaFree(d_species_idx_);
    if (d_domain_idx_) cudaFree(d_domain_idx_);
    if (d_active_) cudaFree(d_active_);
    if (d_born_) cudaFree(d_born_);
    cudaMalloc(&d_species_idx_, n * sizeof(uint32_t));
    cudaMalloc(&d_domain_idx_, n * sizeof(int32_t));
    cudaMalloc(&d_active_, n * sizeof(uint8_t));
    cudaMalloc(&d_born_, n * sizeof(uint8_t));
    ion_capacity_ = n;
}

void GPUReactionBackend::ensure_reaction_buffers(size_t reaction_count, size_t offset_count) {
    if (reaction_capacity_ < reaction_count) {
        if (d_reactions_) cudaFree(d_reactions_);
        cudaMalloc(&d_reactions_, reaction_count * sizeof(DeviceReaction));
        reaction_capacity_ = reaction_count;
    }
    if (offset_capacity_ < offset_count) {
        if (d_offsets_) cudaFree(d_offsets_);
        cudaMalloc(&d_offsets_, offset_count * sizeof(int32_t));
        offset_capacity_ = offset_count;
    }
}

GPUReactionBackend::FlattenedReactions GPUReactionBackend::flatten_reactions(
    const ICARION::config::ReactionDatabase& reaction_db,
    const ICARION::config::SpeciesDatabase& species_db,
    const ICARION::core::IonEnsemble& ensemble,
    const std::vector<ICARION::config::DomainConfig>& domains) {
    (void)species_db;
    FlattenedReactions flat;
    flat.n_domains = static_cast<int>(domains.size());
    flat.reaction_offsets.resize(flat.n_domains + 1, 0);

    std::unordered_map<std::string, int> species_to_idx;
    auto* pool_mut = const_cast<std::vector<std::string>*>(ensemble.species_pool());
    for (size_t i = 0; i < pool_mut->size(); ++i) {
        species_to_idx[pool_mut->at(i)] = static_cast<int>(i);
    }

    int32_t total = 0;
    for (size_t d = 0; d < domains.size(); ++d) {
        const auto& env = domains[d].environment;

        std::unordered_map<std::string, double> concentrations;
        if (!env.gas_mixture.empty()) {
            for (const auto& comp : env.gas_mixture) {
                if (!comp.participates_in_reactions) continue;
                concentrations[comp.species] = comp.density_m3;
            }
        } else {
            concentrations[env.gas_species] = env.particle_density_m_3;
        }

        for (const auto& rxn : reaction_db.reactions) {
            auto it_r = species_to_idx.find(rxn.reactant);
            if (it_r == species_to_idx.end()) {
                continue;
            }
            auto it_p = species_to_idx.find(rxn.product);
            if (it_p == species_to_idx.end() && species_db.has(rxn.product)) {
                int new_idx = static_cast<int>(pool_mut->size());
                pool_mut->push_back(rxn.product);
                species_to_idx[rxn.product] = new_idx;
                it_p = species_to_idx.find(rxn.product);
            }
            if (it_p == species_to_idx.end()) continue;
            double k_eff = rxn.effective_rate_s(env.temperature_K, concentrations);
            DeviceReaction dr;
            dr.reactant_idx = it_r->second;
            dr.product_idx = it_p->second;
            dr.k_eff = k_eff;
            flat.reactions.push_back(dr);
            total++;
        }
        flat.reaction_offsets[d + 1] = total;
    }
    return flat;
}

bool GPUReactionBackend::process_batch(ICARION::core::IonEnsemble& ensemble,
                                       const std::vector<int>& domain_indices,
                                       double dt,
                                       const ICARION::config::ReactionDatabase& reaction_db,
                                       const ICARION::config::SpeciesDatabase& species_db,
                                       const std::vector<ICARION::config::DomainConfig>& domains,
                                       std::vector<ICARION::physics::PhysicsRng>& rng_pool) {
    (void)rng_pool;  // GPU backend keeps its own RNG states; pool unused here.
    if (!cpu_fallback_) {
        return false;
    }
    const size_t n = ensemble.size();
    if (n < config_.gpu_threshold) {
        return false;
    }

    auto flat = flatten_reactions(reaction_db, species_db, ensemble, domains);
    if (flat.reactions.empty()) {
        return false;
    }

    ensure_rng_states(n);
    ensure_ion_buffers(n);

    cudaMemcpy(d_species_idx_, ensemble.species_id_indices(), n * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_domain_idx_, domain_indices.data(), n * sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_active_, ensemble.active_data(), n * sizeof(uint8_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_born_, ensemble.born_data(), n * sizeof(uint8_t), cudaMemcpyHostToDevice);

    ensure_reaction_buffers(flat.reactions.size(), flat.reaction_offsets.size());
    cudaMemcpy(d_reactions_, flat.reactions.data(), flat.reactions.size() * sizeof(DeviceReaction), cudaMemcpyHostToDevice);
    cudaMemcpy(d_offsets_, flat.reaction_offsets.data(), flat.reaction_offsets.size() * sizeof(int32_t), cudaMemcpyHostToDevice);

    auto stream = context_->get_stream();
    launch_reaction_kernel(d_species_idx_,
                           d_domain_idx_,
                           d_active_,
                           d_born_,
                           n,
                           d_reactions_,
                           d_offsets_,
                           flat.n_domains,
                           d_rng_states_,
                           dt,
                           stream);
    cudaStreamSynchronize(stream);

    std::vector<uint32_t> species_idx_host(n);
    cudaMemcpy(species_idx_host.data(), d_species_idx_, n * sizeof(uint32_t), cudaMemcpyDeviceToHost);

    const auto* pool = ensemble.species_pool();
    auto* mass = ensemble.mass_data();
    auto* charge = ensemble.charge_data();
    auto* ccs = ensemble.CCS_data();
    auto* mobility = ensemble.mobility_data();
    auto* species_indices_mut = const_cast<uint32_t*>(ensemble.species_id_indices());

    for (size_t i = 0; i < n; ++i) {
        if (species_idx_host[i] == ensemble.species_id_indices()[i]) {
            continue;
        }
        uint32_t idx = species_idx_host[i];
        if (idx >= pool->size()) {
            continue;
        }
        const std::string& new_species = pool->at(idx);
        if (!species_db.has(new_species)) {
            continue;
        }
        const auto& props = species_db.get(new_species);
        species_indices_mut[i] = idx;
        mass[i] = props.mass_kg;
        charge[i] = props.charge * ELEM_CHARGE_C;
        if (props.CCS_A2) {
            ccs[i] = (*props.CCS_A2) * 1e-20;
        }
        if (props.mobility_cm2Vs) {
            mobility[i] = *props.mobility_cm2Vs;
        }
    }

    return true;
}

}  // namespace icarion::gpu

#endif  // ICARION_USE_GPU
