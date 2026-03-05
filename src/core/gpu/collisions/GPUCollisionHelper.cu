// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
/**
 * @file GPUCollisionHelper.cu
 * @brief GPU collision helper implementation (experimental; EHSS geometry mapping stubbed for v1.0.0)
 */

/**
 * @file GPUCollisionHelper.cu
 * @brief CUDA implementation (kernel launches and device code only)
 * 
 * Factory/constructor/destructor moved to _host.cpp to avoid nvcc issues with
 * GPUContext. This file contains only GPU kernel launches and device memory
 * management. Species-index mapping for EHSS is minimal (single-species unless
 * caller fills indices) and results are not yet validated against CPU.
 */

#include "GPUCollisionHelper.h"
#ifndef __CUDACC__
#include "core/log/Logger.h"
#endif
#include "collision_kernels_gpu.cuh"
#include "utils/constants.h"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_set>

// Minimal C++ utilities (avoid STL headers with nvcc)
#define MAX(a, b) ((a) > (b) ? (a) : (b))

namespace icarion {
namespace gpu {

GPUCollisionHelper::~GPUCollisionHelper() {
    if (d_curand_states_) {
        cudaFree(d_curand_states_);
    }
    
    // Free persistent buffers
    if (d_vx_) cudaFree(d_vx_);
    if (d_vy_) cudaFree(d_vy_);
    if (d_vz_) cudaFree(d_vz_);
    if (d_mass_) cudaFree(d_mass_);
    if (d_ccs_) cudaFree(d_ccs_);
    if (d_active_) cudaFree(d_active_);
    if (d_species_indices_) cudaFree(d_species_indices_);
    
    // Free geometry GPU memory
    if (geometry_gpu_) {
        if (geometry_uploaded_) {
            cudaFree(geometry_gpu_->atom_x);
            cudaFree(geometry_gpu_->atom_y);
            cudaFree(geometry_gpu_->atom_z);
            cudaFree(geometry_gpu_->atom_radii);
            cudaFree(geometry_gpu_->atom_counts);
            cudaFree(geometry_gpu_->atom_offsets);
        }
        delete geometry_gpu_;
    }
}

void GPUCollisionHelper::initialize_curand_states(size_t n_states) {
    if (curand_initialized_ && n_states <= n_curand_states_) {
        return;  // Already initialized with sufficient states
    }
    
    // Free old states if reallocating
    if (d_curand_states_) {
        cudaFree(d_curand_states_);
    }
    
    // Allocate
    cudaMalloc((void**)&d_curand_states_, n_states * sizeof(curandState));
    n_curand_states_ = n_states;
    
    // Initialize states
    constexpr int THREADS_PER_BLOCK = 256;
    int blocks = (n_states + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    
    cudaStream_t stream = (cudaStream_t)cuda_stream_;
    
    init_curand_states<<<blocks, THREADS_PER_BLOCK, 0, stream>>>(
        d_curand_states_, rng_seed_, n_states
    );
    
    cudaStreamSynchronize(stream);
    curand_initialized_ = true;
}

void GPUCollisionHelper::upload_geometry_to_gpu() {
    if (!geometry_map_host_ || geometry_uploaded_) {
        return;
    }
    
    // Allocate GeometryData_GPU struct
    if (!geometry_gpu_) {
        geometry_gpu_ = new GeometryData_GPU();
    }
    
    const auto& geom_map = *geometry_map_host_;
    int n_species = geom_map.size();
    
    if (n_species == 0) {
        throw std::runtime_error("GPUCollisionHelper: Empty geometry map");
    }
    
    // Create deterministic ordering of species to ensure index consistency
    std::vector<std::string> species_ids;
    species_ids.reserve(n_species);
    for (const auto& kv : geom_map) {
        species_ids.push_back(kv.first);
    }
    std::sort(species_ids.begin(), species_ids.end());

    species_index_map_.clear();
    species_index_map_.reserve(n_species);

    // Count total atoms
    size_t total_atoms = 0;
    for (const auto& species_id : species_ids) {
        const auto& geom_data = geom_map.at(species_id);
        total_atoms += geom_data.first.size();
    }
    
    // Allocate host buffers
    std::vector<double> atom_x_host, atom_y_host, atom_z_host, atom_radii_host;
    std::vector<int> atom_counts_host, atom_offsets_host;
    
    atom_x_host.reserve(total_atoms);
    atom_y_host.reserve(total_atoms);
    atom_z_host.reserve(total_atoms);
    atom_radii_host.reserve(total_atoms);
    atom_counts_host.reserve(n_species);
    atom_offsets_host.reserve(n_species);
    
    // Flatten geometry data
    size_t offset = 0;
    int species_idx = 0;
    for (const auto& species_id : species_ids) {
        const auto& geom_data = geom_map.at(species_id);
        const auto& positions = geom_data.first;
        const auto& radii = geom_data.second;
        
        atom_offsets_host.push_back(offset);
        atom_counts_host.push_back(positions.size());
        species_index_map_.emplace(species_id, species_idx++);
        
        for (const auto& pos : positions) {
            atom_x_host.push_back(pos.x);
            atom_y_host.push_back(pos.y);
            atom_z_host.push_back(pos.z);
        }
        
        for (double r : radii) {
            atom_radii_host.push_back(r);
        }
        
        offset += positions.size();
    }
    
    // Allocate GPU memory
    cudaMalloc(&geometry_gpu_->atom_x, total_atoms * sizeof(double));
    cudaMalloc(&geometry_gpu_->atom_y, total_atoms * sizeof(double));
    cudaMalloc(&geometry_gpu_->atom_z, total_atoms * sizeof(double));
    cudaMalloc(&geometry_gpu_->atom_radii, total_atoms * sizeof(double));
    cudaMalloc(&geometry_gpu_->atom_counts, n_species * sizeof(int));
    cudaMalloc(&geometry_gpu_->atom_offsets, n_species * sizeof(int));
    
    // Upload to GPU
    cudaMemcpy(geometry_gpu_->atom_x, atom_x_host.data(), 
               total_atoms * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(geometry_gpu_->atom_y, atom_y_host.data(), 
               total_atoms * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(geometry_gpu_->atom_z, atom_z_host.data(), 
               total_atoms * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(geometry_gpu_->atom_radii, atom_radii_host.data(), 
               total_atoms * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(geometry_gpu_->atom_counts, atom_counts_host.data(), 
               n_species * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(geometry_gpu_->atom_offsets, atom_offsets_host.data(), 
               n_species * sizeof(int), cudaMemcpyHostToDevice);
    
    geometry_gpu_->num_species = n_species;
    geometry_uploaded_ = true;
}

// Static helper function (internal to .cu, no longer in header)
static EnvironmentParams_GPU convert_environment_params(
    const ICARION::config::EnvironmentConfig& env
) {
    EnvironmentParams_GPU params{};
    params.temperature_K = env.temperature_K;
    params.pressure_Pa = env.pressure_Pa;
    params.gas_velocity_x = env.gas_velocity_m_s.x;
    params.gas_velocity_y = env.gas_velocity_m_s.y;
    params.gas_velocity_z = env.gas_velocity_m_s.z;

    if (!env.gas_mixture.empty()) {
        double density_sum = 0.0;
        double mass_sum = 0.0;
        double radius_sum = 0.0;
        int count = 0;

        for (const auto& gas : env.gas_mixture) {
            if (!gas.participates_in_collisions || gas.density_m3 <= 0.0) {
                continue;
            }

            density_sum += gas.density_m3;
            mass_sum += gas.mass_kg * gas.density_m3;
            radius_sum += gas.radius_m * gas.density_m3;

            if (count < MAX_GPU_GAS_COMPONENTS) {
                params.component_density_m3[count] = gas.density_m3;
                params.component_mass_kg[count] = gas.mass_kg;
                params.component_radius_m[count] = gas.radius_m;
                params.component_cross_section_m2[count] = gas.cross_section_m2;
                ++count;
            }
        }

        params.num_components = count;

        if (density_sum > 0.0) {
            params.neutral_mass_kg = mass_sum / density_sum;
            params.neutral_radius_m = radius_sum / density_sum;
        }
    }

    if (params.neutral_mass_kg <= 0.0) {
        if (env.gas_mass_kg > 0.0) {
            params.neutral_mass_kg = env.gas_mass_kg;
            params.neutral_radius_m = env.gas_radius_m;
        } else {
            params.neutral_mass_kg = 28.0 * AMU_TO_KG;
            params.neutral_radius_m = 1.85e-10;
        }
    }

    params.thermal_velocity_m_s = std::sqrt(
        8.0 * BOLTZMANN_CONSTANT * env.temperature_K / (M_PI * params.neutral_mass_kg)
    );

    double number_density = env.pressure_Pa / (BOLTZMANN_CONSTANT * env.temperature_K);
    double sigma_gas = M_PI * params.neutral_radius_m * params.neutral_radius_m;
    params.mean_free_path_m = 1.0 / (1.414 * sigma_gas * number_density);

    return params;
}

void GPUCollisionHelper::set_geometry(const GeometryMap& geometry_map) {
    geometry_map_host_ = &geometry_map;
    
    if (collision_model_ == "EHSS") {
        upload_geometry_to_gpu();
    }
}

bool GPUCollisionHelper::process_collisions_batch(
    std::vector<IonState>& ions,
    double dt,
    const ICARION::config::EnvironmentConfig& env
) {
    size_t n_ions = ions.size();

    if (!warned_mixture_limit_ &&
        !env.gas_mixture.empty() &&
        env.gas_mixture.size() > MAX_GPU_GAS_COMPONENTS) {
#ifndef __CUDACC__
        ICARION::log::Logger::main()->warn(
            "GPUCollisionHandler: Truncating gas mixture to {} components ({} configured)",
            MAX_GPU_GAS_COMPONENTS,
            env.gas_mixture.size());
#endif
        warned_mixture_limit_ = true;
    }
    
    // Check threshold
    if (n_ions < threshold_) {
        return false;  // Too few ions, use CPU
    }
    
    // Initialize cuRAND states if needed
    if (!curand_initialized_) {
        // Allocate enough states for maximum expected ions
        size_t n_states = MAX(n_ions, 10000);
        initialize_curand_states(n_states);
    }
    
    // ====== ALLOCATE/RESIZE PERSISTENT BUFFERS ======
    // Only allocate once, reuse across timesteps
    if (buffer_capacity_ < n_ions) {
        // Need larger buffers - reallocate with 1.5x growth factor
        size_t new_capacity = n_ions + n_ions / 2;
        
        // Free old buffers
        if (d_vx_) cudaFree(d_vx_);
        if (d_vy_) cudaFree(d_vy_);
        if (d_vz_) cudaFree(d_vz_);
        if (d_mass_) cudaFree(d_mass_);
        if (d_ccs_) cudaFree(d_ccs_);
        if (d_active_) cudaFree(d_active_);
        if (d_species_indices_) cudaFree(d_species_indices_);
        
        // Allocate new buffers
        cudaMalloc(&d_vx_, new_capacity * sizeof(double));
        cudaMalloc(&d_vy_, new_capacity * sizeof(double));
        cudaMalloc(&d_vz_, new_capacity * sizeof(double));
        cudaMalloc(&d_mass_, new_capacity * sizeof(double));
        cudaMalloc(&d_ccs_, new_capacity * sizeof(double));
        cudaMalloc(&d_active_, new_capacity * sizeof(uint8_t));
        
        if (collision_model_ == "EHSS") {
            cudaMalloc(&d_species_indices_, new_capacity * sizeof(int));
        }
        
        buffer_capacity_ = new_capacity;
    }
    
    // Note: No try/catch in CUDA code - nvcc has issues with C++ exception handling
    // ====== UPLOAD PHASE ======
        
        // Flatten ion data (CPU -> GPU)
        std::vector<double> vx_host(n_ions), vy_host(n_ions), vz_host(n_ions);
        std::vector<double> mass_host(n_ions), ccs_host(n_ions);
        std::vector<uint8_t> active_host(n_ions);  // Use uint8_t instead of bool
        std::vector<int> species_indices_host(n_ions);
        
        for (size_t i = 0; i < n_ions; ++i) {
            vx_host[i] = ions[i].vel.x;
            vy_host[i] = ions[i].vel.y;
            vz_host[i] = ions[i].vel.z;
            mass_host[i] = ions[i].mass_kg;
            ccs_host[i] = ions[i].CCS_m2;
            active_host[i] = ions[i].active ? 1 : 0;
            if (collision_model_ == "EHSS") {
                auto it = species_index_map_.find(ions[i].species_id);
                if (it != species_index_map_.end()) {
                    species_indices_host[i] = it->second;
                } else {
#ifndef __CUDACC__
                    if (!missing_species_warned_.count(ions[i].species_id)) {
                        ICARION::log::Logger::main()->warn(
                            "GPUCollisionHelper: No geometry uploaded for species '{}', falling back to index 0",
                            ions[i].species_id);
                        missing_species_warned_.insert(ions[i].species_id);
                    }
#endif
                    species_indices_host[i] = 0;
                }
            }
        }
        
        cudaStream_t stream = (cudaStream_t)cuda_stream_;
        
        // Upload to persistent buffers
        cudaMemcpyAsync(d_vx_, vx_host.data(), n_ions * sizeof(double), 
                        cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_vy_, vy_host.data(), n_ions * sizeof(double), 
                        cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_vz_, vz_host.data(), n_ions * sizeof(double), 
                        cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_mass_, mass_host.data(), n_ions * sizeof(double), 
                        cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_ccs_, ccs_host.data(), n_ions * sizeof(double), 
                        cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_active_, active_host.data(), n_ions * sizeof(uint8_t), 
                        cudaMemcpyHostToDevice, stream);
        
        if (collision_model_ == "EHSS") {
            cudaMemcpyAsync(d_species_indices_, species_indices_host.data(), 
                            n_ions * sizeof(int), 
                            cudaMemcpyHostToDevice, stream);
        }
        
        // ====== COMPUTE PHASE ======
        
        EnvironmentParams_GPU env_gpu = convert_environment_params(env);
        
        if (collision_model_ == "HSS") {
            launch_hss_collision_batch(
                d_vx_, d_vy_, d_vz_,
                d_mass_, d_ccs_, d_active_,
                d_curand_states_,
                env_gpu, dt, n_ions,
                stream
            );
        } else {  // EHSS
            if (!geometry_uploaded_) {
                // Fallback to HSS if geometry not uploaded
                launch_hss_collision_batch(
                    d_vx_, d_vy_, d_vz_,
                    d_mass_, d_ccs_, d_active_,
                    d_curand_states_,
                    env_gpu, dt, n_ions,
                    stream
                );
            } else {
                launch_ehss_collision_batch(
                    d_vx_, d_vy_, d_vz_,
                    d_mass_, d_ccs_, d_species_indices_, d_active_,
                    d_curand_states_,
                    env_gpu, *geometry_gpu_, dt, n_ions,
                    stream
                );
            }
        }
        
        // ====== DOWNLOAD PHASE ======
        
        // Download modified velocities from persistent buffers
        cudaMemcpyAsync(vx_host.data(), d_vx_, n_ions * sizeof(double), 
                        cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(vy_host.data(), d_vy_, n_ions * sizeof(double), 
                        cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(vz_host.data(), d_vz_, n_ions * sizeof(double), 
                        cudaMemcpyDeviceToHost, stream);
        
        cudaStreamSynchronize(stream);
        
        // Write back to host ions
        for (size_t i = 0; i < n_ions; ++i) {
            ions[i].vel.x = vx_host[i];
            ions[i].vel.y = vy_host[i];
            ions[i].vel.z = vz_host[i];
        }
        
        // No cleanup needed - buffers are persistent!
        
        // Update statistics
        stats_.total_batches++;
        stats_.total_ions_processed += n_ions;
        
        return true;  // Success
}

} // namespace gpu
} // namespace icarion
