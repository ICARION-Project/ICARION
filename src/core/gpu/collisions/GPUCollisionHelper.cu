// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file GPUCollisionHelper.cu
 * @brief GPU collision helper implementation
 */

/**
 * @file GPUCollisionHelper.cu
 * @brief CUDA implementation (kernel launches and device code only)
 * 
 * Factory/constructor/destructor moved to _host.cpp to avoid nvcc issues with GPUContext.
 * This file contains only GPU kernel launches and device memory management.
 */

#include "GPUCollisionHelper.h"
#include "collision_kernels_gpu.cuh"
#include "utils/constants.h"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdint>

// Minimal C++ utilities (avoid STL headers with nvcc)
#define MAX(a, b) ((a) > (b) ? (a) : (b))

namespace icarion {
namespace gpu {

GPUCollisionHelper::~GPUCollisionHelper() {
    if (d_curand_states_) {
        cudaFree(d_curand_states_);
    }
    
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
    
    // Count total atoms
    size_t total_atoms = 0;
    for (const auto& [species_id, geom_data] : geom_map) {
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
    for (const auto& [species_id, geom_data] : geom_map) {
        const auto& positions = geom_data.first;
        const auto& radii = geom_data.second;
        
        atom_offsets_host.push_back(offset);
        atom_counts_host.push_back(positions.size());
        
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
    EnvironmentParams_GPU params;
    params.temperature_K = env.temperature_K;
    params.pressure_Pa = env.pressure_Pa;
    
    // Use first gas component (TODO: support mixtures)
    if (!env.gas_mixture.empty()) {
        // Compute weighted average mass
        double total_mass = 0.0;
        double total_fraction = 0.0;
        for (const auto& gas : env.gas_mixture) {
            total_mass += gas.mass_kg * gas.mole_fraction;
            total_fraction += gas.mole_fraction;
        }
        params.neutral_mass_kg = total_mass / total_fraction;
        
        // Use first gas radius (simplification)
        params.neutral_radius_m = env.gas_mixture[0].radius_m;
    } else {
        // Fallback to environment-level properties
        params.neutral_mass_kg = 28.0 * AMU_TO_KG;  // N2 default
        params.neutral_radius_m = 1.85e-10;  // N2 default
    }
    
    params.gas_velocity_x = env.gas_velocity_m_s.x;
    params.gas_velocity_y = env.gas_velocity_m_s.y;
    params.gas_velocity_z = env.gas_velocity_m_s.z;
    
    // Compute derived quantities
    params.thermal_velocity_m_s = std::sqrt(
        8.0 * BOLTZMANN_CONSTANT * env.temperature_K / (M_PI * params.neutral_mass_kg)
    );
    
    // Mean free path (simplified)
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
    
    // Note: No try/catch in CUDA code - nvcc has issues with C++ exception handling
    // ====== UPLOAD PHASE ======
        
        // Allocate device memory (could use GPUMemoryPool here)
        double *d_vx, *d_vy, *d_vz;
        double *d_mass, *d_ccs;
        uint8_t *d_active;
        int *d_species_indices = nullptr;
        
        cudaMalloc(&d_vx, n_ions * sizeof(double));
        cudaMalloc(&d_vy, n_ions * sizeof(double));
        cudaMalloc(&d_vz, n_ions * sizeof(double));
        cudaMalloc(&d_mass, n_ions * sizeof(double));
        cudaMalloc(&d_ccs, n_ions * sizeof(double));
        cudaMalloc(&d_active, n_ions * sizeof(uint8_t));
        
        if (collision_model_ == "EHSS") {
            cudaMalloc(&d_species_indices, n_ions * sizeof(int));
        }
        
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
                // TODO: Map species_id to species_index using a species database lookup
                species_indices_host[i] = 0;  // For now, assume single species
            }
        }
        
        cudaStream_t stream = (cudaStream_t)cuda_stream_;
        
        // Upload
        cudaMemcpyAsync(d_vx, vx_host.data(), n_ions * sizeof(double), 
                        cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_vy, vy_host.data(), n_ions * sizeof(double), 
                        cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_vz, vz_host.data(), n_ions * sizeof(double), 
                        cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_mass, mass_host.data(), n_ions * sizeof(double), 
                        cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_ccs, ccs_host.data(), n_ions * sizeof(double), 
                        cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_active, active_host.data(), n_ions * sizeof(uint8_t), 
                        cudaMemcpyHostToDevice, stream);
        
        if (collision_model_ == "EHSS") {
            cudaMemcpyAsync(d_species_indices, species_indices_host.data(), 
                            n_ions * sizeof(int), 
                            cudaMemcpyHostToDevice, stream);
        }
        
        // ====== COMPUTE PHASE ======
        
        EnvironmentParams_GPU env_gpu = convert_environment_params(env);
        
        if (collision_model_ == "HSS") {
            launch_hss_collision_batch(
                d_vx, d_vy, d_vz,
                d_mass, d_ccs, d_active,
                d_curand_states_,
                env_gpu, dt, n_ions,
                stream
            );
        } else {  // EHSS
            if (!geometry_uploaded_) {
                // Fallback to HSS if geometry not uploaded
                launch_hss_collision_batch(
                    d_vx, d_vy, d_vz,
                    d_mass, d_ccs, d_active,
                    d_curand_states_,
                    env_gpu, dt, n_ions,
                    stream
                );
            } else {
                launch_ehss_collision_batch(
                    d_vx, d_vy, d_vz,
                    d_mass, d_ccs, d_species_indices, d_active,
                    d_curand_states_,
                    env_gpu, *geometry_gpu_, dt, n_ions,
                    stream
                );
            }
        }
        
        // ====== DOWNLOAD PHASE ======
        
        // Download modified velocities
        cudaMemcpyAsync(vx_host.data(), d_vx, n_ions * sizeof(double), 
                        cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(vy_host.data(), d_vy, n_ions * sizeof(double), 
                        cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(vz_host.data(), d_vz, n_ions * sizeof(double), 
                        cudaMemcpyDeviceToHost, stream);
        
        cudaStreamSynchronize(stream);
        
        // Write back to host ions
        for (size_t i = 0; i < n_ions; ++i) {
            ions[i].vel.x = vx_host[i];
            ions[i].vel.y = vy_host[i];
            ions[i].vel.z = vz_host[i];
        }
        
        // Cleanup device memory
        cudaFree(d_vx);
        cudaFree(d_vy);
        cudaFree(d_vz);
        cudaFree(d_mass);
        cudaFree(d_ccs);
        cudaFree(d_active);
        if (d_species_indices) cudaFree(d_species_indices);
        
        // Update statistics
        stats_.total_batches++;
        stats_.total_ions_processed += n_ions;
        
        return true;  // Success
}

} // namespace gpu
} // namespace icarion
