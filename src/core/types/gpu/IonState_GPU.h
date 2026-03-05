// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#ifndef ICARION_ION_STATE_GPU_H
#define ICARION_ION_STATE_GPU_H

#include "core/types/IonState.h"
#include <cuda_runtime.h>
#include <vector>

namespace icarion {
namespace gpu {

/**
 * @brief GPU-friendly Structure of Arrays (SoA) for ion state
 * 
 * Converts CPU Array of Structures (AoS) to GPU Structure of Arrays (SoA)
 * for better memory coalescing and GPU performance.
 * 
 * CPU layout (AoS):
 *   ion[0]: {x, y, z, vx, vy, vz, ...}
 *   ion[1]: {x, y, z, vx, vy, vz, ...}
 * 
 * GPU layout (SoA):
 *   x[]:  {ion0_x, ion1_x, ion2_x, ...}
 *   y[]:  {ion0_y, ion1_y, ion2_y, ...}
 *   vx[]: {ion0_vx, ion1_vx, ...}
 */
struct IonStateGPU {
    // Positions (m)
    double* x;
    double* y;
    double* z;
    
    // Velocities (m/s)
    double* vx;
    double* vy;
    double* vz;
    
    // Physical properties
    double* mass;           // kg
    double* charge;         // C
    int* species_id;        // Species identifier
    
    // Simulation state
    int* domain_id;         // Current domain
    bool* active;           // Active flag
    
    size_t count;           // Number of ions
    
    /**
     * @brief Default constructor
     */
    IonStateGPU()
        : x(nullptr), y(nullptr), z(nullptr)
        , vx(nullptr), vy(nullptr), vz(nullptr)
        , mass(nullptr), charge(nullptr), species_id(nullptr)
        , domain_id(nullptr), active(nullptr)
        , count(0)
    {}
    
    /**
     * @brief Allocate device memory for N ions
     * @param N Number of ions
     */
    void allocate(size_t N);
    
    /**
     * @brief Free device memory
     */
    void free();
    
    /**
     * @brief Check if memory is allocated
     */
    bool is_allocated() const { return x != nullptr; }
};

/**
 * @brief Helper functions for IonState <-> IonStateGPU conversion
 */
namespace ion_state_conversion {

/**
 * @brief Upload CPU ion states to GPU (AoS -> SoA)
 * @param ions_cpu CPU ion states (vector of IonState)
 * @param ions_gpu GPU ion states (must be pre-allocated)
 * @param stream CUDA stream for async upload
 */
void upload_ions(
    const std::vector<IonState>& ions_cpu,
    IonStateGPU& ions_gpu,
    cudaStream_t stream = 0
);

/**
 * @brief Download GPU ion states to CPU (SoA -> AoS)
 * @param ions_gpu GPU ion states
 * @param ions_cpu CPU ion states (will be resized)
 * @param stream CUDA stream for async download
 */
void download_ions(
    const IonStateGPU& ions_gpu,
    std::vector<IonState>& ions_cpu,
    cudaStream_t stream = 0
);

/**
 * @brief Upload only positions and velocities (for integration output)
 * @param ions_cpu CPU ion states
 * @param ions_gpu GPU ion states (must be pre-allocated)
 * @param stream CUDA stream for async upload
 */
void upload_positions_velocities(
    const std::vector<IonState>& ions_cpu,
    IonStateGPU& ions_gpu,
    cudaStream_t stream = 0
);

/**
 * @brief Download only positions and velocities (for integration output)
 * @param ions_gpu GPU ion states
 * @param ions_cpu CPU ion states (must have correct size)
 * @param stream CUDA stream for async download
 */
void download_positions_velocities(
    const IonStateGPU& ions_gpu,
    std::vector<IonState>& ions_cpu,
    cudaStream_t stream = 0
);

} // namespace ion_state_conversion

} // namespace gpu
} // namespace icarion

#endif // ICARION_ION_STATE_GPU_H
