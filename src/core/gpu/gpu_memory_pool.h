/**
 * GPU Memory Pool for ICARION
 * 
 * Reduces GPU malloc/free overhead by pre-allocating memory pools and reusing them
 * across simulation runs. Expected 1.2-1.5x speedup for low effort.
 * 
 * Implements Phase 2 Section 9.2 from REVISED_OPTIMIZATION_STRATEGY.md
 * 
 * Date: 2025-01-23
 */

#pragma once

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <vector>
#include <unordered_map>
#include <memory>
#include <iostream>
#include "core/types/IonState.h"
#include "core/types/IonState_GPU.h"
#include "core/gpu/gpu_geometry.h"

namespace ICARION {
namespace gpu {

/**
 * @brief GPU Memory Pool Manager
 * 
 * Pre-allocates GPU memory pools for common data structures to eliminate
 * repeated malloc/free overhead during simulation runs.
 */
class GPUMemoryPool {
public:
    static GPUMemoryPool& getInstance() {
        static GPUMemoryPool instance;
        return instance;
    }

    /**
     * @brief Get pre-allocated ion state memory
     * @param max_ions Maximum number of ions needed
     * @return Device pointer to ion state array
     */
    IonStateGPU* getIonStateMemory(size_t max_ions);

    /**
     * @brief Get pre-allocated RNG state memory  
     * @param max_ions Maximum number of RNG states needed
     * @return Device pointer to curandState array
     */
    curandState* getRNGStateMemory(size_t max_ions);

    /**
     * @brief Get pre-allocated geometry memory
     * @param total_atoms Total atoms across all species
     * @return Device pointer to multi-species geometry
     */
    MultiSpeciesGeometryGPU* getGeometryMemory(size_t total_atoms);

    /**
     * @brief Get temporary work memory for kernels
     * @param size_bytes Size in bytes needed
     * @return Device pointer to temporary memory
     */
    void* getWorkMemory(size_t size_bytes);

    /**
     * @brief Release all memory back to pool for reuse
     * 
     * Call this after each simulation to make memory available for next run.
     * Does NOT free GPU memory - just marks it as available for reuse.
     */
    void releaseAll();

    /**
     * @brief Get memory usage statistics
     */
    void printStats() const;

    /**
     * @brief Free all GPU memory (called at shutdown)
     */
    ~GPUMemoryPool();

private:
    GPUMemoryPool() = default;
    GPUMemoryPool(const GPUMemoryPool&) = delete;
    GPUMemoryPool& operator=(const GPUMemoryPool&) = delete;

    struct MemoryBlock {
        void* ptr;
        size_t size;
        bool in_use;
        
        MemoryBlock(void* p, size_t s) : ptr(p), size(s), in_use(false) {}
    };

    // Memory pools for different data types
    std::vector<std::unique_ptr<MemoryBlock>> ion_state_pool_;
    std::vector<std::unique_ptr<MemoryBlock>> rng_state_pool_;
    std::vector<std::unique_ptr<MemoryBlock>> geometry_pool_;
    std::vector<std::unique_ptr<MemoryBlock>> work_memory_pool_;

    // Statistics
    mutable size_t total_allocations_ = 0;
    mutable size_t pool_hits_ = 0;
    mutable size_t pool_misses_ = 0;

    /**
     * @brief Find or allocate memory from a specific pool
     */
    template<typename T>
    T* getFromPool(std::vector<std::unique_ptr<MemoryBlock>>& pool, 
                   size_t required_size, const std::string& pool_name);

    /**
     * @brief Allocate new GPU memory block
     */
    void* allocateGPUMemory(size_t size_bytes);
};

} // namespace gpu
} // namespace ICARION