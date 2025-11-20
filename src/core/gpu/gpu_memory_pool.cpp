/**
 * GPU Memory Pool Implementation
 * 
 * Date: 2025-01-23
 */

#include "gpu_memory_pool.h"
#include <algorithm>
#include <iomanip>

namespace ICARION {
namespace gpu {

IonStateGPU* GPUMemoryPool::getIonStateMemory(size_t max_ions) {
    size_t required_bytes = max_ions * sizeof(IonStateGPU);
    return getFromPool<IonStateGPU>(ion_state_pool_, required_bytes, "IonState");
}

curandState* GPUMemoryPool::getRNGStateMemory(size_t max_ions) {
    size_t required_bytes = max_ions * sizeof(curandState);
    return getFromPool<curandState>(rng_state_pool_, required_bytes, "RNGState");
}

MultiSpeciesGeometryGPU* GPUMemoryPool::getGeometryMemory(size_t total_atoms) {
    // Estimate size needed for geometry structure
    size_t required_bytes = sizeof(MultiSpeciesGeometryGPU) + 
                           total_atoms * (sizeof(Vec3) + sizeof(double)) + 
                           1024; // Extra buffer for indices/metadata
    return getFromPool<MultiSpeciesGeometryGPU>(geometry_pool_, required_bytes, "Geometry");
}

void* GPUMemoryPool::getWorkMemory(size_t size_bytes) {
    return getFromPool<void>(work_memory_pool_, size_bytes, "WorkMemory");
}

template<typename T>
T* GPUMemoryPool::getFromPool(std::vector<std::unique_ptr<MemoryBlock>>& pool, 
                              size_t required_size, const std::string& pool_name) {
    total_allocations_++;
    
    // Look for existing block that's big enough and not in use
    for (auto& block : pool) {
        if (!block->in_use && block->size >= required_size) {
            block->in_use = true;
            pool_hits_++;
            std::cout << "[GPU Memory Pool] Reusing " << pool_name 
                      << " block: " << block->size << " bytes" << std::endl;
            return static_cast<T*>(block->ptr);
        }
    }
    
    // No suitable block found - allocate new one
    pool_misses_++;
    void* ptr = allocateGPUMemory(required_size);
    
    auto new_block = std::make_unique<MemoryBlock>(ptr, required_size);
    new_block->in_use = true;
    
    T* result = static_cast<T*>(ptr);
    pool.push_back(std::move(new_block));
    
    std::cout << "[GPU Memory Pool] Allocated new " << pool_name 
              << " block: " << required_size << " bytes" << std::endl;
    
    return result;
}

void* GPUMemoryPool::allocateGPUMemory(size_t size_bytes) {
    void* ptr = nullptr;
    cudaError_t err = cudaMalloc(&ptr, size_bytes);
    
    if (err != cudaSuccess) {
        std::cerr << "GPU Memory Pool: cudaMalloc failed for " << size_bytes 
                  << " bytes: " << cudaGetErrorString(err) << std::endl;
        throw std::runtime_error("GPU memory allocation failed");
    }
    
    return ptr;
}

void GPUMemoryPool::releaseAll() {
    // Mark all blocks as not in use (but keep the memory allocated)
    auto release_pool = [](std::vector<std::unique_ptr<MemoryBlock>>& pool) {
        for (auto& block : pool) {
            block->in_use = false;
        }
    };
    
    release_pool(ion_state_pool_);
    release_pool(rng_state_pool_);
    release_pool(geometry_pool_);
    release_pool(work_memory_pool_);
    
    std::cout << "[GPU Memory Pool] Released all memory blocks for reuse" << std::endl;
}

void GPUMemoryPool::printStats() const {
    auto pool_info = [](const std::vector<std::unique_ptr<MemoryBlock>>& pool, 
                        const std::string& name) {
        size_t total_size = 0;
        size_t in_use_count = 0;
        
        for (const auto& block : pool) {
            total_size += block->size;
            if (block->in_use) in_use_count++;
        }
        
        std::cout << "  " << std::setw(12) << name << ": " 
                  << pool.size() << " blocks, "
                  << total_size / (1024.0 * 1024.0) << " MB, "
                  << in_use_count << " in use" << std::endl;
    };
    
    std::cout << "\n=== GPU Memory Pool Statistics ===" << std::endl;
    std::cout << "Total allocations: " << total_allocations_ << std::endl;
    std::cout << "Pool hits: " << pool_hits_ 
              << " (" << (100.0 * pool_hits_ / total_allocations_) << "%)" << std::endl;
    std::cout << "Pool misses: " << pool_misses_ 
              << " (" << (100.0 * pool_misses_ / total_allocations_) << "%)" << std::endl;
    std::cout << "Memory pools:" << std::endl;
    
    pool_info(ion_state_pool_, "IonState");
    pool_info(rng_state_pool_, "RNGState");
    pool_info(geometry_pool_, "Geometry");
    pool_info(work_memory_pool_, "WorkMemory");
    
    std::cout << "===================================" << std::endl;
}

GPUMemoryPool::~GPUMemoryPool() {
    std::cout << "[GPU Memory Pool] Cleaning up..." << std::endl;
    
    auto cleanup_pool = [](std::vector<std::unique_ptr<MemoryBlock>>& pool) {
        for (auto& block : pool) {
            if (block->ptr) {
                cudaFree(block->ptr);
                block->ptr = nullptr;
            }
        }
        pool.clear();
    };
    
    cleanup_pool(ion_state_pool_);
    cleanup_pool(rng_state_pool_);
    cleanup_pool(geometry_pool_);
    cleanup_pool(work_memory_pool_);
    
    std::cout << "[GPU Memory Pool] Cleanup complete" << std::endl;
}

} // namespace gpu
} // namespace ICARION