// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include "GPUMemoryPool.h"

namespace icarion {
namespace gpu {

GPUMemoryPool::GPUMemoryPool(const GPUContext& context)
    : context_(context)
    , stats_{}
{
}

GPUMemoryPool::~GPUMemoryPool() {
    clear();
}

void GPUMemoryPool::clear() {
    // Free all device buffers
    for (auto& [key, buffers] : device_free_list_) {
        for (void* ptr : buffers) {
            cudaFree(ptr);
        }
        buffers.clear();
    }
    device_free_list_.clear();
    
    // Free all pinned buffers
    for (auto& [key, buffers] : pinned_free_list_) {
        for (void* ptr : buffers) {
            cudaFreeHost(ptr);
        }
        buffers.clear();
    }
    pinned_free_list_.clear();
    
    // Reset stats (keep allocated counts, clear cached)
    stats_.device_bytes_cached = 0;
    stats_.device_buffers_cached = 0;
    stats_.pinned_bytes_cached = 0;
    stats_.pinned_buffers_cached = 0;
}

GPUMemoryPool::Stats GPUMemoryPool::get_stats() const {
    return stats_;
}

} // namespace gpu
} // namespace icarion
