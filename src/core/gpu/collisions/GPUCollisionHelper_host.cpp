// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
/**
 * @file GPUCollisionHelper_host.cpp
 * @brief Host-side implementation (non-CUDA)
 * 
 * Separates GPUContext dependencies from CUDA compilation unit.
 * nvcc cannot handle forward declarations properly, so we split:
 * - This file (.cpp): Factory, constructor, destructor (uses GPUContext)
 * - .cu file: Kernel launches and GPU-specific code
 *
 * The helper is experimental: EHSS assumes caller-provided geometry and
 * species indices; there is no validation parity with the CPU path yet.
 */

#include "core/gpu/core/GPUContext.h"
#include "GPUCollisionHelper.h"
#include <cuda_runtime.h>
#include <stdexcept>

namespace icarion {
namespace gpu {

// Helper functions (implemented here, used by factory/constructor)
static void* extract_cuda_stream_impl(const GPUContext& ctx) {
    return (void*)ctx.get_stream();
}

static bool check_cuda_available_impl() {
    return GPUContext::is_cuda_available();
}

std::unique_ptr<GPUCollisionHelper> GPUCollisionHelper::create(
    const gpu::GPUContext& context,
    size_t threshold,
    const std::string& collision_model,
    unsigned long long rng_seed
) {
    if (!check_cuda_available_impl()) {
        return nullptr;
    }
    
    try {
        return std::unique_ptr<GPUCollisionHelper>(
            new GPUCollisionHelper(context, threshold, collision_model, rng_seed)
        );
    } catch (const std::exception&) {
        return nullptr;
    }
}

GPUCollisionHelper::GPUCollisionHelper(
    const gpu::GPUContext& context,
    size_t threshold,
    const std::string& collision_model,
    unsigned long long rng_seed
)
    : cuda_stream_(extract_cuda_stream_impl(context)),
      threshold_(threshold),
      collision_model_(collision_model),
      rng_seed_(rng_seed),
      stats_{}
{
    if (collision_model != "HSS" && collision_model != "EHSS") {
        throw std::invalid_argument(
            "GPUCollisionHelper: Unsupported collision model '" + collision_model +
            "'. Supported: 'HSS', 'EHSS'."
        );
    }
}

// Destructor implemented in .cu file (needs cudaFree)

} // namespace gpu
} // namespace icarion
