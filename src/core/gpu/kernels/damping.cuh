// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include "core/types/gpu/IonState_GPU.h"
#include "core/gpu/damping/DeviceDamping.h"
#include <cuda_runtime.h>

namespace icarion {
namespace gpu {

// Apply linear velocity damping: v <- v * exp(-nu * dt)
// Supports constant nu or per-ion nu (nu_per_ion). No noise/CCS/Teff.
void launch_damping_kernel(
    const IonStateGPU& ions,
    const DeviceDamping& damping,
    double dt,
    cudaStream_t stream = 0
);

}  // namespace gpu
}  // namespace icarion
