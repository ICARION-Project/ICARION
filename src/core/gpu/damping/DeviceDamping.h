// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#pragma once

namespace icarion {
namespace gpu {

struct DeviceDamping {
    int enabled = 0;             // 0 = off, 1 = on
    float nu_const = 0.0f;       // Constant damping rate [1/s]
    const float* nu_per_ion = nullptr;  // Optional per-ion damping rates [1/s] on device
};

}  // namespace gpu
}  // namespace icarion
