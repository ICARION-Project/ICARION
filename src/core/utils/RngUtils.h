// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include "core/types/CollisionTypes.h"
#include "core/types/IonEnsemble.h"
#include <cstddef>
#include <cstdint>
#include <vector>

namespace ICARION {
namespace utils {

uint64_t splitmix64(uint64_t x);
uint64_t double_bits(double value);
uint64_t mix_seed(uint64_t seed, uint64_t value);

uint64_t ion_rng_seed(uint64_t base_seed, size_t ion_index);
std::vector<uint64_t> build_ion_rng_fingerprints(const core::IonEnsemble& ensemble);

void sync_rng_pool_for_ensemble(std::vector<physics::PhysicsRng>& rng_by_ion,
                                const std::vector<uint64_t>& previous_fingerprints,
                                const core::IonEnsemble& ensemble,
                                uint64_t base_seed);

}  // namespace utils
}  // namespace ICARION
