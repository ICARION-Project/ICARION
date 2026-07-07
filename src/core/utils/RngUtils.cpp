// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "core/utils/RngUtils.h"
#include <cstring>
#include <limits>
#include <unordered_map>

namespace ICARION {
namespace utils {
namespace {

uint64_t hash_double_bits(double value) {
    return splitmix64(double_bits(value));
}

uint64_t ion_rng_fingerprint(const core::IonEnsemble& ensemble, size_t idx) {
    uint64_t h = 0x6a09e667f3bcc909ULL;
    h = splitmix64(h ^ hash_double_bits(ensemble.pos_x_data()[idx]));
    h = splitmix64(h ^ hash_double_bits(ensemble.pos_y_data()[idx]));
    h = splitmix64(h ^ hash_double_bits(ensemble.pos_z_data()[idx]));
    h = splitmix64(h ^ hash_double_bits(ensemble.vel_x_data()[idx]));
    h = splitmix64(h ^ hash_double_bits(ensemble.vel_y_data()[idx]));
    h = splitmix64(h ^ hash_double_bits(ensemble.vel_z_data()[idx]));
    h = splitmix64(h ^ hash_double_bits(ensemble.mass_data()[idx]));
    h = splitmix64(h ^ hash_double_bits(ensemble.charge_data()[idx]));
    h = splitmix64(h ^ static_cast<uint64_t>(ensemble.species_id_indices()[idx]));
    return h;
}

}  // namespace

uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

uint64_t double_bits(double value) {
    uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value), "Unexpected double size");
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

uint64_t mix_seed(uint64_t seed, uint64_t value) {
    return splitmix64(seed ^ (value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2)));
}

uint64_t ion_rng_seed(uint64_t base_seed, size_t ion_index) {
    return base_seed + static_cast<uint64_t>(ion_index);
}

std::vector<uint64_t> build_ion_rng_fingerprints(const core::IonEnsemble& ensemble) {
    std::vector<uint64_t> out(ensemble.size());
    for (size_t i = 0; i < ensemble.size(); ++i) {
        out[i] = ion_rng_fingerprint(ensemble, i);
    }
    return out;
}

void sync_rng_pool_for_ensemble(std::vector<physics::PhysicsRng>& rng_by_ion,
                                const std::vector<uint64_t>& previous_fingerprints,
                                const core::IonEnsemble& ensemble,
                                uint64_t base_seed) {
    const size_t n_ions = ensemble.size();
    if (rng_by_ion.size() == n_ions) {
        return;
    }

    if (rng_by_ion.empty() || previous_fingerprints.empty()) {
        rng_by_ion.clear();
        rng_by_ion.reserve(n_ions);
        for (size_t i = 0; i < n_ions; ++i) {
            rng_by_ion.emplace_back(ion_rng_seed(base_seed, i));
        }
        return;
    }

    const auto current_fingerprints = build_ion_rng_fingerprints(ensemble);
    std::unordered_map<uint64_t, std::vector<size_t>> old_slots;
    old_slots.reserve(previous_fingerprints.size() * 2 + 1);
    for (size_t i = 0; i < previous_fingerprints.size(); ++i) {
        old_slots[previous_fingerprints[i]].push_back(i);
    }

    std::vector<uint8_t> old_used(rng_by_ion.size(), 0);
    std::vector<physics::PhysicsRng> remapped;
    remapped.reserve(n_ions);
    for (size_t i = 0; i < n_ions; ++i) {
        size_t match = std::numeric_limits<size_t>::max();
        auto slot_it = old_slots.find(current_fingerprints[i]);
        if (slot_it != old_slots.end()) {
            auto& slots = slot_it->second;
            while (!slots.empty()) {
                const size_t candidate = slots.back();
                slots.pop_back();
                if (candidate < rng_by_ion.size() && !old_used[candidate]) {
                    match = candidate;
                    break;
                }
            }
        }
        if (match == std::numeric_limits<size_t>::max() &&
            i < rng_by_ion.size() &&
            !old_used[i]) {
            match = i;
        }
        if (match != std::numeric_limits<size_t>::max()) {
            old_used[match] = 1;
            remapped.push_back(std::move(rng_by_ion[match]));
        } else {
            remapped.emplace_back(ion_rng_seed(base_seed, i));
        }
    }
    rng_by_ion = std::move(remapped);
}

}  // namespace utils
}  // namespace ICARION
