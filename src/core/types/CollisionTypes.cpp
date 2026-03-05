// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include "core/types/CollisionTypes.h"

namespace ICARION::physics {

PhysicsRng::PhysicsRng(uint64_t seed)
    : eng_(seed)
    , uni_(0.0, 1.0)
    , norm_(0.0, 1.0)
{}

double PhysicsRng::uniform01() {
    return uni_(eng_);
}

double PhysicsRng::normal() {
    return norm_(eng_);
}

} // namespace ICARION::physics
