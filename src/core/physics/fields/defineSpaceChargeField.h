// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once
#include "core/types/IonState.h"
#include <vector>

// Compute Coulombic space-charge field acting on one ion
Vec3 SpaceChargeField(const IonState& ion,
                      const std::vector<IonState>& ions,
                      double eps0);
