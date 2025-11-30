// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once
#include <cmath>

// Shared math helpers for host/device parity.
// Provide consistent constants and thin wrappers so both host and device
// call the same symbols and use the same literal constants.

static constexpr double TWO_PI_CONST = 2.0 * 3.14159265358979323846;

static inline __host__ __device__ double TWO_PI() { return TWO_PI_CONST; }

static inline __host__ __device__ double sin_shared(double x) {
    return sin(x);
}

static inline __host__ __device__ double cos_shared(double x) {
    return cos(x);
}

static inline __host__ __device__ double exp_shared(double x) {
    return exp(x);
}
