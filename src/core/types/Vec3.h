// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#pragma once
#include <cmath>
// --------------------------------------------------------------------
// CUDA host/device annotation helper
// Avoid redefining CUDA built-in macros like __host__/__device__/__global__
// which can trigger redefinition warnings when CUDA headers are present.
// Instead provide a local macro `HD` that expands to the CUDA annotations
// when compiling with NVCC and to nothing for a regular host-only build.
// --------------------------------------------------------------------
#if defined(__CUDACC__)
#define HD __host__ __device__
#else
#define HD
#endif

// Provide empty definitions for CUDA annotation macros if they are missing
// on non-CUDA host compilers. Use individual guards to avoid redefining
// them when CUDA headers are present (prevents redefinition warnings).
#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif
#ifndef __global__
#define __global__
#endif

namespace ICARION {
namespace core {

/**
 * @struct Vec3
 * @brief Simple 3D vector with basic arithmetic operations.
 *
 * Header-only for inlining; CUDA annotations are gated via `HD` when available.
 */
struct Vec3 {
    double x, y, z;

  HD
  Vec3() : x(0), y(0), z(0) {}

  HD
  Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}

    // =============================
    // Vec3 arithmetic
    // =============================

  HD
  Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }

  HD
  Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }

  HD
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }

  HD
  Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }

  HD
  Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }

  HD
  Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
};

}  // namespace core
}  // namespace ICARION

// Bring Vec3 into global namespace for backward compatibility
using ICARION::core::Vec3;
