// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#pragma once
#include "core/types/Vec3.h"
#include <vector>
#include <cmath>
#include <ostream>
#include <array>
#include <initializer_list>

// -----------------------------
// Vec3 operations (GPU-compatible)
// -----------------------------

/** @brief Dot product: a·b = ax*bx + ay*by + az*bz */
__host__ __device__ inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/** @brief Cross product: a×b (right-hand rule) */
__host__ __device__ inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

/** @brief Euclidean norm (magnitude): ||a|| = sqrt(a·a) */
__host__ __device__ inline double norm(const Vec3& a) {
    return sqrt(dot(a, a));
}

/** @brief Normalize vector to unit length: a/||a|| (returns zero if norm=0) */
__host__ __device__ inline Vec3 normalize(const Vec3& a) {
    double n = norm(a);
    return (n > 0.0 ? a / n : Vec3{0, 0, 0});
}

/** @brief Matrix-vector multiplication: R × v (3×3 matrix times 3-vector) */
__host__ __device__ inline Vec3 Rmul(const double R[3][3], const Vec3& v) {
    return {
        R[0][0]*v.x + R[0][1]*v.y + R[0][2]*v.z,
        R[1][0]*v.x + R[1][1]*v.y + R[1][2]*v.z,
        R[2][0]*v.x + R[2][1]*v.y + R[2][2]*v.z
    };
}

// ==============================
// 3x3 matrix
// ==============================

/**
 * @struct Mat3
 * @brief Simple 3×3 matrix.
 *
 * Supports identity, construction from columns, transpose, and multiplication with Vec3.
 */
struct Mat3 {
    double m[3][3];

    Mat3() {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                m[i][j] = (i == j) ? 1.0 : 0.0;
    }

    static Mat3 identity() {
        return Mat3();
    }

    static Mat3 fromColumns(const Vec3& c1, const Vec3& c2, const Vec3& c3) {
        Mat3 R;
        R.m[0][0] = c1.x; R.m[1][0] = c1.y; R.m[2][0] = c1.z;
        R.m[0][1] = c2.x; R.m[1][1] = c2.y; R.m[2][1] = c2.z;
        R.m[0][2] = c3.x; R.m[1][2] = c3.y; R.m[2][2] = c3.z;
        return R;
    }

};

inline Mat3 transpose(const Mat3& M) {
    Mat3 T;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            T.m[i][j] = M.m[j][i];
    return T;
}

// Multiplication with vector
inline Vec3 operator*(const Mat3& R, const Vec3& v) {
    return {
        R.m[0][0]*v.x + R.m[0][1]*v.y + R.m[0][2]*v.z,
        R.m[1][0]*v.x + R.m[1][1]*v.y + R.m[1][2]*v.z,
        R.m[2][0]*v.x + R.m[2][1]*v.y + R.m[2][2]*v.z
    };
}

// Allow scalar * Vec3 in addition to Vec3 * scalar
inline Vec3 operator*(double s, const Vec3& v) {
    return v * s;
}

// Streaming for debug/logging
inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    os << "(" << v.x << "," << v.y << "," << v.z << ")";
    return os;
}

inline Mat3 rotationFromDirection(const Vec3& dir) {
    Vec3 z = normalize(dir);
    Vec3 x = normalize(cross({0,1,0}, z));
    if (norm(x) < 1e-6) x = normalize(cross({1,0,0}, z));
    Vec3 y = cross(z, x);
    return Mat3::fromColumns(x, y, z);
}

// ==============================
// std::vector<double> arithmetic
// ==============================
std::vector<double> add_vectors(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> sub_vectors(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> scale_vector(const std::vector<double>& v, double scalar);
double              get_vector_magnitude(const std::vector<double>& y);
