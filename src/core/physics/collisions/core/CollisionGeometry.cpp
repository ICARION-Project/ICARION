// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "CollisionGeometry.h"
#include "core/utils/mathUtils.h"
#include <cmath>
#include <algorithm>

namespace ICARION::physics::collision_core {

void CollisionGeometry::construct_orthonormal_basis(
    const Vec3& ehat,
    Vec3& t1,
    Vec3& t2
) {
    // PHYSICS: Exact implementation from collisionHelpers.cpp::ortho_basis()
    // NO CHANGES to numerical algorithm - only renamed for clarity
    
    // Pick a vector not parallel to ehat
    // Use x-axis unless ehat is nearly parallel to it
    Vec3 a = (std::fabs(ehat.x) < 0.9) ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    
    // Gram-Schmidt: first tangent perpendicular to ehat
    t1 = normalize(cross(ehat, a));
    
    // Second tangent perpendicular to both ehat and t1
    t2 = normalize(cross(ehat, t1));
}

void CollisionGeometry::generate_random_rotation(
    EhssRng& rng,
    double R[3][3]
) {
    // PHYSICS: Exact implementation from collisionHelpers.cpp::rand_rotation()
    // NO CHANGES to quaternion algorithm - only renamed for clarity
    
    // Uniform unit quaternion sampling (Shoemake 1992)
    double u1 = rng.uniform01();
    double u2 = rng.uniform01();
    double u3 = rng.uniform01();

    double q1 = std::sqrt(1.0 - u1) * std::sin(2.0 * M_PI * u2);
    double q2 = std::sqrt(1.0 - u1) * std::cos(2.0 * M_PI * u2);
    double q3 = std::sqrt(u1) * std::sin(2.0 * M_PI * u3);
    double q4 = std::sqrt(u1) * std::cos(2.0 * M_PI * u3);

    double x = q1, y = q2, z = q3, w = q4;
    
    // Convert quaternion (x,y,z,w) to rotation matrix
    R[0][0] = 1 - 2 * (y * y + z * z);
    R[0][1] = 2 * (x * y - z * w);
    R[0][2] = 2 * (x * z + y * w);
    R[1][0] = 2 * (x * y + z * w);
    R[1][1] = 1 - 2 * (x * x + z * z);
    R[1][2] = 2 * (y * z - x * w);
    R[2][0] = 2 * (x * z - y * w);
    R[2][1] = 2 * (y * z + x * w);
    R[2][2] = 1 - 2 * (x * x + y * y);
}

Vec3 CollisionGeometry::rotate_vector(
    const Vec3& v,
    const double R[3][3]
) {
    // Matrix-vector multiplication: R·v
    return Vec3{
        R[0][0]*v.x + R[0][1]*v.y + R[0][2]*v.z,
        R[1][0]*v.x + R[1][1]*v.y + R[1][2]*v.z,
        R[2][0]*v.x + R[2][1]*v.y + R[2][2]*v.z
    };
}

bool CollisionGeometry::is_valid_rotation(
    const double R[3][3],
    double tolerance
) {
    // Check orthogonality: R^T R = I
    // For each pair (i,j), compute dot product of columns i and j
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double dot = 0.0;
            for (int k = 0; k < 3; ++k) {
                dot += R[k][i] * R[k][j];
            }
            double expected = (i == j) ? 1.0 : 0.0;
            if (std::abs(dot - expected) > tolerance) {
                return false;
            }
        }
    }
    
    // Check det(R) = 1 (proper rotation, no reflection)
    double det = R[0][0]*(R[1][1]*R[2][2] - R[1][2]*R[2][1])
               - R[0][1]*(R[1][0]*R[2][2] - R[1][2]*R[2][0])
               + R[0][2]*(R[1][0]*R[2][1] - R[1][1]*R[2][0]);
    
    return std::abs(det - 1.0) < tolerance;
}

} // namespace ICARION::physics::collision_core
