// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
/**
 * @file CollisionGeometry.h
 * @brief Low-level geometric operations for collision calculations
 * 
 * Stateless geometric helpers used by collision kernels. Random rotations are
 * sampled uniformly on SO(3) via the quaternion method.
 */

#pragma once

#include "core/types/Vec3.h"
#include "core/types/CollisionTypes.h"  // For PhysicsRng
#include <array>

namespace ICARION::physics::collision_core {

/**
 * @brief Collision geometry utilities
 * 
 * Low-level geometric operations for collision calculations.
 * Pure static functions (no state, no side effects).
 */
class CollisionGeometry {
public:
    /**
     * @brief Construct orthonormal basis perpendicular to given vector
     * 
     * Given a unit vector @p ehat, computes two orthonormal vectors @p t1 and @p t2 such that:
     * - @f$ t_1 \perp \hat{e} @f$
     * - @f$ t_2 \perp \hat{e},\ t_2 \perp t_1 @f$
     * - @f$ \{ \hat{e}, t_1, t_2 \} @f$ forms a right-handed orthonormal basis
     * 
     * Uses Gram-Schmidt orthogonalization to construct tangent vectors.
     * Useful for defining impact planes or scattering coordinate frames.
     * 
     * @param ehat Input unit vector (should be normalized)
     * @param t1 Output: first tangent vector (perpendicular to ehat)
     * @param t2 Output: second tangent vector (perpendicular to both)
     * 
     * @pre ehat.norm() ≈ 1.0 (within numerical precision)
     * @post t1·ehat = 0, t2·ehat = 0, t1·t2 = 0
     * @post |t1| = 1, |t2| = 1
     * 
     * @note Equivalent to old `ortho_basis()` function from collisionHelpers.h
     */
    static void construct_orthonormal_basis(
        const Vec3& ehat,
        Vec3& t1,
        Vec3& t2
    );
    
    /**
     * @brief Generate uniform random rotation matrix
     * 
     * Samples a uniformly distributed unit quaternion using three random numbers
     * and converts it to a rotation matrix. The resulting matrix is orthogonal
     * with determinant +1.
     * 
     * Uses quaternion method (Shoemake 1992) for unbiased SO(3) sampling.
     * Ensures isotropic orientation distribution for collision scattering.
     * 
     * @param rng Random number generator
     * @param R Output: 3x3 rotation matrix R[i][j]
     * 
     * @post R is orthogonal (R^T R = I)
     * @post det(R) = 1 (proper rotation, no reflection)
     * 
     * @note Equivalent to old `rand_rotation()` function from collisionHelpers.h
     */
    static void generate_random_rotation(
        PhysicsRng& rng,
        double R[3][3]
    );
    
    /**
     * @brief Rotate vector by rotation matrix
     * 
     * Applies rotation matrix R to vector v: result = R·v
     * 
     * @param v Input vector
     * @param R Rotation matrix R[i][j] (row-major indexing)
     * @return Rotated vector R·v
     * 
     * @post result.norm() == v.norm() (if R is orthogonal)
     */
    static Vec3 rotate_vector(
        const Vec3& v,
        const double R[3][3]
    );

    /**
     * @brief Convert quaternion to rotation matrix
     *
     * Quaternion is expected as [w, x, y, z]. The result is a proper rotation
     * matrix suitable for rotating vectors in the lab frame.
     *
     * @param quat Input quaternion [w, x, y, z]
     * @param R Output: 3x3 rotation matrix R[i][j]
     */
    static void quaternion_to_rotation(
        const std::array<double, 4>& quat,
        double R[3][3]
    );
    
    /**
     * @brief Check if rotation matrix is valid
     * 
     * Validates that R is orthogonal (R^T R = I) and proper (det(R) = 1).
     * Used for debugging and testing.
     * 
     * @param R Rotation matrix to validate
     * @param tolerance Numerical tolerance (default: 1e-6)
     * @return true if R is a valid proper rotation matrix
     */
    static bool is_valid_rotation(
        const double R[3][3],
        double tolerance = 1e-6
    );
};

} // namespace ICARION::physics::collision_core
