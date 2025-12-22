// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file test_collision_geometry.cpp
 * @brief Unit tests for CollisionGeometry module
 * 
 * Tests geometric operations used in collision calculations:
 * - Orthonormal basis construction
 * - Random rotation generation
 * - Vector rotation
 * - Rotation matrix validation
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/physics/collisions/core/CollisionGeometry.h"
#include "core/utils/mathUtils.h"  // For norm(), dot(), cross(), normalize()
#include <cmath>

using namespace ICARION::physics::collision_core;
using namespace ICARION::core;  // For Vec3
using ICARION::physics::PhysicsRng;
using Catch::Matchers::WithinAbs;

// Numerical tolerance for floating-point comparisons
constexpr double TOL = 1e-10;

TEST_CASE("CollisionGeometry: Orthonormal basis construction", "[collision][geometry]") {
    SECTION("Basis perpendicular to z-axis") {
        Vec3 ehat{0.0, 0.0, 1.0};
        Vec3 t1, t2;
        
        CollisionGeometry::construct_orthonormal_basis(ehat, t1, t2);
        
        // Check orthogonality
        REQUIRE_THAT(dot(t1, ehat), WithinAbs(0.0, TOL));
        REQUIRE_THAT(dot(t2, ehat), WithinAbs(0.0, TOL));
        REQUIRE_THAT(dot(t1, t2), WithinAbs(0.0, TOL));
        
        // Check normalization
        REQUIRE_THAT(norm(t1), WithinAbs(1.0, TOL));
        REQUIRE_THAT(norm(t2), WithinAbs(1.0, TOL));
        
        // Check right-handed coordinate system (ehat·(t1×t2) = 1)
        Vec3 cross_product = cross(t1, t2);
        REQUIRE_THAT(dot(ehat, cross_product), WithinAbs(1.0, TOL));
    }
    
    SECTION("Basis perpendicular to x-axis") {
        Vec3 ehat{1.0, 0.0, 0.0};
        Vec3 t1, t2;
        
        CollisionGeometry::construct_orthonormal_basis(ehat, t1, t2);
        
        // Check orthogonality
        REQUIRE_THAT(dot(t1, ehat), WithinAbs(0.0, TOL));
        REQUIRE_THAT(dot(t2, ehat), WithinAbs(0.0, TOL));
        REQUIRE_THAT(dot(t1, t2), WithinAbs(0.0, TOL));
        
        // Check normalization
        REQUIRE_THAT(norm(t1), WithinAbs(1.0, TOL));
        REQUIRE_THAT(norm(t2), WithinAbs(1.0, TOL));
    }
    
    SECTION("Basis perpendicular to arbitrary direction") {
        Vec3 ehat = normalize(Vec3{1.0, 2.0, 3.0});
        Vec3 t1, t2;
        
        CollisionGeometry::construct_orthonormal_basis(ehat, t1, t2);
        
        // Check orthogonality
        REQUIRE_THAT(dot(t1, ehat), WithinAbs(0.0, TOL));
        REQUIRE_THAT(dot(t2, ehat), WithinAbs(0.0, TOL));
        REQUIRE_THAT(dot(t1, t2), WithinAbs(0.0, TOL));
        
        // Check normalization
        REQUIRE_THAT(norm(t1), WithinAbs(1.0, TOL));
        REQUIRE_THAT(norm(t2), WithinAbs(1.0, TOL));
    }
}

TEST_CASE("CollisionGeometry: Random rotation matrix generation", "[collision][geometry]") {
    SECTION("Generated rotation is valid") {
        PhysicsRng rng(42);
        double R[3][3];
        
        CollisionGeometry::generate_random_rotation(rng, R);
        
        REQUIRE(CollisionGeometry::is_valid_rotation(R));
    }
    
    SECTION("Multiple rotations are all valid") {
        PhysicsRng rng(123);
        
        for (int i = 0; i < 100; ++i) {
            double R[3][3];
            CollisionGeometry::generate_random_rotation(rng, R);
            
            REQUIRE(CollisionGeometry::is_valid_rotation(R));
        }
    }
    
    SECTION("Rotation matrix is orthogonal (R^T R = I)") {
        PhysicsRng rng(42);
        double R[3][3];
        
        CollisionGeometry::generate_random_rotation(rng, R);
        
        // Compute R^T R
        double RtR[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                RtR[i][j] = 0.0;
                for (int k = 0; k < 3; ++k) {
                    RtR[i][j] += R[k][i] * R[k][j];
                }
            }
        }
        
        // Check diagonal elements (should be 1)
        REQUIRE_THAT(RtR[0][0], WithinAbs(1.0, TOL));
        REQUIRE_THAT(RtR[1][1], WithinAbs(1.0, TOL));
        REQUIRE_THAT(RtR[2][2], WithinAbs(1.0, TOL));
        
        // Check off-diagonal elements (should be 0)
        REQUIRE_THAT(RtR[0][1], WithinAbs(0.0, TOL));
        REQUIRE_THAT(RtR[0][2], WithinAbs(0.0, TOL));
        REQUIRE_THAT(RtR[1][0], WithinAbs(0.0, TOL));
        REQUIRE_THAT(RtR[1][2], WithinAbs(0.0, TOL));
        REQUIRE_THAT(RtR[2][0], WithinAbs(0.0, TOL));
        REQUIRE_THAT(RtR[2][1], WithinAbs(0.0, TOL));
    }
    
    SECTION("Determinant is +1 (proper rotation)") {
        PhysicsRng rng(42);
        double R[3][3];
        
        CollisionGeometry::generate_random_rotation(rng, R);
        
        double det = R[0][0]*(R[1][1]*R[2][2] - R[1][2]*R[2][1])
                   - R[0][1]*(R[1][0]*R[2][2] - R[1][2]*R[2][0])
                   + R[0][2]*(R[1][0]*R[2][1] - R[1][1]*R[2][0]);
        
        REQUIRE_THAT(det, WithinAbs(1.0, TOL));
    }
}

TEST_CASE("CollisionGeometry: Vector rotation", "[collision][geometry]") {
    SECTION("Rotation preserves vector norm") {
        PhysicsRng rng(42);
        double R[3][3];
        CollisionGeometry::generate_random_rotation(rng, R);
        
        Vec3 v{1.0, 2.0, 3.0};
        Vec3 v_rot = CollisionGeometry::rotate_vector(v, R);
        
        REQUIRE_THAT(norm(v), WithinAbs(norm(v_rot), TOL));
    }
    
    SECTION("Identity rotation leaves vector unchanged") {
        double I[3][3] = {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        };
        
        Vec3 v{1.0, 2.0, 3.0};
        Vec3 v_rot = CollisionGeometry::rotate_vector(v, I);
        
        REQUIRE_THAT(v.x, WithinAbs(v_rot.x, TOL));
        REQUIRE_THAT(v.y, WithinAbs(v_rot.y, TOL));
        REQUIRE_THAT(v.z, WithinAbs(v_rot.z, TOL));
    }
    
    SECTION("90-degree rotation around z-axis") {
        // Rotation matrix for 90° around z-axis
        double Rz[3][3] = {
            {0, -1, 0},
            {1,  0, 0},
            {0,  0, 1}
        };
        
        Vec3 v{1.0, 0.0, 0.0};  // x-axis
        Vec3 v_rot = CollisionGeometry::rotate_vector(v, Rz);
        
        // Should rotate to y-axis
        REQUIRE_THAT(v_rot.x, WithinAbs(0.0, TOL));
        REQUIRE_THAT(v_rot.y, WithinAbs(1.0, TOL));
        REQUIRE_THAT(v_rot.z, WithinAbs(0.0, TOL));
    }
    
    SECTION("Random rotation preserves norm for multiple vectors") {
        PhysicsRng rng(123);
        double R[3][3];
        CollisionGeometry::generate_random_rotation(rng, R);
        
        std::vector<Vec3> test_vectors = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {1.0, 2.0, 3.0},
            {-1.5, 2.7, -0.3}
        };
        
        for (const auto& v : test_vectors) {
            Vec3 v_rot = CollisionGeometry::rotate_vector(v, R);
            REQUIRE_THAT(norm(v), WithinAbs(norm(v_rot), TOL));
        }
    }
}

TEST_CASE("CollisionGeometry: Rotation matrix validation", "[collision][geometry]") {
    SECTION("Identity matrix is valid") {
        double I[3][3] = {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        };
        
        REQUIRE(CollisionGeometry::is_valid_rotation(I));
    }
    
    SECTION("90-degree rotation around z-axis is valid") {
        double Rz[3][3] = {
            {0, -1, 0},
            {1,  0, 0},
            {0,  0, 1}
        };
        
        REQUIRE(CollisionGeometry::is_valid_rotation(Rz));
    }
    
    SECTION("Scaling matrix is not a valid rotation") {
        double S[3][3] = {
            {2, 0, 0},
            {0, 2, 0},
            {0, 0, 2}
        };
        
        REQUIRE_FALSE(CollisionGeometry::is_valid_rotation(S));
    }
    
    SECTION("Reflection matrix is not a valid rotation (det = -1)") {
        double Ref[3][3] = {
            {-1, 0, 0},
            {0,  1, 0},
            {0,  0, 1}
        };
        
        REQUIRE_FALSE(CollisionGeometry::is_valid_rotation(Ref));
    }
    
    SECTION("Non-orthogonal matrix is not valid") {
        double M[3][3] = {
            {1, 0.1, 0},
            {0, 1,   0},
            {0, 0,   1}
        };
        
        REQUIRE_FALSE(CollisionGeometry::is_valid_rotation(M));
    }
}

TEST_CASE("CollisionGeometry: Regression test against old implementation", "[collision][geometry][regression]") {
    SECTION("Orthonormal basis matches old ortho_basis() behavior") {
        // Test specific cases that should produce identical results
        Vec3 ehat_z{0.0, 0.0, 1.0};
        Vec3 t1_z, t2_z;
        CollisionGeometry::construct_orthonormal_basis(ehat_z, t1_z, t2_z);
        
        // Expected behavior from old implementation:
        // For ehat = (0,0,1), arbitrary vector is (1,0,0)
        // t1 = normalize((0,0,1) × (1,0,0)) = normalize((0,1,0)) = (0,1,0)
        // t2 = normalize((0,0,1) × (0,1,0)) = normalize((-1,0,0)) = (-1,0,0)
        
        REQUIRE_THAT(t1_z.x, WithinAbs(0.0, TOL));
        REQUIRE_THAT(t1_z.y, WithinAbs(1.0, TOL));
        REQUIRE_THAT(t1_z.z, WithinAbs(0.0, TOL));
        
        REQUIRE_THAT(t2_z.x, WithinAbs(-1.0, TOL));
        REQUIRE_THAT(t2_z.y, WithinAbs(0.0, TOL));
        REQUIRE_THAT(t2_z.z, WithinAbs(0.0, TOL));
    }
    
    SECTION("Random rotation with fixed seed produces consistent results") {
        // Same seed should produce identical rotation matrix
        PhysicsRng rng1(42);
        PhysicsRng rng2(42);
        
        double R1[3][3];
        double R2[3][3];
        
        CollisionGeometry::generate_random_rotation(rng1, R1);
        CollisionGeometry::generate_random_rotation(rng2, R2);
        
        // Matrices should be identical
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                REQUIRE_THAT(R1[i][j], WithinAbs(R2[i][j], TOL));
            }
        }
    }
}
