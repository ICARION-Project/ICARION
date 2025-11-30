// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file ReflectionAction.h
 * @brief Boundary action: reflect ion velocity at boundary
 * 
 * Three types:
 * - Specular: elastic bounce (v' = v - 2(v·n)n)
 * - Diffuse: random cosine-weighted direction
 * - Thermal: Maxwell-Boltzmann re-emission at wall temperature
 */

#pragma once

#include "BoundaryAction.h"
#include "utils/constants.h"
#include <random>
#include <cmath>

namespace ICARION {
namespace integrator {

/**
 * @brief Reflection action with multiple modes
 */
class ReflectionAction : public BoundaryAction {
public:
    enum class Type {
        SPECULAR,   ///< Mirror reflection: v' = v - 2(v·n)n
        DIFFUSE,    ///< Cosine-weighted random direction
        THERMAL     ///< Maxwell-Boltzmann re-emission
    };
    
    /**
     * @brief Construct reflection action
     * @param type Reflection type (specular/diffuse/thermal)
     * @param accommodation_coeff Energy accommodation coefficient [0,1]
     *        0 = fully specular, 1 = fully thermal
     * @param rng Random number generator (for diffuse/thermal)
     */
    ReflectionAction(
        Type type,
        double accommodation_coeff = 1.0,
        std::mt19937* rng = nullptr
    ) : type_(type),
        accommodation_coeff_(accommodation_coeff),
        rng_(rng)
    {
        if ((type == Type::DIFFUSE || type == Type::THERMAL) && !rng) {
            throw std::invalid_argument("Diffuse/Thermal reflection requires RNG");
        }
    }
    
    void apply(
        IonState& ion,
        const Vec3& normal,
        const Vec3& boundary_pos,
        double temperature_K
    ) override {
        // Set position to boundary
        ion.pos = boundary_pos;
        
        // Apply reflection based on type
        switch (type_) {
            case Type::SPECULAR:
                apply_specular(ion, normal);
                break;
            case Type::DIFFUSE:
                apply_diffuse(ion, normal, temperature_K);
                break;
            case Type::THERMAL:
                apply_thermal(ion, normal, temperature_K);
                break;
        }
    }
    
    std::string name() const override {
        switch (type_) {
            case Type::SPECULAR: return "Specular Reflection";
            case Type::DIFFUSE: return "Diffuse Reflection";
            case Type::THERMAL: return "Thermal Reflection";
            default: return "Unknown Reflection";
        }
    }
    
private:
    Type type_;
    double accommodation_coeff_;
    std::mt19937* rng_;
    
    /**
     * @brief Specular reflection: v' = v - 2(v·n)n
     */
    void apply_specular(IonState& ion, const Vec3& normal) {
        double v_dot_n = ion.vel.x * normal.x + ion.vel.y * normal.y + ion.vel.z * normal.z;
        ion.vel = ion.vel - normal * (2.0 * v_dot_n);
    }
    
    /**
     * @brief Diffuse reflection: cosine-weighted random direction
     * 
     * Physics:
     * - Outgoing direction sampled from cos(θ) distribution
     * - Speed partially thermalized based on accommodation coefficient
     */
    void apply_diffuse(IonState& ion, const Vec3& normal, double temperature_K) {
        // Sample outgoing direction (cosine-weighted)
        Vec3 outgoing_dir = sample_cosine_hemisphere(normal);
        
        // Compute thermal speed (RMS velocity)
        double v_thermal = std::sqrt(2.0 * BOLTZMANN_CONSTANT * temperature_K / ion.mass_kg);
        
        // Blend between specular and thermal speed
        double v_in = std::sqrt(ion.vel.x*ion.vel.x + ion.vel.y*ion.vel.y + ion.vel.z*ion.vel.z);
        double v_out = (1.0 - accommodation_coeff_) * v_in + accommodation_coeff_ * v_thermal;
        
        ion.vel = outgoing_dir * v_out;
    }
    
    /**
     * @brief Thermal reflection: Maxwell-Boltzmann re-emission
     * 
     * Physics:
     * - Velocity magnitude sampled from Maxwell-Boltzmann at wall temperature
     * - Direction sampled from cosine-weighted hemisphere
     */
    void apply_thermal(IonState& ion, const Vec3& normal, double temperature_K) {
        // Sample velocity magnitude from Maxwell-Boltzmann distribution
        // Using Box-Muller transform for 3D Gaussian → Maxwell-Boltzmann speed
        std::normal_distribution<double> normal_dist(0.0, 1.0);
        double sigma = std::sqrt(BOLTZMANN_CONSTANT * temperature_K / ion.mass_kg);
        
        double vx = sigma * normal_dist(*rng_);
        double vy = sigma * normal_dist(*rng_);
        double vz = sigma * normal_dist(*rng_);
        double v_mag = std::sqrt(vx*vx + vy*vy + vz*vz);
        
        // Sample outgoing direction (cosine-weighted)
        Vec3 outgoing_dir = sample_cosine_hemisphere(normal);
        
        ion.vel = outgoing_dir * v_mag;
    }
    
    /**
     * @brief Sample direction from cosine-weighted hemisphere
     * 
     * Uses Malley's method:
     * 1. Sample (x,y) from unit disk
     * 2. z = sqrt(1 - x^2 - y^2)
     * 3. Transform to local frame aligned with normal
     * 
     * @param normal Surface normal (pointing inward, unit length)
     * @return Unit vector in outgoing hemisphere
     */
    Vec3 sample_cosine_hemisphere(const Vec3& normal) {
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        
        // Sample unit disk (Malley's method)
        double r = std::sqrt(uniform(*rng_));
        double theta = 2.0 * M_PI * uniform(*rng_);
        double x = r * std::cos(theta);
        double y = r * std::sin(theta);
        double z = std::sqrt(1.0 - x*x - y*y);
        
        // Build orthonormal basis (n, u, v)
        Vec3 u, v;
        build_orthonormal_basis(normal, u, v);
        
        // Transform to global frame
        return u * x + v * y + normal * z;
    }
    
    /**
     * @brief Build orthonormal basis from normal vector
     * 
     * Given normal n, find two perpendicular vectors u and v
     * such that (u, v, n) form a right-handed orthonormal basis.
     * 
     * @param n Normal vector (assumed unit length)
     * @param u Output: first tangent vector
     * @param v Output: second tangent vector
     */
    void build_orthonormal_basis(const Vec3& n, Vec3& u, Vec3& v) {
        // Choose initial vector not parallel to n
        Vec3 initial = (std::fabs(n.x) < 0.9) ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
        
        // Gram-Schmidt
        double initial_dot_n = initial.x*n.x + initial.y*n.y + initial.z*n.z;
        u = initial - n * initial_dot_n;
        double u_mag = std::sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
        u = u * (1.0 / u_mag);  // Normalize
        
        // Cross product: n × u
        v = Vec3{n.y*u.z - n.z*u.y, n.z*u.x - n.x*u.z, n.x*u.y - n.y*u.x};
    }
};

}  // namespace integrator
}  // namespace ICARION
