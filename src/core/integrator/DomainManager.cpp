// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "DomainManager.h"
#include "core/config/types/InstrumentTypes.h"  // Instrument enum
#include <stdexcept>
#include <cmath>

namespace ICARION {
namespace integrator {

namespace {
/**
 * @brief Compute the radial coordinate r(z) on an Orbitrap hyperbolic electrode.
 * 
 * The Orbitrap electrode shape is defined by the implicit equation:
 *   z² = 0.5·(r² - R²) + R_m² · ln(R/r)
 * 
 * This function numerically solves for r given z, R (electrode radius), and R_m (characteristic radius).
 * Uses bisection method for robust convergence.
 * 
 * @param z Axial position (can be negative, will use |z|)
 * @param R Reference electrode radius (Rin or Rout)
 * @param Rm Characteristic radius parameter
 * @return Radial coordinate r at position z
 */
double orbitrap_r_for_z(double z, double R, double Rm) {
    // Solve: z² = 0.5·(r² - R²) + R_m² · ln(R/r)
    // 
    // WARNING: This implementation is ported from legacy paramUtils.cpp and has
    // known issues with bisection convergence when the solution lies outside [R, 2R].
    // The bisection may not converge to the true root if:
    //   - f(R) and f(2R) have the same sign (no bracketing)
    //   - Rm is not in a "reasonable" range relative to R
    // 
    // This produces incorrect boundaries for certain parameter combinations!
    // TODO(Phase 6): Implement robust root-finding with adaptive bracketing
    const double z2 = z * z;
    const double eps = 1e-12;
    const int max_iter = 100;
    
    double r_lo = R;
    double r_hi = R * 2.0;
    double r_mid = R;
    
    for (int i = 0; i < max_iter; ++i) {
        r_mid = 0.5 * (r_lo + r_hi);
        double f = 0.5 * (r_mid * r_mid - R * R) + Rm * Rm * std::log(R / r_mid) - z2;
        
        if (std::fabs(f) < eps) break;
        if (f > 0) r_lo = r_mid;
        else       r_hi = r_mid;
    }
    
    return r_mid;
}
}  // anonymous namespace

DomainManager::DomainManager(const std::vector<config::DomainConfig>& domains)
    : domains_(domains)
{
    if (domains.empty()) {
        throw std::invalid_argument("DomainManager: domains vector is empty");
    }
}

int DomainManager::find_domain_index(const Vec3& pos) const {
    for (size_t i = 0; i < domains_.size(); ++i) {
        if (is_inside_domain(domains_[i], pos)) {
            return static_cast<int>(i);
        }
    }
    return -1;  // Outside all domains
}

const config::DomainConfig& DomainManager::get_domain(int idx) const {
    if (idx < 0 || idx >= static_cast<int>(domains_.size())) {
        throw std::out_of_range("DomainManager::get_domain: invalid domain index " + 
                                std::to_string(idx));
    }
    return domains_[idx];
}

Vec3 DomainManager::global_to_local_pos(const Vec3& pos, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    return dom.rotation_global_to_local * (pos - dom.geometry.origin_m);
}

Vec3 DomainManager::global_to_local_vel(const Vec3& vel, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    return dom.rotation_global_to_local * vel;
}

Vec3 DomainManager::local_to_global_pos(const Vec3& pos_local, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    return dom.rotation_local_to_global * pos_local + dom.geometry.origin_m;
}

Vec3 DomainManager::local_to_global_vel(const Vec3& vel_local, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    return dom.rotation_local_to_global * vel_local;
}

void DomainManager::check_aperture_crossing(IonState& ion, int domain_idx,
                                             const Vec3& pos_before, const Vec3& pos_after) const {
    const auto& dom = get_domain(domain_idx);
    
    // No aperture constraint
    if (dom.geometry.end_aperture_m <= 0.0) {
        return;  // Ion can pass freely
    }
    
    const double z_ap = dom.geometry.length_m;
    
    // Check if ion crossed the aperture plane (forward or backward)
    const bool crossed_forward = (pos_before.z < z_ap && pos_after.z >= z_ap);
    const bool crossed_backward = (pos_before.z > z_ap && pos_after.z <= z_ap);
    
    if (crossed_forward || crossed_backward) {
        // Interpolate crossing point
        const double alpha = (z_ap - pos_before.z) / (pos_after.z - pos_before.z);
        const Vec3 cross_point = pos_before + (pos_after - pos_before) * alpha;
        
        // Check radial distance at crossing
        const double r_cross = std::sqrt(cross_point.x * cross_point.x + 
                                         cross_point.y * cross_point.y);
        
        if (r_cross > dom.geometry.end_aperture_m) {
            // Ion blocked by aperture
            ion.active = false;
        }
    }
}

void DomainManager::update_domain_properties(IonState& ion, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    
    ion.domain_neutral_mass_kg = dom.environment.gas_mass_kg;
    ion.domain_temperature_K = dom.environment.temperature_K;
    ion.domain_particle_density_m3 = dom.environment.particle_density_m_3;
    ion.domain_gas_velocity_m_s = dom.environment.gas_velocity_m_s;
    ion.current_domain_index = dom.domain_index;
}

void DomainManager::terminate_ion_at_boundary(IonState& ion, int domain_idx,
                                                const Vec3& pos_before_local,
                                                const Vec3& pos_after_local) const {
    const auto& dom = get_domain(domain_idx);
    
    // Orbitrap: Hyperbolic boundary ray tracing
    // Hyperboloid equation: z² - r²/2 = C
    // Ray: p(t) = p0 + t*d, where d = (p_after - p_before)
    if (dom.instrument == config::Instrument::Orbitrap) {
        Vec3 dir = pos_after_local - pos_before_local;
        double dir_len = std::sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
        
        if (dir_len < NUMERICAL_ZERO) {
            // Degenerate case: no movement
            ion.pos = local_to_global_pos(pos_before_local, domain_idx);
            ion.vel = {0.0, 0.0, 0.0};
            ion.active = false;
            return;
        }
        
        // Normalize direction
        dir = dir * (1.0 / dir_len);
        
        // Solve for intersection with hyperboloids: (z0 + t*dz)² - 0.5*(r0 + t*dr)² = C
        // This gives a quadratic equation: a*t² + b*t + c = 0
        // where r² = x² + y², dr = (x*dx + y*dy)/r (radial direction component)
        
        double t_min = dir_len;  // Default: full step
        Vec3 intersection = pos_after_local;
        
        // Test both inner and outer hyperboloids
        for (int surface = 0; surface < 2; ++surface) {
            const double C = (surface == 0) ? dom.geometry.orbitrap_C_in : dom.geometry.orbitrap_C_out;
            
            // Quadratic coefficients for z² - 0.5*r² = C
            // Expanding: (z0 + t*dz)² - 0.5*((x0 + t*dx)² + (y0 + t*dy)²) = C
            const double a = dir.z * dir.z - 0.5 * (dir.x * dir.x + dir.y * dir.y);
            const double b = 2.0 * (pos_before_local.z * dir.z 
                                    - 0.5 * (pos_before_local.x * dir.x + pos_before_local.y * dir.y));
            const double r0_sq = pos_before_local.x * pos_before_local.x + pos_before_local.y * pos_before_local.y;
            const double c = pos_before_local.z * pos_before_local.z - 0.5 * r0_sq - C;
            
            const double discriminant = b*b - 4.0*a*c;
            if (discriminant < 0.0) continue;  // No intersection
            
            const double sqrt_disc = std::sqrt(discriminant);
            const double t1 = (-b - sqrt_disc) / (2.0 * a);
            const double t2 = (-b + sqrt_disc) / (2.0 * a);
            
            // Take smallest positive t in range [0, dir_len]
            for (double t : {t1, t2}) {
                if (t > 0.0 && t < t_min) {
                    t_min = t;
                    intersection = pos_before_local + dir * t;
                }
            }
        }
        
        ion.pos = local_to_global_pos(intersection, domain_idx);
        ion.vel = {0.0, 0.0, 0.0};
        ion.active = false;
        return;
    }
    
    // Cylindrical geometry: Ray tracing with analytical intersection
    // Find intersection of line (pos_before → pos_after) with domain boundary
    // We test all boundaries and take the nearest intersection
    
    Vec3 dir = pos_after_local - pos_before_local;
    double dir_len = std::sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
    if (dir_len < NUMERICAL_ZERO) {
        // Degenerate case: no movement, just place at pos_before
        ion.pos = local_to_global_pos(pos_before_local, domain_idx);
        ion.vel = {0.0, 0.0, 0.0};
        ion.active = false;
        return;
    }
    dir = dir * (1.0 / dir_len);  // Normalize
    
    double t_min = dir_len;  // Default: full step (if no intersection found)
    Vec3 intersection = pos_after_local;
    
    // 1. Check radial boundary (cylindrical wall): r = radius_m
    //    Solve: |pos_before + t*dir|_xy = radius_m
    //    (x0 + t*dx)^2 + (y0 + t*dy)^2 = R^2
    const double R = dom.geometry.radius_m;
    const double a = dir.x*dir.x + dir.y*dir.y;
    const double b = 2.0*(pos_before_local.x*dir.x + pos_before_local.y*dir.y);
    const double c = pos_before_local.x*pos_before_local.x + 
                     pos_before_local.y*pos_before_local.y - R*R;
    
    if (a > NUMERICAL_ZERO) {  // Ray has radial component
        const double discriminant = b*b - 4.0*a*c;
        if (discriminant >= 0.0) {
            const double sqrt_disc = std::sqrt(discriminant);
            const double t1 = (-b + sqrt_disc) / (2.0*a);
            const double t2 = (-b - sqrt_disc) / (2.0*a);
            
            // Take smallest positive t (first intersection)
            for (double t : {t1, t2}) {
                if (t > 0.0 && t < t_min) {
                    Vec3 candidate = pos_before_local + dir * t;
                    // Check if within axial bounds
                    if (candidate.z >= -DOMAIN_BOUNDARY_EPSILON && candidate.z <= dom.geometry.length_m) {
                        t_min = t;
                        intersection = candidate;
                    }
                }
            }
        }
    }
    
    // 2. Check entrance plane: z = 0
    if (std::abs(dir.z) > NUMERICAL_ZERO) {
        double t = (0.0 - pos_before_local.z) / dir.z;
        if (t > 0.0 && t < t_min) {
            Vec3 candidate = pos_before_local + dir * t;
            double r = std::sqrt(candidate.x*candidate.x + candidate.y*candidate.y);
            if (r <= R) {
                t_min = t;
                intersection = candidate;
            }
        }
    }
    
    // 3. Check exit plane: z = length_m
    if (std::abs(dir.z) > NUMERICAL_ZERO) {
        double t = (dom.geometry.length_m - pos_before_local.z) / dir.z;
        if (t > 0.0 && t < t_min) {
            Vec3 candidate = pos_before_local + dir * t;
            double r = std::sqrt(candidate.x*candidate.x + candidate.y*candidate.y);
            // Check aperture constraint if defined
            double aperture_limit = (dom.geometry.end_aperture_m > 0.0) ? 
                                    dom.geometry.end_aperture_m : R;
            if (r <= aperture_limit) {
                t_min = t;
                intersection = candidate;
            }
        }
    }
    
    // Set ion to intersection point (in global coordinates)
    ion.pos = local_to_global_pos(intersection, domain_idx);
    ion.vel = {0.0, 0.0, 0.0};  // Ion absorbed/stopped
    ion.active = false;
}

bool DomainManager::is_inside_domain(const config::DomainConfig& dom, const Vec3& globalPos) const {
    // Transform to local coordinates
    Vec3 local = dom.rotation_global_to_local * (globalPos - dom.geometry.origin_m);
    double r = std::sqrt(local.x * local.x + local.y * local.y);
    
    // Cylindrical geometry (most instruments)
    if (dom.instrument != config::Instrument::Orbitrap) {
        // Ion must be inside domain:
        // - z >= -ε (at or past entrance, with FP tolerance for start position)
        // - z < length_m (before exit)
        // - r < radius_m (within radial boundary)
        // Note: Negative epsilon critical for IMS/SIFDT where ions start at z=0.0
        //       but may have z=-1e-16 due to floating-point roundoff
        return (local.z >= -DOMAIN_BOUNDARY_EPSILON && local.z < dom.geometry.length_m) && (r < dom.geometry.radius_m);
    }
    
    // Orbitrap: Hyperbolic electrode geometry (equipotential surfaces)
    // Hyperboloid equation: z² - r²/2 = C
    // Ion is inside if: C_out <= C <= C_in, where:
    //   C_in  = -0.5 * R_in²  (inner, smaller radius → less negative)
    //   C_out = -0.5 * R_out² (outer, larger radius → more negative)
    // Note: C_out < C_in since R_out > R_in
    const double r2 = local.x * local.x + local.y * local.y;
    const double C = local.z * local.z - 0.5 * r2;
    constexpr double eps = 1e-12;  // Numerical tolerance
    
    return (C >= dom.geometry.orbitrap_C_out - eps) && (C <= dom.geometry.orbitrap_C_in + eps);
}

}  // namespace integrator
}  // namespace ICARION
