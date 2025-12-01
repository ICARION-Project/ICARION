// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "DomainManager.h"
#include "core/config/types/InstrumentTypes.h"
#include <stdexcept>
#include <cmath>

namespace ICARION {
namespace integrator {

// ==================== ORBITRAP HELPERS (from tmp/) ====================

namespace {

/**
 * @brief Residual function for Orbitrap hyperlogarithmic electrode surface
 * 
 * **Hyperlogarithmic Surface Equation (Makarov 2000):**
 * 
 *   z² = 0.5(r² - R²) + R_m² × ln(R/r)
 * 
 * Rearranged as implicit equation F(r,z) = 0:
 * 
 *   F(r,z) = z² - 0.5(r² - R²) - R_m² × ln(R/r) = 0
 * 
 * **Parameters:**
 * @param r Current radial coordinate [m] (ion position)
 * @param z Current axial coordinate [m] (ion position, symmetric: |z|)
 * @param R Electrode radius at z=0 [m] (R_in for inner, R_out for outer)
 * @param R_m Characteristic radius [m] (field shaping parameter, R_m > R_out)
 * 
 * **Return:**
 * @return Residual F(r,z):
 *   - F < 0: Point (r,z) is inside electrode (closer to axis than surface)
 *   - F = 0: Point (r,z) is ON electrode surface
 *   - F > 0: Point (r,z) is outside electrode (further from axis)
 * 
 * **Usage in Bisection:**
 * - Given z, find r where F(r,z) = 0 → r is surface radius at height z
 * - Bracket: r ∈ [r_lo, r_hi] where F(r_lo) × F(r_hi) < 0
 * 
 * **Physical Meaning:**
 * - At z=0: r = R (electrode passes through (R, 0))
 * - As |z| increases: r decreases (electrodes curve inward)
 * - Creates quadro-logarithmic potential U(r,z) for mass analysis
 */
inline double orbitrap_surface_residual(double r, double z, double R, double R_m) {
    // Hyperlogarithmic equation: z² = 0.5(r² - R²) + R_m² × ln(R/r)
    // Rearranged: z² - 0.5(r² - R²) - R_m² × ln(R/r) = 0
    const double term1 = z * z;
    const double term2 = 0.5 * (r * r - R * R);
    const double term3 = R_m * R_m * std::log(R / r);
    
    return term1 - term2 - term3;
}

/**
 * @brief Residual for 3D point (convenience wrapper)
 * 
 * @param p_local 3D position in local coordinates [m]
 * @param R Electrode radius at z=0 [m]
 * @param R_m Characteristic radius [m]
 * @return Residual F(r,z) where r = sqrt(x² + y²), z = |z_local|
 */
inline double orbitrap_surface_residual(const Vec3& p_local, double R, double R_m) {
    const double r = std::sqrt(p_local.x * p_local.x + p_local.y * p_local.y);
    return orbitrap_surface_residual(r, std::fabs(p_local.z), R, R_m);
}

}  // namespace

// ==================== DOMAINMANAGER IMPLEMENTATION ====================

DomainManager::DomainManager(
    const std::vector<config::DomainConfig>& domains,
    unsigned int rng_seed
) : domains_(domains), rng_(rng_seed)
{
    if (domains.empty()) {
        throw std::invalid_argument("DomainManager: domains vector is empty");
    }
    
    // Create boundary actions for each domain
    boundary_actions_.reserve(domains.size());
    for (const auto& domain : domains) {
        boundary_actions_.push_back(
            BoundaryActionFactory::create(domain.boundary, &rng_)
        );
    }
}

int DomainManager::find_domain_index(const Vec3& pos) const {
    for (size_t i = 0; i < domains_.size(); ++i) {
        if (is_inside_domain(domains_[i], pos)) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

const config::DomainConfig& DomainManager::get_domain(int idx) const {
    if (idx < 0 || idx >= static_cast<int>(domains_.size())) {
        throw std::out_of_range("DomainManager::get_domain: invalid index " + 
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
    return dom.geometry.origin_m + dom.rotation_local_to_global * pos_local;
}

Vec3 DomainManager::local_to_global_vel(const Vec3& vel_local, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    return dom.rotation_local_to_global * vel_local;
}

void DomainManager::check_aperture_crossing(IonState& ion, int domain_idx,
                                             const Vec3& pos_before, const Vec3& pos_after) const {
    const auto& dom = get_domain(domain_idx);
    
    if (dom.geometry.end_aperture_m <= 0.0) {
        return;  // No aperture constraint
    }
    
    Vec3 pb_local = global_to_local_pos(pos_before, domain_idx);
    Vec3 pa_local = global_to_local_pos(pos_after, domain_idx);
    
    bool before_exit = (pb_local.z < dom.geometry.length_m);
    bool after_exit = (pa_local.z >= dom.geometry.length_m);
    
    if (before_exit && after_exit) {
        double t = (dom.geometry.length_m - pb_local.z) / (pa_local.z - pb_local.z);
        Vec3 crossing = pb_local + (pa_local - pb_local) * t;
        double r_cross = std::sqrt(crossing.x*crossing.x + crossing.y*crossing.y);
        
        if (r_cross > dom.geometry.end_aperture_m) {
            ion.active = false;
        }
    }
}

void DomainManager::update_domain_properties(IonState& ion, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    
    // Domain cache fields removed from IonState - only update domain index
    ion.current_domain_index = dom.domain_index;
}

void DomainManager::terminate_ion_at_boundary(IonState& ion, int domain_idx,
                                                const Vec3& pos_before_local,
                                                const Vec3& pos_after_local) const {
    const auto& dom = get_domain(domain_idx);
    
    // ========== ORBITRAP: Hyperlogarithmic Electrode Boundary ==========
    if (dom.instrument == config::Instrument::Orbitrap) {
        Vec3 dir = pos_after_local - pos_before_local;
        double step_len = std::sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
        
        if (step_len < NUMERICAL_ZERO) {
            ion.pos = local_to_global_pos(pos_before_local, domain_idx);
            ion.vel = {0.0, 0.0, 0.0};
            ion.active = false;
            return;
        }
        
        dir = dir * (1.0 / step_len);  // Normalize
        
        const double Rin = dom.geometry.radius_in_m;
        const double Rout = dom.geometry.radius_out_m;
        const double R_m = dom.geometry.radius_char_m;
        const double z_max = 0.5 * dom.geometry.length_m;
        
        double t_min = 1.0;  // Parameter in [0,1]
        Vec3 intersection = pos_after_local;
        bool hit = false;
        
        auto point_on_segment = [&](double t) {
            return pos_before_local + dir * (t * step_len);
        };
        
        // Helper: intersect with electrode surface (R, R_m)
        auto intersect_surface = [&](double R) {
            const double f0 = orbitrap_surface_residual(pos_before_local, R, R_m);
            const double f1 = orbitrap_surface_residual(pos_after_local, R, R_m);
            
            if (f0 * f1 > 0.0) return;  // No crossing
            
            // Bisection
            double t_lo = 0.0;
            double t_hi = 1.0;
            double f_lo = f0;
            double f_hi = f1;
            const double eps = 1e-8;
            const int max_iter = 40;
            
            for (int i = 0; i < max_iter; ++i) {
                double t_mid = 0.5 * (t_lo + t_hi);
                Vec3 p_mid = point_on_segment(t_mid);
                double f_mid = orbitrap_surface_residual(p_mid, R, R_m);
                
                if (std::fabs(f_mid) < eps) {
                    if (t_mid < t_min) {
                        t_min = t_mid;
                        intersection = p_mid;
                        hit = true;
                    }
                    return;
                }
                
                if (f_mid * f_lo > 0.0) {
                    t_lo = t_mid;
                    f_lo = f_mid;
                } else {
                    t_hi = t_mid;
                    f_hi = f_mid;
                }
            }
            
            // Fallback: use midpoint
            double t_mid = 0.5 * (t_lo + t_hi);
            if (t_mid < t_min) {
                t_min = t_mid;
                intersection = point_on_segment(t_mid);
                hit = true;
            }
        };
        
        // Check both electrodes
        intersect_surface(Rin);
        intersect_surface(Rout);
        
        // Axial clipping
        auto intersect_axial_plane = [&](double z_plane) {
            if (std::fabs(dir.z) < NUMERICAL_ZERO) return;
            
            double t_plane = (z_plane - pos_before_local.z) / (dir.z * step_len);
            if (t_plane <= 0.0 || t_plane >= 1.0) return;
            
            Vec3 p = point_on_segment(t_plane);
            double r = std::sqrt(p.x * p.x + p.y * p.y);
            
            double r_in = orbitrap_r_for_z(std::fabs(z_plane), Rin, R_m);
            double r_out = orbitrap_r_for_z(std::fabs(z_plane), Rout, R_m);
            
            if (r >= r_in - DOMAIN_BOUNDARY_EPSILON &&
                r <= r_out + DOMAIN_BOUNDARY_EPSILON) {
                if (t_plane < t_min) {
                    t_min = t_plane;
                    intersection = p;
                    hit = true;
                }
            }
        };
        
        intersect_axial_plane(+z_max);
        intersect_axial_plane(-z_max);
        
        if (!hit) {
            intersection = pos_after_local;
        }
        
        ion.pos = local_to_global_pos(intersection, domain_idx);
        ion.vel = {0.0, 0.0, 0.0};
        ion.death_time_s = ion.t;
        ion.active = false;
        return;
    }
    
    // ========== CYLINDRICAL GEOMETRY (Existing Code) ==========
    Vec3 dir = pos_after_local - pos_before_local;
    double dir_len = std::sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
    
    if (dir_len < NUMERICAL_ZERO) {
        ion.pos = local_to_global_pos(pos_before_local, domain_idx);
        ion.vel = {0.0, 0.0, 0.0};
        ion.death_time_s = ion.t;
        ion.active = false;
        return;
    }
    
    dir = dir * (1.0 / dir_len);
    
    double t_min = dir_len;
    Vec3 intersection = pos_after_local;
    
    // 1. Radial boundary
    const double R = dom.geometry.radius_m;
    const double a = dir.x*dir.x + dir.y*dir.y;
    const double b = 2.0*(pos_before_local.x*dir.x + pos_before_local.y*dir.y);
    const double c = pos_before_local.x*pos_before_local.x + 
                     pos_before_local.y*pos_before_local.y - R*R;
    
    if (a > NUMERICAL_ZERO) {
        const double discriminant = b*b - 4.0*a*c;
        if (discriminant >= 0.0) {
            const double sqrt_disc = std::sqrt(discriminant);
            const double t1 = (-b + sqrt_disc) / (2.0*a);
            const double t2 = (-b - sqrt_disc) / (2.0*a);
            
            for (double t : {t1, t2}) {
                if (t > 0.0 && t < t_min) {
                    Vec3 candidate = pos_before_local + dir * t;
                    if (candidate.z >= -DOMAIN_BOUNDARY_EPSILON && 
                        candidate.z <= dom.geometry.length_m) {
                        t_min = t;
                        intersection = candidate;
                    }
                }
            }
        }
    }
    
    // 2. Entrance plane (z=0)
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
    
    // 3. Exit plane (z=length_m)
    if (std::abs(dir.z) > NUMERICAL_ZERO) {
        double t = (dom.geometry.length_m - pos_before_local.z) / dir.z;
        if (t > 0.0 && t < t_min) {
            Vec3 candidate = pos_before_local + dir * t;
            double r = std::sqrt(candidate.x*candidate.x + candidate.y*candidate.y);
            double aperture_limit = (dom.geometry.end_aperture_m > 0.0) ? 
                                    dom.geometry.end_aperture_m : R;
            if (r <= aperture_limit) {
                t_min = t;
                intersection = candidate;
            }
        }
    }
    
    // Compute surface normal at intersection (pointing inward)
    Vec3 normal = compute_surface_normal(intersection, domain_idx);
    
    // Apply boundary action
    const auto& boundary_action = boundary_actions_[domain_idx];
    boundary_action->apply(
        ion,
        normal,
        local_to_global_pos(intersection, domain_idx),
        dom.environment.temperature_K,
        ion.t  // Current simulation time
    );
}

Vec3 DomainManager::compute_surface_normal(const Vec3& pos_local, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    const double R = dom.geometry.radius_m;
    const double EPSILON = 1e-6;
    
    // Determine which surface was hit
    double r = std::sqrt(pos_local.x * pos_local.x + pos_local.y * pos_local.y);
    
    // Radial boundary (cylindrical wall)
    if (std::abs(r - R) < EPSILON) {
        // Normal points radially inward: -r_hat
        return Vec3{-pos_local.x / r, -pos_local.y / r, 0.0};
    }
    
    // Entrance plane (z=0)
    if (std::abs(pos_local.z) < EPSILON) {
        // Normal points into domain: +z direction
        return Vec3{0.0, 0.0, 1.0};
    }
    
    // Exit plane (z=length_m)
    if (std::abs(pos_local.z - dom.geometry.length_m) < EPSILON) {
        // Normal points into domain: -z direction
        return Vec3{0.0, 0.0, -1.0};
    }
    
    // Default: radial inward (safest guess)
    if (r > NUMERICAL_ZERO) {
        return Vec3{-pos_local.x / r, -pos_local.y / r, 0.0};
    } else {
        return Vec3{0.0, 0.0, 1.0};  // Fallback: axial
    }
}

bool DomainManager::is_inside_domain(const config::DomainConfig& dom, const Vec3& globalPos) const {
    Vec3 local = dom.rotation_global_to_local * (globalPos - dom.geometry.origin_m);
    double r = std::sqrt(local.x * local.x + local.y * local.y);
    
    // ========== CYLINDRICAL GEOMETRY ==========
    if (dom.instrument != config::Instrument::Orbitrap) {
        return (local.z >= -DOMAIN_BOUNDARY_EPSILON && 
                local.z < dom.geometry.length_m) && 
               (r < dom.geometry.radius_m);
    }
    
    // ========== ORBITRAP: Hyperlogarithmic Electrodes ==========
    const double Rin = dom.geometry.radius_in_m;
    const double Rout = dom.geometry.radius_out_m;
    const double R_m = dom.geometry.radius_char_m;
    
    const double z_max = 0.5 * dom.geometry.length_m;
    const double z_abs = std::fabs(local.z);
    
    // DEBUG: First few calls
    static int orbitrap_debug_count = 0;
    bool do_debug = (orbitrap_debug_count < 3);
    if (do_debug) {
        std::cerr << "\n=== ORBITRAP is_inside_domain DEBUG ===\n";
        std::cerr << "  Global pos: (" << globalPos.x*1000 << ", " << globalPos.y*1000 << ", " << globalPos.z*1000 << ") mm\n";
        std::cerr << "  Local pos: (" << local.x*1000 << ", " << local.y*1000 << ", " << local.z*1000 << ") mm\n";
        std::cerr << "  r=" << r*1000 << " mm, z_abs=" << z_abs*1000 << " mm\n";
        std::cerr << "  Rin=" << Rin*1000 << " mm, Rout=" << Rout*1000 << " mm, R_m=" << R_m*1000 << " mm\n";
        std::cerr << "  z_max=" << z_max*1000 << " mm, length_m=" << dom.geometry.length_m*1000 << " mm\n";
        orbitrap_debug_count++;
    }
    
    // Check axial bounds
    if (z_abs > z_max + DOMAIN_BOUNDARY_EPSILON) {
        if (do_debug) std::cerr << "  FAIL: z_abs > z_max\n";
        return false;
    }
    
    // Compute allowed radial range at this z
    const double r_in_allowed = orbitrap_r_for_z(z_abs, Rin, R_m);
    const double r_out_allowed = orbitrap_r_for_z(z_abs, Rout, R_m);
    
    if (!(r_in_allowed > 0.0 && r_out_allowed > r_in_allowed)) {
        return false;
    }
    
    // Check if ion is in allowed corridor
    bool inside = (r >= r_in_allowed - DOMAIN_BOUNDARY_EPSILON) &&
                  (r <= r_out_allowed + DOMAIN_BOUNDARY_EPSILON);
    
    return inside;
}

double DomainManager::orbitrap_r_for_z(double z, double R, double R_m) const {
    const double z_abs = std::fabs(z);
    const double eps = 1e-10;
    const int max_iter = 80;
    
    // Special case: z=0 → r=R
    if (z_abs < eps) {
        return R;
    }
    
    // DEBUG
    static int bisect_debug_count = 0;
    bool do_debug = (bisect_debug_count < 2);
    
    // Initial bracket: At z>0, surface curves INWARD, so r < R
    // Start with bracket [0.1*R, R] to search for r where surface exists
    double r_lo = 0.1 * R;
    double r_hi = R;
    
    double f_lo = orbitrap_surface_residual(r_lo, z_abs, R, R_m);
    double f_hi = orbitrap_surface_residual(r_hi, z_abs, R, R_m);
    
    if (do_debug) {
        std::cerr << "\n=== orbitrap_r_for_z DEBUG ===\n";
        std::cerr << "  z=" << z_abs*1000 << " mm, R=" << R*1000 << " mm, R_m=" << R_m*1000 << " mm\n";
        std::cerr << "  Initial bracket: r_lo=" << r_lo*1000 << " mm, r_hi=" << r_hi*1000 << " mm\n";
        std::cerr << "  f_lo=" << f_lo << ", f_hi=" << f_hi << "\n";
        bisect_debug_count++;
    }
    
    // Expand bracket if needed
    int expand_iter = 0;
    while (f_lo * f_hi > 0.0 && expand_iter < 10) {
        r_lo *= 0.5;
        r_hi *= 1.5;
        f_lo = orbitrap_surface_residual(r_lo, z_abs, R, R_m);
        f_hi = orbitrap_surface_residual(r_hi, z_abs, R, R_m);
        expand_iter++;
        if (do_debug) {
            std::cerr << "  Expand iter " << expand_iter << ": r_lo=" << r_lo*1000 << " mm, r_hi=" << r_hi*1000 << " mm\n";
            std::cerr << "    f_lo=" << f_lo << ", f_hi=" << f_hi << "\n";
        }
    }
    
    if (do_debug) {
        std::cerr << "  Final bracket after expansion: r_lo=" << r_lo*1000 << " mm, r_hi=" << r_hi*1000 << " mm\n";
    }
    
    // Bisection
    double r_mid = R;
    for (int i = 0; i < max_iter; ++i) {
        r_mid = 0.5 * (r_lo + r_hi);
        double f_mid = orbitrap_surface_residual(r_mid, z_abs, R, R_m);
        
        if (std::fabs(f_mid) < eps) {
            break;
        }
        
        if (f_mid * f_lo > 0.0) {
            r_lo = r_mid;
            f_lo = f_mid;
        } else {
            r_hi = r_mid;
            f_hi = f_mid;
        }
    }
    
    if (do_debug) {
        std::cerr << "  Result: r_mid=" << r_mid*1000 << " mm\n";
    }
    
    return r_mid;
}

}  // namespace integrator
}  // namespace ICARION
