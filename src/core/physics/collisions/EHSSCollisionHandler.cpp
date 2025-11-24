// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "EHSSCollisionHandler.h"
#include "collisionHelpers.h"
#include "utils/constants.h"
#include "core/log/Logger.h"
#include <cmath>
#include <stdexcept>

namespace {
    constexpr double MIN_RELATIVE_VELOCITY = 1e-10;  ///< Minimum relative velocity to consider collision [m/s]
}

namespace ICARION::physics {

EHSSCollisionHandler::EHSSCollisionHandler(
    const GeometryMap& geometry_map,
    bool enable_logging,
    const config::SpeciesDatabase* species_db
)
    : geometry_map_(geometry_map)  // Store reference (no copy!)
    , enable_logging_(enable_logging)
    , species_db_(species_db)
{
    if (geometry_map_.empty()) {
        throw std::invalid_argument("EHSSCollisionHandler: geometry_map cannot be empty!");
    }
}

bool EHSSCollisionHandler::handle_collision(
    IonState& ion,
    double dt,
    EhssRng& rng,
    const config::EnvironmentConfig& env
) {
    // Mixture-aware handling
    auto do_collision = [&](double n, double T_K, double m_neutral, const Vec3& v_gas, double neutral_radius) -> bool {
        const double sigma_eff = compute_effective_ccs(ion, neutral_radius, "");

        EHSSParams p;
        p.n = n;
        p.dt = dt;
        p.mi = ion.mass_kg;
        p.mn = m_neutral;
        p.kB = BOLTZMANN_CONSTANT;
        p.Tn = T_K;
        p.ubx = v_gas.x;
        p.uby = v_gas.y;
        p.ubz = v_gas.z;
        p.Rn = neutral_radius;
        p.sigma_eff = sigma_eff;

        const Vec3 v_neutral = sample_neutral_velocity(p, rng);
        const Vec3 v_rel = ion.vel - v_neutral;
        const double v_rel_mag = norm(v_rel);
        if (v_rel_mag < MIN_RELATIVE_VELOCITY) {
            return false;
        }

        const double P = 1.0 - std::exp(-n * sigma_eff * v_rel_mag * dt);
        if (rng.uniform01() >= P) {
            return false;
        }

        auto it = geometry_map_.find(ion.species_id);
        Vec3 v_post;
        if (it != geometry_map_.end() && !it->second.first.empty()) {
            const auto& [centers, radii] = it->second;
            v_post = collide_ehss_cpu_geometry_given_neutral(
                ion.vel, v_neutral, p, centers, radii, rng
            );
        } else {
            throw std::runtime_error("[EHSSCollisionHandler] No geometry for species '" + ion.species_id + "'");
        }

        ion.vel = v_post;
        stats_.total_collisions++;
        return true;
    };

    if (!env.gas_mixture.empty()) {
        // Compute rates per component and select gas
        Vec3 v_rel_bulk = ion.vel - env.gas_velocity_m_s;
        double v_rel_mag = norm(v_rel_bulk);
        if (v_rel_mag < MIN_RELATIVE_VELOCITY) {
            return false;
        }
        std::vector<double> k_values;
        k_values.reserve(env.gas_mixture.size());
        double k_total = 0.0;
        for (const auto& comp : env.gas_mixture) {
            double sigma_i = compute_effective_ccs(ion, comp.radius_m > 0.0 ? comp.radius_m : env.gas_radius_m, comp.species);
            double n_i = comp.density_m3;
            if (sigma_i <= 0.0 || n_i <= 0.0) {
                k_values.push_back(0.0);
                continue;
            }
            double k_i = n_i * sigma_i * v_rel_mag;
            k_values.push_back(k_i);
            k_total += k_i;
        }
        if (k_total <= 0.0) {
            return false;
        }
        double P_total = 1.0;
        if (k_total * dt <= 50.0) {
            P_total = 1.0 - std::exp(-k_total * dt);
        }
        if (rng.uniform01() >= P_total) {
            return false;
        }
        double r = rng.uniform01() * k_total;
        size_t idx = 0;
        double cum = 0.0;
        for (; idx < k_values.size(); ++idx) {
            cum += k_values[idx];
            if (r < cum) break;
        }
        if (idx >= env.gas_mixture.size()) {
            idx = env.gas_mixture.size() - 1;
        }
        const auto& comp = env.gas_mixture[idx];
        return do_collision(comp.density_m3, env.temperature_K, comp.mass_kg, env.gas_velocity_m_s,
                            comp.radius_m > 0.0 ? comp.radius_m : env.gas_radius_m);
    }

    // Single gas path
    return do_collision(
        env.particle_density_m_3,
        env.temperature_K,
        env.gas_mass_kg,
        env.gas_velocity_m_s,
        env.gas_radius_m
    );
}

double EHSSCollisionHandler::compute_effective_ccs(
    const IonState& ion,
    double neutral_radius,
    const std::string& gas_id
) const {
    // Prefer precomputed EHSS CCS map if available
    if (species_db_) {
        auto it = species_db_->species.find(ion.species_id);
        if (it != species_db_->species.end()) {
            const auto& map = it->second.ccs_ehss_m2;
            auto itg = map.find(gas_id);
            if (itg != map.end() && itg->second > 0.0) {
                return itg->second;
            } else if (!gas_id.empty() && ion.CCS_m2 > 0.0) {
                std::string key = ion.species_id + ":" + gas_id;
                if (!warned_missing_sigma_.count(key)) {
                    ICARION::log::debug_log(
                        "[EHSSCollisionHandler] Warning: No CCS_EHSS for gas '" + gas_id +
                        "' and species '" + ion.species_id + "'; using geometry/CCS fallback");
                    warned_missing_sigma_.insert(key);
                }
            }
        }
    }

    // Try to find geometry for this species
    auto it = geometry_map_.find(ion.species_id);
    
    if (it != geometry_map_.end() && !it->second.first.empty()) {
        const auto& [centers, radii] = it->second;
        // Simple orientation-averaged projection (same approach as ccs_precompute)
        EhssRng rng(12345);
        double A_sum = 0.0;
        double Rmat[3][3];
        const int n_orientations = 64;
        for (int k = 0; k < n_orientations; ++k) {
            rand_rotation(rng, Rmat);
            double Rmax = 0.0;
            for (size_t i = 0; i < radii.size(); ++i) {
                Vec3 c{
                    Rmat[0][0] * centers[i].x + Rmat[0][1] * centers[i].y + Rmat[0][2] * centers[i].z,
                    Rmat[1][0] * centers[i].x + Rmat[1][1] * centers[i].y + Rmat[1][2] * centers[i].z,
                    Rmat[2][0] * centers[i].x + Rmat[2][1] * centers[i].y + Rmat[2][2] * centers[i].z
                };
                double r_eff = std::sqrt(c.x * c.x + c.y * c.y) + radii[i] + neutral_radius;
                Rmax = std::max(Rmax, r_eff);
            }
            A_sum += M_PI * Rmax * Rmax;
        }
        return A_sum / static_cast<double>(n_orientations);
    } else {
        // No geometry available - use stored CCS from ion state
        return ion.CCS_m2;
    }
}

} // namespace ICARION::physics
