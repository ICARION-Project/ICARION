// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "EHSSCollisionHandler.h"
#include "core/physics/collisions/core/CollisionKernels.h"
#include "core/physics/collisions/core/VelocitySampling.h"
#include "core/physics/collisions/core/CollisionGeometry.h"
#include "utils/constants.h"
#include "core/log/Logger.h"
#include <cmath>
#include <stdexcept>

namespace {
    constexpr double MIN_RELATIVE_VELOCITY = 1e-10;  ///< Minimum relative velocity to consider collision [m/s]
}

namespace ICARION::physics {

EHSSCollisionHandler::EHSSCollisionHandler(
    GeometryMap geometry_map,
    bool enable_logging,
    const config::SpeciesDatabase* species_db
)
    : geometry_map_(std::move(geometry_map))  // Store copy (move to avoid unnecessary copying)
    , enable_logging_(enable_logging)
    , species_db_(species_db)
{
    if (geometry_map_.empty()) {
        throw std::invalid_argument("EHSSCollisionHandler: geometry_map cannot be empty!");
    }

    if (species_db_) {
        for (const auto& [id, props] : species_db_->species) {
            if (!props.ehss_samples_file) {
                continue;
            }
            EHSSOrientationSamples samples;
            std::string error_msg;
            if (!load_ehss_samples_file(*props.ehss_samples_file, samples, &error_msg)) {
                ICARION::log::Logger::get("collision")->warn(
                    "[EHSS] Failed to load orientation samples for '{}': {}",
                    id, error_msg);
                continue;
            }
            orientation_samples_.emplace(id, std::move(samples));
        }
        for (const auto& [id, props] : species_db_->species) {
            if (!props.ehss_offline_samples_file) {
                continue;
            }
            EHSSOfflineSampleSet samples;
            std::string error_msg;
            if (!load_ehss_offline_sample_set_file(*props.ehss_offline_samples_file, samples, &error_msg)) {
                ICARION::log::Logger::get("collision")->warn(
                    "[EHSS] Failed to load offline collision samples for '{}': {}",
                    id, error_msg);
                continue;
            }
            if (!samples.valid()) {
                ICARION::log::Logger::get("collision")->warn(
                    "[EHSS] Offline collision samples for '{}' are invalid (empty or malformed).",
                    id);
                continue;
            }
            offline_samples_.emplace(id, std::move(samples));
        }
    }
}

bool EHSSCollisionHandler::handle_collision(
    core::IonCollisionData& view,
    double dt,
    PhysicsRng& rng,
    const config::EnvironmentConfig& env
) {
    auto mark_once = [&](std::unordered_set<std::string>& set, const std::string& key) {
        std::lock_guard<std::mutex> lock(state_mutex_);
        return set.insert(key).second;
    };

    auto record_collision = [&]() {
        std::lock_guard<std::mutex> lock(state_mutex_);
        stats_.total_collisions++;
    };

    IonState ion;
    ion.vel = view.kin.vel();
    ion.mass_kg = view.kin.get_mass();
    ion.ion_charge_C = view.kin.get_charge();
    ion.CCS_m2 = view.get_CCS();
    ion.species_id = view.species_id();

    const EHSSOrientationSamples* samples = nullptr;
    size_t sample_idx = 0;
    bool use_samples = false;
    const EHSSOfflineSampleSet* offline = nullptr;
    size_t offline_idx = 0;
    bool use_offline = false;

    if (!offline_samples_.empty() && env.gas_mixture.empty()) {
        auto it_off = offline_samples_.find(ion.species_id);
        if (it_off != offline_samples_.end() && it_off->second.valid()) {
            if (!it_off->second.gas.empty() && it_off->second.gas != env.gas_species) {
                const std::string key = ion.species_id + ":" + env.gas_species;
                if (mark_once(warned_offline_gas_, key)) {
                    ICARION::log::Logger::get("collision")->warn(
                        "[EHSS] Offline samples for '{}' are tagged for gas '{}' (requested '{}'); "
                        "falling back to online EHSS.",
                        ion.species_id, it_off->second.gas, env.gas_species);
                }
            } else {
                const size_t n_orientations = it_off->second.sigma_eff_m2.size();
                double r = rng.uniform01() * static_cast<double>(n_orientations);
                offline_idx = static_cast<size_t>(r);
                if (offline_idx >= n_orientations) {
                    offline_idx = n_orientations - 1;
                }
                offline = &it_off->second;
                use_offline = true;
            }
        }
    }

    if (use_offline && offline) {
        if (mark_once(logged_offline_use_, ion.species_id)) {
            ICARION::log::Logger::get("collision")->info(
                "[EHSS] Using offline collision samples for '{}' (gas '{}', N={}, mu_stride={})",
                ion.species_id,
                offline->gas.empty() ? env.gas_species : offline->gas,
                offline->sigma_eff_m2.size(),
                offline->mu_stride
            );
        }
    }

    if (!use_offline && !orientation_samples_.empty()) {
        auto it = orientation_samples_.find(ion.species_id);
        if (it != orientation_samples_.end() && !it->second.orientations_quat.empty()) {
            samples = &it->second;
            const size_t n_orientations = samples->orientations_quat.size();
            bool ready = true;
            if (!env.gas_mixture.empty()) {
                for (const auto& comp : env.gas_mixture) {
                    auto itg = samples->areas_by_gas_m2.find(comp.species);
                    if (itg == samples->areas_by_gas_m2.end() || itg->second.size() != n_orientations) {
                        ready = false;
                        break;
                    }
                }
            } else {
                auto itg = samples->areas_by_gas_m2.find(env.gas_species);
                if (itg == samples->areas_by_gas_m2.end() || itg->second.size() != n_orientations) {
                    ready = false;
                }
            }

            if (ready) {
                double r = rng.uniform01() * static_cast<double>(n_orientations);
                sample_idx = static_cast<size_t>(r);
                if (sample_idx >= n_orientations) {
                    sample_idx = n_orientations - 1;
                }
                use_samples = true;
            }
        }
    }

    auto apply_collision = [&](double n,
                               bool check_probability,
                               double T_K,
                               double m_neutral,
                               const Vec3& v_gas,
                               double sigma_eff,
                               double gas_radius_m,
                               const EHSSOfflineSampleSet* offline_samples,
                               size_t offline_sample_idx,
                               const Vec3* sampled_neutral) -> bool {
        Vec3 v_neutral;
        if (sampled_neutral) {
            v_neutral = *sampled_neutral;
        } else {
            v_neutral = collision_core::VelocitySampling::sample_neutral_velocity(
                T_K, m_neutral, v_gas, rng
            );
        }

        const Vec3 v_rel = ion.vel - v_neutral;
        const double v_rel_mag = norm(v_rel);
        if (v_rel_mag < MIN_RELATIVE_VELOCITY) {
            return false;
        }

        if (check_probability) {
            const double P = 1.0 - std::exp(-n * sigma_eff * v_rel_mag * dt);
            if (rng.uniform01() >= P) {
                return false;
            }
        }

        Vec3 v_post;
        if (offline_samples && offline_samples->valid()) {
            const size_t mu_stride = offline_samples->mu_stride;
            size_t mu_idx = static_cast<size_t>(rng.uniform01() * static_cast<double>(mu_stride));
            if (mu_idx >= mu_stride) {
                mu_idx = mu_stride - 1;
            }
            double mu = offline_samples->mu_at(offline_sample_idx, mu_idx);
            mu = std::min(1.0, std::max(-1.0, mu));
            const double sin_theta = std::sqrt(std::max(0.0, 1.0 - mu * mu));
            const double phi = 2.0 * M_PI * rng.uniform01();

            const Vec3 collision_axis = v_rel / v_rel_mag;
            Vec3 t1, t2;
            collision_core::CollisionGeometry::construct_orthonormal_basis(collision_axis, t1, t2);

            const Vec3 n_contact = t1 * (sin_theta * std::cos(phi))
                                 + t2 * (sin_theta * std::sin(phi))
                                 + collision_axis * (-mu);

            const double vdotn = v_rel.x * n_contact.x + v_rel.y * n_contact.y + v_rel.z * n_contact.z;
            const Vec3 v_rel_reflected = v_rel - n_contact * (2.0 * vdotn);
            const double mt = ion.mass_kg + m_neutral;
            const double inv_mt = 1.0 / mt;
            const Vec3 Vcom = (ion.vel * ion.mass_kg + v_neutral * m_neutral) * inv_mt;
            v_post = Vcom + v_rel_reflected * (m_neutral * inv_mt);
        } else {
            auto it = geometry_map_.find(ion.species_id);
            if (it != geometry_map_.end() && !it->second.first.empty()) {
                const auto& [centers, radii] = it->second;
                // Use gas radius for contact geometry to match EHSS rate definition
                const double contact_radius = gas_radius_m;
                if (use_samples && samples) {
                    double Rori[3][3];
                    collision_core::CollisionGeometry::quaternion_to_rotation(
                        samples->orientations_quat[sample_idx], Rori
                    );
                    v_post = collision_core::CollisionKernels::ehss_collision_with_orientation(
                        ion.vel, v_neutral, ion.mass_kg, m_neutral, contact_radius,
                        centers, radii, Rori, rng, 256, sigma_eff, true
                    );
                } else {
                    v_post = collision_core::CollisionKernels::ehss_collision(
                        ion.vel, v_neutral, ion.mass_kg, m_neutral, contact_radius, centers, radii, rng
                    );
                }
            } else {
                throw std::runtime_error("[EHSSCollisionHandler] No geometry for species '" + ion.species_id + "'");
            }
        }

        ion.vel = v_post;
        view.kin.set_vel(ion.vel);
        record_collision();
        return true;
    };

    if (!env.gas_mixture.empty()) {
        // Compute rates with per-component sampled neutral velocities.
        std::vector<double> k_values;
        k_values.reserve(env.gas_mixture.size());
        std::vector<Vec3> sampled_neutrals;
        sampled_neutrals.reserve(env.gas_mixture.size());

        double k_total = 0.0;
        bool has_effective_component = false;
        for (const auto& comp : env.gas_mixture) {
            if (!comp.participates_in_collisions) {
                k_values.push_back(0.0);
                sampled_neutrals.push_back(Vec3{0.0, 0.0, 0.0});
                continue;
            }
            if (comp.density_m3 > 0.0) {
                has_effective_component = true;
            }

            const Vec3 v_neutral_i = collision_core::VelocitySampling::sample_neutral_velocity(
                env.temperature_K,
                comp.mass_kg,
                env.gas_velocity_m_s,
                rng
            );
            sampled_neutrals.push_back(v_neutral_i);

            double sigma_i = 0.0;
            // Note: offline-samples path is not active here (use_offline is always false
            // when gas_mixture is non-empty; the guard above ensures this).
            if (use_samples && samples) {
                sigma_i = samples->areas_by_gas_m2.at(comp.species)[sample_idx];
            } else {
                sigma_i = compute_effective_ccs(
                    ion,
                    comp.radius_m > 0.0 ? comp.radius_m : env.gas_radius_m,
                    comp.species
                );
            }

            const double n_i = comp.density_m3;
            const double v_rel_i = norm(ion.vel - v_neutral_i);
            if (sigma_i <= 0.0 || n_i <= 0.0 || v_rel_i < MIN_RELATIVE_VELOCITY) {
                k_values.push_back(0.0);
                continue;
            }

            const double k_i = n_i * sigma_i * v_rel_i;
            k_values.push_back(k_i);
            k_total += k_i;
        }

        if (k_total <= 0.0) {
            if (!has_effective_component) {
                static std::once_flag no_active_mixture_warn_once;
                std::call_once(no_active_mixture_warn_once, []() {
                    ICARION::log::Logger::main()->warn(
                        "[EHSS] gas_mixture is configured but has no active collision components; "
                        "collisions are disabled for this domain until mixture flags/densities are fixed.");
                });
            }
            return false;
        }

        double P_total = 1.0;
        if (k_total * dt <= 50.0) {
            P_total = 1.0 - std::exp(-k_total * dt);
        }
        if (rng.uniform01() >= P_total) {
            return false;
        }

        const double r = rng.uniform01() * k_total;
        size_t idx = 0;
        double cum = 0.0;
        for (; idx < k_values.size(); ++idx) {
            cum += k_values[idx];
            if (r < cum) {
                break;
            }
        }
        if (idx >= env.gas_mixture.size()) {
            idx = env.gas_mixture.size() - 1;
        }

        const auto& comp = env.gas_mixture[idx];
        double sigma_eff = 0.0;
        // Note: offline-samples path is not active here (use_offline is always false
        // when gas_mixture is non-empty; the guard above ensures this).
        if (use_samples && samples) {
            sigma_eff = samples->areas_by_gas_m2.at(comp.species)[sample_idx];
        } else {
            sigma_eff = compute_effective_ccs(
                ion,
                comp.radius_m > 0.0 ? comp.radius_m : env.gas_radius_m,
                comp.species
            );
        }
        if (sigma_eff <= 0.0) {
            return false;
        }

        const double gas_radius_m = (comp.radius_m > 0.0) ? comp.radius_m : env.gas_radius_m;
        return apply_collision(
            comp.density_m3,
            false,
            env.temperature_K,
            comp.mass_kg,
            env.gas_velocity_m_s,
            sigma_eff,
            gas_radius_m,
            nullptr,  // offline not used in mixture path
            offline_idx,
            &sampled_neutrals[idx]
        );
    }

    // Single gas path
    double sigma_eff = 0.0;
    if (use_offline && offline) {
        sigma_eff = offline->sigma_eff_m2[offline_idx];
    } else if (use_samples && samples) {
        sigma_eff = samples->areas_by_gas_m2.at(env.gas_species)[sample_idx];
    } else {
        sigma_eff = compute_effective_ccs(ion, env.gas_radius_m, env.gas_species);
    }
    if (sigma_eff <= 0.0) {
        return false;
    }

    return apply_collision(
        env.particle_density_m_3,
        true,
        env.temperature_K,
        env.gas_mass_kg,
        env.gas_velocity_m_s,
        sigma_eff,
        env.gas_radius_m,
        use_offline ? offline : nullptr,
        offline_idx,
        nullptr
    );
}

double EHSSCollisionHandler::compute_effective_ccs(
    const IonState& ion,
    double neutral_radius,
    const std::string& gas_id
) const {
    // 1. Try precomputed EHSS CCS map (BEST: most accurate)
    if (species_db_) {
        auto it = species_db_->species.find(ion.species_id);
        if (it != species_db_->species.end()) {
            const auto& map = it->second.ccs_ehss_m2;
            auto itg = map.find(gas_id);
            if (itg != map.end() && itg->second > 0.0) {
                return itg->second;  // Precomputed CCS (from ccs_precompute tool)
            }
        }
    }

    // 2. Try geometry-based computation (GOOD: accurate but slow)
    auto it_geom = geometry_map_.find(ion.species_id);
    if (it_geom != geometry_map_.end() && !it_geom->second.first.empty()) {
        const auto& [centers, radii] = it_geom->second;
        // Orientation-averaged projection approximation (OAPA)
        PhysicsRng rng(12345);
        double A_sum = 0.0;
        double Rmat[3][3];
        const int n_orientations = 64;
        for (int k = 0; k < n_orientations; ++k) {
            collision_core::CollisionGeometry::generate_random_rotation(rng, Rmat);
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
        return A_sum / static_cast<double>(n_orientations);  // OAPA average over orientations
    }

    // 3. Try automatic CCS derivation from reference (ACCEPTABLE: ~10% error)
    if (species_db_ && !gas_id.empty()) {
        auto it_spec = species_db_->species.find(ion.species_id);
        if (it_spec != species_db_->species.end()) {
            const auto& species = it_spec->second;
            
            if (species.CCS_m2 > 0.0 && species.ccs_reference_gas.has_value()) {
                double derived_ccs = derive_ccs_for_target_gas(
                    species.CCS_m2,
                    *species.ccs_reference_gas,
                    gas_id
                );
                
                if (derived_ccs > 0.0) {
                    // Log once per species:gas combination
                    std::string key = ion.species_id + ":" + gas_id;
                    bool should_log = false;
                    {
                        std::lock_guard<std::mutex> lock(state_mutex_);
                        should_log = warned_missing_sigma_.insert(key).second;
                    }
                    if (should_log) {
                        ICARION::log::Logger::get("collision")->info(
                            "[EHSS] Derived CCS for {}:{} = {:.2f} Å² (from {} reference). "
                            "For better accuracy: ccs_precompute --species {} --ref-gas {} --ref-ccs-A2 {:.2f}",
                            ion.species_id, gas_id, derived_ccs * 1e20, *species.ccs_reference_gas,
                            ion.species_id, *species.ccs_reference_gas, species.CCS_m2 * 1e20
                        );
                    }
                    return derived_ccs;  // Derived CCS (HSS approximation, ~10% error)
                }
            }
        }
    }

    // 4. Last resort: use reference CCS (WARNING: may be wrong gas!)
    if (ion.CCS_m2 > 0.0) {
        std::string key = ion.species_id + ":" + gas_id;
        bool should_log = false;
        {
            std::lock_guard<std::mutex> lock(state_mutex_);
            should_log = warned_missing_sigma_.insert(key).second;
        }
        if (should_log) {
            ICARION::log::Logger::get("collision")->warn(
                "[EHSS] Using reference CCS ({:.2f} Å²) for gas {} - may be inaccurate! "
                "Run: ccs_precompute --species {} --ref-gas <gas> --ref-ccs-A2 <ccs>",
                ion.CCS_m2 * 1e20, gas_id, ion.species_id
            );
        }
        return ion.CCS_m2;  // Last-resort fallback: reference CCS (may be wrong gas)
    }

    // 5. No CCS data available at all!
    throw std::runtime_error(
        "[EHSS] No CCS data for species '" + ion.species_id + "' and gas '" + gas_id + "'"
    );
}

CollisionStats EHSSCollisionHandler::get_stats() const {
    std::lock_guard<std::mutex> lock(state_mutex_);
    return stats_;
}

void EHSSCollisionHandler::reset_stats() {
    std::lock_guard<std::mutex> lock(state_mutex_);
    stats_ = {};
}

double EHSSCollisionHandler::derive_ccs_for_target_gas(
    double sigma_ref_m2,
    const std::string& gas_ref,
    const std::string& gas_target
) const {
    // Gas radii lookup (from utils/constants.h)
    static const std::unordered_map<std::string, double> GAS_RADII = {
        {"He", RADIUS_HE_M},
        {"N2", RADIUS_N2_M},
        {"O2", RADIUS_O2_M},
        {"Ar", RADIUS_AR_M},
        {"CO2", RADIUS_CO2_M},
        {"Ne", RADIUS_NE_M},
        {"H2O", RADIUS_H2O_M}
    };
    
    auto it_ref = GAS_RADII.find(gas_ref);
    auto it_tgt = GAS_RADII.find(gas_target);
    
    if (it_ref == GAS_RADII.end() || it_tgt == GAS_RADII.end()) {
        return 0.0;  // Unknown gas
    }
    
    // Extract ion radius from reference CCS
    // sigma_ref = π(r_ion + r_ref)² → r_ion = sqrt(sigma_ref/π) - r_ref
    double r_ref = it_ref->second;
    double r_ion = std::max(0.0, std::sqrt(sigma_ref_m2 / M_PI) - r_ref);
    
    // Compute CCS for target gas
    // sigma_target = π(r_ion + r_target)²
    double r_target = it_tgt->second;
    return M_PI * (r_ion + r_target) * (r_ion + r_target);
}

} // namespace ICARION::physics
