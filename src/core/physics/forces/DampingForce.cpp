// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "DampingForce.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/config/types/SpeciesConfig.h"
#include "core/config/types/DomainConfig.h"
#include "core/log/Logger.h"
#include "utils/constants.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace ICARION {
namespace physics {

// ============================================================================
// Constructor
// ============================================================================

DampingForce::DampingForce(
    const ICARION::config::EnvironmentConfig& env,
    DampingModel model,
    const ICARION::config::SpeciesDatabase* species_db
)
    : env_(&env)
    , model_(model)
    , species_db_(species_db)
{
    // No validation - all params can be zero (disabled force)
}

// ============================================================================
// IForce Interface Implementation
// ============================================================================

Vec3 DampingForce::compute(const core::IonEnsemble& ensemble, size_t ion_idx, double t,
                           const ForceContext& ctx) const {
    ForceState state{};
    state.pos = ensemble.get_pos(ion_idx);
    state.vel = ensemble.get_vel(ion_idx);
    state.mass_kg = ensemble.mass_data()[ion_idx];
    state.ion_charge_C = ensemble.charge_data()[ion_idx];
    state.CCS_m2 = ensemble.CCS(ion_idx);
    state.reduced_mobility_cm2_Vs = ensemble.mobility(ion_idx);
    state.species_id = ensemble.species_id(ion_idx);
    state.active = ensemble.active_data()[ion_idx] != 0;
    state.born = ensemble.born_data()[ion_idx] != 0;
    state.current_domain_index = ensemble.domain_index(ion_idx);
    state.birth_time_s = ensemble.birth_time(ion_idx);
    state.ensemble_index = ion_idx;

    return compute_soa(state, t, ctx);
}

Vec3 DampingForce::compute_soa(const ForceState& state, double t,
                               const ForceContext& ctx) const {
    (void)t;  // Time-independent (deterministic damping)

    if (model_ == DampingModel::None) {
        return Vec3{0.0, 0.0, 0.0};
    }

    IonState ion = state.to_ion_state();

    double gamma = calculate_gamma(ion, ctx);

    if (gamma <= 0.0) {
        return Vec3{0.0, 0.0, 0.0};
    }

    const auto* active_env = env_;
    const config::GeometryConfig* geometry = nullptr;
    if (ctx.domain) {
        active_env = &ctx.domain->environment;
        geometry = &ctx.domain->geometry;
    }
    const Vec3 gas_velocity = active_env->gas_velocity_at(state.pos, geometry);
    const Vec3 relative_velocity = ion.vel - gas_velocity;
    return relative_velocity * (-gamma * ion.mass_kg);
}

std::string DampingForce::name() const {
    switch (model_) {
        case DampingModel::HardSphere: return "Damping(HardSphere)";
        case DampingModel::Langevin:   return "Damping(Langevin)";
        case DampingModel::Friction:   return "Damping(Friction)";
        default:                       return "Damping(None)";
    }
}

// ============================================================================
// Damping Coefficient Calculation
// ============================================================================

double DampingForce::calculate_gamma(const IonState& ion, const ForceContext& ctx) const {
    const auto* active_env = env_;
    if (ctx.domain) {
        active_env = &ctx.domain->environment;
    }

    // Prefer per-ion environment cache when available (keeps damping in sync with
    // time-varying environment updates from SimulationEngine).
    double gas_density = active_env->particle_density_m_3;
    double temperature_K = active_env->temperature_K;
    double m_neutral = active_env->gas_mass_kg;
    if (ctx.ion_ensemble && ctx.ion_index < ctx.ion_ensemble->size()) {
        gas_density = ctx.ion_ensemble->gas_density(ctx.ion_index);
        temperature_K = ctx.ion_ensemble->temperature(ctx.ion_index);
        m_neutral = ctx.ion_ensemble->neutral_mass(ctx.ion_index);
    }
    const double v_th = (temperature_K > 0.0 && m_neutral > 0.0)
        ? std::sqrt(8.0 * BOLTZMANN_CONSTANT * temperature_K / (M_PI * m_neutral))
        : 0.0;
    
    // ========================================================================
    // Model-specific damping coefficient calculation
    // ========================================================================
    
    switch (model_) {
        case DampingModel::HardSphere: {
            // Hard-sphere collision frequency:
            // γ = ν_collision = n·σ·v_th
            // (Legacy had mass factor, but this contradicts Mason-Schamp mobility theory)
            
            const double CCS = ion.CCS_m2;
            
            if (CCS <= 0.0 || gas_density <= 0.0 || v_th <= 0.0 || m_neutral <= 0.0) {
                return 0.0;  // Missing parameters
            }
            
            // Collision frequency (without mass factor - matches mobility theory)
            return gas_density * CCS * v_th;
        }
        
        case DampingModel::Langevin: {
            // Langevin collision frequency (velocity-dependent):
            // γ = ν_Langevin = n·σ_L(v)·v_th·m_reduced/m_ion
            // where σ_L = π·q·√(α/(4πε₀·m_reduced))/|v|
            
            const double alpha = active_env->gas_polarizability_m3;
            
            const double m_ion = ion.mass_kg;
            const double q = ion.ion_charge_C;
            
            if (alpha <= 0.0 || gas_density <= 0.0 || v_th <= 0.0 || m_neutral <= 0.0) {
                return 0.0;  // Missing parameters
            }
            
            // Ion velocity magnitude (guard against division by zero)
            const double v_mag = std::max(1e-6, norm(ion.vel));
            
            // Reduced mass
            const double m_reduced = (m_ion * m_neutral) / (m_ion + m_neutral);
            
            // Langevin cross-section (velocity-dependent)
            // σ_L = π·q·√(α/(4πε₀·m_reduced))/v
            const double epsilon_factor = 4.0 * M_PI * EPSILON_0;
            const double cs = M_PI * q * std::sqrt(alpha / (epsilon_factor * m_reduced)) / v_mag;
            
            // Collision frequency
            return gas_density * cs * v_th * m_neutral / (m_neutral + m_ion);
        }
        
        case DampingModel::Friction: {
            // Mobility-based friction:
            // γ = q/(K·m_ion) where K = K₀·(n₀/n) is actual mobility
            const double K0_cm2_Vs = ion.reduced_mobility_cm2_Vs;
            const double m_ion = ion.mass_kg;
            const double q = ion.ion_charge_C;
            if (K0_cm2_Vs <= 0.0 || m_ion <= 0.0 || gas_density <= 0.0) {
                return 0.0;  // Missing parameters
            }
            const double K0_m2_Vs = K0_cm2_Vs * 1e-4;
            const double N_ratio = LOSCHMIDT_CONSTANT / gas_density;
            // Temperature correction (Mason-Schamp): K ∝ sqrt(T0 / T)
            const double teff_scale = std::sqrt(STP_TEMP / std::max(temperature_K, 1.0));
            const double K0_teff_m2_Vs = K0_m2_Vs * teff_scale;

            // DEBUG: Print what we're actually using
            static bool debug_once = false;
            if (!debug_once) {
                std::cerr << "\n[DampingForce::Friction DEBUG]\n";
                std::cerr << "  K₀ (from ion) = " << K0_cm2_Vs << " cm²/(V·s)\n";
                std::cerr << "  gas_density = " << gas_density << " m⁻³\n";
                std::cerr << "  LOSCHMIDT = " << LOSCHMIDT_CONSTANT << " m⁻³\n";
                debug_once = true;
            }

            auto lookup_sigma = [&](const config::GasMixtureComponent& comp) -> double {
                // Prefer DB CCS_HSS per gas if available
                if (species_db_) {
                    auto it_spec = species_db_->species.find(ion.species_id);
                    if (it_spec != species_db_->species.end()) {
                        const auto& map = it_spec->second.ccs_hss_m2;
                        auto it_g = map.find(comp.species);
                        if (it_g != map.end() && it_g->second > 0.0) {
                            return it_g->second;
                        }
                    }
                }
                // Config override per gas
                if (comp.cross_section_m2 > 0.0) {
                    return comp.cross_section_m2;
                }
                // Fallback to ion CCS (warn once to avoid silent fallback)
                if (species_db_) {
                    std::string key = ion.species_id + ":" + comp.species;
                    if (!warned_missing_sigma_.count(key)) {
                        ICARION::log::debug_log(
                            "[DampingForce] Warning: No CCS_HSS for gas '" + comp.species +
                            "' and species '" + ion.species_id + "'; using ion.CCS_m2 fallback");
                        warned_missing_sigma_.insert(key);
                    }
                }
                return ion.CCS_m2;
            };

            // If mixture present: use Blanc's law with per-gas mobility scaling
            if (!active_env->gas_mixture.empty()) {
                const double sigma_ref = ion.CCS_m2;
                if (sigma_ref <= 0.0) {
                    throw std::runtime_error("[DampingForce] Missing CCS for species '" + ion.species_id + "'");
                }
                if (active_env->gas_mass_kg <= 0.0) {
                    throw std::runtime_error("[DampingForce] Invalid mixture-averaged gas mass");
                }
                double frac_sum = 0.0;
                double inv_K_sum = 0.0;
                for (const auto& comp : active_env->gas_mixture) {
                    if (!comp.participates_in_collisions) continue;
                    frac_sum += comp.mole_fraction;
                }
                if (frac_sum <= 0.0) {
                    throw std::runtime_error("[DampingForce] Gas mixture has zero participating fraction");
                }
                for (const auto& comp : active_env->gas_mixture) {
                    if (!comp.participates_in_collisions) continue;
                    double f = comp.mole_fraction / frac_sum;
                    double sigma_i = lookup_sigma(comp);
                    if (sigma_i <= 0.0) {
                        throw std::runtime_error("[DampingForce] Missing CCS for species '" + ion.species_id +
                                                 "' in gas '" + comp.species + "'");
                    }
                    if (comp.mass_kg <= 0.0) {
                        throw std::runtime_error("[DampingForce] Invalid mass for gas '" + comp.species + "'");
                    }

                    // Heuristic scaling: mobility ∝ 1/σ · sqrt(m_ref / m_gas)
                    double scale_sigma = sigma_ref / sigma_i;
                    double scale_mass = std::sqrt(active_env->gas_mass_kg / comp.mass_kg);
                    double Ki = K0_teff_m2_Vs * scale_sigma * scale_mass * N_ratio;
                    if (Ki <= 0.0) {
                        continue;
                    }
                    inv_K_sum += f / Ki;
                }
                if (inv_K_sum > 0.0) {
                    double K_mix = 1.0 / inv_K_sum;
                    return q / (K_mix * m_ion);
                }
                throw std::runtime_error("[DampingForce] No valid mobility contributions for gas mixture");
            }

            // Single-gas fallback: K = K₀·(n₀/n)
            const double ion_mobility = K0_teff_m2_Vs * N_ratio;
            const double gamma_result = q / (ion_mobility * m_ion);
            
            // DEBUG: Print calculated values
            static bool debug_calc = false;
            if (!debug_calc) {
                std::cerr << "  K (actual) = " << ion_mobility << " m²/(V·s)\n";
                std::cerr << "  γ = " << gamma_result << " Hz\n";
                std::cerr << "  Expected v_term = q·E/(γ·m) = " << (q * 10000.0 / (gamma_result * m_ion)) << " m/s @ E=10kV/m\n\n";
                debug_calc = true;
            }
            
            return gamma_result;
        }
        
        default:
            return 0.0;
    }
}

} // namespace physics
} // namespace ICARION
