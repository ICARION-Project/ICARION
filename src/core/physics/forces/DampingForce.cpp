// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "DampingForce.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/config/types/SpeciesConfig.h"
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
    (void)t;  // Time-independent (deterministic damping)

    if (model_ == DampingModel::None) {
        return Vec3{0.0, 0.0, 0.0};
    }

    IonState ion;
    ion.pos = ensemble.get_pos(ion_idx);
    ion.vel = ensemble.get_vel(ion_idx);
    ion.mass_kg = ensemble.mass_data()[ion_idx];
    ion.ion_charge_C = ensemble.charge_data()[ion_idx];
    ion.CCS_m2 = ensemble.CCS(ion_idx);
    ion.reduced_mobility_cm2_Vs = ensemble.mobility(ion_idx);
    ion.species_id = ensemble.species_id(ion_idx);
    ion.active = ensemble.active_data()[ion_idx] != 0;
    ion.born = ensemble.born_data()[ion_idx] != 0;
    ion.current_domain_index = ensemble.domain_index(ion_idx);
    ion.birth_time_s = ensemble.birth_time(ion_idx);

    double gamma = calculate_gamma(ion, ctx);

    if (gamma <= 0.0) {
        return Vec3{0.0, 0.0, 0.0};
    }

    return ion.vel * (-gamma * ion.mass_kg);
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
    (void)ctx;  // SSOT: Read from env_ config directly
    
    // SSOT: Read gas properties from config
    const double gas_density = env_->particle_density_m_3;
    const double v_th = env_->mean_thermal_velocity_m_s;
    const double m_neutral = env_->gas_mass_kg;
    
    // ========================================================================
    // Model-specific damping coefficient calculation
    // ========================================================================
    
    switch (model_) {
        case DampingModel::HardSphere: {
            // Hard-sphere collision frequency:
            // γ = ν_collision = n·σ·v_th
            // (Legacy had mass factor, but this contradicts Mason-Schamp mobility theory)
            
            const double CCS = ion.CCS_m2;
            const double m_ion = ion.mass_kg;
            
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
            
            const double alpha = env_->gas_polarizability_m3;
            
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
            if (!env_->gas_mixture.empty()) {
                const double sigma_ref = ion.CCS_m2;
                if (sigma_ref <= 0.0) {
                    throw std::runtime_error("[DampingForce] Missing CCS for species '" + ion.species_id + "'");
                }
                if (env_->gas_mass_kg <= 0.0) {
                    throw std::runtime_error("[DampingForce] Invalid mixture-averaged gas mass");
                }
                double frac_sum = 0.0;
                double inv_K_sum = 0.0;
                for (const auto& comp : env_->gas_mixture) {
                    if (!comp.participates_in_collisions) continue;
                    frac_sum += comp.mole_fraction;
                }
                if (frac_sum <= 0.0) {
                    throw std::runtime_error("[DampingForce] Gas mixture has zero participating fraction");
                }
                for (const auto& comp : env_->gas_mixture) {
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
                    double scale_mass = std::sqrt(env_->gas_mass_kg / comp.mass_kg);
                    double Ki = K0_m2_Vs * scale_sigma * scale_mass * LOSCHMIDT_CONSTANT / gas_density;
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
            const double ion_mobility = K0_m2_Vs * LOSCHMIDT_CONSTANT / gas_density;
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
