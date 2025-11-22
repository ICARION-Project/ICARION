// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "ReactionConfig.h"
#include "utils/constants.h"
#include <cmath>
#include <stdexcept>

namespace ICARION::config {

double Reaction::compute_rate_constant(double temperature_K) const {
    switch (rate_model) {
        case RateModel::Constant:
            // k(T) = k₀ (no temperature dependence)
            return rate_constant_m3s;
        
        case RateModel::Arrhenius: {
            // k(T) = A × exp(-Eₐ / (kB·T))
            
            if (activation_energy_eV <= 0.0) {
                // No barrier → constant rate
                return rate_constant_m3s;
            }
            
            // Convert Eₐ from eV to Joules
            const double Ea_J = activation_energy_eV * ELEM_CHARGE_C;
            
            // Compute exponential term: exp(-Eₐ / (kB·T))
            const double exponent = -Ea_J / (BOLTZMANN_CONSTANT * temperature_K);
            
            // Numerical safety: if exponent > 50, exp(50) overflows
            // if exponent < -50, exp(-50) ≈ 0
            if (exponent < -50.0) {
                return 0.0;  // Rate effectively zero (huge barrier)
            }
            if (exponent > 50.0) {
                // This shouldn't happen (negative barrier?), but be safe
                return rate_constant_m3s * std::exp(50.0);
            }
            
            return rate_constant_m3s * std::exp(exponent);
        }
        
        case RateModel::ModifiedArrhenius: {
            // k(T) = A × (T/T₀)ⁿ × exp(-Eₐ / (kB·T))
            
            // Temperature power term: (T/T₀)ⁿ
            double T_ratio = temperature_K / reference_temperature_K;
            double T_factor = std::pow(T_ratio, temperature_exponent);
            
            // Arrhenius exponential term: exp(-Eₐ / (kB·T))
            double exp_factor = 1.0;
            if (activation_energy_eV > 0.0) {
                const double Ea_J = activation_energy_eV * ELEM_CHARGE_C;
                const double exponent = -Ea_J / (BOLTZMANN_CONSTANT * temperature_K);
                
                // Numerical safety (same as Arrhenius)
                if (exponent < -50.0) {
                    exp_factor = 0.0;
                } else if (exponent > 50.0) {
                    exp_factor = std::exp(50.0);
                } else {
                    exp_factor = std::exp(exponent);
                }
            }
            
            return rate_constant_m3s * T_factor * exp_factor;
        }
        
        default:
            throw std::runtime_error(
                "Reaction '" + id + "': unknown rate model (internal error)"
            );
    }
}

} // namespace ICARION::config
