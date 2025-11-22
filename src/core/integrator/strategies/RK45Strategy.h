// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   ------------------------------------------------
 *   Modular framework for simulating ion trajectories in custom
 *   electric fields and background gas environments.
 *
 *   @file       RK45Strategy.h
 *   @brief      Adaptive 4/5 Runge-Kutta integration (Dormand-Prince)
 *
 *   @details
 *   Embedded Runge-Kutta method with automatic step size control.
 *   SSOT-compliant implementation using ForceRegistry.
 *
 *   **Algorithm (Dormand-Prince):**
 *   - 7 stages with shared evaluations
 *   - 4th-order solution for propagation
 *   - 5th-order solution for error estimation
 *   - Error = ||y5 - y4||
 *   - Step control: dt_new = dt * safety * (tol/error)^(1/5)
 *
 *   **Properties:**
 *   - Order: 4(5) embedded pair
 *   - Timestep: Adaptive (error-based)
 *   - Stability: Good for smooth, non-stiff systems
 *   - Cost: 6 force evaluations per step (FSAL property)
 *   - Overhead: ~50% more than RK4 for same accuracy
 *   - Advantage: Automatic error control, optimal efficiency
 *
 *   **When to Use:**
 *   - Variable dynamics (fast/slow phases)
 *   - Tight error tolerance requirements
 *   - Unknown optimal timestep a priori
 *   - Long simulations where efficiency matters
 *
 *   **When NOT to Use:**
 *   - Stiff systems (eigenvalue spread > 10^6)
 *   - Discontinuous forces (collisions, field boundaries)
 *   - Fixed timestep required (data output synchronization)
 *
 *   @date       2025-11-22
 *   @version    1.0.0
 *   @authors    ICARION Development Team
 *
 * =====================================================================
 */
#pragma once

#include "IIntegrationStrategy.h"
#include <limits>

namespace ICARION {
namespace integrator {

/**
 * @brief Adaptive RK4/5 (Dormand-Prince) integration strategy
 * 
 * Embedded Runge-Kutta method with automatic step size adjustment
 * based on local truncation error estimation.
 * 
 * **Dormand-Prince Coefficients:**
 * ```
 *    0   |
 *   1/5  | 1/5
 *   3/10 | 3/40      9/40
 *   4/5  | 44/45    -56/15     32/9
 *   8/9  | 19372/6561  -25360/2187  64448/6561  -212/729
 *    1   | 9017/3168   -355/33   46732/5247   49/176  -5103/18656
 *    1   | 35/384       0        500/1113     125/192  -2187/6784   11/84
 * ```
 * 
 * **Error Estimation:**
 * - e_i = |y5_i - y4_i| for each component i
 * - error = max_i(e_i / (atol + rtol * |y_i|))
 * - Target: error ≈ 1.0
 * 
 * **Step Size Control (PI Controller):**
 * - Proportional: k_p = 0.7 / 5 = 0.14
 * - Integral: k_i = 0.4 / 5 = 0.08
 * - dt_new = dt * safety * (1/error)^k_p * (error_prev/error)^k_i
 * - Safety factor: 0.9 (conservative)
 * - Min step: dt_initial * 1e-6
 * - Max step: dt_initial * 1e3
 * 
 * **FSAL Property:**
 * First-Same-As-Last: k7(n) = k1(n+1) → 6 evaluations/step
 * 
 * **Thread Safety:**
 * - Each ion needs independent error tracking
 * - Cannot share state between ions
 * - Parallelize with per-ion RK45Strategy instances
 */
class RK45Strategy : public IIntegrationStrategy {
public:
    /**
     * @brief Configuration for adaptive timestep control
     */
    struct AdaptiveConfig {
        double atol = 1e-8;           ///< Absolute error tolerance [m, m/s]
        double rtol = 1e-6;           ///< Relative error tolerance [dimensionless]
        double safety_factor = 0.9;   ///< Step size safety factor (0.8-0.95)
        double min_step_factor = 1e-6; ///< Min step = dt_initial * min_step_factor
        double max_step_factor = 1e3;  ///< Max step = dt_initial * max_step_factor
        double max_step_increase = 5.0; ///< Limit growth per step
        double max_step_decrease = 0.1; ///< Limit shrinkage per step
    };

    /**
     * @brief Default constructor with standard tolerances
     */
    RK45Strategy();
    
    /**
     * @brief Constructor with custom adaptive parameters
     * 
     * @param config Adaptive timestep configuration
     * 
     * **Default Tolerances:**
     * - atol = 1e-8 m, m/s (ion position precision)
     * - rtol = 1e-6 (0.0001% relative error)
     * 
     * **Typical Use Cases:**
     * - Tight: atol=1e-10, rtol=1e-8 (high precision)
     * - Standard: atol=1e-8, rtol=1e-6 (default)
     * - Loose: atol=1e-6, rtol=1e-4 (fast exploration)
     */
    explicit RK45Strategy(const AdaptiveConfig& config);
    
    /**
     * @brief Advance ion by one timestep using RK45 (base interface)
     * 
     * @param ion Ion state (updated in-place)
     * @param t Current time [s]
     * @param dt Timestep [s] (fixed from config)
     * @param force_registry Force computation engine
     * @param domain Domain configuration (SSOT!)
     * @param all_ions All ions (for space charge)
     * 
     * This overrides the base interface with fixed dt.
     * For adaptive timestep control, use step_adaptive().
     */
    void step(
        IonState& ion,
        double t,
        double dt,
        const physics::ForceRegistry& force_registry,
        const config::DomainConfig& domain,
        const std::vector<IonState>& all_ions
    ) override;
    
    /**
     * @brief Advance ion by one adaptive timestep using RK45
     * 
     * @param ion Ion state (updated in-place)
     * @param t Current time [s]
     * @param dt_inout Timestep [s] (input: suggested, output: actual used)
     * @param force_registry Force computation engine
     * @param domain Domain configuration (SSOT!)
     * @param all_ions All ions (for space charge)
     * 
     * **Algorithm:**
     * 1. Compute 7 RK stages (k1...k7) using Dormand-Prince coefficients
     * 2. Calculate 4th-order solution: y4 = y + Σ(b4_i * k_i)
     * 3. Calculate 5th-order solution: y5 = y + Σ(b5_i * k_i)
     * 4. Estimate error: e = ||y5 - y4|| (scaled by tolerance)
     * 5. Accept step if e <= 1.0, otherwise reject and retry
     * 6. Update step size: dt_new based on PI controller
     * 7. Use FSAL: store k7 for next step's k1
     * 
     * **Error Handling:**
     * - If error > 1: reject step, reduce dt, retry (up to 10 attempts)
     * - If error < 0.1: increase dt for next step
     * - If dt < dt_min: accept step with warning
     * 
     * **Step Size Update:**
     * - dt_new = dt * S * (1/e)^(1/5) where S = safety_factor
     * - Clamped to [dt*max_decrease, dt*max_increase]
     * - Clamped to [dt_min, dt_max]
     * 
     * **Output:**
     * - ion.pos, ion.vel: Updated to 4th-order solution
     * - dt_inout: Set to actual timestep used (may differ from input)
     * 
     * **Performance:**
     * - Successful step: 6 force evaluations (FSAL)
     * - Rejected step: 6 evaluations (wasted)
     * - Typical acceptance rate: 95-99%
     * 
     * **Thread Safety:**
     * - Must use separate RK45Strategy instance per ion
     * - Internal state: last_error_, k1_stored_ (not thread-safe)
     */
    void step_adaptive(
        IonState& ion,
        double t,
        double& dt_inout,
        const physics::ForceRegistry& force_registry,
        const config::DomainConfig& domain,
        const std::vector<IonState>& all_ions
    );
    
    /**
     * @brief Get method name
     */
    std::string name() const override { return "RK45"; }
    
    /**
     * @brief Check if adaptive
     */
    bool is_adaptive() const override { return true; }
    
    /**
     * @brief Get statistics for adaptive stepping
     */
    struct StepStats {
        size_t accepted_steps = 0;
        size_t rejected_steps = 0;
        double min_step_used = std::numeric_limits<double>::max();
        double max_step_used = 0.0;
        double avg_error = 0.0;
    };
    
    const StepStats& get_stats() const { return stats_; }
    void reset_stats() { stats_ = StepStats{}; }

private:
    AdaptiveConfig config_;
    StepStats stats_;
    double last_error_ = 1.0;  ///< Previous error for PI controller
    bool fsal_available_ = false; ///< Flag for FSAL k1 reuse
    
    // FSAL storage: k1 for next step = k7 from previous step
    struct {
        double ax = 0.0, ay = 0.0, az = 0.0;  ///< Acceleration [m/s²]
    } k1_stored_;
    
    /**
     * @brief Compute acceleration at given state
     * 
     * Wrapper for ForceRegistry to match RK interface.
     */
    void compute_acceleration(
        double& ax, double& ay, double& az,
        const IonState& ion,
        double t,
        const physics::ForceRegistry& force_registry,
        const config::DomainConfig& domain,
        const std::vector<IonState>& all_ions
    );
    
    /**
     * @brief Estimate local truncation error
     * 
     * @param y4 4th-order solution
     * @param y5 5th-order solution
     * @param y_current Current state (for scaling)
     * @return Scaled error (target: ~1.0)
     * 
     * Error metric: max_i |e_i| / (atol + rtol * |y_i|)
     */
    double estimate_error(
        const IonState& y4,
        const IonState& y5,
        const IonState& y_current
    ) const;
    
    /**
     * @brief Compute new step size using PI controller
     * 
     * @param current_dt Current timestep
     * @param error Current scaled error
     * @param dt_min Minimum allowed step
     * @param dt_max Maximum allowed step
     * @return New timestep
     */
    double compute_new_step(
        double current_dt,
        double error,
        double dt_min,
        double dt_max
    );
};

} // namespace integrator
} // namespace ICARION
