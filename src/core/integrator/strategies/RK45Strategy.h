// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "IIntegrationStrategy.h"
#include "core/types/IonState.h"
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
 * **FSAL Property:** First-Same-As-Last: k7(n) = k1(n+1) → 6 evaluations/step
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
        double absolute_min_step_s = 0.0; ///< Optional absolute dt floor (0 = disabled)
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
     * @brief Advance ion by one timestep using RK45 (SoA interface)
     */
    void step(
        core::IonEnsemble& ensemble,
        size_t ion_idx,
        double t,
        double dt,
        const physics::ForceRegistry& force_registry
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
        double sum_step_used = 0.0;  ///< Accumulated dt of accepted steps
    };
    
    const StepStats& get_stats() const { return stats_; }
    void reset_stats() { stats_ = StepStats{}; }
    
    /**
     * @brief Get adaptive configuration (for GPU integration)
     */
    const AdaptiveConfig& get_config() const { return config_; }
    double last_dt_used() const { return last_dt_used_; }
    double last_dt_suggested() const { return last_dt_suggested_; }
    void enable_stats(bool enabled) { stats_enabled_ = enabled; }

private:
    struct RK45State {
        double last_error = 1.0;
        bool fsal_available = false;
        double last_dt_used = 0.0;
        double last_dt_suggested = 0.0;
        double k1_ax = 0.0, k1_ay = 0.0, k1_az = 0.0;
    };

    std::vector<RK45State> per_ion_state_;

    AdaptiveConfig config_;
    StepStats stats_;
    double last_dt_used_ = 0.0;      ///< Actual dt used in last step (per call)
    double last_dt_suggested_ = 0.0; ///< Suggested dt for next step (per call)
    bool stats_enabled_ = true;
    
    // FSAL storage: k1 for next step = k7 from previous step
    RK45State& state_for(size_t ion_idx);
    
    /**
     * @brief Compute acceleration at given state
     * 
     * Wrapper for ForceRegistry to match RK interface.
     */
    void compute_acceleration_batch(
        double& ax, double& ay, double& az,
        const core::IonEnsemble& ensemble,
        size_t ion_idx,
        double t,
        const physics::ForceRegistry& force_registry
    );

    void compute_acceleration_state(
        double& ax, double& ay, double& az,
        const IonState& state,
        double t,
        const physics::ForceRegistry& force_registry
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
